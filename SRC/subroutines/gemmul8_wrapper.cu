#include "gemmul8.hpp"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <type_traits>
#include <unordered_map>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

// One workspace for all calls, grown on demand.  A cudaMalloc/cudaFree pair per GEMM
// synchronizes the device each time; Sigma_c issues thousands of GEMMs per q point.
static void*  g_work      = nullptr;
static size_t g_work_size = 0;

static void* gemmul8_workspace(size_t need) {
    if (need <= g_work_size) return g_work;
    if (g_work != nullptr) cudaFree(g_work);   // cudaFree waits for the GEMMs that used it
    g_work = nullptr;
    g_work_size = 0;
    cudaError_t err = cudaMalloc(&g_work, need);
    if (err != cudaSuccess || g_work == nullptr) {
        printf("gemmul8 workspace: cudaMalloc failed: %s (requested %zu bytes)\n",
               cudaGetErrorString(err), need);
        g_work = nullptr;
        return nullptr;
    }
    g_work_size = need;
    return g_work;
}

// Kept scaled/split A (GEMMul8 skip_scalA) per caller key.  The caller promises that one key means one
// matrix A until gemmul8_cache_reset_.  A kept A is reused only with the same type, op(A), m, k, moduli and
// fast mode (a key may be used with both rows of a table, which can differ in the moduli); otherwise the call
// runs without keeping A.  Budget: ECALJ_LA_CACHE_GB (default 4) and at most a quarter of the free device memory.
// (Bug fixed 2026-09-27 22:00: only the key, type and op were compared, so a key used with rows of other moduli
// could reuse an A split for a different number of moduli.)
struct KeptA { void* ptr; size_t size; size_t m, k; unsigned int moduli; int fast; };
static std::unordered_map<long long, KeptA> g_kept;
static size_t g_kept_bytes = 0;
static size_t g_kept_cap   = 0;
static bool   g_kept_init  = false;

// Keeping A is OFF unless ECALJ_GEMMUL8_KEEP=1 (Bug found 2026-09-29 07:35, not yet understood): with the kept A the fp64
// Sigma_c of LiTi2O4 6^3 was eV off (SEc at iteration 1 from LDA: max 2.79 eV, rms 1.50 eV against fp32), while GEMMul8
// without keeping A (ECALJ_LA_CACHE_GB=0) and cuBLAS agree with fp32 within 1 meV (kt1, research log 2026-09-29 07:40).
// The conditions of GEMMul8 for skip_scalA (same A, op, m, k, moduli, fast mode, INT8) look met, so the cause is still
// open.  The reuse saved little: iteration 1 from LDA took 1076 s with it and 1109 s without (fp64, kt1, 2 GPUs).
// realsgemm/realhgemm keep their A' as before (fp32/tf32 are right).
static size_t kept_cap() {
    if (!g_kept_init) {
        g_kept_init = true;
        const char* keep = std::getenv("ECALJ_GEMMUL8_KEEP");
        if (keep == nullptr || std::atoi(keep) != 1) { g_kept_cap = 0; return g_kept_cap; }
        double gb = 4.0;
        if (const char* s = std::getenv("ECALJ_LA_CACHE_GB")) gb = std::atof(s);
        size_t free_b = 0, total_b = 0;
        cudaMemGetInfo(&free_b, &total_b);
        g_kept_cap = (size_t)std::min(gb * 1e9, 0.25 * (double)free_b);
    }
    return g_kept_cap;
}

template <typename T> constexpr int typecode() {
    return std::is_same<T, double>::value ? 1 : std::is_same<T, cuDoubleComplex>::value ? 2 :
           std::is_same<T, cuComplex>::value ? 3 : 4;
}

template <typename T>
void gemmul8_gemm_impl(
    void* handle_ptr,
    int transa,
    int transb,
    int m_int,
    int n_int,
    int k_int,
    T alpha,
    const void* devA,
    int lda_int,
    const void* devB,
    int ldb_int,
    T beta,
    void* devC,
    int ldc_int,
    unsigned int num_moduli,
    int fastmode,
    int key
) {
    size_t m   = m_int;
    size_t n   = n_int;
    size_t k   = k_int;
    size_t lda = lda_int;
    size_t ldb = ldb_int;
    size_t ldc = ldc_int;

    cublasHandle_t cublas_handle = *reinterpret_cast<cublasHandle_t*>(handle_ptr);
    cublasOperation_t op_a = static_cast<cublasOperation_t>(transa);
    cublasOperation_t op_b = static_cast<cublasOperation_t>(transb);

    constexpr bool is_complex =
        std::is_same<T, cuComplex>::value ||
        std::is_same<T, cuDoubleComplex>::value;

    // Keep A when a key is given and the budget allows; reuse it when the key was seen.
    void* keptA  = nullptr;
    bool  reuseA = false;
    size_t worksizeA, worksizeB;
    if (key >= 0) {
        gemmul8::workSize<is_complex, gemmul8::Backend::INT8>(m, n, k, num_moduli, true, false,
                                                               &worksizeA, &worksizeB, fastmode != 0);
        long long kk = ((long long)key << 4) | (long long)(typecode<T>() * 4 + transa);
        auto it = g_kept.find(kk);
        if (it != g_kept.end()) {
            const KeptA& e = it->second;
            if (e.size >= worksizeA && e.m == m && e.k == k && e.moduli == num_moduli && e.fast == fastmode) {
                keptA  = e.ptr;
                reuseA = true;
            }
        } else if (g_kept_bytes + worksizeA <= kept_cap()) {
            if (cudaMalloc(&keptA, worksizeA) == cudaSuccess) {
                g_kept[kk] = KeptA{keptA, worksizeA, m, k, num_moduli, fastmode};
                g_kept_bytes += worksizeA;
            } else {
                keptA = nullptr;
            }
        }
    }
    const bool enable_skip_A = (keptA != nullptr);

    const size_t worksize = gemmul8::workSize<is_complex, gemmul8::Backend::INT8>(
        m, n, k, num_moduli, enable_skip_A, false, &worksizeA, &worksizeB, fastmode != 0);

    // Shared workspace: B and the rest, plus A when it is not kept.
    const size_t shared = enable_skip_A ? worksize - worksizeA : worksize;
    void* work_total = gemmul8_workspace(shared);
    if (work_total == nullptr) {
        printf("gemmul8 workSize inputs: m=%zu n=%zu k=%zu num_moduli=%u\n",
               m, n, k, num_moduli);
        return;
    }
    void* p_workA;
    void* p_workB;
    if (enable_skip_A) {
        p_workA = keptA;
        p_workB = work_total;
    } else {
        p_workA = work_total;
        p_workB = static_cast<char*>(work_total) + worksizeA;
    }
    void* p_work_rem = static_cast<char*>(p_workB) + worksizeB;

    gemmul8::gemm<T, gemmul8::Backend::INT8>(
        cublas_handle, op_a, op_b, m, n, k,
        &alpha,
        static_cast<const T*>(devA), lda,
        static_cast<const T*>(devB), ldb,
        &beta,
        static_cast<T*>(devC), ldc,
        num_moduli,
        (fastmode != 0),
        p_work_rem,
        p_workA,
        p_workB,
        enable_skip_A,
        false,
        reuseA,
        false
    );
}

extern "C" {

void gemmul8_dgemm_(void* handle_ptr, int transa, int transb, int m_int, int n_int, int k_int,
                    double alpha, const void* devA, int lda_int, const void* devB, int ldb_int,
                    double beta, void* devC, int ldc_int, unsigned int num_moduli, int fastmode, int key) {
    gemmul8_gemm_impl<double>(handle_ptr, transa, transb, m_int, n_int, k_int, alpha, devA, lda_int,
                              devB, ldb_int, beta, devC, ldc_int, num_moduli, fastmode, key);
}

void gemmul8_zgemm_(void* handle_ptr, int transa, int transb, int m_int, int n_int, int k_int,
                    cuDoubleComplex alpha, const void* devA, int lda_int, const void* devB, int ldb_int,
                    cuDoubleComplex beta, void* devC, int ldc_int, unsigned int num_moduli, int fastmode, int key) {
    gemmul8_gemm_impl<cuDoubleComplex>(handle_ptr, transa, transb, m_int, n_int, k_int, alpha, devA, lda_int,
                                       devB, ldb_int, beta, devC, ldc_int, num_moduli, fastmode, key);
}

void gemmul8_cgemm_(void* handle_ptr, int transa, int transb, int m_int, int n_int, int k_int,
                    cuComplex alpha, const void* devA, int lda_int, const void* devB, int ldb_int,
                    cuComplex beta, void* devC, int ldc_int, unsigned int num_moduli, int fastmode, int key) {
    gemmul8_gemm_impl<cuComplex>(handle_ptr, transa, transb, m_int, n_int, k_int, alpha, devA, lda_int,
                                 devB, ldb_int, beta, devC, ldc_int, num_moduli, fastmode, key);
}

// Forget the kept A (call when the matrices behind the keys change, e.g. at the next q point).
void gemmul8_cache_reset_() {
    if (g_kept.empty()) return;
    cudaDeviceSynchronize();
    for (auto& e : g_kept) cudaFree(e.second.ptr);
    g_kept.clear();
    g_kept_bytes = 0;
}

} // extern "C"

extern "C" void gemmul8_init_handle_(void** handle_out) {
    auto h = new cublasHandle_t;
    cublasCreate(h);
    *handle_out = h;
}

// GEMMul8 launches its kernels and INT8 GEMMs on the stream of its handle; m_blas sets it to the stream of
// its own cuBLAS handle before each call (cublas_set_stream: Sigma_c runs on OpenACC queue 1).
extern "C" void gemmul8_set_stream_(void* handle_ptr, cudaStream_t stream) {
    cublasSetStream(*reinterpret_cast<cublasHandle_t*>(handle_ptr), stream);
}

extern "C" void gemmul8_finalize_handle_(void* handle_ptr) {
    auto h = reinterpret_cast<cublasHandle_t*>(handle_ptr);
    cublasDestroy(*h);
    delete h;
    gemmul8_cache_reset_();
    if (g_work != nullptr) cudaFree(g_work);
    g_work = nullptr;
    g_work_size = 0;
}
