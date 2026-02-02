#include "gemmul8.hpp"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cuComplex.h>

extern "C" {
void gemmul8_dgemm_(
    void* handle_ptr,
    int transa,
    int transb,
    int m_int,
    int n_int,
    int k_int,
    double alpha,
    const void* devA,
    int lda_int,
    const void* devB,
    int ldb_int,
    double beta,
    void* devC,
    int ldc_int,
    unsigned int num_moduli,
    int fastmode,
    int enable_skip_A,
    int enable_skip_B,
    int skip_scalA,
    int skip_scalB
) {
    size_t m = m_int;
    size_t n = n_int;
    size_t k = k_int;
    size_t lda = lda_int;
    size_t ldb = ldb_int;
    size_t ldc = ldc_int;

    cublasHandle_t cublas_handle = *reinterpret_cast<cublasHandle_t*>(handle_ptr);
    cublasOperation_t op_a = static_cast<cublasOperation_t>(transa);
    cublasOperation_t op_b = static_cast<cublasOperation_t>(transb);

    bool enable_skip_A_bool = (enable_skip_A != 0);
    bool enable_skip_B_bool = (enable_skip_B != 0);

    size_t worksizeA, worksizeB;
    const size_t worksize = gemmul8::workSize<true, false>( // UseExtraWorkspace = false
        m, n, k, num_moduli,
        enable_skip_A_bool,
        enable_skip_B_bool,
        &worksizeA,
        &worksizeB
    );

    void* work_total = nullptr;
    cudaMalloc(&work_total, worksize);

    void* p_workA = work_total;
    void* p_workB = static_cast<char*>(p_workA) + worksizeA;
    void* p_work_rem = static_cast<char*>(p_workB) + worksizeB;

    gemmul8::gemm<double, false>(cublas_handle, op_a, op_b, m, n, k, // UseExtraWorkspace = false
        &alpha,
        static_cast<const double*>(devA), lda,
        static_cast<const double*>(devB), ldb, 
        &beta,
        static_cast<double*>(devC), ldc, 
        num_moduli, 
        (fastmode != 0), 
        p_work_rem,
        p_workA,
        p_workB,
        enable_skip_A_bool,
        enable_skip_B_bool,
        (skip_scalA != 0), 
        (skip_scalB != 0)
    );
    cudaFree(work_total);
}
}

extern "C" {
void gemmul8_zgemm_(
    void* handle_ptr,
    int transa,
    int transb,
    int m_int,
    int n_int,
    int k_int,
    cuDoubleComplex alpha,
    const void* devA,
    int lda_int,
    const void* devB,
    int ldb_int,
    cuDoubleComplex beta,
    void* devC,
    int ldc_int,
    unsigned int num_moduli,
    int fastmode,
    int enable_skip_A,
    int enable_skip_B,
    int skip_scalA,
    int skip_scalB
) {
    size_t m = m_int;
    size_t n = n_int;
    size_t k = k_int;
    size_t lda = lda_int;
    size_t ldb = ldb_int;
    size_t ldc = ldc_int;

    cublasHandle_t cublas_handle = *reinterpret_cast<cublasHandle_t*>(handle_ptr);
    cublasOperation_t op_a = static_cast<cublasOperation_t>(transa);
    cublasOperation_t op_b = static_cast<cublasOperation_t>(transb);

    bool enable_skip_A_bool = (enable_skip_A != 0);
    bool enable_skip_B_bool = (enable_skip_B != 0);

    size_t worksizeA, worksizeB;
    const size_t worksize = gemmul8::workSize<true, false>( // UseExtraWorkspace = false
        m, n, k, num_moduli,
        enable_skip_A_bool,
        enable_skip_B_bool,
        &worksizeA,
        &worksizeB
    );

    void* work_total = nullptr;
    cudaMalloc(&work_total, worksize);

    void* p_workA = work_total;
    void* p_workB = static_cast<char*>(p_workA) + worksizeA;
    void* p_work_rem = static_cast<char*>(p_workB) + worksizeB;

    gemmul8::gemm<cuDoubleComplex, false>(cublas_handle, op_a, op_b, m, n, k, // UseExtraWorkspace = false
        &alpha,
        static_cast<const cuDoubleComplex*>(devA), lda,
        static_cast<const cuDoubleComplex*>(devB), ldb, 
        &beta,
        static_cast<cuDoubleComplex*>(devC), ldc, 
        num_moduli, 
        (fastmode != 0), 
        p_work_rem,
        p_workA,
        p_workB,
        enable_skip_A_bool,
        enable_skip_B_bool,
        (skip_scalA != 0), 
        (skip_scalB != 0)
    );
    cudaFree(work_total);
}
} // extern "C"
