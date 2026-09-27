// fp16acc.cu: separate the input rounding from the accumulation in the FP16 route (report 2.1 (a)).
// The same FP16-rounded inputs go through DGEMM (reference), CUDA-core SGEMM (pedantic, round-to-nearest sums),
// the tensor-core FP16 product (inputs 16F, CUBLAS_COMPUTE_32F, as realhgemm) and the tensor-core TF32 product
// (FP16 values are exact in TF32).  a = sum(C C_ref)/sum(C_ref^2) - 1 is the uniform shrink.
// usage: fp16acc m n k posonly(0|1) [seed]
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cublas_v2.h>
#define CK(x) do { auto e_ = (x); if ((int)e_ != 0) { printf("error %d at %s:%d\n", (int)e_, __FILE__, __LINE__); exit(1); } } while (0)

static double pow2scale(double amax) {            // power of 2 that puts amax in [2^13, 2^14)
  int e; std::frexp(amax, &e);                     // amax = f 2^e, f in [0.5,1)
  return std::ldexp(1.0, 14 - e);
}

int main(int argc, char** argv) {
  int m = atoi(argv[1]), n = atoi(argv[2]), k = atoi(argv[3]), pos = atoi(argv[4]);
  unsigned seed = argc > 5 ? atoi(argv[5]) : 1;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  // values over about 3 decades (log-uniform magnitude), positive or random sign
  auto val = [&]() { double v = std::pow(10.0, -3.0*u(rng)); return (pos || u(rng) < 0.5) ? v : -v; };
  std::vector<double> A(size_t(k)*m), B(size_t(k)*n);   // A: k x m (used as A^T), B: k x n, column major
  for (auto& x : A) x = val();
  for (auto& x : B) x = val();
  double amax = 0, bmax = 0;
  for (auto x : A) amax = std::fmax(amax, std::fabs(x));
  for (auto x : B) bmax = std::fmax(bmax, std::fabs(x));
  double fa = pow2scale(amax), fb = pow2scale(bmax);
  std::vector<__half> A16(A.size()), B16(B.size());
  std::vector<float> Af(A.size()), Bf(B.size());
  std::vector<double> Ad(A.size()), Bd(B.size());
  size_t nsubA = 0, nsubB = 0, nzA = 0, nzB = 0;
  for (size_t i = 0; i < A.size(); i++) {
    A16[i] = __double2half(A[i]*fa); float f = __half2float(A16[i]); Af[i] = f; Ad[i] = f;
    if (f == 0) nzA++; else if (std::fabs(f) < 6.103515625e-05) nsubA++;
  }
  for (size_t i = 0; i < B.size(); i++) {
    B16[i] = __double2half(B[i]*fb); float f = __half2float(B16[i]); Bf[i] = f; Bd[i] = f;
    if (f == 0) nzB++; else if (std::fabs(f) < 6.103515625e-05) nsubB++;
  }
  cublasHandle_t h; CK(cublasCreate(&h));
  __half *dA16, *dB16; float *dAf, *dBf, *dC; double *dAd, *dBd, *dCd;
  CK(cudaMalloc(&dA16, A.size()*2)); CK(cudaMalloc(&dB16, B.size()*2));
  CK(cudaMalloc(&dAf, A.size()*4)); CK(cudaMalloc(&dBf, B.size()*4));
  CK(cudaMalloc(&dAd, A.size()*8)); CK(cudaMalloc(&dBd, B.size()*8));
  CK(cudaMalloc(&dC, size_t(m)*n*4)); CK(cudaMalloc(&dCd, size_t(m)*n*8));
  CK(cudaMemcpy(dA16, A16.data(), A.size()*2, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dB16, B16.data(), B.size()*2, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dAf, Af.data(), A.size()*4, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dBf, Bf.data(), B.size()*4, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dAd, Ad.data(), A.size()*8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dBd, Bd.data(), B.size()*8, cudaMemcpyHostToDevice));
  std::vector<double> Cref(size_t(m)*n), C(size_t(m)*n);
  std::vector<float> Cf(size_t(m)*n);
  double one = 1, zero = 0; float onef = 1, zerof = 0;
  // reference: DGEMM of the same rounded inputs (C = A^T B)
  CK(cublasDgemm(h, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k, &one, dAd, k, dBd, k, &zero, dCd, m));
  CK(cudaMemcpy(Cref.data(), dCd, Cref.size()*8, cudaMemcpyDeviceToHost));
  auto report = [&](const char* name) {
    double s1 = 0, s2 = 0, maxrel = 0, se = 0;
    for (size_t i = 0; i < Cref.size(); i++) { s1 += (double)Cf[i]*Cref[i]; s2 += Cref[i]*Cref[i]; }
    double a = s1/s2 - 1;
    double cmax = 0; for (auto x : Cref) cmax = std::fmax(cmax, std::fabs(x));
    for (size_t i = 0; i < Cref.size(); i++) {
      double d = (double)Cf[i] - Cref[i]; se += d*d;
      maxrel = std::fmax(maxrel, std::fabs(d)/cmax);
    }
    printf("  %-22s a = %+.3e   rms/max = %.2e   max/max = %.2e\n", name, a, std::sqrt(se/Cref.size())/cmax, maxrel);
  };
  // CUDA-core SGEMM (no tensor cores)
  CK(cublasSetMathMode(h, CUBLAS_PEDANTIC_MATH));
  CK(cublasSgemm(h, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k, &onef, dAf, k, dBf, k, &zerof, dC, m));
  CK(cudaMemcpy(Cf.data(), dC, Cf.size()*4, cudaMemcpyDeviceToHost)); report("SGEMM (CUDA cores)");
  // tensor cores, FP16 inputs, FP32 accumulation (the realhgemm call)
  CK(cublasSetMathMode(h, CUBLAS_DEFAULT_MATH));
  CK(cublasGemmEx(h, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k, &onef, dA16, CUDA_R_16F, k, dB16, CUDA_R_16F, k, &zerof,
                  dC, CUDA_R_32F, m, CUBLAS_COMPUTE_32F, CUBLAS_GEMM_DEFAULT));
  CK(cudaMemcpy(Cf.data(), dC, Cf.size()*4, cudaMemcpyDeviceToHost)); report("tensor cores FP16");
  // tensor cores, TF32 (FP16 values are exact in TF32, so again only the accumulation differs)
  CK(cublasGemmEx(h, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k, &onef, dAf, CUDA_R_32F, k, dBf, CUDA_R_32F, k, &zerof,
                  dC, CUDA_R_32F, m, CUBLAS_COMPUTE_32F_FAST_TF32, CUBLAS_GEMM_DEFAULT));
  CK(cudaMemcpy(Cf.data(), dC, Cf.size()*4, cudaMemcpyDeviceToHost)); report("tensor cores TF32");
  printf("  (k=%d, %s values over 3 decades; k*2^-24 = %.2e; FP16 subnormal/zero: A %.1e/%.1e, B %.1e/%.1e)\n", k,
         pos ? "positive" : "random-sign", k*std::ldexp(1.0,-24),
         double(nsubA)/A.size(), double(nzA)/A.size(), double(nsubB)/B.size(), double(nzB)/B.size());
  return 0;
}
