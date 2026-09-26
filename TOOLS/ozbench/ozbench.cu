// ozbench: speed and accuracy of complex GEMM on one GPU --
//   cuBLAS FP32 (CUBLAS_COMPUTE_32F), cuBLAS TF32, cuBLAS FP64 (zgemm), and GEMMul8 (Ozaki-II, INT8)
//   emulating cgemm / zgemm with a given number of moduli.
// The shapes are those of the hgw kernels (C = op(A) B with opA = C):
//   Sigma_c imaginary axis: m=k=ngb, n = (#middle states)*(#bands of Sigma)
//   Sigma_c real axis     : m=k=ngb, n = nttp (pairs at one frequency; small)
//   final contraction     : m=n=ntq, k = ngb*(#middle states)
// Usage: ozbench m n k [reps] [rowscale_decades]
//   rowscale_decades > 0 multiplies row i of A and column j of B by 10^(u*d), u in [-1/2,1/2]
//   (a crude stand-in for the dynamic range of W and zmel).
// Error = max_ij |C - C_ref| / max_ij |C_ref|, C_ref from FP64 zgemm of the same (FP32-rounded) inputs.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuComplex.h>
#include "gemmul8.hpp"

#define CK(x) do{cudaError_t e=(x); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__); exit(1);} }while(0)
#define CB(x) do{cublasStatus_t s=(x); if(s!=CUBLAS_STATUS_SUCCESS){printf("cuBLAS %d @%d\n",(int)s,__LINE__); exit(1);} }while(0)

static double now(){ return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

template<class F> double timeit(F f, int reps){
  f(); CK(cudaDeviceSynchronize());
  double t0=now();
  for(int r=0;r<reps;r++) f();
  CK(cudaDeviceSynchronize());
  return (now()-t0)/reps;
}

static double relerr_c(const std::vector<cuComplex>& c, const std::vector<cuDoubleComplex>& r){
  double emax=0, rmax=0;
  for(size_t i=0;i<r.size();i++){
    double dx=c[i].x-r[i].x, dy=c[i].y-r[i].y;
    emax=fmax(emax, sqrt(dx*dx+dy*dy)); rmax=fmax(rmax, sqrt(r[i].x*r[i].x+r[i].y*r[i].y));
  }
  return emax/rmax;
}
static double relerr_z(const std::vector<cuDoubleComplex>& c, const std::vector<cuDoubleComplex>& r){
  double emax=0, rmax=0;
  for(size_t i=0;i<r.size();i++){
    double dx=c[i].x-r[i].x, dy=c[i].y-r[i].y;
    emax=fmax(emax, sqrt(dx*dx+dy*dy)); rmax=fmax(rmax, sqrt(r[i].x*r[i].x+r[i].y*r[i].y));
  }
  return emax/rmax;
}

int main(int argc, char** argv){
  if(argc<4){ printf("usage: ozbench m n k [reps] [rowscale_decades]\n"); return 1; }
  size_t m=atol(argv[1]), n=atol(argv[2]), k=atol(argv[3]);
  int reps = argc>4 ? atoi(argv[4]) : 5;
  double dec = argc>5 ? atof(argv[5]) : 0.0;
  const cublasOperation_t opA=CUBLAS_OP_C, opB=CUBLAS_OP_N;   // A is k x m (stored), op(A) = A^H is m x k
  size_t lda=k, ldb=k, ldc=m;
  double flop = 8.0*m*n*k;
  printf("m n k = %zu %zu %zu   flop/GEMM = %.3g   reps=%d  row/col scale decades=%.1f\n", m,n,k,flop,reps,dec);

  std::mt19937_64 rng(12345); std::normal_distribution<double> g(0,1); std::uniform_real_distribution<double> u(-0.5,0.5);
  std::vector<double> sa(m), sb(n);
  for(auto& s: sa) s=pow(10.0, u(rng)*dec);
  for(auto& s: sb) s=pow(10.0, u(rng)*dec);
  std::vector<cuComplex> hA(k*m), hB(k*n);
  for(size_t j=0;j<m;j++) for(size_t i=0;i<k;i++) hA[i+j*lda]=make_cuComplex(g(rng)*sa[j], g(rng)*sa[j]);   // column j of A = row j of A^H
  for(size_t j=0;j<n;j++) for(size_t i=0;i<k;i++) hB[i+j*ldb]=make_cuComplex(g(rng)*sb[j], g(rng)*sb[j]);
  std::vector<cuDoubleComplex> hAz(k*m), hBz(k*n);                 // the same (FP32-rounded) numbers in FP64
  for(size_t i=0;i<hA.size();i++) hAz[i]=make_cuDoubleComplex(hA[i].x,hA[i].y);
  for(size_t i=0;i<hB.size();i++) hBz[i]=make_cuDoubleComplex(hB[i].x,hB[i].y);

  cuComplex *dA,*dB,*dC; cuDoubleComplex *dAz,*dBz,*dCz;
  CK(cudaMalloc(&dA,sizeof(cuComplex)*k*m)); CK(cudaMalloc(&dB,sizeof(cuComplex)*k*n)); CK(cudaMalloc(&dC,sizeof(cuComplex)*m*n));
  CK(cudaMalloc(&dAz,sizeof(cuDoubleComplex)*k*m)); CK(cudaMalloc(&dBz,sizeof(cuDoubleComplex)*k*n)); CK(cudaMalloc(&dCz,sizeof(cuDoubleComplex)*m*n));
  CK(cudaMemcpy(dA,hA.data(),sizeof(cuComplex)*k*m,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dB,hB.data(),sizeof(cuComplex)*k*n,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dAz,hAz.data(),sizeof(cuDoubleComplex)*k*m,cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dBz,hBz.data(),sizeof(cuDoubleComplex)*k*n,cudaMemcpyHostToDevice));
  cublasHandle_t h; CB(cublasCreate(&h));
  cuComplex a1=make_cuComplex(1,0), b0=make_cuComplex(0,0);
  cuDoubleComplex a1z=make_cuDoubleComplex(1,0), b0z=make_cuDoubleComplex(0,0);

  // reference: FP64 zgemm
  auto zg=[&]{ CB(cublasZgemm(h,opA,opB,m,n,k,&a1z,dAz,lda,dBz,ldb,&b0z,dCz,ldc)); };
  double tz=timeit(zg, reps>2?2:reps);
  std::vector<cuDoubleComplex> ref(m*n); CK(cudaMemcpy(ref.data(),dCz,sizeof(cuDoubleComplex)*m*n,cudaMemcpyDeviceToHost));
  printf("%-28s %9.3f ms %8.2f TFLOPS  err --\n","cuBLAS zgemm (FP64)", tz*1e3, flop/tz/1e12);

  std::vector<cuComplex> hc(m*n); std::vector<cuDoubleComplex> hcz(m*n);
  auto run_c=[&](const char* name, cublasComputeType_t ct){
    auto f=[&]{ CB(cublasGemmEx(h,opA,opB,m,n,k,&a1,dA,CUDA_C_32F,lda,dB,CUDA_C_32F,ldb,&b0,dC,CUDA_C_32F,ldc,ct,CUBLAS_GEMM_DEFAULT)); };
    double t=timeit(f,reps);
    CK(cudaMemcpy(hc.data(),dC,sizeof(cuComplex)*m*n,cudaMemcpyDeviceToHost));
    printf("%-28s %9.3f ms %8.2f TFLOPS  err %.2e\n",name,t*1e3,flop/t/1e12,relerr_c(hc,ref));
  };
  run_c("cuBLAS cgemm FP32", CUBLAS_COMPUTE_32F);
  run_c("cuBLAS cgemm TF32", CUBLAS_COMPUTE_32F_FAST_TF32);
  {
    auto f=[&]{ CB(cublasCgemm3m(h,opA,opB,m,n,k,&a1,dA,lda,dB,ldb,&b0,dC,ldc)); };
    double t=timeit(f,reps);
    CK(cudaMemcpy(hc.data(),dC,sizeof(cuComplex)*m*n,cudaMemcpyDeviceToHost));
    printf("%-28s %9.3f ms %8.2f TFLOPS  err %.2e\n","cuBLAS cgemm3m",t*1e3,flop/t/1e12,relerr_c(hc,ref));
  }

  // cuBLAS built-in emulation (CUDA 13): FP64 via fixed point (Ozaki-type, INT8 tensor cores),
  // FP32 via 3xBF16 (9 BF16 products)
  {
    auto f=[&]{ CB(cublasGemmEx(h,opA,opB,m,n,k,&a1z,dAz,CUDA_C_64F,lda,dBz,CUDA_C_64F,ldb,&b0z,dCz,CUDA_C_64F,ldc,
                                CUBLAS_COMPUTE_64F_EMULATED_FIXEDPOINT,CUBLAS_GEMM_DEFAULT)); };
    for(int strat: {1,2}){   // PERFORMANT, EAGER
      CB(cublasSetEmulationStrategy(h,(cublasEmulationStrategy_t)strat));
      double t=timeit(f,reps);
      CK(cudaMemcpy(hcz.data(),dCz,sizeof(cuDoubleComplex)*m*n,cudaMemcpyDeviceToHost));
      char name[64]; snprintf(name,64,"cuBLAS zgemm EMU-FP %s",strat==1?"perf":"eager");
      printf("%-28s %9.3f ms %8.2f TFLOPS  err %.2e\n",name,t*1e3,flop/t/1e12,relerr_z(hcz,ref));
    }
    CB(cublasSetEmulationStrategy(h,CUBLAS_EMULATION_STRATEGY_DEFAULT));
  }
  run_c("cuBLAS cgemm EMU 3xBF16", CUBLAS_COMPUTE_32F_EMULATED_16BFX9);

  // GEMMul8: persistent workspace (the hgw wrapper mallocs per call; that cost is NOT included here)
  for(int fast=0; fast<=1; fast++){
    for(int nm: {5,6,7,8,9,10}){
      size_t wa,wb; size_t ws=gemmul8::workSize<true,gemmul8::Backend::INT8>(m,n,k,nm,false,false,&wa,&wb);
      void* w; CK(cudaMalloc(&w,ws));
      auto f=[&]{ gemmul8::gemm<cuComplex,gemmul8::Backend::INT8>(h,opA,opB,m,n,k,&a1,dA,lda,dB,ldb,&b0,dC,ldc,nm,fast!=0,w); };
      double t=timeit(f,reps);
      CK(cudaMemcpy(hc.data(),dC,sizeof(cuComplex)*m*n,cudaMemcpyDeviceToHost));
      char name[64]; snprintf(name,64,"Oz cgemm INT8 nm=%d %s",nm,fast?"fast":"acc");
      printf("%-28s %9.3f ms %8.2f TFLOPS  err %.2e  (work %.2f GB)\n",name,t*1e3,flop/t/1e12,relerr_c(hc,ref),ws/1e9);
      CK(cudaFree(w));
    }
  }
  for(int fast=0; fast<=1; fast++){
    for(int nm: {12,14,16,18}){
      size_t wa,wb; size_t ws=gemmul8::workSize<true,gemmul8::Backend::INT8>(m,n,k,nm,false,false,&wa,&wb);
      void* w; CK(cudaMalloc(&w,ws));
      auto f=[&]{ gemmul8::gemm<cuDoubleComplex,gemmul8::Backend::INT8>(h,opA,opB,m,n,k,&a1z,dAz,lda,dBz,ldb,&b0z,dCz,ldc,nm,fast!=0,w); };
      double t=timeit(f,reps);
      CK(cudaMemcpy(hcz.data(),dCz,sizeof(cuDoubleComplex)*m*n,cudaMemcpyDeviceToHost));
      char name[64]; snprintf(name,64,"Oz zgemm INT8 nm=%d %s",nm,fast?"fast":"acc");
      printf("%-28s %9.3f ms %8.2f TFLOPS  err %.2e  (work %.2f GB)\n",name,t*1e3,flop/t/1e12,relerr_z(hcz,ref),ws/1e9);
      CK(cudaFree(w));
    }
  }
  // cudaMalloc+cudaFree per call, as gemmul8_wrapper.cu does now
  {
    int nm=7; size_t wa,wb;
    auto f=[&]{ size_t ws=gemmul8::workSize<true,gemmul8::Backend::INT8>(m,n,k,nm,false,false,&wa,&wb); void* w; CK(cudaMalloc(&w,ws));
                gemmul8::gemm<cuComplex,gemmul8::Backend::INT8>(h,opA,opB,m,n,k,&a1,dA,lda,dB,ldb,&b0,dC,ldc,nm,false,w); CK(cudaFree(w)); };
    double t=timeit(f,reps);
    printf("%-28s %9.3f ms %8.2f TFLOPS  (malloc/free per call)\n","Oz cgemm nm=7 acc +malloc",t*1e3,flop/t/1e12);
  }
  return 0;
}
