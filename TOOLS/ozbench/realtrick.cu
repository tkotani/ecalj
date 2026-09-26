// realtrick: C = A^H B (complex float) by cuBLAS cgemm vs ONE real SGEMM on the same memory.
//   B (k x n complex) is read as a real (2k x n) matrix, C (m x n complex) is written as a real (2m x n) one.
//   Ap (2k x 2m real): column 2j-1 = A(:,j) as (re,im) pairs, column 2j = (-im, re) pairs.
//   Then Ap^T B = [Re C(j,:); Im C(j,:)] interleaved, i.e. C itself.  8mnk flops, as the complex product.
// Usage: realtrick m n k [reps]      (small footprint: A, Ap, B, C only)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuComplex.h>
#define CK(x) do{cudaError_t e=(x); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__); exit(1);} }while(0)
#define CB(x) do{cublasStatus_t s=(x); if(s!=CUBLAS_STATUS_SUCCESS){printf("cuBLAS %d @%d\n",(int)s,__LINE__); exit(1);} }while(0)
static double now(){ return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
__global__ void mkap(const cuComplex* a, int lda, float* ap, int k, int m){
  int i = blockIdx.x*blockDim.x + threadIdx.x, j = blockIdx.y;
  if(i>=k || j>=m) return;
  cuComplex v = a[i + (size_t)j*lda];
  size_t c1 = (size_t)(2*j)*(2*k), c2 = (size_t)(2*j+1)*(2*k);
  ap[c1 + 2*i] = v.x;  ap[c1 + 2*i+1] = v.y;
  ap[c2 + 2*i] = -v.y; ap[c2 + 2*i+1] = v.x;
}
int main(int argc,char**argv){
  size_t m=atol(argv[1]), n=atol(argv[2]), k=atol(argv[3]); int reps=argc>4?atoi(argv[4]):5;
  double flop=8.0*m*n*k;
  std::mt19937_64 rng(1); std::normal_distribution<float> g(0,1);
  std::vector<cuComplex> hA(k*m), hB(k*n);
  for(auto& v:hA) v=make_cuComplex(g(rng),g(rng));
  for(auto& v:hB) v=make_cuComplex(g(rng),g(rng));
  cuComplex *dA,*dB,*dC; float* dAp;
  CK(cudaMalloc(&dA,8*k*m)); CK(cudaMalloc(&dB,8*k*n)); CK(cudaMalloc(&dC,8*m*n)); CK(cudaMalloc(&dAp,4*4*k*m));
  CK(cudaMemcpy(dA,hA.data(),8*k*m,cudaMemcpyHostToDevice)); CK(cudaMemcpy(dB,hB.data(),8*k*n,cudaMemcpyHostToDevice));
  cublasHandle_t h; CB(cublasCreate(&h));
  cuComplex a1=make_cuComplex(1,0), b0=make_cuComplex(0,0); float one=1, zero=0;
  auto cg=[&]{ CB(cublasGemmEx(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1,dA,CUDA_C_32F,k,dB,CUDA_C_32F,k,&b0,dC,CUDA_C_32F,m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT)); };
  auto rt=[&]{ dim3 bl(256), gr((k+255)/256, m); mkap<<<gr,bl>>>(dA,k,dAp,k,m);
               CB(cublasGemmEx(h,CUBLAS_OP_T,CUBLAS_OP_N,2*m,n,2*k,&one,dAp,CUDA_R_32F,2*k,(float*)dB,CUDA_R_32F,2*k,&zero,(float*)dC,CUDA_R_32F,2*m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT)); };
  std::vector<cuComplex> c1(m*n), c2(m*n);
  cg(); CK(cudaDeviceSynchronize()); double t0=now(); for(int r=0;r<reps;r++) cg(); CK(cudaDeviceSynchronize()); double tc=(now()-t0)/reps;
  CK(cudaMemcpy(c1.data(),dC,8*m*n,cudaMemcpyDeviceToHost));
  rt(); CK(cudaDeviceSynchronize()); t0=now(); for(int r=0;r<reps;r++) rt(); CK(cudaDeviceSynchronize()); double tr=(now()-t0)/reps;
  CK(cudaMemcpy(c2.data(),dC,8*m*n,cudaMemcpyDeviceToHost));
  double emax=0,rmax=0; for(size_t i=0;i<m*n;i++){ double dx=c1[i].x-c2[i].x, dy=c1[i].y-c2[i].y; emax=fmax(emax,sqrt(dx*dx+dy*dy)); rmax=fmax(rmax,sqrt(c1[i].x*c1[i].x+c1[i].y*c1[i].y)); }
  printf("m n k = %zu %zu %zu : cgemm FP32 %.3f ms (%.1f TFLOPS) | real SGEMM trick %.3f ms (%.1f TFLOPS) | max rel diff %.2e\n",
         m,n,k, tc*1e3, flop/tc/1e12, tr*1e3, flop/tr/1e12, emax/rmax);
  return 0;
}
