// matinvbench: inverse of one n x n complex matrix like epstilde (W-build matinv, 318 of them per q).
//  ref  : cuSOLVER zgetrf + zgetrs(RHS = I) in FP64
//  mixed: cgetrf + cgetrs in FP32, then Newton-Schulz X <- X (2I - A X) twice with FP64 GEMMs,
//         the GEMMs either native FP64 or cuBLAS's FP64 emulation (fixed point on INT8 tensor cores).
// Error: max|X - X_ref| / max|X_ref| and max|I - A X|.
// Usage: matinvbench n [cond_decades] [reps]   (A = I - 0.5 G/sqrt(n) with singular values stretched by cond)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuComplex.h>
#define CK(x) do{cudaError_t e=(x); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__); exit(1);} }while(0)
#define CB(x) do{cublasStatus_t s=(x); if(s!=CUBLAS_STATUS_SUCCESS){printf("cuBLAS %d @%d\n",(int)s,__LINE__); exit(1);} }while(0)
#define CS(x) do{cusolverStatus_t s=(x); if(s!=CUSOLVER_STATUS_SUCCESS){printf("cuSOLVER %d @%d\n",(int)s,__LINE__); exit(1);} }while(0)
static double now(){ return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
__global__ void z2c(const cuDoubleComplex* a, cuComplex* b, size_t n){ size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x; if(i<n) b[i]=make_cuComplex((float)a[i].x,(float)a[i].y); }
__global__ void c2z(const cuComplex* a, cuDoubleComplex* b, size_t n){ size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x; if(i<n) b[i]=make_cuDoubleComplex(a[i].x,a[i].y); }
__global__ void eyec(cuComplex* b, int n){ size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x; if(i<(size_t)n*n) b[i]=make_cuComplex((i%n)==(i/n)?1.f:0.f,0.f); }
__global__ void eyez(cuDoubleComplex* b, int n){ size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x; if(i<(size_t)n*n) b[i]=make_cuDoubleComplex((i%n)==(i/n)?1.:0.,0.); }
__global__ void twoIminus(cuDoubleComplex* t, int n){ size_t i=blockIdx.x*(size_t)blockDim.x+threadIdx.x; if(i<(size_t)n*n){ t[i].x=-t[i].x; t[i].y=-t[i].y; if((i%n)==(i/n)) t[i].x+=2.; } }
int main(int argc,char**argv){
  int n=atoi(argv[1]); double dec=argc>2?atof(argv[2]):0; int reps=argc>3?atoi(argv[3]):5;
  size_t nn=(size_t)n*n;
  std::mt19937_64 rng(3); std::normal_distribution<double> g(0,1);
  std::vector<cuDoubleComplex> hA(nn);
  // A = I - 0.5 G/sqrt(n), then scale column j by 10^(-dec*j/n) to stretch the condition number
  for(int j=0;j<n;j++){ double s=pow(10.0,-dec*j/n); for(int i=0;i<n;i++){ double re=-0.5*g(rng)/sqrt(n), im=-0.5*g(rng)/sqrt(n); if(i==j) re+=1; hA[i+(size_t)j*n]=make_cuDoubleComplex(re*s,im*s);} }
  cuDoubleComplex *dA,*dLU,*dX,*dXr,*dT; cuComplex *dA32,*dX32; int *ipiv,*info;
  CK(cudaMalloc(&dA,16*nn)); CK(cudaMalloc(&dLU,16*nn)); CK(cudaMalloc(&dX,16*nn)); CK(cudaMalloc(&dXr,16*nn)); CK(cudaMalloc(&dT,16*nn));
  CK(cudaMalloc(&dA32,8*nn)); CK(cudaMalloc(&dX32,8*nn)); CK(cudaMalloc(&ipiv,4*n)); CK(cudaMalloc(&info,4));
  CK(cudaMemcpy(dA,hA.data(),16*nn,cudaMemcpyHostToDevice));
  cusolverDnHandle_t sh; CS(cusolverDnCreate(&sh)); cublasHandle_t bh; CB(cublasCreate(&bh));
  int lwz, lwc; CS(cusolverDnZgetrf_bufferSize(sh,n,n,dLU,n,&lwz)); CS(cusolverDnCgetrf_bufferSize(sh,n,n,dA32,n,&lwc));
  cuDoubleComplex* wz; cuComplex* wc; CK(cudaMalloc(&wz,16*(size_t)lwz)); CK(cudaMalloc(&wc,8*(size_t)lwc));
  int nb=(int)((nn+255)/256);
  auto ref=[&]{ CK(cudaMemcpy(dLU,dA,16*nn,cudaMemcpyDeviceToDevice)); CS(cusolverDnZgetrf(sh,n,n,dLU,n,wz,ipiv,info));
                eyez<<<nb,256>>>(dXr,n); CS(cusolverDnZgetrs(sh,CUBLAS_OP_N,n,n,dLU,n,ipiv,dXr,n,info)); };
  cuDoubleComplex one=make_cuDoubleComplex(1,0), zero=make_cuDoubleComplex(0,0);
  auto mixed=[&](cublasComputeType_t ct, int nit){
    z2c<<<nb,256>>>(dA,dA32,nn); CS(cusolverDnCgetrf(sh,n,n,dA32,n,wc,ipiv,info));
    eyec<<<nb,256>>>(dX32,n); CS(cusolverDnCgetrs(sh,CUBLAS_OP_N,n,n,dA32,n,ipiv,dX32,n,info));
    c2z<<<nb,256>>>(dX32,dX,nn);
    for(int it=0;it<nit;it++){
      CB(cublasGemmEx(bh,CUBLAS_OP_N,CUBLAS_OP_N,n,n,n,&one,dA,CUDA_C_64F,n,dX,CUDA_C_64F,n,&zero,dT,CUDA_C_64F,n,ct,CUBLAS_GEMM_DEFAULT));
      twoIminus<<<nb,256>>>(dT,n);
      CB(cublasGemmEx(bh,CUBLAS_OP_N,CUBLAS_OP_N,n,n,n,&one,dX,CUDA_C_64F,n,dT,CUDA_C_64F,n,&zero,dLU,CUDA_C_64F,n,ct,CUBLAS_GEMM_DEFAULT));
      CK(cudaMemcpy(dX,dLU,16*nn,cudaMemcpyDeviceToDevice));
    } };
  auto err=[&](cuDoubleComplex* x, const char* name, double t){
    std::vector<cuDoubleComplex> hx(nn), hr(nn), ht(nn);
    CK(cudaMemcpy(hx.data(),x,16*nn,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(hr.data(),dXr,16*nn,cudaMemcpyDeviceToHost));
    double e=0,r=0; for(size_t i=0;i<nn;i++){ double dx=hx[i].x-hr[i].x, dy=hx[i].y-hr[i].y; e=fmax(e,sqrt(dx*dx+dy*dy)); r=fmax(r,sqrt(hr[i].x*hr[i].x+hr[i].y*hr[i].y)); }
    // residual I - A X with a native FP64 product
    CB(cublasSetEmulationStrategy(bh,CUBLAS_EMULATION_STRATEGY_DEFAULT));
    CB(cublasZgemm(bh,CUBLAS_OP_N,CUBLAS_OP_N,n,n,n,&one,dA,n,x,n,&zero,dT,n));
    CK(cudaMemcpy(ht.data(),dT,16*nn,cudaMemcpyDeviceToHost));
    double res=0; for(size_t i=0;i<nn;i++){ double re=ht[i].x-((i%n)==(i/n)?1.:0.), im=ht[i].y; res=fmax(res,sqrt(re*re+im*im)); }
    printf("%-34s %8.3f ms   err %.2e   max|I-AX| %.2e\n", name, t*1e3, e/r, res);
  };
  ref(); CK(cudaDeviceSynchronize()); double t0=now(); for(int r=0;r<reps;r++) ref(); CK(cudaDeviceSynchronize()); double tref=(now()-t0)/reps;
  printf("n=%d cond decades=%.1f\n", n, dec);
  err(dXr,"zgetrf+zgetrs FP64 (reference)",tref);
  struct V{ const char* name; cublasComputeType_t ct; int strat; int nit; } vs[]={
    {"cgetrf+cgetrs FP32 only",CUBLAS_COMPUTE_64F,0,0},
    {"FP32 + 1 Newton, native FP64",CUBLAS_COMPUTE_64F,0,1},
    {"FP32 + 2 Newton, native FP64",CUBLAS_COMPUTE_64F,0,2},
    {"FP32 + 1 Newton, emulated FP64",CUBLAS_COMPUTE_64F_EMULATED_FIXEDPOINT,2,1},
    {"FP32 + 2 Newton, emulated FP64",CUBLAS_COMPUTE_64F_EMULATED_FIXEDPOINT,2,2}};
  for(auto& v: vs){
    CB(cublasSetEmulationStrategy(bh,(cublasEmulationStrategy_t)v.strat));
    mixed(v.ct,v.nit); CK(cudaDeviceSynchronize()); t0=now(); for(int r=0;r<reps;r++) mixed(v.ct,v.nit); CK(cudaDeviceSynchronize());
    double t=(now()-t0)/reps;
    CB(cublasSetEmulationStrategy(bh,(cublasEmulationStrategy_t)v.strat));
    err(dX,v.name,t);
  }
  return 0;
}
