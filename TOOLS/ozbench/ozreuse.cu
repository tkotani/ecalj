// ozreuse: GEMMul8 with the scaled/split form of a fixed matrix reused (skip_scalA / skip_scalB), against
// cuBLAS cgemm FP32 and the real-SGEMM route, for the repeated products of hgw (C = A^H B, complex float):
//   case R: Sigma_c real axis  -- A = W(omega) fixed (k x m), B changes (k x n, n small), many calls
//   case I: Sigma_c imag. axis -- B = zmel fixed (k x n, n large), A = W(i omega) changes, 11 calls
//   case M: A fixed, B large and changing (zmel M2E; also W(i omega) across KXloops)
// Usage: ozreuse m k nR nI [reps] [moduli]      e.g. ozreuse 1053 1053 400 49928 40 7
// Time per call (ms), effective TFLOPS (8mnk), max rel. error against FP64 zgemm for one call.
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
__global__ void mkap(const cuComplex* a, int lda, float* ap, int k, int m){
  int i = blockIdx.x*blockDim.x + threadIdx.x, j = blockIdx.y;
  if(i>=k || j>=m) return;
  cuComplex v = a[i + (size_t)j*lda];
  size_t c1 = (size_t)(2*j)*(2*k), c2 = (size_t)(2*j+1)*(2*k);
  ap[c1 + 2*i] = v.x;  ap[c1 + 2*i+1] = v.y;
  ap[c2 + 2*i] = -v.y; ap[c2 + 2*i+1] = v.x;
}
struct Mat { cuComplex* d; cuDoubleComplex* z; size_t r, c; };
static Mat mk(size_t r, size_t c, std::mt19937_64& rng){
  std::normal_distribution<float> g(0,1); Mat M{nullptr,nullptr,r,c};
  std::vector<cuComplex> h(r*c); for(auto& v:h) v=make_cuComplex(g(rng),g(rng));
  std::vector<cuDoubleComplex> hz(r*c); for(size_t i=0;i<r*c;i++) hz[i]=make_cuDoubleComplex(h[i].x,h[i].y);
  CK(cudaMalloc(&M.d,8*r*c)); CK(cudaMemcpy(M.d,h.data(),8*r*c,cudaMemcpyHostToDevice));
  CK(cudaMalloc(&M.z,16*r*c)); CK(cudaMemcpy(M.z,hz.data(),16*r*c,cudaMemcpyHostToDevice));
  return M;
}
int main(int argc,char**argv){
  size_t m=atol(argv[1]), k=atol(argv[2]), nR=atol(argv[3]), nI=atol(argv[4]);
  int reps=argc>5?atoi(argv[5]):40, nm=argc>6?atoi(argv[6]):7;
  std::mt19937_64 rng(7);
  cublasHandle_t h; CB(cublasCreate(&h));
  cuComplex a1=make_cuComplex(1,0), b0=make_cuComplex(0,0); float one=1, zero=0;
  cuDoubleComplex a1z=make_cuDoubleComplex(1,0), b0z=make_cuDoubleComplex(0,0);
  float* ap; CK(cudaMalloc(&ap,16*k*m));
  auto err=[&](cuComplex* c, cuDoubleComplex* cref, size_t n){
    std::vector<cuComplex> x(m*n); std::vector<cuDoubleComplex> r(m*n);
    CK(cudaMemcpy(x.data(),c,8*m*n,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(r.data(),cref,16*m*n,cudaMemcpyDeviceToHost));
    double e=0,s=0; for(size_t i=0;i<m*n;i++){ double dx=x[i].x-r[i].x, dy=x[i].y-r[i].y; e=fmax(e,sqrt(dx*dx+dy*dy)); s=fmax(s,sqrt(r[i].x*r[i].x+r[i].y*r[i].y)); }
    return e/s; };
  auto line=[&](const char* nm_, double t, double flop, double e){ printf("  %-34s %8.3f ms %7.1f TFLOPS  err %.1e\n", nm_, t*1e3, flop/t/1e12, e); };
  // generic GEMMul8 call with persistent workA/workB
  auto oz=[&](cuComplex* A, cuComplex* B, cuComplex* C, size_t n, void* w, void* wA, void* wB, bool eA, bool eB, bool sA, bool sB){
    gemmul8::gemm<cuComplex,gemmul8::Backend::INT8>(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1,A,k,B,k,&b0,C,m,nm,true,w,wA,wB,eA,eB,sA,sB); };
  auto ozws=[&](size_t n, bool eA, bool eB, void** w, void** wA, void** wB){
    size_t sa,sb; size_t tot=gemmul8::workSize<true,gemmul8::Backend::INT8>(m,n,k,nm,eA,eB,&sa,&sb,true);
    CK(cudaMalloc(wA,sa?sa:1)); CK(cudaMalloc(wB,sb?sb:1)); CK(cudaMalloc(w,tot-sa-sb>0?tot-sa-sb:1)); return tot; };

  // ---------- case R and M: A fixed, B changes ----------
  for(size_t n : {nR, 2*nR+200, nI}){
    const char* tag = (n==nI) ? "case M: A fixed, B large and changing" : "case R: A = W(omega) fixed, B small and changing";
    int nb=4, rr = (n==nI) ? 6 : reps;
    Mat A=mk(k,m,rng); std::vector<Mat> Bs; for(int i=0;i<nb;i++) Bs.push_back(mk(k,n,rng));
    cuComplex* C; CK(cudaMalloc(&C,8*m*n)); cuDoubleComplex* Cz; CK(cudaMalloc(&Cz,16*m*n));
    CB(cublasZgemm(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1z,A.z,k,Bs[0].z,k,&b0z,Cz,m));
    double flop=8.0*m*n*k; printf("%s   m n k = %zu %zu %zu  (%d calls)\n", tag, m,n,k, rr);
    double t0,t; CK(cudaDeviceSynchronize());
    t0=now(); for(int r=0;r<rr;r++) CB(cublasGemmEx(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1,A.d,CUDA_C_32F,k,Bs[r%nb].d,CUDA_C_32F,k,&b0,C,CUDA_C_32F,m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT));
    CK(cudaDeviceSynchronize()); t=(now()-t0)/rr;
    CB(cublasGemmEx(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1,A.d,CUDA_C_32F,k,Bs[0].d,CUDA_C_32F,k,&b0,C,CUDA_C_32F,m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT));
    line("cuBLAS cgemm FP32",t,flop,err(C,Cz,n));
    dim3 bl(256), gr((k+255)/256, m);
    t0=now(); for(int r=0;r<rr;r++){ mkap<<<gr,bl>>>(A.d,k,ap,k,m);
      CB(cublasGemmEx(h,CUBLAS_OP_T,CUBLAS_OP_N,2*m,n,2*k,&one,ap,CUDA_R_32F,2*k,(float*)Bs[r%nb].d,CUDA_R_32F,2*k,&zero,(float*)C,CUDA_R_32F,2*m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT)); }
    CK(cudaDeviceSynchronize()); t=(now()-t0)/rr; line("real SGEMM route (A' each call)",t,flop,-1);
    mkap<<<gr,bl>>>(A.d,k,ap,k,m);
    t0=now(); for(int r=0;r<rr;r++)
      CB(cublasGemmEx(h,CUBLAS_OP_T,CUBLAS_OP_N,2*m,n,2*k,&one,ap,CUDA_R_32F,2*k,(float*)Bs[r%nb].d,CUDA_R_32F,2*k,&zero,(float*)C,CUDA_R_32F,2*m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT));
    CK(cudaDeviceSynchronize()); t=(now()-t0)/rr;
    CB(cublasGemmEx(h,CUBLAS_OP_T,CUBLAS_OP_N,2*m,n,2*k,&one,ap,CUDA_R_32F,2*k,(float*)Bs[0].d,CUDA_R_32F,2*k,&zero,(float*)C,CUDA_R_32F,2*m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT));
    line("real SGEMM route (A' kept)",t,flop,err(C,Cz,n));
    { void *w,*wA,*wB; ozws(n,false,false,&w,&wA,&wB);
      t0=now(); for(int r=0;r<rr;r++) oz(A.d,Bs[r%nb].d,C,n,w,wA,wB,false,false,false,false); CK(cudaDeviceSynchronize()); t=(now()-t0)/rr;
      oz(A.d,Bs[0].d,C,n,w,wA,wB,false,false,false,false); CK(cudaDeviceSynchronize());
      char s[64]; snprintf(s,64,"Ozaki %d moduli, no reuse",nm); line(s,t,flop,err(C,Cz,n)); cudaFree(w); cudaFree(wA); cudaFree(wB); }
    { void *w,*wA,*wB; ozws(n,true,false,&w,&wA,&wB);
      oz(A.d,Bs[0].d,C,n,w,wA,wB,true,false,false,false); CK(cudaDeviceSynchronize());       // first call scales A into wA
      t0=now(); for(int r=0;r<rr;r++) oz(A.d,Bs[r%nb].d,C,n,w,wA,wB,true,false,true,false); CK(cudaDeviceSynchronize()); t=(now()-t0)/rr;
      oz(A.d,Bs[0].d,C,n,w,wA,wB,true,false,true,false); CK(cudaDeviceSynchronize());
      char s[64]; snprintf(s,64,"Ozaki %d moduli, A reused",nm); line(s,t,flop,err(C,Cz,n)); cudaFree(w); cudaFree(wA); cudaFree(wB); }
    cudaFree(A.d); cudaFree(A.z); for(auto& B:Bs){ cudaFree(B.d); cudaFree(B.z);} cudaFree(C); cudaFree(Cz);
  }
  // ---------- case I: B = zmel fixed, 11 different A ----------
  {
    size_t n=nI; int na=11; Mat B=mk(k,n,rng); std::vector<Mat> As; for(int i=0;i<na;i++) As.push_back(mk(k,m,rng));
    cuComplex* C; CK(cudaMalloc(&C,8*m*n)); cuDoubleComplex* Cz; CK(cudaMalloc(&Cz,16*m*n));
    CB(cublasZgemm(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1z,As[na-1].z,k,B.z,k,&b0z,Cz,m));
    double flop=8.0*m*n*k; printf("case I: B = zmel fixed, 11 different W(i omega)   m n k = %zu %zu %zu\n", m,n,k);
    double t0,t; CK(cudaDeviceSynchronize());
    t0=now(); for(int i=0;i<na;i++) CB(cublasGemmEx(h,CUBLAS_OP_C,CUBLAS_OP_N,m,n,k,&a1,As[i].d,CUDA_C_32F,k,B.d,CUDA_C_32F,k,&b0,C,CUDA_C_32F,m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT));
    CK(cudaDeviceSynchronize()); t=(now()-t0)/na; line("cuBLAS cgemm FP32",t,flop,err(C,Cz,n));
    dim3 bl(256), gr((k+255)/256, m);
    t0=now(); for(int i=0;i<na;i++){ mkap<<<gr,bl>>>(As[i].d,k,ap,k,m);
      CB(cublasGemmEx(h,CUBLAS_OP_T,CUBLAS_OP_N,2*m,n,2*k,&one,ap,CUDA_R_32F,2*k,(float*)B.d,CUDA_R_32F,2*k,&zero,(float*)C,CUDA_R_32F,2*m,CUBLAS_COMPUTE_32F,CUBLAS_GEMM_DEFAULT)); }
    CK(cudaDeviceSynchronize()); t=(now()-t0)/na; line("real SGEMM route",t,flop,err(C,Cz,n));
    { void *w,*wA,*wB; ozws(n,false,false,&w,&wA,&wB);
      t0=now(); for(int i=0;i<na;i++) oz(As[i].d,B.d,C,n,w,wA,wB,false,false,false,false); CK(cudaDeviceSynchronize()); t=(now()-t0)/na;
      char s[64]; snprintf(s,64,"Ozaki %d moduli, no reuse",nm); line(s,t,flop,err(C,Cz,n)); cudaFree(w); cudaFree(wA); cudaFree(wB); }
    { void *w,*wA,*wB; ozws(n,false,true,&w,&wA,&wB);
      t0=now(); for(int i=0;i<na;i++) oz(As[i].d,B.d,C,n,w,wA,wB,false,true,false,i>0); CK(cudaDeviceSynchronize()); t=(now()-t0)/na;
      char s[64]; snprintf(s,64,"Ozaki %d moduli, B reused",nm); line(s,t,flop,err(C,Cz,n)); cudaFree(w); cudaFree(wA); cudaFree(wB); }
  }
  return 0;
}
