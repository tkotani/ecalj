module m_zhev_gpu_handles
#ifdef __GPU
  use cusolverdn
  use cublas_v2
  type(cusolverDnHandle) :: zhev_cusolver_handle
  type(cublasHandle) :: zhev_cublas_handle
  logical :: zhev_gpu_handles_init = .false.
#endif
end module m_zhev_gpu_handles
