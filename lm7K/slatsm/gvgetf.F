      subroutine gvgetf(ng,n,kv,k1,k2,k3,c,c0)
C- Gathers Fourier coefficients from 3D array c into list c0.
C     implicit none
      integer ng,n,k1,k2,k3,kv(ng,3)
      double complex c0(ng,n),c(k1,k2,k3,n)
      integer ig,i,j1,j2,j3
      do  10  i = 1, n
        do  10  ig = 1, ng
          j1 = kv(ig,1)
          j2 = kv(ig,2)
          j3 = kv(ig,3)
          c0(ig,i) = c(j1,j2,j3,i)
   10 continue
      end
      subroutine gvputf(ng,n,kv,k1,k2,k3,c0,c)
C- Pokes Fourier coefficients from gathered list c0 into 3D array c.
C     implicit none
      integer ng,n,k1,k2,k3,kv(ng,3)
      double complex c0(ng,n),c(k1,k2,k3,n)
      integer ig,i,j1,j2,j3
      call dpzero(c,2*k1*k2*k3*n)
      do  10  i = 1, n
        do  10  ig = 1, ng
          j1 = kv(ig,1)
          j2 = kv(ig,2)
          j3 = kv(ig,3)
          c(j1,j2,j3,i) = c0(ig,i)
   10 continue
      end
      subroutine gvaddf(ng,kv,k1,k2,k3,c0,c)
C-  Adds Fourier coefficients from gathered list c0 into 3D array c.
C     implicit none
      integer ng,k1,k2,k3,kv(ng,3)
      double complex c0(ng),c(k1,k2,k3)
      integer ig,j1,j2,j3
      do  10  ig = 1, ng
        j1 = kv(ig,1)
        j2 = kv(ig,2)
        j3 = kv(ig,3)
        c(j1,j2,j3) = c(j1,j2,j3) + c0(ig)
   10 continue
      end

