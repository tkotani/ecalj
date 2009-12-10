      subroutine yhmlso(ldim,sk,so3c,hk)
C- 3-center terms of the L-S+ hamiltonian, noncollinear case.
C ----------------------------------------------------------------
Ci Inputs
Ci   sk
Ci   so3c: quasidiagonal matrix in the (1,2) block (see mksod)
Co Outputs
Co   hk:  accumulate sk * so3c * sk into hk
Cr Remarks
Cr   matrix multiplication done by blocks for efficiency
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ldim
      double precision sk(ldim,2,ldim,2,2),hk(ldim,2,ldim,2,2),so3c(1)
C Local parameters
      integer nrow,owk,owk2
      parameter (nrow=48)
      integer w(1)
      common /w/ w

      call defdr(owk,nrow*ldim*2)
      call defdr(owk2,nrow*ldim*2*2)
      call xyhmls(nrow,ldim,sk,so3c,w(owk),w(owk2),hk)
      call rlse(owk)
      end
      subroutine xyhmls(nrow,ldim,sk,so3c,wk,wk2,hk)
C     implicit none
      integer ldim,mrow,irow,nrow
      double precision so3c(1),wk(nrow,ldim,2),wk2(nrow,ldim,2,2),
     .  sk(ldim,2,ldim,2,2),hk(ldim,2,ldim,2,2)
      integer i,j

C --- Matrix multiplication by blocks of size nrow ---
      do  10  irow = 1, ldim*2, nrow
        mrow = min(2*ldim-irow+1,nrow)

C ...   sk * soc
        do  20  j = 1, ldim
        do  20  i = 1, mrow
          wk(i,j,1) = sk(i+irow-1,1,j,1,1)*so3c(j)
          wk(i,j,2) = sk(i+irow-1,1,j,1,2)*so3c(j)
   20   continue

C       call yprm('s.sod.4',2,wk,nrow*ldim,nrow,mrow,ldim)

C ...   sk * soc * sk
        call yygemm('N','N',mrow,2*ldim,ldim-1,1d0,wk,wk(1,1,2),nrow,
     .  sk(2,2,1,1,1),sk(2,2,1,1,2),ldim*2,0d0,wk2,wk2(1,1,1,2),nrow)

C       call yprm('s.sod.4.s',2,wk2,nrow*ldim*2,nrow,mrow,ldim*2)

C ...   Add into hk and Hermitian congugate into hk for L+.S- block
        do  30  j = 1, ldim*2
        do  30  i = 1, mrow
          hk(i+irow-1,1,j,1,1) = hk(i+irow-1,1,j,1,1) + wk2(i,j,1,1)
          hk(i+irow-1,1,j,1,2) = hk(i+irow-1,1,j,1,2) + wk2(i,j,1,2)
          hk(j,1,i+irow-1,1,1) = hk(j,1,i+irow-1,1,1) + wk2(i,j,1,1)
          hk(j,1,i+irow-1,1,2) = hk(j,1,i+irow-1,1,2) - wk2(i,j,1,2)
   30   continue

   10 continue

      end

