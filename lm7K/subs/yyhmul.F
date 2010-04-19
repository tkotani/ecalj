#define BLAS3 1
      subroutine yyhmul(jdim,idim,ndim,ldim,sk,ns,d,hk,nh)
C- Multiplies matrix hk = sk*d*sk, where d is a diagonal matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   ldim,idim:   row dimension of product sk*sk, and length inner
Ci                product for sk*sk (sk may be rectangular)
Ci   sk,jdim,ns:  matrix sk, its row dimension, and spacing between
Ci                real and imaginary parts
Ci   hk,ndim,nh:  matrix hk, its row dimension, and spacing between
Ci                real and imaginary parts
Ci   d: diagonal matrix
Co Outputs
Co   hk:  accumulate sk * d * sk into hk
Cr Remarks
Cr   This version uses vectorizable BLAS-style daxpy loops.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer jdim,ldim,idim,ndim,ns,nh
      double precision sk(jdim,ldim),hk(ndim,ldim),d(1,idim)
C Local parameters
      integer i,j,k,nid
      double precision xikr,xiki
      integer nrow,owk,irow,icol,mrow,ncol
      parameter (nrow=48)
      integer w(100)
      common /w/ w
C ... nid is skip length of d; set as a parameter to 1 now.
      nid = 1

C      open(55, file='snot')
C      rewind 55
C      write(55,*) idim,ldim,ndim,nid
C      call dfdump(sk,idim*ldim*2,-55)
C      call dfdump(d, nid*idim,-55)
C      call dfdump(hk,ndim*ldim*2,-55)
C      close(55)

      call defdr(owk,nrow*idim*2)
      do  10  irow = 1, ldim, nrow
        mrow = min(ldim-irow+1,nrow)
        ncol = mrow+irow-1
        call xyhmul(mrow,jdim,idim,nid,ns,d,sk(1,irow),w(owk))
        call zampy(w(owk),mrow,1,idim*mrow,sk,jdim,1,ns,
     .    hk(irow,1),ndim,1,nh,mrow,ncol,idim)
   10 continue
      do  12  irow = 1, ldim
        do  12  icol = irow+1, ldim
          hk(irow,icol) = hk(icol,irow)
          hk(nh+irow,icol) = -hk(nh+icol,irow)
   12 continue
      call rlse(owk)
      end

      subroutine xyhmul(mrow,jdim,idim,nid,ns,d,sk,w)
C Kernel called by yhmul
C     implicit none
      integer mrow,jdim,idim,nid,ns
      double precision d(nid,1),w(mrow,idim,2),sk(jdim,1)
      integer i,j

      call dpzero(w(1,1,2),mrow*idim)
      do  10  j = 1, mrow
        call dcopy(idim,sk(1,j),1,   w(j,1,1),mrow)
        call daxpy(idim,-1d0,sk(ns+1,j),1,w(j,1,2),mrow)
   10 continue

      do  20  j = 1, idim
        do  30  i = 1, mrow
          w(i,j,1) = w(i,j,1)*d(1,j)
          w(i,j,2) = w(i,j,2)*d(1,j)
   30   continue
   20 continue
      end

C ...  test yyhmul
C      subroutine fmain
C      implicit none
C      double precision sk(100*100),hk(100*100),d(100),h2(100*100),
C     .  s2(100*100)
C      common /static/ sk,hk,d,s2
C      integer jdim,idim,ldim,ndim,nid,wksize
C      parameter (wksize=10000)
C      integer w(wksize)
C      common /w/ w
C
C      call wkinit(wksize)
C      open(55, file='snot')
C      rewind 55
C      read(55,*) idim,ldim,ndim,nid
C      call dfdump(sk,idim*ldim*2,55)
C      jdim = idim+03
C      call snit(jdim,idim,ldim,sk,s2)
C      call dfdump(d, nid*idim,55)
C      call dfdump(hk,ndim*ldim*2,55)
C      close(55)
CC     call yyhmul(idim,idim,ndim+1,ldim,sk,d,nid,hk)
C      call yyhmul(jdim,idim,ndim+1,ldim,s2,d,nid,hk)
C      call prmx('output of yyhmul',hk,ndim+1,ldim,ldim)
C      end
C      subroutine snit(jdim,idim,ldim,sk,s2)
C      implicit none
C      integer jdim,idim,ldim,i,j
C      double precision sk(idim,ldim),s2(jdim,ldim)
C
C      do  10  i = 1, ldim
C      do  10  j = 1, idim
C      s2(j+jdim*ldim,i) = sk(j+idim*ldim,i)
C   10 s2(j,i) = sk(j,i)
C
C      end

