      subroutine zhmul(ndim,hk,d,nid,ck)
C- Multiplies matrix ck = hk*d*hk, where d is a diagonal matrix.
C- Only the upper triangle of ck is calculated.
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim,hk,d
Ci   nid: skip length between elements of d
Co Outputs
Co   ck
Cr Remarks
Cr   Hk and ck are dimensioned (ndim*ndim*2) (real followed by imag)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,nid
      double precision hk(ndim*ndim*2),ck(ndim*ndim*2),d(nid*ndim)
C Local parameters
      integer n2,i,j,k,jndim,j2ndim,kndim,k2ndim
      double precision xikr,xiki,xikim

      n2 = ndim**2
      do  100  j = 1, ndim
        jndim = (j-1)*ndim
        j2ndim = jndim+n2
      do  100  k = 1, ndim
        kndim = (k-1)*ndim
        k2ndim = kndim+n2
        xikr = d(k)*hk(jndim+k)
        xiki = d(k)*hk(j2ndim+k)
        xikim = -xiki
        do  110  i = 1, j
          ck(jndim+i) = ck(jndim+i) + xikr * hk(kndim+i) +
     .                                 xikim * hk(k2ndim+i)
          ck(j2ndim+i) = ck(j2ndim+i) + xiki * hk(kndim+i) +
     .                                   xikr * hk(k2ndim+i)
  110   continue
  100 continue
      end
C#elseC
C      subroutine zhmul(ndim,hk,d,nid,ck)
CC- Multiplies matrix ck = hk*d*hk, where d is a diagonal matrix
CC ----------------------------------------------------------------
CCi Inputs
CCi   ndim,hk,d
CCi   nid: skip length between elements of d
CCo Outputs
CCo   ck
CCr Remarks
CCr   Hk and ck are dimensioned (ndim,ndim,2) (real followed by imag)
CCr   The looping order is chosen so as to vectorize the inner loop
CC ----------------------------------------------------------------
CC     implicit none
CC Passed parameters
C      integer ndim,nid
C      double precision hk(ndim,ndim),ck(ndim,ndim),d(nid,1)
CC Local parameters
C      integer n2,i,j,k
C      double precision xikr,xiki
CC
C      n2 = ndim**2
C      do  100  k = 1, ndim
C        do  100  i = 1, ndim
C        xikr = d(1,k)*hk(i,k)
C        xiki = d(1,k)*hk(n2+i,k)
C#ifdefC BLAS
C        call daxpy(ndim, xikr,hk(   k,1),ndim,ck(i,1),   ndim)
C        call daxpy(ndim,-xiki,hk(n2+k,1),ndim,ck(i,1),   ndim)
C        call daxpy(ndim, xikr,hk(n2+k,1),ndim,ck(n2+i,1),ndim)
C        call daxpy(ndim, xiki,hk(   k,1),ndim,ck(n2+i,1),ndim)
C#elseC
C        do  110  j = 1, ndim
C          ck(i,j)    = ck(i,j)    + xikr*hk(k,j) - xiki*hk(n2+k,j)
C          ck(n2+i,j) = ck(n2+i,j) + xikr*hk(n2+k,j) + xiki*hk(k,j)
C  110 continue
C#endifC
C  100 continue
C      end
C#endif
      subroutine zmul0(ndim,a,b,c)
C- Multiplies complex matrix c = a*b,
C- All matrices have the same dimension
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim,a,b
Co Outputs
Co   c
Cr Remarks
Cr   a,b, and c are dimensioned (ndim*ndim*2) (real followed by imag)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim
      double precision a(ndim*ndim*2),b(ndim*ndim*2),c(ndim*ndim*2)
C Local parameters
      integer n2,i,j,k,jndim,j2ndim,kndim,k2ndim
      double precision bkj,bk2j,bk2jm

      n2 = ndim**2

C#ifdefC CRAY
C      call mxma(a,1,ndim,b,1,ndim,c,1,ndim,ndim,ndim,ndim)
C      call mxma(a,1,ndim,b(n2+1),1,ndim,c(n2+1),1,ndim,ndim,ndim,ndim)
C      call dinit(a,n2)
C      call mxma(a(n2+1),1,ndim,b(n2+1),1,ndim,a,1,ndim,ndim,ndim,ndim)
C      do 10 i = 1, n2
C        c(i) = c(i) - a(i)
C   10 continue
C      call dinit(a,n2)
C      call mxma(a(n2+1),1,ndim,b,1,ndim,a,1,ndim,ndim,ndim,ndim)
C      do 20 i = 1, n2
C        c(n2+i) = c(n2+i) + a(i)
C   20 continue
C#else
      do  100  j = 1, ndim
        jndim = (j-1)*ndim
        j2ndim = jndim+n2
      do  100  k = 1, ndim
        kndim = (k-1)*ndim
        k2ndim = kndim+n2
        bkj = b(jndim+k)
        bk2j = b(j2ndim+k)
        bk2jm = -bk2j
        do  110  i = 1, ndim
          c(jndim+i) = c(jndim+i) + bkj*a(kndim+i)
     .                                + bk2jm*a(k2ndim+i)
          c(j2ndim+i) = c(j2ndim+i) + bk2j*a(kndim+i)
     .                                  + bkj*a(k2ndim+i)
  110   continue
  100 continue
C#endif
      end
      subroutine zmul(a,nleft,b,nmid,c,nright)
C- General complex matrix copy c = a*b, cij = \sum aik*bkj
C ----------------------------------------------------------------
Ci Inputs
Ci   nleft,nmid,nright,a,b
Co Outputs
Co   c
Cr Remarks
Cr   c(nleft,nright) = a(nleft,nmid) * b(nmid,nright)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nleft,nmid,nright
      double precision a(nleft*nmid*2),b(nmid*nright*2),
     .                 c(nleft*nright*2)
C Local parameters
      integer nlm,nmr,nlr,i,j,k,jb,j2b,jc,j2c,ka,k2a
      double precision bkj,bk2j,bk2jm

      nlm = nleft * nmid
      nmr = nmid * nright
      nlr = nleft * nright

      do  100  j = 1, nright
        jc = (j-1) * nleft
        j2c = jc + nlr
        jb = (j-1) * nmid
        j2b = jb + nmr
      do  100  k = 1, nmid
        ka = (k-1) * nleft
        k2a = ka + nlm
        bkj = b(jb+k)
        bk2j = b(j2b+k)
        bk2jm = -bk2j
        do  110  i = 1, nleft
          c(jc+i)  = c(jc+i)  + bkj*a(ka+i)  + bk2jm*a(k2a+i)
          c(j2c+i) = c(j2c+i) + bk2j*a(ka+i) + bkj*a(k2a+i)
  110   continue
  100 continue
      end
      subroutine zmult(a,nleft,b,nmid,c,nright)
C- General complex matrix copy c = at*b,
C-                             cij = \sum atik*bkj = \sum aki*bkj
C- at is the transposed of a
C ----------------------------------------------------------------
Ci Inputs
Ci   nleft,nmid,nright,a,b
Co Outputs
Co   c
Cr Remarks
Cr   c(nleft,nright) = at(nleft,nmid) * b(nmid,nright)
Cr   c(nleft,nright) = a (nmid,nleft) * b(nmid,nright)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nleft,nmid,nright
      double precision a(nmid*nleft*2),b(nmid*nright*2),
     .                 c(nleft*nright*2)
C Local parameters
      integer nml,nmr,nlr,i,j,k,ia,i2a,jb,j2b,jc,j2c
      double precision cic,ci2c

      nml = nmid * nleft
      nmr = nmid * nright
      nlr = nleft * nright

      do  100  j = 1, nright
        jc = (j-1) * nleft
        j2c = jc + nlr
        jb = (j-1) * nmid
        j2b = jb + nmr
      do  100  i = 1, nleft
        ia = (i-1) * nmid
        i2a = ia + nml
        cic = 0
        ci2c = 0
        do  110  k = 1, nmid
          cic  = cic  + a(ia+k)*b(jb+k) - a(i2a+k)*b(j2b+k)
          ci2c = ci2c + a(i2a+k)*b(jb+k) + a(ia+k)*b(j2b+k)
  110   continue
        c(jc+i) = c(jc+i) + cic
        c(j2c+i) = c(j2c+i) + ci2c
  100 continue
      end
      SUBROUTINE DVCPY(A,NCA,B,NCB,N)
C- general vector copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca is the source vector and the number of
Ci      elements separating each element in the vector
Ci   b,ncb,nrb is the destination vector and the number of
Ci      elements separating each element in the vector
Ci   n: the number of elements to calculate
Co Outputs:
Co   result matrix stored in b
Cr Remarks:
Cr   This is a general-purpose vectore copy routine
Cr   Example: Set all elements of 3-by-2 matrix c to -1.d0
Cr     call dvcpy(-1.d0,0,c,1,3*2)
Cr   Example: Set block (n,m) of array a(p,m)=0, and a(i,i)=1, i=1,m
Cr     call dmcpy(0.d0,0,0,a,p,1,n,m)
Cr     call dvcpy(1.d0,0,a,p+1,m)
C ----------------------------------------------------------------
      INTEGER NCA,NCB,N
      DOUBLE PRECISION A(0:*), B(0:*)
      INTEGER I,IA,IB

      IA = N*NCA
      IB = N*NCB
      DO  200  I = N-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        B(IB) = A(IA)
  200 CONTINUE
      return
      end
      
      
      subroutine zmpy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- Double precision complex matrix multiplication
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of a (in real words)
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of b (in real words)
Ci   c,ncc,nrc is the product matrix and respectively the number of
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of c (in real words)
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of complex matrix c = a*b  (arrays b,c
Cr     dimensioned complex*16; a real*8 with imag following real)
Cr     call zmpy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b
Cr     call zmpy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   Warning: this routine has not been thoroughly checked!
Cr   This version suitable for Cray
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nrci,nccj,ncbj,nrcicj
      double precision ar,ai

C --- Initialize array to zero ---
      do  10  i = n-1, 0, -1
        nrci = nrc*i
        nccj = -ncc
      do  10  j = m-1, 0, -1
        nccj = nccj + ncc
        nrcicj = nrci + nccj
        c(nrcicj)     = 0
        c(nic+nrcicj) = 0
   10 continue

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
      do  20  i = n-1, 0, -1
        ar = a(      nra*i + nca*k)
        ai = a(nia + nra*i + nca*k)
        nrci = nrc*i
C#ifdefC BLAS
C        nrbk = nrb*k
C        call daxpy(m, ar,b(    nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m,-ai,b(nib+nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m, ar,b(nib+nrbk),ncb,c(nic+nrci),ncc)
C        call daxpy(m, ai,b(    nrbk),ncb,c(nic+nrci),ncc)
C#else
        nccj = -ncc
        ncbj = -ncb + nrb*k
      do  20  j = m-1, 0, -1
        nccj = nccj + ncc
        ncbj = ncbj + ncb
        nrcicj = nrci + nccj
        c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
        c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
C#endif
   20 continue
      end

!!#############################################################################
      integer(4) function llxxx(ilm)
      integer(4),parameter :: lmx=50
      integer(4),save:: lla((lmx+1)**2)
      logical:: init=.true.
      integer:: ilm,l,lend,lini
      if(ilm>(lmx+1)**2) call rx( 'll: ilm too large')
      if(init) then
         do l=0,lmx
            lini= l**2 + 1
            lend=(l+1)**2
            lla(lini:lend)=l
         enddo
      endif
      llxxx = lla(ilm)
      return
      end
