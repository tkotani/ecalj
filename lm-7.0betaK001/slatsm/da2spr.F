      subroutine da2spr(a,nda,n,m,thresh,nmax,sa,ija)
C- Copy a double precision array into sparse format
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :array in normal format
Ci   nda   :leading dimension of a
Ci   n     :number of rows to copy
Ci   m     :number of columns to copy
Ci   thresh:threshold, below which value is deemed to be zero
Ci   nmax  :maximum size of sparse array
Co Outputs
Co   sa    :nonzero elements of sparse array
Co          sa(i), i=1..n diagonal elements of a
Co          sa(i), i=n+2  off-diagonal elements of a, grouped by rows
Co          Must be dimensioned at least= 1+n-m + (# nonzero elements).
Co   ija   :pointer index, with:
Co          ija(1)= n+2
Co          ija(i), i = 1,..n+1 = points to first entry in sa for row i
Co          ija(k), k = n+2... column index of element sa(k).  Thus
Co          for row i, k ranges from ija(i) <= k < ija(i+1), and
Co          sum_j a_ij x_j  = sum sa_k x_(ija(k))
Cr Remarks
Cr   This routine was adapted from Numerical Recipes, which see for
Cr   a description of the storage conventions.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,m,nmax,nda,ija(nmax)
      double precision thresh,a(nda,m),sa(nmax)
C ... Local parameters
      integer i,j,k

C ... Assign diagonal elements first.  sa(min(n,m)+1,n) are not used.
      k = min(n,m)
      do  10  j = 1, k
        sa(j) = a(j,j)
   10 continue

      ija(1) = n+2
      k = n+1
      do  13  i = 1, n
        do  12  j = 1, m
          if (abs(a(i,j)) .ge. thresh) then
            if (i .ne. j) then
              k = k+1
              if (k .gt. nmax) call rx('nmax too small in da2spr')
              sa(k) = a(i,j)
              ija(k) = j
            endif
          endif
   12   continue
        ija(i+1) = k+1
   13 continue
      call awrit2(' da2spr: copied %i nonzero elements out of %i',
     .  ' ',80,6,k,n*m)
      end
C      subroutine fmain
CC     Test da2spr
C      implicit none
C      integer nda,nmax,mda
C      parameter(nda=5,mda=7,nmax=2*nda*mda+1)
C      integer i,j,msize,ija(nmax),m
C      double precision a(nda,mda),aa(nda,mda),sa(nmax)
C      character*7 strn
C      data a/3d0,0d0,0d0,0d0,0d0,
C     .       0d0,4d0,7d0,0d0,0d0,
C     .       1d0,0d0,5d0,0d0,0d0,
C     .       0d0,0d0,9d0,0d0,6d0,
C     .       0d0,0d0,0d0,2d0,5d0,
C     .      -1d0,-0d0,-3d0,-4d0,-5d0,
C     .      -0d0,-12d0,-13d0,-14d0,-0d0/
CC      data a/1.1d0,1.2d0,1.3d0,1.4d0,1.5d0,
CC     .       2.1d0,2.2d0,2.3d0,2.4d0,2.5d0,
CC     .       3.1d0,3.2d0,3.3d0,3.4d0,3.5d0,
CC     .       4.1d0,4.2d0,4.3d0,4.4d0,4.5d0,
CC     .       5.1d0,5.2d0,5.3d0,5.4d0,5.5d0,
CC     .      -1d0,-0d0,-3d0,-4d0,-5d0,
CC     .      -0d0,-12d0,-13d0,-14d0,-0d0/
C
C      strn = '(7f7.1)'
C      m = mda
C   10 continue
C      call da2spr(a,nda,nda,m,0.1d0,nmax,sa,ija)
C      msize = ija(ija(1)-1)-1
C      sa(nda+1) = 0d0
C      write(*,'(t4,a,t18,a,t24,a)') 'index','ija','sa'
C      do  11  i = 1, msize
C        write(*,'(t2,i4,t16,i4,t20,f12.6)') i,ija(i),sa(i)
C   11 continue
C      do  13  i = 1, nda
C        do  12  j = 1, m
C          aa(i,j) = 0d0
C   12   continue
C   13 continue
C      do  15  i = 1, nda
C        aa(i,i) = sa(i)
C        do  14  j = ija(i), ija(i+1)-1
C          aa(i,ija(j)) = sa(j)
C   14   continue
C   15 continue
C      write(*,*) 'original matrix:'
C      write(*,strn) ((a(i,j),j=1,m),i = 1,nda)
C      write(*,*) 'reconstructed matrix:'
C      write(*,strn) ((aa(i,j),j=1,m),i = 1,nda)
C
C      if (m .gt. nda) then
C        m = nda-1
C        strn = '(4f7.1)'
C        print *, 'redo, with n > m'
C        goto 10
C      endif
C
C      end

