      subroutine dsprmv(sa,ija,x,n,b)
C- Double-precision sparse-matrix-vector multiply
C ----------------------------------------------------------------------
Ci Inputs
Ci   sa    :sparse array (see da2spr)
Ci   ija   :pointer index (see da2spr)
Ci   x     :vector
Ci   n     :dimension of result vector b
Co Outputs
Ci   b     :result vector
Cr Remarks
Cr   The column dimension is implicit in the ija.
Cu Updates
C ----------------------------------------------------------------------
      integer n,ija(*)
      double precision b(n),sa(*),x(n)
      integer i,k
      if (ija(1).ne.n+2) call rx('dsprmv: mismatched vector and matrix')
      do  12  i = 1, n
        b(i) = sa(i)*x(i)
        do  11  k = ija(i), ija(i+1)-1
          b(i) = b(i) + sa(k)*x(ija(k))
   11   continue
   12 continue
      end
C      subroutine fmain
Cc      subroutine fmain
Cc     driver for routine dsprmv
C      integer nda,nmax,mda
C      parameter(nda=5,mda=7,nmax=2*nda*mda+1)
C      integer i,j,n,ija(nmax)
C      double precision a(nda,mda),sa(nmax),ax(nda),b(nda),x(mda)
C      data a/3d0,0d0,0d0,0d0,0d0,
C     .       0d0,4d0,7d0,0d0,0d0,
C     .       1d0,0d0,5d0,0d0,0d0,
C     .       0d0,0d0,9d0,0d0,6d0,
C     .       0d0,0d0,0d0,2d0,5d0,
C     .      -1d0,-0d0,-3d0,-4d0,-5d0,
C     .      -0d0,-12d0,-13d0,-14d0,-0d0/
CC      data a/3d0,0d0,0d0,0d0,0d0,
CC     .     0d0,4d0,7d0,0d0,0d0,
CC     .     1d0,0d0,5d0,0d0,0d0,
CC     .     0d0,0d0,9d0,0d0,6d0,
CC     .     0d0,0d0,0d0,2d0,5d0/
C      data x/1d0,2d0,3d0,4d0,5d0,-6d0,-7d0/
C
C      n = nda
C      m = mda
C   10 continue
C      print 345, n,m
C  345 format(' multiply a x with a dimensioned',i2,',',i2)
C      call da2spr(a,nda,n,m,1d0,nmax,sa,ija)
C      n = ija(1)-2
C      call dsprmv(sa,ija,x,n,b)
C      do  12  i = 1, n
C        ax(i) = 0d0
C        do  11  j = 1, m
C          ax(i) = ax(i) + a(i,j)*x(j)
C   11   continue
C   12 continue
C
C      write(*,'(t4,a,t18,a,t34,a)') 'reference','dsprmv result','diff'
C      do  13  i = 1, n
C        write(*,'(t4,f6.2,t22,f6.2,t32,f6.2)') ax(i),b(i),ax(i)-b(i)
C   13 continue
C
C      if (m .gt. nda) then
C        m = nda-1
C        print *, 'redo test, now with n > m'
C        goto 10
C      endif
C
C      end

