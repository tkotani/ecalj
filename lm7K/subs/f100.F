      double precision function f100(fac,j1,j2,j3,m1,m2,m3)
C-
C ----------------------------------------------------------------
Ci Inputs
Ci   FAC,J1,J2,J3,M1,M2,M3
Co Outputs
Co   F100
Cr Remarks
Cr
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer j1,j2,j3,m1,m2,m3
      double precision fac(50)
C Local parameters
      integer m,n,n1,n2
      double precision t,t1

      if (m3 .ne. m1+m2) goto 2
      t = (2*j3+1)*fac(j1+j2-j3+1)*fac(j3+j1-j2+1)*fac(j3+j2-j1+1)/
     .             fac(j1+j2+j3+2)
      t = dsqrt(t*fac(j1+m1+1)*fac(j1-m1+1)*fac(j2+m2+1)*fac(j2-m2+1)*
     .            fac(j3+m3+1)*fac(j3-m3+1))
      n1 = max0(j2-j3-m1,j1-j3+m2,0) + 1
      n2 = min0(j1+j2-j3,j1-m1,j2+m2) + 1
      if (n1 .gt. n2) goto 2
      t1 = 0.d0
      do  1  m = n1, n2
        n = m-1
        t1 = t1 +
     .    dble(1+4*(n/2)-2*n)/(fac(m)*fac(j1+j2-j3-n+1)*fac(j1-m1-n+1)*
     .    fac(j2+m2-n+1)*fac(j3-j2+m1+n+1)*fac(j3-j1-m2+n+1))
    1 continue
      f100 = t*t1
      return
    2 f100 = 0.d0
      return
      end

