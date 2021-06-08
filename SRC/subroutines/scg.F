      subroutine scg(lmax,c,cindx,js)
C- Computes Clebsch-Gordan coefficients
C ----------------------------------------------------------------
Ci Inputs
Ci   lmax
Co Outputs
Co   c,cindx,js
Cr Remarks
Cr   (FORMERLY S104 IN ASW)
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer lmax
      integer cindx(1),js(1)
      double precision c(1)
C Local parameters
      integer i,i1,i2,i3,i31,i32,ic,j1,j1s,j2,j2s,k2,l1,l2,l3,lmindx,
     .m1,m2,m3,mb,n1,n2,n3,nl,nm3,s1,s2,s3,t1,t2,t3
      double precision q1,sr2,t,srpi,fs,f100,f102
      double precision fac(161)
      external f100,f102
      data srpi /1.772453851d0/
      fs(i) = 1 + 4*(i/2) - 2*i

      mb = 999999
      nl = lmax+1
      sr2 = dsqrt(2.d0)
      fac(1) = 1.d0
      do  11  i = 1, 160
        fac(i+1) = i*fac(i)
   11 continue
      ic = 0
      lmindx = 0
      do  104  i1 = 1, nl
        l1 = i1-1
        j1s = 2*l1+1
      do  103  j1 = 1, j1s
        m1 = j1-i1
        n1 = iabs(m1)
        s1 = 0
        if (m1 .lt. 0) s1 = 1
        t1 = 0
        if (m1 .eq. 0) t1 = 1
      do  102  i2 = 1, i1
        l2 = i2-1
        i31 = l1 - l2 + 1
        i32 = l1 + l2 + 1
        j2s = 2*l2 + 1
        k2 = j1s*j2s
        if (i2 .eq. i1) j2s = j1
      do  101  j2 = 1, j2s
        lmindx = lmindx + 1
        cindx(lmindx) = ic+1
        m2 = j2-i2
        n2 = iabs(m2)
        s2 = 0
        if (m2 .lt. 0) s2 = 1
        t2 = 0
        if (m2 .eq. 0) t2 = 1
        if (m1*m2) 2,3,4
    2   m3 = -n1 - n2
        mb = -iabs(n1-n2)
        if (mb .eq. 0) goto 21
        nm3 = 2
        goto 5
    3   m3 = m1+m2
   21   nm3 = 1
        goto 5
    4   m3 = n1+n2
        mb = iabs(n1-n2)
        nm3 = 2
    5   n3 = iabs(m3)
        s3 = 0
        if (m3 .lt. 0) s3 = 1
        t3 = 0
        if (m3 .eq. 0) t3 = 1
        q1 = dsqrt(dble(k2))*fs(n3+(s1+s2+s3)/2)/(2*sr2**(1+t1+t2+t3))
        do  6  i3 = i31, i32,2
          l3 = i3-1
          if (n3 .gt. l3) goto 6
          t = 0.d0
          if (n1+n2 .eq. -n3) t = t + f102(fac,l1,l2,l3)
          if (n1+n2 .eq. n3)
     .          t = t + f100(fac,l1,l2,l3,n1,n2,n3)*fs(n3+s3)
          if (n1-n2 .eq. -n3)
     .          t = t + f100(fac,l1,l2,l3,n1,-n2,-n3)*fs(n2+s2)
          if (n1-n2 .eq. n3)
     .          t = t + f100(fac,l1,l2,l3,-n1,n2,-n3)*fs(n1+s1)
          ic = ic+1
          c(ic) = q1*t*f102(fac,l1,l2,l3)/(srpi*dsqrt(dble(2*l3+1)))
          js(ic) = l3*(l3+1) + m3 + 1
    6   continue
        nm3 = nm3-1
        m3 = mb
        if (nm3 .gt. 0) goto 5
 101  continue
 102  continue
 103  continue
 104  continue
      cindx(lmindx+1) = ic+1
      return
      end

      subroutine scg_sizechk(lmax,lnjcg,lnxcg) !(lmax,c,cindx,js)
c  computes clebsch-gordan coefficients (formerly s104 in asw)
c  but here all is doubleprecision
      implicit real*8 (a-h,p-z)
      implicit integer(o)
      integer s1,s2,s3,t1,t2,t3,lmax,lnjcg,lnxcg,i1,i2,i3,i31,i32,j1,j1s,j2,k1,l1,l2,j2s,l3
     & ,lmindx,m1,m2,m3,mb,ic,k2,n2,n3,nl,nm3,n1
c ,cindx(1),js(1)
c      doubleprecision fac(50),c(1)
c      print *,' scg_sizechk:'
      lnjcg=0
      lnxcg=0
c      fs(i)=dfloat(1+4*(i/2)-2*i)
      mb = 999999
      srpi = dsqrt(4*datan(1d0))
      nl=lmax+1
      sr2=dsqrt(2.d0)
c      fac(1)=1.d0
c      do 11 i=1,49
c  11  fac(i+1)=dfloat(i)*fac(i)
      ic=0
      lmindx=0
      do 1111 i1=1,nl
        l1=i1-1
        j1s=2*l1+1
      do 111 j1=1,j1s
        m1=j1-i1
        n1=iabs(m1)
        s1=0
        if(m1.lt.0) s1=1
        t1=0
        if(m1.eq.0) t1=1
      do 11 i2=1,i1
        l2=i2-1
        i31=l1-l2+1
        i32=l1+l2+1
        j2s=2*l2+1
        k2=j1s*j2s
        if(i2.eq.i1) j2s=j1
      do 1 j2=1,j2s
        lmindx=lmindx+1
c      cindx(lmindx)=ic+1
        lnxcg=max(lnxcg,lmindx)
        m2=j2-i2
        n2=iabs(m2)
        s2=0
        if(m2.lt.0) s2=1
        t2=0
        if(m2.eq.0) t2=1
        if(m1*m2) 2,3,4
 2      continue
        m3=-n1-n2
        mb=-iabs(n1-n2)
        if(mb.eq.0) goto 21
        nm3=2
        goto 5
 3      continue
        m3=m1+m2
  21    nm3=1
        goto 5
 4      continue
        m3=n1+n2
        mb=iabs(n1-n2)
        nm3=2
 5      continue
        n3=iabs(m3)
        s3=0
        if(m3.lt.0) s3=1
        t3=0
        if(m3.eq.0) t3=1
c      q1=dsqrt(dfloat(k2))*fs(n3+(s1+s2+s3)/2)/(2.d0*sr2**(1+t1+t2+t3))
        do 6 i3=i31,i32,2
          l3=i3-1
          if(n3.gt.l3) goto 6
          t=0.d0
c      if(n1+n2.eq.-n3) t=t+f102(fac,l1,l2,l3)
c      if(n1+n2.eq.n3)  t=t+f100(fac,l1,l2,l3,n1,n2,n3)*fs(n3+s3)
c      if(n1-n2.eq.-n3) t=t+f100(fac,l1,l2,l3,n1,-n2,-n3)*fs(n2+s2)
c      if(n1-n2.eq.n3)  t=t+f100(fac,l1,l2,l3,-n1,n2,-n3)*fs(n1+s1)
          ic=ic+1
          lnjcg=max(ic,lnjcg)
c      c(ic)=q1*t*f102(fac,l1,l2,l3)/(srpi*dsqrt(dfloat(2*l3+1)))
c      js(ic)=l3*(l3+1)+m3+1
  6     continue
        nm3=nm3-1
        m3=mb
        if(nm3.gt.0) goto 5
  1   continue
 11   continue
 111  continue
 1111 continue
c      cindx(lmindx+1)=ic+1
      lnxcg=max(lnxcg,lmindx+1)
      return
      end

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
     .fac(j1+j2+j3+2)
      t = dsqrt(t*fac(j1+m1+1)*fac(j1-m1+1)*fac(j2+m2+1)*fac(j2-m2+1)*
     .fac(j3+m3+1)*fac(j3-m3+1))
      n1 = max0(j2-j3-m1,j1-j3+m2,0) + 1
      n2 = min0(j1+j2-j3,j1-m1,j2+m2) + 1
      if (n1 .gt. n2) goto 2
      t1 = 0.d0
      do  1  m = n1, n2
        n = m-1
        t1 = t1 +
     .  dble(1+4*(n/2)-2*n)/(fac(m)*fac(j1+j2-j3-n+1)*fac(j1-m1-n+1)*
     .  fac(j2+m2-n+1)*fac(j3-j2+m1+n+1)*fac(j3-j1-m2+n+1))
    1 continue
      f100 = t*t1
      return
    2 f100 = 0.d0
      return
      end

      double precision function f102(fac,l1,l2,l3)
C-
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer l1,l2,l3
      double precision fac(50)
C Local parameters
      integer lt,p,x

      lt = l1 + l2 + l3
      p = lt/2
      if (2*p .ne. lt) goto 1
      f102 = dsqrt(dble(2*l3+1)/dble(lt+1))
      f102 = f102*fac(p+1)/dsqrt(fac(2*p+1))
      x = p-l1
      f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
      x = p-l2
      f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
      x = p-l3
      f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
      if (x .gt. 2*(x/2)) f102 = -f102
      return
    1 f102 = 0.d0
      end

      
