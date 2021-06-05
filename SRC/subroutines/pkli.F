      subroutine qmpkl(rmax,rsm,k0,kmax,lmax,pmax,aklm,qmp,qval)
C- Integrals of pkl r^(2+2*l+p) from 0 to rmax, and value of pkl r^(l+p)
C  at rmax, p=0..pmax.
C     implicit none
      integer k0,kmax,lmax,pmax,l,p
      double precision rmax,rsm,aklm(0:k0,0:k0,0:lmax),
     .qmp(0:k0,0:lmax,0:pmax),qval(0:k0,0:lmax,0:pmax),
     .sumv,sumq,rq,rv,qq
      integer m,k

C --- Decomposition of pkl into simple polynomials ---
      call pkl2r(rsm,k0,kmax,lmax,aklm)

C --- Integrate the pkl term by term ---
      do  101  p = 0, pmax
      do  10  l = 0, lmax
        do  20  k = 0, kmax
          sumq = 0
          sumv = 0
          qq = 2*l+p+3-2
          rq = rmax**qq*rsm**2/2
          rv = rq/rmax**(l+3)
          do  30  m = 0, k
            qq = qq+2
            rq = rq*2*(rmax/rsm)**2
            rv = rv*2*(rmax/rsm)**2
            sumq = sumq + rq/qq*aklm(m,k,l)
            sumv = sumv + rv*aklm(m,k,l)
   30     continue
          qmp(k,l,p) = sumq
          qval(k,l,p) = sumv
   20   continue
   10 continue
 101  continue
      end

      subroutine spkli(rmax,rsm,k0,kmax,lmax,pmax,aklm,s)
C- Integrals of pkl pk'l r^(2+2l+p) from 0 to rmax, p=0..pmax
C     implicit none
      integer k0,kmax,lmax,pmax,l,p
      double precision rmax,rsm,aklm(0:k0,0:k0,0:lmax),
     .s(0:k0,0:k0,0:lmax,0:pmax),sum,qq
      integer m1,m2,k1,k2

C --- Decomposition of pkl into simple polynomials ---
      call pkl2r(rsm,k0,kmax,lmax,aklm)

C --- Integrate the pkl pk'l product term by term ---
      do  101  p = 0, pmax
      do  10  l = 0, lmax
        do  201  k1 = 0, kmax
        do  20  k2 = 0, kmax
          sum = 0d0
          do  301  m1 = 0, k1
          do  30  m2 = 0, k2
            qq = 2*(m1+m2+l)+p+3
            sum = sum + (2d0/rsm**2)**(m1+m2)/qq*rmax**qq*
     .            aklm(m1,k1,l)*aklm(m2,k2,l)
   30     continue
 301      continue
          s(k1,k2,l,p) = sum
   20   continue
 201    continue
   10 continue
 101  continue
      end
      subroutine tpkli(rmax,rsm,k0,kmax,lmax,aklm,s)
C- Integrals of pkl (-nabla^2 pk'l) r^(2+2l) from 0 to rmax
C  Integrates (pkl 1/r d^2(r*pk'l) - l(l+1)/r^2*pkl*pk'l) r^(2+2l)
C  from r=0 to rmax
C     implicit none
      integer k0,kmax,lmax,l
      double precision rmax,rsm,aklm(0:k0,0:k0,0:lmax),
     .s(0:k0,0:k0,0:lmax),sum,qq
      integer m1,m2,k1,k2

C --- Decomposition of pkl into simple polynomials ---
      call pkl2r(rsm,k0,kmax,lmax,aklm)

C --- Integrate -(-l(l+1)r**-2 pkl pk'l) r^2 ---
      do  10  l = 0, lmax
        do  121  k1 = 0, kmax
        do  12  k2 = 0, kmax
          sum = 0d0
          do  141  m1 = 0, k1
          do  14  m2 = 0, k2
            qq = 2*(m1+m2+l)+1
            sum = sum + l*(l+1)*(2d0/rsm**2)**(m1+m2)/qq*rmax**qq*
     .          aklm(m1,k1,l)*aklm(m2,k2,l)
   14     continue
 141      continue
          s(k1,k2,l) = sum
   12   continue
 121    continue
   10 continue

C --- Integrate pkl -1/r d^2(r*pk'l)/dr^2 r^2 ---
      do  20  l = 0, lmax
        do  221  k1 = 0, kmax
        do  22  k2 = 0, kmax
          sum = 0d0
          do  241  m1 = 0, k1
          do  24  m2 = 0, k2
            qq = 2*(m1+m2+l)+1
            sum = sum - (2d0/rsm**2)**(m1+m2)/qq*rmax**qq*
     .          aklm(m1,k1,l)*aklm(m2,k2,l)*(2*m2+l+1)*(2*m2+l)
   24     continue
 241      continue
          s(k1,k2,l) = s(k1,k2,l) + sum
   22   continue
 221    continue
   20 continue

      end
      subroutine s3pkl(rmax,rsma,k0,l0,km1,km2,km3,lm1,lm2,lm3,aklm,s)
C- Integrals to rmax of products of three radial p_kl's
C  Returns s(k1l1,k2l2,k3l3) where 0<k_j<km_j and 0<l_j<lm_j
C     implicit none
      integer k0,l0,km1,km2,km3,lm1,lm2,lm3
      double precision rmax,rsma,aklm(0:k0,0:k0,0:l0),
     .s(0:k0,0:k0,0:k0,0:l0,0:l0,0:l0),sum,qq
      integer m1,m2,m3,k1,k2,k3,l1,l2,l3
      double precision tar2

C --- Decomposition of pkl into simple polynomials ---
      if (max(km1,km2,km3) .gt. k0) call rx('s3pkl: kmax gt k0')
      if (max(lm1,lm2,lm3) .gt. l0) call rx('s3pkl: lmax gt l0')
      call pkl2r(rsma,k0,max(km1,km2,km3),max(lm1,lm2,lm3),aklm)

C --- Integrate the pkl products term by term ---
      tar2 = 2d0/rsma**2
      do  12  l1 = 0, lm1
      do  11  l2 = 0, lm2
      do  10  l3 = 0, lm3
        do  22  k1 = 0, km1
        do  21  k2 = 0, km2
        do  20  k3 = 0, km3
          sum = 0d0
          do  32  m1 = 0, k1
          do  31  m2 = 0, k2
          do  30  m3 = 0, k3
            qq = 2*(m1+m2+m3)+l1+l2+l3+3
            sum = sum + tar2**(m1+m2+m3)/qq*rmax**qq*
     .                  aklm(m1,k1,l1)*aklm(m2,k2,l2)*aklm(m3,k3,l3)
   30     continue
   31     continue
   32     continue
          s(k1,k2,k3,l1,l2,l3) = sum
   20   continue
   21   continue
   22   continue
   10 continue
   11 continue
   12 continue
      end
      subroutine s3pkl0(rmax,rsm,k0,kmax,kmax3,lmax,aklm,s)
C- Integrals from r=0 to rmax of p_k1,l p_k2,l p_k3,0
C  Integrals are: pkl pk'l p_k3,0 r^(2+2l), 0<k3<kmax3
C     implicit none
      integer k0,kmax,kmax3,lmax,l
      double precision rmax,rsm,aklm(0:k0,0:k0,0:lmax),
     .s(0:k0,0:k0,0:k0,0:*),sum,qq
      integer m1,m2,m3,k1,k2,k3

C --- Decomposition of pkl into simple polynomials ---
      if (max(kmax,kmax3) .gt. k0) call rx('s3pkl: kmax gt k0')
      call pkl2r(rsm,k0,max(kmax,kmax3),lmax,aklm)

C --- Integrate the pkl pk'l product term by term ---
      do  10  l = 0, lmax
        do  22  k1 = 0, kmax
        do  21  k2 = 0, kmax
        do  20  k3 = 0, kmax3
          sum = 0d0
          do  32  m1 = 0, k1
          do  31  m2 = 0, k2
          do  30  m3 = 0, k3
            qq = 2*(m1+m2+m3+l)+3
            sum = sum + (2d0/rsm**2)**(m1+m2+m3)/qq*rmax**qq*
     .              aklm(m1,k1,l)*aklm(m2,k2,l)*aklm(m3,k3,0)
   30     continue
   31     continue
   32     continue
          s(k1,k2,k3,l) = sum
   20   continue
   21   continue
   22   continue
   10 continue
      end
      subroutine s3pklx(rmax,rsm,rsm3,k0,kmax,kmax3,lmax,aklm,aklm3,s)
C- Integrals from r=0 to rmax of p_k1,l p_k2,l p_k3,0
C  Integrals are: pkl pk'l p_k3,0 r^(2+2l), 0<k3<kmax3
C  Pkl for k=1,2 defined by rsm, for k=3 defined by rsm3
C     implicit none
      integer k0,kmax,kmax3,lmax
      double precision rmax,rsm,rsm3,aklm(0:kmax,0:kmax,0:lmax),
     .s(0:k0,0:k0,0:k0,0:*),aklm3(0:kmax3,0:kmax3,0:0)
      integer l,k1,k2,k3,m1,m2,m3
      double precision sum,qq
c|    double precision p1(0:20,0:10),p2(0:20,0:10)


C --- Decomposition of pkl into simple polynomials ---
      call pkl2r(rsm, kmax, kmax, lmax, aklm)
      call pkl2r(rsm3,kmax3,kmax3, 0,   aklm3)

C --- Integrate the pkl pk'l product term by term ---
      do  10  l = 0, lmax
        do  22  k1 = 0, kmax
        do  21  k2 = 0, kmax
        do  20  k3 = 0, kmax3
          sum = 0d0
          do  32  m1 = 0, k1
          do  31  m2 = 0, k2
          do  30  m3 = 0, k3
            qq = 2*(m1+m2+m3+l)+3
            sum = sum + (2d0/rsm**2)**(m1+m2)*(2d0/rsm3**2)**m3
     .              /qq*rmax**qq*aklm(m1,k1,l)*aklm(m2,k2,l)*aklm3(m3,k3,0)
   30     continue
   31     continue
   32     continue
          s(k1,k2,k3,l) = sum
   20   continue
   21   continue
   22   continue
   10 continue

c ... do numerically here
c|      n=801
c|      dr=rmax/(n-1)
c|      do l=0,lmax
c|        do k1=0,kmax
c|          do k2=0,kmax
c|            do k3=0,kmax3
c|              sum=0d0
c|              do i=2,n
c|                r=(i-1)*dr
c|                wgt=2*(mod(i+1,2)+1)/3d0
c|                if (i.eq.n) wgt=1d0/3d0
c|                call radpkl(r,rsm,kmax,lmax,20,p1)
c|                call radpkl(r,rsm3,kmax3,0,20,p2)
c|                sum=sum+dr*wgt*p1(k1,l)*p1(k2,l)*p2(k3,0)*r**(2*l+2)
c|              enddo
c|              write (6,890) l,k1,k2,k3,s(k1,k2,k3,l),sum
c|  890         format(4i4,2f14.8)
c|              s(k1,k2,k3,l)=sum
c|            enddo
c|          enddo
c|        enddo
c|      enddo

      end

      subroutine pkl2r(rsm,k0,kmax,lmax,aklm)
C- Coefficients that decompose pkl in simple polynomials
C  pkl = sum_m=0..k akl(m,k,l) x^m, where x=2*(r/rsm)^2
C     implicit none
      integer k0,kmax,lmax,m,k,l
      double precision a,aklm(0:k0,0:k0,0:lmax),tk2lp1,rsm

      a = 1/rsm
      call dpzero(aklm,(k0+1)**2*(lmax+1))

C --- do explicitly for k=0,1 ---
      do  1  l = 0, lmax
        aklm(0,0,l) = a**l
    1 continue
      if (kmax .eq. 0) return
      do  10  l = 0, lmax
        aklm(0,1,l) = -a**l
        aklm(1,1,l) = a**l/(2*l+3)
   10 continue

C --- do by recursion for higher k ---
      do  21  l = 0, lmax
      do  20  k = 2, kmax
        tk2lp1 = 2*k+2*l+1
        do  24  m = 0, k-1
          aklm(m,k,l) = aklm(m,k,l) - 2*(k-1)*aklm(m,k-2,l)/tk2lp1
          aklm(m,k,l) = aklm(m,k,l) - (4*k+2*l-1)*aklm(m,k-1,l)/tk2lp1
          aklm(m+1,k,l) = aklm(m+1,k,l) + aklm(m,k-1,l)/tk2lp1
   24   continue
   20 continue
   21 continue

C      do  30  l = 0, lmax
C      do  30  k = 0, kmax
C   30 print 333, l,k,(aklm(m,k,l), m=0,kmax)
C  333 format(2i4,10f12.6)

      end

      subroutine phkl2r(rsm,k0,kmax,lmax,aklm)
C- Coefficients that decompose phi_kl in simple polynomials
C  phi_kl = sum_m=0..k akl(m,k,l) x^m, where x=2*(r/rsm)^2

C     implicit none
      integer k0,kmax,lmax,m,k,l
      double precision a,ta2,aklm(0:k0,0:k0,0:lmax),tk2lm1,rsm

      a = 1/rsm
      ta2 = 2*a**2
      call dpzero(aklm,(k0+1)**2*(lmax+1))

C --- do explicitly for k=0,1 ---
      do  1  l = 0, lmax
        aklm(0,0,l) = ta2**l
    1 continue
      if (kmax .eq. 0) return
      do  10  l = 0, lmax
        aklm(0,1,l) = -ta2**(l+1)*(2*l+3)
        aklm(1,1,l) = ta2**(l+1)
   10 continue

C --- do by recursion for higher k ---
      do  201  l = 0, lmax
      do  20  k = 2, kmax
        tk2lm1 = 2*k+2*l-1
        do  24  m = 0, k-1
          aklm(m,k,l) = aklm(m,k,l) - 2*ta2**2*(k-1)*tk2lm1*aklm(m,k-2,l)
          aklm(m,k,l) = aklm(m,k,l) - ta2*(4*k+2*l-1)*aklm(m,k-1,l)
          aklm(m+1,k,l) = aklm(m+1,k,l) + ta2*aklm(m,k-1,l)
   24   continue
   20 continue
 201  continue

      end

      subroutine p2pkl(rsm,k0,mmax,lmax,aklm,cmkl)
C- Decomposition of simple polynomials into P_kl
C  Let F_ml = r**(2m+l)Y_L.  Then F_mL = sum_k=0..m c_pkl p_kl r^l Y_L
C     implicit none
      integer k0,mmax,lmax,p,k,l
      double precision aklm(0:k0,0:k0,0:lmax),cmkl(0:k0,0:k0,0:lmax),
     .a,ta2,dfac,pi,srpi,rsm,f(0:200),fack,fac,sum,fpi
      integer m

C --- Decomposition of phi_kl into simple polynomials ---
      call phkl2r(rsm,k0,mmax,lmax,aklm)

      pi = 4*datan(1d0)
      srpi = dsqrt(4*datan(1d0))
      fpi = 16*datan(1d0)
      a = 1/rsm
      ta2 = 2*a**2

C --- Make f(m) = int r^(2m) exp(-a^2r^2) ---
      dfac = 1
      f(0) = srpi/(2*a)
      if (lmax+2*mmax+1 .gt. 200) call rx('p2pkl: f too small')
      do  5  m = 1, (lmax+2*mmax+1)
        f(m) = (2*m-1)*f(m-1)/ta2
    5 continue

C --- Integrate r^(2m+l) Y_L phi_kl r^l Y_L g0 term by term ---
C     and fold in scaling to make pkl expansion
      do  10  m = 0, mmax
        dfac = 1d0
        do  12  l = 0, lmax
          fack = 1d0
          do  20  k = 0, m
            fac = (a**2/pi)**1.5d0/ta2
            sum = 0
            do  31  p = 0, k
              fac = fac*ta2
              sum = sum + fac*aklm(p,k,l)*f(l+m+p+1)
   31       continue
            fac = (4*a*a)**k * a**l * fack * dfac
            cmkl(m,k,l) = sum*(fpi/fac)
            fack = fack*(k+1)
   20     continue
          dfac = dfac*(2*l+3)
   12   continue
   10 continue

C      do  30  l = 0, lmax
C      do  30  m = 0, mmax
C   30 print 333, l,m,(cmkl(m,k,l), k=0,mmax)
C  333 format(2i4,10f12.6)

      end

