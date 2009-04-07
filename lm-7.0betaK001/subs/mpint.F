      subroutine mpint(f1,f2,lmaxs,lmax,lmaxv,nr,nsp,ri,wi,p)
C- Multipole moments of function products f1_l1(r) f2_l2(r) r**l
C ----------------------------------------------------------------------
Ci Inputs
Ci   f1,f2  functions
Ci   lmxs   leading dimensions of p are 0:lmaxs
Ci   lmax   l-cutoff for f1,f2
Ci   lmaxv  l-cutoff for multipole moments
Ci   nr     number of radial mesh points
Ci   ri,wi  radial mesh, and integration weights; see radwgt.
Co Outputs
Co   p
Co   wi     is DESTROYED on output
Cu Updates
Cu   05 Jan 04 leading dimensions of p distinct from lmax
C ----------------------------------------------------------------------
C     implicit none
      integer lmaxs,lmax,lmaxv,nr,nsp
      double precision f1(nr,lmax+1,nsp),f2(nr,lmax+1,nsp),
     .  ri(nr),wi(nr),p(0:lmaxs,0:lmaxs,0:lmaxv,3,nsp)
      integer ir,l1p1,l2p1,lf,isp
      double precision sum

      do  10  lf = 0, lmaxv
C   ... Integral f1 f1
        do  20  l1p1 = 1, lmax+1
        do  20  l2p1 = 1, l1p1
        do  20  isp = 1, nsp
        call xxint0(f1(1,l1p1,isp),f1(1,l2p1,isp),wi,nr,sum)
        p(l1p1-1,l2p1-1,lf,1,isp) = sum
   20   p(l2p1-1,l1p1-1,lf,1,isp) = sum
C   ... Integral f1 f2
        do  21  l1p1 = 1, lmax+1
        do  21  l2p1 = 1, lmax+1
        do  21  isp = 1, nsp
        call xxint0(f1(1,l1p1,isp),f2(1,l2p1,isp),wi,nr,sum)
   21   p(l1p1-1,l2p1-1,lf,2,isp) = sum
C   ... Integral f2 f2
        do  22  l1p1 = 1, lmax+1
        do  22  l2p1 = 1, l1p1
        do  22  isp = 1, nsp
        call xxint0(f2(1,l1p1,isp),f2(1,l2p1,isp),wi,nr,sum)
        p(l1p1-1,l2p1-1,lf,3,isp) = sum
   22   p(l2p1-1,l1p1-1,lf,3,isp) = sum

        do  25  ir = 1, nr
   25   wi(ir) = wi(ir)*ri(ir)
   10 continue

C      call xxint1(p(0,0,0,1,1),lmax,lmaxv)
C      call xxint1(p(0,0,0,2,1),lmax,lmaxv)
C      call xxint1(p(0,0,0,3,1),lmax,lmaxv)
C      call xxint1(p(0,0,0,1,nsp),lmax,lmaxv)
C      call xxint1(p(0,0,0,2,nsp),lmax,lmaxv)
C      call xxint1(p(0,0,0,3,nsp),lmax,lmaxv)


      end

      subroutine xxint0(a1,a2,h,nr,sum)
C     implicit none
      double precision a1(1),a2(1),h(1),sum
      integer ir,nr
      sum = 0
      do  10  ir = 1, nr
   10 sum = sum + a1(ir)*a2(ir)*h(ir)
      end

C      subroutine xxint1(p,lmax,lmaxv)
C      implicit none
C      integer lmax,lmaxv,l1p1,l2p1,lf
C      double precision p(lmax+1,lmax+1,0:lmaxv)
C      print 996, lmax,lmaxv
C  996 format(' xxint1:  lmax=',i1,'   lmaxv=',i1)
C      do  10  l1p1 = 1, lmax+1
C      do  10  l2p1 = 1, lmax+1
C   10 print 200, l1p1-1,l2p1-1,(p(l1p1,l2p1,lf),lf=0,lmaxv)
C  200 format(2i4,6f12.6)
C      end

