      subroutine hsmml(p,rsm,e,lmax,hsm,hsmp)
C- Solid smooth hankel function and energy derivative
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   e     :energy of smoothed Hankel
Ci   lmax  :l-cutoff for hsm
Co Outputs
Co   hsm   :Bloch-summed smoothed Hankels
Co   hsmp  :Energy derivatives of hsm
Cl Local variables
Cl         :
Cr Remarks
Cr   Correctly handles the p=0 limit, since yl(ilm>1) = 0
Cu Updates
Cu   24 Jan 07 New routine
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmax
      double precision p(3),rsm,e,yl(121)
      double precision hsm(1),hsmp(1)
C ... Local parameters
      integer ilm,l,m
      double precision r2,r,rl
      double precision wk(-1:20),wkd(0:20),wk2(0:20)
C     double precision eps,wd1,wd2,wd
C ... External calls
      external hansmr,ropyln

C  ... Debugging check: xidot by brute force
C      call ropyln(1,p(1),p(2),p(3),lmax,1,yl,r2)
C      r = dsqrt(r2)
C      call hansmr(r,e,1/rsm,wk(0),lmax)
C      eps = .001d0
C      call hansmr(dsqrt(r2),e+2*eps,1/rsm,wkd(1),0)
C      call hansmr(dsqrt(r2),e+1*eps,1/rsm,wkd(2),0)
C      call hansmr(dsqrt(r2),e-1*eps,1/rsm,wkd(3),0)
C      call hansmr(dsqrt(r2),e-2*eps,1/rsm,wkd(4),0)
C      wd1 = (wkd(1)-wkd(4))/4/eps
C      wd2 = (wkd(2)-wkd(3))/2/eps
C      wd  = (4*wd2-wd1)/3
C      wk(-1) = 2*wd
C      print '(4f20.15)', wk(-1:2)

C ... Make YL*r^l
      call ropyln(1,p(1),p(2),p(3),lmax,1,yl,r2)
      r = dsqrt(r2)
C ... Radial part of smooth Hankel, and energy derivative
      call hansmd(10,r,e,rsm,lmax,wk(0),wk2,wk2,wkd(0),wk2,wk2)
      wk(-1) = wkd(0)*2
C ... Convert to xi(l) = h(l)/r^l
      if (r .ne. 0) then
        rl = r
        do  l = 1, lmax
          wk(l) = wk(l)/rl
          rl = rl*r
        enddo
      endif
C     print '(4f20.15)', wk(-1:2)

C ... Solid Hankel is product of radial part and YL*r^l
      ilm = 0
      do  l = 0, lmax
        do  m  =  -l, l
          ilm  =  ilm+1
          hsm(ilm)  = yl(ilm)*wk(l)
          hsmp(ilm) = yl(ilm)*wk(l-1)/2
        enddo
      enddo
      end

