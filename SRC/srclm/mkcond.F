      subroutine mkcond(mode,eval,qp,iq1,iq2,iq3,iq4,cvec,wt)
C- Make conductivity or other weight in one tetrahedron
C-----------------------------------------------------------------------
Ci Input
Ci   mode  :selects which weight is to be estimated; see Outputs
Ci   eval  :eigenvalues at four corners of tetrahedron
Ci   qp    :k-points
Ci   iq1   :index to first corner
Ci   iq2   :index to second corner
Ci   iq3   :index to third corner
Ci   iq4   :index to fourth corner
Ci   cvec  :vector into which to project grad_k E (mode=1)
Ci         :vectors 1,2 defining grad_1 E grad_2 E  (mode=2)
Co Outputs
Co   wt    :some function of eval is interpolated from the four corners
Co         :of the tetrahedron, depending on mode.
Co         :  mode  wt
Co         :   0 : 1 (dos)
Co         :   1 : 1/2 | grad_k E(k) . cvec | (ballistic conductivity)
Co         :       see PRL 74, 586 (1995)).
Co         :   2 : grad_1 E(k) . grad_2 E(k)
Cr Remarks
Cr   For calculating gradient, use:
Cr     dE(1,2,3) = eval(2,2,3)-eval(1)
Cr     dk(1,2,3) = qp(2,3,4) - qp(1)
Cr   The gradient is obtained by the inversion of a:
Cr
Cr       (dE1)                   (dk1x  dk1y  dk1z) (dE/dx)
Cr       (dE2)  = a . grad E  =  (dk2x  dk2y  dk2z) (dE/dy)
Cr       (dE3)                   (dk3x  dk3y  dk3z) (dE/dz)
Cr
Cr   Test:  compare to analytic results for energy bands E(k) = k^2.
Cr     Mode 0:
Cr        N(E) = int d^3k delta(E(k)-E)
Cr             = 4 pi int k^2 delta(k-k0)/2k0
Cr             = 2 pi k0 = 2 pi E^1/2
Cr     Mode 1 (vec=vec2=0,0,1) :
Cr       v1(E) = int d^3k delta(E(k)-E) 1/2 (2k.cvec)
Cr             = 2 pi int k^2 delta(k-k0)/2k0 int 2k |cos| dcos
Cr             = 2 pi k0^2 1/2 2/2 = 1/2 E^1/2 N(E)
Cr     Mode 2 (vec=vec2=0,0,1) :
Cr       v2(E) = int d^3k delta(E(k)-E) (2k.cvec)^2
Cr             = 2 pi int k^2 delta(k-k0)/2k0 int 4k^2 cos^2 dcos
Cr             = 2 pi 2 k0^3 2/3 = 4/3 E N(E)
Cu Updates
Cu   01 Feb 01  Added mode 2
C-----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,iq1,iq2,iq3,iq4
      double precision eval(4),qp(3,1),cvec(3,2),wt
C ... Local parameters
      integer i
      double precision a(3,3),b(3),grad(3),ainv(3,3),det,dabs,ddot
C     Special purpose
C     double precision qbar(3),plat(3,3),qlat(3,3)
C     data plat /1d0,0d0,0d0,-.5d0,0.8660254d0,0d0,0d0,0d0,1.142d0/
C     data qlat /1d0,0.57735027d0,0d0,0d0,1.15470054d0,0d0,
C    .            0d0,0d0,0.87565674d0/
C ... External calls
      external dinv33,dmpy,rx1

      if (mode .eq. 0) then
        wt = 1
      elseif (mode .eq. 1 .or. mode .eq. 2) then
        do  i = 1, 3
          a(1,i) = qp(i,iq2) - qp(i,iq1)
          a(2,i) = qp(i,iq3) - qp(i,iq1)
          a(3,i) = qp(i,iq4) - qp(i,iq1)
          b(i) = eval(i+1) - eval(1)
        enddo
        call dinv33(a,0,ainv,det)
        call dmpy(ainv,3,1,b,3,1,grad,3,1,3,1,3)
        if (mode .eq. 1) then
          wt = 0.5d0*dabs(ddot(3,grad,1,cvec,1))
        else
          wt = ddot(3,grad,1,cvec,1)*ddot(3,grad,1,cvec(1,2),1)
        endif
      else
        call rx1('mkcond: mode %i not implemented',mode)
      endif

C     Special-purpose
C      call dpzero(qbar,3)
C      call dpadd(qbar,qp(1,iq1),1,3,.25d0)
C      call dpadd(qbar,qp(1,iq2),1,3,.25d0)
C      call dpadd(qbar,qp(1,iq3),1,3,.25d0)
C      call dpadd(qbar,qp(1,iq4),1,3,.25d0)
C      call dpcopy(qp(1,iq1),qbar,1,3,1d0)
C      call shorbz(qbar,qbar,qlat,plat)
C      if (sqrt(ddot(2,qbar,1,qbar,1)) .gt. .25d0) wt=0

      end

