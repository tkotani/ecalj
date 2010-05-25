      subroutine hansr4(rsq,lmax,nrx,nr,e,rsm,wk,wk2,chi)
C- Difference between true,smoothed hankels for l=-1...lmax
C  for e nonzero.  Vectorizes for negative e.
C ---------------------------------------------------------------
Ci Inputs
Ci   rsq,nr:vector of points r**2, and number of points.
Ci          Only the first value of rsq may be zero (See Remarks)
Ci   lmax:  highest l for which to evaluate xi.
Ci   e,rsm: smoothing radius and energy
Ci   nrx:   leading dimension of chi
Ci   wk:    array containing y0*dexp(-(r/rsm)**2)
Ci   wk2:   a work array of length nr.
Co Outputs:
Co   chi:  Difference between smoothed and true Hankel for l=-1..lmax.
Cr Remarks
Cr   Smooth Hankels for e<=0 defined in J. Math. Phys. 39, 3393 (1998).
Cr   Notes on conventions:
Cr  *  akap^2 = -e with akap>0
Cr   JMP defines akap in contradistinction to usual convention for kappa:
Cr     kappa = sqrt(e), Im(kappa) >= 0.
Cr   It is related to kappa (defined according to usual conventions) as
Cr     akap = -i kappa (real and positive for e<0)
Cr  *JMP defines
Cr      u(+/-) = exp(-/+ akap*r) erfc(akap*rsm/2 -/+ r/rsm)
Cr             = exp(+/- i kappa r) erfc(-i*kappa*rsm/2 -/+ r/rsm)
Cr   Note: old codes (e.g. rcnsl,hsmbld) define
Cr      uplus = u-  and umins = u+
Cr   The smoothed Hankels for l=0,-1 are:
Cr      h^s_0 (r) = (u+ - u-) / 2r    = (umins - uplus) / 2r
Cr      h^s_-1(r) = (u+ + u-) / 2akap = (umins + uplus) / 2akap
Cr
Cr   We can a new set of functions
Cr      U+/- = 1/2 exp(+/- i kappa r) erfc(r/rsm +/- i*kappa*rsm/2)
Cr   It is convenient to use
Cr      erfc(-x^*) = 2-erfc^*(x).
Cr
Cr  *For e<0, i*kappa is real and U+/- are also real:
Cr      U+   = exp( i kappa r ) - u+/2
Cr      U-   = u-/2
Cr   The difference in smoothed, unsmoothed Hankels is
Cr      h_0  - h^s_0  = (exp(-akap r) - umins/2 + uplus/2) /r
Cr                    = (exp(i kappa r ) - u+/2 + u-/2) /r
Cr                    = (U+ + U-) /r
Cr      h_-1 - h^s_-1 = (exp(-akap r) - umins/2 - uplus/2) /akap
Cr                    = (exp(i kappa r ) - u+/2 - u-/2) /akap
Cr                    = (U+ - U-) /(-i kappa)
Cr
Cr  *For e>0, kappa is real; thus U+ = U-*  since
Cr     erfc(x*) = erfc*(x)
Cr   Keeping the same conventions in h-h^s for e>0 we get
Cr      h_0  - h^s_0  = (U+ + U-) /r = 2 Re (U+) /r
Cr      h_-1 - h^s_-1 = (U+ - U-) / (-i kappa) = -2 Im (U+) / kappa
Cr   The differences are thus real.
Cr
Cr  *The first point may have r->0.  In this case, chi(1,l)=0 except
Cr   chi(1,0), chi(1,-1) which are returned as the (-) value of the
Cr   smoothed hankels, i.e. the infinite unsmoothed part is discarded.
Cr
Cr   In the limit r->0
Cr     h_0  -> 1/r + i*kappa
C
Cr     U+/- -> 1/2 erfc(+/- i*kappa*rsm/2) = 1 - erfc(i*kappa*rsm/2)
Cr     U+   -> 1/2 erfc(i*kappa*rsm/2)
Cr     U-   -> 1/2 erfc(-i*kappa*rsm/2)
Cr          =  1 - 1/2 erfc(i*kappa*rsm/2) = 1 - U+
Cr     U+ + U- -> 1 + terms of order r . Calculate:
Cr     d/dr (erfc(r/rsm+/-i*kappa*rsm/2)) -> 2/srpi/rsm *
Cr                                           exp(-(i*kappa*rsm/2)^2)
Cr     dU+/- /dr |r=0 -> +/- i*kappa U+/-(r->0)
Cr                        - 1/srpi/rsm * exp(-(i*kappa*rsm/2)^2)
Cr     U+ + U- -> 1 + r*i*kappa( U+(0) - U-(0))
Cr                  - r/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
Cr              = 1 + r*i*kappa( 2*U+(0) - 1)
Cr                  - r/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
Cr     h^s_0  = h0 - U+ - U-
Cr           -> 1/r + i*kappa - (1/r + i*kappa( 2*U+(0) - 1))
Cr                  + 1/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
Cr            = i*kappa*erfc(i*kappa*rsm/2) +
Cr              1/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
Cr     h_-1  = exp(i kappa r) /(-i kappa) -> 1/(-i kappa)
Cr     h^s_-1 = 1/(-i kappa) - (U+(0) - U-(0))/(-i kappa)
Cr            = 1/(-i kappa) - (2*U+(0) - 1)/(-i kappa)
Cr            = 1/(-i kappa) - (1 - 2*U-(0))/(-i kappa)
Cr            = 2*U-(0)/(-i kappa)
Cr   For e<0   i*kappa = -akap (real)
Cr     h^s_0  = akap*erfc(akap*rsm/2) + 1/srpi/rsm*exp(-(akap*rsm/2)^2))
Cr     h^s_-1 = 2*U-(0) = erfc(akap*rsm/2)
Cr   For e>0   Use the expressions above
Cr             Warning! hs is NOT REAL ... returns only real part
Cr             h^s_-1 is discontinuous there:
Cr             As e->0, for e>0, 2*U-(0)/(-i kappa) -> i/kappa
Cr             As e->0, for e<0, 2*U-(0)/(-i kappa) -> 1/akap
Cr
Cr  Limiting behavior e<0:
Cr    If both (akap*rsm/2 -/+ r/rsm) >> 1, we have
Cr   -log u(+/-) -> (akap*rsm/2 -/+ r/rsm)^2 -/+ akap*r
Cr               =  (akap*rsm/2)^2 + (r/rsm)^2 >> 1
Cr    u(+/-)     -> exp[-(akap*rsm/2)^2 - (r/rsm)^2]
Cr                  is negligible compared to 1
Cr  Also, if akap*r >> 1,
Cr    h_0 -> exp(-akap*r)/r > h_s^0  is negligible compared to 1
Cu Updates
Cu  15 Aug 00 extended to e>0 (but erfc not vectorized)
C ---------------------------------------------------------------
C     implicit none
      integer nrx,nr,lmax
      double precision rsq(1),e,chi(nrx,-1:0),rsm,wk(nr),wk2(nr)
C Local variables
      logical lpos
      integer l,ir,ir1
      double precision sre,r2,xx,ra,h0,arsm,earsm,kappa,ta
      double precision akap,a,r,um,up,x,facgl
      double complex zikap,zerfc,eikr,zuplus

C ... erfc(x) is evaluated inline as a ratio of polynomials,
C     to a relative precision of <10^-15 for x<5.
C     Different polynomials are used for x<1.3 and x>1.3.
C     Numerators and denominators are t,b respectively.
      double precision w,f1,f2,
     .t10,t11,t12,t13,t14,t15,t16,t17,b11,b12,b13,b14,b15,b16,b17,b18,
     .t20,t21,t22,t23,t24,t25,t26,t27,b21,b22,b23,b24,b25,b26,b27,b28
      parameter (
     .t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0,
     .t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0,
     .t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0,
     .t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0,
     .t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0,
     .t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0,
     .t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0,
     .t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0,
     .b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0,
     .b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0,
     .b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0,
     .b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0,
     .b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0,
     .b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0,
     .b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0,
     .b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0)
C ... f1(w=x-1/2) is erfc(x) for x<1.3, if xx is y0*dexp(-x*x)
      f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)*
     .w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)*
     .w+b13)*w+b12)*w+b11)*w+1)
C ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
      f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)*
     .w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)*
     .w+b23)*w+b22)*w+b21)*w+1)

C --- Setup ---
      if (lmax .lt. 0 .or. nr .eq. 0) return
      a = 1/rsm
      lpos = e .gt. 0d0
      if (lpos) then
        kappa  =  dsqrt(e)
        zikap  = dcmplx(0d0,1d0)*kappa
        ta = 2*a
C       facgl = 4*a*exp(e/(2*a)**2)
        facgl = 4*a*exp(e/ta**2)
      else
        akap = dsqrt(-e)
        arsm = akap*rsm/2
        earsm = dexp(-arsm**2)/2
C       facgl = 4*a*exp(e/(2*a)**2)
        facgl = 8*a*earsm
      endif

C --- chi(*,-1), chi(*,0) = true - sm hankel for each r ---
C ... chi(r,-1) = (h0 - um - up)*r/akap   and
C     chi(r,0)  = (h0 - um + up)          with
C     h0 = exp(-akap*r)/r
C     up = erfc(akap*rsm/2+r/rsm)*exp(akap*r)/(2r)
C     um = erfc(akap*rsm/2-r/rsm)/exp(akap*r)/(2r)
C     um,up correspond to umins/(2r),uplus/(2r) in hansmr
      ir1 = 1
C ... See Remarks for derivation of expressions in this case
      if (rsq(1) .lt. 1d-12) then
        ir1 = 2
        if (lpos) then
          chi(1,0)  = zerfc(zikap/ta)*zikap-facgl/dsqrt(16d0*datan(1d0))
          chi(1,-1) = -dimag(zerfc(zikap/ta))/kappa
        else
C     ... make h0 = erfc(arsm) = erfc(akap*rsm/2)
          xx = earsm/dsqrt(4d0*datan(1d0))
          if (arsm .gt. 1.3d0) then
            h0 = f2(arsm-2d0)
          else
            h0 = f1(arsm-.5d0)
          endif
C         chi(-1) -> erfc(akap/2a)/akap for r=0
          chi(1,0)  = akap*h0 - 4*a*xx
          chi(1,-1) = -h0/akap
        endif
      endif

C ... Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
      if (lpos) then
        do  22  ir = ir1, nr
          r2 = rsq(ir)
          r = dsqrt(r2)
          ra = r*a
          eikr = exp(zikap*r)
C         h_0  - h^s_0  =  Re (U+) /r;  h_-1 - h^s_-1 = -Im (U+) / kappa
          zuplus  = zerfc(r*a + zikap/ta)*eikr

          chi(ir,0)  =  dble(zuplus)/r
          chi(ir,-1) = -dimag(zuplus)/kappa
          wk2(ir) = facgl*wk(ir)
   22   continue

      else
        do  20  ir = ir1, nr
          r2 = rsq(ir)
          r = dsqrt(r2)
          ra = r*a
          sre = akap*r
          h0 = dexp(-sre)/r
          xx = earsm*wk(ir)/r
C     ... Evaluate um; use (akap/2a-ra)^2 = (akap/2a)^2+(ra)^2-akap*r
          x = ra - arsm
          if (x .gt. 1.3d0) then
            um = h0 - f2(x-2d0)
          elseif (x .gt. 0) then
            um = h0 - f1(x-.5d0)
          elseif (x .gt. -1.3d0) then
            um = f1(-x-.5d0)
          else
            um = f2(-x-2d0)
          endif

C     ... Evaluation of up assumes x gt 0
          x = ra + arsm
          if (x .gt. 1.3d0) then
            up = f2(x-2d0)
          else
            up = f1(x-.5d0)
          endif

C     ... Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
          chi(ir,-1) = (h0 - um - up)*r/akap
          chi(ir,0)  =  h0 - um + up
          wk2(ir) = facgl*wk(ir)

   20   continue
      endif

C --- chi(ir,l) for l>1 by upward recursion ---
      facgl = 2*a**2
      do  31  l = 1, lmax
        chi(1,l) = 0
   31 continue
      do  30  l = 1, lmax
        xx = 2*l-1
        do  7  ir = ir1, nr
          chi(ir,l) = (xx*chi(ir,l-1) - e*chi(ir,l-2) + wk2(ir))/rsq(ir)
          wk2(ir) = facgl*wk2(ir)
    7   continue
   30 continue

      end

