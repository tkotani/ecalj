      subroutine hansr4(rsq,lmax,nrx,nr,e,rsm,wk,wk2,chi)
C- Difference between true,smoothed hankels for l=-1...lmax, negative e.
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
Cr   The first point may have r->0.  In this case, chi(1,l)=0 except
Cr   chi(1,0), chi(1,-1) which are returned as the (-) value of the
Cr   smoothed hankels, i.e. the infinite unsmoothed part is discarded.
C ---------------------------------------------------------------
C     implicit none
      integer nrx,nr,lmax
c      double precision rsq(1),e,chi(nrx,-1:0),rsm,wk(nr),wk2(nr)
      double precision rsq(1),e,chi(nrx,-1:lmax),rsm,wk(nr),wk2(nr)
C Local variables
      integer l,ir,ir1
      double precision sre,r2,xx,ra,h0,arsm,earsm
      double precision akap,a,r,um,up,x,facgl
C ... erfc(x) is evaluated as a ratio of polynomials,
C     to a relative precision of <10^-15 for x<5.
C     Different polynomials are used for x<1.3 and x>1.3.
C     Numerators and denominators are t,b respectively.
      double precision w,f1,f2,
     .  t10,t11,t12,t13,t14,t15,t16,t17,b11,b12,b13,b14,b15,b16,b17,b18,
     .  t20,t21,t22,t23,t24,t25,t26,t27,b21,b22,b23,b24,b25,b26,b27,b28
      parameter (
     .  t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0,
     .  t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0,
     .  t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0,
     .  t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0,
     .  t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0,
     .  t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0,
     .  t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0,
     .  t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0,
     .  b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0,
     .  b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0,
     .  b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0,
     .  b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0,
     .  b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0,
     .  b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0,
     .  b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0,
     .  b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0)
C ... f1(w=x-1/2) is erfc(x) for x<1.3, if xx is y0*dexp(-x*x)
      f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)*
     .  w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)*
     .  w+b13)*w+b12)*w+b11)*w+1)
C ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
      f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)*
     .  w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)*
     .  w+b23)*w+b22)*w+b21)*w+1)

C --- Setup ---
      if (lmax .lt. 0 .or. nr .eq. 0) return
      a = 1/rsm
      akap = dsqrt(-e)
      arsm = akap*rsm/2
      earsm = dexp(-arsm**2)/2
      facgl = 8*a*earsm

C --- chi(*,-1), chi(*,0) ---
C ... chi(r,-1) = (um+up)*r/akap and chi(r,0) = (um-up)/(2r) with
C     up = erfc(akap*rsm/2+r/rsm)*exp(akap*r)/(2r)
C     um = erfc(akap*rsm/2-r/rsm)/exp(akap*r)/(2r)
C     um,up correspond to umins/(2r),uplus/(2r) in hansmr
      ir1 = 1
C ... chi1(-1,rsm,r->0) = erfc(arsm)/akap
C ... chi1( 0,rsm,r->0) = 4*y0/rsm*exp(-arsm**2) - akap*erfc(arsm)
      if (rsq(1) .lt. 1d-6) then
        ir1 = 2
C   ... make h0 = erfc(arsm) = erfc(akap*rsm/2)
        xx = earsm/dsqrt(4d0*datan(1d0))
        if (arsm .gt. 1.3d0) then
          h0 = f2(arsm-2d0)
        else
          h0 = f1(arsm-.5d0)
        endif
C   ... chi(-1) -> erfc(akap/2a)/akap for r=0
        chi(1,0)  = akap*h0 - 4*a*xx
        chi(1,-1) = -h0/akap
      endif
      do  20  ir = ir1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        ra = r*a
        sre = akap*r
        h0 = dexp(-sre)/r
        xx = earsm*wk(ir)/r
C ...   Evaluate um
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
C ...   Evaluation of up assumes x gt 0
        x = ra + arsm
        if (x .gt. 1.3d0) then
          up = f2(x-2d0)
        else
          up = f1(x-.5d0)
        endif

C   ... Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
        chi(ir,-1) = (h0 - um - up)*r/akap
        chi(ir,0)  =  h0 - um + up
        wk2(ir) = facgl*wk(ir)

   20 continue

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
