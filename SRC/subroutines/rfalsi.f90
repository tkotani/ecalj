subroutine rfalsi(xn,fn,xtol,ftol,dxmn,dxmx,isw,wk,ir)
  !- Find root of a function by a modified regular falsi method.
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i  xn,fn: (point,function) pair; xn is also output; see below
  !i  xtol,ftol: tolerances in x or f.  1s digit of isw determines
  !i         which tolerance is used; see isw, below.
  !i  dxmn:  smallest change in x for which a change in f can be detected.
  !i  dxmx:  maximum allowed extrapolation in x in a single step.  Can
  !i         be either positive or negative, which sets the step direction
  !i         when rfalsi doesn't have enough information to decide internally
  !i         or 1's digit of isw has '4' bit set.
  !i  isw:   compound of one-digit switches.
  !i    Ones digit (imposes contraints on slope at root)
  !i      0: no constraints imposed
  !i      1: require slope >0 at root (pos curvature of antiderivative)
  !i      2: require slope <0 at root (neg curvature of antiderivative)
  !i      4: while root not bracketed, set suggested x = current x + dxmx
  !i    10's digit: controls tolerance for root
  !i      0: convergence when |f| < ftol
  !i      1: convergence when |bracket| < xtol
  !i      2: convergence when |f| < ftol and |bracket| < xtol
  !i      3: convergence when |f| < ftol or |bracket| < xtol
  !i      4  Set this bit if rfalsi is to return, with ir=0 or 4,
  !i         xn corresponding to some xn actually passed and fn closest to
  !i         zero, rather than estimate the root (or minimum).
  !o  Co Outputs:
  !o  xn:    next suggested x for finding root, OR:
  !o         estimate for root if rfalsi returns ir = 0
  !o         estimate for extremum if rfalsi returns ir = 4
  !o  wk:    work array of length 12.  Should be passed through unaltered
  !o         between minimization steps.  wk holds:
  !o         1..3  x from prior calls
  !o         4..6  f from prior calls
  !o         7..9  scaled f
  !o         10,11 smallest and largest x passed to rfalsi
  !o         12    describes how prior points are reordered.
  !o               0 => new x2 = old x1; new x1 = old x0; new x0 = xn.
  !o               If nonzero, the following bits mean:
  !o               1 (bit 1) => permute points x0,x1
  !o               2 (bit 2) => permute points x0,x2
  !o               4 (bit 3) => permute points x1,x2
  ! o Inputs/Outputs
  ! o ir:    program flow control.
  ! o        Input:  To start a new line minimization, set ir to zero.
  ! o        In each iteration rfalsi will set ir, and that value should
  ! o        be preserved while iterating within a line minimization.
  ! o        On exit, rfalsi sets ir to a positive or negative number
  ! o        as described below.
  ! o        Thus if ir<0, rfalsi expects a new (xn,fn) pair.
  ! o        If ir is returned >0, rfalsi has finished.
  ! o        So, you should not call rfalsi with ir>0.
  ! o        On exit, rfalsi sets ir to one of the following,
  ! o        and suggests caller uses the returned value of xn.
  ! o    >=0: no more function calls expected:
  ! o      0: function has converged to within tolerance.
  ! o      1: input fn equals 0.  rfalsi does nothing.
  ! o      2: nonsensical tolerance or isw.  rfalsi does nothing.
  ! o      3: input xn equals x1 or x2.  rfalsi does nothing.
  ! o      4: extremum bracketed with no root.  rfalsi returns
  ! o         an estimate of the minimum point.  NOT IMPLEMENTED
  ! o     <0: expects new function call, suggesting new point xn:
  ! o     -1: no information yet: require new point at suggested xn
  ! o         to improve the estimate of the root.
  ! o     -2: a root is bracketed and rfalsi is attempting linear
  ! o         interpolation to improve the estimate of the root.
  ! o     -3: a root is bracketed and rfalsi is attempting quadratic
  ! o         interpolation to improve the estimate of the root.
  ! o     -4: a root has not been bracketed.
  ! o         rfalsi will suggests xn from a linear extrapolation.
  ! o     -5: a root has not been bracketed.  rfalsi will try to
  ! o         guess a new point xmin+dxmx or xmax+dxmx depending
  ! o         on the sign of dxmx. Here xmin and xmax are the
  ! o         smallest and largest points passed to rfalsi.
  ! o     -6: a constraint was not fufilled.  rfalsi will try to
  ! o         guess a new point either xmin-|dxmx| or xmax+|dxmx|
  ! o         depending on the information it has.
  !r Remarks:
  !r  If a root is bracketed (some pair of f1*f2 < 0), a root is guaranteed
  !r  for an analytic function.  When not bracketed, rfalsi extrapolates,
  !r  either with a linear extrapolation of prior points or by a fixed
  !r  step size  dxmx  until a root is bracketed or something goes wrong.
  !u Updates
  !u   06 Aug 07 Handle interpolation case when xn=~x0 or xn=~x1
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed Parameters
  double precision :: xtol,ftol,dxmn,dxmx,xn,fn,wk(12)
  integer :: ir,isw
  ! Local Variables
  integer :: i1mach,ic,ipr,i,itol,ir0,lqua,ic4,scrwid,stdo
  parameter (scrwid=100)
  character*(scrwid) outs
  double precision :: f2x,xtl,x0,f0,f0x,x1,f1,x2,f2,f1x,xmin,xmax, &
       x(0:2),f(0:2),fx(0:2),dum,dx2,dx1,df2,df1,fpp,fp,den,dirdx,xtoll
  equivalence (x0,x(0)),(x1,x(1)),(x2,x(2))
  equivalence (f0,f(0)),(f1,f(1)),(f2,f(2))
  equivalence (f0x,fx(0)),(f1x,fx(1)),(f2x,fx(2))
  logical :: ltmp,cnvgx,cnvgf,lextr,root1,root2,cst1,cst2
  double precision :: d1mach

  ! ... Recover local variables from work array
  stdo = i1mach(2)
  if (ir == 0) call dpzero(wk,11)
  call dcopy(3,wk(1),1,x,1)
  call dcopy(3,wk(4),1,f,1)
  call dcopy(3,wk(7),1,fx,1)
  xmin = wk(10)
  xmax = wk(11)
  ! ... Point-independent setup
  wk(12) = 0
  ir0 = ir
  call getpr(ipr)
  !     ipr = 55
  ic = mod(mod(isw,10),4)
  ic4 = mod(isw,10) - ic
  itol = mod(isw/10,4)
  if (ic == 2) ic = -1
  xtl = max(xtol,dxmn)
  if (ipr > 40 .AND. ir == 0) &
       write(stdo,*)'rfalsi newstart ic,itol,xtl,ftol,dxmn,dxmx',ic,itol,xtl,ftol,dxmn,dxmx
  !     .call awrit8('  rfalsi: new start%?#n#  (c)##'//
  !     .'%?#n# xtol=%1;3g#%j#%?#n==2# and##%-1j%?#n==3# or##'//
  !     .'%?#n# ftol=%1;3g#%j#  dxmn=%1;3g  dxmx=%1;3g',
  !     .' ',scrwid,stdo,ic,itol,xtl,itol,itol-1,ftol,dxmn,dxmx)
  if (ftol <= 0 .AND. itol == 0 .OR. &
       xtl <= 0 .AND. itol == 1 .OR. &
       (ftol <= 0 .OR. xtl <= 0) .AND. itol == 2 .OR. &
       (ftol <= 0 .AND. xtl <= 0) .AND. itol == 3 .OR. &
       mod(mod(isw,10),4) == 3) then
     ir = 2
     goto 999
  endif

  ! --- Treat next point ---
  !  10 continue
  !      call awrit3('  rfalsi (ir=%i)  xn=%1,8;8d  fn=%1;4g',
  !     .  ' ',scrwid,stdo,ir,xn,fn)

  !     (x1,f1) -> (x2,f2); (x0,f0) -> (x1,f1); (xn,fn) -> (x0,f0)
  lqua = 0
  do  12  i = 1, 0, -1
     fx(i+1) = fx(i)
     f(i+1) = f(i)
     x(i+1) = x(i)
12 enddo
  fx(0) = fn
  f(0) = fn
  x(0) = xn
  if (ir == -1) then
     x2 = x1
     f2 = f1
  endif

  ! ... Special treatments
  if (fn == 0) then
     ir = 1
     call info2(31,0,0, &
          ' rfalsi:  fn=0 at xn=%1,8;8g; return ir=%i',xn,ir)
     goto 999
  elseif (ir == 0) then
     ir = -1
     xmin = xn
     xmax = xn
     xn = x0+dxmx
     goto 998
  elseif (x1 == xn .OR. x2 == xn) then
     ir = 3
     goto 999
  elseif (ir >  0) then
     goto 999
  endif

  xmin = min(xmin,xn)
  xmax = max(xmax,xn)

  ! ... T if root bracketed, with constraint satisfied
  cst1  = iabs(ic) .ne. 1 .or. (f1-fn)*(x1-xn)*ic .gt. 0
  cst2 =  iabs(ic) .ne. 1 .or. (f2-fn)*(x2-xn)*ic .gt. 0
  ! ... if only 1 cst satisfied, x0 must be in the middle
  if ((cst1 .neqv. cst2) .AND. (x0-x1)*(x0-x2) > 0) then
     if (cst1 .AND. abs(x0-x2) < abs(x0-x1)) cst1 = .FALSE. 
     if (cst2 .AND. abs(x0-x2) > abs(x0-x1)) cst2 = .FALSE. 
  endif
  root1 = fn*f1 .lt. 0 .and. cst1
  root2 = fn*f2 .lt. 0 .and. cst2
  ! ... Swap x1,x2 if (x0,x2) a closer fit
  if (( .NOT. root1 .AND. abs(f0) < abs(f1) .OR. root2) &
       .AND. abs(x1-x0) > abs(x2-x0) .OR. &
       .NOT. root1 .AND. root2 .OR. .NOT. cst1 .AND. cst2) then
     dum = x2
     x2 = x1
     x1 = dum
     dum = f2
     f2 = f1
     f1 = dum
     dum = f2x
     f2x = f1x
     f1x = dum
     ltmp = cst1
     cst1 = cst2
     cst2 = ltmp
     wk(12) = 4
  endif

  ! ... Make first and second derivatives from three points, if available
  dx2 = x2-x0
  dx1 = x1-x0
  df2 = f2-f0
  df1 = f1-f0
  den = dx1*(dx1-dx2)*dx2
  fpp = 0
  if (den /= 0) then
     fp = -((-(df2*dx1**2) + df1*dx2**2)/den)
     fpp = (2*(df1*dx2-df2*dx1))/den
  endif

  ! ... T if local extremum bracketed
  lextr = (f0-f1)*(x0-x1)*(f1-f2)*(x1-x2) .lt. 0 .and. fpp*f1 .ge. 0
  if (lextr .AND. .NOT. (root1 .OR. root2) .AND. ipr >= 30) then
     write (stdo,*) 'rfalsi: local extremum bracketed without root:'
     !        call awrit6('  x0,x1,x2 = %1;6g %1;6g %1;6g  f0,f1,f2 = '//
     !     .  '%1;6g %1;6g %1;6g',' ',scrwid,stdo,x0,x1,x2,f0,f1,f2)
     write(stdo,"(a,3d13.5,2x,3d13.5)")' x0,x1,x2: f0,f1,f2 = ',x0,x1,x2,f0,f1,f2
  endif

  ! --- CASE I: root bracketed with possible constraint satisfied ---
  if (root1 .OR. root2) then

     !   ... Linear estimate for xn
     xn = (f1x*x0-f0*x1)/(f1x-f0)

     ! ... If conditions allow, improve estimate for x.  Conditions:
     !     a sufficiently distinct third point exists; f1<=f1x*2;
     !     corrected slope of same sign; new x still within bracket.
     if (dabs(x2-x0) > dxmn .AND. &
          dabs(x2-x1) > dxmn .AND. abs(f1) <= abs(f1x*2)) then

        if (fp*df1/dx1 > 0) then
           !       ... Estimate for quadratic correction
           dum = x0 - f0/fp - (f0**2*fpp)/(2*fp**3)
           if (dum > min(x0,x1) .AND.  dum < max(x0,x1)) then
              lqua = 2
              xn = dum
              !       ... Estimate for linear term, with improved slope
           else
              dum = x0 - f0/fp
              if (dum > min(x0,x1) .AND.  dum < max(x0,x1)) then
                 lqua = 1
                 xn = dum
              endif
           endif
        endif
     endif
     !   ... Artificially reduce f1 to avoid ill-conditioned cases
     if (min(dabs(xn-x0),dabs(xn-x1)) < .1d0*dabs(x0-x1))f1x=f1x/2

     !   ... Ask for another point, or if set ir if converged to within tol
     ir = -2
     if (lqua /= 0) ir = -3
     cnvgx = abs(x1-x0) .le. xtol
     if (wk(12) == 4) cnvgx = abs(x2-x0) < xtol
     cnvgf = abs(fn) .le. ftol
     if (itol == 0 .AND. cnvgf .OR. itol == 1 .AND. cnvgx .OR. &
          itol == 2 .AND. (cnvgf .AND. cnvgx) .OR. &
          itol == 3 .AND. (cnvgf .OR. cnvgx)) ir = 0

     !   ... If xn near machine precision of x0 or x1:
     !       call dswap(1,x0,1,x1,1)
     xtoll = 256*d1mach(3)
     if (abs(xn-x0) <= abs(xn)*xtoll &
          .AND. (ir == -2 .OR. ir == -3)) then
        !          print *, '1',xn-x0, xn-x1 !, sign(x0*xtoll,x1-x0)
        xn = x0 + sign(x0*xtoll,x1-x0)
        if (abs(xn-x0) > 0.9d0*abs(x1-x0)) then
           !            print *, 'branch 1aa ...',xn-x0
           xn = (x0*7+x1)/8
           !            print *, 'branch 1aa ...',xn-x0
           !           pause
        endif
        !          print *, '1a',xn-x0, xn-x1, (xn-x0)*(xn-x1).gt.0
        !         pause
        ir = -2
     elseif (abs(xn-x1) <= abs(xn)*xtoll &
          .AND. (ir == -2 .OR. ir == -3)) then
        !         print *, '2',xn-x1, xn-x0   !, sign(x1*xtoll,x0-x1)
        xn = x1 + sign(x1*xtoll,x0-x1)
        if (abs(xn-x1) > 0.9d0*abs(x0-x1)) then
           !            print *, 'branch 2aa ...',xn-x1
           xn = (x1*7+x0)/8
           !            print *, 'branch 2aa ...',xn-x1
           !           pause
        endif
        !         print *, '2a',xn-x1, xn-x0, (xn-x0)*(xn-x1).gt.0
        ir = -2
     endif
     !       call dswap(1,x0,1,x1,1)

     ! --- CASE III: root not found, step prescribed ---
  elseif (ic4 == 4) then
     goto 80

     ! --- CASE IV: both constraints violated ---
     !     Use whatever info we have to try xmin-|dxmx| or xmax+|dxmx|
     !                                            *
     !            Illustration for ic=1:        -----------------  ic=1
     !            Slope the wrong sign.               *
     !                                                       *
  elseif (ic /= 0 .AND. .NOT. (cst1 .OR. cst2)) then
     !   ... If have second derivative, choose dir to ameloriate constraint
     if (fpp*ic /= 0) then
        dirdx = sign(1d0,fpp*ic)
        !   ... Otherwise, choose dir to reduce the gradient
     else
        dirdx = sign(1d0,df1*dx1*ic)
     endif
     if (dirdx <= 0) xn = xmin - abs(dxmx)
     if (dirdx >= 0) xn = xmax + abs(dxmx)
     ir = -6
     goto 998

     ! --- CASE II: extremum bracketed without a root ---
     !                *                                      *
     !      *                  handle left case          *
     !           *             excluding right case                 *
     !     ------------------                          -----------------
  elseif (lextr .AND. fpp*f1 >= 0) then
     !    ... This is one possibility: return minimum point.
     !        xn = x0 - fp/fpp
     !        ir = 4
     !   ... This one returns a more distant point
     if (abs(xn-xmin) <= abs(xn-xmax)) then
        xn = xmax + abs(dxmx)
     else
        xn = xmin - abs(dxmx)
     endif
     ir = -6
     goto 998

     ! --- CASE V: no root bracketed: try linear extrapolation ---
  else
     if (f0 == f1) goto 80
     xn = x0 - f0*(x1-x0)/(f1-f0)
     if (xn < min(x0,x1) .AND. dabs(xn-min(x0,x1)) < dxmn) &
          xn = min(x0,x1) - dxmn
     if (xn > max(x0,x1) .AND. dabs(xn-max(x0,x1)) < dxmn) &
          xn = max(x0,x1) + dxmn
     !   ... do not allow step lengths exceeding dxmx
     xn = max(min(x0,x1)-abs(dxmx),xn)
     xn = min(max(x0,x1)+abs(dxmx),xn)
     !   ... But if inside range already seen, try a new point outside
     if (xn > xmin .AND. xn < xmax) goto 80
     ir = -4
  endif
  ! --- End of cases ---

  ! ... Convergence achieved or extremum bracketed with no root:
  !     Force xn on a prior point if 4's bit in 10's digit set.
  if ((ir == 0 .OR. ir == 4) .AND. mod(isw/10,10) > itol) then
     !   ... Look for smallest f: if not x0, swap with x0
     do  74  i = 2, 1, -1
        if (ir0 == -1 .AND. i == 2) goto 74
        if (abs(f(i)) < abs(f(0))) then
           call dswap(1,x(i),1,x,1)
           call dswap(1,f(i),1,f,1)
           call dswap(1,fx(i),1,fx,1)
           wk(12) = wk(12) + i
        endif
74   enddo
     xn = x0
  endif
  if (ir == 4) goto 999
  goto 998

  ! ... Exit ir=-5: suggest new point = max (min) given x +(-) abs(dxmx)
80 continue
  ir = -5
  if (dxmx < 0) xn = xmin - abs(dxmx)
  if (dxmx > 0) xn = xmax + abs(dxmx)

  ! ... Get a new point
998 continue
  if (ir0 == 0 .AND. ipr > 30) then
     !        call awrit2('  rfalsi ir=%,2i: seek xn=%1,8;8g',' ',scrwid,stdo,ir,xn)
     write(stdo,"(a,i5,d13.5)")'  rfalsi ir: seek xn=',ir,xn
  elseif (ipr > 30 .OR. ipr >= 30 .AND. ir >= 0) then
     !        i = awrite('%x  rfalsi ir=%,2i x1=%1;4g f1=%1;4g x2=%1;4g'//
     !     .  ' f2=%1;4g: %?#n#seek#est#',outs,scrwid,0,
     !     .  ir,x0,f0,x1,f1,ir,0,0)
     write(stdo,"(a,i0,2d13.5,2x,2d13.5)")'rfalsi ir,x0,f0,x1,f1=', ir,x0,f0,x1,f1
     !         i = awrite('%x  rfalsi ir=%,2i x1=%1;4g f1=%1;4g x2=%1;4g'//
     !     .  ' f2=%1;4g: %?#n#seek#est#',outs,scrwid,0,
     !     .  ,0,0)
     !        i = min(12,scrwid-3-i)
     !        call awrit2('%a x=%;nF',outs,scrwid,-stdo,i,xn)
     !       continue
  endif

  ! ... Exit: save data
999 continue
  wk(10) = xmin
  wk(11) = xmax
  call dcopy(3,x,1,wk(1),1)
  call dcopy(3,f,1,wk(4),1)
  call dcopy(3,fx,1,wk(7),1)

  !      if (ir .eq. -2 .or. ir .eq. -3) then
  !        if ((xn-x0)*(xn-x1) .gt. 0) then
  !          print *, 'exit', xn-x0,xn-x1
  !          stop 'oops'
  !        endif
  !      endif

end subroutine rfalsi

!     Testing:
!     rfalsi should return xn matching table xnj; output wk should match wkj
!      subroutine fmain

!      double precision xnj(10),fnj(10),xtol,ftol,dxmn,dxmx,wk(12),fn,xn
!      double precision wkj(12,10)
!      integer irj(10),isw,ir,j
!      data irj /0,-1,-2,-3,-3,-3,-2,-3,-3,-2/
!      data xnj /
!     .  -10.0000000000000d0,
!     .  -0.02d0,
!     .  -5.01000000000000d0,
!     .  -3.92391198839667d0,
!     .  -3.55881578946834d0,
!     .  -3.49752901570077d0,
!     .  -3.48943419141083d0,
!     .  -3.49366185972306d0,
!     .  -3.49366451990721d0,
!     .  -3.49366451990776d0/

!      data fnj /
!     .  -0.786470963100663d0,
!     .   1.27119252795989d0,
!     .  -0.254259440060593d0,
!     .   -7.92505718082347D-002,
!     .   -1.24254410725922D-002,
!     .   -7.41481145397956D-004,
!     .    8.12323636572358D-004,
!     .    5.10605076779703D-007,
!     .    1.06474508029808D-013,
!     .    2.58353017798529D-013/

!      data wkj /
!     .  -10.0000000000000D0,0.00000000000000D0,0.00000000000000D0,
!     . -0.786470963100663D0,0.00000000000000D0,0.00000000000000D0,
!     . -0.786470963100663D0,0.00000000000000D0,0.00000000000000D0,
!     .  -10.0000000000000D0,-10.0000000000000D0,0.00000000000000D0,

!     .  -1.99999999999996D-2,-10.0000000000000D0,-10.0000000000000D0,
!     .   1.27119252795989D0,-0.786470963100663D0,-0.786470963100663D0,
!     .   1.27119252795989D0,-0.786470963100663D0,0.00000000000000D0,
!     .  -10.0000000000000D0,-1.99999999999996D-2,0.00000000000000D0,

!     .  -5.01000000000000D0,-1.99999999999996D-2,-10.0000000000000D0,
!     . -0.254259440060593D0,1.27119252795989D0,-0.786470963100663D0,
!     . -0.254259440060593D0,1.27119252795989D0,-0.786470963100663D0,
!     .  -10.0000000000000D0,-1.99999999999996D-2,0.00000000000000D0,

!     .  -3.92391198839667D0,-1.99999999999996D-2,-5.01000000000000D0,
!     .  -7.92505718082347D-2,1.27119252795989D0,-0.254259440060593D0,
!     .  -7.92505718082347D-2,0.635596263979943D0,-0.254259440060593D0,
!     .  -10.0000000000000D0,-1.99999999999996D-2,4.00000000000000D0,

!     .  -3.55881578946834D0,-1.99999999999996D-2,-3.92391198839667D0,
!     .  -1.24254410725922D-2,1.27119252795989D0,-7.92505718082347D-2,
!     .  -1.24254410725922D-2,0.317798131989971D0,-7.92505718082347D-2,
!     .  -10.0000000000000D0,-1.99999999999996D-2,4.00000000000000D0,

!     .  -3.49752901570077D0,-1.99999999999996D-2,-3.55881578946834D0,
!     .  -7.41481145397956D-4,1.27119252795989D0,-1.24254410725922D-2,
!     .  -7.41481145397956D-4,0.158899065994986D0,-1.24254410725922D-2,
!     .  -10.0000000000000D0,-1.99999999999996D-2,4.00000000000000D0,

!     .  -3.48943419141083D0,-3.49752901570077D0,-1.99999999999996D-2,
!     .   8.12323636572358D-4,-7.41481145397956D-4,1.27119252795989D0,
!     .   8.12323636572358D-4,-7.41481145397956D-4,0.158899065994986D0,
!     .  -10.0000000000000D0,-1.99999999999996D-2,0.00000000000000D0,

!     .  -3.49366185972306D0,-3.49752901570077D0,-3.48943419141083D0,
!     .   5.10605076779703D-7,-7.41481145397956D-4,8.12323636572358D-4,
!     .   5.10605076779703D-7,-3.70740572698978D-4,8.12323636572358D-4,
!     .  -10.0000000000000D0,-1.99999999999996D-2,4.00000000000000D0,

!     .  -3.49366451990721D0,-3.49752901570077D0,-3.49366185972306D0,
!     .   1.06474508029808D-13,-7.41481145397956D-4,5.10605076779703D-7,
!     .   1.06474508029808D-13,-1.85370286349489D-4,5.10605076779703D-7,
!     .  -10.0000000000000D0,-1.99999999999996D-2,4.00000000000000D0,

!     .  -3.49366451990776d0,-3.49752901570077D0,-3.49366451990721D0,
!     .   2.58353017798529D-13,-7.41481145397956D-4,1.06474508029808D-13,
!     .   2.58353017798529D-13,-9.26851431747444D-5,1.06474508029808D-13,
!     .  -10.0000000000000D0,-1.99999999999996D-2,4.00000000000000D0/

!      data xtol/1D-12/,ftol/1D-12/,dxmn/1D-13/,dxmx/9.98d0/,isw/10/

!      call pshpr(55)

!      do  j = 1, 10
!        ir = irj(j)
!C       if (j .gt. 1) print *, 'j',j,xn-xnj(j),xn.lt.wk(1),xn.lt.wk(2)
!        if (ir .ne. 0) print *, xn,xn-xnj(j)
!        xn = xnj(j)
!        fn = fnj(j)
!        call rfalsi(xn,fn,xtol,ftol,dxmn,dxmx,isw,wk,ir)

!        print *, wk(:)-wkj(:,j)
!C       print *, wk(:)

!      enddo

!      end

