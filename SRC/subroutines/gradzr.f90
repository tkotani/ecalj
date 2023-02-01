module m_gradzr
  use m_ftox
  use m_lmfinit,only:stdo
  public:: gradzr
  private
contains
  subroutine gradzr(n,p,hess,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
    !- Find zero in gradients of a multivariate function
    ! ----------------------------------------------------------------------
    !i Inputs
    !i  n:    number of independent variables
    !i  hess: inverse Hessian matrix for Fletcher-Powell and Broyden
    !i  xtoll:x tolerance for line minimizations.  See Remarks.
    !i  dxmx: Maximum step length in any one component of x
    !i        Also, for C.G., initial step length; see Remarks.
    !i        A positive sign of dxmx makes line mins tend 'downhill'
    !i        A negative sign of dxmx makes line mins tend 'uphill'
    !i  xtol: tolerance in delta x for global minimization.
    !i        xtol > 0 : tolerance compared to |delta x|
    !i        xtol < 0 : tolerance compared largest component, delta max(|x|)
    !i  gtol: tolerance in gradient for global minimization
    !i        gtol > 0 : tolerance compared to |delta g|
    !i        gtol < 0 : tolerance compared largest component, delta max(|g|)
    !i  grfac:Extrapolation growth factor for line minimizations.
    !i        Whenever a root is not bracketed, the effective value of dxmx
    !i        is scaled by grfac.
    !i  isw:  compound of one-digit switches.
    !i    1's digit: handles constraints and step directions
    !i       *For the line minimization algorithms, digit is passed through
    !i        to the line min. routine rfalsi (which see).
    !i      0 no constraints imposed
    !i      1 require line min grad >0 at root (minimum for each line min)
    !i      2 require line min grad <0 at root (maximum for each line min)
    !i      4 while root not bracketed, set suggested x = current x + dxmx
    !i       *For the Broyden scheme, the digit is passed through to routine
    !i        brmin; the meaning is somewhat different (see brmin.f)
    !i        See Remarks for some further description.
    !i   10's digit: convergence criteria; see also Remarks.
    !i      0 convergence when <g> < gtol
    !i      1 convergence when <dx> < xtol
    !i      2 convergence when <g> < gtol and <dx> < xtol
    !i      3 convergence when <g> < gtol or <dx> < xtol
    !i        Here, <g> and <dx> are the max component of gradient and change
    !i        in x between iterations
    !i      4 Set this bit to specify a line minimization tolerance in the
    !i        gradient; see Remarks.
    !i  100's digit: governs which minimization  to use:
    !i      0 Congugate gradients, adapted from Numerical Recipes, 10.6.
    !i      1 Fletcher-Powell: M.J.Norgett & R.Fletcher, J.Phys.C3,L190(1972)
    !i      2 Broyden
    !i 1000's digit: Hessian matrix switches
    !i      0  no tests are made for Hessian
    !i      1  If Hessian is not positive definite, it is not updated.
    !i      2  project out parts of Hessian corresponding to negative
    !i         eval (not implemented).
    !i      4  On the second iteration, Hessian is globally scaled
    !i         once a first estimate is known (D. Novikov)
    !i 10000's digit: line-minimization-specific switches.
    !i      1 specifies tolerance for change in line minimization gradient;
    !i        see Remarks.
    ! o Inputs/Outputs
    ! o ir:   program flow control.  To start a new minimization,
    ! o       set ir to zero.  On exit, gradzr sets ir to one of the
    ! o       following.  For ir<0, gradzr expects a new set of points
    ! o       and gradients at positions p(1..n).
    ! o      *For the line minimization schemes F.P. and C.G:
    ! o    >=0: no more function calls expected:
    ! o      0: gradzr finds function converged to within tolerance
    ! o      1: input grad equals 0  gradzr does nothing.
    ! o      2: nonsensical tolerance.  gradzr does nothing.
    ! o      3: input xn equals x1 or x2.  gradzr does nothing.
    ! o      4: line min extremum bracketed but with no root.  gradzr returns
    ! o         suggested xn = largest (smallest) of given x so far + dxmx
    ! o      5: line minimization failed to bracket root;
    ! o         see 100's digit of isw
    ! o     <0: expects new function call, suggesting xn:
    ! o     -1: no information yet: require new point at suggested xn
    ! o         to improve the estimate of the root.
    ! o     -2: a root is bracketed and gradzr is attempting linear
    ! o         interpolation to improve the estimate of the root.
    ! o     -3: a root is bracketed and gradzr is attempting quadratic
    ! o         interpolation to improve the estimate of the root.
    ! o     -4: a root has not been bracketed or constraints not fufilled.
    ! o         gradzr will suggests xn from a linear extrapolation.
    ! o     -5: same as -4 above, but gradzr tries xn = current x + dmmx
    ! o         is increasing along the step direction.  If ir is set
    ! o         to -4 for a new point, gradzr will continue along that
    ! o         step direction
    ! o     -6: gradzr is having trouble; asks for new point
    ! o     -9: Line minimization converged; new line minization sought.
    ! o         For now, ir=-9 used just internally.
    ! o         NB: gradzr is organized so that this need not happen
    ! o    -10: ir input -10 works just as input ir=0, but flags that
    ! o         the conjugate direction p(nd) is already specified by
    ! o         the caller.  Valid only for line minimization schemes.
    ! o      *For Broyden (meaning same as in brmin.f but sign has changed)
    ! o      0: brmin has converged to specified tolerance
    ! o     <0: brmin needs gradients (preferably at the output x)
    ! o     >0: like <0, but the Hessian was found to be not positive
    ! o         definite, and was not updated.
    ! o p:    array of dimension at least 6*n for F.P. and C.G. and
    ! o       8*n for Broyden.  If the Hessian eigenvalues are calculated,
    ! o       dimensions must be n*(6+n) and n*(11+2n) for F.P. and Broyden
    ! o       p(*,nx=1) current position vector. (Input on first call)
    ! o                 Internally, p(nx) = p (start of current line min)
    ! o       p(*,ng=2) gradient at current point (Input on each call)
    ! o      *For line minimization schemes:
    ! o       p(*,nm=3) position vector where gradient was found to be minimum
    ! o                 data is kept for informational purposes only
    ! o       p(*,nd=4) conjugate direction for current line min.
    ! o       p(*,n1=5) gradient at the point prior to this one
    ! o       p(*,n2=6) gradient at second prior point
    ! o       p(*,n0=7) gradient at start of line min (C.G. only)
    ! o      *For Broyden, p(*,nd...8) are passed as w(*,1..6); see brmin.f
    ! o wk    is a work array of dimension 0:27.  It should not be altered
    ! o       between succesive calls.  The elements of wk are:
    ! o wk(0) = xn is a measure of the shift in coordinates
    ! o       between successive function calls.
    ! o      *For line minimizations, xn is the amount of direction vector
    ! o       that was added to coordinate positions; see Remarks.
    ! o      *For Broyden, xn is the largest shift in a component of x
    ! o wk(1..12) see rfalsi.
    ! o wk(13) gam = gamma for conjugate gradients; cf Numerical Recipes.
    ! o wk(14) dxmxl = maximum step length along direction vector
    ! o wk(15) dxmxlx = step size that makes shift in largest p(nx)= dxmx
    ! o wk(16) dxlast = Largest change in one component of x;
    ! o          used as estimate for dx in new C.G. line minimization.
    ! o wk(17) xmax = largest value of x along a given line min; used
    ! o          when line minimizations have trouble finding a root
    ! o          satisifying constraints.
    ! o wk(18) lminn = number of line minimizations so far
    ! o wk(19) grad0 = (grad.dir-vec) at start of line minimization
    ! o wk(20) growl = Current extrapolation growth factor
    ! o wk(21..23) x0h,x1h,x2h values of xn in prior iterations
    ! o wk(24..26) g0h,g1h,g2h values of grad.dir-vec corresponding to x*h
    ! o wk(27) global minimum gradient, corrsponding to p(nm)
    ! o        See 10's digit isw for criterion that measures gradient
    ! o        sign(wk(27)) used as a flag to indicate:
    ! o        >0 most recent |g| is also global minimum
    ! o        <0 most recent |g| is not global minimum
    ! o           In this case, p(*,nm) contains positions for global min
    !r Remarks
    !r *gradzr attempts to find a a minimum, a maximum,
    !r  or just a zero in the gradient using one of three schemes:
    !r  Broyden, Fletcher-Powell, and congugate gradients.  There are
    !r  strengths and weaknesses in each approach.  The F.P and C.G. proceed
    !r  by successive line minimizations, in which the zero projection of the
    !r  gradient along each direction vector is sought.  In the program,
    !r  p(*,nd) is the direction vector, and p(*,nx) are the coordinates at
    !r  the start of the line minimization.  xn is varied so that the for
    !r  positions p(nx) + xn * p(nd), grad.p(nd) is zero to within a specified
    !r  tolerance, after which the inverse Hessian matrix is updated and new
    !r  conjugate direction is calculated.  The F.P Hessian is explicit; in
    !r  the C.G. case it is implicit.  Having the Hessian is useful
    !r  but can be expensive if there are many variables.  In practice, F.P.
    !r  seems to converge a little faster than C.G..  Broyden is essentially a
    !r  Newton-Raphson scheme for several variables, and the Hessian is
    !r  explicit.  In well conditioned cases, it tends to get to the root more
    !r  rapidly than the line minimization approaches, except that it is more
    !r  sensitive to roundoff errors, which after a few iterations cause it
    !r  to converge more slowly.
    !r
    !r *gradzr works by accepting a new vector x and its gradient g, in
    !r  p(1..n) and p(n+1..2n) respectively.  It returns with a new
    !r  suggested value for p, for which the caller is expected to generate
    !r  g.  (For the line minimization routines, it is more natural for
    !r  gradzr to return a direction vector and the distance along the
    !r  projection vector, but this is not done here to make the calling
    !r  interface as consistent as possible.  The distance xn is returned in
    !r  case the caller wants to change it.)  Variable ir supplies some
    !r  information as to what gradzr is looking for in the next iteration.
    !r
    !r *Convergence tolerances.
    !r  gradzr continues its minimization procedure until it satisfies some
    !r  combination of:
    !r    (the change in position between iterations) < |xtol|
    !r    (the value  of the gradient)                < |gtol|
    !r  The l.h.s is evaluated either as the largest component
    !r  (xtol>0 or gtol>0) or as the length of the vector (xtol<0, gtol<0)
    !r  The 10's digit of isw specifies which combination is selected.
    !r
    !r  For Broyden, the change in x is the change between each function call;
    !r  for the line minimization algorithms C.G. and F.P., it is the change
    !r  between each line minimization.  In this latter case, there is a
    !r  second, independent tolerance for the individual line minimizations,
    !r  xtoll.  There is also a separate convergence criterion for gradient in
    !r  the individual line minimizations, namely the line minimization stops
    !r  when (grad.dir-vec) has dropped to a specified fraction f of its
    !r  starting value.  Turn this option on by adding 40 to isw.  You can
    !r  specify f using i=10000's digit of isw; f=i/20.
    !r  Using i=0 sets f to a default of 0.25.  NB: when no x tolerance is
    !r  specified, the line gradient tolerance is automatically turned on.
    !r
    !r *Initial step lengths for a new line minimization: when the hessian
    !r  is explicit, it specifies the step length.  You may specify a
    !r  maximum step length dxmx, which imposes a maximum value for a change
    !r  in any component of x.  The C.G. uses dxmx to specify the initial
    !r  step size for the first line minimization.  For subsequent line
    !r  mins, the initial step length is taken to be the total change in
    !r  position from the prior line min.
    !r
    !r *gradzr can be used to seek a minimum, a maximum,
    !r  or just a zero in the gradient; see one's digit of isw.
    !r  When a constraint is imposed, caller is strongly advised to
    !r  also set the 4's digit of isw, which will cause rfalsi to continue
    !r  in the same search direction even if the gradient increases.
    !r  Also, if a maximum is sought (isw=2), caller is strongly advised
    !r  to set dxmx negative.  Otherwise, strange behavior may result.
    !r  For now these switches are automatically set internally.
    !r
    !r *recommended values for line minimizations.
    !r  At least in well behaved cases, it is recommended that
    !r  the caller choose an x tolerance for the line minimizations.
    !r  A 'safe' default for the line tolerance xtoll is the same value
    !r  as the global xtol, though gradzr tends to convergence a little
    !r  faster if it is set somewhat larger.
    !r  Also, convergence seems to be aided by changing line minimizations
    !r  when the gradient has dropped by ~0.25 (add 40 to isw).
    !r
    !r *It is usually simpler and more convenient to call gradzr using
    !r  the driver routine drgrzr.
    !r
    !u  Updates
    !u   08 Mar 06 some improvements to convergence criteria; saves p(gmin)
    !u   07 Sep 03 gradzr returns after first iter
    !u             when g<gtol is satisfied and is sole conv. criterion
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed Parameters
    integer :: n,ir,isw
    double precision :: hess(n,n),p(n,10),xtoll,dxmx,xtol,gtol,wk(0:*), &
         grfac
    ! Local Variables
    logical :: ltmp,cnvgx,cnvgg,nwlin
    integer :: i,j,ipr,npr,lminn,isw0,isw1,isw2,isw3,isw4, &
         isw31,idamax,iswl,itol,iswb
    integer :: nx,ng,nd,n1,n2,n0,nm,nhs
    double precision :: gg,dgg,gam,ddot,dasum,x0h,x1h,x2h,evtol, &
         evmin,dxmxl,dxmxlx,xtll,dxmn,xmax,gfn,grad0,g0h,g1h,g2h, &
         dxm,growl,gammax,xh(0:2),gh(0:2),xn, &
         gd1,gd2,dg,alpha,r,s,gll
    equivalence (xh(0),x0h), (xh(1),x1h), (xh(2),x2h)
    equivalence (gh(0),g0h), (gh(1),g1h), (gh(2),g2h)
    parameter (nx=1,ng=2,nm=3,nd=4,n1=5,n2=6,n0=7,nhs=8, evtol=1d-10)
    double precision :: gwg,q,dxlast,dxtop,gtop,resx,resg,dum(1)
    integer :: scrwid
    parameter (scrwid=80)
    character*(100) outs

    ! --- Iteration independent setup ---
    if (ir == 0) call dpzero(wk,27+1)
    xn =     wk(0)
    gam =    wk(13)
    dxmxl =  wk(14)
    dxmxlx = wk(15)
    dxlast = wk(16)
    xmax   = wk(17)
    lminn  = wk(18)
    grad0  = wk(19)
    growl  = wk(20)
    call dcopy(3,wk(21),1,xh,1)
    call dcopy(3,wk(24),1,gh,1)
    call getpr(ipr)
    npr = min(n,6)
    isw0  = mod(isw,10)
    isw1  = mod(isw/10,10)
    isw2  = mod(isw/100,10)
    isw3  = mod(isw/1000,10)
    isw31 = mod(isw3,4)
    if (isw31 /= 0) call rx('gradzr: this branch is not used now...')
    isw4  = mod(isw/10000,10)
    itol = mod(isw1,4)
    !     Save min gradient, p(nm) using initial p and gradient
    if (ir == 0) then
       call pgradz(0,n,0d0,p,p(1,ng),0d0,0d0, &
            dxtop,resx,gtop,resg,wk(27),cnvgx,cnvgg)
       call dcopy(n,p(1,nx),1,p(1,nm),1)
    endif

    ! ... Update p if gradient to minimum gradient
    call pgradz(0,n,0d0,p,p(1,ng),0d0,0d0, &
         dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
    if (r < abs(wk(27))) then
       wk(27) = r
       call dcopy(n,p(1,nx),1,p(1,nm),1)
    else if (r > abs(wk(27))) then
       wk(27) = -abs(wk(27))
    endif

    ! ... These lines impose sanity for conditions described in Remarks
    if (mod(isw0,4) /= 0) isw0 = mod(isw0,4) + 4
    dxm = dxmx
    j = mod(isw0,4)
    if (dxmx > 0 .AND. j == 2 .OR. dxmx < 0 .AND. j == 1) &
         dxm = -dxmx

    ! ... Re-entry point for new line minimization
1   continue
    ! ... If the gradient is zero, nothing to calculate.
    gg = abs(p(idamax(n,p(1,ng),1),ng))
    if (gg == 0) then
       ir = 0
       goto 999
    endif

    ! --- First iteration ---
    if (ir == 0 .OR. (ir == -10 .AND. isw2 <= 1)) then
       !   ... This the largest dx from the last iteration
       !       For first line min, set it so largest dx is dxm
       dxlast = dxm
       !   ... Number of line minimizations
       if (ir == 0) lminn = 0
       !   ... Congugate gradients gamma (see Numerical Recipes)
       gam = 0
       !   ... Size of initial gradient (and direction) vector
       j = mod(isw0,4)
       ltmp = j .eq. isw0 .or. &
            dxm.gt.0 .and. j.eq.2 .or. dxm.lt.0 .and. j.eq.1
       if (ipr >= 30 .OR. ipr >= 20 .AND. ltmp) then
          !$$$          call awrit7(' gradzr: begin %?#n# xtol=|%1;3g|#%j#'
          !$$$     .    //'%?#n==2# and##%-1j%?#n==3# or##%-1j'
          !$$$     .    //'%?#n<>1# gtol=|%1;3g|#%j#'
          !$$$     .    //'%?#n==0|n>=4#  gtll=%1;3g#%j#'
          !$$$     .    //'  dxmx=%1;3g',outs,
          !$$$     .    scrwid,0,itol,xtol,itol,gtol,
          !$$$     .    isw1,pgrada(isw1,isw4),dxm)
          write(stdo,ftox)'itol,xtol,gtol,isw1,pgrada(isw1,isw4),dxm',itol,xtol,gtol,isw1,pgrada(isw1,isw4),dxm
          !          call info2(-25,0,0,
          !     .    outs//'%a  '//
          !     .    '%?#n==0#C.G.##%-1j'//
          !     .    '%?#n==1#F.P.##%-1j'//
          !     .    '%?#n==2#Broy##'//
          !     .    '  isw=%i',isw2,isw)
          write(stdo,ftox)'cg fp broy',isw2,isw
          !     ... Sanity checks
          !          if (j .eq. isw0)
          !     .    call awrit0(' gradzr (warning) seek extremum but 4''s'
          !     .      //' bit of isw is not set',' ',scrwid,i1mach(2))
          !          if (dxm.gt.0 .and. j.eq.2 .or. dxm.lt.0 .and. j.eq.1)
          !     .    call awrit0(' gradzr (warning) seek extremum but dxmx'//
          !     .      ' is of the wrong sign',' ',scrwid,i1mach(2))
          if (iand(isw1,1) == 0 .AND. gtol == 0) ir = 2
          if (iand(isw1+3,2) == 0 .AND. xtol == 0) ir = 2
          if (xtol == 0 .AND. gtol == 0) ir = 2
          !         if (ir .eq. 2) call rx('OOPS')
          if (ir == 2) goto 999
       endif
    endif

    ! --- Broyden minimization ---
    if (isw2 == 2) then
       ir = -ir
       iswb = isw0 + 10*iand(isw1,3) + 1000*isw3
       call brmin(n,p,p(1,ng),iswb,ipr,dxmx,xtol,gtol,dum,p(1,nd),xn, &
            hess,ir)
       ir = -ir
       goto 999
    endif

    ! ... Case continue with current line minimization:
    if (ir >= -6 .AND. ir <= -1) then
       !       Undo direction vector added to position vector,
       !       restoring position vector to starting point of current line min
       if (ir < 0) call daxpy(n,-xn,p(1,nd),1,p(1,nx),1)
       goto 20
    endif

    ! --- New line minimization ---
    lminn = lminn+1
    ! ... Congugate gradients
    if (isw2 == 0) then
       gg = ddot(n,p(1,ng),1,p(1,ng),1)
       if (gam == 0 .AND. ir /= -10) call dpzero(p(1,nd),n)
       dgg = ddot(n,p(1,ng),1,p(1,nd),1)
       resg = gam*dgg-gg
       gammax = 2
       if (resg > 0) then
          write(stdo,ftox) ' gradzr encountered g.h =' &
               //'  ... set gam to zero',ftod(resg)
          gam = 0
       elseif (gam > gammax) then
          write(stdo,ftox)' gradzr cap gam = ',' to max gam =',ftod(gam),ftod(gammax)
          gam = gammax
       endif
       call dpcopy(p(1,ng),p(1,n0),1,n,-1d0)
       if (ir == -10) then
          ir = 0
       else
          do    j = 1, n
             p(j,nd) = p(j,n0) + gam*p(j,nd)
          enddo
       endif
       dxmxl = sign(dxlast/p(idamax(n,p(1,nd),1),nd),dxm)
       ! ... Fletcher-Powell
    elseif (isw2 == 1) then
       !       If starting Hessian is zero, set to unity.
       if (dasum(n*n,hess,1) == 0) call dcopy(n,1d0,0,hess,n+1)
       !       F.P doesn't use p(n0), but save for compatibility with C. G.
       call dpcopy(p(1,ng),p(1,n0),1,n,-1d0)
       !       New positions generated from inverse Hessian matrix
       !       Shift in conjugate dir p(nd) = -hess . grad
       if (ir == -10) then
          ir = 0
       else
          do  6  i = 1, n
             p(i,nd) = 0
             do    j = 1, n
                p(i,nd) = p(i,nd) - hess(i,j)*p(j,ng)
             enddo
6         enddo
          dxmxl = 1
       endif
    else
       call rxi('gradzr not ready for isw2=',isw2)
    endif

    ! ... This is the 'best guess' for step size, new line min.
    dxmxl = dxmxl*min(1d0,dabs(dxm/(dxmxl*p(idamax(n,p(1,nd),1),nd))))
    ! ... Step size that makes shift in largest component = dxmx:
    dxmxlx = sign(dxmx/p(idamax(n,p(1,nd),1),nd),dxm)
    ! ... Initial extrapolation growth factor
    growl = 1
    ! ... (grad.dir-vec) at start of line minimization
    grad0 = ddot(n,p(1,ng),1,p(1,nd),1)
    if (ir == -9 .OR. ir == 0) then
       ir = 0
       xn = 0d0
    endif
    xmax = 0d0
    x1h = -9d9
    x2h = -9d9
    ! ... Printout of current p,g,h
    if (ipr > 30) then
       gg = ddot(n,p(1,ng),1,p(1,ng),1)
       dgg = ddot(n,p(1,ng),1,p(1,nd),1)
       resg = dsqrt(ddot(n,p(1,ng),1,p(1,ng),1))
       gtop = p(idamax(n,p(1,ng),1),ng)
       ! call info5(-30,0,0,' gradzr new line %i:  g.h=%;4g'// &
       !      '  g.(h-g)=%;4g  max g=%1;3g  |grad|=%1;3g  ', &
       !      lminn,dgg,dgg+gg,gtop,resg)
    endif
    ! call info2(-40,0,0,'  p=%n:;10F',npr,p(1,nx))
    ! call info2(-40,0,0,'  g=%n:;10F',npr,p(1,ng))
    ! call info2(-40,0,0,'  h=%n:;10F',npr,p(1,nd))

    ! --- Begin or continue current line minimization ---
20  continue
    xmax = max(xmax,xn)
    x0h = xn
    gfn = ddot(n,p(1,ng),1,p(1,nd),1)
    g0h = dsqrt(ddot(n,p(1,ng),1,p(1,ng),1))

    ! ... Decide on global convergence before calling rfalsi,
    !     since criteria might not be identical to that in rfalsi

    !     Global convergence criteria: shift relative to last position
    !     call pgradz(itol,n,xn-x0h,p(1,nd),p(1,ng),
    !     .  dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
    !     Global convergence criteria: shift relative to start of line min
    nwlin = ir .eq. 0
    call pgradz(itol,n,xn,p(1,nd),p(1,ng),xtol,gtol, &
         dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
    !     Case gradzr has converged globally
    if (resg == 0 .OR. (cnvgg .AND. cnvgx)) goto 998

    ! ... New line min if g < line g tol
    if (dabs(gfn/grad0) < pgrada(isw1,isw4)) then
       ir = 0

       ! ... Otherwise, call rfalsi to check for convergence and/or new xn
    else
       !       Requires rfalsi to return xn on an existing point when ir=0
       iswl = isw0 + 40
       !       Also suppress rfalsi using gtol; we do that above
       if (mod(isw1,4) /= 0) iswl = iswl + 10
       !       Tolerance should be less than and distinct from dxmx
       xtll = min(xtoll/dabs(p(idamax(n,p(1,nd),1),nd)), &
            abs(dxmxl*growl*.999999d0))
       !       Minimimum step size should be less than, distinct from tolerance
       dxmn=xtll/2
       gll = abs(gfn/10)
       call pshpr(ipr-10)
       call rfalsi(xn,gfn,xtll,gll,dxmn,dxmxl*growl,iswl,wk(1),ir)
       call poppr
       !       Exit if 1st step of new line satisfies 'global convergence'
       if (nwlin) then
          call pgradz(itol,n,xn,p(1,nd),p(1,ng),xtol,gtol, &
               dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
          !         Case gradzr has converged globally
          if (resg == 0 .OR. (cnvgg .AND. cnvgx)) then
             xn = 0
             goto 998
          endif
       endif
       !   ... On the first line min movement, reset dxmx to maximum allowed
       if (mod(iswl,10) == 0 .AND. nwlin) dxmxl = dxmxlx
       !   ... rfalsi wants to extrapolate
       if (ir >= -6 .AND. ir <= -4) then
          !     ... This accelerates linear extrapolation
          if (ir == -4) xn = x0h + (xn-x0h)*growl
          growl = growl * grfac
          !        elseif (ir .eq. -1) then
          !          growl = grfac
       else
          growl = 1
       endif
       !   ... Swap points in the same way rfalsi swapped them.
       if (nint(wk(12)) >= 4) then
          call dswap(n,p(1,n1),1,p(1,n2),1)
          call dswap(1,xh(1),1,xh(2),1)
          call dswap(1,gh(1),1,gh(2),1)
       endif
       if (mod(nint(wk(12)),4) >= 2) then
          call dswap(n,p(1,ng),1,p(1,n2),1)
          call dswap(1,xh(0),1,xh(2),1)
          call dswap(1,gh(0),1,gh(2),1)
       endif
       if (mod(nint(wk(12)),2) >= 1) then
          call dswap(n,p(1,ng),1,p(1,n1),1)
          call dswap(1,xh(0),1,xh(1),1)
          call dswap(1,gh(0),1,gh(1),1)
       endif
       ! ... Handle special case ir=0, xn=0
       if (ir == 0 .AND. xn == 0) then
          ir = -2
          xn = (wk(5)*wk(1)-wk(4)*wk(2))/(wk(5)-wk(4))
       endif
    endif

    ! ... Case line minimization not converged
    if (ir < 0) then

       !   ... After clean, replace xh with wk.  for now:
       if (x0h /= wk(1) .OR. x1h /= wk(2) .AND. x1h /= -9d9 &
            .OR. x2h /= wk(3) .AND. x2h /= -9d9) then
          print *, x0h,wk(1)
          print *, x1h,wk(2)
          print *, x2h,wk(3)
          call rx('bug in gradzr or rfalsi')
       endif

       !   ... preserve gradient for last two points
       x2h = x1h
       g2h = g1h
       call dcopy(n,p(1,n1),1,p(1,n2),1)
       x1h = x0h
       g1h = g0h
       call dcopy(n,p(1,ng),1,p(1,n1),1)

    elseif (ir == 4) then
       xn = xmax + dxmxl
!       call info2(-10,0,0,'%x gradzr: found extremum without root:'// &
!            '  attempt new xn=%1;4g',xn,0)
       ir = -5
       goto 991

    elseif (ir /= 0) then
!       call info2(-1,0,0,' gradzr (abort) line %i encountered ir=%i '// &
!            'from rfalsi',lminn,ir)
       goto 991

       ! --- Line minimization has converged ---
    else
       !       Sanity check
       if (xn /= x0h) call rx('bug in gradzr')
       if (isw2 == 0) then
          !     ... Make gam for this line
          gg = 0d0
          dgg = 0d0
          do  21  j = 1, n
             gg = gg + p(j,n0)**2
             !           Use the following line for Fletcher-Reeves:
             !           dggfp = dgg + p(j,ng)**2
             !           or the following line for Polak-Ribiere:
             dgg = dgg + (p(j,ng)+p(j,n0))*p(j,ng)
21        enddo
          if (gg == 0d0) then
          else
             gam = dgg/gg
             !           gamfp = dggfp/gg
          endif
       elseif (isw2 == 1) then

          !         Local copy of the original Hessian, in case used later
          if (isw31 /= 0) call dcopy(n*n,hess,1,p(1,nhs),1)

          !     ... Update inverse of Hessian
          alpha = x0h - x1h
          gd2 = ddot(n,p(1,ng),1,p(1,nd),1)
          gd1 = ddot(n,p(1,n1),1,p(1,nd),1)
          dg = alpha * (gd2 - gd1)
          gwg = 0
          do  33  i = 1, n
             p(i,n2) = 0
             do   j = 1, n
                p(i,n2) = p(i,n2) + hess(i,j) * (p(j,ng) - p(j,n1))
             enddo
             gwg = gwg + (p(i,ng) - p(i,n1)) * p(i,n2)
33        enddo
          q = 1 + gwg/dg
          gd1 = 0
          do  36  i = 1, n
             r = - alpha*p(i,nd)/dg
             s = (-p(i,n2) + q * alpha*p(i,nd))/dg
             do   j = 1, n
                hess(i,j) = hess(i,j) + r * p(j,n2) + s * alpha*p(j,nd)
             enddo
36        enddo

          !     ... Check that Hessian is positive definite
          if (isw31 /= 0) then
             call rx('this branch is not used now...')
             ! !           Keep a local copy of hessian, since dsev1 destroys it.
             ! call dcopy(n*n,hess,1,p(1,nhs+n),1)
             ! j =chkhss(p(1,nhs+n),n,p(1,n2),evtol,isw31,p(1,nhs),p(1,n1))
             ! evmin = p(1,n1)
             ! if (j > 0 .AND. (isw31 == 1)) then
             !    call dcopy(n*n,p(1,nhs),1,hess,1)
             !    ir = -ir
             ! elseif (j > 0) then
             ! endif
          endif

       endif

       !       Internally flag to start new line min
       ir = -9
       !       Restore position vector to current position
       call daxpy(n,xn,p(1,nd),1,p(1,nx),1)
       !       Largest change in component of x relative to start of line min
       !       Note: used by C.G. line min.
       dxlast = xn*p(idamax(n,p(1,nd),1),nd)

       !       Printout for line min convergence
       ! call info8(-30,0,0,'%x gradzr cvg line %i:'// &
       !      '%?#n==0#  gam=%1;3g#%j#  x=%1;8d'// &
       !      '  |g.h|=%1;3g  dxmax=%1;3g'// &
       !      '%?#n#  evmin=%1;2g#%j#', &
       !      lminn,isw2,gam,xn,gfn,dxlast,isw31,evmin)
    endif

    ! --- Cleanup ---
991 continue
    !     Add direction vector to position vector
    if (ir >= -6 .AND. ir <= -1) &
         call daxpy(n,xn,p(1,nd),1,p(1,nx),1)
    if (ir == -9) goto 1

    ! --- Restore local variables needed to preserve; exit ---
999 continue
    wk(0) =  xn
    wk(13) = gam
    wk(14) = dxmxl
    wk(15) = dxmxlx
    wk(16) = dxlast
    wk(17) = xmax
    wk(18) = lminn
    wk(19) = grad0
    wk(20) = growl
    call dcopy(3,xh,1,wk(21),1)
    call dcopy(3,gh,1,wk(24),1)

    !     Use p(minimum g), if specified and different from p
    !      if (isw5 .eq. 1 .and. ir .ge. 0) then
    !        call info2(-20,0,0,'restore vector for minimum g, g=%1;3g',
    !     .    wk(27),0)
    !        call dcopy(n,p(1,nm),1,p,1)
    !      endif
    return

    ! --- Global convergence FP or CG ---
998 continue
    ! call info5(-20,0,0, &
    !      '%x gradzr converged to dxmax=%1;3g, |dx|=%1;3g,'// &
    !      '  gmax=%1;3g, |grad|=%1;3g in %i line min', &
    !      dxtop,resx,gtop,resg,lminn)
    ! call info2(-30,0,0,' p=%n:;10F',npr,p(1,nx))
    ! call info2(-30,0,0,' g=%n:;10F',npr,p(1,ng))

    !     Restore position vector to current position
    call daxpy(n,xn,p(1,nd),1,p(1,nx),1)

    !     Clean up and exit
    ir = 0
    goto 999
  end subroutine gradzr
!   subroutine drgrzr(n,pnew,gnew,p,hess,xtoll,dxmx,xtol,gtol,grfac, &
!        wk,copt,isw,ir)
!     !- Driver routine for gradzr
!     ! ----------------------------------------------------------------------
!     !i Inputs
!     !i   n      :size of vector; see gradzr
!     !i   pnew   :new positions for which gradient was obtained
!     !i   gnew   :gradient corresponding to new positions
!     !i  The following inputs are passed directly through to gradzr.
!     !i  See that routine for further description.
!     !i   hess   :inverse Hessian matrix for Fletcher-Powell and Broyden
!     !i   xtoll  :x tolerance for line minimizations.  See gradzr
!     !i   dxmx   :Maximum step length in any one component of x.
!     !i   xtol   :tolerance in x for global minimization.
!     !i   gtol   :tolerance in gradient for global minimization
!     !i   grfac  :Extrapolation growth factor for line minimizations.
!     !i   copt   :character string for a convenient specification of special
!     !i           options.  drgrzr converts the following strings in copt
!     !i           into the corresponding switches in isw.  These may be
!     !i           strung together, separated by spaces.  isw is altered
!     !i           only when input ir=0.
!     !i             'def' Set default isw
!     !i             'cg'  congugate gradients
!     !i             'fp'  Fletcher-Powell
!     !i             'br'  Broyden
!     !i             'min' specify minimization
!     !i             'max' specify maximization
!     !i           Options are read left-to-right, so when incompatible
!     !i           switches are set, the last takes precedence.
!     !i   isw    :compound of one-digit switches.  It will be altered
!     !i           if copt is set.
!     !o Inputs/Outputs
!     ! o  ir     :flow control passed to gradzr.
!     ! o          To start a new minimization, set ir to zero.
!     ! o          On exit, gradzr sets ir to a value as described in gradzr.
!     ! o          A return of ir=0 => convergence achieved to prescribed tol
!     ! o          A return of ir>0 => gradzr had trouble and is aborting.
!     ! o  p      :position and work array for input to gradzr.
!     ! o          Array is dimensioned at least 6*n for F.P. and C.G. and 8*n
!     ! o          for Broyden.  If the Hessian eigenvalues are calculated,
!     ! o          it increases to n*(6+n) and n*(11+2n) for F.P. and Broyden.
!     ! o          p should remain untouched between successive calls
!     ! o          drgrzr replaces p(*,nx):replaced by pnew and
!     ! o          p(*,nd) changed to p(nd) = xn * (h(xn)-h(xn=0))
!     ! o          NB: input ir=0 => p(nd) not defined and not changed
!     ! o  wk     :work array passed for input to gradzr.
!     ! o         :wk should remain untouched between successive calls
!     !r Remarks
!     !r   This is a driver routine for gradzr, taking as input positions
!     !r   and gradients, and internally updating the gradzr matrix p.
!     !r   For line minimizations, it also resets the conjugate direction
!     !r   vector h = p(1,nd) in the event the positions passed to drgrzr
!     !r   do not correspond to those kept by p.
!     !u Updates
!     ! ----------------------------------------------------------------------
!     !     implicit none
!     ! ... Passed parameters
!     character*(*) copt
!     integer :: n,ir,isw
!     double precision :: pnew(n),gnew(n)
!     double precision :: hess(n,n),p(n,10),xtoll,dxmx,xtol,gtol,wk(0:26), &
!          grfac
!     ! Local variables
!     integer :: idx,nx,ng,nd,n0,idamax,j1,j2,getdig
!     double precision :: dhmax,hold,hmax,tol,xn
!     parameter (nx=1,ng=2,nd=3,n0=6,tol=1d-4)

!     ! --- Set switches based on copt ---
!     if (ir == 0 .AND. copt /= ' ') then
!        j1 = 1
! 10     continue
!        call nword(copt,1,j1,j2)
!        if (j2 >= j1) then
!           if (copt(j1:j2) == 'cg') then
!              isw = isw + 100*(0-getdig(isw,2,10))
!           elseif (copt(j1:j2) == 'fp') then
!              isw = isw + 100*(1-getdig(isw,2,10))
!           elseif (copt(j1:j2) == 'br') then
!              isw = isw + 100*(2-getdig(isw,2,10))
!           elseif (copt(j1:j2) == 'def') then
!              isw = 40
!           elseif (copt(j1:j2) == 'min') then
!              isw = isw + 1*(1-mod(getdig(isw,0,10),4))
!           elseif (copt(j1:j2) == 'max') then
!              isw = isw + 1*(2-mod(getdig(isw,0,10),4))
!           else
!              call rxs2('drgrzr:  bad option, "',copt(j1:j2),'"')
!           endif
!           j1 = j2+1
!           goto 10
!        endif
!     endif

!     xn = wk(0)
!     ! ... hnew = hold + 1/xn (pnew - pold)
!     if (ir /= 0 .AND. xn /= 0d0) then
!        idx = idamax(n,p(1,nd),1)
!        hmax = abs(p(idx,nd))
!        !       This destroys p(nx) but it isn't needed anymore
!        call daxpy(n,-1d0,pnew,1,p(1,nx),1)
!        idx = idamax(n,p(1,nx),1)
!        dhmax = abs(p(idx,nx)/xn)
!        hold = abs(p(idx,nd))
!        call daxpy(n,-1d0/xn,p(1,nx),1,p(1,nd),1)
!        !   ... Change in direction exceeds tolerance; reset ir
!        if (dhmax > tol*hmax) then
!           call info2(-30,0,0,' gradzr: reset conjugate'// &
!                ' direction.  dhmax = %1,3;3g  hold = %1,3;3g',dhmax,hold)

!           !     ... For ir=-10, need p(nx) = p(xn=0), p(ng) = -p(n0)
!           call dcopy(n,pnew,1,p(1,nx),1)
!           call daxpy(n,-xn,p(1,nd),1,p(1,nx),1)
!           call dpcopy(p(1,n0),p(1,ng),1,n,-1d0)

!           !     ... Call gradzr to reset C.G. or F.P. for new line min
!           ir = -10
!           call gradzr(n,p,hess,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
!           wk(0) = xn

!        endif
!     endif

!     ! ... Set p(nx) to pnew and p(ng) to gnew
!     call dcopy(n,pnew,1,p(1,nx),1)
!     call dcopy(n,gnew,1,p(1,ng),1)

!     ! ... Call gradzr for next step in minimization
!     call gradzr(n,p,hess,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
!   end subroutine drgrzr

  real(8) function pgrada(isw1,isw4)
    !- Line minimization gradient tolerance
    !     implicit none
    integer :: isw1,isw4
    double precision :: gtll

    gtll = 0
    !          requested   req'd when xtol not used
    if (isw1 >= 4 .OR. mod(isw1,4) == 0) then
       gtll = dble(isw4)/20
       !       If zero, use a default value
       if (gtll == 0) gtll = .25d0
    endif

    pgrada = gtll
  END function pgrada
  subroutine pgradz(itol,n,amp,p,g,xtol,gtol,dxmx,resx,dgmx,resg, gx,cnvgx,cnvgg)
    !- Evaluates convergence criteria in x and g
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   itol  :flags which criterion to use
    !i         :0 => convergence when <g> < gtol
    !i         :1 => convergence when <p> < xtol and amp ne 0
    !i         :2 => convergence when <g> < gtol and (<p> < xtol and amp ne 0)
    !i         :3 => convergence when <g> < gtol or  (<p> < xtol and amp ne 0)
    !i   n     :number of elements
    !i   amp   :scaling factor for p
    !i   p     :unscaled vector for position shift
    !i   g     :gradient vector
    !o Outputs
    !o   dxmx  :largest change in single component of p
    !o   resx  :change in length of p
    !o   dgmx  :largest change in single component of g
    !o   resg  :change in length of g
    !o   gx    :whichever of dgmx or resg used pgradz uses as measure of g
    !o         :gx is always positive
    !o   cnvgx :T if convergence condition met in x
    !o   cnvgg :T if convergence condition met in g
    !l Local variables
    !l         :
    !r Remarks
    !u Updates
    !u   08 Mar 06 adapted from gradzr
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: itol,n
    double precision :: amp,p(n),g(n),xtol,gtol,dxmx,resx,dgmx,resg,gx
    logical :: cnvgx,cnvgg
    ! ... Local parameters
    integer :: idamax
    double precision :: ddot

    !     Largest component of p
    dxmx = amp*p(idamax(n,p,1))
    !     Change in length of p
    resx = amp*dsqrt(ddot(n,p,1,p,1))
    !     Largest component of g
    dgmx = g(idamax(n,g,1))
    !     Residual of gradient at new position and printout
    resg = dsqrt(ddot(n,g,1,g,1))

    !     Find cnvgx
    if (xtol >= 0) gx = abs(resx)
    if (xtol < 0) gx = abs(dxmx)
    cnvgx = gx .le. abs(xtol)
    if (amp == 0)  cnvgx = .FALSE. 
    if (itol == 0) cnvgx = .TRUE. 

    !     Find cnvgg
    if (gtol >= 0) gx = abs(resg)
    if (gtol < 0) gx = abs(dgmx)
    cnvgg = abs(gx) .le. abs(gtol)
    if (itol == 1) cnvgg = .TRUE. 

    !     If either is permissible and one is satisfied, set both true
    if (itol == 3 .AND. (cnvgx .OR. cnvgg)) then
       cnvgx = .true.
       cnvgg = .true.
    endif
  end subroutine pgradz

subroutine brmin(n,x,g,isw,ipr,dxmx,xtol,gtol,hmax,w,diff,hess,ir)
  !use m_gradzr,only:pgradz
  !use m_mathlib,only:chkhss
  !- One Broyden step in finding the root of an n-dimensional function
  ! ----------------------------------------------------------------
  !i Inputs
  !i   n:     number of variables
  !i   x      current for which gradients are specified.
  !i   g      gradient at current position
  !i   isw    1s digit (not implemented)
  !i          0  find minimum
  !i          1  find maximum
  !i         10s digit governs convergence criterion:
  !i          0 convergence when <g> < gtol
  !i          1 convergence when <dx> < xtol
  !i          2 convergence when <g> < gtol and <dx> < xtol
  !i          3 convergence when <g> < gtol or <dx> < xtol
  !i            Here, <g> and <dx> are the max component of gradient
  !i            and change in x between iterations
  !i        100s digit
  !i          1  Set this bit if there is an x-dependent hmax.
  !i       1000s digit: Hessian matrix.  4's bit is independent from 1,2
  !i          0  no tests are made for Hessian
  !i          1  If Hessian is not positive definite, it is not updated.
  !i          2  project out parts of Hessian corresponding to negative
  !i             eval (not implemented).
  !i          4  On the second iteration, Hessian is globally scaled
  !i             once a first estimate is known (D. Novikov)
  !i   ipr    verbosity
  !i   dxmx   maximum allowed change in x for this step
  !i   xtol  tolerance in delta x for global minimization.
  !i         xtol > 0 : tolerance compared to |delta x|
  !i         xtol < 0 : tolerance compared largest component, delta max(|x|)
  !i   gtol  tolerance in gradient for global minimization
  !i         gtol > 0 : tolerance compared to |delta g|
  !i         gtol < 0 : tolerance compared largest component, delta max(|g|)
  !i   hmax   a variable-specific maximum; shift is scaled so that
  !i          no variable exceeds maximum (not tested).  Only used for
  !i          those hmax>0.
  !i   w      work array of dimension at least w(n,6), or if Hessian
  !i          tested for positive definiteness, w(n,9+2*n), and if
  !i          negative evals projected out, w(n,9+2*n).  Pieces signify:
  !i          (*,nx)  prior position x (must be preserved between calls)
  !i          (*,ng)  input gradient g (must be preserved between calls)
  !i          (*,ns)  shift = x(output) - x(input)
  !i          (*,nd)  change in g
  !i          (*,nh)  a help array
  !i          (*,nh2) a help array
  !i          (*,nhs) a local copy of the Hessian (needed only if
  !i                  Hessian is tested for positive definiteness)
  ! o  ir:    iteration number.  Should be zero on first iteration.
  ! o         It is updated internally, and returns the following:
  ! o         0  brmin has converged to specified tolerance
  ! o        >0  brmin needs gradients (preferably at the output x)
  ! o        <0  like >0, but the Hessian was found to be not positive
  ! o            definite, and was not updated.
  ! o  hess inverse Hessian.  brmin updates hess as minimization
  ! o         proceeds.  Initial Hessian may be specified; it is set to
  ! o         unity if initially zero.
  !o Outputs
  !o   x      suggested new position at which to evaluate gradients
  !o   diff   the largest shift
  !r Remarks
  !r   Adapted from D. Novikov.
  !r   This routine generates the inverse of the Hessian matrix through
  !r   the Broyden scheme, and supplies an update of the coordinates.
  !u Updates
  !u   08 Mar 06 Modified convergence checks for consistency w/ gradzr
  ! ----------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n,ir,isw,ipr
  double precision :: diff,dxmx,xtol,gtol,x(n),g(n),hmax(n),w(n,10), &
       hess(n,n)
  ! ... Local parameters
  integer :: nx,ng,ns,nd,nh,nh2,nhs,itol
  parameter (nx=1,ng=2,ns=3,nd=4,nh=5,nh2=6,nhs=nd+5)
  double precision :: dxtop,resx,gtop,resg,dum
  logical :: cnvgx,cnvgg
  double precision :: beta,one,evtol
  parameter (beta=1d0, one=1d0, evtol=1d-10)
  double precision :: factor,fmin,a,b,ddot,dasum,gmax,evmin,scale
  !     double precision wk(n)
  integer :: i,j,isw1,isw32,isw2,isw31,npr

  isw1 = mod(isw/10,10)
  isw2 = mod(isw/100,10)
  isw31 = mod(mod(isw/1000,10),4)
  isw32 = mod(isw/1000,10)/4
  scale = 1
  evmin = 0
  if (isw31 /= 0) call rx('brmin: isw31/=0 is not used now')
  call pshpr(ipr)

  ! --- First iteration or restart ---
  if (ir == 0) then
     !   ... If starting Hessian is zero, set to unity.
     if (dasum(n*n,hess,1) == 0) call dcopy(n,1d0,0,hess,n+1)
     ir = 1
     do  6  i = 1, n
        w(i,nx) = x(i)
        w(i,ng) = g(i)
        do    j = 1, n
           x(i) = x(i) - hess(i,j)*g(j)
        enddo
6    enddo

     ! call info8(-25,0,0,' brmin: start'// &
     !      '%?#(n==1|n==2|n==3)#  xtol=|%1;2g|#%j#'// &
     !      '%?#(n==0|n==2|n==3)#  gtol=|%1;2g|#%j#'// &
     !      '  isw=%i  |g|=%1;3g', &
     !      isw1,xtol,isw1,gtol,isw,dsqrt(ddot(n,g,1,g,1)),0,0)

     ! --- Subsequent iterations  ---
  else

     ir = iabs(ir)+1

     !   ... Global scaling of the Hessian matrix
     if (ir == 2 .AND. isw32 /= 0) then
        do    i = 1, n
           w(i,ns) = x(i) - w(i,nx)
           w(i,nd) = g(i) - w(i,ng)
        enddo
        do    i = 1, n
           w(i,nh) = 0d0
           do    j = 1, n
              w(i,nh) = w(i,nh) + hess(i,j)*w(j,nd)
           enddo
        enddo
        a = ddot(n,w(1,nd),1,w(1,nh),1)
        b = ddot(n,w(1,nd),1,w(1,ns),1)
        scale = dabs(b/a)
        call dscal(n*n,scale,hess,1)
     endif

     !   --- Make dx,dg; save Hessian,x,g; make some work arrays ---
     do  20  i = 1, n
        w(i,ns) = x(i) - w(i,nx)
        w(i,nd) = g(i) - w(i,ng)
20   enddo
     !       Local copy of the original Hessian
!     if (isw31 /= 0) call dcopy(n*n,hess,1,w(1,nhs),1)
     do  30  i = 1, n
        w(i,nx) = x(i)
        w(i,ng) = g(i)
        w(i,nh) = 0d0
        do   j = 1, n
           w(i,nh) = w(i,nh) + hess(i,j)*w(j,nd)
        enddo
30   enddo
     a = ddot(n,w(1,nd),1,w(1,nh),1)
     b = ddot(n,w(1,nd),1,w(1,ns),1)
     do    i = 1, n
        w(i,nh2) = 1/a*w(i,nh) - 1/b*w(i,ns)
     enddo

     !   --- Update inverse of Hessian ---
     do    i = 1, n
        do   j = i, n
           hess(j,i) = hess(j,i) &
                - w(j,nh)*w(i,nh)/a + w(j,ns)*w(i,ns)/b &
                + beta*a*b*w(j,nh2)*w(i,nh2)
        enddo
     enddo
     !   ... Make it symmetric
     do    i = 2, n
        do   j = 1, i-1
           hess(j,i) = hess(i,j)
        enddo
     enddo

     !   --- Check that Hessian is positive definite ---
     ! if (isw31 /= 0) then
     !    !         Keep a local copy of hessian, since dsev1 destroys it.
     !    call dcopy(n*n,hess,1,w(1,nhs+n),1)
     !    j = chkhss(w(1,nhs+n),n,w(1,nd),evtol,isw31,w(1,nhs),w(1,ns))
     !    evmin = w(1,ns)
     !    if (j > 0 .AND. (isw31 == 1)) then
     !       call dcopy(n*n,w(1,nhs),1,hess,1)
     !       ir = -ir
     !    elseif (j > 0) then
     !    endif
     ! endif

     !   --- Update x from Hessian, g ---
     do  i = 1, n
        do  j = 1, n
           x(i) = x(i) - hess(i,j)*g(j)
        enddo
     enddo

     !   ... End of subsequent Hessian update
  endif

  ! --- Limit the shift, depending on dxmx and hmax ---
  do  i = 1, n
     w(i,ns) = x(i) - w(i,nx)
  enddo
  diff = 0
  do  i = 1, n
     diff = max(diff,dabs(w(i,ns)))
  enddo
  if (diff == 0) goto 100
  if (diff > dxmx) then
     factor = dxmx/diff
     print 334, diff,factor
334  format(' brmin: max shift =',f9.5, &
          ' is larger than dxmx.  Scale by',f8.4)
     do  84  i = 1, n
        w(i,ns) = factor * w(i,ns)
        x(i) = w(i,ns) + w(i,nx)
84   enddo
     diff = factor*diff
  endif
  ! ... Check local shifts
  if (isw2 /= 0) then
     fmin = 1
     do  85  i = 1, n
        if (hmax(i) > 0) then
           factor = min(one,abs(hmax(i)/x(i)))
           if (factor < fmin) then
              fmin = factor
              print '('' Variable No.'',i3,'' is out of range'',f8.6)', &
                   i,x(i)
           endif
        endif
85   enddo
     if (fmin < one) then
        print '('' Rescale all shifts by'',f8.6)', fmin
        do  86  i = 1, n
           x(i) = fmin*x(i)
           w(i,ns) = x(i) - w(i,nx)
86      enddo
     endif
  endif

  ! --- Cleanup ---
100 continue
  diff = 0
  gmax = 0
  do  i = 1, n
     gmax = max(gmax,dabs(g(i)))
     diff = max(diff,dabs(w(i,ns)))
  enddo

  ! if (scale /= 1) then
  !    call info8(-30,0,0,' brmin: iter %i  max shift=%1;4g  '// &
  !         '|g|=%1;4g  max g=%1;4g  scale H by %1;3g', &
  !         iabs(ir),diff,dsqrt(ddot(n,g,1,g,1)),gmax,scale,0,0,0)
  ! else
  !    call info8(-30,0,0,' brmin: iter %i  max shift=%1;4g  '// &
  !         '|g|=%1;4g  max g=%1;4g%?#n#  evmin=%1;2g#%j#',iabs(ir), &
  !         diff,dsqrt(ddot(n,g,1,g,1)),gmax,isw31,evmin,0,0)
  ! endif
  npr = min(n,8)
  ! call info0(-60,0,0,' Updated x: ')
  ! do  i = 1, n, npr
  !    call info(-60,0,0,' %n;10,6D',min(npr,n-i+1),x(i))
  ! enddo

  !     Global convergence criteria: shift relative to start of line min
  itol = mod(isw1,4)
  call pgradz(itol,n,1d0,w(1,ns),g,xtol,gtol, &
       dxtop,resx,gtop,resg,dum,cnvgx,cnvgg)

  ! ... Case brmin has converged globally
  if (resg == 0 .OR. (cnvgg .AND. cnvgx)) then
     ! call info5(-20,0,0, &
     !      '%x brmin converged to dxmax=%1;3g, |dx|=%1;3g,'// &
     !      '  gmax=%1;3g, |grad|=%1;3g in %i iter', &
     !      dxtop,resx,gtop,resg,ir)
     npr = min(n,6)
     ! call info(-30,0,0,' p=%n:;10F',npr,x)
     ! call info(-30,0,0,' g=%n:;10F',npr,g)
     ir = 0
  endif

  call poppr


end subroutine brmin

end module m_gradzr
