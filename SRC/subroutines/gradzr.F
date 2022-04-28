      module m_gradzr
      use m_ftox
      use m_lmfinit,only:stdo
      public:: gradzr,drgrzr,pgradz,chkhss
      private
      contains
      subroutine gradzr(n,p,hess,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
C- Find zero in gradients of a multivariate function
C ----------------------------------------------------------------------
Ci Inputs
Ci  n:    number of independent variables
Ci  hess: inverse Hessian matrix for Fletcher-Powell and Broyden
Ci  xtoll:x tolerance for line minimizations.  See Remarks.
Ci  dxmx: Maximum step length in any one component of x
Ci        Also, for C.G., initial step length; see Remarks.
Ci        A positive sign of dxmx makes line mins tend 'downhill'
Ci        A negative sign of dxmx makes line mins tend 'uphill'
Ci  xtol: tolerance in delta x for global minimization.
Ci        xtol > 0 : tolerance compared to |delta x|
Ci        xtol < 0 : tolerance compared largest component, delta max(|x|)
Ci  gtol: tolerance in gradient for global minimization
Ci        gtol > 0 : tolerance compared to |delta g|
Ci        gtol < 0 : tolerance compared largest component, delta max(|g|)
Ci  grfac:Extrapolation growth factor for line minimizations.
Ci        Whenever a root is not bracketed, the effective value of dxmx
Ci        is scaled by grfac.
Ci  isw:  compound of one-digit switches.
Ci    1's digit: handles constraints and step directions
Ci       *For the line minimization algorithms, digit is passed through
Ci        to the line min. routine rfalsi (which see).
Ci      0 no constraints imposed
Ci      1 require line min grad >0 at root (minimum for each line min)
Ci      2 require line min grad <0 at root (maximum for each line min)
Ci      4 while root not bracketed, set suggested x = current x + dxmx
Ci       *For the Broyden scheme, the digit is passed through to routine
Ci        brmin; the meaning is somewhat different (see brmin.f)
Ci        See Remarks for some further description.
Ci   10's digit: convergence criteria; see also Remarks.
Ci      0 convergence when <g> < gtol
Ci      1 convergence when <dx> < xtol
Ci      2 convergence when <g> < gtol and <dx> < xtol
Ci      3 convergence when <g> < gtol or <dx> < xtol
Ci        Here, <g> and <dx> are the max component of gradient and change
Ci        in x between iterations
Ci      4 Set this bit to specify a line minimization tolerance in the
Ci        gradient; see Remarks.
Ci  100's digit: governs which minimization  to use:
Ci      0 Congugate gradients, adapted from Numerical Recipes, 10.6.
Ci      1 Fletcher-Powell: M.J.Norgett & R.Fletcher, J.Phys.C3,L190(1972)
Ci      2 Broyden
Ci 1000's digit: Hessian matrix switches
Ci      0  no tests are made for Hessian
Ci      1  If Hessian is not positive definite, it is not updated.
Ci      2  project out parts of Hessian corresponding to negative
Ci         eval (not implemented).
Ci      4  On the second iteration, Hessian is globally scaled
Ci         once a first estimate is known (D. Novikov)
Ci 10000's digit: line-minimization-specific switches.
Ci      1 specifies tolerance for change in line minimization gradient;
Ci        see Remarks.
Cio Inputs/Outputs
Cio ir:   program flow control.  To start a new minimization,
Cio       set ir to zero.  On exit, gradzr sets ir to one of the
Cio       following.  For ir<0, gradzr expects a new set of points
Cio       and gradients at positions p(1..n).
Cio      *For the line minimization schemes F.P. and C.G:
Cio    >=0: no more function calls expected:
Cio      0: gradzr finds function converged to within tolerance
Cio      1: input grad equals 0  gradzr does nothing.
Cio      2: nonsensical tolerance.  gradzr does nothing.
Cio      3: input xn equals x1 or x2.  gradzr does nothing.
Cio      4: line min extremum bracketed but with no root.  gradzr returns
Cio         suggested xn = largest (smallest) of given x so far + dxmx
Cio      5: line minimization failed to bracket root;
Cio         see 100's digit of isw
Cio     <0: expects new function call, suggesting xn:
Cio     -1: no information yet: require new point at suggested xn
Cio         to improve the estimate of the root.
Cio     -2: a root is bracketed and gradzr is attempting linear
Cio         interpolation to improve the estimate of the root.
Cio     -3: a root is bracketed and gradzr is attempting quadratic
Cio         interpolation to improve the estimate of the root.
Cio     -4: a root has not been bracketed or constraints not fufilled.
Cio         gradzr will suggests xn from a linear extrapolation.
Cio     -5: same as -4 above, but gradzr tries xn = current x + dmmx
Cio         is increasing along the step direction.  If ir is set
Cio         to -4 for a new point, gradzr will continue along that
Cio         step direction
Cio     -6: gradzr is having trouble; asks for new point
Cio     -9: Line minimization converged; new line minization sought.
Cio         For now, ir=-9 used just internally.
Cio         NB: gradzr is organized so that this need not happen
Cio    -10: ir input -10 works just as input ir=0, but flags that
Cio         the conjugate direction p(nd) is already specified by
Cio         the caller.  Valid only for line minimization schemes.
Cio      *For Broyden (meaning same as in brmin.f but sign has changed)
Cio      0: brmin has converged to specified tolerance
Cio     <0: brmin needs gradients (preferably at the output x)
Cio     >0: like <0, but the Hessian was found to be not positive
Cio         definite, and was not updated.
Cio p:    array of dimension at least 6*n for F.P. and C.G. and
Cio       8*n for Broyden.  If the Hessian eigenvalues are calculated,
Cio       dimensions must be n*(6+n) and n*(11+2n) for F.P. and Broyden
Cio       p(*,nx=1) current position vector. (Input on first call)
Cio                 Internally, p(nx) = p (start of current line min)
Cio       p(*,ng=2) gradient at current point (Input on each call)
Cio      *For line minimization schemes:
Cio       p(*,nm=3) position vector where gradient was found to be minimum
Cio                 data is kept for informational purposes only
Cio       p(*,nd=4) conjugate direction for current line min.
Cio       p(*,n1=5) gradient at the point prior to this one
Cio       p(*,n2=6) gradient at second prior point
Cio       p(*,n0=7) gradient at start of line min (C.G. only)
Cio      *For Broyden, p(*,nd...8) are passed as w(*,1..6); see brmin.f
Cio wk    is a work array of dimension 0:27.  It should not be altered
Cio       between succesive calls.  The elements of wk are:
Cio wk(0) = xn is a measure of the shift in coordinates
Cio       between successive function calls.
Cio      *For line minimizations, xn is the amount of direction vector
Cio       that was added to coordinate positions; see Remarks.
Cio      *For Broyden, xn is the largest shift in a component of x
Cio wk(1..12) see rfalsi.
Cio wk(13) gam = gamma for conjugate gradients; cf Numerical Recipes.
Cio wk(14) dxmxl = maximum step length along direction vector
Cio wk(15) dxmxlx = step size that makes shift in largest p(nx)= dxmx
Cio wk(16) dxlast = Largest change in one component of x;
Cio          used as estimate for dx in new C.G. line minimization.
Cio wk(17) xmax = largest value of x along a given line min; used
Cio          when line minimizations have trouble finding a root
Cio          satisifying constraints.
Cio wk(18) lminn = number of line minimizations so far
Cio wk(19) grad0 = (grad.dir-vec) at start of line minimization
Cio wk(20) growl = Current extrapolation growth factor
Cio wk(21..23) x0h,x1h,x2h values of xn in prior iterations
Cio wk(24..26) g0h,g1h,g2h values of grad.dir-vec corresponding to x*h
Cio wk(27) global minimum gradient, corrsponding to p(nm)
Cio        See 10's digit isw for criterion that measures gradient
Cio        sign(wk(27)) used as a flag to indicate:
Cio        >0 most recent |g| is also global minimum
Cio        <0 most recent |g| is not global minimum
Cio           In this case, p(*,nm) contains positions for global min
Cr Remarks
Cr *gradzr attempts to find a a minimum, a maximum,
Cr  or just a zero in the gradient using one of three schemes:
Cr  Broyden, Fletcher-Powell, and congugate gradients.  There are
Cr  strengths and weaknesses in each approach.  The F.P and C.G. proceed
Cr  by successive line minimizations, in which the zero projection of the
Cr  gradient along each direction vector is sought.  In the program,
Cr  p(*,nd) is the direction vector, and p(*,nx) are the coordinates at
Cr  the start of the line minimization.  xn is varied so that the for
Cr  positions p(nx) + xn * p(nd), grad.p(nd) is zero to within a specified
Cr  tolerance, after which the inverse Hessian matrix is updated and new
Cr  conjugate direction is calculated.  The F.P Hessian is explicit; in
Cr  the C.G. case it is implicit.  Having the Hessian is useful
Cr  but can be expensive if there are many variables.  In practice, F.P.
Cr  seems to converge a little faster than C.G..  Broyden is essentially a
Cr  Newton-Raphson scheme for several variables, and the Hessian is
Cr  explicit.  In well conditioned cases, it tends to get to the root more
Cr  rapidly than the line minimization approaches, except that it is more
Cr  sensitive to roundoff errors, which after a few iterations cause it
Cr  to converge more slowly.
Cr
Cr *gradzr works by accepting a new vector x and its gradient g, in
Cr  p(1..n) and p(n+1..2n) respectively.  It returns with a new
Cr  suggested value for p, for which the caller is expected to generate
Cr  g.  (For the line minimization routines, it is more natural for
Cr  gradzr to return a direction vector and the distance along the
Cr  projection vector, but this is not done here to make the calling
Cr  interface as consistent as possible.  The distance xn is returned in
Cr  case the caller wants to change it.)  Variable ir supplies some
Cr  information as to what gradzr is looking for in the next iteration.
Cr
Cr *Convergence tolerances.
Cr  gradzr continues its minimization procedure until it satisfies some
Cr  combination of:
Cr    (the change in position between iterations) < |xtol|
Cr    (the value  of the gradient)                < |gtol|
Cr  The l.h.s is evaluated either as the largest component
Cr  (xtol>0 or gtol>0) or as the length of the vector (xtol<0, gtol<0)
Cr  The 10's digit of isw specifies which combination is selected.
Cr
Cr  For Broyden, the change in x is the change between each function call;
Cr  for the line minimization algorithms C.G. and F.P., it is the change
Cr  between each line minimization.  In this latter case, there is a
Cr  second, independent tolerance for the individual line minimizations,
Cr  xtoll.  There is also a separate convergence criterion for gradient in
Cr  the individual line minimizations, namely the line minimization stops
Cr  when (grad.dir-vec) has dropped to a specified fraction f of its
Cr  starting value.  Turn this option on by adding 40 to isw.  You can
Cr  specify f using i=10000's digit of isw; f=i/20.
Cr  Using i=0 sets f to a default of 0.25.  NB: when no x tolerance is
Cr  specified, the line gradient tolerance is automatically turned on.
Cr
Cr *Initial step lengths for a new line minimization: when the hessian
Cr  is explicit, it specifies the step length.  You may specify a
Cr  maximum step length dxmx, which imposes a maximum value for a change
Cr  in any component of x.  The C.G. uses dxmx to specify the initial
Cr  step size for the first line minimization.  For subsequent line
Cr  mins, the initial step length is taken to be the total change in
Cr  position from the prior line min.
Cr
Cr *gradzr can be used to seek a minimum, a maximum,
Cr  or just a zero in the gradient; see one's digit of isw.
Cr  When a constraint is imposed, caller is strongly advised to
Cr  also set the 4's digit of isw, which will cause rfalsi to continue
Cr  in the same search direction even if the gradient increases.
Cr  Also, if a maximum is sought (isw=2), caller is strongly advised
Cr  to set dxmx negative.  Otherwise, strange behavior may result.
Cr  For now these switches are automatically set internally.
Cr
Cr *recommended values for line minimizations.
Cr  At least in well behaved cases, it is recommended that
Cr  the caller choose an x tolerance for the line minimizations.
Cr  A 'safe' default for the line tolerance xtoll is the same value
Cr  as the global xtol, though gradzr tends to convergence a little
Cr  faster if it is set somewhat larger.
Cr  Also, convergence seems to be aided by changing line minimizations
Cr  when the gradient has dropped by ~0.25 (add 40 to isw).
Cr
Cr *It is usually simpler and more convenient to call gradzr using
Cr  the driver routine drgrzr.
Cr
Cu  Updates
Cu   08 Mar 06 some improvements to convergence criteria; saves p(gmin)
Cu   07 Sep 03 gradzr returns after first iter
Cu             when g<gtol is satisfied and is sole conv. criterion
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer n,ir,isw
      double precision hess(n,n),p(n,10),xtoll,dxmx,xtol,gtol,wk(0:*),
     .grfac
C Local Variables
      logical ltmp,cnvgx,cnvgg,nwlin
      integer i,j,bitand,ipr,npr,lminn,isw0,isw1,isw2,isw3,isw4,
     .isw31,idamax,iswl,itol,iswb
      integer nx,ng,nd,n1,n2,n0,nm,nhs
      double precision gg,dgg,gam,ddot,dasum,x0h,x1h,x2h,evtol,
     .evmin,dxmxl,dxmxlx,xtll,dxmn,xmax,gfn,grad0,g0h,g1h,g2h,
     .dxm,growl,gammax,xh(0:2),gh(0:2),xn,
     .gd1,gd2,dg,alpha,r,s,gll
      equivalence (xh(0),x0h), (xh(1),x1h), (xh(2),x2h)
      equivalence (gh(0),g0h), (gh(1),g1h), (gh(2),g2h)
      parameter (nx=1,ng=2,nm=3,nd=4,n1=5,n2=6,n0=7,nhs=8, evtol=1d-10)
      double precision gwg,q,dxlast,dxtop,gtop,resx,resg,dum
      integer scrwid
      parameter (scrwid=80)
      character*(100) outs

C --- Iteration independent setup ---
      if (ir .eq. 0) call dpzero(wk,27+1)
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
      isw4  = mod(isw/10000,10)
      itol = mod(isw1,4)
C     Save min gradient, p(nm) using initial p and gradient
      if (ir .eq. 0) then
        call pgradz(0,n,0d0,p,p(1,ng),0d0,0d0,
     .  dxtop,resx,gtop,resg,wk(27),cnvgx,cnvgg)
        call dcopy(n,p(1,nx),1,p(1,nm),1)
      endif

C ... Update p if gradient to minimum gradient
      call pgradz(0,n,0d0,p,p(1,ng),0d0,0d0,
     .dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
      if (r .lt. abs(wk(27))) then
        wk(27) = r
        call dcopy(n,p(1,nx),1,p(1,nm),1)
      else if (r .gt. abs(wk(27))) then
        wk(27) = -abs(wk(27))
      endif

C ... These lines impose sanity for conditions described in Remarks
      if (mod(isw0,4) .ne. 0) isw0 = mod(isw0,4) + 4
      dxm = dxmx
      j = mod(isw0,4)
      if (dxmx.gt.0 .and. j.eq.2 .or. dxmx.lt.0 .and. j.eq.1)
     .dxm = -dxmx

C ... Re-entry point for new line minimization
    1 continue
C ... If the gradient is zero, nothing to calculate.
      gg = abs(p(idamax(n,p(1,ng),1),ng))
      if (gg .eq. 0) then
        ir = 0
        goto 999
      endif

C --- First iteration ---
      if (ir .eq. 0 .or. (ir .eq. -10 .and. isw2 .le. 1)) then
C   ... This the largest dx from the last iteration
C       For first line min, set it so largest dx is dxm
        dxlast = dxm
C   ... Number of line minimizations
        if (ir .eq. 0) lminn = 0
C   ... Congugate gradients gamma (see Numerical Recipes)
        gam = 0
C   ... Size of initial gradient (and direction) vector
        j = mod(isw0,4)
        ltmp = j .eq. isw0 .or.
     .  dxm.gt.0 .and. j.eq.2 .or. dxm.lt.0 .and. j.eq.1
        if (ipr .ge. 30 .or. ipr .ge. 20 .and. ltmp) then
c$$$          call awrit7(' gradzr: begin %?#n# xtol=|%1;3g|#%j#'
c$$$     .    //'%?#n==2# and##%-1j%?#n==3# or##%-1j'
c$$$     .    //'%?#n<>1# gtol=|%1;3g|#%j#'
c$$$     .    //'%?#n==0|n>=4#  gtll=%1;3g#%j#'
c$$$     .    //'  dxmx=%1;3g',outs,
c$$$     .    scrwid,0,itol,xtol,itol,gtol,
c$$$     .    isw1,pgrada(isw1,isw4),dxm)
           write(stdo,ftox)'itol,xtol,gtol,isw1,pgrada(isw1,isw4),dxm',itol,xtol,gtol,isw1,pgrada(isw1,isw4),dxm
c          call info2(-25,0,0,
c     .    outs//'%a  '//
c     .    '%?#n==0#C.G.##%-1j'//
c     .    '%?#n==1#F.P.##%-1j'//
c     .    '%?#n==2#Broy##'//
c     .    '  isw=%i',isw2,isw)
          write(stdo,ftox)'cg fp broy',isw2,isw
C     ... Sanity checks
C          if (j .eq. isw0)
C     .    call awrit0(' gradzr (warning) seek extremum but 4''s'
C     .      //' bit of isw is not set',' ',scrwid,i1mach(2))
C          if (dxm.gt.0 .and. j.eq.2 .or. dxm.lt.0 .and. j.eq.1)
C     .    call awrit0(' gradzr (warning) seek extremum but dxmx'//
C     .      ' is of the wrong sign',' ',scrwid,i1mach(2))
          if (bitand(isw1,1) .eq. 0 .and. gtol .eq. 0) ir = 2
          if (bitand(isw1+3,2) .eq. 0 .and. xtol .eq. 0) ir = 2
          if (xtol .eq. 0 .and. gtol .eq. 0) ir = 2
C         if (ir .eq. 2) call rx('OOPS')
          if (ir .eq. 2) goto 999
        endif
      endif

C --- Broyden minimization ---
      if (isw2 .eq. 2) then
        ir = -ir
        iswb = isw0 + 10*bitand(isw1,3) + 1000*isw3
        call brmin(n,p,p(1,ng),iswb,ipr,dxmx,xtol,gtol,dum,p(1,nd),xn,
     .  hess,ir)
        ir = -ir
        goto 999
      endif

C ... Case continue with current line minimization:
      if (ir .ge. -6 .and. ir .le. -1) then
C       Undo direction vector added to position vector,
C       restoring position vector to starting point of current line min
        if (ir .lt. 0) call daxpy(n,-xn,p(1,nd),1,p(1,nx),1)
        goto 20
      endif

C --- New line minimization ---
      lminn = lminn+1
C ... Congugate gradients
      if (isw2 .eq. 0) then
        gg = ddot(n,p(1,ng),1,p(1,ng),1)
        if (gam .eq. 0 .and. ir .ne. -10) call dpzero(p(1,nd),n)
        dgg = ddot(n,p(1,ng),1,p(1,nd),1)
        resg = gam*dgg-gg
        gammax = 2
        if (resg .gt. 0) then
          write(stdo,ftox) ' gradzr encountered g.h ='
     .    //'  ... set gam to zero',ftod(resg)
          gam = 0
        elseif (gam .gt. gammax) then
          write(stdo,ftox)' gradzr cap gam = ',' to max gam =',ftod(gam),ftod(gammax)
          gam = gammax
        endif
        call dpcopy(p(1,ng),p(1,n0),1,n,-1d0)
        if (ir .eq. -10) then
          ir = 0
        else
          do    j = 1, n
             p(j,nd) = p(j,n0) + gam*p(j,nd)
          enddo   
        endif
        dxmxl = sign(dxlast/p(idamax(n,p(1,nd),1),nd),dxm)
C ... Fletcher-Powell
      elseif (isw2 .eq. 1) then
C       If starting Hessian is zero, set to unity.
        if (dasum(n*n,hess,1) .eq. 0) call dcopy(n,1d0,0,hess,n+1)
C       F.P doesn't use p(n0), but save for compatibility with C. G.
        call dpcopy(p(1,ng),p(1,n0),1,n,-1d0)
C       New positions generated from inverse Hessian matrix
C       Shift in conjugate dir p(nd) = -hess . grad
        if (ir .eq. -10) then
          ir = 0
        else
          do  6  i = 1, n
            p(i,nd) = 0
            do    j = 1, n
               p(i,nd) = p(i,nd) - hess(i,j)*p(j,ng)
            enddo   
    6     continue
          dxmxl = 1
        endif
      else
        call rxi('gradzr not ready for isw2=',isw2)
      endif

C ... This is the 'best guess' for step size, new line min.
      dxmxl = dxmxl*min(1d0,dabs(dxm/(dxmxl*p(idamax(n,p(1,nd),1),nd))))
C ... Step size that makes shift in largest component = dxmx:
      dxmxlx = sign(dxmx/p(idamax(n,p(1,nd),1),nd),dxm)
C ... Initial extrapolation growth factor
      growl = 1
C ... (grad.dir-vec) at start of line minimization
      grad0 = ddot(n,p(1,ng),1,p(1,nd),1)
      if (ir .eq. -9 .or. ir .eq. 0) then
        ir = 0
        xn = 0d0
      endif
      xmax = 0d0
      x1h = -9d9
      x2h = -9d9
C ... Printout of current p,g,h
      if (ipr .gt. 30) then
        gg = ddot(n,p(1,ng),1,p(1,ng),1)
        dgg = ddot(n,p(1,ng),1,p(1,nd),1)
        resg = dsqrt(ddot(n,p(1,ng),1,p(1,ng),1))
        gtop = p(idamax(n,p(1,ng),1),ng)
        call info5(-30,0,0,' gradzr new line %i:  g.h=%;4g'//
     .  '  g.(h-g)=%;4g  max g=%1;3g  |grad|=%1;3g  ',
     .  lminn,dgg,dgg+gg,gtop,resg)
      endif
      call info2(-40,0,0,'  p=%n:;10F',npr,p(1,nx))
      call info2(-40,0,0,'  g=%n:;10F',npr,p(1,ng))
      call info2(-40,0,0,'  h=%n:;10F',npr,p(1,nd))

C --- Begin or continue current line minimization ---
   20 continue
      xmax = max(xmax,xn)
      x0h = xn
      gfn = ddot(n,p(1,ng),1,p(1,nd),1)
      g0h = dsqrt(ddot(n,p(1,ng),1,p(1,ng),1))

C ... Decide on global convergence before calling rfalsi,
C     since criteria might not be identical to that in rfalsi

C     Global convergence criteria: shift relative to last position
C     call pgradz(itol,n,xn-x0h,p(1,nd),p(1,ng),
C     .  dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
C     Global convergence criteria: shift relative to start of line min
      nwlin = ir .eq. 0
      call pgradz(itol,n,xn,p(1,nd),p(1,ng),xtol,gtol,
     .dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
C     Case gradzr has converged globally
      if (resg .eq. 0 .or. (cnvgg .and. cnvgx)) goto 998

C ... New line min if g < line g tol
      if (dabs(gfn/grad0) .lt. pgrada(isw1,isw4)) then
        ir = 0

C ... Otherwise, call rfalsi to check for convergence and/or new xn
      else
C       Requires rfalsi to return xn on an existing point when ir=0
        iswl = isw0 + 40
C       Also suppress rfalsi using gtol; we do that above
        if (mod(isw1,4) .ne. 0) iswl = iswl + 10
C       Tolerance should be less than and distinct from dxmx
        xtll = min(xtoll/dabs(p(idamax(n,p(1,nd),1),nd)),
     .  abs(dxmxl*growl*.999999d0))
C       Minimimum step size should be less than, distinct from tolerance
        dxmn=xtll/2
        gll = abs(gfn/10)
        call pshpr(ipr-10)
        call rfalsi(xn,gfn,xtll,gll,dxmn,dxmxl*growl,iswl,wk(1),ir)
        call poppr
C       Exit if 1st step of new line satisfies 'global convergence'
        if (nwlin) then
          call pgradz(itol,n,xn,p(1,nd),p(1,ng),xtol,gtol,
     .    dxtop,resx,gtop,resg,r,cnvgx,cnvgg)
C         Case gradzr has converged globally
          if (resg .eq. 0 .or. (cnvgg .and. cnvgx)) then
            xn = 0
            goto 998
          endif
        endif
C   ... On the first line min movement, reset dxmx to maximum allowed
        if (mod(iswl,10) .eq. 0 .and. nwlin) dxmxl = dxmxlx
C   ... rfalsi wants to extrapolate
        if (ir .ge. -6 .and. ir .le. -4) then
C     ... This accelerates linear extrapolation
          if (ir .eq. -4) xn = x0h + (xn-x0h)*growl
          growl = growl * grfac
C        elseif (ir .eq. -1) then
C          growl = grfac
        else
          growl = 1
        endif
C   ... Swap points in the same way rfalsi swapped them.
        if (nint(wk(12)) .ge. 4) then
          call dswap(n,p(1,n1),1,p(1,n2),1)
          call dswap(1,xh(1),1,xh(2),1)
          call dswap(1,gh(1),1,gh(2),1)
        endif
        if (mod(nint(wk(12)),4) .ge. 2) then
          call dswap(n,p(1,ng),1,p(1,n2),1)
          call dswap(1,xh(0),1,xh(2),1)
          call dswap(1,gh(0),1,gh(2),1)
        endif
        if (mod(nint(wk(12)),2) .ge. 1) then
          call dswap(n,p(1,ng),1,p(1,n1),1)
          call dswap(1,xh(0),1,xh(1),1)
          call dswap(1,gh(0),1,gh(1),1)
        endif
C ... Handle special case ir=0, xn=0
        if (ir .eq. 0 .and. xn .eq. 0) then
          ir = -2
          xn = (wk(5)*wk(1)-wk(4)*wk(2))/(wk(5)-wk(4))
        endif
      endif

C ... Case line minimization not converged
      if (ir .lt. 0) then

C   ... After clean, replace xh with wk.  for now:
        if (x0h .ne. wk(1) .or. x1h .ne. wk(2) .and. x1h .ne. -9d9
     .  .or. x2h .ne. wk(3) .and. x2h .ne. -9d9) then
          print *, x0h,wk(1)
          print *, x1h,wk(2)
          print *, x2h,wk(3)
          call rx('bug in gradzr or rfalsi')
        endif

C   ... preserve gradient for last two points
        x2h = x1h
        g2h = g1h
        call dcopy(n,p(1,n1),1,p(1,n2),1)
        x1h = x0h
        g1h = g0h
        call dcopy(n,p(1,ng),1,p(1,n1),1)

      elseif (ir .eq. 4) then
        xn = xmax + dxmxl
        call info2(-10,0,0,'%x gradzr: found extremum without root:'//
     .  '  attempt new xn=%1;4g',xn,0)
        ir = -5
        goto 991

      elseif (ir .ne. 0) then
        call info2(-1,0,0,' gradzr (abort) line %i encountered ir=%i '//
     .  'from rfalsi',lminn,ir)
        goto 991

C --- Line minimization has converged ---
      else
C       Sanity check
        if (xn .ne. x0h) call rx('bug in gradzr')
        if (isw2 .eq. 0) then
C     ... Make gam for this line
          gg = 0d0
          dgg = 0d0
          do  21  j = 1, n
            gg = gg + p(j,n0)**2
C           Use the following line for Fletcher-Reeves:
C           dggfp = dgg + p(j,ng)**2
C           or the following line for Polak-Ribiere:
            dgg = dgg + (p(j,ng)+p(j,n0))*p(j,ng)
   21     continue
          if (gg .eq. 0d0) then
          else
            gam = dgg/gg
C           gamfp = dggfp/gg
          endif
        elseif (isw2 .eq. 1) then

C         Local copy of the original Hessian, in case used later
          if (isw31 .ne. 0) call dcopy(n*n,hess,1,p(1,nhs),1)

C     ... Update inverse of Hessian
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
   33     continue
          q = 1 + gwg/dg
          gd1 = 0
          do  36  i = 1, n
            r = - alpha*p(i,nd)/dg
            s = (-p(i,n2) + q * alpha*p(i,nd))/dg
            do   j = 1, n
               hess(i,j) = hess(i,j) + r * p(j,n2) + s * alpha*p(j,nd)
            enddo   
   36     continue

C     ... Check that Hessian is positive definite
          if (isw31 .ne. 0) then
C           Keep a local copy of hessian, since dsev1 destroys it.
            call dcopy(n*n,hess,1,p(1,nhs+n),1)
            j =chkhss(p(1,nhs+n),n,p(1,n2),evtol,isw31,p(1,nhs),p(1,n1))
            evmin = p(1,n1)
            if (j .gt. 0 .and. (isw31 .eq. 1)) then
              call dcopy(n*n,p(1,nhs),1,hess,1)
              ir = -ir
            elseif (j .gt. 0) then
            endif
          endif

        endif

C       Internally flag to start new line min
        ir = -9
C       Restore position vector to current position
        call daxpy(n,xn,p(1,nd),1,p(1,nx),1)
C       Largest change in component of x relative to start of line min
C       Note: used by C.G. line min.
        dxlast = xn*p(idamax(n,p(1,nd),1),nd)

C       Printout for line min convergence
        call info8(-30,0,0,'%x gradzr cvg line %i:'//
     .  '%?#n==0#  gam=%1;3g#%j#  x=%1;8d'//
     .  '  |g.h|=%1;3g  dxmax=%1;3g'//
     .  '%?#n#  evmin=%1;2g#%j#',
     .  lminn,isw2,gam,xn,gfn,dxlast,isw31,evmin)
      endif

C --- Cleanup ---
  991 continue
C     Add direction vector to position vector
      if (ir .ge. -6 .and. ir .le. -1)
     .call daxpy(n,xn,p(1,nd),1,p(1,nx),1)
      if (ir .eq. -9) goto 1

C --- Restore local variables needed to preserve; exit ---
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

C     Use p(minimum g), if specified and different from p
C      if (isw5 .eq. 1 .and. ir .ge. 0) then
C        call info2(-20,0,0,'restore vector for minimum g, g=%1;3g',
C     .    wk(27),0)
C        call dcopy(n,p(1,nm),1,p,1)
C      endif
      return

C --- Global convergence FP or CG ---
  998 continue
      call info5(-20,0,0,
     .'%x gradzr converged to dxmax=%1;3g, |dx|=%1;3g,'//
     .'  gmax=%1;3g, |grad|=%1;3g in %i line min',
     .dxtop,resx,gtop,resg,lminn)
      call info2(-30,0,0,' p=%n:;10F',npr,p(1,nx))
      call info2(-30,0,0,' g=%n:;10F',npr,p(1,ng))

C     Restore position vector to current position
      call daxpy(n,xn,p(1,nd),1,p(1,nx),1)

C     Clean up and exit
      ir = 0
      goto 999
      end
      subroutine drgrzr(n,pnew,gnew,p,hess,xtoll,dxmx,xtol,gtol,grfac,
     .wk,copt,isw,ir)
C- Driver routine for gradzr
C ----------------------------------------------------------------------
Ci Inputs
Ci   n      :size of vector; see gradzr
Ci   pnew   :new positions for which gradient was obtained
Ci   gnew   :gradient corresponding to new positions
Ci  The following inputs are passed directly through to gradzr.
Ci  See that routine for further description.
Ci   hess   :inverse Hessian matrix for Fletcher-Powell and Broyden
Ci   xtoll  :x tolerance for line minimizations.  See gradzr
Ci   dxmx   :Maximum step length in any one component of x.
Ci   xtol   :tolerance in x for global minimization.
Ci   gtol   :tolerance in gradient for global minimization
Ci   grfac  :Extrapolation growth factor for line minimizations.
Ci   copt   :character string for a convenient specification of special
Ci           options.  drgrzr converts the following strings in copt
Ci           into the corresponding switches in isw.  These may be
Ci           strung together, separated by spaces.  isw is altered
Ci           only when input ir=0.
Ci             'def' Set default isw
Ci             'cg'  congugate gradients
Ci             'fp'  Fletcher-Powell
Ci             'br'  Broyden
Ci             'min' specify minimization
Ci             'max' specify maximization
Ci           Options are read left-to-right, so when incompatible
Ci           switches are set, the last takes precedence.
Ci   isw    :compound of one-digit switches.  It will be altered
Ci           if copt is set.
Co Inputs/Outputs
Cio  ir     :flow control passed to gradzr.
Cio          To start a new minimization, set ir to zero.
Cio          On exit, gradzr sets ir to a value as described in gradzr.
Cio          A return of ir=0 => convergence achieved to prescribed tol
Cio          A return of ir>0 => gradzr had trouble and is aborting.
Cio  p      :position and work array for input to gradzr.
Cio          Array is dimensioned at least 6*n for F.P. and C.G. and 8*n
Cio          for Broyden.  If the Hessian eigenvalues are calculated,
Cio          it increases to n*(6+n) and n*(11+2n) for F.P. and Broyden.
Cio          p should remain untouched between successive calls
Cio          drgrzr replaces p(*,nx):replaced by pnew and
Cio          p(*,nd) changed to p(nd) = xn * (h(xn)-h(xn=0))
Cio          NB: input ir=0 => p(nd) not defined and not changed
Cio  wk     :work array passed for input to gradzr.
Cio         :wk should remain untouched between successive calls
Cr Remarks
Cr   This is a driver routine for gradzr, taking as input positions
Cr   and gradients, and internally updating the gradzr matrix p.
Cr   For line minimizations, it also resets the conjugate direction
Cr   vector h = p(1,nd) in the event the positions passed to drgrzr
Cr   do not correspond to those kept by p.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) copt
      integer n,ir,isw
      double precision pnew(n),gnew(n)
      double precision hess(n,n),p(n,10),xtoll,dxmx,xtol,gtol,wk(0:26),
     .grfac
C Local variables
      integer idx,nx,ng,nd,n0,idamax,j1,j2,getdig
      double precision dhmax,hold,hmax,tol,xn
      parameter (nx=1,ng=2,nd=3,n0=6,tol=1d-4)

C --- Set switches based on copt ---
      if (ir .eq. 0 .and. copt .ne. ' ') then
        j1 = 1
   10   continue
        call nword(copt,1,j1,j2)
        if (j2 .ge. j1) then
          if (copt(j1:j2) .eq. 'cg') then
            isw = isw + 100*(0-getdig(isw,2,10))
          elseif (copt(j1:j2) .eq. 'fp') then
            isw = isw + 100*(1-getdig(isw,2,10))
          elseif (copt(j1:j2) .eq. 'br') then
            isw = isw + 100*(2-getdig(isw,2,10))
          elseif (copt(j1:j2) .eq. 'def') then
            isw = 40
          elseif (copt(j1:j2) .eq. 'min') then
            isw = isw + 1*(1-mod(getdig(isw,0,10),4))
          elseif (copt(j1:j2) .eq. 'max') then
            isw = isw + 1*(2-mod(getdig(isw,0,10),4))
          else
            call rxs2('drgrzr:  bad option, "',copt(j1:j2),'"')
          endif
          j1 = j2+1
          goto 10
        endif
      endif

      xn = wk(0)
C ... hnew = hold + 1/xn (pnew - pold)
      if (ir .ne. 0 .and. xn .ne. 0d0) then
        idx = idamax(n,p(1,nd),1)
        hmax = abs(p(idx,nd))
C       This destroys p(nx) but it isn't needed anymore
        call daxpy(n,-1d0,pnew,1,p(1,nx),1)
        idx = idamax(n,p(1,nx),1)
        dhmax = abs(p(idx,nx)/xn)
        hold = abs(p(idx,nd))
        call daxpy(n,-1d0/xn,p(1,nx),1,p(1,nd),1)
C   ... Change in direction exceeds tolerance; reset ir
        if (dhmax .gt. tol*hmax) then
          call info2(-30,0,0,' gradzr: reset conjugate'//
     .    ' direction.  dhmax = %1,3;3g  hold = %1,3;3g',dhmax,hold)

C     ... For ir=-10, need p(nx) = p(xn=0), p(ng) = -p(n0)
          call dcopy(n,pnew,1,p(1,nx),1)
          call daxpy(n,-xn,p(1,nd),1,p(1,nx),1)
          call dpcopy(p(1,n0),p(1,ng),1,n,-1d0)

C     ... Call gradzr to reset C.G. or F.P. for new line min
          ir = -10
          call gradzr(n,p,hess,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
          wk(0) = xn

        endif
      endif

C ... Set p(nx) to pnew and p(ng) to gnew
      call dcopy(n,pnew,1,p(1,nx),1)
      call dcopy(n,gnew,1,p(1,ng),1)

C ... Call gradzr for next step in minimization
      call gradzr(n,p,hess,xtoll,dxmx,xtol,gtol,grfac,wk,isw,ir)
      end

      double precision function pgrada(isw1,isw4)
C- Line minimization gradient tolerance
C     implicit none
      integer isw1,isw4
      double precision gtll

      gtll = 0
C          requested   req'd when xtol not used
      if (isw1 .ge. 4 .or. mod(isw1,4) .eq. 0) then
        gtll = dble(isw4)/20
C       If zero, use a default value
        if (gtll .eq. 0) gtll = .25d0
      endif

      pgrada = gtll
      end

      subroutine pgradz(itol,n,amp,p,g,xtol,gtol,dxmx,resx,dgmx,resg,
     .gx,cnvgx,cnvgg)
C- Evaluates convergence criteria in x and g
C ----------------------------------------------------------------------
Ci Inputs
Ci   itol  :flags which criterion to use
Ci         :0 => convergence when <g> < gtol
Ci         :1 => convergence when <p> < xtol and amp ne 0
Ci         :2 => convergence when <g> < gtol and (<p> < xtol and amp ne 0)
Ci         :3 => convergence when <g> < gtol or  (<p> < xtol and amp ne 0)
Ci   n     :number of elements
Ci   amp   :scaling factor for p
Ci   p     :unscaled vector for position shift
Ci   g     :gradient vector
Co Outputs
Co   dxmx  :largest change in single component of p
Co   resx  :change in length of p
Co   dgmx  :largest change in single component of g
Co   resg  :change in length of g
Co   gx    :whichever of dgmx or resg used pgradz uses as measure of g
Co         :gx is always positive
Co   cnvgx :T if convergence condition met in x
Co   cnvgg :T if convergence condition met in g
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   08 Mar 06 adapted from gradzr
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer itol,n
      double precision amp,p(n),g(n),xtol,gtol,dxmx,resx,dgmx,resg,gx
      logical cnvgx,cnvgg
C ... Local parameters
      integer idamax
      double precision ddot

C     Largest component of p
      dxmx = amp*p(idamax(n,p,1))
C     Change in length of p
      resx = amp*dsqrt(ddot(n,p,1,p,1))
C     Largest component of g
      dgmx = g(idamax(n,g,1))
C     Residual of gradient at new position and printout
      resg = dsqrt(ddot(n,g,1,g,1))

C     Find cnvgx
      if (xtol .ge. 0) gx = abs(resx)
      if (xtol .lt. 0) gx = abs(dxmx)
      cnvgx = gx .le. abs(xtol)
      if (amp .eq. 0)  cnvgx = .false.
      if (itol .eq. 0) cnvgx = .true.

C     Find cnvgg
      if (gtol .ge. 0) gx = abs(resg)
      if (gtol .lt. 0) gx = abs(dgmx)
      cnvgg = abs(gx) .le. abs(gtol)
      if (itol .eq. 1) cnvgg = .true.

C     If either is permissible and one is satisfied, set both true
      if (itol .eq. 3 .and. (cnvgx .or. cnvgg)) then
        cnvgx = .true.
        cnvgg = .true.
      endif
      end

      integer function chkhss(hess,n,w5,evtol,iopt,z,e)
C- Returns the number of eigenvalues greater than a specified value
Co z: eigenvectors if iopt > 1; otherwise, z unchanged
Co e: eigenvalues
Co hess: DESTROYED on output
C     implicit none
      integer n,iopt
      double precision hess(n,n),w5(n,5),evtol,z(n,n),e(n)
      double precision dum(1,1)
      integer i,j

      j = 0
      chkhss = j
      if (iopt .eq. 0) return
      if (iopt .gt. 1) j = n
      call pshpr(0)
C     call prm(.false.,1,6,' ',w5,n,5)

      call dsev1(n,hess,dum,w5,0,.false.,.false.,0,
     .j,9d9,i,z,e)
      call poppr
      j = 0
      do  i = 1, n
         if (e(i) .lt. evtol) j = j+1
      enddo   
      if (j .gt. 0 .and. iopt .eq. 1) then
        print 333, j,e(1)
  333   format(/' chkhss: Hessian has',i3,' negative evals; e(1)=',
     .  f8.4,'.  Not updated.')
      elseif (j .gt. 0) then
        print 332, j,e(1)
  332   format(/' chkhss: Hessian has',i3,' negative evals; e(1)=',
     .  f8.4,'.  Project unwanted part.')
C        call prmx('e',e(1),n,n,1)
C        call prmx('z',z,n,n,n)
        call rx('not implemented')
      endif
      chkhss = j
      end
      
      subroutine dsev1(n,h,o,wk,ipr,lx,lov,linv,nmx,emx,nev,z,e)
C- Diagonalize secular equation with overlap 
C----------------------------------------------------------------------
Ci Inputs
Ci    n:    dimension of hamiltonian
Ci    h,o:  hamiltonian, overlap matrices
Ci    wk:   work array length at least (5,11)*ndim (linv = (0,>0))
Ci    ipr:  verbosity
Ci    nmx:  maximum number of eigenvectors to be found
Ci    emx:  eigenvalue limit for eigenvectors to be found
Ci    lov:  true: overlap matrix, false, no overlap matrix
Ci   linv:  1 for inverse iteration
Ci    lx:   true, call x version for overlap handling
Co Outputs
Co    z:    eigenvectors; e, eigenvalues
Co    nev:  number of eigenvectors found
Cu Updates
Cu   07 Apr 07 Bug fix: returns true nev, not nmx in tql2 branch
Cr Remarks
Cr    h,o,z are dimensioned (n,n)
Cr    h,o are OVERWRITTEN in this routine
C----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lov,lx
      integer n,ipr,nmx,nev,linv
      double precision h(n,n),o(n,n),z(n,n),e(n),wk(n,11),emx
C Local variables
      integer ierr,j,iprint

      call tcn('dsev1')
      nev = 0

C --- Eigenvalues of O^-1/2  H  O^-1/2 ---
      if (lov) then
        call dschd(n,n,o,wk,lx,ierr)
        call rxx(ierr.ne.0,'DSEV1: error in dschd')
        if (lx) then
          call dsredx(n,n,h,o,z)
        else
          call dsred(n,n,h,o)
        endif
      endif

      if (linv .eq. 1 .and. nmx .gt. 0) then
        call dtridx(n,n,h,wk,wk(1,4),wk(1,5),wk(1,2))
      endif

      if (nmx .le. 0) then
        call dtridx(n,n,h,e,wk,wk(1,4),wk(1,2))
        do   j = 1, n
           wk(j,1) = wk(j,1)**2
        enddo   
        call tqlrat(n,e,wk,ierr)
        call rxx(ierr.ne.0,'DSEV1: tqlrat cannot find all evals') 
        goto 100
      else if (linv .eq. 1) then
        call imtqlv(n,wk,wk(1,4),wk(1,5),e,wk(1,11),ierr,wk(1,6))
        call rxx(ierr.ne.0,'DSEV1: imtqlv cannot find all evals') 
C   ... Determine number of eigenvectors to be calculated
        nev = 1
        do   j = 2, n
           if (j .le. nmx .and. e(j-1) .le. emx) nev = j
        enddo   
        call tinvit(n,n,wk(1,1),wk(1,4),wk(1,5),nev,e,wk(1,11),z,
     .  ierr,wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10))
        call rxx(ierr.ne.0,'DSEV1: tinvit cannot find all evecs') 
        call dtribx(n,n,h,wk(1,2),nev,z)
      else
        call tred2(n,n,h,e,wk(1,2),z)
        call tql2(n,n,e,wk(1,2),z,ierr)
        call rxx(ierr.ne.0,'DSEV1: tql2 cannot find all evecs') 
        nev = n
      endif

C --- Get the eigenvectors of H - E O ---
      if (.not. lov .or. nmx .le. 0) goto 100
      if (lx) then
        call dcopy(n*n,z,1,h,1)
        call dmpy(o,n,1,h,n,1,z,n,1,n,nev,n)
      else
        call dsbak(n,n,o,nev,z)
      endif

C --- Exit ---
  100 continue
      if(iprint() .ge. 60) print 600, e
  600 format(' evl='/(1x,8f10.5))
      call tcx('dsev1')
      end



      subroutine dtribx(nm,n,ar,tau,m,zr)
c
C     implicit none
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),tau(n,2),zr(nm,m)
      double precision h,s
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  htridi  in their
c          full lower triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c     Aug 1990 MvS altered into daxpy-style loops.  Use this
c     version with htridx, for unit strides
c     ------------------------------------------------------------------
c
      if (m .eq. 0) return
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the full symmetric
c                tridiagonal matrix. ..........
      do  k = 1, n
        do j = 1, m
           zr(k,j) = zr(k,j) * tau(k,1)
        enddo
      enddo

      if (n .eq. 1) return
      do 140 i = 2, n
        l = i - 1
        h = -tau(i,2)
        if (h .eq. 0.0d0) go to 140
        do 130 j = 1, m
          s = 0.0d0
          do 110 k = 1, l
             s = s + ar(k,i) * zr(k,j)
 110      continue
c     .......... double divisions avoid possible underflow ..........
          s = (s / h) / h
          do k = 1, l
             zr(k,j) = zr(k,j) - s * ar(k,i)
          enddo
  130   continue
  140 continue
      end

      
      subroutine dtridx(nm,n,ar,d,e,e2,tau)
c
C     implicit none
C     implicit real*8 (v)
      integer i,j,k,l,n,nm
      double precision ar(nm,n),d(n),e(n),e2(n),tau(n,2)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version adapted from august 1983 htridi.  Differences
C     with htridi are that indices of tau and a are permuted (uses
C     upper triangle of a)
c     ------------------------------------------------------------------
c
      call tcn('dtridx')
      tau(n,1) = 1.0d0
c
      do  i = 1, n
         d(i) = ar(i,i)
      enddo   
      do  300  i = n, 1, -1
        l = i - 1
        h = 0.0d0
        scale = 0.0d0
        if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
        do   k = 1, l
           scale = scale + dabs(ar(k,i))
        enddo   
        if (scale .ne. 0.0d0) go to 140
        tau(l,1) = 1.0d0
  130   e(i) = 0.0d0
        e2(i) = 0.0d0
        go to 290
c
  140   do  150  k = 1, l
          ar(k,i) = ar(k,i) / scale
          h = h + ar(k,i)**2
  150   continue
c
        e2(i) = scale * scale * h
        g = dsqrt(h)
        e(i) = scale * g
        f = dabs(ar(l,i))
c     .......... form next diagonal element of matrix t ..........
        if (f .eq. 0.0d0) go to 160
        tau(l,1) = ( - ar(l,i) * tau(i,1)) / f
        h = h + f * g
        g = 1.0d0 + g / f
        ar(l,i) = g * ar(l,i)
        if (l .eq. 1) go to 270
        go to 170
  160   tau(l,1) = -tau(i,1)
        ar(l,i) = g
  170   f = 0.0d0
        do  240  j = 1, l
c     .......... form element of a*u ..........
          g = 0.d0
          do  180  k = 1, j
            g = g + ar(k,j) * ar(k,i)
  180     continue
          do  200  k = j+1, l
            g = g + ar(j,k) * ar(k,i)
  200     continue
c     .......... form element of p ..........
  220     e(j) = g / h
          f = f + e(j) * ar(j,i)
  240   continue
c
        hh = f / (h + h)
c     .......... form reduced a ..........
        do  j = 1, l
          f = ar(j,i)
          g = e(j) - hh * f
          e(j) = g
          do   k = 1, j
             ar(k,j) = ar(k,j) - f * e(k) - g * ar(k,i)
          enddo
        enddo
  270   do  280  k = 1, l
          ar(k,i) = scale * ar(k,i)
  280   continue
  290   hh = d(i)
        d(i) = ar(i,i)
        tau(i,2) = -scale * dsqrt(h)
  300 continue
      call tcx('dtridx')
      end

      end module m_gradzr
