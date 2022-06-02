      subroutine brmin(n,x,g,isw,ipr,dxmx,xtol,gtol,hmax,w,diff,hess,ir)
      use m_gradzr,only:pgradz,chkhss
C- One Broyden step in finding the root of an n-dimensional function
C ----------------------------------------------------------------
Ci Inputs
Ci   n:     number of variables
Ci   x      current for which gradients are specified.
Ci   g      gradient at current position
Ci   isw    1s digit (not implemented)
Ci          0  find minimum
Ci          1  find maximum
Ci         10s digit governs convergence criterion:
Ci          0 convergence when <g> < gtol
Ci          1 convergence when <dx> < xtol
Ci          2 convergence when <g> < gtol and <dx> < xtol
Ci          3 convergence when <g> < gtol or <dx> < xtol
Ci            Here, <g> and <dx> are the max component of gradient
Ci            and change in x between iterations
Ci        100s digit
Ci          1  Set this bit if there is an x-dependent hmax.
Ci       1000s digit: Hessian matrix.  4's bit is independent from 1,2
Ci          0  no tests are made for Hessian
Ci          1  If Hessian is not positive definite, it is not updated.
Ci          2  project out parts of Hessian corresponding to negative
Ci             eval (not implemented).
Ci          4  On the second iteration, Hessian is globally scaled
Ci             once a first estimate is known (D. Novikov)
Ci   ipr    verbosity
Ci   dxmx   maximum allowed change in x for this step
Ci   xtol  tolerance in delta x for global minimization.
Ci         xtol > 0 : tolerance compared to |delta x|
Ci         xtol < 0 : tolerance compared largest component, delta max(|x|)
Ci   gtol  tolerance in gradient for global minimization
Ci         gtol > 0 : tolerance compared to |delta g|
Ci         gtol < 0 : tolerance compared largest component, delta max(|g|)
Ci   hmax   a variable-specific maximum; shift is scaled so that
Ci          no variable exceeds maximum (not tested).  Only used for
Ci          those hmax>0.
Ci   w      work array of dimension at least w(n,6), or if Hessian
Ci          tested for positive definiteness, w(n,9+2*n), and if
Ci          negative evals projected out, w(n,9+2*n).  Pieces signify:
Ci          (*,nx)  prior position x (must be preserved between calls)
Ci          (*,ng)  input gradient g (must be preserved between calls)
Ci          (*,ns)  shift = x(output) - x(input)
Ci          (*,nd)  change in g
Ci          (*,nh)  a help array
Ci          (*,nh2) a help array
Ci          (*,nhs) a local copy of the Hessian (needed only if
Ci                  Hessian is tested for positive definiteness)
Cio  ir:    iteration number.  Should be zero on first iteration.
Cio         It is updated internally, and returns the following:
Cio         0  brmin has converged to specified tolerance
Cio        >0  brmin needs gradients (preferably at the output x)
Cio        <0  like >0, but the Hessian was found to be not positive
Cio            definite, and was not updated.
Cio  hess inverse Hessian.  brmin updates hess as minimization
Cio         proceeds.  Initial Hessian may be specified; it is set to
Cio         unity if initially zero.
Co Outputs
Co   x      suggested new position at which to evaluate gradients
Co   diff   the largest shift
Cr Remarks
Cr   Adapted from D. Novikov.
Cr   This routine generates the inverse of the Hessian matrix through
Cr   the Broyden scheme, and supplies an update of the coordinates.
Cu Updates
Cu   08 Mar 06 Modified convergence checks for consistency w/ gradzr
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,ir,isw,ipr
      double precision diff,dxmx,xtol,gtol,x(n),g(n),hmax(n),w(n,10),
     .hess(n,n)
C ... Local parameters
      integer nx,ng,ns,nd,nh,nh2,nhs,itol
      parameter (nx=1,ng=2,ns=3,nd=4,nh=5,nh2=6,nhs=nd+5)
      double precision dxtop,resx,gtop,resg,dum
      logical cnvgx,cnvgg
      double precision beta,one,evtol
      parameter (beta=1d0, one=1d0, evtol=1d-10)
      double precision factor,fmin,a,b,ddot,dasum,gmax,evmin,scale
C     double precision wk(n)
      integer i,j,isw1,isw32,isw2,isw31,npr

      isw1 = mod(isw/10,10)
      isw2 = mod(isw/100,10)
      isw31 = mod(mod(isw/1000,10),4)
      isw32 = mod(isw/1000,10)/4
      scale = 1
      evmin = 0
      call pshpr(ipr)

C --- First iteration or restart ---
      if (ir .eq. 0) then
C   ... If starting Hessian is zero, set to unity.
        if (dasum(n*n,hess,1) .eq. 0) call dcopy(n,1d0,0,hess,n+1)
        ir = 1
        do  6  i = 1, n
          w(i,nx) = x(i)
          w(i,ng) = g(i)
          do    j = 1, n
             x(i) = x(i) - hess(i,j)*g(j)
          enddo   
    6   continue

        call info8(-25,0,0,' brmin: start'//
     .  '%?#(n==1|n==2|n==3)#  xtol=|%1;2g|#%j#'//
     .  '%?#(n==0|n==2|n==3)#  gtol=|%1;2g|#%j#'//
     .  '  isw=%i  |g|=%1;3g',
     .  isw1,xtol,isw1,gtol,isw,dsqrt(ddot(n,g,1,g,1)),0,0)

C --- Subsequent iterations  ---
      else

        ir = iabs(ir)+1

C   ... Global scaling of the Hessian matrix
        if (ir .eq. 2 .and. isw32 .ne. 0) then
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

C   --- Make dx,dg; save Hessian,x,g; make some work arrays ---
        do  20  i = 1, n
          w(i,ns) = x(i) - w(i,nx)
          w(i,nd) = g(i) - w(i,ng)
   20   continue
C       Local copy of the original Hessian
        if (isw31 .ne. 0) call dcopy(n*n,hess,1,w(1,nhs),1)
        do  30  i = 1, n
          w(i,nx) = x(i)
          w(i,ng) = g(i)
          w(i,nh) = 0d0
          do   j = 1, n
             w(i,nh) = w(i,nh) + hess(i,j)*w(j,nd)
          enddo   
   30   continue
        a = ddot(n,w(1,nd),1,w(1,nh),1)
        b = ddot(n,w(1,nd),1,w(1,ns),1)
        do    i = 1, n
           w(i,nh2) = 1/a*w(i,nh) - 1/b*w(i,ns)
        enddo   

C   --- Update inverse of Hessian ---
        do    i = 1, n
        do   j = i, n
            hess(j,i) = hess(j,i)
     .      - w(j,nh)*w(i,nh)/a + w(j,ns)*w(i,ns)/b
     .            + beta*a*b*w(j,nh2)*w(i,nh2)
         enddo
         enddo
C   ... Make it symmetric
        do    i = 2, n
        do   j = 1, i-1
           hess(j,i) = hess(i,j)
        enddo
        enddo

C   --- Check that Hessian is positive definite ---
        if (isw31 .ne. 0) then
C         Keep a local copy of hessian, since dsev1 destroys it.
          call dcopy(n*n,hess,1,w(1,nhs+n),1)
          j = chkhss(w(1,nhs+n),n,w(1,nd),evtol,isw31,w(1,nhs),w(1,ns))
          evmin = w(1,ns)
          if (j .gt. 0 .and. (isw31 .eq. 1)) then
            call dcopy(n*n,w(1,nhs),1,hess,1)
            ir = -ir
          elseif (j .gt. 0) then
          endif
        endif

C   --- Update x from Hessian, g ---
        do  i = 1, n
        do  j = 1, n
           x(i) = x(i) - hess(i,j)*g(j)
        enddo
        enddo

C   ... End of subsequent Hessian update
      endif

C --- Limit the shift, depending on dxmx and hmax ---
      do  i = 1, n
         w(i,ns) = x(i) - w(i,nx)
      enddo   
      diff = 0
      do  i = 1, n
         diff = max(diff,dabs(w(i,ns)))
      enddo   
      if (diff .eq. 0) goto 100
      if (diff .gt. dxmx) then
        factor = dxmx/diff
        print 334, diff,factor
  334   format(' brmin: max shift =',f9.5,
     .  ' is larger than dxmx.  Scale by',f8.4)
        do  84  i = 1, n
          w(i,ns) = factor * w(i,ns)
          x(i) = w(i,ns) + w(i,nx)
   84   continue
        diff = factor*diff
      endif
C ... Check local shifts
      if (isw2 .ne. 0) then
        fmin = 1
        do  85  i = 1, n
          if (hmax(i) .gt. 0) then
            factor = min(one,abs(hmax(i)/x(i)))
            if (factor .lt. fmin) then
              fmin = factor
              print '('' Variable No.'',i3,'' is out of range'',f8.6)',
     .        i,x(i)
            endif
          endif
   85   continue
        if (fmin .lt. one) then
          print '('' Rescale all shifts by'',f8.6)', fmin
          do  86  i = 1, n
            x(i) = fmin*x(i)
            w(i,ns) = x(i) - w(i,nx)
   86     continue
        endif
      endif

C --- Cleanup ---
  100 continue
      diff = 0
      gmax = 0
      do  i = 1, n
        gmax = max(gmax,dabs(g(i)))
        diff = max(diff,dabs(w(i,ns)))
      enddo  

      if (scale .ne. 1) then
        call info8(-30,0,0,' brmin: iter %i  max shift=%1;4g  '//
     .  '|g|=%1;4g  max g=%1;4g  scale H by %1;3g',
     .  iabs(ir),diff,dsqrt(ddot(n,g,1,g,1)),gmax,scale,0,0,0)
      else
        call info8(-30,0,0,' brmin: iter %i  max shift=%1;4g  '//
     .  '|g|=%1;4g  max g=%1;4g%?#n#  evmin=%1;2g#%j#',iabs(ir),
     .  diff,dsqrt(ddot(n,g,1,g,1)),gmax,isw31,evmin,0,0)
      endif
      npr = min(n,8)
      call info0(-60,0,0,' Updated x: ')
      do  i = 1, n, npr
        call info(-60,0,0,' %n;10,6D',min(npr,n-i+1),x(i))
      enddo

C     Global convergence criteria: shift relative to start of line min
      itol = mod(isw1,4)
      call pgradz(itol,n,1d0,w(1,ns),g,xtol,gtol,
     .dxtop,resx,gtop,resg,dum,cnvgx,cnvgg)

C ... Case brmin has converged globally
      if (resg .eq. 0 .or. (cnvgg .and. cnvgx)) then
        call info5(-20,0,0,
     .  '%x brmin converged to dxmax=%1;3g, |dx|=%1;3g,'//
     .  '  gmax=%1;3g, |grad|=%1;3g in %i iter',
     .  dxtop,resx,gtop,resg,ir)
        npr = min(n,6)
        call info(-30,0,0,' p=%n:;10F',npr,x)
        call info(-30,0,0,' g=%n:;10F',npr,g)
        ir = 0
      endif

      call poppr


      end

