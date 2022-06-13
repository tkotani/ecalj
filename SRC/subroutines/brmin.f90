subroutine brmin(n,x,g,isw,ipr,dxmx,xtol,gtol,hmax,w,diff,hess,ir)
  use m_gradzr,only:pgradz
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

     call info8(-25,0,0,' brmin: start'// &
          '%?#(n==1|n==2|n==3)#  xtol=|%1;2g|#%j#'// &
          '%?#(n==0|n==2|n==3)#  gtol=|%1;2g|#%j#'// &
          '  isw=%i  |g|=%1;3g', &
          isw1,xtol,isw1,gtol,isw,dsqrt(ddot(n,g,1,g,1)),0,0)

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

  if (scale /= 1) then
     call info8(-30,0,0,' brmin: iter %i  max shift=%1;4g  '// &
          '|g|=%1;4g  max g=%1;4g  scale H by %1;3g', &
          iabs(ir),diff,dsqrt(ddot(n,g,1,g,1)),gmax,scale,0,0,0)
  else
     call info8(-30,0,0,' brmin: iter %i  max shift=%1;4g  '// &
          '|g|=%1;4g  max g=%1;4g%?#n#  evmin=%1;2g#%j#',iabs(ir), &
          diff,dsqrt(ddot(n,g,1,g,1)),gmax,isw31,evmin,0,0)
  endif
  npr = min(n,8)
  call info0(-60,0,0,' Updated x: ')
  do  i = 1, n, npr
     call info(-60,0,0,' %n;10,6D',min(npr,n-i+1),x(i))
  enddo

  !     Global convergence criteria: shift relative to start of line min
  itol = mod(isw1,4)
  call pgradz(itol,n,1d0,w(1,ns),g,xtol,gtol, &
       dxtop,resx,gtop,resg,dum,cnvgx,cnvgg)

  ! ... Case brmin has converged globally
  if (resg == 0 .OR. (cnvgg .AND. cnvgx)) then
     call info5(-20,0,0, &
          '%x brmin converged to dxmax=%1;3g, |dx|=%1;3g,'// &
          '  gmax=%1;3g, |grad|=%1;3g in %i iter', &
          dxtop,resx,gtop,resg,ir)
     npr = min(n,6)
     call info(-30,0,0,' p=%n:;10F',npr,x)
     call info(-30,0,0,' g=%n:;10F',npr,g)
     ir = 0
  endif

  call poppr


end subroutine brmin

