integer function broyj(n,xin,gin,ir,isw,ipr,beta,dxmx,xtol,gtol, wc,wk,ndw,xnew)
  !- One Broyden step in finding gin = f[xin]-xin = 0
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n:     number of variables
  !i   ir:    Number of iterations of x and g.
  !i          1 initiates a new sequence of mixing;
  !i          broyj uses linear mixing for this iteration.
  !i   isw    1s digit (not implemented)
  !i          0  find minimum
  !i          1  find maximum
  !i         10s digit not used
  !i        100s digit not used
  !i       1000s digit governs convergence criterion:
  !i          1 return when |grad| < gtol
  !i          2 return when max dx < xtol
  !i          3 return when either (1) or (2) satisfied
  !i          4 return when both (1) and (2) are satisfied
  !i   beta:  linear mixing parameter (ir=1 only)
  !i   xin:   input vector, this iteration
  !i   gin:   output-input vector, f[xin]-xin, this iteration
  !i   wc:    weighting for this iteration
  ! o  wk     workspace of 2*ndw*(ir+2), ndw>=n
  ! o         wk must be preserved between calls to broyj.
  ! o         (*,1,0) x of the prior iteration.
  ! o         (*,2,0) g of the prior iteration.
  ! o         (*,1..2,1..ir-1) u and vt of this and prior iterations
  ! o         (*,1,ir) g(this iter) - g (prior iter).
  !o Outputs
  !o   xnew   estimate of x
  !o   broyj
  !r Remarks
  !r   Adapted from Duane Johnson
  ! ----------------------------------------------------------------------
  implicit none
  integer :: isw,ir,n,ipr,ndw
  double precision :: beta,dxmx,wc,xin(n),gin(n),xnew(n),xtol,gtol, wk(ndw,2,0:ir)
  ! Local variables
  integer :: i,ip,j,k,irm1,irm2,lm,ln,nn,i1mach,isw1,isw2,isw3 !dinv
  integer :: ierr
  double precision :: aij,cmj,dfnorm,fac1,fac2,gmi,one,zero,ddot,w0
  parameter (zero=0d0,one=1d0,nn=20)
  double precision :: a(nn,nn),cm(nn),w(nn),d(nn,nn)
  double precision :: betx,diff,gmax,xmax
  !     double precision wl(nn,3),u(nn,nn),v(nn,nn)
  save w,cm,a,w0

  isw1 = mod(isw/10,10)
  isw2 = mod(isw/100,10)
  isw3 = mod(isw/1000,10)
  if (ir > nn) call rxi('broyj: increase nn, need',ir)

  ! --- First iteration: simple mixing ---
  if (ir == 1) then
     w0 = wc
     betx = beta
     gmax = 0
     do  k = 1, n
        gmax = max(gmax,abs(gin(k)))
     enddo

     !$$$c#if AWRITE
     !$$$        if (ipr .ge. 30) then
     !$$$          j = isw3
     !$$$          call awrit8(' broyj:  start'//
     !$$$     .    '%?#(n==2|n==3|n==4)#  xtol=%1;2g#%j#'//
     !$$$     .    '%?#(n==1|n==3|n==4)#  gtol=%1;2g#%j#  beta=%1;2g'//
     !$$$     .    '  w0=%1;2g  isw=%i  gmax=%1;2g',' ',80,i1mach(2),
     !$$$     .    j,xtol,j,gtol,beta,w0,isw,gmax)
     !$$$        endif
     !$$$c#endif

     If (dxmx > 0d0 .AND. gmax > dxmx) then
        betx = beta*dxmx/gmax
        !$$$c#if AWRITE
        !$$$          call awrit3(' broyj:  max shift = %1;3g'//
        !$$$     .    ' is larger than dxmx = %1;3g.  Scale by %1;3g',
        !$$$     .    ' ',80,i1mach(2),gmax,dxmx,dxmx/gmax)
        !$$$c#endif
     endif
     do    k = 1, n
        xnew(k) = xin(k) + betx*gin(k)
     enddo

     ! --- Subsequent iterations: Broyden mixing ---
  else

     !   ... Make xold, gold
     do  20  k = 1, n
        wk(k,1,0) = xin(k) - wk(k,1,0)
        wk(k,1,ir) = gin(k) - wk(k,2,0)
20   enddo

     !   --- Coefficient matrices and the sum for corrections ---
     !   ... dfnorm = |g(i)-g(i-1)|, used for normalization
     dfnorm = dsqrt(ddot(n,wk(1,1,ir),1,wk(1,1,ir),1))
     fac2 = one/dfnorm
     fac1 = beta*fac2
     !   ... Shuffle each prior u,vt to prior+1 iteration
     irm1 = ir-1
     irm2 = ir-2
     do  30  j = irm2, 1, -1
        call dcopy(n,wk(1,1,j),1,wk(1,1,j+1),1)
        call dcopy(n,wk(1,2,j),1,wk(1,2,j+1),1)
30   enddo
     !   ... Make u,vt for this iteration
     do    k = 1, n
        wk(k,1,1) = fac1*wk(k,1,ir) + fac2*wk(k,1,0)
        wk(k,2,1) = fac2*wk(k,1,ir)
     enddo

     !   --- Make  a and b = ( w0**2 I + a )^-1 (symmetric) ---
     do  42  j = 1, irm2
        aij = zero
        cmj = zero
        do    k = 1, n
           cmj = cmj + wk(k,2,ir-j)*gin(k)
           aij = aij + wk(k,2,ir-j)*wk(k,2,1)
        enddo
        a(irm1,j) = aij
        a(j,irm1) = aij
        cm(j) = cmj
42   enddo
     aij = zero
     cmj = zero
     do  k = 1, n
        cmj = cmj + wk(k,2,1)*gin(k)
        aij = aij + wk(k,2,1)*wk(k,2,1)
     enddo
     a(irm1,irm1) = aij
     cm(irm1) = cmj
     w(irm1) = wc

     !   ... Set up and calculate beta matrix
     do   lm = 1, irm1
        do  ln = 1, irm1
           d(ln,lm) = a(ln,lm)*w(ln)*w(lm)
        enddo
        d(lm,lm) = w0**2 + a(lm,lm)*w(lm)*w(lm)
     enddo

     !   --- Invert to make d ---
     !        if (dinv(' ',irm1,nn,d) .ne. 0) then
     !          call rx('broyj: matrix singular')
     !       endif
     call matinv2(irm1,d(1:irm1,1:irm1),ierr)

     !   ... Invert with singular value decomposition
     !        call svd(nn,irm1,irm1,d,wl(1,2),.true.,u,.true.,v,ierr,wl)
     !        call dpzero(d,nn**2)
     !        do  60  ln = 1, irm1
     !   60   d(ln,ln) = 1
     !        call svbksb(nn,irm1,irm1,irm1,wl(1,2),u,v,d,d,wl)
     !    ... This one sometimes hangs up
     !        call rs(nn,irm1,d,wl(1,3),1,v,wl,wl(1,2),ierr)
     !        call dpzero(d,nn**2)
     !        do  60  ln = 1, irm1
     !          print *, 'evl',ln, wl(ln,3)
     !   60   d(ln,ln) = 1
     !        call svbksb(nn,irm1,irm1,irm1,wl(1,3),v,v,d,d,wl)

     !   --- xnew <- vector for the new iteration ---
     do   k = 1, n
        xnew(k) = xin(k) + beta*gin(k)
     enddo
     do  70  i = 1, irm1
        gmi = zero
        do  ip = 1, irm1
           gmi = gmi + cm(ip)*d(ip,i)*w(ip)
        enddo
        do  k = 1, n
           xnew(k) = xnew(k) - gmi*wk(k,1,ir-i)*w(i)
        enddo
70   enddo

     !   ... Cap to maximum allowed shift xnew-xin
     if (dxmx > 0d0) then
        diff = 0
        do  k = 1, n
           diff = max(diff,abs(xnew(k)-xin(k)))
        enddo
        if (diff > dxmx) then
           betx = dxmx/diff
           !$$$c#if AWRITE
           !$$$            call awrit3(' broyj:  max shift = %1;3g'//
           !$$$     .      ' is larger than dxmx = %1;3g.  Scale by %1;3g',
           !$$$     .      ' ',80,i1mach(2),diff,dxmx,dxmx/diff)
           !$$$c#endif
           do   k = 1, n
              xnew(k) = xin(k) + betx*(xnew(k)-xin(k))
           enddo
        endif
     endif

  endif

  ! --- Cleanup, setup for next call ---
  xmax = 0
  gmax = 0
  diff = 0
  do  110  k = 1, n
     xmax = max(xmax,abs(xnew(k)-xin(k)))
     gmax = max(gmax,dabs(gin(k)))
     diff = diff + (xnew(k)-xin(k))**2
     wk(k,2,0) = gin(k)
     wk(k,1,0) = xin(k)
110 enddo
  diff = dsqrt(diff/n)

  j = ir+1
  if (isw3 /= 0 .AND. (gmax == 0 .OR. &
       gmax < gtol .AND. xmax < xtol .AND. isw3 == 4 .OR. &
       gmax < gtol .AND. (isw3 == 1 .OR. isw3 == 3)  .OR. &
       xmax < xtol .AND. (isw3 == 2 .OR. isw3 == 3)  .OR. &
       gmax < gtol .AND. (isw3 == 1 .OR. isw3 == 3))) j = 0
  !$$$c#if AWRITE
  !$$$      if (ipr .ge. 20 .and. j .eq. 0) then
  !$$$        call awrit3(' broyj: converged to max dx'//
  !$$$     .  '=%1;2g, gmax=%1;2g using %i iterations',' ',80,
  !$$$     .  i1mach(2),xmax,gmax,ir)
  !$$$      elseif (ipr .ge. 30 .and. ir .ne. 1) then
  !$$$        call awrit4(' broyj:  ir=%i  dxmax=%1;2g  gmax=%1;2g'//
  !$$$     .  '  wc=%1;2g',' ',80,i1mach(2),ir,xmax,gmax,wc)
  !$$$      endif
  !$$$c#endif
  broyj = j
end function broyj

