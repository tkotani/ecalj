module m_mtchae !- Matches augmentation function to envelope function
  public mtchre
contains
  subroutine mtchae(mode,rsm,eh,l,r,phi,dphi,alfa,beta)!- Matches augmentation function to envelope function
    use m_lgunit,only:stdo !    use m_hansr,only:hansmd
    !i Inputs
    !i   mode  :0 match phi,dphi to h at r
    !i         :1 match phi to h,hdot at r
    !i         :2 return log deriv of h in alfa, K.E in beta
    !i            phi is not used
    !i         :3 match phi,dphi to h,hdot  return K.E in beta
    !i   rsm   :smoothing radius of basis function
    !i   eh    :energy of basis function
    !i   l     :l quantum number
    !i   r     :radius at which to perform matching
    !i   phi   :wave function at r
    !i   dphi  :logarithmic derivative of phi at r
    !l Local variables
    !l   hs    :sm hankel at r
    !l   dhs   :radial derivative of hs
    !l   hsp   :energy derivative of hs
    !l   dhsp  :radial derivative of hsp
    !o Outputs
    !o         (mode 0)
    !o   alfa  :alfa*phi + beta*phidot matches differentiably onto hs
    !o   beta  :alfa*phi + beta*phidot matches differentiably onto hs
    !o         (mode 1)
    !o   alfa  :alfa*hs + beta*hsdot matches differentiably onto phi
    !o   beta  :alfa*hs + beta*hsdot matches differentiably onto phi
    !o         (mode 2)
    !o   alfa  :log deriv of hs (phi not used)
    !o   beta  :K.E. of hs
    implicit none
    integer :: mode,l
    double precision :: dphi,dphip,eh,phi,phip,r,rsm,alfa,beta
    double precision :: det,a,b,val,slo
    double precision :: hs(0:l),dhs(0:l),ddhs(0:l)
    double precision :: hsp(0:l),dhsp(0:l),ddhsp(0:l)
    logical:: isanrg,l_dummy_isanrg
    ! --- Radial part of smooth hankels and its derivatives ---
    call hansmd(12,r,eh,rsm,l,hs,dhs,ddhs,hsp,dhsp,ddhsp)
    ! --- Match hs,dhs to a linear combination of phi,phidot ---
    !     Use  phi=phi(R); phip=phidot(R) dphi=phi'(R); dphip=phidot'(R)
    !     (phi  phip ) (alfa)   (hs )    (alfa)    1  (dphip -phip) (hs )
    !     (          ) (    ) = (   ) -> (    ) = --- (           ) (   )
    !     (dphi dphip) (beta)   (dhs)    (beta)   det (-dphi  phi ) (dhs)
    if (mode == 2) then! --- Match phi to hs; return dhs/hs in alfa, K.E. in beta ---
       alfa = dhs(l)/hs(l)
       beta = -ddhs(l)/hs(l)
       ! --- Match phi,dphi to a linear combination of hs,hsdot; return log derivative in alfa, K.E. in beta
    elseif (mode == 3) then
       det  = hs(l)*dhsp(l) - dhs(l)*hsp(l)
       a    = (phi*dhsp(l) - dphi*hsp(l))/det
       b    = (dphi*hs(l) - phi*dhs(l))/det
       val  = a*hs(l) + b*hsp(l)
       slo  = a*dhs(l) + b*dhsp(l)
       alfa = slo/val
       slo  = -(a*ddhs(l) + b*ddhsp(l))
       beta = slo/val
    endif
  end subroutine mtchae
  subroutine mtchre(l,rsmin,rsmax,emin,emax,r1,r2,phi1,dphi1, & !- Finds envelope function parameters that match conditions on sphere
       phi2,dphi2,rsm,eh,ekin,ir)
    use m_lgunit,only:stdo
    use m_ftox
    !i Inputs
    !i   mode  :controls how matching is done
    !i         :10s digit :
    !i         :0  vary rsm to match h(rsm,eh) to phi1,dphi1 at r1
    !i         :1  Vary eh  to match h(rsm,eh) to phi1,dphi1 at r1
    !i         :2  Vary rsm to match K.E. at r1 (no matching of slope)
    !i         :   hankel energy is input.
    !i         :3  Vary rsm and eh to match both slope and K.E. at r1
    !i         :4  Vary rsm to match K.E. and slope at r1, using
    !i         :   linear combination of h and hdot to match slope.
    !i         :   hankel energy is input.
    !i         :10 Try mode 0.  If unsuccessful, switch to mode 1
    !i         :11 Try mode 1.  If unsuccessful, switch to mode 0
    !i         :100s digit
    !i         :0 always return
    !i         :1 abort if match in log. derivative unsuccessful
    !i         :2 abort if match in K.E. unsuccessful
    !i   l     :l quantum number
    !i   rsmin :when rsm allowed to vary, lower bound for rsm
    !i   rsmax :when rsm allowed to vary, upper bound for rsm
    !i   emin  :when eh allowed to vary, lower bound for eh
    !i   emax  :when eh allowed to vary, upper bound for eh
    !i   r1    :first radius at which to perform matching
    !i   r2    :second radius at which to perform matching (not used now)
    !i   phi1  :wave function at r1
    !i   dphi  :derivative of phi1 at r1
    !i   phi2  :(mode=2,3) K.E. at r1
    !i         :(not used now) wave function at r2
    !i   dphi2 :derivative of phi2 at r2 (not used now)
    ! o Inputs/Outputs
    ! o  rsm   :envelope smoothing radius for smoothed hankel
    ! o        :mode 0: rsm is output
    ! o        :mode 1: rsm is input
    ! o  eh    :envelope energy
    ! o        :mode 0: eh is input
    ! o        :mode 1: eh is output
    !o Outputs
    !o  ir     :information description
    !o         :0  successful match first mode choice
    !o         :1  successful match second mode choice
    !o         :-1 failed to match log derivative
    !o         :-2 failed to match K.E.
    !o         :-3 maximum iterations exceeded
    !o  ekin   :(modes 2,3) kinetic energy at MT boundary
    !l Local variables
    !l   dxrsm :rsm step length used for two-variable search
    !r Remarks
    !r   Properties of h(rsm,eh):
    !r   dh/h is a monotonically increasing function of eh for fixed rsm;
    !r   it is also a monotonically increasing function of rsm for fixed eh.
    !r   However, the kinetic energy is not.
    !r
    !r   Procedure for mode=3 (Vary rsm and eh to match dh/h and K.E. at r1)
    !r   1. Setup.  Find some (rsm,eh) pair that satisfies constraint
    !r         dphi1/phi1=dh/h.  K.E. is ignored in this step.
    !r      a. Begin with rsm=rsmin.
    !r      b. Determine eh that matches dh/h.
    !r         If matching fails, increment rsm by dxrsm/2 and repeat step
    !r         until match found or rsm exceeds rsmax.
    !r      c. mtchre aborts or exits unless a match is found
    !r   2. Iteratively search for K.E. match
    !r      a. Begin with rsm from step 1.  step 1 guarantees that there
    !r         is at least one (rsm,eh) pair that satisfies slope constraint
    !r      b. Iteratively search for rsm that matches K.E.
    !r         Iteration proceeds by regula falsi, which accepts as input
    !r         some (rsm,K.E) pair, returning a new rsm for next iteration.
    !r         In each step, rsm is given.
    !r         1  determine eh by matching slope constraint.
    !r         2  if slope constraint fails, return with (rsm,eh) pair that
    !r            best satisfies constraint (see 4)
    !r         3  compute K.E. for given (rsm,eh)
    !r         4  Follow the (rsm,eh) pair that best satisfies constraint.
    !r            (prepares for possible failure in a future iteration)
    !r         5  call rfalsi for a new estimate for rsm, or until
    !r            iteration converges.
    !r
    !     implicit none
    ! ... Passed parameters
    integer :: mode,l,ir
    double precision :: dphi1,dphi2,eh,phi1,phi2,r1,r2,rsm,rsmin,rsmax,    emin,emax,ekin
    ! ... Local parameters
    integer :: iter,IPRT1,IPRTW,ipr,maxit,mode0,mode2,ir1,         ipass
    double precision :: tol,eh0,dxrsm
    parameter (tol=1d-12,IPRT1=100/1,IPRTW=50/1,maxit=50,dxrsm=.05d0)
    double precision :: xnow,alfa,beta,xclose,bclose,eclose,tclose,wk(12)
    logical:: isanrg, l_dummy_isanrg
    character(8):: xt
    call getpr(ipr)
    mode=103
    mode0 = mod(mode,100)
    mode2 = mod(mode/100,10)
    l_dummy_isanrg=isanrg(l,0,8,'mtchre','l',.true.)
    ! ... Start of matching in current mode
10  continue
    ! --- Vary rsm to match phi to h ---
    if (mode0 == 0 .OR. mode0 == 10) then
       call mtchr2(0,l,rsmin,rsmax,(rsmin+rsmax)/2,r1,phi1,dphi1,rsm,  eh,ekin,ir)
       !   ... rfalsi unable to bound the search : handle error
       if (ir < 0) then
          if (ipr >= IPRTW) write (stdo,'('' mtchre (warning) matching failed varying rsm'')')
          if (mode0 == 10) then
             mode0 = 1
             goto 10
          endif
       endif
       ! --- Vary eh to match phi to h ---
    elseif (mode0 == 1 .OR. mode0 == 11) then
       call mtchr2(1,l,emin,emax,(emin+emax)/2,r1,phi1,dphi1,rsm,eh, ekin,ir)
       !   ... rfalsi unable to bound the search : handle error
       if (ir < 0) then
          if (ipr >= IPRTW) write (stdo, '('' mtchre (warning) matching failed varying eh'')')
          if (mode0 == 10) then
             mode0 = 0
             goto 10
          endif
       endif
       ! --- Vary rsm (and eh) to match K.E. (and log der) to h ---
       !     See Remarks for description of procedure for mode3
    elseif (mode0 >= 2 .AND. mode0 <= 4) then
       xnow = min(rsmin,rsmax)
       if (xnow < 0) call rx1('mtchr2: bad range, rsm = ',xnow)
       xclose = 0
       bclose = 9d9
       !       Not used now, but keep anyway
       ipass = 1
       iter = 0
       !       Initial value for rfalsi
       ir = 0
       !       call pshpr(min(ipr,1))
       call pshpr(ipr-10)
       if (ipr >= IPRT1) write(stdo,261)
261    format(' l  it  ir',6x,'Rsm',9x,'Eh',7x,'slope/val',6x,'K.E.',5x, 'target K.E.')
       !   ... mode 3-specific setup: find some (rsm,eh) satisfying slope cond.
       !       Loop through rsm(=xnow) in uniform increments;
       !       Find first rsm which can match dh/h to dphi1/phi1
       !       K.E. is ignored in this first step
       if (mode0 == 3) then
110       continue
          eh0 = (emin+emax)/2
          call mtchr2(1,l,emin,emax,eh0,r1,phi1,dphi1,xnow,eh,ekin,ir1)
          !         No match found for this rsm ... increment rsm
          if (ir1 < 0) then
             !           No match for any rsm ... give up
             if (xnow < min(rsmin,rsmax)-1d-6 .OR. &
                  xnow > max(rsmin,rsmax)+1d-6 .OR. &
                  rsmin == rsmax) then
                ir = -1
                if(mode2>=1)call rx('mtchre:err to match phi to envelope. l dphi1/phi1='//trim(xt(l))//' '//trim(ftof(dphi1/phi1)))
                return
             endif
             xnow = xnow + dxrsm
             goto 110
          endif
       endif
       !       End of mode 3-specific setup
       !   ... Iteratively try to match K.E. with mode-specific constraints
120    continue       !       mode 2-specific : No matching of slope
       if (mode0 == 2) then
       elseif (mode0 == 3) then !eh determined from dphi1/phi1=dh/h
          eh0 = eh
          call mtchr2(1,l,emin,emax,eh0,r1,phi1,dphi1,xnow,eh,ekin,ir1)
          !         Failed to match log derivative.  Resort to best prior case
          if (ir1 < 0) then
             if (iter == 0) then
                ir = -1
                if (mode2 >= 1) call rxi('mtchre: '// &
                     'failed to match phi to envelope, mode',mod(mode,100))
                return
             endif
             !           Flags that energy is near boundary point
             !           ipass = 2
             ir = -1
             if (iter /= 0) ir = -2
             goto 80
          endif
       endif
       !   ... K.E. for xnow = current guess for rsm.
       if (mode0 /= 4) then          !         phi1 only dummy here
          call mtchae(2,xnow,eh,l,r1,phi1,phi1,alfa,ekin)
       else
          call mtchae(3,xnow,eh,l,r1,phi1,dphi1,alfa,ekin)
       endif
       beta = ekin - phi2
       rsm = xnow
       !       Keep running track of closest approach in case no match found
       if (abs(beta) < bclose) then
          xclose = xnow
          bclose = abs(beta)
          eclose = eh
          tclose = ekin
       endif
       call pshpr(0)
       call rfalsi(xnow,beta,tol,tol,tol/10,dxrsm,34,wk,ir)
       call poppr
       if (ir > 1) call rxi('bug in mtchre, ir=',ir)
       iter = iter+1
       if (ipr >= IPRT1)  write (stdo,60) l,iter,ir,wk(1),eh,alfa,beta+phi2,phi2
60     format(i2,i4,i4,5f12.6)
       !   ... rfalsi unable to find K.E... Use closest point
       if (xnow < min(rsmin,rsmax)-1d-6 .OR. &
            xnow > max(rsmin,rsmax)+1d-6) then
          ir = -2
          goto 80
          !   ... rfalsi either has converged or requires another iteration
       else
          if (iter < maxit .AND. ir < 0) goto 120
          if (iter >= maxit) ir = -3
          if (iter >= maxit) goto 80
          rsm = xnow
       endif
    else
       call rxi('mtchre: bad mode,',mode)
    endif
    ! ... Cleanup and exit
    if (ir < 0 .AND. mode2 == 1) call rxi('mtchre: failed to match phi to envelope, mode',mod(mode,100))
    ir = 0
    if (mode0 /= mod(mode,100)) ir = 1
    call poppr
    return
    ! --- Handle mode2 when matching failed ---     !     ir should be set before jumping here
80  continue
    xnow = xclose
    eh0  = eclose
    ekin = tclose
    ! ... Closest point at boundary point in rsm or eh; no further search
    if (xclose == min(rsmin,rsmax) .OR. xclose == max(rsmin,rsmax) .OR. ipass == 2 .OR. .TRUE. ) then
       !       get K.E. at closest point.  phi1 only dummy here
       rsm  = xclose
       eh   = eclose
       ekin = tclose
       call mtchae(2,rsm,eh,l,r1,phi1,phi1,alfa,beta)
       ! ... Not implemented: find K.E. closest to target
    else
       print *, xclose
       call rx('mtchre : not implemented')
    endif
    call poppr
    if (ir == -1) call rxi('mtchre: failed to match phi to envelope, mode',mod(mode,100))
  end subroutine mtchre
  subroutine mtchr2(mode,l,x1,x2,x0,r,phi,dphi,rsm,eh,ekin,info) ! Match rsm or eh to value and slope of phi at surface of sphere
    use m_lgunit,only:stdo
    !i   mode  :controls how matching is done
    !i         :0 vary rsm to match h(rsm,eh) to phi,dphi at r
    !i         :  where phi=val, dphi=slope
    !i         :1 Vary eh  to match h(rsm,eh) to phi,dphi at r
    !i         :  where phi=val, dphi=slope
    !i         :2 Vary rsm to match h(rsm,eh) to phi,dphi at r
    !i         :  where phi=val, dphi=K.E.
    !i   l     :l quantum number
    !i   x1    :lower bound for rsm or eh
    !i   x2    :upper upper for rsm or eh
    !i   x0    :initial guess for rsm or eh (-99 => not used)
    !i   r     :first radius at which to perform matching
    !i   phi   :wave function at r
    !i   dphi  :derivative of phi at r
    !l Local variables
    ! o Inputs/Outputs
    ! o  rsm   :envelope smoothing radius for smoothed hankel
    ! o        :mode 0: rsm is output
    ! o        :mode 1: rsm is input
    ! o  eh    :envelope energy
    ! o        :mode 0: eh is input
    ! o        :mode 1: eh is output
    ! o  ir    :error
    !o Outputs
    !o  ekin   :Kinetic energy of h(rsm,eh,r)
    !o  info   :information description
    !o         :0 successful match first mode choice
    !o         :1 successful match second mode choice
    !o         :-1 (min,emax) do not bound search
    !o         :-2 root finding unsuccessful after maxit iterations
    implicit none
    integer :: mode,l,info
    double precision :: dphi,eh,phi,r,rsm,x1,x2,x0,ekin
    integer :: ir,iter,IPRT,ipr,maxit
    double precision :: xnow,dxmx,alfa,beta,tol,wk(12)
    parameter (tol=1d-12,IPRT=100/1,maxit=50)
    call getpr(ipr)
    iter = 0
    ir = 0
    info = 0
    if (ipr >= IPRT) write(stdo,261) mode,l,r
261 format(' mtchr2 mode',i2,'  l =',i2,'  r =',f10.6,/' l  it  ir',6x,'Rsm',9x,'Eh',8x,'phi''/phi     target')
    ! ... Vary rsm to match phi to h.  Do iteratively:
    if (mode == 0) then
       xnow = min(x1,x2)
       dxmx = max(x1,x2) - xnow
       if (xnow < 0) call rx1('mtchr2: bad lower rsm : ',xnow)
100    continue
       !       Match phi,dphi to hs,hsdot; find point where hsdot=0.
       !       call mtchae(1,xnow,eh,l,r,phi,dphi,0d0,0d0,alfa,beta)
       !       alfa = alfa/phi
       !       beta = beta/phi
       !       Match dphi/phi to dhs/hs. dhs/hs is monotonic in rsm for rsm<r
       call mtchae(2,xnow,eh,l,r,phi,dphi,alfa,ekin)
       beta = alfa - dphi/phi
       !     mode=2 : match phi to hs; match dphi to K.E. of hs
       !     not ready
       rsm = xnow
       call pshpr(0)
       call rfalsi(xnow,beta,tol,tol,tol/10,dxmx,10,wk,ir)
       call poppr
       if (ir > 1) call rxi('bug in mtchr2, ir=',ir)
       iter = iter+1
       if (ipr >= IPRT)write (stdo,60) l,iter,ir,wk(1),eh,alfa,dphi/phi
60     format(i2,i4,i4,2f12.6,1x,1p,4e12.3)
       !  ... after 2rd iteration try estimate x0, if supplied
       if (iter == 2 .AND. x0 > min(x1,x2) .AND. x0 < x2) xnow = x0
       !   ... rfalsi unable to bound the search : handle error
       if (xnow < min(x1,x2)-1d-6 .OR. xnow > max(x1,x2)+1d-6) then
          info = -1 
       else           !   ... rfalsi either has converged or requires another iteration
          if (iter < maxit .AND. ir < 0) goto 100
          if (iter >= maxit) info = -2
          rsm = xnow
       endif
    elseif (mode == 1) then ! ... Vary eh iteratively to match log deriv of phi to that of h.
       xnow = min(x1,x2)
       dxmx = max(x1,x2) - xnow
       if (xnow >= 0) call rx1('mtchr2: bad upper eh : ',xnow)
200    continue
       !       Try and match phi to hs,hsdot, and find point where hsdot=0.
       !        call mtchae(1,rsm,xnow,l,r,phi,dphi,0d0,0d0,alfa,beta)
       !        alfa = alfa/phi
       !        beta = beta/phi

       !       Try and match dphi/phi to alfa=dhs/hs.
       !       NB: dhs/hs is monotonic in eh for rsm<r
       call mtchae(2,rsm,xnow,l,r,phi,dphi,alfa,ekin)
       beta = alfa - dphi/phi
       eh = xnow
       call pshpr(0)
       call rfalsi(xnow,beta,tol,tol,tol/10,dxmx,10,wk,ir)
       call poppr
       if (ir > 1) call rxi('bug in mtchr2, ir=',ir)
       iter = iter+1
       if (ipr >= IPRT) write(stdo,60)l,iter,ir,rsm,wk(1),alfa,dphi/phi
       !        if (ipr.ge.IPRT) write(stdo,60)l,iter,ir,rsm,wk(1),alfa,dphi/phi,xnow-wk(1),alfa-dphi/phi
       if (iter == 2 .AND. x0 > min(x1,x2) .AND. x0 < x2) xnow = x0 !after 2rd iteration try estimate x0, if supplied
       !   ... rfalsi unable to bound the search : handle error
       if (xnow < min(x1,x2)-1d-6 .OR. xnow > max(x1,x2)+1d-6) then
          info = -1           !   ... rfalsi either has converged or requires another iteration
       else
          if (iter < maxit .AND. ir < 0) goto 200
          if (iter >= maxit) info = -2
          eh = xnow
       endif
       if (ipr >= IPRT) print '('' exit mtchr2 info ='',i2)',info
       ! ... Vary rs to match K.E. of h to dphi.  Do iteratively, matching
       !     phi to hs; find point where K.E. matches.
    elseif (mode == 2) then
       continue
    elseif (mode == 2) then
       continue
    else
       call rxi('mtchr2 : bad mode',mode)
    endif
  end subroutine mtchr2
  subroutine hansmd(mode,r,e,rsm,lmax,hs,dhs,ddhs,hsp,dhsp,ddhsp) !Value and some derivatives of smoothed radial Hankel functions
    use m_hansmr,only: hansmr,hansmronly
    use m_hansr,only:  hansr
    !i Inputs
    !i   r     :radius
    !i   e     :hankel energy
    !i   rsm   :hankel smoothing radius
    !i   lmax  :make function values for l between (0:lmax)
    !o  Outputs:
    !o   hs    : function values of the radial sm. hankel h(e,r,0:lmax)
    !o         : A solid hankel is H=h*Y_L, where Y_L are the spherical harmonics for unit radius (no scaling by r**l)
    !o   dhs   : radial derivative of hs
    !o   ddhs  : radial part of Laplacian of hs, i.e. 1/r d^2 (r h) /dr^2 OR some other second derivative (see mode)
    !o   hsp   : energy derivative of hs
    !o   dhsp  : mixed energy + radial derivative of hs
    !o   ddhsp : 3-d laplacian of hsp YL
    !r Remarks
    !i        ddhs = 1/r d^2 (r*h_l) / dr^2  - l(l+1)/r^2 h_l :  NB: ddhs = laplacian of 3-dimensional hs_l YL
    !r  See J. Math. Phys. 39, 3393 (1998).
    !r    For radial derivative, see JMP 39, 3393, Eq. 4.7   h'  = l/r h_l - h_l+1
    !r    Second radial derivative:                          h'' = l(l-1)/r^2 xi_l - (2l+1)/r xi_l+1 + xi_l+2 1/r d^2/dr^2 (r*h_l)
    !r                                                           = l(l+1)/r^2 h_l - (2l+3)/r h_l+1 + h_l+2
    !r    Energy derivative:  see JMP 39, 3393, Eq. 7.5
    !r      hp_l = r/2 h_l-1  Special case l=0: hp_0 = 1/2 h_-1 ?
    !r    Mixed energy + radial derivative: hp'_l = h(l-1)*l/2 - h(l)*r/2
    !r    Mixed energy + kinetic energy     hp'' = -(2l+3)/2 h_l + h_l+1*r/2
    !r
    !r  Note connection with hansmr, which makes xi(l) = h(l) / r^l
    implicit none
    integer :: mode,lmax,idx,l,mode0,mode1,i
    real(8) :: r,e,rsm, hs(0:lmax),dhs(0:lmax),ddhs(0:lmax),hsp(0:lmax),dhsp(0:lmax),ddhsp(0:lmax), xi(-1:lmax+2)
    real(8):: ra,sre,akap,arsm,earsm,x,xx,erfcee,um,up
    real(8)::     y0 = 1/dsqrt(16*datan(1d0))
    if (lmax < 0) return
    mode0 = 2 !mod(mode,10)
    mode1 = 1 !mod(mode/10,10)
    if(hansmronly) then
        call hansmr(r,e,1d0/rsm,xi(0:lmax+2),lmax+2) 
       ra = r/rsm
       akap = dsqrt(-e)
       sre = akap*r
       arsm = akap*rsm/2
       earsm = dexp(-arsm**2)/2
       xx = earsm*y0*dexp(-(r/rsm)**2)/r
       x = ra-arsm
       if(x>0 ) um=dexp(-sre)/r-xx*erfcee(x) ! ---   Evaluate um,up ---
       if(x<=0) um=xx*erfcee(x)
       up= xx*erfcee(ra + arsm) !assumes x gt 0
       xi(-1) = (um + up)*r/akap
       xi(1:lmax+2)=[(xi(i)*r**i,i=1,lmax+2)]
    else       
       call hansr(rsm,-1,lmax+2,1,[lmax+2],[e],[r**2],1,1,[idx],11,xi)
    endif
    hs=xi(0:lmax)
    dhs(:)  = [(xi(l)*l/r - xi(l+1),l=0,lmax)]
    ddhs = [(                   - (2*l+3)/r*xi(l+1) + xi(l+2), l=0,lmax)]
    hsp   = [(xi(l-1)*r/2,                    l=0,lmax)]
    dhsp  = [((xi(l-1)*l - xi(l)*r)/2,        l=0,lmax)]
    ddhsp = [(- (2*l+3)*xi(l)/2 + xi(l+1)*r/2,l=0,lmax)]
    hsp(0) = xi(-1)/2
  end subroutine hansmd
endmodule m_mtchae
