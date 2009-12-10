      subroutine mtchae(mode,rsm,eh,l,r,phi,dphi,phip,dphip,alfa,beta)
C- Matches augmentation function to envelope function
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 match phi,dphi to h at r
Ci         :1 match phi to h,hdot at r
Ci         :2 return log deriv of h in alfa, K.E in beta
Ci            phi is not used
Ci         :3 match phi,dphi to h,hdot  return K.E in beta
Ci   rsm   :smoothing radius of basis function
Ci   eh    :energy of basis function
Ci   l     :l quantum number
Ci   r     :radius at which to perform matching
Ci   phi   :wave function at r
Ci   dphi  :logarithmic derivative of phi at r
Ci   phip  :radial derivative of phi  (used when mode=0)
Ci   dphip :radial derivative of dphi (used when mode=0)
Cl Local variables
Cl   hs    :sm hankel at r
Cl   dhs   :radial derivative of hs
Cl   hsp   :energy derivative of hs
Cl   dhsp  :radial derivative of hsp
Co Outputs
Co         (mode 0)
Co   alfa  :alfa*phi + beta*phidot matches differentiably onto hs
Co   beta  :alfa*phi + beta*phidot matches differentiably onto hs
Co         (mode 1)
Co   alfa  :alfa*hs + beta*hsdot matches differentiably onto phi
Co   beta  :alfa*hs + beta*hsdot matches differentiably onto phi
Co         (mode 2)
Co   alfa  :log deriv of hs (phi not used)
Co   beta  :K.E. of hs
Cr Remarks
Cu Updates
Cu   16 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,l
      double precision dphi,dphip,eh,phi,phip,r,rsm,alfa,beta
C ... Local parameters
      double precision det,a,b,val,slo
      double precision hs(0:l),dhs(0:l),ddhs(0:l)
      double precision hsp(0:l),dhsp(0:l),ddhsp(0:l)

C --- Radial part of smooth hankels and its derivatives ---
      call hansmd(12,r,eh,rsm,l,hs,dhs,ddhs,hsp,dhsp,ddhsp)

C --- Match hs,dhs to a linear combination of phi,phidot ---
C     Use  phi=phi(R); phip=phidot(R) dphi=phi'(R); dphip=phidot'(R)
C     (phi  phip ) (alfa)   (hs )    (alfa)    1  (dphip -phip) (hs )
C     (          ) (    ) = (   ) -> (    ) = --- (           ) (   )
C     (dphi dphip) (beta)   (dhs)    (beta)   det (-dphi  phi ) (dhs)
      if (mode .eq. 0) then
        det  = phi*dphip - dphi*phip
        alfa = (hs(l)*dphip - dhs(l)*phip)/det
        beta = (dhs(l)*phi - hs(l)*dphi)/det

C --- Match phi,dphi to a linear combination of hs,hsdot ---
C     (hs  hsp ) (alfa)   (phi )    (alfa)    1  (dhsp -hsp) (phi )
C     (        ) (    ) = (    ) -> (    ) = --- (         ) (    )
C     (dhs dhsp) (beta)   (dphi)    (beta)   det (-dhs  hs ) (dphi)
      elseif (mode .eq. 1) then
        det  = hs(l)*dhsp(l) - dhs(l)*hsp(l)
        alfa = (phi*dhsp(l) - dphi*hsp(l))/det
        beta = (dphi*hs(l) - phi*dhs(l))/det

C --- Match phi to hs; return dhs/hs in alfa, K.E. in beta ---
      elseif (mode .eq. 2) then
        alfa = dhs(l)/hs(l)
        beta = -ddhs(l)/hs(l)

C --- Match phi,dphi to a linear combination of hs,hsdot; ---
C     return log derivative in alfa, K.E. in beta
      elseif (mode .eq. 3) then
        det  = hs(l)*dhsp(l) - dhs(l)*hsp(l)
        a    = (phi*dhsp(l) - dphi*hsp(l))/det
        b    = (dphi*hs(l) - phi*dhs(l))/det
        val  = a*hs(l) + b*hsp(l)
        slo  = a*dhs(l) + b*dhsp(l)
        alfa = slo/val
        slo  = -(a*ddhs(l) + b*ddhsp(l))
        beta = slo/val

      else
        call rxi('mtchae: bad mode',mode)
      endif

      end

      subroutine mtchre(mode,l,rsmin,rsmax,emin,emax,r1,r2,phi1,dphi1,
     .  phi2,dphi2,rsm,eh,ekin,ir)
C- Finds envelope function parameters that match conditions on sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :controls how matching is done
Ci         :10s digit :
Ci         :0  vary rsm to match h(rsm,eh) to phi1,dphi1 at r1
Ci         :1  Vary eh  to match h(rsm,eh) to phi1,dphi1 at r1
Ci         :2  Vary rsm to match K.E. at r1 (no matching of slope)
Ci         :   hankel energy is input.
Ci         :3  Vary rsm and eh to match both slope and K.E. at r1
Ci         :4  Vary rsm to match K.E. and slope at r1, using
Ci         :   linear combination of h and hdot to match slope.
Ci         :   hankel energy is input.
Ci         :10 Try mode 0.  If unsuccessful, switch to mode 1
Ci         :11 Try mode 1.  If unsuccessful, switch to mode 0
Ci         :100s digit
Ci         :0 always return
Ci         :1 abort if match in log. derivative unsuccessful
Ci         :2 abort if match in K.E. unsuccessful
Ci   l     :l quantum number
Ci   rsmin :when rsm allowed to vary, lower bound for rsm
Ci   rsmax :when rsm allowed to vary, upper bound for rsm
Ci   emin  :when eh allowed to vary, lower bound for eh
Ci   emax  :when eh allowed to vary, upper bound for eh
Ci   r1    :first radius at which to perform matching
Ci   r2    :second radius at which to perform matching (not used now)
Ci   phi1  :wave function at r1
Ci   dphi  :derivative of phi1 at r1
Ci   phi2  :(mode=2,3) K.E. at r1
Ci         :(not used now) wave function at r2
Ci   dphi2 :derivative of phi2 at r2 (not used now)
Cio Inputs/Outputs
Cio  rsm   :envelope smoothing radius for smoothed hankel
Cio        :mode 0: rsm is output
Cio        :mode 1: rsm is input
Cio  eh    :envelope energy
Cio        :mode 0: eh is input
Cio        :mode 1: eh is output
Co Outputs
Co  ir     :information description
Co         :0  successful match first mode choice
Co         :1  successful match second mode choice
Co         :-1 failed to match log derivative
Co         :-2 failed to match K.E.
Co         :-3 maximum iterations exceeded
Co  ekin   :(modes 2,3) kinetic energy at MT boundary
Cl Local variables
Cl   dxrsm :rsm step length used for two-variable search
Cr Remarks
Cr   Properties of h(rsm,eh):
Cr   dh/h is a monotonically increasing function of eh for fixed rsm;
Cr   it is also a monotonically increasing function of rsm for fixed eh.
Cr   However, the kinetic energy is not.
Cr
Cr   Procedure for mode=3 (Vary rsm and eh to match dh/h and K.E. at r1)
Cr   1. Setup.  Find some (rsm,eh) pair that satisfies constraint
Cr         dphi1/phi1=dh/h.  K.E. is ignored in this step.
Cr      a. Begin with rsm=rsmin.
Cr      b. Determine eh that matches dh/h.
Cr         If matching fails, increment rsm by dxrsm/2 and repeat step
Cr         until match found or rsm exceeds rsmax.
Cr      c. mtchre aborts or exits unless a match is found
Cr   2. Iteratively search for K.E. match
Cr      a. Begin with rsm from step 1.  step 1 guarantees that there
Cr         is at least one (rsm,eh) pair that satisfies slope constraint
Cr      b. Iteratively search for rsm that matches K.E.
Cr         Iteration proceeds by regula falsi, which accepts as input
Cr         some (rsm,K.E) pair, returning a new rsm for next iteration.
Cr         In each step, rsm is given.
Cr         1  determine eh by matching slope constraint.
Cr         2  if slope constraint fails, return with (rsm,eh) pair that
Cr            best satisfies constraint (see 4)
Cr         3  compute K.E. for given (rsm,eh)
Cr         4  Follow the (rsm,eh) pair that best satisfies constraint.
Cr            (prepares for possible failure in a future iteration)
Cr         5  call rfalsi for a new estimate for rsm, or until
Cr            iteration converges.
Cr
Cu Updates
Cu   16 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,l,ir
      double precision dphi1,dphi2,eh,phi1,phi2,r1,r2,rsm,rsmin,rsmax,
     .  emin,emax,ekin
C ... Local parameters
      integer iter,stdo,lgunit,IPRT1,IPRTW,ipr,maxit,mode0,mode2,ir1,
     .  ipass
      double precision tol,eh0,dxrsm
      parameter (tol=1d-12,IPRT1=100/1,IPRTW=50/1,maxit=50,dxrsm=.05d0)
      double precision xnow,alfa,beta,xclose,bclose,eclose,tclose,wk(12)

      stdo = lgunit(1)
      call getpr(ipr)
      mode0 = mod(mode,100)
      mode2 = mod(mode/100,10)
      call isanrg(l,0,8,'mtchre','l',.true.)

C ... Start of matching in current mode
   10 continue

C --- Vary rsm to match phi to h ---
      if (mode0 .eq. 0 .or. mode0 .eq. 10) then
        call info2(IPRT1,0,-1,
     .  ' mtchre: match sm H(l=%i) to phi''/phi=%,1;3d varying rsm ...',
     .    l,dphi1/phi1)
        call mtchr2(0,l,rsmin,rsmax,(rsmin+rsmax)/2,r1,phi1,dphi1,rsm,
     .    eh,ekin,ir)
C   ... rfalsi unable to bound the search : handle error
        if (ir .lt. 0) then
          if (ipr .ge. IPRTW) write (stdo,
     .      '('' mtchre (warning) matching failed varying rsm'')')
          if (mode0 .eq. 10) then
            mode0 = 1
            goto 10
          endif
        else
          call info2(IPRT1,0,0,' found rsm=%;4d',rsm,0)
        endif

C --- Vary eh to match phi to h ---
      elseif (mode0 .eq. 1 .or. mode0 .eq. 11) then
        call info2(IPRT1,0,-1,
     .    ' mtchre: match sm H(l=%i) to phi varying eh ...',l,0)
        call mtchr2(1,l,emin,emax,(emin+emax)/2,r1,phi1,dphi1,rsm,eh,
     .    ekin,ir)
C   ... rfalsi unable to bound the search : handle error
        if (ir .lt. 0) then
          if (ipr .ge. IPRTW) write (stdo,
     .      '('' mtchre (warning) matching failed varying eh'')')
          if (mode0 .eq. 10) then
            mode0 = 0
            goto 10
          endif
        else
          call info2(IPRT1,0,0,' found eh=%;4d',eh,0)
        endif

C --- Vary rsm (and eh) to match K.E. (and log der) to h ---
C     See Remarks for description of procedure for mode3
      elseif (mode0 .ge. 2 .and. mode0 .le. 4) then

C   ... debugging ...
C       open(66,file='out')
C       do  eh = -5d0, -.02d0, .25d0
C       do  xnow = .8d0, 3.5d0, .1d0
C         call mtchae(2,xnow,eh,l,r1,phi1,phi1,0d0,0d0,alfa,ekin)
C         write(66,'(3f15.8)') xnow, eh, alfa
C       enddo
C       enddo
C       close(66)
C       stop

C   ... Setup: xnow is current guess for rsm.  Initialize w/ min val.
        xnow = min(rsmin,rsmax)
C       dxmx = max(rsmin,rsmax) - xnow
        if (xnow .lt. 0) call rx1('mtchr2: bad range, rsm = ',xnow)
        xclose = 0
        bclose = 9d9
C       Not used now, but keep anyway
        ipass = 1
        iter = 0
C       Initial value for rfalsi
        ir = 0
C       call pshpr(min(ipr,1))
        call pshpr(ipr-10)
      if (ipr .ge. IPRT1) write(stdo,261)
  261 format(' l  it  ir',6x,'Rsm',9x,'Eh',7x,'slope/val',6x,'K.E.',5x,
     .  'target K.E.')

C   ... mode 3-specific setup: find some (rsm,eh) satisfying slope cond.
C       Loop through rsm(=xnow) in uniform increments;
C       Find first rsm which can match dh/h to dphi1/phi1
C       K.E. is ignored in this first step
        if (mode0 .eq. 3) then
  110     continue
          eh0 = (emin+emax)/2
          call mtchr2(1,l,emin,emax,eh0,r1,phi1,dphi1,xnow,eh,ekin,ir1)
C         No match found for this rsm ... increment rsm
          if (ir1 .lt. 0) then
C           No match for any rsm ... give up
            if (xnow .lt. min(rsmin,rsmax)-1d-6 .or.
     .          xnow .gt. max(rsmin,rsmax)+1d-6 .or.
     .          rsmin .eq. rsmax) then
              ir = -1
              if (mode2 .ge. 1) call fexit2(-1,111,
     .          ' Exit -1 : mtchre : failed to match phi''/phi=%,1;3d'//
     .          ' to envelope, l=%i',dphi1/phi1,l)
              return
            endif
            xnow = xnow + dxrsm
            goto 110
          endif
        endif
C       End of mode 3-specific setup

C   ... Iteratively try to match K.E. with mode-specific constraints
  120   continue
C       mode 2-specific : No matching of slope
        if (mode0 .eq. 2) then

C       mode 3-specific : eh determined from dphi1/phi1=dh/h
        elseif (mode0 .eq. 3) then
          eh0 = eh
          call mtchr2(1,l,emin,emax,eh0,r1,phi1,dphi1,xnow,eh,ekin,ir1)
C         Failed to match log derivative.  Resort to best prior case
          if (ir1 .lt. 0) then
            if (iter .eq. 0) then
              ir = -1
              if (mode2 .ge. 1) call rxi('mtchre: '//
     .          'failed to match phi to envelope, mode',mod(mode,100))
              return
            endif
C           Flags that energy is near boundary point
C           ipass = 2
            ir = -1
            if (iter .ne. 0) ir = -2
            goto 80
          endif
        endif

C   ... K.E. for xnow = current guess for rsm.
        if (mode0 .ne. 4) then
C         phi1 only dummy here
          call mtchae(2,xnow,eh,l,r1,phi1,phi1,0d0,0d0,alfa,ekin)
        else
          call mtchae(3,xnow,eh,l,r1,phi1,dphi1,0d0,0d0,alfa,ekin)
        endif
        beta = ekin - phi2
C        print *, '!!'
C        beta = xnow
        rsm = xnow
C       Keep running track of closest approach in case no match found
        if (abs(beta) .lt. bclose) then
          xclose = xnow
          bclose = abs(beta)
          eclose = eh
          tclose = ekin
        endif
        call pshpr(0)
        call rfalsi(xnow,beta,tol,tol,tol/10,dxrsm,34,wk,ir)
        call poppr
        if (ir .gt. 1) call rxi('bug in mtchre, ir=',ir)
        iter = iter+1
        if (ipr .ge. IPRT1)
     .    write (stdo,60) l,iter,ir,wk(1),eh,alfa,beta+phi2,phi2
   60   format(i2,i4,i4,5f12.6)

C   ... rfalsi unable to find K.E... Use closest point
        if (xnow .lt. min(rsmin,rsmax)-1d-6 .or.
     .      xnow .gt. max(rsmin,rsmax)+1d-6) then
          ir = -2
          goto 80

C   ... rfalsi either has converged or requires another iteration
        else
          if (iter .lt. maxit .and. ir .lt. 0) goto 120
          if (iter .ge. maxit) ir = -3
          if (iter .ge. maxit) goto 80
          rsm = xnow
        endif
      else
        call rxi('mtchre: bad mode,',mode)
      endif

C ... Cleanup and exit
      if (ir .lt. 0 .and. mode2 .eq. 1) call rxi(
     .  'mtchre: failed to match phi to envelope, mode',mod(mode,100))
      ir = 0
      if (mode0 .ne. mod(mode,100)) ir = 1
      call poppr
      return

C --- Handle mode2 when matching failed ---
C     ir should be set before jumping here
   80 continue
      xnow = xclose
      eh0  = eclose
      ekin = tclose

C      call poppr
C      call mtchr2(1,l,emin,emax,eh0,r1,phi1,dphi1,xnow,eh,ir1)
C      if (ir1 .lt. 0) then
C        call rx('bug in mtchre: cannot match log deriv')
C      endif

C ... Closest point at boundary point in rsm or eh; no further search
      if (xclose .eq. min(rsmin,rsmax) .or.
     .    xclose.eq.max(rsmin,rsmax) .or. ipass .eq. 2
     .   .or. .true.) then

C       get K.E. at closest point.  phi1 only dummy here
        rsm  = xclose
        eh   = eclose
        ekin = tclose
        call mtchae(2,rsm,eh,l,r1,phi1,phi1,0d0,0d0,alfa,beta)

C ... Not implemented: find K.E. closest to target
      else
        print *, xclose
        call rx('mtchre : not implemented')
      endif

      call poppr

      if (ir .eq. -2)
     . call info5(IPRTW,0,0,' mtchre (warning) failed to match K.E. '//
     .  'of sm H(l=%i):  sought %;3d; best value %;3d',l,phi2,ekin,0,0)
      if (ir .eq. -1) call rxi(
     .  'mtchre: failed to match phi to envelope, mode',mod(mode,100))

      end

      subroutine mtchr2(mode,l,x1,x2,x0,r,phi,dphi,rsm,eh,ekin,info)
C- Match rsm or eh to value and slope of phi at surface of sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :controls how matching is done
Ci         :0 vary rsm to match h(rsm,eh) to phi,dphi at r
Ci         :  where phi=val, dphi=slope
Ci         :1 Vary eh  to match h(rsm,eh) to phi,dphi at r
Ci         :  where phi=val, dphi=slope
Ci         :2 Vary rsm to match h(rsm,eh) to phi,dphi at r
Ci         :  where phi=val, dphi=K.E.
Ci   l     :l quantum number
Ci   x1    :lower bound for rsm or eh
Ci   x2    :upper upper for rsm or eh
Ci   x0    :initial guess for rsm or eh (-99 => not used)
Ci   r     :first radius at which to perform matching
Ci   phi   :wave function at r
Ci   dphi  :derivative of phi at r
Cl Local variables
Cio Inputs/Outputs
Cio  rsm   :envelope smoothing radius for smoothed hankel
Cio        :mode 0: rsm is output
Cio        :mode 1: rsm is input
Cio  eh    :envelope energy
Cio        :mode 0: eh is input
Cio        :mode 1: eh is output
Cio  ir    :error
Co Outputs
Co  ekin   :Kinetic energy of h(rsm,eh,r)
Co  info   :information description
Co         :0 successful match first mode choice
Co         :1 successful match second mode choice
Co         :-1 (min,emax) do not bound search
Co         :-2 root finding unsuccessful after maxit iterations
Cr Remarks
Cr
Cu Updates
Cu   16 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,l,info
      double precision dphi,eh,phi,r,rsm,x1,x2,x0,ekin
C ... Local parameters
      integer ir,iter,stdo,lgunit,IPRT,ipr,maxit
      double precision xnow,dxmx,alfa,beta,tol,wk(12)
      parameter (tol=1d-12,IPRT=100/1,maxit=50)

      stdo = lgunit(1)
      call getpr(ipr)
C     call pshpr(ipr-10)
C     ipr=100

C     Mode-indpendent setup
      iter = 0
      ir = 0
      info = 0

      if (ipr .ge. IPRT) write(stdo,261) mode,l,r
  261 format(' mtchr2 mode',i2,'  l =',i2,'  r =',f10.6,
     .  /' l  it  ir',6x,'Rsm',9x,'Eh',8x,'phi''/phi     target')

C ... Vary rsm to match phi to h.  Do iteratively:
      if (mode .eq. 0 .or. mode .eq. 2) then
C       Setup
        xnow = min(x1,x2)
        dxmx = max(x1,x2) - xnow
        if (xnow .lt. 0) call rx1('mtchr2: bad lower rsm : ',xnow)
  100   continue
C
C       Match phi,dphi to hs,hsdot; find point where hsdot=0.
C       call mtchae(1,xnow,eh,l,r,phi,dphi,0d0,0d0,alfa,beta)
C       alfa = alfa/phi
C       beta = beta/phi
C       Match dphi/phi to dhs/hs. dhs/hs is monotonic in rsm for rsm<r
        call mtchae(2,xnow,eh,l,r,phi,dphi,0d0,0d0,alfa,ekin)
        beta = alfa - dphi/phi

C     mode=2 : match phi to hs; match dphi to K.E. of hs
C     not ready

        rsm = xnow
        call pshpr(0)
        call rfalsi(xnow,beta,tol,tol,tol/10,dxmx,10,wk,ir)
        call poppr
        if (ir .gt. 1) call rxi('bug in mtchr2, ir=',ir)
        iter = iter+1
        if (ipr.ge.IPRT)write (stdo,60) l,iter,ir,wk(1),eh,alfa,dphi/phi
   60   format(i2,i4,i4,2f12.6,1x,1p,4e12.3)

C  ... after 2rd iteration try estimate x0, if supplied
        if (iter .eq. 2 .and. x0 .gt. min(x1,x2)
     .      .and. x0 .lt. x2) xnow = x0

C   ... rfalsi unable to bound the search : handle error
        if (xnow .lt. min(x1,x2)-1d-6 .or.
     .      xnow .gt. max(x1,x2)+1d-6) then
          info = -1
C   ... rfalsi either has converged or requires another iteration
        else
          if (iter .lt. maxit .and. ir .lt. 0) goto 100
          if (iter .ge. maxit) info = -2
          rsm = xnow
        endif

C ... Vary eh iteratively to match log deriv of phi to that of h.
      elseif (mode .eq. 1) then

C       Setup:
        xnow = min(x1,x2)
        dxmx = max(x1,x2) - xnow
        if (xnow .ge. 0) call rx1('mtchr2: bad upper eh : ',xnow)
  200   continue

C       Try and match phi to hs,hsdot, and find point where hsdot=0.
C        call mtchae(1,rsm,xnow,l,r,phi,dphi,0d0,0d0,alfa,beta)
C        alfa = alfa/phi
C        beta = beta/phi
C
C       Try and match dphi/phi to alfa=dhs/hs.
C       NB: dhs/hs is monotonic in eh for rsm<r
        call mtchae(2,rsm,xnow,l,r,phi,dphi,0d0,0d0,alfa,ekin)
        beta = alfa - dphi/phi

        eh = xnow
C       call pshpr(100)
C       if (iter .eq. 26) then
C         print *, iter
C       endif
        call pshpr(0)
        call rfalsi(xnow,beta,tol,tol,tol/10,dxmx,10,wk,ir)
        call poppr
        if (ir .gt. 1) call rxi('bug in mtchr2, ir=',ir)

        iter = iter+1
        if (ipr.ge.IPRT) write(stdo,60)l,iter,ir,rsm,wk(1),alfa,dphi/phi
C        if (ipr.ge.IPRT)
C     .    write(stdo,60)l,iter,ir,rsm,wk(1),alfa,dphi/phi,
C     .    xnow-wk(1),alfa-dphi/phi

C  ... after 2rd iteration try estimate x0, if supplied
        if (iter .eq. 2 .and. x0 .gt. min(x1,x2)
     .      .and. x0 .lt. x2) xnow = x0

C   ... rfalsi unable to bound the search : handle error
        if (xnow .lt. min(x1,x2)-1d-6 .or.
     .      xnow .gt. max(x1,x2)+1d-6) then
          info = -1
C   ... rfalsi either has converged or requires another iteration
        else
          if (iter .lt. maxit .and. ir .lt. 0) goto 200
          if (iter .ge. maxit) info = -2
          eh = xnow
        endif

        if (ipr.ge.IPRT) print '('' exit mtchr2 info ='',i2)',info
C        if (info .eq. -1) then
C        print *, 'exit mtchr2 info=',info
C        endif

C ... Vary rs to match K.E. of h to dphi.  Do iteratively, matching
C     phi to hs; find point where K.E. matches.
      elseif (mode .eq. 2) then

      else
        call rxi('mtchr2 : bad mode',mode)
      endif


C     call poppr

      end


