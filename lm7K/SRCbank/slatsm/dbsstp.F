      subroutine dbsstp(x,y,dydx,yerr,neq,tol,hmin,h,yscal,mx,nmpseq,
     .ir,wk,ipr)
C- Bulirsch-Stoer integration step with monitoring of local accuracy
C  ---------------------------------------------------
Ci Inputs:
Ci   x:     initial value for independent variable x
Ci   y:     function values at current x
Ci   dydx:  y'(x,y) This array is DESTROYED on output but y'
Ci          at the starting x,y is saved in wk.
Ci   neq:   number of simultaneous equations
Ci   h:     step size algorithm will attempt
Ci   hmin:  smallest step size algorithm will use before giving up
Ci   tol:   tolerance of result, calculated from max(yerr/yscal)
Ci   yscal: scale each y relative to yscal to calculate error
Ci          and compare to tolerance.  Not used if yscal(1)=0.
Ci   mx:    maximum order of rational extrapolation polynomial
Ci          m=1 is just midpoint rule.
Ci   nmpseq:the sequence of divisions used by the midpoint rule.
Ci          dbsstp will take its default values if caller passes 0.
Ci          nmpseq(0): the maximum number of midpoint rule integrations
Ci          that can be invoked before giving up and reducing h.
Ci          nmpseq(1..): the number of divisions for each rule.
Ci   ir:    a flag telling dbsstp where it is in the current
Ci          integration step.  When starting a new integration
Ci          step, ir should be 0 (See Outputs, Remarks.)
Ci   wk:    a work array of dimensions 5+(4+mx)*neq.  On each call
Ci          wk(1..neq) holds dy/dx at current (x,y)
Ci          nint(wk(5)) contains the number of function calls requested
Ci          so far in this integration step
Ci   ipr:   verbosity
Co Outputs:
Co   ir,    a flag indicating what action to take (see Remarks)
Co    ir>0  dbsstp wants a new evaluation of dy/dx for
Co          the (x,y) it returns.  Caller must evaluate dy/dx
Co          and call dbsstp again, leaving ir untouched.
Co    ir=0  dbsstp has completed integration from xs to new x
Co    ir=-1 dbsstp did not converge to tol for any h>hmin.
Co   x,y:   (ir>0) is some intermediate values which caller uses
Co          to make a new dy/dx, and to call dbsstp again.
Co          (ir=0) x,y are replaced by x+h and y(x+h)
Co   yerr:  estimated errors in y(x+h)
Co   h:     a new suggested step size
C  Local variables (used in midpoint rule):
Cv   xs:    initial value of independent variable x
Cv   ir:    0, beginning of midpoint rule, or dbsstp converged
Cv         -1, returned when dbsstp did not converge even for h=hmin
Cv        <-1, dbsstp called with a bad wk array.
Cv          1  dbsstp is at the imp'th step of a midpoint rule.
Cv   imp:   index to current midpoint-rule, containing nseq(imp) points.
Cv   ip:    index to which point in current midpoint-rule
Cv   dydx:  as the midpoint rule progresses, dy/dx;
Cv          at close of midpoint rule, midpoint estimate of y(x+h)
Cv   wk:    a work array dimensioned 5+(3+mx)*neq, broken into
Cv          5 elements holding scalars and vectors of length neq.
Cv          Offset:
Cv          nv0: starting value of dy/dx
Cv          nv1: y at xs
Cv          nv2: y_k-1 (midpoint rule)
Cv          nv3: coffs for rational extrapolation
Cv              (equivalent to matrix d in Num. Recip. rzextr)
Cr Remarks:
Cr   bsstp integrates the equations from x to x+h with several midpoint
Cr   rules combined with rational function extrapolation to accelerate
Cr   convergence.  Midpoint rules with successively finer meshes are
Cr   invoked until the estimated error in the integrated y's is less
Cr   than abs(tol).
Cr   Case tol>0:
Cr     If it fails to converge to tol with the finest available mesh, 
Cr     dbsstp reduces the step size h and starts over.
Cr   Case tol=0:
Cr     dbsstp will use all the midpoint rules available, and return 
Cr     normally without adjusting step size.
Cr   Case tol<0:
Cr     dbsstp exits with a warning, returns ir=-2 and a suggested new
Cr     step size
Cb Bugs/improvements
Cb   mi,nseq0 are the default for nmpseq, taken from Numerical Recipes.
Cb   Improvements: copy work arrays when h is reduced.
C  ---------------------------------------------------
C     implicit none
      integer neq,mx,ipr,ir,nmpseq(0:1)
      double precision y(neq),dydx(neq),yscal(neq),yerr(neq),
     .wk((2+mx)*neq)
      double precision tol,hmin,h,x
      double precision b,b1,c,ddy,hi,swap,yy,v,h2,errmax,
     .hgrow,hshrnk
C Local parameters
      integer mi,mi0,m
      parameter(mi0=11,hgrow=1.1d0,hshrnk=.9d0)
      integer nseq0(0:mi0),nseq
      double precision xs,hfac,xx
      integer imp,ip,j,k,nv0,nv1,nv2,nv3,i1mach,nseq1,ncall
      data nseq0 /mi0,2,4,6,8,12,16,24,32,48,64,96/

C --- ir-independent setup ---
      mi = nmpseq(0)
      if (mi .le. 0) mi = mi0
      xs =  wk(1)
      imp = wk(2)
      ip  = wk(3)
      nseq = nint(wk(4))
      ncall = nint(wk(5))
      nv0 = 5
      nv1 = nv0+neq
      nv2 = nv1+neq
      nv3 = nv2+neq

C --- ir-dependent setup ---
      if (ir .eq. 0) then
        if (ipr .ge. 30) call awrit6
     .  (' dbsstp:  new step  neq=%i  x=%1;6g  h=%1;6g  tol=%1;6g'/
     .  /'  mx=%i  mi=%i',' ',80,i1mach(2),neq,x,h,tol,mx,mi)
        imp = 0
        xs = x
        ncall = 0
        do  8  j = 1, neq
          wk(nv0+j) = dydx(j)
    8   wk(nv1+j) = y(j)
      elseif (ir .eq. 1) then
C   ... check for bad imp
        if (imp .le. 0 .or. imp .gt. mi) then
          ir = -2
          return
        endif
        hi = h/nseq
        if (ip .eq. 1) goto 901
        goto 902
      else
        ir = -3
        return
      endif

C --- Entry point for a new Stoer-Bulirsch iteration ---
   11 continue
      imp = imp+1
      nseq = nseq0(imp)
      if (nmpseq(0) .ne. 0) nseq = nmpseq(imp)
      if (ipr .ge. 40)
     .call awrit3(' dbsstp:  begin mp rule %i (%i points)'//
     .'  ncall=%i',' ',80,i1mach(2),imp,nseq,ncall+1)

C   --- Midpoint rule ---
      hi = h/nseq
      ip = 0
      ip = ip+1
C   ... First point
      do  111  j = 1, neq
        wk(nv2+j) = wk(nv1+j)
  111 y(j) = wk(nv1+j) + hi*wk(nv0+j)
      x = xs + hi
C   ... Exit to let caller compute dydx(x,y) and jump back here (901)
      ir = 1
      goto 900
  901 continue
C   ... Loop over all succeeding points ip = 2 .. nseq
  113 ip = ip+1
      h2 = 2*hi
      do  112  j = 1, neq
        swap = wk(nv2+j) + h2*dydx(j)
        wk(nv2+j) = y(j)
        y(j) = swap
  112 continue
      x = x + hi
C     ... Exit to let caller compute dydx(x,y) and jump back here (902)
      ir = 1
      goto 900
  902 if (ip .gt. nseq) then
        ir = -3
        return
      endif
C   ... End of 113 loop
      if (ip .lt. nseq) goto 113
C   ... End result of modified midpoint rule
      do  114  j = 1, neq
  114 dydx(j) = (y(j) + wk(nv2+j) + hi*dydx(j))/2

C --- Estimate y, error as a rational function of h*h ---
C     (Equivalent to: call rzextr(imp,hi**2,dydx(1),y,yerr,neq,mx))
c       write(*,'(1x,a,i4,5f14.8)') 'm',imp,(dydx(j),j=1,min(neq,5))
C   ... First diagonal in tableau
      if (imp .eq. 1) then
        do  211  j = 1, neq
          y(j) = dydx(j)
          wk(nv3+j) = dydx(j)
          yerr(j) = dydx(j)
  211   continue
      else
C     ... Use at most mx previous points
        m = min(imp,mx)
C     ... Evaluate next diagonal in tableau
        do  214  j = 1, neq
          yy = dydx(j)
          v = wk(nv3+j)
          c = yy
          wk(nv3+j) = yy
          do  213  k = 2, m
            if (v .gt. 1d10) call rx('dbsstp: bad work array')
            nseq1 = nseq0(imp-k+1)
            if (nmpseq(0) .ne. 0) nseq1 = nmpseq(imp-k+1)
            b1 = v*(dble(nseq)/dble(nseq1))**2
            b = b1 - c
            if (b .ne. 0) then
              b = (c-v)/b
              ddy = c*b
              c = b1*b
            else
              ddy = v
            endif
C         ... for k<m, v will obtain from a previous i; not used for k=m
            v = wk(nv3+j+(k-1)*neq)
            wk(nv3+j+(k-1)*neq) = ddy
            yy = yy + ddy
  213     continue
          yerr(j) = ddy
          y(j) = yy
  214   continue
      endif
c       write(*,'(a,i4,5f14.8)') 'a',imp, (y(j),j=1,min(neq,5))
c       write(*,'(a,i4,5f14.8)') 'a',imp,(yerr(j),j=1,min(neq,5))

C ---   Error estimate; check for convergence ---
      errmax = 0
      if (yscal(1) .gt. 0) then
        do  12  j = 1, neq
   12   errmax = max(errmax,dabs(yerr(j)/yscal(j)))
      else
        do  13  j = 1, neq
   13   errmax = max(errmax,dabs(yerr(j)))
      endif
C   ... Normal exit when this step is completed
      if (errmax .lt. abs(tol)) then
        ir = 0
C     ... Tune next h to make last imp come as close to mx as possible
        if (imp .eq. mx) then
          hfac = hshrnk
        elseif (imp .eq. mx-1) then
          hfac = hgrow
        else
          nseq1 = nseq0(mx-1)
          if (nmpseq(0) .ne. 0) nseq1 = nmpseq(mx-1)
          hfac = dble(nseq1)/nseq
        endif
        x = xs + h
        if (ipr .ge. 20)
     .  call awrit5(' dbsstp: completed xs=%1;5d to %1;5d (%i'
     .  //' calls, %i mp)  hnext=%1;5d',
     .  ' ',80,i1mach(2),xs,x,ncall,imp,h*hfac)
        h = hfac * h
        return
      endif
      if (imp .lt. mi) goto 11

C --- Step failed for this h:  reduce h and try again ---
      if (tol .gt. 0) then
        xx = h/4/2**((mi-mx)/2)
        if (ipr .ge. 30)
     .  call awrit4(' dbsstp: err=%1;3g > tol=%1;3g.  Reduce h from'//
     .  ' %1;6g to %1;6g',' ',80,i1mach(2),errmax,tol,h,xx)
        h = xx
        imp = 0
        if (h .gt. hmin) goto 11
        ir = -1
        return
      endif
C ... Handle case tol <= 0
      ir = 0
      x = xs + h
      if (tol .eq. 0) then
        if (ipr .ge. 30)
     .  call awrit5(' dbsstp: completed xs=%1;5d to %1;5d (%i calls,'
     .  //' %i mp)  err=%1;3g',' ',80,i1mach(2),xs,x,ncall,imp,errmax)
      else
        xx = h/4/2**((mi-mx)/2)
        if (ipr .ge. 20)
     .  call awrit7(
     .  ' dbsstp: err(%1;2g)>tol(%1;2g) for x=%1;5d to %1;5d'//
     .  ' (%i calls, %i mp) hnext=%1;3g',
     .  ' ',80,i1mach(2),errmax,dabs(tol),xs,x,ncall,imp,xx)
        h = xx
        ir = -1
      endif
      return

C --- Save internal parms for a new call ---
  900 continue
      wk(1) = xs
      wk(2) = imp
      wk(3) = ip
      wk(4) = nseq
      wk(5) = ncall+1
      if (ipr .gt. 50)
     .call awrit3(' dbsstp: continuing mp rule %i, point %i  x=%1;6g',
     .' ',80,i1mach(2),imp,ip,x)

      end

