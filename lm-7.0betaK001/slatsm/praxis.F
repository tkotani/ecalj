      subroutine praxis(f,n,x,t,h,eps,scbd,ktm,illc,ipm,lq1,lq2,lc,v,d,
     .                  y,z,q0,q1,tvec,fx,ir)
C- Minimisation of a function of n variables without derivatives
C-----------------------------------------------------------------------
Cr Praxis is a FORTRAN translation of the ALGOL algorithm given by
Cr R.P.Brent in Chapter 7 of "Algorithms for minimisation without
Cr derivatives", Prentice Hall, 1973 with the following modifications:
Cr 1. allows 1-D minimisation
Cr 2. allows user-initialisation of starting search directions
Cr 3. allow an upper limit to the number of Powell iterations to be done
Cr 4. incorporates switches to reduce the number of function calls for
Cr    a function which is very nearly quadratic. These control the level
Cr    of approximation in the line minimisations.
Cr 5. Includes an early termination switch, ir, passed back from the
Cr    function evaluation routine: if f(x,n,ir) ever sets ir to -1
Cr    then praxis immediately terminates returning ir = -1. This feature
Cr    is useful if the function calls involve operations of large
Cr    suites of programs and/or shell scripts and work by building up
Cr    a look-up table of function values.
Ci Inputs (see type declarations for cast):
Ci  f    name of the function to be minimised, called as f(x,n,ir)
Ci  n    the dimension of the function f
Ci  x    a vector containing the argument to f; on input, an
Ci       initial value
Ci  t    tolerance within which to find the minimum. For an analytic
Ci       function suffering only machine rounding error this can be
Ci       typically 10^-5, but for total energy calculations it should
Ci       not be smaller than the error from noise in the total energy.
Ci  h    maximum step length for line minimisations. If lq2 is set (see
Ci       below) then h is the step length used, otherwise the program
Ci       chooses an optimum step length from ||x||, eps, h and t. In
Ci       this case h should be the estimated distance to the minimum
Ci       from the intitial estimate of x.
Ci eps   square root of the machine precision: smallest relative spacing
Ci       This can be set larger than actually if there is noise in the
Ci       total energy.
Ci scbd  normally set to 1. If the algorithm could benefit from scaling
Ci       the axes then set to 10. See Brent p. 131
Ci ktm   number of iterations without improvement to perform before
Ci       terminating. Set to zero for a strictly quadratic function or
Ci       if function calls are very expensive, otherwise 1 is usually
Ci       quite adequate and 4 is very cautious.
Ci illc  normally F on input. Set to T if problem may be ill-conditioned
Ci ipm   a maximum number of basic Powell procedures to do, then
Ci       terminate. For truly quadratic function, can be set equal to n.
Ci       If the conjugate directions are known then ipm=1 will find the
Ci       minimum of a quadratic function in n line minimisations.
Ci lq1   A switch to limit the line searches to no more than 4 function
Ci       calls: 3 to estimate the minimum quadratically and one to
Ci       evaluate the function at the estimated minimum
Ci lq2   van Schilfgaarde's kamikaze switch: limits the line searches
Ci       to 3 function calls, the function value at the minimum being
Ci       estimated by the quadratic interpolation. If the function is
Ci       truly quadratic then setting lq2 and ipm=n will find the
Ci       minimum with the least possible number of function calls.
Ci       There is a strong likelihood that the algorithm will diverge
Ci       otherwise and there is no trap for this. Use this switch at
Ci       your own peril!!
Ci  lc   if set, then a set of conjugate directions (search directions)
Ci       are entered. If F, then v_ij are set to unit orthogonal
Ci       directions in praxis
Ci  v    a matrix whose columns are search directions for the
Ci       minimisation. As the algorithm proceeds these become the
Ci       "conjugate directions" or eigenvectors of the hessian matrix
Ci  d    a vector of second derivatives: eigenvalues of the hessian
Ci       initialised to zero in praxis and then determined by
Ci       successive line minimisations and singular value decompositions
Ci  y,z,q0,q1,tvec are work arrays
Co Outputs:
Co  x    the estimated minimum of f
Co  fx   the value f at the minimum
Co  v_ij and d contain information on second derivates, the hessian
Co       matrix A may be estimated from V(dagger)AV=diag(d)
Co  ir   returned by f: see point 5, above
Cp Procedures called by praxis:
Cp   Translated from Brent's ALGOL and included in this suite:
Cp         flin, lmin, quad, sortdv
Cp   Blas: daxpy, dcopy, dscal
Cp   Eispack: minfit
Cp   Local procedures:
Cp       iprint() returns value of print verbosity:
Cp       iprint() < 20: no printing;
Cp                = 20-40: increasing normal verbosity;
Cp                > 40: debugging
Cp       fexit: handles error exit after call to minfit in praxis
Cp       ran1() returns a single precision random number between 0 and 1
Cp       ran1in(iseed) initialises ran1().
C-----------------------------------------------------------------------
C     implicit none
      integer n,ktm,ipm,ir
C passed arguments
      double precision f,x(n),fx,t,h,eps,scbd,v(n,n),d(n),
     .                 y(n),z(n),q0(n),q1(n),tvec(n)
      external f
      logical illc,lc,lq1,lq2
C local
      integer nl,nf,kl,kt,i,j,k,k2,ipr,ierr,ipl,nits,is,iprint,i1mach
      double precision s,sl,dn,dmin,f1,lds,ldt,sf,df,qf1,qd0,qd1,
     .                 m2,m4,small,vsmall,large,vlarge,ldfac,t2
      double precision dsqrt,dabs,ax,cx,fa,fc,tin,hin
      real ran1
      character st*100

      call getpr(ipr)
      if (n .eq. 1) then
        ax = x(1)
        cx = x(1)+h
        call pshprt(0)
        call togprt
        call brentm(ax,x,cx,fa,fx,fc,f,1,t,eps,ipm,ir)
        if (ir .gt. 0 .and. (ir .le. ipm .or. ipm .eq. 0)) ir = 0
        call togprt
        call popprt
        return
      endif

C initialise machine dependent numbers
      small = eps**2
      vsmall = small**2
      large = 1d0/small
      vlarge = 1d0/vsmall
      m2 = dsqrt(eps)
      m4 = dsqrt(m2)
C seed random number generator
      call ran1in(-1)
C other initialisation
      ir = 0
      ipl = 1
      if (illc) then
        ldfac = 0.1
      else
        ldfac = 0.01
      endif
      nl = 0
      kt = 0
      tin = t
      hin = h
      t2 = small + dabs(t)
      t = t2
      dmin = small
C Brent has ..
C     if (h .lt. 100*t .and. .not.lq2) h = 100*t
      ldt = h
      if (iprint() .ge. 60) then 
        print *,'praxis choosing initial step length ...'
        call awrit5('small=%g,(tol,step-length) Before: (%g,%g)'//
     .              ' Now: (%g,%g)',' ',120,i1mach(2),
     .              small,tin,hin,t,ldt)
      endif
      qd0 = 0
      d(1) = 0
      call dcopy(n,x,1,q1,1)
C initialise conjugate directions to unit vectors
      if (.not. lc) then
        call dcopy(n*n,0d0,0,v,1)
        call dcopy(n,1d0,0,v,n+1)
      endif
C --- First call to function evaluation ---
      if (ipr .ge. 30) then
        is = 0
        st = ' '
        call bin2a('Praxis: start at x=',0,0,0,1,0,80,st,is)
        do  30  i = 1, n
   30   call bin2a(' ',1,0,x(i),4,0,100,st,is)
        print *, st(1:is)
      endif
      nf = 1
      fx = f(x,n,ir)
      if (ir .eq. -1) return
      qf1 = fx
C --- begin main loop ---
10000 continue
      sf = d(1)
      s = 0
      d(1) = 0
C ---  Minimise along the first direction ---
      nits = 2
      if (lq1) nits = 0
      if (n .eq. 1) then
        nits = 10
        lq2 = .false.
      endif
      call lmin(f,n,x,fx,v,1,qd0,qd1,q0,q1,d(1),dmin,.false.,s,fx,lq2,
     .         h,t,ldt,eps,nits,nf,nl,tvec,ir)
      if (ir .eq. -1) return
      if (n .eq. 1) goto 30000
      if (s .le. 0) call dscal(n,-1d0,v,1)
      if (sf.le.0.9*d(1).or.0.9*sf.ge.d(1)) call dcopy(n-1,0d0,0,d(2),1)
C --- Do n basic Powell procedures ---
      do  10  k = 2, n
        call dcopy(n,x,1,y,1)
        sf = fx
        illc = illc .or. (kt .gt. 0)
20000   continue
        kl = k
        df = 0
        if (illc) then
C ...     random step to get off resolution valley
          if (ipr .gt. 30) print *,'praxis: making random step ..'
          do  1  i = 1, n
            s = (0.1*ldt + t2*10**kt)*(ran1() - 0.5d0)
            z(i) = s
            call daxpy(n,s,x,1,v(1,i),1)
    1     continue
          fx = f(x,n,ir)
          if (ir .eq. -1) return
          nf = nf + 1
        endif
C --- Minimise along "non-conjugate directions" ---
        do  2  k2 = k, n
          sl = fx
          s = 0
          nits = 2
          if (lq1) nits = 0
          call lmin(f,n,x,fx,v,k2,qd0,qd1,q0,q1,d(k2),dmin,.false.,s,fx,
     .             lq2,h,t,ldt,eps,nits,nf,nl,tvec,ir)
          if (ir .eq. -1) return
          if (illc) then
            s = d(k2)*(s + z(k2))**2
          else
            s = sl - fx
          endif
          if (df .lt. s) then
            df = s
            kl = k2
          endif
    2   continue
C --- No success with illc = F; try once with illc = T ---
c      if (.not. illc .and. (df .lt. dabs(100*eps*fx))) then
c      if (.not. illc .and. (df .lt. eps)) then
        if (.not. illc .and. (df .lt. dabs(10*eps*fx))) then
          illc = .true.
          if (ipr .gt. 30)
     .      print *, 'praxis: problem seems ill-conditioned ..'
          goto 20000
        endif
        if (k .eq. 2 .and. ipr .gt. 50) then
          print *,'praxis: new vector of second derivatives:'
          write (*,100) (2*d(i), i = 1, n)
        endif
C --- Minimise along "conjugate" directions ---
        do  3  k2 = 1, k-1
          s = 0
          nits = 2
          if (lq1) nits = 0
          call lmin(f,n,x,fx,v,k2,qd0,qd1,q0,q1,d(k2),dmin,.false.,s,fx,
     .             lq2,h,t,ldt,eps,nits,nf,nl,tvec,ir)
          if (ir .eq. -1) return
    3   continue
C --- Get new "conjugate" direction x_n - x_0 ---
        f1 = fx
        fx = sf
        lds = 0
        do  4  i = 1, n
          sl = x(i)
          x(i) = y(i)
          sl = sl - y(i)
          y(i) = sl
          lds = lds + sl*sl
    4   continue
        lds = dsqrt(lds)
        if (lds .gt. small) then
C throw away direction kl and minimise along new "conjugate" direction
          if (ipr .gt. 50) print*,'praxis: shuffling directions ..'
          do  5  i = kl-1, k, -1
            call dcopy(n,v(1,i),1,v(1,i+1),1)
            d(i+1) = d(i)
    5     continue
          d(k) = 0
          call dcopy(n,y,1,v(1,k),1)
          call dscal(n,1d0/lds,v(1,k),1)
          if (ipr .gt. 50) then
            print*,'praxis: new conjugate direction u_n :'
            write (*,100) (v(i,k), i = 1, n)
          endif
          nits = 4
          if (lq1) nits = 0
          call lmin(f,n,x,fx,v,k,qd0,qd1,q0,q1,d(k),dmin,.true.,lds,f1,
     .             lq2,h,t,ldt,eps,nits,nf,nl,tvec,ir)
          if (ir .eq. -1) return
          if (lds .lt. 0) then
            lds = -lds
            call dscal(n,-1d0,v(1,k),1)
          endif
        endif
C --- Check stopping criteria ---
        ldt = ldfac*ldt
        if (ldt .lt. lds) ldt = lds
        t2 = 0
        do  6  i = 1, n
          t2 = t2 + x(i)**2
    6   continue
        t2 = m2*dsqrt(t2) + t
C --- See if step length exceeds half the tolerance ---
        if (ldt .gt. 0.5*t2) then
          kt = 0
        else
          kt = kt + 1
        endif
        if (kt .gt. ktm) goto 30000
C --- Check ipm switch ---
        ipl = ipl + 1
        if (ipm .gt. 0 .and. ipl .ge. ipm) then
          if (ipr .ge. 30) print*,
     .      'praxis: given number of Powell iterations completed.'
          goto 30000
      endif
C .. end basic Powell procedure
   10 continue
      if (ipr .ge. 30) then
        print*,'praxis: done n basic Powell procedures. Now x = '
        write (*,100) (x(i), i = 1, n)
        write (*,200) fx
      endif
C --- Try quadratic extrapolation in case stuck in a curved valley ---
      call quad(f,n,x,fx,dmin,nf,nl,h,t,ldt,eps,qf1,qd0,qd1,q0,q1,tvec,
     .          ir)
      dn = 0
      do  7  i = 1, n
        d(i) = 1d0/dsqrt(d(i))
        if (dn .lt. d(i)) dn = d(i)
    7 continue
      if (ipr .gt. 30) then
        print *,'praxis: new search directions:'
        do  8  j = 1, n
          write (*,100) (v(i,j),i=1,n)
    8   continue
      endif
      do  9  j = 1, n
        s = d(j)/dn
        call dscal(n,s,v(1,j),1)
    9 continue
      if (scbd .gt. 1) then
C --- Scale axes to try and reduce condition number ---
        s = vlarge
        do  12  i = 1, n
          sl = 0
          do  11  j = 1, n
            sl = sl + v(i,j)**2
   11     continue
          z(i) = dsqrt(sl)
          if (z(i) .lt. m4) z(i) = m4
          if (s .gt. z(i)) s = z(i)
   12   continue
        do  13  i = 1, n
          z(i) = z(i) / s
          if (z(i) .gt. scbd) z(i) = scbd
   13   continue
      endif
C                                           T
C --- Singular value decomposition of v - Golub & Reinsch (Eispack) ---
      do  14  i = 2, n
      do  14  j = 1, i-1
        s = v(i,j)
        v(i,j) = v(j,i)
        v(j,i) = s
   14 continue
      call minfit(n,n,n,v,d,0,tvec,ierr,tvec)
      if (ierr .ne. 0) call fexit(-1,1,'minfit in praxis: ierr=%d',ierr)
      if (scbd .gt. 1) then
C unscaling
        do  15  i = 1, n
          call dscal(n,z(i),v(i,1),n)
   15   continue
        do  18  i = 1, n
          s = 0
          do  16  j = 1, n
            s = s + v(i,j)**2
   16     continue
          s = dsqrt(s)
          d(i) = s * d(i)
          s = 1d0 / s
          do  17  j = 1, n
            v(i,j) = s * v(i,j)
   17     continue
   18   continue
      endif
C --- Eigenvalues of the approximate hessian from singular values of v
      do  19  i = 1, n
        if (dn*d(i) .gt. large) then
          d(i) = vsmall
        else if (dn*d(i) .lt. small) then
          d(i) = vlarge
        else
          d(i) = 1d0/(dn*d(i))**2
        endif
   19 continue
C --- Sort new eigenvalues and eigenvectors ---
      call sortdv(n,v,d)
      dmin = d(n)
      if (dmin .lt. small) dmin = small
      illc = m2*d(1) .gt. dmin
      if (ipr .gt. 30) then
        if (scbd .gt. 1) then
          print *,'praxis: scale factors:'
          write (*,100) (z(i), i = 1, n)
        endif
        print *,'praxis: eigenvalues of approximate hessian:'
        write (*,100) (d(i), i = 1, n)
      endif
      if (ipr .gt. 40) then
        print *, 'praxis: eigenvectors of approximate hessian:'
        do  20  j = 1, n
          write (*,100) (v(i,j), i = 1, n)
   20   continue
      endif
      goto 10000
30000 continue
      if (ipr .ge. 20) then
        print *,'praxis exiting:'
        print *,'x = '
        write (*,100) (x(i), i = 1, n)
        write (*,200) fx
        write (*,300) nl, nf
        write (*,400) dble(nf)/dble(nl)
      endif
  100 format(10(1pg10.3))
  200 format('f(x) = ',1pg13.6)
  300 format(/'number of line minimisations :',i6
     .       /'    number of function calls :',i6)
  400 format( '                         f/l :',f6.2)
      end
      subroutine lmin(f,n,x,fx,v,j,qd0,qd1,q0,q1,d2,dmin,fk,x1,f1,lq,
     .               h,t,ldt,eps,nits,nf,nl,tvec,ir)
C- Line minimisation from x along search direction.
C-----------------------------------------------------------------------
Ci  name of function: f;current point: x(1:n); function value at x: fx;
Ci  either: search direction, v(.,j), or (if j<1) points q0(.) and q1(.)
Ci  which define (with x) a space-curve along which which to estimate
Ci  a minimum, see Brent. p. 134,eq. (5.3) q0 and q1 are x' and x" and
Ci  qd0 and qd1 are d_0 and d_1;
Ci  dmin: a minimum estimate of d2 in case d2 near 0;
Ci  logical switch: fk: if T there is already known a value of f in the
Ci  search direction, namely f1 at a distance x1 along the s.d., if F
Ci  then f1 and x1 should be 0 and fx on entry; lq: logical, if T the
Ci  value of the function at the minimum is estimated from the
Ci  quadratic approximation used to find x; t: global tolerance;
Ci  h: max. step length; ldt: current step length; eps:  machine
Ci  precision; nits: maximium number of bisections in line search.
Cio incrementing # of function and line calls: nf,nl; tvec is workspace
Cr  lmin estimates the minimum of the function f along the search
Cr  direction v(.,j) starting at point x. That is, lmin finds the
Cr  minimum of the function (eq. 1.11, Brent) f(x_0 + lambda*v)
Co  On exit, x and fx are updated to the value and argument at the
Co  estimated minimum. The distance along the search direction,
Co  lambda, is returned in x1.
Co  if ir returned as -1, praxis will return to calling segment
C-----------------------------------------------------------------------
C     implicit none
      integer n,ir
      double precision f,x(n),v(n,n),d2,x1,f1,fx,h,t,ldt,dmin,eps,
     .                 qd0,qd1,q0(n),q1(n),tvec(n)
      external f
      integer nf,j,nits,nl
      logical fk,lq
C local
      double precision x2,xm,f0,f2,fm,d1,t2,s,sf1,sx1,f2q,d,m2,m4,small,
     .                 s1,t2fs,flin
      double precision dabs,dsqrt
      integer i,k,iprint,i1mach
      logical dz
C initialise
      small = eps**2
      m2 = dsqrt(eps)
      m4 = dsqrt(m2)
      sf1 = f1
      sx1 = x1
      k = 0
      xm = 0
      fm = fx
      f0 = fm
      dz = d2 .lt. eps
C find optimum step length, x1, along search direction
      s = 0d0
      do  1  i = 1, n
        s = s + x(i)**2
    1 continue
      s = dsqrt(s)
      s1 = s
      d = d2
      if (dz) d = dmin
      t2 = m4*dsqrt(dabs(fx)/d + s*ldt) + m2*ldt
      t2fs = t2
      s = m4*s + t
      if (dz .and. (t2 .gt. s)) t2 = s
      if (t2 .lt. small) t2 = small
      if (t2 .gt. 0.01*h) t2 = 0.01*h
      if (lq) t2 = h
      if (fk .and. (f1 .le. fm)) then
        xm = x1
        fm = f1
      endif
      if (iprint() .ge. 60) then 
        print 200, (v(i,j), i = 1, n)
        print *,'lmin, choosing initial step length ...'
        call awrit7('t=%g,h=%g,ldt=%g,||x||=%g(step=%g),s=%g,small=%g',
     .             ' ',120,i1mach(2),t,h,ldt,s1,t2fs,s,small)
        call awrit3('d2=%g,d=%g,step-length now:%g',' ',120,i1mach(2),
     .              d2,d,t2)
      endif
C if f(x + x1) along search direction not available, then evaluate it
      if (.not. fk .or. (dabs(x1) .lt. t2)) then
        if (x1 .ge. 0d0) then
          x1 = t2
        else
          x1 = -t2
        endif
        if (iprint() .ge. 60) then
          call awrit1('lmin calling flin. New step length:%g',' ',
     .      80,i1mach(2),x1)
        endif
        f1 = flin(f,n,x,v,x1,j,qd0,qd1,q0,q1,nf,tvec,ir)
        if (ir .eq. -1) return
      endif
      if (f1 .le. fm) then
        xm = x1
        fm = f1
      endif
10000 continue
      if (dz) then
C evaluate f at another point, x2, and estimate 2nd derivative
        if (f0 .lt. f1) then
          x2 = -x1
        else
          x2 = 2*x1
        endif
        if (iprint() .ge. 60) then
          call awrit1('lmin calling flin. New step length:%g',' ',
     .      80,i1mach(2),x2)
        endif
        f2 = flin(f,n,x,v,x2,j,qd0,qd1,q0,q1,nf,tvec,ir)
        if (ir .eq. -1) return
        if (f2 .le. fm) then
          xm = x2
          fm = f2
        endif
        d2 = (x2*(f1 - f0) - x1*(f2 - f0))/(x1*x2*(x1 - x2))
      endif
C estimate 1st derivative at x
      d1 = (f1 - f0)/x1 - x1*d2
      dz = .true.
C predict minimum ..
C here, Brent has if (d2 .le. small) then
      if (d2 .le. small*dabs(d1)) then
        if (d1 .lt. 0) then
          x2 = h
        else
          x2 = -h
        endif
      else
        x2 = -0.5d0*d1/d2
      endif
      if (dabs(x2) .gt. h) then
        if (x2 .gt. 0) then
          x2 = h
        else
          x2 = -h
        endif
      endif
20000 continue
C evaluate f at the predicted minimum
      f2q = f0 + d1*x2 + d2*x2*x2
      if (.not. lq) then
        if (iprint() .ge. 60) then
          call awrit1('lmin calling flin. New step length:%g',' ',
     .      80,i1mach(2),x2)
          call awrit3('      h=%g d1=%g d2=%g',' ',70,i1mach(2),h,d1,d2)
        endif
        f2 = flin(f,n,x,v,x2,j,qd0,qd1,q0,q1,nf,tvec,ir)
        if (ir .eq. -1) return
        if ((k .lt. nits) .and. (f2 .gt. f0)) then
C no success so try again ..
          k = k + 1
          if ((f0 .lt. f1) .and. (x1*x2 .gt. 0)) goto 10000
          x2 = 0.5d0*x2
          goto 20000
        endif
        if (iprint() .gt. 50) write (*,300) f2,f2q
      endif
      if (lq) f2 = f2q
C increment one-dimensional search counter
      nl = nl + 1
      if (f2 .gt. fm) then
        x2 = xm
      else
        fm = f2
      endif
C get new estimate of second derivative
      if (dabs(x2*(x2 - x1)) .gt. small) then
        d2 = (x2*(f1 - f0) - x1*(fm - f0))/(x1*x2*(x1 - x2))
      else if (k .gt. 0) then
        d2 = 0
      endif
      if (d2 .le. small) d2 = small
      x1 = x2
      fx = fm
      if (sf1 .lt. fx) then
        fx = sf1
        x1 = sx1
      endif
C Update x from linear search
      if (j .gt. 0) call daxpy(n,x1,v(1,j),1,x,1)
  200 format(' lmin: search direction:',256f8.4)
  300 format(' Lmin: f(x_min) = ',1pg15.8,
     .       '; quadratic estimate :',1pg15.8)
      end
      subroutine quad(f,n,x,fx,dmin,nf,nl,h,t,ldt,eps,qf1,qd0,qd1,q0,
     .                q1,tvec,ir)
C- Extrapolation down a long valley
C-----------------------------------------------------------------------
Ci f,n,x,fx,dmin,nf,nl,h,t,ldt,eps,qd0,qd1,q0,q1: see subroutine lmin
Ci tvec is workspace
Co x,fx,qf1,q0,q1; nl,nf incremented
Co if ir returned as -1, praxis will return to calling segment
Cr See p. 134, Brent. Here, l,qd0 and qd1 are respectively lambda,d_0
Cr and d_1. q0 and q1 are last two positions and x is current value.
Cr quad extrapolates quadratically down the curved valley defined by
Cr these three points, moves x to the minimum and updates q0 and q1.
C-----------------------------------------------------------------------
C     implicit none
      integer n,nf,nl,ir
      double precision f,x(n),fx,dmin,h,t,ldt,eps,qf1,qd0,qd1,q0(n),
     .                 q1(n),tvec(n)
      external f
C local
      double precision qa,qb,qc,s,l,dum,dsqrt
      integer i,iprint
C --- Initialise ---
      s = fx
      fx = qf1
      qf1 = s
      qd1 = 0

      do  1  i = 1, n
        s = x(i)
        l = q1(i)
        x(i) = l
        q1(i) = s
        qd1 = qd1 + (s - l)**2
    1 continue
      qd1 = dsqrt(qd1)
      l = qd1
      s = 0
      if (qd0 .gt. 0 .and. qd1 .gt. 0 .and. nl .ge. 3*n*n) then
        if (iprint() .gt. 30) print*,'quad: trying extrapolation ... '
        call lmin(f,n,x,fx,dum,0,qd0,qd1,q0,q1,s,dmin,.true.,l,qf1,
     .           .false.,h,t,ldt,eps,2,nf,nl,tvec,ir)
        if (ir .eq. -1) return
        qa = l*(l - qd1)/(qd0*(qd0 + qd1))
        qb = (l + qd0)*(qd1 - l)/(qd0*qd1)
        qc = l*(l + qd0)/(qd1*(qd0 + qd1))
      else
        fx = qf1
        qa = 0
        qb = 0
        qc = 1
      endif
      qd0 = qd1
      do  2  i = 1, n
        s = q0(i)
        q0(i) = x(i)
        x(i) = qa*s + qb*x(i) + qc*q1(i)
    2 continue
      end
      double precision function flin(f,n,x,v,l,j,qd0,qd1,q0,q1,nf,t,ir)
C- f(x_0 + l(x*_1 - x_0)) along j'th direction v(.,j)
C-----------------------------------------------------------------------
Ci nf, # of function calls so far; f, function name; n, # of dimensions;
Ci x, current position; v, conjugate directions; l, parameter; j, index
Ci to current search direction; qd0; qd1; q0; q1
Co t, x_0 + l(x*_1 - x_0); flin, f(t); nf incremented by one
Co if ir returned as -1, praxis will return to calling segment
Cr eq (1.11) Brent.
Cr flin places t in x if f returns ir=-1
C-----------------------------------------------------------------------
C     implicit none
      external f
      integer n,j,nf,ir
      double precision f,x(n),v(n,n),l,t(n),qd0,qd1,q0(n),q1(n)
C local
      double precision qa,qb,qc
      integer i

      if (j .gt. 0) then
C --- Linear search -- this branch called by lmin ---
        do  1  i = 1, n
          t(i) = x(i) + l*v(i,j)
    1   continue
      else
C --- Quadratic search -- this branch called by quad ---
        qa = l*(l - qd1)/(qd0*(qd0 + qd1))
        qb = (l + qd0)*(qd1 - l)/(qd0*qd1)
        qc = l*(l + qd0)/(qd1*(qd0 + qd1))
        do  2  i = 1, n
          t(i) = qa*q0(i) + qb*x(i) + qc*q1(i)
    2   continue
      endif
      nf = nf + 1
      flin = f(t,n,ir)
      if (ir .eq. -1) call dcopy(n,t,1,x,1)
      end
      subroutine sortdv(n,v,d)
C- widget to sort eigenvalues into descending order for praxis
C-----------------------------------------------------------------------
Ci n, dimension; v, eigenvectors; d, eigenvalues
Co v and d sorted in order of descending values in d
Cr bubble sort: no need to be fancy here.
C-----------------------------------------------------------------------
C     implicit none
      integer n
      double precision v(n,n),d(n)
C local
      integer i,j,k
      double precision s

      do  3  i = 1, n-1
        k = i
        s = d(i)
        do  1  j = i+1, n
          if (d(j) .gt. s) then
            k = j
            s = d(j)
          endif
    1   continue
        if (k .gt. i) then
          d(k) = d(i)
          d(i) = s
          do  2  j = 1, n
            s = v(j,i)
            v(j,i) = v(j,k)
            v(j,k) = s
    2     continue
        endif
    3 continue
      end

