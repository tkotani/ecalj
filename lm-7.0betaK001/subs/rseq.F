      subroutine rseq(eb1,eb2,e,tol,z,l,nod,val,slo,v,g,q,a,b,rofi,nr,
     .  nre)
C- Solves radial wave equation for given BCs and number of nodes
C ----------------------------------------------------------------------
Ci Inputs:
Ci   eb1,eb2 lower and upper bounds for the eigenvalue
Ci   tol     tolerance: maximum relative error in energy;
Ci           absolute error if |energy| < 1
Ci   z       nuclear charge
Ci   a,b     mesh points given by rofi(i) = b [e^(a(i-1)) -1]
Ci   l       angular momentum
Ci   nod     number of nodes
Ci   val,slo BC'S for large component u(r)=g(r,1) with psi=(u/r)*ylm
Ci           val is u(rmax), i.e. (rmax * radial w.f.)
Ci           slo is radial derivative of g at rmax (rmax * radial w.f.)'
Ci   v       spherical potential (true potential excluding nuclear part)
Ci   rofi,nr radial mesh points, and number
Co Outputs:
Co   e       eigenvalue
Co   g       Wave function times r normalized so that int (g*g) dr = 1
Co           Large component in g(ir,1); small component = g(ir,2)
Co   q       integral of unnormalized product, that is of
Co           int (g*g) dr if g were normalized to input val,slo
Co   nre     index to rofi at which outward and inward solutions join
Cr Remarks:
Cr   Scalar relativistic version
Cr   Output wavefunction normalized to 1 ... so g(nr) .ne. val*wsr
Cr   Note: if r = b(exp(a*z)-1) then Jacobian dr=a(r+b) dz
Cr
Cu Updates
Cu   08 Feb 01 if rseq fails to converge, but the number of nodes is
Cu             correct, rseq returns with a warning instead of aborting
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l,nod,nr,nre
      double precision eb1,eb2,e,tol,z,val,slo,v(nr),g(nr,2),q,rofi(nr),
     .  a,b
C ... Local parameters
      integer ipr,k,k2,kc,nctp,nctp0,nit,nitmax,nod1,nod2,node
      double precision c,de,e1,e2,fac,fllp1,r,ratio,re,rhok,slo1,slo2,
     .  slop,tmcr,val1,val2,valu,wgt
C ... Speed of light, or infinity in nonrelativistic case
      common /cc/ c
C     double complex gc(nr,2),ec,valc,sloc

      call getpr(ipr)
      nitmax = 80
      if (ipr .ge. 110) write(*,815) z,l,nod,val,slo,eb1,eb2
  815 format(' RSEQ:  Z=',f5.1,'  l=',i1,'  nod=',i1,'  bc=',2f7.3,
     .   '  e1,e2=',2f8.2)
      e1 = eb1
      e2 = eb2
      call fctp0(l,nr,rofi,v,z,nctp0)
      if (ipr .ge. 120) print 301
  301 format(' nit l node nod nre kr',
     . '        e1            e              e2            de')
C --- Start iterations to find energy ---
      do  1  nit = 1, nitmax
        if (e .le. e1 .or. e .ge. e2) e = (e1 + e2)/2
C       call fctp(e,nctp,nctp0,xrim,xmin,nsave,l,rofi,v,z,nr,a,b)
        call fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)
        re = 15d0*rofi(nctp)
        nre = int(dlog(re/b + 1d0)/a + 1d0)
        nre = (nre/2)*2 + 1
        nre = max0(35,min0(nre,nr))
        valu = val
        slop = slo
        if (nre .lt. nr) valu = 1d-5
        if (nre .lt. nr) slop = -1d-5
        k2 = min(nre,30)
        if (nod .eq. 0) k2 = nre/3
        if (valu*slop .gt. 0d0 .and. nod .eq. 0) k2 = nre - 10
C       Integrate the scalar relativistic eqn inward from nre to kc
        call rsq2(e,l,z,v,nre,k2,valu,slop,g,val2,slo2,nod2,kc,
     .              a,b,rofi,nr)
C       Integrate the scalar relativistic eqn outward from origin to kc
        call rsq1(0,e,l,z,v,kc,g,val1,slo1,nod1,a,b,rofi,nr)

        node = nod1 + nod2
        if (node .ne. nod) then
          if (ipr .ge. 120 .or. nit .ge. nitmax-5)
     .    write(*,101) nit,l,node,nod,nre,kc,e1,e,e2
  101     format(2i3,4i4,3f15.7,1p,d13.4)
          if (node .gt. nod) e2 = e
          if (node .lt. nod) e1 = e
          e = (e1 + e2)/2
          goto 1
        endif
C   ... Calculate q = norm of wave function from trapezoidal rule and
C       de = estimated correction to eigenvalue
        ratio = val2/val1
        q = 0d0
        do  5  k = 2, kc
          q = q + (rofi(k)+b)*g(k,1)*g(k,1)
    5   continue
        q = q*ratio**2
        do  6  k = kc+1, nre
          q = q + (rofi(k)+b)*g(k,1)*g(k,1)
    6   continue
        q = a*(q - (rofi(nre)+b)*g(nre,1)*g(nre,1)/2)
        de = -val2*(slo2 - ratio*slo1)/q

        if (ipr .ge. 120 .or. nit .ge. nitmax-5)
     .    write(*,101) nit,l,node,nod,nre,kc,e1,e,e2,de
        if (de .gt. 0d0) e1 = e
        if (de .lt. 0d0) e2 = e
        e = e + de
C       Exit loop when de meets tolerance; eval found
        if (dabs(de/dmax1(dabs(e),1d0)) .lt. tol) goto 2
    1 continue
C --- Search for eigenvalue failed ---
      nit = nitmax+1
C     Fatal if node mismatch
      if (nod .ne. node) goto 99

C --- Normalize g ---
    2 continue
      fllp1 = l*(l+1)
      e = e - de
      do  8  k = 1, kc
        g(k,1) = g(k,1)*ratio
        g(k,2) = g(k,2)*ratio
    8 continue
      q = 0d0
      do  10  k = 2, nre
        r = rofi(k)
        wgt = (mod(k+1,2) + 1)*(r + b)
        tmcr = (c - (v(k) - 2d0*z/r - e)/c)*r
        rhok = g(k,1)*g(k,1)*(1d0 + fllp1/(tmcr*tmcr)) + g(k,2)*g(k,2)
        q = q + wgt*rhok
   10 continue
      q = (q - wgt*rhok/2)*a*2d0/3d0
      fac = 1d0/dsqrt(q)
      do  11  k = 1, nre
        g(k,1) = g(k,1)*fac
        g(k,2) = g(k,2)*fac
   11 continue
      do  12  k = nre+1, nr
        g(k,1) = 0d0
        g(k,2) = 0d0
   12 continue

C --- Possible warning or error exit ---
   99 continue
      if (ipr .ge. 110 .or. nit .ge. nitmax)
     .  write(*,701) e,q,nr,nre,kc,de
  701 format(' e=',f13.5,'  q=',f10.5,'   nr/nre/kc=',3i4,
     .  '   de=',1p,d10.2)
      if (nit .gt. nitmax) then
        if (nod .ne. node) then
          write(*,814) nitmax,l,nod,node,e
  814     format(' RSEQ : nit gt',i3,' and bad nodes for l=',i1,
     .      '.  Sought',i2,' but found',i2,'.  e=',1pd11.4)
          call rx('RSEQ: bad nodes')
        else
          write(*,816) l,nod,abs(de),e
  816     format(' RSEQ (warning) eval for l=',i1,' node=',i1,
     .      ' did not converge: de=',1pd8.2,' e=',1pd11.4)
        endif
      endif

C     ec = e
C     valc = val
C     sloc = slo
C     call zseq(ec,tol,z,l,nod,valc,sloc,v,gc,q,a,b,rofi,nr,nre)

      end
      subroutine rsq1(i0,e,l,z,v,kr,g,val,slo,nn,a,b,rofi,nr)
C- Integrate the scalar relativistic eqn outward to rofi(kr)
C ----------------------------------------------------------------
Ci Inputs:
Ci   i0      *If i0<5, integrate from origin to rofi(kr)
Ci           *Else, integrate from rofi(i0) to rofi(kr).
Ci            NB: in this case, g(i0-4..ir) must be input!
Ci   e        Energy
Ci   l,z,v    angular momentum, nuc. charge, potential (v_true)
Ci   a,b      mesh points given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rofi,nr  radial mesh points, and number
Ci   kr       integration from origin to rofi(kr)
Co Outputs:
Co   g        radial wave function times r (but see i0),
Co            normalized to unity.
Co   val,slo  value and slope of g at rofi(kr)
Co   nn       number of nodes encountered
Cr Remarks:
Cr   Boundary condition does not fix value of wave function near the
Cr   origin; integration can be scaled by an arbritrary factor
Cu Updates
Cu   21 Jun 04 Integration can start from point other than origin
C ----------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer i0,l,kr,nn,nr
      double precision e,z,v(nr),g(nr,2),val,slo,a,b,rofi(nr)
C Local parameters:
      integer ir,i1
      double precision d(2,3),zz,c,fllp1,r83sq,r1,r2,r3,h83,g0,s,sf,aa,
     .    f0,dg1,dg2,dg3,df1,df2,df3,r,drdi,phi,u,x,y,det,b1,b2
      equivalence (dg1,d(1,1)),(dg2,d(1,2)),(dg3,d(1,3)),
     .            (df1,d(2,1)),(df2,d(2,2)),(df3,d(2,3))
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

C      double complex gc(nr,2),ec,valc,sloc
C      ir = 0
C      if (ir .gt. 0) then
C        ec = e
C        call zsq1(ec,l,z,v,kr,gc,valc,sloc,a,b,rofi,nr)
C      endif

      nn = 0
      zz = z+z
      fllp1 = l*(l+1)
      r83sq = 64d0/9d0
      r1 = 1d0/9d0
      r2 = -5d0*r1
      r3 = 19d0*r1
      h83 = 8d0/3d0

C --- Approximate g,f by leading term near zero ----
      if (i0 .lt. 5) then
      g0 = 1
      if (z .lt. 0.9d0) then
        s = l+1
        sf = l
        f0 = l/c
      else
        aa = zz/c
        s  = dsqrt(fllp1 + 1d0 - aa*aa)
        sf = s
        f0 = g0*(s - 1d0)/aa
      endif
      g(1,1) = 0d0
      g(1,2) = 0d0
      do  2  ir = 2, 4
        r = rofi(ir)
        drdi = a*(r+b)
        g(ir,1) = (r**s)*g0
        g(ir,2) = (r**sf)*f0
        d(1,ir-1) = drdi*g(ir,1)*s/r
        d(2,ir-1) = drdi*g(ir,2)*sf/r
    2 continue
      endif

C --- Setup to integrate over rest of points ------
      if (i0 .lt. 5) then
        i1 = 5
C       dg1 = d(1,1)
C       dg2 = d(1,2)
C       dg3 = d(1,3)
C       df1 = d(2,1)
C       df2 = d(2,2)
C       df3 = d(2,3)
      else
        i1 = i0
        call dpzero(d,6)
        do  3  ir = i1-3, i1-1
        r     = rofi(ir)
        drdi  = a*(r + b)
        phi   = (e + zz/r - v(ir))*drdi/c
        u     = drdi*c + phi
        x     = -drdi/r
        y     = -fllp1*x*x/u + phi
C       det   = r83sq - x*x + u*y
C       b1    = g(ir-1,1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
C       b2    = g(ir-1,2)*h83 + r1*df1 + r2*df2 + r3*df3
C       g(ir,1) = (b1*(h83-x) + b2*u)/det
C       g(ir,2) = (b2*(h83+x) - b1*y)/det
C       if (g(ir,1)*g(ir-1,1) .lt. 0d0) nn = nn+1
        dg1   = dg2
        df1   = df2
        dg2   = dg3
        df2   = df3
        dg3   = u*g(ir,2) - x*g(ir,1)
        df3   = x*g(ir,2) - y*g(ir,1)
    3   continue
      endif

C --- Integrate over rest of points ------
      do  4  ir = i1, kr
        r     = rofi(ir)
        drdi  = a*(r + b)
        phi   = (e + zz/r - v(ir))*drdi/c
        u     = drdi*c + phi
        x     = -drdi/r
        y     = -fllp1*x*x/u + phi
        det   = r83sq - x*x + u*y
        b1    = g(ir-1,1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b2    = g(ir-1,2)*h83 + r1*df1 + r2*df2 + r3*df3
        g(ir,1) = (b1*(h83-x) + b2*u)/det
        g(ir,2) = (b2*(h83+x) - b1*y)/det
        if (g(ir,1)*g(ir-1,1) .lt. 0d0) nn = nn+1
        dg1   = dg2
        df1   = df2
        dg2   = dg3
        df2   = df3
        dg3   = u*g(ir,2) - x*g(ir,1)
        df3   = x*g(ir,2) - y*g(ir,1)
    4 continue
      val = g(kr,1)
      slo = dg3/(a*(rofi(kr) + b))

C     call prrmsh('g',rofi,g,nr,kr,2)
      end
      subroutine rsq2(e,l,z,v,nre,ncmin,val1,slo1,g,val,slo,nn,nc,
     .  a,b,rofi,nr)
C- Integrate the scalar relativistic eqn inward from nre to nc
C ----------------------------------------------------------------
Ci Inputs:
Ci   e     :energy
Ci   l     :angular momentum
Ci   z     :nuclear charge
Ci   v     :spherical potential = v_true (excluding nuclear 2*Z/r)
Ci   nre   :rofi(nre) = starting radius from which rsq2 integrates
Ci   ncmin: rsq2 integrates to cutoff nc, of which ncmin is lower bound
Ci   val1  :value of large component of g at nre
Ci   slo1  :slope of g at nre
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   rofi  :radial mesh points
Ci   nr    :leading dimension of g
Co Outputs:
Co   g     :radial wave function times r at points nc..nre
Co   nc    :rofi(nc) = radius to which rsq2 integrates
Co         :nc is the larger of ncmin and the position of the
Co         :first maximum encountered
Co   val   :value of g(nc) (large component)
Co   slo   :slope of g(nc)
Co   nn    :number of nodes found between nc and nre
Cu Updates
Cu   21 Jun 04 Attempt to extend integration also to outward.
Cu             Doesn't work.  Use rsq1
Cr Remarks:
Cr   Integrates inward from nre to nc
Cr   Cutoff nc is chosen at first maximum, but nc .ge. ncmin
C ----------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer l,nre,ncmin,nn,nc,nr
      double precision e,z,v(nr),val1,slo1,g(nr,2),val,slo,a,b,rofi(nr)
C Local parameters:
      integer i,ir,irp1,i2,i1,ifac
      double precision d(2,3),zz,c,fllp1,r83sq,r1,r2,r3,h83,ea,rpb,
     .      q,ag1,ag2,ag3,af1,af2,af3,gg,ff,vb,drdi,r,phi,u,x,y,
     .      dg1,dg2,dg3,df1,df2,df3,det,b1,b2
C Speed of light in common to  relativity
      common /cc/ c

C      double complex gc(nr,2),ec,valc,sloc,val1c,slo1c
C      ir = 0
C      if (ir .gt. 0) then
C        ec = e
C        val1c = val1
C        slo1c = slo1
C        call zsq2(ec,l,z,v,nre,ncmin,val1c,slo1c,gc,valc,sloc,nn,nc,
C     .    a,b,rofi,nr)
C      endif

      nn = 0
      zz = z + z
      fllp1 = l*(l + 1)
      r83sq =64d0/9d0
      r1    = 1d0/9d0
      r2    =-5d0/9d0
      r3    =19d0/9d0
      h83   =-8d0/3d0
      ifac = -1
C     if (ncmin .gt. nre) ifac = 1
      if (ncmin .gt. nre)
     .  call rx('rsq2 not implemented for outward integration')

C --- First point ------
      r      = rofi(nre)
      rpb    = r+b
      drdi   = a*rpb
      phi    = (e + zz/r - v(nre))*drdi/c
      u      = drdi*c + phi
      x      = -drdi/r
      y      = -fllp1*x*x/u + phi
      g(nre,1) = val1
      g(nre,2) = (slo1*drdi + x*val1)/u
      ag1    = slo1*drdi
      af1    = x*g(nre,2) - y*g(nre,1)
      ir     = nre
      dg3    = ag1
      if (ncmin .eq. nre) goto 3

C --- Runge-Kutta for next three points -----
      ea = dexp(a)
      q  = 1d0/dsqrt(ea)
      if (ifac .eq. 1) q = dsqrt(ea)
      do  1  i = 1, 3
        irp1 = ir
        ir   = ir+ifac
        rpb  = rpb*q
        drdi = rpb*a
        r    = rpb - b
        gg   = g(irp1,1)-.5d0*ag1
        ff   = g(irp1,2)-.5d0*af1
        vb   = (3d0*v(irp1) + 6d0*v(ir) - v(ir-1))*.125d0
        phi  = (e + zz/r - vb)*drdi/c
        u    = drdi*c + phi
        x    = -drdi/r
        y    = -fllp1*x*x/u + phi
        ag2  = u*ff - x*gg
        af2  = x*ff - y*gg
        gg   = g(irp1,1)-.5d0*ag2
        ff   = g(irp1,2)-.5d0*af2
        ag3  = u*ff - x*gg
        af3  = x*ff - y*gg

        rpb  = rpb*q
        drdi = a*rpb
        r    = rpb - b
        phi  = (e + zz/r - v(ir))*drdi/c
        u    = drdi*c + phi
        x    = -drdi/r
        y    = -fllp1*x*x/u + phi
        gg   = g(irp1,1) - ag3
        ff   = g(irp1,2) - af3
        g(ir,1) = g(irp1,1) - (ag1 + 2d0*(ag2 + ag3) + u*ff - x*gg)/6d0
        g(ir,2) = g(irp1,2) - (af1 + 2d0*(af2 + af3) + x*ff - y*gg)/6d0
        if (g(ir,1)*g(irp1,1) .lt. 0d0) nn = nn + 1
        ag1  = u*g(ir,2) - x*g(ir,1)
        af1  = x*g(ir,2) - y*g(ir,1)
        if (ir .eq. ncmin) goto 3
        d(1,i) = ag1
        d(2,i) = af1
    1 continue

C --- All remaining points -----
      q = 1d0/ea
      if (ifac .eq. 1) q = ea
      dg1 = d(1,1)
      dg2 = d(1,2)
      dg3 = d(1,3)
      df1 = d(2,1)
      df2 = d(2,2)
      df3 = d(2,3)
      i2 = nre-4
      i1 = 0
      if (ifac .eq. 1) then
        i2 = nre+4
        i1 = nr
      endif
      do  2  i = i2, i1, ifac
        ir = i
        irp1  = ir-ifac
C        rpb   = rpb*q
C        drdi  = a*rpb
C        r     = rpb - b
C        print *, rofi(ir)-r
C        print *, drdi - a*(r + b)
        r    = rofi(ir)
        drdi = a*(r + b)
        phi  = (e + zz/r - v(ir))*drdi/c
        u    = drdi*c + phi
        x    = -drdi/r
        y    = -fllp1*x*x/u + phi
        det  = r83sq - x*x + u*y
        b1   = g(irp1,1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b2   = g(irp1,2)*h83 + r1*df1 + r2*df2 + r3*df3
        g(ir,1) = (b1*(h83-x) + b2*u)/det
        g(ir,2) = (b2*(h83+x) - b1*y)/det
        if (g(ir,1)*g(irp1,1) .lt. 0d0) nn = nn+1
        dg1  = dg2
        df1  = df2
        dg2  = dg3
        df2  = df3
        dg3  = u*g(ir,2) - x*g(ir,1)
        df3  = x*g(ir,2) - y*g(ir,1)
        if (ifac .eq. -1 .and. mod(ir,2) .ne. 0 .and.
     .     (ir.le.ncmin .or. g(ir,1)*dg3.ge.0d0)) goto 3
    2 continue

C --- Integration done, clean up ---
    3 nc  = ir
      val = g(nc,1)
      drdi= a*(rofi(nc) + b)
      slo = dg3/drdi

C      do  ir  = min(nre+10*ifac,nre),max(nre+10*ifac,nre)
C        write(*,*) rofi(ir), g(ir,1)
C      enddo

      end
      subroutine setcc(lrel)
C- Set speed of light for radial S-eqn in /cc/ common block
C     implicit none
      integer lrel
      double precision c
      common /cc/ c

      if (lrel .ne. 0) then
C       should be:
C       c = 274.072d0
        c = 274.074d0
      else
        c = 1d10
      endif
      end


