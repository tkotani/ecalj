module m_rseq
  use m_lmfinit,only: c=>cc
  public rseq,rsq1
  private
contains
  subroutine rseq(eb1,eb2,e,tol,z,l,nod,val,slo,v,g,q,a,b,rofi,nr,nre) !- Solves radial wave equation for given BCs and number of nodes
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   eb1,eb2 lower and upper bounds for the eigenvalue
    !i   tol     tolerance: maximum relative error in energy;
    !i           absolute error if |energy| < 1
    !i   z       nuclear charge
    !i   a,b     mesh points given by rofi(i) = b [e^(a(i-1)) -1]
    !i   l       angular momentum
    !i   nod     number of nodes
    !i   val,slo BC'S for large component u(r)=g(r,1) with psi=(u/r)*ylm
    !i           val is u(rmax), i.e. (rmax * radial w.f.)
    !i           slo is radial derivative of g at rmax (rmax * radial w.f.)'
    !i   v       spherical potential (true potential excluding nuclear part)
    !i   rofi,nr radial mesh points, and number
    !o Outputs:
    !o   e       eigenvalue
    !o   g       Wave function times r normalized so that int (g*g) dr = 1
    !o           Large component in g(ir,1); small component = g(ir,2)
    !o   q       integral of unnormalized product, that is of
    !o           int (g*g) dr if g were normalized to input val,slo
    !o   nre     index to rofi at which outward and inward solutions join
    !r Remarks:
    !r   Scalar relativistic version
    !r   Output wavefunction normalized to 1 ... so g(nr) .ne. val*wsr
    !r   Note: if r = b(exp(a*z)-1) then Jacobian dr=a(r+b) dz
    !r
    !u Updates
    ! takao       open(7838,file='RSEQ_ERROR') is written when rseq is not converged.
    !u   08 Feb 01 if rseq fails to converge, but the number of nodes is
    !u             correct, rseq returns with a warning instead of aborting
    ! ----------------------------------------------------------------
    implicit none
    integer :: l,nod,nr,nre
    double precision :: eb1,eb2,e,tol,z,val,slo,v(nr),g(nr,2),q,rofi(nr), a,b
    integer :: ipr,k,k2,kc,nctp,nctp0,nit,nitmax,nod1,nod2,node,irsq
    double precision :: de,e1,e2,fac,fllp1,r,ratio=1d99,re,rhok=1d99,slo1,slo2, slop,tmcr,val1,val2,valu,wgt=1d99
    call getpr(ipr)
    nitmax = 80
    if (ipr >= 110) write(*,815) z,l,nod,val,slo,eb1,eb2
815 format(' RSEQ:  Z=',f5.1,'  l=',i1,'  nod=',i1,'  bc=',2f7.3,'  e1,e2=',2f8.2)
    e1 = eb1
    e2 = eb2
    call fctp0(l,nr,rofi,v,z,nctp0)
    if (ipr >= 120)write(*,*)' nit l node nod nre kr','        e1            e              e2            de'
    ! --- Start iterations to find energy ---
    do  1  nit = 1, nitmax
       if (e <= e1 .OR. e >= e2) e = (e1 + e2)/2
       !       call fctp(e,nctp,nctp0,xrim,xmin,nsave,l,rofi,v,z,nr,a,b)
       call fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)
       re = 15d0*rofi(nctp)
       nre = int(dlog(re/b + 1d0)/a + 1d0)
       nre = (nre/2)*2 + 1
       nre = max0(35,min0(nre,nr))
       valu = val
       slop = slo
       if (nre < nr) valu = 1d-5
       if (nre < nr) slop = -1d-5
       k2 = min(nre,30)
       if (nod == 0) k2 = nre/3
       if (valu*slop > 0d0 .AND. nod == 0) k2 = nre - 10
       !       Integrate the scalar relativistic eqn inward from nre to kc
       call rsq2(e,l,z,v,nre,k2,valu,slop,g,val2,slo2,nod2,kc, &
            a,b,rofi,nr)
       !       Integrate the scalar relativistic eqn outward from origin to kc
       call rsq1(0,e,l,z,v,kc,g,val1,slo1,nod1,a,b,rofi,nr)
       node = nod1 + nod2
       if (node /= nod) then
          if (ipr >= 120 .OR. ((nit >= nitmax-5) .AND. ipr>20) ) &
               write(*,101) nit,l,node,nod,nre,kc,e1,e,e2
101       format(2i3,4i4,3f15.7,1p,d13.4)
          if (node > nod) e2 = e
          if (node < nod) e1 = e
          e = (e1 + e2)/2
          goto 1
       endif
1011   continue
       !   ... Calculate q = norm of wave function from trapezoidal rule and
       !       de = estimated correction to eigenvalue
       ratio = val2/val1
       q = 0d0
       do  5  k = 2, kc
          q = q + (rofi(k)+b)*g(k,1)*g(k,1)
5      enddo
       q = q*ratio**2
       do  6  k = kc+1, nre
          q = q + (rofi(k)+b)*g(k,1)*g(k,1)
6      enddo
       q = a*(q - (rofi(nre)+b)*g(nre,1)*g(nre,1)/2)
       de = -val2*(slo2 - ratio*slo1)/q
       if (de > 0d0) e1 = e
       if (de < 0d0) e2 = e
       e = e + de
       !       Exit loop when de meets tolerance; eval found
       if (dabs(de/dmax1(dabs(e),1d0)) < tol) goto 2
1   enddo
    ! --- Search for eigenvalue failed ---
    nit = nitmax+1
    !     Fatal if node mismatch
    if (nod /= node) goto 99
    ! --- Normalize g ---
2   continue
    fllp1 = l*(l+1)
    e = e - de
    do  8  k = 1, kc
       g(k,1) = g(k,1)*ratio
       g(k,2) = g(k,2)*ratio
8   enddo
    q = 0d0
    do  10  k = 2, nre
       r = rofi(k)
       wgt = (mod(k+1,2) + 1)*(r + b)
       tmcr = (c - (v(k) - 2d0*z/r - e)/c)*r
       rhok = g(k,1)*g(k,1)*(1d0 + fllp1/(tmcr*tmcr)) + g(k,2)*g(k,2)
       q = q + wgt*rhok
10  enddo
    q = (q - wgt*rhok/2)*a*2d0/3d0
    fac = 1d0/dsqrt(q)
    do  11  k = 1, nre
       g(k,1) = g(k,1)*fac
       g(k,2) = g(k,2)*fac
11  enddo
    do  12  k = nre+1, nr
       g(k,1) = 0d0
       g(k,2) = 0d0
12  enddo
    ! --- Possible warning or error exit ---
99  continue
    if (ipr >= 110 .OR. nit >= nitmax) &
         write(*,701) e,q,nr,nre,kc,de
701 format(' e=',f13.5,'  q=',f10.5,'   nr/nre/kc=',3i4, &
         '   de=',1p,d10.2)
    if (nit > nitmax .AND. ipr >20) then
       if (nod /= node) then
          write(*,814) nitmax,l,nod,node,e
814       format(' RSEQ : nit gt',i3,' and bad nodes for l=',i1, &
               '.  Sought',i2,' but found',i2,'.  e=',1pd11.4)
          call rx('RSEQ: bad nodes')
       else
          write(*,816) l,nod,abs(de),e
816       format(' RSEQ (warning) eval for l=',i1,' node=',i1, &
               ' did not converge: de=',1pd9.2,' e=',1pd11.4)
          block
            character(256):: aaaa
            write(aaaa,816) l,nod,abs(de),e
            call rx('rseq_error'//trim(aaaa))
          endblock
          !          open(newunit=irsq,file='RSEQ_ERROR')
          !          write(irsq,816) l,nod,abs(de),e
          !          close(irsq)
       endif
    endif
  end subroutine rseq
  subroutine rsq1(i0,e,l,z,v,kr,g,val,slo,nn,a,b,rofi,nr)
    !- Integrate the scalar relativistic eqn outward to rofi(kr)
    ! ----------------------------------------------------------------
    !i Inputs:
    !i   i0      *If i0<5, integrate from origin to rofi(kr)
    !i           *Else, integrate from rofi(i0) to rofi(kr).
    !i            NB: in this case, g(i0-4..ir) must be input!
    !i   e        Energy
    !i   l,z,v    angular momentum, nuc. charge, potential (v_true)
    !i   a,b      mesh points given by rofi(i) = b [e^(a(i-1)) -1]
    !i   rofi,nr  radial mesh points, and number
    !i   kr       integration from origin to rofi(kr)
    !o Outputs:
    !o   g        radial wave function times r (but see i0),
    !o            normalized to unity.
    !o   val,slo  value and slope of g at rofi(kr)
    !o   nn       number of nodes encountered
    !r Remarks:
    !r   Boundary condition does not fix value of wave function near the
    !r   origin; integration can be scaled by an arbritrary factor
    !u Updates
    !u   21 Jun 04 Integration can start from point other than origin
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters:
    integer :: i0,l,kr,nn,nr
    double precision :: e,z,v(nr),g(nr,2),val,slo,a,b,rofi(nr)
    ! Local parameters:
    integer :: ir,i1
    double precision :: d(2,3),zz,fllp1,r83sq,r1,r2,r3,h83,g0,s,sf,aa, &
         f0,dg1,dg2,dg3,df1,df2,df3,r,drdi,phi,u,x,y,det,b1,b2
    equivalence (dg1,d(1,1)),(dg2,d(1,2)),(dg3,d(1,3)), (df1,d(2,1)),(df2,d(2,2)),(df3,d(2,3))
    nn = 0
    zz = z+z
    fllp1 = l*(l+1)
    r83sq = 64d0/9d0
    r1 = 1d0/9d0
    r2 = -5d0*r1
    r3 = 19d0*r1
    h83 = 8d0/3d0
    ! --- Approximate g,f by leading term near zero ----
    if (i0 < 5) then
       g0 = 1
       if (z < 0.9d0) then
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
2      enddo
    endif
    ! --- Setup to integrate over rest of points ------
    if (i0 < 5) then
       i1 = 5
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
          dg1   = dg2
          df1   = df2
          dg2   = dg3
          df2   = df3
          dg3   = u*g(ir,2) - x*g(ir,1)
          df3   = x*g(ir,2) - y*g(ir,1)
3      enddo
    endif

    ! --- Integrate over rest of points ------
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
       if (g(ir,1)*g(ir-1,1) < 0d0) nn = nn+1
       dg1   = dg2
       df1   = df2
       dg2   = dg3
       df2   = df3
       dg3   = u*g(ir,2) - x*g(ir,1)
       df3   = x*g(ir,2) - y*g(ir,1)
4   enddo
    val = g(kr,1)
    slo = dg3/(a*(rofi(kr) + b))
  end subroutine rsq1
  subroutine rsq2(e,l,z,v,nre,ncmin,val1,slo1,g,val,slo,nn,nc, a,b,rofi,nr)
    !- Integrate the scalar relativistic eqn inward from nre to nc
    ! ----------------------------------------------------------------
    !i Inputs:
    !i   e     :energy
    !i   l     :angular momentum
    !i   z     :nuclear charge
    !i   v     :spherical potential = v_true (excluding nuclear 2*Z/r)
    !i   nre   :rofi(nre) = starting radius from which rsq2 integrates
    !i   ncmin: rsq2 integrates to cutoff nc, of which ncmin is lower bound
    !i   val1  :value of large component of g at nre
    !i   slo1  :slope of g at nre
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :                 -//-
    !i   rofi  :radial mesh points
    !i   nr    :leading dimension of g
    !o Outputs:
    !o   g     :radial wave function times r at points nc..nre
    !o   nc    :rofi(nc) = radius to which rsq2 integrates
    !o         :nc is the larger of ncmin and the position of the
    !o         :first maximum encountered
    !o   val   :value of g(nc) (large component)
    !o   slo   :slope of g(nc)
    !o   nn    :number of nodes found between nc and nre
    !u Updates
    !u   21 Jun 04 Attempt to extend integration also to outward.
    !u             Doesn't work.  Use rsq1
    !r Remarks:
    !r   Integrates inward from nre to nc
    !r   Cutoff nc is chosen at first maximum, but nc .ge. ncmin
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters:
    integer :: l,nre,ncmin,nn,nc,nr
    double precision :: e,z,v(nr),val1,slo1,g(nr,2),val,slo,a,b,rofi(nr)
    ! Local parameters:
    integer :: i,ir,irp1,i2,i1,ifac
    double precision :: d(2,3),zz,fllp1,r83sq,r1,r2,r3,h83,ea,rpb, &
         q,ag1,ag2,ag3,af1,af2,af3,gg,ff,vb,drdi,r,phi,u,x,y, &
         dg1,dg2,dg3,df1,df2,df3,det,b1,b2
    nn = 0
    zz = z + z
    fllp1 = l*(l + 1)
    r83sq =64d0/9d0
    r1    = 1d0/9d0
    r2    =-5d0/9d0
    r3    =19d0/9d0
    h83   =-8d0/3d0
    ifac = -1
    !     if (ncmin .gt. nre) ifac = 1
    if (ncmin > nre) &
         call rx('rsq2 not implemented for outward integration')

    ! --- First point ------
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
    if (ncmin == nre) goto 3

    ! --- Runge-Kutta for next three points -----
    ea = dexp(a)
    q  = 1d0/dsqrt(ea)
    if (ifac == 1) q = dsqrt(ea)
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
       if (g(ir,1)*g(irp1,1) < 0d0) nn = nn + 1
       ag1  = u*g(ir,2) - x*g(ir,1)
       af1  = x*g(ir,2) - y*g(ir,1)
       if (ir == ncmin) goto 3
       d(1,i) = ag1
       d(2,i) = af1
1   enddo
    ! --- All remaining points -----
    q = 1d0/ea
    if (ifac == 1) q = ea
    dg1 = d(1,1)
    dg2 = d(1,2)
    dg3 = d(1,3)
    df1 = d(2,1)
    df2 = d(2,2)
    df3 = d(2,3)
    i2 = nre-4
    i1 = 0
    if (ifac == 1) then
       i2 = nre+4
       i1 = nr
    endif
    do  2  i = i2, i1, ifac
       ir = i
       irp1  = ir-ifac
       !        rpb   = rpb*q
       !        drdi  = a*rpb
       !        r     = rpb - b
       !        print *, rofi(ir)-r
       !        print *, drdi - a*(r + b)
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
       if (g(ir,1)*g(irp1,1) < 0d0) nn = nn+1
       dg1  = dg2
       df1  = df2
       dg2  = dg3
       df2  = df3
       dg3  = u*g(ir,2) - x*g(ir,1)
       df3  = x*g(ir,2) - y*g(ir,1)
       if (ifac == -1 .AND. mod(ir,2) /= 0 .AND. &
            (ir <= ncmin .OR. g(ir,1)*dg3 >= 0d0)) goto 3
2   enddo
    ! --- Integration done, clean up ---
3   nc  = ir
    val = g(nc,1)
    drdi= a*(rofi(nc) + b)
    slo = dg3/drdi
    !      do  ir  = min(nre+10*ifac,nre),max(nre+10*ifac,nre)
    !        write(*,*) rofi(ir), g(ir,1)
    !      enddo
  end subroutine rsq2
  subroutine fctp0(l,nr,rofi,v,z,nctp0)
    !- Initialize things for FCTP, which finds classical turning point
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   l     :angular momentum
    !i   rofi  :radial mesh points
    !i   v     :spherical potential
    !i   z     :nuclear charge
    !i   nr    :number of mesh points
    !o Outputs:
    !o   nctp0 :minimum of effective potential
    !o          or nr if v(nr) > vmin + 3 Ry
    ! ----------------------------------------------------------------------
    implicit none
    integer :: l,nr,nctp0
    double precision :: z,v(nr),rofi(nr)
    integer :: ir,irmin
    double precision :: veff,zz,fllp1,vi,vim1
    parameter (irmin=11)
    veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)
    zz = z+z
    fllp1 = l*(l+1)
    ir = irmin
    vi = veff(ir)
    vim1 = veff(ir-1)
    do while (vi <= vim1 .AND. ir < nr) 
       ir = ir+1
       vim1 = vi
       vi = veff(ir)
    enddo
    nctp0 = ir-1
    if (veff(nctp0) >= veff(nr)-3.d0) nctp0 = nr
  end subroutine fctp0
  subroutine fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)
    !- Finds classical turning point
    !  ---------------------------------------------------------------------
    !i Inputs:
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :                 -//-
    !i   e     :energy
    !i   l     :angular momentum
    !i   nctp0 :minimum of effective potential (fctp0)
    !i   nr    :number of mesh points
    !i   rofi  :radial mesh points
    !i   v     :spherical potential
    !i   z     :nuclear charge
    !i Inputs/Outputs:
    ! o  nctp  :rofi(nctp) is classical turning point
    !u Updates
    !u   13 Jun 00 Added safety checks to guard against jumps in potential
    !  ---------------------------------------------------------------------
    implicit none
    integer :: nctp,nctp0,l,nr
    double precision :: e,rofi(nr),v(nr),z,a,b
    integer :: ir,irep,n1,n2,nlast,ntry
    double precision :: r,veff,fllp1,fofr,dfdr,rtry,zz
    veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)
    zz=z+z
    fllp1 = l*(l+1)
    if (nctp0 == nr .OR. e > veff(nr)) then
       nctp = nr
    elseif (e < veff(nctp0)) then
       nctp = 2
    else
       n1 = nctp0
       n2 = nr-1
       nlast = -10
       do  10  irep = 1, 20
          if (nctp > n2 .OR. nctp < n1) nctp = (n1 + n2)/2
          r = rofi(nctp)
          fofr =  veff(nctp)-e
          dfdr =-(fllp1+fllp1)/r/r/r + zz/r/r + &
               (v(nctp+1) - v(nctp-1))/(2.d0*a*(r+b))
          rtry = dmax1(r-fofr/dfdr,rofi(2))
          ntry = dlog(rtry/b+1.d0)/a + 1.5d0
          !         If there was a large change, check for safety
          if (nlast == nctp .AND. iabs(ntry-nctp) > 10) then
             if (fofr > 0) then
                do  ntry = nctp-1, nctp0+1, -1
                   fofr =  veff(ntry)-e
                   if (fofr < 0) goto 20
                   nctp = ntry
                   n2 = ntry
                enddo
20              continue
             else
                do  ntry = nctp+1, nr-1
                   fofr =  veff(ntry)-e
                   if (fofr > 0) goto 21
                   nctp = ntry
                   n1 = ntry
                enddo
21              continue
             endif
          endif
          !         Exit point
          if (nlast == nctp) then
             if (nctp == nctp0+1) nctp = 2
             return
          endif
          if (fofr > 0.d0)  then
             n2 = nctp
          else
             n1 = nctp
          endif
          nlast = nctp
          nctp = min0(ntry,nr-1)
10     enddo
       if (nctp == nctp0+1) nctp = 2
    endif
  end subroutine fctp
endmodule m_rseq

