module m_vxtrap !- Extrapolate potential, radial wave function outside MT boundary
  public vxtrap,rwftai
  private
contains
  subroutine vxtrap(mode,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrmax, lpz,nrbig,rofi,rwgt,vbig)
    !- Extrapolate potential outside MT boundary
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :tells vxtrap what to do.
    !i         :0 do nothing, just return
    !i         :1 determine whether mesh should be extended based on need
    !i         :  for extension by one or more local orbitals.
    !i         :  Find extended radius suitable for extending local orbitals.
    !i         :  Make points and weights for extended mesh (whether or not
    !i         :  mesh need be extended, rofi, rgwt are returned).
    !i         :2 intended for 3-fold OKA-type extrapolation. Not used now
    !i         :  Switches 1 and 2 may be taken in combination
    !i         :4 make potential on extended mesh
    !i   z     :(mode=4 only) nuclear charge
    !i   rmax  :augmentation radius, in a.u.
    !i   lmxa  :(mode=1 only) augmentation L-cutoff
    !i   v     :(mode=4 only) spherical potential to be extrapolated
    !i   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
    !i   nr    :number of radial mesh points
    !i   nsp   :(mode=4 only) number of spins for which to make vbig
    !i   pnz   :boundary conditions for second p.q.n (local orbital).
    !i          10s digit controls how local orbital included in hamiltonian
    !i         :The following four parameters are used to extrapolate
    !i         :quantities outside the MT radius.  Extrapolation is used to
    !i         :correct matrix elements for local orbitals.
    !i   rs3   :(mode=4 only) minimum smoothing radius for extrapolation
    !i         :of MT potential
    !i   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
    !i         :The asymptotic form of V-vmtz is taken to be zero.
    !i   nrmax :(mode=4 only) leading dimenson of vbig
    ! o Inputs/Outputs
    ! o  nrbig :number of points on extended mesh.
    ! o        :nrbig is output if mode=1 or mode=2; nrbig is input for mode=4
    ! o        :In the latter case, nrbig must be consistent with the mesh
    ! o        :points specified by (a,nr,rmax) and also rbig.
    ! o  rofi  :(Output, mode=1 or 2; Input, mode 4)
    ! o        :radial mesh points: rofi(1..nrbig) will be generated
    ! o        :rofi(nrbig) is rmax for extended mesh
    !o Outputs
    !o   lpz   :(mode=1 only)
    !o         : 0 if no local orbitals
    !o         : 1 if no local orbitals requiring extended mesh
    !o         : 2 if no at least one loc. orbital requiring extended mesh
    !o   rwgt  :(mode=1 or 2 only)
    !o         :radial mesh weights: rofi(1..nrbig) will be generated
    !o         :NB: rwgt is actually designed for two integration radii:
    !o         :int(0,rmax) = I(1..nr) and int(rmax,rbig) = I(nr..nrbig).
    !o         :Integral int(1..nrbig) must be done in two steps, by summing
    !o         :I(1..nr) and I(nr..nrbig)
    !o   vbig  :(mode=4 only) potential extrapolated from v(*,isp)
    !o         :NB: vbig is returned for only one spin at present
    !b Bugs
    !b   This routine is not sufficiently robust.  Two problems to solve:
    !b     vmtz exceed v(rmax), causing
    !b     fitting procedure may not work; see rwftai
    !b     Test case, ctrl.srtio3, Sr atom.
    !b     Possible solution:  enhance rwftai so that it can vary sm Hankel
    !b     energy in order to force match of K.E.
    !r Remarks
    !u Updates
    !u   12 Jul 04 First created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: mode,lmxa,nr,nrmax,nsp,lpz,nrbig
    double precision :: z,rmax,a,v(nr,nsp),pnz(1+lmxa),rs3,vmtz, &
         rofi(*),rwgt(*),vbig(nrmax,nsp)
    ! ... Local parameters
    logical :: lmode0,lmode1,lmode2
    integer :: ir,j,k,l,idn,isp
    !     nptirp = number of points to include in fitting potential
    !              for extrapolation to large sphere
    !     ehv    = smooth Hankel energy used to extrapolate V
    !              The deeper ehv, the more rapidly V approaches vmtz.
    !              Anyway, extrapolation is somewhat arbitrary.
    integer :: nrx,nptirp
    double precision :: ehv
    parameter (nrx=1501,nptirp=10,ehv=-1d0)
    double precision :: dasum,tphi,rbig,phi,dphi,xx
    double precision :: gzbig(nrx*2),gbig(nrx*2)
    ! --- Check for local orbital types ---
       lpz = 0
       if (dasum(lmxa+1,pnz,1) /= 0) lpz = 1
       if (lpz == 1) then
          do  2  l = 0, lmxa
             k = l+1
             if (pnz(k) >= 10) lpz = 2
2         enddo
       endif
    ! --- Determine large radius for local orbitals ---
    !     r(nrx)/rmax = (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
    if (lpz > 1) then
       call radext(1,nr,nrx,2d0,a,rmax,nrbig,rbig,xx,xx)
       !        rbig = rmax * (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
       !C       If rbig>2*rmax, estimate from exp((nrbig-nr)a) = 2
       !        if (rbig .gt. 2*rmax) then
       !          idn = dlog(2d0)/a
       !          if (mod(idn,2) .eq. 1) idn = idn-1
       !          nrbig = min(nr+idn,nrx)
       !          rbig = rmax * (dexp(a*nrbig-a)-1d0)/(dexp(a*nr-a)-1d0)
       !        endif
    else
       nrbig = nr
       rbig = rmax
    endif
    ! --- Points and weights on extended mesh ---
    call radext(10,nr,nrx,2d0,a,rmax,nrbig,rbig,rofi,rwgt)
  end subroutine vxtrap
  subroutine rwftai(mode,rmt,a,nrmt,nrbig,ribig,phi,dphi,tphi,l, ehl,rsml,g)
    use m_hansr,only :hansr
    !- Extend radial wave function outside MT boundary
    ! ----------------------------------------------------------------------
    !i Inputs (mode=dummy)
    !i         :Compute radial wave function on extended mesh using
    !i         : rsml,ehl for tail, scale gz so that value
    !i         : matches envelope h(rsm,eh)
    !i   rmt   :augmentation radius, in a.u., by species
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nrmt  :number of radial mesh points between origin and rmt
    !i   l     :quantum number for this wave function
    !i   ehl   :energy of smoothed Hankel tail for local orbital
    !i   rsml  :smoothing radius for smoothed Hankel tail of local orbital
    !i  nrbig  :number of radial mesh points on extended mesh
    !i  ribig  :points for extended radial mesh
    !o Outputs
    !o   g      :radial wave function extended to points nrmt..nrbig
    !l Local variables
    !l   rbig   :rmax of extended mesh
    !l   nxtn   :number of points between rmt and rbig
    !r Remarks
    !r
    !u Updates
    !u   27 Jun 04 First created
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode,nrmt,nrbig,l
    double precision :: a,rmt,phi,dphi,tphi,ribig(*),g(nrbig,2),ehl,rsml
    integer :: nrx,idn,IPRTDB,nxtn,ir,info,mode0
    parameter (nrx=1501,IPRTDB=10)
    double precision :: rbig,rwgt(nrx),xi(nrx,0:l),r2(nrx)
    double precision :: fac1,fac2,rsm,ekin,eh,alfa,beta
    rbig = ribig(nrbig)
    nxtn = nrbig-nrmt+1
    call radwgt(rbig,a,nrbig,rwgt)
    rwgt(nrmt) = rwgt(nrmt)/2
    r2(1:nxtn) = [(ribig(ir-1+nrmt)**2,ir=1,nxtn)]
    call hansr(rsml,0,l,1,[l],[ehl],[r2],nrx,nxtn,[idn],11,xi)
    fac1 = g(nrmt,1)/rmt/xi(1,l) ! Hankel function scaled to match g at nrmt
    fac2 = g(nrmt,2)/rmt/xi(1,l)
    g(1:nrmt,:)=g(1:nrmt,:)/fac1 
    fac2 = fac2 / fac1
    fac1 = 1
    g(nrmt:nrmt+nxtn-1,1) = [(xi(ir,l)*ribig(ir-1+nrmt) * fac1,ir=1,nxtn)]
    g(nrmt:nrmt+nxtn-1,2) = [(xi(ir,l)*ribig(ir-1+nrmt) * fac2,ir=1,nxtn)]
  end subroutine rwftai
end module m_vxtrap
