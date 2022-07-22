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

    logical:: l_dummy_isanrg,isanrg

    ! ... Setup
    if (mode == 0) return
    ! ino isanrg is logical function,       call isanrg(mode,0,8,'vxtrap:','mode',.true.)
    l_dummy_isanrg=isanrg(mode,0,8,'vxtrap:','mode',.true.)
    lmode0 = mod(mode,2) .ne. 0
    lmode1 = mod(mode/2,2) .ne. 0
    lmode2 = mod(mode/4,2) .ne. 0
    if (lmode2 .AND. (lmode0 .OR. lmode1)) &
         call rx('vxtrap: incompatible switches')

    ! --- Check for local orbital types ---
    if (lmode0) then
       lpz = 0
       if (dasum(lmxa+1,pnz,1) /= 0) lpz = 1
       if (lpz == 1) then
          do  2  l = 0, lmxa
             k = l+1
             if (pnz(k) >= 10) lpz = 2
2         enddo
       endif
    endif

    ! --- Determine large radius for local orbitals ---
    !     r(nrx)/rmax = (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
    if (lmode0 .AND. lpz > 1) then
       call radext(1,nr,nrx,2d0,a,rmax,nrbig,rbig,xx,xx)
       !        rbig = rmax * (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
       !C       If rbig>2*rmax, estimate from exp((nrbig-nr)a) = 2
       !        if (rbig .gt. 2*rmax) then
       !          idn = dlog(2d0)/a
       !          if (mod(idn,2) .eq. 1) idn = idn-1
       !          nrbig = min(nr+idn,nrx)
       !          rbig = rmax * (dexp(a*nrbig-a)-1d0)/(dexp(a*nr-a)-1d0)
       !        endif
    elseif(lmode0) then
       nrbig = nr
       rbig = rmax
    elseif (lmode1) then
       call rx('vxtrap not ready for 2s digit mode')
    endif

    ! --- Points and weights on extended mesh ---
    if (lmode0 .OR. lmode1) then
       call radext(10,nr,nrx,2d0,a,rmax,nrbig,rbig,rofi,rwgt)
       !        call radmsh(rbig,a,nrbig,rofi)
       !        call radwgt(rbig,a,nrbig,rwgt)
       !        if (nr .lt. nrbig) rwgt(nr) = rwgt(nr)/2
       return
    endif

    ! --- Potential on extended mesh ---
    if ( .NOT. lmode2) return
    ! ino isanrg is logical function,       call isanrg(nrbig,nr,nrx,'vxtrap:','nrbig',.true.)
    l_dummy_isanrg=isanrg(nrbig,nr,nrx,'vxtrap:','nrbig',.true.)
    if (rbig < rmax) call rx('vxtrap: illegal value for rbig')

    do  isp = 1, nsp
       ! ... Extrapolate potential, if nrbig > nr
       if (nr < nrbig) then

          !       Make r*(true V - vmtz) to facilitate 2nd derivative
          do  ir = 2, nr
             vbig(ir,isp) = (v(ir,isp) - 2*z/rofi(ir) - vmtz)*rofi(ir)
             !         vbig(ir,isp) = exp(-2d0*rofi(ir))
          enddo
          !       call prrmsh('v*r',rofi,vbig,nr,nr,1)

          !       First derivative in gbig; second derivative in gzbig
          k = nr-2*nptirp+1
          if (k < 2) call rxi('vxtrap: nr too small, nr=',nr)
          call poldvm(rofi(k),vbig(k,isp),2*nptirp,nptirp,.false.,1d-8,j, &
               gbig(k))
          call poldvm(rofi(k),gbig(k),2*nptirp,nptirp,.false.,1d-8,j, &
               gzbig(k))
          !       call prrmsh('v 1st d',rofi(k),gbig(k),2*nptirp,2*nptirp,1)
          !       call prrmsh('v 2nd d',rofi(k),gzbig(k),2*nptirp,2*nptirp,1)

          !       Use r*V' = (rV)' - V
          phi  = vbig(nr,isp)/rofi(nr)
          dphi = (gbig(nr) - phi)/rofi(nr)
          tphi = gzbig(nr)/vbig(nr,isp)
          call info5(40/10,0,0,' vxtrap: extrapolate V using '// &
               'vmtz=%;3,3d  v-vmtz=%;3,3d  v''=%;3,3d  (nabla v)/v=%;3,3d', &
               vmtz,phi,dphi,tphi,0)
          call rwftai(3,rmax,a,nr,nrbig,rofi,phi,dphi,tphi,0,ehv,rs3, &
               vbig(1,isp))

          !       Restore vbig to estimate for true V + 2Z/r
          vbig(1,isp) = v(1,isp)
          do  ir = 2, nrbig
             vbig(ir,isp) = vbig(ir,isp)/rofi(ir) + 2*z/rofi(ir) + vmtz
          enddo

          !       For pictures ...
          !        do  ir = 2, nrbig
          !          vbig(ir,isp) = vbig(ir,isp) - 2*z/rofi(ir)
          !        enddo
          !        call prrmsh('vbig',rofi,vbig,nrmax,nrbig,nsp)
          !        stop

          ! ... Just copy potential if nr = nrbig
       else
          call dcopy(nr,v(1,isp),1,vbig(1,isp),1)
       endif
    enddo
  end subroutine vxtrap
  subroutine rwftai(mode,rmt,a,nrmt,nrbig,ribig,phi,dphi,tphi,l, ehl,rsml,g)
    use m_hansr,only :hansr,hansmd
    !- Extend radial wave function outside MT boundary
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :controls approximation in which w.f. is extended
    !i         :1s digit
    !i         :0 Make extended radial mesh (to smaller of 2*rmax or r(nrx))
    !i         :  NO LONGER IMPLEMENTED; use radext
    !i         :1 Use info about sm. Hankel tail to estimate radius to
    !i         :  to which w.f. need be extended to incorporate full charge
    !i         :  (not implemented)
    !i         :2 Compute radial wave function on extended mesh using
    !i         :  rsml,ehl for tail, scaling envelope to match gz.
    !i         :3 Compute radial wave function on extended mesh matching
    !i         :  phi,dphi to h,hdot.  Use ehl for sm. Hankel energy
    !i         :  and find rsm that best matches kinetic energy tphi
    !i         :4 similar to mode 3, except if no KE match, also attempt
    !i         :  to vary ehl.  NB: not implemented.
    !i         :  Used to extrapolate potential.  In this case, g=r*V.
    !i         :5 Similar to mode 2, except scale gz so that value
    !i         :  matches envelope h(rsm,eh), rather than scaling envelope
    !i   rmt   :augmentation radius, in a.u., by species
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nrmt  :number of radial mesh points between origin and rmt
    !i   phi   :value of normalized wave function at rmt, i.e. g/rmt
    !i         :(used for mode = 3)
    !i   dphi  :slope of normalized w.f. at rmt, i.e. (d(g/r)/dr)_rmt
    !i         :(used for mode = 3)
    !i   tphi  :kinetic energy of w.f. at rmt
    !i         :(used for mode = 3)
    !i   l     :quantum number for this wave function
    !i   ehl   :energy of smoothed Hankel tail for local orbital
    !i   rsml  :smoothing radius for smoothed Hankel tail of local orbital
    ! o Inputs/Outputs
    ! o  nrbig  :number of radial mesh points on extended mesh
    ! o         :mode  0 : output
    ! o         :mode >0 : input
    ! o  ribig  :points for extended radial mesh
    ! o         :mode  0 : output
    ! o         :mode >0 : input
    !o Outputs
    !o   g      :radial wave function extended to points nrmt..nrbig
    !o          :(done in modes >0)
    !l Local variables
    !l   rbig   :rmax of extended mesh
    !l   nxtn   :number of points between rmt and rbig
    !r Remarks
    !r
    !u Updates
    !u   27 Jun 04 First created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: mode,nrmt,nrbig,l
    double precision :: a,rmt,phi,dphi,tphi,ribig(*),g(nrbig,2),ehl,rsml
    ! ... Local parameters
    integer :: nrx,idn,IPRTDB,nxtn,ir,info,mode0
    parameter (nrx=1501,IPRTDB=10)
    double precision :: rbig,rwgt(nrx),xi(nrx,0:l),r2(nrx) !,wk(2*nrx)
    double precision :: fac1,fac2,rsm,ekin,eh,alfa,beta
    double precision :: hs(0:l),dhs(0:l),ddhs(0:l)
    double precision :: hsp(0:l),dhsp(0:l),ddhsp(0:l)
    !     double precision dot3,sum1,sum2

    mode0 = mod(mode,10)
    !     mode1 = mod(mode/10,10)

    ! --- Get rbig which corresponds to smaller of r(nrx) and 2*rmt ---
    if (mode0 == 0) then
       call rx('rwftai: mode not implemented')

       !     elseif (mode0 .eq. 1) then
    elseif (mode0 == 2 .OR. mode0 == 5) then
       rbig = ribig(nrbig)
       !       Uncomment following 2 lines to replace entire function with envelope
       !       nrmt = 1
       !       g = 0
       nxtn = nrbig-nrmt+1
       call radwgt(rbig,a,nrbig,rwgt)
       rwgt(nrmt) = rwgt(nrmt)/2
       do  ir = 1, nxtn
          r2(ir) = ribig(ir-1+nrmt)**2
       enddo
       !       Hankel function scaled to match g at nrmt
       !        call hansr(rsml,0,l,1,l,ehl,r2,nrx,nxtn,idn,wk,11,xi)
       call hansr(rsml,0,l,1,[l],[ehl],[r2],nrx,nxtn,[idn],11,xi)
       !       fac1 = phi/xi(1,l)
       fac1 = g(nrmt,1)/rmt/xi(1,l)
       fac2 = g(nrmt,2)/rmt/xi(1,l)
       if (mode0 == 5) then
          call dscal(nrmt,1/fac1,g(1,1),1)
          call dscal(nrmt,1/fac1,g(1,2),1)
          fac2 = fac2 / fac1
          fac1 = 1
       endif
       do  ir = 1, nxtn
          g(ir-1+nrmt,1) = xi(ir,l)*ribig(ir-1+nrmt) * fac1
          g(ir-1+nrmt,2) = xi(ir,l)*ribig(ir-1+nrmt) * fac2
       enddo
       !C       Integral of gg from 0 to rmt, neglecting small component
       !        sum1 = dot3(nrmt,rwgt,g,g)
       !C       Integral of gg from rmt to rbig, neglecting small component
       !        sum2 = dot3(nxtn,rwgt(nrmt),g(nrmt,1),g(nrmt,1))
       !       call info2(IPRTDB,0,0,' rwftai: <gg>_MT=%g  <gg>_big=%g',
       !    .    sum1,sum2)
       !       call prrmsh('gzbig',ribig,g,nrbig,nrbig,2)


    elseif (mode0 == 3) then
       eh = min(ehl,-.05d0)
       call pshpr(110)
       call mtchre(104,l,rsml,rmt,0d0,0d0,rmt,rmt,phi,dphi, &
            tphi,dphi,rsm,eh,ekin,info)
       call poppr
       call mtchae(1,rsm,eh,l,rmt,phi,dphi,0d0,0d0,alfa,beta)
       do  ir = nrmt+1, nrbig
          call hansmd(12,ribig(ir),eh,rsm,l,hs,dhs,ddhs,hsp,dhsp,ddhsp)
          g(ir,1) = (alfa*hs(l) + beta*hsp(l)) * ribig(ir)
       enddo
    else
       call rxi('rwftai: mode not implemented:',mode)
    endif
  end subroutine rwftai
end module m_vxtrap
