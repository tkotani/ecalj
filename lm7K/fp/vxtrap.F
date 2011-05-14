      subroutine vxtrap(mode,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrmax,
     .lpz,nrbig,rofi,rwgt,vbig)
C- Extrapolate potential outside MT boundary
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :tells vxtrap what to do.
Ci         :0 do nothing, just return
Ci         :1 determine whether mesh should be extended based on need
Ci         :  for extension by one or more local orbitals.
Ci         :  Find extended radius suitable for extending local orbitals.
Ci         :  Make points and weights for extended mesh (whether or not
Ci         :  mesh need be extended, rofi, rgwt are returned).
Ci         :2 intended for 3-fold OKA-type extrapolation. Not used now
Ci         :  Switches 1 and 2 may be taken in combination
Ci         :4 make potential on extended mesh
Ci   z     :(mode=4 only) nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   lmxa  :(mode=1 only) augmentation L-cutoff
Ci   v     :(mode=4 only) spherical potential to be extrapolated
Ci   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nsp   :(mode=4 only) number of spins for which to make vbig
Ci   pnz   :boundary conditions for second p.q.n (local orbital).
Ci          10s digit controls how local orbital included in hamiltonian
Ci         :The following four parameters are used to extrapolate
Ci         :quantities outside the MT radius.  Extrapolation is used to
Ci         :correct matrix elements for local orbitals.
Ci   rs3   :(mode=4 only) minimum smoothing radius for extrapolation
Ci         :of MT potential
Ci   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
Ci         :The asymptotic form of V-vmtz is taken to be zero.
Ci   nrmax :(mode=4 only) leading dimenson of vbig
Cio Inputs/Outputs
Cio  nrbig :number of points on extended mesh.
Cio        :nrbig is output if mode=1 or mode=2; nrbig is input for mode=4
Cio        :In the latter case, nrbig must be consistent with the mesh
Cio        :points specified by (a,nr,rmax) and also rbig.
Cio  rofi  :(Output, mode=1 or 2; Input, mode 4)
Cio        :radial mesh points: rofi(1..nrbig) will be generated
Cio        :rofi(nrbig) is rmax for extended mesh
Co Outputs
Co   lpz   :(mode=1 only)
Co         : 0 if no local orbitals
Co         : 1 if no local orbitals requiring extended mesh
Co         : 2 if no at least one loc. orbital requiring extended mesh
Co   rwgt  :(mode=1 or 2 only)
Co         :radial mesh weights: rofi(1..nrbig) will be generated
Co         :NB: rwgt is actually designed for two integration radii:
Co         :int(0,rmax) = I(1..nr) and int(rmax,rbig) = I(nr..nrbig).
Co         :Integral int(1..nrbig) must be done in two steps, by summing
Co         :I(1..nr) and I(nr..nrbig)
Co   vbig  :(mode=4 only) potential extrapolated from v(*,isp)
Co         :NB: vbig is returned for only one spin at present
Cb Bugs
Cb   This routine is not sufficiently robust.  Two problems to solve:
Cb     vmtz exceed v(rmax), causing
Cb     fitting procedure may not work; see rwftai
Cb     Test case, ctrl.srtio3, Sr atom.
Cb     Possible solution:  enhance rwftai so that it can vary sm Hankel
Cb     energy in order to force match of K.E.
Cr Remarks
Cu Updates
Cu   12 Jul 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,lmxa,nr,nrmax,nsp,lpz,nrbig
      double precision z,rmax,a,v(nr,nsp),pnz(1+lmxa),rs3,vmtz,
     .rofi(*),rwgt(*),vbig(nrmax,nsp)
C ... Local parameters
      logical lmode0,lmode1,lmode2
      integer ir,j,k,l,idn,isp
C     nptirp = number of points to include in fitting potential
C              for extrapolation to large sphere
C     ehv    = smooth Hankel energy used to extrapolate V
C              The deeper ehv, the more rapidly V approaches vmtz.
C              Anyway, extrapolation is somewhat arbitrary.
      integer nrx,nptirp
      double precision ehv
      parameter (nrx=1501,nptirp=10,ehv=-1d0)
      double precision dasum,tphi,rbig,phi,dphi,xx
      double precision gzbig(nrx*2),gbig(nrx*2)

C ... Setup
      if (mode .eq. 0) return
Ckino isanrg is logical function,       call isanrg(mode,0,8,'vxtrap:','mode',.true.)
      l_dummy_isanrg=isanrg(mode,0,8,'vxtrap:','mode',.true.)
      lmode0 = mod(mode,2) .ne. 0
      lmode1 = mod(mode/2,2) .ne. 0
      lmode2 = mod(mode/4,2) .ne. 0
      if (lmode2 .and. (lmode0 .or. lmode1))
     .call rx('vxtrap: incompatible switches')

C --- Check for local orbital types ---
      if (lmode0) then
        lpz = 0
        if (dasum(lmxa+1,pnz,1) .ne. 0) lpz = 1
        if (lpz .eq. 1) then
          do  2  l = 0, lmxa
            k = l+1
            if (pnz(k) .ge. 10) lpz = 2
    2     continue
        endif
      endif

C --- Determine large radius for local orbitals ---
C     r(nrx)/rmax = (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
      if (lmode0 .and. lpz .gt. 1) then
        call radext(1,nr,nrx,2d0,a,rmax,nrbig,rbig,xx,xx)
C        rbig = rmax * (dexp(a*nrx-a)-1d0)/(dexp(a*nr-a)-1d0)
CC       If rbig>2*rmax, estimate from exp((nrbig-nr)a) = 2
C        if (rbig .gt. 2*rmax) then
C          idn = dlog(2d0)/a
C          if (mod(idn,2) .eq. 1) idn = idn-1
C          nrbig = min(nr+idn,nrx)
C          rbig = rmax * (dexp(a*nrbig-a)-1d0)/(dexp(a*nr-a)-1d0)
C        endif
      elseif(lmode0) then
        nrbig = nr
        rbig = rmax
      elseif (lmode1) then
        call rx('vxtrap not ready for 2s digit mode')
      endif

C --- Points and weights on extended mesh ---
      if (lmode0 .or. lmode1) then
        call radext(10,nr,nrx,2d0,a,rmax,nrbig,rbig,rofi,rwgt)
C        call radmsh(rbig,a,nrbig,rofi)
C        call radwgt(rbig,a,nrbig,rwgt)
C        if (nr .lt. nrbig) rwgt(nr) = rwgt(nr)/2
        return
      endif

C --- Potential on extended mesh ---
      if (.not. lmode2) return
Ckino isanrg is logical function,       call isanrg(nrbig,nr,nrx,'vxtrap:','nrbig',.true.)
      l_dummy_isanrg=isanrg(nrbig,nr,nrx,'vxtrap:','nrbig',.true.)
      if (rbig .lt. rmax) call rx('vxtrap: illegal value for rbig')

      do  isp = 1, nsp
C ... Extrapolate potential, if nrbig > nr
        if (nr .lt. nrbig) then

C       Make r*(true V - vmtz) to facilitate 2nd derivative
          do  ir = 2, nr
            vbig(ir,isp) = (v(ir,isp) - 2*z/rofi(ir) - vmtz)*rofi(ir)
C         vbig(ir,isp) = exp(-2d0*rofi(ir))
          enddo
C       call prrmsh('v*r',rofi,vbig,nr,nr,1)

C       First derivative in gbig; second derivative in gzbig
          k = nr-2*nptirp+1
          if (k .lt. 2) call rxi('vxtrap: nr too small, nr=',nr)
          call poldvm(rofi(k),vbig(k,isp),2*nptirp,nptirp,.false.,1d-8,j,
     .    gbig(k))
          call poldvm(rofi(k),gbig(k),2*nptirp,nptirp,.false.,1d-8,j,
     .    gzbig(k))
C       call prrmsh('v 1st d',rofi(k),gbig(k),2*nptirp,2*nptirp,1)
C       call prrmsh('v 2nd d',rofi(k),gzbig(k),2*nptirp,2*nptirp,1)

C       Use r*V' = (rV)' - V
          phi  = vbig(nr,isp)/rofi(nr)
          dphi = (gbig(nr) - phi)/rofi(nr)
          tphi = gzbig(nr)/vbig(nr,isp)
          call info5(40/10,0,0,' vxtrap: extrapolate V using '//
     .    'vmtz=%;3,3d  v-vmtz=%;3,3d  v''=%;3,3d  (nabla v)/v=%;3,3d',
     .    vmtz,phi,dphi,tphi,0)
          call rwftai(3,rmax,a,nr,nrbig,rofi,phi,dphi,tphi,0,ehv,rs3,
     .    vbig(1,isp))

C       Restore vbig to estimate for true V + 2Z/r
          vbig(1,isp) = v(1,isp)
          do  ir = 2, nrbig
            vbig(ir,isp) = vbig(ir,isp)/rofi(ir) + 2*z/rofi(ir) + vmtz
          enddo

C       For pictures ...
C        do  ir = 2, nrbig
C          vbig(ir,isp) = vbig(ir,isp) - 2*z/rofi(ir)
C        enddo
C        call prrmsh('vbig',rofi,vbig,nrmax,nrbig,nsp)
C        stop

C ... Just copy potential if nr = nrbig
        else
          call dcopy(nr,v(1,isp),1,vbig(1,isp),1)
        endif
      enddo

      end
      subroutine rwftai(mode,rmt,a,nrmt,nrbig,ribig,phi,dphi,tphi,l,
     .ehl,rsml,g)
C- Extend radial wave function outside MT boundary
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :controls approximation in which w.f. is extended
Ci         :1s digit
Ci         :0 Make extended radial mesh (to smaller of 2*rmax or r(nrx))
Ci         :  NO LONGER IMPLEMENTED; use radext
Ci         :1 Use info about sm. Hankel tail to estimate radius to
Ci         :  to which w.f. need be extended to incorporate full charge
Ci         :  (not implemented)
Ci         :2 Compute radial wave function on extended mesh using
Ci         :  rsml,ehl for tail, scaling envelope to match gz.
Ci         :3 Compute radial wave function on extended mesh matching
Ci         :  phi,dphi to h,hdot.  Use ehl for sm. Hankel energy
Ci         :  and find rsm that best matches kinetic energy tphi
Ci         :4 similar to mode 3, except if no KE match, also attempt
Ci         :  to vary ehl.  NB: not implemented.
Ci         :  Used to extrapolate potential.  In this case, g=r*V.
Ci         :5 Similar to mode 2, except scale gz so that value
Ci         :  matches envelope h(rsm,eh), rather than scaling envelope
Ci   rmt   :augmentation radius, in a.u., by species
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nrmt  :number of radial mesh points between origin and rmt
Ci   phi   :value of normalized wave function at rmt, i.e. g/rmt
Ci         :(used for mode = 3)
Ci   dphi  :slope of normalized w.f. at rmt, i.e. (d(g/r)/dr)_rmt
Ci         :(used for mode = 3)
Ci   tphi  :kinetic energy of w.f. at rmt
Ci         :(used for mode = 3)
Ci   l     :quantum number for this wave function
Ci   ehl   :energy of smoothed Hankel tail for local orbital
Ci   rsml  :smoothing radius for smoothed Hankel tail of local orbital
Cio Inputs/Outputs
Cio  nrbig  :number of radial mesh points on extended mesh
Cio         :mode  0 : output
Cio         :mode >0 : input
Cio  ribig  :points for extended radial mesh
Cio         :mode  0 : output
Cio         :mode >0 : input
Co Outputs
Co   g      :radial wave function extended to points nrmt..nrbig
Co          :(done in modes >0)
Cl Local variables
Cl   rbig   :rmax of extended mesh
Cl   nxtn   :number of points between rmt and rbig
Cr Remarks
Cr
Cu Updates
Cu   27 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nrmt,nrbig,l
      double precision a,rmt,phi,dphi,tphi,ribig(*),g(nrbig,2),ehl,rsml
C ... Local parameters
      integer nrx,idn,IPRTDB,nxtn,ir,info,mode0
      parameter (nrx=1501,IPRTDB=10)
      double precision rbig,rwgt(nrx),xi(nrx,0:l),r2(nrx) !,wk(2*nrx)
      double precision fac1,fac2,rsm,ekin,eh,alfa,beta
      double precision hs(0:l),dhs(0:l),ddhs(0:l)
      double precision hsp(0:l),dhsp(0:l),ddhsp(0:l)
C     double precision dot3,sum1,sum2

      mode0 = mod(mode,10)
C     mode1 = mod(mode/10,10)

C --- Get rbig which corresponds to smaller of r(nrx) and 2*rmt ---
      if (mode0 .eq. 0) then
        call rx('rwftai: mode not implemented')

C     elseif (mode0 .eq. 1) then
      elseif (mode0 .eq. 2 .or. mode0 .eq. 5) then
        rbig = ribig(nrbig)
C       Uncomment following 2 lines to replace entire function with envelope
C       nrmt = 1
C       g = 0
        nxtn = nrbig-nrmt+1
        call radwgt(rbig,a,nrbig,rwgt)
        rwgt(nrmt) = rwgt(nrmt)/2
        do  ir = 1, nxtn
          r2(ir) = ribig(ir-1+nrmt)**2
        enddo
C       Hankel function scaled to match g at nrmt
c        call hansr(rsml,0,l,1,l,ehl,r2,nrx,nxtn,idn,wk,11,xi)
        call hansr(rsml,0,l,1,l,ehl,r2,nrx,nxtn,idn,11,xi)
C       fac1 = phi/xi(1,l)
        fac1 = g(nrmt,1)/rmt/xi(1,l)
        fac2 = g(nrmt,2)/rmt/xi(1,l)
        if (mode0 .eq. 5) then
          call dscal(nrmt,1/fac1,g(1,1),1)
          call dscal(nrmt,1/fac1,g(1,2),1)
          fac2 = fac2 / fac1
          fac1 = 1
        endif
        do  ir = 1, nxtn
          g(ir-1+nrmt,1) = xi(ir,l)*ribig(ir-1+nrmt) * fac1
          g(ir-1+nrmt,2) = xi(ir,l)*ribig(ir-1+nrmt) * fac2
        enddo
CC       Integral of gg from 0 to rmt, neglecting small component
C        sum1 = dot3(nrmt,rwgt,g,g)
CC       Integral of gg from rmt to rbig, neglecting small component
C        sum2 = dot3(nxtn,rwgt(nrmt),g(nrmt,1),g(nrmt,1))
C       call info2(IPRTDB,0,0,' rwftai: <gg>_MT=%g  <gg>_big=%g',
C    .    sum1,sum2)
C       call prrmsh('gzbig',ribig,g,nrbig,nrbig,2)


      elseif (mode0 .eq. 3) then
        eh = min(ehl,-.05d0)
        call pshpr(110)
        call mtchre(104,l,rsml,rmt,0d0,0d0,rmt,rmt,phi,dphi,
     .  tphi,dphi,rsm,eh,ekin,info)
        call poppr
        call mtchae(1,rsm,eh,l,rmt,phi,dphi,0d0,0d0,alfa,beta)
        do  ir = nrmt+1, nrbig
          call hansmd(12,ribig(ir),eh,rsm,l,hs,dhs,ddhs,hsp,dhsp,ddhsp)
          g(ir,1) = (alfa*hs(l) + beta*hsp(l)) * ribig(ir)
        enddo

      else
        call rxi('rwftai: mode not implemented:',mode)
      endif

      end
C      subroutine xxxx
C      print *, 'xxx'
C      end

