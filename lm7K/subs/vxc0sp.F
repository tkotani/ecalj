      subroutine vxc0sp(a,b,rofi,rho,nr,v,rho0,rep,rmu,nsp,exrmx)


C- Adds xc part to spherical potential, makes integrals rmu and rep
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rofi  :radial mesh points
Ci   rho   :density = (true rho)*(4*pi*r**2)
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Co   v     :vxc is added to v
Co   rho0  :density extrapolated to origin
Co   rep   :integral rho * exc.
Co   rmu   :integral rho * vxc.
Co   exrmx :exchange energy density at rmax
Cl Local variables
Cl  lx     :type of local xc potential
Cl           1's digit specifies local potential (evxcv)
Cl           1: Ceperly Alder, 2: Hedin-Lundqvist
Cl           3: Perdew-Wang '91, 3: Burke, Perdew and Ernzerhof
Cl         10's  digit specifies gradient corrections
Cl           1: Langreth-Mehl-Hu 2: PW91,
Cl           3: PBE, 4: PBE with Becke exchange
Cr Remarks
Cu Updates
Cu   18 Jun 04 lxcfun is no longer used and should be deleted from cmd-line
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nsp,nrmx
      double precision a,b,rofi(nr),v(nr,nsp),rho(nr,nsp),
     .  rep(nsp),rmu(nsp),rho0(2),qs(2),exrmx
C ... Local parameters
      parameter (nrmx=1501)
      double precision pi,rho2,rho3,
     .  ub4pi,wgt,exc(nrmx),vxc(nrmx,2),rp(nrmx,2),repl(2),rmul(2)
      integer:: lx , lxcf2 , lxcg2 , i , ir , isp , iprint , nglob
      real(8) ,allocatable :: grh_rv(:)
      real(8) ,allocatable :: agrh_rv(:)
      real(8) ,allocatable :: ggrh_rv(:)
      real(8) ,allocatable :: grgag_rv(:)

      real(8) ,allocatable :: excx_rv(:)
      real(8) ,allocatable :: excc_rv(:)
      real(8) ,allocatable :: vxcx_rv(:)
      real(8) ,allocatable :: vxcc_rv(:)

      character *2 st
C ... Heap
      integer w(1)
      common /w/ w

      lx    = nglob('lxcf')
      lxcf2 = mod(lx,10)
      lxcg2 = lx/10
      if (nr .gt. nrmx) call rx('vxc0sp: nr > nrmx')
      pi = 4d0*datan(1d0)
      ub4pi = 1d0/(4d0*pi)

C --- Add background rho to calculate vxc ---
*      rhobg = 0d0
*      call getsyv('rhobg',rhobg,i)
*      call addzbk(rofi,nr,1,nsp,rho,rhobg,1d0)

C --- Extrapolate rho to origin ---
      do  10  isp = 1, nsp
        rep(isp) = 0d0
        rmu(isp) = 0d0
        rho2 = rho(2,isp)/rofi(2)**2
        rho3 = rho(3,isp)/rofi(3)**2
        rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
   10 continue

C --- Make true rho ---
      do  20  isp = 1, nsp
        rp(1,isp) = rho0(isp)
      do  20  ir = 2, nr
        rp(ir,isp) = rho(ir,isp)*ub4pi/rofi(ir)**2
   20 continue

C --- Generate vxc,exc on a mesh ---
      if (lxcf2 .gt. 2) then
        call evxcp(rp,rp(1,2),nr,nsp,lxcf2,exc,vxc,vxc(1,2))
        do  60  isp = 1, nsp
          vxc(1,isp) =
     .      (vxc(2,isp)*rofi(3)-vxc(3,isp)*rofi(2))/(rofi(3)-rofi(2))
   60   continue
      else
        allocate(excx_rv(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) excx_rv(:)=0.0d0

        allocate(excc_rv(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) excc_rv(:)=0.0d0

        allocate(vxcx_rv(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) vxcx_rv(:)=0.0d0

        allocate(vxcc_rv(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) vxcc_rv(:)=0.0d0

        if (nsp .eq. 1) then
          call evxcv ( rp , rp , nr , nsp , lxcf2 , exc , excx_rv , 
     .     excc_rv , vxc , vxcx_rv , vxcc_rv )

        else
          call dpadd(rp(1,2),rp,1,nr,1d0)
          call evxcv ( rp ( 1 , 2 ) , rp , nr , 2 , lxcf2 , exc , excx_rv 
     .     , excc_rv , vxc , vxcx_rv , vxcc_rv )

          call dpadd(rp(1,2),rp,1,nr,-1d0)
          call dpadd(rp,rp(1,2),1,nr,1d0)
          call evxcv ( rp , rp ( 1 , 2 ) , nr , 2 , lxcf2 , exc , excx_rv 
     .     , excc_rv , vxc ( 1 , 2 ) , vxcx_rv , vxcc_rv )

          call dpadd(rp,rp(1,2),1,nr,-1d0)
        endif
        if (allocated(vxcc_rv)) deallocate(vxcc_rv)
        if (allocated(vxcx_rv)) deallocate(vxcx_rv)
        if (allocated(excc_rv)) deallocate(excc_rv)
        if (allocated(excx_rv)) deallocate(excx_rv)

      endif

C --- Integrals ---
      do  14  i  = 1, nsp
        qs(i)  = 0d0
        rep(i) = 0d0
        rmu(i) = 0d0
        do  12  ir = 1, nr
          wgt = 2*(mod(ir+1,2)+1)/3d0
          if (ir .eq. 1 .or. ir .eq. nr) wgt = 1d0/3d0
          wgt = wgt * a*(rofi(ir)+b)
          qs(i)  = qs(i)  + wgt*rho(ir,i)
          rep(i) = rep(i) + wgt*rho(ir,i)*exc(ir)
          rmu(i) = rmu(i) + wgt*rho(ir,i)*vxc(ir,i)
   12   continue
        repl(i) = rep(i)
        rmul(i) = rmu(i)
   14 continue

C --- Gradient correction ---
      if (lxcg2 .ne. 0) then
        allocate(grh_rv(nrmx*nsp))
Ctakao_ZeroClear_NotRequiered          if (nrmx*nsp<0) grh_rv(:)=0.0d0

        allocate(ggrh_rv(nrmx*nsp))
Ctakao_ZeroClear_NotRequiered          if (nrmx*nsp<0) ggrh_rv(:)=0.0d0

        allocate(agrh_rv(nrmx*(3*nsp-2)))
Ctakao_ZeroClear_NotRequiered          if (nrmx*(3*nsp-2)<0) agrh_rv(:)=0.0d0

        allocate(grgag_rv(nrmx*(2*nsp-1)))
Ctakao_ZeroClear_NotRequiered          if (nrmx*(2*nsp-1)<0) grgag_rv(:)=0.0d0

        call vxcgr2 ( nr , nsp , nrmx , rofi , rp , grh_rv , ggrh_rv 
     .   , agrh_rv , grgag_rv , exc , vxc )

        if (allocated(grgag_rv)) deallocate(grgag_rv)
        if (allocated(agrh_rv)) deallocate(agrh_rv)
        if (allocated(ggrh_rv)) deallocate(ggrh_rv)
        if (allocated(grh_rv)) deallocate(grh_rv)

C ...   Redo integrals, with gradient correction
        do  24  i  = 1, nsp
          repl(i) = rep(i)
          rmul(i) = rmu(i)
          rep(i) = 0d0
          rmu(i) = 0d0
          do  22  ir = 1, nr
            wgt = 2*(mod(ir+1,2)+1)/3d0
            if (ir .eq. 1 .or. ir .eq. nr) wgt = 1d0/3d0
            wgt = wgt * a*(rofi(ir)+b)
C         print *, ir, exc(ir), vxc(ir,1)
            rep(i) = rep(i) + wgt*rho(ir,i)*exc(ir)
            rmu(i) = rmu(i) + wgt*rho(ir,i)*vxc(ir,i)
   22     continue
   24   continue
      endif

C --- Add to V ---
      call dpadd(v,vxc,1,nr,1d0)
      if (nsp .eq. 2) call dpadd(v(1,2),vxc(1,2),1,nr,1d0)
      exrmx = exc(nr)

C --- Undo background rho for purposes of calculating vxc ---
*     call addzbk(rofi,nr,1,nsp,rho,rhobg,-1d0)

      if (iprint() .lt. 80) return
      if (lxcg2 .eq. 0) print 333
      if (lxcg2 .ne. 0) print 334
  333 format(/' vxc0sp: reps(l)     rmu(l)')
  334 format(/' vxc0sp: reps(l)     rmu(l)      reps(nl)    rmu(nl)')
      do  30  i = 1, nsp
        st = ' '
        if (i .lt. nsp) st = 'up'
        if (i .eq. 2)   st = 'dn'
        if (lxcg2 .eq. 0) print 335, st, rep(i),  rmu(i)
        if (lxcg2 .ne. 0) print 335, st, repl(i), rmul(i),
     .  rep(i)-repl(i), rmu(i)-rmul(i)
  335   format(1x,a2,2x,4f12.6)
   30 continue
      if (nsp .eq. 2 .and. lxcg2 .eq. 0)
     .  print 335, '  ', rep(1)+rep(2), rmu(1)+rmu(2)
      if (nsp .eq. 2 .and. lxcg2 .ne. 0)
     .  print 335, '  ', repl(1)+repl(2), rmul(1)+rmul(2),
     .  rep(1)+rep(2)-repl(1)-repl(2), rmu(1)+rmu(2)-rmul(1)-rmul(2)

      end subroutine vxc0sp


