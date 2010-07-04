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
Cl  lxcfun     :type of local xc potential
Cl           1's digit specifies local potential (evxcv)
Cl           1: Ceperly Alder, 2: Hedin-Lundqvist
Cl           3: Perdew-Wang '91, 3: Burke, Perdew and Ernzerhof
Cl         100's  digit specifies gradient corrections
Cl           1: Langreth-Mehl-Hu 2: PW91,      
Cl           3: PBE, 4: PBE with Becke exchange
Cr Remarks
Cu Updates
Cu   18 Jun 04 lxcfun is no longer used and should be deleted from cmd-line
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nr,nsp,nrmx,stdo
      double precision a,b,rofi(nr),v(nr,nsp),rho(nr,nsp),
     .rep(nsp),rmu(nsp),rho0(2),qs(2),exrmx
C ... Local parameters
      parameter (nrmx=1501)
      double precision pi,rho2,rho3,
     .ub4pi,wgt,exc(nrmx),vxc(nrmx,2),rp(nrmx,2),repl(2),rmul(2)
      integer:: lx , lxcf , lxcg , i , ir , isp , iprint , nglob, lxcfun
      real(8) ,allocatable :: grh_rv(:)
      real(8) ,allocatable :: agrh_rv(:)
      real(8) ,allocatable :: ggrh_rv(:)
      real(8) ,allocatable :: grgag_rv(:)

      real(8) ,allocatable :: excx_rv(:)
      real(8) ,allocatable :: excc_rv(:)
      real(8) ,allocatable :: vxcx_rv(:)
      real(8) ,allocatable :: vxcc_rv(:),vx2(:),vc2(:)

      character *2 st
C ... Heap
c      integer w(1)
c      common /w/ w

      stdo = nglob('stdo')

c old version before lm72
c      lx    = nglob('lxcf')
c      lxcf2 = mod(lx,10)
c      lxcg = lx/10
clm72
      lxcfun = nglob('lxcf')
      lxcf = mod(lxcfun,100)
      lxcg = lxcfun/100

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
c      print *,'end of do 10'

C --- Make true rho ---
      do  20  isp = 1, nsp
        rp(1,isp) = rho0(isp)
        do  20  ir = 2, nr
          rp(ir,isp) = rho(ir,isp)*ub4pi/rofi(ir)**2
   20 continue
c      print *,'end of do 20'
C --- Generate vxc,exc on a mesh ---
      allocate(excx_rv(nr))
      allocate(excc_rv(nr))
      allocate(vxcx_rv(nr))
      allocate(vxcc_rv(nr))
c      print *,'end of do 20 xxxx lxcf2=',lxcf2
      if (lxcf > 2) then
        allocate(vx2(nr),vc2(nr))
c        call evxcp(rp,rp(1,2),nr,nsp,lxcf2,exc,vxc,vxc(1,2))
         call evxcp(rp,rp(1,2),nr,nsp,lxcf,excx_rv,excc_rv,exc, !w(oex),w(oec),exc,
c     .    w(ovx),w(ovx2),w(ovc),w(ovc2),vxc,vxc(1,2))
     .    vxcx_rv,vx2, vxcc_rv,vc2,vxc,vxc(1,2))
        deallocate(vx2,vc2)
        do  60  isp = 1, nsp
          vxc(1,isp) =
     .    (vxc(2,isp)*rofi(3)-vxc(3,isp)*rofi(2))/(rofi(3)-rofi(2))
   60   continue
      else
        if (nsp .eq. 1) then
          call evxcv ( rp , rp , nr , nsp , lxcf, 
     .    exc , excx_rv , excc_rv , 
     .    vxc , vxcx_rv , vxcc_rv )

        else
          call dpadd(rp(1,2),rp,1,nr,1d0)
          call evxcv ( rp ( 1 , 2 ) , rp , nr , 2 , lxcf, exc , excx_rv 
     .    , excc_rv , vxc , vxcx_rv , vxcc_rv )
          call dpadd(rp(1,2),rp,1,nr,-1d0)
          call dpadd(rp,rp(1,2),1,nr,1d0)
          call evxcv ( rp , rp ( 1 , 2 ) , nr , 2 , lxcf, exc , excx_rv 
     .    , excc_rv , vxc ( 1 , 2 ) , vxcx_rv , vxcc_rv )

          call dpadd(rp,rp(1,2),1,nr,-1d0)
        endif

      endif
      deallocate(vxcc_rv,vxcx_rv,excc_rv,excx_rv)

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
      if (lxcg /= 0) then
c        allocate(grh_rv(nrmx*nsp))
c        allocate(ggrh_rv(nrmx*nsp))
c        allocate(agrh_rv(nrmx*(3*nsp-2)))
c        allocate(grgag_rv(nrmx*(2*nsp-1)))
c        call vxcgr2 ( nr , nsp , nrmx , rofi , rp , grh_rv , ggrh_rv 
c     .  , agrh_rv , grgag_rv , exc , vxc )
c        deallocate(grgag_rv,agrh_rv,ggrh_rv,grh_rv)
        call vxcgr2 ( nr , nsp , nrmx , rofi , rp , exc , vxc )

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
      if (lxcg .eq. 0) write(stdo,333)
      if (lxcg .ne. 0) write(stdo,334)
  333 format(/' vxc0sp: reps(l)     rmu(l)')
  334 format(/' vxc0sp: reps(l)     rmu(l)      reps(nl)    rmu(nl)')
      do  30  i = 1, nsp
      st = ' '
      if (i .lt. nsp) st = 'up'
      if (i .eq. 2)   st = 'dn'
      if (lxcg .eq. 0) write(stdo,335) st, rep(i),  rmu(i)
      if (lxcg .ne. 0) write(stdo,335) st, repl(i), rmul(i),
     .  rep(i)-repl(i), rmu(i)-rmul(i)
  335 format(1x,a2,2x,4f12.6)
   30 continue
      if (nsp .eq. 2 .and. lxcg .eq. 0)
     .write(stdo,335) '  ', rep(1)+rep(2), rmu(1)+rmu(2)
      if (nsp .eq. 2 .and. lxcg .ne. 0)
     .write(stdo,335) '  ', repl(1)+repl(2), rmul(1)+rmul(2),
     .rep(1)+rep(2)-repl(1)-repl(2), rmu(1)+rmu(2)-rmul(1)-rmul(2)
C$$$      if (lxcg .eq. 0) print 333
C$$$      if (lxcg .ne. 0) print 334
C$$$  333 format(/' vxc0sp: reps(l)     rmu(l)')
C$$$  334 format(/' vxc0sp: reps(l)     rmu(l)      reps(nl)    rmu(nl)')
C$$$      do  30  i = 1, nsp
C$$$        st = ' '
C$$$        if (i .lt. nsp) st = 'up'
C$$$        if (i .eq. 2)   st = 'dn'
C$$$        if (lxcg .eq. 0) print 335, st, rep(i),  rmu(i)
C$$$        if (lxcg .ne. 0) print 335, st, repl(i), rmul(i),
C$$$     .  rep(i)-repl(i), rmu(i)-rmul(i)
C$$$  335   format(1x,a2,2x,4f12.6)
C$$$   30 continue
C$$$      if (nsp .eq. 2 .and. lxcg .eq. 0)
C$$$     .print 335, '  ', rep(1)+rep(2), rmu(1)+rmu(2)
C$$$      if (nsp .eq. 2 .and. lxcg .ne. 0)
C$$$     .print 335, '  ', repl(1)+repl(2), rmul(1)+rmul(2),
C$$$     .rep(1)+rep(2)-repl(1)-repl(2), rmu(1)+rmu(2)-rmul(1)-rmul(2)
      end subroutine vxc0sp


