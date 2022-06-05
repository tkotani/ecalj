subroutine vxc0sp(a,b,rofi,rho,nr,v,rho0,rep,rmu,nsp,exrmx)
  use m_lmfinit,only: lxcf_g=>lxcf, stdo
  use m_ropyln,only: ropyln
  !!- Adds xc part to spherical potential, makes integrals rmu and rep
  !! ----------------------------------------------------------------------
  ! i Inputs
  ! i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  ! i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  ! i   rofi  :radial mesh points
  ! i   rho   :density = (true rho)*(4*pi*r**2)
  ! i   nr    :number of radial mesh points
  ! i   nsp   :2 for spin-polarized case, otherwise 1
  ! i  globalvariables%lxcf  :type of xc potential
  !!    1 : VWN
  !!    2 : Barth-Hedin
  !!   103: GGA PBE
  ! o Outputs
  ! o   v     :vxc is added to v
  ! o   rho0  :density extrapolated to origin
  ! o   rep   :integral rho * exc.
  ! o   rmu   :integral rho * vxc.
  ! o   exrmx :exchange energy density at rmax
  !!   vxc0sp contained in lm7.0 beta compiled by M.van Schilfgaarde.
  !! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nr,nsp
  double precision :: a,b,rofi(nr),v(nr,nsp),rho(nr,nsp), &
       rep(nsp),rmu(nsp),rho0(2),qs(2),exrmx
  ! ... Local parameters
  !      parameter (nrmx=1501)
  double precision :: pi,rho2,rho3, &
       ub4pi,wgt,exc(nr),vxc(nr,2),rp(nr,2),repl(2),rmul(2) !exc(nrmx)-->exc(nr) automatic array.
  integer:: lx ,  i , ir , isp , iprint , nglob, lxcfun,lxcf
  real(8) ,allocatable :: grh_rv(:)
  real(8) ,allocatable :: agrh_rv(:)
  real(8) ,allocatable :: ggrh_rv(:)
  real(8) ,allocatable :: grgag_rv(:)
  real(8) ,allocatable :: excx_rv(:)
  real(8) ,allocatable :: excc_rv(:)
  real(8) ,allocatable :: vxcx_rv(:)
  real(8) ,allocatable :: vxcc_rv(:),vx2(:),vc2(:)
  character (2) :: st

  integer:: nth,nph,lxcg,np,lmax,nlm
  integer,parameter:: nnn=122
  real(8):: p(3,nnn),wp(nnn),rwgt(nr) ,p2(nnn,3),r2(nnn)
  real(8),allocatable:: gyl(:,:),yl(:),rll(:,:),vxc_(:,:)
  logical:: newv=.true.
  pi = 4d0*datan(1d0)
  ub4pi = 1d0/(4d0*pi)
  lxcfun = lxcf_g
  lxcf = mod(lxcfun,100)
  ! --- Add background rho to calculate vxc ---
  !      rhobg = 0d0
  !      call getsyv('rhobg',rhobg,i)
  !      call addzbk(rofi,nr,1,nsp,rho,rhobg,1d0)

  !! === Extrapolate rho to origin ===
  do  10  isp = 1, nsp
     rep(isp) = 0d0
     rmu(isp) = 0d0
     rho2 = rho(2,isp)/rofi(2)**2
     rho3 = rho(3,isp)/rofi(3)**2
     rho0(isp) = ub4pi*(rho2*rofi(3)-rho3*rofi(2))/(rofi(3)-rofi(2))
10 enddo
  !! === Make true rho ===
  do   isp = 1, nsp
     rp(1,isp) = rho0(isp)
     do   ir = 2, nr
        rp(ir,isp) = rho(ir,isp)*ub4pi/rofi(ir)**2
     enddo
  enddo

  !      print *,'end of do 20'
  !! === Generate vxc,exc on a mesh ===
  !      print *,'end of do 20 xxxx lxcf2=',lxcf2
  if (lxcfun==103) then
     if( .NOT. newv) then
        allocate(excx_rv(nr))
        allocate(excc_rv(nr))
        allocate(vxcx_rv(nr))
        allocate(vxcc_rv(nr))
        allocate(vx2(nr),vc2(nr))
        call evxcp(rp,rp(1,2),nr,nsp,lxcf,excx_rv,excc_rv,exc, &
             vxcx_rv,vx2, vxcc_rv,vc2,vxc,vxc(1,2))
        deallocate(vx2,vc2)
        do isp = 1, nsp
           vxc(1,isp) = (vxc(2,isp)*rofi(3)-vxc(3,isp)*rofi(2))/(rofi(3)-rofi(2))
        enddo
        deallocate(vxcc_rv,vxcx_rv,excc_rv,excx_rv)
     else

        !         np=1
        !         wp(1)=4d0*pi

        nph=0
        call fpiint(-4,nph,np,p,wp)
        p2=transpose(p)
        call radwgt(rofi(nr),a,nr,rwgt)
        lmax=0
        nlm=(lmax+1)**2 !yl(l=1) is needed for gradfl called in vxcnls
        allocate (yl(np*(lmax+2)**2),gyl(np*nlm,3)) !nlm=1 Only s channel. Sphelical atom.
        call ropyln(np,p2(1,1),p2(1,2),p2(1,3),lmax+1,np,yl,r2)
        !         call ropylg(1,lmax,nlm,np,np,p2,p2(1+np),p2(1+2*np),r2,yl,gyl)
        gyl=0d0
        !         print *,'goto vxcnls only s channel',sum(abs(gyl))
        !         nlm=1
        !         allocate(rll(nr,nsp))
        !         rll(:,1:nsp)=
        !         print *,'goto vxcnls',sum(rll)
        call vxcnls(a,rofi,0,nr,np,nlm,nsp, &
             yl,gyl,rwgt,wp,sqrt(4d0*pi)*rp(:,1:nsp),lxcfun-100,  vxc,rep,rmu)
        deallocate(yl,gyl)
        vxc= vxc/sqrt(4d0*pi)
        !         deallocate(rll)

     endif
  elseif(lxcfun==1 .OR. lxcfun==2) then
     allocate(excx_rv(nr))
     allocate(excc_rv(nr))
     allocate(vxcx_rv(nr))
     allocate(vxcc_rv(nr))
     if (nsp == 1) then
        call evxcv ( rp , rp , nr , nsp , lxcf, &
             exc , excx_rv , excc_rv , &
             vxc , vxcx_rv , vxcc_rv )
     else
        rp(1:nr,2)= rp(1:nr,1)+rp(1:nr,2)
        !          call dpadd(rp(1,2),rp,1,nr,1d0)
        call evxcv ( rp ( 1 , 2 ) , rp , nr , 2 , lxcf, exc , excx_rv &
             , excc_rv , vxc , vxcx_rv , vxcc_rv )
        rp(1:nr,2) = rp(1:nr,2) - rp(1:nr,1)
        rp(1:nr,1) = rp(1:nr,1) + rp(1:nr,2)
        !          call dpadd(rp(1,2),rp,1,nr,-1d0)
        !          call dpadd(rp,rp(1,2),1,nr,1d0)
        call evxcv ( rp , rp ( 1 , 2 ) , nr , 2 , lxcf, exc , excx_rv &
             , excc_rv , vxc ( 1 , 2 ) , vxcx_rv , vxcc_rv )
        rp(1:nr,1) = rp(1:nr,1) - rp(1:nr,2)
        !          call dpadd(rp,rp(1,2),1,nr,-1d0)
     endif
     deallocate(vxcc_rv,vxcx_rv,excc_rv,excx_rv)
     !! === Integrals ===
     do  14  i  = 1, nsp
        !        qs(i)  = 0d0
        rep(i) = 0d0
        rmu(i) = 0d0
        do  12  ir = 1, nr
           wgt = 2*(mod(ir+1,2)+1)/3d0
           if (ir == 1 .OR. ir == nr) wgt = 1d0/3d0
           wgt = wgt * a*(rofi(ir)+b)
           !          qs(i)  = qs(i)  + wgt*rho(ir,i)
           rep(i) = rep(i) + wgt*rho(ir,i)*exc(ir)
           rmu(i) = rmu(i) + wgt*rho(ir,i)*vxc(ir,i)
12      enddo
        !        repl(i) = rep(i)
        !        rmul(i) = rmu(i)
14   enddo
  endif

  !$$$C --- Gradient correction ---
  !      if(.false.) then
  if (lxcfun==103 .AND. ( .NOT. newv)) then
     !      if (lxcg /= 0) then
     !        allocate(grh_rv(nrmx*nsp))
     !        allocate(ggrh_rv(nrmx*nsp))
     !        allocate(agrh_rv(nrmx*(3*nsp-2)))
     !        allocate(grgag_rv(nrmx*(2*nsp-1)))
     !        call vxcgr2 ( nr , nsp , nrmx , rofi , rp , grh_rv , ggrh_rv
     !     .  , agrh_rv , grgag_rv , exc , vxc )
     !        deallocate(grgag_rv,agrh_rv,ggrh_rv,grh_rv)
     !        call vxcgr2 ( nr , nsp , nrmx , rofi , rp , exc , vxc )
     call vxcgr2 ( nr , nsp , nr , rofi , rp , exc , vxc )
     ! ...   Redo integrals, with gradient correction
     do  24  i  = 1, nsp
        repl(i) = rep(i)
        rmul(i) = rmu(i)
        rep(i) = 0d0
        rmu(i) = 0d0
        do  22  ir = 1, nr
           wgt = 2*(mod(ir+1,2)+1)/3d0
           if (ir == 1 .OR. ir == nr) wgt = 1d0/3d0
           wgt = wgt * a*(rofi(ir)+b)
           rep(i) = rep(i) + wgt*rho(ir,i)*exc(ir)
           rmu(i) = rmu(i) + wgt*rho(ir,i)*vxc(ir,i)
22      enddo
24   enddo
  endif
  !      endif

  ! ccccccccccccccccccccccccccccccccccc
  !      do   i = 1, nsp
  !      do  ir = 1, nr
  !         print *, 'pppp ', ir,i,  vxc(ir,i)
  !      enddo
  !      enddo
  !      stop 'xxxxxxxxxxxx vxc0sp xxxxxxxxxxxxxxx'
  ! ccccccccccccccccccccccccccccccccc

  !! --- Add to V ---
  !      call dpadd(v,vxc,1,nr,1d0)
  !      if (nsp .eq. 2) call dpadd(v(1,2),vxc(1,2),1,nr,1d0)
  v(:,1:nsp)=v(:,1:nsp)+vxc(:,1:nsp)
  exrmx = exc(nr)

  !! --- Undo background rho for purposes of calculating vxc ---
  !     call addzbk(rofi,nr,1,nsp,rho,rhobg,-1d0)
  if (iprint() < 80) return
  !      if (lxcg .eq. 0) write(stdo,333)
  !  333 format(/' vxc0sp: reps(l)     rmu(l)')
  !      if (lxcg .ne. 0) write(stdo,334)
  !  334 format(/' vxc0sp: reps(l)     rmu(l)      reps(nl)    rmu(nl)')
  write(stdo,"(/' vxc0sp: reps        rmu   ')")
  do i = 1, nsp
     st = ' '
     if (i < nsp) st = 'up'
     if (i == 2)   st = 'dn'
     write(stdo,335) st, rep(i),  rmu(i)
  enddo
  if(nsp==2) write(stdo,335) '  ', rep(1)+rep(2), rmu(1)+rmu(2)
335 format(1x,a2,2x,4f12.6)
  !      if (lxcg .ne. 0) write(stdo,335) st, repl(i), rmul(i),
  !     .  rep(i)-repl(i), rmu(i)-rmul(i)
  !      if (lxcg .eq. 0) write(stdo,335) st, rep(i),  rmu(i)
  !      if (lxcg .ne. 0) write(stdo,335) st, repl(i), rmul(i),
  !     .  rep(i)-repl(i), rmu(i)-rmul(i)
  !      if (nsp .eq. 2 .and. lxcg .eq. 0)
  !     .write(stdo,335) '  ', rep(1)+rep(2), rmu(1)+rmu(2)
  !      if (nsp .eq. 2 .and. lxcg .ne. 0)
  !     .write(stdo,335) '  ', repl(1)+repl(2), rmul(1)+rmul(2),
  !     .rep(1)+rep(2)-repl(1)-repl(2), rmu(1)+rmu(2)-rmul(1)-rmul(2)
  !$$$      if (lxcg .eq. 0) print 333
  !$$$      if (lxcg .ne. 0) print 334
  !$$$  333 format(/' vxc0sp: reps(l)     rmu(l)')
  !$$$  334 format(/' vxc0sp: reps(l)     rmu(l)      reps(nl)    rmu(nl)')
  !$$$      do  30  i = 1, nsp
  !$$$        st = ' '
  !$$$        if (i .lt. nsp) st = 'up'
  !$$$        if (i .eq. 2)   st = 'dn'
  !$$$        if (lxcg .eq. 0) print 335, st, rep(i),  rmu(i)
  !$$$        if (lxcg .ne. 0) print 335, st, repl(i), rmul(i),
  !$$$     .  rep(i)-repl(i), rmu(i)-rmul(i)
  !$$$  335   format(1x,a2,2x,4f12.6)
  !$$$   30 continue
  !$$$      if (nsp .eq. 2 .and. lxcg .eq. 0)
  !$$$     .print 335, '  ', rep(1)+rep(2), rmu(1)+rmu(2)
  !$$$      if (nsp .eq. 2 .and. lxcg .ne. 0)
  !$$$     .print 335, '  ', repl(1)+repl(2), rmul(1)+rmul(2),
  !$$$     .rep(1)+rep(2)-repl(1)-repl(2), rmu(1)+rmu(2)-rmul(1)-rmul(2)
end subroutine vxc0sp
