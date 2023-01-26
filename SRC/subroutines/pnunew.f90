subroutine pnunew(eferm)
  use m_ftox
  use m_MPItk,only: master_mpi
  use m_lmfinit,only:nbas,nsp,ispec,sspec=>v_sspec,ham_frzwf,idmodis=>idmod,slabl,&
       pmin=>ham_pmin,pmax=>ham_pmax,n0,mxcst4
  use m_mkrout,only: hbyl=>hbyl_rv,qbyl=>qbyl_rv
  use m_lgunit,only:stdo
  use m_phidx,only: phidx

  use m_density,only: v0pot,pnuall,pnzall !output
  
  !!= Makes new boundary conditions pnu for phi,phidot =
  !!* takao sep2010. P is setted to that at efermi if PZ is for semicore.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   lfrzw :0, float pnu to band !G, provided IDMOD=0 for that channel
  !i         :1 freeze all pnu for all species.
  !i         :  NB: pnu are also frozen for specific species
  !i         :      that have nonzero 4's bit of species->mxcst.
  !i   pmin  :lower bound for fractional part of P
  !i   pmax  :upper bound for fractional part of P
  !i   hab   :<u,s | H | u,s> for each pair uu, us, su, ss; see Remarks
  !i   sab   :<u,s | 1 | u,s>
  !i         :NB: hab and sab are no longer used.
  !i   qbyl  :l-decomposed charge
  !i   hbyl  :l-decomposed eigenvalue sum
  !l Local variables
  !l   lfrzv  if T, freeze valence pnu
  !l   lfrzz  if T, freeze local orbital pnz
  !r Remarks
  !u Updates
  !u   28 Jun 06 Handles case idmod=3; New constraints pmin,pmax
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u    9 May 02 Added species-specific freezing of pnu
  !u   22 Dec 01 Adjustments to accomodate changes in phidx
  !u   20 Sep 01 Some patches to work in pathological cases
  !u   17 Sep 01 When local orbital present, allow semiore pnz to float
  !u   28 Aug 01 Extended to local orbitals.  For now, freeze pnu
  ! ----------------------------------------------------------------------
  implicit none
  logical lpz,lfrzv,lfrzz
  integer:: idmod(n0) ,ipr , ib , is , lmxa , l , ipqn , m , isp , nr , nn , mxcst ! 
  real(8) ,allocatable :: g_rv(:)
  real(8) ,allocatable :: gp_rv(:)
  real(8) ,allocatable :: rofi_rv(:)
  real(8) ,allocatable :: v0i_rv(:)
  double precision pi,rmt,p1,ebar,a,d0l,pfree,pold,ptry,z,val(5),slo(5),dl,phi,dphi
!  real(8):: pnu(n0,2),pnz(n0,2)
  real(8),pointer:: pnu(:,:),pnz(:,:)
  double precision ez,umegam,phip,dphip,dlphi,dlphip,cz
  double precision pznew,fi(0:10),gi(0:10),xx,dnz
  character spid*8
  integer ::iwdummy ,i_copy_size,nnz,nnv
  real(8):: eferm,eee
  logical:: lsemicorepz,phispinsym,cmdopt0
  real(8):: pmean
  character strn*120
  call tcn('pnunew')
  pi = 4d0*datan(1d0)
  call getpr(ipr)
  if(ipr>20)write(stdo,ftox)' Make new boundary conditions for phi,phidot..'
  if(ipr>40) then
     print *, ' pnunew: ebar: '
     print *, '   without lo    : ebar = center of gravity of occupied states'
     print *, '   with lo & PZ>P: ebar for lo is meaningless(zero is shown). Use empty-sphere PZ.'
     print *, '                   ebar for valence is at the center of gravity of occ. states.'
     print *, '   with lo & PZ<P: ebar for lo is by atomic calculation phidx for given PZ'
     print *, '                   ebar for valence is at the Fermi energy.'
     print *, ' idmod=0-> If P=prty(fractional part is log.-derivative)<pfeee, we use pfree.'
  endif
  ! --- For each site, do ---
  do  ib = 1, nbas
     is=ispec(ib) !ssite(ib)%spec
     pnu=>pnuall(:,:,ib)
     pnz=>pnzall(:,:,ib)
!     pnu=ssite(ib)%pnu
!     pnz=ssite(ib)%pz
     lmxa=sspec(is)%lmxa
     rmt=sspec(is)%rmt
     idmod=idmodis(:,is)
     if (lmxa .eq. -1) goto 10
     spid=slabl(is) !sspec(is)%name
     if (mxcst4(is)) idmod=1 !call ivset(idmod,1,n0,1)
     if (ipr >40) write(stdo,320) ib,is,spid
320  format(/' site',i5,'   species',i4,':',a)
     if (ipr >40) write(stdo,311)
     do  l = 0, lmxa
        do  isp = 1, nsp
           m = l+1
           p1 = 2d10
           pznew = 2d10
           lfrzv = (mod(idmod(m),10).ne.0 .and. mod(idmod(m),10).ne.3).or. ham_frzwf !lfrzw .ne. 0
           lfrzz = lfrzv !freezing swiches
           lpz = pnz(m,1) .ne. 0
           lsemicorepz = .false.
           if(lpz)then
              nnv = int(mod(pnu(m,1),10d0))-l-1
              nnz = int(mod(pnz(m,1),10d0))-l-1
              if(nnz<nnv) then !this means semicore for local orbital.
                 lsemicorepz = .true.
              endif
           endif
           if (dabs(qbyl(m,isp,ib)) .gt. 1d-8) then
              ebar = hbyl(m,isp,ib)/qbyl(m,isp,ib)
              is=ispec(ib) !ssite(ib)%spec
              z=sspec(is)%z
              a=sspec(is)%a
              nr=sspec(is)%nr
              allocate(g_rv(nr*2))
              allocate(rofi_rv(nr*2))
              call radmsh ( rmt , a , nr , rofi_rv )
              allocate(v0i_rv(nr))
              call dpscop ( v0pot(ib)%v , v0i_rv , nr , 1 + nr * ( isp - 1 ) , 1 , 1d0 )
              if (mod(idmod(m),10) .eq. 3) then
                 val(1) = rmt
                 dl = dtan(pi*(0.5d0-mod(pnu(m,isp),10d0)))
                 slo(1) = dl + 1d0
                 nn = int(mod(pnu(m,1),10d0))-l-1
                 allocate(gp_rv(8*nr))
                 call phidx ( 0 , z , l , v0i_rv , 0d0 , 0d0 , rofi_rv , &
                      nr , 2 , 1d-12 , ebar , val , slo , nn , g_rv , gp_rv &
                      , phi , dphi , phip , dphip , xx , xx , [xx] , xx , [xx] )
                 !         ... cz = estimate for energy of orbital with b.c. connecting
                 !             to Hankel of energy 0
                 dlphi = rmt*dphi/phi
                 dlphip = rmt*dphip/phip
                 umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
                 cz = ebar + umegam
                 ebar = cz
              endif
              if(lsemicorepz)then
                 ebar=eferm
                 if(ipr>40) write(6,"(' pnunew: valence with semicore ebar=efermi=',f12.6)")ebar
              endif
              call phidx ( 2 , z , l , v0i_rv , 0d0 , 0d0 , rofi_rv , nr ,  &
                   0 , 1d-12 , ebar , val , slo , nn , g_rv , gp_rv , phi , &
                   dphi , 0d0 , 0d0 , 0d0 , 0d0 , [0d0] , 0d0 , [0d0] )
              nnv=nn
              if (nn .eq. int(pnu(m,1))-l-1) then
                 dl = rmt*slo(1)/val(1) - 1
                 p1 = 0.5d0 - datan(dl)/pi
              elseif (ipr>=10.and.(.not.lsemicorepz)) then
                 write(stdo,ftox)'node=',nn, &
                 '(warning,probably no problem) not expecting node count for l=',l,'ebar=',ftof(ebar)
              endif
              if(lsemicorepz) then
                 dl = rmt*slo(1)/val(1) - 1
                 p1 = 0.5d0 - datan(dl)/pi
              endif
              !       ... Estimate new pnz for semicore state
              if (lpz .and. int(mod(pnz(m,1),10d0)).lt.int(pnu(m,1))) then !local obital is semicore
                 val(1) = rmt
                 dnz = dtan(pi*(0.5d0-mod(pnz(m,isp),10d0)))
                 slo(1) = dnz + 1d0
                 nn = int(mod(pnz(m,1),10d0))-l-1
                 allocate(gp_rv(8*nr))
                 call phidx ( 0 , z , l , v0i_rv , 0d0 , 0d0 , rofi_rv , &
                      nr , 2 , 1d-12 , ez , val , slo , nn , g_rv , gp_rv , &
                      phi , dphi , phip , dphip , xx , xx , [xx] , xx , [xx] )
                 dlphi = rmt*dphi/phi
                 dlphip = rmt*dphip/phip
                 umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
                 !         ... cz = estimate for energy of orbital with b.c. connecting
                 !             to Hankel of energy 0 (=> solution for interstitial
                 !             is constant potential, value C
                 cz = ez + umegam
                 !         ... estimate for log der. of wf for constant pot of value C
                 !             Maybe better to recalc. val,slo for this cz like above?
                 call bessl(cz*rmt**2,m,fi,gi)
                 !             dh_l/dr = l*h_l/r - h_l+1, h=g_l/r**(l+1)
                 !             when cz -> 0,  dl -> -l-1
                 dl = (l*gi(l) - gi(l+1))/gi(l)
                 pznew = 0.5d0 - datan(dl)/pi
                 !                lfrzv = .true. !takao commentout Sep2010
              else
                 lfrzz = .true.
                 pznew = 0.5d0 - datan(dble(l))/pi
                 ez = 0
              endif
              if (allocated(gp_rv)) deallocate(gp_rv)
              if (allocated(v0i_rv)) deallocate(v0i_rv)
              if (allocated(rofi_rv)) deallocate(rofi_rv)
              if (allocated(g_rv)) deallocate(g_rv)
           endif
           !     ... Free-electron value for pnu
           ipqn = pnu(m,isp)
           d0l = l
           pfree = ipqn + 0.5d0 - datan(d0l)/pi
           !     --- Set the new pnu ---
           pold = pnu(m,isp)
           ipqn = pold
           ptry = pold
           if (dabs(p1) .lt. 1d10) ptry = nnv + l+1 +p1 !takao
           if (.not. lfrzv) then
              pnu(m,isp) = ptry
              !       ... Permit pnu no lower than free electron value or pmin
              if (ptry .lt. pfree) pnu(m,isp) = pfree
              if (pmin(m) .gt. 0 .and. pmin(m) .lt. 1) then
                 if (ptry .lt. ipqn+pmin(m)) pnu(m,isp) = ipqn+pmin(m)
              endif
              !       ... Permit pnu no higher than pmax
              if (pmax(m) .gt. 0 .and. pmax(m) .lt. 1) then
                 if (ptry .gt. ipqn+pmax(m)) pnu(m,isp) = ipqn+pmax(m)
              endif
           endif
           if (ipr>40) write(stdo,310) &
                l,isp,idmod(m),qbyl(m,isp,ib),ebar,pold,ptry,pfree,pnu(m,isp)
!           if (ipr>40 .and. isp .eq. 2) write(stdo,410) &
!                idmod(m),qbyl(m,isp,ib),ebar,pold,ptry,pfree,pnu(m,isp)
310        format(i2,i2,i6,6f12.6,l)
!410        format(' spn 2',i2,6f12.6,l)
311        format(' l isp idmod     ql',9x,'ebar',7x,' pold',8x,'ptry',8x,'pfree',8x,'pnew',8x)
           !     --- Set the new pnz ---
           if (lpz) then
              pold = mod(pnz(m,isp),10d0)
              ipqn = pold
              ptry = pold
              pfree = ipqn + 0.5d0 - datan(d0l)/pi
              if (dabs(pznew) .lt. 1d10) ptry = ipqn+pznew
              if (.not. lfrzz) then
                 pnz(m,isp) = ptry + (pnz(m,isp)-mod(pnz(m,isp),10d0))
                 !         ... Permit pnu no lower than free electron value
                 d0l = l
                 if (ptry .lt. pfree) pnz(m,isp) = pfree + (pnz(m,isp)-mod(pnz(m,isp),10d0))
              endif
              if (ipr>40) write(stdo,520)l,isp,idmod(m),ez,pold,ptry,pfree,pnz(m,isp)
!              if (ipr>40.and. isp .eq. 2) write(stdo,520)l,isp,idmod(m),ez,pold,ptry,pfree,pnz(m,isp)
520           format(i2,i2,i6,'      ---   ',6f12.6)
!620           format(' spn 2',i2,'      ---   ',6f12.6)
           elseif (lpz) then
           endif
        enddo !end of spin loop
        phispinsym= cmdopt0('--phispinsym') !! spin averaged pnu takaoAug2019
        if(phispinsym) then
           if(ipr>0.and.m==lmxa+1) write(6,*)'pnunew: --phispinsym enforces spin-averaged pnu' 
           pmean = sum(pnu(m,1:nsp))/nsp
           pnu(m,1:nsp) = pmean
           if (lpz) then
              pmean = sum(pnz(m,1:nsp))/nsp
              pnz(m,1:nsp) = pmean
           endif
        endif
     enddo 
     !ssite(ib)%pnu=pnu
     !ssite(ib)%pz=pnz
     pnuall(:,1:nsp,ib)=pnu(:,1:nsp)
     pnzall(:,1:nsp,ib)=pnz(:,1:nsp)
10   continue
  enddo
  call tcx('pnunew')
end subroutine pnunew
