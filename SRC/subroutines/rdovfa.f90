subroutine rdovfa()
  use m_density,only: zv_a_osmrho=>osmrho,sv_p_orhoat=>orhoat

  use m_supot,only: lat_nabc,lat_ng,rv_a_ogv,iv_a_okv,rv_a_ogv
  use m_lmfinit,only:lat_alat,nsp,nbas,nspec,ssite=>v_ssite,sspec=>v_sspec,qbg=>zbak
  use m_lattic,only: lat_plat,lat_vol
  use m_struc_def,only: s_rv1
  use m_struc_func, only: mpibc1_s_spec
  use m_ext,only: sname
  use m_lgunit,only:stdo,stdl
  !!- Read and overlap free atom densities.
  !  allocates orhoca with free-atom core density.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   nspec :number of species
  !i   ssite :struct containing site-specific information
  !i   sspec :struct containing species-specific information
  !i   qbg  :constant background charge  qcore+qval-qbg = \sum_i Zi
  !o Outputs
  !o   orhoat: vector of offsets containing site density, in standard
  !o           3-component form (true rho, smoothed rho, core rho)
  !o   smrho :smoothed interstitial density
  !o         :* for smrho = smoothed mesh density, smrho is complex and
  !o         :  smrho = smrho(k1,k2,k3)
  ! ----------------------------------------------------------------------
  implicit none
  integer :: procid, master, mpipid
  integer :: nrmx, n0
  parameter ( nrmx=1501, n0=10 )
  integer:: nxi(nspec)
  type(s_rv1) :: rv_a_orhofa(nspec)
  type(s_rv1) :: rv_a_ov0a(nspec)
  double precision :: rsmfa(nspec),pnu(n0,2),exi(n0,nspec), &
       hfc(n0,2,nspec),hfct(n0,2,nspec)
  character(8) :: spid(nspec),spidr
  integer:: ipr , iprint , ngabc(3) , n1 , n2 , n3 , k1 , k2 , &
       k3 , i , ifi , is &
       , nr , lfoc , nr0 , i1 , nch , ib , igetss , lmxl , nlml , ng &
       , iofa , kcor , lcor
  real(8) ,allocatable :: rwgt_rv(:)
  complex(8) ,allocatable :: cv_zv(:)
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  double precision :: alat,plat(3,3),a,rmt,z,rfoc,z0,rmt0,a0,qc,ccof, &
       ceh,stc,ztot,ctot,corm,sum,fac,sum1,sum2,sqloc,dq,vol,smom, &
       slmom,qcor(2)
  character msg*23, strn*120
  logical :: mlog,cmdopt,lfail
  logical:: l_dummy_isanrg,isanrg
  integer:: i_copy_size,i_spec,ifile_handle
  call tcn('rdovfa')
  ipr   = iprint()
  !      stdo  = lgunit(1)
  !      stdl  = lgunit(2)
  msg   = '         File mismatch:'
  procid = mpipid(1)
  master = 0
  mlog = cmdopt('--mlog',6,0,strn)
  if (ipr >= 10) write(stdo,700)
700 format(/' rdovfa: read and overlap free-atom densities',' (mesh density) ...')
  alat=lat_alat
  plat=lat_plat
  ngabc=lat_nabc
  vol=lat_vol
  call fftz30(n1,n2,n3,k1,k2,k3)
  call dpzero(hfc,n0*2*nspec)
  call dpzero(pnu,n0*2)
  !! Read free-atom density for all species ---
  if (procid == master) then
     ifi=ifile_handle()
     open(ifi,file='atm.'//trim(sname))
  endif
  do  10  is = 1, nspec
     allocate(rv_a_orhofa(is)%v(abs(-nrmx*nsp)))
     rv_a_orhofa(is)%v(:)=0.0d0
     if (allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
     allocate(sspec(is)%rv_a_orhoc(   nrmx*nsp) )
     sspec(is)%rv_a_orhoc=0.0d0
     allocate(rv_a_ov0a(is)%v(nrmx*nsp))
     rv_a_ov0a(is)%v(:)=0.0d0
     spid(is)=sspec(is)%name
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     z=sspec(is)%z
     lfoc=sspec(is)%lfoca
     rfoc=sspec(is)%rfoca
     lfail = .false.
     if (procid == master) then
        if (z == 0 .AND. rmt == 0) then
           nxi(is) = 0
           call dpzero(exi(1,is),n0)
           call dpzero(hfc(1,1,is),2*n0)
           call dpzero(hfct(1,1,is),2*n0)
           rsmfa(is) = 0
           z0=0
           rmt0=0
           a0=0
           nr0=0
           qc=0
           ccof=0
           ceh=0
           stc=0
           if (allocated(rv_a_ov0a(is)%v)) deallocate(rv_a_ov0a(is)%v)
           !            deallocate(sspec(is)%rv_a_orhoc)
           if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
        else
           nr0=nrmx ! %rv_a_orhoc(nrmx*nsp)
           lfail = ( iofa ( spidr , n0 , nxi ( is ) , exi ( 1 , is ) , hfc &
                ( 1 , 1 , is ) , hfct ( 1 , 1 , is ) , rsmfa ( is ) , z0 , rmt0 &
                , a0 , nr0 , qc , ccof , ceh , stc , rv_a_orhofa( is )%v , sspec &
                ( is ) %rv_a_orhoc , rv_a_ov0a ( is ) %v , ifi ) .lt. 0 )
        endif
     endif
     call mpibc1(nr0,1,2,mlog,'rdovfa','nr0')
     call mpibc1(lfail,1,1,mlog,'rdovfa','read error')
     if (lfail) call rxs('missing species data, species ',spid(is))
     !       i = mpipid(3)
     !   ... Broadcast file data
     call mpibc1(nxi(is),1,2,mlog,'rdovfa','nxi')
     call mpibc1(exi(1,is),nxi(is),4,mlog,'rdovfa','exi')
     call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'rdovfa','hfc')
     call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'rdovfa','hfct')
     call mpibc1(rsmfa(is),1,4,mlog,'rdovfa','rsmfa')
     call mpibc1(a0,1,4,mlog,'rdovfa','a0')
     call mpibc1 ( rv_a_orhofa( is )%v , nr0 * nsp , 4 , mlog , 'rdovfa' &
          , 'rhofa' )
     ! ino Dec.28.2011:          call mpibc1(sspec(is)%rv_p_orhoc,nr0*nsp,4,mlog,'rdovfa','rhoca')
     call mpibc1(sspec(is)%rv_a_orhoc,nr0*nsp,4,mlog,'rdovfa','rhoca')
     call mpibc1 ( rv_a_ov0a( is )%v , nr0 * nsp , 4 , mlog , 'rdovfa', 'v0a' )
     i = mpipid(3)
     ! ...   Defaults
     if (procid == master) then
        call strip(spid(is),i1,nch)
        if (ipr >= 30 .AND. rmt0 /= 0) &
             write(stdo,400) spid(is)(1:nch),spidr,rmt0,nr0,a0
400     format(' rdovfa: expected ',a,',',T27,' read ',a, &
             ' with rmt=',f8.4,'  mesh',i6,f7.3)
     endif
     if (nr <= 0)   nr = nr0
     if (a <= 1d-6) a = a0
     if (z == 0 .AND. rmt == 0) then
        a = 0
        nr = 0
     endif

     ! ...   Sanity checks
     if (procid == master) then
        call fsanrg(z0,z,z,0d-9,msg,'z',.true.)
        call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.true.)
        call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
        ! ino isanrg is logical function,           call isanrg(nr0,nr,nr,msg,'nr',.true.)
        l_dummy_isanrg=isanrg(nr0,nr,nr,msg,'nr',.true.)
     endif
     sspec(is)%a=a
     sspec(is)%nr=nr
     sspec(is)%qc=qc
     sspec(is)%nxi=nxi(is)
     i_copy_size=size(sspec(is)%exi)
     call dcopy(i_copy_size,exi(1,is),1,sspec(is)%exi,1)
     i_copy_size=size(sspec(is)%chfa)
     call dcopy(i_copy_size,hfc(1,1,is),1,sspec(is)%chfa,1)
     sspec(is)%rsmfa=rsmfa(is)
     sspec(is)%ctail=ccof
     sspec(is)%etail=ceh
     sspec(is)%stc=stc
10 enddo
  !     Wait for all proccesses to synchronize
  i = mpipid(3)
  !     Re-broadcast entire species structure, and arrays used below
  ! i      call mpibc1(sspec,nspec*nint(sspec(1)%size),4,0,'rdovfa','sspec')
  do i_spec=1,nspec
     call mpibc1_s_spec(sspec(i_spec),'rdovfa_sspec')
  enddo

  if (procid == master) close(ifi)

  ! --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
  ztot = 0d0
  ctot = 0d0
  corm = 0d0
  allocate(sv_p_orhoat(3,nbas))
  do  20  ib = 1, nbas
     is = int(ssite(ib)%spec)
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     lmxl=sspec(is)%lmxl
     z=sspec(is)%z
     qc=sspec(is)%qc
     lfoc=sspec(is)%lfoca
     nlml = (lmxl+1)**2
     !!
     allocate(sv_p_orhoat(1,ib)%v(nr*nlml*nsp))
     allocate(sv_p_orhoat(2,ib)%v(nr*nlml*nsp))
     allocate(sv_p_orhoat(3,ib)%v(nr*nsp))
     if (allocated(ssite(ib)%rv_a_ov0)) deallocate(ssite(ib)%rv_a_ov0)
     allocate(ssite(ib)%rv_a_ov0(abs(nr*nsp)))
     if (allocated(ssite(ib)%rv_a_ov1)) deallocate(ssite(ib)%rv_a_ov1)
     allocate(ssite(ib)%rv_a_ov1(abs(nr*nsp)))

     !       Core magnetic moment (possible if magnetized core hole)
     if (nsp == 2 .AND. lmxl > -1) then
        allocate(rwgt_rv(nr))
        call radwgt ( rmt , a , nr , rwgt_rv )
        call radsum ( nr , nr , 1 , nsp , rwgt_rv , sspec(is)%rv_a_orhoc , sum )
        call radsum ( nr , nr , 1 , 1 , rwgt_rv , sspec(is)%rv_a_orhoc , sum1 )
        sum2 = sum - sum1
        call gtpcor(sspec,is,kcor,lcor,qcor)
        if (dabs(qcor(2)-(sum1-sum2)) > 0.01d0) then
           call info5(10,0,0,' (warning) core moment mismatch spec %i:' &
                //'  input file=%;6d  atom file=%;6d', is,qcor(2),sum1-sum2,0,0)
        endif
        corm = corm + qcor(2)
        if (allocated(rwgt_rv)) deallocate(rwgt_rv)
     endif
     if (lmxl > -1) then
        call dpcopy ( rv_a_ov0a( is )%v , ssite(ib)%rv_a_ov0 , 1 , nr * nsp , 1d0  )
        call dpcopy ( rv_a_ov0a( is )%v , ssite(ib)%rv_a_ov1 , 1 , nr * nsp , 1d0  )
        call dpcopy ( sspec ( is ) %rv_a_orhoc , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , 1d0 )
        if (lfoc == 0) then
           allocate(rwgt_rv(nr))
           call radwgt ( rmt , a , nr , rwgt_rv )
           call radsum ( nr , nr , 1 , nsp , rwgt_rv , sv_p_orhoat( 3 , ib )%v , sum )
           fac = 1d0
           if(dabs(sum) > 1d-7) fac = qc/sum
           if (ipr >= 40) write(stdo,787) is,qc,sum,fac
787        format(' scale foca=0 core species',i2,': qc,sum,scale=', &
                3f12.6,f12.6)
           call dpcopy ( sv_p_orhoat( 3 , ib )%v , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , fac )
           if (allocated(rwgt_rv)) deallocate(rwgt_rv)
        endif
     endif
     ztot = ztot+z
     ctot = ctot+qc
     !     end loop over sites
20 enddo

  !! allocate array for iteration
  allocate(zv_a_osmrho(k1*k2*k3*nsp))
  zv_a_osmrho(:)=0d0

  ! --- Overlap smooth hankels to get smooth interstitial density ---
  ng=lat_ng
  allocate(cv_zv(ng*nsp))
  call ovlpfa ( ssite , nbas , nxi , n0 , exi , hfc , rsmfa, ng , ng , rv_a_ogv , cv_zv )
  call gvputf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , cv_zv , zv_a_osmrho )
  if (allocated(cv_zv)) deallocate(cv_zv)
  ! ... FFT to real-space mesh
  call fftz3 ( zv_a_osmrho , n1 , n2 , n3 , k1 , k2 , k3 , nsp, 0 , 1 )
  ! ... Add compensating uniform electron density to compensate background
  call addbkgsm ( zv_a_osmrho , k1 , k2 , k3 , nsp , qbg , vol, - 1d0 )
  ! ... integrate
  call mshint ( vol , nsp , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho, sum1 , sum2 )
  if (nsp == 2) then
     call mshint ( vol , 1 , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho, smom , sum2 )
     smom = 2*smom - sum1
  endif
  ! --- Set up local densities using rmt from atm file ---
  call ovlocr ( nbas , ssite , sspec ,  n0 , nxi , exi ,&
       hfc , rsmfa , rv_a_orhofa , sv_p_orhoat , sqloc , slmom )
  ! --- Add compensating uniform electron density to compensate background
  call adbkql ( sv_p_orhoat , nbas , nsp , qbg , vol , - 1d0 , sspec , ssite )
  if (abs(qbg) /= 0) call info(10,0,0,' Uniform '// &
       'density added to neutralize background, q=%;6,6d',qbg,0)
  ! --- Print charges ---
  dq = sum1+sqloc+ctot-ztot+qbg
  if (nsp == 1) then
     if (ipr >= 10) &
          write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
895  format(/' Smooth charge on mesh:    ',f16.6 &
          /    ' Sum of local charges:     ',f16.6 &
          /    ' Total valence charge:     ',f16.6 &
          /    ' Sum of core charges:      ',f16.6 &
          /    ' Sum of nuclear charges:   ',f16.6 &
          /    ' Homogeneous background:   ',f16.6 &
          /    ' Deviation from neutrality:',f16.6)
     if (ipr >= 10) write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
710  format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6, &
          '   bg',f10.6,'  dQ',f10.6)
  else
     if (ipr >= 10) write(stdo,896) sum1,smom,sqloc,slmom, &
          sum1+sqloc,smom+slmom,ctot,corm,-ztot,qbg,dq
896  format(/' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6, &
          /    ' Sum of local charges:     ',f16.6,4x,'moments',f11.6, &
          /    ' Total valence charge:     ',f16.6,4x,'moment', f12.6, &
          /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6, &
          /    ' Sum of nuclear charges:   ',f16.6 &
          /    ' Homogeneous background:   ',f16.6 &
          /    ' Deviation from neutrality:',f16.6)
     if (ipr >= 10) &
          write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom+slmom
711  format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6, &
          '   bg',f11.6,' mm',f11.6)
  endif
  if (dabs(dq) > 1d-4 .AND. ipr > 0) &
       write(stdo,"(' rdovfa (warning) overlapped' &
       //' density not neutral'//', dq= ',d13.5)") dq
  do is=1,nspec
     if (allocated(rv_a_ov0a(is)%v)) deallocate(rv_a_ov0a(is)%v)
     if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
  enddo
  call tcx('rdovfa')
end subroutine rdovfa

