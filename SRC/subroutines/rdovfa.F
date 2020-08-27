      subroutine rdovfa ( nbas , nspec , ssite , sspec  !, spot !, slat
     . , qbg , sv_p_orhoat )
      use m_work,only: rv_a_ogv,iv_a_okv,rv_a_ogv,zv_a_osmrho
      use m_struc_def  !Cgetarg
      use m_struc_func , only: mpibc1_s_spec
      use m_globalvariables
      use m_lmfinit,only:lat_alat,lat_plat,lat_nabc,lat_vol,lat_ng
C- Read and overlap free atom densities.
C  Also allocates orhoca with free-atom core density.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   slat  :struct containing information about the lattice
Ci   qbg   :constant background charge
Co Outputs
Co   orhoat:vector of offsets containing site density, in standard
Co          3-component form (true rho, smoothed rho, core rho)
Co   smrho :smoothed interstitial density
Co         :* for smrho = smoothed mesh density, smrho is complex and
Co         :  smrho = smrho(k1,k2,k3)
Cl Local variables
Ci   k1,k2,k3 dimensions of smrho for smoothed mesh density
Cr Remarks
Cu Updates
Cu   12 May 07 package mpi-specific calls
Cu   02 Jan 06 generates core magnetic moment, checks against spec->qcor
Cu   01 Jul 05 Zero-radius empty spheres treated as local orbitals
Cu   12 Apr 03 (WRL) Added constant charge background
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   13 Jun 00 spin polarized
Cu   21 Apr 00 Adapted from nfp rdovfa.f
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer procid, master, mpipid
      integer:: nbas , nspec
      type(s_rv1) :: sv_p_orhoat(3,nbas)
      real(8):: qbg
      type(s_site)::ssite(nbas)
      type(s_spec)::sspec(nspec)
c      type(s_lat)::slat
c      type(s_pot)::spot

C ... Local parameters
      integer nrmx, n0
      parameter ( nrmx=1501, n0=10 )
      integer:: nxi(nspec)
      type(s_rv1) :: rv_a_orhofa(nspec)
      type(s_rv1) :: rv_a_ov0a(nspec)

      double precision rsmfa(nspec),pnu(n0,2),exi(n0,nspec),
     .hfc(n0,2,nspec),hfct(n0,2,nspec)
      character*8 spid(nspec),spidr
      integer:: ipr , iprint , ngabc(3) , n1 , n2 , n3 , k1 , k2 , 
     .     k3 , nsp , nglob , lgunit , stdo , stdl , i , ifi , fopna , is 
     .     , nr , lfoc , nr0 , i1 , nch , ib , igetss , lmxl , nlml , ng 
     .     , iofa , kcor , lcor
      real(8) ,allocatable :: rwgt_rv(:)
      complex(8) ,allocatable :: cv_zv(:)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,plat(3,3),a,rmt,z,rfoc,z0,rmt0,a0,qc,ccof,
     .ceh,stc,ztot,ctot,corm,sum,fac,sum1,sum2,sqloc,dq,vol,smom,
     .slmom,qcor(2)
      character msg*23, strn*120
      logical mlog,cmdopt,lfail
      logical:: l_dummy_isanrg,isanrg
      integer:: i_copy_size,i_spec
      call tcn('rdovfa')
      ipr   = iprint()
      stdo  = lgunit(1)
      stdl  = lgunit(2)
      nsp   = globalvariables%nsp
      msg   = '         File mismatch:'
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      if (ipr .ge. 10) write(stdo,700)
  700 format(/' rdovfa: read and overlap free-atom densities',' (mesh density) ...')
      alat=lat_alat
c      i_copy_size=size(lat_plat)
c      call dcopy(i_copy_size,lat_plat,1,plat,1)
c      i_copy_size=size(lat_nabc)
c     call icopy(i_copy_size,lat_nabc,1,ngabc,1)
      plat=lat_plat
      ngabc=lat_nabc
      vol=lat_vol
      call fftz30(n1,n2,n3,k1,k2,k3)
      call dpzero(hfc,n0*2*nspec)
      call dpzero(pnu,n0*2)

C  --- Read free-atom density for all species ---
C      if (ipr .ge. 30) write(stdo,*) ' '
      if (procid .eq. master) then
        ifi = fopna('atm',-1,0)
        rewind ifi
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
        if (procid .eq. master) then
          if (z .eq. 0 .and. rmt .eq. 0) then
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
            deallocate(sspec(is)%rv_a_orhoc)
            if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
          else
            nr0=nrmx ! %rv_a_orhoc(nrmx*nsp)
            lfail = ( iofa ( spidr , n0 , nxi ( is ) , exi ( 1 , is ) , hfc 
     .        ( 1 , 1 , is ) , hfct ( 1 , 1 , is ) , rsmfa ( is ) , z0 , rmt0 
     .        , a0 , nr0 , qc , ccof , ceh , stc , rv_a_orhofa( is )%v , sspec
     .        ( is ) %rv_a_orhoc , rv_a_ov0a ( is ) %v , ifi ) .lt. 0 )
          endif
        endif
        call mpibc1(nr0,1,2,mlog,'rdovfa','nr0')
        call mpibc1(lfail,1,1,mlog,'rdovfa','read error')
        if (lfail) call rxs('missing species data, species ',spid(is))
C       i = mpipid(3)
C   ... Broadcast file data
        call mpibc1(nxi(is),1,2,mlog,'rdovfa','nxi')
        call mpibc1(exi(1,is),nxi(is),4,mlog,'rdovfa','exi')
        call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'rdovfa','hfc')
        call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'rdovfa','hfct')
        call mpibc1(rsmfa(is),1,4,mlog,'rdovfa','rsmfa')
        call mpibc1(a0,1,4,mlog,'rdovfa','a0')
        call mpibc1 ( rv_a_orhofa( is )%v , nr0 * nsp , 4 , mlog , 'rdovfa' 
     .   , 'rhofa' )
ckino Dec.28.2011:          call mpibc1(sspec(is)%rv_p_orhoc,nr0*nsp,4,mlog,'rdovfa','rhoca')
        call mpibc1(sspec(is)%rv_a_orhoc,nr0*nsp,4,mlog,'rdovfa','rhoca')
        call mpibc1 ( rv_a_ov0a( is )%v , nr0 * nsp , 4 , mlog , 'rdovfa', 'v0a' )
        i = mpipid(3)
C ...   Defaults
        if (procid .eq. master) then
          call strip(spid(is),i1,nch)
          if (ipr .ge. 30 .and. rmt0 .ne. 0)
     .    write(stdo,400) spid(is)(1:nch),spidr,rmt0,nr0,a0
  400     format(' rdovfa: expected ',a,',',T27,' read ',a,
     .    ' with rmt=',f8.4,'  mesh',i6,f7.3)
        endif
        if (nr .le. 0)   nr = nr0
        if (a .le. 1d-6) a = a0
        if (z .eq. 0 .and. rmt .eq. 0) then
          a = 0
          nr = 0
        endif

C ...   Sanity checks
        if (procid .eq. master) then
          call fsanrg(z0,z,z,0d-9,msg,'z',.true.)
          call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.true.)
          call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
Ckino isanrg is logical function,           call isanrg(nr0,nr,nr,msg,'nr',.true.)
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
   10 continue
C     Wait for all proccesses to synchronize
      i = mpipid(3)
C     Re-broadcast entire species structure, and arrays used below
cki      call mpibc1(sspec,nspec*nint(sspec(1)%size),4,0,'rdovfa','sspec')
      do i_spec=1,nspec
        call mpibc1_s_spec(sspec(i_spec),.false.,'rdovfa','sspec')
      enddo

      if (procid .eq. master) then
        call fclose(ifi)
      endif

C --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
      ztot = 0d0
      ctot = 0d0
      corm = 0d0
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

        allocate(sv_p_orhoat(1,ib)%v(abs(nr*nlml*nsp)))
        if (nr*nlml*nsp<0) sv_p_orhoat(1,ib)%v(:)=0.0d0
        allocate(sv_p_orhoat(2,ib)%v(abs(nr*nlml*nsp)))
        if (nr*nlml*nsp<0) sv_p_orhoat(2,ib)%v(:)=0.0d0
        allocate(sv_p_orhoat(3,ib)%v(abs(nr*nsp)))
        if (nr*nsp<0) sv_p_orhoat(3,ib)%v(:)=0.0d0
        if (allocated(ssite(ib)%rv_a_ov0)) deallocate(ssite(ib)%rv_a_ov0)
        allocate(ssite(ib)%rv_a_ov0(abs(nr*nsp)))
        if (nr*nsp<0) ssite(ib)%rv_a_ov0(:)=0.0d0
        if (allocated(ssite(ib)%rv_a_ov1)) deallocate(ssite(ib)%rv_a_ov1)
        allocate(ssite(ib)%rv_a_ov1(abs(nr*nsp)))
c        if (nr*nsp<0) ssite(ib)%rv_a_ov1(:)=0.0d0

C       Core magnetic moment (possible if magnetized core hole)
        if (nsp .eq. 2 .and. lmxl .gt. -1) then
          allocate(rwgt_rv(nr))
          call radwgt ( rmt , a , nr , rwgt_rv )
          call radsum ( nr , nr , 1 , nsp , rwgt_rv , sspec(is)%rv_a_orhoc , sum )
          call radsum ( nr , nr , 1 , 1 , rwgt_rv , sspec(is)%rv_a_orhoc , sum1 )
          sum2 = sum - sum1
          call gtpcor(sspec,is,kcor,lcor,qcor)
          if (dabs(qcor(2)-(sum1-sum2)) .gt. 0.01d0) then
            call info5(10,0,0,' (warning) core moment mismatch spec %i:'
     .      //'  input file=%;6d  atom file=%;6d', is,qcor(2),sum1-sum2,0,0)
          endif
          corm = corm + qcor(2)
          if (allocated(rwgt_rv)) deallocate(rwgt_rv)
        endif
        if (lmxl .gt. -1) then
          call dpcopy ( rv_a_ov0a( is )%v , ssite(ib)%rv_a_ov0 , 1 , nr * nsp , 1d0  )
          call dpcopy ( rv_a_ov0a( is )%v , ssite(ib)%rv_a_ov1 , 1 , nr * nsp , 1d0  )
          call dpcopy ( sspec ( is ) %rv_a_orhoc , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , 1d0 )
          if (lfoc .eq. 0) then
            allocate(rwgt_rv(nr))
            call radwgt ( rmt , a , nr , rwgt_rv )
            call radsum ( nr , nr , 1 , nsp , rwgt_rv , sv_p_orhoat( 3 , ib )%v , sum )
            fac = 1d0
            if(dabs(sum) .gt. 1d-7) fac = qc/sum
            if (ipr .ge. 40) write(stdo,787) is,qc,sum,fac
  787       format(' scale foca=0 core species',i2,': qc,sum,scale=',
     .      3f12.6,f12.6)
            call dpcopy ( sv_p_orhoat( 3 , ib )%v , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , fac )
            if (allocated(rwgt_rv)) deallocate(rwgt_rv)
          endif
        endif
        ztot = ztot+z
        ctot = ctot+qc
C     end loop over sites
   20 continue

C --- Overlap smooth hankels to get smooth interstitial density ---
      ng=lat_ng
      allocate(cv_zv(ng*nsp))
      call ovlpfa ( ssite , nbas , nxi , n0 , exi , hfc , rsmfa !slat , 
     . , ng , ng , rv_a_ogv , cv_zv )
      call gvputf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , cv_zv , zv_a_osmrho )
      if (allocated(cv_zv)) deallocate(cv_zv)
C ... FFT to real-space mesh
      call fftz3 ( zv_a_osmrho , n1 , n2 , n3 , k1 , k2 , k3 , nsp, 0 , 1 )
C ... Add compensating uniform electron density to compensate background
      call addbkgsm ( zv_a_osmrho , k1 , k2 , k3 , nsp , qbg , vol, - 1d0 )
C ... integrate
      call mshint ( vol , nsp , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho, sum1 , sum2 )
      if (nsp .eq. 2) then
        call mshint ( vol , 1 , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho, smom , sum2 )
        smom = 2*smom - sum1
      endif
C --- Set up local densities using rmt from atm file ---
      call ovlocr ( nbas , ssite , sspec ,  n0 , nxi , exi , !slat ,
     .  hfc , rsmfa , rv_a_orhofa , sv_p_orhoat , sqloc , slmom )
C --- Add compensating uniform electron density to compensate background
      call adbkql ( sv_p_orhoat , nbas , nsp , qbg , vol , - 1d0 , sspec , ssite )
      if (abs(qbg) .ne. 0) call info(10,0,0,' Uniform '//
     .'density added to neutralize background, q=%;6,6d',qbg,0)
C --- Print charges ---
      dq = sum1+sqloc+ctot-ztot+qbg
      if (nsp .eq. 1) then
        if (ipr .ge. 10)
     .  write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
  895   format(/' Smooth charge on mesh:    ',f16.6
     .  /    ' Sum of local charges:     ',f16.6
     .  /    ' Total valence charge:     ',f16.6
     .  /    ' Sum of core charges:      ',f16.6
     .  /    ' Sum of nuclear charges:   ',f16.6
     .  /    ' Homogeneous background:   ',f16.6
     .  /    ' Deviation from neutrality:',f16.6)
        if (ipr .ge. 10) write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
  710   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .  '   bg',f10.6,'  dQ',f10.6)
      else
        if (ipr .ge. 10) write(stdo,896) sum1,smom,sqloc,slmom,
     .  sum1+sqloc,smom+slmom,ctot,corm,-ztot,qbg,dq
  896   format(/' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6,
     .  /    ' Sum of local charges:     ',f16.6,4x,'moments',f11.6,
     .  /    ' Total valence charge:     ',f16.6,4x,'moment', f12.6,
     .  /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6,
     .  /    ' Sum of nuclear charges:   ',f16.6
     .  /    ' Homogeneous background:   ',f16.6
     .  /    ' Deviation from neutrality:',f16.6)
        if (ipr .ge. 10)
     .  write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom+slmom
  711   format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,
     .  '   bg',f11.6,' mm',f11.6)
      endif
      if (dabs(dq) .gt. 1d-4 .and. ipr .gt. 0)
     .call awrit1(' rdovfa (warning) overlapped'
     .//' density not neutral'//', dq= %d',' ',80,stdo,dq)
      do is=1,nspec
        if (allocated(rv_a_ov0a(is)%v)) deallocate(rv_a_ov0a(is)%v)
        if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
      enddo
      call tcx('rdovfa')
      end subroutine rdovfa


