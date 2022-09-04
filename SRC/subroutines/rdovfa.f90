module m_rdovfa
  public rdovfa
contains
  subroutine rdovfa() !- Read atm files and overlap free atom densities.
    use m_density,only: zv_a_osmrho=>osmrho,sv_p_orhoat=>orhoat,v1pot,v0pot !Outputs. allocated

    use m_supot,only: lat_nabc,lat_ng,rv_a_ogv,iv_a_okv,rv_a_ogv
    use m_lmfinit,only:lat_alat,nsp,nbas,nspec,ispec,sspec=>v_sspec,qbg=>zbak,slabl
    use m_lattic,only: lat_plat,lat_vol
    use m_struc_def,only: s_rv1
    use m_struc_func, only: mpibc1_s_spec
    use m_ext,only: sname
    use m_lgunit,only:stdo,stdl
    use m_ftox
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
    integer :: procid, master, mpipid, nrmx, n0,i_spec,ifile_handle
    parameter ( nrmx=1501, n0=10 )
    integer:: nxi(nspec)
    type(s_rv1) :: rv_a_orhofa(nspec)
    type(s_rv1) :: rv_a_ov0a(nspec)
    real(8):: rsmfa(nspec),exi(n0,nspec), hfc(n0,2,nspec),hfct(n0,2,nspec),&
         alat,plat(3,3),a,rmt,z,rfoc,z0,rmt0,a0,qc,ccof, &
         ceh,stc,ztot,ctot,corm,ssum,fac,sum1,sum2,sqloc,dq,vol,smom, slmom,qcor(2)
    character(8) :: spid(nspec),spidr
    integer:: ipr , iprint , ngabc(3) , n1 , n2 , n3 , k1 , k2 , iofa , kcor , lcor,&
         k3 , i , ifi , is, nr , lfoc , nr0 , i1 , nch , ib , igetss , lmxl , nlml , ng,ierr 
    real(8) ,allocatable :: rwgt_rv(:)
    complex(8) ,allocatable :: cv_zv(:)
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    character msg*23, strn*120
    logical :: mlog,cmdopt,lfail, l_dummy_isanrg,isanrg
    call tcn('rdovfa')
    ipr   = iprint()
    msg   = '         File mismatch:'
    procid = mpipid(1)
    master = 0
    mlog = cmdopt('--mlog',6,0,strn)
    if(ipr>=10)write(stdo,"(/'rdovfa: read and overlap free-atom densities',' (mesh density) ...')")
    alat=lat_alat
    plat=lat_plat
    ngabc=lat_nabc
    vol=lat_vol
    call fftz30(n1,n2,n3,k1,k2,k3)
    hfc=0d0
    exi=0d0
    hfc=0d0
    hfct=0d0
    if (procid == master) then
       ifi=ifile_handle()
       open(ifi,file='atm.'//trim(sname))  !! Read free-atom density for all species ---
    endif
    isloop: do  10  is = 1, nspec
       spid(is)=slabl(is)
       a=sspec(is)%a
       nr=sspec(is)%nr
       allocate(rv_a_orhofa(is)%v(nr*nsp))
       rv_a_orhofa(is)%v=0d0
       if(allocated(sspec(is)%rv_a_orhoc)) deallocate(sspec(is)%rv_a_orhoc)
       allocate(sspec(is)%rv_a_orhoc(nr*nsp) )
       sspec(is)%rv_a_orhoc=0d0
       allocate(rv_a_ov0a(is)%v(nr*nsp))
       rv_a_ov0a(is)%v(:)=0d0
       rmt=sspec(is)%rmt
       z=sspec(is)%z
       lfoc=sspec(is)%lfoca
       rfoc=sspec(is)%rfoca
       if (procid == master) then
          if (z == 0 .AND. rmt == 0) then
             nxi(is) = 0
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
             if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
          else
             nr0=nrmx 
             lfail = .false.
             lfail = iofa ( spidr , n0 , nxi ( is ) , exi ( 1 , is ) , hfc &
                  ( 1 , 1 , is ) , hfct ( 1 , 1 , is ) , rsmfa ( is ) , z0 , rmt0 &
                  , a0 , nr0 , qc , ccof , ceh , stc , rv_a_orhofa( is )%v , sspec &
                  ( is ) %rv_a_orhoc , rv_a_ov0a ( is ) %v , ifi )    < 0
             if (lfail) call rxs('missing species data, species ',spid(is))
          endif
       endif
       call mpibc1(nr0,1,2,mlog,'rdovfa','nr0')
       call mpibc1(nxi(is),1,2,mlog,'rdovfa','nxi')
       call mpibc1(exi(1,is),nxi(is),4,mlog,'rdovfa','exi')
       call mpibc1(hfc(1,1,is),nsp*n0,4,mlog,'rdovfa','hfc')
       call mpibc1(hfct(1,1,is),nsp*n0,4,mlog,'rdovfa','hfct')
       call mpibc1(rsmfa(is),1,4,mlog,'rdovfa','rsmfa')
       call mpibc1(a0,1,4,mlog,'rdovfa','a0')
       call mpibc1(rv_a_orhofa( is )%v , nr0 * nsp , 4 , mlog , 'rdovfa'  , 'rhofa' )
       call mpibc1(sspec(is)%rv_a_orhoc,nr0*nsp,4,mlog,'rdovfa','rhoca')
       call mpibc1( rv_a_ov0a( is )%v , nr0 * nsp , 4 , mlog , 'rdovfa', 'v0a' )
       i = mpipid(3)
       if (procid == master) then
          !call strip(spid(is),i1,nch)
          if (ipr >= 30 .AND. rmt0 /= 0) write(stdo,400) trim(spid(is)),spidr,rmt0,nr0,a0
400       format(' rdovfa: expected ',a,',',T27,' read ',a, ' with rmt=',f8.4,'  mesh',i6,f7.3)
       endif
       if (nr <= 0)   nr = nr0
       if (a <= 1d-6) a = a0
       if (z == 0 .AND. rmt == 0) then
          a = 0
          nr = 0
       endif
       if (procid == master) then
          print *,'zzzzzzzzz',z,z-z0
          call fsanrg(z0,z,z,0d-9,msg,'z',.true.)
          call fsanrg(rmt0,rmt,rmt,1d-6,msg,'rmt',.true.)
          call fsanrg(a0,a,a,0d-9,msg,'a',.true.)
          l_dummy_isanrg=isanrg(nr0,nr,nr,msg,'nr',.true.)
       endif
       sspec(is)%a=a
       sspec(is)%nr=nr
       sspec(is)%qc=qc
       sspec(is)%nxi=nxi(is)
       sspec(is)%exi=exi(:,is)
       sspec(is)%chfa=hfc(:,:,is)
       sspec(is)%rsmfa=rsmfa(is)
       sspec(is)%ctail=ccof
       sspec(is)%etail=ceh
       sspec(is)%stc=stc
10  enddo isloop
    i = mpipid(3)
    !     Re-broadcast entire species structure, and arrays used below
    do i_spec=1,nspec
       call mpibc1_s_spec(sspec(i_spec),'rdovfa_sspec')
    enddo
    if (procid == master) close(ifi)
    ! --- Define arrays for local densities rho1,rho2,rhoc and v0,v1 ---
    ztot = 0d0
    ctot = 0d0
    corm = 0d0
    if(allocated(v0pot)) deallocate(v0pot,v1pot) !this may cause mem leak?(v0pot%v is not deallocated).
    allocate(v1pot(nbas),v0pot(nbas))
    if(allocated(sv_p_orhoat)) deallocate(sv_p_orhoat)
    allocate(sv_p_orhoat(3,nbas))
    ibloop: do  20  ib = 1, nbas
       is = ispec(ib) !int(ssite(ib)%spec)
       a=sspec(is)%a
       nr=sspec(is)%nr
       rmt=sspec(is)%rmt
       lmxl=sspec(is)%lmxl
       z=sspec(is)%z
       qc=sspec(is)%qc
       lfoc=sspec(is)%lfoca
       nlml = (lmxl+1)**2
       allocate(sv_p_orhoat(1,ib)%v(nr*nlml*nsp))
       allocate(sv_p_orhoat(2,ib)%v(nr*nlml*nsp))
       allocate(sv_p_orhoat(3,ib)%v(nr*nsp))
       if (nsp == 2 .AND. lmxl > -1) then
          allocate(rwgt_rv(nr))
          call radwgt ( rmt , a , nr , rwgt_rv )
          call radsum ( nr , nr , 1 , nsp , rwgt_rv , sspec(is)%rv_a_orhoc , ssum )
          call radsum ( nr , nr , 1 , 1 , rwgt_rv , sspec(is)%rv_a_orhoc , sum1 )
          sum2 = ssum - sum1
          call gtpcor(is,kcor,lcor,qcor)
          if (dabs(qcor(2)-(sum1-sum2)) > 0.01d0) then
             if(ipr>=10) write(stdo,ftox)' (warning) core moment mismatch spec ',is, &
                  'input file=',ftof(qcor(2)),'atom file=',ftof(sum1-sum2)
          endif
          corm = corm + qcor(2)
          if (allocated(rwgt_rv)) deallocate(rwgt_rv)
       endif
       if (lmxl > -1) then
          allocate(v0pot(ib)%v(nr*nsp))
          allocate(v1pot(ib)%v(nr*nsp))
          v0pot(ib)%v=rv_a_ov0a( is )%v
          v1pot(ib)%v=rv_a_ov0a( is )%v
          call dpcopy ( sspec ( is ) %rv_a_orhoc , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , 1d0 )
          if (lfoc == 0) then
             allocate(rwgt_rv(nr))
             call radwgt ( rmt , a , nr , rwgt_rv )
             call radsum ( nr , nr , 1 , nsp , rwgt_rv , sv_p_orhoat( 3 , ib )%v , ssum )
             fac = 1d0
             if(dabs(ssum) > 1d-7) fac = qc/ssum
             if (ipr >= 40) write(stdo,787) is,qc,ssum,fac
787          format(' scale foca=0 core species',i2,': qc,sum,scale=', 3f12.6,f12.6)
             call dpcopy ( sv_p_orhoat( 3 , ib )%v , sv_p_orhoat( 3 , ib )%v, 1 , nr * nsp , fac )
             if (allocated(rwgt_rv)) deallocate(rwgt_rv)
          endif
       endif
       ztot = ztot+z
       ctot = ctot+qc
20  enddo ibloop
    
    if(procid == master) then
       v0wrireblock:block
         real(8):: ov0mean
         integer:: ir,isp
         logical:: cmdopt0,v0write
         character(8):: charext
         !v0write=cmdopt0('--v0write')
         !if(v0write) then
         do ib=1,nbas
            do ir=1,nr
               ov0mean = 0d0
               do isp=1,nsp
                  ov0mean = ov0mean + v0pot(ib)%v( ir + nr*(isp-1) )
               enddo
               ov0mean = ov0mean/nsp !spin averaged
               do isp=1,nsp
                  v0pot(ib)%v(ir + nr*(isp-1))= ov0mean
               enddo
            enddo
            open(newunit=ifi,file='v0pot.'//trim(charext(ib)),form='unformatted')
            write(ifi) v0pot(ib)%v(1:nr)
            close(ifi)
         enddo
         !endif
       endblock v0wrireblock
    endif

    if(allocated(zv_a_osmrho)) deallocate(zv_a_osmrho)
    allocate(zv_a_osmrho(k1*k2*k3,nsp))
    zv_a_osmrho=0d0
    ! --- Overlap smooth hankels to get smooth interstitial density ---
    ng=lat_ng
    allocate(cv_zv(ng*nsp))
    call ovlpfa( nbas , nxi , n0 , exi , hfc , rsmfa, ng , ng , rv_a_ogv , cv_zv )
    call gvputf( ng , nsp , iv_a_okv , k1 , k2 , k3 , cv_zv , zv_a_osmrho )
    if (allocated(cv_zv)) deallocate(cv_zv)
    ! ... FFT to real-space mesh
    call fftz3 ( zv_a_osmrho , n1 , n2 , n3 , k1 , k2 , k3 , nsp, 0 , 1 )
    ! ... Add compensating uniform electron density to compensate background
!    call addbkgsm ( zv_a_osmrho , k1 , k2 , k3 , nsp , qbg , vol, - 1d0 )
    zv_a_osmrho=zv_a_osmrho-qbg/vol/nsp
    ! ... integrate
    !call mshint ( vol , nsp , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho, sum1 , sum2 )
    !if (nsp == 2) then
    !  call mshint ( vol , 1 , n1 , n2 , n3 , k1 , k2 , k3 , zv_a_osmrho, smom , sum2 )
    !endif
    sum1 = dreal(sum(zv_a_osmrho(:,:)))*vol/(n1*n2*n3)
    if(nsp==2) smom = 2d0*dreal(sum(zv_a_osmrho(:,1)))*vol/(n1*n2*n3) - sum1
    ! --- Set up local densities using rmt from atm file ---
    call ovlocr(nbas,n0,nxi,exi,hfc,rsmfa,rv_a_orhofa,sv_p_orhoat,sqloc,slmom)
    ! --- Add compensating uniform electron density to compensate background
    call adbkql ( sv_p_orhoat , nbas , nsp , qbg , vol , - 1d0 )!, sspec )!, ssite )
    if (abs(qbg)/=0d0.and. ipr>=10) write(stdo,ftox) ' Uniform '// &
         'density added to neutralize background q=',ftof(qbg)
    dq = sum1+sqloc+ctot-ztot+qbg !charge
    if (nsp == 1) then
       if (ipr >= 10) write(stdo,895) sum1,sqloc,sum1+sqloc,ctot,-ztot,qbg,dq
895    format(/' Smooth charge on mesh:    ',f16.6 &
            /    ' Sum of local charges:     ',f16.6 &
            /    ' Total valence charge:     ',f16.6 &
            /    ' Sum of core charges:      ',f16.6 &
            /    ' Sum of nuclear charges:   ',f16.6 &
            /    ' Homogeneous background:   ',f16.6 &
            /    ' Deviation from neutrality:',f16.6)
       if (ipr >= 10) write (stdl,710) sum1+sqloc,sum1,sqloc,qbg,dq
710    format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,'   bg',f10.6,'  dQ',f10.6)
    else
       if (ipr >= 10) write(stdo,896) sum1,smom,sqloc,slmom, &
            sum1+sqloc,smom+slmom,ctot,corm,-ztot,qbg,dq
896    format(/' Smooth charge on mesh:    ',f16.6,4x,'moment', f12.6, &
            /    ' Sum of local charges:     ',f16.6,4x,'moments',f11.6, &
            /    ' Total valence charge:     ',f16.6,4x,'moment', f12.6, &
            /    ' Sum of core charges:      ',f16.6,4x,'moment', f12.6, &
            /    ' Sum of nuclear charges:   ',f16.6 &
            /    ' Homogeneous background:   ',f16.6 &
            /    ' Deviation from neutrality:',f16.6)
       if (ipr >= 10) write (stdl,711) sum1+sqloc,sum1,sqloc,qbg,smom+slmom
711    format('ov qvl',f11.6,'  sm',f11.6,'  loc',f11.6,'   bg',f11.6,' mm',f11.6)
    endif
    if (dabs(dq) > 1d-4 .AND. ipr>0) write(stdo,"(' rdovfa (warning) overlapped' &
         //' density not neutral'//', dq= ',d13.5)") dq
    do is=1,nspec
       if (allocated(rv_a_ov0a(is)%v)) deallocate(rv_a_ov0a(is)%v)
       if (allocated(rv_a_orhofa(is)%v)) deallocate(rv_a_orhofa(is)%v)
    enddo
    call tcx('rdovfa')
  end subroutine rdovfa

  subroutine ovlocr(nbas,nxi0,nxi,exi,hfc,rsmfa,rv_a_orhofa, sv_p_orhoat , sqloc, slmom )
    use m_lmfinit,only: nsp,ispec,sspec=>v_sspec
    use m_struc_def
    use m_lgunit,only:stdo
    use m_smhankel,only: hxpbl
    use m_lattic,only: rv_a_opos
    use m_smhankel,only: hxpos
    use m_hansr,only:corprm
    !- Makes the site densities for overlapped free atoms.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   ssite :struct containing site-specific information
    !i   sspec :struct containing species-specific information
    !i   slat  :struct containing information about the lattice
    !i   nxi   :number of Hankels
    !i   nxi0  :leading dimension of hfc
    !i   exi   :smoothed Hankel energies; see Remarks
    !i   hfc   :coefficients to smoothed Hankels
    !i   rsmfa :Hankel smoothing radius
    !i   orhofa:free-atom density, by species
    !o Outputs
    !   orhoat :local density, given by true and smooth densities
    !    sqloc :sum of local charges (integral over rho1-rho2)
    !    slmom :sum of local magnetic moments
    !r Remarks
    !u Updates
    !u   xxx 12 May 07 parallelized (MPI)
    !u   01 Jul 05 Zero-radius empty spheres skip as having no local part
    !u   14 Jun 00 spin polarized
    !u   24 Apr 00 Adapted from nfp ovlocr.f
    ! ----------------------------------------------------------------------
    implicit none
    integer:: kmxv=15 !Hardwired taken from original code.

    integer:: nbas , nxi(1) , nxi0
    type(s_rv1) :: sv_p_orhoat(3,nbas)
    type(s_rv1) :: rv_a_orhofa(nbas)
    real(8):: rsmfa(1) , exi(nxi0,1) , hfc(nxi0,2,1) , sqloc , slmom
    !  type(s_site)::ssite(*)
    !  type(s_spec)::sspec(*)
    integer:: ib , ipr , iprint , is , jb , je , js , lfoca &
         , lmxl , nlmh , nlml , nr , i
    double precision :: ceh,cofg,cofh,eh,qcorg,qcorh,qsc,qcsm,qloc,rfoca,rmt,rsmh,rsmv,z,amom
    double precision :: a,p1(3), p2(3),q(3) !,b0(ktop0+1),acof((ktop0+1),nlmx,2)
    real(8),allocatable:: acof(:,:,:),b0(:,:),rofi(:),rwgt(:)
    complex(8),allocatable:: b(:,:)
    data q /0d0,0d0,0d0/
    integer:: ibini,ibend

    call tcn('ovlocr')
    ipr  = iprint()
    sqloc = 0
    slmom = 0
    if (ipr >= 30) write (stdo,300)
300 format(/' Free atom and overlapped crystal site charges:' &
         /'   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local')
    ibini= 1
    ibend= nbas
    do  ib = ibini,ibend
       is=ispec(ib) !ssite(ib)%spec
       p1(:)=rv_a_opos(:,ib) !ssite(ib)%pos(:)
       lmxl=sspec(is)%lmxl
       !kmxv=sspec(is)%kmxv
       rsmv=sspec(is)%rsmv
       nlml = (lmxl+1)**2
       allocate(acof(0:kmxv,nlml,nsp),b(0:kmxv,nlml))
       acof=0d0
       b=0d0
       a=sspec(is)%a
       nr=sspec(is)%nr
       rmt=sspec(is)%rmt
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoca,rfoca, z)
       qcsm = qcorg+qcorh
       if (lmxl == -1) goto 10
       allocate(rofi(nr),rwgt(nr))
       call radmsh(rmt,a,nr,rofi)
       call radwgt(rmt,a,nr,rwgt)
       !   ... Loop over other sites, add up tail expansion
       do  jb = 1, nbas
          js=ispec(jb) !ssite(jb)%spec
          p2(:)=rv_a_opos(:,jb) !ssite(jb)%pos(:)
          do  je = 1, nxi(js)
             rsmh = rsmfa(js)
             eh   = exi(je,js)
             nlmh = 1
             call hxpbl ( p2,p1,q,[rsmh], rsmv,[eh],kmxv,nlmh,nlml, kmxv,nlml,  b ) 
             allocate(b0(0:kmxv,nlmh))
             b0=0d0
             if (ib == jb) call hxpos([rsmh],rsmv,[eh],kmxv,nlmh,kmxv,b0)
             do  i = 1, nsp
                call p1ovlc(kmxv,nlml,hfc(je,i,js),b,b0,acof(0,1,i))
             enddo
             deallocate(b0)
          enddo
       enddo
       call p2ovlc ( ib,nsp,rsmv,kmxv,nr,nlml,acof,rofi &
           ,rwgt,nxi0,nxi(is),exi(1,is),hfc(1,1,is),rsmfa(is),rv_a_orhofa(is)%v,sv_p_orhoat(3,ib)%v &
           ,lfoca,qcsm,qloc,amom,sv_p_orhoat(1,ib)%v,sv_p_orhoat( 2,ib )%v )
       sqloc = sqloc + qloc
       slmom = slmom + amom
       deallocate(rofi,rwgt)
10     continue
       deallocate(acof,b)
    enddo
    call tcx('ovlocr')
  end subroutine ovlocr
  subroutine p1ovlc(kmxv,nlml,hfc,b,b0,a)
    !- Adds contribution to P_kl expansion of density from one basis function
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nlml  :density expanded to nlml
    !i   kmxv  :k-cutoff for P_kl expansion
    !i   hfc   :coefficient to basis function
    !i   b     :P_kl expansion of density from one basis function
    !i   b0    :P_kl expansion of on-site density
    !o Outputs
    !o   a     :cumulative P_kl expansion of density for this site
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nlml,kmxv
    double precision :: a(0:kmxv,nlml),b0(0:kmxv,1),hfc
    double complex b(0:kmxv,nlml)
    integer :: k,ilm
    do  10  k = 0, kmxv
       do  12  ilm = 1, nlml
          a(k,ilm) = a(k,ilm) + hfc*dble(b(k,ilm))
12     enddo
       a(k,1) = a(k,1) - hfc*b0(k,1)
10  enddo
  end subroutine p1ovlc


  subroutine p2ovlc(ib,nsp,rsmv,kmxv,nr,nlml,acof,rofi,rwgt, &
       nxi0,nxi,exi,hfc,rsmfa,rhofa,rhoc,lfoca,qcsm,qloc,amom,rho1,rho2)
    use m_lgunit,only:stdo
    use m_hansr,only: hansmr
    !- Assemble local density from P_kl expansion for one site
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ib    :site for which to assemble local density
    !i   nsp   :number of spin channels
    !i   rsmv  :smoothing radius for P_kl expansion
    !i   kmxv  :k-cutoff for P_kl expansion
    !i   nr    :number of radial mesh points
    !i   nlml  :L-cutoff for P_kl expansion
    !i   acof  :coefficients to P_kl expansion
    !i   rofi  :radial mesh points for tabulation on a radial mesh
    !i   rwgt  :radial mesh weights for integration on a radial mesh
    !o   rhohd :work array (holds on-site smoothed head density)
    !i   nxi0  :leading dimension of hfc
    !i   nxi   :number of smoothed Hankel energies in head expansion
    !i   exi   :smoothed Hankel energies in head expansion
    !i   hfc   :coefficients to Hankel energies in head expansion
    !i   rsmfa :Hankel smoothing radius in head expansion
    !i   rhofa :head free-atom density
    !i   rhoc  :core density --- used to integrate core charge
    !i   lfoca :switch specifying treatment of core density.
    !i          0 => val,slo = 0 at sphere boundary
    !i          1 => core tails included explicitly with valence
    !i          2 => tails included perturbatively
    !i   qcsm  :smoothed core density
    !o Outputs
    !i   rho1  :local true density, tabulated on a radial mesh
    !i   rho2  :local smoothed density, tabulated on a radial mesh
    !o   qloc  :sphere charge
    !o   amom  :sphere magnetic moment
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: nr,nxi0,ib,nsp,kmxv,nlml,nxi,lfoca
    double precision :: qcsm,qloc,rhofa(nr,nsp),rho1(nr,nlml,nsp), &
         rho2(nr,nlml,nsp),rofi(nr),rwgt(nr),rhohd(nr,nsp),exi(1), &
         hfc(nxi0,nsp),rhoc(nr,nsp),acof(0:kmxv,nlml,nsp),rsmv,rsmfa,amom
    ! ... Local parameters
    !      integer stdo !kmx,lmx,
    !      parameter (kmx=20, lmx=6)
    integer :: i,ie,ilm,ipr,iprint,k,l,ll,lmax,lmxl,isp
    double precision :: asm,gam,pi,qall,qexa,qin,qlc,qnum,qout,qsmo,qut, &
         r,rl,rmt,srfpi,sum,sumfa,sumhd,sumsm,sumtr,y0, &
         xi(0:10),x0(0:2),ddot !pkl(0:kmx,0:lmx)
    real(8),allocatable:: pkl(:,:)

    ipr   = iprint()
    !      stdo  = lgunit(1)
    pi    = 4d0*datan(1d0)
    srfpi = dsqrt(4*pi)
    y0    = 1d0/srfpi
    lmxl  = ll(nlml)
    allocate(pkl(0:kmxv,0:lmxl))
    !      if (lmxl .gt. lmx) call rxi('ovlocr: increase lmx, need',lmxl)

    !     do  ilm = 1, nlml
    !       do  k = 0, kmxv
    !         if (dabs(acof(k,ilm,1)).gt.1d-6)
    !    .      write(stdo,780) ilm,k,acof(k,ilm,1),acof(k,ilm,nsp)
    ! 780     format('ilm,k',2i5,2f14.8)
    !       enddo
    !     enddo

    ! --- Assemble smooth on-site head density in rhohd ---
    qnum = 0d0
    qexa = 0d0
    qsmo = 0d0
    qut = 0d0
    call dpzero(rhohd, nr*nsp)
    asm = 1d0/rsmfa
    lmax = 0
    do  ie = 1, nxi
       sum = 0d0
       do  i = 1, nr
          r = rofi(i)
          call hansmr(r,exi(ie),asm,xi,lmax)
          sum = sum + srfpi*rwgt(i)*xi(0)*r*r
          do  isp = 1, nsp
             rhohd(i,isp) = rhohd(i,isp) + srfpi*hfc(ie,isp)*xi(0)*r*r
          enddo
       enddo
       gam = 0.25d0*rsmfa**2
       qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
       rmt = rofi(nr)
       call hansmr(rmt,0d0,1/rsmfa,x0,1)
       call hansmr(rmt,exi(ie),1/rsmfa,xi,1)
       qout = srfpi/exi(ie)*(-dexp(rsmfa**2/4*exi(ie)) &
            - rmt**3*(xi(1)-dexp(rsmfa**2/4*exi(ie))*x0(1)))
       qin = qall-qout
       do  isp = 1, nsp
          qnum = qnum + hfc(ie,isp)*sum
          qexa = qexa + hfc(ie,isp)*qin
          qsmo = qsmo + hfc(ie,isp)*qall
          qut  = qut  + hfc(ie,isp)*qout
       enddo
    enddo

    !|      write(stdo,917) qnum,qexa,qsmo,qut
    !|  917 format('summed smooth charge:  num',f14.8,'   exact',f14.8
    !|     .   /' total smooth q',f14.8,'  outside',f14.8)

    ! --- Assemble overlapped tail density in rho2 ---
    !      if (kmxv .gt. kmx) call rx('ovlocr: increase kmx')
    call dpzero(rho2,  nr*nlml*nsp)
    do  i = 1, nr
       r = rofi(i)
       call radpkl(r,rsmv,kmxv,lmxl,kmxv,pkl)
       do  isp = 1, nsp
          do  ilm = 1, nlml
             l = ll(ilm)
             rl = 0.d0
             if ( r > 0.d0 ) rl = r**l
             do  k = 0, kmxv
                rho2(i,ilm,isp) = rho2(i,ilm,isp) + &
                     acof(k,ilm,isp)*pkl(k,l)*r*r*rl
             enddo
          enddo
       enddo
    enddo
    ! ... Make the true density in rho1, smooth density in rho2
    call dpcopy(rho2,rho1,1,nr*nlml*nsp,1d0)
    do   isp = 1, nsp
       do   i = 1, nr
          rho1(i,1,isp) = rho1(i,1,isp) + y0*rhofa(i,isp)
          rho2(i,1,isp) = rho2(i,1,isp) + y0*rhohd(i,isp)
       enddo
    enddo
    ! ... Do some integrals
    sumfa = 0d0
    sumsm = 0d0
    sumhd = 0d0
    sumtr = 0d0
    qlc = 0d0
    do   isp = 1, nsp
       do   i = 1, nr
          sumfa = sumfa + rwgt(i)*rhofa(i,isp)
          sumhd = sumhd + rwgt(i)*rhohd(i,isp)
          sumtr = sumtr + rwgt(i)*rho1(i,1,isp)
          sumsm = sumsm + rwgt(i)*rho2(i,1,isp)
          qlc = qlc + rwgt(i)*rhoc(i,isp)
       enddo
    enddo
    sumsm = sumsm*srfpi
    sumtr = sumtr*srfpi
    qloc = sumtr-sumsm
    amom = -srfpi* &
         (ddot(nr,rwgt,1,rho1(1,1,nsp),1)-ddot(nr,rwgt,1,rho1(1,1,1),1) &
         -ddot(nr,rwgt,1,rho2(1,1,nsp),1)+ddot(nr,rwgt,1,rho2(1,1,1),1))
    if (lfoca == 0) qloc = qloc + qlc - qcsm
    if (ipr >= 30) then
       write(stdo,810) ib,sumfa,sumhd,sumtr,sumsm,qloc
       if (nsp == 2) write(stdo,811) &
            ddot(nr,rwgt,1,rhofa,1)-ddot(nr,rwgt,1,rhofa(1,2),1), &
            ddot(nr,rwgt,1,rhohd,1)-ddot(nr,rwgt,1,rhohd(1,2),1), &
            srfpi*(ddot(nr,rwgt,1,rho1,1)-ddot(nr,rwgt,1,rho1(1,1,2),1)), &
            srfpi*(ddot(nr,rwgt,1,rho2,1)-ddot(nr,rwgt,1,rho2(1,1,2),1)), &
            amom
    endif
810 format(i5,6f12.6)
811 format(' amom',6f12.6)
  end subroutine p2ovlc
  subroutine adbkql( sv_p_orhoat , nbas , nsp , qbg , vol , fac )
    use m_struc_def
    use m_lmfinit,only: ispec, sspec=>v_sspec
    !- Add uniform bkg charge density to local smooth rho
    !i orhoat: pointers to local density in spheres
    !i nbas: number of atoms in basis
    !i qbg: background charge
    !i sspec: species structure
    !i nsp: spins
    !i vol: vol of cell
    !i fac: fac * backg density is added
    !u Updates
    !u   01 Jul 05 Zero-radius sites skipped over
    !----------------------------------------
    implicit none
    integer :: nrmx,nlmx,nlml,lmxl,nbas
    parameter (nrmx=1501,nlmx=64)
    integer:: nsp
    type(s_rv1) :: sv_p_orhoat(3,nbas)
    real(8):: qbg , fac
    integer :: ib,nr,is
    double precision :: rhobkg,vol,a,rmt,rofi(nrmx)
    rhobkg = fac*qbg/vol
    do  ib = 1, nbas
       is=ispec(ib)
       a=sspec(is)%a
       nr=sspec(is)%nr
       rmt=sspec(is)%rmt
       lmxl=sspec(is)%lmxl
       if (lmxl == -1) goto 10
       nlml=(lmxl+1)**2
       call rxx(nr .gt. nrmx,  'addbkgloc: increase nrmx')
       call rxx(nlml .gt. nlmx,'addbkgloc: increase nlmx')
       call radmsh(rmt,a,nr,rofi)
       call addbkgl(sv_p_orhoat(1,ib )%v,sv_p_orhoat(2,ib)%v, rhobkg , nr , nsp , rofi , nlml )
10     continue
    enddo
  end subroutine adbkql
  subroutine addbkgl(rho1,rho2,rhobkg,nr,nsp,rofi,nlml)
    ! adds uniform background to local smooth density at this site for l=0 component (ilm=1) 
    implicit none
    integer :: nsp,is,nr,nlml,i
    real(8):: rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rofi(nr),rhobkg
    real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4*pi)
    do is = 1, nsp
       rho1(:,1,is) = rho1(:,1,is)+srfpi*rofi(:)**2*rhobkg/nsp
       rho2(:,1,is) = rho2(:,1,is)+srfpi*rofi(:)**2*rhobkg/nsp
    enddo
  end subroutine addbkgl
end module m_rdovfa
