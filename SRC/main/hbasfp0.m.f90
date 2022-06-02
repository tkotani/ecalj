program hbasfp0
  !-- Generates orthonormal optimal product basis and required radial integrals in each MT.
  ! input files
  !  GWinput : input data for GW
  !  PHICV  : radial functions Valence and Core

  ! output files
  !  BASFP//ibas :: product basis for each ibas
  !  PPBRD_V2_//ibas :: radial <ppb> integrals. Note indexing of ppbrd

  ! The main part of this routine is in the subroutine basnfp_v2
  !     use m_rgwinf_v3,only:rgwinf_v3,
  use m_genallcf_v3,only: Genallcf_v3, &
       alat,nclass,natom,nspin,nl,nnv,nnc,nrx, cutbase,lcutmx,nindxc, &
       nindxv,occv,unoccv,occc,unoccc,iclass
  use m_keyvalue,only: getkeyvalue
  use m_anf,only: ibasf,laf,anfcond !may2015takao
  !      use m_lgunit,only: m_lgunit_init
  !      use m_mpi,only: MPI__Initialize
  implicit none
  integer:: &
       ifphiv(2),ifphic(2), iphiv(2),iphivd(2),iphic(2),iphi(2),iphidot(2), &
       ifev(2),ifevf(2),ibas,ibas1,ic,icx,ifaln,iflmto,ifphi, &
       ii,ir,irad,isp,ix,lmx,lmx2,n,nbas,ncoremx,l,nn,icore,ifianf,nphi,nradmx,nsp,maxnn, &
       idummy,ifinin,incwfin,ret,ndat
  integer,allocatable:: lcutmxa(:),nrofi(:), nocc(:,:),nunocc(:,:),nindx(:,:), &
       ncindx(:,:),lcindx(:,:), nrad(:), nindx_r(:,:), lindx_r(:,:), nc_max(:,:),ncore(:)
  real(8):: qbas(3,3),ginv(3,3)
  real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:),phitotr(:,:,:,:,:)
  real(8),allocatable:: cutbasex(:), bb(:),zz(:), phic(:,:)
  logical :: checkdid ,anfexist, cmdopt2, ptest=.false. !See ptest in hvccfp0.f
  character(12) :: aaa
  character(11) :: ffaln
  character(20):: outs=''
  !      call M_lgunit_init()
  write(6,'(a)') ' --- Input normal(=0); coremode(=3);'// &
       ' ptest(=4); Excore(=5); for core-valence Ex(=6);'// &
       ' val-val Ex(7);  normal+<rho_spin|B> (8); version(-9999) ?'
  if(cmdopt2('--job=',outs)) then
     read(outs,*) ix
  else
     read(5,*) ix
  endif

  selectcase(ix)
  case(3)
     write(6,*)'  ### coremode; Product basis for SEXcore ### '
     incwfin = -2
  case(0)
     write(6,*)'  ### usual mode use occ and unocc for core ### '
     incwfin = 0
  case(4)
     write(6,*) ' ### ptest mode. now special for Q0P. GWIN_V2 is neglected ### '
     write(6,*) '  See basnfp.f of ptest section.'
     incwfin = 0
  case(5)
     write(6,*) '  ### calculate core exchange energy ### ix==5'
     incwfin = 0
  case(6)
     write(6,*) '  ### calculate p-basis for core-valence Ex ix==6'
     write(6,*) ' occ=1:unocc=0 for all core'
     incwfin = -3
  case(7)
     write(6,*) '  ### calculate p-basis for val-val Ex ix==7'
     write(6,*) ' occ=0:unocc=0 for all core'
     incwfin = -4
  case(8)
     write(6,"('  ### usual mode use occ and unocc for core', &
          ' and <rho_spin |B(I)> ### ')")
     incwfin = 0
  case default
     write(6,*)' hbasfp: input is out of range'
     call rx( ' hbasfp: input is out of range')
  endselect

  !! incwfin gives the setting of Product basis
  call Genallcf_v3(incwfin)
  nsp=nspin
  !! readin lcutmxa ------------
  call Getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifinin,status=ret)
  allocate(lcutmxa(1:natom))
  do
     read(ifinin,*,err=980) aaa
     if(aaa=='lcutmx(atom)') then
        read(ifinin,*) lcutmxa(1:natom)
        !          write(6, '(" lcutmxa=",20i3)' ) lcutmxa(1:natom)
        goto 990
     endif
  enddo
980 continue
  lcutmxa=lcutmx
990 continue
  close(ifinin)
  if(ix==8) then
     write(6,*)' Enfoece lcutmx=0 for all atoms'
     lcutmxa=0
  endif
  write(6,"(' lcutmxa=',$)")
  write(6,'(20i3)') lcutmxa(1:natom)
  lmx        = 2*(nl-1)
  lmx2       = (lmx+1)**2
  nn         = maxnn (nindxv,nindxc,nl,nclass)
  nphi       = nrx*nl*nn*nclass

  ! -optimal orthonormal product basis
  !> reindex nocc,nunocc,nindx
  ! For valence  from GWIN_V2
  ! occv   : occ    switch
  ! unoccv : unocc  switch
  ! nindexv: n index
  !-----------------------------
  ! For core  from GWIN_V2
  ! occc   : occ  switch
  ! unoccc : unocc switch
  ! nindexc: n index
  !-----------------------------
  ! For valence+core
  ! nocc
  ! nunocc
  ! nindx
  allocate( nocc(nl*nn,nclass), nunocc(nl*nn,nclass), nindx(nl,nclass) )
  call reindx  (occv,unoccv,nindxv,   occc,unoccc,nindxc, &
       nl,nn,nnv,nnc,nclass, &
       nocc,nunocc,nindx)
  write(6,*)' --- end of reindx ---'
  !-----------
  ! read PHIVC  and reserve it to phitot
  !----------
  open(newunit=ifphi,file='PHIVC',form='unformatted')
  read(ifphi) nbas, nradmx, ncoremx
  allocate(  ncindx(ncoremx,nbas), &
       lcindx(ncoremx,nbas), &
       nrad(nbas), &
       nindx_r(1:nradmx,1:nbas), &
       lindx_r(1:nradmx,1:nbas), &
       aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) , &
       phitoto(nrx,0:nl-1,nn,nbas,nsp), &
       phitotr(nrx,0:nl-1,nn,nbas,nsp), &
       nc_max(0:nl-1,nbas),ncore(nbas) )
  read(ifphi) nrad(1:nbas)
  read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
  nc_max=0
  do ibas=1,nbas
     write(6,*)' --- read PHIVC of ibas=',ibas
     ic = ibas
     read(ifphi) ncore(ic), ncoremx                            !core
     read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
     read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
     if(ic/=icx) then
        write(6,*) 'ic icx=',ic,icx
        call rx( 'hbasfp0: ic/=icx')
     endif
     read(ifphi) rr(1:nrofi(ic),ic)
     do isp = 1, nsp
        write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
        do icore = 1, ncore(ic)
           l =  lcindx(icore,ic)
           n =  ncindx(icore,ic)
           read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
           phitotr(1:nrofi(ic),l,n, ic,isp)=      &       !  core raw= core orthgonal
           phitoto(1:nrofi(ic),l,n, ic,isp)               !
           if(n>nc_max(l,ic)) nc_max(l,ic)=n
        enddo
        do irad = 1, nrad(ic)
           l = lindx_r (irad,ic)
           n = nindx_r (irad,ic) + nc_max(l,ic)
           read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
           read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
        enddo
     enddo
  enddo
  close(ifphi)
  !-----------
  !! check write
  open(newunit=ifaln,file='PHIV.chk')
  do ibas = 1,nbas
     ic = ibas
     do irad = 1, nrad(ic)
        l = lindx_r (irad,ic)
        n = nindx_r (irad,ic) + nc_max(l,ic)
        write(ifaln,"(a,5i5)")'------- ibas l n =',ibas,l,n
        do ir=1,nrofi(ic)
           write(ifaln,"(3d24.15)")rr(ir,ic), phitotr(ir,l,n,ic,1:nsp)
        enddo
     enddo
  enddo
  close(ifaln)
  !!  excore mode ---------
  if(ix==5 ) then
     call excore(nrx,nl,nnc,nclass,nsp,natom, &
          phitotr(1:nrx,0:nl-1,1:nnc,1:nclass,1:nsp), nindxc,iclass, aa,bb,nrofi,rr)
     goto 998
  endif
  !! antiferro or not.
  !! For AF case, we have laf=.true. and we have data set for 'call anfsig', stored in m_anf.
  call anfcond()
  if(laf) then ! Check iclass =ibas ; CLASS file contains true classs information.
     write(6,*) '--- Antiferro mode --- '
     do ibas=1,natom
        if(iclass(ibas)/=ibas) call rx( ' iclass(ibas)/=ibas: ')
     enddo
     ii=0
     do ic=1,nclass
        ibas=ic
        if( ibasf(ibas)>0 ) then
           phitotr(:,:,:,ibasf(ibas), :)=phitotr(:,:,:,ibas, :)
           write(6,"(a,2i4)") &
                '  radial functions: phi(ibasf)=phi(ibas): ibasf ibas=',ibasf(ibas),ibas
        endif
     enddo
  endif
  !! override cutbase to make epsPP_lmfh safer. may2013takao
  allocate( cutbasex(0:2*(nl-1)) )
  cutbasex=cutbase
  if(ix==4) then
     write(6,*)' !!! set tolerance for PB to be 1d-6 ---'
     cutbasex=1d-6
  endif
  do ic = 1,nclass
     call basnfp_v2(nocc(1,ic),nunocc(1,ic),nindx(1,ic),nl,nn,nrx, nrofi(ic),rr(1,ic),aa(ic),bb(ic),ic, &
          phitoto,phitotr,nsp,nclass, cutbasex, lcutmxa(ic),ix,alat,nc_max(0,ic) )
  end do
  if(ix==0) call rx0s( ' OK! hbasfp0 ix=0 normal mode ')
  if(ix==3) call rx0s( ' OK! hbasfp0 ix=3 core mode ')
  if(ix==4) call rx0s( ' OK! hbasfp0 ix=4 ptest mode  ')
  if(ix==6) call rx0s( ' OK! hbasfp0 ix=6 Exx core-val mode  ')
  if(ix==7) call rx0s( ' OK! hbasfp0 ix=7 Exx val-val mode  ')
  if(ix==8) call rx0s( ' OK! hbasfp0 ix=8 normal(ix==0) + <B|spin den>. Enforce lcutmx=0.')
998 if(ix==5) call rx0s( ' OK! hbasfp0 ix=5 ex core mode  ')
END PROGRAM hbasfp0
