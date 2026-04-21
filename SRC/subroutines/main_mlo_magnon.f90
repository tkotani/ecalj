!>  Calculate Chi^+-, spin susceptibility. 
module m_mlo_magnon 
  public :: mlo_magnon
  contains
subroutine mlo_magnon() bind(C)
  use m_mlo_ham, only: read_ham_rs, calc_ham_eigen, nwf => ndimMTO, nsite, ib_tableM, nspx
  use m_mlo_scrw, only: nnwf, scrw, mlo_pairs, trace, trace_onsite, nnwf_init, scrw_init,  &
                        contract_to_site, pair_site, pair_lorb, nnwf_mask, nnwf2_mask, trace2
  use m_mlo_formfactor, only: get_formfactor_q
  use m_mlo_ovlppair, only: get_ovlppair_q
  use m_HamPMT,only: ReadHamPMTInfo
  use m_ReadEfermi, only: readefermi
  use m_read_bzdata, only: nqbz, qbz
  use m_genallcf_v3, only: genallcf_v3, nspin
  use m_keyvalue, only: getkeyvalue
  use m_freq, only: getfreq, freq_r, nwhis, nw_i, nw, npm
  use m_tetwt, only: tetdeallocate, gettetwt, whw, ihw, nhw, jhw, n1b, n2b, nbnb
  use m_readgwinput, only: ReadGWinputKeys
  use m_lgunit, only: m_lgunit_init, stdo
  use m_dpsion, only: dpsion_init, dpsion_chiq
  use m_mpi, only: MPI__Initialize, MPI__consoleout, MPI__SplitXq
  use m_mpi, only: mpi__rank, mpi__size, mpi__root, comm, comm_k, mpi__rank_k, mpi__size_k, mpi__root_k, ipr
  use m_mpiio, only: openm, closem, writem, readm, mpiio_buf, buf_put, buf_get, writem_buf, readm_buf
  use m_blas, only: m_op_C, zmm => zmm_h, int_split
  use m_lapack, only: zminv => zminv_h, zhev => zhev_h, zgev => zgev_h, zhgv_lindep => zhgv_lindep_h, zggv => zggv_h
  use m_mem, only: writemem
  use m_ftox, only: ftox
  implicit none
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  xxx removed--> eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  integer:: iwf, jwf, inwf, jnwf
  integer:: file_tr_kr, file_jq_site, file_jq_full
  integer:: iqxini, iqxend, i, iw, iq, kx, istat, nqcalc
  real(8):: q(3), rydberg, hartree, delta, eta, delta_dos, freq_ratio, freq_dw
  real(8), allocatable:: qibze(:,:)
  complex(8), pointer:: zxq(:,:,:) => null()
  complex(8), allocatable, target :: kmat(:,:,:)
  complex(8), allocatable:: imat(:,:)
  complex(8), parameter :: img=(0d0,1d0)
  logical:: cmdopt0
  logical:: realomega, imagomega, epsmode
  logical, allocatable :: mpi__task(:)
  character(8):: charext
  character(len=128) :: msg
  integer, parameter :: is=1, isf=2  !K_down up = Kpm
  real(8), parameter :: pi = 4d0*datan(1d0), eta_default =-1d0
  logical, parameter :: nnwf_size_reduction = .true.
  logical :: w_onsite_dddd, geteta, negative_cut, ganmma_only, gettetwt_split, calcdos
  !For dos calculation
  integer :: nqibz_dos, ntetf_dos, nteti_dos
  integer, allocatable :: idteti_dos(:,:)
  real(8), allocatable :: qibz_dos(:,:), rho(:,:,:), sz(:), sz_site(:)
  real(8), allocatable :: freq(:)
  complex(8), allocatable :: ovlppair(:,:), formfactor(:), ovlppair_inv(:,:)

  hartree = 2d0*rydberg()
  geteta = cmdopt0('--geteta')
  calcdos  = cmdopt0('--dos')
  ganmma_only = geteta  !GammaPoint only calculation

  call m_lgunit_init()
  call MPI__Initialize()
  msg ='mlo_magnon'
  if(geteta) msg = trim(msg)//'_geteta_mode'
  call MPI__consoleout(trim(msg))
  realomega = .true.
  imagomega = .false.
  epsmode   = .true.
  call genallcf_v3(incwfx=0) !!incwfin=0 =>ForX0 for core in GWIN. in module m_genallcf_v3 Readin by genallcf. Set basic data for crystal
  if(nspin < 2) call rx(' mlo_magnon: nspin<2: not supported. exit.')
  call ReadGWinputKeys() ! jun2020 new routint to read all inputs
  call getkeyvalue("GWinput","magnon_w_onsite_dddd",w_onsite_dddd,default=.true.)
  call getkeyvalue("GWinput","magnon_delta", delta, default=0d0) !1d-6 for Insulator case
  call getkeyvalue("GWinput","magnon_delta_dos", delta_dos, default=1d-6)
  call getkeyvalue("GWinput","HistBin_ratio",freq_ratio, default=1.03d0)
  call getkeyvalue("GWinput","HistBin_dw",freq_dw, default=1d-5)
  call getkeyvalue("GWinput","magnon_HistBin_ratio",freq_ratio, default=freq_ratio)
  call getkeyvalue("GWinput","magnon_HistBin_dw", freq_dw, default=freq_dw)
  call getkeyvalue("GWinput","magnon_negative_cut",negative_cut,default=.false.)
  if(ipr) write(stdo,ftox) "magnon_w_onsite_dddd", w_onsite_dddd
  if(ipr) write(stdo,ftox) "magnon_geteta", geteta
  if(ipr) write(stdo,ftox) "magnon_delta", delta
  if(ipr) write(stdo,ftox) "magnon_delta_dos", delta_dos
  if(ipr) write(stdo,ftox) "magnon_negative_cut", negative_cut
  if(ipr) write(stdo,ftox) "HistBin_ratio/HistBin_dw", freq_ratio, freq_dw
  if(calcdos) then
    delta = delta_dos
    if(ipr) write(stdo,ftox) "dos calculation: delta_dos is used", delta
  endif

  SetBZDATAandQveclist: block
    use m_read_bzdata, only:read_BZDATA, idteti, nteti, ntetf, qibz, nqibz, nq0i, q0i
    if(calcdos) then
      call read_BZDATA(dosmesh=.true.) !read __BZDATA.DOS
      nqibz_dos = nqibz
      ntetf_dos = ntetf
      nteti_dos = nteti
      allocate(qibz_dos(3,nqibz_dos), source=qibz(1:3,1:nqibz_dos))
      allocate(idteti_dos(0:4,nteti_dos), source=idteti(0:4,1:nteti_dos))
    endif
    call read_BZDATA() !overwrite by __BZDATA
    if(mpi__root) then
      do i=1, nqbz
        if(i<10 .OR. i>nqbz-10) write(stdo,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
        if(i==10 .AND. nqbz>18) write(stdo,"('... ')")
      enddo
      write(stdo,*)' !!nqbz =',nqbz
    endif
    if(ganmma_only) then
      allocate(qibze(3,1), source = 0d0)
      iqxini = 1
      iqxend = 1
    elseif(calcdos) then
      allocate(qibze(3,nqibz_dos), source = qibz_dos(1:3,1:nqibz_dos))
      iqxini = 1
      iqxend = nqibz_dos
    else
      allocate(qibze(3,nq0i), source = q0i(1:3,1:nq0i))
      iqxini = 1
      iqxend = nq0i
    endif
    if(ipr) write(stdo,"('iqxini,iqxend ',2I8)") iqxini, iqxend
    do iq = iqxini, iqxend
      if(ipr) write(stdo,"('iq, qibze:',I8,3f9.4)") iq, qibze(:,iq)
    enddo
  endblock SetBZDATAandQveclist

  SetMPI_Rankdivider: block
    integer :: n_bpara, n_kpara, worker_inQtask
    logical:: cmdopt2
    character(20):: outs
    nqcalc = iqxend-iqxini+1
    n_bpara = 1
    n_kpara = max(mpi__size/(n_bpara*nqcalc), 1)  !Default setting of parallelization. b-parallel is 1.
    if(cmdopt2('--nk=', outs)) read(outs,*) n_kpara
    worker_inQtask = n_bpara * n_kpara
    ! gettetwt_split = n_kpara > 2
    gettetwt_split = .true. ! We always split gettetwt by k to prevent memory exhaustion.
    if(ipr) write(stdo,ftox) 'MPI: worker_inQtask', worker_inQtask
    allocate(mpi__task(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=iqxini,iqxend)])
    if(ipr) write(stdo,ftox) 'mpi_rank',mpi__rank,'mpi__Qtask=',mpi__task
    call MPI__SplitXq(n_bpara, n_kpara)
  endblock SetMPI_Rankdivider

  SetFreqencyMesh:block
    use m_freq, only: niw
    integer :: niw_in
    real(8) :: wemax, omg2max
    real(8), parameter :: ua = 1d0
    ! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
    ! wemax   = 5d0 !max value for plot
    ! omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
    ! narrow energy range: 2026-03-01
    wemax = 2d0 !max value for plot
    omg2max = wemax*.5d0+.5d0 ! (in Hartree) covers all relevant omega, +.2 for margin
    !! NOTE: npmtwo=T sets npm=2   !! optional npmtwo is added aug2017   !! 20190604 Im[K]
    if( .NOT. imagomega) niw_in=1  !dummy
    call Getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw_in,ua, npmtwo=.true.,dw=freq_dw, ratio=freq_ratio)!,tetra
    if(ipr) write(6,"(' nw_i nw niw npm=',4i5)") nw_i,nw,niw,npm
    allocate(freq(nw_i:nw), source = ([-freq_r(nw:1:-1), freq_r(0:nw)]))
    ! write(stdo,ftox) freq
  endblock SetFreqencyMesh

  call readefermi() !!! ef:     Fermi energy at 0 K

  SetMLOAndScreendCoulombData: block
    logical :: cmdopt2
    character(20):: Wtype, opts
    call ReadHamPMTInfo()  ! Read info from PMTHamiltonianInfo (lattice structures and index of basis).
    call read_ham_rs()
    call nnwf_init(nnwf_size_reduction) !set nnwf ~ # of RiRj (onsite_approx = .true.), RiR'j (onsite_approx = .flase. ), wan_pair_index
    if(ipr) write(stdo,ftox) '# nwf, nnwf:', nwf, nnwf
    Wtype = 'up' !options: up, down, up_down, down_up
    if(cmdopt2('--Wtype=', opts)) Wtype = trim(opts)
    call scrw_init(w_onsite_dddd, Wtype=Wtype, enforce_Hermite=.false.) !set scrw
  endblock SetMLOAndScreendCoulombData

  allocate(rho(nwf,nwf,nspin), source = 0d0)
  ReadEta: if(.not. geteta) then
    block
      logical :: exist_etafile
      integer :: ios, iunit, site
      logical, allocatable :: mask(:)
      eta = eta_default
      inquire(file='__EtaMagnon',exist=exist_etafile)
      if(exist_etafile) then
        open(newunit=iunit, file='__EtaMagnon',status='old',form='unformatted',action='read')
        read(iunit, iostat=ios) eta
        if(ios == 0 .and. ipr) write(stdo, ftox) "# Read eta from __EtaMagnon"
        if(ios /= 0) eta = eta_default
        read(iunit, iostat=ios) rho(:,:,:)
        close(iunit)
        allocate(sz(nnwf))
        sz(1:nnwf) = pack(reshape((rho(:,:,1)-rho(:,:,2)), shape=[nwf*nwf]), mask = nnwf_mask)
        allocate(sz_site(nsite))
        do site=1, nsite
          mask = [((pair_site(inwf,1) == site .and. pair_site(inwf,2) == site), inwf=1,nnwf)]
          sz_site(site) = sum(pack(sz, mask=mask))
        enddo
        where(abs(sz_site) < 0.01d0) sz_site = 0d0
        write(stdo,ftox) 'sz_site (threshold = 0.01):', sz_site(1:nsite)
      endif
      if(ipr) write(stdo, ftox) "# WK is scalled to  eta WK: eta =",eta
    endblock
  endif ReadEta

  SetTemporaryFiles: if(.not. geteta) then
    istat = openm(newunit=file_tr_kr, file='__TrKR', recl=(nw-nw_i+1)*16*7)
    istat = openm(newunit=file_jq_site, file='__JqSite', recl=nsite*(nw-nw_i+1)*16)
    istat = openm(newunit=file_jq_full, file='__JqFull', recl=nnwf*(nw-nw_i+1)*16)
  endif SetTemporaryFiles

  allocate(imat(1:nnwf,1:nnwf),source=(0d0,0d0))
  forall(iwf=1:nnwf) imat(iwf,iwf) = 1d0 + img*delta !check the sign
  allocate(kmat(1:nnwf,1:nnwf,(1-npm)*nwhis:nwhis))
  allocate(ovlppair(nnwf,nnwf), formfactor(nnwf), ovlppair_inv(nnwf,nnwf))
  BIGiqloop: do iq = iqxini,iqxend
    if(.NOT. MPI__task(iq)) cycle BIGiqloop
    q = qibze(:,iq)
    if(ipr) write(6,"('===== do : iq wibz(iq) q=',i6,f13.6,3f9.4,' ========')") iq,q !,wibz(iqlist(iq)),qshort !qq
    GETzxq: block ! zxq and zxqi are the main output after Hilbert transformation, ! zxqi is not used in hmagnon (imagomega=.false.)
      use m_mpi,only: MPI__AllreduceSumReal, MPI__AllreduceSum
      real(8) :: evkx_w1(nwf,nspx), evkx_w2(nwf,nspx) !dummy
      complex(8) :: zxqi(1,1,1), evc_w1(nwf,nwf), evc_w2(nwf,nwf)
      complex(8) :: ov_evc_w1(nwf,nwf), ov_evc_w2(nwf,nwf)
      complex(8) :: ovlp_w1(nwf,nwf), ovlp_w2(nwf,nwf), oovlp4(nwf,nwf,nwf,nwf)
      complex(8), allocatable :: evc_w1_kx(:,:,:), evc_w2_kx(:,:,:)
      integer, allocatable :: nttp(:),  itw(:,:), itpw(:,:), ik(:,:)
      integer :: nttp_max, ittp, jpm, it, itp, ibib, isdummy, kx_ini, kx_fin, kx_num, kx_start,kx_end
      integer :: kwf, lwf
      real(8), allocatable :: whwc(:,:), ev_w1(:,:), ev_w2(:,:)
      real(8), parameter:: schi = 1d0
      integer, parameter:: nkblock = 1024
      complex(8), allocatable :: zw(:,:), wzw(:,:)
      allocate(ev_w1(nwf,nqbz), ev_w2(nwf,nqbz), source=0d0)

      if(ipr) call writemem('mlo_magnon start gettetwt')
      call int_split(nqbz, mpi__size_k, mpi__rank_k, kx_ini, kx_fin, kx_num)
      CalcEigenEnergy: do kx = kx_ini, kx_fin !!! ev_w1, ev_w2 unit: [Ry]
        call calc_ham_eigen( is,  is,   qbz(:,kx), evkx_w1(:,:), evec=evc_w1, ovlp_evec=ov_evc_w1)
        call calc_ham_eigen(isf, isf, q+qbz(:,kx), evkx_w2(:,:), evec=evc_w2, ovlp_evec=ov_evc_w2)
        ev_w1(:,kx) = evkx_w1(:, is)
        ev_w2(:,kx) = evkx_w2(:,isf)
        if(ganmma_only) then
          block
          use m_ReadEfermi,only: ef
          real(8):: occ1, occ2
          integer :: iband
          do iband = 1, nwf
            occ1 = merge(1d0,0d0,ev_w1(iband,kx)<ef)/dble(nqbz)
            occ2 = merge(1d0,0d0,ev_w2(iband,kx)<ef)/dble(nqbz)
            do jwf = 1, nwf
              do iwf = 1, nwf
                rho(iwf,jwf, is) = rho(iwf,jwf, is) + dble(dconjg(ov_evc_w1(iwf,iband))*evc_w1(jwf,iband))*occ1
                rho(iwf,jwf,isf) = rho(iwf,jwf,isf) + dble(dconjg(ov_evc_w2(iwf,iband))*evc_w2(jwf,iband))*occ2
              enddo
            enddo
          enddo
          endblock
        endif
      enddo CalcEigenEnergy
      call MPI__AllreduceSumReal(ev_w1, nwf*nqbz, communicator=comm_k)
      call MPI__AllreduceSumReal(ev_w2, nwf*nqbz, communicator=comm_k)
      if(ganmma_only) call MPI__AllreduceSumReal(rho, nwf*nwf*nspin, communicator=comm_k)
      if(.not.gettetwt_split) call gettetwt(q,iq,isdummy,isdummy,ev_w1,ev_w2,nwf,.true.) !! tetrahedron weight. iq is dummy index
        !!     ihw(ibjb,kx): omega index, to specify the section of the histogram., ibjb=1,nbnb
        !!     nhw(ibjb,kx): the number of histogram sections
        !!     jhw(ibjb,kx): pointer to whw
        !!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
        !!     : histogram weights for given ib,jb,kx for histogram sections
        !!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.

      if(ipr) call writemem('mlo_magnon start Im kmat')
      kmat(:,:,:) = 0d0
      kxblock_loop: do kx_start = kx_ini, kx_fin, nkblock
        kx_end = min(kx_fin, kx_start + nkblock -1)
        allocate(evc_w1_kx(nwf,nwf,kx_start:kx_end), evc_w2_kx(nwf,nwf,kx_start:kx_end))
        if(gettetwt_split) call gettetwt(q,iq,isdummy,isdummy,ev_w1,ev_w2,nwf,.true.,ikbz_in=kx_start,fkbz_in=kx_end)
        CalcEigenFunction: do kx = kx_start, kx_end
        call calc_ham_eigen( is, is,  qbz(:,kx), evkx_w1(:,:), evec=evc_w1_kx(:,:,kx)) !evkx_w1  is dummy
        call calc_ham_eigen(isf,isf,q+qbz(:,kx), evkx_w2(:,:), evec=evc_w2_kx(:,:,kx)) !evkx_w2  is dummy
        enddo CalcEigenFunction

        jpmloop:do jpm=1, npm ! jpm=2: negative frequency
!           ibibloop: do 2013 ibib=1,nbnb(kx,jpm) !! n,n' pair band index loop
!             it=n1b(ibib,kx,jpm)  !index for n  for q   ! n1b(ibib,k,jpm) = n :band index for k (occupied),   
!             itp=n2b(ibib,kx,jpm) !index for n' for q+k ! n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
!             wanmat=0d0
!             do concurrent(iwf=1:nwf,jwf=1:nwf, kwf=1:nwf,lwf=1:nwf)
!               ijwf=(iwf-1)*nwf+jwf
!               klwf=(kwf-1)*nwf+lwf 
!               if ( ijwf > klwf ) cycle         ! calculate numerator of Kmatrix
!               wan_j=dconjg(evc_w1(jwf,it,kx))  !a_{Rj   beta}^{kn}*
!               wan_i=evc_w2(iwf,itp,kx)         !a_{Ri  alpha}^{(k+q)n'}
!               wan_l=evc_w1(lwf,it,kx)          !a_{R'l  beta}^{kn}
!               wan_k=dconjg(evc_w2(kwf,itp,kx)) !a_{R'k alpha}^{(k+q)n'}*
!               wanijkl=wan_j*wan_i*wan_k*wan_l
!               wanmat(ijwf,klwf)=wanijkl
!               wanmat(klwf,ijwf)=dconjg(wanijkl) !!! Suppose Hermite Kmat ! wanmat (dimension:nnwf)
!             enddo
!             do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
!               imagweight=whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
!               kmat(:,:,iw,jpm)=kmat(:,:,iw,jpm)+imagweight*wanmat(:,:) !accumulate Im[K]
!             enddo
! 2013      enddo ibibloop
         ! 2025-12-02 MO optimize calculation of kmat same as in x0gemm
         allocate(nttp(nwhis), source = 0)
         do kx = kx_start, kx_end
           do ibib = 1, nbnb(kx,jpm)
             do iw = ihw(ibib,kx,jpm), ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
               nttp(iw) = nttp(iw) + 1
             enddo
           enddo
         enddo
         nttp_max = maxval(nttp(1:nwhis))
         allocate(itw(nttp_max,nwhis), itpw(nttp_max,nwhis), whwc(nttp_max,nwhis), ik(nttp_max,nwhis))
         nttp(:) = 0
         do kx = kx_start, kx_end
           do ibib = 1, nbnb(kx,jpm) !! n,n' pair band index loop
             ! n1b, n2b has opposite meaning (occ <-> unocc) in jpm =2
             it = n1b(ibib,kx,jpm)  !index for n  for q   ! n1b(ibib,k,jpm) = n :band index for k (occupied),
             itp = n2b(ibib,kx,jpm) !index for n' for q+k ! n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
             do iw = ihw(ibib,kx,jpm), ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
               nttp(iw) = nttp(iw) + 1
               ittp = nttp(iw)
               itw(ittp,iw) = it
               itpw(ittp,iw) = itp
               ik(ittp,iw) = kx
               whwc(ittp,iw) = whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
             enddo
           enddo
         enddo
         allocate(zw(nttp_max,nnwf), wzw(nttp_max,nnwf))
         do iw = 1, nwhis
           if (nttp(iw) < 1) cycle
           do concurrent(inwf = 1:nnwf, ittp = 1:nttp(iw))
             iwf = mlo_pairs(inwf,1)
             jwf = mlo_pairs(inwf,2)
             it = itw(ittp,iw)
             itp = itpw(ittp,iw)
             kx = ik(ittp,iw)
             !12: Dual
             ! wzw(ittp,inwf) = whwc(ittp,iw)*dconjg(ov_evc_w2_kx(iwf,itp,kx))*ov_evc_w1_kx(jwf,it,kx)
             !  zw(ittp,inwf) =               dconjg(   evc_w2_kx(iwf,itp,kx))*   evc_w1_kx(jwf,it,kx)
             !13: Dual
             ! wzw(ittp,inwf) = whwc(ittp,iw)*dconjg(   evc_w2_kx(iwf,itp,kx))*ov_evc_w1_kx(jwf,it,kx)
             !  zw(ittp,inwf) =               dconjg(ov_evc_w2_kx(iwf,itp,kx))*   evc_w1_kx(jwf,it,kx)
             zw(ittp, inwf) = dconjg(evc_w2_kx(iwf,itp,kx))*evc_w1_kx(jwf,it,kx)
             wzw(ittp,inwf) = whwc(ittp,iw)*zw(ittp,inwf)
           enddo
           istat = zmm(zw, wzw, kmat(1,1,iw*(3-2*jpm)), nnwf, nnwf, nttp(iw), opA=m_op_C, beta=(1d0,0d0), ldA=nttp_max, ldB=nttp_max)
         enddo
         deallocate(nttp, itw, itpw, whwc, zw, wzw, ik)
        enddo jpmloop
        deallocate(evc_w1_kx, evc_w2_kx)
        if(gettetwt_split) call tetdeallocate()
      enddo kxblock_loop
      if(.not.gettetwt_split) call tetdeallocate()      ! --> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
      if(negative_cut) kmat(:,:,-nwhis:-1) = 0d0
      if(ipr) call writemem('mlo_magnon start dpsion')
      mpi_k_accumulate: block
        use m_mpi,only: MPI__reduceSum
        do jpm=1, npm
          do iw=1, nwhis
            call MPI__reduceSum(0, kmat(1,1,iw*(3-2*jpm)), nnwf*nnwf, communicator=comm_k)
          enddo
        enddo
      endblock mpi_k_accumulate
      if(mpi__root_k) then
        call dpsion_init(realomega, imagomega, .false.)
        call dpsion_chiq(realomega, imagomega, .false., kmat, zxqi, nnwf, nnwf, schi, 1, 1d99) !Inplace routine: kmat is overwritten by zxq
      endif
    endblock GETzxq

    if(geteta .and. (.not. mpi__root_k)) exit BIGiqloop
    if(.not. mpi__root_k) cycle BIGiqloop

    !Below lines are executed only by root of mpi__rank_k
    if(ipr) call writemem('mlo_magnon start getting R')
    if(associated(zxq)) nullify(zxq)
    zxq(1:,1:,nw_i:) => kmat(1:nnwf,1:nnwf,nw_i:nw)
    where(abs(dimag(zxq))<1d-15) zxq = dreal(zxq) ! threshold for Im[K] (zxq)

    formfactor = get_formfactor_q(q, isf, spinflip=.true.)
    ovlppair = get_ovlppair_q(q, isf, spinflip=.true.)
    ovlppair_inv = ovlppair
    istat = zminv(ovlppair_inv, n=nnwf)

    IfGetEta: if(geteta) then
      block
        use m_intg, only: intg_trapezoidal_nonuniform
        integer :: iunit, nev
        complex(8) ::wkmat(1:nnwf,1:nnwf), chi0(nnwf,nnwf), eval(nnwf), chiH(nnwf,nnwf)
        real(8) :: www, evl(nnwf), evlw(nw_i:nw)
        chi0(:,:) = zxq(:,:,0)
        istat = zmm(scrw, chi0(:,:), wkmat, nnwf, nnwf, nnwf)
        write(stdo,ftox) 'sum wkmat', sum(wkmat)
        istat = zgev(wkmat, n=nnwf, evl=eval)
        eta = -1d0/maxval(abs(eval))
        write(stdo,ftox) 'sum evl', sum(eval)
        write(stdo,ftox) "now eigenvalue WK",eval(1),"is inversed"
        write(stdo,ftox) "wkmat calculated eta:", eta !negative value
        open(newunit=iunit,file='__EtaMagnon',status='replace',form='unformatted',action='write')
        write(iunit) eta
        write(iunit) rho(:,:,:)
        close(iunit)
        open(newunit=iunit, file='Kpmdiag_q0.dat', status='replace', action='write')
        write(iunit,ftox) "# iw omega(eV) Tr K TrdiagK"
        do inwf=1,nnwf
          write(stdo,ftox) 'inwf, formfactor, diag_ovlppair', inwf, formfactor(inwf), ovlppair(inwf,inwf)
        enddo
        do iw = nw_i,nw
          www = merge(-freq_r(-iw),freq_r(iw),iw<0)
          chi0(:,:) = zxq(:,:,iw)
          chiH = (chi0 - transpose(dconjg(chi0)))*0.5d0/img
          istat = zhgv_lindep(chiH, ovlppair, n=nnwf, evl=evl, nev=nev)
          evlw(iw) = sum(evl(1:nev))
          write(iunit,"(e14.6,30e17.9)") www*hartree,  dot_product(formfactor, matmul(chi0(:,:), formfactor)), &
                                        trace_onsite(matmul(ovlppair_inv, chi0(:,:)))/hartree, evlw(iw)
        enddo
        close(iunit)
        write(stdo,ftox) 'Int dw sum of evl: Im K', intg_trapezoidal_nonuniform(freq, evlw(nw_i:nw))
        write(stdo,*) sum(rho(:,:,1)), sum(rho(:,:,2)), sum(rho(:,:,1) -rho(:,:,2))
        do i = 1, nwf
          write(stdo,*) rho(i,i,1), rho(i,i,2), rho(i,i,1) - rho(i,i,2)
        enddo
      endblock
      exit Bigiqloop
    endif IfGetEta

    OmegaLoop: block
      complex(8) :: r_tr(nw_i:nw), k_tr(nw_i:nw), jq_w_full(nw_i:nw,1:nnwf), &
                    k_tr_onsite(nw_i:nw), r_tr_onsite(nw_i:nw), &
                    jq_w_site(nw_i:nw,1:nsite),  wkmat(nnwf,nnwf), rmat(nnwf,nnwf), &
                    chi0(nnwf,nnwf), rmat_site(nsite,nsite,nw_i:nw), kmat_site(nsite,nsite,nw_i:nw), &
                    r_uovlp(nw_i:nw), k_uovlp(nw_i:nw), chiH(nnwf,nnwf), OchiH(nnwf,nnwf)
      real(8) :: evl(nnwf)
      real(8) :: spectrum_k(nw_i:nw), spectrum_r(nw_i:nw)
      integer :: nev
      complex(8) :: tmp(nnwf,nnwf)
      iwloop: do iw = nw_i, nw
        chi0(:,:) = zxq(:,:,iw)
        istat = zmm(scrw, chi0(:,:), wkmat, nnwf, nnwf, nnwf, alpha=dcmplx(eta,0d0)) !wkmat = etaWK
        wkmat(1:nnwf,1:nnwf) = imat(1:nnwf,1:nnwf) - wkmat(1:nnwf,1:nnwf) ! wkamt = 1 - etaWK
        istat = zminv(wkmat, n=nnwf) ! wkmat = (1- eta WK)^-1
        istat = zmm(chi0, wkmat, rmat, nnwf, nnwf, nnwf) !rmat = K (1-eta WK)^-1

        !Retarted
        chi0(:,:) = merge(conjg(transpose(chi0)), chi0, iw <0)
        rmat(:,:) = merge(conjg(transpose(rmat)), rmat, iw <0)

        chiH = -(chi0 - transpose(conjg(chi0)))*0.5d0*img
        istat = zhgv_lindep(chiH, ovlppair, n=nnwf, evl=evl, nev=nev, nkeep=nwf)
        spectrum_k(iw) = sum(evl(1:nev))

        ! OchiH = matmul(chi0, transpose(ovlppair))
        ! OchiH = -(OchiH- transpose(conjg(OchiH)))*0.5d0*img
        ! istat = zhev(OchiH, n=nnwf, evl=evl)
        ! spectrum_k(iw) = sum(evl(1:nnwf))

        chiH = -(rmat- transpose(conjg(rmat)))*0.5d0*img
        istat = zhgv_lindep(chiH, ovlppair, n=nnwf, evl=evl, nev=nev, nkeep=nwf)
        spectrum_r(iw) = sum(evl(1:nev))

        ! OchiH = matmul(rmat, transpose(ovlppair))
        ! OchiH = -(OchiH- transpose(conjg(OchiH)))*0.5d0*img
        ! istat = zhev(OchiH, n=nnwf, evl=evl)
        ! spectrum_r(iw) = sum(evl(1:nnwf))

        k_tr(iw) = trace2(matmul(ovlppair,chi0))
        r_tr(iw) = trace2(matmul(ovlppair,rmat))

        k_tr_onsite(iw) = trace_onsite(chi0)
        r_tr_onsite(iw) = trace_onsite(rmat)

        k_uovlp(iw) = dot_product(formfactor, matmul(chi0, formfactor))
        r_uovlp(iw) = dot_product(formfactor, matmul(rmat, formfactor))
        ! kmat_site(1,1,iw) = sum(ovlppair*chi0)
        ! rmat_site(1,1,iw) = sum(ovlppair*rmat)
        kmat_site(:,:,iw) = contract_to_site(chi0)
        rmat_site(:,:,iw) = contract_to_site(rmat)

        istat = zminv(rmat(:,:), n=nnwf)
        do concurrent(inwf=1:nnwf,jnwf=1:nnwf)
          rmat(inwf,jnwf) = sz(inwf)*rmat(inwf,jnwf)
        enddo
        istat = zgev(rmat(:,:), n=nnwf, evl=jq_w_full(iw,:))
      enddo iwloop

      CalcJqSite:block
        use m_intg, only: intg_trapezoidal_nonuniform
        use m_intg, only: intg_pade_nonuniform
        real(8) :: sz_kmat(nsite), sz_rmat(nsite), correction_factor
        integer :: isite, isite1, isite2
        write(stdo,ftox) '     q Int ev SK/pi, Int ev SR/pi:',q, &
           -intg_pade_nonuniform(freq, spectrum_k(nw_i:nw))/pi, -intg_pade_nonuniform(freq, spectrum_r(nw_i:nw))/pi
        write(stdo,ftox) '     q Int tr SK/pi, Int tr SR/pi:',q, &
           -intg_pade_nonuniform(freq, imag(k_tr(nw_i:nw)))/pi, -intg_pade_nonuniform(freq, imag(r_tr(nw_i:nw)))/pi
        write(stdo,ftox) '     q Int ff K/pi, Int  ff R/pi:',q, &
           -intg_pade_nonuniform(freq, imag(k_uovlp(nw_i:nw)))/pi, -intg_pade_nonuniform(freq, imag(r_uovlp(nw_i:nw)))/pi
        do isite = 1, nsite
          sz_kmat(isite) = -intg_pade_nonuniform(freq, imag(kmat_site(isite,isite,nw_i:nw)))/pi
          sz_rmat(isite) = -intg_pade_nonuniform(freq, imag(rmat_site(isite,isite,nw_i:nw)))/pi
          write(stdo,ftox) 'Site q Int K/pi, Int R/pi:',q, isite, sz_kmat(isite), sz_rmat(isite), trace(ovlppair)
        enddo
        do iw = nw_i, nw
          istat = zminv(rmat_site(:,:,iw), n=nsite)
          do concurrent(isite1=1:nsite,isite2=1:nsite)
            rmat_site(isite1,isite2,iw) = sz_site(isite1)*rmat_site(isite1,isite2,iw)
          enddo
          istat = zgev(rmat_site(:,:,iw), n=nsite, evl=jq_w_site(iw,:))
        enddo
      endblock CalcJqSite

      SaveBufferFile:block
        type(mpiio_buf) :: buf
        call buf_put(buf, k_uovlp);     call buf_put(buf, k_tr)
        call buf_put(buf, k_tr_onsite)
        call buf_put(buf, r_uovlp);     call buf_put(buf, r_tr)
        call buf_put(buf, r_tr_onsite)
        call buf_put(buf, spectrum_k)
        call buf_put(buf, spectrum_r)
        istat = writem_buf(file_tr_kr, rec=iq, buf=buf)
        istat = writem(file_jq_full, rec=iq, data=jq_w_full(nw_i:nw,1:nnwf))
        istat = writem(file_jq_site, rec=iq, data=jq_w_site(nw_i:nw,1:nsite))
      endblock SaveBufferFile
    endblock OmegaLoop

    if(ipr) call writemem('mlo_magnon end iq='//trim(charext(iq)))
  enddo BIGiqloop
  deallocate(formfactor, ovlppair, ovlppair_inv)

  if(geteta) call rx0( ' OK! mlo_magnon get eta')
  call mpi_barrier(comm, istat)

  if(ipr) write(stdo,"('eta for 1-eta*WK:',f13.8)") eta
  ReformatOutputFilesForDOS:if(mpi__root .and. calcdos) then
    block
      use m_bz_integ, only: ibz_integ
      integer :: file_tr_kpm_out, file_tr_rpm_out
      complex(8) :: dos_tr_kpm, dos_tr_onsite_kpm, dos_tr_rpm, dos_tr_onsite_rpm
      complex(8) :: k_uovlp_ibz(nqibz_dos,nw_i:nw), k_tr_ibz(nqibz_dos,nw_i:nw), k_tr_onsite_ibz(nqibz_dos,nw_i:nw), &
                    r_uovlp_ibz(nqibz_dos,nw_i:nw), r_tr_ibz(nqibz_dos,nw_i:nw), r_tr_onsite_ibz(nqibz_dos,nw_i:nw), &
                    jq_site_ibz(nqibz_dos,nw_i:nw,nsite)
      real(8) :: www, omega, tetra_vol
      integer :: tetra_nodes(4,nteti_dos), tetra_weight(nteti_dos)
      open(newunit=file_tr_rpm_out, file='TrRpm.dos', status='replace', form='formatted', action='write')
      open(newunit=file_tr_kpm_out, file='TrKpm.dos', status='replace', form='formatted', action='write')
      write(file_tr_kpm_out, '(A)')' # omega(eV) Real_Tr_K/eV Imag_Tr_K/eV !omega=0 is commented'
      write(file_tr_rpm_out, '(A)')' # omega(eV) Real_Tr_R/eV Imag_Tr_R/eV !omega=0 is commented'
      tetra_nodes(1:4,1:nteti_dos) = idteti_dos(1:4,1:nteti_dos)
      tetra_weight(1:nteti_dos) = idteti_dos(0,1:nteti_dos)
      do iq=1, nqcalc
        ReadBufferFileDOS:block
          type(mpiio_buf) :: buf
          istat = readm_buf(file_tr_kr, rec=iq, buf=buf)
          call buf_get(buf, k_uovlp_ibz(iq,nw_i:nw));     call buf_get(buf, k_tr_ibz(iq,nw_i:nw))
          call buf_get(buf, k_tr_onsite_ibz(iq,nw_i:nw))
          call buf_get(buf, r_uovlp_ibz(iq,nw_i:nw));     call buf_get(buf, r_tr_ibz(iq,nw_i:nw))
          call buf_get(buf, r_tr_onsite_ibz(iq,nw_i:nw))
          istat = readm(file_jq_site, rec=iq, data=jq_site_ibz(iq,nw_i:nw,1:nsite))
        endblock ReadBufferFileDOS
      enddo

      tetra_vol = 1d0/ntetf_dos ! ntetf was =6*n1*n2*n3
      do iw = nw_i, nw
        dos_tr_kpm = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, k_tr_ibz(:,iw))
        dos_tr_rpm = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, r_tr_ibz(:,iw))
        dos_tr_onsite_kpm = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, k_tr_onsite_ibz(:,iw))
        dos_tr_onsite_rpm = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, r_tr_onsite_ibz(:,iw))
        www = merge(-freq_r(-iw),freq_r(iw),iw<0)
        omega = www*hartree
        if(iw==0) then
          write(file_tr_kpm_out, "(A)", advance = "no") "#"
          write(file_tr_rpm_out, "(A)", advance = "no") "#"
        endif
        write(file_tr_kpm_out, "(e14.6,4e17.9)") omega, dos_tr_kpm/hartree, dos_tr_onsite_kpm/hartree
        write(file_tr_rpm_out, "(e14.6,4e17.9)") omega, dos_tr_rpm/hartree, dos_tr_onsite_rpm/hartree
      enddo
      close(file_tr_rpm_out)
      close(file_tr_kpm_out)
    endblock
  endif ReformatOutputFilesForDOS

  ReformatOutputFiles:if(mpi__root .and. .not.calcdos) then
    block
      use m_read_bzdata, only: epslgroup
      use m_sort, only: sort_index, lower_bound, upper_bound
      integer :: epslgroup_old, file_tr_kpm_out, file_tr_rpm_out, file_jq_out, file_jq_site_out, file_jq_full_out
      integer :: file_spec_out
      character(3) :: charnum3
      real(8) :: q_old(3), q_position, dq, omega, www
      real(8), parameter :: jq_cut_min =  1d-3, jq_cut_max=1
      integer :: idx_min, idx_max
      logical :: opened
      integer, allocatable, target :: idx_sort(:)
      complex(8) :: r_tr(nw_i:nw), k_tr(nw_i:nw), jq_w_site(nw_i:nw,1:nsite), &
                    jq_w_full(nw_i:nw,1:nnwf), r_tr_onsite(nw_i:nw), k_tr_onsite(nw_i:nw), r_uovlp(nw_i:nw), k_uovlp(nw_i:nw)
      real(8) :: k_spec(nw_i:nw), r_spec(nw_i:nw)
      epslgroup_old = -1 !epslgroup starts from 1
      q_old(:) = qibze(:,1)
      q_position = 0d0
      do iq=1, nqcalc
        q(:) = qibze(:,iq)
        if(epslgroup(iq) /= epslgroup_old) then
          if(any(q_old(:) /= q(:))) q_old(:) = q(:) + ([0.05d0, 0d0, 0d0])
          inquire(unit=file_tr_kpm_out, opened=opened)
          if(opened) close(file_tr_kpm_out)
          inquire(unit=file_tr_rpm_out, opened=opened)
          if(opened) close(file_tr_rpm_out)
          inquire(unit=file_jq_out, opened=opened)
          if(opened) close(file_jq_out)
          inquire(unit=file_jq_site_out, opened=opened)
          if(opened) close(file_jq_site_out)
          open(newunit=file_tr_kpm_out, file='TrKpm.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=file_tr_rpm_out, file='TrRpm.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=file_spec_out, file='Spec.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=file_jq_site_out, file='JqSite.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=file_jq_full_out, file='JqFull.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          write(file_tr_kpm_out, '(A)')' # qx qy qz q_pos omega(eV) Real_Tr_K/eV Imag_Tr_K/eV'
          write(file_tr_rpm_out, '(A)')' # qx qy qz q_pos omega(eV) Real_Tr_R/eV Imag_Tr_R/eV'
        endif
        ReadBufferFile:block
          type(mpiio_buf) :: buf
          istat = readm_buf(file_tr_kr, rec=iq, buf=buf)
          call buf_get(buf, k_uovlp);     call buf_get(buf, k_tr)
          call buf_get(buf, k_tr_onsite)
          call buf_get(buf, r_uovlp);     call buf_get(buf, r_tr)
          call buf_get(buf, r_tr_onsite)
          call buf_get(buf, k_spec)
          call buf_get(buf, r_spec)
          istat = readm(file_jq_site, rec=iq, data=jq_w_site(:,:))
          istat = readm(file_jq_full, rec=iq, data=jq_w_full(:,:))
        endblock ReadBufferFile
        dq = sqrt(sum((q(:)-q_old(:))**2))
        q_position = q_position + dq
        do iw = nw_i, nw
          www = merge(-freq_r(-iw),freq_r(iw),iw<0)
          omega = www*hartree
          write(file_tr_kpm_out, "(4f9.5,e14.6,8e17.9)") q(1:3), q_position, omega, k_uovlp(iw)/hartree, k_tr(iw)/hartree, k_tr_onsite(iw)/hartree
          write(file_tr_rpm_out, "(4f9.5,e14.6,8e17.9)") q(1:3), q_position, omega, r_uovlp(iw)/hartree, r_tr(iw)/hartree, r_tr_onsite(iw)/hartree
          write(file_spec_out, "(4f9.5,e14.6,2e17.9)") q(1:3), q_position, omega, k_spec(iw)/hartree, r_spec(iw)/hartree
          idx_sort = sort_index(abs(jq_w_site(iw,:)))
          write(file_jq_site_out, "(4f9.5,e14.6,20e17.9)") q(1:3), q_position, omega, ((jq_w_site(iw,idx_sort(i))*hartree),i=1,nsite)
          idx_sort = sort_index(abs(jq_w_full(iw,:)))
          idx_min = lower_bound(abs(jq_w_full(iw,:)), jq_cut_min, idx_sort)
          idx_max = upper_bound(abs(jq_w_full(iw,:)), jq_cut_max, idx_sort)
          write(file_jq_full_out, "(4f9.5,e14.6,81e17.9)") q(1:3), q_position, omega, ((dble(jq_w_full(iw,idx_sort(i)))*hartree),i=idx_min,idx_max)
        enddo
        write(file_tr_kpm_out, *)
        write(file_tr_rpm_out, *)
        write(file_jq_site_out, *)
        write(file_jq_full_out, *)
        write(file_spec_out, *)
        epslgroup_old = epslgroup(iq)
        q_old = q(:)
      enddo
      inquire(unit=file_tr_kpm_out, opened=opened)
      if(opened) close(file_tr_kpm_out)
      inquire(unit=file_tr_rpm_out, opened=opened)
      if(opened) close(file_tr_rpm_out)
    endblock
  endif ReformatOutputFiles
  istat = closem(file_tr_kr)
  istat = closem(file_jq_site)
  call rx0( ' OK! mlo_magnon mode')
END subroutine mlo_magnon
end module m_mlo_magnon
