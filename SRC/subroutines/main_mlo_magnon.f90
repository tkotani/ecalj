!>  Calculate Chi^+-, spin susceptibility. 
module m_mlo_magnon 
  use m_cmdopt_registry, only: c0_dos, c0_geteta, c2_nk, c2_sp1, c2_sp2
  implicit none
  public :: mlo_magnon
  contains
subroutine mlo_magnon() bind(C)
  use m_mlo_ham, only: read_ham_rs, calc_ham_eigen, nwf => ndimMTO, nsite, ib_tableM, nspx
  use m_mlo_scrw, only: nnwf, scrw, mlo_pairs, trace, nnwf_init, scrw_init,  &
                        pair_site, pair_lorb, nnwf_mask, nnwf2_mask, trace2
  use m_mlo_formfactor, only: get_formfactor_q
  use m_mlo_ovlppair, only: get_ovlppair_q
  use m_HamPMT,only: ReadHamPMTInfo
  use m_ReadEfermi, only: readefermi
  use m_read_bzdata, only: nqbz, qbz
  use m_genallcf_v3, only: genallcf_v3; use m_struct_from_lmf, only: nspin
  use m_keyvalue, only: getkeyvalue
  use m_GWinput, only: gwinput_init, gwinput_loaded, &
                       tg_magnon_w_onsite_dddd => magnon_w_onsite_dddd, &
                       tg_magnon_delta         => magnon_delta, &
                       tg_magnon_delta_dos     => magnon_delta_dos, &
                       tg_HistBin_ratio        => HistBin_ratio, &
                       tg_HistBin_dw           => HistBin_dw, &
                       tg_magnon_HistBin_ratio => magnon_HistBin_ratio, &
                       tg_magnon_HistBin_dw    => magnon_HistBin_dw, &
                       tg_magnon_negative_cut  => magnon_negative_cut
  use m_freq, only: getfreq, freq_r, nwhis, nw_i, nw, npm
  use m_tetwt, only: tetdeallocate, gettetwt, whw, ihw, nhw, jhw, n1b, n2b, nbnb
  use m_readgwinput, only: ReadGWinputKeys
  use m_lgunit, only: m_lgunit_init, stdo
  use m_dpsion, only: dpsion_init, dpsion_chiq => dpsion_chiq_h
  use m_mpi, only: MPI__Initialize, MPI__consoleout, MPI__InitQgroups, MPI__SplitXq
  use m_mpi, only: mpi__rank, mpi__size, mpi__root, comm, comm_k => comm_k_xq, &
                   mpi__rank_k => mpi__rank_k_xq, mpi__size_k => mpi__size_k_xq, &
                   mpi__root_k => mpi__root_k_xq, ipr
  use m_mpiio, only: openm, closem, writem, readm, mpiio_buf, buf_put, buf_get, writem_buf, readm_buf, buf_reset
  use m_blas, only: m_op_C, zmm => zmm_h, int_split
  use m_lapack, only: zminv => zminv_h, zhev => zhev_h, zgev => zgev_h,  zggv => zggv_h
  use m_mem, only: writemem
  use m_sort, only: sort_index, lower_bound, upper_bound
  use m_ftox, only: ftox
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  xxx removed--> eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  integer:: iwf, jwf, inwf, jnwf
  integer :: file_magnon
  integer:: iqxini, iqxend, i, iw, iq, kx, istat, nqcalc
  real(8):: q(3), rydberg, hartree, delta, eta, delta_dos, freq_ratio, freq_dw
  real(8), allocatable:: qibze(:,:)
  complex(8), pointer:: zxq(:,:,:) => null()
  complex(8), allocatable, target :: kmat(:,:,:)
  complex(8), allocatable:: imat(:,:)
  complex(8), parameter :: img=(0d0,1d0)
  integer :: isp1, isp2, is, isf
  logical:: realomega, imagomega, epsmode
  logical, allocatable :: mpi__task(:)
  character(8):: charext
  character(len=128) :: msg
  real(8), parameter :: pi = 4d0*datan(1d0), eta_default =1d0
  logical, parameter :: nnwf_size_reduction = .true.
  logical :: w_onsite_dddd, geteta, negative_cut, ganmma_only, gettetwt_split, calcdos
  !For dos calculation
  integer :: nqibz_dos, ntetf_dos, nteti_dos
  integer, allocatable :: idteti_dos(:,:)
  real(8), allocatable :: qibz_dos(:,:), rho(:,:,:), sz(:), sz_site(:)
  real(8), allocatable :: freq(:)

  hartree = 2d0*rydberg()
  isp1 = 2; isp2 = 1  ! default DNUP
  if (c2_sp1 >= 0) isp1 = c2_sp1
  if (c2_sp2 >= 0) isp2 = c2_sp2
  is = isp2; isf = isp1
  geteta = c0_geteta
  calcdos = c0_dos
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
  call gwinput_init()
  if (gwinput_loaded) then
     w_onsite_dddd = tg_magnon_w_onsite_dddd
     delta         = tg_magnon_delta
     delta_dos     = tg_magnon_delta_dos
     freq_ratio    = tg_HistBin_ratio
     freq_dw       = tg_HistBin_dw
     if (tg_magnon_HistBin_ratio /= 1.03d0) freq_ratio = tg_magnon_HistBin_ratio
     if (tg_magnon_HistBin_dw    /= 1d-5)   freq_dw    = tg_magnon_HistBin_dw
     negative_cut  = tg_magnon_negative_cut
  else
     call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!     call getkeyvalue("GWinput","magnon_w_onsite_dddd",w_onsite_dddd,default=.true.)
!     call getkeyvalue("GWinput","magnon_delta", delta, default=0d0) !1d-6 for Insulator case
!     call getkeyvalue("GWinput","magnon_delta_dos", delta_dos, default=1d-6)
!     call getkeyvalue("GWinput","HistBin_ratio",freq_ratio, default=1.03d0)
!     call getkeyvalue("GWinput","HistBin_dw",freq_dw, default=1d-5)
!     call getkeyvalue("GWinput","magnon_HistBin_ratio",freq_ratio, default=freq_ratio)
!     call getkeyvalue("GWinput","magnon_HistBin_dw", freq_dw, default=freq_dw)
!     call getkeyvalue("GWinput","magnon_negative_cut",negative_cut,default=.false.)
  endif
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
    nqcalc = iqxend-iqxini+1
    n_bpara = 1
    n_kpara = max(mpi__size/(n_bpara*nqcalc), 1)  !Default setting of parallelization. b-parallel is 1.
    if (c2_nk >= 0) n_kpara = c2_nk
    worker_inQtask = n_bpara * n_kpara
    ! gettetwt_split = n_kpara > 2
    gettetwt_split = .true. ! We always split gettetwt by k to prevent memory exhaustion.
    if(ipr) write(stdo,ftox) 'MPI: worker_inQtask', worker_inQtask
    allocate(mpi__task(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=iqxini,iqxend)])
    if(ipr) write(stdo,ftox) 'mpi_rank',mpi__rank,'mpi__Qtask=',mpi__task
    call MPI__InitQgroups(worker_inQtask)
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
    call ReadHamPMTInfo()  ! Read info from PMTHamiltonianInfo (lattice structures and index of basis).
    call read_ham_rs()
    call nnwf_init(nnwf_size_reduction) !set nnwf ~ # of RiRj (onsite_approx = .true.), RiR'j (onsite_approx = .flase. ), wan_pair_index
    if(ipr) write(stdo,ftox) '# nwf, nnwf:', nwf, nnwf
    call scrw_init(w_onsite_dddd, isp1, isp2, enforce_Hermite=.false.) !set scrw
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
    istat = openm(newunit=file_magnon, file='__MagnonData', recl=(nw-nw_i+1)*(16*(nsite*nsite+1)*2*2 + 8*(nsite+1)*2))
  endif SetTemporaryFiles

  allocate(imat(1:nnwf,1:nnwf),source=(0d0,0d0))
  forall(iwf=1:nnwf) imat(iwf,iwf) = 1d0 + img*delta !check the sign
  allocate(kmat(1:nnwf,1:nnwf,(1-npm)*nwhis:nwhis))
  BIGiqloop: do iq = iqxini,iqxend
    if(.NOT. MPI__task(iq)) cycle BIGiqloop
    q = qibze(:,iq)
    if(ipr) write(6,"('===== do : iq wibz(iq) q=',i6,f13.6,3f9.4,' ========')") iq,q !,wibz(iqlist(iq)),qshort !qq
    GETzxq: block ! zxq and zxqi are the main output after Hilbert transformation, ! zxqi is not used in hmagnon (imagomega=.false.)
      use m_mpi,only: MPI__AllreduceSumReal, MPI__AllreduceSum
      real(8) :: evkx_w1(nwf,nspx), evkx_w2(nwf,nspx) !dummy
      complex(8) :: zxqi(1,1,1), evc_w1(nwf,nwf), evc_w2(nwf,nwf)
      complex(8) :: ov_evc_w1(nwf,nwf), ov_evc_w2(nwf,nwf)
      complex(8), allocatable :: evc_w1_kx(:,:,:), evc_w2_kx(:,:,:)
      integer, allocatable :: nttp(:),  itw(:,:), itpw(:,:), ik(:,:)
      integer :: nttp_max, ittp, jpm, it, itp, ibib, isdummy, kx_ini, kx_fin, kx_start,kx_end, kx_num
      real(8), allocatable :: whwc(:,:), ev_w1(:,:), ev_w2(:,:)
      integer, parameter:: nkblock = 1024
      complex(8), allocatable :: zw(:,:), wzw(:,:)
      allocate(ev_w1(nwf,nqbz), ev_w2(nwf,nqbz), source=0d0)

      if(ipr) call writemem('mlo_magnon start gettetwt')
      call int_split(nqbz, mpi__size_k, mpi__rank_k, kx_ini, kx_fin, kx_num)
      CalcEigenEnergy: do kx = kx_ini, kx_fin !!! ev_w1, ev_w2 unit: [Ry]
        call calc_ham_eigen( is,   qbz(:,kx), evkx_w1(:, is), evec=evc_w1, ovlp_evec=ov_evc_w1)
        call calc_ham_eigen(isf, q+qbz(:,kx), evkx_w2(:,isf), evec=evc_w2, ovlp_evec=ov_evc_w2)
        ev_w1(:,kx) = evkx_w1(: ,is)
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
          call calc_ham_eigen( is,   qbz(:,kx), evkx_w1(:, is), evec=evc_w1_kx(:,:,kx)) !evkx_w1  is dummy
          call calc_ham_eigen(isf, q+qbz(:,kx), evkx_w2(:,isf), evec=evc_w2_kx(:,:,kx)) !evkx_w2  is dummy
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
        call dpsion_chiq(realomega, imagomega, .false., kmat, zxqi, nnwf, nnwf, 1d0, 1, 1d99) !Inplace routine: kmat is overwritten by zxq
      endif
    endblock GETzxq

    if(geteta .and. (.not. mpi__root_k)) exit BIGiqloop
    if(.not. mpi__root_k) cycle BIGiqloop

    !Below lines are executed only by root of mpi__rank_k
    if(ipr) call writemem('mlo_magnon start getting R')
    if(associated(zxq)) nullify(zxq)
    zxq(1:,1:,nw_i:) => kmat(1:nnwf,1:nnwf,nw_i:nw)
    zxq(:,:,:) = -zxq(:,:,:)
    where(abs(dimag(zxq))<1d-15) zxq = dreal(zxq) ! threshold for Im[K] (zxq)

    IfGetEta: if(geteta) then
      block
        integer :: iunit
        complex(8) :: wkmat(1:nnwf,1:nnwf), chi0(nnwf,nnwf), eval(nnwf)
        chi0(:,:) = zxq(:,:,0)
        istat = zmm(scrw, chi0(:,:), wkmat, nnwf, nnwf, nnwf)
        write(stdo,ftox) 'sum wkmat', sum(wkmat)
        istat = zgev(wkmat, n=nnwf, evl=eval)
        eta = 1d0/maxval(abs(eval))
        write(stdo,ftox) 'sum evl', sum(eval)
        write(stdo,ftox) "now eigenvalue WK",eval(1),"is inversed"
        write(stdo,ftox) "wkmat calculated eta:", eta !negative value
        open(newunit=iunit,file='__EtaMagnon',status='replace',form='unformatted',action='write')
        write(iunit) eta
        write(iunit) rho(:,:,:)
        close(iunit)
        write(stdo,ftox) sum(rho(:,:,1)), sum(rho(:,:,2)), sum(rho(:,:,1) -rho(:,:,2))
        do i = 1, nwf
          write(stdo,ftox) rho(i,i,1), rho(i,i,2), rho(i,i,1) - rho(i,i,2)
        enddo
      endblock
      exit Bigiqloop
    endif IfGetEta

    OmegaLoop: block
      complex(8) :: wkmat(nnwf,nnwf), rmat(nnwf,nnwf), chi0(nnwf,nnwf), OchiH(nnwf,nnwf), weight_site(nnwf,nsite)
      real(8) :: evl(nnwf), k_spec(nw_i:nw), r_spec(nw_i:nw), k_spec_site(nw_i:nw,nsite), r_spec_site(nw_i:nw,nsite)
      complex(8) :: k_chi(nw_i:nw), r_chi(nw_i:nw), k_chi_site(nw_i:nw,nsite,nsite), r_chi_site(nw_i:nw,nsite,nsite)
      complex(8) :: k_dms(nw_i:nw), r_dms(nw_i:nw), k_dms_site(nw_i:nw,nsite,nsite), r_dms_site(nw_i:nw,nsite,nsite)
      complex(8) :: ovlppair(nnwf,nnwf), formfactor(nnwf), ovlp(nnwf), ovlp_site(nnwf, nsite), formfactor_site(nnwf, nsite)
      integer, allocatable, target :: idx_sort(:)
      integer :: isite, jsite

      formfactor = get_formfactor_q(q, isp1, isp2)
      ovlp = get_formfactor_q([0d0,0d0,0d0], isp1, isp2)
      ovlppair = get_ovlppair_q(q, isp1, isp2)
      do isite=1, nsite
        ovlp_site(:,isite) = merge(ovlp, (0d0,0d0), pair_site(:,1) == isite)
        formfactor_site(:,isite) = merge(formfactor, (0d0,0d0), pair_site(:,1) == isite)
      enddo
      iwloop: do iw = nw_i, nw
        chi0(:,:) = zxq(:,:,iw)
        istat = zmm(scrw, chi0(:,:), wkmat, nnwf, nnwf, nnwf, alpha=dcmplx(eta,0d0)) !wkmat = etaWK
        wkmat(1:nnwf,1:nnwf) = imat(1:nnwf,1:nnwf) - wkmat(1:nnwf,1:nnwf) ! wkamt = 1 - etaWK
        istat = zminv(wkmat, n=nnwf) ! wkmat = (1- eta WK)^-1
        istat = zmm(chi0, wkmat, rmat, nnwf, nnwf, nnwf) !rmat = K (1-eta WK)^-1

        !Retarted
        chi0(:,:) = merge(conjg(transpose(chi0)), chi0, iw <0)
        rmat(:,:) = merge(conjg(transpose(rmat)), rmat, iw <0)

        !overlap
        k_chi(iw) = dot_product(ovlp, matmul(chi0, ovlp))
        r_chi(iw) = dot_product(ovlp, matmul(rmat, ovlp))
        do concurrent(isite=1:nsite, jsite=1:nsite)
          k_chi_site(iw,isite,jsite) = dot_product(ovlp_site(:,isite), matmul(chi0, ovlp_site(:,jsite)))
          r_chi_site(iw,isite,jsite) = dot_product(ovlp_site(:,isite), matmul(rmat, ovlp_site(:,jsite)))
        enddo

        !formfactor
        k_dms(iw) = dot_product(formfactor, matmul(chi0, formfactor))
        r_dms(iw) = dot_product(formfactor, matmul(rmat, formfactor))
        do concurrent(isite=1:nsite, jsite=1:nsite)
          k_dms_site(iw,isite,jsite) = dot_product(formfactor_site(:,isite), matmul(chi0, formfactor_site(:,jsite)))
          r_dms_site(iw,isite,jsite) = dot_product(formfactor_site(:,isite), matmul(rmat, formfactor_site(:,jsite)))
        enddo

        !mode
        OchiH = matmul(chi0, ovlppair)
        OchiH = -(OchiH- transpose(conjg(OchiH)))*0.5d0*img
        istat = zhev(OchiH, n=nnwf, evl=evl)
        k_spec(iw) = sum(evl(1:nnwf))
        do concurrent(i=1:nnwf, isite=1:nsite)
          weight_site(i,isite) = sum(pack(OchiH(:,i), mask=(pair_site(:,1)==isite)))
        enddo
        do isite = 1, nsite
          do i = 1, nnwf
            k_spec_site(iw,isite) = k_spec_site(iw,isite) + evl(i)*sum(abs(OchiH(:,i))**2, mask=(pair_site(:,1)==isite))
          enddo
        enddo

        OchiH = matmul(rmat, ovlppair)
        OchiH = -(OchiH- transpose(conjg(OchiH)))*0.5d0*img
        istat = zhev(OchiH, n=nnwf, evl=evl)
        r_spec(iw) = sum(evl(1:nnwf))
        do concurrent(i=1:nnwf, isite=1:nsite)
          weight_site(i,isite) = sum(pack(OchiH(:,i), mask=(pair_site(:,1)==isite)))
        enddo
        do isite = 1, nsite
          do i = 1, nnwf
            r_spec_site(iw,isite) = r_spec_site(iw,isite) + evl(i)*sum(abs(OchiH(:,i))**2, mask=(pair_site(:,1)==isite))
          enddo
        enddo
      enddo iwloop

      SaveBufferFile:block
        type(mpiio_buf) :: buf
        call buf_put(buf, k_chi);      call buf_put(buf, r_chi)
        call buf_put(buf, k_chi_site); call buf_put(buf, r_chi_site)
        call buf_put(buf, k_dms);      call buf_put(buf, r_dms)
        call buf_put(buf, k_dms_site); call buf_put(buf, r_dms_site)
        call buf_put(buf, k_spec);     call buf_put(buf, r_spec)
        call buf_put(buf, k_spec_site); call buf_put(buf, r_spec_site)
        istat = writem_buf(file_magnon, rec=iq, buf=buf)
      endblock SaveBufferFile
    endblock OmegaLoop

    if(ipr) call writemem('mlo_magnon end iq='//trim(charext(iq)))
  enddo BIGiqloop

  if(geteta) call rx0( ' OK! mlo_magnon get eta')
  call mpi_barrier(comm, istat)

  if(ipr) write(stdo,"('eta for 1-eta*WK:',f13.8)") eta
  ! ReformatOutputFilesForDOS:if(mpi__root .and. calcdos) then
  !   block
  !     use m_bz_integ, only: ibz_integ
  !     integer :: file_dos_out
  !     real(8) :: dos_k_spec, dos_r_spec, dos_k_dms, dos_r_dms
  !     ! real(8) :: k_dms_ibz(nqibz_dos,nw_i:nw), k_spec_ibz(nqibz_dos,nw_i:nw), &
  !     !            r_dms_ibz(nqibz_dos,nw_i:nw), r_spec_ibz(nqibz_dos,nw_i:nw)
  !     complex(8) :: k_dms_ibz(nqibz_dos,nw_i:nw), k_spec_ibz(nqibz_dos,nw_i:nw), &
  !                r_dms_ibz(nqibz_dos,nw_i:nw), r_spec_ibz(nqibz_dos,nw_i:nw)
  !     real(8) :: www, omega, tetra_vol
  !     integer :: tetra_nodes(4,nteti_dos), tetra_weight(nteti_dos)
  !     open(newunit=file_dos_out, file='MagSpecDMS.dos', status='replace', form='formatted', action='write')
  !     write(file_dos_out, '(A)')' # omega(eV) Real_Tr_R/eV Imag_Tr_R/eV !omega=0 is commented'
  !     tetra_nodes(1:4,1:nteti_dos) = idteti_dos(1:4,1:nteti_dos)
  !     tetra_weight(1:nteti_dos) = idteti_dos(0,1:nteti_dos)
  !     do iq=1, nqcalc
  !       ReadBufferFileDOS:block
  !         type(mpiio_buf) :: buf
  !         istat = readm_buf(file_spec, rec=iq, buf=buf)
  !         call buf_get(buf, k_spec_ibz(iq,nw_i:nw))
  !         call buf_get(buf, r_spec_ibz(iq,nw_i:nw))
  !         istat = readm_buf(file_dms, rec=iq, buf=buf)
  !         call buf_get(buf, k_dms_ibz(iq,nw_i:nw))
  !         call buf_get(buf, r_dms_ibz(iq,nw_i:nw))
  !       endblock ReadBufferFileDOS
  !     enddo
  !
  !     tetra_vol = 1d0/ntetf_dos ! ntetf was =6*n1*n2*n3
  !     do iw = nw_i, nw
  !       dos_k_spec = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, k_spec_ibz(:,iw))
  !       dos_r_spec = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, r_spec_ibz(:,iw))
  !       dos_k_dms = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, k_dms_ibz(:,iw))
  !       dos_r_dms = ibz_integ(tetra_vol, nteti_dos, tetra_nodes, tetra_weight, qibz_dos, nqibz_dos, r_dms_ibz(:,iw))
  !       www = merge(-freq_r(-iw),freq_r(iw),iw<0)
  !       omega = www*hartree
  !       if(iw==0) then
  !         write(file_dos_out, "(A)", advance = "no") "#"
  !       endif
  !       write(file_dos_out, "(e14.6,4e17.9)") omega, dos_k_spec/hartree, dos_r_spec/hartree, dos_r_dms/hartree, dos_r_dms/hartree
  !     enddo
  !     close(file_dos_out)
  !   endblock
  ! endif ReformatOutputFilesForDOS

  ReformatOutputFiles:if(mpi__root .and. .not.calcdos) then
    block
      use m_read_bzdata, only: epslgroup
      use m_intg, only: intg_pade_nonuniform
      integer :: epslgroup_old
      integer :: out_file_spec, out_file_chi, out_file_dms, out_file_spec_site(nsite), out_file_chi_site(nsite), out_file_dms_site(nsite), out_file_inv_chi
      character(3) :: charnum3
      real(8) :: q_old(3), q_position, dq, omega, www
      integer :: idx_min, idx_max, isite, jsite
      integer, allocatable, target :: idx_sort(:)
      complex(8) :: r_dms(nw_i:nw), k_dms(nw_i:nw), k_dms_site(nw_i:nw,nsite,nsite), r_dms_site(nw_i:nw,nsite,nsite)
      complex(8) :: k_chi(nw_i:nw), r_chi(nw_i:nw), k_chi_site(nw_i:nw,nsite,nsite), r_chi_site(nw_i:nw,nsite,nsite)
      real(8) :: k_spec(nw_i:nw), r_spec(nw_i:nw), k_spec_site(nw_i:nw,nsite), r_spec_site(nw_i:nw,nsite)
      epslgroup_old = -1 !epslgroup starts from 1
      q_old(:) = qibze(:,1)
      q_position = 0d0
      write(stdo,'(A)') 'sum-check q: chi_k, chi_r, mds_k, mds_r, mode_k, mode_r'
      do iq=1, nqcalc
        q(:) = qibze(:,iq)

        OpenOutFiles:if(epslgroup(iq) /= epslgroup_old) then
          if(any(q_old(:) /= q(:))) q_old(:) = q(:) + ([0.05d0, 0d0, 0d0])
          open(newunit=out_file_spec, file='MagSpec.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=out_file_chi, file='MagSuscep.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=out_file_dms, file='DynMagSuscep.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          open(newunit=out_file_inv_chi, file='InvMagSuscep.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          do isite = 1, nsite
            open(newunit=out_file_spec_site(isite), file='MagSpecSite'//charnum3(isite)//'.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
            open(newunit=out_file_chi_site(isite), file='MagSuscepSite'//charnum3(isite)//'.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
            open(newunit=out_file_dms_site(isite), file='DynMagSuscepSite'//charnum3(isite)//'.syml'//charnum3(epslgroup(iq)), status='replace', form='formatted', action='write')
          enddo
        endif OpenOutFiles

        ReadBufferFile:block
          type(mpiio_buf) :: buf
          istat = readm_buf(file_magnon, rec=iq, buf=buf)
          call buf_get(buf, k_chi);      call buf_get(buf, r_chi)
          call buf_get(buf, k_chi_site); call buf_get(buf, r_chi_site)
          call buf_get(buf, k_dms);      call buf_get(buf, r_dms)
          call buf_get(buf, k_dms_site); call buf_get(buf, r_dms_site)
          call buf_get(buf, k_spec);     call buf_get(buf, r_spec)
          call buf_get(buf, k_spec_site); call buf_get(buf, r_spec_site)
        endblock ReadBufferFile

        dq = sqrt(sum((q(:)-q_old(:))**2))
        q_position = q_position + dq
        write(stdo,'(3f8.5,6f18.5)') q,&
           intg_pade_nonuniform(freq, imag(k_chi(nw_i:nw)))/pi, intg_pade_nonuniform(freq, imag(r_chi(nw_i:nw)))/pi ,&
           intg_pade_nonuniform(freq, imag(k_dms(nw_i:nw)))/pi, intg_pade_nonuniform(freq, imag(r_dms(nw_i:nw)))/pi ,&
           intg_pade_nonuniform(freq, k_spec(nw_i:nw))/pi, intg_pade_nonuniform(freq, r_spec(nw_i:nw))/pi

        GetLLGparameters:block
          complex(8) :: r_inv_site(nsite,nsite)
          character(40) :: fmt_inv
          write(fmt_inv,'("(4f9.5,e14.6,",I0,"e17.9)")') 2*nsite
          do iw = nw_i, nw
            r_inv_site(:,:) = r_chi_site(iw,:,:)
            istat = zminv(r_inv_site, n=nsite)
            www = merge(-freq_r(-iw),freq_r(iw),iw<0)
            omega = www*hartree
            write(out_file_inv_chi, fmt_inv) q(1:3), q_position, omega, (r_inv_site(isite,isite)*hartree, isite=1,nsite)
          enddo
          write(out_file_inv_chi, *)
        endblock GetLLGparameters

        WriteOutFiles:block
          do iw = nw_i, nw
            www = merge(-freq_r(-iw),freq_r(iw),iw<0)
            omega = www*hartree
            write(out_file_spec, "(4f9.5,e14.6,2e17.9)") q(1:3), q_position, omega, k_spec(iw)/hartree, r_spec(iw)/hartree
            write(out_file_chi,  "(4f9.5,e14.6,4e17.9)") q(1:3), q_position, omega, k_chi(iw)/hartree,  r_chi(iw)/hartree
            write(out_file_dms,  "(4f9.5,e14.6,4e17.9)") q(1:3), q_position, omega, k_dms(iw)/hartree,  r_dms(iw)/hartree
            do isite = 1, nsite
              write(out_file_spec_site(isite), "(4f9.5,e14.6,2e17.9)") q(1:3), q_position, omega, &
                k_spec_site(iw,isite)/hartree, r_spec_site(iw,isite)/hartree
              write(out_file_chi_site(isite),  "(4f9.5,e14.6,4e17.9)") q(1:3), q_position, omega, &
                k_chi_site(iw,isite,isite)/hartree, r_chi_site(iw,isite,isite)/hartree
              write(out_file_dms_site(isite),  "(4f9.5,e14.6,4e17.9)") q(1:3), q_position, omega, &
                k_dms_site(iw,isite,isite)/hartree, r_dms_site(iw,isite,isite)/hartree
            enddo
          enddo
          write(out_file_spec, *)
          write(out_file_chi, *)
          write(out_file_dms, *)
          do isite = 1, nsite
            write(out_file_spec_site(isite), *)
            write(out_file_chi_site(isite), *)
            write(out_file_dms_site(isite), *)
          enddo
        endblock WriteOutFiles

        epslgroup_old = epslgroup(iq)
        q_old = q(:)
      enddo
    endblock
  endif ReformatOutputFiles
  call rx0( ' OK! mlo_magnon mode')
END subroutine mlo_magnon
end module m_mlo_magnon
