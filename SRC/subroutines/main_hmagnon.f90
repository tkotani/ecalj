!>  Calculate Chi^+-, spin susceptibility. 
module m_hmagnon 
  contains
subroutine hmagnon() bind(C)
  use m_readwan,only: wan_readeval2, read_wandata, nwf, tr_mat_onsite, tr_mat_onsite_diag, &
                    & set_wan_nnwf, nnwf, set_wan_scrw, scrw, wan_pair_index
  use m_ReadEfermi, only: readefermi
  use m_read_bzdata, only: read_bzdata, nqbz, nqibz, qbz, wqt=>wt, epslgroup, nq0i, q0i, neps
  use m_genallcf_v3, only: genallcf_v3, nspin
  use m_keyvalue, only: getkeyvalue
  use m_freq, only: getfreq, freq_r, nwhis, nw_i, nw, npm
  use m_tetwt, only: tetdeallocate, gettetwt, whw, ihw, nhw, jhw, n1b, n2b, nbnb
  use m_readgwinput, only: ReadGWinputKeys
  use m_lgunit, only: m_lgunit_init, stdo
  use m_dpsion, only: dpsion_init, dpsion_chiq
  use m_mpi, only: MPI__Initialize, MPI__consoleout, MPI__SplitXq
  use m_mpi, only: mpi__rank, mpi__size, mpi__root ,comm, ipr, comm_k, mpi__rank_k, mpi__size_k, mpi__root_k
  use m_mpiio, only: openm, closem, writem, readm
  use m_blas, only: m_op_C, zmm => zmm_h
  use m_lapack, only: zminv => zminv_h
  use m_mem, only: writemem
  use m_ftox
  implicit none
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  xxx removed--> eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  integer:: iwf, jwf, inwf, kwf, lwf, ijwf, klwf 
  integer:: file_tr_kpm, file_tr_rpm, file_tr_diag_kpm, file_tr_diag_rpm
  integer:: iww, iqxini, iqxend, i, iw, iq, kx, ik, istat, nqcalc
  real(8):: q(3), omg2max, wemax, rydberg, hartree
  real(8), parameter:: schi=1d0, ua = 1d0, eta_default =-1d0
  real(8):: nms_delta, www
  real(8), allocatable:: qibze(:,:)
  complex(8), pointer:: zxq(:,:,:) => null()
  complex(8), allocatable, target :: kmat(:,:,:)
  complex(8), allocatable:: wkmat(:,:), rmat(:,:), r_tr(:), r_diag(:), k_tr(:), k_diag(:)
  complex(8), parameter :: img=(0d0,1d0)
  complex(8), allocatable::imat(:,:) !unit matrix for 1-WK
  logical:: cmdopt0
  logical:: realomega, imagomega, epsmode, wan, nms !, lhm, lsvd
  logical, allocatable :: mpi__task(:)
  character(8):: charext
  character(len=128) :: msg
  real(8) :: eta
  integer, parameter :: is=1, isf=2  !K_down up = Kpm
  real(8), parameter :: pi = 4d0*datan(1d0), znorm=-1d0*pi ! normalization of Im[K]:
  logical :: onsite_approx, w_onsite_dddd, geteta, negative_cut, ganmma_only
!!! q on symline
  
  ! MO cma mode is commented out 2025-12-06. cma mode is no longer maintained. For CMA mode, use old version
  ! logical:: cma_mode !cma_mode for Cu2MnAl only 2019/09/27
  ! real(8):: cma_up_shift, cma_dn_shift, cma_wshift
  ! integer(4):: cma_iwf_s, cma_iwf_e
  ! complex(8), allocatable:: scrw_original(:,:) !screening W

  hartree  = 2d0*rydberg()
  geteta = cmdopt0('--geteta')
  ganmma_only = geteta  !GammaPoint only calculation

  call m_lgunit_init()
  call MPI__Initialize()
  msg ='hmagnon'
  if(geteta) msg = trim(msg)//'_geteta_mode'
  call MPI__consoleout(trim(msg))
  call cputid(0)
  realomega = .true.
  imagomega = .false.
  epsmode   = .true.
  wan       = .true.
  call genallcf_v3(incwfx=0) !!incwfin=0 =>ForX0 for core in GWIN. in module m_genallcf_v3 Readin by genallcf. Set basic data for crystal
  if(nspin < 2) call rx(' hmagnon: nspin<2: not supported. exit.')
  !! Prof.Naraga said " write(6,*)'Timereversal=',Timereversal()" here caused a stop in ifort ver.1x.x. Why? May be a compilar bug, and fixed now.
  !! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
  !! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
  call read_BZDATA() !  !! Read electron gas mode or not.
  call ReadGWinputKeys() ! jun2020 new routint to read all inputs
  ! W is enforced as on site regardless onsite_approx, onsite_approx specifies whether nnwf is set as onsite or not.
  call getkeyvalue("GWinput","magnon_onsite_approximation",onsite_approx,default=.true.)
  call getkeyvalue("GWinput","magnon_w_onsite_dddd",w_onsite_dddd,default=.true.)
  call getkeyvalue("GWinput","nms",nms,default=.false.)  !!! For NiMnSb
  call getkeyvalue("GWinput","nms_delta",nms_delta,default=1d-6)
  call getkeyvalue("GWinput","negative_cut",negative_cut,default=.false.)
  write(6,*) "negative_cut",negative_cut

  SetQvecList: block
    if(mpi__root) then
      do i=1,nqbz
        if(i<10 .OR. i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
        if(i==10 .AND. nqbz>18) write(6,"('... ')")
      enddo
      write(6,*)' !!nqbz nqibz =',nqbz, nqibz
    endif
    if(ganmma_only) then
      allocate(qibze(3,1), source = 0d0)
      iqxini = 1
      iqxend = 1
    else
      write(6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      allocate(qibze(3,nq0i), source = q0i(1:3,1:nq0i))
      iqxini = 1
      iqxend = nq0i
    endif
    write(6,"('iqxini,iqxend ',2I8)") iqxini,iqxend
    do iq = iqxini,iqxend
      write(6,"('iq, qibze:',I8,3f9.4)") iq, qibze(:,iq)
    enddo
  endblock SetQvecList

  SetMPI_Rankdivider: block
    integer :: n_bpara, n_kpara, worker_inQtask
    integer, allocatable :: mpi__ranktab(:)
    logical:: cmdopt2
    character(20):: outs
    nqcalc = iqxend-iqxini+1
    n_bpara = 1
    n_kpara = max(mpi__size/(n_bpara*nqcalc), 1)  !Default setting of parallelization. b-parallel is 1.
    if(cmdopt2('--nk=', outs)) read(outs,*) n_kpara
    worker_inQtask = n_bpara * n_kpara
    write(stdo,ftox) 'MPI: worker_inQtask', worker_inQtask
    allocate(mpi__ranktab(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)*worker_inQtask           ,iq=iqxini,iqxend)])
    allocate(mpi__task(iqxini:iqxend),    source=[(mod(iq-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=iqxini,iqxend)])
    write(stdo,ftox) 'mpi_rank',mpi__rank,'mpi__Qtask=',mpi__task
    write(stdo,ftox) 'mpi_qrank', mpi__ranktab
    call MPI__SplitXq(n_bpara, n_kpara)
  endblock SetMPI_Rankdivider

  SetFreqencyMesh:block
    use m_freq, only: niw
    integer :: niw_in
    ! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
    wemax   = 5d0 !max value for plot
    omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
    !! NOTE: npmtwo=T sets npm=2   !! optional npmtwo is added aug2017   !! 20190604 Im[K]
    if( .NOT. imagomega) niw_in=1  !dummy
    call Getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw_in,ua, npmtwo=.true.)!,tetra
    if(mpi__root) write(6,"(' nw_i nw niw npm=',4i5)") nw_i,nw,niw,npm
  endblock SetFreqencyMesh

  call readefermi() !!! ef:     Fermi energy at 0 K

  SetWannierAndScreendCoulombData: block
    logical :: cmdopt2
    character(20):: Wtype, opts
    call read_wandata()    ! nwf, nsp_w,nqtt_w ! --- okumura Read dimensions of hamiltonian_wannier, spin, nqtt
    call set_wan_nnwf(onsite_approx) !set nnwf ~ # of RiRj (onsite_approx = .true.), RiR'j (onsite_approx = .flase. ), wan_pair_index
    Wtype = 'up' !options: up, down, up_down, down_up
    if(cmdopt2('--Wtype=', opts)) Wtype = trim(opts)
    call set_wan_scrw(onsite_approx, w_onsite_dddd, Wtype=Wtype) !set scrw
    if(mpi__root) write(stdo,ftox) '# nwf, nnwf:', nwf, nnwf
  endblock SetWannierAndScreendCoulombData

  ReadEta: if(.not. geteta) then
    block
      logical :: exist_etafile
      integer :: ios, iunit
      eta = eta_default
      inquire(file='__EtaMagnon',exist=exist_etafile)
      if(exist_etafile) then
        open(newunit=iunit, file='__EtaMagnon',status='old',form='formatted',action='read')
        read(iunit, *, iostat=ios) eta
        close(iunit)
        if(ios == 0) write(stdo, ftox) "# Read eta from __EtaMagnon"
        if(ios /= 0) eta = eta_default
      endif
      write(stdo, ftox) "# WK is scalled to  eta WK: eta =",eta
    endblock
  endif ReadEta

  SetTemporaryFiles: if(.not. geteta) then
    istat = openm(newunit=file_tr_kpm, file='__TrKpm', recl=(nw-nw_i+1)*16)
    istat = openm(newunit=file_tr_rpm, file='__TrRpm', recl=(nw-nw_i+1)*16)
    istat = openm(newunit=file_tr_diag_kpm, file='__TrDiagKpm', recl=(nw-nw_i+1)*16)
    istat = openm(newunit=file_tr_diag_rpm, file='__TrDiagRpm', recl=(nw-nw_i+1)*16)
  endif SetTemporaryFiles

  allocate(imat(1:nnwf,1:nnwf),source=(0d0,0d0))
  forall(iwf=1:nnwf) imat(iwf,iwf)=1d0+img*merge(nms_delta,0d0,nms) !identical matrix
  allocate(kmat(1:nnwf,1:nnwf,(1-npm)*nwhis:nwhis))
  BIGiqloop: do iq = iqxini,iqxend
    if(.NOT. MPI__task(iq)) cycle BIGiqloop
    q = qibze(:,iq)
    write(6,"('===== do : iq wibz(iq) q=',i6,f13.6,3f9.4,' ========')") iq,q !,wibz(iqlist(iq)),qshort !qq
    call writemem('hmagnon start gettetwt')
    GETtet: block
      integer:: isdummy
      real(8) :: ev_w1(nwf,nqbz), ev_w2(nwf,nqbz)
      complex(8):: evc_w1(nwf,nwf), evc_w2(nwf,nwf)
      readeigen: do kx=1,nqbz      !!! ev_w1, ev_w2 unit: [Ry]
        call wan_readeval2(  qbz(:,kx), is,  ev_w1(1:nwf,kx), evc_w1) !eigenvalue eigenfunciton
        call wan_readeval2(q+qbz(:,kx), isf, ev_w2(1:nwf,kx), evc_w2)
      enddo readeigen
      call gettetwt(q,iq,isdummy,isdummy,ev_w1,ev_w2,nwf,wan) !! tetrahedron weight. iq is dummy index
      !!     ihw(ibjb,kx): omega index, to specify the section of the histogram., ibjb=1,nbnb
      !!     nhw(ibjb,kx): the number of histogram sections
      !!     jhw(ibjb,kx): pointer to whw
      !!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
      !!     : histogram weights for given ib,jb,kx for histogram sections
      !!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
    endblock GETtet
    call writemem('hmagnon start Im kmat')
    GETzxq: block ! zxq and zxqi are the main output after Hilbert transformation, ! zxqi is not used in hmagnon (imagomega=.false.)
      real(8) :: ev_w1(nwf), ev_w2(nwf) !dummy
      complex(8) :: zxqi(1,1,1), evc_w1(nwf,nwf), evc_w2(nwf,nwf)
      integer, allocatable :: nttp(:),  itw(:,:), itpw(:,:)
      integer :: nttp_max, ittp, jpm, it, itp, ibib
      real(8), allocatable :: whwc(:,:)
      complex(8), allocatable :: zw(:,:), wzw(:,:)
      kmat(:,:,:) = 0d0
      kxloop: do kx=1, nqbz
        if(mod(kx-1, mpi__size_k) /= mpi__rank_k)  cycle kxloop
        call wan_readeval2(  qbz(:,kx), is,  ev_w1, evc_w1) !eigenvalue eigenfunciton
        call wan_readeval2(q+qbz(:,kx), isf, ev_w2, evc_w2)
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
         do ibib = 1, nbnb(kx,jpm)
           do iw = ihw(ibib,kx,jpm), ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
             nttp(iw) = nttp(iw) + 1
           enddo
         enddo
         nttp_max = maxval(nttp(1:nwhis))
         allocate(itw(nttp_max,nwhis), itpw(nttp_max,nwhis), whwc(nttp_max,nwhis))
         nttp(:) = 0
         do ibib = 1, nbnb(kx,jpm) !! n,n' pair band index loop
           it = n1b(ibib,kx,jpm)  !index for n  for q   ! n1b(ibib,k,jpm) = n :band index for k (occupied),
           itp = n2b(ibib,kx,jpm) !index for n' for q+k ! n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
           do iw = ihw(ibib,kx,jpm), ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
             nttp(iw) = nttp(iw) + 1
             ittp = nttp(iw)
             itw(ittp,iw) = it
             itpw(ittp,iw) = itp
             whwc(ittp,iw) = whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
           enddo
         enddo
         allocate(zw(nttp_max,nnwf), wzw(nttp_max,nnwf))
         do iw = 1, nwhis
           if (nttp(iw) < 1) cycle
           do ittp = 1, nttp(iw)
             it = itw(ittp,iw)
             itp = itpw(ittp,iw)
             do inwf =1, nnwf
               iwf = wan_pair_index(inwf,1)
               jwf = wan_pair_index(inwf,2)
               zw(ittp, inwf) = dconjg(evc_w2(jwf,itp))*evc_w1(iwf,it) !a_{Rk alpha}^{(k+q)n'}* a_{Rl beta}^{kn}
               wzw(ittp,inwf) = whwc(ittp,iw)*zw(ittp,inwf)
             enddo
           enddo
           istat = zmm(zw, wzw, kmat(1,1,iw*(3-2*jpm)), nnwf, nnwf, nttp(iw), opA=m_op_C, beta=(1d0,0d0), ldA=nttp_max, ldB=nttp_max)
         enddo
         deallocate(nttp, itw, itpw, whwc, zw, wzw)
        enddo jpmloop
      enddo kxloop
      call tetdeallocate()      ! --> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
      if (negative_cut) kmat(:,:,-nwhis:-1) = 0d0
      call writemem('hmagnon start dpsion')
      mpi_k_accumulate: block
        use m_mpi,only: MPI__reduceSum
        do jpm=1, npm
          do iw=1, nwhis
            call MPI__reduceSum(0, kmat(1,1,iw*(3-2*jpm)), nnwf*nnwf, communicator=comm_k)
          enddo
        enddo
      end block mpi_k_accumulate
      if(mpi__root_k) then
        call dpsion_init(realomega, imagomega, .false.)
        call dpsion_chiq(realomega, imagomega, .false., kmat, zxqi, nnwf, nnwf, schi, 1, 1d99) !Inplace routine: kmat is overwritten by zxq
      endif
    endblock GETzxq

    if(geteta .and. (.not. mpi__root_k)) exit BIGiqloop
    if(.not. mpi__root_k) cycle BIGiqloop

    !Below lines are executed only by root of mpi__rank_k
    call writemem('hmagnon start getting R')
    if(associated(zxq)) nullify(zxq)
    zxq(1:,1:,nw_i:) => kmat(1:nnwf,1:nnwf,nw_i:nw)
    where(abs(dimag(zxq))<1d-15) zxq = dreal(zxq) ! threshold for Im[K] (zxq)
    allocate(wkmat(1:nnwf,1:nnwf), rmat(1:nnwf,1:nnwf)) !WKmatrix, WKmatrix_inv

    IfGetEta: if(geteta) then
      block
        integer :: iunit
        complex(8) ::eval_wk(nnwf)
        istat = zmm(scrw, zxq(:,:,0), wkmat, nnwf, nnwf, nnwf)
        call diagcvuh3(wkmat(:,:),nnwf,eval_wk) !!   eval_wk is complex array because of Non-Hermite WK
        eta = -1d0/maxval(abs(eval_wk))
        write(stdo,ftox) "now eigenvalue abs(WK)",abs(eval_wk(1)),"is inversed"
        write(stdo,ftox) "check eigenvalue Re(WK)",dreal(eval_wk(1))
        write(stdo,ftox) "check eigenvalue Im(WK)",dimag(eval_wk(1))
        write(stdo,ftox) "wkmat calculated eta:", eta !negative value
        open(newunit=iunit,file='__EtaMagnon',status='replace',form='formatted',action='write')
        write(iunit,*) eta
        close(iunit)
        open(newunit=iunit, file='Kpmdiag_q0.dat', status='replace', action='write')
        write(iunit,ftox) "# iw omega(eV) Tr K/znorm TrdiagK/znorm"
        do iw = nw_i,nw
          www = merge(-freq_r(-iw),freq_r(iw),iw<0)
          write(iunit,"(e14.6,4e17.9)") www*hartree, tr_mat_onsite(zxq(:,:,iw))/hartree/znorm, &
                                      & hartree*tr_mat_onsite_diag(zxq(:,:,iw))/hartree/znorm
        enddo
        close(iunit)
      endblock
      exit Bigiqloop
    endif IfGetEta

    allocate(r_tr(nw_i:nw), r_diag(nw_i:nw), k_tr(nw_i:nw), k_diag(nw_i:nw))
    iwloop: do iw = nw_i,nw
      istat = zmm(scrw, zxq(:,:,iw), wkmat, nnwf, nnwf, nnwf, alpha=dcmplx(eta,0d0)) !wkmat = etaWK
      wkmat(1:nnwf,1:nnwf) = imat(1:nnwf,1:nnwf) - wkmat(1:nnwf,1:nnwf) ! wkamt = 1 - etaWK
      istat = zminv(wkmat, n=nnwf) ! wkmat = (1- eta WK)^-1
      istat = zmm(zxq(:,:,iw), wkmat, rmat, nnwf, nnwf, nnwf) !rmat = K (1-eta WK)^-1
      k_tr(iw) = tr_mat_onsite(zxq(:,:,iw))/znorm
      r_tr(iw) = tr_mat_onsite(rmat)
      k_diag(iw) = tr_mat_onsite_diag(zxq(:,:,iw))/znorm
      r_diag(iw) = tr_mat_onsite_diag(rmat)
    enddo iwloop

    istat = writem(file_tr_kpm, rec=iq, data=k_tr(nw_i:nw))
    istat = writem(file_tr_rpm, rec=iq, data=r_tr(nw_i:nw))
    istat = writem(file_tr_diag_kpm, rec=iq, data=k_diag(nw_i:nw))
    istat = writem(file_tr_diag_rpm, rec=iq, data=r_diag(nw_i:nw))
    deallocate(rmat, wkmat)
    deallocate(r_tr, r_diag, k_tr, k_diag)
    call writemem('hmagnon end iq'//trim(charext(iq)))
  enddo BIGiqloop

  if(geteta) call rx0( ' OK! hmagnon get eta')
  call mpi_barrier(comm, istat)

  write(stdo,"('eta for 1-eta*WK:',f13.8)") eta
  ReformatOutputFiles:if(mpi__root) then
    block
      integer :: epslgroup_old, file_tr_kpm_out, file_tr_rpm_out
      character(128):: filename_kmat, filename_rmat
      character(3) :: charnum3
      real(8) :: q_old(3), q_position, dq, omega
      logical :: opened
      allocate(r_tr(nw_i:nw), r_diag(nw_i:nw), k_tr(nw_i:nw), k_diag(nw_i:nw))
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
          filename_kmat = 'TrKpm.syml'//charnum3(epslgroup(iq))
          filename_rmat = 'TrRpm.syml'//charnum3(epslgroup(iq))
          open(newunit=file_tr_kpm_out, file=filename_kmat, status='replace', form='formatted', action='write')
          open(newunit=file_tr_rpm_out, file=filename_rmat, status='replace', form='formatted', action='write')
          write(file_tr_kpm_out, '(A)')' # qx qy qz q_pos omega(eV) Real_Tr_K/eV Imag_Tr_K/eV Real_Tr_Diag_K/eV Imag_Tr_Diag_K/eV'
          write(file_tr_rpm_out, '(A)')' # qx qy qz q_pos omega(eV) Real_Tr_R/eV Imag_Tr_R/eV Real_Tr_Diag_R/eV Imag_Tr_Diag_R/eV'
        endif
        istat = readm(file_tr_kpm, rec=iq, data=k_tr(:))
        istat = readm(file_tr_rpm, rec=iq, data=r_tr(:))
        istat = readm(file_tr_diag_kpm, rec=iq, data=k_diag(:))
        istat = readm(file_tr_diag_rpm, rec=iq, data=r_diag(:))
        dq = sqrt(sum((q(:)-q_old(:))**2))
        q_position = q_position + dq
        do iw = nw_i, nw
          www = merge(-freq_r(-iw),freq_r(iw),iw<0)
          omega = www*hartree
          write(file_tr_kpm_out, "(4f9.5,e14.6,4e17.9)") q(1:3), q_position, omega, k_tr(iw)/hartree, k_diag(iw)/hartree
          write(file_tr_rpm_out, "(4f9.5,e14.6,4e17.9)") q(1:3), q_position, omega, r_tr(iw)/hartree, r_diag(iw)/hartree
        enddo
        write(file_tr_kpm_out, *)
        write(file_tr_rpm_out, *)
        epslgroup_old = epslgroup(iq)
        q_old = q(:)
      enddo
      inquire(unit=file_tr_kpm_out, opened=opened)
      if(opened) close(file_tr_kpm_out)
      inquire(unit=file_tr_rpm_out, opened=opened)
      if(opened) close(file_tr_rpm_out)
      deallocate(r_tr, r_diag, k_tr, k_diag)
    endblock
  endif ReformatOutputFiles
  istat = closem(file_tr_kpm)
  istat = closem(file_tr_rpm)
  istat = closem(file_tr_diag_kpm)
  istat = closem(file_tr_diag_rpm)
  call rx0( ' OK! hmagnon mode')
END subroutine hmagnon
end module m_hmagnon
