!>  Calculate Chi^+-, spin susceptibility. 
module m_hmagnon 
  contains
subroutine hmagnon() bind(C)
  use m_readwan,only: write_qdata, wan_readeval2, readscr, read_wandata, nwf, tr_mat_onsite, tr_mat_onsite_diag, &
                    & set_wan_nnwf, nnwf, set_wan_scrw, scrw, wan_pair_index
  use m_ReadEfermi,only: readefermi
  use m_read_bzdata,only: read_bzdata, nqbz, nqibz, ginv, qbz, qibz, wibz, nstbz, wqt=>wt, q0i, nq0i ,nq0iadd, epslgroup, neps
  use m_genallcf_v3,only: genallcf_v3, nspin, niw_in=>niw, plat
  use m_keyvalue,only: getkeyvalue
  use m_freq,only: getfreq, frhis, freq_r, freq_i, nwhis, nw_i, nw, npm
  use m_tetwt,only: tetdeallocate, gettetwt, whw, ihw, nhw, jhw, ibjb, nbnbx, nhwtot, n1b, n2b, nbnb
  use m_readgwinput,only: ReadGWinputKeys
  use m_lgunit,only: m_lgunit_init, stdo
  use m_dpsion,only: dpsion_init, dpsion_chiq
  use m_mpi,only: MPI__Initialize_magnon, MPI__consoleout_magnon, MPI__AllreduceSum
  use m_mpi,only: MPI__rank=>mpi__rankMG,MPI__size=>mpi__sizeMG, MPI__root ,comm
  use m_blas, only: m_op_C, zmm => zmm_h
  use m_lapack, only: zminv => zminv_h
  use m_ftox
  implicit none
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  xxx removed--> eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  integer, parameter :: ndble=8
  integer:: iqbz, iqindx, iww, iqq
  integer:: iwf, jwf, inwf, kwf, lwf, ijwf, klwf, ijwf_j, imaximr
  integer:: ifchipmz_wan, ifchipmr_wan
  integer:: niw, ifif, ierr, MPI__MEq
  integer:: iqxini, iqxend, iqxendx, i, ini, ix, is ,iw, iq, nspinmx, kx, isf, ik, ibib, verbose, istat
  integer:: nqbze, nqibze, isdummy
  real(8):: q(3), omg2max, wemax, rydberg, hartree
  real(8), parameter:: schi=1d0, ua = 1d0
  real(8):: maximr, w_maximr, nms_delta, www
  real(8), allocatable:: qbze(:,:), qibze(:,:)
  real(8), allocatable:: rpa_maximr(:),mf_maximr(:) !! for MAX(Im[R])
  complex(8), pointer:: zxq(:,:,:) => null()
  complex(8), allocatable, target :: kmat(:,:,:)
  complex(8), allocatable:: evc_w1(:,:), evc_w2(:,:), wkmat(:,:), rmat(:,:)
  complex(8), allocatable:: scrw_original(:,:) !screening W
  complex(8), allocatable:: eval_wk(:), eval_wk2(:),trmat22(:)
  complex(8), parameter :: img=(0d0,1d0)
  complex(8)::trmat,trmatt,trmat1,trmat2
  complex(8),allocatable::imat(:,:) !unit matrix for 1-WK
  logical:: debug=.false.
  logical:: realomega=.true., imagomega=.true.,omitqbz=.false. !, noq0p
  logical:: chipm=.false., nolfco=.false., epsmode=.false., normalm=.false. ,autogamma=.false., addgamma=.false.
  logical:: wan=.true., lhm, lsvd, nms 
  logical, allocatable :: mpi__task(:)
  character(4):: charnum4
  real(8)::qlat(3,3), eta
  real(8), parameter :: pi = 4d0*datan(1d0)
  real(8), parameter :: znorm=-1d0*pi ! normalization of Im[K]:
  integer, parameter :: size_lim=999
  logical, parameter :: onsite_approx = .false.
!!! q on symline
  integer:: nqsym
  logical:: negative_cut, write_hmat, output_ddmat
  logical:: cma_mode !cma_mode for Cu2MnAl only 2019/09/27
  real(8):: cma_up_shift, cma_dn_shift, cma_wshift
  integer(4):: cma_iwf_s, cma_iwf_e
  complex(8):: sumrpa_maximr(1), summf_maximr(1)
  hartree  = 2d0*rydberg()
  call m_lgunit_init()
  call MPI__Initialize_magnon()
  call MPI__consoleout_magnon('hmagnon',size_lim) ! size_lim for saving memory (avoid swapping)
  call cputid(0)
  imagomega =.false.
  omitqbz =.true. ! for QforEPS (20Feb, 2020)
  epsmode = .true.
  chipm   =.true.
  nolfco  =.true.
  wan     =.true.
  call genallcf_v3(incwfx=0) !!incwfin=0 =>ForX0 for core in GWIN. in module m_genallcf_v3 Readin by genallcf. Set basic data for crystal
  write(6,"(' nqbz nqibz =',2i5)") nqbz,nqibz
  if(chipm .AND. nspin==1) call rx( 'chipm mode is for nspin=2')  ! We fix newaniso2=T now.
  !! Prof.Naraga said " write(6,*)'Timereversal=',Timereversal()" here caused a stop in ifort ver.1x.x. Why? May be a compilar bug, and fixed now.
  !! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
  !! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
  call read_BZDATA() !  !! Read electron gas mode or not.
  call ReadGWinputKeys() ! jun2020 new routint to read all inputs
  call getkeyvalue("GWinput","lHermite",lhm,default=.false.)
  call getkeyvalue("GWinput","lsvd",lsvd,default=.false.)
  call getkeyvalue("GWinput","nms",nms,default=.false.)  !!! For NiMnSb
  call getkeyvalue("GWinput","nms_delta",nms_delta,default=1d-6)
  call getkeyvalue("GWinput","negative_cut",negative_cut,default=.false.)
  call getkeyvalue("GWinput","write_hmat",write_hmat,default=.false.)
  call getkeyvalue("GWinput","output_ddmat",output_ddmat,default=.false.) 
  call getkeyvalue("GWinput","cma",cma_mode,default=.false.)   ! 2019/09/27 cma_mode for Cu2MnAl
  call getkeyvalue("GWinput","cma_dn_shift",cma_dn_shift,default=0d0)
  call getkeyvalue("GWinput","cma_up_shift",cma_up_shift,default=0d0)
  call getkeyvalue("GWinput","cma_wshift",cma_wshift,default=1d-6)
  call getkeyvalue("GWinput","cma_iwf_start",cma_iwf_s,default=999)
  call getkeyvalue("GWinput","cma_iwf_end"  ,cma_iwf_e,default=999)
  write(6,*) "lsvd, lhm, nms, nms_delta",lsvd,lhm,nms,nms_delta,cma_mode
  write(6,*) "negative_cut",negative_cut
!  write(6,*) "reduce mpi_size for saving memory (default =999) ",size_lim
  if(MPI__root) then
    do i=1,nqbz
      if(i<10 .OR. i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
      if(i==10 .AND. nqbz>18) write(6,"('... ')")
    enddo
    write(6,*)' !!nqbz nqibz =',nqbz,nqibz
  endif
  call readefermi() !!! ef:     Fermi energy at 0 K
  write( 6,*) ' num of zero weight q0p=',neps
  write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  !! Readin q+G. nqbze and nqibze are for adding Q0P related points to nqbz and nqibz.
  nqbze  = nqbz *(1 + nq0i)
  nqibze = nqibz + nq0i
  allocate( qbze(3, nqbze), qibze(3, nqibze))
  qbze(:,1:nqbz)   = qbz(:,1:nqbz)
  qibze(:,1:nqibz) = qibz(:,1:nqibz)
  do i = 1,nq0i
    qibze(:,nqibz+i)  = q0i(:,i)
    ini = nqbz*(1 + i -1)
    do ix=1,nqbz
      qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
      if( abs(qbze(1,ini+ix)+0.1d0)+abs(qbze(2,ini+ix)+0.1d0)<1d-6 ) then
        write(6,"('hx0fp0 qbze q0i=',i8,3f18.14,2x,3f14.10)") ini+ix,qbze(:,ini+ix),q0i(:,i)
      endif
    enddo
  enddo
  call read_wandata()    ! nwf, nsp_w,nqtt_w ! --- okumura Read dimensions of hamiltonian_wannier, spin, nqtt
  call set_wan_nnwf(onsite_approx) !set nnwf ~ # of RiRj (onsite_approx = .true.), RiR'j (onsite_approx = .flase. ), wan_pair_index
  call set_wan_scrw(onsite_approx) !set scrw
  if(mpi__root) write(stdo,ftox) '# nwf, nnwf:', nwf, nnwf
  !Weight for irreducible q-point (qibz); do iq=1,nqibz; write(6,"('wibz',4f9.4)") wibz(iq),qibz(:,iq); enddo
  iqxend = nqibz !+ nq0i
  !! Shift W by hand (Oct.02, 2019)
  if (cma_mode) allocate(scrw_original(nnwf,nnwf), source = scrw)
  call minv33tp(plat,qlat)
  if(verbose()>50) print *,'eeee exit of init_readeigen2'
  do iq=1,nqibz
    iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)      !     nstibz(iq) = nstbz(iqbz)
    write(6,"(' iq qibz nstibz=',2i5,3f9.4,i5)")iq,iqbz,qibz(:,iq) !,nstibz(iq)
  enddo
  ! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
  wemax   = 5d0 !max value for plot
  omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
  !! NOTE: npmtwo=T sets npm=2   !! optional npmtwo is added aug2017   !! 20190604 Im[K]
  niw=niw_in
  if( .NOT. imagomega) niw=1  !dummy
  call Getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw,ua, npmtwo=.true.)!,tetra
  if(realomega .AND. mpi__root) then
    open(newunit=ifif,file='freq_r') !write number of frequency points nwp and frequensies in 'freq_r' file
    write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
    do iw= nw_i,-1
      write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw
    enddo
    do iw= 0,nw
      write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw
    enddo
    close(ifif)
  endif
  if(MPI__root) write(6,"(' nw_i nw niw npm=',4i5)") nw_i,nw,niw,npm
  nspinmx = nspin
  iqxini = merge(nqibz + 1,1,omitqbz)
  iqxend = nqibz + nq0i
  write(6,"('iqxini,iqxend ',2I8)") iqxini,iqxend
  do iq = iqxini,iqxend
    write(6,"('iq, qibze:',I8,3f9.4)") iq-iqxini,qibze(:,iq)
    autogamma=.true.
  enddo
  iqxendx=iqxend
  call write_qdata(ginv,nqbz,qbz(:,:))
  !! Gamma point is automatically added if qibze does not include.
  if ( .NOT. sum(abs(qibze(:,iqxini)**2)) == 0d0) then
    write(6,*) "Gamma point is automatically added"
    iqxini=iqxini-1 !! Gamma: iqxini-1
    addgamma=.true.
  endif
  nqsym=iqxend-iqxini+1
  write(6,"('nqsym:',I4)") nqsym
  rankdivider: block 
    use m_mpi,only: mpi__sizeMG,mpi__rankMG
    integer :: mpi__ranktab(1:nqsym)
    if(mpi__rankMG==0) write(6,*) "MPI_hmagnon_rankdivider:"
    allocate( mpi__task(1:nqsym) )
    mpi__task(:) = .false.
    mpi__ranktab(1:nqsym)=999999
    mpi__MEq=1+nqsym/mpi__sizeMG
    if( mpi__sizeMG == 1 ) then
      mpi__task(:) = .true.
      mpi__ranktab(:) = mpi__rankMG
    else   
      do iq=1,nqsym
        mpi__ranktab(iq) = mod(iq-1,mpi__sizeMG)  !rank_table for given iq. iq=1 must give rank=0
        if(mpi__ranktab(iq) == mpi__rankMG) mpi__task(iq) = .true.         !mpi__task is nodeID-dependent.
        if(mpi__rankMG==0) write(6,"('  iq irank=',2i5)")iq,mpi__ranktab(iq)
      enddo
    endif
    write(6,*) "mpi__sizeMG nqsym mpi_MEq(:): ",mpi__sizeMG,nqsym,mpi__MEq
  endblock rankdivider
  allocate(rpa_maximr(mpi__MEq),mf_maximr(mpi__MEq)) !list of w(MAX(Im[R])): magnon peak
  rpa_maximr=0d0
  mf_maximr=0d0
  do iq=1,nqsym
    if (MPI__task(iq)) write(6,'("iq,MPI_rank",3I8)') iq,MPI__rank!,mpi__ranktab(iq)
  enddo
  call MPI_barrier(comm,ierr)
  allocate(imat(1:nnwf,1:nnwf),source=(0d0,0d0))
  forall(iwf=1:nnwf) imat(iwf,iwf)=1d0+img*merge(nms_delta,0d0,nms) !identical matrix
  allocate(evc_w1(nwf,nwf),evc_w2(nwf,nwf)) !, zxq_d(nnwf,nnwf) )
  allocate(eval_wk(nnwf),eval_wk2(nnwf),trmat22(nw_i:nw))
  allocate(kmat(1:nnwf,1:nnwf,(1-npm)*nwhis:nwhis))
  imaximr = 0
  BIGiqqloop: do 1001 iqq = iqxini,iqxend      ! NOTE: q=(0,0,0) is active iqq=iqxini (see autogamma)
!   if(MPI__rank > size_lim) cycle !reduce mpi-size for test (skip 21-32)
    iq = iqq-iqxini+1 !! start with iq=1 for convenience
    if( .NOT. MPI__task(iq) .AND. iq /= 1) cycle
    q = merge([0d0,0d0,0d0], qibze(:,iqq),(autogamma .AND. iq==1)) !! automatically set q=(0 0 0)
    imaximr = imaximr+1      !!! imaximr=1,2,...,mpi__MEq
    if(iq==1.and.sum(q**2)>1d-10) call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
    write(6,"('===== do 1001: iq wibz(iq) q=',i6,f13.6,3f9.4,' ========')") iq,q !,wibz(iqlist(iq)),qshort !qq
    if(lhm) cycle
    GETtet: block
      integer,parameter:: is=1,isf=2
      real(8)::ev_w1(nwf,nqbz),ev_w2(nwf,nqbz)
      readeigen: do 5001 kx=1,nqbz      !!! ev_w1, ev_w2 unit: [Ry]
        call wan_readeval2(  qbz(:,kx), is,  ev_w1(1:nwf,kx), evc_w1) !eigenvalue eigenfunciton
        call wan_readeval2(q+qbz(:,kx), isf, ev_w2(1:nwf,kx), evc_w2)
        onlyCu2MnAl: if (cma_mode) then !! only Cu2MnAl (cma)         !! Energy of Mn3d(dn) is moved by cma_shitf
          !$$$     if (iq==1) write(6,"('cma_mode: iwf_s, iwf_e',2i4)") cma_iwf_s,cma_iwf_e
          !$$$     if (iq==1) write(6,"('cma_mode: cma_up_shift, cma_dn_shift',2E13.5,' [eV]')") cma_up_shift,cma_dn_shift
          do iwf = 1, nwf
            if ( cma_iwf_s <= iwf .AND. iwf <= cma_iwf_e ) then
              ev_w1(iwf,kx) = ev_w1(iwf,kx) + cma_up_shift/rydberg()
              ev_w2(iwf,kx) = ev_w2(iwf,kx) + cma_dn_shift/rydberg()
            endif
          enddo
        endif onlyCu2MnAl
5001  enddo readeigen
      write(6,"(' = start wan_gettetwt =',2i6,3f9.4)") nwf,iq,q
      call gettetwt(q,iq,isdummy,isdummy,ev_w1,ev_w2,nwf,wan) !! tetrahedron weight. 
      !!     ihw(ibjb,kx): omega index, to specify the section of the histogram., ibjb=1,nbnb
      !!     nhw(ibjb,kx): the number of histogram sections
      !!     jhw(ibjb,kx): pointer to whw
      !!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
      !!     : histogram weights for given ib,jb,kx for histogram sections
      !!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
    endblock GETtet
    GETzxq: block ! zxq and zxqi are the main output after Hilbert transformation, ! zxqi is not used in hmagnon (imagomega=.false.)
      integer,parameter:: is=1,isf=2
      complex(8)::zxqi(1,1,1) !,wanmat(1:nnwf,1:nnwf)
      real(8) ::ev_w1(nwf), ev_w2(nwf) !dummy
      integer, allocatable :: nttp(:),  itw(:,:), itpw(:,:)
      integer :: nttp_max, ittp, jpm, it, itp
      real(8), allocatable :: whwc(:,:)
      complex(8), allocatable :: zw(:,:), wzw(:,:)
      ! zxq=0d0
      kmat=0d0
      kxloop:       do 2011 kx=1,nqbz 
        call wan_readeval2(  qbz(:,kx), is,  ev_w1, evc_w1) !eigenvalue eigenfunciton
        call wan_readeval2(q+qbz(:,kx), isf, ev_w2, evc_w2)
        jpmloop:    do 2012 jpm=1,npm ! jpm=2: negative frequency
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
             it = itw(ittp,iw); itp = itpw(ittp,iw)
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
2012    enddo jpmloop
2011  enddo kxloop
      call tetdeallocate()      ! --> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
      if (negative_cut) kmat(:,:,-nwhis:-1)=0d0
      call dpsion_init(realomega, imagomega, .false.)
      call dpsion_chiq(realomega, imagomega, .false., kmat, zxqi, nnwf, nnwf, schi, 1, 1d99) !! Inplace routine: kmat is overwritten by zxq
      ! call dpsion5( realomega, imagomega, kmat, nnwf,nnwf, zxq, zxqi,.false., schi,1,1d99,1d99)  !! kmat ---> zxq
    endblock GETzxq
    if(associated(zxq)) nullify(zxq)
    zxq(1:nnwf,1:nnwf,nw_i:nw) => kmat(1:nnwf,1:nnwf,nw_i:nw)
    if(lhm) then !Enforce zxq Hermitian 
      ZxqHermitian: block
      complex(8), allocatable:: zxq2(:,:,:)
      allocate(zxq2(nnwf,nnwf,nw_i:nw),source=zxq(1:nnwf,1:nnwf,nw_i:nw))
      zxq=0d0
      ijwf=0
      do 3006 iwf=1,nwf
        do 3007 jwf=1,nwf
          ijwf=ijwf+1
          klwf=0
          do 3008 kwf=1,nwf
            do 3009 lwf=1,nwf
              klwf=klwf+1
              if (ijwf == klwf) then
                zxq(ijwf,ijwf,:)= zxq2(ijwf,ijwf,:)
              elseif(ijwf>klwf) then             !!! ijwf > klwf
                zxq(ijwf,klwf,:)=( zxq2(ijwf,klwf,:) + dconjg(zxq2(klwf,ijwf,:)) )/2.0
                zxq(klwf,ijwf,:)=dconjg( zxq(ijwf,klwf,:) )
              endif
3009        enddo
3008      enddo
3007    enddo
3006  enddo
      ijwf=0; klwf=0
      deallocate(zxq2)
      endblock ZxqHermitian
    endif
    where(abs(dimag(zxq))<1d-15) zxq=dreal(zxq) ! threshold for Im[K] (zxq)
    allocate(wkmat(1:nnwf,1:nnwf), rmat(1:nnwf,1:nnwf)) !WKmatrix, WKmatrix_inv
    GetEta: if (iq==1) then ! (1-eta*WK)
      ! wkmat(1:nnwf,1:nnwf) =matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,0)) !omega=0
      istat = zmm(scrw, zxq(:,:,0), wkmat, nnwf, nnwf, nnwf)
      call diagcvuh3(wkmat(:,:),nnwf,eval_wk) !!   eval_wk is complex array because of Non-Hermite WK
      eta=-1d0/maxval(abs(eval_wk))
      write(6,*) "now eigenvalue abs(WK)",abs(eval_wk(1)),"is inversed"
      write(6,*) "check eigenvalue Re(WK)",real(eval_wk(1))
      write(6,*) "check eigenvalue Im(WK)",aimag(eval_wk(1))
      write(6,*) "wkmat calculated eta:", eta !negative value
    endif GetEta                     !iq==1
    InitalWritewan_ChiPM: if (MPI__task(iq)) then
      open(newunit=ifchipmz_wan,file="wan_ChiPMz.mat"//charnum4(iq))
      open(newunit=ifchipmr_wan,file="wan_ChiPMr.mat"//charnum4(iq))
      print *,'ifchipm=',ifchipmz_wan,ifchipmr_wan
      if (iq==1) then
        write(ifchipmz_wan,*) "# syml: Gamma"
        write(ifchipmr_wan,*) "# syml: Gamma"
      else
        if (addgamma) then
          write(ifchipmz_wan,*) "# syml: ",epslgroup(iq-1)," "
          write(ifchipmr_wan,*) "# syml: ",epslgroup(iq-1)," "
        else
          write(ifchipmz_wan,*) "# syml: ",epslgroup(iq)," "
          write(ifchipmr_wan,*) "# syml: ",epslgroup(iq)," "
        endif
      endif
    endif InitalWritewan_ChiPM
    maximr=0d0; w_maximr=0d0 !!! search for Im[R] peak
    iwloop: do 2050 iw = nw_i,nw
      ! call diagcvuh3(zxq(:,:,iw),nnwf,eval_wk)
!!! Hermite matrix diagonization        ! all diagcvh2(zxq(:,:,iw),nnwf,eval_wk)
      ! where(abs(dimag(eval_wk)) < 1d-16) eval_wk=dreal(eval_wk)
      www = merge(-freq_r(-iw),freq_r(iw),iw<0)
      if(MPI__task(iq)) then ! Check for KK relatioin ! Some Error: remove this section or modified ... (Okumura, Oct02,2019)
        if(iw==nw_i) then !only first line
          if(debug) then !make it debug mode
            write(6,"(' --check iww sum(zxq)',i5,2E13.5)") iw,sum((zxq(:,:,:)))
            trmat22=0d0
            do iww =nw_i,nw
              if (iww==0) cycle !skip w=0 (Cauthy principle integral)
              call diagcvuh3(zxq(:,:,iww),nnwf,eval_wk2)
              where(abs(dimag(eval_wk2)) < 1d-16) eval_wk2=dreal(eval_wk)
              trmat22(iww) =sum(eval_wk2)/(znorm)* merge(-freq_r(iww),freq_r(iww),iww<0) *abs(freq_r(abs(iww))-freq_r(abs(iww)-1))
            enddo
            trmat2=sum(trmat22)
            call diagcvuh3(zxq(:,:,0),nnwf,eval_wk2) !for Re[K(w=0)]
            write(ifchipmz_wan,"('# int ImK/omega dw, Re[K(0)]=',2E13.5)") dimag(trmat2),sum(dreal(eval_wk2))
          else
            write(ifchipmz_wan,'(A)') "# Skip calculation of int ImK/omega dw"
          endif
        endif
        ! write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,2E17.9)") q,iw,www*2d0,hartree*sum(eval_wk)/(znorm)
        write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,4E17.9)") q,iw,www*2d0,hartree*tr_mat_onsite(zxq(:,:,iw))/(znorm), &
                                                               & hartree*tr_mat_onsite_diag(zxq(:,:,iw))/(znorm)
        !! Im[K] [1/Ry] ? write d-d (diagonal) and d-other (non-diagonal)
        ! if (output_ddmat) then
        !   if (iw==1000) then !!! Note: writeddmat(matrix, nwf, nw_i, nw, filename, diagonal or non-diagnal)
        !     call writeddmat(zxq(:,:,iw),nwf,"wan_ChiPMr.mat.dd",.true.,zxq_d(:,:)) 
        !     call writehmat(zxq_d(:,:),nwf,"zxqdmat_check.dat")
        !     write(6,*) "PASS for zxqmat_check"
        !     deallocate(zxq_d)
        !   endif
        ! endif
      endif
      if(cma_mode) then
        scrw = scrw_original
        if (sum(q**2)> 0.5**2) then
          do jwf=1,nwf
            if ( .NOT. cma_iwf_s <= jwf .AND. jwf <= cma_iwf_e) cycle
            ijwf_j=(jwf-1)*nwf+iwf
            scrw(ijwf_j,ijwf_j)=scrw_original(ijwf_j,ijwf_j) + cmplx(dble(cma_wshift/hartree),0d0,kind(0d0))
          enddo
        endif
      endif
      !!  K/(1-WK) = K(1-WK)^(-1)  !W shift (q is far from Gamma) ! W shift for Cu2MnAl
      ! wkmat(1:nnwf,1:nnwf) =eta*matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,iw)) !! WK matrix
      ! wkmat2(1:nnwf,1:nnwf)=imat(1:nnwf,1:nnwf)-wkmat(1:nnwf,1:nnwf)        ! c Hermite check for 1-eWK for developing code
      ! call matcinv(nnwf,wkmat2(1:nnwf,1:nnwf)) ! inv(1-WK)
      ! rmat(1:nnwf,1:nnwf,iw)=matmul(zxq(1:nnwf,1:nnwf,iw), wkmat2(1:nnwf,1:nnwf))
      ! call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
      ! call diagwan(rmat(:,:,iw),eval_wk) !! diagonalization for R
      ! trmat = sum(eval_wk(1:nnwf))
      ! call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
      !MO replace above lines with zmm calls
      istat = zmm(scrw, zxq(:,:,iw), wkmat, nnwf, nnwf, nnwf, alpha=dcmplx(eta,0d0)) !wkmat = etaWK
      wkmat(1:nnwf,1:nnwf) = imat(1:nnwf,1:nnwf) - wkmat(1:nnwf,1:nnwf) ! wkamt = 1 - etaWK
      istat = zminv(wkmat, n=nnwf) ! wkmat = (1- eta WK)^-1
      istat = zmm(zxq(:,:,iw), wkmat, rmat, nnwf, nnwf, nnwf) !rmat = K (1-eta WK)^-1
      trmat = tr_mat_onsite(rmat)
      if(MPI__task(iq)) then
        write(ifchipmr_wan,"(3f9.4,i6,E12.4,x,4E12.4)")q,iw,www*2d0,hartree*trmat, hartree*tr_mat_onsite_diag(rmat)
        if (0d0 < www.and. maximr < -1d0*aimag(trmat)) then !search for MAX(Im[R]) 20180706 (0 - 1500 meV) if (0d0 <www*hartree .and. www*hartree < 1.5)
          maximr=-1d0*aimag(trmat)
          w_maximr=www
        endif
      endif
2050 enddo iwloop
    deallocate(rmat, wkmat)
    if(iq/=1) then         ! c(MF)        ! BZweight*omega[eV] for sum(E(q))/N
      mf_maximr(imaximr) = wibz(iq)*w_maximr*hartree         ! c(RPA)         ! BZweight*omega[eV] for sum(1/E(q))/N
      rpa_maximr(imaximr)= wibz(iq)/(w_maximr*hartree)
    endif
    if(MPI__task(iq)) then  !magnon peak
      write(6,"(' AAAA',I4,3f9.4)") iq,q
      write(6,"(' AAAA, w(MAX(im[R])), Im[R]',f13.5,E12.4)") w_maximr*2d0*rydberg()*1000,maximr/hartree
      write(ifchipmz_wan,*)
      write(ifchipmr_wan,*)    
      close(ifchipmz_wan)
      close(ifchipmr_wan)
    endif
1001 enddo BIGiqqloop
  write(6,*) "maximr RPA",rpa_maximr
  write(6,*) "maximr MFA",mf_maximr
  sumrpa_maximr=cmplx(sum(rpa_maximr(:)),0d0,kind(0d0))
  summf_maximr =cmplx(sum( mf_maximr(:)),0d0,kind(0d0))
  call MPI_barrier(comm,ierr)
  write(6,*) "MPIcheck",MPI__size,MPI__rank,sumrpa_maximr
  call MPI__AllreduceSum(sumrpa_maximr(1),1)
  call MPI__AllreduceSum( summf_maximr(1),1)
  write(6,"('sum(E(q)) for MFA:',f9.4)") real(summf_maximr(1))
  write(6,"('[sum(1/E(q))]inv for RPA',f9.5)") 1d0/real(sumrpa_maximr(1))
  call cputid(0)   !      call MPI__Finalize
  write(6,"('eta for 1-eta*WK:',f13.8)") eta
  call rx0( ' OK! hmagnon mode')
END subroutine hmagnon
end module m_hmagnon
