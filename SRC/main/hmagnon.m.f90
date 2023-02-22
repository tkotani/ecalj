!! hx0fp0.m.F in order to develop by Okumura (2017/03/23)
!!  Calculate x0, \epsilon, spin susceptibility.
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_magnon_chipm mode works).
program hmagnon
  use m_readwan,only: write_qdata, wan_readeigen, wan_readeval, wan_readeval2, &
       readscr, checkorb, checkorb2, diagwan, diagwan_tr, wan_imat, &
       writehmat, writeddmat, &
       read_wandata, nwf, nsp_w, nqtt_w
  use m_ReadEfermi,only: readefermi,ef_read=>ef, readefermi_kbt,ef_kbt
  use m_readeigen,only: readeval,init_readeigen,init_readeigen2
  use m_hamindex,only:qtt,nqtt
  use m_read_bzdata,only: read_bzdata, &
       nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,ginv, &
       qbz,wbz,qibz,wibz,qbzw, &
       idtetf,ib1bz,idteti, &
       nstar,irk,nstbz &
       ,wqt=>wt,q0i,nq0i ,nq0iadd,ixyz,epslgroup,nq0ix,neps
  use m_genallcf_v3,only: genallcf_v3, &
       nclass,natom,nspin,nl,nn, &
       nlmto,nlnmx, nctot,niw_in=>niw, &
       alat, delta,deltaw,esmr,clabl,iclass, &
       il,in,im,nlnm, &
       plat, pos,ecore
  use m_keyvalue,only: getkeyvalue
  use m_mpi,only: MPI__task,MPI__Initialize,MPI__Finalize,MPI__root, &
       MPI__hx0fp0_rankdivider2, &
       MPI__hmagnon_rankdivider, MPI__Initialize_magnon, MPI__consoleout_magnon, &
       MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank=>mpi__rankMG,MPI__size=>mpi__sizeMG, &
       MPI__ranktab,MPI__consoleout,MPI__barrier,MPI__MEq,MPI__REAL8send,MPI__REAL8recv,MPI__AllreduceSum
  !! frequency
  use m_freq,only: getfreq, frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw 
  use m_tetwt,only: tetdeallocate,gettetwt, whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  !! w0 and w0i (head part at Gamma point)
  use m_w0w0i,only: w0w0i, w0,w0i
  use m_readgwinput,only: ReadGWinputKeys
  use m_lgunit,only:m_lgunit_init

  implicit none
  !! ------------------------------------------------
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  !!-------------------------------------------------

  ! cccc this may be wrong or correct cccccccccc
  !r Be careful for the indexing...
  !r      A routine idxlnmc(nindxv,nindxc,...  in index.f
  !r      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
  !r      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
  !r      The indexing starts with core first and then valence on top of core
  !r      So n-index in "in" for valence electron is different from "inv".
  ! cccccccccccccccccccccccccccccccccccccccccccccccc
  real(8):: q(3),  qgbin(3),qx(3)
  real(8):: ua=1d0 ! this is a dummy.
  integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
  integer:: ndble=8
  integer:: nwordr
  real(8),allocatable:: vxcfp(:,:), &
       wgt0(:,:) !,q0i(:,:) wqt(:),
  complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:),zxq2(:,:,:), &
       zxq_d(:,:)

  !  tetrahedron method
  logical :: tetra=.true., usetetrakbt
  ! so as to reduce the memory usage.
  complex(8) :: fff,img=(0d0,1d0)
  logical :: debug
  logical :: realomega=.true., imagomega=.true.

  real(8) :: omg2max, wemax

  integer::maxocc2,iclose, &
       ixc,iqxini,iqxend,iqxendx, &
       i,mxx,ini,ix,is &
       ,iw,noccxv,noccx,iq,ngb &
       ,nprecx,nblochpmx,mrecl,ifwd,nspinmx,ibas &
       ,kx,isf,job &
       ,ihis,ik,ibib,ib1,ib2,j

  integer:: incwfin,  verbose
  logical:: omitqbz=.false. !, noq0p

  real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree &
       ,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
  integer:: iqbz,iqindx,iflegas,nmx &
       ,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1
  real(8) :: eclda_bh,eclda_pz,wk4ec,faca

  integer:: ifv,lxxilmx,ilm_r,nx_r,lb,nb,mb
  integer,allocatable:: nxx_r(:)
  real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
  character*3:: charnum3
  character*4:: charnum4
  character*5:: charnum5

  real(8)::schi=1d0
  integer:: new,nmxx,ii,iy,ipl1,ixx

  real(8),allocatable::SS(:),rwork(:),ss0(:)
  complex(8),allocatable:: UU(:,:),VT(:,:),work(:),ddd(:,:) &
       ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
  integer::lwork,info,imin,ifzxq
  complex(8),allocatable:: UU0(:,:),VT0(:,:)

  logical ::  chipm=.false.,nolfco=.false.,epsmode=.false.,normalm=.false., crpa=.false. &
       ,autogamma=.false., l1wkout=.false., addgamma=.false.
  real(8):: qs,qt,ww,muu, ddq(3)

  ! Feb2006 time-reversal=off case
  logical :: timereversal, testtimer,onceww
  integer:: jpm,ncc

  integer:: istat,imb,imb1,imb2,nmbas_in
  integer,allocatable:: imbas(:), imbas_s(:),iibas(:)
  real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq

!  integer,allocatable:: nwgt(:,:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
  integer:: ificlass,k !,ifile_handle
  real(8)::ebmx
  integer:: nbmx,mtet(3)

  ! Feb2020 QforEPS mode
  integer:: nqbze, nqibze!, ngcmx, ngpmx ifq0p,
  real(8), allocatable:: qbze(:,:), qibze(:,:)
  integer:: idammy
  !      integer, allocatable:: epslgroup(:)

  !     okumura
  integer::iwf,jwf,inwf,kwf,lwf,ijwf,klwf,nnwf,ijwf_j
  integer::ifgas,ifchipmz_wan,ifchipmr_wan,ifchipmrk_wan
  real(8),allocatable::ev_w1(:,:),ev_w2(:,:)
  complex(8),allocatable::evc_w1(:,:,:),evc_w2(:,:,:)
  logical::wan=.true.,lhm,lsvd,nms!,eibzmode
  logical(8):: lqsh_all
  real(8):: nms_delta
  logical(8):: ijklmag !for checkorb2
  complex(8),allocatable::kmat(:,:,:,:),wanmat(:,:) &
       ,wkmat(:,:),wkmat2(:,:),rmat(:,:,:),rmat2(:,:) &
       !     &     ,swkwmat(:,:),swkwmat2(:,:)
       !     &     ,sqemat(:,:),plmat(:,:),plpmat(:,:),plpmat_inv(:,:)
       ,wmat_check(:,:)
  complex(8)::trmat,trmatt,trmat1,trmat2
  complex(8),allocatable::imat(:,:) !unit matrix for 1-WK
  real(8):: wan_ecore(1), imagweight, www, www2
  logical::npmtwo,diag=.false.,t2g
  complex(8)::wan_i,wan_j,wan_k,wan_l,wanijkl
  integer::igv
  real(8)::qlat(3,3),qsh1(3),qsh2(3),znorm,rnqbz, eta
  ! or screening W
  complex(8),allocatable:: scrw(:,:) ,cmat2(:,:), scrw_original(:,:)
  complex(8),allocatable:: eval_wk(:),eval_k(:),eval_wk2(:)
  integer:: it,itp,isdummy,lorb, size_lim
!!! q on symline
  real(8)::qrot(3),cr=2d-5
  logical(8)::init, init2=.true.,gskip
  real(8)::rlatp(3,3),xmx2(3),qqin(3),qshort(3),qshort2
  integer:: nlatout(3,48),nout,iout,nqsym
  integer,allocatable:: iqlist(:)
  real(8),allocatable:: qshort_all(:,:)
  real(8),allocatable:: qshort_dammy(:,:)
!!! WK eigenvalue check
  integer::ifwkeigen,ifwkeigen2,ifwkeigen3
!!! temp variable
  logical(8)::lsmo,threshold=.true.
  logical:: negative_cut, write_hmat, output_ddmat
  logical:: cma_mode !cma_mode for Cu2MnAl only 2019/09/27
  real(8):: cma_up_shift, cma_dn_shift, cma_wshift
  integer(4):: cma_iwf_s,cma_iwf_e
  integer::ifwanmat, npm2
!!! for MAX(Im[R])
  integer::imaximr=0,niw,ifif
  real(8)::maximr,w_maximr
  complex(8)::sumrpa_maximr(1),summf_maximr(1)
  real(8),allocatable:: rpa_maximr(:),mf_maximr(:)

  real(8):: ef
  !! -------------------------------------------------------------------
  call getkeyvalue("GWinput","mpi_size_lim",size_lim,default=999)

  call MPI__Initialize_magnon()
  call M_lgunit_init()
  call MPI__consoleout_magnon('hmagnon',size_lim) ! size_lim for saving memory (avoid swapping)

  call cputid (0)
  hartree  = 2d0*rydberg()
  pi       = 4d0*datan(1d0)
  fourpi   = 4d0*pi
  sqfourpi = sqrt(fourpi)

  call cputid(0)

  write(6,*) " OK ixc=223    chipm_wannier sergey's "
  imagomega =.false.
  omitqbz=.true. ! for QforEPS (20Feb, 2020)
  epsmode = .true.
  chipm=.true.
  nolfco=.true.
  wan=.true.


  !! === Readin by genallcf. Set basic data for crystal
  !! See "use m_genallcf_v3" at the begining of this routine
  !!
  incwfin=0  !use ForX0 for core in GWIN
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  write(6,"(' nqbz nqibz =',2i5)") nqbz,nqibz
  if(chipm .AND. nspin==1) call rx( 'chipm mode is for nspin=2')

  !     ! We fix newaniso2=T now.

  !! Prof.Naraga said " write(6,*)'Timereversal=',Timereversal()"
  !! here caused a stop in ifort ver.1x.x. Why? May be a compilar bug, and fixed now.

  !! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
  !! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
  call read_BZDATA()

  !! Read electron gas mode or not.
  call ReadGWinputKeys() ! jun2020 new routint to read all inputs

  call getkeyvalue("GWinput","lHermite",lhm,default=.false.)
  call getkeyvalue("GWinput","lsvd",lsvd,default=.false.)
  call getkeyvalue("GWinput","nms",nms,default=.false.)  !!! For NiMnSb
  call getkeyvalue("GWinput","nms_delta",nms_delta,default=1d-6)
  call getkeyvalue("GWinput","negative_cut",negative_cut,default=.false.)
  call getkeyvalue("GWinput","write_hmat",write_hmat,default=.false.)
  call getkeyvalue("GWinput","output_ddmat",output_ddmat,default=.false.)
  ! 2019/09/27 cma_mode for Cu2MnAl
  call getkeyvalue("GWinput","cma",cma_mode,default=.false.)
  call getkeyvalue("GWinput","cma_dn_shift",cma_dn_shift,default=0d0)
  call getkeyvalue("GWinput","cma_up_shift",cma_up_shift,default=0d0)
  call getkeyvalue("GWinput","cma_wshift",cma_wshift,default=1d-6)
  call getkeyvalue("GWinput","cma_iwf_start",cma_iwf_s,default=999)
  call getkeyvalue("GWinput","cma_iwf_end"  ,cma_iwf_e,default=999)
  !! How to use "call getkeyvalue("IN","n1n2n3",ivkey,nsize,status=ret ) "
  write(6,*) "lsvd, lhm, nms, nms_delta",lsvd,lhm,nms,nms_delta,cma_mode
  write(6,*) "negative_cut",negative_cut
  write(6,*) "reduce mpi_size for saving memory (default =999) ",size_lim

  !$$$!! check write
  if(MPI__root) then
     do i=1,nqbz
        if(i<10 .OR. i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
        if(i==10 .AND. nqbz>18) write(6,"('... ')")
     enddo
     write(6,*)' !!nqbz nqibz =',nqbz,nqibz
  endif

  !! EFERMI
  !c need for gettetwt; automatically read file(okumura) OK
  !! tetrakbt mode (usetetrakbt)
  call getkeyvalue("GWinput","tetrakbt",usetetrakbt,default=.false.)
  if (usetetrakbt) then
     call readefermi_kbt()  !!! ef_kbt: Fermi energy at finite temperature
     ef = ef_kbt
  else
     call readefermi()      !!! ef:     Fermi energy at 0 K
     ef = ef_read
  endif

  tpioa=2d0*pi/alat

  !! QforEPS mode: 19Feb, 2020
  ! ccccccccccccccccccccccccccccccccccccccccccccccccc
  !! --- Readin Offset Gamma --------
  !$$$      if(debug) write(6,*) 'reading QOP'
  !$$$      open (newunit=ifq0p,file='Q0P')
  !$$$      read (ifq0p,"(i5)") nq0i
  !$$$      write(6,*) ' ### nqibz nq0i=', nqibz,nq0i
  !$$$      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i),epslgroup(nq0i) )
  !$$$      do i=1,nq0i
  !$$$        read (ifq0p, * ) wqt(i),q0i(1:3,i),idammy,epslgroup(i)
  !$$$      enddo
  !$$$      close(ifq0p)
  !$$$      nq0ix = nq0i
  !$$$      do i=1,nq0i
  !$$$        if(wqt(i)==0d0 ) then
  !$$$          nq0ix = i-1
  !$$$          exit
  !$$$        endif
  !$$$      enddo
  !$$$      neps = nq0i - nq0ix  ! number of zero weight q0p
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
  ! ccccccccccccccccccccccccccccccccccccccccccccccccc

  !! --- okumura Read dimensions of hamiltonian_wannier, spin, nqtt
  call read_wandata() ! nwf, nsp_w,nqtt_w
  call readscr(nwf,scrw)
  nnwf=nwf*nwf
  !      endif
  !$$$      !Weight for irreducible q-point (qibz)
  !$$$      do iq=1,nqibz
  !$$$         write(6,"('wibz',4f9.4)") wibz(iq),qibz(:,iq)
  !$$$      enddo


  !! Screening W for magnon
  scrw(:,:)=scrw(:,:)/hartree
  iqxend = nqibz !+ nq0i

  !! Shift W by hand (Oct.02, 2019)
  if (cma_mode) then
     allocate(scrw_original(nnwf,nnwf))
     scrw_original = scrw
  endif

  !! Initialization of readEigen !readin m_hamindex
  !      ginv = transpose(plat)
  call minv33tp(plat,qlat)
  if(verbose()>50) print *,'eeee exit of init_readeigen2'

  do iq=1,nqibz
     iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
     !     nstibz(iq) = nstbz(iqbz)
     write(6,"(' iq qibz nstibz=',2i5,3f9.4,i5)")iq,iqbz,qibz(:,iq) !,nstibz(iq)
  enddo

  !! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
  wemax=    5d0 !max value for plot
  omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
  !      if(MPI__root) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max

  !! getfreq returun date given at " use m_freq,only:".
  !      if(.not.epsmode) call getwemax(lqall,wemax) !wemax is to determine nw !real axis divisions

  !! NOTE: npmtwo=T sets npm=2
  !! optional npmtwo is added aug2017
  !! 20190604 Im[K]
  !$$$      if (negative_cut) then
  !$$$         call getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,MPI__root, npmtwo=.false.)
  !$$$      else
  !$$$         call getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,MPI__root, npmtwo=.true.)
  !$$$      endif
  niw=niw_in
  if( .NOT. imagomega) niw=1  !dummy
  call Getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw,ua, npmtwo=.true.)!,tetra
  !! Write freq_r
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

  if(MPI__root) write(6,"(' nw npm=',2i5)") nw,npm

  !! Tetra hedron initialization
  ! f (wan) then
  noccxv = nwf
  print *,"noccxv ",noccxv
  noccx  = noccxv + nctot
  nprecx = ndble  !We use double precision arrays only.
  mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
  nspinmx = nspin
  iqxini=1
  mtet=(/1,1,1/) !dummy
  !$$$      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/))
  ! kumura

  !! Set iqxini !omitqbz means skip loopf for iq=1,nqibz !19Feb,2020
  if(omitqbz) then
     iqxini= nqibz + 1
  else
     iqxini= 1
  endif
  iqxend = nqibz + nq0i

  write(6,"('iqxini,iqxend',2I8)") iqxini,iqxend
  do iq = iqxini,iqxend
     write(6,"('iq, qibze:',I8,3f9.4)") iq-iqxini,qibze(:,iq)
     autogamma=.true.
  enddo

  !! === Use of symmetry. EIBZ procedure PRB81,125102 ===
  !! EIBZ mode memo for nolfco (right?)
  !! If eibzmode=T, it is efficient but can slightly break crystal symmetry.(how much?)
  !! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
  !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
  !     okumura
  !eibzmode=.false.          !simple symmetry
  !$$$c     end okumura
  iqxendx=iqxend
!  allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx) &
!       ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
!  nwgt=1

  !! Wannier eigenvalues are read
  do is=1,nspinmx
     ! all wan_readeigen(qbz(:,:),nqbz,is) !QforEPSL ver. (20Feb, 2020)
     call write_qdata(ginv,nqbz,qbz(:,:))
  enddo

  !! Gamma point is automatically added if qibze does not include.
  if ( .NOT. sum(abs(qibze(:,iqxini)**2)) == 0d0) then
     write(6,*) "Gamma point is automatically added"
     iqxini=iqxini-1 !! Gamma: iqxini-1
     addgamma=.true.
  endif
  nqsym=iqxend-iqxini+1
  write(6,"('nqsym:',I4)") nqsym

  !! mpi processing for QforEPSL
  call MPI__hmagnon_rankdivider(nqsym)
  write(6,*) "nqsym mpi_MEq(:)",nqsym,mpi__MEq
  allocate(rpa_maximr(mpi__MEq),mf_maximr(mpi__MEq)) !list of w(MAX(Im[R])): magnon peak
  rpa_maximr=0d0; mf_maximr=0d0

!!! for MPI (finally nqsym ---> nqibz)
  do iq=1,nqsym
     if (MPI__task(iq)) write(6,'("iq,MPI_rank, mpi__ranktab",3I8)') iq,MPI__rank,mpi__ranktab(iq)
  enddo
  call MPI__barrier()

  !! ======== Loop over iq ================================
  do 1001 iqq = iqxini,iqxend      ! NOTE: q=(0,0,0) is active iqq=iqxini (see autogamma)
     iq = iqq-iqxini+1 !! start with iq=1 for convenience
     if(MPI__rank > size_lim) cycle !reduce mpi-size for test (skip 21-32)
     if( .NOT. MPI__task(iq) .AND. iq /= 1) cycle

     !! automatically set q=(0 0 0)
     if (autogamma .AND. iq==1) then
        q = 0d0
     else
        q = qibze(:,iqq)
     endif
     imaximr=imaximr+1      !!! imaximr=1,2,...,mpi__MEq

     !! Caution : confusing point
     !!  ngc by QGcou is shown at the bottom of lqg4gw.
     !!  ngc read from PPOVL are given by rdata4gw---> ngc(iq>nqibz )=ngc for q=0
     if( iq==1 ) then       ! *sanity check
        if(sum(q**2)>1d-10) then
           call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
        endif
     endif
!!! write qshort instead of q
     ! ibz: BZ weight
     !         wibz=0d0 !! not valid this version 20Feb,2020
     write(6,"('===== do 1001: iq wibz(iq) q=',i6,f13.6,3f9.4,' ========')") &
          iq,q!,wibz(iqlist(iq)),qshort !qq
     if (lhm) cycle

     !     ! zxq and zxqi are the main output after Hilbert transformation
     !     ! zxqi is not used in hmagnon (imagomega=.false.)
     if ( .NOT. allocated(zxq) )  allocate( zxq (nnwf,nnwf,nw_i:nw), zxqi (nnwf,nnwf,niw))
     if ( .NOT. allocated(zxq_d)) allocate( zxq_d(nnwf,nnwf) )
     zxq=0d0;  zxqi=0d0
     if(debug) write(6,*)' niw nw=',niw,nw

     !     ! ==== spin chi_charge or chi_+- ====
     is=1
     isf=2
     !! Tetrahedron weight.
     !! output
     !!     nbnbx
     !!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
     !!     nhw(ibjb,kx): the number of histogram sections
     !!     jhw(ibjb,kx): pointer to whw
     !!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
     !!     : histogram weights for given ib,jb,kx for histogram sections
     !!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
     !            write(6,*) ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2

     ! ccccccccccccc simple band for check
     if( .NOT. allocated(ev_w1)) allocate(ev_w1(nwf,nqbz))
     if( .NOT. allocated(ev_w2)) allocate(ev_w2(nwf,nqbz))
     if( .NOT. allocated(evc_w1)) allocate(evc_w1(nwf,nwf,nqbz))
     if( .NOT. allocated(evc_w2)) allocate(evc_w2(nwf,nwf,nqbz))

     write(6,*) "mpm nctot",npm,nctot

     do 5001 kx=1,nqbz      !!! ev_w1, ev_w2 unit: [Ry]
        call wan_readeval2(  qbz(:,kx), is, &
             ev_w1(1:nwf,kx), evc_w1(1:nwf,1:nwf,kx))
        call wan_readeval2(q+qbz(:,kx), isf, &
             ev_w2(1:nwf,kx), evc_w2(1:nwf,1:nwf,kx))

        !! only Cu2MnAl (cma)
        !! Energy of Mn3d(dn) is moved by cma_shitf
        if (cma_mode) then
           !$$$               if (iq==1) write(6,*) "cma_mode: True"
           !$$$               if (iq==1) write(6,"('cma_mode: iwf_s, iwf_e',2i4)") cma_iwf_s,cma_iwf_e
           !$$$               if (iq==1) write(6,"('cma_mode: cma_up_shift, cma_dn_shift',2E13.5,' [eV]')") cma_up_shift,cma_dn_shift
           do iwf = 1, nwf
              if ( cma_iwf_s <= iwf .AND. iwf <= cma_iwf_e ) then
                 ev_w1(iwf,kx) = ev_w1(iwf,kx) + cma_up_shift/rydberg()
                 ev_w2(iwf,kx) = ev_w2(iwf,kx) + cma_dn_shift/rydberg()
              endif
           enddo
        endif
        !! end cma
5001 enddo
     write(6,"('= start wan_gettetwt =',2i6,3f9.4)") nwf,iq,q

     !! tetrahedron weight
!call gettetwt(q,iq,isdummy,isdummy,nwgt(:,iq),ev_w1,ev_w2,nwf,eibzmode,wan)
     call gettetwt(q,iq,isdummy,isdummy,ev_w1,ev_w2,nwf,wan)
     write(6,*) "=== end gettetwt"
     deallocate(ev_w1,ev_w2)
!!!! normalization of Im[K]:
     rnqbz=1/dble(nqbz)
     znorm=-1d0*pi
!!!! generate wanmat
     if ( .NOT. allocated(kmat)) allocate(kmat(1:nnwf,1:nnwf,1:nwhis,1:npm), &
          wanmat(1:nnwf,1:nnwf)) !! (nwf^2*nwf^2)
     !           kmat=(0d0,1d-15)
     kmat=(0d0,0d0)       !orbtest
     do 2011 kx=1,nqbz    !! discrete k-point loop
        do 2012 jpm=1,npm !! jpm=2: negative frequency
           do 2013 ibib=1,nbnb(kx,jpm) !! n,n' band loop
!!!     n1b(ibib,k,jpm) = n :band index for k (occupied)
!!!     n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
              if (debug) then
                 write(6,*) 'jpm,ibib',jpm,ibib
                 write(6,*) 'kx,nbnb',kx,nbnb(kx,jpm)
                 if (jpm==1) then
                    write(6,*) 'n1b(occ),n2b(unocc):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                 else
                    write(6,*) 'n1b(unocc),n2b(occ):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                 endif
              endif

              !==========================caliculate M(:,:) : wanmat(nwf,nwf)
              ! c   matrix elements from wannier eigenvector
              !$$$  !!! evc_w1: upspin ; evc_w2: downspin
              it=n1b(ibib,kx,jpm) !index for n  for q
              itp=n2b(ibib,kx,jpm) !index for n' for q+k
              wanmat=(0d0,0d0)
              ijwf=0
              !$$$ O(nwf^4)
              do 3002 iwf=1,nwf
                 do 3003 jwf=1,nwf
                    ijwf=ijwf+1
                    klwf=0
                    do 3004 kwf=1,nwf
                       do 3005 lwf=1,nwf
                          klwf=klwf+1
                          if ( ijwf <= klwf ) then !!! Herimite
                             !                   call checkorb2(iwf,jwf,kwf,lwf,ijklmag)
                             !                   ijklmag=.true. !test for FeFe
                             !                   print *,"checkorb2:",iwf,jwf,kwf,lwf,ijklmag !! ijklmag=T: {ijkl} in same atom
                             !                   if (ijklmag) then !!! same atom
                             !                 elseif (.true.) then !!! same atom
!!! calculate numerator of Kmatrix

                             wan_j=dconjg(evc_w1(jwf,it,kx)) !a_{Rj   beta}^{kn}*
                             wan_i=evc_w2(iwf,itp,kx) !a_{Ri  alpha}^{(k+q)n'}
                             wan_l=evc_w1(lwf,it,kx) !a_{R'l  beta}^{kn}
                             wan_k=dconjg(evc_w2(kwf,itp,kx)) !a_{R'k alpha}^{(k+q)n'}*
                             wanijkl=wan_j*wan_i*wan_k*wan_l

                             wanmat(ijwf,klwf)=wanijkl
                             wanmat(klwf,ijwf)=dconjg(wanijkl) !!! Suppose Hermite Kmat

                             !$$$                      write(ifwanmat,"('i j k l ijkltag abs(wan) Re(wan) Im(wan)',4I3,L4,5E13.5)") iwf,jwf,
                             !$$$     &                     kwf,lwf,ijklmag,abs(wanijkl),wanijkl
!!! wanmat (dimension:nnwf)
                          else
                             wanijkl=(0d0,0d0)
                          endif
                          !$$$                 if (abs(aimag(wanijkl))<1d-16) then
                          !$$$                    wanijkl=cmplx(dble(wanijkl),0d0,kind(0d0))
                          !$$$                 endif
                          !                 wanmat(ijwf,klwf)=wanijkl+(0d0,1d-9) ! +delta (~ 1d-9?)
                          ! delta makes a width of Im[R] peaks

                          continue     !lwf
3005                   enddo
                       continue     !kwf
3004                enddo
                    continue     !jwf
3003             enddo
                 continue     !iwf
3002          enddo
              debug=.False.
              !========================== end caliculate M(:,:) : wanmat(nwf,nwf)
              ! c   print *,"ihw,nhw",ihw(ibib,kx,jpm),nhw(ibib,kx,jpm)
              do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
                 imagweight=whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
!!!   accumulate Im[K]
                 kmat(:,:,iw,jpm)=kmat(:,:,iw,jpm)+imagweight*wanmat(:,:)
              enddo

              continue                  !ibib-loop
2013       enddo
           continue                  !jpm-loop
2012    enddo
        continue                  !k-loop
2011 enddo
     deallocate(evc_w1,evc_w2,wanmat)
     call tetdeallocate()      !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)

     !! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
     !     ! npm is set by npmtwo=T
     schi=1                    ! flip over maj/min spins (<0 flip??).
     if (negative_cut) kmat(:,:,:,2)=(0d0,0d0)
     !! kmat ---> zxq
     call dpsion5( realomega, imagomega, kmat, nnwf,nnwf, zxq, zxqi, &
          .false., schi,1,1d99,1d99) 
     deallocate(kmat)
     deallocate(zxqi) !!! not used in hmagnon (2018/06/20)
!!!Hermite for Kmat (zxq)
     if (lhm) then
        allocate(zxq2(nnwf,nnwf,nw_i:nw))
        zxq2=zxq
        zxq=0d0
        ijwf=0
        do 3006 iwf=1,nwf
           do 3007 jwf=1,nwf
              ijwf=ijwf+1
              klwf=0
              do 3008 kwf=1,nwf
                 do 3009 lwf=1,nwf
                    klwf=klwf+1
                    if ( ijwf < klwf ) then !!! Herimite
                       continue

                    elseif (ijwf == klwf) then
                       zxq(ijwf,ijwf,:)= zxq2(ijwf,ijwf,:)
                    else             !!! ijwf > klwf
                       zxq(ijwf,klwf,:)=( zxq2(ijwf,klwf,:) + dconjg(zxq2(klwf,ijwf,:)) )/2.0
                       zxq(klwf,ijwf,:)=dconjg( zxq(ijwf,klwf,:) )
                    endif
3009             enddo
3008          enddo
3007       enddo
3006    enddo
        ijwf=0; klwf=0
        deallocate(zxq2)
     endif

!!! 20190604 negative cut for Im[K]
!$$$      if (negative_cut) then
!$$$         do iw=nw_i,0
!$$$            do iwf=1,nnwf
!$$$               do jwf=1,nnwf
!$$$                  zxq(iwf,jwf,iw)=cmplx(dble(zxq(iwf,jwf,iw)),0d0,kind(0d0))
!$$$               enddo
!$$$            enddo
!$$$         enddo
!$$$      endif


!!! threshold for Im[K] (zxq)
threshold=.True.
if (threshold) then
do iw=nw_i,nw
do iwf=1,nnwf
  do jwf=1,nnwf
     if (log10(abs(aimag(zxq(iwf,jwf,iw)))) < -15) then
        zxq(iwf,jwf,iw)=cmplx(dble(zxq(iwf,jwf,iw)),0d0,kind(0d0))
        !$$$  if ( aimag(zxq(iwf,jwf,iw)) > 0d0 ) then
        !$$$                     zxq(iwf,jwf,iw)=cmplx(dble(zxq(iwf,jwf,iw)),-1d0*aimag(zxq(iwf,jwf,iw)),kind(0d0))
     endif
  enddo
enddo
enddo
endif
threshold=.False.
allocate(wkmat(1:nnwf,1:nnwf),wkmat2(1:nnwf,1:nnwf)) !WKmatrix, WKmatrix_inv
allocate(rmat(1:nnwf,1:nnwf,nw_i:nw)) !Rmatrix
wkmat=0d0;wkmat2=0d0;rmat=0d0

!!! get eta for (1-eta*WK)
if (iq==1) then

if ( .NOT. allocated(eval_wk)) allocate(eval_wk(nnwf))
if ( .NOT. allocated(eval_wk2)) allocate(eval_wk2(nnwf))
if ( .NOT. allocated(eval_k)) allocate(eval_k(nnwf))

wkmat(1:nnwf,1:nnwf) &
 =matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,0)) !omega=0

if (write_hmat) call writehmat(scrw,nwf,"wmat_check.dat")

!!!   eval_wk is complex array because of Non-Hermite WK
call diagcvuh3(wkmat(:,:),nnwf,eval_wk)
! all diagwan(wkmat(:,:),eval_wk)

if (debug) then
write(6,"('eigenvalue of WK real&imag =',2E15.5)" &
    ,advance='NO') eval_wk
write(6,*)
endif

!         eta= 1d0/minval(dble(eval_wk)) !moderate peak
eta=-1d0/maxval(abs(eval_wk))
!         eta=  -1d0/abs(eval_wk(1)) !sharp peak
write(6,*) "now eigenvalue abs(WK)",abs(eval_wk(1)),"is inversed"
write(6,*) "check eigenvalue Re(WK)",real(eval_wk(1))
write(6,*) "check eigenvalue Im(WK)",aimag(eval_wk(1))
write(6,*) "wkmat calculated eta:", eta !negative value
endif                     !iq==1

! cc open file each iq
if (MPI__task(iq)) then
open(newunit=ifchipmz_wan,file="wan_ChiPMz.mat"//charnum4(iq))
open(newunit=ifchipmr_wan,file="wan_ChiPMr.mat"//charnum4(iq))
print *,'ifchipm=',ifchipmz_wan,ifchipmr_wan
!!!
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
!$$$         ifchipmrk_wan=ifile_handle()
!$$$         open(ifchipmrk_wan,file="wan_ChiPMr-k.mat"//charnum4(iq))
endif

!!! unit matrix (dimension: nwf*nwf), need for 1-WK matrix
if ( .NOT. allocated(imat)) then
allocate(imat(1:nnwf,1:nnwf))
imat=0d0
if (nms) then
do iwf=1,nnwf
  imat(iwf,iwf)=cmplx(1d0,nms_delta,kind(0d0))
enddo
write(6,*) "imat(iwf,iwf)=",imat(1,1)
else
do iwf=1,nnwf
  imat(iwf,iwf)=(1d0,0d0)
enddo
endif
!$$$         imat=0d0
!$$$         call wan_imat(nwf,imat) !

endif
! ccccc
if (l1wkout) then
if (MPI__task(iq)) then
   !$$$  ifwkeigen2=ifile_handle()
   !$$$  open(ifwkeigen2,file="wk_eval_list.dat"//charnum4(iq))
open(newunit=ifwkeigen3,file="1wk_eval_list.dat"//charnum4(iq))
endif
endif

! cccccdiagonalization for K
print *,"nw_i,nw",nw_i,nw
maximr=0d0; w_maximr=0d0 !!! search for Im[R] peak
do 2050 iw = nw_i,nw
trmat=(0d0,0d0)

!!! diag
call diagcvuh3(zxq(:,:,iw),nnwf,eval_wk)
!!! Hermite matrix diagonization
! all diagcvh2(zxq(:,:,iw),nnwf,eval_wk)
!!!

do iwf=1,nnwf
if (abs(aimag(eval_wk(iwf))) < 1d-16) eval_wk(iwf)=cmplx(dble(eval_wk(iwf)),0d0,kind(0d0))
enddo
!         call diagwan_tr(zxq(:,:,iw),trmat)
!!!
www=freq_r(iw)
if(iw<0) www=-freq_r(-iw)
if( .NOT. gskip .AND. MPI__task(iq)) then
   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Check for Kramers-Kronig relatioin
   ! Some Error: remove this section or modified ... (Okumura, Oct02,2019)
if (iw==nw_i) then !only first line
  trmat2=0d0
  write(6,"('--check iww sum(zxq)',i5,2E13.5)") &
       iw,sum((zxq(:,:,:)))
  do iww=nw_i,nw
     if (iww==0) cycle !skip w=0 (Cauthy principle integral)
     call diagcvuh3(zxq(:,:,iww),nnwf,eval_wk2)
     do iwf=1,nnwf
        if (abs(aimag(eval_wk2(iwf))) < 1d-16) eval_wk2(iwf)=cmplx(dble(eval_wk(iwf)),0d0,kind(0d0))
     enddo

     www2 = freq_r(iww)
     if (iww<0) www2 = -freq_r(iww)
     trmat2 = trmat2 + &
          sum(eval_wk2)/(znorm)*www2*abs(freq_r(abs(iww))-freq_r(abs(iww)-1))
  enddo
  call diagcvuh3(zxq(:,:,0),nnwf,eval_wk2) !for Re[K(w=0)]
  write(ifchipmz_wan,"('# int ImK/omega dw, Re[K(0)]=',2E13.5)") aimag(trmat2),sum(real(eval_wk2))
endif
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,2E17.9)") q, &
                                !     &        iw,www*2d0,trmat/hartree
    iw,www*2d0,hartree*sum(eval_wk)/(znorm)
!!! Im[K] [1/Ry] ?
!!! write d-d (diagonal) and d-other (non-diagonal)
if (output_ddmat) then
  if (iw==1000) then
!!! Note: writeddmat(matrix, nwf, nw_i, nw, filename, diagonal or non-diagnal)
     call writeddmat(zxq(:,:,iw),nwf,"wan_ChiPMr.mat.dd",.true.,zxq_d(:,:))
     ! all writeddmat(zxq(:,:,:),nwf,nw_i,nw,"wan_ChiPMr.mat.do",.false.)
     call writehmat(zxq_d(:,:),nwf,"zxqdmat_check.dat")
     write(6,*) "PASS for zxqmat_check"
     deallocate(zxq_d)
  endif
endif
endif

!!! K/(1-WK) = K(1-WK)^(-1)
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c W shift (q is far from Gamma)
!c W shift for Cu2MnAl
if (cma_mode) then
scrw = scrw_original
if (sqrt(dot_product(q,q)) > 0.5) then
  do jwf=1,nwf
     if ( .NOT. cma_iwf_s <= jwf .AND. jwf <= cma_iwf_e) cycle
     ijwf_j=(jwf-1)*nwf+iwf
     scrw(ijwf_j,ijwf_j)=scrw_original(ijwf_j,ijwf_j) + cmplx(dble(cma_wshift/hartree),0d0,kind(0d0))
     ! write (6,"('wshift (iwf,scrw):'i6,f9.4)",advance="no") ijwf_j,abs(scrw(ijwf_j,ijwf_j))
  enddo
endif
endif
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!! WK matrix
wkmat(1:nnwf,1:nnwf) &
 =eta*matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,iw))

!??     &        =eta*matmul(scrw(1:nnwf,1:nnwf),abs(zxq(1:nnwf,1:nnwf,iw)))
! cc Hermite check (WKmat)
! cc  For developing program
if (write_hmat) then
if (iq==4 .AND. iw==205) then
  eval_wk=0d0
  write(6,*) "wkmat check"
  call writehmat(wkmat,nwf,"wkmat_check.dat")
  if (lsvd) then
     call zgesvdnn(nnwf,zxq(:,:,iw),SS,UU,VT)
     do iwf=1,nnwf
        write(6,'("wkmat eigenvalue:",f9.4)') SS(iwf)
     enddo
  endif
endif
endif
! cc
if (debug) then
call diagcvuh3(wkmat(:,:),nnwf,eval_wk)
if (iw==45) then
  write(6,"('check sign eval WK:omega=45 case, iq=',i4)") iq
  write(6,"('check_sign eval WK:',2E13.5)" &
       ,advance='NO') eval_wk(1:10)
endif
call diagwan(wkmat(:,:),eval_wk)
if (iw==45) then
  write(6,"('check sign eval WK:omega=45 case, iq=',i4)") iq
  write(6,"('check_sign eval WK:',2E13.5)" &
       ,advance='NO') eval_wk(1:10)
endif
endif
!!! 1-eWK ===> 1-WK
wkmat2(1:nnwf,1:nnwf)=imat(1:nnwf,1:nnwf)-wkmat(1:nnwf,1:nnwf)
! c Hermite check for 1-eWK for developing code
if (write_hmat) then
if (iq==4 .AND. iw==205) then
  call writehmat(wkmat2,nwf,"1wkmat_check.dat")
endif
endif

!!! eigenvalue of 1-WK
if (debug) then
if (iq==1 .AND. iw==0) then
  call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
  !               call diagwan(wkmat2(:,:),eval_wk)
  write(6,*) "q=0 and omega=0 case:"
  write(6,"('eigenvalue of (1-eta*WK) real&imag=',2E13.5)" &
       ,advance='NO') eval_wk
  write(6,*)
endif
endif
! cc
if (debug) then
call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
if (iw==45) then
  write(6,"('check sign eval 1-WK:omega=45 case, iq=',i4)") iq
  write(6,"('check_sign eval 1-WK:',2E13.5)" &
       ,advance='NO') eval_wk(1:10)
endif
call diagwan(wkmat2(:,:),eval_wk)
if (iw==45) then
  write(6,"('check sign eval 1-WK:omega=45 case, iq=',i4)") iq
  write(6,"('check_sign eval 1-WK:',2E13.5)" &
       ,advance='NO') eval_wk(1:10)
endif
endif
! cc

! c   WK eigenvalue writing (wk_eval_list.dat)
!$$$         call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
!$$$         if (MPI__task(iq)) write(ifwkeigen2,"(3f9.5,I6,10E13.5)") qshort,iw,www*2d0,
!$$$     &        eval_wk(1),abs(eval_wk(1))!,eval_wk(2),abs(eval_wk(2)),eval_wk(3),abs(eval_wk(3))
! c

call matcinv(nnwf,wkmat2(1:nnwf,1:nnwf)) ! inv(1-WK)
rmat(1:nnwf,1:nnwf,iw)=matmul(zxq(1:nnwf,1:nnwf,iw), &
 wkmat2(1:nnwf,1:nnwf))

!!! eigenvalue of (1-WK)inv
if (debug) then
if (iq==1 .AND. iw==0) then
  call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
  ! c
  !               call diagwan(wkmat2(:,:),eval_wk)
  ! c
  write(6,*) "q=0 and omega=0 case:"
  write(6,"('eigenvalue of (1-eta*WK)inv real&imag=',2E13.5)" &
       ,advance='NO') eval_wk
  write(6,*)
endif
! c   find eigenvalue of (1-WK) at all omega (2018/04/24)
if (iw==0) then
   !               call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
  write(6,"('WKWK: eig 1-WK inv Re',3f9.4,f15.8)") q, dble(eval_wk(1))
  write(6,"('WKWK: eig 1-WK inv Im',3f9.4,f15.8)") q,aimag(eval_wk(1))
endif
endif
! c 1-WK eigenvalue writing (1wk_eval_list.dat)
!          if (iq==1 .and. iw .ge. 0) then
!          if (iw .ge. 0) then ! omega >= 0

call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
!         call diagwan(wkmat2(:,:),eval_wk)

!!! eigenvalue of rmat
if (debug) then
if (iq==1 .AND. iw==0) then
   !               call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
  write(6,*) "q=0 and omega=0 case:"
  write(6,"('eigenvalue of rmat',i5,2f15.8)") iw,eval_wk(1)
  write(6,*)
endif
endif

! cc
if (debug) then
if (iw==45) then
  write(6,"('check sign eval Rmat: omega=45 case, iq=',i4)") iq
  write(6,"('check_sign eval Rmat:',2E13.5)" &
       ,advance='NO') eval_wk(1:10)
  call diagwan(wkmat2(:,:),eval_wk)
  write(6,"('check sign eval Rmat: omega=45 case, iq=',i4)") iq
  write(6,"('check_sign eval Rmat:',2E13.5)" &
       ,advance='NO') eval_wk(1:10)
endif
endif
! cc
if (l1wkout) then
if (MPI__task(iq)) write(ifwkeigen3,"(3f9.5,I6,6E13.5)") q,iw,www*2d0, &
    eval_wk(1)     !,abs(eval_wk(1)),eval_wk(2),abs(eval_wk(2)),eval_wk(3),abs(eval_wk(3)) !iwf=1,2,3
endif


!!! diagonalization for R
!         trmat=0d0
!         trmatt=0d0
!$$$         do iwf=1,nwf
!$$$            ijwf=(1+nwf)*iwf-nwf !!! iwf=jwf
!$$$            do jwf=1,nwf
!$$$               klwf=(1+nwf)*jwf-nwf !!! kwf=lwf
!$$$c$$$               call checkorb2(iwf,iwf,jwf,jwf,ijklmag) ; ijklmag=.true.
!$$$c$$$               if(ijklmag) then
!$$$               trmat=trmat+rmat(ijwf,klwf,iw)
!$$$c               trmatt=trmatt+(rmat(ijwf,klwf,iw)-zxq(ijwf,klwf,iw))
!$$$               !endif
!$$$            enddo
!$$$         enddo

!         write(6,'("trmat_check3",81E12.4)') aimag(eval_wk)
call diagwan(rmat(:,:,iw),eval_wk)
trmat=sum(eval_wk(1:nnwf))

call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
!$$$         write(6,'("trmat_check diagwan and diagcvuh",3f9.4,i6,3E12.4)') qshort,iw,
!$$$     &        aimag(trmat),aimag(sum(eval_wk(1:nnwf))),aimag(sum(eval_wk(:)))
! cc
if (MPI__task(iq)) then
if( .NOT. gskip) write(ifchipmr_wan,"(3f9.4,i6,E12.4,x,12E12.4)")q,iw, &
    www*2d0,hartree*trmat!, wibz(iqlist(iq))!,
!     &           www*2d0,sum(eval_wk(1:nnwf))/hartree
!            if(.not. gskip) write(ifchipmrk_wan,"(3f9.4,i6,E12.4,x,12E12.4)")qshort,iw,
!     &           www*2d0,trmatt/hartree

! cc search for MAX(Im[R]) 20180706 (0 - 1500 meV)
!            if (0d0 <www*hartree .and. www*hartree < 1.5) then
if (0d0 < www) then
  if (maximr < -1d0*aimag(trmat)) then
     maximr=-1d0*aimag(trmat)
     w_maximr=www
  endif
endif
! cc
endif
! c END K/(1-WK)

gskip=.false.          !!! Skip Gamma point
! cc Hermite check (zxq and rmat) for developing code
if (write_hmat) then
if (iq==4 .AND. iw==205) then
  call writehmat(zxq(:,:,iw),nwf,"kmat_check.dat")
  call writehmat(rmat(:,:,iw),nwf,"rmat_check.dat")
endif
endif
! cc

2050 enddo
deallocate(zxq, rmat, wkmat, wkmat2)

if ( .NOT. iq==1) then
   ! c(MF)
   ! BZweight*omega[eV] for sum(E(q))/N
mf_maximr(imaximr)=wibz(iq)*w_maximr*hartree
! c(RPA)
! BZweight*omega[eV] for sum(1/E(q))/N
rpa_maximr(imaximr)=wibz(iq)/(w_maximr*hartree)
endif

if (MPI__task(iq)) then  !magnon peak
write(6,"('AAAA',I4,3f9.4)") iq,q
write(6,"('AAAA, w(MAX(im[R])), Im[R]',f13.5,E12.4)") w_maximr*2d0*rydberg()*1000,maximr/hartree
endif

if (MPI__task(iq)) then
   !$$$         write(ifwkeigen2,*)
if (l1wkout) write(ifwkeigen3,*)
write(ifchipmz_wan,*)
write(ifchipmr_wan,*)
!         write(ifchipmrk_wan,*)
endif

!$$$  !! ImagOmega end =================
debug=.False.
if (MPI__task(iq)) then
   !$$$         close(ifwkeigen2)
if (l1wkout) close(ifwkeigen3)
close(ifchipmz_wan)
close(ifchipmr_wan)
!         close(ifchipmrk_wan)
endif
continue                  !q point loop
1001 enddo
write(6,*) "maximr RPA",rpa_maximr
write(6,*) "maximr MFA",mf_maximr
sumrpa_maximr=cmplx(sum(rpa_maximr(:)),0d0,kind(0d0))
summf_maximr =cmplx(sum( mf_maximr(:)),0d0,kind(0d0))
!$$$!! =================== end of loop 1001 for q point ========================
call MPI__barrier()

write(6,*) "MPIcheck",MPI__size,MPI__rank,sumrpa_maximr
if(MPI__size/=1) then
call MPI__AllreduceSum(sumrpa_maximr(1),1)
call MPI__AllreduceSum( summf_maximr(1),1)
endif

write(6,"('sum(E(q)) for MFA:',f9.4)") real(summf_maximr(1))
write(6,"('[sum(1/E(q))]inv for RPA',f9.5)") 1d0/real(sumrpa_maximr(1))

call cputid(0)
!      call MPI__Finalize
write(6,"('eta for 1-eta*WK:',f13.8)") eta
call rx0( ' OK! hmagnon mode')
END PROGRAM



