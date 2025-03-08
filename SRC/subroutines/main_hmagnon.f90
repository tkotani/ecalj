!>  Calculate Chi^+-, spin susceptibility. 
module m_hmagnon 
  contains
subroutine hmagnon() bind(C)
  use m_readwan,only: write_qdata, wan_readeigen, wan_readeval, wan_readeval2, &
       readscr, checkorb, checkorb2, diagwan, diagwan_tr, wan_imat, &
       writehmat, writeddmat, read_wandata, nwf, nsp_w, nqtt_w
  use m_ReadEfermi,only: readefermi,ef_read=>ef, readefermi_kbt,ef_kbt
  use m_readeigen,only: readeval,init_readeigen,init_readeigen2
  use m_hamindex,only:qtt,nqtt
  use m_read_bzdata,only: read_bzdata,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,ginv, &
       qbz,wbz,qibz,wibz,qbzw, idtetf,ib1bz,idteti,  nstar,irk,nstbz &
       ,wqt=>wt,q0i,nq0i ,nq0iadd,ixyz,epslgroup,nq0ix,neps
  use m_genallcf_v3,only: genallcf_v3,natom,nspin,nl,nn, &
       ndima,nlnmx, nctot,niw_in=>niw, alat, delta,deltaw,esmr,clabl,iclass, &
       il,in,im,nlnm, plat, pos,ecore
  use m_keyvalue,only: getkeyvalue
  use m_freq,only: getfreq, frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw 
  use m_tetwt,only: tetdeallocate,gettetwt, whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb
  use m_w0w0i,only: w0w0i, w0,w0i ! w0 and w0i (head part at Gamma point)
  use m_readgwinput,only: ReadGWinputKeys
  use m_lgunit,only:m_lgunit_init,stdo
  use m_dpsion,only: dpsion5
  use m_kind,only:kindrcxq
  use m_mpi,only: MPI__Initialize_magnon, MPI__consoleout_magnon, MPI__AllreduceSum
  use m_mpi,only: MPI__rank=>mpi__rankMG,MPI__size=>mpi__sizeMG, MPI__root ,comm
  implicit none
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  xxx removed--> eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2), ndble=8,MPI__MEq
  integer:: iqbz,iqindx,iflegas,nmx,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1
  integer:: ifv,lxxilmx,ilm_r,nx_r,lb,nb,mb, ii,iy,ipl1,ixx
  integer:: jpm,ncc, istat,imb,imb1,imb2,nmbas_in
  integer::iwf,jwf,inwf,kwf,lwf,ijwf,klwf,nnwf,ijwf_j
  integer::ifgas,ifchipmz_wan,ifchipmr_wan,ifchipmrk_wan
  integer::imaximr=0,niw,ifif,ierr
  integer::maxocc2,ixc,iqxini,iqxend,iqxendx, i,mxx,ini,ix,is &
       ,iw,noccxv,noccx,iq,ngb,nprecx,nblochpmx,ifwd,nspinmx,ibas &
       ,kx,isf,job,ihis,ik,ibib,ib1,ib2,j, incwfin,  verbose
  integer:: ificlass,k , nbmx, nqbze, nqibze
  integer,allocatable:: imbas(:), imbas_s(:),iibas(:), nxx_r(:)
  real(8):: q(3),  qgbin(3),qx(3), ua=1d0 ! ua is a dummy.
  real(8) :: omg2max, wemax
  real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
  real(8)::schi=1d0
  real(8) :: eclda_bh,eclda_pz,wk4ec,faca, qs,qt,ww,muu, ddq(3)
  real(8)::qrot(3),cr=2d-5
  real(8)::maximr,w_maximr, ef,tpioa,absq,vcou1,vcou1sq, ebmx, nms_delta
  real(8):: wan_ecore(1), imagweight, www, www2
  real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
  real(8),allocatable::SS(:),rwork(:),ss0(:)
  real(8),allocatable:: vxcfp(:,:),     wgt0(:,:) 
  real(8), allocatable:: qbze(:,:), qibze(:,:)
  real(8),allocatable::ev_w1(:,:),ev_w2(:,:)
  real(8),allocatable:: rpa_maximr(:),mf_maximr(:) !! for MAX(Im[R])
  complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:),zxq2(:,:,:), zxq_d(:,:)
  complex(8) :: fff,img=(0d0,1d0)
  complex(8),allocatable:: UU(:,:),VT(:,:),work(:),ddd(:,:) &
       ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
  complex(8),allocatable:: UU0(:,:),VT0(:,:)
  complex(8),allocatable::evc_w1(:,:,:),evc_w2(:,:,:)
  complex(8),allocatable::wanmat(:,:) ,wkmat(:,:),wkmat2(:,:),rmat(:,:,:),rmat2(:,:),wmat_check(:,:)
  complex(8),allocatable:: scrw(:,:) ,cmat2(:,:), scrw_original(:,:) !screening W
  complex(8),allocatable:: eval_wk(:),eval_k(:),eval_wk2(:),trmat22(:)
  complex(8)::trmat,trmatt,trmat1,trmat2
  complex(8),allocatable::imat(:,:) !unit matrix for 1-WK
  complex(kindrcxq),allocatable::kmat(:,:,:,:)
  logical :: tetra=.true., usetetrakbt !  tetrahedron method
!  logical :: debug=.false.,  
  logical:: realomega=.true., imagomega=.true.,omitqbz=.false. !, noq0p
  logical ::  chipm=.false.,nolfco=.false.,epsmode=.false.,normalm=.false., crpa=.false. &
       ,autogamma=.false., l1wkout=.false., addgamma=.false.
  logical :: timereversal, testtimer,onceww ! Feb2006 time-reversal=off case
  logical::wan=.true.,lhm,lsvd,nms , ijklmag 
  logical, allocatable :: mpi__task(:)
  logical(8):: gskip=.false.
  character*3:: charnum3
  character*4:: charnum4
  character*5:: charnum5
  logical::npmtwo,diag=.false.,t2g
  complex(8)::wan_i,wan_j,wan_k,wan_l,wanijkl
  integer::igv
  real(8)::qlat(3,3),qsh1(3),qsh2(3),znorm,rnqbz, eta 
  integer:: it,itp,isdummy,lorb, size_lim=999
!!! q on symline
  real(8)::rlatp(3,3),xmx2(3),qqin(3)
  integer:: nlatout(3,48),nout,iout,nqsym
  integer,allocatable:: iqlist(:)
  integer::ifwkeigen,ifwkeigen2,ifwkeigen3 ! WK eigenvalue check
  logical(8)::lsmo,threshold=.true.
  logical:: negative_cut, write_hmat, output_ddmat
  logical:: cma_mode !cma_mode for Cu2MnAl only 2019/09/27
  real(8):: cma_up_shift, cma_dn_shift, cma_wshift
  integer(4):: cma_iwf_s,cma_iwf_e
  integer::ifwanmat, npm2
  complex(8)::sumrpa_maximr(1),summf_maximr(1)
  real(8),parameter::   pi = 4d0*datan(1d0),  fourpi   = 4d0*pi,  sqfourpi = sqrt(fourpi)
  hartree  = 2d0*rydberg()
  call m_lgunit_init()
!  call getkeyvalue("GWinput","mpi_size_lim",size_lim,default=999)
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
! need for gettetwt; automatically read file(okumura) OK   !! tetrakbt mode (usetetrakbt)
!  call getkeyvalue("GWinput","tetrakbt",usetetrakbt,default=.false.) !tetrakbt is not yet tested
!  if (usetetrakbt) then
!    call readefermi_kbt()  !!! ef_kbt: Fermi energy at finite temperature
!    ef = ef_kbt
!  else
    call readefermi()      !!! ef:     Fermi energy at 0 K
    ef = ef_read
!  endif
  tpioa=2d0*pi/alat
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
  call readscr(nwf,scrw)
  nnwf=nwf*nwf
  !Weight for irreducible q-point (qibz); do iq=1,nqibz; write(6,"('wibz',4f9.4)") wibz(iq),qibz(:,iq); enddo
  scrw(:,:)=scrw(:,:)/hartree !! Screening W for magnon
  iqxend = nqibz !+ nq0i
  !! Shift W by hand (Oct.02, 2019)
  if (cma_mode) then
    allocate(scrw_original(nnwf,nnwf))
    scrw_original = scrw
  endif
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
  noccxv = nwf
  noccx  = noccxv + nctot
  nprecx = ndble  !We use double precision arrays only.
  nspinmx = nspin
  iqxini = merge(nqibz + 1,1,omitqbz)
  iqxend = nqibz + nq0i
  write(6,"('iqxini,iqxend noccxv',3I8)") iqxini,iqxend,noccxv
  do iq = iqxini,iqxend
    write(6,"('iq, qibze:',I8,3f9.4)") iq-iqxini,qibze(:,iq)
    autogamma=.true.
  enddo
  iqxendx=iqxend
  do is=1,nspinmx  
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
  rankdivider: block 
    use m_mpi,only: mpi__sizeMG,mpi__rankMG
    integer :: iq,i,mpi__ranktab(1:nqsym)
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
  allocate(zxq (nnwf,nnwf,nw_i:nw),evc_w1(nwf,nwf,nqbz),evc_w2(nwf,nwf,nqbz)) !, zxq_d(nnwf,nnwf) )
  allocate(eval_wk(nnwf),eval_wk2(nnwf),eval_k(nnwf),trmat22(nw_i:nw))
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
      write(6,*) "mpm nctot",npm,nctot
      readeigen: do 5001 kx=1,nqbz      !!! ev_w1, ev_w2 unit: [Ry]
        call wan_readeval2(  qbz(:,kx), is,    ev_w1(1:nwf,kx), evc_w1(1:nwf,1:nwf,kx)) !eigenvalue eigenfunciton
        call wan_readeval2(q+qbz(:,kx), isf,   ev_w2(1:nwf,kx), evc_w2(1:nwf,1:nwf,kx))
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
      complex(8):: kmat(1:nnwf,1:nnwf,1:nwhis,1:npm), wanmat(1:nnwf,1:nnwf), zxqi(1,1,1)
      zxq=0d0
      rnqbz=1/dble(nqbz)
      znorm=-1d0*pi ! normalization of Im[K]:
      kmat=0d0
      kxloop:       do 2011 kx=1,nqbz 
        jpmloop:    do 2012 jpm=1,npm ! jpm=2: negative frequency
          ibibloop: do 2013 ibib=1,nbnb(kx,jpm) !! n,n' pair band index loop
            it=n1b(ibib,kx,jpm)  !index for n  for q   ! n1b(ibib,k,jpm) = n :band index for k (occupied),   
            itp=n2b(ibib,kx,jpm) !index for n' for q+k ! n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
            wanmat=0d0
            do concurrent(iwf=1:nwf,jwf=1:nwf, kwf=1:nwf,lwf=1:nwf)
              ijwf=(iwf-1)*nwf+jwf
              klwf=(kwf-1)*nwf+lwf 
              if ( ijwf > klwf ) cycle         ! calculate numerator of Kmatrix
              wan_j=dconjg(evc_w1(jwf,it,kx))  !a_{Rj   beta}^{kn}*
              wan_i=evc_w2(iwf,itp,kx)         !a_{Ri  alpha}^{(k+q)n'}
              wan_l=evc_w1(lwf,it,kx)          !a_{R'l  beta}^{kn}
              wan_k=dconjg(evc_w2(kwf,itp,kx)) !a_{R'k alpha}^{(k+q)n'}*
              wanijkl=wan_j*wan_i*wan_k*wan_l
              wanmat(ijwf,klwf)=wanijkl
              wanmat(klwf,ijwf)=dconjg(wanijkl) !!! Suppose Hermite Kmat ! wanmat (dimension:nnwf)
            enddo
            do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
              imagweight=whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
              kmat(:,:,iw,jpm)=kmat(:,:,iw,jpm)+imagweight*wanmat(:,:) !accumulate Im[K]
            enddo
2013      enddo ibibloop
2012    enddo jpmloop
2011  enddo kxloop
      call tetdeallocate()      ! --> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
      schi=1                    ! flip over maj/min spins (<0 flip??).
      if (negative_cut) kmat(:,:,:,2)=0d0
      call dpsion5( realomega, imagomega, kmat, nnwf,nnwf, zxq, zxqi,.false., schi,1,1d99,1d99)  !! kmat ---> zxq
    endblock GETzxq
    if(lhm) then !Enforce zxq Hermitian 
      allocate(zxq2(nnwf,nnwf,nw_i:nw),source=zxq)
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
    endif
    where(abs(dimag(zxq))<1d-15) zxq=dreal(zxq) ! threshold for Im[K] (zxq)
    allocate(wkmat(1:nnwf,1:nnwf),wkmat2(1:nnwf,1:nnwf),rmat(1:nnwf,1:nnwf,nw_i:nw),source=(0d0,0d0)) !WKmatrix, WKmatrix_inv
    GetEta: if (iq==1) then ! (1-eta*WK)
      wkmat(1:nnwf,1:nnwf) =matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,0)) !omega=0
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
      !print *,'ifchipm=',ifchipmz_wan,ifchipmr_wan
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
      call diagcvuh3(zxq(:,:,iw),nnwf,eval_wk)
!!! Hermite matrix diagonization        ! all diagcvh2(zxq(:,:,iw),nnwf,eval_wk)
      where(abs(dimag(eval_wk)) < 1d-16) eval_wk=dreal(eval_wk)
      www = merge(-freq_r(-iw),freq_r(iw),iw<0)
      if(MPI__task(iq)) then ! Check for KK relatioin ! Some Error: remove this section or modified ... (Okumura, Oct02,2019)
        if(iw==nw_i) then !only first line
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
        endif
        write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,2E17.9)") q,iw,www*2d0,hartree*sum(eval_wk)/(znorm)
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
      wkmat(1:nnwf,1:nnwf) =eta*matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,iw)) !! WK matrix
      wkmat2(1:nnwf,1:nnwf)=imat(1:nnwf,1:nnwf)-wkmat(1:nnwf,1:nnwf)        ! c Hermite check for 1-eWK for developing code
      call matcinv(nnwf,wkmat2(1:nnwf,1:nnwf)) ! inv(1-WK)
      rmat(1:nnwf,1:nnwf,iw)=matmul(zxq(1:nnwf,1:nnwf,iw), wkmat2(1:nnwf,1:nnwf))
      call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
      call diagwan(rmat(:,:,iw),eval_wk) !! diagonalization for R
      trmat = sum(eval_wk(1:nnwf))
      call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
      if(MPI__task(iq)) then
        write(ifchipmr_wan,"(3f9.4,i6,E12.4,x,12E12.4)")q,iw,www*2d0,hartree*trmat
        if (0d0 < www.and. maximr < -1d0*aimag(trmat)) then !search for MAX(Im[R]) 20180706 (0 - 1500 meV) if (0d0 <www*hartree .and. www*hartree < 1.5)
          maximr=-1d0*aimag(trmat)
          w_maximr=www
        endif
      endif
2050 enddo iwloop
    deallocate(rmat, wkmat, wkmat2)
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
