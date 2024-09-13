!! hhomogas, originally included in hmagnon.m.F by Okumura 
subroutine hhomogas()
  use m_readefermi,only: setefermi
  use m_genallcf_v3,only: Genallcf_v3
  use m_shortn3_qlat,only: shortn3_qlat
  use m_ftox
  use m_lgunit,only:stdo,m_lgunit_init
  use m_shortvec,only:shortvec,shortvecinitialize
  use m_dpsion,only: dpsion5
  use m_hamindex0,only:   readhamindex0, symgg=>symops, ngrp,alat,plat,qlat
  use m_read_bzdata,only: read_bzdata, nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,ginv,&
       dq_,qbz,wbz,qibz,wibz,qbzw, idtetf,ib1bz,idteti, nstar,irk,nstbz
  use m_mpi,only: MPI__Initialize,MPI__root,&!MPI__Finalize,
       MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank, & !MPI__size,&
       MPI__consoleout!,MPI__barrier!MPI__ranktab,
  !! frequency
  use m_freq,only: getfreq,& !NOTE: call getfreq generate following data.
       frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw !, frhis0,nwhis0 !output of getfreq
  !! tetwt
  use m_tetwt,only: tetdeallocate,gettetwt, whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
  use m_homoelectron,only: read_qgband, efermi_egas !for gsq 
  use m_readgwinput,only: ReadGwinputKeys!,egauss,ecut,ecuts,mtet,ebmx,nbmx,imbas
!  use m_shortn3,only: shortn3_initialize,shortn3
  implicit none
  integer::nctot=0,nspin=1,niw
  !! ------------------------------------------------
  !! We calculate chi0 by the follwoing three steps.
  !!  gettetwt: tetrahedron weights
  !!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
  !!  dpsion5: calculate real part by the Hilbert transformation from the Im part
  !!  eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
  !!-------------------------------------------------
  real(8):: q(3),  qgbin(3),qx(3)
  real(8):: ua=1d0 ! this is a dummy.
  integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
  integer:: ndble=8
  integer:: nwordr
  real(8),allocatable:: vxcfp(:,:), wqt(:), wgt0(:,:),q0i(:,:)
  integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
  complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
  real(8),allocatable :: eqt(:) !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:), qbze(:,:),qibze(:,:)
  !,ecore(:,:)  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
  complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:)
  !  tetrahedron method
  logical :: tetra=.true. !,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
  ! so as to reduce the memory usage.
  complex(8) :: fff,img=(0d0,1d0)
  complex(8),allocatable :: wwk(:,:,:)
  integer,allocatable ::    noccxvv(:),n2bminimum(:,:,:)
  real(8) ::qbzx(3),anfvec(3)
  logical :: debug
  integer,allocatable:: ibasf(:)
  real(8),allocatable :: transaf(:,:)
  logical :: realomega=.true., imagomega=.false.
  complex(8),allocatable:: epsi(:,:),gbvec(:),zzr(:,:),x0mean(:,:,:),zzr0(:)
  complex(8) :: epxxx,vcmean, vcmmmm
  complex(8),allocatable:: vcmmm(:)
  character*11 fileps
  character*11 fileps23
  character*16 filepsnolfc
  character*11  filele
  character(5) :: charnum5
  character(20):: xxt

  real(8) :: Emin, Emax,emin2,emax2
  real(8) :: omg2max,omg1max,wemax
  real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)

  integer::iopen,maxocc2,iclose, ixc,iqxini,iqxend,iqxendx, ifhbe, nprecb,mrecb,mrece,nlmtot,nqbzt, & !nband,
       nq0i,i,nq0ix,neps,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx,nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, &
       ifepscond,nxx,nw0,iw,ifinin,iw0,ifwwk,noccxv,noccx,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas&
       ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0dummy&
       ,kx,isf,kqxx,kp,job,noccxvx(2)=-9999,nwmax,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j
  real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm, dwry,dwh,omg_c,omg2

  integer:: incwfin,  verbose

  integer:: ngc,mrecg !bzcase, 
  real(8):: quu(3), deltaq(3)!,qq(3) !,qqq(3)=0d0
  logical,allocatable :: iwgt(:,:,:,:)
  complex(8),allocatable:: wgt(:,:,:)

  real(8),allocatable:: qbz2(:,:)
  logical :: qbzreg !if true, we use off-gamma mesh.
  integer:: nbcut,nbcut2

  integer,allocatable:: nstibz(:) !Nov2004 Miyake's tote
  real(8),allocatable:: ecqw(:,:) !,wiw(:)
  real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree,pi,efz,qfermi,alpha,voltot,ecelgas,efx,valn
  integer:: iqbz,iqindx,iflegas,nmx,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1 
  real(8) :: eclda_bh,eclda_pz,wk4ec,faca
  real(8),allocatable::    evall(:)
  complex(8),allocatable:: ovlpc(:,:),evecc(:,:)
  integer:: nev !,  ifdpin

  real(8),allocatable:: ecut(:),ecuts(:) ,totexc(:), trpv(:),trlog(:)
  integer:: necut,iecut

  integer:: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
  integer,allocatable:: nxx_r(:)
  real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
  character*3:: charnum3
  character*4:: charnum4
  complex(8),allocatable:: jcoup(:,:), mcm(:,:,:)
  real(8)::chg1,chg2,spinmom,schi=1d0

  complex(8),allocatable:: ovlp(:,:),evec(:,:),ovlpi(:,:)
  real(8),allocatable::eval(:)
  integer:: new,nmxx,ii,iy,ipl1,ixx

  complex(8),allocatable :: ppovl(:,:),oo(:,:),x0meanx(:,:),x0inv(:,:),ppovlzinv(:,:)
  real(8)::qxx(3),ssm
  ! svd. not used now
  real(8),allocatable::SS(:),rwork(:),ss0(:)
  complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:),ddd(:,:),vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
  integer::lwork,info,imin,ifzxq
  complex(8)::x0mx
  complex(8),allocatable:: UU0(:,:),VT0(:,:)

  logical ::  chipm=.false.,nolfco=.false. ,epsmode=.false.,normalm=.false., crpa=.false. 
  integer::  ife, idum4 !ifchipmn,ifchipm,
  real(8):: qs,qt,ww,muu, ddq(3)
  character*11 ::ttt
  integer:: nnmx,nomx

  ! Feb2006 time-reversal=off case
  logical :: timereversal, testtimer,onceww
  integer:: jpm,ncc
  real(8):: frr

  integer:: ipm,nrecoff

  real(8),allocatable:: ebb(:)
  logical :: evaltest !for a debug test
  character*300:: aline
  integer:: istat,nmbas,imb,imb1,imb2,nmbas_in
  integer,allocatable:: imbas(:), imbas_s(:),iibas(:)
  !...
  complex(8),allocatable:: am1(:),am2(:),mmat(:,:), x0mat(:,:),x0matinv(:,:),eiqrm(:)
  integer:: ifchipmn_mat, ifchipm_fmat !,ifchipm_mat
  integer::ifstoner,ifx,i1
  real(8):: Istoner,zz1,zz2,zz3,zz4,Istoner0,jzero2,dumm1,dumm2
  complex(8):: trr,trr0,trr1     , zzzx(4,4), zzzy(4,4),trrx,mmatx(4,4),denom(4,4)
  real(8),allocatable:: eee(:),mmnorm(:), asvec(:,:),ssv(:,:),sproj(:,:),sprojx(:,:), momsite(:)
  real(8):: eex(4),eey(4),qvv(3)
  integer :: ngb0,ifvcoud,idummy,ifepstinv,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x,ifiss,iq0x
  complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:)
  real(8),allocatable:: vcousq(:)
  real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq

  !! Eq.(40) in PRB81 125102
  complex(8),allocatable::sk(:),sks(:),skI(:),sksI(:), w_k(:),w_ks(:),w_kI(:), w_ksI(:), s_vc(:),vw_k(:),vw_ks(:)
  complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
  integer:: lxklm,nlxklm,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
  complex(8):: vc1vc2
!  integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)

  real(8),allocatable:: aik(:,:,:,:)
  integer,allocatable:: aiktimer(:,:)
  integer:: l2nl
!  logical:: eibz4x0,tiii,iprintx,symmetrize,eibzmode
  real(8):: qread(3),imagweight,tot_imagweight

  character(128):: vcoudfile,aaax
  integer:: src,dest
  integer,allocatable :: iclasst(:), invgx(:)
  integer:: ificlass,ifile_handle,k
  complex(8),allocatable:: ppovl_(:,:)

  logical:: readw0w0itest=.false.

  real(8)::ebmx
  integer:: nbmx,mtet(3),ifq0p
  real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

  !     okumura
  integer::nwf,nsp_w,nqtt_w,iwf,jwf,inwf,iexc,kwf,lwf,ijwf,klwf,nnwf,ijwf_j
  integer::ifhamw,ifchipm_wan,ifgas,ifchipmz_wan,ifchipmr_wan,ifchipmrk_wan
  real(8),allocatable::ev_w1(:,:),ev_w2(:,:)
  complex(8),allocatable::evc_w1(:,:,:),evc_w2(:,:,:)
  logical(8):: ijklmag !for checkorb2
  complex(8),allocatable::kmat(:,:,:,:),wanmat(:,:)  ,wkmat(:,:),wkmat2(:,:),rmat(:,:,:),rmat2(:,:),&
       swkwmat(:,:),swkwmat2(:,:) ,sqemat(:,:),plmat(:,:),plpmat(:,:),plpmat_inv(:,:),wmat_check(:,:),kmat_check(:,:,:,:)
  complex(8)::trmat,trmatt
  real(8),allocatable::trmat2(:)
  real(8),allocatable::imat(:,:) !unit matrix for 1-WK
  real(8):: wan_ecore(1)
  logical::npmtwo,diag=.false.,t2g
  complex(8)::wan_i,wan_j,wan_k,wan_l,tmpwan,wanijkl
  integer::igv!,bsimple=3
  real(8):: qsh1(3),qsh2(3),znorm,rnqbz
  !For electron gas
  real(8),allocatable::qgsq1(:,:),qgsq2(:,:)
  real(8)::ntot_r,www,eta,deta,eta2
  !For screening W
  complex(8),allocatable:: scrw(:,:) ,cmat2(:,:)
  complex(8),allocatable:: eval_wk(:),vr_w(:,:),eval_k(:),eval_swkw(:)
  complex(8),allocatable:: eval_sqw(:)
  integer:: it,itp,isdummy,lorb
!!! q on symline
  real(8)::qrot(3),cr=2d-5
  logical(8)::init2=.true.,llsym=.true.,gskip
  logical(8)::d100,d110,d111,d1xx,dhpb,dxwf,dhnb,dnpb
  real(8)::rlatp(3,3),xmx2(3),qqin(3),qshort(3),qshort2,ppin(3),qlength,rs
  integer:: nlatout(3,48),nout,iout
  integer,parameter:: noutmx=48
  logical:: initiq,cmdopt0
  integer:: ifz,ifi,ifif
  real(8):: ef
  integer:: comm
  include "mpif.h"
! Pay attension to the following bootstrap sequence to fill data to modules!  
  debug = cmdopt0('--debug')
  comm = MPI_COMM_WORLD
  call m_lgunit_init() 
  hartree  = 2d0*rydberg()
  pi       = 4d0*datan(1d0)
  fourpi   = 4d0*pi
  sqfourpi = sqrt(fourpi)
  call MPI__Initialize()
  call MPI__consoleout('hhomogas')
  call cputid (0)
  call cputid(0)
  call read_BZDATA()
  call readhamindex0()
  call Genallcf_v3(incwfx=0) !Basic data. incwfin= 0 takes 'ForX0 for core' in GWinput
  call ReadGWinputKeys() !Readin dataset in GWinput   !      call Readq0p()    !Readin Offset Gamma
  tpioa=2d0*pi/alat
  iqxend = nqibz !+ nq0i
  do iq=1,nqibz
    iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
    write(6,"(' iq qibz nstibz=',2i5,3f9.4,i5)")iq,iqbz,qibz(:,iq) !,nstibz(iq)
 enddo
 
!!!!!!!!!!!! ! test for Homogenious gas corresponding to the density of Li 
  wemax=    3d0 !max value for plot
  omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call Getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw,ua,npmtwo=.false.)
  if(realomega .and. mpi__root) then
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
  niw=1 !no data along imag axis now
  noccxv = 1
  nmbas = 1
  print *,"noccxv ",noccxv
  noccx  = noccxv !+ nctot
  nprecx = ndble  !We use double precision arrays only.
  mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
  nspinmx = nspin
  iqxini=1
  iqxendx=iqxend
  ntot=1  ! We get ef=efermi for given ntot (number of electrons). We show rs parameter is shown in efermi_egas.
  ntot_r=ntot/2.0
  voltot = abs(alat**3*det33(plat)) !cell volume (bohr unit. We use Rydberg unit m=1/2 e^2=2 hbar=1)
  call efermi_egas(ntot_r,alat,plat,efz)
  ef=efz                    !ef for electron gas
  write(stdo,ftox)' alat=',ftof(alat)
  write(stdo,ftox)' plat(*,1)=',ftof(plat(:,1))
  write(stdo,ftox)' plat(*,2)=',ftof(plat(:,2))
  write(stdo,ftox)' plat(*,3)=',ftof(plat(:,3))
  write(stdo,ftox)" voltot ntot ef npm,ngc::",ftof(voltot),ntot,ftof(ef),npm,ngc
  call setefermi(ef) !set Fermi energy
  !! ======== Loop over iq ================================
  initiq=.true.
  iqloop: do 1001 iq = 1,10 !2,iqxend !iqxini,iqxend ! NOTE: q=(0,0,0) is omitted when iqxini=2
    print *
    print *
    q = [0.001d0,0d0,0d0]*iq !qibz(:,iq)          !qibze ! you can spefify any q, which you like.
    write(6,"('===== do 1001: iq q=',i7,3f9.4,' ========')")iq,q !qq
    !! get q in 1st BZ =shortest q. nout is the number of shortest (the same size of) data.
    !! When q is at the BZ boundary, nout>1 can be.
    ppin = matmul(transpose(plat),q)
    call shortn3_qlat(ppin) 
    iout=1
    qlength = tpioa * (sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2))**.5 !one of shortest q dueto the periodicicty of BZ
    do iout=1,nout
      write(*,"(a,3i5,f10.4,3f8.4)")'rrrr: nlat qlength q =',nlatout(:,iout),&
           (sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2))**.5, matmul(qlat(:,:),ppin+nlatout(:,iout))
    enddo
!!! |q+G|
    ngc=1 !only single MPB now
    if(.not.allocated(ngveccB)) allocate(ngveccB(3,1))
    ngveccB=0d0
    nwf=1
    nmbas=1 !dummy. matrix dimension.
    nmbas1=nmbas
    nmbas2=nmbas
    !! rcxq: imaginary part after x0kf_v4h and symmetrization. 
    !! zxq and zxqi are the main output after Hilbert transformation
    allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
    allocate( zxq (nmbas1,nmbas2,nw_i:nw), zxqi (nmbas1,nmbas2,niw))
    zxq=0d0;  zxqi=0d0;  rcxq = 0d0
    if(debug) write(6,*)' niw nw=',niw,nw
    !! ==== spin chi_charge or chi_+- ====
    is=1
    isf=1 !2 for magnetic
!!! qgsq: |q+G|^2
    allocate(ev_w1(nwf,nqbz),ev_w2(nwf,nqbz),source=0d0)
    allocate(qgsq1(ngc,nqbz),qgsq2(ngc,nqbz),source=0d0)
    write(stdo,ftox) 'bands: nwf=',nwf
    do kx=1,nqbz
      call read_qgband(alat,plat,  qbz(:,kx),ngc,ngveccB,is,qgsq1(1:ngc,kx),qsh1)
      call read_qgband(alat,plat,q+qbz(:,kx),ngc,ngveccB,isf,qgsq2(1:ngc,kx),qsh2)
      !c                   write(6,"('qqq111=',3f9.4,2x,3f9.4)") qbz(:,kx), qsh1
      !c                   write(6,"('qqq222=',3f9.4,2x,3f9.4)")q+qbz(:,kx),qsh2
!!! replace eigenvalue , egasmode nwf=1
      ev_w1(1:nwf,kx) = qgsq1(:,kx) !Rydberg
      ev_w2(1:nwf,kx) = qgsq2(:,kx) !Rydberg
      if(debug) then
      if(ev_w1(1,kx)<ef.and.ev_w2(1,kx)-ev_w1(1,kx)>0) then
        write(stdo,ftox)'qbz=',ftof(qbz(:,kx)),'ev1 (eV)=',ftof((ev_w1(1,kx)-ef)*rydberg()),&
             'de=',ftof((ev_w2(1,kx)-ev_w1(1,kx))*rydberg())
      endif
      endif
    enddo
    deallocate(qgsq1,qgsq2)
    !! Get tetrahedron weight. ! See m_tetwt. The tetrahedron weight is stored in a complicated manner.
    !! No core, thus nctot=0, wan_ecore is dummy.
    wan_ecore=0d0
    is=1
    isf=1
    write(stdo,ftox) ' === goto gettetwt iq =',iq,'ef(eV)=',ftof(ef*rydberg())
    call gettetwt(q,iq,is,isf,ev_w1,ev_w2,nwf) !Tetrahedron weight whh
    write(stdo,ftox) " === end gettetwt. we now have the tetrahedron weight whw"
    jpmloop: do 2011 jpm=1,npm     !! jpm=2: negative frequency
      kloop: do 2012 kx=1,nqbz     !! discrete k-point loop
        ibibloop: do 2013 ibib=1,nbnb(kx,jpm) !! n,n' band loop
!!!     n1b(ibib,k,jpm) = n :band index for k (occupied)
!!!     n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
          if(debug) then
            write(6,*) 'jpm,ibib',jpm,ibib
            write(6,*) 'kx,nbnb',kx,nbnb(kx,jpm)
            if (jpm==1) then
              write(6,*) 'n1b(occ),n2b(unocc):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
            else
              write(6,*) 'n1b(unocc),n2b(occ):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
            endif
          endif
          
!          i1=n1b(ibib,kx,jpm) 
!          i2=n2b(ibib,kx,jpm)
!          imgweight=merge(1d0,0d0,ev_w1(i1)>ev_w2(i2))
              
!!!   print *,"ihw,nhw",ihw(ibib,kx,jpm),nhw(ibib,kx,jpm)
          do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
            imagweight = whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm)) ! imagweight is the tetrahedron weight for (k,n1b), (q+k,n2b) .
            rcxq(1:nmbas1,1:nmbas2,iw,jpm) = rcxq(1:nmbas1,1:nmbas2,iw,jpm) + imagweight !in a.u.
          enddo
2013    enddo ibibloop
2012  enddo kloop
2011 enddo jpmloop
    deallocate(ev_w1,ev_w2)
    rcxq=rcxq*2 ! Spin factor  !    write(6,"('rcxq/nbnb 4:',3E13.5)") sum(abs(rcxq(:,:,:,:)))
    call tetdeallocate() !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
    write(6,"('  nmbas1 nmbas2 npm=',3i8)") nmbas1,nmbas2,npm
    schi=1                   ! flip over maj/min spins.
    ! Hilbert transformation.  We get zxq (complex(8), along real axis), and zxqi (complex(8), along imag axis).
    call dpsion5( realomega, imagomega, rcxq, nmbas1,nmbas2, zxq, zxqi,.false., schi,1,1d99,1d99) 
    if(initiq) then
      initiq=.false.
      open(newunit=ifz,file="x0homo.dat")
   endif
   !NOTE: zxq = chi0 \times voltot, where chi0 is defined in Fetter-Walecka.
   !      The demominator in zxq is in a.u. Thus zxq/2.0 has the dimension 1/Ry.
    do iw=nw_i,nw
      write(ifz, "(2i4,3f9.4,x,f10.5,E13.5,x,2E13.5,x,2E13.5)") iq,iw,q,qlength,2*freq_r(iw), zxq(1,1,iw) !vol/energy(a.u.) 
    enddo
    write(ifz,*)
    write(ifz,*)
    deallocate(rcxq,zxq,zxqi)
1001 enddo iqloop !continue !q point loop
  close(ifz)
  call cputid(0)
  call rx0( ' OK! homogas test')
  contains
  pure real(8) function det33(am)
    implicit none
    real(8),intent(in) :: am(3,3)
    det33= am(1,1)*am(2,2)*am(3,3) &
         -am(1,1)*am(3,2)*am(2,3) &
         -am(2,1)*am(1,2)*am(3,3) &
         +am(2,1)*am(3,2)*am(1,3) &
         +am(3,1)*am(1,2)*am(2,3) &
         -am(3,1)*am(2,2)*am(1,3)
  END function det33
end subroutine hhomogas
  
  ! c$$$!!  electron gas bare exchange (exact)
  ! c$$$        if (legas.and.exchange) then
  ! c$$$          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
  ! c$$$          pi         = 4.d0*datan(1.d0)
  ! c$$$          tpia       = 2.d0*pi/alat
  ! c$$$          qfermi= dsqrt(efz)
  ! c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
  ! c$$$          write (6,*)' --- exact electron gas bare exchange --- '
  ! c$$$          write (6,*)' density parameter rs= ', alpha/qfermi
  ! c$$$          write (6,*)' kf= ',qfermi
  ! c$$$          do      ip = 1,nq
  ! c$$$            qreal =  tpia*q(1:3,ip)
  ! c$$$            qm    = dsqrt ( sum(qreal**2) )
  ! c$$$            xsex  = hartree * egex (qm,efz)
  ! c$$$            write (6,*)
  ! c$$$            write (6,"(' True qm-ef Sx=',2f14.6,' q/qf=',f14.6)")
  ! c$$$     &       rydberg()*(qm**2-efz), xsex, qm/qfermi
  ! c$$$            write (6,"(' Num  qm-ef Sx=',2f14.6)") 
  ! c$$$     &       eqx(1,ip,is),        hartree*dreal(zsec(1,1,ip)) !sf 21May02
  ! c$$$            write (6,"(' === diff     =',2f14.6)") 
  ! c$$$     &       rydberg()*(qm**2-efz)-eqx(1,ip,is)
  ! c$$$     &       , xsex - hartree*dreal(zsec(1,1,ip)) !sf 21May02
  ! c$$$            write (661,"(' qm True qm-ef Sx=',3f14.6)") 
  ! c$$$     &       qm,rydberg()*(qm**2-efz), xsex
  ! c$$$            write (662,"(' qm Num  qm-ef Sx=',3f14.6)") 
  ! c$$$     &       qm,eqx(1,ip,is),     hartree*dreal(zsec(1,1,ip)) !sf 21May02
  ! c$$$ccc   write (ifsex(is),6600) qreal(1),qreal(2),qreal(3),xsex
  ! c$$$ccc   write (6,6600) qreal(1),qreal(2),qreal(3),xsex
  ! c$$$ccc   6600   format (' qreal =',3f8.4,'   SEx(q) =',d13.5)
  ! c$$$            write (663,"(2f14.6)") qm/qfermi, qfermi
  ! c$$$          end do
  ! c$$$        endif
  
  ! c$$$!! ef is taken as rs for the empty-sphere test case of legas=T case 
  ! c$$$!! HOMOGENIOUS GAS code. Usually not used. Need fixing if necessary.
  ! c$$$!! Keep this just as a memo.
  ! c$$$      legas = .false.
  ! c$$$      if(.false.) then
  ! c$$$        INQUIRE (FILE = 'LEGAS', EXIST = legas)
  ! c$$$        if(legas) then          !!! test for electron gas case.
  ! c$$$          write(6,*)' find LEGAS. legas =',legas
  ! c$$$          iflegas = 2101
  ! c$$$          open (iflegas,file='LEGAS')
  ! c$$$          read(iflegas,*)rs
  ! c$$$          close(iflegas)
  ! c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
  ! c$$$          qfermi = alpha/rs
  ! c$$$          efx  = qfermi**2
  ! c$$$          valn = efx**1.5d0*voltot/3d0/pi**2
  ! c$$$          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
  ! c$$$          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
  ! c$$$          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
  ! c$$$          if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
  ! c$$$        endif
  ! c$$$      endif
