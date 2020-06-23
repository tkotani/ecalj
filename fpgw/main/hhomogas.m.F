!! hhomogas, originally included in hmagnon.m.F by Okumura 
      program hhomogas
c      use m_ReadEfermi,only: readefermi,ef_read
c      use m_readqg,only: readqg,readngmx
c      use m_readeigen,only: readeval,init_readeigen,init_readeigen2
c      use m_hamindex,only:qtt,nqtt

      use m_read_bzdata,only: read_bzdata,
     & ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,
     & dq_,qbz,wbz,qibz,wibz,qbzw,
     & idtetf,ib1bz,idteti,
     & nstar,irk,nstbz

c      use m_genallcf_v3,only: genallcf_v3,
c     & alat,plat,ngrp,symgg

c     & nclass,natom,nspin,nl,nn, ngrp,
c     & nlmto,nlnmx, nctot,niw, !nw_input=>nw,
c     & alat, delta,deltaw,esmr,symgrp,clabl,iclass, !diw,dw,
c     & invg, il,in,im,nlnm, 
c     & plat, pos,ecore, symgg 
c      use m_keyvalue,only: getkeyvalue
c      use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
c      use m_readqgcou,only: readqgcou
      use m_mpi,only: MPI__hx0fp0_rankdivider2,MPI__task,MPI__Initialize,MPI__Finalize,MPI__root,
     & MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size,
     & MPI__ranktab,MPI__consoleout,MPI__barrier
c!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
c      use m_rdpp,only: rdpp,    !NOTE: "call rdpp" generate following data.
c     & nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr
c!! Generate matrix element for "call get_zmelt".
      use m_zmel,only:       !NOTE: these data set are stored in this module, and used  
     & nband
c    & ,itq,ngcmx,ngpmx, ppovlz,
c     & ppbir,shtvg, miat,tiat , ntq
!! frequency
      use m_freq,only: getfreq, !NOTE: call getfreq generate following data.
     & frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw !, frhis0,nwhis0 !output of getfreq
!! tetwt
      use m_tetwt,only: tetdeallocate,gettetwt, !followings are output of 'L871:call gettetwt')
     &  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
c!! w0 and w0i (head part at Gamma point)
c      use m_w0w0i,only: w0w0i,
c     & w0,w0i
c
c      use m_lldata,only: ll
      use m_homoelectron,only: read_qgband, efermi_egas !for gsq 
      implicit none

      integer::nctot=0,nspin=1,niw,ngrp
      real(8):: alat,plat(3,3),det33
      real(8),allocatable:: symgg(:,:,:)
!! ------------------------------------------------
!! We calculate chi0 by the follwoing three steps.
!!  gettetwt: tetrahedron weights
!!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
!!  dpsion5: calculate real part by the Hilbert transformation from the Im part
!!  eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
!!-------------------------------------------------

cccccc this may be wrong or correct cccccccccc
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
Cr      The indexing starts with core first and then valence on top of core
Cr      So n-index in "in" for valence electron is different from "inv".
cccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8):: q(3),  qgbin(3),qx(3)
      real(8):: ua=1d0 ! this is a dummy.
      integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
      integer:: ndble=8
      integer:: nword
      real(8),allocatable:: vxcfp(:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:)
      integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
c     &           , ngvecp(:,:), ngvecc(:,:), !,ngpn(:),ngcni(:),iqib(:),
c     &   ifppb(:)   !ongveccBr(:,:,:),nx(:,:),nblocha(:),lx(:),
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &  zw(:,:),zw0(:,:),
     &  zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: eqt(:), !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),
c     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
c     &  ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), 
     &  qbze(:,:),qibze(:,:)  !,ecore(:,:)  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:)
c     & ,rcxqmean(:,:,:,:),rcxqmeanc(:,:,:,:) !now rcxqmean is treated as a case of rcxq(nmbas,nmbas)

!  tetrahedron method
      logical :: tetra=.true. !,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
      ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer,allocatable :: 
     &     noccxvv(:),n2bminimum(:,:,:)
!     &     n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),
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

c      logical::imagonly=.false.,realonly=.false. !,readgwinput
      integer::iopen,maxocc2,iclose,
     & ixc,iqxini,iqxend,iqxendx,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,!nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx !ngcmx,
     &   ,nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, !ifvcfpout,,mdimx,nbloch
     &   ifepscond,nxx !,ifvxcpout,ifgb0vec
     &   ,nw0,iw,ifinin,iw0,ifwwk,noccxv,noccx
     &   ,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0dummy
     &   ,kx,isf,kqxx,kp,job,noccxvx(2)=-9999,nwmax  !,ifev1,ifev2 nbnbx,nhwtot,
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j
c     &   ,ngpmx !,  ifchipmlog

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      integer:: incwfin,  verbose

      integer:: ngc,mrecg !bzcase, 
      real(8):: quu(3), deltaq(3)!,qq(3) !,qqq(3)=0d0
c      logical:: omitqbz=.false., noq0p

      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)

      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg !if true, we use off-gamma mesh.
      integer:: nbcut,nbcut2

      integer,allocatable:: nstibz(:) !Nov2004 Miyake's tote
      real(8),allocatable:: ecqw(:,:) !,wiw(:)
      real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree
     &   ,pi,efz,qfermi,alpha,voltot,ecelgas,efx,valn
      integer:: iqbz,iqindx,iflegas,nmx
     &   ,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1 
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
      complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:),ddd(:,:)
     & ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
      integer::lwork,info,imin,ifzxq
      complex(8)::x0mx
      complex(8),allocatable:: UU0(:,:),VT0(:,:)

      logical ::  chipm=.false.,nolfco=.false. !sergeyv only ngczero=.false.,
     &  ,epsmode=.false.,normalm=.false., crpa=.false. 
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
      complex(8),allocatable:: am1(:),am2(:),mmat(:,:),
     &     x0mat(:,:),x0matinv(:,:),eiqrm(:)
      integer:: ifchipmn_mat, ifchipm_fmat !,ifchipm_mat
      integer::ifstoner,ifx,i1
      real(8):: Istoner,zz1,zz2,zz3,zz4,Istoner0,jzero2,dumm1,dumm2
      complex(8):: trr,trr0,trr1     , zzzx(4,4), zzzy(4,4),trrx,mmatx(4,4),denom(4,4)
      real(8),allocatable:: eee(:),mmnorm(:),
     &    asvec(:,:),ssv(:,:),sproj(:,:),sprojx(:,:), momsite(:)
      real(8):: eex(4),eey(4),qvv(3)
!!
c      logical :: newaniso,newaniso2,newanisox !,z1offd
      integer :: ngb0,ifvcoud,idummy,ifepstinv,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x,ifiss,iq0x
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:)
      real(8),allocatable:: vcousq(:)
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq

!! Eq.(40) in PRB81 125102
c      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:),sks(:),skI(:),sksI(:), 
     &  w_k(:),w_ks(:),w_kI(:), w_ksI(:), s_vc(:),vw_k(:),vw_ks(:)
      complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
      integer:: lxklm,nlxklm,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
      complex(8):: vc1vc2
      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)

      real(8),allocatable:: aik(:,:,:,:)
      integer,allocatable:: aiktimer(:,:)
      integer:: l2nl
      logical:: eibz4x0,tiii,iprintx,symmetrize,eibzmode
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
c      logical::wan=.true.,egasmode
      logical(8):: ijklmag !for checkorb2
      complex(8),allocatable::kmat(:,:,:,:),wanmat(:,:)
     &     ,wkmat(:,:),wkmat2(:,:),rmat(:,:,:),rmat2(:,:)
     &     ,swkwmat(:,:),swkwmat2(:,:)
     &     ,sqemat(:,:),plmat(:,:),plpmat(:,:),plpmat_inv(:,:)
     &     ,wmat_check(:,:),kmat_check(:,:,:,:)
      complex(8)::trmat,trmatt
      real(8),allocatable::trmat2(:)
      real(8),allocatable::imat(:,:) !unit matrix for 1-WK
!      real,allocatable::kmat(:,:,:,:),wanmat(:,:)
      real(8):: wan_ecore(1)
      logical::npmtwo,diag=.false.,t2g
      complex(8)::wan_i,wan_j,wan_k,wan_l,tmpwan,wanijkl
      integer::igv!,bsimple=3
      real(8)::qlat(3,3),qsh1(3),qsh2(3),znorm,rnqbz
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

      logical:: initiq
      integer:: ifz,ifi
      real(8):: ef
!! -------------------------------------------------------------------
      hartree  = 2d0*rydberg()
      pi       = 4d0*datan(1d0)
      fourpi   = 4d0*pi
      sqfourpi = sqrt(fourpi)
      call MPI__Initialize()
      call MPI__consoleout('hhomogas')
      call cputid (0)
      call cputid(0)

!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
!! As described at the top of this routine (see  use m_read_bzdata,only: read_bzdata,...),
!! we have many data set prepared after we finish read_BZDATA (all data successive to read_bzdata in the use statement).
      call read_BZDATA()

!! EFERMI is given by  call efermi_egas(ntot_r,alat,plat,efz) afterwards.
cc need for gettetwt; automatically read file(okumura) OK
c      call readefermi()


!! === Readin by genallcf. Set basic data for crystal
!! See "use m_genallcf_v3" at the begining of this routine
!!
c      incwfin=0  !use ForX0 for core in GWIN
c      call genallcf_v3(incwfin) !in module m_genallcf_v3
c      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
c      if(ngrp/=ngrp2 ) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')

      ifi = ifile_handle()
      open (ifi, file='SYMOPS')
      read(ifi,*) ngrp
      allocate(symgg(3,3,ngrp))
      do ig = 1,ngrp
        read(ifi,*)
        do i=1,3
          read(ifi,*) symgg(i,1:3,ig)
        enddo
c        print *,'ssssssss',ig,symgg(:,:,ig)
      enddo
      close(ifi)
      ifi = ifile_handle()
      open(ifi,file='LATTC')
      read(ifi,*) alat
      read(ifi,*) plat(1:3,1)
      read(ifi,*) plat(1:3,2)
      read(ifi,*) plat(1:3,3)
      close(ifi)

      tpioa=2d0*pi/alat
      iqxend = nqibz !+ nq0i
c      ginv = transpose(plat)   !ginv is inverse of plat.
      call minv33tp(plat,qlat) !plat -> qlat (confusing. qlat is reciprocal vector). 
      do iq=1,nqibz
         iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
         write(6,"(' iq qibz nstibz=',2i5,3f9.4,i5)")iq,iqbz,qibz(:,iq) !,nstibz(iq)
      enddo
      call shortn3_initialize(qlat,rlatp,xmx2)


!! We get 
!!   frhis (histgram bins along real axis), 
!!   freq_r (data point along real axis),
!!   freq_i (data point along img axis)
!!   nwhis,nw,npm,wiw  
!!   by getfreq
c      call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
c      omg2max = (Emax-Emin)*.5d0+.2d0 

ccccccccccccccccccccccccccccccccccc
c test for Li
      wemax=    3d0 !max value for plot
      omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
cccccccccccccccccccccccccccccccccccc


!! npmtwo=F means only positive energy for real axis.
      call getfreq(epsmode,realomega,imagomega,omg2max,wemax,niw,ua,MPI__root, npmtwo=.false.) !tetra,

      if(MPI__root) write(6,"(' nw npm=',2i5)") nw,npm
c      nwp = nw+1
c      if(.not.imagomega) niw=1
      niw=1 !no data along imag axis.

      noccxv = 1
      nmbas = 1
      print *,"noccxv ",noccxv
      noccx  = noccxv !+ nctot
      nprecx = ndble  !We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
      nspinmx = nspin
      iqxini=1
      mtet=(/1,1,1/) !dummy
      eibzmode=.false.          !simple symmetry
      iqxendx=iqxend
      allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &     ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      nwgt=1

!! Calculate x0(q,iw) and W == main loop 1001 for iq. 
 !! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
 !! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
 !! (or need to modify do 1001 loop).
 !! iq>nqibz for ixc=11 is not time-consuming (right???)
      call MPI__hx0fp0_rankdivider2(iqxini,iqxend)

      print *,' alat=',alat
      print *,' plat(*,1)=',plat(:,1)
      print *,' plat(*,2)=',plat(:,2)
      print *,' plat(*,3)=',plat(:,3)
      print *," voltot ntot ef npm,ngc::",voltot,ntot,ef,npm,ngc

!! efermi for given ntot. The rs parameter is shown in efermi_egas.
      ntot=1
      ntot_r=ntot/2.0
      voltot = abs(alat**3*det33(plat))
      call efermi_egas(ntot_r,alat,plat,efz)
      ef=efz                    !ef for electron gas


!! ======== Loop over iq ================================
      initiq=.true.
      do 1001 iq = iqxini,iqxend ! NOTE: q=(0,0,0) is omitted when iqxini=2
        print *
        print *
        print *
        q = qibz(:,iq)          !qibze ! you can spefify any q, which you like.
        write(6,"('===== do 1001: iq q=',i7,3f9.4,' ========')")iq,q !qq
         !!! symmetry check

!! get q in 1st BZ =shortest q. nout is the number of shortest (the same size of) data.
!! When q is at the BZ boundary, nout>1 can be.
        ppin = matmul(transpose(plat),q) 
        call shortn3(rlatp,xmx2,ppin, noutmx, nout,nlatout)
        iout=1
        qlength = tpioa * (sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2))**.5
        do iout=1,nout
          write(*,"(a,3i5,f10.4,3f8.4)")'rrrr: nlat qlength q =',nlatout(:,iout),
     &     (sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2))**.5,
     &     matmul(qlat(:,:),ppin+nlatout(:,iout))
         enddo

!!! |q+G|
        ngc=1
        if(.not.allocated(ngveccB)) allocate(ngveccB(3,1))
        ngveccB=0d0
        nwf=1
        nmbas=1 !dummy. matrix dimension.
        nmbas1=nmbas
        nmbas2=nmbas

!! rcxq: imaginary part after x0kf_v4h and symmetrization. 
!! zxq and zxqi are the main output after Hilbert transformation
        if(allocated(zxq) )  deallocate(zxq)
        if(allocated(zxqi) ) deallocate(zxqi)
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
        allocate( zxq (nmbas1,nmbas2,nw_i:nw), zxqi (nmbas1,nmbas2,niw))
        zxq=0d0;  zxqi=0d0;  rcxq = 0d0
        if(debug) write(6,*)' niw nw=',niw,nw
          
!! ==== spin chi_charge or chi_+- ====
        is=1
        isf=1 !2 for magnetic
!! Tetrahedron weight.
!! output
!!     nbnbx
!!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
!!     nhw(ibjb,kx): the number of histogram sections
!!     jhw(ibjb,kx): pointer to whw
!!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
!!     : histogram weights for given ib,jb,kx for histogram sections
!!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c            write(6,*) ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2
!! input
!!     ekxx1 for   rk,is
!!     ekxx2 for q+rk,isf 

c          !! get Efermi, qfermi, rs of electron gas
c          if (egasmode) then
!!
c        print *,' rs=',rrs
!!! qgsq: |q+G|^2
        if(.not. allocated(ev_w1)) allocate(ev_w1(nwf,nqbz))
        if(.not. allocated(ev_w2)) allocate(ev_w2(nwf,nqbz))
        allocate(qgsq1(ngc,nqbz),qgsq2(ngc,nqbz))
        qgsq1=0d0;qgsq2=0d0
        ev_w1=0d0;ev_w2=0d0

c$$$                ifgas=ifile_handle()
c$$$                if (is==1 .and. iq==54) open(ifgas,file="band_electrongas.dat")                
        do kx=1,nqbz
          call read_qgband(alat,plat,  qbz(:,kx),ngc,ngveccB,is,qgsq1(1:ngc,kx),qsh1)
          call read_qgband(alat,plat,q+qbz(:,kx),ngc,ngveccB,isf,qgsq2(1:ngc,kx),qsh2)
c                   write(6,"('qqq111=',3f9.4,2x,3f9.4)") qbz(:,kx), qsh1
c                   write(6,"('qqq222=',3f9.4,2x,3f9.4)")q+qbz(:,kx),qsh2
!!! replace eigenvalue , egasmode nwf=1
          ev_w1(1:nwf,kx) = qgsq1(:,kx) !Rydberg
          ev_w2(1:nwf,kx) = qgsq2(:,kx) !Rydberg
        enddo
        if (allocated(qgsq1)) deallocate(qgsq1)
        if (allocated(qgsq2)) deallocate(qgsq2)

!! Get tetrahedron weight.
!! See m_tetwt. The tetrahedron weight is stored in a complicated manner.
!! No core, thus nctot=0, wan_ecore is dummy.
        wan_ecore=0d0
        is=1
        isf=1
        call gettetwt(q,iq, is,isf,nwgt(:,iq),frhis,nwhis,npm,
     i            qbas,ginv, ef, nqibz, nwf,ev_w1,ev_w2, nctot,wan_ecore,
     i            nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
     i            nwf,ebmx,mtet,eibzmode) !nov2016
        write(6,*) "=== end gettetwt. we now have the tetrahedron weight whw"
        deallocate(ev_w1,ev_w2)

        rnqbz=1/real(nqbz)
        do 2011 jpm=1,npm     !! jpm=2: negative frequency
        do 2012 kx=1,nqbz     !! discrete k-point loop
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
!!!   print *,"ihw,nhw",ihw(ibib,kx,jpm),nhw(ibib,kx,jpm)
          do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
            imagweight = whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm)) ! imagweight is the tetrahedron weight for (k,n1b), (q+k,n2b) .
            rcxq(1:nmbas1,1:nmbas2,iw,jpm) = rcxq(1:nmbas1,1:nmbas2,iw,jpm) + imagweight
          enddo
 2013   continue            !ibib-loop
 2012   continue              !k-loop
 2011   continue                !jpm-loop

!! normalization check
       write(6,*) "rcxq/nbnb rnqbz:",rnqbz
       write(6,"('rcxq/nbnb 4:',3E13.5)") sum(abs(rcxq(:,:,:,:)))
       call tetdeallocate()     !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)


       write(6,"('  nmbas1 nmbas2 npm=',3i8)") nmbas1,nmbas2,npm
       schi=1                   ! flip over maj/min spins.

!! Hilbert transformation. zxq (complex(8), along real axis), and zxqi (complex(8), along imag axis).
       call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i         rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o         zxq, zxqi,
     i         .false., schi,1,1d99,1d99) !ecut(iecut),ecuts(iecut))


!! Write final results. zxq is just 1x1 matrix for homogenious gas.
       if(initiq) then
         ifz = ifile_handle()
         initiq=.false.
         open(ifz,file="x0homo.dat")
       endif
       do iw=nw_i,nw
         write(ifz, "(2i4,3f9.4,x,f10.5,E13.5,x,2E13.5)") iq,iw,q,qlength,freq_r(iw), zxq(1,1,iw)*hartree
       enddo
!!       write(ifz,*)
       write(ifz,*)

ccccccccccccccccccccccccccccccccccccccccccccc 20190818
       zxqi=0d0
       write(6,*) "AAAAA, nw_i,nw:",nw_i,nw
       do iw=nw_i,nw
          if (iw/=0) zxqi(1,1,:) = zxqi(1,1,:) + zxq(1,1,:)/freq_r(iw)*(freq_r(iw)-freq_r(iw-1))/2
       enddo
       write(6,'("sum(aimag(zxq(1,1,:)))=",3E13.5)') sum(aimag(zxq(1,1,:))),sum(aimag(zxqi(1,1,:))),real(zxq(1,1,0))
       !write(6,*) "sum(aimag(zxq(1,1,:)))=",sum(aimag(zxq(1,1,:))),real(zxq(1,1,0))
ccccccccccccccccccccccccccccccccccccccccccccc




       if(allocated(rcxq) ) deallocate(rcxq)
       if(allocated(zw0)) deallocate(zw0)
       if(allocated(zxq )) deallocate(zxq)
       if(allocated(zxqi)) deallocate(zxqi)
 1001 continue !q point loop
      close(ifz)
      call MPI__barrier()
      call cputid(0)
      call MPI__Finalize
      call rx0( ' OK! homogas')
      end



