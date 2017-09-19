!! cp hx0fp0.m.F in order to develop by Okumura (2017/03/23)
!!  Calculate x0, \epsilon, spin susceptibility. 
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_lmf_chipm mode works).
      program hmagnon
      use m_readwan,only: write_qdata, wan_readeigen, wan_readeval
!     &     evecc_w, eval_w  !! wannier eigen
      use m_ReadEfermi,only: readefermi,ef
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: readeval,init_readeigen,init_readeigen2
      use m_hamindex,only:qtt,nqtt
      use m_read_bzdata,only: read_bzdata,
     & ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,
     & dq_,qbz,wbz,qibz,wibz,qbzw,
     & idtetf,ib1bz,idteti,
     & nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     & nclass,natom,nspin,nl,nn, ngrp,
     & nlmto,nlnmx, nctot,niw, !nw_input=>nw,
     & alat, delta,deltaw,esmr,symgrp,clabl,iclass, !diw,dw,
     & invg, il,in,im,nlnm, 
     & plat, pos,ecore, symgg 
      use m_keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
      use m_readqgcou,only: readqgcou
      use m_mpi,only: MPI__hx0fp0_rankdivider2,MPI__task,MPI__Initialize,MPI__Finalize,MPI__root,
     & MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size,
     & MPI__ranktab,MPI__consoleout,MPI__barrier
!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: rdpp,    !NOTE: "call rdpp" generate following data.
     & nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr
!! Generate matrix element for "call get_zmelt".
      use m_zmel,only:       !NOTE: these data set are stored in this module, and used  
     & nband,itq,ngcmx,ngpmx, ppovlz,
     & ppbir,shtvg, miat,tiat , ntq
!! frequency
      use m_freq,only: getfreq, !NOTE: call getfreq generate following data.
     & frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw !, frhis0,nwhis0 !output of getfreq
!! tetwt
      use m_tetwt,only: tetdeallocate,gettetwt, !followings are output of 'L871:call gettetwt')
     &  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
!! w0 and w0i (head part at Gamma point)
      use m_w0w0i,only: w0w0i,
     & w0,w0i

      use m_lldata,only: ll

      use m_homoelectron,only: read_qgband, efermi_egas !for gsq 
      implicit none
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
      logical :: realomega=.true., imagomega=.true.
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
      logical:: omitqbz=.false., noq0p

      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)

      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg !if true, we use off-gamma mesh.
      integer:: nbcut,nbcut2

      integer,allocatable:: nstibz(:) !Nov2004 Miyake's tote
      real(8),allocatable:: ecqw(:,:) !,wiw(:)
      real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree
     &   ,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
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
      real(8):: qread(3),imagweight

      character(128):: vcoudfile,aaax
      integer:: src,dest
      logical:: lqall
      integer,allocatable :: iclasst(:), invgx(:)
      integer:: ificlass,ifile_handle,k
      complex(8),allocatable:: ppovl_(:,:)

      logical:: readw0w0itest=.false.

      real(8)::ebmx
      integer:: nbmx,mtet(3),ifq0p
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

!     okumura
      integer::nwf,nsp_w,nqtt_w,iwf,jwf,inwf,iexc
      integer::ifhamw,ifchipm_wan,ifgas,ifchipmz_wan
      real(8),allocatable::ev_w1(:,:),ev_w2(:,:)
      complex(8),allocatable::evc_w1(:,:,:),evc_w2(:,:,:)
      logical::wan=.true.,egasmode
      complex(8),allocatable::kmat(:,:,:,:),wanmat(:,:),wanij(:,:)
!      real,allocatable::kmat(:,:,:,:),wanmat(:,:)
      real(8):: wan_ecore(1)
      logical::npmtwo,diag=.false.,t2g
      complex(8)::wan_i,wan_j,wan_k,wan_l,tmpwan,wanijkl
      integer::igv,bsimple=3
      real(8)::qlat(3,3),qsh1(3),qsh2(3)
      complex(8):: wsum(5) !normalization check
      !For electron gas
      real(8),allocatable::qgsq1(:,:),qgsq2(:,:)
      real(8)::ntot_r,www
      
      integer:: it,itp,isdummy
!! -------------------------------------------------------------------
      call MPI__Initialize()
      call MPI__consoleout('hx0fp0')
      call cputid (0)
      hartree  = 2d0*rydberg()
      pi       = 4d0*datan(1d0)
      fourpi   = 4d0*pi
      sqfourpi = sqrt(fourpi)

c$$$!! computational mode select ! takao keeps only the Sergey mode.
c$$$      write(6,"(a)") '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'
c$$$      write(6,"(a)") ' #1:run mode'
c$$$      write(6,"(a)") '    11  : normal    '
c$$$      write(6,"(a)") '    111 : normal    fullband '
c$$$      write(6,"(a)") '    10111 : normal  crpa     '
c$$$      write(6,"(a)") '    202 : epsNoLFC  '
c$$$      write(6,"(a)") '    203 : eps       '
c$$$      write(6,"(a)") '    222 : chi^+- NoLFC'
c$$$      write(6,"(a)") '    223 : chi^+- wannier'
c$$$c      write(6,"(a)") '    12  : total energy Miyake: developint'
c$$$      write(6,"(a)")  '-------------------------------------------------------'
c$$$      if( MPI__root ) read(5,*) ixc
c$$$      call MPI__Broadcast(ixc)
      call cputid(0)

!! Set switches: ---
!!  normalm: normal eps mode
!!    crpa: crpa mode
!!  epsmode: (normalm or epsmode)
!!    omitqbz: qbz>nqbz+1 are calulated
!!  realomega:
!!   \chi on real axis
!!  imagomega:
!!   \chi on imag omega
!!  lqall:
!!   limited range of \chi on real axis (mainly for memory reduction
!!  chipm:
!!    \Chi_pm mode (nspin=2)
!!  nolfco:
!!    no local field correction
!!
      lqall=.true.
c$$$      if(ixc==11) then; write(6,*) " OK ixc=11 normal mode "
c$$$        normalm=.true.
c$$$        lqall=.false.
c$$$      elseif(ixc==111) then; write(6,*) " OK ixc=111 normal mode. fullband"
c$$$        normalm=.true.
c$$$      elseif(ixc==10011) then; write(6,*) " OK ixc=10011 crpa mode "
c$$$        normalm=.true.
c$$$        crpa=.true.
c$$$!! -- eps mode NoLFC
c$$$      elseif(ixc==202) then
c$$$        write(6,*) " OK ixc=202  sergey's eps mode Only NoLFC "
c$$$        imagomega=.false.
c$$$        omitqbz=.true.
c$$$        epsmode = .true.
c$$$        nolfco=.true.
c$$$!!    eps mode with LFC
c$$$      elseif(ixc==203) then
c$$$        write(6,*) " OK ixc=203 sergey's eps mode with LFC "
c$$$        imagomega=.false.
c$$$        omitqbz=.true.
c$$$        epsmode = .true.
c$$$!! chipm mode NoLFC
c$$$      elseif(ixc==222) then
c$$$        write(6,*) " OK ixc=222    chipm sergey's "
c$$$        imagomega =.false.
c$$$        omitqbz=.true.
c$$$        epsmode = .true.
c$$$        chipm=.true.
c$$$        nolfco=.true.
c$$$      elseif(ixc==223) then

      write(6,*) " OK ixc=223    chipm_wannier sergey's "
        imagomega =.false.
        omitqbz=.false.
        epsmode = .true.
        chipm=.true.
        nolfco=.true.
        wan=.true.
c$$$      else
c$$$        call rx( ' hx0fp0: given mode ixc is not appropriate')
c$$$      endif

!! We fix newaniso2=T now.

!! Prof.Naraga said " write(6,*)'Timereversal=',Timereversal()" 
!! here caused a stop in ifort ver.1x.x. Why? May be a compilar bug, and fixed now.

!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
      call read_BZDATA()

!! Read electron gas mode or not.
      call getkeyvalue("GWinput","EGASmode",egasmode,default=.false.)
      print *,"electron gas mode (egasmode): ",egasmode
      if (egasmode) bsimple=2
      
      
c$$$!! Use off-regular mesh (not including Gamma) for qbzreg()=F See hx0fp0.m.sc.F also.
c$$$!! This must be consistent with qg4gw.F-mkqg.F
c$$$      if(.not.qbzreg()) then 
c$$$         deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
c$$$         do i=1,nqbz
c$$$            qbz(:,i) = qbz(:,i) - deltaq/2d0
c$$$            write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
c$$$         enddo
c$$$      endif
      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp

c$$$!! check write
c$$$      if(MPI__root) then
c$$$         do i=1,nqbz
c$$$            if(i<10.or.i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
c$$$            if(i==10.and.nqbz>18) write(6,"('... ')")
c$$$         enddo
c$$$         write(6,*)' nqbz nqibz =',nqbz,nqibz
c$$$      endif

!! EFERMI
cc need for gettetwt; automatically read file(okumura) OK
      call readefermi()
      write(6,"(a,f12.6)")' --- READIN aa ef from EFERMI. ef=',ef



!! === Readin by genallcf. Set basic data for crystal
!! See "use m_genallcf_v3" at the begining of this routine
!!
      incwfin=0  !use ForX0 for core in GWIN
      call genallcf_v3(incwfin) !in module m_genallcf_v3
!     ! incwfin is for core handling
c     write(6,*)' ncore nw_input delta=', ncore,nw_input,delta
      if(ngrp/=ngrp2       ) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      if(chipm.and.nspin==1) call rx( 'chipm mode is for nspin=2')

c      debug=.false.
c      if(verbose()>=100) debug=.true.
c      if(debug) write(6,*)' end of genallc'

      tpioa=2d0*pi/alat


!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
c      if(nclass /= natom) call rx( ' nclass /= natom ')


c$$$!! --- Only tetrahedron method is allowed now.
c$$$      if(delta > 0d0) then
c$$$        call rx(' hx0fp0: only tetra=T is supported. delta in GWinput should be negative')
c$$$      else
c$$$        delta = -delta
c$$$        write(6,"(' hx0fp0: tetrahedron mode delta=',d13.6)") delta
c$$$      endif

c$$$!! --- Read dimensions of h,hb
c$$$      ifhbe  = iopen('hbe.d',1,0,0)
c$$$      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg !Note: nband is in m_zmel
c$$$      is = iclose('hbe.d')
c$$$      if(nlmto/=nlmtot) call rx('hx0fp0: nlmto/=nlmtot in hbe.d')

!     ! --- okumura Read dimensions of hamiltonian_wannier, spin, nqtt
c      if(ixc==223) then
      ifhamw=ifile_handle()
      open(ifhamw,file="wan4chi.d",form='unformatted')
      read(ifhamw) nwf,nsp_w,nqtt_w
      close(ifhamw)
      write(6,*) "nwf,nsp_w,nqtt",nwf,nsp_w,nqtt_w
c      endif

      
c$$$!! --- Readin Offset Gamma --------
c$$$      if(debug) write(6,*) 'reading QOP'
c$$$      ifq0p=ifile_handle()
c$$$      open (ifq0p,file='Q0P')
c$$$      read (ifq0p,"(i5)") nq0i
c$$$      write(6,*) ' ### nqibz nq0i=', nqibz,nq0i
c$$$      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
c$$$      do i=1,nq0i
c$$$        read (ifq0p, * ) wqt(i),q0i(1:3,i)
c$$$      enddo
c$$$      close(ifq0p)
c$$$      nq0ix = nq0i
c$$$      do i=1,nq0i
c$$$        if(wqt(i)==0d0 ) then
c$$$          nq0ix = i-1
c$$$          exit
c$$$        endif
c$$$      enddo
c$$$      neps = nq0i - nq0ix  ! number of zero weight q0p 
c$$$      write( 6,*) ' num of zero weight q0p=',neps
c$$$      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
c$$$
c$$$!! Readin q+G. nqbze and nqibze are for adding Q0P related points to nqbz and nqibz.
c$$$      call getsrdpp2( nclass,nl,nxx)
c$$$      call readngmx('QGpsi',ngpmx)
c$$$      call readngmx('QGcou',ngcmx)
c$$$

c      nq0i=0
c      write(6,"(' ngcmx ngpmx= ',2i8)") ngcmx,ngpmx 
c      nqbze  = nqbz *(1 + nq0i)
c      nqibze = nqibz + nq0i
c      allocate( qbze(3, nqbze)), qibze(3, nqibze))
c      qbze(:,1:nqbz)   = qbz(:,1:nqbz)
c      qibze(:,1:nqibz) = qibz(:,1:nqibz)
c$$$      do i = 1,nq0i
c$$$        qibze(:,nqibz+i)  = q0i(:,i)
c$$$        ini = nqbz*(1 + i -1)
c$$$        do ix=1,nqbz
c$$$          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
c$$$          if( abs(qbze(1,ini+ix)+0.1d0)+abs(qbze(2,ini+ix)+0.1d0)<1d-6 ) then
c$$$            write(6,"('hx0fp0 qbze q0i=',i8,3f18.14,2x,3f14.10)") ini+ix,qbze(:,ini+ix),q0i(:,i)
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$!! ngrpx=1 (no symmetry operation in hx0fp0), whereas we use ngrp in eibzmode=T.
c$$$      ngrpx = 1
c$$$      l2nl=2*(nl-1)
c$$$      allocate(symope(3,3))
c$$$      symope(1:3,1) = (/1d0,0d0,0d0/)
c$$$      symope(1:3,2) = (/0d0,1d0,0d0/)
c$$$      symope(1:3,3) = (/0d0,0d0,1d0/)
c$$$      ificlass=ifile_handle()
c$$$      open (ificlass,file='CLASS')
c$$$      allocate(iclasst(natom),invgx(ngrp)
c$$$     & ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
c$$$      write(6,*)'  --- Readingin CLASS info ---'
c$$$      do ibas = 1,natom
c$$$        read(ificlass,*) ibasx, iclasst(ibas)
c$$$        write(6, "(2i10)") ibasx, iclasst(ibas)
c$$$      enddo
c$$$      close(ificlass)
c$$$!! Get space-group transformation information. See header of mptaouof. 
c$$$!! Here we use ngrpx=1 ==> "no symmetry operation in hx0fp0", c.f. hsfp0.sc.m.F case.
c$$$      call mptauof(symope,ngrpx,plat,natom,pos,iclasst
c$$$     o ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
c$$$      if(verbose()>=40) write (*,*)' hsfp0.sc.m.F: end of mptauof'
c$$$
c$$$!!  ppbrd = radial integrals,  cgr = rotated cg coeffecients (no rotatio here since nrgpx=1 for identity matrix)
c$$$      call rdpp(nxx, nl, ngrpx, nn, nclass, nspin, symope,qbas) !see m_rdpp
c$$$      ntq = nband
c$$$      allocate(itq(ntq)) !itq=i for i=1,ntq. a dummy. c.f. hsfp0.sc.F
c$$$      do i=1,ntq
c$$$        itq(i)=i
c$$$      enddo
c$$$     
c$$$!! Pointer to optimal product basis
c$$$      nblochpmx = nbloch + ngcmx
c$$$      allocate(ngveccB(3,ngcmx)) ! work arry

      iqxend = nqibz !+ nq0i
c$$$      write(6,*) ' nqibz nqibze=',nqibz,nqibze

!! Initialization of readEigen !readin m_hamindex
      ginv = transpose(plat)
c      call init_readeigen(ginv,nspin,nband,mrece) !EVU EVD are read in init_readeigen
c      call init_readeigen2(mrecb,nlmto,mrecg)
      if(verbose()>50) print *,'eeee exit of init_readeigen2'







ccccccccccccccccccccccccccccccccccccc
!! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
c      call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
c      if (nctot > 0) Emin=minval(ecore(:,1:nspin))
c      omg2max = (Emax-Emin)*.5d0+.2d0 

      wemax=    5d0 !max value for plot
      omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
c      if(MPI__root) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
c
!! getfreq returun date given at " use m_freq,only:". 
c      if(.not.epsmode) call getwemax(lqall,wemax) !wemax is to determine nw !real axis divisions

!! NOTE: npmtwo=T sets npm=2
!! optional npmtwo is added aug2017
      call getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,MPI__root, npmtwo=.true.)
      if(MPI__root) write(6,"(' nw npm=',2i5)") nw,npm
c      nwp = nw+1
      if(.not.imagomega) niw=1

!! Get eigenvector corresponds to exp(iqr) (q is almost zero).
c      if(epsmode) then 
c        allocate(epsi(nw_i:nw,neps)) 
c      endif


cccccccccccccccccccccccccccccccccc
!! Tetra hedron initialization
      if (egasmode) then
         noccxv = 1
         nmbas = 1
      elseif (wan) then
         noccxv = nwf
         nmbas = nwf
      else
         call rx("noccxv can not be defined") 
      endif
      print *,"noccxv ",noccxv
c      noccxv = maxocc2 (nspin,ef, nband, qbze,nqbze) ! maximum no. occupied valence states
ccccccccccccccccccccccccccccccc

c      if(noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
ccccccccccc
      noccx  = noccxv + nctot
cccccccccc
      nprecx = ndble  !We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
c$$$      if (MPI__root) then
c$$$        ifwd   = iopen('WV.d',1,-1,0)
c$$$        write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
c$$$        ifwd = iclose('WV.d'); ifwd=0
c$$$      endif
c$$$      allocate(  zw(nblochpmx,nblochpmx) )
      nspinmx = nspin

c$$$!!... these are used x0k
c$$$      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
c$$$      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
c$$$      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2

c$$$!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,R(r))> 
c$$$!! This is general for rotated CG coefficient; but hx0fp0 mode is only for  ngrpx=1 (not rotated).
c$$$!! Compare usage in hsfp0 modes.
c$$$      irot=1
c$$$      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,irot,nspin))
c$$$      do is = 1,nspin
c$$$        call ppbafp_v2 (irot,ngrpx,is,nspin,
c$$$     i   il,in,im,nlnm,      
c$$$     i   nl,nn,nclass,nlnmx,
c$$$     i   mdimx,lx,nx,nxx,       !Bloch wave    
c$$$     i   cgr, nl-1,             !rotated CG
c$$$     i   ppbrd,                 !radial integrals
c$$$     o   ppbir(:,irot,is))      ! this is in m_zmel, used to generate <phi|phi B>
c$$$      enddo  
c$$$      if(debug)write(6,*) ' end of ppbafp_v2'

!! Set iqxini !omitqbz means skip loopf for iq=1,nqibz
c      if(omitqbz) then
c        iqxini= nqibz + 1
c      else
        iqxini= 1
c      endif

!!!!!
c$$$      necut=1
c$$$      allocate(ecut(necut),ecuts(necut))
c$$$      call getkeyvalue("GWinput","ecut_p" ,ecut, necut,default=(/1d10/) )
c$$$      call getkeyvalue("GWinput","ecuts_p",ecuts,necut,default=(/1d10/) )
!!
c$$$      if( chipm ) then
c$$$        nmbas=natom
c$$$        allocate(imbas(nmbas),imbas_s(nmbas))
c$$$        istat=-9999 ! istat=-9999 means noumber of readin arguments is returened in istat.
c$$$        call getkeyvalue("GWinput","MagAtom", imbas,nmbas,status=istat)
c$$$        nmbas = istat
c$$$        write(6,"('Readin MagAtom nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)
c$$$        imbas_s(1:nmbas) = imbas(1:nmbas)
c$$$        imbas(1:nmbas)   = abs(imbas(1:nmbas))
c$$$        allocate(jcoup(nw_i:nw,neps) )
c$$$        allocate( svec(nbloch,nmbas) )  !sep2006
c$$$        svec=0d0
c$$$        allocate( cvec(nbloch,nmbas),momsite(nmbas), mmnorm(nmbas)) !May2007
c$$$        cvec=0d0
c$$$        do imb=1,nmbas
c$$$          ibas= imbas(imb)
c$$$          ifv = iopen ('MixSpin.'//charnum3(ibas),1,3,0)
c$$$          read(ifv,*) ibasx,lxx
c$$$          allocate(nxx_r(0:lxx))
c$$$          do i=0,lxx
c$$$            read(ifv,*) nxx_r(i)   !   write(6,"(2i5,d13.6)") nxx_r(i)
c$$$          enddo
c$$$          allocate(spinvec((lxx+1)**2,maxval(nxx_r)))
c$$$          allocate(consvec((lxx+1)**2,maxval(nxx_r)))
c$$$          spinvec=0d0
c$$$          do ilmx = 1, (lxx+1)**2
c$$$            lb = ll(ilmx )         !  write(6,*)' lb=',lb,lxx,ilmx
c$$$            do ixx = 1, nxx_r(lb)  !  write(6,*)' nn=',nn,nxx_r(lb)
c$$$              if(ilmx==1) then
c$$$                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2 ,consvec(ilmx,ixx)
c$$$              else
c$$$                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),dumm1,dumm2 ,consvec(ilmx,ixx)
c$$$              endif
c$$$!           write(6,"(2i5,d13.6)") ilmx, ixx, spinvec(ilmx,ixx)
c$$$            enddo
c$$$          enddo
c$$$!! Calculate ChiPM. So sign of omega should be correct.
c$$$          if(imb==1) then !determine spin direction with respect to ibas=imbas(imb=1)
c$$$            spinmom=(chg1-chg2)
c$$$            schi=1d0
c$$$            if(spinmom<0d0) then
c$$$              schi  = -1d0   ! This affects to dpsion. Obtained results
c$$$                             ! should be the same in both mode.
c$$$            endif
c$$$          endif
c$$$!!  ReOrdering of spinvec in natom ordering...
c$$$          i=0
c$$$          if(ibas>1) i= sum(nblocha(1:ibas-1))
c$$$          do lb  = 0, lx (ibas)
c$$$            do nb  = 1, nx (lb,ibas)
c$$$              do mb  = -lb, lb
c$$$                i = i+1
c$$$                ilmx = lb**2+ lb+ mb +1
c$$$                svec(i,imb) = spinvec(ilmx,nb)
c$$$                cvec(i,imb) = consvec(ilmx,nb)
c$$$                write(6,"(' i lb mb svec svec**2=',3i4,2d13.5)") 
c$$$     &      i,lb,mb,svec(i,imb),svec(i,imb)**2
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          deallocate(nxx_r,spinvec,consvec)
c$$$          close(ifv)
c$$$          mmnorm (imb) = sqrt(sum(svec(:,imb)**2))
c$$$          momsite(imb) = chg1-chg2
c$$$c        write(6,"(' svecsum=',e23.15)") sum(svec(:,imb)**2
c$$$c        write(ifchipmlog,"(2e23.15,' ! mmom mmnorm')")momsite(imb),mmnorm(imb)
c$$$          write(6,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
c$$$        enddo
c$$$      endif

c$$$!! nov2016 moved from tetwt5 --> here
c$$$      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband )
c$$$      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
      mtet=(/1,1,1/) !dummy
c$$$      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
      ! multitet=T ==> micro tetrahedron method (divided-tetrahedron). Not used so much now...
cccccccccccccccccccccccccccccccccccccccc
!     allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))
cccccccccccccccccccccccccccccccccccccccc


!okumura
      
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
 !! EIBZ mode memo for nolfco (right?) 
 !! If eibzmode=T, it is efficient but can slightly break crystal symmetry.(how much?)
 !! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
 !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
c     okumura
c      eibzmode=eibz4x0()
      eibzmode=.false. !simple symmetry
c$$$c     end okumura
c$$$      if(eibzmode) then
c$$$ !! commentout block inversion Use iqxendx=iqxend because of full inversion
c$$$        call cputid(0)
c$$$        write(6,*)' ---goto eibzmode block ---'
c$$$        iqxendx=iqxend
c$$$        if(epsmode) iqxendx=iqxend
c$$$        allocate( nwgt(nqbz,iqxini:iqxendx), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
c$$$     &   igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
c$$$     &   eibzsym(ngrp,-1:1,iqxini:iqxendx))
c$$$ !! Check timereversal is required for symmetrization operation or not. 
c$$$ !! If tiii=timereversal=F is enforced,
c$$$ !!   the symmetrization procedure in x0kf_v4h becomes a little time-consuming.
c$$$        write(6,*)
c$$$        write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal() 
c$$$        if(MPI__root) iprintx=.true.
c$$$        call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,
c$$$     i   timereversal(),ginv,iprintx,
c$$$     o   nwgt,igx,igxt,eibzsym,tiii)
c$$$        write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
c$$$        call cputid(0)
c$$$ !! PBindex: index for product basis. 
c$$$ !! We will unify this system; still similar is used in ppbafp_v2.
c$$$        call PBindex(natom,lx,l2nl,nx)!all input. Returns requied index stored in arrays in m_pbindex.
c$$$        call cputid(0)
c$$$        call readqgcou() !no input. Read QGcou and store date into variables.
c$$$      else     !dummy allocation to overlaid -check bound !sep2014

        iqxendx=iqxend
        allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &   ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
        nwgt=1
c      endif

c!! llw, and llwI are for L(omega) for Q0P in PRB81,125102
c      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )


!! Calculate x0(q,iw) and W == main loop 1001 for iq. 
 !! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
 !! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
 !! (or need to modify do 1001 loop).
 !! iq>nqibz for ixc=11 is not time-consuming (right???)
      call MPI__hx0fp0_rankdivider2(iqxini,iqxend)


cccc Okumura 2017/06/13 for gnuplot (3D)
      ifchipm_wan =ifile_handle()
      open(ifchipm_wan, file="wan_ChiPM.mat")          
      ifchipmz_wan=ifile_handle()
      open(ifchipmz_wan,file="wan_ChiPMz.mat")          
      print *,'ifchipm=',ifchipm_wan,ifchipmz_wan
      write(ifchipm_wan,*) '# nwf:',nwf
      write(ifchipm_wan,*) '# nqbz:',nqbz

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do is=1,nspinmx
         call wan_readeigen(qbz(:,:),nqbz,is,nsp_w)
         call write_qdata(ginv,nqbz,qbz(:,:))
      enddo      
!      call rx("=== OK! test end ===")  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          
!! ======== Loop over iq ================================
      do 1001 iq = iqxini,iqxend ! NOTE: q=(0,0,0) is omitted when iqxini=2
c        if(.not. MPI__task(iq) ) cycle
c        if(normalm) then
c          ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
c          ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
c        endif

!        call cputid (0)
        q = qibz(:,iq) !qibze ! you can spefify any q, which you like.

!!!   For Fe (x,0,0)
c$$$        if (q(2)/=0 .and. q(2)/=1) cycle
c$$$        if (q(3)/=0 .and. q(3)/=1) cycle

!!!  For Fe (x,x,0) !! G to N
        if (q(1)/=0) cycle
        if (q(2)/=q(3)) cycle
!     if (q(1)>1.0) q(1)=q(1)-1d0

c$$$        !!! only Gamma
c$$$        if (q(1)/=0) then
c$$$           cycle
c$$$        endif
        
!!! |q+G|
c$$$        if (wan) then !!!
c$$$           if (iq==1) then !only Gamma
c$$$              call readqg('QGcou', q, ginv,  quu,ngc,ngveccB)
c$$$           endif
c$$$        else
c$$$           call readqg('QGcou', q, ginv,  quu,ngc,ngveccB)
c$$$        endif
        ngc=1
        if (egasmode) then
           if(.not.allocated(ngveccB)) allocate(ngveccB(3,1))
           ngveccB=0        
           nwf=1
        endif
        

c     endif
        
c        do iwf=1,ngc
c           print *,"iq,igc,ngveccB:",iq,iwf,ngveccB(:,iwf)
c        enddo
        write(6,*)'iq=',iq,q
        
        
 !! Caution : confusing point
 !!  ngc by QGcou is shown at the bottom of lqg4gw.
 !!  ngc read from PPOVL are given by rdata4gw---> ngc(iq>nqibz )=ngc for q=0
        if( iq==1 ) then ! *sanity check
          if(sum(q**2)>1d-10) then
            call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
          endif
        endif

!! Readin Coulomb matrix ====
c        ngb = nbloch + ngc !ngb is readin from vcoul 25jan2006
c        write(6,*)
        write(6,"('===== do 1001: iq q=',i7,3f9.4,' ========')")iq,q !qq
c        write(6,"('  nbloch ngb ngc=',3i10)") nbloch,ngb,ngc

        
c$$$!! Readin diagonalized Coulomb interaction ===
c$$$ !! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
c$$$        if((.not.chipm)) then
c$$$          vcoudfile='Vcoud.'//charnum5(iq)  !this is closed at the end of do 1001.  iq was iqqv
c$$$          ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$          read(ifvcoud) ngb0
c$$$          if( ngb0/=ngb ) call rx( 'hx0fp0.m.f:ngb0/=ngb')
c$$$          read(ifvcoud) qvv
c$$$          if(sum(abs(qvv-q))>1d-10) then
c$$$            write(6,*)'qvv =',qvv
c$$$            call rx( 'hx0fp0: qvv/=0 hvcc is not consistent')
c$$$          endif
c$$$          if(allocated(zcousq)) deallocate( zcousq,vcousq )
c$$$          allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
c$$$          read(ifvcoud) vcousq
c$$$          read(ifvcoud) zcousq
c$$$          idummy=iclose(trim(vcoudfile))
c$$$          vcousq=sqrt(vcousq)
c$$$          if(allocated(zzr)) deallocate(zzr)
c$$$          allocate(zzr(1,1)) !dummy
c$$$          zzr=0d0
c$$$        endif


        nmbas=1 !dummy. matrix dimension.
        nmbas1=nmbas
        nmbas2=nmbas
c        print *,'bbbbbbbb1111'
c        if(epsmode.and.nolfco) then 
c          continue
c        else
c          write(6,*) "rcxq alloc ngb nwhis npm ---",ngb,nwhis,npm
c          allocate( rcxq(ngb,ngb,nwhis,npm) )
c        endif


!! rcxq: imaginary part after x0kf_v4h and symmetrization. 
!! zxq and zxqi are the main output after Hilbert transformation
c        if(nolfco) then
c          if(allocated(rcxq)) deallocate(rcxq)
          if(allocated(zxq) ) deallocate(zxq)
          if(allocated(zxqi) ) deallocate(zxqi)
          allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
          allocate( zxq (nmbas1,nmbas2,nw_i:nw), zxqi (nmbas1,nmbas2,niw))
c        else
c          allocate( zw0(ngb,ngb), zxq (ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
c        endif
          zxq=0d0;  zxqi=0d0;  rcxq = 0d0


          if(debug) write(6,*)' niw nw=',niw,nw

          
!! ==== spin chi_charge or chi_+- ====
          is=1
          isf=2
c          isf = is
c          if(chipm) then
c            write(6,*)" chi_+- mode nolfc=",nolfco
c            if(is==1) isf=2
c            if(is==2) isf=1
c          endif

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
c$$$

!! okumura
          !evud_w(nwf,nwf,nqbz,nsp)
!     wan=.false.

ccccccccccccccc simple band for check

          if(.not. allocated(ev_w1)) allocate(ev_w1(nwf,nqbz))
          if(.not. allocated(ev_w2)) allocate(ev_w2(nwf,nqbz))
          if(.not. allocated(evc_w1)) allocate(evc_w1(nwf,nwf,nqbz))
          if(.not. allocated(evc_w2)) allocate(evc_w2(nwf,nwf,nqbz))
          
          
          !! get Efermi, qfermi, rs of electron gas
          if (egasmode) then
             ntot=1
             ntot_r=ntot/2.0
             call efermi_egas(ntot_r,alat,plat,efz)
             ef=efz             !ef for electron gas
             diag=.true.
          endif
!!
          print *,"npm,ngc::",npm,ngc
          write(6,*) "nctot,ecore",nctot,ecore
          if (wan) then
c$$$             !!! avoid n1b=0 or n2b=0 check
c$$$             evc_w1(1:nwf,0,1:nqbz)=0d0
c$$$             evc_w2(1:nwf,0,1:nqbz)=0d0
             !!! simple band for chec 1:simple, 2:gas, 3:else
             print *,"wan: EFERMI, bsimple",ef,bsimple
             if (bsimple==1) then
                ev_w1(1:nwf,1:nqbz)=3.0+ef
                ev_w2(1:nwf,1:nqbz)=4.0+ef
                ev_w1(1,1:nqbz)=-2.0+ef
                ev_w2(1,1:nqbz)=-1.0+ef
                ev_w1(2,1:nqbz)=1.0+ef
                ev_w2(2,1:nqbz)=2.0+ef
                ev_w1(3,1:nqbz)=2.0+ef
                ev_w2(3,1:nqbz)=3.0+ef
                print *,"bsimple pass"

             elseif (bsimple==2) then
            !!! qgsq: |q+G|^2
                allocate(qgsq1(ngc,nqbz),qgsq2(ngc,nqbz))
                qgsq1=0d0;qgsq2=0d0
                ev_w1=0d0;ev_w2=0d0
c$$$                ifgas=ifile_handle()
c$$$                if (is==1 .and. iq==54) open(ifgas,file="band_electrongas.dat")                
                do kx=1,nqbz
                   call read_qgband(alat,plat,qbz(:,kx),ngc,ngveccB,is,qgsq1(1:ngc,kx),qsh1)
                   call read_qgband(alat,plat,q+qbz(:,kx),ngc,ngveccB,isf,qgsq2(1:ngc,kx),qsh2)
c                   write(6,"('qqq111=',3f9.4,2x,3f9.4)") qbz(:,kx), qsh1
c                   write(6,"('qqq222=',3f9.4,2x,3f9.4)")q+qbz(:,kx),qsh2

                   !!! replace eigenvalue
                   ev_w1(1:nwf,kx)=qgsq1(:,kx) !Rydberg
                   ev_w2(1:nwf,kx)=qgsq2(:,kx) !Rydberg

                enddo
                if (allocated(qgsq1)) deallocate(qgsq1)
                if (allocated(qgsq2)) deallocate(qgsq2)
                
             else               !!! else
                iwf=ifile_handle()

                print *,"wan_readeval start"
                
                do kx=1,nqbz
                   call wan_readeval(  qbz(:,kx), is,
     o                  ev_w1(1:nwf,kx), evc_w1(1:nwf,1:nwf,kx))
                   call wan_readeval(q+qbz(:,kx), isf,
     o                  ev_w2(1:nwf,kx), evc_w2(1:nwf,1:nwf,kx))

                   if(iq==1 .and. kx==1) open(iwf,file="waneval_check.data")
!     if (kx==1 .and. iq > 100 ) open(iwf,file="wan_eval_check.data")
                   if (qbz(2,kx)==1.0) then
                      if (qbz(3,kx)==1.0) then
!                   if (iq > 100) then

                         write(iwf,"('q(1) ev_w1(:,kx)',6f9.4)") qbz(1,kx),ev_w1(:,kx)
                      endif
                   endif
                   if (iq==1 .and. kx==nqbz) close(iwf)
!!!                   if (iq==1 .and. kx==nqbz) call rx("check end")
c$$$  !! band shift for normalization check
c$$$                   ev_w1(:,kx)=ev_w1(:,kx)+2d0
c$$$                   ev_w2(:,kx)=ev_w2(:,kx)-2d0
                   
!!!   for wannier eigenvector (evc_w1, evc_w2) for M(:)

c$$$                   if (kx>10 .and. kx<30) then
c$$$                      !sumcheck
                   do iwf=1,nwf
                      do jwf=1,nwf
                         if(iwf==jwf) cycle
                         write(6,"('sum evc1,evc2',2I4,4f9.4)") iwf,jwf,
     &                        abs(sum((evc_w1(:,iwf,kx))*dconjg((evc_w1(:,jwf,kx))))),
     &                        abs(sum((evc_w2(:,iwf,kx))*dconjg((evc_w2(:,jwf,kx)))))
c     &                        sum(abs(evc_w2(:,iwf,kx)))
                      enddo
                   enddo
c$$$  endif
                enddo

             endif
             
c$$$             !!!!!!!
c$$$             do iwf=1,nwf
c$$$                write(6,"('max[ev_w1(',I1,',:)]',f13.5) =") iwf,minval(ev_w1(iwf,:))
c$$$                write(6,"('min[ev_w2(',I1,',:)]',f13.5) =") iwf,maxval(ev_w2(iwf,:))
c$$$             enddo
c$$$             !!!!!!!
             
             write(6,*) "===start wan_gettetwt",nwf,iq,q
c             if(iq==54) then
c             do kx=1,nqbz
c               write(6,"('evwwwww1',3f10.4)") ev_w1(1,kx)-ef, ev_w2(1,kx)-ef, ev_w1(1,kx)-ev_w2(1,kx)
c             enddo
c             endif
!             wan_ecore=0.0;nbmx=nwf

!!!!  shift band (2017, Sep, 13)
             if (.false.) then
                ev_w1=ev_w1-(maxval(ev_w1)-minval(ev_w2))-0.001
                ev_w2=ev_w2+(maxval(ev_w1)-minval(ev_w2))+0.001
             endif
             write(6,*) "minvalmaxval:",maxval(ev_w1),minval(ev_w2)

             iwf=ifile_handle()
             open(iwf,file="ev_band1.data")
             jwf=ifile_handle()
             open(jwf,file="ev_band2.data")
             do kx=1,10
                write(iwf,"(3f9.4,5f9.4)") qbz(:,kx),ev_w1(:,kx)
                write(jwf,"(3f9.4,5f9.4)") qbz(:,kx),ev_w2(:,kx)             
             enddo
             close(jwf)
             close(iwf)

!             call rx("test end properly")
             
cccccccccccccccccccccccccccccccccccc
             call gettetwt(q,iq,isdummy,isdummy,nwgt(:,iq),frhis,nwhis,npm,
     i            qbas,ginv, ef, nqibz, nwf,ev_w1,ev_w2, nctot,wan_ecore,
     i            nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
     i            nwf,ebmx,mtet,eibzmode) !nov2016
             write(6,*) "=== end gettetwt"
             deallocate(ev_w1,ev_w2)
cccccccccccccccccccccccccccccccccccc

c$$$          else
c$$$             do kx = 1, nqbz
c$$$c                call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) ) 
c$$$c                call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
c$$$             enddo
c$$$             call gettetwt(q,iq,isdummy,isdummy,nwgt(:,iq),frhis,nwhis,npm,
c$$$     i            qbas,ginv, ef, nqibz, nband,ekxx1,ekxx2, nctot,ecore,
c$$$     i            nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
c$$$     i            nbmx,ebmx,mtet,eibzmode) !nov2016
           endif
          
!! end okumura
c$$$  !!! test code to print tetrahedron weight !!! okumura
!     do k=1,nqbz !k-vector loop?

          allocate(kmat(1:nwf,1:nwf,1:nwhis,1:npm),wanmat(1:nwf,1:nwf))
          kmat=0d0; wanmat=0d0
          
c          print *,'sizen1b =',size(n1b),size(n2b),shape(n1b),shape(n2b)
c          print *,'sizen1b =',sum(n1b),sum(n2b)
          allocate(wanij(1:nwf,1:nwf))
          do 2011 jpm=1,npm     !! jpm=2: negative frequency
          wsum=0d0
          do 2012 kx=1,nqbz !! discrete k-point loop
          wanij=(0d0,0d0)
          do 2013 ibib=1,nbnb(kx,jpm) !! n,n' band loop
!     n1b(ibib,k,jpm) = n :band index for k (occupied)
!     n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
                if (debug) then
                   write(6,*) 'jpm,ibib',jpm,ibib
                   write(6,*) 'kx,nbnb',kx,nbnb(kx,jpm)
                   if (jpm==1) then
                      write(6,*) 'n1b(occ),n2b(unocc):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                   else
                      write(6,*) 'n1b(unocc),n2b(occ):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                   endif
                endif

!!! for M(:,:) ???
c                print *,'ibib kx jpm=',ibib,kx,jpm,' n1b(occ),n2b(unocc):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                it=n1b(ibib,kx,jpm)  !index for n  for q
                itp=n2b(ibib,kx,jpm) !index for n' for q+k

c========================== caliculate M(:,:) : wanmat(nwf,nwf)
c matrix elements from wannier eigenvector

                diag=.false.
                if (diag) then
                   do iwf=1,nwf
                      wanmat(iwf,iwf)=(1d-1,0d0)
                   enddo
                else
                   
c$$$  !!! evc_w1: upspin ; evc_w2: downspin
c$$$  !!! from down(beta) to upp(alpha)
c$$$  !!! module ?
 !                  print *,"create M (jpm,kx,ibib)=",jpm,kx,ibib
 !                  print *,"it,itp =",it,itp

c$$$                   do iwf=1,nwf
c$$$                      write(6,"('wan_i:',5f9.4)") evc_w1(iwf,:,kx)
c$$$                   enddo


                   wanmat=0d0
                   do 3001 iexc=1,1 !Direct and exchange
                      do 3002 iwf=1,nwf
                         do 3003 jwf=1,nwf
                            wan_j=dconjg(evc_w2(jwf,it,kx)) !a_{Rjbeta}^{kn}*
                            wan_i=evc_w1(iwf,itp,kx) !a_{Rialpha}^{(k+q)n'}
                            if (iexc==1) then !!!   direct
                               wan_l=evc_w2(jwf,it,kx) !a_{R'lbeta}^{kn}
                               wan_k=dconjg(evc_w1(iwf,itp,kx)) !a_{R'kalpha}^{(k+q)n'}*  
                            elseif (iexc==2) then      !!! exchange
                               wan_l=evc_w2(iwf,it,kx) !a_{R'lbeta}^{kn}
                               wan_k=dconjg(evc_w1(jwf,itp,kx)) !a_{R'kalpha}^{(k+q)n'}*
                            endif
                            wanijkl=wan_j*wan_l*wan_i*wan_k

!                            write(6,"('wansize:',4f9.4)") abs(wan_i),abs(wan_j),abs(wan_k),abs(wan_l)
                           
cccccccccccccccc
!                            if (log(abs(wanijkl))<-0.5) wanijkl=0d0
!     if (log(aimag(wanijkl))<-5) wanijkl=0d0
                            
                            if (iexc==2) wanij(iwf,jwf)=wanij(iwf,jwf)+wan_i*wan_j
                            if (iexc==2 .and. iwf==jwf) cycle
                            
                            wanmat(iwf,jwf)=wanmat(iwf,jwf)+wanijkl
!                            write(6,"('iwf,jwf,wanijkl',2I3,2E13.4)") iwf,jwf,wanijkl
cccccccccccccccc
                            
c$$$  if (iq==1 .and. ibib==1) then
c$$$  write(6,*) "wanmat for Magnon(log-scale):",log(wanmat(nwf,1))
c$$$  endif

 3003                    continue  !jwf
 3002                 continue     !iwf
 3001              continue        !iexc                     
                endif
                !!! writing wanmat (Re and Im)
                if (jpm==1. and. kx==1) then         
                   write(6,*) "Re[wanmat(iwf,jwf)]"
                   do iwf=1,nwf
                      write(6,"(5E11.3)") real(wanmat(iwf,:))
                   enddo
                   write(6,*) "Im[wanmat(iwf,jwf)]"
                   do iwf=1,nwf
                      write(6,"(5E11.3)") aimag(wanmat(iwf,:))
                   enddo
                   
                endif

c========================== end caliculate M(:,:) : wanmat(nwf,nwf)

!!!   print *,"ihw,nhw",ihw(ibib,kx,jpm),nhw(ibib,kx,jpm)
                do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
                   imagweight=whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
                   
!!!   accumulate Im[K]
                   kmat(:,:,iw,jpm) = kmat(:,:,iw,jpm)+imagweight*wanmat(:,:) 
!     rcxq(1,1,iw,jpm) = rcxq(1,1,iw,jpm) + imagweight
!                   rcxq=kmat

!                   wsum(1)=wsum(1)+imagweight !!! normalization check (wsum=pi)
                enddo
 2013        continue !ibib-loop
             if (jpm==1) then
                write(6,*) "abswan(i,j),ik",kx
                do iwf=1,nwf
                   write(6,"(5E11.3)") abs(wanij(iwf,:))
                enddo
             endif
             write(6,*) "sum abswan(i,j)",sum(abs(wanij(:,:)))
             
 2012     continue !k-loop
!          write(6,"('jpm =',I3,' wsum(:)==',f9.4)") jpm,wsum(1) !normalization
          write(6,"('jpm kmatrix_sum =',I3,2E14.5)") jpm,sum(kmat(:,:,:,jpm))
ccccc rcxq = tr(kmat)
          t2g=.true. !!! switch t2g or eg
          do iw=1,nwhis
             do iwf=1,nwf
                do jwf=1,nwf
                   if (t2g) then !!! t2g only : iwf=1,2,4
                      if (iwf/=1 .and. iwf/=2 .and. iwf/=4) cycle
                      if (jwf/=1 .and. jwf/=2 .and. jwf/=4) cycle
                   else          !!!  eg only : iwf=3,5
                      if (iwf/=3 .and. iwf/=5) cycle
                      if (jwf/=3 .and. jwf/=5) cycle
                   endif
                   if (iwf/=jwf) cycle
!     if (iwf<7 .and. iwf==8) cycle
c                rcxq(1,1,iw,jpm) = rcxq(1,1,iw,jpm) + kmat(iwf,iwf,iw,jpm)
                      rcxq(1,1,iw,jpm) = rcxq(1,1,iw,jpm) + kmat(iwf,jwf,iw,jpm)
                   enddo
             enddo
c             rcxq(1,1,iw,jpm) = rcxq(1,1,iw,jpm) + sum(kmat(:,:,iw,jpm))
          enddo
          write(6,"('jpm rcxq_sum    =',I3,2E14.5)") jpm,sum(rcxq(1,1,:,jpm))

c$$$          do iw=1,nwhis
c$$$             print *,"--- iw",iw
c$$$             if (real(rcxq(1,1,iw,jpm))>100) then
c$$$                print *,"Re rcxq",real(rcxq(1,1,iw,jpm))
c$$$             endif
c$$$             if (aimag(rcxq(1,1,iw,jpm))/=0) then
c$$$                print *,"Im rcxq",aimag(rcxq(1,1,iw,jpm))
c$$$             endif
c$$$          enddo
      
c$$$  do iwf=1,nwf
c$$$             do jwf=1,nwf
c$$$!                if (iwf/=jwf) cycle
c$$$                write(6,"('pppp iwf,jwf,sum(kmat)',2I2,2E13.5)") iwf,jwf,sum(kmat(iwf,jwf,:,jpm))
c$$$             enddo
c$$$          enddo

 2011  continue
       
       Deallocate(evc_w1,evc_w2)

!       call rx("kmatsum end")      
!!-------------------------------------       
       if (wan) then            !writing Kmatrix(kmat) 
          if (iq==1 .and. is==1) write(ifchipm_wan,*) '# nwhis:',nwhis
!     if (iq<=4) then    !check X to Gamma
          if (iq>=0) then       !check X to Gamma
             do iw=1,nwhis
                write(ifchipm_wan,"(f9.4,4E13.5)") q(1),frhis(iw),
     &               sum(kmat(:,:,iw,1)), (frhis(iw+1)-frhis(iw))
!     write(ifchipm_wan,"(f9.4,2E13.5)") q(1),-1d0*frhis(iw)*2.0,sum(kmat(:,:,iw,2)) !Rydberg unit?
             enddo
             write(ifchipm_wan,*) ""
          endif
       endif

       
!     if (iq>100) call rx("test end")          

c$$$  !! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 into rcxq ==
c$$$  call cputid(0)
c$$$  if(npm==1) then
c$$$  ncc=0
c$$$  else
c$$$  ncc=nctot
c$$$  endif
c$$$  call x0kf_v4hz(npm,ncc,   
c$$$  i        ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
c$$$  i        n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,

c$$$  i              q,  
c$$$  i              nspin,is,isf, !symmetrize, !
c$$$  i              qbas,ginv,  qbz,wbz, 
c$$$  d            nlmto,nqbz,nctot, !noccx,noccxv,
c$$$  d            nbloch,  nwhis, !nlnmx,mdimx,
c$$$  i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
c$$$  i    nqbze,nband,nqibz, 
c$$$  o    rcxq,   !  rcxq is the accumulating variable for spins 
c$$$  i    nolfco,zzr,nmbas_in, zcousq,   !ppovl,nmbas1,nmbas2, is removed ppovlz, 
c$$$  i    chipm,eibzmode, !z1offd,!for nolfco Add nmbas Sep2006
c$$$  i    nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),ngrp, eibzsym(:,:,iq),crpa)
       
!     ! ------Question, Apr2015takao. -------------------
!     ! ???  we may need 
!     ! ???  "if(is==nspinmx.or.chipm) then" for chipm mode.
!     ! ???  really OK ??? Need check more... Compare with old code...
!     ! ---------------------------------------------------

!     ! == Symmetrizer for crystal symmetry (and also for spin)
!     !    Symmetrize and convert to Enu basis by dconjg(tranpsoce(zcousq)*rcxq8zcousq if eibzmode
c$$$  if(is==nspinmx.or.chipm) then !Apr2015. TK think ".or.chipm" is required for chipm mode
c$$$  !Because rcxq is calculated for each is, symmetrized and its contribution
c$$$  !is added to zxq in dpsion5.
c$$$  call x0kf_v4hz_symmetrize(npm, !ncc,
c$$$  i       q,
c$$$  i       nspin,is,isf, !symmetrize, !
c$$$  i       qbas,ginv,  !qbz,wbz,
c$$$  d       nbloch,  nwhis,    ! nlnmx,mdimx,
c$$$  i       iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
c$$$  i       nqbze,nband,nqibz,
c$$$  o       rcxq,              !  rcxq is the accumulating variable for spins
c$$$  i       nolfco,zzr,nmbas_in, zcousq, !ppovl,nmbas1,nmbas2, is removed ppovlz, 
c$$$  i       chipm,eibzmode,    !z1offd,!for nolfco Add nmbas Sep2006
c$$$  i       ngrp, eibzsym(:,:,iq))
c$$$  endif
       call tetdeallocate()     !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
c$$$  iecut=1
c$$$  if(debug) write(6,"(a)") ' --- goto dpsion5 --- '
c$$$  if(is==nspinmx.or.chipm) then
c$$$  write(6,"('  nmbas1,nmbas2=',2i10)") nmbas1,nmbas2
c$$$  call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
c$$$  i       rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
c$$$  o       zxq, zxqi,
c$$$  i       chipm, schi,is,  ecut(iecut),ecuts(iecut))
c$$$  if(nolfco.and.epsmode) then
c$$$  do iw=nw_i,nw
c$$$  x0mean(iw,:,:)=zxq(:,:,iw)
c$$$  enddo
c$$$  endif
c$$$  write(6,*)' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
c$$$  endif

!     ! dpsion5. zxq, zxqi are accumulating for spins.

ccccccccccccccccccccccccccccccccccccccccccccccccc
!     ! npm is set by npmtwo=T
       write(6,"('  nmbas1 nmbas2 npm=',3i8)") nmbas1,nmbas2,npm
       schi=1                   ! flip over maj/min spins.
       call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i      rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o      zxq, zxqi,
     i      .false., schi,1,1d99,1d99) !ecut(iecut),ecuts(iecut))
!     !
c     1003   continue  !end of spin loop =====

       do iw = nw_i,nw
          www=freq_r(iw)
          if(iw<0) www=-freq_r(-iw)
          write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,2E13.4)")q,iw,www*2d0,zxq(1,1,iw)
       enddo
       write(ifchipmz_wan,*)
       
ccccccccccccccccccccccccccccccccccccccccccccccccc

       if(allocated(kmat)) deallocate(kmat)
       if(allocated(wanmat)) deallocate(wanmat)

       if(allocated(wanij)) deallocate(wanij)
       
       if(allocated(rcxq) ) deallocate(rcxq)
c$$$  
c$$$  
c$$$  !! ===  RealOmega ====================================
c$$$  if (realomega) then
c$$$  if(epsmode.and.nolfco.and.(.not.chipm)) then 
c$$$  if (nspin==1) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
c$$$  if (nspin==1) zxq = 2d0*zxq       !if paramagnetic, multiply x0 by 2
c$$$  else
c$$$  if (nspin == 1) zxq = 2d0*zxq   !if paramagnetic, multiply x0 by 2
c$$$  endif
c$$$  
c$$$  if(epsmode) then
c$$$  if(nolfco) then
c$$$  ttt='without LFC'
c$$$  else
c$$$  ttt='with LFC'
c$$$  endif
c$$$  if(chipm) then
c$$$  write(6,*) '--- chi0_{+-}}^{-1}      --- '//ttt
c$$$  else
c$$$  write(6,*) '--- dielectric constant --- '//ttt
c$$$  write(6, *)" trace check for W-V"
c$$$  endif
c$$$  endif
c$$$  
c$$$  !! prepare for iq0.
c$$$  iq0 = iq - nqibz
c$$$  if(allocated(epstilde)) deallocate(epstilde,epstinv)
c$$$  allocate(epstilde(ngb,ngb),epstinv(ngb,ngb))
c$$$  
c$$$  !! === iw loop for real axiw ===
c$$$  do 1015 iw  = nw_i,nw  !Feb2006. Before it was 1:nwp (nwp=nw+1).
c$$$  !  So freq_r(iw-1) is shifted to freq_r(iw).
c$$$  frr= dsign(freq_r(abs(iw)),dble(iw))
c$$$  if(.not.epsmode) then  
c$$$  !! === wcf: W= (1-v zxq)^{-1} v ===
c$$$  c              if(newaniso2.and.iq<=nqibz) then !for mmmw
c$$$  if(iq<=nqibz) then !for mmmw
c$$$  if(iq==1) then
c$$$  ix=1
c$$$  zw0(:,1)=0d0
c$$$  zw0(1,:)=0d0
c$$$  else
c$$$  ix=0
c$$$  endif
c$$$  !!  Eqs.(37),(38) in PRB81 125102
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  
c$$$  
c$$$  do igb1=1+ix,ngb
c$$$  do igb2=1+ix,ngb
c$$$  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$  enddo
c$$$  enddo
c$$$  endif
c$$$  
c$$$  !! for iq>nqibz Full inversion to calculalte eps with LFC.
c$$$  if(iq>nqibz) then
c$$$  ix=0
c$$$  vcou1 = fourpi/sum(q**2*tpioa**2) !--> vcousq(1)**2! !fourpi/sum(q**2*tpioa**2-eee)
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  if(igb1==1.and.igb2==1) then
c$$$  epstilde(igb1,igb2)= 1d0 - vcou1*zxq(1,1,iw)
c$$$  cycle
c$$$  endif
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) then
c$$$  epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
c$$$  endif   
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  llw(iw,iq0)= 1d0/epstinv(1,1)
c$$$  write(6,"('iq iw_real eps(withLFC) eps(woLFC) ',2i5,d13.6,x,d13.6,2x,d13.6,x,d13.6)")
c$$$  &           iq,iw,llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
c$$$  endif
c$$$  if(iq>nqibz) then
c$$$  continue
c$$$  else
c$$$  zw(1:ngb,1:ngb) = zw0
c$$$  write(ifrcw, rec= iw-nw_i+1) zw   !  WP = vsc-v
c$$$  call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
c$$$  endif
c$$$  !! epsmode
c$$$  elseif(epsmode.and.(.not.chipm)) then 
c$$$  if(debug)write(6,*) 'xxx2 epsmode iq,iw=',iq,iw
c$$$  !! there is difference of two vcmean below since we use (sligthy) screened Coulomb (screenfac() in switch.F)
c$$$  !!   NOTE that we use vcoul with screening (screenfac() is used in hvccfp0.F
c$$$  vcmean=vcousq(1)**2 !fourpi/sum(q**2*tpioa**2) !aug2012
c$$$  epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxq(1,1,iw))
c$$$  write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2e23.15,3x, 2e23.15,
c$$$  &          " vcmean x0mean =", 2e23.15,3x, 2e23.15)') iqixc2,iw,2*frr,
c$$$  &          1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean, zxq(1,1,iw)!x0mean(iw,1,1)
c$$$  write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
c$$$  &          q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$  if(.not.nolfco) then
c$$$  ix=0
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  if(igb1==1.and.igb2==1) then
c$$$  epstilde(igb1,igb2)= -vcmean*zxq(igb1,igb2,iw) !aug2012
c$$$  else   
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  endif  
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  epsi(iw,iqixc2)= epstinv(1,1)
c$$$  write(6,'( " iq iw omega eps epsi  wLFC="
c$$$  &           ,2i6,f8.3,2e23.15,3x, 2e23.15)') 
c$$$  &            iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$  write(6,*)
c$$$  write(ifepsdat,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
c$$$  &            q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$  endif
c$$$  !! ChiPM mode
c$$$  elseif(epsmode.and.chipm) then
c$$$  allocate( x0meanx(nmbas,nmbas) )
c$$$  if(nolfco) then  ! ChiPM mode without LFC
c$$$  c$$$c! --- three lines below may work for test purpose for legas. But not sure.
c$$$  c$$$c       vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
c$$$  c$$$c       write(ifchipmn,'(3f12.8,2x,f8.5,2x,2e23.15)')
c$$$  c$$$c     & q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw,1,1)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
c$$$  x0meanx = x0mean(iw,:,:)/2d0 !in Ry unit.
c$$$  else
c$$$  !! ChiPM mode with LFC... NoLFC part
c$$$  zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)/2d0  ! in Ry.
c$$$  do imb1=1,nmbas
c$$$  do imb2=1,nmbas
c$$$  x0meanx(imb1,imb2)= 
c$$$  &              sum( svec(1:nbloch,imb1)*
c$$$  &              matmul(zxq(1:nbloch,1:nbloch,iw),svec(1:nbloch,imb2))) 
c$$$  enddo
c$$$  enddo
c$$$  !  x0meanx= <m|chi^+-(\omega)|m>/<m|m>**2
c$$$  endif
c$$$  do imb1=1,nmbas
c$$$  do imb2=1,nmbas
c$$$  x0meanx(imb1,imb2) = 
c$$$  &            x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
c$$$  enddo
c$$$  enddo
c$$$  !!! Check!!! okumura 
c$$$  write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,255e23.15)')q, 2*schi*frr, x0meanx(:,:)
c$$$  if(.not.nolfco) write(ifchipm_fmat) q, 2*schi*frr, zxq(1:nbloch,1:nbloch,iw)
c$$$  deallocate(x0meanx)
c$$$  endif
c$$$  1015     continue
       if( allocated(zzr)   ) deallocate(zzr)
       if( allocated(x0mean)) deallocate(x0mean)
       if( allocated(gbvec) ) deallocate(gbvec)

c$$$  endif
c$$$  !! ... Close files for epsmode
c$$$  if(epsmode .and. .not. wan) then 
c$$$  if(chipm) then
c$$$  ifchipmn_mat=iclose('ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
c$$$  !             ifchipmn_mat=iclose('ChiPM1111.nlfc.mat')
c$$$  if(.not.nolfco) then
c$$$  ifchipm_fmat=iclose( 'ChiPM'//charnum4(iqixc2)//'.fmat')
c$$$  endif
c$$$  else
c$$$  filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
c$$$  ifepsdatnolfc = iclose( filepsnolfc)
c$$$  if(.not.nolfco) then
c$$$  fileps = 'EPS'//charnum4(iqixc2)//'.dat'
c$$$  ifepsdat  = iclose(fileps)
c$$$  endif
c$$$  endif
c$$$  endif
c$$$  !! RealOmega end ==============================
c$$$  
c$$$  
c$$$  
c$$$  !! ImagOmega start ============================
c$$$  if (imagomega) then
c$$$  write(6,*)' goto imag omega'
c$$$  if (nspin == 1) zxqi = 2d0*zxqi   ! if paramagnetic, multiply x0 by 2
c$$$  !! === iw loop for imag axiw ===
c$$$  do 1016 iw  = 1,niw
c$$$  if( iq<=nqibz ) then
c$$$  !!  Eqs.(37),(38) in PRB81 125102
c$$$  if(iq==1) then
c$$$  ix=1
c$$$  zw0(:,1)=0d0
c$$$  zw0(1,:)=0d0
c$$$  else
c$$$  ix=0
c$$$  endif
c$$$  !!  Eqs.(37),(38) in PRB81 125102
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv=epstilde
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$  enddo
c$$$  enddo
c$$$  endif
c$$$  if(iq>nqibz) then
c$$$  !! Full inversion to calculalte eps with LFC.
c$$$  ix=0
c$$$  vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  if(igb1==1.and.igb2==1) then
c$$$  epstilde(igb1,igb2)= 1d0 - vcou1*zxqi(1,1,iw)
c$$$  cycle
c$$$  endif
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) then
c$$$  epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
c$$$  endif   
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  llwI(iw,iq0)= 1d0/epstinv(1,1)
c$$$  write(6,*) 'iq iw_img  eps(withLFC) eps(woLFC)',iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
c$$$  endif
c$$$  if(iq>nqibz) then
c$$$  continue
c$$$  else
c$$$  zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
c$$$  write(ifrcwi, rec=iw)  zw    !  WP = vsc-v
c$$$  call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
c$$$  endif
c$$$  
c$$$  !!!
c$$$  1016     continue
c$$$  deallocate(epstinv)
       if(allocated(epstilde)) deallocate(epstilde)
c$$$  endif
c$$$  !! ImagOmega end =================
       if(allocated(vcoul)) deallocate(vcoul)
       if(allocated(zw0)) deallocate(zw0)
       if(allocated(zxq )) deallocate(zxq)
       if(allocated(zxqi)) deallocate(zxqi)
c$$$  if (normalm) then
c$$$  ifrcwi = iclose('WVI.'//charnum5(iq))
c$$$  ifrcw  = iclose('WVR.'//charnum5(iq))
c$$$  endif
       
 1001 continue !q point loop

      close(ifchipm_wan)
      close(ifchipmz_wan)

      
      !!! okumura
      
c$$$!! =================== end of loop 1001 for q point ========================
!!! okumura comment out below
      call MPI__barrier()

c$$$!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
c$$$      if(MPI__size/=1) then
c$$$        do iq=nqibz+1,iqxend
c$$$          iq0 = iq - nqibz
c$$$          if(MPI__ranktab(iq)/=0) then !jan2012
c$$$            if(MPI__ranktab(iq) == MPI__rank) then
c$$$              dest=0
c$$$              call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
c$$$              call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
c$$$            elseif(MPI__root) then
c$$$              src=MPI__ranktab(iq)
c$$$              call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
c$$$              call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
c$$$            endif
c$$$          endif
c$$$        enddo  
c$$$      endif

c$$$!! == W(0) divergent part and W(0) non-analytic constant part.==
c$$$!!   Note that this is only for q=0 -->iq=1
c$$$      if(normalm.and.MPI__rank==0) then
c$$$!! get w0 and w0i (diagonal element at Gamma point
c$$$!! This return w0, and w0i
c$$$        call w0w0i(llw,llwI,nw_i,nw,nq0i,niw,q0i)
c$$$!! === w0,w0i are stored to zw for q=0 ===
c$$$!! === w_ks*wk are stored to zw for iq >nqibz ===
c$$$        do iq = 1,1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
c$$$          q = qibze(:,iq)
c$$$          do ircw=1,2
c$$$            if    (ircw==1) then;  nini=nw_i;   nend=nw;
c$$$              ifrcwx = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
c$$$            elseif(ircw==2) then;  nini=1;     nend=niw;
c$$$              ifrcwx = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
c$$$            endif
c$$$            do iw=nini,nend
c$$$              read(ifrcwx, rec= iw-nini+1) zw !(1:ngb,1:ngb)
c$$$              if(ircw==1) zw(1,1) = w0(iw)
c$$$              if(ircw==2) zw(1,1) = w0i(iw)
c$$$c              if(ircw==1) then
c$$$c                 write(6,"('ffffrrr:', f13.6,2x,f13.6,x,f13.6)") hartree*freq_r(iw),w0(iw)
c$$$c              endif
c$$$              write(ifrcwx,rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
c$$$            enddo
c$$$            if    (ircw==1) then
c$$$              ifrcwx = iclose('WVR.'//charnum5(iq))
c$$$            elseif(ircw==2) then
c$$$              ifrcwx = iclose('WVI.'//charnum5(iq))
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$      endif
          
      call cputid(0)
      call MPI__Finalize
      call rx0( ' OK! hmagnon mode')

c$$$      if(ixc==11)   call rx0( ' OK! hx0fp0 mode=11     read <Q0P> normal sergeyv')
c$$$      if(ixc==111)  call rx0( ' OK! hx0fp0 mode=111    normal sergeyv')
c$$$      if(ixc==10011)call rx0( ' OK! hx0fp0 mode=10011  crpa normal sergeyv')
c$$$      if(ixc==12)   call rx0( ' OK! hx0fp0 mode=12  Ecor sergeyv mode')
c$$$      if(ixc==101)  call rx0( ' OK! hx0fp0 mode=101 Ecor ')
c$$$      if(ixc==202)  call rx0( ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC')
c$$$      if(ixc==203)  call rx0( ' OK! hx0fp0 mode=203 sergeyv eps LFC ')
c$$$      if(ixc==222)  call rx0( ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv')
c$$$      if(ixc==223)  call rx0( ' OK! hx0fp0 mode=223 chi+- wannier sergeyv')

      end

c--------------------------------------------------------------------
      real*8 function eclda_bh(rs)
      real(8) :: rs,cp,rp,z
      cp       = 0.0504d0*0.5d0 ! 0.5 changes unit from Ry to Hartree
      rp       = 30.d0
      z        = rs / rp
      eclda_bh = -cp * ( (1.d0+z**3)*dlog(1.d0+1.d0/z) 
     .                 + 0.5d0*z - z**2 - 0.33333333d0 )
      end
c--------------------------------------------------------------------
      real*8 function eclda_pz(rs)
      real(8) :: rs
      if (rs.ge.1.d0) then
        eclda_pz = -0.1423d0 / (1.d0 + 1.0529d0*dsqrt(rs) + 0.334d0*rs)
      else
        eclda_pz = -0.0480d0 + 0.0311d0*dlog(rs) - 0.0116d0 * rs 
     .         + 0.0020d0*rs*dlog(rs)
      endif
      end
c--------------------------------------------------------------------
      subroutine wecqw(ifcor,
     d                 nqibz,nqbz,nq0i,nqitot,niw,
     o                 wibz,wqt,wx,freqx,ecqw)

      implicit double precision (a-h,o-z)
      dimension   wibz(nqibz),wqt(nq0i),wx(niw),
     .            freqx(niw),ecqw(nqitot,niw)
      real(8):: rydberg
      write(ifcor,*)'### ecqw(q,w) ###'
      write(ifcor,*)'nqibz =',nqibz
      write(ifcor,*)'nq0i  =',nq0i
      write(ifcor,*)'niw   =',niw
      do ip = 2,nqitot
        if (ip <= nqibz) then
          wk = wibz(ip)*0.5d0 ! 0.5 for the normalization of wibz
        else
c        wk = wqt(ip-nqibz)*wibz(1)*0.5d0 ! 0.5 for the normalization of wibz
          wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
        endif
        write(ifcor,*)'### iq,wq = ',ip,wk
        sume=0d0
        do ix = 1,niw
          write(ifcor,*)freqx(ix),ecqw(ip,ix),wx(ix)
          sume=sume+  wx(ix)/(freqx(ix)*freqx(ix)) * ecqw(ip,ix)
        enddo
        write(ifcor,*) '  sum ecqw*wx=', wk*sume*2d0*rydberg()
! end of ip-loop
      enddo
      return
      end
c---------------------------------------------------------------------
      subroutine getsqovlp(q,ngc,ngb,sqovlp)
!! == Get sqrt of ppovl ==
      implicit none
      real(8)::q(3)
      integer:: ngc,ngb,nbloch,i,nmxx,ix,iy,nev
      complex(8):: sqovlp(ngb,ngb)
      complex(8),allocatable:: ooo(:,:),ppo(:,:),sqovlpi(:,:),ppovl(:,:)
      complex(8),allocatable:: ovlp(:,:),evec(:,:)
      real(8),allocatable:: eval(:)
      nbloch = ngb-ngc
      if(ngc==0) goto 888

      allocate(ppovl(1:ngc,1:ngc))
      call readppovl0(q,ngc,ppovl)
      allocate(ooo(ngc,ngc),ppo(ngc,ngc),evec(ngc,ngc),eval(ngc))
      ooo= 0d0
      do ix=1,ngc
        ooo(ix,ix)=1d0
      enddo
      ppo = ppovl
      deallocate(ppovl)
      nmxx = ngc
      evec = 0d0
      eval = 0d0
      call diagcv(ooo, ppo,
     &     evec, ngc, eval, nmxx, 1d99, nev)
      write(6,*)' diagcv overlap ngc nev=',ngc,nev
      deallocate(ooo,ppo)
c
 888  continue
      sqovlp=0d0
      do i=1,nbloch
        sqovlp(i,i)=1d0
      enddo
      do i=1,ngc
        if(eval(i)<0d0) then
          call rx( 'getsqovlp:  eval(i) <0d0')
        endif
        do ix=1,ngc;  do iy=1,ngc
          sqovlp(ix+nbloch,iy+nbloch)= 
     &     sqovlp(ix+nbloch,iy+nbloch)
     &      + evec(ix,i)* sqrt(eval(i))* dconjg(evec(iy,i))
        enddo ;      enddo
      enddo
      if(allocated(evec)) deallocate(evec)
      if(allocated(eval)) deallocate(eval)
      write(6,*)' end of getsqovlp'
c         sqovlpi = sqovlp
c         call matcinv(ngb,sqovlp)     !  inverse
c         ovlpi=ovlp
c         deallocate(ppovl,ovlp)
      end
c--------------------------------------------------------------------
      subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
!! == check write for zw, no output == !!
      implicit none
      integer:: nblochpmx,nbloch,ngb,iw,i,iq
      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
      real(8):: freqq
c      logical :: smbasis
      character*(*)::tagname
      trwv=0d0
c      if(.not.smbasis()) then
        do i = 1,nbloch
          trwv = trwv + zw(i,i)
        enddo
c      endif
      trwv2 = 0d0
      do i = 1,ngb
        trwv2 = trwv2 + zw(i,i)
      enddo !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(6,'(a,f10.6,2i5,4d20.12)')trim(adjustl(tagname)),freqq,iq,iw,trwv,trwv2
c      do i = 1,ngb
c        write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c      enddo
      end

c--------------------------------------------------------------------
C      subroutine test_xxx(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
C      implicit none
C      integer:: nblochpmx,nbloch,ngb,iw,i,iq
C      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
C      real(8):: freqq
C      logical :: smbasis
C      character*(*)::tagname
C      trwv2 = 0d0
C      forall( i = 1:ngb)
C        trwv2 = trwv2 + zw(i,i)
C      end forall
C      end
c--------------------------------------------------------------------

      function matcinvf(a) result(b)
!!== Test routine for Inversion ==
      implicit none
      integer :: info,n,n2(2)
      integer,allocatable :: ipiv(:)
      complex(8):: a(:,:), b(1)
      complex(8),allocatable:: work(:)
      n2= SHAPE(a)
      n=n2(1)
      call zcopy(n,b,1,a,1)
      call zgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        write(6,*)' matcinv: zegtrf info=',info
        call rx( ' matcinv: zegtrf ')
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        write(6,*)'matcinv: zegtri info=',info
        call rx( 'matcinv: zegtri ')
      endif
      end

c--------------------------------------------------------------------
      subroutine diagno00(nbloch,wpvc,eval)
!! == ontain eigenvalue only for input complex matrix wpvc(nbloch,nbloch)
      implicit none
      integer:: nbloch,nmx,nev,i
      complex(8),allocatable:: ovlpc(:,:),evecc(:,:),wpvcc(:,:)
      real(8)::emx,eval(nbloch)
      complex(8):: wpvc(nbloch,nbloch)
      allocate( ovlpc(nbloch,nbloch),evecc(nbloch,nbloch),wpvcc(nbloch,nbloch))
      wpvcc= wpvc
      ovlpc= 0d0
      do i=1,nbloch
        ovlpc(i,i)=1d0
      enddo
      eval=0d0
      nev  = nbloch
      nmx  = nbloch
      call diagcv(ovlpc,wpvcc, evecc, nbloch, eval, nmx, 1d99, nev)
      deallocate(ovlpc,evecc,wpvcc)
      end
