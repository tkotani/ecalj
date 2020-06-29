!!  Calculate x0, \epsilon, spin susceptibility. 
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_lmf_chipm mode works).
      program hx0fp0
      use m_ReadEfermi,only: Readefermi,ef
      use m_readqg,only: Readqg,Readngmx2,ngpmx,ngcmx
      use m_readeigen,only: Readeval,Init_readeigen,Init_readeigen2
      use m_read_bzdata,only: Read_bzdata,
     &     ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,
     &     dq_,qbz,wbz,qibz,wibz,qbzw,
     &     idtetf,ib1bz,idteti,
     &     nstar,irk,nstbz
      use m_genallcf_v3,only: Genallcf_v3,
     &     nclass,natom,nspin,nl,nn, ngrp,
     &     nlmto,nlnmx, nctot, !niw_in=>niw, !nw_input=>nw,
     &     alat, deltaw,esmr,symgrp,clabl,iclass, !diw,dw,delta,
     &     invg, il,in,im,nlnm, 
     &     plat, pos,ecore, symgg 
c      use m_keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
      use m_readqgcou,only: readqgcou
      use m_mpi,only: MPI__hx0fp0_rankdivider2,MPI__task,MPI__Initialize,MPI__Finalize,MPI__root,
     &     MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size,
     &     MPI__ranktab,MPI__consoleout,MPI__barrier
!     ! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: Rdpp,    !NOTE: "call rdpp" generate following data.
     &     nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
!     ! Generate matrix element for "call get_zmelt".
      use m_zmel,only:          !NOTE: these data set are stored in this module, and used  
     &      Mptauof_zmel, Setppovlz,Setppovlz_chipm !Ppbafp_v2_zmel,Setppovlz_ee,
c     & nband,itq,ngcmx,ngpmx, ppovlz, ppbir,shtvg, miat,tiat , ntq
      
      use m_itq,only: Setitq !set itq,ntq,nband,ngcmx,ngpmx to m_itq
      
!! frequency
      use m_freq,only: Getfreq3, !NOTE: call getfreq generate following data.
     & frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw,niw !, frhis0,nwhis0 !output of getfreq
!! tetwt
      use m_tetwt,only: Tetdeallocate,Gettetwt, !followings are output of 'L871:call gettetwt')
     &  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
!! w0 and w0i (head part at Gamma point)
      use m_w0w0i,only: W0w0i,
     & w0,w0i
      use m_lldata,only: ll
      use m_readgwinput,only: ReadGwinputKeys,
     &     ecut,ecuts,nbcut,nbcut2,mtet,ebmx,nbmx,nmbas,imbas
      use m_qbze,only: Setqbze,
     &     nqbze,nqibze,qbze,qibze
      use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
!! q0p
      use m_readq0p,only: Readq0p,
     &     wqt,q0i,nq0i ,nq0iadd,ixyz,nq0ix,neps
      use m_readVcoud,only: Readvcoud,vcousq,zcousq
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
     &   wgt0(:,:)!wqt(:),,q0i(:,:)
      integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     &           , ngvecp(:,:), ngvecc(:,:), !,ngpn(:),ngcni(:),iqib(:),
     &   ifppb(:)   !ongveccBr(:,:,:),nx(:,:),nblocha(:),lx(:),
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &  zw(:,:),zw0(:,:),
     &  zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: eqt(:), !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
     &  ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:) 
c     &  qbze(:,:),qibze(:,:)  !,ecore(:,:)  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:)
c     & ,rcxqmean(:,:,:,:),rcxqmeanc(:,:,:,:) !now rcxqmean is treated as a case of rcxq(nmbas,nmbas)

!  tetrahedron method
c      logical :: tetra=.true. !,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
      ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer,allocatable :: 
     &         noccxvv(:) !n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),
      real(8) ::qbzx(3),anfvec(3)
      logical :: debug=.false.
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
c     &   ifhbe,
c     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   i,ngrpmx,mxx,ini,ix,ngrpx,!ngcmx,ngpmx nq0i,nq0ix,
     &   nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, !ifvcfpout,,mdimx,nbloch
     &   ifepscond !,ifvxcpout,ifgb0vec,nxx
     &   ,nw0,iw,ifinin,iw0,ifwwk,noccxv,noccx
     &   ,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0dummy
     &   ,kx,isf,kqxx,kp,job,noccxvx(2)=-9999,nwmax  !,ifev1,ifev2 nbnbx,nhwtot,
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j
c     &   ,ngpmx !,  ifchipmlog

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      integer:: incwfin,  verbose

      integer:: ngc!,mrecg !bzcase, 
      real(8):: quu(3), deltaq(3)!,qq(3) !,qqq(3)=0d0
      logical:: omitqbz=.false., noq0p

      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)

      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg !if true, we use off-gamma mesh.
c      integer:: nbcut,nbcut2

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

      real(8),allocatable:: totexc(:), trpv(:),trlog(:) !ecut(:),ecuts(:) ,
      integer:: necut,iecut

      integer:: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
      integer,allocatable:: nxx_r(:)
      real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
      character*3:: charnum3
      character*4:: charnum4
      complex(8),allocatable:: jcoup(:,:), mcm(:,:,:)
      real(8)::chg1,chg2,spinmom,schi=1d0
c$$$#ifdef COMMONLL
c$$$      integer::ll(51**2)
c$$$      common/llblock/ll
c$$$#else
c$$$      integer :: ll
c$$$      external ll
c$$$#endif
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
      integer:: imb,imb1,imb2,nmbas_in !nmbas,
      integer,allocatable:: aimbas(:), iibas(:) !imbas_s(:),
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
      complex(8),allocatable:: epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:)!zcousq(:,:),
c      real(8),allocatable:: vcousq(:) 
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq

!! Eq.(40) in PRB81 125102
c      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:),sks(:),skI(:),sksI(:), 
     &  w_k(:),w_ks(:),w_kI(:), w_ksI(:), s_vc(:),vw_k(:),vw_ks(:)
      complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
      integer:: lxklm,nlxklm,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,iwx,icc1,icc2!,niw,niwxxx,
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
      integer,allocatable ::  invgx(:) !iclasst(:),
      integer:: ificlass,ifile_handle,k
      complex(8),allocatable:: ppovl_(:,:)

      logical:: readw0w0itest=.false.,hx0

c      real(8)::ebmx
      integer:: ifq0p,ifwc,ifif !nbmx,mtet(3),
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)
c      integer:: nqbze,nqibze
c      real(8),allocatable:: qbze(:,:),qibze(:,:)

!! -------------------------------------------------------------------
      call MPI__Initialize()
      call MPI__consoleout('hx0fp0')
      call cputid (0)
      if(verbose()>=100) debug=.true.
      hartree  = 2d0*rydberg()
      pi       = 4d0*datan(1d0)
      fourpi   = 4d0*pi
      sqfourpi = sqrt(fourpi)

!! computational mode select ! takao keeps only the Sergey mode.
      write(6,"(a)") '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'
      write(6,"(a)") ' #1:run mode'
      write(6,"(a)") '    11  : normal    '
      write(6,"(a)") '    111 : normal    fullband '
      write(6,"(a)") '    10111 : normal  crpa     '
      write(6,"(a)") '    202 : epsNoLFC  '
      write(6,"(a)") '    203 : eps       '
      write(6,"(a)") '    222 : chi^+- NoLFC'
c      write(6,"(a)") '    12  : total energy Miyake: developint'
      write(6,"(a)")  '-------------------------------------------------------'
      if( MPI__root ) read(5,*) ixc
      call MPI__Broadcast(ixc)
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
      if(ixc==11) then; write(6,*) " OK ixc=11 normal mode "
        normalm=.true.
        lqall=.false.
      elseif(ixc==111) then; write(6,*) " OK ixc=111 normal mode. fullband"
        normalm=.true.
      elseif(ixc==10011) then; write(6,*) " OK ixc=10011 crpa mode "
        normalm=.true.
        crpa=.true.
!! -- eps mode NoLFC
      elseif(ixc==202) then
        write(6,*) " OK ixc=202  sergey's eps mode Only NoLFC "
        imagomega=.false.
        omitqbz=.true.
        epsmode = .true.
        nolfco=.true.
!!    eps mode with LFC
      elseif(ixc==203) then
        write(6,*) " OK ixc=203 sergey's eps mode with LFC "
        imagomega=.false.
        omitqbz=.true.
        epsmode = .true.
!! chipm mode NoLFC
      elseif(ixc==222) then
        write(6,*) " OK ixc=222    chipm sergey's "
        imagomega =.false.
        omitqbz=.true.
        epsmode = .true.
        chipm=.true.
        nolfco=.true.
c$$$!! Total energy test mode --> need fixing
c$$$      elseif(ixc==12) then
c$$$        write(6,*) " ixc=12 Miyake's total energy Sergey--->need to fix this mode"
c$$$        call rx( " ixc=12 Miyake's total energy Sergey--->need to fix this mode")
c$$$        realomega=.false.
c$$$        ecorr_on=901
      else
        call rx( ' hx0fp0: given mode ixc is not appropriate')
      endif

c$$$!! files for RPA correlation energy mode.
c$$$      if(ecorr_on > 0) then
c$$$        ieclog = 8155
c$$$        if(ecorr_on==901) then
c$$$          ieceig=8156
c$$$          open(ieceig,file='rpa_eigen.chk')
c$$$          close(ieceig,status='delete')
c$$$        endif
c$$$        open(ieclog, file='ecorr.chk')
c$$$      endif

!! We fix newaniso2=T now.
!!  write(6,*)'Timereversal=',Timereversal() works ok?
      
!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
      call Read_BZDATA(hx0)
      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp

!! check write
      if(MPI__root) then
         do i=1,nqbz
            if(i<10.or.i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
            if(i==10.and.nqbz>18) write(6,"('... ')")
         enddo
         write(6,*)' nqbz nqibz =',nqbz,nqibz
      endif

!! EFERMI
      call Readefermi()
      write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef

!! === Readin by genallcf. Set basic data for crystal
!! See "use m_genallcf_v3" at the begining of this routine
!!
      incwfin=0  !use ForX0 for core in GWIN
      call genallcf_v3(incwfin) !in module m_genallcf_v3
!! incwfin is for core handling
c     write(6,*)' ncore nw_input delta=', ncore,nw_input,delta
      if(ngrp/=ngrp2       ) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      if(chipm.and.nspin==1) call rx( 'chipm mode is for nspin=2')
      if(debug) write(6,*)' end of genallc'
      tpioa=2d0*pi/alat
!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call rx( ' nclass /= natom ')

!! --- Read dimensions 
c      open(newunit=ifhbe,file='hbe.d')
c      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg !Note: nband is in m_zmel
c      close(ifhbe)
      call Readhbe()   
      if(nqbz /=nqbzt ) call rx( ' hx0fp0_sc: nqbz /=nqbzt  in hbe.d')
      if(nlmto/=nlmtot) call rx('hx0fp0: nlmto/=nlmtot in hbe.d')
!!
      call ReadGWinputKeys()    !Readin GWinput
      
!! --- Readin Offset Gamma --------
      call ReadQ0P()
      write( 6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)

!! Readin q+G. nqbze and nqibze are for adding Q0P related points to nqbz and nqibz.
      call Readngmx2() !return ngpmx and ngcmx in m_readqg
      write(6,"(' ngcmx ngpmx= ',2i8)") ngcmx,ngpmx 
      call Setqbze()    ! extented BZ points list

!! This check write was required to pass check for commit 299baa Tue Jun 23 21:16:02 2020 +0900
!!    gfortran gcc version 7.5.0 (Ubuntu 7.5.0-3ubuntu1~18.04) 
      do i = 1,nq0i+1
        ini = nqbz*(i-1)
        do ix=1,nqbz
            write(6,"('hx0fp0 qbze q0i=',i8,3f10.4,2x,3f10.4)") ini+ix,qbze(:,ini+ix)!,q0i(:,i)
         enddo
         write(6,*)
      enddo

c$$$  do i = 1,nq0i+1
c$$$        ini = nqbz*(i-1)
c$$$        do ix=1,nqbz
c$$$            write(6,"('hx0fp0 qbze q0i=',i8,3f10.4,2x,3f10.4)") ini+ix,qbze(:,ini+ix)!,q0i(:,i)
c$$$         enddo
c$$$         write(6,*)
c$$$      enddo

!! Get space-group transformation information. See header of mptaouof. 
!! Here we use ngrpx=1 ==> "no symmetry operation in hx0fp0", c.f. hsfp0.sc.m.F case.
!! ngrpx=1 (no symmetry operation in hx0fp0), whereas we use ngrp in eibzmode=T.
      ngrpx = 1
      l2nl=2*(nl-1)
      allocate(symope(3,3))
      symope(1:3,1) = [1d0,0d0,0d0]
      symope(1:3,2) = [0d0,1d0,0d0]
      symope(1:3,3) = [0d0,0d0,1d0]
      call Mptauof_zmel(symope,ngrpx) !we set thing in m_zmel for matrix elemenets generator
!! ppbrd = radial integrals, cgr = rotated cg coeffecients (no rotatio here since nrgpx=1 for identity matrix)
c      call Rdpp(nl, ngrpx, nn, nclass, nspin, symope,qbas) !nxx, see m_rdpp
c      call Rdpp(ngrpx,symope)
!! Set nband, ngcmx, nqpmx and itq in m_zmel
      call Setitq()

!! Pointer to optimal product basis
      nblochpmx = nbloch + ngcmx
      allocate(ngveccB(3,ngcmx)) ! work arry
      iqxend = nqibz + nq0i
      write(6,*) ' nqibz nqibze=',nqibz,nqibze

!! Initialization of readEigen !readin m_hamindex
      call init_readeigen()!nband,mrece)!EVU EVD are read in init_readeigen
      call init_readeigen2()!(mrecb,nlmto,mrecg)
      if(verbose()>50) print *,'eeee exit of init_readeigen2'
!! Frequency      
      call Getfreq3(lqall,epsmode,realomega,imagomega,ua,mpi__root)
!! Write freq_r
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
      if(MPI__root) write(6,"(' nw=',i5)") nw
      nwp = nw+1
      
!! Get eigenvector corresponds to exp(iqr) (q is almost zero).
      if(epsmode) allocate(epsi(nw_i:nw,neps)) 

!! Tetrahedron initialization
      noccxv = maxocc2 (nspin,ef, nband, qbze,nqbze) ! maximum no. occupied valence states
      if(noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
      noccx  = noccxv + nctot
      nprecx = ndble  !We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
      if (MPI__root) then
        open(newunit=ifwd,file='WV.d')
        write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
        close(ifwd)
      endif
      allocate( zw(nblochpmx,nblochpmx) )
      nspinmx = nspin

!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,R(r))> 
!! This is general for rotated CG coefficient; but hx0fp0 mode is only for  ngrpx=1 (not rotated).
!! Compare usage in hsfp0 modes.
c      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,ngrpx,nspin))
c      call ppbafp_v2_zmel (ngrpx)!,nspin, !all inputs. This is in m_zmel
c     i   il,in,im,nlnm,         !w(i_mnl),
c     i   nl,nn,nclass,nlnmx,
c     i   mdimx,lx,nx,nxx,       !Bloch wave    
c     i   cgr, nl-1,             !rotated CG
c     i   ppbrd) !,              !radial integrals
      if(debug)write(6,*) ' end of ppbafp_v2'

!! Set iqxini !omitqbz means skip loopf for iq=1,nqibz
      if(omitqbz) then
        iqxini= nqibz + 1
      else
        iqxini= 1
      endif
!!
      if( chipm ) then
c        nmbas=natom
c        allocate(imbas_s(nmbas))
c        imbas_s(1:nmbas) = imbas(1:nmbas)
        allocate(aimbas(nmbas))
        aimbas(1:nmbas)   = abs(imbas(1:nmbas))
        allocate(jcoup(nw_i:nw,neps) )
        allocate( svec(nbloch,nmbas) )  !sep2006
        svec=0d0
        allocate( cvec(nbloch,nmbas),momsite(nmbas), mmnorm(nmbas)) !May2007
        cvec=0d0
        do imb=1,nmbas
          ibas= aimbas(imb)
          ifv = iopen ('MixSpin.'//charnum3(ibas),1,3,0)
          read(ifv,*) ibasx,lxx
          allocate(nxx_r(0:lxx))
          do i=0,lxx
            read(ifv,*) nxx_r(i)   !   write(6,"(2i5,d13.6)") nxx_r(i)
          enddo
          allocate(spinvec((lxx+1)**2,maxval(nxx_r)))
          allocate(consvec((lxx+1)**2,maxval(nxx_r)))
          spinvec=0d0
          do ilmx = 1, (lxx+1)**2
            lb = ll(ilmx )         !  write(6,*)' lb=',lb,lxx,ilmx
            do ixx = 1, nxx_r(lb)  !  write(6,*)' nn=',nn,nxx_r(lb)
              if(ilmx==1) then
                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2 ,consvec(ilmx,ixx)
              else
                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),dumm1,dumm2 ,consvec(ilmx,ixx)
              endif
!           write(6,"(2i5,d13.6)") ilmx, ixx, spinvec(ilmx,ixx)
            enddo
          enddo
!! Calculate ChiPM. So sign of omega should be correct.
          if(imb==1) then !determine spin direction with respect to ibas=aimbas(imb=1)
            spinmom=(chg1-chg2)
            schi=1d0
            if(spinmom<0d0) then
              schi  = -1d0   ! This affects to dpsion. Obtained results
                             ! should be the same in both mode.
            endif
          endif
!!  ReOrdering of spinvec in natom ordering...
          i=0
          if(ibas>1) i= sum(nblocha(1:ibas-1))
          do lb  = 0, lx (ibas)
            do nb  = 1, nx (lb,ibas)
              do mb  = -lb, lb
                i = i+1
                ilmx = lb**2+ lb+ mb +1
                svec(i,imb) = spinvec(ilmx,nb)
                cvec(i,imb) = consvec(ilmx,nb)
                write(6,"(' i lb mb svec svec**2=',3i4,2d13.5)") 
     &      i,lb,mb,svec(i,imb),svec(i,imb)**2
              enddo
            enddo
          enddo
          deallocate(nxx_r,spinvec,consvec)
          close(ifv)
          mmnorm (imb) = sqrt(sum(svec(:,imb)**2))
          momsite(imb) = chg1-chg2
c        write(6,"(' svecsum=',e23.15)") sum(svec(:,imb)**2
c        write(ifchipmlog,"(2e23.15,' ! mmom mmnorm')")momsite(imb),mmnorm(imb)
          write(6,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
        enddo
      endif
!!
      allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))
      
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
 !! EIBZ mode memo for nolfco (right?) 
 !! If eibzmode=T, it is efficient but can slightly break crystal symmetry.(how much?)
 !! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
 !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      eibzmode=eibz4x0() 
      if(eibzmode) then
 !! comment-out block inversion. Use iqxendx=iqxend because of full inversion
        call cputid(0)
        write(6,*)' ---goto eibzmode block ---'
        iqxendx=iqxend
        if(epsmode) iqxendx=iqxend
        allocate( nwgt(nqbz,iqxini:iqxendx), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     &   igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     &   eibzsym(ngrp,-1:1,iqxini:iqxendx))
 !! Check timereversal is required for symmetrization operation or not. 
 !! If tiii=timereversal=F is enforced,
 !!   the symmetrization procedure in x0kf_v4h becomes a little time-consuming.
        write(6,*)
        write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal() 
        if(MPI__root) iprintx=.true.
        call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,
     i   timereversal(),ginv,iprintx,
     o   nwgt,igx,igxt,eibzsym,tiii)
        write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
        call cputid(0)
 !! PBindex: index for product basis. 
 !! We will unify this system; still similar is used in ppbafp_v2.
        call PBindex(natom,lx,l2nl,nx)!all input. Returns requied index stored in arrays in m_pbindex.
        call cputid(0)
        call readqgcou() !no input. Read QGcou and store date into variables.
      else     !dummy allocation to overlaid -check bound !sep2014
        iqxendx=iqxend
        allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &   ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
        nwgt=1
      endif
!! Calculate x0(q,iw) and W == main loop 1001 for iq. 
 !! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
 !! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
 !! (or need to modify do 1001 loop).
 !! iq>nqibz for ixc=11 is not time-consuming (right???)
      call MPI__hx0fp0_rankdivider2(iqxini,iqxend)
      

!! llw, and llwI are for L(omega) for Q0P in PRB81,125102
      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
!! ======== Loop over iq ================================
      do 1001 iq = iqxini,iqxend  ! NOTE: q=(0,0,0) is omitted when iqxini=2
        if(.not. MPI__task(iq) ) cycle
        if(normalm) then
          ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
          ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
        endif
        call cputid (0)
        q = qibze(:,iq)
        call readqg('QGcou', q,  quu,ngc,ngveccB) 

 !! Caution : confusing point
 !!  ngc by QGcou is shown at the bottom of lqg4gw.
 !!  ngc read from PPOVL are given by rdata4gw---> ngc(iq>nqibz )=ngc for q=0
        if( iq==1 ) then ! *sanity check
          if(sum(q**2)>1d-10) then
            call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
          endif
        endif

!! Readin Coulomb matrix ====
        ngb = nbloch + ngc !ngb is readin from vcoul 25jan2006
        write(6,*)
        write(6,"('===== do 1001: iq q=',i7,3f9.4,' ========')")iq,q !qq
        write(6,"('  nbloch ngb ngc=',3i10)") nbloch,ngb,ngc

!! Readin diagonalized Coulomb interaction ===
 !! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
        if((.not.chipm)) then
          call Readvcoud(q,iq,ngb) ! Readin vcousq,zcousq !for the Coulomb matrix
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
          if(allocated(zzr)) deallocate(zzr)
          allocate(zzr(1,1)) !dummy
          zzr=0d0
        endif
      
        if(chipm .and. nolfco) then 
          nmbas_in = nmbas
        elseif(nolfco) then
          nmbas_in = 1
        else
          nmbas_in = ngb
        endif
        nmbas1 = nmbas_in
        nmbas2 = nmbas1

!! Open output files for epsilon mode =====
        if(epsmode) then 
          iqixc2 = iq- (nqibz+nq0ix)
          if((.not.chipm).and.nolfco) then
            allocate( x0mean(nw_i:nw,1,1) )
            x0mean=0d0
          endif  
          if((.not.chipm).and.wqt(iq-nqibz)==0d0) then
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nlfc.dat'
            ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
            write(ifepsdatnolfc,"(a)")' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
            if(.not.nolfco) then 
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              ifepsdat = iopen ( fileps,1,3,0)
              write(ifepsdat,"(a)") ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
            endif
          endif
!! zzr is only for chipm.and.nolfco mode
          if(chipm .and.nolfco) then
            allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,nmbas))
            x0mean=0d0
            zzr   =0d0
            zzr(1:nbloch,1:nmbas) = svec(1:nbloch,1:nmbas)
          endif
!! ... Open ChiPM* files for \Chi_+-
          if(chipm.and. wqt(iq-nqibz)==0d0) then
            ifchipmn_mat=iopen('ChiPM'//charnum4(iqixc2)//'.nlfc.mat',1,3,0)
            write(ifchipmn_mat,"(255i5)") nmbas
            write(ifchipmn_mat,"(255i5)") aimbas(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)") momsite(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)")  mmnorm(1:nmbas)
            write(ifchipmn_mat,"( ' Here was eiqrm: If needed, need to fix hx0fp0')") 
            if(.not.nolfco) then
              ifchipm_fmat=iopen('ChiPM'//charnum4(iqixc2)//'.fmat',0,3,0)
              write(ifchipm_fmat) nbloch, natom,nmbas, iqxini,iqxend, nw_i,nw
              write(ifchipm_fmat) aimbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
              write(ifchipm_fmat) nblocha(1:natom),svec(1:nbloch,1:nmbas)
              write(ifchipm_fmat) zzr0(1:nbloch) !zzr(1:nbloch,1)
            endif
          endif
        endif

        if(epsmode.and.nolfco) then 
          continue
        else
          write(6,*) "rcxq alloc ngb nwhis npm ---",ngb,nwhis,npm
          allocate( rcxq(ngb,ngb,nwhis,npm) )
        endif

c       print *,' iiiiiiiiii of setppovlz ...',nmbas1,nolfco
!! Set ppovlz: zmelt conversion. Non-orthogonality and Enu(coulomb-diagonal) related.
!!         ppovlz is used in get_zmelt2 in m_zmel (called in x0kf_v4h).
        if(chipm.and.nolfco) then
          call setppovlz_chipm(zzr,ngb,nmbas1)
        elseif(nolfco .and. nmbas1==1) then !for <e^iqr|x0|e^iqr>
          call Setppovlz(q,ngc,ngb,matz=.true.)
        else                     !may2013  this removes O^-1 factor from zmelt
          call Setppovlz(q,ngc,ngb,matz=.not.eibz4x0())
       endif
c       print *,' end of setppovlz ...'

!! rcxq: imaginary part after x0kf_v4h and symmetrization. 
!! zxq ans zxqi are the main output after Hilbert transformation
        if(nolfco) then
          if(allocated(rcxq)) deallocate(rcxq)
          if(allocated(zxq) ) deallocate(zxq)
          if(allocated(zxqi) ) deallocate(zxqi)
          allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
          allocate( zxq (nmbas1,nmbas2,nw_i:nw), zxqi (nmbas1,nmbas2,niw))
        else
          allocate( zw0(ngb,ngb), zxq (ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
        endif
        zxq=0d0;  zxqi=0d0;  rcxq = 0d0

!! === loop over spin=== =========================================
        do 1003 is = 1,nspinmx
          write(6,"(' ##### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ##### ')")iq, is, nqibz + nq0i,nspin
          if(debug) write(6,*)' niw nw=',niw,nw
!! ==== spin chi_charge or chi_+- ====
          isf=is
          if(chipm) then
            write(6,*)" chi_+- mode nolfc=",nolfco
            if(is==1) isf=2
            if(is==2) isf=1
            rcxq=0d0
          endif

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
          do kx = 1, nqbz
            call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) ) 
            call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
          enddo
          call gettetwt(q,iq,is,isf,nwgt(:,iq),ekxx1,ekxx2,nband,eibzmode) 
c          call gettetwt(q,iq,is,isf,nwgt(:,iq),frhis,nwhis,npm,
c     i     qbas,ginv, ef, nqibz, nband,ekxx1,ekxx2, nctot,ecore,
c     i     nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
c     i     nbmx,ebmx,mtet,eibzmode) !nov2016

c$$$ !!! test code to print tetrahedron weight !!!
c$$$        jpm=1
c$$$        do k=1,nqbz
c$$$c           nkqmin= 999999
c$$$c           do ibib = 1, nbnb(k,jpm)
c$$$c              nkqmin = min(n2b(ibib,k,jpm),nkqmin)
c$$$c           enddo
c$$$           do ibib = 1, nbnb(k,jpm) !---  ibib loop
c$$$c             print *,' k ibib=',k,ibib
c$$$c              it =  n1b(ibib,k,jpm)    !valence
c$$$c              itp = n2b(ibib,k,jpm) - itps + 1 !val
c$$$c              if(n1b(ibib,k,jpm)==n2b(ibib,k,jpm)) then
c$$$              do iw = ihw(ibib,k,jpm),ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 !iiww=iw+ihw(ibib,k)-1
c$$$c                 if(iw<20) then
c$$$                 imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
c$$$c                 write(*,'("eeee ",4i5,2x,d13.5,x,d13.5)') k, n1b(ibib,k,jpm),n2b(ibib,k,jpm), iw, imagweight
c$$$c                 endif
c$$$              enddo             ! iw
c$$$c              endif
c$$$           enddo
c$$$        enddo
c$$$        stop 'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv'


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
!! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 into rcxq ==
c          call x0kf_v4hz(q, is,isf,  iq,ngb,ngc, 
c     o    rcxq,     !  rcxq is the accumulating variable for spins 
c     i    nmbas_in, 
c     i    eibzmode,nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq), eibzsym(:,:,iq), crpa)
!! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 into rcxq ==
          call cputid(0)
c          if(npm==1) then
c            ncc=0
c          else
c            ncc=nctot
c          endif
          call x0kf_v4hz(!npm,ncc,   
c     i        ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
c     i        n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf, !symmetrize, !
     i              qbas,ginv,  qbz,wbz, 
     d            nlmto,nqbz,nctot, !noccx,noccxv,
     d            nbloch,  nwhis, !nlnmx,mdimx,
     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
     i    nqbze,nband,nqibz, 
     o    rcxq,   !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,nmbas_in, zcousq,   !ppovl,nmbas1,nmbas2, is removed ppovlz, 
     i    chipm,eibzmode, !z1offd,!for nolfco Add nmbas Sep2006
     i    nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),ngrp, eibzsym(:,:,iq),crpa)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          
!! == Symmetrizer for crystal symmetry (and also for spin)
!!    Symmetrize and convert to Enu basis by dconjg(tranpsoce(zcousq)*rcxq8zcousq if eibzmode
          if(is==nspinmx.or.chipm) then !Apr2015. TK think ".or.chipm" is required for chipm mode
                                        !Because rcxq is calculated for each is, symmetrized and its contribution
                                        !is added to zxq in dpsion5.
            call x0kf_v4hz_symmetrize(npm, !ncc,
     i       q,
     i       nspin,is,isf, !symmetrize, !
     i       qbas,ginv,  !qbz,wbz,
     d       nbloch,  nwhis,    ! nlnmx,mdimx,
     i       iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
     i       nqbze,nband,nqibz,
     o       rcxq,              !  rcxq is the accumulating variable for spins
     i       nolfco,zzr,nmbas_in, zcousq, !ppovl,nmbas1,nmbas2, is removed ppovlz, 
     i       chipm,eibzmode,    !z1offd,!for nolfco Add nmbas Sep2006
     i       ngrp, eibzsym(:,:,iq))
          endif
          call tetdeallocate() !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
c          iecut=1
          if(debug) write(6,"(a)") ' --- goto dpsion5 --- '
          if(is==nspinmx.or.chipm) then
            write(6,"('  nmbas1,nmbas2=',2i10)") nmbas1,nmbas2
            call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i       rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o       zxq, zxqi,
     i       chipm, schi,is,  ecut,ecuts)!ecut(iecut),ecuts(iecut))
            if(nolfco.and.epsmode) then
              do iw=nw_i,nw
                x0mean(iw,:,:)=zxq(:,:,iw)
              enddo
            endif
            write(6,*)' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
          endif
 1003   continue  !end of spin loop =====
        if(allocated(rcxq) ) deallocate(rcxq)


!! ===  RealOmega ====================================
        if (realomega) then
          if(epsmode.and.nolfco.and.(.not.chipm)) then 
            if (nspin==1) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
            if (nspin==1) zxq = 2d0*zxq       !if paramagnetic, multiply x0 by 2
          else
            if (nspin == 1) zxq = 2d0*zxq   !if paramagnetic, multiply x0 by 2
          endif

          if(epsmode) then
            if(nolfco) then
              ttt='without LFC'
            else
              ttt='with LFC'
            endif
            if(chipm) then
              write(6,*) '--- chi0_{+-}}^{-1}      --- '//ttt
            else
              write(6,*) '--- dielectric constant --- '//ttt
              write(6, *)" trace check for W-V"
            endif
          endif

!! prepare for iq0.
          iq0 = iq - nqibz
          if(allocated(epstilde)) deallocate(epstilde,epstinv)
          allocate(epstilde(ngb,ngb),epstinv(ngb,ngb))

!! === iw loop for real axiw ===
          do 1015 iw  = nw_i,nw  !Feb2006. Before it was 1:nwp (nwp=nw+1).
            !  So freq_r(iw-1) is shifted to freq_r(iw).
            frr= dsign(freq_r(abs(iw)),dble(iw))
            if(.not.epsmode) then  
!! === wcf: W= (1-v zxq)^{-1} v ===
c              if(newaniso2.and.iq<=nqibz) then !for mmmw
              if(iq<=nqibz) then !for mmmw
                if(iq==1) then
                  ix=1
                  zw0(:,1)=0d0
                  zw0(1,:)=0d0
                else
                  ix=0
                endif
!!  Eqs.(37),(38) in PRB81 125102
                do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                    if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cmmm direct inversion vs. block inversion
c$$$                if(iq>nqibz) then
c$$$c direct inversion
c$$$                ix=0
c$$$                do igb1=ix+1,ngb
c$$$                do igb2=ix+1,ngb
c$$$                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$                enddo
c$$$                enddo
c$$$                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$                do igb1=1+ix,ngb
c$$$                do igb2=1+ix,ngb
c$$$                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$                enddo
c$$$                enddo
c$$$c block inversion
c$$$                ix=1
c$$$                do igb1=ix+1,ngb
c$$$                do igb2=ix+1,ngb
c$$$                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$                enddo
c$$$                enddo
c$$$                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$                absq=sqrt(sum(q**2*tpioa**2))
c$$$                sk(  1:ngb)= zxq(1,1:ngb,iw)
c$$$                sks( 1:ngb)= zxq(1:ngb,1,iw)
c$$$                w_k(1) =0d0
c$$$                w_ks(1)=0d0
c$$$                w_k( 2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(vcousq(1)*sk(2:ngb)*vcousq(2:ngb),epstinv(2:ngb,2:ngb))
c$$$                w_ks(2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(epstinv(2:ngb,2:ngb),vcousq(1)*sks(2:ngb)*vcousq(2:ngb))
c$$$                llw(iw,iq0)=
c$$$     &             1d0
c$$$     &            -vcousq(1)*sk(1)*vcousq(1) ! sk(1,1,iw)=sks(1,1,iw)=H of Eq.(40).
c$$$     &            -vcousq(1)*vcousq(1)* sum( vcousq(2:ngb)*sk(2:ngb) * matmul(epstinv(2:ngb,2:ngb),sks(2:ngb)*vcousq(2:ngb)))
c$$$                endif
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                do igb1=1+ix,ngb
                  do igb2=1+ix,ngb
                    zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                    if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                  enddo
                enddo
              endif

!! for iq>nqibz Full inversion to calculalte eps with LFC.
              if(iq>nqibz) then
                ix=0
                vcou1 = fourpi/sum(q**2*tpioa**2) !--> vcousq(1)**2! !fourpi/sum(q**2*tpioa**2-eee)
                do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    if(igb1==1.and.igb2==1) then
                      epstilde(igb1,igb2)= 1d0 - vcou1*zxq(1,1,iw)
                      cycle
                    endif
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                    if(igb1==igb2) then
                      epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                    endif   
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                llw(iw,iq0)= 1d0/epstinv(1,1)
                write(6,"('iq iw_real eps(withLFC) eps(woLFC) ',2i5,d13.6,x,d13.6,2x,d13.6,x,d13.6)")
     &           iq,iw,llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
              endif
              if(iq>nqibz) then
                continue
              else
                zw(1:ngb,1:ngb) = zw0
                write(ifrcw, rec= iw-nw_i+1) zw   !  WP = vsc-v
                call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
              endif
!! epsmode
            elseif(epsmode.and.(.not.chipm)) then 
               if(debug)write(6,*) 'xxx2 epsmode iq,iw=',iq,iw
!! there is difference of two vcmean below since we use (sligthy) screened Coulomb (screenfac() in switch.F)
!!   NOTE that we use vcoul with screening (screenfac() is used in hvccfp0.F
                vcmean=vcousq(1)**2 !fourpi/sum(q**2*tpioa**2) !aug2012
                epsi(iw,iqixc2)= 1d0/(1d0 - vcmean*zxq(1,1,iw))
                write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2e23.15,3x, 2e23.15,
     &          " vcmean x0mean =", 2e23.15,3x, 2e23.15)') iqixc2,iw,2*frr,
     &          1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean, zxq(1,1,iw)!x0mean(iw,1,1)
                write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
     &          q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                if(.not.nolfco) then
                  ix=0
                  do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    if(igb1==1.and.igb2==1) then
                      epstilde(igb1,igb2)= -vcmean*zxq(igb1,igb2,iw) !aug2012
                    else   
                      epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                    endif  
                    if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                  enddo
                  enddo
                  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                  epsi(iw,iqixc2)= epstinv(1,1)
                  write(6,'( " iq iw omega eps epsi  wLFC="
     &           ,2i6,f8.3,2e23.15,3x, 2e23.15)') 
     &            iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                  write(6,*)
                  write(ifepsdat,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
     &            q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                endif
!! ChiPM mode
            elseif(epsmode.and.chipm) then
              allocate( x0meanx(nmbas,nmbas) )
              if(nolfco) then  ! ChiPM mode without LFC
c$$$c! --- three lines below may work for test purpose for legas. But not sure.
c$$$c       vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
c$$$c       write(ifchipmn,'(3f12.8,2x,f8.5,2x,2e23.15)')
c$$$c     & q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw,1,1)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
                x0meanx = x0mean(iw,:,:)/2d0 !in Ry unit.
              else
!! ChiPM mode with LFC... NoLFC part
                zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)/2d0  ! in Ry.
                do imb1=1,nmbas
                  do imb2=1,nmbas
                    x0meanx(imb1,imb2)= 
     &              sum( svec(1:nbloch,imb1)*
     &              matmul(zxq(1:nbloch,1:nbloch,iw),svec(1:nbloch,imb2))) 
                  enddo
                enddo
!  x0meanx= <m|chi^+-(\omega)|m>/<m|m>**2
              endif
              do imb1=1,nmbas
                do imb2=1,nmbas
                  x0meanx(imb1,imb2) = 
     &            x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
                enddo
              enddo
              write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,255e23.15)')q, 2*schi*frr, x0meanx(:,:)
              if(.not.nolfco) write(ifchipm_fmat) q, 2*schi*frr, zxq(1:nbloch,1:nbloch,iw)
              deallocate(x0meanx)
            endif
 1015     continue
          if( allocated(zzr)   ) deallocate(zzr)
          if( allocated(x0mean)) deallocate(x0mean)
          if( allocated(gbvec) ) deallocate(gbvec)
        endif
!! ... Close files for epsmode
        if(epsmode) then 
          if(chipm) then
            ifchipmn_mat=iclose('ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
            if(.not.nolfco) then
              ifchipm_fmat=iclose( 'ChiPM'//charnum4(iqixc2)//'.fmat')
            endif
          else
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
            ifepsdatnolfc = iclose( filepsnolfc)
            if(.not.nolfco) then
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              ifepsdat  = iclose(fileps)
            endif
          endif
        endif
!! RealOmega end ==============================



!! ImagOmega start ============================
        if (imagomega) then
          write(6,*)' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi   ! if paramagnetic, multiply x0 by 2
!! === iw loop for imag axiw ===
          do 1016 iw  = 1,niw
            if( iq<=nqibz ) then
!!  Eqs.(37),(38) in PRB81 125102
              if(iq==1) then
                ix=1
                zw0(:,1)=0d0
                zw0(1,:)=0d0
              else
                ix=0
              endif
!!  Eqs.(37),(38) in PRB81 125102
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv=epstilde
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                enddo
              enddo
            endif
            if(iq>nqibz) then
!! Full inversion to calculalte eps with LFC.
              ix=0
              vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  if(igb1==1.and.igb2==1) then
                    epstilde(igb1,igb2)= 1d0 - vcou1*zxqi(1,1,iw)
                    cycle
                  endif
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) then
                    epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                  endif   
                enddo
              enddo
              epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              llwI(iw,iq0)= 1d0/epstinv(1,1)
              write(6,*) 'iq iw_img  eps(withLFC) eps(woLFC)',iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
            endif
            if(iq>nqibz) then
              continue
            else
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
              write(ifrcwi, rec=iw)  zw    !  WP = vsc-v
              call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
            endif
c$$$!! --- Miyake's total energy branch !Nov2004. not maintained now... need to fix this maybe(2012takao) -----------
c$$$            if(.false.) then
c$$$              if (debug) write(6,*)'ip,ix=',iq,iw,'  niw=',niw
c$$$              call getwk(iq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)
c$$$              call ecorq2 (vcoul, zw0, ngb, iq,iw,ieceig,
c$$$     o        erpaqw, trpvqw, trlogqw)
c$$$C --- Integration along imaginary axis.
c$$$! omit k and basis index for simplicity
c$$$! wint = -(i/4pi) < [w'=-inf,inf] Q(w') >
c$$$!
c$$$! When w' ==> iw', w' is now real,
c$$$!   wint =  (1/2pi) < [w'=0,inf] Q(iw') >
c$$$!
c$$$! transform: x = 1/(1+w')
c$$$! this leads to a denser mesh in w' around 0 for equal mesh x
c$$$! which is desirable since Q is peaked around w'=0
c$$$!    wint =  (1/2pi) < [x=0,1] Q(iw') / x^2 >
c$$$              faca   =  wk4ec* wiw(iw)
c$$$              trpv(iecut)    = trpv(iecut)  + faca* trpvqw
c$$$              trlog(iecut)   = trlog(iecut) + faca* trlogqw
c$$$              totexc(iecut)  = totexc(iecut)+ faca* erpaqw !  = trpv+ trlog
c$$$c            ecqw(iq,iw) = erpaqw
c$$$              if(iw==1) then
c$$$                write(ieclog,*)
c$$$              endif
c$$$              if(iw==1.and.iq==iqxini) then
c$$$                write(ieclog,
c$$$     &         "('   iq   iw   omega/i(Ry)        IntWgt',
c$$$     &         '    trpvqw(eV)      ecqw(eV)   ecqw*IntWgt',  
c$$$     &         ' :  ecut   ecuts')")
c$$$              endif
c$$$              write(ieclog,"( 2i5,3f14.6,3f14.6,2f8.3)") 
c$$$     &       iq,iw, 2d0*freq_i(iw), faca, trpvqw*hartree, erpaqw*hartree, 
c$$$     &       faca*erpaqw*hartree, ecut(iecut),ecuts(iecut)
c$$$              close(ieclog)
c$$$              open(ieclog,file="ecorr.chk",access='append')
c$$$            endif
 1016     continue
            deallocate(epstinv)
            if(allocated(epstilde)) deallocate(epstilde)
            if(allocated(epstilde)) deallocate(epstilde)
        endif
!! ImagOmega end =================
        if(allocated(vcoul)) deallocate(vcoul)
        if(allocated(zw0)) deallocate(zw0)
        if(allocated(zxq )) deallocate(zxq)
        if(allocated(zxqi)) deallocate(zxqi)
        if (normalm) then
          ifrcwi = iclose('WVI.'//charnum5(iq))
          ifrcw  = iclose('WVR.'//charnum5(iq))
        endif
 1001 continue 
!! =================== end of loop 1001 for q point ========================
      call MPI__barrier()

!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
      if(MPI__size/=1) then
        do iq=nqibz+1,iqxend
          iq0 = iq - nqibz
          if(MPI__ranktab(iq)/=0) then !jan2012
            if(MPI__ranktab(iq) == MPI__rank) then
              dest=0
              call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
              call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
            elseif(MPI__root) then
              src=MPI__ranktab(iq)
              call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
              call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
            endif
          endif
        enddo  
      endif

!! == W(0) divergent part and W(0) non-analytic constant part.==
!!   Note that this is only for q=0 -->iq=1
      if(normalm.and.MPI__rank==0) then
!! get w0 and w0i (diagonal element at Gamma point
!! This return w0, and w0i
        call w0w0i(llw,llwI,nw_i,nw,nq0i,niw,q0i)
!! === w0,w0i are stored to zw for q=0 ===
!! === w_ks*wk are stored to zw for iq >nqibz ===
        do iq = 1,1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
          q = qibze(:,iq)
          do ircw=1,2
            if    (ircw==1) then;  nini=nw_i;   nend=nw;
              ifrcwx = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
            elseif(ircw==2) then;  nini=1;     nend=niw;
              ifrcwx = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
            endif
            do iw=nini,nend
              read(ifrcwx, rec= iw-nini+1) zw !(1:ngb,1:ngb)
              if(ircw==1) zw(1,1) = w0(iw)
              if(ircw==2) zw(1,1) = w0i(iw)
c              if(ircw==1) then
c                 write(6,"('ffffrrr:', f13.6,2x,f13.6,x,f13.6)") hartree*freq_r(iw),w0(iw)
c              endif
              write(ifrcwx,rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
            enddo
            if    (ircw==1) then
              ifrcwx = iclose('WVR.'//charnum5(iq))
            elseif(ircw==2) then
              ifrcwx = iclose('WVI.'//charnum5(iq))
            endif
          enddo
        enddo
      endif

c$$$!! --- legas mode is not working now. Need fixing... voltot ntot are not given.
c$$$      if(epsmode.and.legas) then
c$$$        call rx( ' LEGAS mode is not maintained well. Need some fixing.')
c$$$        voltot=0d0
c$$$        ntot=0d0
c$$$        write(6,*)' Find LEGAS. legas =',legas
c$$$        iflegas = 2101
c$$$        open (iflegas,file='LEGAS')
c$$$        read(iflegas,*)rs
c$$$        close(iflegas)
c$$$        alpha  = (9*pi/4d0)**(1d0/3d0)
c$$$        qfermi = alpha/rs
c$$$        efx  = qfermi**2
c$$$        valn = efx**1.5d0*voltot/3d0/pi**2
c$$$        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
c$$$        write (6,*)'     Exact Fermi momentum  qf  =', qfermi
c$$$        write (6,*)'     Exact Fermi energy    Ef  =', efx
c$$$        do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
c$$$          if(iq<=nqibz) cycle
c$$$          write(6,*)' iq=',iq
c$$$          iqixc2 = iq- (nqibz+nq0ix)
c$$$          filele ='EPSEG'//charnum4(iqixc2)//'.dat'
c$$$          ife = iopen ( filele,1,3,0)
c$$$          write(ife,"(a)")
c$$$     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
c$$$          q = qibze(:,iq)
c$$$          qt= sqrt(sum(qibze(1:,iq)**2))*2d0*pi/alat
c$$$          qs= qt/qfermi
c$$$          write(6,"(' qs qfermi=',2d13.5)"    ) qs,qfermi
c$$$          write(6,"(' q-q^2/2 q+q^2=',2d13.5)") qs-qs**2/2d0,qs+qs**2/2d0
c$$$          do iw  = nw_i,nw
c$$$            ww  = freq_r(iw)
c$$$            muu = ww/qfermi**2
c$$$            if(     qs<2d0 .and. muu < qs-qs**2/2d0) then
c$$$              x0mx= -img*qfermi/(4*pi*qs)*2*muu
c$$$            elseif( qs<2d0 .and. muu < qs+qs**2/2d0) then
c$$$              x0mx= -img*qfermi/(4*pi*qs)*( 1d0-(muu/qs-.5d0*qs)**2 )
c$$$            else
c$$$              x0mx=0d0
c$$$            endif
c$$$            vcmmmm= 4*pi/qt**2
c$$$            epsi(iw,iqixc2) = 1d0/(1- vcmmmm * x0mx)
c$$$c            epsi(iw,iqixc2) = 1d0/(1- vcmmm(iq) * x0meanx)
c$$$            write(ife,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
c$$$     &        q, 2*ww,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$          enddo
c$$$        enddo
c$$$        write(6,*)' ----------legas end--------'
c$$$      endif

c$$$!! Write TEECOR ecorr_on mode
c$$$      if(imagomega.and.ecorr_on>0) then
c$$$        hartree=2d0*rydberg()
c$$$        ifcor   = iopen('TEECORR2',1,-1,0) ! output files
c$$$        do iecut=1,necut
c$$$          write(6,"( ' RPA Ec =' 3f23.15,'   ecut ecuts (Ry)=',2d12.4)") 
c$$$     &   totexc(iecut)*hartree,trpv(iecut)*hartree, trlog(iecut)*hartree 
c$$$     &    ,ecut(iecut),ecuts(iecut)
c$$$          write(ifcor,*) '============================'
c$$$          write(ifcor,*) 'Correlation energy Erpa (eV)'
c$$$          write(ifcor,*) '============================'
c$$$          write(ifcor,*)' ### '
c$$$          write(ifcor,"(5e23.15)") 
c$$$     &     totexc(iecut)*hartree,trpv(iecut)*hartree,trlog(iecut)*hartree
c$$$     &    ,ecut(iecut),ecuts(iecut)
c$$$        enddo
c$$$!! output ecqw !    write(ifcor,*)'### ecqw(q,w) ###'
c$$$        write(ifcor,*)' nqibz =',nqibz
c$$$        write(ifcor,*)' nq0i  =',nq0i
c$$$        write(ifcor,*)' niw   =',niw
c$$$        write(ifcor,*)' --- See details of Ec in ecor.chk ---'
c$$$C... Write electron gas correlation energy
c$$$c$$$        legas = .false.
c$$$c$$$        INQUIRE (FILE = 'LEGAS', EXIST = legas)
c$$$c$$$        if(legas) then !!! test for electron gas case.
c$$$c$$$          call rx( ' LEGAS mode is not maintained well. Need some fixing.')
c$$$c$$$          voltot=0d0
c$$$c$$$          ntot=0d0
c$$$c$$$          write(6,*)' find LEGAS. legas =',legas
c$$$c$$$          iflegas = 2101
c$$$c$$$          open (iflegas,file='LEGAS')
c$$$c$$$          read(iflegas,*)rs
c$$$c$$$          close(iflegas)
c$$$c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
c$$$c$$$          qfermi = alpha/rs
c$$$c$$$          efx  = qfermi**2
c$$$c$$$          valn = efx**1.5d0*voltot/3d0/pi**2
c$$$c$$$          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
c$$$c$$$          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
c$$$c$$$          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
c$$$c$$$          if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
c$$$c$$$          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
c$$$c$$$          qfermi= dsqrt(efz)
c$$$c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
c$$$c$$$          rs    = alpha/qfermi
c$$$c$$$          write (ifcor,*)' --- electron gas ---'
c$$$c$$$          write (ifcor,*)' density parameter rs= ', rs
c$$$c$$$          write (ifcor,*)' kf= ',qfermi
c$$$c$$$          write (ifcor,*)' ### Barth-Hedin formula'
c$$$c$$$          ecelgas = eclda_bh(rs) * hartree * ntot
c$$$c$$$          write (ifcor,*)ecelgas
c$$$c$$$          write (ifcor,*)' ### Perdew-Zunger formula'
c$$$c$$$          ecelgas = eclda_pz(rs) * hartree * ntot
c$$$c$$$          write (ifcor,*)ecelgas
c$$$c$$$          write (ifcor,*)' ### Gell-Mann and Brueckner formula'
c$$$c$$$          ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
c$$$c$$$          write (ifcor,*)ecelgas
c$$$c$$$        endif
c$$$      endif
      call cputid(0)
      call MPI__Finalize
      if(ixc==11)   call rx0( ' OK! hx0fp0 mode=11     read <Q0P> normal sergeyv')
      if(ixc==111)  call rx0( ' OK! hx0fp0 mode=111    normal sergeyv')
      if(ixc==10011)call rx0( ' OK! hx0fp0 mode=10011  crpa normal sergeyv')
      if(ixc==12)   call rx0( ' OK! hx0fp0 mode=12  Ecor sergeyv mode')
      if(ixc==101)  call rx0( ' OK! hx0fp0 mode=101 Ecor ')
      if(ixc==202)  call rx0( ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC')
      if(ixc==203)  call rx0( ' OK! hx0fp0 mode=203 sergeyv eps LFC ')
      if(ixc==222)  call rx0( ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv')
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

