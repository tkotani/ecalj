!!  Calculate x0, \epsilon, spin susceptibility. 
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_lmf_chipm mode works).
      program hx0fp0
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: readeval,init_readeigen,init_readeigen2
      use m_read_bzdata,only: read_bzdata,
     & ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,qbasmc,
     & dq_,qbz,wbz,qibz,wibz,qbzw,
     & idtetf,ib1bz,idteti,
     & nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     & nclass,natom,nspin,nl,nn, ngrp,
     & nlmto,nlnmx, nctot,niw, !nw_input=>nw,
     & alat, diw,dw,delta,deltaw,esmr,symgrp,clabl,iclass,
     & invg, il,in,im,nlnm, 
     & plat, pos,ecore, symgg 
      use keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
      use m_readqgcou,only: readqgcou
      use m_mpi,only: MPI__hx0fp0_rankdivider2,MPI__task,MPI__Initialize,MPI__Finalize,MPI__root,
     & MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size,
     & MPI__ranktab,MPI__consoleout,MPI__barrier
!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: rdpp, !NOTE: "call rdpp" generate following data.
     & nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr
!! Generate matrix element for "call get_zmelt".
      use m_zmel,only:       !NOTE: these data set are stored in this module, and used when 
     & nband,itq,ngcmx,ngpmx, ppovlz,
     & ppbir,shtvg, miat,tiat , ntq
!! frequency
      use m_freq,only: getfreq, !NOTE: call getfreq generate following data.
     & frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw !output of getfreq
!! tetwt
      use m_tetwt,only: gettetwt, !followings are output of 'L871:call gettetwt')
     &  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
      use m_lldata,only: ll
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
     &           , ngvecp(:,:), ngvecc(:,:), !,ngpn(:),ngcni(:),iqib(:),
     &   ifppb(:)   !ongveccBr(:,:,:),nx(:,:),nblocha(:),lx(:),
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &  zw(:,:),zw0(:,:),
     &  zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: eqt(:), !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
     &  ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), 
     &  qbze(:,:),qibze(:,:)  !,ecore(:,:)  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:)
c     & ,rcxqmean(:,:,:,:),rcxqmeanc(:,:,:,:) !now rcxqmean is treated as a case of rcxq(nmbas,nmbas)

!  tetrahedron method
      logical :: tetra !,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
      ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer,allocatable :: 
     &         nbnbtt(:,:),noccxvv(:) !n1b(:,:,:),n2b(:,:,:),nbnb(:,:),
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)!,ekxx3(:,:) !,qbzw(:,: )
      real(8) ::qbzx(3),anfvec(3)
      logical :: debug
      integer,allocatable:: ibasf(:)
      real(8),allocatable :: transaf(:,:)
      logical :: realomega, imagomega
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

      logical::imagonly=.false.,realonly=.false. !,readgwinput
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
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
c     &   ,ngpmx !,  ifchipmlog

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      integer::nwin, incwfin,  verbose
      real(8)::ef

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
     &  ,epsmode=.false.,normalm=.false., eiqr=.false.
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
      complex(8),allocatable:: llw(:,:), llwI(:,:),w0(:),w0i(:),aaamat(:,:)
      real(8),allocatable:: dmlx(:,:),epinvq0i(:,:),epinv(:,:,:),epinvq0i_m1(:,:),wklm(:),qeibz(:,:,:)
      integer:: lxklm,nlxklm,ifidmlx,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
      complex(8):: vc1vc2
      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)

      real(8),allocatable:: aik(:,:,:,:)
      integer,allocatable:: aiktimer(:,:)
      integer:: l2nl
      logical:: eibz4x0,tiii,iprintx,symmetrize,eibzmode
      real(8):: qread(3)

      character(128):: vcoudfile
      integer:: src,dest
      logical:: crpa
      integer,allocatable :: iclasst(:), invgx(:)
      integer:: ificlass,ifile_handle,ifief
      complex(8),allocatable:: ppovl_(:,:)

!! -------------------------------------------------------------------
      call MPI__Initialize()
      call MPI__consoleout('hx0fp0')
      call cputid (0)
      hartree  = 2d0*rydberg()
      pi       = 4d0*datan(1d0)
      fourpi   = 4d0*pi
      sqfourpi = sqrt(fourpi)
!! computational mode select
c takao keeps only the Sergey mode.
      write(6,"(a)") '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'
      write(6,"(a)") ' #1:run mode'
      write(6,"(a)") '    11  : normal    Sergey'
      write(6,"(a)") '    202 : epsNoLFC  Sergey'
      write(6,"(a)") '    203 : eps       Sergey'
      write(6,"(a)") '    222 : chi^+- NoLFC Sergey'
      write(6,"(a)") '    223 : chi^+- Sergey'
      write(6,"(a)") '    12  : total energy Miyake Sergey'
      write(6,"(a)") '    -9999: just show version num'
c      write(6,"(a)") ' #2=iqxini   #3=iqxend' '
c      write(6,"(a)") '  10222 : <e^{iqr}|chi^+-|e^{iqr}> NoLFC'
      write(6,"(a)")  '-------------------------------------------------------'
      if( MPI__root ) then
         read(5,*) ixc
      endif
      call MPI__Broadcast(ixc)
      call headver('hx0fp0',ixc)
      call cputid(0)
      crpa=.false.
      if(ixc<=6.or.ixc==22.or.ixc==23) then
        write(6,*)'these modes are removed now'
        call rx( 'these modes are not supported')
! Sergey (Hilbert-transformation) modes
      elseif(ixc==11) then; write(6,*) " OK ixc=11 normal mode "
        normalm=.true.
      elseif(ixc==10011) then; write(6,*) " OK ixc=10011 crpa mode "
        normalm=.true.
        crpa=.true.
!     -- eps mode NoLFC
      elseif(ixc==202) then
        write(6,*) " OK ixc=202  sergey's eps mode Only NoLFC "
        realonly=.true.
c        iepsmode=202
        omitqbz=.true.
!     -- eps mode with LFC
      elseif(ixc==203) then
        write(6,*) " OK ixc=203 sergey's eps mode with LFC "
        realonly=.true.
c        iepsmode=203
        omitqbz=.true.
! Total energy modes
      elseif(ixc==12) then
        write(6,*) " ixc=12 Miyake's total energy Sergey--->need to fix this mode"
        call rx( " ixc=12 Miyake's total energy Sergey--->need to fix this mode")
        imagonly=.true.
        ecorr_on=901
!     -- chipm mode NoLFC
      elseif(ixc==222) then
        write(6,*) " OK ixc=222    chipm sergey's "
        realonly=.true.
        omitqbz=.true.
        eiqr =.false. ! .true. aug2012
!     -- chipm mode NoLFC
c      elseif(ixc==10222) then
c        write(6,*) " OK ixc=10222  <q|chipm_0|q> sergey"
c        sergeyv=.true.
c        realonly=.true.
c        omitqbz=.true.
c        eiqr =.true.
!     -- eps mode with LFC
      elseif(ixc==223) then
        write(6,*) " ixc=223  chipm with LFC sergey's -->commented out not. need to fix this mode if necessary."
        call rx( " ixc=223  chipm with LFC sergey's -->commented out not. need to fix this mode if necessary.")
        realonly=.true.
        omitqbz=.true.
        eiqr =.true.
      else
        call rx( ' hx0fp0: mode ixc is not appropriate')
      endif

      if(ixc==202.or.ixc==203.or.ixc==222.or.ixc==223) then
        epsmode = .true.
        if(mod(ixc,200)==22.or.mod(ixc,200)==23) chipm =.true.
        if(mod(ixc,10)==2)                       nolfco=.true.
      endif

C ... files for RPA correlation energy mode.
      if(ecorr_on > 0) then
        ieclog = 8155
        if(ecorr_on==901) then
          ieceig=8156
          open(ieceig,file='rpa_eigen.chk')
          close(ieceig,status='delete')
        endif
        open(ieclog, file='ecorr.chk')
      endif
!! ====newaniso2====
c$$$      newaniso2=.false.
c$$$      if(newaniso()) then
c$$$        newaniso2=.true.
c$$$      endif

!! naraga says this cause a stop in ifort --->why???
c      write(6,*)'Timereversal=',Timereversal()

!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
      call read_BZDATA()

!! read bzdata; See use m_read_bzdata,only:
!! Use off-regular mesh for qbzreg()=F See hx0fp0.m.sc.F also.
!! This must be consistent with qg4gw.F-mkqg.F
      if(.not.qbzreg()) then 
         deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
         do i=1,nqbz
            qbz(:,i) = qbz(:,i) - deltaq/2d0
            write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
         enddo
      endif
      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp

C --- Use regular mesh even for bzcase==2 and qbzreg()=T
! A little confusing...
c      ddq = 0d0
c      if(bzcase()==2) ddq= dq_
c      do iq = 1, nqbz
c        qbz(1:3,iq) = qbz(1:3,iq) + ddq
c        ! This new qbz is regular mesh, which are identical in the both bzcase.
c      enddo
c      if(.not.qbzreg()) then ! off-regular mesh case
c        do i=1,nqbz
c          qbz(:,i) = qbz(:,i) - dq_
c        enddo
c      endif
      if(MPI__root) then
         do i=1,nqbz
            if(i<10.or.i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
            if(i==10.and.nqbz>18) write(6,"('... ')")
         enddo
         write(6,*)' nqbz nqibz =',nqbz,nqibz
      endif

c$$$!!- oct2005 not implimented cases.
c$$$      if(smbasis().and.chipm) then
c$$$        write(6,*)' smbasis=T & chipm=T is not implimented yet.'//
c$$$     &       ' Supply consistent MixSpin for smbasis!'//
c$$$     &       ' MixSpin should be converted at the end of hvccfp0.'
c$$$        call rx( ' smbasis=T & chipm=T is not implimented yet.')
c$$$      endif

c      call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
c      write(6,"(' ScaledGapX0=',f5.3)") sciss

!! === Readin by genallcf ===
!! See "use m_genallcf_v3" at the begining of this routine
!! We set basic data.

c$$$      if(epsmode) then
c$$$        nwin = -999
c$$$      else
c$$$        nwin = 0      !Readin nw from NW file
c$$$      endif
      incwfin= 0  !use ForX0 for core in GWIN
c      efin =  0d0 !readin EFERMI
c--- EFERMI
      ifief=ifile_handle()
      open(ifief,file='EFERMI')
      read(ifief,*) ef
      close(ifief)
      write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
      call genallcf_v3(incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      tpioa=2d0*pi/alat

      if(chipm.and.nspin==1) call rx( 'chipm mode is for nspin=2')
      debug=.false.; if(verbose()>=100) debug=.true.
      if(debug) write(6,*)' end of genallc'
c      write(6,"(' ncore=',i4)") ncore
c      write(6,*) 'nw_input delta=',nw_input,delta

!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call rx( ' nclass /= natom ')

!! --- tetra or not
      if(delta <= 0d0) then
c        tetra =  .true.
        delta = -delta
        write(6,*)' hx0fp0: tetrahedron mode delta=',delta
      else
c        tetra = .false. ! switch for tetrahedron method for dielectric functions
        call rx(' hx0fp0: only tetra=T support')
      endif

!! --- read dimensions of h,hb
      ifhbe  = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg !warn nband is in m_zmel
      is = iclose('hbe.d')
      if(nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
      if(nqbz /=nqbzt ) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')

!! --- Readin Offset Gamma --------
      if(debug) write(6,*) 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      write(6,*) ' ### nqibz nq0i=', nqibz,nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      do i=1,nq0i
        read (101, * ) wqt(i),q0i(1:3,i)
      enddo
      nq0ix = nq0i
      do i=1,nq0i
        if(wqt(i)==0d0 ) then
          nq0ix = i-1
          exit
        endif
      enddo
      neps = nq0i - nq0ix  ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      write( 6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)
c$$$      if(.not.newaniso2) then
c$$$        wqtsum = sum(abs(wqt(1:nq0i)))
c$$$        call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
c$$$        if(normalm.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) 
c$$$     &  call rx( ' wqtsum of Q0P /=1 ')
c$$$      endif

C --- readin by rdpp ; Radial integrals ppbrd and plane wave part
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx
! qibze(3,nqbze) qbze(3,nqibze)
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
        enddo
      enddo
      ngrpx = 1
      l2nl=2*(nl-1)
      allocate(symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      ificlass=ifile_handle()
      open (ificlass,file='CLASS')
      allocate(iclasst(natom),invgx(ngrp)
     & ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
      write(6,*)'  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(ificlass,*) ibasx, iclasst(ibas)
        write(6, "(2i10)") ibasx, iclasst(ibas)
      enddo
      close(ificlass)
!! Get space-group transformation information. See header of mptaouof.
      call mptauof(symope,ngrpx,plat,natom,pos,iclasst
     o ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
      if(verbose()>=40) write (*,*)' hsfp0.sc.m.F: end of mptauof'
!!  ppbrd = radial integrals,  cgr = rotated cg coeffecients.
      call rdpp(nxx, nl, ngrpx, nn, nclass, nspin, symope,qbas)
      ntq = nband
      allocate(itq(ntq)) !itq=i for i=1,ntq. a dummy. c.f. hsfp0.sc.F
      do i=1,ntq
        itq(i)=i
      enddo
!! Pointer to optimal product basis
c      allocate(imdim(natom))
c      call indxmdm (nblocha,nclass,iclass,natom,
c     o imdim )                  !use in m_zmel
c      if(.not.smbasis()) nblochpmx = nbloch + ngcmx
      nblochpmx = nbloch + ngcmx
      allocate(ngveccB(3,ngcmx)) ! work arry

!! ... for legas test (not used so often. To compare homogeneos electron gas).
c      legas = .false.
c      INQUIRE (FILE = 'LEGAS', EXIST = legas)
!!
      iqxend = nqibz + nq0i
      write(6,*) ' nqibz nqibze=',nqibz,nqibze

!! Initialization of readEigen !readin m_hamindex
      call init_readeigen(ginv,nspin,nband,mrece)!EVU EVD are read in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)

!! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
      call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
      if (nctot > 0) Emin=minval(ecore(:,1:nspin))
      omg2max = (Emax-Emin)*.5d0+.2d0 
             ! (in Hartree) covers all relevant omega, +.2 for margin
      if(MPI__root) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
      realomega = .true.
      imagomega = .true.
      tetra     = .true.
      if(imagonly) then !WVI only for imagonly for ixc==12
        realomega =.false.
        imagomega =.true.
      endif
      if(realonly) then !epsPP noLFC mode for ixc==13
        realomega =.true.
        imagomega =.false.
      endif
!! getfreq returun date given at " use m_freq,only:". 
      call getwemax(.false.,wemax) !wemax is to determine nw !real axis divisions
      if(MPI__root) write(6,"(' wemax=  ',f13.4)") wemax
      call getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,MPI__root)
      nwp = nw+1
      if(.not.imagomega) niw=1
!! ... get eigenvector corresponds to exp(iqr) (q is almost zero).
      if(epsmode) then !iepsmode/=0) then ;  write(6,*) ' read in Mix0vec'
        allocate(epsi(nw_i:nw,neps)) !5July2005 nwp should be used after it is defined!
      endif

!! Miyake tote mode Nov2004. Need fixing.
      if(ecorr_on>0) then  !it was bzcase()==2 Was it bug?
        allocate(nstibz(nqibz))
        do iq=1,nqibz
          iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
          nstibz(iq) = nstbz(iqbz)
c          write(6,"(' iq qibz nstibz=',i5,3f9.4,i5)")iq,qibz(:,iq),nstibz(iq)
        enddo
      endif

!! tetra init
c      call getkeyvalue("GWinput","TmpWWK",tmpwwk,default=.false.)
c      if(tetra) then
        allocate( !wgt(nband+nctot,nband,nqbz), !noccxvv(nw+niw),
     &    nbnb(nqbz,npm),nbnbtt(nqbz,npm), ekxx1(nband,nqbz), ekxx2(nband,nqbz)) !!! nband=nlmto 
c        if(tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
c      endif ;  if(debug) write(6,*)' xxx1:'
      noccxv     = maxocc2 (nspin,ef, nband, qbze,nqbze) ! maximum no. occupied valence states
      if(noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
      noccx      = noccxv + nctot

C allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )
c$$$C ... This is just to get nblochpmx
c$$$      if(smbasis()) then
c$$$        call getngbpomat(nqibz+nq0i,  nnmx,nomx)
c$$$        nblochpmx = nnmx
c$$$      endif

      nprecx = ndble  !We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
      if (MPI__root) then
        ifwd   = iopen('WV.d',1,-1,0)
        write (ifwd,"(1x,10i14)") nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
        ifwd = iclose('WV.d'); ifwd=0
      endif
      allocate(  zw(nblochpmx,nblochpmx) )
      nspinmx = nspin

!!... these are used x0k
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2

!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,R(r))> 
!! This is general for rotated CG coefficient; but hx0fp0 mode is only for  ngrpx=1 (not rotated).
!! Compare usage in hsfp0 modes.
      irot=1
      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,irot,nspin))
      do is = 1,nspin
        call ppbafp_v2 (irot,ngrpx,is,nspin,
     i   il,in,im,nlnm,      
     i   nl,nn,nclass,nlnmx,
     i   mdimx,lx,nx,nxx,       !Bloch wave    
     i   cgr, nl-1,             !rotated CG
     i   ppbrd,                 !radial integrals
     o   ppbir(:,irot,is))      ! this is in m_zmel, used to generate <phi|phi B>
      enddo  
      if(debug)write(6,*) ' end of ppbafp_v2'

!! Set iqxini
      if(omitqbz) then
        iqxini= nqibz + 1
      else
        iqxini= 1
      endif

!! check write 1st part for Ec mode to ecorr.chk Nov2004
      if(ecorr_on>0) then
!!!!!!!!!!!!! this path is under developing. !!!!!!!!!!!
        call rx(' ! hx0fp0: need to fix this path. check subroutine getwk and so on in this path')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(ecorr_on >0.and. (.not.imagomega) ) 
     &   call rx( ' hx0fp0: ecorr_on .and. (.not.imagomega)')
        write(ieclog, "('   iq                q                       wk')")
        do iqq = iqxini,iqxend
          call getwk(iqq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)
          write(ieclog,"(i5,3x,3f12.8, f15.5)") iqq, q, wk4ec
        enddo
        write(ieclog,*)
        write(ieclog,"('   iw omega(Ry)       wiw')")
        do iww=1,niw
          write(ieclog,"(i5, f10.5, f10.5)") iww,2d0*freq_i(iww),wiw(iww)
        enddo
        write(ieclog,*)
        write(ieclog,"(' Note:IntWgt=wk*wiw.',
     &'  Ec =\sum_{k,iw) IntWgt(k,iw)*ecqw(k,iw)')")
        close(ieclog)
        open(ieclog,file="ecorr.chk",access='append')

        call getkeyvalue("GWinput","necut_p",necut, default=1 )
        allocate(totexc(necut),trpv(necut),trlog(necut))
        totexc = 0d0
        trpv   = 0d0
        trlog  = 0d0
      else
        necut=1
      endif
!!
      allocate(ecut(necut),ecuts(necut))
      call getkeyvalue("GWinput","ecut_p" ,ecut, necut,default=(/1d10/) )
      call getkeyvalue("GWinput","ecuts_p",ecuts,necut,default=(/1d10/) )
!!
      if( chipm ) then
        nmbas=natom
        allocate(imbas(nmbas),imbas_s(nmbas))
        istat=-9999 ! istat=-9999 means noumber of readin arguments is returened in istat.
        call getkeyvalue("GWinput","MagAtom", 
     &        imbas,nmbas,status=istat)
        nmbas = istat
        write(6,*)
        write(6,"('Readin MagAtom nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)
        imbas_s(1:nmbas) = imbas(1:nmbas)
        imbas(1:nmbas)   = abs(imbas(1:nmbas))
        allocate(jcoup(nw_i:nw,neps) )
        allocate( svec(nbloch,nmbas) )  !sep2006
        svec=0d0
        allocate( cvec(nbloch,nmbas),momsite(nmbas),
     &    mmnorm(nmbas))              !May2007
        cvec=0d0
        do imb=1,nmbas
          ibas= imbas(imb)
          ifv = iopen ('MixSpin.'//charnum3(ibas),1,3,0)
          read(ifv,*) ibasx,lxx
          allocate(nxx_r(0:lxx))
          do i=0,lxx
            read(ifv,*) nxx_r(i)  !   write(6,"(2i5,d13.6)") nxx_r(i)
          enddo
          allocate(spinvec((lxx+1)**2,maxval(nxx_r)))
          allocate(consvec((lxx+1)**2,maxval(nxx_r)))
          spinvec=0d0
          do ilmx = 1, (lxx+1)**2
            lb = ll(ilmx )       !  write(6,*)' lb=',lb,lxx,ilmx
            do ixx = 1, nxx_r(lb)  !  write(6,*)' nn=',nn,nxx_r(lb)
              if(ilmx==1) then
                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2
     &       ,consvec(ilmx,ixx)
              else
                read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),dumm1,dumm2
     &       ,consvec(ilmx,ixx)
              endif
!           write(6,"(2i5,d13.6)") ilmx, ixx, spinvec(ilmx,ixx)
            enddo
          enddo
!! Calculate ChiPM. So sign of omega should be correct.
          if(imb==1) then !determine spin direction with respect to ibas=imbas(imb=1)
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
c         write(ifchipmlog,"(2e23.15,' ! mmom mmnorm')")momsite(imb),mmnorm(imb)
          write(6,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
        enddo
      endif

! I assume 1 is for majority for eiqr case.
c      if(ix==10222) then
cc        schi=1d0 !1d0 means Majority is isp=1. If Majority is isp=2, use schi=-1d0.
cc        allocate(jcoup(nw_i:nw,neps))
c        mmnorm=1d0
c      endif
c
c! nmbas_in is for rcxqmean
c      if(chipm .and. nolfco) then !.and. ix/=10222) then
c        nmbas_in = nmbas
c      else
c        nmbas_in = 1
c      endif
c      if(epsmode.and.nolfco) then
c        allocate( rcxqmean(nwhis,npm,nmbas_in,nmbas_in))
c        if(debug) write(6,"('fff:',3i5)") nwhis,npm,nmbas_in
c      else
c        allocate( rcxqmean(1,1,1,1)) !dummy
c      endif
c
c      if(chipm) allocate(eiqrm(nmbas)) 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cctakao
c$$$      allocate( x0meanx(nmbas,nmbas) )
c$$$      allocate( x0mat(nmbas,nmbas),x0matinv(nmbas,nmbas) )
c$$$      do 1101 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
c$$$         if(iq==iqxini+2) exit
c$$$         q = qibze(:,iq)
c$$$         write(6,*)'aaaaaaaaaa q=',q
c$$$          read(ifgb0vec,*) qgbin(1:3),ngbin,igc0,dnorm
c$$$          if(sum(abs(q))==0d0)then
c$$$            if(sum(qgbin**2) >1d-7)stop'qgbin=0 xxx See hvccfp0'
c$$$          elseif(sum(abs(qgbin(1:3)-q)) >1d-8)then
c$$$            stop'qgbin inconsistent'
c$$$          endif
c$$$          write(6,"(' --- Readin Mix0vec: ',3d13.6,2i5,d18.8)")
c$$$     &      qgbin(1:3),ngbin,igc0,dnorm
c$$$c          if(ngb/=ngbin) stop 'hx0fp0: ngb/=ngbgin'
c$$$          ngb=ngbin
c$$$          write(6,"(' ngb nwp niw=',3i8)")ngb,nwp,niw
c$$$          nmbas_in=1
c$$$          allocate( gbvec(ngb),zzr(ngb,1),x0mean(nw_i:nw,1,1))
c$$$          x0mean=0d0
c$$$          do i=1,ngb
c$$$            read(ifgb0vec,"(4d24.15)") zz1,zz2,zz3,zz4
c$$$            gbvec(i)= dcmplx(zz1,zz2)
c$$$            zzr(i,1)= dcmplx(zz3,zz4)
c$$$          enddo
c$$$          write(6,"(' normchk=',255e23.15)") sum( dconjg(gbvec)*zzr(:,1) )
c$$$     &        ,sum(abs(gbvec(:))), sum(abs(zzr(:,1)))
c$$$            allocate(eiqrm(nmbas))
c$$$            do imb=1,nmbas
c$$$            eiqrm(imb)= sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )
c$$$            write(6,"(' <eiqr|m> =',255e23.15)") eiqrm(imb)
c$$$            if( imbas_s(imb)<-1) eiqrm(imb)= -eiqrm(imb)
c$$$            enddo
c$$$            write(6,"('<eiqr|m>:Set \pm in GWinput(for stuggard chi)')")
c$$$              iqixc2 = iq- (nqibz+nq0ix)
c$$$              ifx = iopen ('StonerNLFC.dat',1,3,0)
c$$$              read(ifx,*) jzero2
c$$$              ifx= iclose('StonerNLFC.dat')
c$$$
c$$$              ifchipm2=iopen(
c$$$     &                 'ChiPM'//charnum4(iqixc2)//'.nolfc.mat',1,3,0)
c$$$              do iw=1,10
c$$$              read(ifchipm2,
c$$$     &        '(36x,2x,20x,2x,255e23.15)') x0meanx(:,:)
c$$$              write(6,'("xxx x0mat=",255d13.5)') x0meanx
c$$$              x0matinv=x0meanx
c$$$              call matcinv(nmbas,x0matinv)
c$$$              do i=1,nmbas
c$$$                x0matinv(i,i)= x0matinv(i,i) - jzero2 ! (chipm_0^+-)^-1 - I
c$$$              enddo
c$$$              x0mat = x0matinv
c$$$              do i=1,nmbas
c$$$              x0mat(i,i) = x0mat(i,i)+ img*1d-30 ! to avoid inversion error.
c$$$              enddo
c$$$
c$$$              call matcinv(nmbas,x0mat) !this is full x0_+-
c$$$              trr = sum( eiqrm*matmul(x0mat,dconjg(eiqrm)) ) !*mmnorm
c$$$              write(6,
c$$$     &        '("ttt",3f12.8,2x,f10.5,2x,2e23.15,2x,2e23.15)') q, 2*schi*frr, trr,1d0/trr
c$$$              enddo
c$$$            ifx=iclose(ifchipm2)
c$$$ 1101  continue
c$$$          stop 'xxxxxxxxxxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!! -- EIBZ mode for nolfco -----------------------
      eibzmode=eibz4x0()
!! If eibzmode=T, it is efficient but can slightly break crystal symmetry.

!! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
c      if(nolfco) then
c          eibzmode = .false. 
c      endif    
c!! ----------------------------------------------

!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      if(eibzmode) then
!! commentout block inversion Use iqxendx=iqxend because of full inversion
        call cputid(0)
        write(6,*)' ---goto eibzmode block ---'
        iqxendx=iqxend
        if(epsmode) iqxendx=iqxend
        allocate( nwgt(nqbz,iqxini:iqxendx), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     &   igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     &   eibzsym(ngrp,-1:1,iqxini:iqxendx))
!! Check timereversal is required for symmetrization operation or not. If tiii=timereversal=F is enforced,
!! the symmetrization procedure in x0kf_v4h becomes a little time-consuming.
        write(6,*)
        write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal() 
        if(MPI__root) iprintx=.true.
        call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,
     i   timereversal(),ginv,iprintx,
     o   nwgt,igx,igxt,eibzsym,tiii)
        write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
        call cputid(0)
! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
        call PBindex(natom,lx,l2nl,nx) !all input. Returns requied index stored in arrays in m_pbindex.
        call cputid(0)
        call readqgcou()        !no input. Read QGcou and store date into variables.
c     call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.
      else                      !dummy allocation to overlaid -check bound !sep2014
        iqxendx=iqxend
        allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &   ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      endif

      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
!! == Calculate x0(q,iw) and W == main loop 1001 for iq. 
!! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
!! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
!! (or need to modify do 1001 loop).
!! iq>nqibz for ixc=11 is not time-consuming.
      call MPI__hx0fp0_rankdivider2(iqxini,iqxend)

!! ---------------------------------------------------------------
!! === loop over iq ============================================
!! ---------------------------------------------------------------
      do 1001 iq = iqxini,iqxend  ! NOTE: q=(0,0,0) is omitted when iqxini=2
        if( .not. MPI__task(iq) ) cycle
        if(ixc==101.or.normalm) then
          ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
        endif  
        if (normalm) then
          ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
        endif
        call cputid (0)

        q = qibze(:,iq)
        call readqg('QGcou', q, ginv,  quu,ngc,ngveccB) 

!! Caution : confusing point
!!  ngc by QGcou is shown at the bottom of lqg4gw.
!!  ngc read from PPOVL are given by rdata4gw---> ngc(iq>nqibz )=ngc for q=0
!!
c        if( newaniso2.and.iq==1 ) then ! *sanity check
        if( iq==1 ) then ! *sanity check
          if(sum(q**2)>1d-10) then
            call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
          endif
        endif

!! ==== readin Coulomb matrix ====
        ngb = nbloch + ngc !ngb is readin from vcoul 25jan2006
        write(6,*)
        write(6,"('===== do 1001: iq q=',i7,3f9.4,' ========')")iq,q !qq
        write(6,"('  nbloch ngb ngc=',3i10)") nbloch,ngb,ngc

!! === readin diagonalized Coulomb interaction ===
!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
c        if(newaniso2.and.(.not.chipm)) then
        if((.not.chipm)) then
          vcoudfile='Vcoud.'//charnum5(iq)  !this is closed at the end of do 1001.  iq was iqqv
          ifvcoud = iopen(trim(vcoudfile),0,-1,0)
          read(ifvcoud) ngb0
          if( ngb0/=ngb ) call rx( 'hx0fp0.m.f:ngb0/=ngb')
          read(ifvcoud) qvv
          if(sum(abs(qvv-q))>1d-10) then
            write(6,*)'qvv =',qvv
            call rx( 'hx0fp0: qvv/=0 hvcc is not consistent')
          endif
          if(allocated(zcousq)) deallocate( zcousq,vcousq )
          allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
          read(ifvcoud) vcousq
          read(ifvcoud) zcousq
          idummy=iclose(trim(vcoudfile))
          vcousq=sqrt(vcousq)
          if(allocated(zzr)) deallocate(zzr)
          allocate(zzr(1,1)) !dummy
          zzr=0d0
        endif
      
        if(chipm .and. nolfco) then !.and. ix/=10222) then
          nmbas_in = nmbas
        elseif(nolfco) then
          nmbas_in = 1
        else
          nmbas_in = ngb
        endif
        nmbas1 = nmbas_in
        nmbas2 = nmbas1

!! ==== set up for epsilon mode =====
        if(epsmode) then 
          iqixc2 = iq- (nqibz+nq0ix)
          if((.not.chipm).and.nolfco) then
            allocate( x0mean(nw_i:nw,1,1) )
            x0mean=0d0
          endif  
!! zzr is only for chipm.and.nolfco mode
          if(chipm .and.nolfco) then
            allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,nmbas))
            x0mean=0d0
            zzr   =0d0
            zzr(1:nbloch,1:nmbas) = svec(1:nbloch,1:nmbas)
          endif
!! ... Open ChiPM*.nolfc_mat
          if( wqt(iq-nqibz)==0d0.and.chipm ) then
            ifchipmn_mat=iopen('ChiPM'//charnum4(iqixc2)//'.nlfc.mat',1,3,0)
            write(ifchipmn_mat,"(255i5)") nmbas
            write(ifchipmn_mat,"(255i5)") imbas(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)") momsite(1:nmbas)
            write(ifchipmn_mat,"(255e23.15)")  mmnorm(1:nmbas)
c            write(ifchipmn_mat,"(255e23.15)") eiqrm(1:nmbas)!if necessary, fix code to give eiqrm. takaoAug2012
            write(ifchipmn_mat,"( ' Here was eiqrm: If needed, need to fix hx0fp0')") 
            if(.not.nolfco) then
              ifchipm_fmat=iopen('ChiPM'//charnum4(iqixc2)//'.fmat',0,3,0)
              write(ifchipm_fmat) nbloch, natom,nmbas, iqxini,iqxend, nw_i,nw
              write(ifchipm_fmat) imbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
              write(ifchipm_fmat) nblocha(1:natom),svec(1:nbloch,1:nmbas)
              write(ifchipm_fmat) zzr0(1:nbloch) !zzr(1:nbloch,1)
            endif
          elseif(wqt(iq-nqibz)==0d0.and.(.not.chipm)) then
!! ... Open EPS* file
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nlfc.dat'
            ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
            write(ifepsdatnolfc,"(a)")' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
            if(.not.nolfco) then 
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              ifepsdat = iopen ( fileps,1,3,0)
              write(ifepsdat,"(a)") ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
            endif
          endif
        endif

        if(epsmode.and.nolfco) then !iepsmode==202) then
c          rcxqmean=0d0
        else
          write(6,*) "rcxq alloc ngb nwhis npm ---",ngb,nwhis,npm
          allocate( rcxq(ngb,ngb,nwhis,npm) )
        endif

!! === zmelt conversion on different basis. ppovlz is used in get_zmelt2 in m_zmel (called in x0kf_v4h).
        if(chipm.and.nolfco) then
          if(allocated(ppovlz)) deallocate(ppovlz)
          allocate(ppovlz(ngb,nmbas1))
          ppovlz= zzr
        elseif(nolfco .and. nmbas1==1) then !for <e^iqr|x0|e^iqr>
c          if(allocated(ppovlzinv)) deallocate(ppovlzinv)
          if(allocated(ppovlz)) deallocate(ppovlz)
          if(allocated(ppovl)) deallocate(ppovl)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb)) !,   ppovlzinv(ngb,ngb))
          call readppovl0(q,ngc,ppovl)
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
          ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
          write(6,*)'nnnnn',nbloch+ngc,ngb
        else                    !may2013  this removes O^-1 factor from zmelt
c          if(allocated(ppovlzinv)) deallocate(ppovlzinv)
          if(allocated(ppovlz)) deallocate(ppovlz)
          if(allocated(ppovl)) deallocate(ppovl)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb)) !,   ppovlzinv(ngb,ngb))
          call readppovl0(q,ngc,ppovl)
          allocate(ppovl_(ngb,ngb))
          ppovl_=0d0
          do i=1,nbloch
            ppovl_(i,i)=1d0
          enddo
          ppovl_(nbloch+1:nbloch+ngc,nbloch+1:nbloch+ngc)=ppovl
          if(.not.eibz4x0()) then !sep2014 added for eibz4x0=F
            ppovl_= matmul(ppovl_,zcousq)
          endif
          ppovlz = ppovl_
          deallocate(ppovl_,ppovl)
        endif

!! takao apr2012
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
!! ---------------------------------------------------------------
!! === loop over spin=== =========================================
!! ---------------------------------------------------------------
        do 1003 is = 1,nspinmx
          write(6,"(' ##### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ##### ')")iq, is, nqibz + nq0i,nspin
          if(debug) write(6,*)' niw nw=',niw,nw
!! ==== spin chi_charge or chi_+- ====
          isf=is
          if(chipm) then
            write(6,*)" chi_+- mode ixc=",ixc
            if(is==1) isf=2
            if(is==2) isf=1
            rcxq=0d0
          endif

!! Tetrahedron weight.
!!     nbnbx
!!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
!!     nhw(ibjb,kx): the number of histogram sections
!!     jhw(ibjb,kx): pointer to whw
!!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
!!     : histogram weights for given ib,jb,kx for histogram sections
!!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c            write(6,*) ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2
          call gettetwt(q,iq,is,isf,nwgt(:,iq))

!! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 ==
          call cputid(0)
          if(npm==1) then
            ncc=0
          else
            ncc=nctot
          endif
          symmetrize=.false.
          call x0kf_v4hz(npm,ncc,   
     i        ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
     i        n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf, symmetrize, !
     i              qbas,ginv,  qbz,wbz, 
     d            nlmto,nqbz,nctot, !noccx,noccxv,
     d            nbloch,  nwhis, !nlnmx,mdimx,
     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
     i    nqbze,nband,nqibz, 
     o    rcxq,   !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,nmbas_in, zcousq,   !ppovl,nmbas1,nmbas2, is removed ppovlz, 
     i    chipm,eibzmode, !z1offd,!for nolfco Add nmbas Sep2006
     i    nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),ngrp, eibzsym(:,:,iq),crpa)
!! ------Question, Apr2015takao. -------------------
!! ???  we may need 
!! ???  "if(is==nspinmx.or.chipm) then" for chipm mode.
!! ???  really OK ??? Need check more... Compare with old code...
!! ---------------------------------------------------
!kino 2014-08-19 add

          if (is==nspinmx.or.chipm) then ! bugfix? at apr2015. Takao think ".or.chipm" is required for chipm mode
                                         ! Because rcxq is calculated for each is, symmetrized and its contribution
                                         ! is added to zxq in dpsion5.
            symmetrize=.true.
            call x0kf_v4hz(npm,ncc,
     i       ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i       n1b,n2b,nbnbx,nbnb, ! use whw by tetwt5 ,
     i       q,
     i       nspin,is,isf, symmetrize, !
     i       qbas,ginv,  qbz,wbz,
c     i       nblocha,           !nlnm,nlnmv,nlnmc,iclass,
c     i              ppb(1,is),
c     i              icore,ncore,
     d       nlmto,nqbz,nctot,  !noccx,noccxv,
c     d       natom,             !nl,nclass,natom,nnc,
     d       nbloch,  nwhis,    ! nlnmx,mdimx,
     i       iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
     i       nqbze,nband,nqibz,
     o       rcxq,              !  rcxq is the accumulating variable for spins
     i       nolfco,zzr,nmbas_in, zcousq, !ppovl,nmbas1,nmbas2, is removed ppovlz, 
     i       chipm,eibzmode,    !z1offd,!for nolfco Add nmbas Sep2006
     i       nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),ngrp, eibzsym(:,:,iq),crpa)
          endif
          deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)
          iecut=1
          if(debug) write(6,"(a)") ' --- goto dpsion5 --- '
          if(is==nspinmx.or.chipm) then
            write(6,"('  nmbas1,nmbas2=',2i10)") nmbas1,nmbas2
            call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i       rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o       zxq, zxqi,
     i       chipm, schi,is,  ecut(iecut),ecuts(iecut))
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
          if(chipm) then !ixc==22.or.ixc==23) then
            if (nspin==1) call rx( 'chipm modes are for nspin==2')
          elseif(epsmode.and.nolfco.and.(.not.chipm)) then !ixc==2.or.iepsmode==202) then
            if (nspin==1) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
            if (nspin==1) zxq = 2d0*zxq       !if paramagnetic, multiply x0 by 2
          else
            if (nspin == 1) zxq = 2d0*zxq   !if paramagnetic, multiply x0 by 2
          endif

c         write (ifxd,"(1x,3f10.4)") q(1),q(2),q(3)
c         write (ifrx) rxq,cxq
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
            if(.not.epsmode) then  !if(ixc==1.or.sergeyv) then
              imode = 1
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
c$$$c              write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$              write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,1d0/epstinv(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$c             write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(1,2:10:3),zw0(1,63:70:3)
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
c$$$                write(*,"('mmmmzwp99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0), !(1d0/llw(iw,iq0)-1d0)*vcousq(1)**2,
c$$$c     &                  w_k(2:10:3)/llw(iw,iq0), w_k(63:70:3)/llw(iw,iq0)
c$$$     &                  w_ks(2:10:3)/llw(iw,iq0), w_ks(63:70:3)/llw(iw,iq0)
c$$$                write(*,"('mmmmzwp99x ')")
c$$$                endif
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                do igb1=1+ix,ngb
                  do igb2=1+ix,ngb
                    zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                    if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                  enddo
                enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmmmm
c       if(iq>nqibz) then
c         write(*,"('mmmmz99x ',i3,10(2d13.5,2x))") iw,zw0(1,1)+vcousq(1)**2,zw0(2:10:3,1),zw0(63:70:3,1)
c       endif
c       if(iq==1.or.iq>nqibz) then
c                write(*,"('mmmz0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c                write(*,"('mmmz0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2:10:3,1,iw),zxq(63:70:3,1,iw)
c         write(*,"('mmmmz99x ',i3,10(2d13.5,2x))") iw,zw0(1,1)+vcousq(1)**2,zw0(1,2:10:3),zw0(1,63:70:3)
c       write(*,"('mmmzx  ',2i3,10(2d13.5,2x))") iq,iw,zxq(2,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c       write(*,"('mmmzx  ',2i3,10(2d13.5,2x))") iq,iw,zxq(3,1,iw),zxq(3,2:10:3,iw),zxq(3,63:70:3,iw)
c       write(*,"('mmmzxs ',2i3,10(2d13.5,2x))") iq,iw,zxq(1,1,iw),zxq(2:10:3,1,iw),zxq(63:70:3,1,iw)
c       write(*,"('mmmzxs ',2i3,10(2d13.5,2x))") iq,iw,zxq(1,2,iw),zxq(2:10:3,2,iw),zxq(63:70:3,2,iw)
c       write(*,"('mmmmzee',2i3,10(2d13.5,2x))")iq,iw,epstilde(2,2),epstilde(2,2:10:3),epstilde(2,63:70:3)
c       write(*,"('mmmmzee',2i3,10(2d13.5,2x))")iq,iw,epstilde(3,2),epstilde(3,2:10:3),epstilde(3,63:70:3)
c       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              endif

!! for iq>nqibz
c              if(newaniso2.and.iq>nqibz) then
              if(iq>nqibz) then
!! Full inversion to calculalte eps with LFC.
                ix=0
                vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
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
                write(6,*) 'iq iw_real  eps(withLFC) eps(woLFC)',iq,iw,llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
c$$$c                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$c                write(6,*)'sumcheck epstinv=',sum(abs(epstinv(2:ngb,2:ngb)))
c$$$c                if(iw/=iwx) stop 'hx0fp0: iw/=iwx'  !sanity check
c$$$c                sk  (1:ngb)= zxq(1,1:ngb,iw)
c$$$c                sks (1:ngb)= zxq(2,1:ngb,iw) ! WARNING: zxq(2,1:ngb,iw) contains zxq(1:ngb,1,iw).
c$$$                                             ! A little confusing. See nolfco=T case in x0kf_v4h.F.
c$$$                                             !   
c$$$cc                sk(  1:ngb)= zxq(1,1:ngb,iw)
c$$$cc                sks( 1:ngb)= zxq(1:ngb,1,iw)
c$$$                vcou1= fourpi/sum(q**2*tpioa**2) !test --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
c$$$                vcou1sq= sqrt(vcou1)   ! only vcousq(1) should be replaced.
c$$$c                write(ifiss)iw,iq0,ngb,q
c$$$c                write(ifiss)vcou1,vcou1sq,vcousq(2:ngb),sk(1:ngb),sks(1:ngb)
c$$$c                w_k(1) =0d0
c$$$c                w_ks(1)=0d0
c$$$c                w_k( 2:ngb)= vcou1sq*matmul( sk(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
c$$$c                w_ks(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), sks(2:ngb)*vcousq(2:ngb))
c$$$cmmm epsPP mode - vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
c$$$c                llw(iw,iq0)=  1d0 -vcou1*sk(1) !- vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
c$$$
c$$$                llw(iw,iq0)=  1d0 - vcou1*zxq(1,1,iw) !- vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
c$$$                write(6,*) 'epsPP iq iw',iq,iw, 1d0 - fourpi* zxq(1,1,iw)/sum(q**2*tpioa**2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                write(*,"('mmmw0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c                write(*,"('mmmw0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c                write(*,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,fourpi/sum(q**2*tpioa**2)/llw(iw,iq0),
c     &                  w_k(2:10:3)/llw(iw,iq0),w_k(63:70:3)/llw(iw,iq0)
c                write(*,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0),
c     &                  w_ks(2:10:3)/llw(iw,iq0),w_ks(63:70:3)/llw(iw,iq0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                write(6,*) 'uuuu iq iw llw vc x0mean',iq,iw,fourpi/sum(q**2*tpioa**2),sk(1)
c                write(ifisk) iw,iq0,q
c                write(ifisk) vcousq(2:ngb)*w_k( 2:ngb),vcousq(2:ngb)*w_ks( 2:ngb)
              endif

c$$$              if(.not.newaniso2) then               ! Original mode
c$$$                call rx( 'not checked here')
c$$$c               call wcf( ngb, vcoul, zxq(1,1,iw), imode, zw0)
c$$$              endif

c$$$!!... a debug mode
c$$$                write(6,"('hhh --- EigenValues for Im( W) --------')")
c$$$                allocate(ebb(ngb))
c$$$                call diagcvh2( (zw0-transpose(dconjg(zw0)))/2d0/img, ngb, ebb)
c$$$                do ii=1,ngb
c$$$                  if( abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
c$$$                    write(6, "('hhhIWq : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
c$$$                  else
c$$$                    write(6, "('hhhIWqxxx : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
c$$$                  endif
c$$$                enddo
c$$$                deallocate(ebb)

c              if(newaniso2.and.iq>nqibz) then
              if(iq>nqibz) then
c                zw(1:ngb,1:ngb) = 0d0
c                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
              else
                zw(1:ngb,1:ngb) = zw0
c                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
                write(ifrcw, rec= iw-nw_i+1) zw   !  WP = vsc-v
                call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
              endif
!! epsmode
            elseif(epsmode.and.(.not.chipm)) then !.and.ixc/=23) then ! No LFC (local field correction). It's better to use echo 4| hbasfp0.
               if(debug)write(6,*) 'xxx2 epsmode iq,iw=',iq,iw
c              write(6,*)'ppppp sumcheck zxq=',sum(abs(zxq)),sum(abs(zzr)),sum(abs(vcoul)),sum(abs(gbvec))
c              if(newaniso2) then
!! there is difference of two vcmean below since we use (sligthy) screened Coulomb (screenfac() in switch.F)
!!   NOTE that we use vcoul with screening (screenfac() is used in hvccfp0.F
c               vcmean = fourpi/sum(q**2*tpioa**2) !aug2012
                vcmean=vcousq(1)**2
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
c$$$              else   
c$$$                 write(6,*)'not support epsmode and newaniso=F mode now.'
c$$$c$$$                vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
c$$$c$$$                if(.not.nolfco) then
c$$$c$$$                 x0mean(iw,1,1) = sum( dconjg(zzr(:,1))* matmul(zxq(:,:,iw),zzr(:,1)))
c$$$c$$$                endif
c$$$c$$$                epsi(iw,iqixc2) = 1d0/(1- vcmean * x0mean(iw,1,1))
c$$$c$$$                write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2e23.15,3x, 2e23.15,
c$$$c$$$     &          " vcmean x0mean =", 2e23.15,3x, 2e23.15)') iqixc2,iw,2*frr,
c$$$c$$$     &          1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean,x0mean(iw,1,1)
c$$$c$$$                write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
c$$$c$$$     &          q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$c$$$                if(.not.nolfco) then ! With LFC
c$$$c$$$                  imode=2
c$$$c$$$                  call wcf( ngb, vcoul, zxq(1,1,iw), imode,
c$$$c$$$     &                zw0) !  write(6,"('ssschk1=',3d13.5)") sum(abs(zw0)) sum(abs(gbvec))
c$$$c$$$                  epsi(iw,iqixc2)= sum( dconjg(gbvec) * matmul(zw0,zzr(:,1)) )
c$$$c$$$                  write(6,'( " iq iw omega eps epsi  wLFC="
c$$$c$$$     &            ,2i6,f8.3,2e23.15,3x, 2e23.15)') 
c$$$c$$$     &            iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$c$$$                  write(6,*)
c$$$c$$$                  write(ifepsdat,'(3f12.8,2x,d12.4,2e23.15,2x,2e23.15)') 
c$$$c$$$     &            q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
c$$$c$$$                endif   
c$$$              endif
C --- ChiPM mode
            elseif(epsmode.and.chipm) then
              allocate( x0meanx(nmbas,nmbas) )
              if(nolfco) then  ! ChiPM mode without LFC
c$$$                if(legas) then
c$$$                   call rx( ' hx0fp0.m.F need to implement thigs here if required')
c$$$! --- three lines below may work for test purpose for legas. But not sure.
c$$$c                  vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
c$$$c                  write(ifchipmn,'(3f12.8,2x,f8.5,2x,2e23.15)')
c$$$c     &            q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw,1,1)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
c$$$                else
c$$$                  x0meanx = x0mean(iw,:,:)/2d0  !in Ry unit.
c$$$                endif
                x0meanx = x0mean(iw,:,:)/2d0 !in Ry unit.
              else
C ... ChiPM mode with LFC... NoLFC part
                zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)/2d0  ! in Ry.
                do imb1=1,nmbas
                  do imb2=1,nmbas
                    x0meanx(imb1,imb2)= 
     &              sum( svec(1:nbloch,imb1)*
     &              matmul(zxq(1:nbloch,1:nbloch,iw),svec(1:nbloch,imb2)))  !/ mmnorm**2  I removed mmnorm may2007
                  enddo
                enddo
!       x0meanx= <m|chi^+-(\omega)|m>/<m|m>**2
              endif
              do imb1=1,nmbas
                do imb2=1,nmbas
                  x0meanx(imb1,imb2) = 
     &            x0meanx(imb1,imb2)/mmnorm(imb1)/mmnorm(imb2)
                enddo
              enddo
              write(ifchipmn_mat,'(3f12.8,2x,f20.15,2x,255e23.15)')q, 2*schi*frr, x0meanx(:,:)
              if(.not.nolfco) write(ifchipm_fmat) q, 2*schi*frr, zxq(1:nbloch,1:nbloch,iw)

c! These lines commented by "c! ' are histories ---> For Takao's memo. Maybe not so useful for others.
c! ! for NoLFC, Get I from q=0, and calculate Tr(Chipm)
c!               allocate( x0mat(nmbas,nmbas),x0matinv(nmbas,nmbas) )
c!               ifx = iopen ('StonerNLFC.dat',1,3,0)
c!               if(iw==0 .and. sum(q**2) <1d-13) then
c!                 x0mat = x0meanx
c!                 x0mat(:,:)= x0mat +transpose(dconjg(x0mat))
c!                 x0matinv= 0.5d0*x0mat
c!                 call matcinv(nmbas,x0matinv)
c!                 write(6,*) ' q=',q
c!                 write(6,*) ' nmbas ifx=',nmbas,ifx
c!                 write(6,*) ' x0matinv=',x0matinv
c!                 allocate(evall(nmbas))
c!                 call diagno00(nmbas,x0matinv,evall)
c! ! Note that x0matinv at omega=0 is negative definite matrix (by definition).
c!                 do i1=1,nmbas
c!                   write(6,'(" eval(iw=0)=",i5,f15.5)') i1, -evall(i1)
c!                 enddo
c!                 jzero2 = minval(-evall)
c!                 deallocate(evall)
c!                 write(ifx,"(e23.15)")  jzero2
c!                 do imb=1,nmbas !temporary
c!                   write(ifx,"(e23.15,' ! tttt temporary... U_mm in eV')")
c!      &            rydberg()*jzero2*mmnorm(imb)**2/momsite(imb)**2
c!                 enddo
c!               elseif(iw==0) then
c!                 read(ifx,*,end=1013,err=1013) jzero2
c!                 goto 1014
c!  1013           continue
c!                 stop " i/o error StonerNLFC.dat"
c!  1014           continue
c!               endif
c!               ifx= iclose('StonerNLFC.dat')
c!               if(onceww(6)) write(6,*)' i/o end: StonerNLFC.dat'
c! !
c!               x0matinv = x0meanx
c!               call matcinv(nmbas,x0matinv)
c!               do i=1,nmbas
c!                 x0matinv(i,i)= x0matinv(i,i) + jzero2 ! (chipm_0^+-)^-1 + I
c!               enddo
c!               x0mat= x0matinv
c!               do i=1,nmbas
c!                 x0mat(i,i) = x0mat(i,i)+ img*1d-30 ! to avoid inversion error.
c!               enddo
c!               call matcinv(nmbas,x0mat) !this is full x0_+-
c!               trr = sum( eiqrm*matmul(x0mat,dconjg(eiqrm)) )
c!               write(ifchipmn,
c!      &        '(3f12.8,2x,f20.15,2x,2e23.15,2x,2e23.15)') q, 2*schi*frr, trr,1d0/trr
c!               deallocate( x0mat,x0matinv)
c!
c! C--- With LFC ! save or read Istoner
c!               if(.not.nolfco) then
c!                 zzz = zxq(1:nbloch,1:nbloch,iw)
c!                 ifstoner = iopen ('Stoner.dat',1,3,0)
c!                 if( sum(q**2) < 1d-10 .and. iw==0 ) then
c!                   call diagno00(nbloch,zzz,ss0)
c!                   ! zzz is negative definite at omegw=0 if the ground state is stable.
c!                   ! minval(ss0) is for the largest negative value (softest mode).
c!                   Istoner = -1d0/minval(ss0)
c!                   do ii= 1,nbloch
c!                     if(verbose()>50.or.iw<=2) then
c!                       write(6,"(' eig chi^0_+- =',
c!      &                      i4,d13.5,256d13.5 )" ) ii, ss0(ii)
c!                     endif
c!                   enddo
c!
c! cxxxx thisa SVD procedure is not used now.
c! c!      SVD of chi^-1: !now only look for lowest eigenvalue problem... So rather eigenvalue problem instead of SVD
c! c                  write(6,"(a,i5)")' ----SVD: chiinv --- iw=',iw
c! c                  zxq(1:nbloch,1:nbloch,iw)=zzz
c! c                  call zgesvdnn(
c! c     i            nbloch, zxq(1:nbloch,1:nbloch,iw),
c! c     o            SS0,UU0,VT0)
c! c                  Istoner = -sum(UU0(:,1)*VT0(1,:))/ss0(1)
c! c!                  do ii= 1,nbloch
c! c!                    write(ifstoner,'(4e23.15)') UU0(ii,1),VT0(1,ii)
c! c!                  enddo
c!                   write(ifstoner,"(e23.15)")
c!      &             Istoner
c!                   do imb=1,nmbas
c!                     write(ifstoner,"(e23.15,'!tttt temporary U_mm in eV')")
c!      &            Istoner*rydberg()*mmnorm(imb)**2/momsite(imb)**2
c!                   enddo
c!                 elseif(iw==0) then
c!                   read(ifstoner,*) Istoner
c!                 endif
c!                 ifstoner = iclose('Stoner.dat')
c! C...  <eqir| 1/(1 + I chi^0_+-) | eiqr>
c!                 mmat = + Istoner * zzz
c!                 do i = 1, nbloch
c!                   mmat(i,i) = mmat(i,i) + 1d0
c!                 enddo
c! c                trr0 = sum( dconjg(zzr(1:nbloch,1))*
c! c     &                      matmul( mmat,zzr(1:nbloch,1) )  )
c! c                write(6,"(' <eiqr| 1 + I chi0^+-|eiqr> =',255e23.15)") trr0
c!                 do i=1,nbloch
c!                   mmat(i,i) = mmat(i,i)+ img*1d-30 ! to avoid inversion error.
c!                 enddo
c!
c! c$$$c prtest for NiO with 4 bloch basis
c! c$$$              zzzx = mmat   !matmul(sproj,matmul(mmat,sproj))
c! c$$$                  call zgesvdnn(
c! c$$$     i            nbloch, zzzx,
c! c$$$     o            eex,UU0,VT0)
c! c$$$! projected denominator
c! c$$$              denom = matmul( sproj,matmul(Istoner * zzz,sproj))
c! c$$$              do i = 1, nbloch
c! c$$$                denom(i,i) = denom(i,i) + 1d0
c! c$$$              enddo
c! c$$$!
c! c$$$              zzzx=denom
c! c$$$                  call zgesvdnn(
c! c$$$     i            nbloch, zzzx,
c! c$$$     o            eey,UU0,VT0)
c! c$$$              write(ifchipm2,
c! c$$$     &        '(3f12.8,2x,f20.15,2x,4f11.5,3x,4f11.5)') q, 2*schi*frr, eex,eey
c! c$$$!
c! c$$$              zzzy=denom
c! c$$$              do i=1,nbloch
c! c$$$                zzzy(i,i) = zzzy(i,i)+ img*1d-30 ! to avoid inversion error.
c! c$$$              enddo
c! c$$$              call matcinv(nbloch, zzzy)
c! c$$$Ctest --- another inversion procedure  ! zzzy is the inverse of denom
c! c$$$c              VT = dconjg(transpose(UU0))
c! c$$$c              UU = dconjg(transpose(VT0))
c! c$$$c              zzzy=0d0
c! c$$$c              do i=1,nbloch
c! c$$$c                do ix=1,nbloch
c! c$$$c                do iy=1,nbloch
c! c$$$c                zzzy(ix,iy) = zzzy(ix,iy) + UU(ix,i)*VT(i,iy)/eey(i)
c! c$$$c                enddo
c! c$$$c                enddo
c! c$$$c              enddo
c! c$$$c              zzzx = matmul(denom,zzzy)
c! c$$$c              do i=1,nbloch
c! c$$$c              do j=1,nbloch
c! c$$$c                 write(6,"('zzzx=',2i5,2d13.6)")i,j,zzzx(i,j)
c! c$$$c              enddo
c! c$$$c              enddo
c! c$$$              zzzx = matmul( sproj,matmul(zzz,sproj) )
c! c$$$              mmatx = matmul(zzzx, zzzy)
c! c$$$              trrx  = sum( dconjg(zzr(1:nbloch,1)) *
c! c$$$     &                     matmul(mmatx,zzr(1:nbloch,1)) )
c! c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!                 call matcinv(nbloch, mmat)
c!                 mmat = matmul(zzz,mmat)
c! c--- prtest I found that This makes the difference at high energy part!!! Nov-9-2006
c! c---            mmat = matmul( sproj,matmul(mmat,sproj) )
c!                 trr  = sum( dconjg(zzr(1:nbloch,1)) *
c!      &                     matmul(mmat,zzr(1:nbloch,1)) )
c!                 write(ifchipm,
c!      &         '(3f12.8,2x,f20.15,2x,2e23.15,2x,4e23.15)')
c!      &              q,   2*schi*frr,  trr, 1d0/trr
c!               endif
              deallocate(x0meanx)
            endif
c            write(6,*)'tttt aaa iw=',iw
 1015     continue
c          write(6,*)'tttt  end of  do 1015 loop'
c          if(newaniso2) then
c            if(allocated(sk)) deallocate(sk)!,sks,w_k,w_ks)
c          endif


c
c          if(chipm.and.(.not.nolfco))
c     &          deallocate(sqovlp,sqovlpi,UU0,VT0,SS0,mmat,zzz)
c$$$          if( ixc==5.or.ixc==6 ) then
c$$$            jpm=1
c$$$c            nwmax = nw
c$$$c            if(ixc==5) nwmax =nw
c$$$            allocate(trwv(nw_i:nw),trwv2(nw_i:nw))
c$$$            do iw = nw_i,nw !max ! trace check
c$$$              trwv(iw) = zxq(6,7,iw)
c$$$              trwv2(iw) = 0d0
c$$$              do i = 1,ngb
c$$$                trwv2(iw) = trwv2(iw) + zxq(i,i,iw)
c$$$              enddo
c$$$            enddo
c$$$            do iw= nw_i,nw-1
c$$$              if(ixc==5)
c$$$     &         write(6,"('iq iw[min_max]=',2i5,2f7.4,' trwv by wwk*h= ',
c$$$     &         12d13.5)") iq, iw, freq_r(iw), freq_r(iw+1),
c$$$     &         (trwv2(iw)+trwv2(iw+1))/2d0*(freq_r(iw)-freq_r(iw+1)),
c$$$     &         (trwv(iw)+trwv(iw+1))  /2d0*(freq_r(iw)-freq_r(iw+1))
c$$$              !weight for the histgram range. by tetwt5
c$$$              if(ixc==6)
c$$$     &         write(6,"('iq iw[min_max]=',2i5,2f7.4,' trwv by whw  = ',
c$$$     &         12d13.5)") iq,iw, freq_r(iw), freq_r(iw+1),
c$$$     &            trwv2(iw),trwv(iw) !weight for the histgram range. by tetwt5
c$$$            enddo
c$$$            deallocate(trwv,trwv2)
c$$$          endif
          if( allocated(zzr)   ) deallocate(zzr)
          if( allocated(x0mean)) deallocate(x0mean)
          if( allocated(gbvec) ) deallocate(gbvec)
        endif


C ... Close files for epsmode
        if(epsmode) then !iepsmode/=0) then      ! only calculate iq>nqibz
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
C --- RealOmega end ==============================



C --- ImagOmega ==================================
        if (imagomega) then
          write(6,*)' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi   ! if paramagnetic, multiply x0 by 2
          if (ecorr_on>0)then !ixc==101.or.(sergeyv.and.imagonly)) then
            imode=0
          else
            imode=1
          endif

!! === iw loop for imag axiw ===
          do 1016 iw  = 1,niw
c            if(newaniso2 .and. iq<=nqibz ) then
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
c              if(iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
            endif
c            if(newaniso2.and.iq>nqibz) then
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

c$$$c              read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$c              if(iw/=iwx) then
c$$$c                write(6,*)'iw iwx=',iw,iwx
c$$$c                stop 'hx0fp0: iw/=iwx' !sanity check
c$$$c              endif
c$$$cmmm3
c$$$cc                skI  (1:ngb)= zxqi(1,1:ngb,iw)
c$$$cc                sksI (1:ngb)= zxqi(1:ngb,1,iw)
c$$$c              skI  (1:ngb)= zxqi(1,1:ngb,iw)
c$$$c              sksI (1:ngb)= zxqi(2,1:ngb,iw)
c$$$              vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
c$$$              vcou1sq= sqrt(vcou1)
c$$$c              write(ifiss) iw,iq0,ngb,q
c$$$c              write(ifiss) vcou1,vcou1sq,vcousq(2:ngb),skI(1:ngb),sksI(1:ngb)
c$$$c              w_kI(1) = 0d0
c$$$c              w_ksI(1)= 0d0
c$$$c              w_kI( 2:ngb)= vcou1sq*matmul( skI(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
c$$$c              w_ksI(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), sksI(2:ngb)*vcousq(2:ngb))
c$$$cmmm epsPP mode ---> no - vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
c$$$c              llwI(iw,iq0)=  1d0 -vcou1*skI(1) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
c$$$
c$$$              llwI(iw,iq0)=  1d0 - vcou1*zxqi(1,1,iw)  !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
c$$$              write(6,*) 'iq iw llwI',iq,iw,llwI(iw,iq0)
c$$$c              write(ifisk) iw,iq0,q
c$$$c              write(ifisk) vcousq(2:ngb)*w_kI(2:ngb),vcousq(2:ngb)*w_ksI( 2:ngb)
            endif

c$$$            if(.not.newaniso2) then               ! original mode
c$$$               call rx( 'not checked here')
c$$$c              call wcf( ngb, vcoul,zxqi(1,1,iw),imode,  zw0)
c$$$            endif

c            if(newaniso2.and.iq>nqibz) then
            if(iq>nqibz) then
c              zw(1:ngb,1:ngb) = 0d0 ! zw(nblochpmx,nblochpmx)
c              write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v
            else
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
c              write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v
              write(ifrcwi, rec=iw)  zw    !  WP = vsc-v
              call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
            endif


!! --- Miyake's total energy branch !Nov2004. not maintained now... need to fix this maybe(2012takao) -----------
            if(.false.) then
c            if(ecorr_on>0 .and. (.not.newaniso2)) then !I did not modified this for newaniso2 2012takao
              if (debug) write(6,*)'ip,ix=',iq,iw,'  niw=',niw
              call getwk(iq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)
              call ecorq2 (vcoul, zw0, ngb, iq,iw,ieceig,
     o        erpaqw, trpvqw, trlogqw)
C --- Integration along imaginary axis.
! omit k and basis index for simplicity
! wint = -(i/4pi) < [w'=-inf,inf] Q(w') >
!
! When w' ==> iw', w' is now real,
!   wint =  (1/2pi) < [w'=0,inf] Q(iw') >
!
! transform: x = 1/(1+w')
! this leads to a denser mesh in w' around 0 for equal mesh x
! which is desirable since Q is peaked around w'=0
!    wint =  (1/2pi) < [x=0,1] Q(iw') / x^2 >
              faca   =  wk4ec* wiw(iw)
              trpv(iecut)    = trpv(iecut)  + faca* trpvqw
              trlog(iecut)   = trlog(iecut) + faca* trlogqw
              totexc(iecut)  = totexc(iecut)+ faca* erpaqw !  = trpv+ trlog
c            ecqw(iq,iw) = erpaqw
              if(iw==1) then
                write(ieclog,*)
              endif
              if(iw==1.and.iq==iqxini) then
                write(ieclog,
     &         "('   iq   iw   omega/i(Ry)        IntWgt',
     &         '    trpvqw(eV)      ecqw(eV)   ecqw*IntWgt',  
     &         ' :  ecut   ecuts')")
              endif
              write(ieclog,"( 2i5,3f14.6,3f14.6,2f8.3)") 
     &       iq,iw, 2d0*freq_i(iw), faca, trpvqw*hartree, erpaqw*hartree, 
     &       faca*erpaqw*hartree, ecut(iecut),ecuts(iecut)
              close(ieclog)
              open(ieclog,file="ecorr.chk",access='append')
ccccccccccccccccccccccccccccccccccccccccccccccccc
c           allocate( ovlpc(ngb,ngb),evall(ngb),
c     &          evecc(ngb,ngb))
c           evall=0d0
c           ovlpc=0d0
c           do i=1,ngb
c             ovlpc(i,i)=1d0
c           enddo
c           nmx=ngb
cc1           call diagcv(ovlpc,zw0/2d0+transpose(dconjg(zw0))/2d0,evecc,ngb, evall,nmx,1d99, nev)
c           call diagcv(ovlpc,zw0,evecc,ngb, evall,nmx,1d99, nev)
c           write(6,"('ngb nev=',2i5)") ngb,nev
c           write(6,"('chk eigen of zw0 Max Min=',2d13.6)")maxval(evall),minval(evall)
c           do i=1,3
c             write(6,*) i, evall(i)
c           enddo
c           do i=ngb-3,ngb
c            write(6,*) i, evall(i)
c           enddo
c           deallocate( ovlpc,evall,evecc)
ccccccccccccccccccccccccccccccccccccccccccccccc
            endif
 1016     continue
c          if(newaniso2) then
            deallocate(epstinv)
            if(allocated(epstilde)) deallocate(epstilde)
            if(allocated(epstilde)) deallocate(epstilde)
c          endif
c$$$          enddo
        endif
C... ImagOmega end ===============================
        if(allocated(vcoul)) deallocate(vcoul)
        if(allocated(zw0)) deallocate(zw0)
        if(allocated(zxq )) deallocate(zxq)
        if(allocated(zxqi)) deallocate(zxqi)
        if    (ixc==101.or.normalm) then
          ifrcwi = iclose('WVI.'//charnum5(iq))
        endif  
        if (normalm) then
          ifrcw  = iclose('WVR.'//charnum5(iq))
        endif
 1001 continue 
!! == end of loop 1001 for q point ==
      call MPI__barrier()


!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
      if(MPI__size/=1) then
        do iq=nqibz+1,iqxend
          iq0 = iq - nqibz
c     write(6,*)' iq iq0 mpi_rank mpi_ranktab(iq)=',iq, iq0,MPI__rank,MPI__ranktab(iq),MPI__root,nw,nw_i,niw
          if(MPI__ranktab(iq)/=0) then !jan2012
            if(MPI__ranktab(iq) == MPI__rank) then
              dest=0
              call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
              call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
            elseif(MPI__root) then
c     write(6,*)' mpi_recv iq from',iq,MPI__ranktab(iq),nw,nw_i,niw
              src=MPI__ranktab(iq)
              call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
              call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
c     do i=nw_i,nw
c     write(6,*)'recivxxx',i,llw(i,iq0)
c     enddo   
c     write(6,*)' recv llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
c     write(6,*)' recv llwI sum=',sum(abs(llwI(:,iq0))),niw
            endif
          endif
        enddo  
      endif

c$$$      deallocate( llw, llwI )
c$$$!! == generate llw and llwI ==
c$$$      ifiss=iopen('SkSks',0,-1,0)
c$$$      ifepstinv = iopen('EPS0inv',0,0,0)
c$$$      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
c$$$      read(ifepstinv) ngb
c$$$      do 1501 iq0=1,nq0i
c$$$        rewind ifepstinv
c$$$        read(ifepstinv) ngb
c$$$        allocate(vcousq(2:ngb),sk(ngb),sks(ngb),w_k(ngb),w_ks(ngb))
c$$$        do iw=nw_i,nw 
c$$$          read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$          read(ifiss) iwx,iq0x,ngb,q
c$$$          if(iw/=iwx) stop 'hx0fp0:1501 iw/=iwx'
c$$$          read(ifiss)vcou1,vcou1sq,vcousq(2:ngb),sk(1:ngb),sks(1:ngb)
c$$$          w_k(1) =0d0
c$$$          w_ks(1)=0d0
c$$$          w_k( 2:ngb)= vcou1sq*matmul( sk(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
c$$$          w_ks(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), sks(2:ngb)*vcousq(2:ngb))
c$$$          !! epsPP mode - vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
c$$$          llw(iw,iq0)=  1d0 -vcou1*sk(1) !- vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
c$$$          write(6,*) 'epsPP iq iw',iq,iw, 1d0 - fourpi* sk(1)/sum(q**2*tpioa**2)
c$$$        enddo  
c$$$        deallocate(vcousq,sk,sks,w_k,w_ks)
c$$$        allocate(vcousq(2:ngb),skI(ngb),sksI(ngb),w_kI(ngb),w_ksI(ngb))
c$$$        do iw=1,niw
c$$$          read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$          read(ifiss) iwx,iq0x,ngb,q
c$$$          if(iw/=iwx) stop 'hx0fp0:1501 iw/=iwx'
c$$$          read(ifiss) vcou1,vcou1sq,vcousq(2:ngb),skI(1:ngb),sksI(1:ngb)
c$$$          w_kI(1) = 0d0
c$$$          w_ksI(1)= 0d0
c$$$          w_kI( 2:ngb)= vcou1sq*matmul( skI(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
c$$$          w_ksI(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), sksI(2:ngb)*vcousq(2:ngb))
c$$$          !! epsPP mode ---> no - vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
c$$$          llwI(iw,iq0)=  1d0 -vcou1*skI(1) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
c$$$          write(6,*) 'iq iw llwI',iq,iw,llwI(iw,iq0)
c$$$        enddo
c$$$        deallocate(vcousq,skI,sksI,w_kI,w_ksI)
c$$$ 1501   continue

!! == W(0) divergent part and W(0) non-analytic constant part.==
!!   Note that this is only for q=0 -->iq=1
c      if(newaniso2.and.ixc==11.and.MPI__rank==0) then
      if(ixc==11.and.MPI__rank==0) then
        write(6,*)
        write(6,*)' ==== newaniso mode W(0) divergent part ==== '
!! == W(0) divergent part ==
!! getw0 routine...
        ifidmlx = iopen('EPSwklm',0,0,0)
        read(ifidmlx) nq0i,lxklm
        allocate( dmlx(nq0i,9))
        allocate( epinvq0i(nq0i,nq0i),epinv(3,3,nq0i))
        nlxklm=(lxklm+1)**2
        allocate( wklm(nlxklm))
        read(ifidmlx) dmlx, epinv,epinvq0i
        read(ifidmlx) wklm
        ifidmlx = iclose('EPSwklm') !jan2013 iclose(ifidmlx)
!! starting from llw(iw,iq0),llwI(iw,iq0)
!! === <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
        allocate(w0(nw_i:nw),w0i(niw))
c       write(6,*)' goto getw0 nq0i epinvq0i=',nq0i,epinvq0i
!! wbz(1) is the weight for q=0 = 1/(n1*n2*n3)
        call getw0(llw, nw_i,nw,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0)
        call getw0(llwI,1,niw  ,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0i)

c       write(6,*)'sumcheck w0,w0i=',sum(abs(w0)),sum(abs(w0i))
!! === w0,w0i are stored to zw for q=0 ===
!! === w_ks*wk are stored to zw for iq >nqibz ===
cmmm
        do iq = 1,1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
          q = qibze(:,iq)
c          if(iq>nqibz) then
c            iq0 = iq - nqibz
c            read(ifisk) ngb,nw_ixxx,nwxxx,niwxxx
c            allocate(vw_k(ngb),vw_ks(ngb))
c          endif
          do ircw=1,2
            if    (ircw==1) then;  nini=nw_i;   nend=nw;
              ifrcwx = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
            elseif(ircw==2) then;  nini=1;     nend=niw;
              ifrcwx = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
            endif
            do iw=nini,nend
c             if(iq<=nqibz) read(ifrcwx, rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
              read(ifrcwx, rec= iw-nini+1) zw !(1:ngb,1:ngb)
c              if( iq==1 ) then
              if(ircw==1) zw(1,1) = w0(iw)
              if(ircw==2) zw(1,1) = w0i(iw)
c$$$cmmm3 cccccccccccccccccccccccccccccccccccccccccc
c$$$          elseif( iq>nqibz ) then !-->In future, we store sperical average of zw below to zw(at q=0)===
c$$$            write(6,*)'ddd skip readin ifisk ddddddddd'
c$$$            read(ifisk) iwxx,iq0xx,qxx
c$$$            if(iwxx /=iw) stop 'iwxx/=iw'
c$$$           if(iq0xx /=iq-nqibz) stop 'iq0xx /=iq'
c$$$            if(sum(abs(qibze(:,iq)-qxx))>1d-8) stop 'sum(abs(qq-qxx))>1d-8'
c$$$            read(ifisk) vw_k(2:ngb),vw_ks(2:ngb)
c$$$            zw=0d0
c$$$            do igb1=1+1,ngb
c$$$            do igb2=1+1,ngb
c$$$              vc1vc2 = vw_ks(igb1)*vw_k(igb2)
c$$$              if(ircw==1) zw(igb1,igb2)=vc1vc2/llw(iw,iq0)
c$$$              if(ircw==2) zw(igb1,igb2)=vc1vc2/llwI(iw,iq0)
c$$$            enddo
c$$$            enddo
c$$$ccccccccccccccccccccccccccccccccccccccccccccccc
c              endif
c          if(iq==1.or.iq>nqibz) write(ifrcwx,rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
c             write(ifrcwx,rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
              write(ifrcwx,rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccmmm3
c          if(iq<=nqibz) then
c             zw=0d0
c             write(ifrcwx, rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 )  ) zw !(1:ngb,1:ngb)
c          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cooo  check write
c          if(mod(iw,3)==1) then
c           do igb1=1,ngb,23
c           do igb2=1,ngb,23
c             if(ircw==1) write(*,"('zzzwr:',4i4,2d13.5)")iq,iw,igb1,igb2,zw(igb1,igb2)
c             if(ircw==2) write(*,"('zzzwi:',4i4,2d13.5)")iq,iw,igb1,igb2,zw(igb1,igb2)
c           enddo
c           enddo
c          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            enddo
            if    (ircw==1) then
              ifrcwx = iclose('WVR.'//charnum5(iq))
            elseif(ircw==2) then
              ifrcwx = iclose('WVI.'//charnum5(iq))
            endif
          enddo
c          if(iq>nqibz) deallocate(vw_k,vw_ks)
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

!! Write TEECOR ecorr_on mode
      if(imagomega.and.ecorr_on>0) then
        hartree=2d0*rydberg()
        ifcor   = iopen('TEECORR2',1,-1,0) ! output files
        do iecut=1,necut
          write(6,"( ' RPA Ec =' 3f23.15,'   ecut ecuts (Ry)=',2d12.4)") 
     &   totexc(iecut)*hartree,trpv(iecut)*hartree, trlog(iecut)*hartree 
     &    ,ecut(iecut),ecuts(iecut)
          write(ifcor,*) '============================'
          write(ifcor,*) 'Correlation energy Erpa (eV)'
          write(ifcor,*) '============================'
          write(ifcor,*)' ### '
          write(ifcor,"(5e23.15)") 
     &     totexc(iecut)*hartree,trpv(iecut)*hartree,trlog(iecut)*hartree
     &    ,ecut(iecut),ecuts(iecut)
        enddo
C... output ecqw !    write(ifcor,*)'### ecqw(q,w) ###'
        write(ifcor,*)' nqibz =',nqibz
        write(ifcor,*)' nq0i  =',nq0i
        write(ifcor,*)' niw   =',niw
        write(ifcor,*)' --- See details of Ec in ecor.chk ---'
c        nqitot = nqibz + nq0i
c        call wecqw(ifcor,
c     d           nqibz,nqbz,nq0i,nqitot,niw,
c     o           wibz,wqt,wx,freqx,ecqw)
C... Write electron gas correlation energy
c$$$        legas = .false.
c$$$        INQUIRE (FILE = 'LEGAS', EXIST = legas)
c$$$        if(legas) then !!! test for electron gas case.
c$$$          call rx( ' LEGAS mode is not maintained well. Need some fixing.')
c$$$          voltot=0d0
c$$$          ntot=0d0
c$$$          write(6,*)' find LEGAS. legas =',legas
c$$$          iflegas = 2101
c$$$          open (iflegas,file='LEGAS')
c$$$          read(iflegas,*)rs
c$$$          close(iflegas)
c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
c$$$          qfermi = alpha/rs
c$$$          efx  = qfermi**2
c$$$          valn = efx**1.5d0*voltot/3d0/pi**2
c$$$          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
c$$$          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
c$$$          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
c$$$          if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
c$$$          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
c$$$          qfermi= dsqrt(efz)
c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
c$$$          rs    = alpha/qfermi
c$$$          write (ifcor,*)' --- electron gas ---'
c$$$          write (ifcor,*)' density parameter rs= ', rs
c$$$          write (ifcor,*)' kf= ',qfermi
c$$$          write (ifcor,*)' ### Barth-Hedin formula'
c$$$          ecelgas = eclda_bh(rs) * hartree * ntot
c$$$          write (ifcor,*)ecelgas
c$$$          write (ifcor,*)' ### Perdew-Zunger formula'
c$$$          ecelgas = eclda_pz(rs) * hartree * ntot
c$$$          write (ifcor,*)ecelgas
c$$$          write (ifcor,*)' ### Gell-Mann and Brueckner formula'
c$$$          ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
c$$$          write (ifcor,*)ecelgas
c$$$        endif
      endif
      call cputid(0)
      call MPI__Finalize
      if(ixc==11)  call rx0( ' OK! hx0fp0 mode=11 normal sergeyv')
      if(ixc==12)  call rx0( ' OK! hx0fp0 mode=12 Ecor sergeyv mode')
      if(ixc==101) call rx0( ' OK! hx0fp0 mode=101 Ecor ')
      if(ixc==202) call rx0( ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC')
      if(ixc==203) call rx0( ' OK! hx0fp0 mode=203 sergeyv eps LFC ')
      if(ixc==222) call rx0( ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv')
      if(ixc==223) call rx0( ' OK! hx0fp0 mode=223 chi+- LFC sergeyv')
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
      write(6,'(a,f8.4,2i5,4d20.12)')trim(adjustl(tagname)),freqq,iq,iw,trwv,trwv2
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

