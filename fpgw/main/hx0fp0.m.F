      program hx0fp0
      use m_readqg
      use m_readeigen,only:readeval,init_readeigen,init_readeigen2
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
C-------------------------------------------------
C  Calculate x0, W-V, \epsilon
Cr
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
Cr      The indexing starts with core first and then valence on top of core
Cr      So n-index in in for valence electron is different from inv.
C
Cr the key routine is x0k_v2.
C
C --- subrouine variables. but not perfect yet. In future, we need to make it automatic
C    i input
C    o output
C    m input and output
C
Cr comment rule
Cr C --- Chapter 
Cr C ... section
Cr ! note for what it is
Cr c Commnetout a line (this can be memo for development).
C-------------------------------------------------
      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      real(8):: ua=1d0 ! this is a dummy.
      integer(4):: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
      integer(4):: ndble=8
      integer(4):: nword

      real(8),allocatable:: vxcfp(:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     &           , ngvecp(:,:), ngvecc(:,:), !,ngpn(:),ngcni(:),iqib(:),
     &   nx(:,:),nblocha(:),lx(:),ifppb(:)   !ongveccBr(:,:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &  zw(:,:),zw0(:,:),
     &  zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
     &  ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), 
     &  qbze(:,:),qibze(:,:) !,ecore(:,:)!,rcxq(:,:,:) !sf..13June
c     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:)
     & ,rcxqmean(:,:,:,:) !takao sep2006

!  tetrahedron method
      logical :: tetra,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
      ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer(4),allocatable :: 
     &         n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),noccxvv(:)
c     &         idtetf(:,:),ib1bz(:)

      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)!,ekxx3(:,:) !,qbzw(:,: )
      real(8) ::qbzx(3),anfvec(3)
      logical :: anfexist
      logical :: debug
      integer(4),allocatable:: ibasf(:)
      real(8),allocatable :: transaf(:,:)

      real(8),allocatable :: freq_r(:),freq_i(:),freqx(:),wx(:),expa(:)
     &                      ,frhis(:)
      logical :: realomega, imagomega
c      real(8),allocatable:: ppovl(:,:)
      complex(8),allocatable:: epsi(:,:),gbvec(:),zzr(:,:),x0mean(:,:,:)
      complex(8) :: epxxx,vcmean, vcmmmm
      complex(8),allocatable:: vcmmm(:)
      character*11 fileps
      character*11 fileps23
      character*16 filepsnolfc
      character*11  filele
      logical :: paralellx0=.true. !, hist
c      character(5) :: charnum5
c      character(8) :: xt
      character(20):: xxt

! faleev
      real(8) :: Emin, Emax,emin2,emax2
c      integer :: nw2
      real(8) :: omg2max,omg1max
      real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)

! for tetwt5
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)

      integer(4),allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)

      logical::imagonly=.false.,realonly=.false. !,readgwinput
      integer(4)::iopen,maxocc2,iclose,
     & ixc,iqxini,iqxend,
     &   nw_input,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0dummy
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx(2)=-9999,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx !,  ifchipmlog

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      integer(4)::nwin, incwfin,  verbose
      real(8)::efin

      integer(4):: bzcase, ngc,mrecg
      real(8):: qq(3),quu(3), deltaq(3),qqq(3)=0d0
      logical:: omitqbz=.false., noq0p

      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)

      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg
      integer(4):: nbcut,nbcut2

      integer(4),allocatable:: nstibz(:) !Nov2004 Miyake's tote
      real(8),allocatable:: ecqw(:,:),wiw(:)
      real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree
     &   ,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
      integer(4):: iqbz,iqindx,iflegas,nmx
     &   ,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1
      real(8) :: eclda_bh,eclda_pz,wk4ec,faca
      logical ::legas

      real(8),allocatable::    evall(:)
      complex(8),allocatable:: ovlpc(:,:),evecc(:,:)
      integer(4):: nev !,  ifdpin


      real(8),allocatable:: ecut(:),ecuts(:) ,totexc(:), trpv(:),trlog(:)
      integer(4):: necut,iecut

! mode23
      integer(4):: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
      integer(4),allocatable:: nxx_r(:)
      real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
      character*3:: charnum3
      character*4:: charnum4
      complex(8),allocatable:: jcoup(:,:), mcm(:,:,:)
      complex(8):: aaax,aaax0,addx !bbbx,bbbx25,bbbx35
      real(8)::chg1,chg2,spinmom,schi=1d0
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4) :: ll
      external ll
#endif
      complex(8),allocatable:: ovlp(:,:),evec(:,:),ovlpi(:,:)
      real(8),allocatable::eval(:)
      integer(4):: new,nmxx,ii,iy,ipl1,ixx

      complex(8),allocatable :: ppovl(:,:),oo(:,:),x0meanx(:,:),x0inv(:,:)
      real(8)::qxx(3),ssm
! svd. not used now
      real(8),allocatable::SS(:),rwork(:),ss0(:)
      complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:),ddd(:,:)
     & ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
      integer(4)::lwork,info,imin,ifzxq
      complex(8)::x0mx
      complex(8),allocatable:: UU0(:,:),VT0(:,:)

      logical :: ngczero=.false., sergeyv=.true., chipm=.false.,nolfco=.false. !sergeyv only
     &  ,epsmode=.false.,normalm=.false., eiqr=.false.
      integer(4)::  ife, idum4 !ingczero,ifchipmn,ifchipm,
      real(8):: qs,qt,ww,muu, ddq(3)
      character*11 ::ttt

      logical ::smbasis
      integer(4):: nnmx,nomx

! Feb2006 time-reversal=off case
      logical :: timereversal, testtimer,onceww
      integer(4):: npm ,jpm,ncc,nw_i
      real(8):: frr

      real(8):: sciss
      integer(4):: ipm,nrecoff

      real(8),allocatable:: ebb(:)
      logical :: evaltest
      character*300:: aline
      integer(4):: istat,nmbas,imb,imb1,imb2,nmbas_in
      integer(4),allocatable:: imbas(:), imbas_s(:),iibas(:)
!...
      complex(8),allocatable:: am1(:),am2(:),mmat(:,:),
     &     x0mat(:,:),x0matinv(:,:),eiqrm(:)
      integer(4):: ifchipmn_mat, ifchipm_fmat !,ifchipm_mat
      integer(4)::ifstoner,ifx,i1
      real(8):: Istoner,zz1,zz2,zz3,zz4,Istoner0,jzero2,dumm1,dumm2
      complex(8):: trr,trr0,trr1     , zzzx(4,4), zzzy(4,4),trrx,mmatx(4,4),denom(4,4)
      real(8),allocatable:: eee(:),mmnorm(:),
     &    asvec(:,:),ssv(:,:),sproj(:,:),sprojx(:,:), momsite(:)
      real(8):: eex(4),eey(4)

      integer:: ngb_in

      hartree=2d0*rydberg()
      pi         = 4d0*datan(1d0)

!-------------------------------------------------------------------------
C --- computational mode select
      write(6,"(a)") 
     &      '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'
      write(6,"(a)") ' #1:run mode'
c      write(6,"(a)") '    1   : normal'
c      write(6,"(a)") '    2   : epsNolfc'
c      write(6,"(a)") '    3   : eps'
c      write(6,"(a)") '    22  : chi^+- NoLFC'
c      write(6,"(a)") '    23  : chi^+- '
c      write(6,"(a)") '    101 : total energy Miyake'
      write(6,"(a)") '    11  : normal    Sergey'
      write(6,"(a)") '    202 : epsNoLFC  Sergey'
      write(6,"(a)") '    203 : eps       Sergey'
      write(6,"(a)") '    222 : chi^+- NoLFC Sergey'
      write(6,"(a)") '    223 : chi^+- Sergey'
      write(6,"(a)") '    12  : total energy Miyake Sergey'
      write(6,"(a)") '    -9999: just show version num'
      write(6,"(a)") ' #2=iqxini   #3=iqxend' !    #4=ngczero or not'
c      write(6,"(a)") '  10222 : <e^{iqr}|chi^+-|e^{iqr}> NoLFC'
      write(6,"(a)") 
     &      '-------------------------------------------------------'
! Test mode   write(6,"(a)") '   tettest(4) Realtetwt4(5) Realtetwt5(6) '
      call readin6(ixc,iqxini,iqxend, idum4 ) !ingczero)
      call headver('hx0fp0',ixc)
c      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      if(iqxini ==0 ) then
        paralellx0=.false.
        write(6,"(' ixc = ',i4)") ixc
      else
        paralellx0=.true.
        write(6,"(' PARALELL mode: ixc iqxini iqxend =',3i5)")
     &  ixc, iqxini, iqxend
      endif

      if(ixc<=5.or.ixc==22.or.ixc==23.or.ixc==101) then
         write(*,*)'these modes are removed now'
         stop 'these modes are not supported'
! test mode
      elseif(ixc==6) then
        write(6,*) ' OK ixc=6. Real-axis tetwt5 mode'
! Sergey (Hilbert-transformation) modes
      elseif(ixc==11) then; write(6,*) " OK ixc=11 s.f.'s normal mode "
        normalm=.true.
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
        write(6,*) " OK ixc=12 Miyake's total energy Sergey"
        imagonly=.true.
        ecorr_on=901
!     -- chipm mode NoLFC
      elseif(ixc==222) then
        write(6,*) " OK ixc=222    chipm sergey's "
        realonly=.true.
        omitqbz=.true.
        eiqr =.true.
!     -- chipm mode NoLFC
c      elseif(ixc==10222) then
c        write(6,*) " OK ixc=10222  <q|chipm_0|q> sergey"
c        sergeyv=.true.
c        realonly=.true.
c        omitqbz=.true.
c        eiqr =.true.
!     -- eps mode with LFC
      elseif(ixc==223) then
        write(6,*) " OK ixc=223  chipm with LFC sergey's  "
        realonly=.true.
        omitqbz=.true.
        eiqr =.true.
      else
        stop' hx0fp0: mode ixc is not appropriate'
      endif

      if(ixc==202.or.ixc==203.or.ixc==222.or.ixc==223) then
        epsmode=.true. 
        if(mod(ixc,200)==22.or.mod(ixc,200)==23) chipm =.true.
        if(mod(ixc,10)==2)                       nolfco=.true.
      endif

c      if(epsmode) then
c        if(chipm) ifchipmlog = iopen ('ChiPM.log',1,3,0)
c      endif

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

C --- Readin BZDATA. See gwsrc/rwbzdata.f
!--------readin data set when you call read_BZDATA ---------------
!       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
!     &   n_index_qbz
!       integer(4):: n1,n2,n3
!       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
!       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
!     &    ,wibz(:),qbzw(:,:)
!       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
!     &    ,nstar(:),irk(:,:),index_qbz(:,:,:),nstbz(:)
!-----------------------------------------------------------------
      call read_BZDATA()


C --- Use regular mesh even for bzcase==2 and qbzreg()=T
! A little confusing...
      ddq = 0d0
      if(bzcase()==2) ddq= dq_bzcase2
      do iq = 1, nqbz
        qbz(1:3,iq) = qbz(1:3,iq) + ddq
        ! This new qbz is regular mesh, which are identical in the both bzcase.
      enddo
      if(.not.qbzreg()) then ! off-regular mesh case
        do i=1,nqbz
          qbz(:,i) = qbz(:,i) - dq_bzcase2
        enddo
      endif
      do i=1,nqbz
        write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
      enddo
      print *,' nqbz nqibz =',nqbz,nqibz
c      print *,  qbz
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

!- oct2005 not implimented cases.
      if(smbasis().and.chipm) then
        print *,' smbasis=T & chipm=T is not implimented yet.'//
     &       ' Supply consistent MixSpin for smbasis!'//
     &       ' MixSpin should be converted at the end of hvccfp0.'
        stop ' smbasis=T & chipm=T is not implimented yet.'
      endif
!
      call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
      write(6,"(' ScaledGapX0=',f5.3)") sciss

C --- readin GWIN and LMTO, then allocate and set datas.
      if(epsmode) then
        nwin = -999
      else
        nwin = 0      !Readin nw from NW file
      endif

C ... Readin by genallcf
      incwfin= 0  !use ForX0 for core in GWIN
      efin =  0d0 !readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
!...  These are allocated and setted by genallcf_v3
!      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
!     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
!      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
!      character(120):: symgrp
!      character(6),allocatable :: clabl(:)
!      integer(4),allocatable:: iclass(:)
!     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
!     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
!     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
!     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
!     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
!     &    occv(:,:,:),unoccv(:,:,:)
!     &   ,occc(:,:,:),unoccc(:,:,:),
!     o    nocc(:,:,:),nunocc(:,:,:)
!      real(8), allocatable::
!     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
!-----------------------------------------------------------------------
!      call genallcf_v2x(
!     o                   plat,alat,natom,nclass,ipos,   !structure
!     o                   symgrp,gen,ngnmax,ngrp,igrp,   !symmetry
!     o                   iinvg,ef,
!     o                   clabl,nspin,nl,nn,nnv,nnc,                  !l,n, dimensions
!     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions
!     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
!     o                   izdummy,
!     o   iil,  iin,iim,   iilnm,i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
!     o   iilv,iinv,iimv,iilnmv, i_mnlv,! l,n,m for Phi
!     o   iilc,iinc,iimc,iilnmc, i_mnlc,! l,n,m for Phi
!     o   incwf,iecore,ikonf,iicore,incore,nctot,            !core
!     o   imagw_dummy,niw,dum1,nw,dw,delta,dum2,dum3,idummy) ! frequency

      print *,' ncore=',ncore
      print *,' iclass=',iclass

! ...
      if(chipm.and.nspin==1) stop 'chipm mode is for nspin=2'
      nw_input = nw ; print *, 'nw delta=',nw_input,delta
      debug=.false.; if(verbose()>=100) debug=.true.
      if(debug) print *,' end of genallc'
!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) stop ' nclass /= natom '

C ... tetra or not
      if(delta <= 0d0) then
        tetra =  .true.
        delta = -delta
        print *,' hx0fp0: tetrahedron mode delta=',delta
      else
        tetra = .false. ! switch for tetrahedron method for dielectric functions
      endif

C ... read dimensions of h,hb
      ifhbe  = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'

C ... q near zero
      if(debug) print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      write(6,*) ' *** nqibz nq0i=', nqibz,nq0i
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
      wqtsum = sum(abs(wqt(1:nq0i)))
      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      if(normalm.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) 
     &  stop ' wqtsum of Q0P /=1 '

C --- readin by rdpp ; Radial integrals ppbrd and plane wave part
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
! qibze(3,nqbze) qbze(3,nqibze)
      nqbze  = nqbz *(1 + nq0i)
      nqibze = nqibz + nq0i
      allocate( qbze(3, nqbze), qibze(3, nqibze))
      call dcopy(3*nqbz, qbz,  1, qbze,1)
      call dcopy(3*nqibz,qibz, 1, qibze,1)
      do i = 1,nq0i
        qibze(:,nqibz+i)  = q0i(:,i)
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
        enddo
      enddo
      ngrpx = 1
      allocate( nx(0:2*(nl-1),nclass),  nblocha(nclass) ,lx(nclass),
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,ngrpx), symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      if(.not.smbasis()) nblochpmx = nbloch + ngcmx
      allocate(ngveccB(3,ngcmx)) ! work arry

C ... for legas test
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)

C --- Check PARALELL.X0. Recently, we have not checked parerell mode. So be careful to use this.
      if(chipm.and.(.not.legas).and.(.not.eiqr) ) then
        iqxend = nqibz + nq0i
      elseif(paralellx0) then
        if(iqxini >99999) stop ' hx0fp0: iqxini >99999'
        if(iqxend >99999) stop ' hx0fp0: iqxend >99999'
        ifvcfpout = iopen('VCCFP.'
     &     //xxt(iqxini,iqxend)
     &  ,0,-1,0)
        read(ifvcfpout) ndummy1, ndummy2
      else
        iqxend = nqibz + nq0i
        ifvcfpout = iopen('VCCFP',0,-1,0)
        read(ifvcfpout) ndummy1, ndummy2
      endif
      print *, ' nqibz nqibze=',nqibz,nqibze

C ... initialization of readEigen
      call init_readeigen(ginv,nspin,nband,mrece) !
      ! EVU EVD are readin in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)

C --- set realomega, imagomega tetra nw niw nwp ifgb0vec ------------
c$$$      if( ixc==101 ) then
c$$$        realomega =.false.
c$$$        imagomega =.true.
c$$$        nw = 0
c$$$        nwp= 0
c$$$      elseif( ixc==1 ) then !old imagw = 2 case
c$$$        realomega =.true.
c$$$        imagomega =.true.
c$$$      elseif( epsmode.and.(.not.sergeyv)) then !ixc==2.or.ixc==3.or.ixc==22.or.ixc==23 ) then
c$$$        print *, 'mode=',ixc
c$$$        realomega =.true.
c$$$        imagomega =.false.
c$$$        niw = 0
c$$$        call getkeyvalue("GWinput","EPSrange",epsrng )
c$$$        call getkeyvalue("GWinput","EPSdw",  dwry )
c$$$        epsrng=epsrng/2d0
c$$$        dw = dwry/2d0
c$$$        if(dw==0d0) then
c$$$          nw = 1
c$$$        else
c$$$          nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
c$$$        endif
      if(ixc==6) then
! test modes. ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.
        realomega = .true.
        imagomega = .false.
        tetra     = .true.
        niw = 0
! ---  For tetwt5 ---  the tetrahedron weight for spectrum function (imaginary part)
!   Histogram bins are specified by freq_r(1:nwp)
!     nwp=nw+1; frhis(1)=0
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].

! ... These parameters specifies a test histgram bins;Sergey's mesh just for test modes.
        nw0 = 200    !100    800
        dwh = 0.01d0 !0.02d0 0.0025d0 !in hartree
        call findemaxmin(nband,qbze,nqbze,nspin,emax,emin)
c$$$ Even if we neglect this section, result is not so bad---anyway we may need to
c$$$ fix this block if necessary. Oct2005
c$$$        if(bzcase()==2) then
c$$$          allocate(qbz2(3,nqbz))
c$$$          do iq=1,nqbz
c$$$           qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
c$$$          enddo
c$$$          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
c$$$          emax=max(emax,emax2)
c$$$          emin=min(emin,emin2)
c$$$          deallocate(qbz2)
c$$$        endif
        if (nctot .gt. 0) Emin = minval(ecore(:,1:nspin))
        omg2max = (Emax-Emin)*.5+.2d0    !(in Hartree) covers  all relevant omega, +.2 for margin
        omg1max = dwh*(nw0-1)
        nwp = int(sqrt(omg2max*(2*nw0-1d0)/dwh-(nw0**2-3*nw0+1d0)))+1 ! + 1 for margin
        nw  = nwp-1
        write(*,*) Emax,Emin,nw0,nw  ! nwp is new max number in frequency array
        write(6,'(a32,2i7,2d15.3)')'hx0fp0: nw0,nw,omg1max,omg2max='
     &              , nw0,nw,  omg1max,omg2max
        if (nw <= nw0) stop 'hx0fp0:ixc==[456] nw2 <= nw'
        allocate(freq_r(0:nw))
        do iw=0,nw  !This is a test mesh by Sergey.Faleev
          if(iw+1<=nw0) then;  freq_r(iw)=dwh*iw
          else;  freq_r(iw)=dwh*((iw+1)**2+nw0**2-3*nw0+1)/(2*nw0-1d0)
          endif
        enddo 
! freq_r(iw) is linear for iw<=nw and quadratic for nw<iw<=nw2
! freq_r(iw) chosen in such a way that it is continues with
!   nw nwp=nw+1 freq_r(1:nwp) are used after here.
        allocate(frhis(nwp))
        frhis=freq_r(0:nw)
        nwhis=nw

C ... Sergey Faleev's normal mode (Hilbert transformation for \chi.
      else
        realomega = .true.
        imagomega = .true.
        tetra     = .true.
        call getkeyvalue("GWinput","omg_c",omg_c )
        write(*,*)'dw, omg_c= ',dw, omg_c
C ... frhis
        if(bzcase()==2) then
          write(6,'(" dq_bzcase2=",3f9.4)')dq_bzcase2
        endif
c        call findemaxmin(ifev,nband,nqbz,nspin
        call findemaxmin(nband,qbze,nqbze,nspin
     o  ,emax,emin)
c$$$ Even if we neglect this section, result is not so bad---anyway we may need to
c$$$ fix this block if necessary. Oct2005
c$$$        if(bzcase()==2) then
c$$$          allocate(qbz2(3,nqbz))
c$$$          do iq=1,nqbz
c$$$           qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
c$$$          enddo
c$$$          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
c$$$          emax=max(emax,emax2)
c$$$          emin=min(emin,emin2)
c$$$          deallocate(qbz2)
c$$$        endif
        if (nctot > 0) Emin=minval(ecore(:,1:nspin))
        omg2max=(Emax-Emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
        print *,' emin emax omg2max=',emin, emax,omg2max
        nw2=int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
        allocate(freqr2(nw2))                                  !+1 b/c (iw-1)
        do iw=1,nw2
          freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
        enddo !linear + quadratic term
        if (nw2 < 2 ) stop "hx0fp0: nw2 < 2"
        if (dw*(nw-2) > freqr2(nw2-1)) stop
     &                           "hx0fp0: dw*(nw-2) > freqr2(nw2-1)"
        nwhis = nw2-1
        allocate(frhis(1:nwhis+1))
        frhis = freqr2(1:nwhis+1)
        print *,' we set frhis nwhis=',nwhis
C ... for frhis_m
        nw=nw2-1   ! nw+1 is how many points of real omega we use
        ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
        ! maximum nw=nw2-1 because nwhis=nw2-1
        do iw=3,nw2-1
          !nw is chosen from condition that frhis_m(nw-3)<dw*(nw-3) <frhis_m(nw-2).
          !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
          !nw was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
          ! for correlation Coulomb Wc(omg),
          ! and one more point omg=dw*(nw-1) needed for extrapolation.
          ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
          ! and one more point omg=frhis_m(nw) needed for extropolation
          ! used in subroutine alagr3z in  sxcf.f.
          omg2 = (freqr2(iw-2)+freqr2(iw-1))/2d0
          if (omg2 > dw*(nw_input-3)) then
            nw=iw
            exit   ! 'nw_input' is only used to get maximum frequency for
            !  dressed coulomb line
          endif
        enddo
        if(epsmode) then 
          nw  = nwhis-1
          niw = 0
        endif
        allocate(freq_r(0:nw))
        freq_r(0)=0d0
        do iw=1,nw
          freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
        enddo
        nwp = nw+1
        delta=0d0 !only for writeemesh
c        write(6,"(' ixc=11: niw nw nwhis=',3i5)") niw,nw,nwhis
      endif

! Timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
        write(6,"('TimeReversal off mode')")
        npm=2
        nw_i=-nw
        if(.not.tetra)      stop ' tetra=T for timereversal=off'
        if(chipm)           stop ' chipm=F for timereversal=off'
      endif

! ...
c$$$      if(realomega.and.(.not.sergeyv)) then !     elseif (realomega.and.(ixc==22.or.ixc==23)) then
c$$$        nwp = nw+1 !Feb2006 I want to unify the meaning of nwp as nwp=nw+1
c$$$        !With this change x0kf_v2z calculate one more freq_r along real axis.
c$$$        ! at freq_r(nw). --->It WAS not calculated. It wat up to freq_r(nw-1).
c$$$        ! But not tested well.
c$$$        allocate( freq_r(0:nw) )
c$$$        do iw = 0,nw
c$$$          freq_r(iw) = dw*iw ! freq_r
c$$$        enddo
c$$$      endif

C ... get eigenvector corresponds to exp(iqr) (q is almost zero).
      if(epsmode) then !iepsmode/=0) then ;  print *, ' read in Mix0vec'
        allocate(epsi(nw_i:nw,neps)) !5July2005 nwp should be used after it is defined!
        if(paralellx0) then
          ifgb0vec = iopen ( "Mix0vec."//xxt(iqxini,iqxend),1,3,0)
        else
          ifgb0vec = iopen ( "Mix0vec",1,3,0)
        endif
        print *, ' read end of Mix0vec'
      endif
! ...
      if(realomega) then
        open(UNIT=3111,file='freq_r') !write number of frequency
        !points nwp and frequensies in 'freq_r' file
        write(3111,"(2i8,'  !(a.u.=2Ry)')") nwp, nw_i
        do iw= nw_i,-1
          write(3111,"(d23.15,2x,i6)") -freq_r(-iw),iw !This file is reffere by hsfp0 and so.
        enddo
        do iw= 0,nw
          write(3111,"(d23.15,2x,i6)") freq_r(iw),iw !This file is reffere by hsfp0 and so.
        enddo
        close(3111)
      endif

C --- WVI only for imagonly for ixc==12
      if(imagonly) then
        realomega =.false.
        imagomega =.true.
      endif

C --- epsPP noLFC mode for ixc==13
      if(realonly) then
        realomega =.true.
        imagomega =.false.
      endif

C --- set freq_i
      if (imagomega) then
        print *,' freqimg: niw =',niw
        allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
        call freq01 (niw,ua,  ! Generate gaussian frequencies x between (0,1) and w=(1-x)/x      
     o                freqx,freq_i,wx,expa)
        allocate(wiw(niw))
        do iw=1,niw
          wiw(iw)=wx(iw)/(2d0*pi*freqx(iw)*freqx(iw))
        enddo
c        deallocate(freqx,wx,expa)! freqx,wx is used for Miyake's tote mode !Nov2004
      endif

c ... for Miyake tote mode Nov2004
      if(ecorr_on>0) then  !it was bzcase()==2 Was it bug?
        allocate(nstibz(nqibz))
        do iq=1,nqibz
          iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
          nstibz(iq) = nstbz(iqbz)
c          write(6,"(' iq qibz nstibz=',i5,3f9.4,i5)")iq,qibz(:,iq),nstibz(iq)
        enddo
      endif


C --- tetra init
      call getkeyvalue("GWinput","TmpWWK",tmpwwk,default=.false.)
      if(tetra) then
        allocate( !wgt(nband+nctot,nband,nqbz), !noccxvv(nw+niw),
     &    nbnb(nqbz,npm),nbnbtt(nqbz,npm)
     &   ,ekxx1(nband,nqbz),ekxx2(nband,nqbz) !!! nband=nlmto 
     &   )
        if(tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
      endif ;  if(debug) print *,' xxx1:'
C ... maximum no. occupied valence states
!      noccxv     = maxocc (ifev,nspin,ef,  nband,nqbze)
      noccxv     = maxocc2 (nspin,ef, nband, qbze,nqbze)
      if(noccxv>nband) stop 'hx0fp0: all the bands filled! too large Ef'
      noccx      = noccxv + nctot

C --- write energy mesh ---
      if(paralellx0) then
        ifemesh = iopen('emesh.hx0fp0.'//xxt(iqxini,iqxend),1,-1,0)
      else
        ifemesh = iopen('emesh.hx0fp0',1,-1,0)
      endif
      call writeemesh(ifemesh,freq_i,niw,freq_r,nwp,delta)
      if(debug) print *,' xxx1aax:'

C ... ppb
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )
c      nwblk  = nw /nwdiv
c      niwblk = niw/nwdiv
c      if ( nwdiv*nwblk  .ne. nw ) nwblk = nwblk + 1
c      if ( nwdiv*niwblk .ne. niw) niwblk = niwblk + 1
c      allocate(  freqr(nwblk), freqi(niwblk)   )
      nprecx = ndble  !We use double precision arrays only.

C ... This is just to get nblochpmx
      if(smbasis()) then
        call getngbpomat(nqibz+nq0i,
     o nnmx,nomx)
        nblochpmx = nnmx
      endif

      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
C --- open WV* files
      if(paralellx0) then
        ifwd   = iopen('WV.d.'//xxt(iqxini,iqxend),1,-1,0)
      else
        ifwd   = iopen('WV.d',1,-1,0)
      endif

! July 2003. now WV.d contains true size of WVR.
      write (ifwd,"(1x,10i14)") !change format oct2005
     &   nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
      allocate(  zw(nblochpmx,nblochpmx) )
      if(debug) print *,' xxx1a2:'
      if    (ixc==101) then
        if(paralellx0) then
          ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend) ,0,-1,mrecl)
        else
          ifrcwi = iopen('WVI',0,-1,mrecl)
        endif
      elseif (normalm) then !ixc==1.or.sergeyv) then
        if(paralellx0) then
          ifrcw  = iopen('WVR.'//xxt(iqxini,iqxend),0,-1,mrecl)
          ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend),0,-1,mrecl)
        else
          ifrcw  = iopen('WVR',0,-1,mrecl)
          ifrcwi = iopen('WVI',0,-1,mrecl)
        endif
c      elseif(chipm) then
c        ifrcw  = iopen('WVR',0,-1,mrecl)
      endif
      nspinmx = nspin

C --- READIN ANFcond------------------------------------
      inquire(file='ANFcond',EXIST=anfexist)
      print *,' ANFcond exist? =',anfexist
      if(anfexist) then
        allocate(ibasf(natom),transaf(1:3,natom))
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*)
        read(ifianf,*)
        read(ifianf,*) anfvec(1:3)
        write(6,'(" ANFcond=",3f13.6)') anfvec
        do ibas=1,natom
          read(ifianf,*) ibas1, ibasf(ibas)
          if(ibas1/=ibas) stop ' readin ANFcond error:ibas1/=ibas'
          write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
c     &     ,transaf(1:3,ibas)
        enddo
        close(ifianf)
        nspinmx =1
      endif ;   if(debug) print *,' xxx1b:'
C ... oct2005
      if(smbasis() .and. anfexist) then
        print *,' smbais=T & anfexist=T is not yet.'//
     &       'At leaset anfx0k shold be changed'
        stop ' smbais=T & anfexist=T is not yet.'
      endif
      if((.not.timereversal()) .and. anfexist) then
        print *,' timereversal=F & anfexist=T is not yet examined.'//
     &       'At leaset you have to consider about anfx0k and so.'
        stop ' timereversal=F & anfexist=T is not yet.'
      endif

C... these are used x0k
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2

C --- ppb ---
!      allocate(ecore(nctot,nspin)) !core energies
      do  is = 1,nspin
! <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin
        if(debug) print *," goto ppbafp is=",is
        irot = 1
        call ppbafp_v2 (irot,ngrpx,is,nspin,
     i              il,in,im, nlnm,  
     i                 nl,nn,nclass,nlnmx,
     i                 mdimx,lx,nx,nxx, 
     i                 cgr, nl-1,       
     i                 ppbrd,           
     o              ppb(1,is))
      enddo
      if(debug) print *,' xxx1c:'

C ... Set iqxini
      if(omitqbz) then
        iqxini= nqibz + 1
      elseif(paralellx0) then
     &  !skip
      elseif(bzcase()==1) then
        iqxini = 2
      else
        iqxini = 1
      endif

C... check write 1st part for Ec mode to ecorr.chk Nov2004
      if(ecorr_on>0) then
        if(ecorr_on >0.and. (.not.imagomega) ) 
     &    stop ' hx0fp0: ecorr_on .and. (.not.imagomega)'
        write(ieclog,
     &          "('   iq                q                       wk')")
        do iqq = iqxini,iqxend
          if(iqq<=nqibz) then
            qq = qibze(:,iqq)
          else
            qq = 0d0
          endif
          call getwk(iqq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)
          write(ieclog,"(i5,3x,3f12.8, f15.5)")
     &          iqq, qq, wk4ec
        enddo
        write(ieclog,*)
        write(ieclog,"('   iw omega(Ry)       wiw')")
        do iww=1,niw
          write(ieclog,"(i5, f10.5, f10.5)") iww,2d0*freq_i(iww),wiw(iww)
        enddo
        write(ieclog,*)
        write(ieclog,"(' Note:IntWgt=wk*wiw.',
     &     '  Ec =\sum_{k,iw) IntWgt(k,iw)*ecqw(k,iw)')")
        close(ieclog)
        open(ieclog,file="ecorr.chk",access='append')

        call getkeyvalue("GWinput","necut_p",necut, default=1 )
        allocate(totexc(necut),trpv(necut),trlog(necut))
        totexc = 0d0
        trpv   = 0d0
        trlog  = 0d0
c        if(.not.sergeyv) necut=1
      else
        necut=1
      endif
      allocate(ecut(necut),ecuts(necut))
      call getkeyvalue("GWinput","ecut_p" ,ecut, necut,default=(/1d10/) )
      call getkeyvalue("GWinput","ecuts_p",ecuts,necut,default=(/1d10/) )
      if(necut/=1) stop 'Sorry! this version is not for necut=1'

css
css   if(chipm.and.(.not.legas).and.(.not.eiqr)) then
      if(chipm.and.(.not.legas) ) then
c--- sanity check
        if(smbasis()) then
          print *,'This version is not for spin sus for smbais=T.'//
     &         'You need to modify MixSpin.'
          stop  'This version is not for spin sus for smbais=T'
        endif

!Sep2006
c        nmbas=1
c        allocate(imbas(nmbas),imbas_s(nmbas))
c        imbas(1)=1
c        imbas_s(1)=1

c        allocate(imbas(natom),imbas_s(natom),iibas(natom))
c        print *,' goto MagAtom read natom=',natom
c        nmbas=natom
c        istat=-9999 ! istat=-9999 means noumber of readin arguments is returened in istat.
c        call getkeyvalue("GWinput","MagAtom",
c     &        imbas,nmbas,status=istat)
c        print *,' end of MagAtom read istat=',istat
c       nmbas = natom
c        do i=1, natom
c          iibas(i) = i
c        enddo
c        do i=1, istat
c          ibas = imbas(i)
c          if(ibas<0) iibas( abs(ibas) ) = ibas
c        enddo
c        write(6,"('Readin MagAtom nmbas =',i3,' imbas_s= ',10i3)")
c     &  nmbas, imbas(1:nmbas)
c        imbas_s(1:natom) =     iibas(1:natom)
c        imbas(1:natom)   = abs(iibas(1:natom))
c        write(6,"('   Use MagAtom nmbas =',i3,' imbas_s= ',10i3)")
c     &  nmbas, imbas_s(1:nmbas)

        nmbas=natom
        allocate(imbas(nmbas),imbas_s(nmbas))
        istat=-9999 ! istat=-9999 means noumber of readin arguments is returened in istat.
        call getkeyvalue("GWinput","MagAtom", 
     &        imbas,nmbas,status=istat)
        nmbas = istat
        write(6,"('Readin MagAtom nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)
        imbas_s(1:nmbas) = imbas(1:nmbas)
        imbas(1:nmbas)   = abs(imbas(1:nmbas))
c
        allocate(jcoup(nw_i:nw,neps) )
        allocate( svec(nbloch,nmbas) )  !sep2006
        svec=0d0
        allocate( cvec(nbloch,nmbas),momsite(nmbas),
     &    mmnorm(nmbas))              !May2007
        cvec=0d0
        do imb=1,nmbas
          ibas= imbas(imb)
c          open(ifv,file='MixSpin.'//charnum3(ibas))
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
            lb = ll(ilmx )       !  print *,' lb=',lb,lxx,ilmx
            do ixx = 1, nxx_r(lb)  !  print *,' nn=',nn,nxx_r(lb)
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

! Calculate ChiPM. So sign of omega should be correct.
          if(imb==1) then !determine spin direction with respect to ibas=imbas(imb=1)
            spinmom=(chg1-chg2)
            schi=1d0
            if(spinmom<0d0) then
c              if(.not.sergeyv)  then
c                freq_r= -freq_r
c              else
                schi  = -1d0   ! This affects to dpsion. Obtained results
                ! should be the same in both mode.
c              endif
            endif
          endif

C ... ReOrdering of spinvec in nbas ordering...
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
!       enddo

c         write(ifchipmlog,"(a)")
c     &   ' --- Tabel of mcm v.s. ss cutoff.'//
c     &   '  : mcm means <m|chi^{-1}|m>. '//
c     &   '  ss is the SVD values of chi^{-1}.'
c         if(.not.nolfco) then
c          write(ifchipmlog,"(13x,a)")'q              omega'//
c     &    '   mcm(ss<ssmin+2Ry)    mcm(ss<ssmin+10Ry)'//
c     &    '    mcm(ss<ssmin+5Ry)   mcm(ss<ssmin+100Ry)'//
c     &    '    mcm(all)'
c         endif
          mmnorm (imb) = sqrt(sum(svec(:,imb)**2))
          momsite(imb) = chg1-chg2
c        write(6,"(' svecsum=',d23.15)") sum(svec(:,imb)**2
c          write(ifchipmlog,"(2d23.15,' ! mmom mmnorm')")momsite(imb),mmnorm(imb)
          write(6,"( 'mmom mmnorm= ',2f14.10)")  momsite(imb),mmnorm(imb)
        enddo
c        ifx= iclose('ChiPM.log')
      endif

! I assume 1 is for majority for eiqr case.
c      if(ix==10222) then
cc        schi=1d0 !1d0 means Majority is isp=1. If Majority is isp=2, use schi=-1d0.
cc        allocate(jcoup(nw_i:nw,neps))
c        mmnorm=1d0
c      endif

! nmbas_in is for rcxqmean
      if(chipm .and. nolfco) then !.and. ix/=10222) then
        nmbas_in = nmbas
      else
        nmbas_in = 1
      endif
      if(epsmode.and.nolfco) then
        allocate( rcxqmean(nwhis,npm,nmbas_in,nmbas_in))
        if(debug) write(6,"('fff:',3i5)") nwhis,npm,nmbas_in
      else
        allocate( rcxqmean(1,1,1,1)) !dummy
      endif

      if(chipm) allocate(eiqrm(nmbas))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cctakao
c$$$      allocate( x0meanx(nmbas,nmbas) )
c$$$      allocate( x0mat(nmbas,nmbas),x0matinv(nmbas,nmbas) )
c$$$      do 1101 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
c$$$         if(iq==iqxini+2) exit
c$$$         q = qibze(:,iq)
c$$$         print *,'aaaaaaaaaa q=',q
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
c$$$          write(6,"(' normchk=',255d23.15)") sum( dconjg(gbvec)*zzr(:,1) )
c$$$     &        ,sum(abs(gbvec(:))), sum(abs(zzr(:,1)))
c$$$            allocate(eiqrm(nmbas))
c$$$            do imb=1,nmbas
c$$$            eiqrm(imb)= sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )
c$$$            write(6,"(' <eiqr|m> =',255d23.15)") eiqrm(imb)
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
c$$$     &        '(36x,2x,20x,2x,255d23.15)') x0meanx(:,:)
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
c$$$     &        '("ttt",3f12.8,2x,f10.5,2x,2d23.15,2x,2d23.15)') q, 2*schi*frr, trr,1d0/trr
c$$$              enddo
c$$$            ifx=iclose(ifchipm2)
c$$$ 1101  continue
c$$$          stop 'xxxxxxxxxxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





C --- Calculate x0(q,iw) and W ------------------------------------------
!============ loop over q point ========================================
!=======================================================================
!      if(iepsmode==202) iqxini = nqibz+1
      do 1001 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
        call cputid (0)
        q = qibze(:,iq)
        if(iq<=nqibz.or.ngczero) then ! iqx = iq
          qq=q
        else  !    iqx = 1            ! corresponding q=0
          qq=0d0
        endif

        call readqg('QGcou', qq, ginv,  quu,ngc,ngveccB)

C --- readin Coulomb matrix
        if(chipm.and.(.not.legas).and.(.not.eiqr)) then
          allocate(vcoul(1,1))!dummy
          ngb = nbloch + ngc !ngb is readin from vcoul 25jan2006
        else
          read(ifvcfpout) ngb       ! oct2005 !I added ngb to VCCFP
          allocate(vcoul(ngb,ngb))
          read(ifvcfpout) vcoul,qx  ! coulomb matrix
css
          if(sum(q**2)<1d-10) then
            if(sum(qx**2) >1d-7)stop'qx=0 xxx See hvccfp0'
          elseif(sum(abs(q-qx))>1d-10) then
            write(6,"('iq q=',i5,3d15.7,' qx(VCCFP)=',3d15.7)")iq,q,qx
            stop 'hx0fp0: VCCFP is not compatible'
          endif
ccss
c         if(sum(abs(q-qx))>1d-10) then
c         write(6,"('iq q=',i5,3f9.4,' qx(VCCFP)=',3f9.4)")iq,q,qx
c         stop 'hx0fp0: VCCFP is not compatible'
c         endif
        endif
        write(6,"('do 1001: iq q=',i5,3f9.4,' qq=',3f9.4)")iq,q,qq

C --- eps mode
        if(epsmode) then !iepsmode/=0) then      ! only calculate iq>nqibz
          if(iq<=nqibz) then
            deallocate(vcoul)
            cycle  !Cycle here
          endif
          iqixc2 = iq- (nqibz+nq0ix)

! Readin Mix0vec for any eps mode.
!          if((.not.chipm).or.legas.or.eiqr) then
          read(ifgb0vec,*) qgbin(1:3),ngbin,igc0dummy,dnorm
          if(sum(abs(q))==0d0) then
            if(sum(qgbin**2) >1d-7)stop'qgbin=0 xxx See hvccfp0'
          elseif(sum(abs(qgbin(1:3)-q)) >1d-8)then
            stop 'qgbin inconsistent'
          endif
! zzr and gbvec
          write(6,"(' --- Readin Mix0vec: ',3d13.6,2i5,d18.8)") qgbin(1:3),ngbin
          if(ngb/=ngbin) stop 'hx0fp0: ngb/=ngbgin'
          write(6,"(' ngb nwp niw=',3i8)")ngb,nwp,niw
          allocate( gbvec(ngb),zzr(ngb,1),x0mean(nw_i:nw,1,1))
          x0mean=0d0
          do i=1,ngb
            read(ifgb0vec,"(4d24.15)") zz1,zz2,zz3,zz4
! /dnorm at 14June2008 
! gbvec and zzr do not change because we avoid /dnorm in hvccfp0.
! This measn that we changed the definition of gbvec and zzr. See hvccfp0
! dcmplx(zz3,zz4) is zzr defined in hvccfp0; this equals <eiqr|B_I>.
! NOTE: dnorm = \sum_I <eiqr|I><I|eiqr> = volume of unit cell (or sum of MT volume if no IPW).
            gbvec(i)= dcmplx(zz1,zz2)/dnorm
            zzr(i,1)= dcmplx(zz3,zz4)/dnorm
          enddo
          write(6,"(' normchk=',255d23.15)") sum( dconjg(gbvec)*zzr(:,1) )
     &        ,sum(abs(gbvec(:))), sum(abs(zzr(:,1)))

! this is for chipm
          if(chipm) then
            ii=0
            do imb=1,nmbas
              ibas= imbas(imb) !dec26 2006 bug: imb--->ibas
              if(ibas>1) i= sum(nblocha(1:ibas-1))
! *dnorm at 14June2008 
!  Before it, eiqrm is a little wrong when no PW is included (normlzimation error).
              eiqrm(imb)= sum( dconjg(gbvec(1:nbloch))*svec(1:nbloch,imb) )/mmnorm(imb) *dnorm
!  mmnorm is moved to trr
              write(6,"(' <eiqr|m>/|m| =',255d23.15)") eiqrm(imb)
              if( imbas_s(imb)<-1) eiqrm(imb)= -eiqrm(imb)
c              if( imbas_s(imb)<-1) zzr(i+1:i+nblocha(ibas), 1)= -zzr( i+1:i+nblocha(ibas), 1)
            enddo
          endif

! This override zzr
          if(chipm .and.nolfco) then
            deallocate(zzr,x0mean)
            allocate(zzr(ngb,nmbas),x0mean(nw_i:nw,nmbas,nmbas))
            x0mean=0d0
            zzr   =0d0
            zzr(1:nbloch,1:nmbas) = svec(1:nbloch,1:nmbas)
cc            print *,' zzr=cvec mode'
cc            zzr(1:nbloch,1:nmbas) = cvec(1:nbloch,1:nmbas)
          endif

          if( wqt(iq-nqibz)==0d0.and.chipm ) then
c            ifchipmn=iopen('ChiPM'//charnum4(iqixc2)//'.nlfc.dat',1,3,0)
C ... Open ChiPM*.nolfc_mat
            ifchipmn_mat=iopen('ChiPM'//charnum4(iqixc2)//'.nlfc.mat',1,3,0)
            write(ifchipmn_mat,"(255i5)") nmbas
            write(ifchipmn_mat,"(255i5)") imbas(1:nmbas)
            write(ifchipmn_mat,"(255d23.15)") momsite(1:nmbas)
            write(ifchipmn_mat,"(255d23.15)")  mmnorm(1:nmbas)
            write(ifchipmn_mat,"(255d23.15)")  eiqrm(1:nmbas)
            if(.not.nolfco) then
c              ifchipm = iopen ('ChiPM'//charnum4(iqixc2)//'.dat',1,3,0)
c              write(ifchipm,"(10x,'q(1:3)',12x,'  w(Ry)',24x
c     &        ,' <eiqr|chipm|eiqr>  ',30x
c     &        ,' <eiqr|chipm|eiqr>^{-1}   ')")
              ifchipm_fmat=iopen('ChiPM'//charnum4(iqixc2)//'.fmat',0,3,0)
              write(ifchipm_fmat) nbloch, natom,nmbas, iqxini,iqxend, nw_i,nw
              write(ifchipm_fmat) imbas(1:nmbas),momsite(1:nmbas),mmnorm(1:nmbas)
              write(ifchipm_fmat) nblocha(1:natom),svec(1:nbloch,1:nmbas)
              write(ifchipm_fmat) zzr(1:nbloch,1)
c            if(ix==10222) then
c                  write(ifchipmn,"(
c     &            'mode=10222: <eiqr|chipm0|eiqr> is calculated',
c     &            ' ... See ChiPM*.mat>')")
c            else
c                  write(ifchipmn,"(10x,'q(1:3)',12x,'  w(Ry)',24x
c     &            ,' <eiqr|chipm|eiqr>  ',30x
c     &            ,' <eiqr|chipm|eiqr>^{-1}   ')")
c            endif
            endif
          elseif(wqt(iq-nqibz)==0d0.and.(.not.chipm)) then
C ... Open EPS* file
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nlfc.dat'
            ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
            write(ifepsdatnolfc,"(a)")
     &      ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
            if(.not.nolfco) then !ixc==3.or.iepsmode==203) then
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
              ifepsdat = iopen ( fileps,1,3,0)
              write(ifepsdat,"(a)") 
     &        ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
            endif
          endif
        endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prtest projector sproj . Tested for NiO case with 4 bloch basis.
c              allocate(ssv(nbloch,nbloch),sproj(nbloch,nbloch)) !,sprojx(nbloch,nbloch))
c              ssv=0d0
c              do ix=1,nmbas
c              ssv(1:nbloch,ix)=  svec(1:nbloch,ix)/mmnorm
c              enddo
c              sproj = matmul (ssv, transpose(ssv))
cc              sprojx=-sproj
cc             do i=1,nbloch
cc                 sprojx(i,i)= 1d0 + sprojx(i,i)
cc              enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        allocate( zw0(ngb,ngb),
     &    zxq (ngb,ngb,nw_i:nw), 
     &    zxqi(ngb,ngb,niw)) !,npm))
c        write(6,*) "---  end of allocation zw0 zxq zxqi"
        if(ixc==6 ) then
          allocate( rcxq(ngb,ngb,nwp,npm) )
        elseif(epsmode.and.nolfco) then !iepsmode==202) then
          allocate( rcxq(1,1,1,npm) )
          rcxqmean=0d0
        else
          write(6,*) "rcxq alloc ngb nwhis npm ---",ngb,nwhis,npm
          allocate( rcxq(ngb,ngb,nwhis,npm) )
        endif

!! takao apr2012
        if(nolfco) then
            if(allocated(rcxq)) deallocate(rcxq)
            if(allocated(zxq) ) deallocate(zxq)
            if(allocated(zxqi) ) deallocate(zxqi)
            allocate( rcxq(nmbas_in,nmbas_in,nwhis,npm) )
            allocate( zxq (nmbas_in,nmbas_in,nw_i:nw), zxqi (nmbas_in,nmbas_in,niw))
        endif  

        zxq=0d0;  zxqi=0d0;  rcxq = 0d0
c============ loop over spin============================================
!=======================================================================
        do 1003 is = 1,nspinmx
          write(6,"(' *** ',2i4,' out of nqibz+n0qi nsp=',2i4,' *** ')") 
     &    iq, is, nqibz + nq0i,nspin
          if(debug) print *,' niw nw=',niw,nw

C ... spin chi_charge or chi_+-
          isf=is
          if(chipm) then
            write(6,*)" chi_+- mode ixc=",ixc
            if(is==1) isf=2
            if(is==2) isf=1
            rcxq=0d0
            if(epsmode.and.nolfco) rcxqmean=0d0
          endif

C ---  tetraini block tetra==.true.===============================1ini
! --- tetrahedron method --- preparatory part
!  nbnb    : total number of pairs (ib,jb) with non-zero weight.
!  ib=n1b  : band index for occ.   1\ge n1b \ge nband+nctot.
!         Valence-core order(Core index follows valence index).
!  jb=n2b  : band index for unocc. 1\ge n2b \ge nband
!  wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
          if(tetra) then
            write(6,"(' tetra mode niw nw nwhis=',3i5)") niw,nw,nwhis
C ... kqxx(kx) ekxx ---
!    ekxx1  for  rk
!    ekxx2  for q+rk.  See tetwt4
            do kx = 1, nqbz
              call readeval(    qbz(:,kx)-qqq, is, ekxx1(1:nband, kx) )
              call readeval(q + qbz(:,kx)-qqq,isf, ekxx2(1:nband, kx) )
            enddo
C ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
!     the pairs are not dependent on the energy omega
!     in the denominator of the dielectric function.
            write(6,"(' -- First tetwt5 is to get size of array --')")
            job = 0
            if(npm==1) then
              ncc=0
            else
              ncc=nctot
            endif
            allocate( demin(nband+nctot,nband+ncc,nqbz,npm),
     &                demax(nband+nctot,nband+ncc,nqbz,npm) )
            allocate( iwgt (nband+nctot,nband+ncc,nqbz,npm) )
            if(nctot==0) then
              deallocate(ecore)
              allocate(ecore(1,2)) !this is dummry
            endif
            call tetwt5x_dtet4(  npm,ncc,
     i           q, ekxx1, ekxx2, qbas,ginv,ef, 
     d          ntetf,nqbzw,  nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i          job,
     o          iwgt,nbnb,   !job=0
     o          demin,demax, !job=0
     i          frhis, nwhis,  ! job=1    not-used
     i          nbnbx,ibjb,nhwtot,  ! job=1    not-used !bug:kino-fix  nhwtot was missing though it was dummy
     i          ihw,nhw,jhw,   ! job=1    not-used
     o          whw,           ! job=1    not-used
     i        iq,is,isf,nqibz)
c            if(ixc/=4.and.ixc/=6.and.(.not.sergeyv)) deallocate(demin,demax)
            nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
            if(debug) print *,' nbnbx=',nbnbx
            allocate(  n1b(nbnbx,nqbz,npm)
     &                ,n2b(nbnbx,nqbz,npm))
            do jpm=1,npm
              call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx, !2006
     o        n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
            enddo
            if(sum(abs(nbnb-nbnbtt))/=0)stop'hx0fp0:sum(nbnb-nbnbtt)/=0'
            noccxv = maxval(noccxvx)
            noccx  = nctot + noccxv
            print *,' Tetra mode: nctot noccxv= ',nctot,noccxv
            deallocate(iwgt)
          endif
C --- end of tetraini block==========================================1end

c sergeyv=T only now feb2012
c$$$C --- TetrahedronWeight_4 block. tettwt4  ixc<=5 ==============2ini
c$$$          if(tetra.and.(.not.sergeyv)) then !(ixc==101.or.ixc<=5.or.ixc==22.or.ixc==23)) then
c$$$             allocate( wgt(nband+nctot,nband,nqbz) )
c$$$             if(.not.tmpwwk) allocate( wwk(nbnbx,nqbz,niw+nw+1)) ! +1 is for nwp case
c$$$             if(tmpwwk) rewind(ifwwk)
c$$$             nwmax = niw +nw +1
c$$$             if(ixc==4.or.ixc==5) nwmax = nw+1
c$$$             do iw = 1, nwmax
c$$$                if(ixc==4.or.ixc==5.or.ixc==6) then
c$$$                   fff = freq_r(iw-1)+img*delta
c$$$                   write(6,"(' *** iw on real= ',i3,' omega(Ry)=',2d15.6)")iw-1,fff*2d0
c$$$                elseif(iw<=niw) then
c$$$                   fff = img * freq_i(iw) ! along img axis
c$$$                   write(6,"(' *** iw on imag= ',i3,' omega(Ry)=',2d15.6)")iw,fff*2d0
c$$$                endif
c$$$             enddo
c$$$
c$$$             do 1010 iw = 1,nwmax
c$$$                if(ixc==4.or.ixc==5.or.ixc==6) then
c$$$                   fff = freq_r(iw-1)+img*delta
c$$$!     write(6,"(' *** iw on real= ',i3,' omega(Ry)=',2d15.6)")
c$$$!     &          iw-1,fff*2d0
c$$$                else
c$$$                   if(iw<=niw) then
c$$$                      fff = img * freq_i(iw) ! along img axis
c$$$!     write(6,"(' *** iw on imag= ',i3,' omega(Ry)=',2d15.6)")
c$$$!     &           iw,fff*2d0
c$$$                   else
c$$$!     nw iw-niw ---> iw-niw-1
c$$$                      fff = freq_r(iw-niw-1) + img*delta ! delta is in a.u.  ! along real axis
c$$$!     write(6,"(' *** iw on real= ',i3,' omega(Ry)=',2d15.6)") iw,fff*2d0
c$$$                   endif
c$$$                endif
c$$$
c$$$                call tetwt4x_dev(fff,  
c$$$     i               q, ekxx1, ekxx2, qbas,ginv,ef,
c$$$     d               ntetf, nqbzw, nband,nqbz,
c$$$     i               nctot,ecore(1,is),idtetf,qbzw,ib1bz, chipm, !ixc,
c$$$     o               wgt, nbnbtt,
c$$$     i       iq,is,isf,nqibz) ! new input for devided-tetrahedron method.
c$$$
c$$$!-------------------
c$$$
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$C ... Be careful about sign of spin direction.
c$$$c This may be confusing.
c$$$c Maybe better to make it clean with clearer notation.
c$$$              if(chipm.and.spinmom<0) wgt= dconjg(wgt)
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$              if(dreal(fff)==0d0) wgt = dreal(wgt) ! imag part zero clear for omg=0 !this is not necessary.
c$$$
c$$$              if(tmpwwk) then
c$$$                print *,' keyx rewind ifwwk nwmax=',iw,nbnbx,nqbz
c$$$                call rsvwwk4_w( wgt,nqbz,nband,nctot,nbnbx,
c$$$     i          n1b, n2b, nbnb, ifwwk)
c$$$              else
c$$$                call rsvwwk4( wgt,nqbz,nband,nctot,nbnbx,
c$$$     i          n1b, n2b, nbnb,
c$$$     o          wwk(1,1,iw) )
c$$$              endif
c$$$ccccccccccccccccccccccccccccc
c$$$c             write(996,"(' --- iq iw=',2i3)" ),iq,iw
c$$$c             do ik  = 1,nqbz
c$$$c             do ibib= 1,nbnb(ik)
c$$$c               write(996,"(4i4,2d16.8)" )
c$$$c     &         ik,iw,n1b(ibib,ik),n2b(ibib,ik),wwk(ibib,ik,iw)
c$$$c             enddo
c$$$c             enddo
c$$$cccccccccccccccccccccccccccc
c$$$c--------------------------
c$$$c              if( maxval(nbnb(1:nqbz,iw)) > nbnbx) then ! security check
c$$$c               print *,' hx0fp0: maxval(nbnb(1:nqbz,iw)) > nbnbx ',
c$$$c     &          nbnbx,nbnb(1:nqbz,iw)
c$$$c                stop ' hx0fp0: maxval(nbnb(1:nqbz,iw)) > nbnbx '
c$$$c              endif
c$$$c--------------------------
c$$$ccccccccccccccccccccccccccccccc
c$$$c            if(iq==2.and.iw==nwmax) write(996,*)' --- end of iq=2 ---'
c$$$c            if(iq==3.and.iw==1) then
c$$$c             write(996,"(' --- iq iw=',2i3)" ),iq,iw
c$$$c             do ik  = 1,nqbz
c$$$c             do ibib= 1,nbnb(ik)
c$$$c               write(996,"(' k ib jb=',3i4,' wwk=',2d16.8)" )
c$$$c     &         ik,n1b(ibib,ik),n2b(ibib,ik),wwk(ibib,ik,iw)
c$$$c             enddo
c$$$c             enddo
c$$$c             stop '----- test end ------------'
c$$$c           endif
c$$$ccccccccccccccccccccccccccccccc
c$$$ 1010       continue
c$$$            print *,'end of 1010 loop'
c$$$            deallocate(wgt)
c$$$            if(tmpwwk) then
c$$$              allocate( wwk(nbnbx,nqbz,nwmax) )
c$$$              print *,' wwk allocated'
c$$$              ifwwk=iclose('tmp.wwk')
c$$$c              rewind(ifwwk)
c$$$              ifwwk = iopen('tmp.wwk',0,-1,0)
c$$$              print *,' keyx rewind ifwwk nwmax=',nwmax,nbnbx,nqbz,niw,nw
c$$$              do iw=1,nwmax
c$$$                read(ifwwk) wwk(1:nbnbx,1:nqbz,iw)
c$$$              enddo
c$$$              ifwwk=iclose('tmp.wwk')
c$$$              print *,' readend of wwk'
c$$$            endif
c$$$          endif
c$$$C --- TetrahedronWeight_4 block end.========================2end
c$$$
c$$$
c$$$c==============x0kf_v2 block original(.not.sergeyv) case =====================3ini
c$$$          if(tetra.and.((.not.sergeyv).and.ixc/=4) )then !ixc<=3.or.ixc==5.or.ixc==101.or.ixc==22.or.ixc==23)) then
c$$$            if(ixc==5) then
c$$$              wwk = dimag(wwk)
c$$$            endif
c$$$            print *,' --- goto x0kf_v2z'
c$$$            if(npm==2) stop 'x0kf_v2z may not work for npm==2'
c$$$            call x0kf_v2z(   ! lessl,
c$$$     i            wwk,n1b,n2b,nbnbx,nbnb, tetra, ! for tetra=T
c$$$     i              q, is,isf, nspin, !is,  !,ifcphi ifev1,ifev2 ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
c$$$     i              qbas,ginv,ef,qbz,wbz,   !iindxk,
c$$$c     i            index_qbz,n_index_qbz, jan2004
c$$$     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, !mdim=nblocha
c$$$     i              ppb(1,is), !,pdb(1,is),dpb(1,is),ddb(1,is),
c$$$     i              freq_r,freq_i,delta,realomega,imagomega, 
c$$$     i              icore,ncore,
c$$$     &              ecore(1,is),
c$$$     d            nlmto,nqbz,noccx,noccxv,
c$$$     d            nl,nclass,natom,nnc,
c$$$     d            nlnmx,mdimx,nbloch,nwp,  niw,
c$$$c     i    iq,ngpn,ngcni(iqx),ngpmx,ngcmx  , geigB(1,1,1,is), Jan2004
c$$$c     i    ngvecpB,ngveccB, ngveccB(1,1,iqx),
c$$$     i     iq,ngc,ngpmx,ngcmx  ,
c$$$c     i    ngveccB,  !ngveccB(1,1,iqx),
c$$$     i     nqbze,nband,nqibz, 
c$$$     o     zxq,zxqi,            ! zxq... are accumulated variable for spins 
c$$$     i     nolfco,zzr,    !for iepsmode
c$$$     o     x0mean, !   )  !for iepsmode
c$$$     i     vcoul,ngb ) ! for te
c$$$          endif
c$$$          if(epsmode.and.nolfco) 
c$$$     &      print *,' nolfco --- end of x0kf_v2 sum x0mean',sum(abs(x0mean))
c$$$c=============x0kf_v2 block end===============================3end


c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
c sergeyv=T only now feb2012
c          if(ixc==4.or.sergeyv) then
C --- METHOD (tetwt5) for the tetrahedron weight
!   Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
! takao-feb/2002
            if(frhis(1)/=0d0) stop ' hx0fp0: we assume frhis(1)=0d0'
            write(6,*)' --- Frequency  histogram  sections (a.u.)---- '
            if(onceww(1)) then
              do ihis= 1, nwhis
                write(6,"(' ihis Init  End=', i4,2f13.6)")
     &        ihis,frhis(ihis),frhis(ihis+1)
              enddo
            endif
            write(6,*)' ----------------nbnbx nqbz= ',nbnbx,nqbz
C ... make index sets
            allocate(ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm))
            jhwtot = 1
            do jpm =1,npm
              do ik   = 1,nqbz
                do ibib = 1,nbnb(ik,jpm)
c             print *,' ik ibib = ',ik,ibib
cc              ib1=  n1b(ibib,ik,jpm)
cc              ib2 = n2b(ibib,ik,jpm)
c              print *,' goto histrange=',ib1,ib2
c     &         ,demin(n1b(ibib,ik),n2b(ibib,ik),ik)
c     &         ,demax(n1b(ibib,ik),n2b(ibib,ik),ik)
                  call hisrange( frhis, nwhis,  
     i         demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     i         demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     o         ihw(ibib,ik,jpm),nhw(ibib,ik,jpm))
                  jhw(ibib,ik,jpm)= jhwtot
                  jhwtot = jhwtot + nhw(ibib,ik,jpm)
                enddo
              enddo
            enddo

            nhwtot = jhwtot-1
            print *,' nhwtot=',nhwtot
            deallocate(demin,demax)
            allocate( whw(nhwtot),   ! histo-weight
     &                ibjb(nctot+nband,nband+ncc,nqbz,npm) )

            whw=0d0
            ibjb = 0
            do jpm=1,npm
              do ik   = 1,nqbz
                do ibib = 1,nbnb(ik,jpm)
                  ib1  = n1b(ibib,ik,jpm)
                  ib2  = n2b(ibib,ik,jpm)
                  ibjb(ib1,ib2,ik,jpm) = ibib
                enddo
              enddo
            enddo
C ... Generate the histogram weights whw
            job=1
c            write(6,*) 'goto tetwt5x_dtet4 job=',job
            call tetwt5x_dtet4(  npm,ncc,
     i           q, ekxx1, ekxx2, qbas,ginv,ef, 
     d          ntetf,nqbzw, nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i          job,
     o          iwgt,nbnb,      ! job=0
     o          demin,demax,   ! job=0
     i          frhis,nwhis,       ! job=1
     i          nbnbx,ibjb,nhwtot, ! job=1
     i          ihw,nhw,jhw,       ! job=1
     o          whw,               ! job=1
     i   iq,is,isf,nqibz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$            print *,' nqbz=',nqbz
c$$$            do ik   = 1,nqbz
c$$$            print *,' nbnb=',ik,nbnb(ik)
c$$$            do jpm  = 1,npm
c$$$            do ibib = 1,nbnb(ik,jpm)
c$$$             ib1  = n1b(ibib,ik,jpm)
c$$$             ib2  = n2b(ibib,ik,jpm)
c$$$             write(6,"('-sumchk whw-- jpm ik ib jb sum=',4i4,d13.6)")
c$$$     &       jpm,ik,ib1,ib2,
c$$$     &       sum( whw(jhw(ibib,ik,jpm):jhw(ibib,ik,jpm)+nhw(ibib,ik,jpm)))
c$$$c              ihww = ihw(ibib,ik)
c$$$c              do iw = 1, nhw(ibib,ik)
c$$$c               write(3666,"('jpm ihis [a b] whw =',i3,2x,i4,2f8.5,f18.12,d14.6)")
c$$$c     &         jpm, ihww+iw-1, frhis(ihww+iw-1),frhis(ihww+iw),
c$$$c     &         whw(jhw(ibib,ik)+iw-1), whw(jhw(ibib,ik)+iw-1)
c$$$c              enddo
c$$$            enddo
c$$$            enddo
c$$$            enddo
c$$$c            stop 'hx0fp0.f ---  testend xxx2'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            write(1236,"('-----',5i8)") noccxvx
c            do ik  = 1,nqbz
c              write(1236,"('--- ',5i8)") ik, nbnb(ik),nbnbtt(ik)
c            do ibib= 1,nbnb(ik)
c              write(1236,"(2i6,2f12.6)") n1b(ibib,ik),n2b(ibib,ik)
c     &         ,demin(n1b(ibib,ik),n2b(ibib,ik),ik)
c     &         ,demax(n1b(ibib,ik),n2b(ibib,ik),ik)
c            enddo
c            enddo
c            stop 'hx0fp0.m.fal.f ---  testend'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          endif
c=======TetrahedronWeight_5 block end ==================4end


c==============x0kf_v4hz block.  =========================5ini
c Instead of wwk, you have to use whw(nwhtot).
c You will need to pass these arrays to x0kf
c   ihw(ibjb,kx): omega index, to specify the section of the histogram.
c   nhw(ibjb,kx): the number of histogram sections
c   jhw(ibjb,kx): pointer to whw
c   whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
c     : histogram weights for given ib,jb,kx for histogram sections
c       from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c------------------------------------------------------------------
c          if(sergeyv) then ;  
            print *, ' --- goto x0kf_v4hz '
c- new version only along real axis for
            call x0kf_v4hz(npm,ncc,   
     i          ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf,
     i              qbas,ginv,  qbz,wbz, 
     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, 
     i              ppb(1,is), 
     i              icore,ncore,
     d            nlmto,nqbz,noccx,noccxv,
     d            nl,nclass,natom,nnc,
     d            nlnmx,mdimx,nbloch,  nwhis,
     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
     i    nqbze,nband,nqibz, 
     o    rcxq, ! See rcxq_zxq in x0kf_v3h.f; rcxq is equivalent with zxq;
            !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,nmbas_in,  !for nolfco Add nmbas Sep2006
     o    rcxqmean   ) !for nolfco
c          endif
c          print *,'ttttt: rcxqmean sum =',sum(abs(rcxqmean))

c          stop '--- test end of x0kf_v2hz ---------'
c$$$c=========== Check write mode4 block ixc==4 ============6ini
c$$$          if(ixc==4) then
c$$$C --- check write ---
c$$$            ichkhis = iopen ('tethis.chk',1,-1,0)
c$$$            do jpm   = 1,npm
c$$$              do ik   = 1,nqbz
c$$$                do ibib = 1,nbnb(ik,jpm)
c$$$                  ib1  = n1b(ibib,ik,jpm)
c$$$                  ib2  = n2b(ibib,ik,jpm)
c$$$                  ihww = ihw(ibib,ik,jpm)
c$$$                  write(ichkhis,"('--- ik ib jb jpm =',4i4)") ik,ib1,ib2,jpm
c$$$                  do iw = 1, nhw(ibib,ik,jpm)
c$$$                    write(ichkhis, "('his [min max]=',i6,2f7.4,
c$$$     &             ' whw wwk*h= ',2d12.5,'( =',d12.5,'+',d12.5,')')")
c$$$     &          ihww+iw-1, freq_r(ihww+iw-1), freq_r(ihww+iw),
c$$$!-------------------------------------------------------
c$$$     &          whw(jhw(ibib,ik,jpm)+iw-1),                     !weight for the histgram range. by tetwt5
c$$$!------------ v.s. -------------
c$$$     &          dimag(wwk(ibib,ik,ihww+iw-1)+wwk(ibib,ik,ihww+iw))/2d0
c$$$     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight=weight1+weight2 for each bin of the histgram range by tetwt4.fal
c$$$!-------------------------------------------------------
c$$$     &          dimag(wwk(ibib,ik,ihww+iw-1))/2d0
c$$$     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight1 for the histgram range by tetwt4
c$$$     &          dimag(wwk(ibib,ik,ihww+iw))/2d0
c$$$     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw))     !weight2 for the histgram range by tetwt4
c$$$                  enddo
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$            ichkhis = iclose('tethis.chk')
c$$$            print *,' OK! hx0fp1_tetwt5_test end. Histgram method.'
c$$$            print *,' The final reulst wwk vs.whw is in tethis.chk!'
c$$$            print *,' This code is a routine to test tetwt5.'
c$$$            stop    ' ixc==4  Weight Histgram check. wwk vs.whw'
c$$$          endif  ! end of if(tetra)
c$$$c=========== Check write mode4 block end  ==============6end
         if(ixc==4.or.sergeyv)         deallocate(ihw,nhw,jhw, whw,ibjb)
c$$$     if(tetra.and.(.not.sergeyv) ) deallocate(wwk) ! !ixc<=5.or.ixc==101.or.ixc==22.or.ixc==23))
         if(tetra)                     deallocate( n1b,n2b)

          iecut=1

c sergeyv=T only feb2012
c          if(sergeyv) then !This was in 1003 loop as if(ixc==11.and.is==nspinmx) then
            if(debug) write(6,"(a)") ' --- goto dpsion5 --- '
            if(is==nspinmx.or.chipm) then
              if(ixc==6) then ! ... convert from rcxq to zxq for ixc==6
                stop "--- !This path is not tested yet !!!!---"
                zxq(:,:,0)=0d0
                do jpm =1,npm
                  do iw  = 1,nw
                    do i=1,ngb
                      do j=1,i
                        ixx=1
                        if(jpm==2) ixx=-1
                        zxq(j,i,iw*ixx )=rcxq(j,i,iw,jpm)
                        zxq(i,j,iw*ixx )=dconjg(zxq(j,i,iw*ixx))
                      enddo
                    enddo
                  enddo
                enddo
              else
                write(6,'("goto dpsion5: "
     &       //"nwhis nw_i niw nw_w nmbas_in=",5i5)') nwhis,nw_i,nw,niw,nmbas_in
                write(6,*)'sumchk rcxq goto dpsion5=', sum(abs(rcxq))
c$$$             call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega,
c$$$     i        ngb, rcxq, npm,nw_i,  ! rcxq is alterd---used as work
c$$$     o        zxq, zxqi,
c$$$     i        nolfco,chipm, schi,is,  rcxqmean,ecut(iecut),ecuts(iecut), !iepsmode,
c$$$     o        x0mean)  !zxq x0mean is accumlating veriable for chipm mode


c rcxqmean --->rcxq index conversion
               if(nolfco) then
                do iw=1,nw
                do jpm=1,npm
                  rcxq(:,:,iw,jpm)=rcxqmean(iw,jpm,:,:)
                enddo
                enddo
               endif
               ngb_in=ngb
               if(nolfco) ngb_in=nmbas_in

               call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i         ngb_in, rcxq, npm,nw_i,  ! rcxq is alterd---used as work
     o         zxq, zxqi, 
c    i         nolfco,chipm, schi,is,  rcxqmean,nmbas_in,ecut(iecut),ecuts(iecut),
     i         chipm, schi,is,  ecut(iecut),ecuts(iecut))
c     o          x0mean)  !zxq x0mean is accumlating veriable for chipm mode
               if(nolfco) then
                 do iw=nw_i,nw
                  x0mean(iw,:,:)=zxq(:,:,iw)
                 enddo  
               endif
               
              endif 
            endif  ;  print *,' --- end of dpsion5 ----'
            if(epsmode.and.nolfco) print *,' end of dpsion5 sum x0mean',sum(abs(x0mean))
c          endif
 1003   continue  !end of spin loop =====
        if(allocated(rcxq) ) deallocate(rcxq)
 1033   continue
        if(debug)  print *," end of 1033 xxxxxxxxx"
!=======================================================================



c==============x0kf_v3h block end  =========================
c$$$- 13July2005
c$$$    Takao needed to remove do 2005 loop so as to put dpsion2 into do 1003
c$$$    loop so as to make chipm mode works ok(do dpsion2 for eachs spin).
c$$$    If necessary we may need to recover do 2005 mode---this was for e_rpa test...
c$$$C... this loop is to calculate Ec for different ecut, ecuts.
c$$$        do 2005 iecut=1,necut !In usual mode, necut=1
c$$$c-----------------
c$$$        if(sergeyv) then !This was in 1003 loop as if(ixc==11.and.is==nspinmx) then
c$$$          if(debug) print *, ' --- goto dpsion2'
c$$$          call dpsion2(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega,
c$$$     i      ngb, rcxq,
c$$$     o      zxq, zxqi, ! zxq... are Not the accumulating variable for spins
c$$$     i      nolfco,chipm,is,  rcxqmean,ecut(iecut),ecuts(iecut), !iepsmode,
c$$$     o      x0mean)
c$$$          print *,' --- end of dpsion2'
c$$$          if(epsmode) print *,' end of dpsion sum x0mean',sum(abs(x0mean))
c$$$c          do iw=1,nwp
c$$$c            write(6,"(a,i5,2d13.5)")' sumcheck zxq=',iw,sum(abs(zxq(:,:,iw)))
c$$$c          enddo
c$$$        endif

c       write(ifdpin)  nwhis,nw,niw, ngb, frhis,freq_r, freq_i,rcxq


C ---  RealOmega ====================================
        if (realomega) then
          if(chipm) then !ixc==22.or.ixc==23) then
            if (nspin==1) stop 'chipm modes are for nspin==2'
            if (anfexist) then
! Not tested yet
!          if(anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
              write(6,*) 'I have not tested anfexist=t and ixc==22 23 '//
     &      '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
              stop 'hx0fp0: tail output! Not yet for anfexist&ixc==22 23'
            endif
          elseif(epsmode.and.nolfco.and.(.not.chipm)) then !ixc==2.or.iepsmode==202) then
            if (nspin==1) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
            if (anfexist) then
! Not tested yet
!          if(anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
              write(6,*) 'I have not tested anfexist=t and ixc==2 '//
     &     '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
              stop 'hx0fp0: tail output! Not yet for anfexist&ixc==2'
            endif
          else
            if (nspin == 1) zxq = 2d0*zxq   !if paramagnetic, multiply x0 by 2
            if (anfexist) then   !Antiferro accelaration.
c           nwmax = nw
c           if(ixc==4.or.ixc==5.or.sergeyv) nwmax=nw
              do iw  = nw_i,nw !max
                call anfx0k(natom,nclass,nblocha,iclass,pos
     i     ,nbloch,ngc, !ngcni(iqx),
c     i    q, ngveccB(1,1,iqx),qbas, ! for q+G
     i     q, ngveccB,qbas, ! for q+G
     i     anfvec,ibasf,             ! these are antiferro informations.
     i     zxq(1,1,iw))              ! i/o
              enddo
            endif
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
            endif
c          elseif(ixc==1.or.sergeyv) then
            write(6, *)" trace check for W-V"
          endif

c C ... Get O^{1/2}= sqrt(ovlp)
c          if(chipm.and.(.not.nolfco)) then
c            allocate(sqovlp(ngb,ngb),sqovlpi(ngb,ngb))
c            write(6,*)"goto getsqovlp"
c            call getsqovlp(q,ngc,ngb, sqovlp)
c            sqovlpi=sqovlp
c            write(6,*)"goto matcinv"
c            call matcinv(ngb,sqovlpi)
c            allocate(UU0(nbloch,nbloch),VT0(nbloch,nbloch)
c     &     ,SS0(nbloch),mmat(nbloch,nbloch)
c     &     ,zzz(nbloch,nbloch))
c          endif
c          write(6,*)"zzzz nw_i nw=",nw_i,nw

          do iw  = nw_i,nw  !Feb2006. Before it was 1:nwp (nwp=nw+1).
            !  So freq_r(iw-1) is shifted to freq_r(iw).
cnw all freq is from zero freq(0:nwp-1); nwp is the number of real flequencies.
            frr= dsign(freq_r(abs(iw)),dble(iw))
            if(.not.epsmode) then  !if(ixc==1.or.sergeyv) then
              imode = 1
              if(debug) write(6,"(' rrr sum zxq=',i5,4d23.14)") iw,sum(zxq(:,:,iw)),sum(abs(zxq(:,:,iw)))

! a debug mode.
              if(evaltest()) then
                write(6,"('hhh --- EigenValues for Im( xq) --------')")
                allocate(ebb(ngb))
                call diagcvh2((zxq(:,:,iw)-transpose(dconjg(zxq(:,:,iw))))/2d0/img, ngb, ebb)
                do ii=1,ngb
                  if( abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
                    write(6, "('hhhIxq : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
                  else
                    write(6, "('hhhIxqxxx : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
                  endif
                enddo
                deallocate(ebb)
                write(6,*) "test end of h0000000000000000000"
              endif

C ... wcf
              call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                zw0)

!... a debug mode
              if(evaltest()) then
                write(6,"('hhh --- EigenValues for Im( W) --------')")
                allocate(ebb(ngb))
                call diagcvh2( (zw0-transpose(dconjg(zw0)))/2d0/img, ngb, ebb)
                do ii=1,ngb
                  if( abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
                    write(6, "('hhhIWq : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
                  else
                    write(6, "('hhhIWqxxx : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
                  endif
                enddo
                deallocate(ebb)
                write(6,*) '---Stop here: end of evaletes---'
                if(iw==nw)  stop  '---Stop here: end of evaletes---'
              endif

              if(debug) write(6,"(' rrr sum zw0=',i5,4d23.14)") iw,sum(zw0),sum(abs(zw0))
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
              write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v

              if(debug) print *,' wvr3 nrec=',(iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1,
     &                nblochpmx,iq-iqxini,iw,nw
              !----------------------------
              call tr_chkwrite("freq_r iq iw realomg trwv=", 
     &               zw, iw, frr,nblochpmx, nbloch,ngb,iq)
              !----------------------------
            elseif(epsmode.and.(.not.chipm)) then !.and.ixc/=23) then ! No LFC (local field correction). It's better to use echo 4| hbasfp0.
              vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
              if(.not.nolfco) then
                x0mean(iw,1,1) = sum( dconjg(zzr(:,1))* matmul(zxq(:,:,iw),zzr(:,1)))
              endif
              epsi(iw,iqixc2) = 1d0/(1- vcmean * x0mean(iw,1,1))
              write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2d23.15,3x, 2d23.15,
     &        " vcmean x0mean =", 2d23.15,3x, 2d23.15)') iqixc2,iw,2*frr,
     &        1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean,x0mean(iw,1,1)
              write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &        q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              if(.not.nolfco) then ! With LFC
                imode=2
                call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                zw0) !  write(6,"('ssschk1=',3d13.5)") sum(abs(zw0)) sum(abs(gbvec))
                epsi(iw,iqixc2)= sum( dconjg(gbvec) * matmul(zw0,zzr(:,1)) )
                write(6,'( " iq iw omega eps epsi  wLFC="
     &          ,2i6,f8.3,2d23.15,3x, 2d23.15)') 
     &          iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                write(6,*)
                write(ifepsdat,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &          q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              endif

C --- ChiPM mode
            elseif(epsmode.and.chipm) then
              allocate( x0meanx(nmbas,nmbas) )

              if(nolfco) then  ! ChiPM mode without LFC
                if(legas) then
! --- three lines below may work for test purpose for legas. But not sure.
c                  vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) )
c                  write(ifchipmn,'(3f12.8,2x,f8.5,2x,2d23.15)') 
c     &            q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw,1,1)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
                else
                  x0meanx = x0mean(iw,:,:)/2d0  !in Ry unit.
                endif
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

              write(ifchipmn_mat,
     &        '(3f12.8,2x,f20.15,2x,255d23.15)')q, 2*schi*frr, x0meanx(:,:)
              if(.not.nolfco) 
     &          write(ifchipm_fmat) q, 2*schi*frr, zxq(1:nbloch,1:nbloch,iw)

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
c!                 write(ifx,"(d23.15)")  jzero2
c!                 do imb=1,nmbas !temporary
c!                   write(ifx,"(d23.15,' ! tttt temporary... U_mm in eV')") 
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
c!      &        '(3f12.8,2x,f20.15,2x,2d23.15,2x,2d23.15)') q, 2*schi*frr, trr,1d0/trr
c!               deallocate( x0mat,x0matinv)
c! 
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
c! cxxxx this SVD procedure is not used now.
c! c!      SVD of chi^-1: !now only look for lowest eigenvalue problem... So rather eigenvalue problem instead of SVD
c! c                  write(6,"(a,i5)")' ----SVD: chiinv --- iw=',iw
c! c                  zxq(1:nbloch,1:nbloch,iw)=zzz
c! c                  call zgesvdnn(
c! c     i            nbloch, zxq(1:nbloch,1:nbloch,iw),
c! c     o            SS0,UU0,VT0)
c! c                  Istoner = -sum(UU0(:,1)*VT0(1,:))/ss0(1)
c! c!                  do ii= 1,nbloch
c! c!                    write(ifstoner,'(4d23.15)') UU0(ii,1),VT0(1,ii)
c! c!                  enddo
c!                   write(ifstoner,"(d23.15)") 
c!      &             Istoner
c!                   do imb=1,nmbas
c!                     write(ifstoner,"(d23.15,'!tttt temporary U_mm in eV')") 
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
c! c                write(6,"(' <eiqr| 1 + I chi0^+-|eiqr> =',255d23.15)") trr0
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
c!      &         '(3f12.8,2x,f20.15,2x,2d23.15,2x,4d23.15)')
c!      &              q,   2*schi*frr,  trr, 1d0/trr
c!               endif

              deallocate(x0meanx)
            endif

          enddo

c          if(chipm.and.(.not.nolfco)) 
c     &          deallocate(sqovlp,sqovlpi,UU0,VT0,SS0,mmat,zzz)
          if( allocated(zzr)   ) deallocate(zzr)
          if( allocated(x0mean)) deallocate(x0mean)
          if( allocated(gbvec) ) deallocate(gbvec)
          if( ixc==5.or.ixc==6 ) then
            jpm=1
c            nwmax = nw
c            if(ixc==5) nwmax =nw
            allocate(trwv(nw_i:nw),trwv2(nw_i:nw))
            do iw = nw_i,nw !max ! trace check
              trwv(iw) = zxq(6,7,iw)
              trwv2(iw) = 0d0
              do i = 1,ngb
                trwv2(iw) = trwv2(iw) + zxq(i,i,iw)
              enddo
            enddo
            do iw= nw_i,nw-1
              if(ixc==5) 
     &         write(6,"('iq iw[min_max]=',2i5,2f7.4,' trwv by wwk*h= ',
     &         12d13.5)") iq, iw, freq_r(iw), freq_r(iw+1), 
     &         (trwv2(iw)+trwv2(iw+1))/2d0*(freq_r(iw)-freq_r(iw+1)),
     &         (trwv(iw)+trwv(iw+1))  /2d0*(freq_r(iw)-freq_r(iw+1))
              !weight for the histgram range. by tetwt5
              if(ixc==6) 
     &         write(6,"('iq iw[min_max]=',2i5,2f7.4,' trwv by whw  = ',
     &         12d13.5)") iq,iw, freq_r(iw), freq_r(iw+1),
     &            trwv2(iw),trwv(iw) !weight for the histgram range. by tetwt5
            enddo
            deallocate(trwv,trwv2)
          endif
        endif

C ... Close files for epsmode
        if(epsmode) then !iepsmode/=0) then      ! only calculate iq>nqibz
          if(chipm) then
c            ifchipmn    =iclose('ChiPM'//charnum4(iqixc2)//'.nlfc.dat')
            ifchipmn_mat=iclose('ChiPM'//charnum4(iqixc2)//'.nlfc.mat')
            if(.not.nolfco) then
c              ifchipm    = iclose ('ChiPM'//charnum4(iqixc2)//'.dat')
              ifchipm_fmat=iclose( 'ChiPM'//charnum4(iqixc2)//'.fmat')
c             ifchipm_mat= iclose( 'ChiPM'//charnum4(iqixc2)//'.mat')
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



C --- ImagOmega ===================================
        if (imagomega) then
          print *,' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi   ! if paramagnetic, multiply x0 by 2
c ... Antiferro accelaration.
          if (anfexist) then
c           do ipm  = 1,npm
            do iw  = 1,niw
              call anfx0k(natom,nclass,nblocha,iclass,pos
     i       ,nbloch,ngc, !ngcni(iqx),
c     i      q, ngveccB(1,1,iqx),qbas,    ! for q+G
     i       q, ngveccB,qbas,    ! for q+G
     i       anfvec,ibasf,  ! these are antiferro informations.
     m       zxqi(1,1,iw)) !,ipm))  ! i/o
            enddo
c           enddo
          endif
          if (ecorr_on>0)then !ixc==101.or.(sergeyv.and.imagonly)) then
            imode=0
          else
            imode=1
          endif

c$$$          do ipm=1,npm
          do iw  = 1,niw
c            print *,'sss sumchk zxqi=',sum(abs(zxqi(:,:,iw))) !,ipm)))
            call wcf( ngb, vcoul,zxqi(1,1,iw),imode,  !,ipm),imode, 
     o              zw0)
c            print *,'sss sumchk zw0=',sum(abs(zw0))
            zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
            write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v

C --- for Miyake's total energy !Nov2004 ----------------
            if(ecorr_on>0) then
              if (debug) write(*,*)'ip,ix=',iq,iw,'  niw=',niw
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
            call tr_chkwrite("freq_i iq iw imgomg  trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
          enddo
c$$$          enddo
        endif
C... ImagOmega end ===============================
c$$$ 2005 continue ! iecut loop


c 1002 continue  ! end of frequency block-loop
        deallocate(vcoul, zw0)
        if(allocated(zxq )) deallocate(zxq)
        if(allocated(zxqi)) deallocate(zxqi)
 1001 continue
C--- end of loop over q point =================================
c=======================================================================
 1100 continue


!--- legas mode is not working now. Need fixing... voltot ntot are not given.
      if(epsmode.and.legas) then
        stop ' LEGAS mode is not maintained well. Need some fixing.'
        voltot=0d0
        ntot=0d0
        print *,' Find LEGAS. legas =',legas
        iflegas = 2101
        open (iflegas,file='LEGAS')
        read(iflegas,*)rs
        close(iflegas)
        alpha  = (9*pi/4d0)**(1d0/3d0)
        qfermi = alpha/rs
        efx  = qfermi**2
        valn = efx**1.5d0*voltot/3d0/pi**2
        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
        write (6,*)'     Exact Fermi momentum  qf  =', qfermi
        write (6,*)'     Exact Fermi energy    Ef  =', efx
        do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
          if(iq<=nqibz) cycle
          print *,' iq=',iq
          iqixc2 = iq- (nqibz+nq0ix)
          filele ='EPSEG'//charnum4(iqixc2)//'.dat'
          ife = iopen ( filele,1,3,0)
          write(ife,"(a)")
     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
          q = qibze(:,iq)
          qt= sqrt(sum(qibze(1:,iq)**2))*2d0*pi/alat
          qs= qt/qfermi
          write(6,"(' qs qfermi=',2d13.5)"    ) qs,qfermi
          write(6,"(' q-q^2/2 q+q^2=',2d13.5)") qs-qs**2/2d0,qs+qs**2/2d0
          do iw  = nw_i,nw
            ww  = freq_r(iw)
            muu = ww/qfermi**2
            if(     qs<2d0 .and. muu < qs-qs**2/2d0) then
              x0mx= -img*qfermi/(4*pi*qs)*2*muu
            elseif( qs<2d0 .and. muu < qs+qs**2/2d0) then
              x0mx= -img*qfermi/(4*pi*qs)*( 1d0-(muu/qs-.5d0*qs)**2 )
            else
              x0mx=0d0
            endif
            vcmmmm= 4*pi/qt**2
            epsi(iw,iqixc2) = 1d0/(1- vcmmmm * x0mx)
c            epsi(iw,iqixc2) = 1d0/(1- vcmmm(iq) * x0meanx)
            write(ife,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &        q, 2*ww,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
          enddo
        enddo
        print *,' ----------legas end--------'
      endif

C--- Write TEECOR ecorr_on mode
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
          write(ifcor,*)' *** '
          write(ifcor,"(5d23.15)") 
     &     totexc(iecut)*hartree,trpv(iecut)*hartree,trlog(iecut)*hartree
     &    ,ecut(iecut),ecuts(iecut)
        enddo
C... output ecqw !    write(ifcor,*)'*** ecqw(q,w) ***'
        write(ifcor,*)' nqibz =',nqibz
        write(ifcor,*)' nq0i  =',nq0i
        write(ifcor,*)' niw   =',niw
        write(ifcor,*)' --- See details of Ec in ecor.chk ---'
c        nqitot = nqibz + nq0i
c        call wecqw(ifcor,
c     d           nqibz,nqbz,nq0i,nqitot,niw,
c     o           wibz,wqt,wx,freqx,ecqw)
C... Write electron gas correlation energy
        legas = .false.
        INQUIRE (FILE = 'LEGAS', EXIST = legas)
        if(legas) then !!! test for electron gas case.
        stop ' LEGAS mode is not maintained well. Need some fixing.'
        voltot=0d0
        ntot=0d0
          print *,' find LEGAS. legas =',legas
          iflegas = 2101
          open (iflegas,file='LEGAS')
          read(iflegas,*)rs
          close(iflegas)
          alpha = (9*pi/4d0)**(1d0/3d0)
          qfermi = alpha/rs
          efx  = qfermi**2
          valn = efx**1.5d0*voltot/3d0/pi**2
          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
          if(tetra) stop 'legas You have to give ef of  tetrahedron'
          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
          qfermi= dsqrt(efz)
          alpha = (9*pi/4d0)**(1d0/3d0)
          rs    = alpha/qfermi
          write (ifcor,*)' --- electron gas ---'
          write (ifcor,*)' density parameter rs= ', rs
          write (ifcor,*)' kf= ',qfermi
          write (ifcor,*)' *** Barth-Hedin formula'
          ecelgas = eclda_bh(rs) * hartree * ntot
          write (ifcor,*)ecelgas
          write (ifcor,*)' *** Perdew-Zunger formula'
          ecelgas = eclda_pz(rs) * hartree * ntot
          write (ifcor,*)ecelgas
          write (ifcor,*)' *** Gell-Mann and Brueckner formula'
          ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
          write (ifcor,*)ecelgas
        endif
      endif
C...close files
c        close(ifdpin)
      is = iclose('hbe.d')
      call cputid(0)
      if(ixc==101) stop ' OK! hx0fp0 mode=101 Ecor '
c      if(ixc==1)   stop ' OK! hx0fp0 mode=1 normal'
c      if(ixc==2)   stop ' OK! hx0fp0 mode=2  eps NoLFC'
c      if(ixc==22)  stop ' OK! hx0fp0 mode=22 chi+- NoLFC'
c      if(ixc==23)  stop ' OK! hx0fp0 mode=23 chi+- LFC  '
c      if(ixc==3)   stop ' OK! hx0fp0 mode=3 eps with LFC'
      if(ixc==202) stop ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC'
      if(ixc==203) stop ' OK! hx0fp0 mode=203 sergeyv eps LFC '
      if(ixc==222) stop ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv'
c      if(ixc==10222)stop' OK! hx0fp0 mode=10222 <q|chi+-|q> NoLFC sergeyv'
      if(ixc==223) stop ' OK! hx0fp0 mode=223 chi+- LFC sergeyv'
      if(ixc==11)  stop ' OK! hx0fp0 mode=11 normal sergeyv'
      if(ixc==12)  stop ' OK! hx0fp0 mode=12 Ecor sergeyv mode'
c      if(ixc==4)   stop ' OK! hx0fp0 ixc=4 test'
c      if(ixc==5)   stop ' OK! hx0fp0 ixc=5 test'
c      if(ixc==6)   stop ' OK! hx0fp0 ixc=6 test'
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
      write(ifcor,*)'*** ecqw(q,w) ***'
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
        write(ifcor,*)'*** iq,wq = ',ip,wk

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
C- Get ppovl
      implicit none
      real(8)::q(3)
      integer(4):: ngc,ngb,nbloch,i,nmxx,ix,iy,nev
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
      print *,' diagcv overlap ngc nev=',ngc,nev
      deallocate(ooo,ppo)
c
 888  continue
      sqovlp=0d0
      do i=1,nbloch
        sqovlp(i,i)=1d0
      enddo
      do i=1,ngc
        if(eval(i)<0d0) then
          stop 'getsqovlp:  eval(i) <0d0'
        endif
        do ix=1,ngc;  do iy=1,ngc
          sqovlp(ix+nbloch,iy+nbloch)= 
     &     sqovlp(ix+nbloch,iy+nbloch)
     &      + evec(ix,i)* sqrt(eval(i))* dconjg(evec(iy,i))
        enddo ;      enddo
      enddo
      if(allocated(evec)) deallocate(evec)
      if(allocated(eval)) deallocate(eval)
      print *,' end of getsqovlp'
c         sqovlpi = sqovlp
c         call matcinv(ngb,sqovlp)     !  inverse
c         ovlpi=ovlp
c         deallocate(ppovl,ovlp)
      end
c---------------------------------------------------------------------
      subroutine readppovl0(q,ngc,ppovl)
      implicit none
      integer(4):: ngc_r,ippovl0,ixx,ngc
      complex(8):: ppovl(ngc,ngc)
      real(8):: q(3),qx(3)
      ippovl0=2301
      open(ippovl0,file='PPOVL0',form='unformatted')
      do
        read(ippovl0) qx,ngc_r
        if(sum(abs(qx-q))<1d-6) then
c           print *,' q qx ngc ngc_r=',q, qx, ngc, ngc_r
          if(ngc_r/=ngc) stop 'readin ppovl: ngc_r/=ngc'
          read(ippovl0) ppovl
c           print *,' readin PPOVL0 OK q=',q
          exit
        endif
      enddo
      close(ippovl0)
      end
c--------------------------------------------------------------------
      subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
      implicit none
      integer(4):: nblochpmx,nbloch,ngb,iw,i,iq
      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
      real(8):: freqq
      logical :: smbasis
      character*(*)::tagname
      trwv=0d0
      if(.not.smbasis()) then
        do i = 1,nbloch
          trwv = trwv + zw(i,i)
        enddo
      endif
      trwv2 = 0d0
      do i = 1,ngb
        trwv2 = trwv2 + zw(i,i)
      enddo !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
c      do i = 1,ngb
c        write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c      enddo
      end
c--------------------------------------------------------------------
C      subroutine test_xxx(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
C      implicit none
C      integer(4):: nblochpmx,nbloch,ngb,iw,i,iq
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
C- Test routine for Inversion
      implicit none
      integer(4) :: info,n,n2(2)
      integer(4),allocatable :: ipiv(:)
      complex(8):: a(:,:), b(1)
      complex(8),allocatable:: work(:)
      n2= SHAPE(a)
      n=n2(1)
      call zcopy(n,b,1,a,1)
      call zgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matcinv: zegtrf info=',info
        stop    ' matcinv: zegtrf '
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matcinv: zegtri info=',info
        stop    'matcinv: zegtri '
      endif
      end
c--------------------------------------------------------------------
      subroutine diagno00(nbloch,wpvc,eval)
      implicit none
      integer(4):: nbloch,nmx,nev,i
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
