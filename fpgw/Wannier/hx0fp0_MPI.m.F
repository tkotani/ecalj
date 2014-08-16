      program hx0fp0_MPI
      use rsmpi
      use rsmpi_qkgroup
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
c-------------------------------------------------
ctm, Takashi Miyake, 0810;  cRPA for entangled bands
c  Calculate x0, W-V, \epsilon
Cr
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic). 
Cr      The indexing starts with core first and then valence on top of core 
Cr      So n-index in in for valence electron is different from inv. 
C
Cr the key routine is x0k_v2.
c
c-------------------------------------------------
      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      real(8):: ua=1d0 ! this is a dummy.

c worksize in megabytes (1 word = 4 bytes)
c      integer(4),parameter :: mbytes=50, mwords=mbytes/4,
c     &  iwksize=mwords * 1000 *1000
c      integer w
c      common /w/ w(iwksize)
c      integer(4),parameter:: ngnmax=10,mxclass=100
c      real(8):: gen(9,ngnmax)
c      character*6 clabl(mxclass)
c      character*120 symgrp

      integer(4):: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
      integer(4):: ndble=8
      integer(4):: nword

c> for genallc.f
c      real(8):: plat(3,3) !,qbas(3,3),ginv(3,3)
c     l          ifphiv(2),ifphic(2),
c     l          iphi(2),iphidot(2),ipgs(2),ipgsdot(2),
c     l          iphiv(2),iphivd(2),ipgv(2),ipgvd(2),
c     l          iphic(2),ipgc(2)
c
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
      complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:),rcxqmean(:,:) !sf..13June

c  tetra 
      logical :: tetra,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
                                    ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer(4),allocatable :: 
     &         n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),noccxvv(:)
c     &         idtetf(:,:),ib1bz(:) 
c 
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
      complex(8),allocatable:: epsi(:,:),gbvec(:),zzr(:),x0mean(:)
      complex(8) :: epxxx,vcmean, vcmmmm
      complex(8),allocatable:: vcmmm(:)
      character*9 fileps
      character*11 fileps23
      character*15 filepsnolfc
      character*11  filele
      logical :: paralellx0=.true. !, hist
c      character(5) :: charnum5
c      character(8) :: xt
      character(20):: xxt

c faleev 
      real(8) :: Emin, Emax,emin2,emax2          
c      integer :: nw2
      real(8) :: omg2max,omg1max
      real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)
ctm, 080213
      real(8) :: emax_freqr,emin_freqr

c for tetwt5
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)

      integer(4),allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      
      logical::imagonly=.false.,realonly=.false. !,readgwinput
      integer(4)::iopen,maxocc2,iclose,
     & ixc,iqxini,iqxend,
c,nw,incwf,natom,nclass,ipos,igrp,
c     & iinvg,nspin,nl,nn,nnv,nnc,
c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions 
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c     o                   izdummy, 
c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
c     o   iecore,ikonf,iicore,incore,nctot,             !core
c     o   imagw_dummy,niw,idummy,
     &   nw_input,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx(2)=-9999,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx,  ifchipmlog

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      integer(4)::nwin, incwfin,  verbose
      real(8)::efin

c for ixc==13
c      integer(4) ::iepsmode=0
c      logical :: ifgeigb

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

c mode23
      integer(4):: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
      integer(4),allocatable:: nxx_r(:)
      real(8),allocatable:: svec(:),spinvec(:,:)
      character*3:: charnum3
      character*4:: charnum4
      complex(8),allocatable:: jcoup(:,:), mcm(:)
      complex(8):: aaax,aaax0,addx !bbbx,bbbx25,bbbx35
      real(8)::chg1,chg2,spinmom,schi
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

      complex(8),allocatable :: ppovl(:,:),oo(:,:)
      real(8)::qxx(3),mmnorm,ssm
c svd
      real(8),allocatable::SS(:),rwork(:),ss0(:)
      complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:)
     & ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:),sqovlp(:,:),sqovlpi(:,:)
      integer(4)::lwork,info,imin,ifzxq
      complex(8)::ddd, x0meanx
      complex(8),allocatable:: UU0(:,:),VT0(:,:)

      logical :: ngczero=.false., sergeyv=.false., chipm=.false.,nolfco=.false.
     &  ,epsmode=.false.,normalm=.false., eiqr=.false.
      integer(4)::  ife, ifchipmn,ifchipm,isvd,idum4 !ingczero,
      real(8):: qs,qt,ww,muu, ddq(3)
      character*11 ::ttt

      logical ::smbasis
      integer(4):: nnmx,nomx

C... Feb2006 time-reversal=off case
      logical :: timereversal, testtimer,onceww
      integer(4):: npm ,jpm,ncc,nw_i
      real(8):: frr


c RS: variables for MPI
      integer(4) :: input4(4),iq_local

      complex(8),allocatable :: rcxq_mpi(:,:,:,:,:),rcxqmean_mpi(:,:,:)

ctm, 060517
      logical :: lmlw=.false.
      integer(4):: nwf,ifwd_full,ifrcw_full,ifrcwi_full,
     &             nprecx2,mrecl2,nblochpmx2,nwp2,niw2,nqibze2,
     &             nw_i2,nrec
      complex(8),allocatable:: zw_full(:,:),zw_full2(:,:)
ctm, 080311
      logical :: lwbnd=.false.
      integer(4),allocatable:: iwbnd(:,:,:)
ctm, 080721
      logical :: lent=.false.
      integer(4),allocatable:: ib_ent(:,:,:)

c RS: program begins here
      call RSMPI_Init()

      hartree=2d0*rydberg()
      pi         = 4d0*datan(1d0)

      if (Is_IO_Root_RSMPI()) then
!-------------------------------------------------------------------------
C --- computational mode select
         write(6,"(a)") 
     &        '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'  
         write(6,"(a)") ' #1:run mode'
         write(6,"(a)") '    1   : normal' 
         write(6,"(a)") '    2   : epsNolfc'
         write(6,"(a)") '    3   : eps'
         write(6,"(a)") '    22  : chi^+- NoLFC'
         write(6,"(a)") '    23  : chi^+- '
         write(6,"(a)") '    101 : total energy Miyake'
         write(6,"(a)") '    11  : normal    Sergey'
         write(6,"(a)") '    202 : eps       Sergey'
         write(6,"(a)") '    203 : eps NoLFC Sergey'
         write(6,"(a)") '    222 : chi^+- NoLFC Sergey'
         write(6,"(a)") '    223 : chi^+- Sergey'
         write(6,"(a)") '    12  : total energy Miyake Sergey'
         write(6,"(a)") '    -9999: just show version num'
         write(6,"(a)") ' #2=iqxini   #3=iqxend' !    #4=ngczero or not'
         write(6,"(a)") '  10222 : <e^{iqr}|chi^+-|e^{iqr}> NoLFC'
         write(6,"(a)") 
     &        '-------------------------------------------------------'
c     Test mode   write(6,"(a)") '   tettest(4) Realtetwt4(5) Realtetwt5(6) '
         call readin6(ixc,iqxini,iqxend, idum4 ) !ingczero)
         call headver('hx0fp0_RSMPI',ixc)
         input4(1)=ixc
         input4(2)=iqxini
         input4(3)=iqxend
         input4(4)=idum4
      endif
      call MPI_Bcast(input4,4,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(input4)",ierror_rsmpi)

      ixc    = input4(1)
      iqxini = input4(2)
      iqxend = input4(3)
      idum4  = input4(4)
      

c      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      if(ixc==0) then
         call RSMPI_Stop( ' --- ixc=0 --- Choose computational mode!')
      endif

      if(iqxini ==0 ) then
        paralellx0=.false.
        if (Is_IO_Root_RSMPI()) write(6,"(' ixc = ',i4)") ixc
      else
c RS: paralellx0=.true. is not implemented...
c        paralellx0=.true.
c        write(6,"(' PARALELL mode: ixc iqxini iqxend =',3i5)")
c     &  ixc, iqxini, iqxend
         call RSMPI_Stop("iqxini!=0 :
     & This PARALELL mode is not implemented")
      endif

! original mode
c RS: only ixc==11,202,203 are available now...
! Sergey (Hilbert-transformation) modes
      if(ixc==11) then
         if (Is_IO_Root_RSMPI()) then
            write(6,*) " OK ixc=11 s.f.'s normal mode " 
         endif
         sergeyv=.true.         !ixc=11
         normalm=.true.
!     -- eps mode NoLFC
      elseif(ixc==202) then
         if (Is_IO_Root_RSMPI()) then
            write(6,*) " OK ixc=202  sergey's eps mode Only NoLFC "
         endif
        sergeyv=.true. !ixc=11
        realonly=.true.
c        iepsmode=202
        omitqbz=.true.
!     -- eps mode with LFC 
      elseif(ixc==203) then
         if (Is_IO_Root_RSMPI()) then
            write(6,*) " OK ixc=203 sergey's eps mode with LFC "
         endif
         sergeyv=.true.         !ixc=11
         realonly=.true.
c        iepsmode=203
        omitqbz=.true.
! Total energy modes
      else
         call RSMPI_Stop("hx0fp0_RSMPI: Only ixc=11,202,203 is allowed")
      endif
!
      if(    ixc==  2.or.ixc==  3.or.ixc== 22.or.ixc== 23
     &   .or.ixc==202.or.ixc==203.or.ixc==222.or.ixc==223
     &   .or.ixc==10222) then
        epsmode=.true. !        iepsmode = ixc
      endif


      if(epsmode) then
c RS: chipm=.true. is not implemented...
c        if(mod(ixc,200)==22.or.mod(ixc,200)==23) chipm =.true.
        if(mod(ixc,10)==2)                       nolfco=.true.
c        if(chipm) ifchipmlog = iopen ('ChiPM.log',1,3,0)
      endif
c
c      call getkeyvalue("GWinput","Ecorr",ecorr_on,default=0)


c     RS: commented out the following block because ecorre_on = -1
c      if(ecorr_on > 0) then
c
c         ieclog = 8155
c         if(ecorr_on==901) then
c           ieceig=8156
c           open(ieceig,file='rpa_eigen.chk') 
c           close(ieceig,status='delete')
c         endif
c         open(ieclog, file='ecorr.chk')
c      endif


c        ifdpin = iopen('DPSIONin',0,-1,0) 

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz  
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:),nstbz(:)
c-----------------------------------------------------------------
      call read_BZDATA()

c--- Use regular mesh even for bzcase==2 and qbzreg()=T
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
      if (Is_IO_Root_RSMPI()) then
         do i=1,nqbz
            write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
         enddo
         write(6,*) ' nqbz nqibz =',nqbz,nqibz
      endif

c      print *,  qbz 
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

C- oct2005 not implimented cases.
c RS: smbasis=T is not implemented...
      if(smbasis()) then
         call RSMPI_Stop('smbasis=T is not implemented yet!')
c        if(.not.sergeyv) 
c     &  stop ' smbasis=T & sergeyv=F is not implimented yet'
c        if(chipm) then
c        print *,' smbasis=T & chipm=T is not implimented yet.'//
c     &       ' Supply consistent MixSpin for smbasis!'//
c     &       ' MixSpin should be converted at the end of hvccfp0.'
c        stop ' smbasis=T & chipm=T is not implimented yet.'
c        endif
      endif

C--- readin GWIN and LMTO, then allocate and set datas.
c      call wkinit (iwksize)
      if(epsmode.and.sergeyv) then
         nwin = -999
      elseif(ixc==101.or.ixc==1.or.sergeyv) then
        nwin = 0      !Readin nw from NW file
      else
        nwin = -999   !Not readin NW file
      endif
c      if(iepsmode==202.or.iepsmode==203) nwin = -999

      incwfin= 0  !use ForX0 for core in GWIN
      efin =  0d0 !readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
c      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
      if(ngrp/= ngrp2) then
         call RSMPI_Stop('ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      endif

c---  These are allocated and setted by genallcf_v3
c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr   
c      character(120):: symgrp
c      character(6),allocatable :: clabl(:)
c      integer(4),allocatable:: iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)   
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c      real(8), allocatable:: 
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
c-----------------------------------------------------------------------

ctm, 060517
        if (niw==0) realonly = .true.

c RS: chipm == F... 
c      if(chipm.and.nspin==1) stop 'chipm mode is for nspin=2'

c      nw_input = nw ; print *, 'nw delta=',nw_input,delta
      nw_input = nw
      if (Is_IO_Root_RSMPI()) write(6,*) 'nw delta=',nw_input,delta

cc   if (nclass > mxclass) stop ' hx0fp0: increase mxclass'
      debug=.false.; if(verbose()>=100) debug=.true.

c RS: debug-mode is disabled...
      if(debug) print *,' end of genallc'
      if(debug) call RSMPI_Stop("sorry, debug mode is disabled.")

!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
c      if(nclass /= natom) stop ' nclass /= natom '
      if(nclass /= natom) call RSMPI_Stop(' nclass /= natom ')

c ... tetra or not
      if(delta <= 0d0) then
        tetra =  .true.
        delta = -delta
c        print *,' hx0fp0: tetrahedron mode delta=',delta
        if (Is_IO_Root_RSMPI()) write(6,*) 
     &       ' hx0fp0: tetrahedron mode delta=',delta
      else
        tetra = .false. ! switch for tetrahedron method for dielectric functions
      endif


c ... read dimensions of h,hb
      ifhbe      = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
c RS: close hbe.d
      is = iclose('hbe.d')
      if(nlmto/=nlmtot) call RSMPI_Stop(' hx0fp0: nlmto/=nlmtot in hbe.d')
      if(nqbz /=nqbzt ) call RSMPI_Stop(' hx0fp0: nqbz /=nqbzt  in hbe.d')

c --- q near zero --------
c      if(debug) print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      if (Is_IO_Root_RSMPI()) write(6,*) ' *** nqibz nq0i=', nqibz,nq0i

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
      close(101)
      if (Is_IO_Root_RSMPI()) then
         write( 6,*) ' num of zero weight q0p=',neps
         write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      endif

      wqtsum = sum(abs(wqt(1:nq0i))) 
      if (Is_IO_Root_RSMPI())
     & call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      call MPI_Bcast(noq0p,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(noq0p)",ierror_rsmpi)
c      if(normalm.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) 
c     &  stop ' wqtsum of Q0P /=1 '
      if(normalm.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p))  then
         call RSMPI_Stop(' wqtsum of Q0P /=1 ')
      endif

c --- readin by rdpp ; Radial integrals ppbrd and plane wave part
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
! --- qibze(3,nqbze) qbze(3,nqibze)
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
     
c... for legas test 
      legas = .false.
c RS: legas = .true. is not implemented..
c      INQUIRE (FILE = 'LEGAS', EXIST = legas)

c --- Check PARALELL.X0
      if(chipm.and.(.not.legas).and.(.not.eiqr) ) then
        iqxend = nqibz + nq0i
c RS:
        call RSMPI_Stop("This cannot happen")
      elseif(paralellx0) then
        if(iqxini >99999) stop ' hx0fp0: iqxini >99999'
        if(iqxend >99999) stop ' hx0fp0: iqxend >99999'
        ifvcfpout = iopen('VCCFP.'
     &     //xxt(iqxini,iqxend)
     &  ,0,-1,0) 
        read(ifvcfpout) ndummy1, ndummy2
c RS:
        call RSMPI_Stop("This cannot happen")
      else
        iqxend = nqibz + nq0i
        ifvcfpout = iopen('VCCFP',0,-1,0) 
        read(ifvcfpout) ndummy1, ndummy2
      endif
      if (Is_IO_Root_RSMPI()) write(6,*) ' nqibz nqibze=',nqibz,nqibze

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
                                                  ! EVU EVD are readin in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)

c --- set realomega, imagomega tetra nw niw nwp ifgb0vec ------------
      if( ixc==101 ) then
         call RSMPI_Stop("This cannot happen")
         realomega =.false.
         imagomega =.true.
         nw = 0
         nwp= 0
      elseif( ixc==1 ) then !old imagw = 2 case
         call RSMPI_Stop("ixc==1, this cannot happen")
        realomega =.true.
        imagomega =.true.
      elseif( epsmode.and.(.not.sergeyv)) then
c!ixc==2.or.ixc==3.or.ixc==22.or.ixc==23 ) then
         call RSMPI_Stop("epsmode=T and sergeyv=F, this cannot happen")
         print *, 'mode=',ixc
         realomega =.true.
         imagomega =.false.
         niw = 0
         if (Is_IO_Root_RSMPI()) then
         call getkeyvalue("GWinput","EPSrange",epsrng )
         call getkeyvalue("GWinput","EPSdw",  dwry )
         endif
      call MPI_Bcast(epsrng,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(epsrng)",ierror_rsmpi)
      call MPI_Bcast(dwry,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(dwry)",ierror_rsmpi)
         epsrng=epsrng/2d0
         dw = dwry/2d0
         if(dw==0d0) then
            nw = 1
         else
            nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
         endif
      elseif(ixc==4.or.ixc==5.or.ixc==6) then 
         call RSMPI_Stop("ixc==4 or 5 or 6, this cannot happen")

! ... These are test modes.
! ixc=4 tetrahedren weight test. tetwt5.vs.tetwt5. Write tethis.chk
! ixc=5 Spectrum function (Img part) along the Real axis with tetwt4
! ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.

c RS: rest of this block is removed..


C ... Sergey Faleev's normal  mode
      elseif(sergeyv) then 
        realomega = .true.
        imagomega = .true.
        tetra     = .true.
        if (Is_IO_Root_RSMPI()) then
        call getkeyvalue("GWinput","omg_c",omg_c )
        write(6,*)'dw, omg_c= ',dw, omg_c
        endif ! Is_IO_Root_RSMPI
        call MPI_Bcast(omg_c,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_Bcast(omg_c)",ierror_rsmpi)
C ... frhis        
        if(bzcase()==2) then
           if (Is_IO_Root_RSMPI()) write(6,'(" dq_bzcase2=",3f9.4)')dq_bzcase2
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
ctm, 080213 >
        if (Is_IO_Root_RSMPI()) then
        call getkeyvalue("GWinput","emax_freqr",emax_freqr,default=-999d0)
        call getkeyvalue("GWinput","emin_freqr",emin_freqr,default=999d0)
        if (emax_freqr.gt.emin_freqr) then
           if (emax_freqr.gt.emax) emax = emax_freqr
           if (emin_freqr.lt.emin) emin = emin_freqr
           write(*,*)'hx0fp0: emax/emin from GWinput'
           write(*,*)'emax=',emax
           write(*,*)'emin=',emin
        endif
        endif ! Is_IO_Root_RSMPI
        call MPI_Bcast(emax,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_Bcast(emax)",ierror_rsmpi)
        call MPI_Bcast(emin,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_Bcast(emin)",ierror_rsmpi)
ctm, 080213 <
        omg2max=(Emax-Emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
        if (Is_IO_Root_RSMPI())
     &    write(6,*) ' emin emax omg2max=',emin, emax,omg2max
        nw2=int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
        allocate(freqr2(nw2))                                  !+1 b/c (iw-1)
        do iw=1,nw2
          freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
        enddo !linear + quadratic term
c        if (nw2 < 2 ) stop "hx0fp0: nw2 < 2"
        if (nw2 < 2 ) call RSMPI_Stop( "hx0fp0: nw2 < 2")
c        if (dw*(nw-2) > freqr2(nw2-1)) stop
c     &                           "hx0fp0: dw*(nw-2) > freqr2(nw2-1)"
        if (dw*(nw-2) > freqr2(nw2-1)) then
           if (Is_IO_Root_RSMPI()) then
             write(*,*)'dw=',dw 
             write(*,*)'nw=',nw 
             write(*,*)'nw2=',nw2
             write(*,*)'omg_c=',omg_c
             write(*,*)'omg2max,',omg2max
             write(*,*)'freqr2(nw2-1)',freqr2(nw2-1)
             write(*,*)'recommended nw<',int(freqr2(nw2-1)/dw) +1
!             do iw=1,nw2
!               write(*,*)iw,freqr2(iw)
!             enddo
           endif ! Is_IO_Root_RSMPI
           call RSMPI_Stop("hx0fp0: dw*(nw-2) > freqr2(nw2-1)")
        endif ! dw*(nw-2)

        nwhis = nw2-1
        allocate(frhis(1:nwhis+1)) 
        frhis = freqr2(1:nwhis+1)
c        print *,' we set frhis nwhis=',nwhis
        if (Is_IO_Root_RSMPI()) write(6,*) ' we set frhis nwhis=',nwhis
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
        if(epsmode.and.sergeyv) then !epsmode==202.or.iepsmode==203) then
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
      endif

C ... timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
        write(6,"('TimeReversal off mode')")
        npm=2
        nw_i=-nw
c        if(.not.tetra)      stop ' tetra=T for timereversal=off'
        if(.not.tetra) call RSMPI_Stop(' tetra=T for timereversal=off')
c        if(.not.sergeyv)    stop ' sergeyv=T for timereversal=off'
        if(.not.sergeyv) call RSMPI_Stop(' sergeyv=T for timereversal=off')
        if(chipm)  call RSMPI_Stop(' chipm=F for timereversal=off')
      endif

! ...
      if(realomega.and.(.not.sergeyv)) then !     elseif (realomega.and.(ixc==22.or.ixc==23)) then
        call RSMPI_Stop("realomega=T and sergeyv=F: this cannot happen")
c RS: removed
      endif

!...get eigenvector corresponds to exp(iqr) (q is almost zero).
      if(epsmode) then !iepsmode/=0) then ;  print *, ' read in Mix0vec'
        allocate(epsi(nw_i:nw,neps)) !5July2005 nwp should be used after it is defined!
        if(paralellx0) then
          call RSMPI_Stop("paralellx0=T: this cannot happen")
          ifgb0vec = iopen ( "Mix0vec."//xxt(iqxini,iqxend),1,3,0)
        else
          ifgb0vec = iopen ( "Mix0vec",1,3,0) 
        endif
c        print *, ' read end of Mix0vec'
        if (Is_IO_Root_RSMPI()) write(6,*) ' read end of Mix0vec'
      endif
! ...
      if(realomega) then
c RS: only IO_Root (rank=0) opens output file "freq_r".
         if (Is_IO_Root_RSMPI()) then
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
      endif  

C --- WVI only for imagonly for ixc==12
      if(imagonly.and.sergeyv) then
        call RSMPI_Stop("imagonly=T and sergeyv=T: this cannot happen")
        realomega =.false.
        imagomega =.true.
      endif

C --- epsPP noLFC mode for ixc==13
      if(realonly.and.sergeyv) then
        realomega =.true.
        imagomega =.false.
      endif

C --- set freq_i 
      if (imagomega) then
         if (Is_IO_Root_RSMPI()) write(6,*)' freqimg: niw =',niw

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
         call RSMPI_Stop("ecorr_on > 0 : this cannot happen")
c RS: removed
      endif

ctm, 080721
        nwf =  1 ! temporary
        if (Is_IO_Root_RSMPI())
     &   call getkeyvalue("GWinput","crpa_ent",lent,default=.false.)
        call MPI_Bcast(lent,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_Bcast(lent)",ierror_rsmpi)
        if (lent) then 
          if (Is_IO_Root_RSMPI())write(*,*)'cRPA_ent mode on'
          allocate(ib_ent(nband,nqbze,nspin))
          ib_ent = 0
          call init_readeigen_crpa_ent(ib_ent,nwf)
        else ! lent
ctm, 080311
        if (Is_IO_Root_RSMPI())
     &   call getkeyvalue("GWinput","crpa_band",lwbnd,default=.false.)
        call MPI_Bcast(lwbnd,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_Bcast(lwbnd)",ierror_rsmpi)
          if (lwbnd) then
            if (Is_IO_Root_RSMPI())write(*,*)'cRPA_band mode on'
            call get_nwf_MPI(nwf)
            allocate(iwbnd(nwf,nqbze,nspin))
            call choose_wanband_MPI(iwbnd,nwf,nqbze,nspin)
          else ! lwbnd  
            if (Is_IO_Root_RSMPI())write(*,*)'cRPA_window mode on'
          endif ! lwbnd
        endif ! lent
      if (.not.lwbnd) allocate(iwbnd(nwf,nqbze,nspin))
      if (.not.lent) allocate(ib_ent(nband,nqbze,nspin))

C --- tetra init
      if (Is_IO_Root_RSMPI())
     & call getkeyvalue("GWinput","TmpWWK",tmpwwk,default=.false.)
       call MPI_Bcast(tmpwwk,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
       call RSMPI_Check("MPI_Bcast(tmpwwk)",ierror_rsmpi)
      if(tetra) then
        allocate( !wgt(nband+nctot,nband,nqbz), !noccxvv(nw+niw),
     &    nbnb(nqbz,npm),nbnbtt(nqbz,npm)
     &   ,ekxx1(nband,nqbz),ekxx2(nband,nqbz) !!! nband=nlmto 
     &   )
        if(tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
      endif ;  if(debug) print *,' xxx1:'
c ... maximum no. occupied valence states
c      noccxv     = maxocc (ifev,nspin,ef,  nband,nqbze)
      noccxv     = maxocc2 (nspin,ef, nband, qbze,nqbze)
      if (Is_IO_Root_RSMPI())write(*,*)'ef=',ef
c      if(noccxv>nband) stop 'hx0fp0: all the bands filled! too large Ef'
      if(noccxv>nband) then
         call RSMPI_Stop('hx0fp0: all the bands filled! too large Ef')
      endif
      noccx      = noccxv + nctot

c --- write energy mesh ---
      if(paralellx0) then
        call RSMPI_Stop("paralellx0==T: this cannot happen")
        ifemesh = iopen('emesh.hx0fp0.'//xxt(iqxini,iqxend),1,-1,0)
      else
c RS: Only IO root opens emesh.hx0fp0
         if (Is_IO_Root_RSMPI()) then
            ifemesh = iopen('emesh.hx0fp0',1,-1,0)
         endif
      endif
      if (Is_IO_Root_RSMPI()) then
c RS: Only IO root opens emesh.hx0fp0
         call writeemesh(ifemesh,freq_i,niw,freq_r,nwp,delta)
         is = iclose("emesh.hx0fp0")
      endif

c      if(debug) print *,' xxx1aax:'

c ... ppb
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )
c      nwblk  = nw /nwdiv
c      niwblk = niw/nwdiv
c      if ( nwdiv*nwblk  .ne. nw ) nwblk = nwblk + 1
c      if ( nwdiv*niwblk .ne. niw) niwblk = niwblk + 1
c      allocate(  freqr(nwblk), freqi(niwblk)   )
      nprecx = ndble  !We use double precision arrays only.

C... This is just to get nblochpmx
c RS: removed 
c      if(smbasis()) then
c      endif

      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
C --- open WV* files
c RS: only io_root_rsmpi opens WV.d
      if(paralellx0) then
         call RSMPI_Stop("parallelx0=T,this cannot happen")
         ifwd   = iopen('WV.d.'//xxt(iqxini,iqxend),1,-1,0)
      else
         if (Is_IO_Root_RSMPI()) then
            ifwd   = iopen('WV.d',1,-1,0)
         endif
      endif

c! July 2003. now WV.d contains true size of WVR.
      if (IS_IO_Root_RSMPI()) then
         write (ifwd,"(1x,10i14)") !change format oct2005
     &        nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
c RS: close WV.d
         is = iclose('WV.d')
      endif

      allocate(  zw(nblochpmx,nblochpmx) )



c RS: opens WVI and WVR
c this part is moved (just before do(1001)-loop)

c      if(debug) print *,' xxx1a2:'
c      if    (ixc==101) then
c         if(paralellx0) then
c           ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend) ,0,-1,mrecl)
c         else
c           ifrcwi = iopen('WVI',0,-1,mrecl)
c         endif
c         call RSMPI_Stop("ixc==101, this cannot happen")
c      elseif (normalm) then !ixc==1.or.sergeyv) then
c         if(paralellx0) then
c           ifrcw  = iopen('WVR.'//xxt(iqxini,iqxend),0,-1,mrecl)
c           ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend),0,-1,mrecl)
c           call RSMPI_Stop("paralellx0==T, this cannot happen")
c         else
c           ifrcw  = iopen('WVR',0,-1,mrecl)
c           ifrcwi = iopen('WVI',0,-1,mrecl)
c           ifrcw  = iopen('WVR.RSMPI'//qkgroup_id,0,-1,mrecl)
c           ifrcwi = iopen('WVI.RSMPI'//qkgroup_id,0,-1,mrecl)
c         endif
c      elseif(chipm) then
c        ifrcw  = iopen('WVR',0,-1,mrecl)
c      endif

      nspinmx = nspin 

c --- READIN ANFcond------------------------------------
      inquire(file='ANFcond',EXIST=anfexist)
      if (Is_IO_Root_RSMPI()) then
         write(6,*) ' ANFcond exist? =',anfexist
      endif

      if(anfexist) then
        allocate(ibasf(natom),transaf(1:3,natom))
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*)
        read(ifianf,*)
        read(ifianf,*) anfvec(1:3)
        if (Is_IO_Root_RSMPI()) then
           write(6,'(" ANFcond=",3f13.6)') anfvec
        endif
        do ibas=1,natom
          read(ifianf,*) ibas1, ibasf(ibas) 
c          if(ibas1/=ibas) stop ' readin ANFcond error:ibas1/=ibas'
          if(ibas1/=ibas) 
     &         call RSMPI_Stop(' readin ANFcond error:ibas1/=ibas')

          if (Is_IO_Root_RSMPI()) then
             write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
          endif
c     &     ,transaf(1:3,ibas)
        enddo
        close(ifianf)
        nspinmx =1
      endif ;   if(debug) print *,' xxx1b:'
C... oct2005
      if(smbasis() .and. anfexist) then

         write(6,*) ' smbais=T & anfexist=T is not yet.'//
     &        'At leaset anfx0k shold be changed'
         call RSMPI_Stop("smbasis==T and anfexist==T")
c       stop ' smbais=T & anfexist=T is not yet.'
      endif
      if((.not.timereversal()) .and. anfexist) then
         write(6,*) ' timereversal=F & anfexist=T is not yet examined.'//
     &        'At leaset you have to consider about anfx0k and so.'
         call RSMPI_Stop("!timereversal and anfexist")
c       stop ' timereversal=F & anfexist=T is not yet.'
      endif

C... these are used x0k
      if (Is_IO_Root_RSMPI()) then
         call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
         call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
         write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
      endif
      call MPI_Bcast(nbcut,1,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(nbcut)",ierror_rsmpi)
      call MPI_Bcast(nbcut2,1,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(nbcut2)",ierror_rsmpi)
C --- ppb ---
!      allocate(ecore(nctot,nspin)) !core energies
      do  is = 1,nspin
! <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin 
        if(debug) print *," goto ppbafp is=",is
        irot = 1
        call ppbafp_v2 (irot,ngrpx,is,nspin,
     i              il,in,im, nlnm,  
     d                 nl,nn,nclass,nlnmx,
     i                 mdimx,lx,nx,nxx, 
     i                 cgr, nl-1,       
     i                 ppbrd,           
     o              ppb(1,is)) 
      enddo
      if(debug) print *,' xxx1c:'

c... Set iqxini 
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
c RS: since ecorr_on = -1 ,this part was removed
      if(ecorr_on>0) then
c
      else
        necut=1
      endif 
      allocate(ecut(necut),ecuts(necut))
      if (Is_IO_Root_RSMPI()) then
      call getkeyvalue("GWinput","ecut_p" ,ecut, necut,default=(/1d10/) )
      call getkeyvalue("GWinput","ecuts_p",ecuts,necut,default=(/1d10/) )
      endif ! Is_IO_Root_RSMPI
      call MPI_Bcast(ecut,necut,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(ecut)",ierror_rsmpi)
      call MPI_Bcast(ecuts,necut,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(ecuts)",ierror_rsmpi)
      if(necut/=1) then
         call RSMPI_Stop('Sorry! this version is not for necut=1')
      endif
! I assume 1 is for majority for eiqr case.
c RS: chipm == F.
c      if(chipm.and.eiqr) then
c        schi=1d0 !1d0 means Majority is isp=1. If Majority is isp=2, use schi=-1d0.
c        allocate(jcoup(nw_i:nw,neps))
c        mmnorm=1d0
c        call RSMPI_Stop("chipm.and.eiqr: this cannot happen")
c      endif

c RS: chipm == F. Removed this part.
c      if(chipm.and.(.not.legas).and.(.not.eiqr)) then
c      endif

      if(epsmode.and.sergeyv.and.nolfco) then
        allocate( rcxqmean(nwhis,npm))
! RS:
        allocate( rcxqmean_mpi(nwhis,npm,nspin))
      endif



c RS: open WVI and WVR
c RS: WVI and WVR is divided for q-subgroups (gwsrc/RSMPI_qkgroup_mod.F)
c RS: this block is moved here by Sakuma
c RS: total number of q-points is needed at this point (iqxend-iqxini+1)
c RS: subgroup for double parallelization (q and k points) is created
c RS: see gwsrc/RSMPI_qkgroup_mod.F
c
      call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
      call RSMPI_qkgroup_Init(iqxend-iqxini+1,nqbz)
      call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
      
      if(debug) print *,' xxx1a2:'
      if    (ixc==101) then
         call RSMPI_Stop("ixc==101, this cannot happen")
c RS:removed
      elseif (normalm) then !ixc==1.or.sergeyv) then
         if(paralellx0) then
           call RSMPI_Stop("paralellx0==T, this cannot happen")
c           ifrcw  = iopen('WVR.'//xxt(iqxini,iqxend),0,-1,mrecl)
c           ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend),0,-1,mrecl)
         else
c           ifrcw  = iopen('WVR',0,-1,mrecl)
c           ifrcwi = iopen('WVI',0,-1,mrecl)
c
c RS: For each subgroup, only processes with rank==0 handle file I/O
c RS: "file_io_qkgroup" is defined in gwsrc/RSMPI_qkgroup_mod.F
            if (file_io_qkgroup) then
               ifrcw  = iopen('WVR.RSMPI'//qkgroup_id,0,-1,mrecl)
               ifrcwi = iopen('WVI.RSMPI'//qkgroup_id,0,-1,mrecl)
            endif
         endif
c      elseif(chipm) then
c        ifrcw  = iopen('WVR',0,-1,mrecl)
      endif

c RS: log file for each process
      ifile_rsmpi = iopen ('lx0.RSMPI'//myrank_id_rsmpi,1,3,0) 
c RS: print information in lx0.RSMPIxxx
      write(ifile_rsmpi,*) 'lx0.RSMPI'//myrank_id_rsmpi
      write(ifile_rsmpi,*) "ifile_rsmpi = ",ifile_rsmpi
      write(ifile_rsmpi,*) "file_io_qkgroup = ",file_io_qkgroup
      write(ifile_rsmpi,*) "rank_global iqkgroup nproc_local rank_local"
      write(ifile_rsmpi,*) myrank_id_rsmpi, "  ",qkgroup_id,
     & nproc_qkgroup, myrank_qkgroup
      write(ifile_rsmpi,*) "nq_local : iq_index"
      write(ifile_rsmpi,*) nq_local_qkgroup," : ", iq_index_qkgroup(:)
      write(ifile_rsmpi,*) "nk_local : ik_index"
      write(ifile_rsmpi,*) nk_local_qkgroup," : ", ik_index_qkgroup(:)
      write(ifile_rsmpi,*) 
      write(ifile_rsmpi,*) "--- q-loop begins here... ---"
      write(ifile_rsmpi,*) 

c RS: set pointer of vcoul and Mix0vec to the right position..
      do iq=1,iq_index_qkgroup(1)-1
         read(ifvcfpout) ngb
         allocate(vcoul(ngb,ngb))
         read(ifvcfpout) vcoul,qx ! coulomb matrix
         deallocate(vcoul)
      enddo
      if(epsmode) then
         do iq=1,iq_index_qkgroup(1)-1
            read(ifgb0vec,*) qgbin(1:3),ngbin,igc0,dnorm
            allocate(gbvec(ngb),zzr(ngb))
            do i=1,ngb
             read(ifgb0vec,"(4d24.15)") gbvec(i), zzr(i)
            enddo
            deallocate(gbvec,zzr)
         enddo
      endif
c RS:
      if (Is_IO_Root_RSMPI()) then
         write(6,*) "RS: do 1001 start----- see lx0.RSMPIxxxxx"
      endif
      call RSMPI_Print_WTime()

C --- Calculate x0(q,iw) and W ------------------------------------------
c============ loop over q point ========================================
c=======================================================================
c      if(iepsmode==202) iqxini = nqibz+1
c
c
c      do 1001 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
      do 1001 iq_local = 1,nq_local_qkgroup
c RS:
c RS: set correct index
         iq = iqxini -1 + iq_index_qkgroup(iq_local)
c         call cputid (0)
         call cputid (ifile_rsmpi)
         q = qibze(:,iq)
         if(iq<=nqibz.or.ngczero) then ! iqx = iq
            qq=q 
         else      !    iqx = 1            ! corresponding q=0 
            qq=0d0 
         endif

         call readqg('QGcou', qq, ginv,  quu,ngc,ngveccB)

c ---readin coulomb matrix
         if(chipm.and.(.not.legas).and.(.not.eiqr)) then
            allocate(vcoul(1,1)) !dummy
            ngb = nbloch + ngc  !ngb is readin from vcoul 25jan2006
            call RSMPI_Stop("chipm.and.(.not.legas).and.(.not.eiqr)")
         else
            read(ifvcfpout) ngb ! oct2005 !I added ngb to VCCFP
            allocate(vcoul(ngb,ngb))
            read(ifvcfpout) vcoul,qx ! coulomb matrix
            if(sum(abs(q-qx))>1d-10) then
               write(6,"('iq q=',i5,3f9.4,' qx(VCCFP)=',3f9.4)")iq,q,qx
c     stop 'hx0fp0: VCCFP is not compatible'
               call RSMPI_Stop('hx0fp0: VCCFP is not compatible')
            else
               write(ifile_rsmpi,"('iq q=',i5,3f9.4,' qx(VCCFP)=',3f9.4)")iq,q,qx
               write(ifile_rsmpi,*) "ngb = ",ngb
            endif
         endif
c         write(6,"('do 1001 start----- iq q=',i5,3f9.4,' qq=',3f9.4)")iq,q,qq
         write(ifile_rsmpi,"('do 1001 start----- iq q=',i5,3f9.4,' qq=',3f9.4)")iq,q,qq

c --- eps mode
         if(epsmode) then       !iepsmode/=0) then      ! only calculate iq>nqibz
            if(iq>nqibz) then
               if(wqt(iq-nqibz)==0d0) then
                  iqixc2 = iq- (nqibz+nq0ix)
               if(chipm) then
c RS: chipm==F
                  call RSMPI_Stop("chipm")
               else
                  filepsnolfc =
     &                 'EPS'//charnum4(iqixc2)//'.nolfc.dat'
                  if (file_io_qkgroup) then
                     ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
                     write(ifepsdatnolfc,"(a)")
     &                    ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
                  endif
                  
                  if(.not.nolfco) then !ixc==3.or.iepsmode==203) then
                     fileps = 'EPS'//charnum4(iqixc2)//'.dat'
                     if (file_io_qkgroup) then
                        ifepsdat = iopen ( fileps,1,3,0)
                        write(ifepsdat,"(a)") 
     &                       ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
                     endif
                  endif
               endif
               goto 110
               endif ! if(wqt==0d0)
            endif   ! if(iq>nqibz)
            deallocate(vcoul)
            cycle               !Cycle here
 110        continue

            if((.not.chipm).or.legas.or.eiqr) then
               read(ifgb0vec,*) qgbin(1:3),ngbin,igc0,dnorm
               write(ifile_rsmpi,"(' --- Readin Mix0vec: ',3d13.6,2i5,d18.8)")
     &      qgbin(1:3),ngbin,igc0,dnorm
            if(ngb/=ngbin) call RSMPI_Stop('hx0fp0: ngb/=ngbgin')
            if(sum(abs(qgbin(1:3)-q)) >1d-8) call RSMPI_Stop('hx0fp0: qgbin/=q')
            write(ifile_rsmpi,"(' ngb nwp niw=',3i8)")ngb,nwp,niw
            allocate(gbvec(ngb),zzr(ngb),x0mean(nw_i:nw))
            x0mean=0d0
            do i=1,ngb
             read(ifgb0vec,"(4d24.15)") gbvec(i), zzr(i)
            enddo
         elseif(chipm.and.nolfco) then
            allocate(zzr(ngb),x0mean(nw_i:nw))
            x0mean=0d0
            zzr   =0d0
            zzr(1:nbloch)=svec(1:nbloch)
            call RSMPI_Stop("chipm.and.nolfco")
         endif
      endif ! if (epsmode)
   
      allocate( zw0(ngb,ngb),
     &     zxq (ngb,ngb,nw_i:nw), 
     &     zxqi(ngb,ngb,niw))
      write(ifile_rsmpi,*) "---  end of allocation zw0 zxq zxqi"
      if(ixc==6 ) then
         allocate( rcxq(ngb,ngb,nwp,npm) )
         call RSMPI_Stop("ixc==6")
      elseif(epsmode.and.sergeyv.and.nolfco) then !iepsmode==202) then
         allocate( rcxq(1,1,1,npm) )
!RS:
         allocate( rcxq_mpi(1,1,1,npm,nspin) )
         rcxqmean=0d0
! RS:
         rcxqmean_mpi=0d0
      elseif(sergeyv) then
         allocate( rcxq(ngb,ngb,nwhis,npm) )
! RS:
         allocate( rcxq_mpi(ngb,ngb,nwhis,npm,nspin) )
         if (.not. allocated(rcxqmean_mpi)) then
           allocate(rcxqmean_mpi(1,1,nspin))
         endif
      endif
      zxq=0d0; zxqi=0d0
      if(sergeyv) then
        rcxq=(0d0,0d0)
        rcxq_mpi=(0d0,0d0)
      endif

c============ loop over spin============================================
c=======================================================================
      do 1003 is = 1,nspinmx  
c         write(6,"(' *** ',2i4,' out of nqibz+n0qi nsp=',2i4,' *** ')") 
c     &        iq, is, nqibz + nq0i,nspin
         write(ifile_rsmpi,"(' *** ',2i4,' out of nqibz+n0qi nsp=',2i4,' *** ')") 
     &        iq, is, nqibz + nq0i,nspin
c         if(debug) print *,' niw nw=',niw,nw

C.. spin chi_charge or chi_+-
         isf=is

c RS: removed
c         if(chipm) then
c         endif


c===========tetraini block tetra==.true.===============================1ini
c --- tetrahedron method --- preparatory part
c  nbnb    : total number of pairs (ib,jb) with non-zero weight.
c  ib=n1b  : band index for occ.   1\ge n1b \ge nband+nctot. 
c         Valence-core order(Core index follows valence index).
c  jb=n2b  : band index for unocc. 1\ge n2b \ge nband
c  wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
          if(tetra) then
           write(ifile_rsmpi,"(' tetra mode niw nw=',2i5)") niw,nw  
C ... kqxx(kx) ekxx --- 
c    ekxx1  for  rk
c    ekxx2  for q+rk.  See tetwt4
            do kx = 1, nqbz
              call readeval(    qbz(:,kx)-qqq, is, ekxx1(1:nband, kx) )
              call readeval(q + qbz(:,kx)-qqq,isf, ekxx2(1:nband, kx) )
            enddo
C ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
c     the pairs are not dependent on the energy omega 
c     in the denominator of the dielectric function.
            write(ifile_rsmpi,"(' -- First tetwt5 is to get size of array --')")
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
            if(ixc/=4.and.ixc/=6.and.(.not.sergeyv)) deallocate(demin,demax)
            nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
            if(debug) print *,' nbnbx=',nbnbx
            allocate(  n1b(nbnbx,nqbz,npm)
     &                ,n2b(nbnbx,nqbz,npm))
            do jpm=1,npm
              call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx, !2006
     o        n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
            enddo
c            if(sum(abs(nbnb-nbnbtt))/=0)stop'hx0fp0:sum(nbnb-nbnbtt)/=0'
            if(sum(abs(nbnb-nbnbtt))/=0) then
               call RSMPI_Stop('hx0fp0:sum(nbnb-nbnbtt)/=0')
            endif
            noccxv = maxval(noccxvx)
            noccx  = nctot + noccxv
c            print *,' Tetra mode: nctot noccxv= ',nctot,noccxv
            write(ifile_rsmpi,*) ' Tetra mode: nctot noccxv= ',nctot,noccxv
            deallocate(iwgt)
         endif !if(tetra)
c=========end of tetraini block==========================================1end


c=======TetrahedronWeight_4 block. tettwt4  ixc<=5 ==============2ini
C RS: removed
c=======TetrahedronWeight_4 block end.========================2end


c==============x0kf_v2 block original(.not.sergeyv) case =====================3ini
C RS: removed
c=============x0kf_v2 block end===============================3end
         if(epsmode.and.nolfco) then
            write(ifile_rsmpi,*)
     &           ' nolfco --- end of x0kf_v2 sum x0mean',sum(abs(x0mean))
         endif


c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
        if(ixc==4.or.sergeyv) then
C --- METHOD (tetwt5) for the tetrahedron weight
!   Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
! takao-feb/2002 
            if(frhis(1)/=0d0) call RSMPI_Stop(' hx0fp0: we assume frhis(1)=0d0')
            write(ifile_rsmpi,*)' --- Frequency  histogram  sections (a.u.)---- '
            if(onceww(1)) then
               do ihis= 1, nwhis
                  write(ifile_rsmpi,"(' ihis Init  End=', i4,2f13.6)")
     &                 ihis,frhis(ihis),frhis(ihis+1)  
               enddo
            endif  
            write(ifile_rsmpi,*)' ----------------nbnbx nqbz= ',nbnbx,nqbz
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
c            print *,' nhwtot=',nhwtot
            write(ifile_rsmpi,*)' nhwtot=',nhwtot
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
            write(ifile_rsmpi,*) 'goto tetwt5x_dtet4 job=',job
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
        endif     
c=======TetrahedronWeight_5 block end ==================4end




c==============x0kf_v3h block. ixc==6,11 =========================5ini
c Instead of wwk, you have to use whw(nwhtot). 
c You will need to pass these arrays to x0kf
c   ihw(ibjb,kx): omega index, to specify the section of the histogram.
c   nhw(ibjb,kx): the number of histogram sections
c   jhw(ibjb,kx): pointer to whw 
c   whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
c     : histogram weights for given ib,jb,kx for histogram sections
c       from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c------------------------------------------------------------------
        if(sergeyv) then ;  write(ifile_rsmpi,*) ' --- goto x0kf_v3hz '
c- new version only along real axis for 
          call x0kf_v3hz_MPI(npm,ncc,   
     i          ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf,
     i              qbas,ginv,  ef,qbz,wbz, 
     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, 
     i              ppb(1,is), 
     i              icore,ncore,
     d            nlmto,nqbz,noccx,noccxv,
     d            nl,nclass,natom,nnc,
     d            nlnmx,mdimx,nbloch,  nwhis,
     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
     i    nqbze,nband,nqibz, 
ctm, 080311
     i    qbze,iwbnd,nwf,lwbnd,
ctm, 080820
     i    ib_ent,lent,
c     o    rcxq, ! See rcxq_zxq in x0kf_v3h.f; rcxq is equivalent with zxq; 
     o    rcxq_mpi(:,:,:,:,is), ! See rcxq_zxq in x0kf_v3h.f; rcxq is equivalent with zxq; 
                !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,  !for nolfco
c     o    rcxqmean   ) !for nolfco
c RS:
     o    rcxqmean_mpi(:,:,is)   ) !for nolfco
          rcxq(:,:,:,:)=rcxq(:,:,:,:)+rcxq_mpi(:,:,:,:,is)
          if (epsmode.and.sergeyv.and.nolfco) then
            rcxqmean(:,:)=rcxqmean(:,:)+rcxqmean_mpi(:,:,is)
          endif
        endif ! if(sergeyv)

cccccccccccccccccccc
c          do iw  = 1,nwhis 
c          do i=1,ngb
c          do j=1,ngb
c             if(i==j) rcxq(j,i,iw)= dreal(rcxq(j,i,iw))
c          enddo; enddo; enddo
cccccccccccccccccccccc

c=========== Check write mode4 block ixc==4 ============6ini
c RS: removed
c         if(ixc==4) then
c         endif  ! end of if(tetra)
c=========== Check write mode4 block end  ==============6end

         if(ixc==4.or.sergeyv)         deallocate(ihw,nhw,jhw, whw,ibjb)
         if(tetra.and.(.not.sergeyv) ) deallocate(wwk) ! !ixc<=5.or.ixc==101.or.ixc==22.or.ixc==23)) 
         if(tetra)                     deallocate( n1b,n2b)

         iecut=1

         if(sergeyv) then !This was in 1003 loop as if(ixc==11.and.is==nspinmx) then
          if(debug) write(6,"(a)") ' --- goto dpsion4 --- '
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
cccccccccccccddddddccccccccccccccccccccccc
c          rcxq(:,:,:,1) =rcxq(:,:,:,2)
cccccccccccccccccccddddddccccccccccccccccc
               write(ifile_rsmpi,*)"goto dpsion4 nwhis nw_i nw_w=",nwhis,nw_i,nw
               call dpsion4(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i        ngb, rcxq, npm,nw_i,  ! rcxq is alterd---used as work
     o        zxq, zxqi, 
     i        nolfco,chipm, schi,is,  rcxqmean,ecut(iecut),ecuts(iecut), !iepsmode,
     o        x0mean)  !zxq x0mean is accumlating veriable for chipm mode
            endif ! else
          endif  !if(is==nspinmx)
          write(ifile_rsmpi,*)' --- end of dpsion4 ----'
          if(epsmode.and.nolfco) write(ifile_rsmpi,*) 
     &         ' end of dpsion4 sum x0mean',sum(abs(x0mean))
       endif !if(sergeyv)
 1003   continue  !end of spin loop =====
        if(allocated(rcxq) ) deallocate(rcxq)
! RS:
        if(allocated(rcxq_mpi) ) deallocate(rcxq_mpi)

 1033    continue
         write(ifile_rsmpi,*) " end of 1033 xxxxxxxxx"

c=======================================================================


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
            if(chipm) then      !ixc==22.or.ixc==23) then
c RS: removed
               call RSMPI_Stop("chipm==T")
            elseif(epsmode.and.nolfco.and.(.not.chipm)) then !ixc==2.or.iepsmode==202) then
               if (nspin==1) x0mean= 2d0*x0mean !if paramagnetic, multiply x0 by 2
               if (anfexist) then 
! Not tested yet
!          if(anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
                  write(6,*) 'I have not tested anfexist=t and ixc==2 '//
     &                 '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
                  call RSMPI_Stop('hx0fp0: tail output! Not yet for anfexist&ixc==2')
               endif            ! if(anfexist)
            else 
               if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
               if (anfexist) then !Antiferro accelaration.
c           nwmax = nw   
c           if(ixc==4.or.ixc==5.or.sergeyv) nwmax=nw
                  do iw  = nw_i,nw !max 
                     call anfx0k(natom,nclass,nblocha,iclass,pos
     i                    ,nbloch,ngc, !ngcni(iqx),
c     i    q, ngveccB(1,1,iqx),qbas, ! for q+G
     i                    q, ngveccB,qbas, ! for q+G
     i                    anfvec,ibasf, ! these are antiferro informations.
     i                    zxq(1,1,iw)) ! i/o 
                  enddo
               endif !anfexist
            endif  !if(chipm),else....

c         write (ifxd,"(1x,3f10.4)") q(1),q(2),q(3)
c         write (ifrx) rxq,cxq
            if(epsmode) then
               if(nolfco) then
                  ttt='without LFC'
               else
                  ttt='with LFC'
               endif
               if(chipm) then
                  write(6,*) '--- chi_{+-}}^{-1}      --- '//ttt
               else
                  write(6,*) '--- dielectric constant --- '//ttt
               endif
            elseif(ixc==1.or.sergeyv) then
               write(ifile_rsmpi, *)" trace check for W-V"
            endif

C... Get O^{1/2}= sqrt(ovlp)
            if(chipm.and.(.not.nolfco)) then
               call  RSMPI_Stop("chipm.and.(.not.nolfco)")
c RS:removed
            endif

            write(ifile_rsmpi,*)"zzzz nw_i nw=",nw_i,nw
            do iw  = nw_i,nw    !Feb2006. Before it was 1:nwp (nwp=nw+1).
                           !  So freq_r(iw-1) is shifted to freq_r(iw).
cnw all freq is from zero freq(0:nwp-1); nwp is the number of real flequencies.
               frr= dsign(freq_r(abs(iw)),dble(iw))
               if(.not.epsmode) then !if(ixc==1.or.sergeyv) then
                  imode = 1
                  if(debug) write(6,"(' rrr sum zxq=',i5,4d23.14)") iw,sum(zxq(:,:,iw)),sum(abs(zxq(:,:,iw)))
                  call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                 zw0)
                  write(ifile_rsmpi,"(' rrr sum zw0=',i5,4d23.14)") iw,sum(zw0),sum(abs(zw0))
                  zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)        
c     write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
                  if (file_io_qkgroup) then
c RS: index changed..
                     write(ifrcw, rec=((iq_local-1)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw !  WP = vsc-v
                  endif
             !----------------------------
c RS: arguments of subroutine tr_chkwrite() is changed
c
c              call tr_chkwrite("freq_r iq iw realomg trwv=", 
c     &               zw, iw, frr,nblochpmx, nbloch,ngb,iq)
c
                  call tr_chkwrite_RSMPI("freq_r iq iw realomg trwv=", 
     &                 zw, iw, frr,nblochpmx, nbloch,ngb,iq,ifile_rsmpi)
             !----------------------------
               elseif(epsmode.and.(.not.chipm)) then !.and.ixc/=23) then ! No LFC (local field correction). It's better to use echo 4| hbasfp0. 
                  vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) ) 
                  if(.not.nolfco) then
                     x0mean(iw) = sum( dconjg(zzr)* matmul(zxq(:,:,iw),zzr)) 
                  endif
                  epsi(iw,iqixc2) = 1d0/(1- vcmean * x0mean(iw)) 
                  write(ifile_rsmpi,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2d23.15,3x, 2d23.15,
     &                 " vcmean x0mean =", 2d23.15,3x, 2d23.15)') iqixc2,iw,2*frr,
     &                 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean,x0mean(iw)
                  if (file_io_qkgroup) then
                     write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &                    q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                  endif

                  if(.not.nolfco) then ! With LFC
                     imode=2
                     call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                    zw0)  !  write(6,"('ssschk1=',3d13.5)") sum(abs(zw0)) sum(abs(gbvec))
                     epsi(iw,iqixc2)= sum( dconjg(gbvec) * matmul(zw0,zzr) ) 
                     write(ifile_rsmpi,'( " iq iw omega eps epsi  wLFC="
     &                    ,2i6,f8.3,2d23.15,3x, 2d23.15)') 
     &                    iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                     write(ifile_rsmpi,*)
                     if (file_io_qkgroup) then
                        write(ifepsdat,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &                       q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                     endif
                  endif

c RS: removed
c            elseif(epsmode.and.chipm) then
c              if(nolfco) then !ixc==22) then 
c              else
c              endif
c            endif
               endif
            enddo               ! do iw=..

c         if(chipm.and.(.not.nolfco)) 
c     &          deallocate(zzz,UU,VT,SS, sqovlp,sqovlpi,UU0,VT0,SS0)
            if( allocated(zzr)   ) deallocate(zzr)
            if( allocated(x0mean)) deallocate(x0mean)
            if( allocated(gbvec) ) deallocate(gbvec)

c RS:removed
c         if( ixc==5.or.ixc==6 ) then
c         endif
         endif                  !if(realomega)

C... Close files for epsmode
         if(epsmode) then       !iepsmode/=0) then      ! only calculate iq>nqibz
            if(chipm) then
               ifchipmn = iclose ( 'ChiPM'//charnum4(iqixc2)//'.nolfc.dat')
               if(.not.nolfco) then
                  ifchipm = iclose ('ChiPM'//charnum4(iqixc2)//'.dat')
               endif 
            else
               if (file_io_qkgroup) then
                  filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
                  ifepsdatnolfc = iclose( filepsnolfc)
               endif
               if(.not.nolfco .and. file_io_qkgroup) then
                  fileps = 'EPS'//charnum4(iqixc2)//'.dat'
                  ifepsdat = iclose( fileps)
               endif
            endif
         endif
C--- RealOmega end ==============================


     
C--- ImagOmega ===================================
         if (imagomega) then
            write(ifile_rsmpi,*) ' goto imag omega'
            if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
c ... Antiferro accelaration.
            if (anfexist) then
               do iw  = 1,niw 
                  call anfx0k(natom,nclass,nblocha,iclass,pos
     i                 ,nbloch,ngc, !ngcni(iqx),
c     i      q, ngveccB(1,1,iqx),qbas,    ! for q+G
     i                 q, ngveccB,qbas, ! for q+G
     i                 anfvec,ibasf, ! these are antiferro informations.
     i                 zxqi(1,1,iw)) ! i/o 
               enddo
            endif
            if (ecorr_on>0)then !ixc==101.or.(sergeyv.and.imagonly)) then 
               imode=0
            else
               imode=1
            endif

            do iw  = 1,niw
               write(ifile_rsmpi,*) 'sss sumchk zxqi=',sum(abs(zxqi(:,:,iw)))
               call wcf( ngb, vcoul,zxqi(1,1,iw),imode, 
     &              zw0)
               write(ifile_rsmpi,*) 'sss sumchk zw0=',sum(abs(zw0))
               zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)    
               if (file_io_qkgroup) then
                  write(ifrcwi, rec=((iq_local-1)*niw+iw) ) zw !  WP = vsc-v 
               endif
C--- for Miyake's total energy !Nov2004 ----------------
c RS: removed
c               if(ecorr_on>0) then
c               endif

               call tr_chkwrite_RSMPI("freq_i iq iw imgomg  trwv=",
     &              zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq,ifile_rsmpi)
            enddo
         endif !if (imagomega)
C... ImagOmega end ===============================
c$$$ 2005 continue ! iecut loop


c 1002 continue  ! end of frequency block-loop
      deallocate(vcoul, zw0, zxq, zxqi)
 1001 continue  
C--- end of loop over q point =================================
c=======================================================================
 1100 continue


c RS: removed
c      if(epsmode.and.legas) then
c      endif

C--- Write TEECOR ecorr_on mode
c RS: removed
c      if(imagomega.and.ecorr_on>0) then
c      endif
C...close files
c        close(ifdpin)

c RS: file "hbe.d" is closed earlier
c      is = iclose('hbe.d')
      call cputid(ifile_rsmpi)


      if(ixc==202) then 
         write(ifile_rsmpi,*) ' OK! hx0fp0_RSMPI
     &        mode=202 sergeyv epsPP NoLFC'
         if (Is_IO_Root_RSMPI()) then
            write(0,*) ' OK! hx0fp0_RSMPI
     &           mode=202 sergeyv epsPP NoLFC'
            write(6,*) ' OK! hx0fp0_RSMPI
     &           mode=202 sergeyv epsPP NoLFC'
         endif
      else if(ixc==203) then
         write(ifile_rsmpi,*) ' OK! hx0fp0_RSMPI
     &        mode=203 sergeyv eps LFC '
         if (Is_IO_Root_RSMPI()) then
            write(0,*) ' OK! hx0fp0_RSMPI
     &           mode=203 sergeyv eps LFC '
            write(6,*) ' OK! hx0fp0_RSMPI
     &           mode=203 sergeyv eps LFC '
         endif
      else if(ixc==11) then
         write(ifile_rsmpi,*) ' OK! hx0fp0_RSMPI
     &        mode=11 normal sergeyv'
         if (Is_IO_Root_RSMPI()) then
            write(0,*) ' OK! hx0fp0_RSMPI
     &           mode=11 normal sergeyv'
            write(6,*) ' OK! hx0fp0_RSMPI
     &           mode=11 normal sergeyv'
         endif
      else
         call RSMPI_Stop("This cannot happen")
      endif
      
      call RSMPI_Finalize()

      write(ifile_rsmpi,*) "MPI_Finalize is called."
      end program hx0fp0_MPI


c----------------------------------------------
      subroutine tr_chkwrite_RSMPI(tagname,zw,iw,freqq,
     &     nblochpmx,nbloch,ngb,iq,ifile) !RS: ifile is added
      implicit none
      integer(4):: nblochpmx,nbloch,ngb,iw,i,iq
      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
      real(8):: freqq
      logical :: smbasis
      character*(*)::tagname
c RS:
      integer(4) :: ifile
      trwv=0d0
      if(.not.smbasis()) then
         do i = 1,nbloch
            trwv = trwv + zw(i,i) 
         enddo
      endif
      trwv2 = 0d0
      do i = 1,ngb
         trwv2 = trwv2 + zw(i,i) 
      enddo                     !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(ifile,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
      end
c----------------------------------------------

