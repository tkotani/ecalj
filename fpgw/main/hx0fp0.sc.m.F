      program hx0fp0_sc
!!  Calculate W-V for QSGW mode.
!! We calculate chi0 by the follwoing three steps.
!!  tetwt5: tetrahedron weights
!!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
!!  dpsion5: calculate real part by the Hilbert transformation from the Im part
      use m_readqg,only: readngmx,readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval
      use m_read_bzdata,only: read_bzdata,
     &     ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,qbasmc,
     &     dq_bzcase2,qbz,wbz,qibz,wibz,qbzw,
     &     idtetf,ib1bz,idteti,
     &     nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     &     nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     &     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw,
     &     alat,ef, diw,dw,delta,deltaw,esmr,symgrp,clabl,iclass,nindxv,nindxc,ncwf,
     &     invg, il, in, im, ilnm, nlnm, ilv,inv,imv,  ilnmv, nlnmv,
     &     ilc,inc,imc,  ilnmc, nlnmc,nindx,konf,icore,ncore,
     &     occv,unoccv ,occc,unoccc, nocc, nunocc, plat, pos,z,ecore, freq, symgg ! symgg=w(igrp)
      use keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
      use m_readqgcou,only: readqgcou
      use m_mpi,only: MPI__hx0fp0_rankdivider2Q,MPI__hx0fp0_rankdivider2S,
     &   mpi__Qtask,MPI__InitializeQSPBM,MPI__Finalize,mpi__root,
     &   MPI__Broadcast,MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,mpi__rank,mpi__size,
     &   mpi__Qranktab,MPI__consoleout,mpi__Ss,mpi__Se, mpi__allreducesumS,
     &   mpi__barrier, mpi__rankQ,mpi__rootQ,mpi__rootS


!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: rdpp,    !"call rdpp" generate following data.
     & nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr
!! Generate matrix element for "call get_zmelt".
      use m_zmel,only:    !these data set are stored in this module, and used when 
     & nband,itq,imdim,ngcmx,ngpmx,    ppovlz,
     & ppbir,shtvg, miat,tiat , ntq

      implicit none
      integer::iopen,maxocc2,iclose, ixc,iqxini,iqxend,
     &     nw_input, ifhbe,  nprecb,mrecb,mrece,nlmtot,nqbzt,!nband,
     &     nq0i,i,nq0ix,neps,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx !ngcmx,
     &     ,nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, !ifvcfpout,,mdimx,nbloch
     &     ifepscond,nxx,ifvxcpout,ifgb0vec
     &     ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &     ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &     ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &     ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx(2)=-9999,nwmax !,ifev1,ifev2
     &     ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &     ,  ifchipmlog ,   nw_w,nwmin  ! ,ngpmx
      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm, dwry,dwh,omg_c,omg2, q(3),  qgbin(3),qx(3) 
      real(8):: ua=1d0          ! this is a dummy.
      integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2), ndble=8, nword
      real(8),allocatable:: vxcfp(:,:),  wqt(:), wgt0(:,:),q0i(:,:) !,nx(:,:),nblocha(:)
      integer,allocatable :: ngveccB(:,:), ngvecc(:,:),iqib(:),ifppb(:) !,lx(:)
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &     zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: eqt(:), !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:)
     &     ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
     &     ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), qbze(:,:),qibze(:,:) !,ecore(:,:)
c     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: rcxq(:,:,:,:)
      logical :: tetra,tmpwwk=.false. ! If tmpwwk=.true., this use a temporary file tmp.wwk so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer,allocatable :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),noccxvv(:) !  &         idtetf(:,:),ib1bz(:)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:) !qbzw(:,: )
      real(8) ::qbzx(3),anfvec(3)
      logical :: anfexist,  debug
      integer,allocatable:: ibasf(:)
      real(8),allocatable :: transaf(:,:), freq_r(:),freq_i(:),freqx(:),wx(:),expa(:),frhis(:) !ppovl(:,:)
      logical :: realomega, imagomega 
      complex(8),allocatable:: zzr(:,:),ppovl(:,:),ppovlzinv(:,:) !,ppovlz(:,:)
      complex(8) :: epxxx,vcmean
      character*9 fileps
      character*15 filepsnolfc
      logical :: paralellx0=.true. !, hist
      character(5) :: charnum5
      character(20):: xxt
      real(8) :: Emin, Emax      ,emax2,emin2
c     integer :: iSigma_en  !sf..21May02  !iSigma_en is integer
                                !parameter stored in GWIN_V2
                                !which determines approximation for  self-energy.
                                !Self-energy should be made hermitian for energies to be real
cxxx  !iSigma_en==0 SE_nn'(ef)+img integral:delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
cxxx  !iSigma_en==1 SE_nn'(ef)+delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
                                !iSigma_en==2 [SE_nn'((e_n+e_n')/2)+h.c.]/2
                                !iSigma_en==3 [(SE_nn'(e_n)+SE_nn'(e_n'))/2+h.c.]/2
      real(8) :: omg2max,omg1max
      real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)
      integer,allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      logical::imagonly=.false. , noq0p !,readgwinput
      integer::nwin, incwfin, verbose,nbcut,nbcut2,ifpomat,nnmx,ikpo,nn_,noo,iqxxx,nomx
      real(8)::efin
      logical :: nolfco=.false.
      complex(8) :: x0meann     !dummy
      integer:: isp1,isp2, ngc,mrecg ! bzcase,
      real(8)::  quu(3),deltaq(3),qqq(3)=0d0 !
      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)
      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg,smbasis
      character(2):: soflag
      real(8):: q_r(3)
      complex(8),allocatable:: pomat(:,:)
      logical   :: timereversal,onceww
      integer:: npm ,jpm,ncc,nw_i
      real(8)   :: frr !, sciss
c      logical :: newaniso,newaniso2 !,z1offd
      integer :: ngb0,ifvcoud,idummy,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x !ifepstinv,
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:),eemat(:,:),zcousq0(:,:)
      real(8),allocatable:: vcousq(:),vcousq0(:),vcoudummy(:)
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq
!! Eq.(40) in PRB81 125102
c     complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:), 
     &     w_k(:),w_ks(:),w_kI(:), w_ksI(:) 
      complex(8),allocatable:: llw(:,:), llwI(:,:),w0(:),w0i(:),aaamat(:,:)
      real(8),allocatable:: dmlx(:,:),epinvq0i(:,:),epinv(:,:,:),epinvq0i_m1(:,:),wklm(:),qeibz(:,:,:)
      integer:: lxklm,nlxklm,ifidmlx,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
      complex(8):: vc1vc2
      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
      real(8),allocatable:: aik(:,:,:,:)
      integer,allocatable:: aiktimer(:,:)
      integer:: l2nl, nmbas_in , iqxendx,imb2 !iqqv,
      logical:: eibz4x0,tiii,iprintx,chipm=.false.,iqinit,localfieldcorrectionllw
      real(8)::qvv(3),ecut,ecuts,hartree,rydberg,pi
      character(128):: vcoudfile
      integer :: iqeibz
      complex(8):: epslfc, axxx(10)
      integer:: src,dest
      integer:: ifw0w0i
      logical :: readw0w0i,symmetrize,eibzmode
      real(8):: schi=-9999 !dummy

      integer:: i_reduction_npm, i_reduction_nwhis,  i_reduction_nmbas2
      logical:: crpa

      integer,allocatable :: iclasst(:), invgx(:)
      integer:: ibasx,ificlass,ifile_handle
      complex(8),allocatable:: ppovl_(:,:)
!-------------------------------------------------------------------------
!TIME0 
!TIME0 
      call MPI__InitializeQSPBM()
      call MPI__consoleout('hx0fp0_sc')
      call cputid (0)
      allocate( zzr(1,1)) !dummy
      hartree= 2d0*rydberg()
      pi     = 4d0*datan(1d0)
      fourpi = 4d0*pi
      sqfourpi=sqrt(fourpi)
      write(6,*) ' --- hx0fp0_sc Choose omodes below ----------------'
      write(6,*) '  ixc= 11,10011,or 1011 '
      write(6,*) ' --- Put number above ! -----------------'
      if( mpi__root ) then
         read(5,*) ixc !c     call readin5(ixc,iqxini,iqxend)
      end if
      call MPI__Broadcast(ixc)
      crpa=.false.
      if(ixc==0) call rx( ' --- ixc=0 --- Choose computational mode!')
      call headver('hx0fp0_sc',ixc)
c      call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
c      write(6,"(' ScaledGapX0=',f8.3)") sciss
      if(ixc==11) then
         write(6,*) " OK ixc=11 normal mode "
      elseif(ixc==10011) then
         write(6,*) " OK ixc=10011 crpa mode "
         crpa=.true.
      elseif(ixc==1011) then
         write(6,*) 'OK ixc=1011 Add W0W0I part at q=0'
      else
         write(6,*)'we only allow ixc==11. given ixc=',ixc
         call rx( 'error:we only allow ixc==11.')
      endif
!! newaniso2 is now fixed to be .true.
      call getkeyvalue("GWinput","ecut_p" ,ecut, default=1d10 )
      call getkeyvalue("GWinput","ecuts_p",ecuts,default=1d10 )
!! naraga says this cause a stop in ifort --->why???
c    write(6,*)'Timereversal=',Timereversal()

!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
      call read_BZDATA()
!TIME1 "after readbzdata"
!TIME0 
!! Use regular mesh even for bzcase==2 and qbzreg()=T
!!     off-regular mesh for bzcase==1 and qbzreg()=F
c      if( ( bzcase()==2.and.qbzreg() )       .or.
c     &     ( bzcase()==1.and.(.not.qbzreg()))      ) then
      if(.not.qbzreg()) then ! set off-gamma mesh
         deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
         do i=1,nqbz
            qbz(:,i) = qbz(:,i) - deltaq/2d0
            write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
         enddo
      endif
      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
!! === Readin by genallcf ===
!! See "use m_genallcf_v3" at the begining of this routine
!! We set basic data.
      nwin   = 0                !Readin nw from NW file
      incwfin= 0                !use ForX0 for core in GWIN
      efin =  0d0               !readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      nw_input = nw ; write(6,*) 'nw delta=',nw_input,delta
      debug=.false.
      if(verbose()>=100) debug=.true.
      if(debug) write(6,*)' end of genallc'
      tpioa=2d0*pi/alat
!!!!  WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call rx( ' hx0fp0_sc: nclass /= natom ')
!! --- tetra or not
      if(delta <= 0d0) then
         tetra =  .true.
         delta = -delta
         write(6,*)' hx0fp0: tetrahedron mode delta=',delta
      else
         tetra = .false. ! switch for tetrahedron method for dielectric functions
      endif
!! --- read dimensions of h,hb
      ifhbe      = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if(nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
      if(nqbz /=nqbzt ) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')
!! --- Readin Offset Gamma --------
      if(debug) write(6,*) 'reading QOP'
      ifiq0p=ifile_handle()
      open (ifiq0p,file='Q0P')
      read (ifiq0p,"(i5)") nq0i
      write(6,*) ' ### nqibz nq0i=', nqibz,nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      do i=1,nq0i
         read (ifiq0p, * ) wqt(i),q0i(1:3,i)
      enddo
      nq0ix = nq0i
      do i=1,nq0i
         if(wqt(i)==0d0 ) then
            nq0ix = i-1
            exit
         endif
      enddo
      neps=nq0i-nq0ix ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      write(6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(ifiq0p)
c$$$      if(.not.newaniso2) then
c$$$         wqtsum = sum(abs(wqt(1:nq0i)))
c$$$         call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
c$$$      endif 

!TIME1 "before readngmx"
!TIME0
      call getsrdpp2(nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx
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
!! ---------------- dummy ngrpx=1 -------------------
      ngrpx = 1
      l2nl=2*(nl-1)
      allocate(symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
!!  dummy. Get space-group transformation information. See header of mptaouof.
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
      call mptauof(symope,ngrpx,plat,natom,pos,iclasst
     o ,miat,tiat,invgx,shtvg ) !note: miat,tiat,shtvg are defined in m_zmel.
      if(verbose()>=40) write (*,*)' hsfp0.sc.m.F: end of mptauof'
!!  call rdpp gives ppbrd = radial integrals and cgr = rotated cg coeffecients.
      call rdpp(nxx, nl, ngrpx, nn, nclass, nspin, symope,qbas)
      ntq=nband
      allocate(itq(ntq))
      do i=1,ntq
        itq(i)=i
      enddo
!! Pointer to optimal product basis
      allocate(imdim(natom))
      call indxmdm (nblocha,nclass,iclass,natom,
     o imdim )                  !in m_zmel
      nblochpmx = nbloch + ngcmx
      allocate(ngveccB(3,ngcmx))
      iqxend = nqibz + nq0i
      write(6,*) ' nqibz nqibze=',nqibz,nqibze
!TIME1 "before initreadeigen"
!TIME0
!!... initialization of readEigen !readin m_hamindex
      call init_readeigen(ginv,nspin,nband,mrece)!EVU EVD are read in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)
c     --- ecore ---
c     allocate(ecore(nctot,nspin)) !core energies
c     do  is = 1,nspin
c     if (nctot .gt. 0) then
c     call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c     write(6,*)' ecore is=',is,ecore(:,is)
c     endif
c     enddo

c     --- set realomega, imagomega tetra nw niw nwp ifgb0vec ---------------------
!     nwp, freq_r,  frhis(1:nwhis+1)
c     if   ( ixc==1 ) then !old imagw = 2 case
c     realomega =.true.
c     imagomega =.true.
c     stop 'hsfp0sc: ixc==1 is not implimented'
ccccccccccccccccccccfaleev 21May02,  use only ixc=1,11 modes cccccccccc
c     elseif( ixc==2.or.ixc==3 ) then
c     realomega =.true.
c     imagomega =.false.
c     niw = 0
c     ifepscond = 2102
c     open (ifepscond,file='EPScond')
c     read (ifepscond,*) epsrng, dwry !epsrng dw in Ry
c     dw = dwry/2d0
c     close(ifepscond)
c     if(dw==0d0) then
c     nw = 1
c     else
c     nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
c     endif
c     allocate(epsi(nw,neps))
c     if(paralellx0) then
c     ifgb0vec = iopen ( "Mix0vec."
c     &     //xxt(iqxini,iqxend),1,3,0)
c     c          "//charnum5(iqxini)//'to'//charnum5(iqxend),1,3,0)
c     else
c     ifgb0vec = iopen ( "Mix0vec",1,3,0)
c     endif
c     elseif(ixc==4.or.ixc==5.or.ixc==6) then
c     ! ... These are test modes.
c     ! ixc=4 tetrahedren weight test. tetwt5.vs.tetwt5. Write tethis.chk
c     ! ixc=5 Spectrum function (Img part) along the Real axis with tetwt4
c     ! ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.
c     realomega = .true.
c     imagomega = .false.
c     tetra     = .true.
c     niw = 0
c     ! ---  For tetwt5 ---  the tetrahedron weight for spectrum function (imaginary part)
c     !   Histogram bins are specified by freq_r(1:nwp)
c     !     nwp=nw+1; frhis(1)=0
c     !     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
c     !     The last  bin  is     [frhis(nw), frhis(nwp)].
c     
c     ! ... These parameters specifies a test histgram bins;Sergey's mesh just for test modes.
c     nw0 = 200    !100    800
c     dwh = 0.01d0 !0.02d0 0.0025d0 !in hartree
c     ! ...
c     call findemaxmin(ifev,nband,nqbz,nspin,emax,emin)
c     if (nctot .gt. 0) Emin = minval(ecore)
c     omg2max = (Emax-Emin)*.5+.2d0    !(in Hartree) covers  all relevant omega, +.2 for margin
c     omg1max = dwh*(nw0-1)
c     nwp = int(sqrt(omg2max*(2*nw0-1d0)/dwh-(nw0**2-3*nw0+1d0)))+1 ! + 1 for margin
c     nw  = nwp-1
c     write(6,*) Emax,Emin,nw0,nw  ! nwp is new max number in frequency array
c     write(6,'(a32,2i7,2d15.3)')'hx0fp1: nw0,nw,omg1max,omg2max='
c     &              , nw0,nw,  omg1max,omg2max
c     if (nw <= nw0) stop 'hx0fp0:ixc==[456] nw2 <= nw'
c     allocate(freq_r(nwp))
c     do iw=1,nwp  !This is a test mesh by Sergey.Faleev
c     if(iw<=nw0) then;  freq_r(iw)=dwh*(iw-1)
c     else;  freq_r(iw)=dwh*(iw**2+nw0**2-3*nw0+1)/(2*nw0-1d0)
c     endif
c     enddo !freq_r(iw) is linear for iw<=nw and quadratic for nw<iw<=nw2
c     !freq_r(iw) chosen in such a way that it is continues with
c!!!  nw nwp=nw+1 freq_r(1:nwp) are used after here.
c     allocate(frhis(nwp))
c     frhis=freq_r(1:nwp)
c     nwhis=nw

!! We first accumulate Imaginary parts. Then do K-K transformation to get real part.
      realomega = .true.
      imagomega = .true.
      tetra     = .true.
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      call freq

C     --- tetra init
      if(tetra) then
         allocate(              !wgt(nband+nctot,nband,nqbz), !noccxvv(nw+niw),
     &        nbnb(nqbz,npm),nbnbtt(nqbz,npm)
     &        ,ekxx1(nband,nqbz),ekxx2(nband,nqbz)) !, !!! nband=nlmto
         if(tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
      endif
      if(debug) write(6,*)' xxx1:'
c     ... maximum no. occupied valence states
      noccxv     = maxocc2 (nspin,ef, nband, qbze,nqbze)
      if(noccxv>nband) call rx( 'hx0fp0: all the bands filled! too large Ef')
      noccx      = noccxv + nctot
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )
      nprecx = ndble            !We use double precision arrays only.

C...  This is just to get nblochpmx
c      if(smbasis()) then
c         call getngbpomat(nqibz+nq0i,
c     o        nnmx,nomx)
c         nblochpmx = nnmx
c      endif
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()

!TIME1 "before open WVd"
!TIME0
C     --- open WV* files
      if (mpi__root)then
      ifwd   = iopen('WV.d',1,-1,0)
      write (ifwd,"(1x,10i14)") !"(1x,i3,i8,i5,5i4)") 
     &     nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
      ifwd = iclose('WV.d'); ifwd=0
      endif
      allocate(  zw(nblochpmx,nblochpmx) )
      nspinmx = nspin
!TIME1 "after write ifwd"
!TIME0

!! --- READIN ANFcond. we don't use this path so much -----------
      inquire(file='ANFcond',EXIST=anfexist)
      write(6,*)' ANFcond exist? =',anfexist
      if(anfexist) then
         allocate(ibasf(natom),transaf(1:3,natom))
         ifianf = 211
         open(ifianf,file='ANFcond')
         read(ifianf,*) soflag
         if(soflag=='SO') then  !SO case is not implimented
            write(6,*)' ANFcond found. but SO case is not implimented. '
     &           ,'Not using ANF symm'
            anfexist=.false.
            goto 889
         endif
         read(ifianf,*)
         read(ifianf,*) anfvec(1:3)
         write(6,'(" ANFcond=",3f13.6)') anfvec
         do ibas=1,natom
            read(ifianf,*) ibas1, ibasf(ibas)
            if(ibas1/=ibas) call rx( ' readin ANFcond error:ibas1/=ibas')
            write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
c     &     ,transaf(1:3,ibas)
         enddo
         close(ifianf)
         nspinmx =1
      endif
 889  continue
c      if(smbasis() .and. anfexist) then
c         write(6,*)' smbais=T & anfexist=T is not yet. '//
c     &        ' At leaset anfx0k shold be changed'
c         call rx( ' smbais=T & anfexist=T is not yet. ')
c      endif
      if((.not.timereversal()) .and. anfexist) then
        write(6,*)' timereversal=F & anfexist=T is not yet examined.'//
     &       'At leaset you have to consider about anfx0k and so.'
        call rx( ' timereversal=F & anfexist=T is not yet.')
      endif
!!... these are used x0k
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2

c$$$!! --- ppb ---
c$$$!!   <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin
c$$$      do  is = 1,nspin
c$$$         if(debug) write(6,*)" goto ppbafp is=",is
c$$$         irot = 1
c$$$         call ppbafp_v2 (irot,ngrpx,is,nspin,
c$$$     i        il,in,im, nlnm,
c$$$     d        nl,nn,nclass,nlnmx,
c$$$     i        mdimx,lx,nx,nxx,  !Bloch wave    
c$$$     i        cgr, nl-1,        !rotated CG
c$$$     i        ppbrd,            !radial integrals
c$$$     o        ppb(1,is))        !,pdb(1,is),dpb(1,is),ddb(1,is) )
c$$$      enddo
!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
      irot=1
      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,irot,nspin))
      do is = 1,nspin
        call ppbafp_v2 (irot,ngrpx,is,nspin,
     i   il,in,im,nlnm,         !w(i_mnl),
     i   nl,nn,nclass,nlnmx,
     i   mdimx,lx,nx,nxx,       !Bloch wave    
     i   cgr, nl-1,             !rotated CG
     i   ppbrd,                 !radial integrals
     o   ppbir(:,irot,is))      !this is in m_zmel
      enddo  
      if(debug) write(6,*) ' end of ppbafp_v2'
!TIME1 "endo of ppbafp_v2"
!TIME0
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
      iqxini=1 !for newaniso
      eibzmode = eibz4x0()

!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      if(eibzmode) then
!! commentout block inversion Use iqxendx=iqxend because of full inversion
         iqxendx=iqxend
         allocate( nwgt(nqbz,iqxini:iqxendx), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     &        igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     &        eibzsym(ngrp,-1:1,iqxini:iqxendx))
         iprintx=.false.

         write(6,*)
         write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal() 
         if(mpi__root) iprintx=.true.
         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,
     i        timereversal(),ginv,iprintx,
     o        nwgt,igx,igxt,eibzsym,tiii)
         write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
         call cputid(0)

c$$$
c$$$         write(6,"('TimeRevesal switch = ',l1)") timereversal()
c$$$         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,timereversal(),ginv,iprintx,
c$$$     o        nwgt,igx,igxt,eibzsym)
c$$$!! Check timereversal is required for symmetrization operation or not. If not tiii=timereversal=F is used.
c$$$!! this is because the symmetrization is a little time-consuming.
c$$$         tiii=timereversal() 
c$$$         if(minval(igxt)==1) tiii=.false.
c$$$         iprintx=.true.
c$$$cccccccccccccccccccc
c$$$c     tiii=.true.
c$$$cccccccccccccccccccc
c$$$         write(6,"('=== goto eibzgen === used timereversal=',l1)")tiii
c$$$         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,tiii,ginv,iprintx,
c$$$     o        nwgt,igx,igxt,eibzsym)

         call PBindex(natom,lx,l2nl,nx) !all input. this returns requied index stored in arrays in m_pbindex.
                                ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
         call readqgcou()       !no input. Read QGcou and store date into variables.
c     call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.

      else !dummy allocation to overlaid -check bound !sep2014
         iqxendx=iqxend
         allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &    ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      endif

      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
      llw=1d99
      llwI=1d99

      if(ixc==1011)then !ixc==11 is a debug mode to test contrib. at \Gamma point.
         goto 1191
      endif

!! rank divider
      call MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
      call MPI__hx0fp0_rankdivider2S(nspinmx)

!! == Calculate x0(q,iw) and W == main loop 1001 for iq. 
!! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
!! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
!! (or need to modify do 1001 loop).
!! ---------------------------------------------------------------
!! === do 1001 loop over iq ============================================
!! ---------------------------------------------------------------
      iqinit=.true.
      write(6,'("irank=",i5," allocated(mpi__qtask)=",L5)')mpi__rank,allocated(mpi__qtask)
      do iq = iqxini,iqxend
         if(mpi__qtask(iq)) write(6,'("irank iq=",i5,i5)') mpi__rank,iq
      enddo
!TIME1 "before 1001 loop"
      do 1001 iq = iqxini,iqxend
         if( .not. mpi__Qtask(iq) ) cycle
!TIME0 
         if (mpi__rootS) then
         ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
         ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
         endif
         call cputid (0)
         q = qibze(:,iq)
         call readqg('QGcou', q, ginv,  quu,ngc,ngveccB) ! q was qq

!! Caution : confusing point
!!  ngc by QGcou is shown at the bottom of lqg4gw.
!!  ngc read from PPOVL are given by rdata4gw.
!!  Note that  ngc(iq>nqibz )=ngc (q=0), because when it is generated in mkqg.F
!!
c         if( newaniso2.and.iq==1 ) then ! *sanity check
         if( iq==1 ) then ! *sanity check
            if(sum(q**2)>1d-10) call rx( ' hx0fp0.sc: sanity check. |q(iqx)| /= 0')
         endif

!! ==== readin Coulomb matrix ====
         ngb = nbloch + ngc
         write(6,"('do 1001: iq q=',i5,3f9.4)")iq,q
         write(6,*)'nbloch ngb ngc=',nbloch,ngb,ngc

!! === readin diagonalized Coulomb interaction ===
!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.
c         if(newaniso2) then
         vcoudfile='Vcoud.'//charnum5(iq) ! iq was iqqv this is closed at the end of do 1001
         ifvcoud = iopen(trim(vcoudfile),0,-1,0)
         read(ifvcoud) ngb0
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
c         endif

c         if(newaniso2.and. iq>nqibz.and.(.not.localfieldcorrectionllw())  ) then 
          if(iq>nqibz.and.(.not.localfieldcorrectionllw())  ) then 
            if( ngb0/=ngb ) then
               call rx( 'hx0fp0.m.f:ngb0/=ngb')
            endif
            nolfco =.true.
            nmbas_in = 1 
c         elseif(newaniso2) then !.and.iq==1) then
          else
            nolfco = .false.
            nmbas_in = ngb
         endif
         nmbas1 = nmbas_in
         nmbas2 = nmbas1

!! newaniso=T case. Used in get_zmelt in m_zmel called in x0kf_v4hz
         if(allocated(ppovlz)) deallocate(ppovlz)
         if(allocated(ppovlzinv)) deallocate(ppovlzinv)
         if(allocated(ppovl)) deallocate(ppovl)
         allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb),   ppovlzinv(ngb,ngb))
         call readppovl0(q,ngc,ppovl) !q was qq
c         ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
c         ppovlz(nbloch+1:nbloch+ngc,:) 
c     &     = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
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
         deallocate(ppovl_)

         
c$$$  if(ixc==11) then
c$$$  write(6,*)" xxx2: memsize 8*ngb*ngb*nwhis=", 8*ngb*ngb*nwhis,' ngb nwhis=',ngb,nwhis
c$$$  allocate( rcxq(ngb,ngb,nwhis,npm) )
c$$$  rcxq=(0d0,0d0)
c$$$  else
c$$$  if(onceww(2)) write(6,*)" xxx2:allocate zxq zxqi memsize 16*ngb*ngb*(nwp+niw)=",
c$$$  &    16*ngb*ngb*(1+nwp+niw),' ngb nwp niw=',ngb,nwp,niw
c$$$  allocate(
c$$$  &    zxq (ngb,ngb,nw_i:nw),   !,nwp) feb2006
c$$$  &    zxqi(ngb,ngb,niw))
c$$$  zxq=0d0; zxqi=0d0
c$$$  endif

!KINO         write(6,'(a,4(1x,i7))')'kino: nmbas1,nmbas2,nwhis,npm',nmbas1,nmbas2,nwhis,npm
         allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
         allocate( zw0(ngb,ngb) ) !, zxq (ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
         rcxq = 0d0
!! ---------------------------------------------------------------
!! === loop over spin=== =========================================
!! ---------------------------------------------------------------
!TIME1 "before 1003 loop"
!         do 1003 is = 1,nspinmx
         do 1003 is = mpi__Ss,mpi__Se
!TIME0 
            write(6,"(' ### ',2i4,' out of nqibz+n0qi nsp=',2i4,' ### ')") 
     &           iq, is, nqibz + nq0i,nspin
            if(debug) write(6,*)' niw nw=',niw,nw
            isf=is
!!===========tetraini block tetra==.true.===============================1ini
!!     --- tetrahedron method --- preparation
!!     nbnb = total number of weight.
!!     n1b  = band index for occ.   1\ge n1b \ge nband+nctot.
!!     Valence-core order(Core index follows valence index).
!!     n2b  = band index for unocc. 1\ge n2b \ge nband
!!     wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
            if(tetra) then
               write(6,"(' tetra mode niw nw nwhis=',3i5)") niw,nw,nwhis
!!     ekxx1  for rk
!!     ekxx2 for q+rk See tetwt4
               do kx = 1, nqbz
                  call readeval(qbz(:,kx),  is, ekxx1(1:nband, kx) ) 
                  call readeval(q+qbz(:,kx), is, ekxx2(1:nband, kx) )
               enddo

c     takao-feb/2002 I replaced tetwt4 (1d30) with tetwt5(job=0) -----
C     ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
c     the pairs are not dependent on the energy omega
c     in the denominator of the dielectric function.
               write(6,"(' -- First tetwt5 is to get size of array --')")
               job = 0
               if(npm==1) then
                  ncc=0
               else
                  ncc=nctot
               endif
               allocate( demin(nband+nctot,nband+ncc,nqbz,npm),
     &              demax(nband+nctot,nband+ncc,nqbz,npm) )
               allocate( iwgt (nband+nctot,nband+ncc,nqbz,npm) )
c     allocate( demin(nband+nctot,nband,nqbz),
c     &                demax(nband+nctot,nband,nqbz) )
c     allocate(iwgt(nband+nctot,nband,nqbz) )
!     wgt, demin, demax may require too much memory in epsilon mode.
!     We will have to remove these memory allocations in future.
!     tetwt5x_dtet2 can be very slow because of these poor memory allocation.
               if(nctot==0) then
                  deallocate(ecore)
                  allocate(ecore(1,2)) !this is dummry
               endif
               allocate(ibjb(1,1,1,1),ihw(1,1,1),jhw(1,1,1),nhw(1,1,1),whw(1))!dummy
               call tetwt5x_dtet4(npm,ncc,  
     i              q, ekxx1, ekxx2, qbas,ginv,ef, 
     d              ntetf,nqbzw, nband,nqbz,
     i              nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i              job,
     o              iwgt,nbnb,  !job=0
     o              demin,demax, !job=0
     i              frhis, nwhis, ! job=1    not-used
     i              nbnbx,ibjb,nhwtot, ! job=1    not-used
     i              ihw,nhw,jhw, ! job=1    not-used
     o              whw,        ! job=1    not-used
     i              iq,is,is,nqibz, eibzmode,nwgt(:,iq)) 
               deallocate(ibjb,ihw,jhw,nhw,whw) !dummy
               nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
               if(debug) write(6,*)' nbnbx=',nbnbx
               allocate(  n1b(nbnbx,nqbz,npm)
     &              ,n2b(nbnbx,nqbz,npm))
               n1b=0; n2b=0 
               do jpm=1,npm
                  call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx,
     o                 n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
               enddo
c     ... check
               if(debug) then
                  do kx  = 1, nqbz
                     do jpm = 1, npm
c     do ix  = 1, nbnb(kx,jpm)
c     write(6,"('jpm kx ixx n1b n2b=',5i7)")jpm,kx,ix,n1b(ix,kx,jpm),n2b(ix,kx,jpm)
c     enddo
                        write(6,"('jpm kx  minval n1b n2b=',4i5)")jpm,kx,
     &                       minval(n1b(1:nbnb(kx,jpm),kx,jpm)),
     &                       minval(n2b(1:nbnb(kx,jpm),kx,jpm))
                     enddo
                  enddo
               endif
c     
               if(sum(abs(nbnb-nbnbtt))/=0)then
                  do ik=1,nqbz
                     write(6,*)
                     write(6,*)"nbnb  =",nbnb(ik,:)
                     write(6,*)"nbnbtt=",nbnbtt(ik,:)
                  enddo
                  call rx( 'hx0fp0:sum(nbnb-nbnbtt)/=0')
               endif
               noccxv = maxval(noccxvx)
               noccx  = nctot + noccxv
               write(6,*)' Tetra mode: nctot noccxv= ',nctot,noccxv
               deallocate(iwgt)
            endif
c=========end of tetraini block==========================================1end



c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
c     if(ixc==11) then !sf 21May02
C     --- METHOD (tetwt5) for the tetrahedron weight
!     Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
!     takao-feb/2002
            if(frhis(1)/=0d0) call rx( ' hx0fp0: we assume frhis(1)=0d0')
            write(6,*)' ----------------nbnbx nqbz= ',nbnbx,nqbz
!!     ... make index sets
            allocate(ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm))
            ihw=0; nhw=0; jhw=0
            jhwtot = 1
            do jpm =1,npm
               do ik   = 1,nqbz
                  do ibib = 1,nbnb(ik,jpm)
                     call hisrange( frhis, nwhis,  
     i                    demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     i                    demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     o                    ihw(ibib,ik,jpm),nhw(ibib,ik,jpm))
                     jhw(ibib,ik,jpm)= jhwtot
                     jhwtot = jhwtot + nhw(ibib,ik,jpm)
                  enddo
               enddo
            enddo
            nhwtot = jhwtot-1
            write(6,*)' nhwtot=',nhwtot
            deallocate(demin,demax)
            allocate( whw(nhwtot), ! histo-weight
     &           ibjb(nctot+nband,nband+ncc,nqbz,npm) )

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
!!     ... Generate the histogram weights whw
            job=1
            write(6,*) 'goto tetwt5x_dtet4 job=',job
            allocate(demin(1,1,1,1),demax(1,1,1,1),iwgt(1,1,1,1)) !dummy
            call tetwt5x_dtet4(  npm,ncc, 
     i           q, ekxx1, ekxx2, qbas,ginv,ef, 
     d           ntetf,nqbzw, nband,nqbz,
     i           nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i           job,
     o           iwgt,nbnb,     ! job=0
     o           demin,demax,   ! job=0
     i           frhis,nwhis,   ! job=1
     i           nbnbx,ibjb,nhwtot, ! job=1
     i           ihw,nhw,jhw,   ! job=1
     o           whw,           ! job=1
     i           iq,is,is,nqibz, eibzmode,nwgt(:,iq)) 
            deallocate(demin,demax,iwgt) !duumy

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     do ik   = 1,nqbz
c     do ibib = 1,nbnb(ik)
c     ib1  = n1b(ibib,ik)
c     ib2  = n2b(ibib,ik)
c     write(3666,"('--- ik ib jb =',3i4)") ik,ib1,ib2
c     ihww = ihw(ibib,ik)
c     do iw = 1, nhw(ibib,ik)
c     write(3666,"('ihis [a b] whw =',i4,2f8.5,f18.12,d14.6)")
c     &         ihww+iw-1, frhis(ihww+iw-1),frhis(ihww+iw),
c     &         whw(jhw(ibib,ik)+iw-1), whw(jhw(ibib,ik)+iw-1)
c     enddo
c     enddo
c     enddo
c     stop 'hx0fp0.f ---  testend xxx2'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write(1236,"('-----',5i8)") noccxvx
c     do ik  = 1,nqbz
c     write(1236,"('--- ',5i8)") ik, nbnb(ik),nbnbtt(ik)
c     do ibib= 1,nbnb(ik)
c     write(1236,"(2i6,2f12.6)") n1b(ibib,ik),n2b(ibib,ik)
c     &         ,demin(n1b(ibib,ik),n2b(ibib,ik),ik)
c     &         ,demax(n1b(ibib,ik),n2b(ibib,ik),ik)
c     enddo
c     enddo
c     stop 'hx0fp0.m.fal.f ---  testend'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                ! endif                        !ixc==11   !sf 21May02
!! ======TetrahedronWeight_5 block end =========


!TIME1 "before x0kf_v4hz"
!TIME0 
!! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 ==
!!     You will need to pass these arrays to x0kf
!!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
!!     nhw(ibjb,kx): the number of histogram sections
!!     jhw(ibjb,kx): pointer to whw
!!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
!!     : histogram weights for given ib,jb,kx for histogram sections
!!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c            write(6,*) ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2
            iqeibz=iq
            !symmetrize=(is==nspinmx)
            symmetrize=.false.

            call x0kf_v4hz(npm,ncc,   
     i           ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i           n1b,n2b,nbnbx,nbnb, ! use whw by tetwt5 ,
     i           q,  
     i           nspin,is,isf, symmetrize, !
     i           qbas,ginv,  qbz,wbz, 
     i           nblocha ,nlnm,nlnmv,nlnmc,iclass, 
     i           ppb(1,is), 
     i           icore,ncore,
     d           nlmto,nqbz,noccx,noccxv,
     d           nl,nclass,natom,nnc,
     d           nlnmx,mdimx,nbloch,  nwhis,
     i           iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
     i           nqbze,nband,nqibz, 
     o           rcxq,          ! rcxq is the accumulating variable for spins 
     i           nolfco,zzr,nmbas_in, ppovl,zcousq, !nmbas2 is removed.,nmbas1 ppovlz, 
     i           chipm,eibzmode,  !z1offd,
     i           nwgt(:,iqeibz),igx(:,:,iqeibz),igxt(:,:,iqeibz),ngrp, eibzsym(:,:,iqeibz),crpa)
            write(6,*)' end of x0kf_v4h sum rcxq=',sum(abs(rcxq))
!TIME1 "x0vk_v4h"

c===========Check write mode4 block ixc==4 ============6ini
ccccccccccccccccccfaleev 21May02 ixc==4 not used
c     if(ixc==4) then
c     C --- check write ---
c     ichkhis = iopen ('tethis.chk',1,-1,0)
c     do ik   = 1,nqbz
c     do ibib = 1,nbnb(ik)
c     ib1  = n1b(ibib,ik)
c     ib2  = n2b(ibib,ik)
c     ihww = ihw(ibib,ik)
c     write(ichkhis,"('--- ik ib jb =',3i4)") ik,ib1,ib2
c     do iw = 1, nhw(ibib,ik)
c     write(ichkhis, "('his [min max]=',i6,2f7.4,
c     &             ' whw wwk*h= ',2d12.5,'( =',d12.5,'+',d12.5,')')")
c     &          ihww+iw-1, freq_r(ihww+iw-1), freq_r(ihww+iw),
c     !-------------------------------------------------------
c     &          whw(jhw(ibib,ik)+iw-1),                     !weight for the histgram range. by tetwt5
c     !------------ v.s. -------------
c     &          dimag(wwk(ibib,ik,ihww+iw-1)+wwk(ibib,ik,ihww+iw))/2d0
c     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight=weight1+weight2 for each bin of the histgram range by tetwt4.fal
c     !-------------------------------------------------------
c     c
c     &          dimag(wwk(ibib,ik,ihww+iw-1))/2d0
c     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight1 for the histgram range by tetwt4
c     &          dimag(wwk(ibib,ik,ihww+iw))/2d0
c     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw))     !weight2 for the histgram range by tetwt4
c     enddo
c     enddo
c     enddo
c     ichkhis = iclose('tethis.chk')
c     write(6,*)' OK! hx0fp1_tetwt5_test end. Histgram method.'
c     write(6,*)' The final reulst wwk vs.whw is in tethis.chk!'
c     write(6,*)' This code is a routine to test tetwt5.'
c     stop ' ixc==4  Weight Histgram check. wwk vs.whw'
c     endif  ! end of if(tetra)
c===========Check write mode4 block end  ==============6end
            deallocate(ihw,nhw,jhw, whw,ibjb )
            if(tetra)            deallocate( n1b,n2b)
 1003       continue;write(6,*) 'end of spin-loop nwp=',nwp !end of spin-loop
c===========end of spin loop============================================
c=======================================================================

!! symmetrize and convert to Enu basis by dconjg(tranpsoce(zcousq)*rcxq8zcousq if eibzmode
            if(eibzmode)  then
               symmetrize= .true.
               is=1             ! dummy
               call x0kf_v4hz(npm,ncc,   
     i              ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i              n1b,n2b,nbnbx,nbnb, ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf, symmetrize, !
     i              qbas,ginv,  qbz,wbz, 
     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, 
     i              ppb(1,is), 
     i              icore,ncore,
     d              nlmto,nqbz,noccx,noccxv,
     d              nl,nclass,natom,nnc,
     d              nlnmx,mdimx,nbloch,  nwhis,
     i              iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
     i              nqbze,nband,nqibz, 
     o              rcxq,       ! rcxq is the accumulating variable for spins 
     i              nolfco,zzr,nmbas_in, ppovl,zcousq, !nmbas2 is removed.,nmbas1 ppovlz, 
     i              chipm,eibzmode, !z1offd,
     i              nwgt(:,iqeibz),igx(:,:,iqeibz),igxt(:,:,iqeibz),ngrp, eibzsym(:,:,iqeibz),crpa)
            endif

!! reduction rcxq in the S-axis
            write(6,*) 'MPI__AllreduceSumS start'
            do i_reduction_npm=1,npm
              do i_reduction_nwhis=1,nwhis
                do i_reduction_nmbas2=1,nmbas2
                    call MPI__AllreduceSumS(
     .              rcxq(1,i_reduction_nmbas2,i_reduction_nwhis,i_reduction_npm), nmbas1)
                enddo
              enddo
            enddo
            write(6,*) 'MPI__AllreduceSumS end'
!TIME0 

!! --- Hilbert transform.  Genrerate Real part from Imaginary part. ======
            if(allocated(zxq) ) deallocate(zxq,zxqi)
            allocate(zxq (nmbas1,nmbas2,nw_i:nw), zxqi(nmbas1,nmbas2,niw))
            write(6,'("goto dpsion5: "
     &           //"nwhis nw_i niw nw_w nmbas1 nmbas2=",6i5)') nwhis,nw_i,nw,niw,nmbas1,nmbas2
            write(6,*)' -------- nmbas1,nmbas2=', nmbas1,nmbas2
            call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i           rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o           zxq, zxqi,
     i           chipm, schi,is,  ecut,ecuts)
            write(6,*)' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
            if(allocated(rcxq) ) deallocate(rcxq)
!TIME1 "before real omega"

!! ===  RealOmega ===
!TIME0 
            if (realomega) then
               if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
               nwmax = nw
               nwmin = nw_i
               if (anfexist) then !Antiferro accelaration.
                  do iw  = nwmin,nwmax
                     call anfx0k(natom,nclass,nblocha,iclass,pos
     i                    ,nbloch,ngc,
     i                    q, ngveccB,qbas, ! for q+G
     i                    anfvec,ibasf, ! these are antiferro informations.
     i                    zxq(1,1,iw)) ! i/o
                  enddo
               endif
!! prepare for iq0.
               iq0 = iq - nqibz
c               if(newaniso2) then
c$$$                 if( iq==1 ) then 
c$$$                    write(6,*)'open EPS0inv mpi=',mpi__rank 
c$$$                   ifepstinv = iopen('EPS0inv',0,-1,0)
c$$$                   write(ifepstinv) ngb
c$$$                 endif  
               if(iqinit) then
                  allocate( sk(ngb,nwmin:nwmax,nq0i),  sks(ngb,nwmin:nwmax,nq0i) )
                  allocate( skI(ngb,niw,nq0i), sksI(ngb,niw,nq0i))
                  iqinit=.false.
               endif
               allocate(epstilde(ngb,ngb))
               allocate(epstinv(ngb,ngb))
c               endif
!KINO               write(6,*)'kino: nwmin,nwmax,ngb=',nwmin,nwmax,ngb
!TIME0 
!TIME0 
               write(6, *)" === trace check for W-V === nwmin nwmax=",nwmin,nwmax
!TIME1 "after write trace msg"
               do 1015 iw  = nwmin,nwmax
                  frr= dsign(freq_r(abs(iw)),dble(iw))
                  imode = 1
c                  if(newaniso2.and.iq<=nqibz) then !for mmmw
                  if(iq<=nqibz) then !for mmmw
!TIME0 
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
!TIME1 "after matcinv"
!TIME0

c$$$  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$  cmmm direct inversion vs. block inversion
c$$$  if(iq>nqibz) then
c$$$  c direct inversion
c$$$  ix=0
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  do igb1=1+ix,ngb
c$$$  do igb2=1+ix,ngb
c$$$  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$  enddo
c$$$  enddo
c$$$  c              write(6,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$  write(6,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,1d0/epstinv(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$  c             write(6,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(1,2:10:3),zw0(1,63:70:3)
c$$$  c block inversion
c$$$  ix=1
c$$$  do igb1=ix+1,ngb
c$$$  do igb2=ix+1,ngb
c$$$  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$  enddo
c$$$  enddo
c$$$  epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$  call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$  absq=sqrt(sum(q**2*tpioa**2))
c$$$  sk(  1:ngb)= zxq(1,1:ngb,iw)
c$$$  sks( 1:ngb)= zxq(1:ngb,1,iw)
c$$$  w_k(1) =0d0
c$$$  w_ks(1)=0d0
c$$$  w_k( 2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(vcousq(1)*sk(2:ngb)*vcousq(2:ngb),epstinv(2:ngb,2:ngb))
c$$$  w_ks(2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(epstinv(2:ngb,2:ngb),vcousq(1)*sks(2:ngb)*vcousq(2:ngb))
c$$$  llw(iw,iq0)=
c$$$  &             1d0
c$$$  &            -vcousq(1)*sk(1)*vcousq(1) ! sk(1,1,iw)=sks(1,1,iw)=H of Eq.(40).
c$$$  &            -vcousq(1)*vcousq(1)* sum( vcousq(2:ngb)*sk(2:ngb) * matmul(epstinv(2:ngb,2:ngb),sks(2:ngb)*vcousq(2:ngb)))
c$$$  write(6,"('mmmmzwp99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0), !(1d0/llw(iw,iq0)-1d0)*vcousq(1)**2,
c$$$  c     &                  w_k(2:10:3)/llw(iw,iq0), w_k(63:70:3)/llw(iw,iq0)
c$$$  &                  w_ks(2:10:3)/llw(iw,iq0), w_ks(63:70:3)/llw(iw,iq0)
c$$$  write(6,"('mmmmzwp99x ')")
c$$$  endif
c$$$  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     do igb1=1+ix,ngb
                        do igb2=1+ix,ngb
                           zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                           if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                        enddo
                     enddo
c$$$                     if(iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
                     zw(1:ngb,1:ngb) = zw0
!TIME1  "after zw=zw0"
!TIME0
                     if (mpi__rootS)then
                     write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
                     endif
                     call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
!TIME1 "end of tr_chkwrite"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmmmm 
c     if(iq>nqibz) then
c     write(6,"('mmmmz99x ',i3,10(2d13.5,2x))") iw,zw0(1,1)+vcousq(1)**2,zw0(2:10:3,1),zw0(63:70:3,1)
c     endif
c     if(iq==1.or.iq>nqibz) then
c     write(6,"('mmmz0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c     write(6,"('mmmz0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2:10:3,1,iw),zxq(63:70:3,1,iw)
c     write(6,"('mmmmz99x ',i3,10(2d13.5,2x))") iw,zw0(1,1)+vcousq(1)**2,zw0(1,2:10:3),zw0(1,63:70:3)
c     write(6,"('mmmzx  ',2i3,10(2d13.5,2x))") iq,iw,zxq(2,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c     write(6,"('mmmzx  ',2i3,10(2d13.5,2x))") iq,iw,zxq(3,1,iw),zxq(3,2:10:3,iw),zxq(3,63:70:3,iw)
c     write(6,"('mmmzxs ',2i3,10(2d13.5,2x))") iq,iw,zxq(1,1,iw),zxq(2:10:3,1,iw),zxq(63:70:3,1,iw)
c     write(6,"('mmmzxs ',2i3,10(2d13.5,2x))") iq,iw,zxq(1,2,iw),zxq(2:10:3,2,iw),zxq(63:70:3,2,iw)
c     write(6,"('mmmmzee',2i3,10(2d13.5,2x))")iq,iw,epstilde(2,2),epstilde(2,2:10:3),epstilde(2,63:70:3)
c     write(6,"('mmmmzee',2i3,10(2d13.5,2x))")iq,iw,epstilde(3,2),epstilde(3,2:10:3),epstilde(3,63:70:3)
c     endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  endif

c                  if(newaniso2.and.iq>nqibz) then
                  if(iq>nqibz) then
!! Full inversion to calculalte eps with LFC.
                   vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
                   if(localfieldcorrectionllw()) then
                     ix=0
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
c !TIME0
                     epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                     call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
                     llw(iw,iq0)= 1d0/epstinv(1,1)
c !TIME1 "end of matcinv_epstinv" !this gives wrong message, probably 
c         because of a bug of !TIME1 processing for MPI mode.
                   else
c commentout block inversion
c$$$                     sk  (1:ngb,iw,iq0)= zxq(1,1:ngb,iw)
c$$$                     sks (1:ngb,iw,iq0)= zxq(1:ngb,1,iw)
c$$$c                     sks (1:ngb,iw,iq0)= zxq(2,1:ngb,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
                     llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw) 
                   endif 
c$$$

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmmmm 
c     write(6,"('mmmw0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c     write(6,"('mmmw0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c     write(6,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,fourpi/sum(q**2*tpioa**2)/llw(iw,iq0),
c     &                  w_k(2:10:3)/llw(iw,iq0),w_k(63:70:3)/llw(iw,iq0)
c     write(6,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0),
c     &                  w_ks(2:10:3)/llw(iw,iq0),w_ks(63:70:3)/llw(iw,iq0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!TIME0 
                   write(6,"('iq iw_real eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)")iq,iw,llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
!TIME1 "write iq iwreal" 

c     write(6,*) 'uuuu iq iw llw vc x0mean',iq,iw,fourpi/sum(q**2*tpioa**2),sk(1)
c     write(ifisk) iw,iq0,q
c     write(ifisk) vcousq(2:ngb)*w_k( 2:ngb),vcousq(2:ngb)*w_ks( 2:ngb)
                  endif

c$$$                  if(.not.newaniso2) then ! Original mode
c$$$                     call rx( 'not checked here')
c$$$c     call wcf( ngb, vcoul, zxq(1,1,iw), imode, zw0)
c$$$                  endif

c$$$  !!... a debug mode
c$$$  write(6,"('hhh --- EigenValues for Im( W) --------')")
c$$$  allocate(ebb(ngb))
c$$$  call diagcvh2( (zw0-transpose(dconjg(zw0)))/2d0/img, ngb, ebb)
c$$$  do ii=1,ngb
c$$$  if( abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
c$$$  write(6, "('hhhIWq : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
c$$$  else
c$$$  write(6, "('hhhIWqxxx : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
c$$$  endif
c$$$  enddo
c$$$  deallocate(ebb)

c     if(newaniso2.and.iq>nqibz) then
c     c                zw(1:ngb,1:ngb) = 0d0
c     c                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
c     else
c     zw(1:ngb,1:ngb) = zw0
c     c                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
c     write(ifrcw, rec= iw-nw_i+1 ) zw   !  WP = vsc-v
c     call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
c     endif
 1015          continue         !iw
!TIME1 "end 1015 loop"

c     if(newaniso2) then
c     if(allocated(sk)) deallocate(sk,sks,w_k,w_ks)
c     endif
               if( allocated(zzr) ) deallocate(zzr)
            endif 
!! === RealOmega end ===
!TIME1 "real omega"


!! === ImagOmega ===
!TIME0 
            if (imagomega) then
               write(6,*)' goto imag omega'
               if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
c     ... Antiferro accelaration.
               if (anfexist) then
                  do iw  = 1,niw
                     call anfx0k(natom,nclass,nblocha,iclass,pos
     i                    ,nbloch,ngc,
c     i      q, ngveccB(1,1,iqx),qbas,    ! for q+G 12jan2004
     i                    q, ngveccB,qbas, ! for q+G
     i                    anfvec,ibasf, ! these are antiferro informations.
     i                    zxqi(1,1,iw)) ! i/o
                  enddo
               endif
               imode=1
               do 1016 iw  = 1,niw
c                  if( newaniso2 .and. iq<=nqibz ) then
                  if( iq<=nqibz ) then
!!  Eqs.(37),(38) in PRB81 125102
                     if(iq==1) then
                        ix=1
                        zw0(:,1)=0d0
                        zw0(1,:)=0d0
                     else
                        ix=0
                     endif
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
c$$$                     if(iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw

                     zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
                     if (mpi__rootS) then
                     write(ifrcwi, rec= iw)  zw !  WP = vsc-v
                     endif
                     call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
                  endif

c                  if( newaniso2.and.iq>nqibz) then
                  if(iq>nqibz) then
!! Full inversion to calculalte eps with LFC.
                      vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
                     if(localfieldcorrectionllw()) then
                      ix=0
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
                     else
c commentout block inversion
c$$$                     skI  (1:ngb,iw,iq0)= zxqi(1,1:ngb,iw)
c$$$c                     sksI (1:ngb,iw,iq0)= zxqi(2,1:ngb,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     sksI (1:ngb,iw,iq0)= zxqi(1:ngb,1,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
c$$$                     vcou1sq= sqrt(vcou1)
c$$$!! llwI without LFC. LFC contribution is added in 
                      llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
                     endif  
                     write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)")iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
                  endif

 1016          continue
c               if(newaniso2) then
c$$$                 if(iq==1) ifepstinv = iclose('EPS0inv') !iq==1 close write mode.
                 deallocate(epstinv)
                 if(allocated(epstilde)) deallocate(epstilde)
c               endif
            endif 
!! === ImagOmega end ===
!TIME1 "imag omega"

c     1002 continue  ! end of frequency block-loop
            if(allocated(vcoul)) deallocate(vcoul)
            if(allocated(zw0)) deallocate(zw0)
            if(allocated(zxq )) deallocate(zxq)
            if(allocated(zxqi)) deallocate(zxqi)

            if (mpi__rootS) then
            ifrcwi = iclose('WVI.'//charnum5(iq))
            ifrcw  = iclose('WVR.'//charnum5(iq))
            endif
!!  
 1001    continue
c============end of loop over q point =================================
c=======================================================================
!TIME0
         call mpi__barrier()
!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
c         write(6,*)'mpi__sizerrr=',mpi__size,mpi__rank,mpi__root,mpi__size,nqibz,iqxend
ccccccccccccccccccccccccccc
CYY!$OMP parallel 
CYY!$OMP master
         write(6,*)' eeeeeeeeeeeeeeeee111 mpi__rank=',mpi__rank
         if(mpi__size/=1) then
            do iq=nqibz+1,iqxend
               iq0 = iq - nqibz
c     write(6,*)' iq iq0 mpi_rank mpi_ranktab(iq)=',iq, iq0,mpi__rank,mpi__ranktab(iq),mpi__root,nw,nw_i,niw
               if(mpi__Qranktab(iq)/=0) then !jan2012
               if(mpi__Qranktab(iq) == mpi__rankQ) then
c     write(6,*)' mpi_send iq from',iq,mpi__ranktab(iq)
c     write(6,*)' send llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
c     do i=nw_i,nw
c     write(6,*)'sendxxx',i,llw(i,iq0)
c     enddo   
c     write(6,*)' send llwI sum=',sum(abs(llwI(:,iq0))),niw
                  dest=0
                  call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
                  call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
               elseif(mpi__rootQ) then
c     write(6,*)' mpi_recv iq from',iq,mpi__ranktab(iq),nw,nw_i,niw
                  src=mpi__Qranktab(iq)
                  call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
                  call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
c     do i=nw_i,nw
c     write(6,*)'recivxxx',i,llw(i,iq0)
c     enddo   
c     write(6,*)' recv llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
c     write(6,*)' recv llwI sum=',sum(abs(llwI(:,iq0))),niw
               endif
               endif
            enddo  
         endif
         write(6,*)' eeeeeeeeeeeeeeeee222 mpi__rank=',mpi__rank
CYY!$OMP end master
CYY!$OMP end parallel 

!TIME1 "end eeeeeee"

c commentout block inversion
c$$$!! Add LFC (local field correction) to llw and llwI
c$$$         if(newaniso2 .and. mpi__rank == 0 ) then ! only on root node
c$$$            iq=1 !for q=0
c$$$            vcoudfile='Vcoud.'//charnum5(iq)
c$$$            ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$            read(ifvcoud) ngb0
c$$$            read(ifvcoud) qvv
c$$$            if(sum(abs(qvv))>1d-10) then
c$$$               write(6,*)'qvv =',qvv
c$$$               stop 'hx0fp0: qvv/=0 hvcc is not consistent'
c$$$            endif
c$$$            if(allocated(zcousq0)) deallocate( zcousq0,vcousq0 )
c$$$            allocate( zcousq0(ngb0,ngb0),vcousq0(ngb0))
c$$$            read(ifvcoud) vcousq0
c$$$            read(ifvcoud) zcousq0
c$$$            idummy=iclose(trim(vcoudfile))
c$$$            vcousq=sqrt(vcousq)
c$$$            allocate(epstinv(ngb0,ngb0),w_k(ngb0),w_ks(ngb0),w_kI(ngb0),w_ksI(ngb0),eemat(ngb0,ngb0))
c$$$
c$$$            do iq0=1,nq0i
c$$$              iq = iq0 + nqibz
c$$$              q = qibze(:,iq)
c$$$
c$$$              vcoudfile='Vcoud.'//charnum5(iq)
c$$$              ifvcoud = iopen(trim(vcoudfile),0,-1,0)
c$$$              read(ifvcoud) ngb
c$$$              read(ifvcoud) qvv
c$$$              if(sum(abs(qvv-q))>1d-10) then
c$$$               write(6,*)'qvv =',qvv
c$$$               stop 'hx0fp0: qvv/=0 hvcc is not consistent'
c$$$              endif
c$$$              if(allocated(zcousq)) deallocate(zcousq)
c$$$              if(allocated(vcousq)) deallocate(vcousq)
c$$$              allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
c$$$              read(ifvcoud) vcousq
c$$$              read(ifvcoud) zcousq
c$$$              idummy=iclose(trim(vcoudfile))
c$$$              vcousq=sqrt(vcousq)
c$$$
c$$$              ifepstinv = iopen('EPS0inv',0,0,0)
c$$$              read(ifepstinv) ngb
c$$$
c$$$               ngc=ngb-nbloch
c$$$               if(allocated(ppovlz)) deallocate(ppovlz)
c$$$               if(allocated(ppovl)) deallocate(ppovl)
c$$$               allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
c$$$               call readppovl0(q,ngc,ppovl) !q was qq
c$$$               ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
c$$$               ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
c$$$
c$$$!  eemat: Z\mu_i(\bfk=0)^* <i|j> Z\nu_j(\bfk) 
c$$$               eemat =matmul(transpose(dconjg(zcousq0)),matmul(ppovlz,zcousq))
c$$$               vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
c$$$               vcou1sq = vcou1**.5
c$$$               write(6,*)
c$$$
c$$$              do iw=nwmin,nwmax
c$$$                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$                epstinv(2:ngb,2:ngb) = matmul( transpose(dconjg(eemat(2:ngb,2:ngb))),
c$$$     &                                matmul(epstinv(2:ngb,2:ngb),eemat(2:ngb,2:ngb)) )
c$$$                if(iw/=iwx) then
c$$$                write(6,*)'iw iwx=',iw,iwx
c$$$                stop 'hx0fp0_sc: iw/=iwx'
c$$$                endif
c$$$                w_k(2:ngb) = vcou1sq*matmul( epstinv(2:ngb,2:ngb), sk(2:ngb,iw,iq0)*vcousq(2:ngb))
c$$$                epslfc = -vcou1sq*sum( sks(2:ngb,iw,iq0) * w_k(2:ngb) *vcousq(2:ngb) )
c$$$                llw(iw,iq0) = llw(iw,iq0)  + epslfc
c$$$                write(6,"('eps(on real) iq iw',2i4,2f9.3,2x,2f9.3)") iq0,iw, llw(iw,iq0)-epslfc,llw(iw,iq0)
c$$$              enddo
c$$$              do iw=1,niw
c$$$                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c$$$                if(iw/=iwx) then
c$$$                 write(6,*)'iw iwx=',iw,iwx
c$$$                 stop 'hx0fp0_sc: iw/=iwx'
c$$$                endif
c$$$                w_kI(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), skI(2:ngb,iw,iq0)*vcousq(2:ngb))
c$$$                epslfc=- vcou1sq*sum( sksI(2:ngb,iw,iq0)* w_kI(2:ngb)*vcousq(2:ngb) )
c$$$                llwI(iw,iq0)= llwI(iw,iq0)+epslfc 
c$$$                write(6,"('eps(on img ) iq iw',2i4,2f9.3,2x,2f9.3)")iq0,iw, llwI(iw,iq0)-epslfc,llwI(iw,iq0)
c$$$              enddo
c$$$              ifepstinv = iclose('EPS0inv')
c$$$           enddo
c$$$         endif


!! == W(0) divergent part and W(0) non-analytic constant part.==
 1191    continue
!TIME0 
c         if(newaniso2 .and. mpi__rank == 0 ) then ! MIZUHO-IR only on root node
         if(mpi__rank == 0 ) then ! MIZUHO-IR only on root node

            if(ixc==1011) then !this is only for test.
            ifw0w0i = iopen('W0W0I',0,-1,0) 
            read(ifw0w0i) nw_i,nw,niw,nq0i
            write(6,*)'w0w0i: n=',nw_i,nw,niw,nq0i
            read(ifw0w0i) llw(nw_i:nw,1:nq0i)
            read(ifw0w0i) llwI(1:niw,1:nq0i)
c            read(ifw0w0i) w0(nw_i:nw)
c            read(ifw0w0i) w0i(1:niw)
            ifw0w0i = iclose('W0W0I')
            endif  
         
            write(6,*)
            write(6,*)' ==== newaniso2 mode W(0) divergent part ==== '
!! == W(0) divergent part ==
!! getw0 routine...
!!NOTE: we usually only use lxklm=1 --> this should be stable.
!!  EPSwklm is generated in gwsrc/mkqg.F
            ifidmlx = iopen('EPSwklm',0,0,0)
            read(ifidmlx) nq0i,lxklm
            allocate( dmlx(nq0i,9))
            allocate( epinvq0i(nq0i,nq0i),epinv(3,3,nq0i))
            nlxklm=(lxklm+1)**2
            allocate( wklm(nlxklm))
            read(ifidmlx) dmlx, epinv,epinvq0i
            read(ifidmlx) wklm
            ifidmlx = iclose('EPSwklm')
!! starting from llw(iw,iq0),llwI(iw,iq0)
!! === <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
            allocate(w0(nw_i:nw),w0i(niw))
c     write(6,*)' goto getw0 nq0i epinvq0i=',nq0i,epinvq0i
!! wbz(1) is the weight for q=0 = 1/(n1*n2*n3)
c            write(6,*)'wbz=',wbz
            call getw0(llw, nw_i,nw,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0)
            call getw0(llwI,1,niw  ,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0i)
            if(ixc/=1011) then
            ifw0w0i = iopen('W0W0I',0,-1,0)
            write(ifw0w0i) nw_i,nw,niw,nq0i
            write(ifw0w0i) llw(nw_i:nw,1:nq0i)
            write(ifw0w0i) llwI(1:niw,1:nq0i)
            write(ifw0w0i) w0(nw_i:nw)
            write(ifw0w0i) w0i(1:niw)
            ifw0w0i = iclose('W0W0I')
            endif

            do i=nw_i,nw
               write(6,"('w0 =',i4,2f13.4)")i,w0(i)
            enddo   
            do i=1,niw
               write(6,"('w0i=',i4,2f13.4)")i,w0i(i)
            enddo   
c     write(6,*)'sumcheck w0,w0i=',sum(abs(w0)),sum(abs(w0i))
!! === w0,w0i are stored to zw for q=0 ===
!! === w_ks*wk are stored to zw for iq >nqibz ===
! We assume iq=1 is for rank=0
            do iq = 1,1         !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
c               if( .not. mpi__task(iq) ) cycle
               q = qibze(:,iq)
               do ircw=1,2
                  if    (ircw==1) then
                     nini=nw_i
                     nend=nw
                     ifrcwx = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
                  elseif(ircw==2) then;  nini=1;      nend=niw;
                     ifrcwx = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
                  endif
                  do iw=nini,nend
c     if(iq<=nqibz) read(ifrcwx, rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
                     read(ifrcwx, rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
                     if( iq==1 ) then
                        if(ircw==1) zw(1,1) = w0(iw)
                        if(ircw==2) zw(1,1) = w0i(iw)
                     endif
c     write(ifrcwx,rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
                     write(ifrcwx,rec=iw-nini+1) zw !(1:ngb,1:ngb)
                  enddo
                  if    (ircw==1) then
                     ifrcwx = iclose('WVR.'//charnum5(iq))
                  elseif(ircw==2) then
                     ifrcwx = iclose('WVI.'//charnum5(iq))
                  endif
               enddo
            end do
         end if
!TIME1 "before closing hbe.d"
      is = iclose('hbe.d')
 999  continue
!TIME1 "end program"
!TIMESHOW
      call cputid(0)
      write(6,*) '--- end of hx0fp0_sc --- irank=',mpi__rank
      call flush(6)
      call MPI__Finalize
      if(ixc==11) call rx0( ' OK! hx0fp0_sc ixc=11 Sergey F. mode')
      if(ixc==1011) call rx0( ' OK! hx0fp0_sc ixc=1011 W0W0Ionly')
      end program hx0fp0_sc 


C===================================================================
      subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
      implicit none
      integer:: nblochpmx,nbloch,ngb,iw,i,iq
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
      enddo                     !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
c     do i = 1,ngb
c     write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c     enddo
      end




!!------------------------------------------------------------------------
      module m_freq
      real(8),allocatable:: frhis(:)

      contains
      subroutine(realomega,imagomega,tetra, nband,qbze,nqbze,nspin,
     o  frhis,freq_r,freq_i)
!!- frequency histgram, real axis points, im axis points.
      implicit none
      logical:: realomega,imagomega,tetra
      integer::
      real(8)::omg_c
      real(8) :: Emin, Emax      ,emax2,emin2
!! We first accumulate Imaginary parts. 
!! Then it is K-K transformed to obtain real part.
      call getkeyvalue("GWinput","omg_c",omg_c )
      call getkeyvalue("GWinput","omg_c",omg_c )
      write(6,*)'dw, omg_c= ',dw, omg_c
      call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
c      if(bzcase()==2.or.(.not.qbzreg())) then
      if(.not.qbzreg()) then
         allocate(qbz2(3,nqbz))
         do iq=1,nqbz
            qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
         enddo
         call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
         emax=max(emax,emax2)
         emin=min(emin,emin2)
         deallocate(qbz2)
      endif
      if (nctot .gt. 0) Emin=minval(ecore(:,1:nspin))
      omg2max=(Emax-Emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
      write(6,*)' emin emax omega2max=',emin, emax, omg2max
      write(6,*)

!TIME1 "before histogram division"
!TIME0 
ccccccccccccccccccccccccccccccccchistogram divisions cccccccccccccccccc

      nw2=int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
      allocate(freqr2(nw2))     !+1 b/c (iw-1)
      do iw=1,nw2
         freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
      enddo                     !linear + quadratic term
      if (nw2 < 2 ) call rx( "hx0fp0: nw2 < 2")
      if (dw*(nw-2) > freqr2(nw2-1)) call rx("hx0fp0: dw*(nw-2) > freqr2(nw2-1)")
      nwhis = nw2-1
      allocate(frhis(1:nwhis+1))
      frhis = freqr2(1:nwhis+1)
      write(6,*)' we set frhis nwhis=',nwhis
!!     ... for frhis_m
      nw=nw2-1                  ! nw+1 is how many points of real omega we use
                                ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
                                ! maximum nw=nw2-1 because nwhis=nw2-1
      do iw=3,nw2-1
                                !nw is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
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
            exit                ! 'nw_input' is only used to get maximum frequency for
                                !  dressed coulomb line
         endif
      enddo
      allocate(freq_r(0:nw))
      freq_r(0)=0d0
      do iw=1,nw
         freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
      enddo
      nwp = nw+1
      delta=0d0                 !only for writeemesh

!! Plot frhis ----------
          if(onceww(1)) then
            write(6,*)' --- Frequency bins to accumulate Im part  (a.u.) are ---- '
            do ihis= 1, min(10,nwhis)
              write(6,"(' ihis Init  End=', i5,2f13.6)") ihis,frhis(ihis),frhis(ihis+1)
            enddo
            write(6,*) 'ihis  ...'
            do ihis= max(min(10,nwhis),nwhis-10), nwhis
              write(6,"(' ihis Init  End=', i5,2f13.6)") ihis,frhis(ihis),frhis(ihis+1)
            enddo
          endif

!!     ... timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
         write(6,"('TimeReversal off mode')")
         npm=2
         nw_i=-nw
         if(.not.tetra)      call rx( ' tetra=T for timereversal=off')
      endif
      write(6,*)'Timereversal=',Timereversal()

!TIME1 "before realomega write3111"
!TIME0
      if(realomega .and. mpi__root) then
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

C     --- set freq_i
      if (imagomega) then
         write(6,*)' freqimg: niw =',niw
         allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
         call freq01 (niw,ua,   ! Generate gaussian frequencies x between (0,1) and w=(1-x)/x      
     o        freqx,freq_i,wx,expa)
         deallocate(freqx,wx,expa)
      endif

