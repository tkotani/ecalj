      program hx0fp0_sc
!!  Calculate W-V for QSGW mode.
!! We calculate chi0 by the follwoing three steps.
!!  tetwt5: tetrahedron weights
!!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
!!  dpsion5: calculate real part by the Hilbert transformation from the Im part

c      use m_readeps,only: read_eps, epsinv, w_mu, llmat2=>llmat,deallocate_eps
!!

      use m_ReadEfermi,only: readefermi,ef
      use m_readqg,only: readngmx,readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval
      use m_read_bzdata,only: read_bzdata, !<--- 'call read_bzdata' sets up following data.
     &   ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv,
     &   dq_,qbz,wbz,qibz,wibz,
     &     ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
c     &     idteti, nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     &     nclass,natom,nspin,nl,nn,ngrp,
     &     nlmto,nlnmx, nctot,niw, !nw_input=>nw,
     &     alat, delta,deltaw,esmr,symgrp,clabl,iclass, !diw,dw,
     &     invg, il, in, im, nlnm, 
     &     plat, pos, ecore, symgg 

      use m_keyvalue,only: getkeyvalue
      use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
      use m_readqgcou,only: readqgcou

!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: rdpp,    !"call rdpp" generate following data.
     & nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr
!! Generate matrix element for "call get_zmelt".
      use m_zmel,only:    !these data set are stored in this module, and used when 
     & nband,itq,ngcmx,ngpmx,    ppovlz,
     & ppbir,shtvg, miat,tiat , ntq
!! frequency
      use m_freq,only: getfreq,
     &   frhis,freq_r,freq_i, nwhis,nw_i,nw,npm !output of getfreq
!! antiferro
c      use m_anf,only: anfcond,
c     & laf,ibasf !,ldima,pos,natom
!! tetwt
      use m_tetwt,only: tetdeallocate,gettetwt, !followings are output of 'L871:call gettetwt')
     &  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
!! w0 and w0i (head part at Gamma point)
      use m_w0w0i,only: w0w0i,
     & w0,w0i,llmat

!! MPI
      use m_mpi,only: MPI__hx0fp0_rankdivider2Q,MPI__hx0fp0_rankdivider2S,
     &   MPI__Qtask,MPI__InitializeQSPBM,MPI__Finalize,MPI__root,
     &   MPI__Broadcast,MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__rank,MPI__size,
     &   MPI__Qranktab,MPI__consoleout,MPI__Ss,MPI__Se, MPI__allreducesumS,
     &   MPI__barrier, MPI__rankQ,MPI__rootQ,MPI__rootS
!! q0p
      use m_readq0p,only: readq0p,
     & wqt,q0i,nq0i ,nq0iadd,ixyz

      implicit none
      integer,allocatable:: nwgt(:,:)
      integer::iopen,maxocc2,iclose, ixc,iqxini,iqxend,
     &     ifhbe,  nprecb,mrecb,mrece,nlmtot,nqbzt,!nband,
     &     i,nq0ix,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx !ngcmx,
     &     ,nblochpmx,ndummy1,ndummy2,ifcphi,is,nwp, !ifvcfpout,,mdimx,nbloch
     &     ifepscond,nxx,ifvxcpout,ifgb0vec
     &     ,nw0,iw,ifinin,iw0,noccxv,noccx
     &     ,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &     ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &     ,kx,isf,kqxx,kp,job,nwmax !,ifev1,ifev2 !,nhwtot
     &     ,ihis,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &     ,  ifchipmlog ,   nw_w,nwmin  ! ,ngpmx
      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm, dwry,dwh,omg2, q(3),  qgbin(3),qx(3) 
      real(8):: ua=1d0          ! this is a dummy.
      integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2), ndble=8, nword
      integer,allocatable :: ngveccB(:,:), iqib(:),ifppb(:) !,lx(:) ngvecc(:,:),
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &     zw(:,:),zw0(:,:), zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: eqt(:), !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:)
     &     ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
     &     ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), qbze(:,:),qibze(:,:) !,ecore(:,:)
c     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: rcxq(:,:,:,:)
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      real(8) ::qbzx(3)
      logical :: debug
c      integer,allocatable:: ibasf(:)
      logical :: realomega, imagomega 
      complex(8),allocatable:: zzr(:,:),ppovl(:,:),ppovlzinv(:,:) !,ppovlz(:,:)
      complex(8) :: epxxx,vcmean
      character*9 fileps
      character*15 filepsnolfc
c      logical :: paralellx0=.true. !, hist
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
      real(8) :: omg2max,omg1max,wemax
      logical::imagonly=.false. , noq0p !,readgwinput
      integer::nwin, incwfin, verbose,nbcut,nbcut2,ifpomat,nnmx,ikpo,nn_,noo,iqxxx,nomx
c      real(8)::efin
      logical :: nolfco=.false.
      integer:: isp1,isp2, ngc,mrecg ! bzcase,
      real(8)::  quu(3),deltaq(3),qqq(3)=0d0 !
      complex(8),allocatable:: wgt(:,:,:)
      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg 
!    logical ::smbasis !smbasis will be implemented in m_zmel.f which generates <phi|phi M>
      real(8):: q_r(3)
      complex(8),allocatable:: pomat(:,:)
      logical   :: timereversal,onceww
      integer :: jpm,ncc
      real(8) :: frr !, sciss
      integer :: ngb0,ifvcoud,idummy,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig,nmbas1x !ifepstinv,
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:),eemat(:,:),zcousq0(:,:)
      real(8),allocatable:: vcousq(:),vcousq0(:),vcoudummy(:)
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq
!! Eq.(40) in PRB81 125102
c     complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:), 
     &     w_k(:),w_ks(:),w_kI(:), w_ksI(:) 
      complex(8),allocatable:: llw(:,:), llwI(:,:),aaamat(:,:)
      integer:: lxklm,nlxklm,ifidmlx,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
      complex(8):: vc1vc2
      integer,allocatable:: neibz(:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
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
      logical :: symmetrize,eibzmode
      real(8):: schi=-9999 !dummy
      integer:: i_reduction_npm, i_reduction_nwhis,  i_reduction_nmbas2
      logical:: crpa
      integer,allocatable :: iclasst(:), invgx(:)
      integer:: ibasx,ificlass,ifile_handle,ifiq0p
      complex(8),allocatable:: ppovl_(:,:)
      logical:: tetra ,readw0w0itest=.false.
      integer::nw_ixx,nwxx

      logical:: w4pmode
      complex(8),allocatable:: wmu(:,:),wmuk(:,:)
      integer:: ifw4p,ngbq0,igb
      real(8):: qv(3,3)

      real(8)::ebmx
      integer:: nbmx,mtet(3)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

      logical:: eginit=.true.
      real(8),allocatable,save:: gfmat(:,:)
      complex(8),allocatable:: rcxqin(:)
      real(8):: egauss
      integer:: imbas1,imbas2,ipm
c      integer:: ifief
c      real(8):: ef
!-------------------------------------------------------------------------
      call MPI__InitializeQSPBM()
       call realtimediff(1,'')
       call realtimediff(3,'')
      call MPI__consoleout('hx0fp0_sc')
      call cputid (0)
C--------------------------
C var 0 zzr(1,1)
C name= zzr
C dim 0 1 1
C dim 1 1 1
C--------------------------
      call add_alloclist("zzr",sizeof(zzr),
     & int((1)*(1),kind=8))  !omplevel0
      allocate( zzr(1,1)) !dummy
      hartree= 2d0*rydberg()
      pi     = 4d0*datan(1d0)
      fourpi = 4d0*pi
      sqfourpi=sqrt(fourpi)
      write(6,*) ' --- hx0fp0_sc Choose omodes below ----------------'
      write(6,*) '  ixc= 11,10011,or 1011 '
      write(6,*) ' --- Put number above ! -----------------'
      if( MPI__root ) then
         read(5,*) ixc !c     call readin5(ixc,iqxini,iqxend)
      end if
      call MPI__Broadcast(ixc)
      crpa=.false.
      if(ixc==0) call rx( ' --- ixc=0 --- Choose computational mode!')
c      call headver('hx0fp0_sc',ixc)
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
c    Prof.Nagara says this cause a stop in ifort --->why???
c    write(6,*)'Timereversal=',Timereversal()

!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
      call read_BZDATA()

       call realtimediff(3,"readbzdata")
       call realtimediff(5,'')
!! Use regular mesh even for bzcase==2 and qbzreg()=T
!!     off-regular mesh for bzcase==1 and qbzreg()=F
c      if( ( bzcase()==2.and.qbzreg() )       .or.
c     &     ( bzcase()==1.and.(.not.qbzreg()))      ) then
!! this mechanism for qbzreg=F is too complicated. We may need to modify difinition of qbz for qbzreg=F.
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
c      nwin   = 0                !Readin nw from NW file
      incwfin= 0                !use ForX0 for core in GWIN
c      efin =  0d0               !readin EFERMI
c--- EFERMI
      call readefermi()
      call genallcf_v3(incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
c      nw_input = nw ; 
c      write(6,*) 'nw delta=',nw_input,delta
      debug=.false.
      if(verbose()>=100) debug=.true.
      if(debug) write(6,*)' end of genallc'
      tpioa=2d0*pi/alat
!!!!  WE ASSUME iclass(iatom)= iatom,  nclass = natom.  !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) call rx( ' hx0fp0_sc: nclass /= natom ')
!! --- tetra or not
c      if(delta <= 0d0) then
      tetra =  .true.
      delta = -delta
      write(6,*)' hx0fp0.sc: tetrahedron mode delta=',delta
c      else
c         tetra = .false. ! switch for tetrahedron method for dielectric functions
c      endif
!! --- read dimensions of h,hb
      ifhbe      = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if(nlmto/=nlmtot) call rx( ' hx0fp0: nlmto/=nlmtot in hbe.d')
      if(nqbz /=nqbzt ) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')

!! --- Readin Offset Gamma --------
      call readq0p()
      write(6,"(' ### nqibz nq0i nq0iadd=', 3i5)")nqibz,nq0i,nq0iadd

c$$$      if(.not.newaniso2) then
c$$$         wqtsum = sum(abs(wqt(1:nq0i)))
c$$$         call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
c$$$      endif 

       call realtimediff(5,"Q0P")
       call realtimediff(7,'')
      call getsrdpp2(nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx
      nqbze  = nqbz *(1 + nq0i+nq0iadd)
      nqibze = nqibz + nq0i+nq0iadd
C--------------------------
C var 0 qbze(3,nqbze)
C name= qbze
C dim 0 3 3
C dim 1 nqbze nqbze
C--------------------------
      call add_alloclist("qbze",sizeof(qbze),
     & int((3)*(nqbze),kind=8))  !omplevel0
C--------------------------
C var 1 qibze(3,nqibze)
C name= qibze
C dim 0 3 3
C dim 1 nqibze nqibze
C--------------------------
      call add_alloclist("qibze",sizeof(qibze),
     & int((3)*(nqibze),kind=8))  !omplevel0
      allocate( qbze(3, nqbze), qibze(3, nqibze))
      call dcopy(3*nqbz, qbz,  1, qbze,1)
      call dcopy(3*nqibz,qibz, 1, qibze,1)
      do i = 1,nq0i+nq0iadd
         qibze(:,nqibz+i)  = q0i(:,i)
         ini = nqbz*(1 + i -1)
         do ix=1,nqbz
            qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
         enddo
      enddo
!! ---------------- dummy ngrpx=1 -------------------
      ngrpx = 1
      l2nl=2*(nl-1)
C--------------------------
C var 0 symope(3,3)
C name= symope
C dim 0 3 3
C dim 1 3 3
C--------------------------
      call add_alloclist("symope",sizeof(symope),
     & int((3)*(3),kind=8))  !omplevel0
      allocate(symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
!!  dummy. Get space-group transformation information. See header of mptaouof.
      ificlass=ifile_handle()
      open (ificlass,file='CLASS')
C--------------------------
C var 0 iclasst(natom)
C name= iclasst
C dim 0 natom natom
C--------------------------
      call add_alloclist("iclasst",sizeof(iclasst),
     & int((natom),kind=8))  !omplevel0
C--------------------------
C var 1 invgx(ngrp)
C name= invgx
C dim 0 ngrp ngrp
C--------------------------
      call add_alloclist("invgx",sizeof(invgx),
     & int((ngrp),kind=8))  !omplevel0
C--------------------------
C var 2 miat(natom,ngrp)
C name= miat
C dim 0 natom natom
C dim 1 ngrp ngrp
C--------------------------
      call add_alloclist("miat",sizeof(miat),
     & int((natom)*(ngrp),kind=8))  !omplevel0
C--------------------------
C var 3 tiat(3,natom,ngrp)
C name= tiat
C dim 0 3 3
C dim 1 natom natom
C dim 2 ngrp ngrp
C--------------------------
      call add_alloclist("tiat",sizeof(tiat),
     & int((3)*(natom)*(ngrp),kind=8))  !omplevel0
C--------------------------
C var 4 shtvg(3,ngrp)
C name= shtvg
C dim 0 3 3
C dim 1 ngrp ngrp
C--------------------------
      call add_alloclist("shtvg",sizeof(shtvg),
     & int((3)*(ngrp),kind=8))  !omplevel0
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
C--------------------------
C var 0 itq(ntq)
C name= itq
C dim 0 ntq ntq
C--------------------------
      call add_alloclist("itq",sizeof(itq),
     & int((ntq),kind=8))  !omplevel0
      allocate(itq(ntq))
      do i=1,ntq
        itq(i)=i
      enddo
!! Pointer to optimal product basis
c      allocate(imdim(natom))
c      call indxmdm (nblocha,nclass,iclass,natom,
c     o imdim )                  !in m_zmel
      nblochpmx = nbloch + ngcmx
C--------------------------
C var 0 ngveccB(3,ngcmx)
C name= ngveccB
C dim 0 3 3
C dim 1 ngcmx ngcmx
C--------------------------
      call add_alloclist("ngveccB",sizeof(ngveccB),
     & int((3)*(ngcmx),kind=8))  !omplevel0
      allocate(ngveccB(3,ngcmx))
      iqxend = nqibz + nq0i + nq0iadd
      write(6,*) ' nqibz nqibze=',nqibz,nqibze
       call realtimediff(7,"mptauof")
       call realtimediff(9,'')
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

!! We get frhis,freq_r,freq_i, nwhis,nw,npm  by getfreq
      realomega = .true.
      imagomega = .true.
      tetra     = .true.
      call findemaxmin(nband,qbze,nqbze,nspin, emax,emin)
      if(.not.qbzreg()) then
C--------------------------
C var 0 qbz2(3,nqbz)
C name= qbz2
C dim 0 3 3
C dim 1 nqbz nqbz
C--------------------------
      call add_alloclist("qbz2",sizeof(qbz2),
     & int((3)*(nqbz),kind=8))  !omplevel0
         allocate(qbz2(3,nqbz))
         do iq=1,nqbz
            qbz2(:,iq)=qbz(:,iq)+dq_
         enddo
         call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
         emax=max(emax,emax2)
         emin=min(emin,emin2)
C--------------------------
C var 0 qbz2
C name= qbz2
C--------------------------
      call del_alloclist( "qbz2") !omplevel0
         deallocate(qbz2)
      endif
      if (nctot > 0) Emin=minval(ecore(:,1:nspin))
      omg2max = (Emax-Emin)*.5d0+.2d0 
             ! (in Hartree) covers all relevant omega, +.2 for margin
      if(MPI__root) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
      call getwemax(.true.,wemax) !wemax is to determine nw !real axis divisions
      if(MPI__root) write(6,"(' wemax=  ',f13.4)") wemax
      call getfreq(.false.,realomega,imagomega,tetra,omg2max,wemax,niw,ua,MPI__root)
      nwp = nw+1

!! We first accumulate Imaginary parts. Then do K-K transformation to get real part.
      noccxv = maxocc2 (nspin,ef, nband, qbze,nqbze) 
        !max no. of occupied valence states
      if(noccxv>nband) call rx( 'hx0fp0_sc: all the bands filled! too large Ef')
      noccx  = noccxv + nctot
      nprecx = ndble        !We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
      if(MPI__root)then
        ifwd   = iopen('WV.d',1,-1,0)
        write (ifwd,"(1x,10i14)") !"(1x,i3,i8,i5,5i4)") 
     &  nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
        ifwd = iclose('WV.d'); ifwd=0
      endif
C--------------------------
C var 0 zw(nblochpmx,nblochpmx)
C name= zw
C dim 0 nblochpmx nblochpmx
C dim 1 nblochpmx nblochpmx
C--------------------------
      call add_alloclist("zw",sizeof(zw),
     & int((nblochpmx)*(nblochpmx),kind=8))  !omplevel0
      allocate(  zw(nblochpmx,nblochpmx) )
      nspinmx = nspin
       call realtimediff(9,"init_readeigen")
       call realtimediff(11,'')

!!... these are used x0k
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
!! This is general for rotated CG coefficient; but hx0fp0 mode is only for  ngrpx=1 (not rotated).
!! Compare usage in hsfp0 modes.
      irot=1
C--------------------------
C var 0 ppbir(nlnmx*nlnmx*mdimx*nclass,irot,nspin)
C name= ppbir
C dim 0 nlnmx*nlnmx*mdimx*nclass nlnmx*nlnmx*mdimx*nclass
C dim 1 irot irot
C dim 2 nspin nspin
C--------------------------
      call add_alloclist("ppbir",sizeof(ppbir),
     & int((nlnmx*nlnmx*mdimx*nclass)*(irot)*(nspin),kind=8))  !omplevel0
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
       call realtimediff(11,"ppbafp_v2")
       call realtimediff(13,'')
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
      iqxini=1 !for newaniso
      eibzmode = eibz4x0()

!! nov2016 moved from tetwt5 --> here
      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband )
      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
      mtet=(/1,1,1/)
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
      ! multitet=T ==> micro tetrahedron method (divided-tetrahedron). Not used so much now...
C--------------------------
C var 0 ekxx1(nband,nqbz)
C name= ekxx1
C dim 0 nband nband
C dim 1 nqbz nqbz
C--------------------------
      call add_alloclist("ekxx1",sizeof(ekxx1),
     & int((nband)*(nqbz),kind=8))  !omplevel0
C--------------------------
C var 1 ekxx2(nband,nqbz)
C name= ekxx2
C dim 0 nband nband
C dim 1 nqbz nqbz
C--------------------------
      call add_alloclist("ekxx2",sizeof(ekxx2),
     & int((nband)*(nqbz),kind=8))  !omplevel0
      allocate(ekxx1(nband,nqbz),ekxx2(nband,nqbz))

!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      if(eibzmode) then
!! commentout block inversion Use iqxendx=iqxend because of full inversion
         iqxendx=iqxend
C--------------------------
C var 0 nwgt(nqbz,iqxini:iqxendx)
C name= nwgt
C dim 0 nqbz nqbz
C dim 1 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("nwgt",sizeof(nwgt),
     & int((nqbz)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
C--------------------------
C var 1 igx(ngrp*2,nqbz,iqxini:iqxendx)
C name= igx
C dim 0 ngrp*2 ngrp*2
C dim 1 nqbz nqbz
C dim 2 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("igx",sizeof(igx),
     & int((ngrp*2)*(nqbz)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
C--------------------------
C var 2 igxt(ngrp*2,nqbz,iqxini:iqxendx)
C name= igxt
C dim 0 ngrp*2 ngrp*2
C dim 1 nqbz nqbz
C dim 2 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("igxt",sizeof(igxt),
     & int((ngrp*2)*(nqbz)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
C--------------------------
C var 3 eibzsym(ngrp,-1:1,iqxini:iqxendx)
C name= eibzsym
C dim 0 ngrp ngrp
C dim 1 -1:1 ((1)-(-1)+1)
C dim 2 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("eibzsym",sizeof(eibzsym),
     & int((ngrp)*(((1)-(-1)+1))*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
         allocate( nwgt(nqbz,iqxini:iqxendx), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     &        igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     &        eibzsym(ngrp,-1:1,iqxini:iqxendx))
         iprintx=.false.

         write(6,*)
         write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal() 
         if(MPI__root) iprintx=.true.
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

!All input. this returns required index stored in arrays in m_pbindex.
         call PBindex(natom,lx,l2nl,nx) 
                          ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
         call readqgcou() ! no input. Read QGcou and store date into variables.
!!  call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.
      else !dummy allocation to overlaid -check bound !sep2014
         iqxendx=iqxend
C--------------------------
C var 0 nwgt(1,iqxini:iqxendx)
C name= nwgt
C dim 0 1 1
C dim 1 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("nwgt",sizeof(nwgt),
     & int((1)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
C--------------------------
C var 1 igx(1,1,iqxini:iqxendx)
C name= igx
C dim 0 1 1
C dim 1 1 1
C dim 2 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("igx",sizeof(igx),
     & int((1)*(1)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
C--------------------------
C var 2 igxt(1,1,iqxini:iqxendx)
C name= igxt
C dim 0 1 1
C dim 1 1 1
C dim 2 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("igxt",sizeof(igxt),
     & int((1)*(1)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
C--------------------------
C var 3 eibzsym(1,1,iqxini:iqxendx)
C name= eibzsym
C dim 0 1 1
C dim 1 1 1
C dim 2 iqxini:iqxendx ((iqxendx)-(iqxini)+1)
C--------------------------
      call add_alloclist("eibzsym",sizeof(eibzsym),
     & int((1)*(1)*(((iqxendx)-(iqxini)+1)),kind=8))  !omplevel0
         allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &    ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      endif

C--------------------------
C var 0 llw(nw_i:nw,nq0i)
C name= llw
C dim 0 nw_i:nw ((nw)-(nw_i)+1)
C dim 1 nq0i nq0i
C--------------------------
      call add_alloclist("llw",sizeof(llw),
     & int((((nw)-(nw_i)+1))*(nq0i),kind=8))  !omplevel0
C--------------------------
C var 1 llwI(niw,nq0i)
C name= llwI
C dim 0 niw niw
C dim 1 nq0i nq0i
C--------------------------
      call add_alloclist("llwI",sizeof(llwI),
     & int((niw)*(nq0i),kind=8))  !omplevel0
      allocate( llw(nw_i:nw,nq0i), llwI(niw,nq0i) )
      llw=1d99
      llwI=1d99
      if(ixc==1011)then !ixc==11 is a debug mode to test contrib. at \Gamma point.
         goto 1191
      endif
!! for w4phonon. all nodes have wmu array.
      w4pmode=.false.
      if(sum(ixyz)/=0) w4pmode=.true.
      if(w4pmode) then
C--------------------------
C var 0 wmuk(2:nblochpmx,3)
C name= wmuk
C dim 0 2:nblochpmx ((nblochpmx)-(2)+1)
C dim 1 3 3
C--------------------------
      call add_alloclist("wmuk",sizeof(wmuk),
     & int((((nblochpmx)-(2)+1))*(3),kind=8))  !omplevel0
        allocate( wmuk(2:nblochpmx,3))
        wmuk=1d99
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
      write(6,'("irank=",i5," allocated(MPI__qtask)=",L5)')MPI__rank,allocated(MPI__qtask)
      do iq = iqxini,iqxend
        if(MPI__qtask(iq)) write(6,'("irank iq=",i5,i5)') MPI__rank,iq
      enddo

!! Get ngbq0 (for q=0) and broadcast for w4p
      if( MPI__root.and. w4pmode ) then
        q = (/0d0,0d0,0d0/)
        call readqg('QGcou', q, ginv,  quu,ngc,ngveccB)
        ngbq0 = nbloch+ngc
      endif
      call MPI__Broadcast(ngbq0)

       call realtimediff(13,"readqgcou")
       call realtimediff(15,'')
      do 1001 iq = iqxini,iqxend
        if( .not. MPI__Qtask(iq) ) cycle
        if (MPI__rootS) then
          ifrcwi = iopen('WVI.'//charnum5(iq),0,-1,mrecl)
          ifrcw  = iopen('WVR.'//charnum5(iq),0,-1,mrecl)
        endif
!!
        call cputid (0)
        q = qibze(:,iq)
        call readqg('QGcou', q, ginv,  quu,ngc,ngveccB) ! q was qq

!! Caution : confusing point
!!  ngc by QGcou is shown at the bottom of lqg4gw.
!!  ngc read from PPOVL are given by rdata4gw.
!!  Note that  ngc(iq>nqibz )=ngc (q=0), because when it is generated in mkqg.F
!!
c         if( newaniso2.and.iq==1 ) then ! *sanity check
        if( iq==1 ) then        ! *sanity check
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
C         if(allocated(zcousq)) deallocate( zcousq,vcousq )
        if(allocated(zcousq))  then
C--------------------------
C var 0 zcousq
C name= zcousq
C--------------------------
      call del_alloclist( "zcousq") !omplevel0
C--------------------------
C var 1 vcousq
C name= vcousq
C--------------------------
      call del_alloclist( "vcousq") !omplevel0
        deallocate( zcousq,vcousq )
        endif
C--------------------------
C var 0 zcousq(ngb0,ngb0)
C name= zcousq
C dim 0 ngb0 ngb0
C dim 1 ngb0 ngb0
C--------------------------
      call add_alloclist("zcousq",sizeof(zcousq),
     & int((ngb0)*(ngb0),kind=8))  !omplevel0
C--------------------------
C var 1 vcousq(ngb0)
C name= vcousq
C dim 0 ngb0 ngb0
C--------------------------
      call add_alloclist("vcousq",sizeof(vcousq),
     & int((ngb0),kind=8))  !omplevel0
        allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
        read(ifvcoud) vcousq
        read(ifvcoud) zcousq
        idummy=iclose(trim(vcoudfile))
        vcousq=sqrt(vcousq)

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
C         if(allocated(ppovlz)) deallocate(ppovlz)
        if(allocated(ppovlz))  then
C--------------------------
C var 0 ppovlz
C name= ppovlz
C--------------------------
      call del_alloclist( "ppovlz") !omplevel0
        deallocate(ppovlz)
        endif
C         if(allocated(ppovlzinv)) deallocate(ppovlzinv)
        if(allocated(ppovlzinv))  then
C--------------------------
C var 0 ppovlzinv
C name= ppovlzinv
C--------------------------
      call del_alloclist( "ppovlzinv") !omplevel0
        deallocate(ppovlzinv)
        endif
C         if(allocated(ppovl)) deallocate(ppovl)
        if(allocated(ppovl))  then
C--------------------------
C var 0 ppovl
C name= ppovl
C--------------------------
      call del_alloclist( "ppovl") !omplevel0
        deallocate(ppovl)
        endif
C--------------------------
C var 0 ppovl(ngc,ngc)
C name= ppovl
C dim 0 ngc ngc
C dim 1 ngc ngc
C--------------------------
      call add_alloclist("ppovl",sizeof(ppovl),
     & int((ngc)*(ngc),kind=8))  !omplevel0
C--------------------------
C var 1 ppovlz(ngb,ngb)
C name= ppovlz
C dim 0 ngb ngb
C dim 1 ngb ngb
C--------------------------
      call add_alloclist("ppovlz",sizeof(ppovlz),
     & int((ngb)*(ngb),kind=8))  !omplevel0
C--------------------------
C var 2 ppovlzinv(ngb,ngb)
C name= ppovlzinv
C dim 0 ngb ngb
C dim 1 ngb ngb
C--------------------------
      call add_alloclist("ppovlzinv",sizeof(ppovlzinv),
     & int((ngb)*(ngb),kind=8))  !omplevel0
        allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb),   ppovlzinv(ngb,ngb))
        call readppovl0(q,ngc,ppovl) !q was qq
c         ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
c         ppovlz(nbloch+1:nbloch+ngc,:) 
c     &     = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
C--------------------------
C var 0 ppovl_(ngb,ngb)
C name= ppovl_
C dim 0 ngb ngb
C dim 1 ngb ngb
C--------------------------
      call add_alloclist("ppovl_",sizeof(ppovl_),
     & int((ngb)*(ngb),kind=8))  !omplevel0
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
C--------------------------
C var 0 ppovl_
C name= ppovl_
C--------------------------
      call del_alloclist( "ppovl_") !omplevel0
C--------------------------
C var 1 ppovl
C name= ppovl
C--------------------------
      call del_alloclist( "ppovl") !omplevel0
        deallocate(ppovl_,ppovl)
        
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

C--------------------------
C var 0 rcxq(nmbas1,nmbas2,nwhis,npm)
C name= rcxq
C dim 0 nmbas1 nmbas1
C dim 1 nmbas2 nmbas2
C dim 2 nwhis nwhis
C dim 3 npm npm
C--------------------------
      call add_alloclist("rcxq",sizeof(rcxq),
     & int((nmbas1)*(nmbas2)*(nwhis)*(npm),kind=8))  !omplevel0
        allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
C--------------------------
C var 0 zw0(ngb,ngb)
C name= zw0
C dim 0 ngb ngb
C dim 1 ngb ngb
C--------------------------
      call add_alloclist("zw0",sizeof(zw0),
     & int((ngb)*(ngb),kind=8))  !omplevel0
        allocate( zw0(ngb,ngb) ) !, zxq (ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
        rcxq = 0d0

!! ---------------------------------------------------------------
!! === loop over spin=== =========================================
!! ---------------------------------------------------------------
       call realtimediff(17,'')
!         do 1003 is = 1,nspinmx
        do 1003 is = MPI__Ss,MPI__Se
          write(6,"(' ### ',2i4,' out of nqibz+n0qi+nq0iadd nsp=',2i4,' ### ')") 
     &     iq, is, nqibz + nq0i+nq0iadd, nspin
          if(debug) write(6,*)' niw nw=',niw,nw
          isf = is

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
          call gettetwt(q,iq,is,isf,nwgt(:,iq),frhis,nwhis,npm,
     i     qbas,ginv, ef, nqibz, nband,ekxx1,ekxx2, nctot,ecore,
     i     nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
     i     nbmx,ebmx,mtet,eibzmode) !nov2016

!! == x0kf_v4hz is the main routine to accumalte imaginary part of x0 ==
          iqeibz=iq
          if(npm==1) then
            ncc=0
          else
            ncc=nctot
          endif
          call x0kf_v4hz(npm,ncc,   
     i     ihw,nhw,jhw,whw,nhwtot, ! tetwt5
     i     n1b,n2b,nbnbx,nbnb,  ! use whw by tetwt5 ,
     i     q,  
     i     nspin,is,isf, !symmetrize, !
     i     qbas,ginv,  qbz,wbz, 
     d     nlmto,nqbz,nctot,    !noccx,noccxv,
     d     nbloch,  nwhis,      !nlnmx,mdimx,
     i     iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
     i     nqbze,nband,nqibz, 
     o     rcxq,                ! rcxq is the accumulating variable for spins 
     i     nolfco,zzr,nmbas_in, zcousq, !ppovl,nmbas2 is removed.,nmbas1 ppovlz, 
     i     chipm,eibzmode,      !z1offd,
     i     nwgt(:,iqeibz),igx(:,:,iqeibz),igxt(:,:,iqeibz),ngrp, eibzsym(:,:,iqeibz),crpa)
          write(6,*)' end of x0kf_v4h sum rcxq=',sum(abs(rcxq))
C--------------------------
C var 0 
C name= 
C--------------------------
      call del_alloclist( "") !omplevel0
          call tetdeallocate() !deallocate(ihw,nhw,jhw, whw,ibjb )
c          if(tetra)  deallocate( n1b,n2b)
 1003   continue;write(6,*) 'end of spin-loop nwp=',nwp !end of spin-loop
       call realtimediff(17,"Do1003")
c===========end of spin loop============================================

!! symmetrize and convert to Enu basis by dconjg(tranpsoce(zcousq)*rcxq8zcousq if eibzmode
       call realtimediff(19,'')
        if(eibzmode)  then
          is=1                  ! dummy
          call x0kf_v4hz_symmetrize(npm,!ncc,   
c     i     ihw,nhw,jhw,whw,nhwtot, ! tetwt5
c     i     n1b,n2b,nbnbx,nbnb,  ! use whw by tetwt5 ,
     i     q,  
     i     nspin,is,isf, !symmetrize, !
     i     qbas,ginv,  !qbz,wbz, 
c     i              nblocha,!nlnm,nlnmv,nlnmc,iclass, 
c     i              ppb(1,is), 
c     i              icore,ncore,
c     d     nlmto,nqbz,nctot,    !noccx,noccxv,
c     d              natom, !nl,nclass,natom,nnc,
     d     nbloch,  nwhis,      !nlnmx,mdimx,
     i     iq,ngb,ngc,ngpmx,ngcmx, !ngb/=ngc+nbloch for smbasis()=T oct2005
     i     nqbze,nband,nqibz, 
     o     rcxq,                ! rcxq is the accumulating variable for spins 
     i     nolfco,zzr,nmbas_in, zcousq, !ppovl,nmbas2 is removed.,nmbas1 ppovlz, 
     i     chipm,eibzmode,      !z1offd,
     i     ngrp, eibzsym(:,:,iqeibz))
        endif

!! reduction rcxq in the S-axis
        write(6,*) 'MPI__AllreduceSumS start'
        do i_reduction_npm=1,npm
          do i_reduction_nwhis=1,nwhis
            do i_reduction_nmbas2=1,nmbas2
              call MPI__AllreduceSumS(
     .         rcxq(1,i_reduction_nmbas2,i_reduction_nwhis,i_reduction_npm), nmbas1)
            enddo
          enddo
        enddo
        write(6,*) 'MPI__AllreduceSumS end'

!! Gaussian filtering of rcxq. Smearging Imag(X0). We may use egauss = 0.05 a.u.\sim 1eV for example.
        if(eginit) then
           call getkeyvalue("GWinput","GaussianFilterX0", egauss, default=0d0 )
           if(MPI__root) write(6,'("GaussianFilterX0= ",d13.6)') egauss
           if(abs(egauss)>1d-15) then
C--------------------------
C var 0 gfmat(nwhis,nwhis)
C name= gfmat
C dim 0 nwhis nwhis
C dim 1 nwhis nwhis
C--------------------------
      call add_alloclist("gfmat",sizeof(gfmat),
     & int((nwhis)*(nwhis),kind=8))  !omplevel0
              allocate(gfmat(nwhis,nwhis))
              call gaussianfilterhis(egauss,frhis,nwhis,gfmat)
           endif
C--------------------------
C var 0 rcxqin(1:nwhis)
C name= rcxqin
C dim 0 1:nwhis ((nwhis)-(1)+1)
C--------------------------
      call add_alloclist("rcxqin",sizeof(rcxqin),
     & int((((nwhis)-(1)+1)),kind=8))  !omplevel0
           allocate(rcxqin(1:nwhis))
           eginit=.false.
c$$$           do i=1,nwhis
c$$$           do j=1,nwhis
c$$$                 write(6,"(2i5,x,f10.8)")i,j,gfmat(i,j)
c$$$           enddo
c$$$           enddo
        endif
        if(abs(egauss)>1d-15) then
           do ipm=1,npm
           do imbas1=1,nmbas1
           do imbas2=1,nmbas2
              rcxqin = rcxq(imbas1,imbas2,1:nwhis,ipm)
              rcxq(imbas1,imbas2,1:nwhis,ipm) = matmul(gfmat,rcxqin)
           enddo
           enddo
           enddo
           write(6,"(' End of Gaussian Filter egauss=',f9.4)") egauss
        endif
        
       call realtimediff(19,"x0kf_sym")
       call realtimediff(21,'')
!! --- Hilbert transform.  Genrerate Real part from Imaginary part. ======
C         if(allocated(zxq) ) deallocate(zxq,zxqi)
        if(allocated(zxq) )  then
C--------------------------
C var 0 zxq
C name= zxq
C--------------------------
      call del_alloclist( "zxq") !omplevel0
C--------------------------
C var 1 zxqi
C name= zxqi
C--------------------------
      call del_alloclist( "zxqi") !omplevel0
        deallocate(zxq,zxqi)
        endif
C--------------------------
C var 0 zxq(nmbas1,nmbas2,nw_i:nw)
C name= zxq
C dim 0 nmbas1 nmbas1
C dim 1 nmbas2 nmbas2
C dim 2 nw_i:nw ((nw)-(nw_i)+1)
C--------------------------
      call add_alloclist("zxq",sizeof(zxq),
     & int((nmbas1)*(nmbas2)*(((nw)-(nw_i)+1)),kind=8))  !omplevel0
C--------------------------
C var 1 zxqi(nmbas1,nmbas2,niw)
C name= zxqi
C dim 0 nmbas1 nmbas1
C dim 1 nmbas2 nmbas2
C dim 2 niw niw
C--------------------------
      call add_alloclist("zxqi",sizeof(zxqi),
     & int((nmbas1)*(nmbas2)*(niw),kind=8))  !omplevel0
        allocate(zxq (nmbas1,nmbas2,nw_i:nw), zxqi(nmbas1,nmbas2,niw))
        write(6,'("goto dpsion5: nwhis nw_i niw nw_w nmbas1 nmbas2=",6i5)') nwhis,nw_i,nw,niw,nmbas1,nmbas2
        write(6,*)' -------- nmbas1,nmbas2=', nmbas1,nmbas2
        call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i   rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o   zxq, zxqi,
     i   chipm, schi,is,  ecut,ecuts)
        write(6,*)' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
C         if(allocated(rcxq) ) deallocate(rcxq)
        if(allocated(rcxq) )  then
C--------------------------
C var 0 rcxq
C name= rcxq
C--------------------------
      call del_alloclist( "rcxq") !omplevel0
        deallocate(rcxq)
        endif
       call realtimediff(21,"HilbertTransformation")

!! ===  RealOmega ===
        if (realomega) then
       call realtimediff(23,'')
          if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
          nwmax = nw
          nwmin = nw_i
!! prepare for iq0.
          iq0 = iq - nqibz
c               if(newaniso2) then
c$$$                 if( iq==1 ) then 
c$$$                    write(6,*)'open EPS0inv mpi=',MPI__rank 
c$$$                   ifepstinv = iopen('EPS0inv',0,-1,0)
c$$$                   write(ifepstinv) ngb
c$$$                 endif  
c$$$
c$$$          if(iqinit) then
c$$$            allocate( sk(ngb,nwmin:nwmax,nq0i),  sks(ngb,nwmin:nwmax,nq0i) )
c$$$            allocate( skI(ngb,niw,nq0i), sksI(ngb,niw,nq0i))
c$$$            iqinit=.false.
c$$$          endif
C--------------------------
C var 0 epstilde(ngb,ngb)
C name= epstilde
C dim 0 ngb ngb
C dim 1 ngb ngb
C--------------------------
      call add_alloclist("epstilde",sizeof(epstilde),
     & int((ngb)*(ngb),kind=8))  !omplevel0
          allocate(epstilde(ngb,ngb))
C--------------------------
C var 0 epstinv(ngb,ngb)
C name= epstinv
C dim 0 ngb ngb
C dim 1 ngb ngb
C--------------------------
      call add_alloclist("epstinv",sizeof(epstinv),
     & int((ngb)*(ngb),kind=8))  !omplevel0
          allocate(epstinv(ngb,ngb))
c               endif
!KINO               write(6,*)'kino: nwmin,nwmax,ngb=',nwmin,nwmax,ngb
          write(6, *)" === trace check for W-V === nwmin nwmax=",nwmin,nwmax
       call realtimediff(23,"ralloc")
       call realtimediff(25,'')
          do 1015 iw  = nwmin,nwmax
            frr= dsign(freq_r(abs(iw)),dble(iw))
            imode = 1
c                  if(newaniso2.and.iq<=nqibz) then !for mmmw
            if(iq<=nqibz) then  !for mmmw
              if(iq==1) then
                ix=1
                zw0(:,1)=0d0
                zw0(1,:)=0d0
              else
                ix=0
              endif

!!  Eqs.(37),(38) in PRB81 125102 (Friedlich)
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))

!! w4p writing eps
              if(iw==0.and.w4pmode) then 
                !static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion.
                             !(ix=1 is set for q=0) 
                ifw4p = ifile_handle()
                open(ifw4p,file='W4PHONON.'//charnum5(iq),form='unformatted')
                write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
                write(ifw4p) epstinv(ix+1:ngb,ix+1:ngb) 
                close(ifw4p)
              endif  
       call realtimediff(27,'')

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
       call realtimediff(27,"zweqzw0")
       call realtimediff(29,'')
              if (MPI__rootS)then
                write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
              endif
              call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
       call realtimediff(29,"tr_chkwrite")
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
                if(iq0<=nq0i) llw(iw,iq0)= 1d0/epstinv(1,1)
c !TIME1 "end of matcinv_epstinv" !this gives wrong message, probably 
c         because of a bug of !TIME1 processing for MPI mode.

!! Wing elements calculation july2016
!! We need check nqb is the same as that of q=0
                if(ixyz(iq0)/=0.and.iw==0) then
                   if(ngb/=ngbq0) then
                      write(6,*)q,iq0,ngb,ngbq0
                      call rx('hx0p0_sc: ngb/=ngbq0')
                   endif
                   wmuk(2:ngb,ixyz(iq0))=epstinv(1,2:ngb)/epstinv(1,1) !this is dot(q(:)*w_mu(:,igb)). See PRB125102(2016) eq.(36)
                endif
              else
c commentout block inversion
c$$$                     sk  (1:ngb,iw,iq0)= zxq(1,1:ngb,iw)
c$$$                     sks (1:ngb,iw,iq0)= zxq(1:ngb,1,iw)
c$$$c                     sks (1:ngb,iw,iq0)= zxq(2,1:ngb,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
                if(iq0<=nq0i) llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw) 
              endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmmmm 
c     write(6,"('mmmw0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c     write(6,"('mmmw0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c     write(6,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,fourpi/sum(q**2*tpioa**2)/llw(iw,iq0),
c     &                  w_k(2:10:3)/llw(iw,iq0),w_k(63:70:3)/llw(iw,iq0)
c     write(6,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0),
c     &                  w_ks(2:10:3)/llw(iw,iq0),w_ks(63:70:3)/llw(iw,iq0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       call realtimediff(31,'')
              if(iq0<=nq0i) write(6,"('iq iw_R omg(iw) eps(wFC) eps(woLFC) ',2i5,x,10(d13.6,2x,d13.6,x,d13.6,2x,d13.6,x,d13.6))")
     &         iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)

       call realtimediff(31,"writeiqiwreal")
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
 1015     continue              !iw
       call realtimediff(25,"do1015")

c     if(newaniso2) then
c     if(allocated(sk)) deallocate(sk,sks,w_k,w_ks)
c     endif
C           if( allocated(zzr) ) deallocate(zzr)
          if( allocated(zzr) )  then
C--------------------------
C var 0 zzr
C name= zzr
C--------------------------
      call del_alloclist( "zzr") !omplevel0
        deallocate(zzr)
        endif
        endif 
!! === RealOmega end ===

!! === ImagOmega ===
       call realtimediff(33,'')
        if (imagomega) then
          write(6,*)' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
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
              if (MPI__rootS) then
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
                if(iq0<=nq0i) llwI(iw,iq0)= 1d0/epstinv(1,1)
              else
c commentout block inversion
c$$$                     skI  (1:ngb,iw,iq0)= zxqi(1,1:ngb,iw)
c$$$c                     sksI (1:ngb,iw,iq0)= zxqi(2,1:ngb,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     sksI (1:ngb,iw,iq0)= zxqi(1:ngb,1,iw) !nmbas1=2 see z1stcol in x0kf_v4h.
c$$$                     vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
c$$$                     vcou1sq= sqrt(vcou1)
c$$$!! llwI without LFC. LFC contribution is added in 
                if(iq0<=nq0i) llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
              endif  
              if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)")
     &         iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
            endif

 1016     continue
c               if(newaniso2) then
c$$$                 if(iq==1) ifepstinv = iclose('EPS0inv') !iq==1 close write mode.
C--------------------------
C var 0 epstinv
C name= epstinv
C--------------------------
      call del_alloclist( "epstinv") !omplevel0
          deallocate(epstinv)
C           if(allocated(epstilde)) deallocate(epstilde)
          if(allocated(epstilde))  then
C--------------------------
C var 0 epstilde
C name= epstilde
C--------------------------
      call del_alloclist( "epstilde") !omplevel0
        deallocate(epstilde)
        endif
c               endif
        endif 
!! === ImagOmega end ===
       call realtimediff(33,"imagomega")

c     1002 continue  ! end of frequency block-loop
C         if(allocated(vcoul)) deallocate(vcoul)
        if(allocated(vcoul))  then
C--------------------------
C var 0 vcoul
C name= vcoul
C--------------------------
      call del_alloclist( "vcoul") !omplevel0
        deallocate(vcoul)
        endif
C         if(allocated(zw0)) deallocate(zw0)
        if(allocated(zw0))  then
C--------------------------
C var 0 zw0
C name= zw0
C--------------------------
      call del_alloclist( "zw0") !omplevel0
        deallocate(zw0)
        endif
C         if(allocated(zxq )) deallocate(zxq)
        if(allocated(zxq ))  then
C--------------------------
C var 0 zxq
C name= zxq
C--------------------------
      call del_alloclist( "zxq") !omplevel0
        deallocate(zxq)
        endif
C         if(allocated(zxqi)) deallocate(zxqi)
        if(allocated(zxqi))  then
C--------------------------
C var 0 zxqi
C name= zxqi
C--------------------------
      call del_alloclist( "zxqi") !omplevel0
        deallocate(zxqi)
        endif

        if (MPI__rootS) then
          ifrcwi = iclose('WVI.'//charnum5(iq))
          ifrcw  = iclose('WVR.'//charnum5(iq))
        endif
!!  
 1001 continue
       call realtimediff(15,"do1001")
c============end of loop over q point =================================
c=======================================================================
      call MPI__barrier()

       call realtimediff(35,'')
!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
      if(MPI__size/=1) then
        do iq=nqibz+1,iqxend
          iq0 = iq - nqibz
c     write(6,*)' iq iq0 mpi_rank mpi_ranktab(iq)=',iq, iq0,MPI__rank,MPI__ranktab(iq),MPI__root,nw,nw_i,niw
          if(MPI__Qranktab(iq)/=0) then !jan2012
            if(MPI__Qranktab(iq) == MPI__rankQ) then
c     write(6,*)' mpi_send iq from',iq,MPI__ranktab(iq)
c     write(6,*)' send llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
c     do i=nw_i,nw
c     write(6,*)'sendxxx',i,llw(i,iq0)
c     enddo   
c     write(6,*)' send llwI sum=',sum(abs(llwI(:,iq0))),niw
              dest=0
              if(iq0<=nq0i) then
                 call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
                 call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
              endif
              if(ixyz(iq0)/=0) then
                call MPI__DbleCOMPLEXsendQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,dest)
              endif
            elseif(MPI__rootQ) then
c     write(6,*)' mpi_recv iq from',iq,MPI__ranktab(iq),nw,nw_i,niw
              src=MPI__Qranktab(iq)
              if(iq0<=nq0i) then
                 call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
                 call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
              endif
              if(ixyz(iq0)/=0) then
                call MPI__DbleCOMPLEXrecvQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,src)
              endif
c     do i=nw_i,nw
c     write(6,*)'recivxxx',i,llw(i,iq0)
c     enddo   
c     write(6,*)' recv llw sum=',sum(abs(llw(:,iq0))),nw,nw_i
c     write(6,*)' recv llwI sum=',sum(abs(llwI(:,iq0))),niw
            endif
          endif
        enddo  
      endif
       call realtimediff(35,"w0mpi")

c commentout block inversion
c$$$!! Add LFC (local field correction) to llw and llwI
c$$$         if(newaniso2 .and. MPI__rank == 0 ) then ! only on root node
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
 1191 continue
       call realtimediff(37,'')
c         if(newaniso2 .and. MPI__rank == 0 ) then ! MIZUHO-IR only on root node
      if(MPI__rank == 0 ) then  ! MIZUHO-IR only on root node


!! ix=1011 is a special mode to overwrite llw and llwI for test purpose
!! A file W0W0I is  generated by call w0w0i, but usually unused at anywhere.
        if(ixc==1011) then  
          ifw0w0i = ifile_handle('W0W0I')
          open(ifw0w0i,form='unformatted')
          read(ifw0w0i) nw_ixx,nwxx,niw,nq0ix
          write(6,*)'w0w0i: n=',nw_ixx,nwxx,niw,nq0ix
          if(nq0i/=nq0ix)  call rx('nq0i/=nq0ix')
          if(nw_i/=nw_ixx) call rx(nw_i/=nw_ixx)
          if(nw/=nwxx) call rx(nw/=nwxx)
          read(ifw0w0i) llw(nw_i:nw,1:nq0i)
          read(ifw0w0i) llwI(1:niw,1:nq0i)
c            read(ifw0w0i) w0(nw_i:nw)
c            read(ifw0w0i) w0i(1:niw)
          close(ifw0w0i)
        endif  

!! get w0 and w0i (diagonal element at Gamma point)
!! This return w0 and w0i. (llw and llwi are input)
!! Outputs w0,w0i,llmat. See use m_w0w0i at the begining of this routine.
        call w0w0i(llw,llwI,nw_i,nw,nq0i,niw,q0i)  !all inputs. get effective W0,W0i, and L(omega=0) matrix.

!! Finalize w4phonon
!! wmuk(ix)= matmul(wmu,qv) ==> wmu= matmul(wmuk,qvinv)
        if(w4pmode) then
          do i=1,3
            qv(:,i)= tpioa*q0i(:,ixyz(i))
            qv(:,i)= qv(:,i)/sqrt(sum(qv(:,i)**2))
          enddo
          call matinv(3,qv)
C--------------------------
C var 0 wmu(2:ngbq0,3)
C name= wmu
C dim 0 2:ngbq0 ((ngbq0)-(2)+1)
C dim 1 3 3
C--------------------------
      call add_alloclist("wmu",sizeof(wmu),
     & int((((ngbq0)-(2)+1))*(3),kind=8))  !omplevel0
          allocate( wmu(2:ngbq0,3) )
          do igb=2,ngbq0
            wmu(igb,:) =matmul(wmuk(igb,:),qv)
          enddo
          ifw4p = ifile_handle()
          open(ifw4p,file='W4PHONON.HeadWing',form='unformatted')
          write(ifw4p) llmat(1:3,1:3),ngbq0 !for q~0
          write(ifw4p) wmu(2:ngbq0,1:3) !for q~0
          close(ifw4p)
ccccccccccccccccccccc
c          write(6,*)'nqbq0=',ngbq0
c          write(6,*)'llmat=',llmat
c          write(6,*)'wmu sum=',sum(abs(wmu(2:ngbq0,1:3)))
c          write(6,*)'wmuksum=',sum(abs(wmuk(2:ngbq0,1:3)))
c          call rx(' test end xxxxxxxxxxxx')
cccccccccccccccccccccc
C--------------------------
C var 0 wmu
C name= wmu
C--------------------------
      call del_alloclist( "wmu") !omplevel0
C--------------------------
C var 1 wmuk
C name= wmuk
C--------------------------
      call del_alloclist( "wmuk") !omplevel0
          deallocate(wmu,wmuk)
        endif

!! Read WVR and WVI at Gamma point, and give correct W(0) (averaged in the Gamma cell, where
!! Gamma cell) is the micro cell of BZ including Gamma point).

c     write(6,*)'sumcheck w0,w0i=',sum(abs(w0)),sum(abs(w0i))
!! === w0,w0i are stored to zw for q=0 ===
!! === w_ks*wk are stored to zw for iq >nqibz ===
! We assume iq=1 is for rank=0
        do iq = 1,1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
c               if( .not. MPI__task(iq) ) cycle
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
      endif
      is = iclose('hbe.d')
       call realtimediff(37,"WVRI")
       call realtimediff(1,"ProgAll")
       call print_realtimediff()
      call cputid(0)
      write(6,*) '--- end of hx0fp0_sc --- irank=',MPI__rank
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
      character*(*)::tagname
      trwv=0d0
      do i = 1,nbloch
        trwv = trwv + zw(i,i)
      enddo
      trwv2 = 0d0
      do i = 1,ngb
         trwv2 = trwv2 + zw(i,i)
      enddo                     !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
c     do i = 1,ngb
c     write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c     enddo
      end


