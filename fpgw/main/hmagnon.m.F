!! cp hx0fp0.m.F in order to develop by Okumura (2017/03/23)
!!  Calculate x0, \epsilon, spin susceptibility. 
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_lmf_chipm mode works).
      program hmagnon
c      use matm,only: matm,matcinv
      use m_readwan,only: write_qdata, wan_readeigen, wan_readeval,
     &     readscr, checkorb, getqsym, checkorb2
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
     & MPI__hmagnon_rankdivider,
     & MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size,
     & MPI__ranktab,MPI__consoleout,MPI__barrier,MPI__MEq,MPI__REAL8send,MPI__REAL8recv,MPI__AllreduceSum
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
     &  zw(:,:),
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
      complex(8),allocatable:: UU(:,:),VT(:,:),work(:),ddd(:,:)
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
      integer:: istat,imb,imb1,imb2,nmbas_in
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
      integer :: ngb0,ifvcoud,idummy,ifepstinv,igb1,igb2,ngb_in,iq0,ifisk,iqx,ig,nmbas1x,ifiss,iq0x
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),zcousqrsum(:,:,:),zcousqr(:,:)
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
      logical:: lqall
      integer,allocatable :: iclasst(:), invgx(:)
      integer:: ificlass,ifile_handle,k
      complex(8),allocatable:: ppovl_(:,:)

      logical:: readw0w0itest=.false.

      real(8)::ebmx
      integer:: nbmx,mtet(3),ifq0p
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)

!     okumura
      integer::nwf,nsp_w,nqtt_w,iwf,jwf,inwf,iexc,kwf,lwf,ijwf,klwf,nnwf,ijwf_j
      integer::ifhamw,ifgas,ifchipmz_wan,ifchipmr_wan,ifchipmrk_wan
      real(8),allocatable::ev_w1(:,:),ev_w2(:,:)
      complex(8),allocatable::evc_w1(:,:,:),evc_w2(:,:,:)
      logical::wan=.true.,egasmode
      logical(8):: ijklmag !for checkorb2
      complex(8),allocatable::kmat(:,:,:,:),wanmat(:,:)
     &     ,wkmat(:,:),wkmat2(:,:),rmat(:,:,:),rmat2(:,:)
     &     ,swkwmat(:,:),swkwmat2(:,:)
     &     ,sqemat(:,:),plmat(:,:),plpmat(:,:),plpmat_inv(:,:)
     &     ,wmat_check(:,:)
      complex(8)::trmat,trmatt
      real(8),allocatable::imat(:,:) !unit matrix for 1-WK
      real(8):: wan_ecore(1),rewanijkl
      logical::npmtwo,diag=.false.,t2g
      complex(8)::wan_i,wan_j,wan_k,wan_l,wanijkl
      integer::igv,bsimple=3
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
      logical(8)::init2=.true.,llsym,gskip
      logical(8)::d100,d110,d111,d1xx,dhpb,dxwf,dhnb,dnpb
      real(8)::rlatp(3,3),xmx2(3),qqin(3),qshort(3),qshort2
      integer:: nlatout(3,48),nout,iout,nqsym
      integer,allocatable:: iqlist(:)
!!! WK eigenvalue check
      integer::ifwkeigen,ifkeigen,ifweigen,ifvecr,ifwkeigen2,ifwkeigen3
!!! temp variable
      logical(8)::lsmo
      integer::ifwanmat
!!! for MAX(Im[R])
      integer::imaximr=0
      real(8)::maximr,w_maximr
      complex(8)::sumrpa_maximr(1),summf_maximr(1)
      real(8),allocatable:: rpa_maximr(:),mf_maximr(:)

!!! for debug (2018/05/23)
!      real(8)::sumkm1,sumkm2,sumkm3,sumkm4,sumkm5,znorm,rnqbz,cellvol

!! -------------------------------------------------------------------
      call MPI__Initialize()
      call MPI__consoleout('hmagnon')
      call cputid (0)
      hartree  = 2d0*rydberg()
      pi       = 4d0*datan(1d0)
      fourpi   = 4d0*pi
      sqfourpi = sqrt(fourpi)

      call cputid(0)

      lqall=.true.
      write(6,*) " OK ixc=223    chipm_wannier sergey's "
      imagomega =.false.
      omitqbz=.false.
      epsmode = .true.
      chipm=.true.
      nolfco=.true.
      wan=.true.

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
      
c$$$!! check write
      if(MPI__root) then
         do i=1,nqbz
            if(i<10.or.i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
            if(i==10.and.nqbz>18) write(6,"('... ')")
         enddo
         write(6,*)' !!nqbz nqibz =',nqbz,nqibz
      endif

!! EFERMI
cc need for gettetwt; automatically read file(okumura) OK
      call readefermi()


!! === Readin by genallcf. Set basic data for crystal
!! See "use m_genallcf_v3" at the begining of this routine
!!
      incwfin=0  !use ForX0 for core in GWIN
      call genallcf_v3(incwfin) !in module m_genallcf_v3
!     ! incwfin is for core handling
c     write(6,*)' ncore nw_input delta=', ncore,nw_input,delta
      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
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

c     Weight for irreducible q-point (qibz)
c      do iq=1,nqibz
c         write(6,"('wibz',4f9.4)") wibz(iq),qibz(:,iq)
c      enddo
c

cccccccccccccccccc Screening W for magnon 
      nnwf=nwf*nwf
      allocate(scrw(nnwf,nnwf))
      call readscr(nwf,scrw)
cccc
      do iexc=1,2
         do iwf=1,nwf
            ijwf=(1+nwf)*iwf-nwf

            if (iexc==1) then
               if (iwf==1) write (6,*) "--- W  (direct) -------------------"
               do jwf=1,nwf
                  if (iwf==1) write (6,"(A6)",advance="no") "     "
                  if (iwf==1) write (6,"(i3)",advance="no") jwf
               enddo
               if (iwf==1) write(6,*) ""

               write (6,"(i3)",advance="no") iwf
               do jwf=1,nwf
                  ijwf_j=(1+nwf)*jwf-nwf
                  write (6,"(f9.4)",advance="no") abs(scrw(ijwf,ijwf_j))
               enddo
            else
               if (iwf==1) write (6,*) "--- W' (exchange) -----------------"
               do jwf=1,nwf
                  if (iwf==1) write (6,"(A4)",advance="no") "    "
                  if (iwf==1) write (6,"(i4)",advance="no") jwf
               enddo
               if (iwf==1) write(6,*) ""

               write (6,"(i3)",advance="no") iwf
               do jwf=1,nwf
                  ijwf_j=(jwf-1)*nwf+iwf
                  write (6,"(f9.4)",advance="no") abs(scrw(ijwf_j,ijwf_j))
               enddo
            endif
            write(6,*) ""
         enddo
      enddo
      scrw(:,:)=scrw(:,:)/hartree
cccccccccccccccccc
      iqxend = nqibz !+ nq0i

!! Initialization of readEigen !readin m_hamindex
      ginv = transpose(plat)
      call minv33tp(plat,qlat)
      if(verbose()>50) print *,'eeee exit of init_readeigen2'

      do iq=1,nqibz
         iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
c     nstibz(iq) = nstbz(iqbz)
         write(6,"(' iq qibz nstibz=',2i5,3f9.4,i5)")iq,iqbz,qibz(:,iq) !,nstibz(iq)
      enddo


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

cccccccccccccccccccccccccccccccccc
!! Tetra hedron initialization
      if (egasmode) then
         noccxv = 1
      elseif (wan) then
         noccxv = nwf
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
      nspinmx = nspin
      iqxini=1
      mtet=(/1,1,1/) !dummy
c$$$      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
!okumura
      
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
 !! EIBZ mode memo for nolfco (right?) 
 !! If eibzmode=T, it is efficient but can slightly break crystal symmetry.(how much?)
 !! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
 !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
c     okumura
c      eibzmode=eibz4x0()
      eibzmode=.false.          !simple symmetry
c$$$c     end okumura
      iqxendx=iqxend
      allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &     ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      nwgt=1

!! Calculate x0(q,iw) and W == main loop 1001 for iq. 
 !! NOTE: iq=1 (q=0,0,0) write 'EPS0inv', which is used for iq>nqibz for ixc=11 mode
 !! Thus it is necessary to do iq=1 in advance to performom iq >nqibz. 
 !! (or need to modify do 1001 loop).
 !! iq>nqibz for ixc=11 is not time-consuming (right???)
 !    call MPI__hx0fp0_rankdivider2(iqxini,iqxend)
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do is=1,nspinmx
         call wan_readeigen(qbz(:,:),nqbz,is,nsp_w)
         call write_qdata(ginv,nqbz,qbz(:,:))
      enddo      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !!! For MPI test (intended symmetric iq)
      nqsym=0
      do iq=1,nqibz
         q = qibz(:,iq)         !qibze ! you can spefify any q, which you like.
         call getqsym(qlat,ngrp,iq,symgg,q,qshort,llsym)
         if (llsym) nqsym=nqsym+1
      enddo
      write(6,*) "nqsym:",nqsym
      call MPI__hmagnon_rankdivider(nqsym)

      write(6,*) "nqsym mpi_MEq(:)",nqsym,mpi__MEq
      allocate(rpa_maximr(mpi__MEq),mf_maximr(mpi__MEq)) !list of w(MAX(Im[R])): magnon peak
      rpa_maximr=0d0; mf_maximr=0d0

      !iqlist for syml
      !qibz(:,iqlist(iq)): intended q (iq=1,nqsym)
      allocate(iqlist(1:nqsym))
      nqsym=0
      do iq=1,nqibz
         q = qibz(:,iq)         !qibze ! you can spefify any q, which you like.
         call getqsym(qlat,ngrp,iq,symgg,q,qshort,llsym)
         if (llsym) then
            nqsym=nqsym+1
            iqlist(nqsym)=iq
         endif
      enddo
c$$$      do iq=1,nqsym
c$$$         write(6,*) "iqlist",iqlist(iq),qibz(:,iqlist(iq))
c$$$      enddo

      !!! for MPI (finally nqsym ---> nqibz)
      do iq=1,nqsym
         if (MPI__task(iq)) write(6,'("iq,MPI_rank, mpi__ranktab",3I8)') iq,MPI__rank,mpi__ranktab(iq)
!         if (MPI__task(iq)) write(6,'("iq,qshort,I8,3f13.5")') iq,qibz(:,iqlist(iq))
      enddo
      call MPI__barrier()

!      call rx("aa")
c
!! ======== Loop over iq ================================
      do 1001 iq = 1,nqsym      ! NOTE: q=(0,0,0) is omitted when iqxini=2
         if(.not. MPI__task(iq) .and. iq .ne. 1) cycle
         q=qibz(:,iqlist(iq))
         if (iq == 1) then
            qshort=0d0; llsym=.true.
         else
            call getqsym(qlat,ngrp,iq,symgg,q,qshort,llsym) !!! exclude Gamma
            imaximr=imaximr+1 !!! imaximr=1,2,...,mpi__MEq
         endif
         if (llsym) write(6,"('qshort (on syml) = ',3f9.5)") qshort            

!!!   |q+G|
         ngc=1
         if (egasmode) then
            if(.not.allocated(ngveccB)) allocate(ngveccB(3,1))
            ngveccB=0        
            nwf=1
            nnwf=1
         endif

!! Caution : confusing point
!!  ngc by QGcou is shown at the bottom of lqg4gw.
!!  ngc read from PPOVL are given by rdata4gw---> ngc(iq>nqibz )=ngc for q=0
         if( iq==1 ) then       ! *sanity check
            if(sum(q**2)>1d-10) then
               call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
            endif
         endif

!!! write qshort instead of q
         write(6,"('===== do 1001: iq q(qshort)=',i7,3f9.4,' ========')")iqlist(iq),qshort !qq
!     ! rcxq: imaginary part after x0kf_v4h and symmetrization. 
!     ! zxq and zxqi are the main output after Hilbert transformation
         if (egasmode) allocate( rcxq(nnwf,nnwf,nwhis,npm) )
         allocate( zxq (nnwf,nnwf,nw_i:nw), zxqi (nnwf,nnwf,niw))
         zxq=0d0;  zxqi=0d0;  rcxq = 0d0
         if(debug) write(6,*)' niw nw=',niw,nw
         
!     ! ==== spin chi_charge or chi_+- ====
         is=1
         isf=2
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
c$$$  !!! avoid n1b=0 or n2b=0 check
c$$$  evc_w1(1:nwf,0,1:nqbz)=0d0
c$$$  evc_w2(1:nwf,0,1:nqbz)=0d0
!!! simple band for check ---> bsimple = 1:simple, 2:gas, 3:else
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
c$$$  ifgas=ifile_handle()
c$$$  if (is==1 .and. iq==54) open(ifgas,file="band_electrongas.dat")                
                do kx=1,nqbz
                   call read_qgband(alat,plat,qbz(:,kx),ngc,ngveccB,is,qgsq1(1:ngc,kx),qsh1)
                   call read_qgband(alat,plat,q+qbz(:,kx),ngc,ngveccB,isf,qgsq2(1:ngc,kx),qsh2)
c     write(6,"('qqq111=',3f9.4,2x,3f9.4)") qbz(:,kx), qsh1
c     write(6,"('qqq222=',3f9.4,2x,3f9.4)")q+qbz(:,kx),qsh2

!!! replace eigenvalue , egasmode nwf=1
                   ev_w1(1:nwf,kx)=qgsq1(:,kx) !Rydberg
                   ev_w2(1:nwf,kx)=qgsq2(:,kx) !Rydberg

                enddo
                if (allocated(qgsq1)) deallocate(qgsq1)
                if (allocated(qgsq2)) deallocate(qgsq2)

c$$$  ev_w1(:,kx)=ev_w1(:,kx)+2d0
c$$$  ev_w2(:,kx)=ev_w2(:,kx)-2d0                

             else               !! bsimple==3
                iwf=ifile_handle()
                print *,"wan_readeval start"
                
                do kx=1,nqbz    !!! ev_w1, ev_w2 unit: [Ry]
                   call wan_readeval(  qbz(:,kx), is,
     o                  ev_w1(1:nwf,kx), evc_w1(1:nwf,1:nwf,kx))
                   call wan_readeval(q+qbz(:,kx), isf,
     o                  ev_w2(1:nwf,kx), evc_w2(1:nwf,1:nwf,kx))

c     compare Friedrich ---> symmetric wannier function (evc_w2=evc_w1) (2018/05/22)
c      evc_w2=evc_w1

c$$$                   if(iq==1 .and. kx==1) open(iwf,file="waneval_check.data")
c$$$!     if (kx==1 .and. iq > 100 ) open(iwf,file="wan_eval_check.data")
c$$$                   if (qbz(2,kx)==1.0) then
c$$$                      if (qbz(3,kx)==1.0) then
c$$$!     if (iq > 100) then
c$$$                         write(iwf,"('q(1) ev_w1(:,kx)',10f9.4)") qbz(1,kx),ev_w1(:,kx)
c$$$                      endif
c$$$                   endif
c$$$                   if (iq==1 .and. kx==nqbz) close(iwf)

c$$$  !! band shift for normalization check
c$$$                   ev_w1(:,kx)=ev_w1(:,kx)+2d0
c$$$                   ev_w2(:,kx)=ev_w2(:,kx)-2d0
                   
!!!   for wannier eigenvector (evc_w1, evc_w2) for M(:)


cccccccc checkorb iwf ---> lorb(1:s, 2:p, 3:d, 4:f, 5:g) 
c$$$                   do iwf=1,nwf
c$$$                      call checkorb(iwf,nwf,lorb)
c$$$                      print *,"checkorb:: ",lorb
c$$$                   enddo
cc                   call rx("test end")
cccccccc
                   
c$$$                   do iwf=1,nwf
c$$$                      do jwf=1,nwf
c$$$                         write(6,"('sum evc1,evc2',2I4,4f9.4)") iwf,jwf,
c$$$     &                        abs(sum(evc_w1(:,iwf,kx)*dconjg((evc_w1(:,jwf,kx))))),
c$$$     &                        abs(sum(evc_w2(:,iwf,kx)*dconjg((evc_w2(:,jwf,kx))))),
c$$$     &                        abs(sum(evc_w1(:,iwf,kx)*dconjg((evc_w2(:,jwf,kx))))),
c$$$     &                        abs(sum(evc_w2(:,iwf,kx)*dconjg((evc_w1(:,jwf,kx)))))
c$$$c$$$!!! band sum
c$$$c$$$     &                        abs(sum(evc_w1(iwf,:,kx)*dconjg((evc_w1(jwf,:,kx))))),
c$$$c$$$     &                        abs(sum(evc_w2(iwf,:,kx)*dconjg((evc_w2(jwf,:,kx))))),
c$$$c$$$     &                        abs(sum(evc_w1(iwf,:,kx)*dconjg((evc_w2(jwf,:,kx))))),
c$$$c$$$     &                        abs(sum(evc_w2(iwf,:,kx)*dconjg((evc_w1(jwf,:,kx)))))
c$$$
c$$$
c$$$c     &                        abs(sum((evc_w2(:,iwf,kx))*dconjg((evc_w2(:,jwf,kx)))))
c$$$c     &                        sum(abs(evc_w2(:,iwf,kx)))
c$$$                      enddo
c$$$                   enddo

c$$$  endif
                enddo
             endif
             
c$$$             !!!!!!!
c$$$             do iwf=1,nwf
c$$$                write(6,"('max[ev_w1(',I1,',:)]',f13.5) =") iwf,minval(ev_w1(iwf,:))
c$$$                write(6,"('min[ev_w2(',I1,',:)]',f13.5) =") iwf,maxval(ev_w2(iwf,:))
c$$$             enddo
c$$$             !!!!!!!
             
             write(6,"('= start wan_gettetwt =',2i6,3f9.4)") nwf,iq,q
!             wan_ecore=0.0;nbmx=nwf

!!!!  shift band (2017, Sep, 13)
c$$$             if (.true.) then !!! ev_w1 < ev_w2
c$$$                ev_w1=ev_w1-(maxval(ev_w1)-minval(ev_w2))-0.001
c$$$                ev_w2=ev_w2+(maxval(ev_w1)-minval(ev_w2))+0.001
c$$$                write(6,*) "minvalmaxval:",maxval(ev_w1),minval(ev_w2)
c$$$             endif

c$$$             if (.true.) then !!! ev_w1 < ef and ev_w1 > ef
c$$$                ev_w1=ev_w1-(maxval(ev_w1)-ef)-0.0001
c$$$                ev_w2=ev_w2-(minval(ev_w2)-ef)+0.0001
c$$$                write(6,*) "minvalmaxval:",maxval(ev_w1),minval(ev_w2)
c$$$             endif

ccccccccccc tetrahedron part             
             call gettetwt(q,iq,isdummy,isdummy,nwgt(:,iq),frhis,nwhis,npm,
     i            qbas,ginv, ef, nqibz, nwf,ev_w1,ev_w2, nctot,wan_ecore,
     i            nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
     i            nwf,ebmx,mtet,eibzmode) !nov2016
             write(6,*) "=== end gettetwt"
             deallocate(ev_w1,ev_w2)
ccccccccc normalization of Im[K]:
             rnqbz=1/real(nqbz)
             znorm=-1d0*pi
!             znorm=sum(nbnb(:,:))*pi*rnqbz
           endif !end if(wan)
ccccccccc generate wanmat         
           if (.not. allocated(kmat)) allocate(kmat(1:nnwf,1:nnwf,1:nwhis,1:npm),
     &          wanmat(1:nnwf,1:nnwf)) !! (nwf^2*nwf^2)
           kmat=(0d0,1d-15)
!           kmat=(0d0,1d-10) !orbtest
           do 2011 kx=1,nqbz    !! discrete k-point loop
           do 2012 jpm=1,npm !! jpm=2: negative frequency
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
              
c==========================caliculate M(:,:) : wanmat(nwf,nwf)
ccc   matrix elements from wannier eigenvector

c$$$                diag=.false.
c$$$                if (diag) then
c$$$                   do iwf=1,nwf
c$$$                      wanmat(iwf,iwf)=(1d-1,0d0)
c$$$                   enddo
c$$$                endif
c                else
                   
c$$$  !!! evc_w1: upspin ; evc_w2: downspin
              it=n1b(ibib,kx,jpm) !index for n  for q
              itp=n2b(ibib,kx,jpm) !index for n' for q+k
              wanmat=(0d0,0d0)
              ijwf=0

              do 3002 iwf=1,nwf
              do 3003 jwf=1,nwf
                 ijwf=ijwf+1
                 klwf=0
              do 3004 kwf=1,nwf
              do 3005 lwf=1,nwf
                 klwf=klwf+1               
c                   call checkorb2(iwf,jwf,kwf,lwf,ijklmag)
c                   ijklmag=.true. !test for FeFe
c                   print *,"checkorb2:",iwf,jwf,kwf,lwf,ijklmag !! ijklmag=T: {ijkl} in same atom
!     if (ijklmag) then !!! same atom
                 if (.true.) then !!! same atom
                   !!! calculate numerator of Kmatrix
                    wan_j=dconjg(evc_w1(jwf,it,kx)) !a_{Rj   beta}^{kn}*
                    wan_i=evc_w2(iwf,itp,kx) !a_{Ri  alpha}^{(k+q)n'}
                    wan_l=evc_w1(lwf,it,kx) !a_{R'l  beta}^{kn}
                    wan_k=dconjg(evc_w2(kwf,itp,kx)) !a_{R'k alpha}^{(k+q)n'}*
                    wanijkl=wan_j*wan_i*wan_k*wan_l

c$$$                      write(ifwanmat,"('i j k l ijkltag abs(wan) Re(wan) Im(wan)',4I3,L4,5E13.5)") iwf,jwf,
c$$$     &                     kwf,lwf,ijklmag,abs(wanijkl),wanijkl
                   !!! wanmat (dimension:nnwf)
                 else
                    wanijkl=(0d0,0d0)
                 endif
c$$$                 if (abs(aimag(wanijkl))<1d-12) then
c$$$                    rewanijkl=real(wanijkl)
c$$$                    wanijkl=cmplx(rewanijkl,0d0,kind(0d0))
c$$$                 endif
                 wanmat(ijwf,klwf)=wanijkl
!                 wanmat(ijwf,klwf)=wanijkl+(0d0,1d-9) ! +delta (~ 1d-9?)
                 ! delta makes a width of Im[R] peaks
                   
 3005  continue     !lwf
 3004  continue     !kwf
 3003  continue     !jwf
 3002  continue     !iwf
       if (egasmode) wanmat=1d0

c$$$c diag wanmat
c$$$       do iwf=1,nnwf
c$$$          wanmat(iwf,iwf)=(1d0,0d0)
c$$$       enddo
c$$$c

c========================== end caliculate M(:,:) : wanmat(nwf,nwf)
ccc   print *,"ihw,nhw",ihw(ibib,kx,jpm),nhw(ibib,kx,jpm)
       do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
          imagweight=whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
!!!   accumulate Im[K]
          kmat(:,:,iw,jpm)=kmat(:,:,iw,jpm)+imagweight*wanmat(:,:)
          if (egasmode) rcxq(1,1,iw,jpm)=rcxq(1,1,iw,jpm)+imagweight*wanmat(1,1)
       enddo
       
 2013 continue                  !ibib-loop
 2012 continue                  !jpm-loop
 2011 continue                  !k-loop

      deallocate(evc_w1,evc_w2)
      deallocate(wanmat)
      call tetdeallocate()      !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)

!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
!     ! npm is set by npmtwo=T
      schi=1                    ! flip over maj/min spins.
      if (egasmode) then
         call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i        rcxq, npm,nw_i, nnwf,nnwf, ! rcxq is alterd---used as work
     o        zxq, zxqi,
     i        .false., schi,1,1d99,1d99) !ecut(iecut),ecuts(iecut))
         deallocate(rcxq)
      else                      !!! Wannier
         call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i        kmat, npm,nw_i, nnwf,nnwf, ! rcxq is alterd---used as work
     o        zxq, zxqi,
     i        .false., schi,1,1d99,1d99) !ecut(iecut),ecuts(iecut))
         deallocate(kmat)
      endif
      deallocate(zxqi) !!! not used in hmagnon (2018/06/20)

      allocate(wkmat(1:nnwf,1:nnwf),wkmat2(1:nnwf,1:nnwf)) !WKmatrix, WKmatrix_inv
c     allocate(swkwmat(1:nnwf,1:nnwf),swkwmat2(1:nnwf,1:nnwf))!sW Kmatrix, sW K sW matrix_inv
      allocate(rmat(1:nnwf,1:nnwf,nw_i:nw)) !Rmatrix
ccc
      debug=.False.              !!01/16
ccc

       !!! get eta for (1-eta*WK)
      if (iq==1) then

         if (.not. allocated(eval_wk)) allocate(eval_wk(nnwf))
         if (.not. allocated(eval_k)) allocate(eval_k(nnwf))

         wkmat(1:nnwf,1:nnwf)
     &        =matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,0)) !omega=0

!!!   eval_wk is complex array because of Non-Hermite WK
         call diagcvuh3(wkmat(:,:),nnwf,eval_wk)
c$$$  call diagcvuh3(swkwmat2(:,:),nnwf,eval_swkw)

         if (debug) then
            write(6,"('eigenvalue of WK real&imag =',2E15.5)"
     &           ,advance='NO') eval_wk
            write(6,*)
         endif

!eta= 1d0/minval(real(eval_wk)) !moderate peak
         eta=  -1d0/abs(eval_wk(1)) !sharp peak
c         eta= -1d0
         write(6,*) "now eigenvalue abs(WK)",abs(eval_wk(1)),"is inversed"
         write(6,*) "check eigenvalue Re(WK)",real(eval_wk(1))
         write(6,*) "check eigenvalue Im(WK)",aimag(eval_wk(1))
C     eta2= -1d0/abs(eval_swkw(1)) !sharp peak
         print *,"wkmat calculated eta:", eta !negative value
      endif                     !iq==1

cccc open file each iq
      if (MPI__task(iq)) then
         ifchipmz_wan=ifile_handle()
         open(ifchipmz_wan,file="wan_ChiPMz.mat"//charnum4(iq))
         ifchipmr_wan=ifile_handle()
         open(ifchipmr_wan,file="wan_ChiPMr.mat"//charnum4(iq))          
         print *,'ifchipm=',ifchipmz_wan,ifchipmr_wan
         ifchipmrk_wan=ifile_handle()
         open(ifchipmrk_wan,file="wan_ChiPMr-k.mat"//charnum4(iq))          
      endif

!!! unit matrix (dimension: nwf*nwf), need for 1-WK matrix
      if (.not. allocated(imat)) then
         allocate(imat(1:nnwf,1:nnwf))
         imat=0d0
c     do iwf=1,nwf
         do iwf=1,nnwf
c     ijwf=(1+nwf)*iwf-nwf !!! iwf=jwf 
c     imat(ijwf,ijwf)=1d0
            imat(iwf,iwf)=1d0
         enddo
      endif
ccccccc
      if (MPI__task(iq)) then
         ifwkeigen2=ifile_handle()
         open(ifwkeigen2,file="wk_eval_list.dat"//charnum4(iq))
         ifwkeigen3=ifile_handle()
         open(ifwkeigen3,file="1wk_eval_list.dat"//charnum4(iq))
      endif
cccccccdiagonalization for K
      print *,"nw_i,nw",nw_i,nw
      maximr=0d0; w_maximr=0d0 !!! search for Im[R] peak
      do 2050 iw = nw_i,nw
         trmat=0d0
         do iwf=1,nwf
            ijwf=(1+nwf)*iwf-nwf !!! iwf=jwf 
            do jwf=1,nwf
               klwf=(1+nwf)*jwf-nwf !!! kwf=lwf
               call checkorb2(iwf,iwf,jwf,jwf,ijklmag); ijklmag=.true.
               if(ijklmag) then
                  trmat=trmat+zxq(ijwf,klwf,iw)/znorm
               endif
            enddo
         enddo
         !!!
         if (aimag(trmat)<0d0 .and. iw > 0 .and. iw < 250) write(6,*) "iq iw (Im(K)<0):",iq,iw
         !!!   
         www=freq_r(iw)
         if(iw<0) www=-freq_r(-iw)
         if(.not. gskip .and. MPI__task(iq)) write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,2E13.4)")qshort,
     &        iw,www*2d0,trmat/hartree
ccccc K/(1-WK)
cccccccccK(1-WK)^(-1)
         call getkeyvalue("GWinput","deta",deta,default=0d0)          
         if (.not. allocated(wkmat)) allocate(wkmat(1:nnwf,1:nnwf))
         wkmat(1:nnwf,1:nnwf)
     &        =matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,iw))
!!! 1-eWK ===> 1-WK
         if (.not. allocated(wkmat2)) allocate(wkmat2(1:nnwf,1:nnwf))
         wkmat2(1:nnwf,1:nnwf)=imat(1:nnwf,1:nnwf)-(eta+deta)*wkmat(1:nnwf,1:nnwf)
         deallocate(wkmat)
!!! eigenvalue of 1-WK
         if (debug) then
            if (iq==1 .and. iw==0) then
               call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
               Write(6,*) "q=0 and omega=0 case:"
               write(6,"('eigenvalue of (1-eta*WK) real&imag=',2f15.8)"
     &              ,advance='NO') eval_wk
               write(6,*)
            endif
         endif
ccc   WK eigenvalue writing (wk_eval_list.dat)
         call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)         
         if (MPI__task(iq)) write(ifwkeigen2,"(3f9.5,I6,10E13.5)") qshort,iw,www*2d0,
     &        eval_wk(1),abs(eval_wk(1)),eval_wk(2),abs(eval_wk(2)),eval_wk(3),abs(eval_wk(3))
ccc   
         call matcinv(nnwf,wkmat2(1:nnwf,1:nnwf)) ! inv(1-WK)
         rmat(1:nnwf,1:nnwf,iw)=matmul(zxq(1:nnwf,1:nnwf,iw),
     &        wkmat2(1:nnwf,1:nnwf))

!!! eigenvalue of (1-WK)inv
         if (debug) then
            if (iq==1 .and. iw==0) then
               call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
               write(6,*) "q=0 and omega=0 case:"
               write(6,"('eigenvalue of (1-eta*WK)inv real&imag=',2f15.8)"
     &              ,advance='NO') eval_wk
               write(6,*)
            endif
ccc   find eigenvalue of (1-WK) at all omega (2018/04/24)
            if (iw==0) then
               call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
               write(6,"('WKWK: eig 1-WK inv Re',3f9.4,f15.8)") qshort, real(eval_wk(1))
               write(6,"('WKWK: eig 1-WK inv Im',3f9.4,f15.8)") qshort,aimag(eval_wk(1))
            endif
         endif

ccc 1-WK eigenvalue writing (1wk_eval_list.dat)
c          if (iq==1 .and. iw .ge. 0) then

c          if (iw .ge. 0) then ! omega >= 0
         call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
         deallocate(wkmat2)
         if (MPI__task(iq)) write(ifwkeigen3,"(3f9.5,I6,10E13.5)") qshort,iw,www*2d0,
     &        eval_wk(1),abs(eval_wk(1)),eval_wk(2),abs(eval_wk(2)),eval_wk(3),abs(eval_wk(3)) !iwf=1,2,3
         
!!! eigenvalue of rmat
         if (debug) then
            if (iq==1 .and. iw==0) then
               call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
               write(6,*) "q=0 and omega=0 case:"
               write(6,"('eigenvalue of rmat',i5,2f15.8)") iw,eval_wk(1)
               write(6,*)
            endif
         endif

!!! diagonalization for R
         trmat=0d0
         trmatt=0d0
         do iwf=1,nwf
            ijwf=(1+nwf)*iwf-nwf !!! iwf=jwf 
            do jwf=1,nwf
               klwf=(1+nwf)*jwf-nwf !!! kwf=lwf
               call checkorb2(iwf,iwf,jwf,jwf,ijklmag) ; ijklmag=.true.
c     if (ijwf/=klwf) cycle
!               if (1<iwf<7) then
               !!! orbtest
!                if (iwf==1)  then
!                if (iwf==4 .or. iwf==6)  then

c$$$               if ( (iwf==2 .or. iwf==3 .or. iwf==5) .and. 
c$$$     &              (jwf==2 .or. jwf==3 .or. jwf==5))  then
c$$$               if ( (iwf==2 .or. iwf==3 .or. iwf==5) .and. 
c$$$     &              (jwf==4 .or. jwf==6))  then
c$$$               if ( (iwf==4 .or. iwf==6)   .and.
c$$$     &              (jwf==4 .or. jwf==6))  then
c$$$               if (iwf==2 .and. jwf==2) then
c$$$               if (iwf==6 .and. jwf==6) then
c$$$               if (iwf==5 .and. jwf==5) then
c$$$               if (iwf==2 .and. jwf==6) then
c$$$               if (iwf==1 .and. jwf==1) then

!               if (jwf==2 .or. jwf==3 .or. jwf==5)  then
               if(ijklmag) then
                  trmat=trmat+rmat(ijwf,klwf,iw)
                  trmatt=trmatt+(rmat(ijwf,klwf,iw)-zxq(ijwf,klwf,iw))
               endif
            enddo
         enddo
         if (MPI__task(iq)) then
            if(.not. gskip) write(ifchipmr_wan,"(3f9.4,i6,E12.4,x,12E12.4)")qshort,iw,
     &           www*2d0,trmat/hartree
            if(.not. gskip) write(ifchipmrk_wan,"(3f9.4,i6,E12.4,x,12E12.4)")qshort,iw,
     &           www*2d0,trmatt/hartree

cccc search for MAX(Im[R]) 20180706 (0 - 1500 meV)
c            if (0d0 <www*hartree .and. www*hartree < 1.5) then
            if (0d0 < www) then
               if (maximr < -1d0*aimag(trmat)) then
                  maximr=-1d0*aimag(trmat)
                  w_maximr=www
               endif
            endif
cccc
         endif
ccc END K/(1-WK)

         gskip=.false.          !!! Skip Gamma point
 2050 continue

      if (.not. iq==1) then 
         !Tc(MF)
         ! BZweight*omega[eV] for sum(E(q))/N
         mf_maximr(imaximr)=wibz(iq)*w_maximr*hartree
         !Tc(RPA)
         ! BZweight*omega[eV] for sum(1/E(q))/N
         rpa_maximr(imaximr)=wibz(iq)/(w_maximr*hartree)
      endif

      if (MPI__task(iq)) then  !magnon peak
         write(6,"('AAAA',I4,3f9.4)") iq,qshort
         write(6,"('AAAA, w(MAX(im[R])), Im[R]',f13.5,E12.4)") w_maximr*2d0*13600,maximr/hartree
      endif

      if (MPI__task(iq)) then
         write(ifwkeigen2,*) 
         write(ifwkeigen3,*) 
         write(ifchipmz_wan,*)
         write(ifchipmr_wan,*)
         write(ifchipmrk_wan,*)
      endif
      
      if(allocated(rmat)) deallocate(rmat)
      
c$$$  !! ImagOmega end =================
      if(allocated(zxq )) deallocate(zxq)
      
      debug=.False.
      
      if (MPI__task(iq)) then
         close(ifwkeigen2)
         close(ifwkeigen3)
         close(ifchipmz_wan)
         close(ifchipmr_wan)
         close(ifchipmrk_wan)
      endif      
 1001 continue                  !q point loop
      write(6,*) "maximr RPA",rpa_maximr
      write(6,*) "maximr MFA",mf_maximr
      sumrpa_maximr=cmplx(sum(rpa_maximr(:)),0d0,kind(0d0))
      summf_maximr =cmplx(sum( mf_maximr(:)),0d0,kind(0d0))
c$$$!! =================== end of loop 1001 for q point ========================

      call MPI__barrier()

!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===      
      write(6,*) "MPIcheck",MPI__size,MPI__rank,sumrpa_maximr
      if(MPI__size/=1) then
         call MPI__AllreduceSum(sumrpa_maximr(1),1)
         call MPI__AllreduceSum( summf_maximr(1),1)
      endif
      
      write(6,*) "mpi__size MPI__size",mpi__size,MPI__size
      write(6,"('sum(E(q)) for MFA:',f9.4)") real(summf_maximr(1))
      write(6,"('[sum(1/E(q))]inv for RPA',f9.5)") 1d0/real(sumrpa_maximr(1))

      call cputid(0)
      call MPI__Finalize
      write(6,"('eta for 1-eta*WK , delta eta:',2f13.8)") eta,deta
      call rx0( ' OK! hmagnon mode')

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
