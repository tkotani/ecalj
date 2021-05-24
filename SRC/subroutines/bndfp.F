!! = Main part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration =
!!
!! Memo july2012takao
!! 1)lwtkb algolism (double path method for METAL=3, lwtkb=-1 is too complicated.)
!!   This makes this code too complicated.
!!   (second path formalism; 1st path for weight, second path for integration for density)
!!   "goto 99" mechanism should be avoided.
!!
!! 2) nevmx is given correctly ---> we need to examine a mechanism to determine nevmx in call subzi.
!!   It will be effective to give nevmx to reduce computational time-->nmx; we need to 
!     !
      module m_bndfp
      real(8),protected:: ham_ehf,ham_ehk,bz_dosw(2)
      real(8),protected:: bz_ef,bz_def ! revised when called.
      logical,private:: binit=.true.,bbinit=.true.
      
      contains
      subroutine m_bndfp_ef_set(bz_ef00,bz_def00) !from iors when rst is read.
      real(8):: bz_ef00,bz_def00
      bz_ef = bz_ef00
      bz_def= bz_def00
      binit=.false.
      end
      
      subroutine bndfp(
     .     leks,lrout,dmxp,iter,maxit, frc, 
     .     dmatu,vorb,llmfgw)
      use m_supot,only: zv_a_osmrho,sv_a_oorhat
      use m_supot,only: ngabc=>lat_nabc
      use m_suham,only: iv_a_oindxo, nchan=>pot_nlma, nvl=>pot_nlml, ndham=>ham_ndham,ham_ldham
      use m_hamindex, only: napwmx,ngrp,norbmto,ibastab,ltab,ktab,offl, symops_af,symops
      use m_lmfinit, only: n0,nab,nppn,ncutovl,lso,bz_ndos,bz_def0,bz_ef0,bz_w,bz_efmax,fsmom=>bz_fsmom,
     &     bz_dosmax,bz_nevmx,lmet=>bz_lmet,bz_fsmommethod,bz_n,
     &     ctrl_nspec,ctrl_pfloat,ctrl_ldos,ctrl_zbak,ctrl_lfrce,
     &     ham_pwmode,pwemin=>ham_pwemin,pwemax=>ham_pwemax,ham_lsig,epsovl=>ham_oveps,ham_eterms,ham_elind,ham_rsrnge,
     &     ham_scaledsigma,ham_pmin,ham_pmax, !
     &     alat=>lat_alat,stdo,stdl,procid,master,bz_doswin,
     &     nsp,nspc,nkaph,nlmax,nl,nbas,nsp, ham_frzwf,
     &     sstrnmix,
     &     lcplxp ,  lekkl,lmaxu,nlibu,lldau,ssite=>v_ssite,sspec=>v_sspec, lpztail
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,bz_nkp,iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
      use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
      use m_seneinterp,only: Senebroadcast,Rdsigm2,ndimsig,nspsigm,Sigmainit,Getsene,sene,Dsene
      use m_seneinterp,only: ham_ndhrs,ham_nqsig
      use m_subzi, only: m_subzi_init, lwtkb,nevmx,numq,lswtk,rv_a_owtkb,rv_a_oswtk,ldos
      use m_rsibl_ev, only : rsibl_ev ! to plot wavefunction in the fat band mode
      use m_struc_def,only : s_rv1
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize
      use m_lmfgw,only: m_lmfgw_init,sv_p_osigx,sv_p_otaux,sv_p_oppix,spotx,jobgw
      use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate,
     o    zv_a_osmpot, qmom_rv, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o     , vval_rv , fes1_rv 
      use m_clsmode,only: m_clsmode_init,m_clsmode_set1,m_clsmode_finalize
c      use m_seneinterp,only:,seneinterp,rv_a_ohrs,iv_a_oiaxs,iv_a_ontabs,rv_p_oqsig
!! ==== followings are just notations in this code (in cases wrong).
Ci   nbas  : size of basis
Ci   nsp   : number of spins
Ci   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
Ci   lmaxu : lmax for U used to dimension vorb and dmatu. lmaxu=2 if d is, but lmaxu=3 if f is included.
Ci   lldau :lldau(ib)=0 => no U on this site otherwise 
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   Following are old. Moved to lmfinit
Ci     ssite :struct for site-specific information; see routine usite
Ci     sspec :struct for species-specific information; see routine uspec
c      
Ci   ndham :dimensioning parameter, at least as large as largest
Ci         :hamiltonian dimension
Ci   leks  :>0 make the Hohnberg-Kohn-Sham energy
Ci         :>1 use the HKS forces
Ci   lrout :>0 generate output density and attendant quantities
Ci   lfrce : 0 suppress generation of forces
Ci   lpnu=1  : =1 means new pnu's
Ci   sstrnmix:mixing
Ci   dmxp  :vector of mixing parameters; see mixrho.f for dmxp(1..25)
Ci         :Additionally:
Ci         :dmxp(33)  is the Lindhard parameter
Ci   iter  :current iteration number
Ci   maxit :maximum number of iterations
Ci   evl   :work array holding band eigenvalues for current spin, qp
Cio LDA+U inputs and outputs
Cio  dmatu :density matrix for LDA+U (changed upon output)
Ci  vorb  :orbital dependent LDA+U potential
Co Outputs
Co   frc   :forces.  Only calculated if lfrce>0.
Co         :If leks<2, forces are HF  forces
Co         :If leks>1, forces are HKS forces
Cl Local variables
Cl   k1,k2,k3: dimensions smrho,smpot
Cl   lpdiag:0 use standard diagonalization (zhev)
ccCl   lwndow:T if to make density in a specified energy window
Cl   lcplxp:0 if ppi is real; 1 if ppi is complex
Cl   jsp   :current spin index.
Cl         :In the collinear case, jsp and isp are equivalent
Cl         :In the noncollinear case, isp loops 1..2 for the
Cl         :purpose of assembling the hamiltonian.
Cl         :Once assembled, isp should not be used; and jsp=1
Cl   ispc  :2 when working on (2,2) block of noncollinear hamiltonian;
Cl         :otherwise 1
Cl   nspx: number of independent spin channels containing eigenvalues
Cl         and total DOS; nspx=nsp unless nspc=2, in which case nspx=1
Cl   onesp :do only one spin branch of isp loop (spec'd by onesp)
Cl         :also used when usual loop order (iq=1..nq, isp=1..2)
Cl         :needs to be reversed, as it does, e.g. when transforming
Cl         :sigma matrix.  Then onesp plays the role of spin index
Cl   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
      
Cr Remarks ---- now improving...
Cr   Band pass consists of:
Cr   (1) make the effective potential,
Cr   (2) generate eigenvalues (and eigenvectors if lrout)
Cr   (3) if lrout, assemble the output density by BZ integration
Cr   (4) evaluate hf (and KS, if leks) energy by BZ integration
Cr   (5) mix the output density to make a new input density.

      
Cu Updates before github 2009. See history github ecalj after 2009.
Cu   05 Jul 08 (T. Kotani) new PW basis
Cu             Option to accumulate energy-weighted output density
Cu   27 Jun 08 Redesigned transformation of sigma to new basis
Cu   09 Jul 07 MPIK enabled to plot bands
Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1
Cu   09 Jun 07 Fixed-spin-moment, noncollinear case
Cu   16 Jan 07 First cut at I/O of sigm transformed to,from LDA basis
Cu   26 May 07 Some preparation for rotation betw/ LDA, GW basis
Cu   17 Jul 06 Some MPI changes, fixes SO case
Cu   27 Jun 06 New constraints for floating pnu
Cu   08 Jun 06 Bug fix total DOS, noncollinear case;
Cu             Mulliken works for noncollinear case
Cu   02 Jan 06 better printout of magnetic moment
Cu   27 Nov 05 LDA+U => complex potential
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   06 Oct 05 (A. Chantis) bug fix dos when nspc=2
Cu   25 Jul 05 bug fix partial dos combined with LDA+U
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   29 Jun 05 (MvS) extended LDA+U to local orbitals
Cu   27 Apr 05 LDA+U (Lambrecht)
Cu   14 Feb 05 fixes for band plot, contour mode
Cu   03 Feb 05 (A. Chantis) implemented spin-orbit coupling by L.S
Cu   11 Jan 05 double-counting term rho*sig subtracted from ehks.
Cu   23 Dec 04 Extended to spin-coupled case
Cu   18 Nov 04 Sampling integration properly handles Fermi distribtion
Cu   25 Sep 04 (ATP) some patches for MPI parallelization
Cu    1 Sep 04 Adapted to handle complex ppi. Spin-orbit folded into ppi
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   29 Jun 04 (A Chantis) First implementation of spin-orbit coupling
Cu             (Lz.Sz only)
Cu   19 Sep 03 (ATP) Modifications of CLS spectroscopy
Cuxxx   24 May 03 New --window switch for accumulating density
Cu             in a specific energy window
Cu   24 May 03 New interpolation mode for sigma
Cu   14 Aug 02 Added file input of self-energy addition to LDA
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   24 Oct 01 Updated mixrho, dfrce
Cu   24 Aug 01 Extended to local orbitals.
Cu   22 Apr 01 Added driver for Kotani's GW
Cu   21 Mar 01 bug fix in call to makdos
Cu   20 Mar 01 (ATP) Added Mulliken analysis, CLS
Cu   15 Feb 01 eliminated smrho, smpot from passed arguments
Cu   23 Jan 01 bug fixes connected with lrout=0
!! -------------------------------------------------------------------------
      implicit none
      logical :: llmfgw
#if MPI | MPIK
      include "mpif.h"
      integer ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      double precision starttime, endtime
#endif
      integer obuf,oqptr
      integer igetss,lmxh,nlmh,nelt(3),kmax !,nlma
      integer, dimension(:),allocatable :: kpproc
      double precision sttime,entime
      integer:: pidorigin, ndimhx_iq, nev_iq,idx
      integer:: work_mpi(10)    ! use only 1 and 2 now
      integer nbx,iter,maxit!,ndham !nbas
      integer:: k1,k2,k3,leks,lrout,lfrce,lpnu=1
      real(8):: dmxp(33),frc(3,*)
      logical lgors,ltet,cmdopt,cmdopt0,cmdopt2!,lwndow
      character strn*120,strn2*120
      integer:: ipr , iprint , isp , jsp , ispc , isum , nglob 
     .  , i , iq , ismidb  ! ldos  !lmet !, nspc  , nsp
     .  , lrep , ndos , ndos0 , nev , nevl 
     .  , nkp  , n1 , n2 , n3 !,  lswtk=0 !, nl lwtkb ,
     .  , lgunit , igets , mpsord ,  isw  !stdo , stdl 
     .  , ifi , fxst ,  lfrzw , i1 , i2 , plbnd 
     .  , nfilem , iobzwt , lrsig ,  nll, lnoxc !, lpdiag lwsig  , 
     .  , iv(10) , parg , ndhamx , nspx , nk1 , nk2 , nk3 , lshft(3) 
     .  , onesp , mpipid , ldham(16) , ldim , ndimh , ndimhx !nfbn(2) , 
      real(8) ,allocatable :: dos_rv(:,:)
      integer:: os , osoptc
      real(8) ,allocatable :: fh_rv(:)
      real(8) ,allocatable :: fes2_rv(:)
      real(8) ,allocatable :: hbyl_rv(:,:,:)
      real(8) ,allocatable :: qbyl_rv(:,:,:)
      complex(8) ,allocatable :: srout_zv(:)
      complex(8) ,allocatable :: t_zv(:,:)
c      real(8) ,allocatable :: qm1_rv(:)
c      real(8) ,allocatable :: qm2_rv(:)
c      real(8) ,allocatable :: ww_rv(:)
c      complex(8) ,allocatable :: ww_zv(:)
c      integer ,allocatable :: ifbls_iv(:)
      integer ,allocatable :: ips_iv(:)
      real(8) ,allocatable :: orbtm_rv(:)
      type(s_rv1) :: sv_p_orhat1(3,nbas)
      type(s_rv1) :: sv_p_oeqkkl(3,nbas)
      type(s_rv1) :: sv_p_oqkkl(3,nbas)
      double precision sigp(10)
      double precision ef00,eferm,emax,emin,esmearm,sev,
     .  sev00,sev1,sumtv,alfa,qbg,ebot,def,ef0,ehar,eks,
     .  dosw(2),efmax,dgets,sumev(2,3),
     .  sumqv(3,2),qp(3),shftqp(3),eterms(20),elind,dum,evtop,
     .  ecbot,rsrnge,xv(20),dosrng,dval !,eomin ,vrmt(nbas)
      equivalence (emin,dosw(1)),(emax,dosw(2))
!!    PW basis
      integer,allocatable:: igv2x(:,:)
      integer ,allocatable :: kv_iv(:)
      integer pwmode,napw
      double precision pwgmin,pwgmax
C     real(8),allocatable :: qpe(:,:)
      logical T, F
!!     for spin-orbit
      complex(8) ,allocatable :: auso_zv(:)
!!     for partial dos
      complex(8) ,allocatable :: ausp_zv(:)
!!     for Mulliken
      integer:: moddos , lsites(nbas) , ng , iomoms 
     .  , nchmx , lmdim , lmxch
      logical:: mullwrite=.false.
      real(8) ,allocatable :: doswt_rv(:)
      integer ,allocatable :: chan_iv(:),chan0_iv(:)
!!     for pzhev
      integer nblk,nprow,npcol,nmx
!!     for LDA+U
      real(8) ,allocatable :: ausu_rv(:)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),
     .  dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      integer:: i_spackv, i_spacks
      integer ::iwdummy
      logical :: debug=.false.  !debug is only for check write. However, it can cause error
!for 'make check',
!because fp/test/test.fp assumes 'something not inserted in console output'.

      data T /.true./ F /.false./
      integer,allocatable:: nevls(:,:)
ctakao
      integer :: ifiogw,ifiese,og,ib,ig,iorb,is,nsgrp,uspec_size,nspec
      logical :: nexist,sigmamode,ok !, newsigmasw
      character spid*8
      integer :: ib1,ib2,ii,ik1,ik2,ik3,j1,j2,ikx
      integer:: nhmax , nknknk , jx , nsp_dummy , ndimh_dummy
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:) ,senex(:,:),ovlmtoi(:,:),ovliovl(:,:) 
      integer:: isize_hamm(3)
      real(8):: rrr(3),rrrr
      complex(8):: pi2i
      real(8),allocatable:: range(:,:)!,pos(:,:)
      integer:: j1x,j2x,ix1,ix2,ix3 !ikc(3,0:1,0:1,0:1), igshift(3,0:1,0:1,0:1),
      integer,allocatable::igx(:,:,:,:)
      real(8):: wgc(0:1,0:1,0:1),qpara(3)
      real(8):: www,qmm(3),qqq(3),sumtvxxx,sevxxx
      integer:: iktm,ikt,igvecc(3),nkk(3),igvecc0(3), nnspc,osmpotbk
      logical::bzwtsf_bisection !noshorbz,
      real(8),allocatable:: qplist(:,:)
      real(8):: qpx(3)          !qpo(3),

      integer:: iqini,iqend
      complex(8),allocatable:: smrhobk(:)
      real(8),allocatable:: evlall(:,:,:),evl(:,:),dwgtall(:,:,:,:,:)
      integer:: idummy
      real(8):: epsnevec
      logical:: onlyonce=.true.
      logical:: isanrg, l_dummy_isanrg,oncewrite, savez,getz,gotosecondpass,initiqisp
      integer:: getef, i_dummy_getef,ifig!, ndimsig
c      integer:: n_listwf
      integer, allocatable:: listwf(:)
      integer:: ntqxx,ispx,ix99 !,nevec unused.
      real(8)::eseavrmean
c      logical:: mtosigmaonly
      real(8)::scaledsigma
c      real(8):: qsfz(3,100)
      integer:: nqps,ii1,ii2,ii4,iq1 !, fxsts
      logical:: lfx

      integer:: unlink
!! for mpi and mpik, parallel write
      integer:: nsizeohrs,nsizeoiaxs,nsizeontabs,noqsig
      integer:: ifis2,ipts,ifid
      real(8):: eee
      logical :: tdos !,pdos
      integer:: ifves,ifvesintatm,ifvesintloc,ispec,ifvbm,ifqbyl,il

      real(8):: rydberg=13.6058d0,evll,qrmtx,qsetsum,esold, vnow
      integer,allocatable:: ndimhx_(:),nev_(:)
      integer:: npri,ibas,ifiqb,ifiwv,lmxa, ifimag
      real(8),parameter::    NULLR =-99999

      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw
      real(8)::ppin(3)
      real(8):: rlatp(3,3),xmx2(3),vadd

      logical:: fsmode,PROCARon,allband !for --fermisurface for xcrysden.
      integer:: nkk1,nkk2,nkk3  !for fermi surface mode for xcrysden.
      integer:: ififm
      real(8):: qshort(3)
! add --fermisurface #1 #2 #3
! only implemented for lmf-MPIK now.
      integer:: ini,iend,nbxx,ival

c mar2014
      integer:: lb1,lb2,onespdummy
      integer:: iprocar,iprocar1,iprocar2,lrsigx
      complex(8),allocatable:: auspp(:,:,:,:,:)
      integer:: l,m,ilm,nlma,iband,hreal,mode,nttabs
      complex(8):: au,as,az
      real(8):: s11,s22,s33,s12,s13,s23,dwgt(100),dwgtt(100),xdat,qold(3)
      character*256:: sss,s222
      character*1000::ccc
      logical:: initbd=.true.

      integer:: ispsigm !nov2014
      logical:: spinoffdiag,    goto99!jan2015
      real(8),allocatable:: xdatt(:)
      logical:: fullmesh

      integer,allocatable:: ipqe(:,:,:),idtete(:,:)
      integer:: ntete,init,ient,nchanp,nevmin,ifip,ifisyml,isyml
      real(8):: qx(3)
      real(8),allocatable:: wkd(:), pdosall(:,:,:)
      character*3::charnum3
      character*100::filenm(2)
      integer::ifbndsp(2)
      real(8):: eminp,emaxp
!! direct read syml
      integer,parameter:: nsymlmax=100
      integer:: nqp_syml(nsymlmax),nqps_syml(nsymlmax),nqpe_syml(nsymlmax)
      integer:: nqp2_syml(nsymlmax),nqp2s_syml(nsymlmax),nqp2e_syml(nsymlmax),nqp2n_syml(nsymlmax)
      logical:: masslineon(nsymlmax)

      integer:: nsyml=0,ikp,ifbnd,nsymln
      real(8):: qps_syml(3,nsymlmax), qpe_syml(3,nsymlmax),rq
      character*20 ::labeli(nsymlmax),labele(nsymlmax)
      real(8)::  totxdatt, dqsyml(nsymlmax),etolv,etolc
      integer:: ifile_handle, fxsts

      character(512),save::schar
      character(512):: aaachar
      logical:: massline , writeham=.false.
      logical,save:: siginit=.true. 
      character(10):: i2char
      character(50)::infoq
      integer::ipos,ndos_,j
      real(8)::dosmax,emaxp_,posp(3)
      integer:: jspend,jspini,jspp
      integer:: iqread,iqindex,ibread(100)=0,ibindex,nibread, job,ist,ip,ni,ix

      integer,allocatable:: npair(:,:), nlat(:,:,:,:),ib_table(:),k_table(:),l_table(:)
      integer::npairmx,ikpd,ifqplist
      integer,allocatable:: nqwgt(:,:,:)
      real(8):: tpioa,pi=4d0*datan(1d0)
      complex(8):: aaaa,img=(0d0,1d0)

      real(8):: pnu(n0,2),pnz(n0,2),pmean,symxx(3,3),qvec(3),dee,dosi(2)
      integer:: i_copy_size, ndummy1,ndummy2,ndummy3,ncount,nbandx,ifii
      logical:: phispinsym,ljobgw
      integer:: nr,lmxl,nlml,nelt1,nelt2,nso,ififft,ifih
      character*8:: charext
      logical:: lfill=.false.,ixx
!! ----------------------------
      call tcn ('bndfp')
!! vbmmode =Get VBM and CBM relative to vaccum (a simple approximaiton. Need fixing if necessary.).
      if(cmdopt0('--vbmonly')) then
         if(master_mpi) call vbmmode()
         call rx0('--vbmonly mode done')
      endif
!! getq mode. Current version is not for spin dependent nor many restrictions.
!! spin symmetic (or nspin=1, not 2 channell binded and so on...
      if(cmdopt0('--getq')) then
         if(master_mpi) call getqmode()
         call rx0('--getq mode done')
      endif
!! checkwrite     
      if(master_mpi) then ! check Write atomic positions
         write(stdo,"(a)")' --- bndfp: Atomic Positions from rst/ctrl. '
         write(stdo,"(a)")'   Note: rst and ctrl may contain different atomic positions when DYN mode'
         do ib=1,nbas
            is  = ssite(ib)%spec   !species index
            spid= ssite(is)%clabel !atom label
            write(stdo,"('ib clabl pos =',i5,a6,3(x,f13.5))") ib, trim(spid),ssite(ib)%pos
         enddo
      endif
!!      
      debug    = cmdopt0('--debugbndfp')
      tdos     = cmdopt0('--tdos')
      PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
      fullmesh = cmdopt0('--fullmesh') !full mesh stop just after do 2010 iq loop.
      fsmode   = cmdopt0('--fermisurface')
      if(fullmesh.and.(.not.procaron)) call rx(' --fullmesh is meaningful only with --mkprocar')
      if(nlibu>0) dmatu=0d0 !dmat initialization
!! fsmode: Fermi surfece mode for xcrysden. It's format is in  http://www.xcrysden.org/doc/XSF.html#2l.16
      if(fsmode) fullmesh=.true.
      nkk1=nkabc(1)
      nkk2=nkabc(2)
      nkk3=nkabc(3)
!!      
      savez=.false.
      getz =.false.
      if (maxit==0) return
      napw = 0
      ipr  = iprint() ! for procid/=master, we set iprint=0 at lmv7.F
      idummy = unlink('RSEQ_ERROR') !the error flag in rseq. Remove RSEQ_ERROR here because
                        ! this error is only fatal when it is finaly converged. takaoNov28_2010
!! Dimension maximum of Hamiltonian is ndhamx (ispx=1,nspx)
!! Note spinoffdiag=T case: we use nspc,nsp,nspx,ndhamx (a little complicated, I think).
      if(nspc==2.and.nsp==2) then !nspc==2
         spinoffdiag=.true. !non-colinear case
         nspx=1
         ndhamx = ndham*2
      elseif(nspc==1.and.nsp==2) then
         spinoffdiag=.false. !spin no-offdiagonal
         nspx=2
         ndhamx= ndham
      elseif(nspc==1.and.nsp==1) then
         spinoffdiag=.false. !paramagnetic
         nspx=1
         ndhamx= ndham
      else
         call rx('bndfp: nspc==2 but nsp=1')
      endif
      allocate(evl(ndhamx,nspx))
      ldham = ham_ldham
      ldim  = ldham(1)
      if(size(ham_ldham)/=size(ldham)) call rx("bndfp:size(ham_ldham)/=size(ldham)")
      pwmode = ham_pwmode
      if (master_mpi) then
         open(newunit=ifiese,file='PWMODE')
         write(ifiese,*)pwmode
         close(ifiese)
      endif 
      onesp = 0
      lrsig = ham_lsig  !12 for read sigma mode
      lfrzw = 0
      if(ham_frzwf) lfrzw = 1  !freeze all augmentation wave
      nspec = ctrl_nspec
      if (master_mpi) then
         inquire(file='sigm.'//trim(sname),exist=nexist)
         if (lrsig/=0 .and. (.not.nexist) ) then
            write(6,*)' bndfp (warning): no sigm file found ... LDA calculation only'
            lrsig = 0
         endif
      endif
      call mpibc1(lrsig,1,2,mlog,'bndfp','lrsig')
      if (master_mpi) then
        aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
        write(stdo,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
      endif
      if (ndham.le.0)  call rx('bndfp: hamiltonian matrix has zero dimension')
      eterms=-99d0
      eterms(19) = 0d0
      ham_eterms=eterms
      eks = 0d0
      n1=ngabc(1)
      n2=ngabc(2)
      n3=ngabc(3)
      k1 = n1
      k2 = n2
      k3 = n3
      nkp  = bz_nkp
      efmax= bz_efmax
      ndos = bz_ndos
      if(bbinit) then !from m_lmfinit
         bz_dosw=bz_doswin
         bbinit=.false.
      endif
      dosw = bz_dosw
      if(binit) then !from m_lmfinit when m_bndfp_ef_set in not called.
         bz_ef =bz_ef0
         bz_def=bz_def0
         binit=.false.
      endif
      ef0  = bz_ef
      def  = bz_def
      ltet = ntet>0
      qbg = ctrl_zbak(1) !homogenious background charge
      alfa = 0
      lfrce= ctrl_lfrce
!!  Switch to plot bands at specified qp
      if(cmdopt0('--band').or.fullmesh) then
         plbnd = 1
         lrout = 0
         lfrce = 0
c         numq = 1
c         allocate(ifbls_iv(ndhamx*2))
c         ifbls_iv(:)=0
      else
         plbnd = 0
      endif
!! PROCAR files 
      if(PROCARon)            open(newunit=iprocar1,file='PROCAR.UP.'//trim(strprocid))
      if(PROCARon.and.nsp==2) open(newunit=iprocar2,file='PROCAR.DN.'//trim(strprocid))     
!! spin-symmetrized pnu. aug2019. See also in pnunew and locpot
      if(cmdopt0('--phispinsym')) call phispinsym_set(ssite) !pnu,pz spin symmetrized      
!! Write out Hamiltonian stated at jan2018kaist
      writeham= cmdopt0('--writeham')
      if(writeham) then
         open(newunit=ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
         if(master_mpi) call writeham_init(ififft) !ififft is file handle returned.
      endif      
!! Make one-particle potential without XC part for GWdriver---
      if(llmfgw) then        
         call m_lmfgw_init() ! Initializer of m_lmfgw. Get sig,tau,ppi,spot without XC(LDA).
#if MPIK
         call mpi_barrier(MPI_comm_world,ierr)
#endif        
      endif
!! Generate one-body potential and energy-related quantities. See use m_mkpot:
      call m_mkpot_init(llmfgw,vorb) 
!! elind for accelarating convergence
      elind=ham_elind
      if(ham_elind<0d0) elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*ham_elind
      dmxp(33) = elind
!!Setup for BZ integration wtkb
      if(cmdopt0('--band').or.fullmesh) nkp = 0 ! dummy for m_subzi_init
      if (plbnd==0) call m_subzi_init(lrout>0,nkp,ef0,def,bz_n,bz_w,eferm)
c      allocate( evlall(ndhamx,nspx,nkp))
c      evlall = 0d0
      if (ldos==0 .and. lmet==0) ndos = 1
      ndos = abs(ndos)
      allocate(dos_rv(ndos,2*nsp))
      if (lrout/=0) then
         call dfqkkl( sv_p_oqkkl )
         if(lekkl==1) call dfqkkl( sv_p_oeqkkl )
         allocate( srout_zv(k1*k2*k3*numq*nsp) )
         allocate( fh_rv(3*nbas)  )
         allocate( fes2_rv(3*nbas))
      endif
!   TK think sspec%ngcut of local orbital is changing when PZ=+10 mode (hankel tail mode)
!   because Hankel's e may be changing. ! lpztail: local orbital of 2nd type(hankel tail) +10 exist or not
c      call suham2(sspec)
      if(lpztail) call sugcut(2, sspec) !sspec%ngcut
!! cls mode
      if(cmdopt0('--cls')) then
         if (lrout .eq. 0) call rx('bndfp: need output density for cls')
         efmax = 1d3
         nevmx = ndhamx
         call m_clsmode_init()
      endif   
      
!! --- get index for hamiltonian for m_hamindex takao june2009
!! these are used in sigm mode(QSGW).
!! memo:
!!  ib = atom index
!!  ltab= L (angular momentum index)
!!  ktab=  =1 for EH, =2 for EH2, =3 for lo
!! SYMOPS
      if(onlyonce) then 
         onlyonce=.false.
         if(master_mpi) then
            call gen_hamindex(ssite, sspec, ldim, .true., jobgw==0)
c           write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
            allocate(ib_table(ldim),l_table(ldim),k_table(ldim))
            do iorb = 1, norbmto !Total number of MTO's
               ib   = ibastab(iorb)
               is   = ssite(ib)%spec
               spid=sspec(is)%name
               ib_table(offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ib
               l_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ltab(iorb)
               k_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ktab(iorb)
            enddo
            write(6,*) ' i  ib l  k trim(spec)'
            do i= 1, ldim
               ib   = ib_table(i)
               is   = ssite(ib)%spec
               spid = sspec(is)%name
               write(6,"(i3,x,3i3,x,a)")i, ib_table(i),l_table(i),k_table(i),trim(spid)
            enddo
            if(writeham) then
               write(ififft)ldim,lso,nsp
               write(ififft)ib_table,l_table,k_table
               close(ififft)
            endif  
         endif
      endif
!!
C ... Case only generate bands at supplied qp: setup
      if (plbnd .ne. 0) then
         if (master_mpi) then !Try and read Fermi level from file
            inquire(file='wkp.'//trim(sname),exist=nexist)
            if(nexist) then
               open(newunit=ifi,file='wkp.'//trim(sname),form='unformatted')
               read(ifi) ndummy1,ndummy2,ndummy3,ef0
               close(ifi)
            endif
         endif   
         call mpibc1(ef0,1,4,.false.,'bndfp','ef0')
      endif   
!! ... Setup for spin-orbit coupling
      if (lso /= 0) then
c         if(allocated(orbtm_rv)) deallocate(orbtm_rv)
         allocate(orbtm_rv(nl*nsp*nbas))
      endif
!     !
!! === Get nkp and qplist =========== feb2015
!! fatband mode ---> --mkprocar
!! plbnd/=0 ---> band plot mode.
      if(debug) write(*,*) 'bndfp: uuuuuuu1111111'
c      if(allocated(qplist)) deallocate(qplist)
      if( master_mpi .and. plbnd==0) then
         allocate(qplist(3,nkp))
         call dcopy(3*nkp, rv_p_oqp,1,qplist,1)
c-----------------
         open(newunit=ifqplist,file='QPLIST.IBZ')
         do iq=1,nkp
            write(ifqplist,"(i5,3f23.15,3x,3f23.15)")iq, qplist(:,iq), rv_a_owtkp(iq)/2d0
         enddo !!I think sbz%rv_p_oqp is not used in the following.
         close(ifqplist)
      elseif(master_mpi) then
         if (cmdopt('--onesp',7,0,strn) .and. nspc .eq. 1) onesp = 1
         if(debug) print *,'vvvvv nkp=',nkp
         if(fullmesh) then
            nkp = nkk1*nkk2*nkk3
            allocate(qplist(3,nkp))
!! Fermi surface version for xcrysden
            if(fsmode) then
               iq=0
               do ik1=1,nkk1    !ordering is differnt from procaron case...
                 do ik2=1,nkk2
                   do ik3=1,nkk3
                     iq=iq+1
                     qplist(:,iq)   =  qlat(:,1)*dble(ik1-1)/(nkk1-1)
     &                 +   qlat(:,2)*dble(ik2-1)/(nkk2-1)
     &                 +   qlat(:,3)*dble(ik3-1)/(nkk3-1)
c                     call dpscop ( qp , sbz%rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 )
                   enddo
                 enddo
               enddo
            endif   
!! pdos (--mkprocar and --fullmesh)
            if(procaron) then !to fit to the tetirr.F
               iq=0
               do ik3=1,nkk3
                 do ik2=1,nkk2
                   do ik1=1,nkk1
                     iq=iq+1
                     qplist(:,iq) =  qlat(:,1)*dble(ik1-1)/nkk1
     &                 +   qlat(:,2)*dble(ik2-1)/nkk2
     &                 +   qlat(:,3)*dble(ik3-1)/nkk3
c                     call dpscop ( qp , sbz%rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 )
                   enddo
                 enddo
               enddo
            endif
         else 
!! syml direct read for plbnd mode. See "call writeband" below. feb2015 
!!   Read efermi.lmf
            write(6,*)' --- Readin efermi.lmf --- '
            ifi=ifile_handle()
            open(unit=ifi,file='efermi.lmf',status='old',err=1012)
            read(ifi,*,err=1012) eferm
            read(ifi,*,err=1012) evtop
            read(ifi,*,err=1012) ecbot
            close(ifi) 
            goto 1013
 1012       continue
            call rx('No efermi.lmf!: copy it, or run lmf-MPIK (sc mode) to get efermi.lmf.')
 1013       continue

!! --- example of syml file ---
! ndiv qleft(1:3) qright(1:3) llabel rlabel  ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry) 
! 5  0 0 0   .5 .5  .5        Gamma  L       1025  1  16     0.1      0.01 
! 5  0 0 0    1.  0  0        Gamma  X    
! 5  0 0 0   .75 .75 0        Gamma  K       1025  8  38     0.1      0.01 
!! As this shows you can add or not add line after ndiv2. These are for highly resolved calculations.
            write(6,*)' --- Readin syml file --- '
            ifisyml=ifile_handle()
            open(ifisyml,status='old',file='syml.'//trim(sname))
            nsyml=0
            nsymln=0
            nqp2_syml=0
c            massline=.false.
            do 
              if(nsyml+1>nsymlmax) call rx('bndfp: Enlarge nsymlmax')
              read(ifisyml,"(a)",end=1015) schar
              if(len(trim(schar))==0.or.schar(1:1)=='#' .or. schar(1:1)=='!' .or. schar(1:1)=='%') cycle !comment line
              read(schar,*,err=1014,end=1014) 
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         labeli(nsyml+1),labele(nsyml+1), 
     &         nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1),etolv,etolc
              masslineon(nsyml+1)=.true.
              nqp2n_syml(nsyml+1)= nqp2e_syml(nsyml+1)-nqp2s_syml(nsyml+1)+1
              nqps_syml(nsyml+1)=1
              nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
              write(6,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a,'  Massl:div,init,end=',3i5)")
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         trim(labeli(nsyml+1)),trim(labele(nsyml+1)),
     &         nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1)
              goto 1025

 1014         continue
              read(schar,*,err=1015,end=1015) 
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         labeli(nsyml+1),labele(nsyml+1) 
              masslineon(nsyml+1)=.false.
              nqps_syml(nsyml+1)=1
              nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
              write(6,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a)")
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         trim(labeli(nsyml+1)),trim(labele(nsyml+1))
              nsymln = nsymln+1
              nqp2n_syml(nsyml+1)= 0
 1025         continue
              if(nqp_syml(nsyml+1)==0) exit
              nsyml = nsyml + 1
            enddo
 1015       continue
            close(ifisyml)
            nkp = sum(nqp_syml(1:nsyml)+nqp2n_syml(1:nsyml))
            allocate(qplist(3,nkp))
            if(allocated(xdatt)) deallocate(xdatt)
            allocate(xdatt(nkp))
            totxdatt=0d0
            ikp=0
            do isyml=1,nsyml
              dqsyml(isyml) = dsqrt(sum((qpe_syml(1:3,isyml) -qps_syml(1:3,isyml))**2))
              do i=1,nqp_syml(isyml)+nqp2n_syml(isyml)
                ikp= ikp+1
                if(i<=nqp_syml(isyml)) then
                  rq = dble(i-1)/(nqp_syml(isyml)-1)
                else  
                  ii= i-nqp_syml(isyml)-1
                  rq = dble(nqp2s_syml(isyml)-1+ii)/(nqp2_syml(isyml)-1)
                endif
                qplist(:,ikp)= (1d0-rq)*qps_syml(1:3,isyml) +rq*qpe_syml(1:3,isyml)
                xdatt(ikp) = totxdatt + dqsyml(isyml)*rq
              enddo  
              totxdatt = totxdatt + dqsyml(isyml)
            enddo  
            write(6,"('nsyml nkp=',3i5)") nsyml,nkp
         endif
         if (nkp <= 0) call rx('bndfp: nkp<=0') ! quit if nkp==0
      endif

!! broadcase nkp and qplist
      call mpibc1( nkp,1,2,mlog,'bndfp','nkp')
      call mpibc1( onesp,1,2,mlog,'bndfp','onesp')
      if(procid/=master) allocate(qplist(3,nkp))
      call mpibc1 ( qplist, 3 * nkp , 4 , .false. , 'bndfp' , 'qp'  )
      if (master_mpi .and.oncewrite(11)) then
         ifqplist=ifile_handle()
         open(ifqplist,file='QPLIST')
         print *,'-------- qplist --------'
         iq=0
         do isyml=1,nsyml
         do i=1,nqp_syml(isyml) + nqp2n_syml(isyml)
            iq=iq+1
            infoq=''
            if(i==1) infoq=' <-- isyml= '//charnum3(isyml)
            if(i==nqp_syml(isyml)+1) infoq=' <-- isyml Mass= '//charnum3(isyml)
            write(6,"(i5,3f8.3,' ',a)")iq,qplist(:,iq),trim(infoq)
            write(ifqplist,"(i5,3f23.15,x,f12.6,' ',a)")iq,qplist(:,iq),xdatt(iq),trim(infoq)
         enddo
           write(ifqplist,*)
         enddo
         close(ifqplist)
      endif 
      if(procaron.and.nsyml==0.and.master_mpi) then !xdatt is dummy
        if(allocated(xdatt)) deallocate(xdatt)
        allocate(xdatt(nkp))
        xdatt=0d0
      endif

      if(PROCARon) then         !if(plbnd/=0) then ! march2015
        if(procid/=master) allocate(xdatt(nkp))
        call mpibc1( xdatt, nkp , 4 , .false. , 'bndfp' , 'xdatt'  )
      endif  
#if MPIK
      call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
      
!! allocation for all eigenvalues  More k-point independent local arrays ---
c      if (allocated(evlall)) deallocate(evlall)
      allocate(evlall(ndhamx,nspx,nkp))
      evlall=0d0
      nchanp=25 !total of s,p,d,f
      if(procaron.and.fullmesh) then
         allocate(dwgtall(nchanp,nbas,ndhamx,nsp,nkp))
         dwgtall=0d0
      endif
      
!! === \Sigma-Vxc to ham_rv_a_ohrs. Neighbor table ham_iv_a_oiaxs ====
!! ndimsig is the dimension of the self-energy. 
!! ndimsig = ldim means we use only projection onto MTO spaces even when PMT. !aug2012takao
!! --- New sigma inerpolation mode. Simplified.
!! OUTPUT in this block is used in
!!   call bloch 
!!  where we use ham_iv_a_oiaxs , ham_rv_a_ohrs, ham_ndhrs and so on.
!!  In future, we should make this block to a subroutine which shows clear i/o.
      sigmamode = mod(lrsig,10) .ne. 0
      if(master_mpi.and.ix99==1) print *,'sigmamode=',sigmamode
      if( llmfgw.and.jobgw==0 ) then
         continue
      elseif( sigmamode.and.siginit ) then
         call sigmainit(ldim,procid,master)
         siginit=.false. !this is equivalent to broad casting.
      endif                     ! end of sigmamode
!! == GW driver ==
      if(llmfgw) then
         call sugw(lrsig,ssite , sspec ,   nbas , ndham , zv_a_osmpot 
     .     , vconst , lcplxp , sv_p_osig , sv_p_otau , sv_p_oppi , ppnl_rv 
     .     , vrmt , spotx    , sv_p_osigx, sv_p_otaux, sv_p_oppix, jobgw 
     .     , nk1,nk2,nk3,sigmamode,ndimsig, numprocs,procid,master)!,lso) 
! NOTE: Sigma-Vxc in real space is give by rv_a_ohrs.
         call tcx('bndfp')
#if MPIK
         call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
         call rx0('sugw mode')  !exit program here normally.
      endif
      if(ham_scaledsigma/=1d0.and.sigmamode) then
        if(master_mpi) write(6,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
      endif
      
#if MPIK
      call info0(30,1,0, ' ... Start MPI k-loop')
      sttime = MPI_WTIME()
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(nkp,numprocs,1,kpproc(0))
      iqini = kpproc(procid)
      iqend = kpproc(procid+1)-1
#else
      iqini=1
      iqend=nkp
#endif
      if(pwmode>0.and.pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)
      allocate(ndimhx_(nkp),nev_(nkp),nevls(nkp,nspx))
      ndimhx_=0
      nev_   =0                
      nevls  =0
      
!! == Start loop over k points; also, re-entry for second band pass ==
!! two pass mechanism, historically. ix99=1 and ix99=2 for tetrahedron method.
!!     1st for eigenvalue (and save eigenfunctions to eigze_* files), resulting tetrahedron weight.
!!     2nd for re-read eigze_* and accumulate density.                           
!!    It will be better to remove this double path mechanism in future.
      gotosecondpass =.false.
      if(plbnd==0) then
        gotosecondpass= (lmet>=0 .and. (lmet/=4 .or. ltet) .and. lrout>0 )!this covers limited cases.
!this judges going through second pass or not.
      endif
      ix99 = 0
      xdat = 0d0
 99   continue
      ix99 = ix99+1               !ix99=2 for second pass
      ebot = 1000d0
      call surho ( nbas , ssite , sspec , lmet , ldos , lrout , lekkl 
     .  , numq , sv_p_oqkkl , sv_p_oeqkkl , k1 , k2 , k3 , srout_zv , 
     .  ndos , dos_rv , sumev , sumqv )
      if (lfrce .gt. 0) then
         frc(1:3,1:nbas*numq) = 0d0
         fh_rv= 0d0              ! 3 * nbas 
      endif
      if (lswtk==1) then
         efmax = 1d3
         nevmx = ndhamx
         call dpzero ( rv_a_oswtk , ndhamx * nkp )
      endif
!! --> Mullikan mode was here
!! ... For insulator, valence band top and conduction band bottom
      evtop = -99999
      ecbot = -evtop
      if(lso/=0) orbtm_rv=0d0
      initiqisp=.true.
C Start k and isp loop (first loop in parallel mode)
C     h,s,evecs are dimensioned ndimh in this loop (may be q dependent)
C     evl is dimensioned evl(ndham,2)
C     evlall, w(owtkb), w(oswtk) are dimensioned (ndham,nsp,nkp)
!! === do loop 2010 for iq and do 2005 isp ===
      if(debug) print *,' --- goto do loop 2010 ---'
      do 2010 iq = iqini, iqend !This is a big iq loop
         qp=qplist(:,iq)
         if(debug) print *,' do 2010 iq procid=',iq,procid
#if MPIK
         if(debug) write(*,'(a,i3,5i5)')'debug iq=',procid,iq,iqini,iqend
         if (iq .eq. kpproc(procid)) then
            if(mlog) then
               call gettime(datim)
               aaachar=' bndfp '//datim//' Process '//
     &              trim(i2char(procid))//' of '//trim(i2char(numprocs))//' on '//
     &              trim(charext(procid))//' starting k-points '//
     &              trim(i2char( kpproc( procid)))//' to '//trim(i2char( kpproc(procid+1)-1))
               write(lgunit(3),"(a)")trim(aaachar)
            endif
         endif
#endif
         if(pwmode>0.and.pwmode<10) then
            ppin=matmul(transpose(plat),qp) !basis on the qlat coordinate. qp in Cartesian.
            call shortn3(rlatp,xmx2,ppin,noutmx, nout,nlatout)
            if(debug) then
               do iout=1,nout
                  write(*,"(a,3i5,f10.4,3f8.4)")'rrrrn1 =',nlatout(:,iout),
     &              sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &              matmul(qlat(:,:),ppin+nlatout(:,iout))
               enddo
            endif
         endif
         if(allocated(igv2x)) deallocate(igv2x)
         if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
            pwgmin = pwemin**.5d0
            pwgmax = pwemax**.5d0
            qqq = 0d0 
            if (mod(pwmode/10,10) == 1) qqq = qp !pwmode 1 in 10th digit means q-dependent nw
            call pshpr(1)
            call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,napw,napw,dum,dum,dum,dum) !get napw
            allocate(igv2x(3,napw), kv_iv(3*napw))
            call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napw,napw,kv_iv,dum,dum,igv2x)
            call poppr
            if (pwmode<10) then
               do iapw=1,napw
                  igv2x(:,iapw)=igv2x(:,iapw)+nlatout(:,1)
               enddo
            endif
            deallocate(kv_iv)   !we only keep 
            ndimh = ldim + napw !ndimh is iq-dependent for pwmode >=10
            if (mod(pwmode,10) .eq. 2) ndimh = napw !APW-only mode
            if (ndimh .gt. ndham) then
               call rx('bndfp:something wrong ndimh>ndham')
            endif
         else
            ndimh = ldim
            allocate(igv2x(1,1)) !dummy
         endif
         ndimhx = ndimh*nspc    !this is iq-dependent.
         nmx    = min(nevmx,ndimhx) 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cC... range check qqq is in [ikc1/nk1, (ikc1+1)/nk1],... qpara is [0,1],[0,1],[0,1]
c      print *,' --- goto qparam--- '
c      call qparam(qp,qsmesh,nk1,nk2,nk3,plat, ikc,qpara,igshift)
c      ! eight corners -->  qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
c      call interparam(qpara, wgc)
c      cycle
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!! ... Loop over majority, minority spins
!!     In the SO case, still a loop isp=1..2
!!     isp=1 generates just the 11 block
c$$$        qpo=qp !moved here on mar2012 only for band plot
c$$$        if(.not.noshorbz()) call shorbz(qp,qp,qlat,plat)
         do 2005 isp = 1, nsp
            if(.not.(onesp==0.or.isp==onesp)) cycle
            if(debug) print *,' bndfp: do 2005 isp loop isp pwemax pwmode=',isp,pwemax,pwmode
!! takao removed shorbz(shortest qp) here. !I think not necessary
!! === For this qp, G vectors for PW basis and hamiltonian dimension ===
            ispc = min(isp,nspc) !nspc=2 for non-colinier case
            jsp = isp
            if (ispc .eq. 2) jsp = 1
ccccccccccccccccccccccccccccccc
            if(ix99==2.and.(nspc==2.and.isp==1)) cycle
            if(ix99==2) goto 1919
ccccccccccccccccccccccccccccccc            
!!   ... Make Hamiltonian and overlap matrices
c            nqsig=ham_nqsig
            if(isp==1) then
               if(allocated(hamm)) deallocate(hamm)
               if(allocated(ovlm)) deallocate(ovlm)
               nnspc=1
               if(nspc==2) nnspc=4
               allocate(hamm(ndimh,ndimh,nnspc),ovlm(ndimh,ndimh,nnspc))
            endif
            if(ispc==1) then
               hamm=0d0
               ovlm=0d0
            endif
!! === Sigma-Vxc ===
!! Generate sene(Sigma-Vxc) for given sfz.
!! Determine interpolated self-energy sene at qp from sfz. 
!! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only) 
!!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
!! Main input  => ham_iv_a_oiaxs,ham_rv_a_ohrs
!! Main output => sene
            if(sigmamode) then
               call getsene(lrsig,qp,isp) !Get sene = Sigma-Vxc 
            endif
!!note kino
!! hambl calls augmbl.  augmbl calculates C_kl(=array b in the program)
!! finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
!! Note that F2=Hankel head at a site + Hankel tail contributions from the other site.

!! == Set up Hamiltonian. ==============
!!    Hamiltonian: hamm(1:ndimh,1:ndimh,ispc). ispc=3 means off-diagonal section when SO=1.
!!    Overlap matrix: ovlm 
!!    senex: Sigma-Vxc
!! ==========================================
            jx = 0
            if (lso == 1) jx = 10
            call hambl ( jx , nbas , ssite , sspec ,  isp ,!sham ,slat , 
     .        qp , k1 , k2 , k3 , zv_a_osmpot , vconst , sv_p_osig , sv_p_otau 
     .        , sv_p_oppi , lcplxp , alfa , ndimh , napw , igv2x , hamm ( 1 
     .           , 1 , ispc ) , ovlm ( 1 , 1 , ispc ) , hamm ( 1 , 1 , 3 ) )
c            print *,'hamble sumcheck =',sum(abs(ovlm)),sum(abs(hamm))
            if(sigmamode) then
               allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
               ovlmtoi = ovlm(1:ndimsig,1:ndimsig,ispc)
               call matcinv(ndimsig,ovlmtoi)
               ovliovl = matmul(ovlmtoi,ovlm(1:ndimsig,1:ndimh,ispc))
               deallocate(ovlmtoi)
               allocate(senex(ndimh,ndimh))
               senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
               call dsene()
               deallocate(ovliovl)!,sene)
               hamm(:,:,ispc) = hamm(:,:,ispc) + ham_scaledsigma * senex
               deallocate(senex)
            endif
C   ... In the noncollinear case (spin-orbit SO=1), this marks the end of looping over
C       separate spins.
C       For the rest of the isp loop,
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
            if (ispc .ne. nspc) cycle !do2005 goto 30
c            jsp = isp
c            if (ispc .eq. 2) jsp = 1
            
!! === Diagonalize and add to density ===
            if(iprint()>=30.and.ix99==1) 
     &        write(*,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4," ndimh = nmto+napw = ",3i5,f13.5)')
     &        iq,nkp,qp,ndimh,ndimh-napw,napw
            nmx = min(nevmx,ndimhx) !lager is safer, but I think nmx can be given in a better manner.
c!!  Calculate all eigenvalues for plotting mode, pdos, tdos modes.
c            if (tdos) then 
c               nmx = ndimhx
c               efmax = 99999
c            endif
            if(debug) write(*,*)'kino, allocate t_zv' ,ndimhx,nmx,ndimhx,ix99
!! == lpdiag=0 only. ===
            if(nspc .eq. 2) then
               call sopert2 ( 0 , ndimh , nspc, hamm , hamm )
               call sopert2 ( 0 , ndimh , nspc,  ovlm , ovlm )
            endif
            if(debug) write(*,*) 'allocate ww_rv', 11*ndimhx
            
!! === diagonalization ===
!! It is better to clean up lwtkb-algolism in this routine. !!!!
!     ! savez save eigenfuncitons for second pass.(second-pass procedure should be removed in future).
ccccccccccccccccccccccccccccccc            
 1919       continue
ccccccccccccccccccccccccccccccc            
!!  Calculate all eigenvalues for plotting mode, pdos, tdos modes.
            if (tdos) then 
               nmx = ndimhx
               efmax = 99999
            endif
            if (allocated(t_zv)) deallocate(t_zv)
            allocate(t_zv(ndimhx,nmx))
            savez=.false.
            getz =.false.
            if(gotosecondpass) then
               if(ix99==1) then
                  savez=.true.  ;  getz=.false.
               else
                  savez=.false. ;  getz=.true.
               endif
               if(initiqisp) then
#if (MPI|MPIK)
           !!  for parallel write, and eigze is temporaly used only in this subroutine. 
                  inquire(file='eigze_'//trim(strprocid),opened=nexist,number=ifig)
                  if(nexist) close(ifig)
                  open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
#else
                  inquire(file='eigze',opened=nexist,number=ifig)
                  if(nexist) close(ifig)
                  ifig = ifile_handle()
                  open(ifig,file='eigze',form='unformatted')
#endif
                  initiqisp=.false.
               endif  
            endif
            
            if(PROCARon) then   !!2014mar for band weight mode
               nmx=ndimhx
               nev=ndimhx
               savez=.false.
               getz=.false.
               nmx = ndimhx
               if(allocated(t_zv)) deallocate(t_zv)
               allocate(t_zv(ndimhx,ndimhx))
            endif

!! == Diagonalize Hamiltonian ==
!! ndimhx: dimension of Hamitonian
!! hamm:Hamiltonian, ovlm: overlap matrix
!! t_zv:eigenfunciton. evl: eigenvalue.
!! nmx: input, number of requested eigenvalues(funcitons)
!! nev: out number of obtained eigenfvalues(funcitons)
!! ---
            if(getz) then   ! second pass
              read(ifig) nev,nmx
              read(ifig) evl(1:nev,jsp)
              read(ifig) t_zv(1:ndimhx,1:nmx)
            else            ! 1st pass
!! save hammovlm jan2018kaist
              if(writeham) then
                write(6,"(a,3f9.5)") "Hamiltonian: Writing hamm and ovlm for qp= ",qp
                write(ifih) qp,ndimhx,lso,epsovl,jsp
                if(lso==1) then !L.S case nspc=2 spin block
                  write(ifih) ovlm(1:ndimhx,1:ndimhx,1:4) 
                  write(ifih) hamm(1:ndimhx,1:ndimhx,1:4)
                else  !nspc=1 spin diagonal case
                  write(ifih) ovlm(1:ndimhx,1:ndimhx,1)
                  write(ifih) hamm(1:ndimhx,1:ndimhx,1)
                endif
              endif
              if(epsovl<1.000001d-14) then
                call zhev_tk2( ndimhx , hamm , ovlm , nmx , nev,
     .           evl(1 , jsp ), t_zv, iq==1,savez,getz,ifig) !stock eigenfunctions z
              else    !nov2015 recover epsovl. Remove space by overlap matrix diagonalization.
                call zhev_tk3( ndimhx , hamm , ovlm , nmx , nev,
     .           evl(1 , jsp ), t_zv, iq==1,savez,getz,ifig,epsovl) !stock eigenfunctions z
              endif
              if(savez) then
                write(ifig) nev,nmx
                write(ifig) evl(1:nev,jsp)
                write(ifig) t_zv(1:ndimhx,1:nmx)
              endif  
#ifdef MPIK            
              if(writeham.and.master_mpi) then
                call prtev(t_zv, ndimhx , evl(1, jsp ) , nmx , efmax , nev ) !efmax is required?
              endif
#endif
            endif
!! It is inefficient to calculate all eigenfunctions and eigenvalues for molecule cases.
!! nevl=nev: number of eigenvalues (if nmx=0, nev=n. See description in zhev_tk2.)
            nevl=nev            !ndimhx    this OK? why nevl and nev???
            evl(nevl+1:ndhamx,jsp)=1d99 !to skip these data
            nevls(iq,jsp) = nevl  !nov2014 isp and jsp is confusing...
            if(debug) print *,'iq jsp nevls=',iq,jsp,nevls(iq,jsp) !nov201
#ifndef MPIK            
            call prtev(t_zv, ndimhx , evl ( 1 , jsp ) , nmx , efmax , nev ) !efmax is required?
#else
            if( master_mpi.and. epsovl>=1.000001d-14.and.plbnd/=0) then
              write(6,"('      : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
            endif  
#endif
!  ! save to evlall
            if(debug) print *,' nmx nev plbnd=',nmx,nev,plbnd
            evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
            
!!   ... Save data for this qp to moments file
            if( plbnd .eq. 0) then
!!   --- Orbital magnetic moment (requires lso) ---
               if (lso .ne. 0 .and. lwtkb .ne. -1) then
                  if (lwtkb .eq. 0) call rx('metal weights required for orb. moment calculation')
                  allocate(auso_zv(nlmax*ndhamx*3*nsp*nbas))
                  if(debug)   write(*,*) 'allocate auso_zv',nlmax*ndhamx*3*nsp*nbas
                  auso_zv(:)=0.0d0
                  call makusq ( 0 , ssite , sspec ,  nbas , nbas , 
     .              0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .              nspc , jsp , 1 , qp , t_zv , ppnl_rv , auso_zv )
                  call mkorbm ( ssite , sspec , jsp , nsp , nspc , nlmax , ndham
     .              , nev , rv_a_owtkb , iq , nbas , ppnl_rv , auso_zv , nl , nkp 
     .              , orbtm_rv )
                  deallocate(auso_zv)
               endif

c$$$!! Mulliken analysis was here
c$$$               if (lwtkb .ne. -1) then
c$$$                  if (mullwrite) then
c$$$                        call rxx(nev.ne.ndimhx,'Mulliken requires all eigenvectors')
c$$$                        call rxx(nevmx.ne.ndhamx,'Mulliken requires nevmx=nhdamx')
c$$$                        allocate(doswt_rv(nchan*nevmx*nspc))
c$$$                        doswt_rv(:)=0.0d0
c$$$                        call mullmf(nbas , ssite , sspec , ham_iv_a_oindxo , t_zv , ndimh 
c$$$     .                    , nspc , iq , jsp , moddos , nsites , lsites , lmxch , nchan 
c$$$     .                    , chan_iv , lmdim , ndham , doswt_rv )
c$$$                        call iomomn(.true.,2,.false.,1,nspc,1,1,i)
c$$$                        i1 = iomoms ( - nfilem , nl , nsp , nspc , nkp , ndimh , i , 
c$$$     .                    nspc , 1 , 1 , ndhamx , nevmx , nchan , nchan , nev , evl ( 1 
c$$$     .                    , jsp ) , 0d0 , doswt_rv , 0d0 , 0d0 , 0d0 )
c$$$                        if (allocated(doswt_rv)) deallocate(doswt_rv)
c$$$                  endif
c$$$               endif   

!! Make new density matrix dmatu for LDA+U (implementation of LDA+U is a little complicated).
               if (lwtkb .ne. -1) then
                  if (nlibu .gt. 0 .and. nev .gt. 0) then
                     if (lwtkb .eq. 0)
     .                 call rx('metal weights required for LDA+U calculation')
c                     nl = int(ctrl_nl)
c                     nlmax = nl*nl
                     allocate(ausu_rv(2*nlmax*ndhamx*3*nsp*nbas))
                     if(debug) write(*,*)'allocate ausu_rv',2*nlmax*ndhamx*3*nsp*nbas
                     ausu_rv(:)=0.0d0
                     call makusq ( 0 , ssite , sspec ,   nbas , nbas , !sham ,slat ,
     .                 0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .                 nspc , jsp , 1 , qp , t_zv , ppnl_rv , ausu_rv )
                     call mkdmtu ( ssite , sspec , rv_a_owtkb , jsp , iq , nsp , nspc 
     .                 , ndham , nlmax , nbas , nev , ppnl_rv , ausu_rv , dmatu , nlibu 
     .                 , lmaxu , lldau )
                     if (allocated(ausu_rv)) deallocate(ausu_rv)
                  endif
               endif
!! Core-level spectroscopy
               if(cmdopt0('--cls') .and. lwtkb .ne. -1) then
                  call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
#if MPI
                  call rx('CLS only k-parallel')
#endif
                  call m_clsmode_set1(napw,ndimh,ndimhx,nmx,jsp,iq,qp,nev,igv2x,t_zv) !all inputs
               elseif (lwtkb.ne.-1 .and. (lrout.ne.0 .or. lmet.eq.4)) then
!! Density mode. wtkb is the weight for iq,ib. With this setting, 
!! we can plot |\psi_{iq,ib}(\bfr)|**2 in the manner of electron density
!! lmf --density,iq=12,ib=3,4,5   (here ib=3,4,5 are superposed).
                 if(cmdopt('--density',9,0,strn)) then
                   iqread=0
                   iqindex=index(strn(10:),'iq=')+2
                   if(iqindex==2) then 
                     iqread=-999
                   else
                     read(strn(10+iqindex:),*) iqread
                   endif
                   ibindex=index(strn(10:),'ib=')+2
                   if(ibindex==2) then
                     ibread=-999
                   else
                     do nibread=1,100
                       read(strn(10+ibindex:),*,err=2019,end=2019) ibread(1:nibread)
c                       print *,'xxx ibread=',ibread(1:nibread)
                     enddo
 2019                continue
                     nibread=nibread-1
                   endif  
                   if(iqread>0) then
                     if(maxval(ibread(1:nibread))<=0) call rx('--density mode: wrong ib=foobar. Try,e.g. --density,iq=12,ib=5')
                     write(6,"('--density bandmode: psi**2 for iq=',i5,' ib=',255i5)") iqread,ibread(1:nibread)
                     rv_a_owtkb(:,:,:)=0d0
                     do ib=1,nibread
                       rv_a_owtkb(ibread(ib),:,iqread)=1d0
                     enddo
                     do ib=1,ndimh
                       if(abs(rv_a_owtkb(ib,jsp,iq))>1d-3) write(6,"('ib wtkb=',i5,2f13.6)") ib,rv_a_owtkb(ib,jsp,iq)
                     enddo
                   endif
                 endif   
!! Accumulate output density and sampling DOS
!! Even if no output rho, still call addrbl to make DOS when lmet=4
                 if(debug) print *,'goto addrbl nevl=',nevl,lwtkb
                 if(lwtkb>=0) then
                   call addrbl ( ssite , sspec ,  jsp , nsp , nspc !sham ,slat , 
     .              , qp , rv_a_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout
     .              , lwtkb , rv_a_owtkb , lswtk , rv_a_oswtk , iq , lfrce , ldos
     .              , lekkl , k1 , k2 , k3 , zv_a_osmpot , vconst , sv_p_osig , sv_p_otau
     .              , sv_p_oppi , lcplxp , numq , qval - qbg , t_zv , evl , nev
     .              , ef0 , def ,  emin , emax , ndos , dos_rv , srout_zv  !esmear ,
     .              , sumqv , sumev , sv_p_oqkkl , sv_p_oeqkkl , frc )
                 endif     
               endif
            endif
            

!!   === Save evals for this qp ===
c            elseif (plbnd .ne. 0) then
c$$$#ifndef MPIK
c$$$               if(debug) write(*,*)  'plbnd/=0 (noMPI)', ndham,jsp 
c$$$               if (nfbn(1)>0 .and. n_listwf>0) then
c$$$                  call rsibl_ev( ssite,sspec,slat,nbas,isp,qp,iq,ndimh,nspc
c$$$     .              ,  napw,igv2x,ham_iv_a_oindxo,nev,t_zv,k1,k2,k3
c$$$     .              ,  n_listwf,listwf    )
c$$$               endif
c$$$               if(debug) print *,'dddddd111'
c$$$               i = nsp
c$$$               if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
c$$$c              print *,' eee1',nevl,jsp,i
c$$$c              print *,' eee2',allocated(evl)
c$$$c              print *,' eee3',evl(1:2,jsp)
c$$$c              print *,' eee3',evl(1:nevl,jsp)
c$$$               call suqlsw2(nevl,jsp,i,evl(1,jsp),qp) !takao
c$$$               if(debug) print *,'ddddd2222'
c$$$               if (nfbn(1) .ne. 0) then
c$$$                  if (ndimhx .ne. nevl)
c$$$     .              call rx('color weights not implemented when '//
c$$$     .              'nevl < hamiltonian dimension')
c$$$                  call suqlse ( ndimhx , jsp , i , ndimhx , 1 , nfbn , ifbls_iv 
c$$$     .              , ndhamx , t_zv , evl )
c$$$               endif
c$$$               if(debug) print *,'ddddd333'
c$$$               if (nfbn(2) .ne. 0) then
c$$$                  if (ndimhx .ne. nevl)
c$$$     .              call rx('color weights not implemented when '//
c$$$     .              'nevl < hamiltonian dimension')
c$$$                  call suqlse ( ndimhx , jsp , i , ndimhx , 2 , nfbn , ifbls_iv 
c$$$     .              , ndhamx , t_zv , evl )
c$$$               endif
c$$$               if(debug) print *,'ddddd4444'
c$$$#else
c               if (procid.eq.master) then
c               write(6,*)  'case plbnd.ne.0 (MPIK),procid.iq=',procid,iq
c               if(debug) write(*,*)  'nnn ', ndham,jsp 
c               endif
c              evlall(1:ndhamx,jsp,iq)= evl(1:ndhamx,jsp)
c$$$#endif
c            endif


!! Decompose to alm components. PROCAR generation.
!! ??? no nsites is used? ! Compare following codes with those of mkpdos.F (PDOS).
            if(PROCARon) then
c               nlmax = nl*nl
               allocate( auspp(nlmax,ndhamx,3,nsp,nbas) )
               initbd=.false.
               qold=qp
               auspp = 0d0
               call makusq( 1 , ssite , sspec ,   nbas , nbas , !sham ,slat ,
     .           0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .           nspc , jsp , 1 , qp , t_zv , ppnl_rv , auspp )
               if(spinoffdiag) then
                  jspini=1
                  jspend=2
               else
                  jspini=isp
                  jspend=isp
               endif
               do jspp= jspini, jspend
                  if(jspp==1) iprocar=iprocar1
                  if(jspp==2) iprocar=iprocar2
                  if(debug) print *,'iiiii =',iprocar1,iprocar2,jspp,iprocar,ef0
                  if(debug) print *,'nnnnn=',nlmax,ndham,nspc,nsp,nbas
                  ccc="ion        s       py       pz       px      dxy      dyz      dz2      dxz   dx2-y2"//
     &                 "      f-3      f-2      f-1       f0       f1       f2       f3"//
     &                 "                                                                                 tot"
                  write(iprocar,*)
                  write(iprocar,*)
                  write(iprocar,"('k-point ',i4,' :    ',3f11.8,'     weight = -------  : x =',f15.8)")iq,qp,xdatt(iq)
                  write(iprocar,*)
               do iband = 1, nev
                  write(iprocar,*)
                  write(iprocar,"('band ',i3,' # energy ',f13.8,' # occ. -----' )")iband,(evl(iband,jsp)-ef0)*rydberg
                  write(iprocar,*)
                  dwgtt=0d0
                  do ib = 1, nbas
                     is  = ssite(ib)%spec
                     ilm = 0
                     dwgt=0d0
                     do  l = 0, sspec(is)%lmxa
                        do  m = -l, l
                           ilm = ilm+1
                           au = auspp(ilm,iband,1,jspp,ib)
                           as = auspp(ilm,iband,2,jspp,ib)
                           az = auspp(ilm,iband,3,jspp,ib)
                           s11 = dconjg(au)*au*ppnl_rv(2,l+1,jspp,ib)
                           s22 = dconjg(as)*as*ppnl_rv(7,l+1,jspp,ib)
                           s33 = dconjg(az)*az*ppnl_rv(8,l+1,jspp,ib)
                           s12 = 0
                           s13 = 2*dconjg(au)*az*ppnl_rv(9,l+1,jspp,ib)
                           s23 = 2*dconjg(as)*az*ppnl_rv(10,l+1,jspp,ib)
                           dwgt(ilm)= s11+s22+s33 + s12+s13+s23
                        enddo
                     enddo
                     dwgtt = dwgtt + dwgt(1:ilm)
                     if(ib==1)  write(iprocar,"(a)") trim(ccc)
                     write(iprocar,"(i3,100(x,f8.5))")ib,(dwgt(i),i=1,nchanp),sum(dwgt)
                     if(ib==nbas) write(iprocar,"('tot',100(x,f8.5))")(dwgtt(i),i=1,nchanp),sum(dwgtt)
                     if(fullmesh) dwgtall(1:nchanp,ib,iband,jspp,iq) = dwgt(1:nchanp)
                  enddo
               enddo
               enddo
               if(debug) print *,'mmmmmmmm procid sum dwgt check=',procid,iq,jsp,sum(dwgtall),nev
               deallocate( auspp )
            endif 
            if(allocated(t_zv)) deallocate(t_zv)
!! == end loop over isp (main loop in parallel mode)==
            if(debug) print *,'bbbbbbbb before 30'
c 30         continue
c            deallocate(igv2x)
 2005    continue
         if(allocated(hamm)) deallocate(hamm,ovlm)
         ndimhx_(iq)=ndimhx
         nev_(iq)   =nev
 2010 continue
      if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid,ix99
!!----------------------
!! fixed at sep2020. nev_ and ndimhx_ are not correctly passed via MPI      
      if(numprocs/=1.and.ix99==1) then
         do iq = 1, nkp
            pidorigin=-999
            do idx=0,numprocs-1
               if (kpproc(idx)<=iq .and. iq<= kpproc(idx+1)-1) then
                  pidorigin=idx
                  exit
               endif
            enddo   
            call MPI_BCAST(ndimhx_(iq),1,MPI_INTEGER, pidorigin,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(nev_(iq),1,MPI_INTEGER,    pidorigin,MPI_COMM_WORLD,ierr)
         enddo   
      endif
      nevmin=minval(nev_(1:nkp))
!! ========== end of iq loop ====================================================

      
!! === pdos mode (--mkprocar and --fullmesh). ===
      if(debug.and.PROCARon) print *,'mmmm procid sum dwgt check=',procid,sum(dwgtall)
      if(PROCARon.and.fullmesh) then
#ifdef MPIK
         call mpi_barrier(MPI_comm_world,ierr)  
         call xmpbnd(kpproc,ndhamx,nkp,nspx,evlall) !ndhamx,nspx nov2014 bugfix
         call xmpbnd(kpproc,nbas*nchanp*ndhamx,nkp,nspx,dwgtall) 
#endif
         if(master_mpi) then  
            if(debug) print *,'mmmmm tot procid sum dwgt check=',procid,sum(dwgtall)
!! tetrahedron without symmetry
            allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            iq=0
            do ik3 = 1, nkk3 
               do ik2 = 1, nkk2
                  do ik1 = 1, nkk1
                     iq = iq+1
                     ipqe(ik1,ik2,ik3)=iq
                     qx = (ik1-1)*qlat(:,1)/nkk1 + (ik2-1)*qlat(:,2)/nkk2+ (ik3-1)*qlat(:,3)/nkk3
                     if(abs(sum(qx-qplist(:,iq)))>1d-6) then !sanity check
                        call rx("bndfp: qx/=qplist something strang")
                     endif  
                  enddo
               enddo
            enddo
            call tetirr(qlat, nkk1,nkk2,nkk3, ipqe, ntete,idtete) !tetrahedron
            write(*,*)" ntete 6*nkk1*nkk2*nkk3 nkp=",ntete,6*nkk1*nkk2*nkk3,nkp
!! doswt(1:nchanp=16,ndhamx,nsp,nkp)
!! dwgtall(1:16,nbas,ndhamx,nsp,nkp)
!!   If we will have a memory problem for dwgtall, we may store dwgtall(i,...) for a file.
!!   and repeat nchanp times.
c            eminp=-25.0/rydberg
c            emaxp=30.0/rydberg
c            ndos=5500
            write(6,"(' pdosmode: ndhamx nsp nspx =',4i7)") ndhamx, nsp, nspx,nevmin !nevmx
c            print *,' pdosmode: nchanp nbas ndos emin emax ef0=', nchanp,nbas,ndos, emin,emax,ef0
            write(6,"(' pdosmode: nchanp nbas emin emax ef0=',2i7,3f12.5)") nchanp,nbas,emin,emax,ef0 !mar2015
            open(newunit=ifip,file='pdosdata.'//trim(sname),form='unformatted')
            write(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp !ndos,nkp mar2015
            write(ifip) idtete,ipqe ! allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            write(ifip) evlall  !  allocate(evlall(ndhamx,nspx,nkp))
            write(ifip) dwgtall !  allocate(dwgtall(nchanp,nbas,ndhamx,nspx,nkp))
            write(ifip) ef0 !eminp,emaxp,ef0
            close(ifip)
            if(nkp/=nkk1*nkk2*nkk3) call rx('pdosmode but nkp/=nkk1*nkk2*nkk3')
            if( cmdopt('--tetraw',8,0,strn)) then
               open(newunit=ifip,file='tetradata.dat',form='unformatted')
               write(ifip) ndhamx,nkp,ntete !ndos,nkp mar2015
               write(ifip) idtete
               close(ifip)
               open(newunit=ifip,file='qplistf.dat')
               do iq=1,nkp
                  write(ifip,"(3f15.8)") qplist(:,iq)
               enddo   
               close(ifip)
            endif               !! we need evlall(ndhamx,nkp) for qplist
            deallocate(idtete,ipqe) !remove this if you uncomment above lines.
         endif  
         if(writeham) call rx0('Done --writeham: --mkprocar --fullmesh are needed now. HamiltonianMTO* genereted')
         call rx0('Done pdos: --mkprocar & --fullmesh. Check by "grep k-point PROCAR.*.*"')
      endif

!! for MPI, collect data of nodes see xmpbnd.
#if MPIK
      entime = MPI_WTIME()
      call info2(30,0,0, ' ... Done MPI k-loop: %;1d seconds elapsed',entime-sttime,0)
      call info0(20,0,-1,' ... Sharing data between processes...')
      sttime = MPI_WTIME()
      call xmpbnd(kpproc,ndhamx,nkp,nspx,evlall) !ndhamx,nspx nov2014 bugfix
      if (cmdopt0('--cls') .and. lwtkb .ne. -1) then
         call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
      endif
      call mpibc2(nevls,nkp*nspx,2,mlog,'bndfp','nevls')
#endif
!!
!! boltztrap data 
      if( cmdopt('--boltztrap',11,0,strn) .and. master_mpi) then
         iqread=0
         iqindex = index(strn(12:),'nb=')+2
         nbandx = ndhamx
         if(iqindex/=2) then 
            read(strn(12+iqindex:),*) nbandx
            nbandx=min(nbandx,ndhamx)
         endif
         open(newunit=ifi,file='efermi.lmf') !readin fermi energy from efermi.lmf
         read(ifi,*)  eferm
         close(ifi)
        open(newunit=ifi,file=trim(sname)//'.intrans_template.boltztrap')
        write(ifi,"(a)")'GENE          # format '
        write(ifi,"(a)")'0 0 0 0       # iskip (not presently used) idebug setgap shiftgap'
        write(ifi,"(f20.16,a,f9.4,a)") eferm,' 0.0005 0.4 ', qval-qbg,
     &  '       # efermi.lmf (Ry), energygrid, energy span around Fermilevel, number of electrons'
        write(ifi,"(a)")'CALC          # CALC (calculate expansion coeff), NOCALC read from file'
        write(ifi,"(a)")'5             # lpfac, number of latt-points per k-point'
        write(ifi,"(a)")'BOLTZ         # run mode (only BOLTZ is supported)'
        write(ifi,"(a)")'0.15          # (efcut) energy range of chemical potential'
        write(ifi,"(a)")'800.0 50.0  # Tmax, temperature grid'
        write(ifi,"(a)")'-1.0  # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)'
        write(ifi,"(a)")'TETRA'
        close(ifi)
        open(newunit=ifi,file=trim(sname)//'.struct.boltztrap')
        write(ifi,"(a)") trim(sname)
        write(ifi,"(3d16.8, ' # plat1  ')") plat(:,1)*alat
        write(ifi,"(3d16.8, ' # plat2  ')") plat(:,2)*alat
        write(ifi,"(3d16.8, ' # plat3  ')") plat(:,3)*alat
        write(ifi,*) ngrp
        do ig=1,ngrp
           symxx = matmul(transpose(qlat), matmul(symops(:,:,ig),plat))
           if(abs(sum( nint(symxx(:,:))-symxx(:,:) )) >1d-6) call rx('bndfp: sym not integers --- probably strange')
           write(ifi,"(9i3)") ((nint(symxx(i,j)),j=1,3),i=1,3)
        enddo
        close(ifi)   
        do jsp=1,nspx         !=isp.  nspx=1 for so=1
           open(newunit=ifi,file=trim(sname)//'.energy.isp'//trim(i2char(jsp))//'.boltztrap')
           write(ifi,"(a)") trim(sname)
           write(ifi,"(i10)") nkp
        do iq=1,nkp
!! true q= 2pi/alat *qplist(:,iq) in cartesian
           qvec= matmul(transpose(plat),qplist(:,iq)) !qvec in qlat unit
           ncount=0
           do iband=1,nbandx
              if(evlall(iband,jsp,iq)>1d98) cycle
              ncount=ncount+1
           enddo
           write(ifi,"(3f15.8,i10)") qvec,ncount
           do iband=1,nbandx
              if(evlall(iband,jsp,iq)>1d98) cycle !evlall=1d99 is dummy
              write(ifi,"(d23.16)") evlall(iband,jsp,iq)
           enddo   
        enddo
           close(ifi)
        enddo
      endif
      if( cmdopt('--boltztrap',11,0,strn)) then
         call rx0(' end of --boltztrap mode: boltztrap.* are generated')
      endif
      
      
!!  Allreduce density-related quantities
      if(debug) print *,'goto all reduce'
      if (lrout .ne. 0) then
         call mpibc2(sumqv,6,4,mlog,'bndfp','sumqv')
         call mpibc2(sumev,6,4,mlog,'bndfp','sumev')
         call mpibc2 ( srout_zv , k1 * k2 * k3 * nsp * numq , 6 , mlog, 'bndfp' , 'smrho' )
         if (lswtk .eq. 1) then
            call mpibc2 ( rv_a_oswtk , ndhamx * nkp , 4 , mlog , 'bndfp' , 'swtk' )
         endif
C     Allreduce qkkl
c         nkaph = globalvariables%nkaph
         do  ib = 1, nbas
            is = int(ssite(ib)%spec)
            lmxa=sspec(is)%lmxa
            lmxh=sspec(is)%lmxb
            kmax=sspec(is)%kmxt
            nlma = (lmxa+1)**2
            nlmh = (lmxh+1)**2
            if (lmxa .gt. -1) then
               nelt(1) = (kmax+1)*(kmax+1)*nlma*nlma
               nelt(2) = (kmax+1)*nkaph*nlma*nlmh
               nelt(3) = nkaph*nkaph*nlmh*nlmh
               do  i = 1, 3
                  call mpibc2 ( sv_p_oqkkl( i , ib )%v , nelt ( i ) * numq * nsp 
     .              * nspc , 4 , mlog , 'bndfp' , 'qkkl' )
                  if ( lekkl .eq. 1 ) call mpibc2 ( sv_p_oeqkkl( i , ib )%v , nelt 
     .              ( i ) * numq * nsp * nspc , 4 , mlog , 'bndfp' , 'eqkkl' )
               enddo
            endif
         enddo
!!    Allreduce DOS, forces, dmatu, orbtm
         if ( ndos .gt. 0 ) call mpibc2 ( dos_rv , ndos * 2 * nsp , 4 , mlog , 'bndfp' , 'dos' )
         if (lfrce .ne. 0)  call mpibc2(frc,3*nbas*numq,4,mlog,'bndfp','frc')
         if (nlibu .gt. 0)  call mpibc2(dmatu,nsp*nlibu*(lmaxu*2+1)**2,6,mlog,'bndfp','dmatu')
         if ( lso .ne. 0 .and. lwtkb .ne. - 1 ) call mpibc2 ( orbtm_rv 
     .     , nl * nsp * nbas , 4 , mlog , 'bndfp' , 'orbtm' )
      endif
!! eterms are alwasy broadcasted. comment out if conditions jan2012
      eterms=ham_eterms ! call dcopy(size(ham_eterms),ham_eterms,1,eterms,1)
      call mpibc2(eterms(19),1,4,mlog,'bndfp','rhosig')
      ham_eterms=eterms ! call dcopy(size(ham_eterms),eterms,1,ham_eterms,1)
#ifdef MPIK      
      entime = MPI_WTIME()
      call info2(20,0,0,' MPI broadcast took %;1d sec', (entime-sttime),0)
#endif

!! ===  Write bands in bands-plotting case: loop over qp getting evals from array ===
      if(plbnd .ne. 0) then
         if (master_mpi) then
            print *,'procid=',procid
!! fermi surface mode together with band data around Ef (Okumura, 2019/11/14)
            if(fsmode) call writefs(ndhamx,nspx,nkp,evlall,qplist,ef0)
!! write bands directly. All inputs
            write(stdo,*)' Writing bands to bands file ...'
            if(nsyml/=0) then
              call writeband(stdo,nkp,nsyml,nspx,nsp,ix99,ndhamx, !,fsmom/=NULLR,vnow
     i        xdatt,nqp_syml,nqp2n_syml,nevls,evlall,qplist,labeli,labele,alat,eferm,evtop,ecbot,
     i        sname,nqps_syml,nqpe_syml,nsyml,dqsyml,
     i        etolv,etolc)
            endif
         endif
         if(fsmom/=NULLR) then
         write(6,"(a)")'NOTE: Be carful! fsmom (fixed moment on). In sc cycle, we use additional bias mag. field  '
         write(6,"(a)")'NOTE: Use mag.field is written in MagField. But it is not used for --band mode! mar2015.'
         endif
         if(fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
         call rx0('plot band mode done')
      endif
      
!! ===   Repeat loop for printout and goto 99 ===
!!   jsp=isp in the collinear case; jsp=1 in the noncollinear
!!     Thus jsp should be used in place of isp
!!     isp serves as a flag for the noncollinear case
      goto99=.false.
      do  iq = 1, nkp
         qp=qplist(:,iq)
         do isp = 1, nsp
            ispc = min(isp,nspc)
            if( .not. ((onesp.eq.0 .or. isp.eq.onesp) .and. (ispc.eq.nspc)) ) cycle
            jsp = isp
            if (ispc .eq. 2) jsp = 1
            if(iprint()>20.and.ix99==1) then
               write(*,'(" bndfp: kpt",i5," of",i5, " k isp=",3f8.4,i2," ndimh nev=",2i5)')
     &           iq,nkp,qp,jsp,ndimhx_(iq),nevls(iq,jsp)
            endif   
            if(iprint()>=35.and.ix99==1) then
               write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
            endif
            evl(1:ndhamx,jsp) = evlall(1:ndhamx,jsp,iq)
            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
            if(debug) call info5(30,0,0,' bndfp 299loop:  kpt %i of %i, k=%3:2,5;5d', iq,nkp,qp,0,0)
            if(debug) print *,'eeeee11111111 ',procid,numprocs
            if(numprocs==1) then
               pidorigin = master
            else
               pidorigin=-999
               do idx=0,numprocs-1
                  if (kpproc(idx)<=iq .and. iq<= kpproc(idx+1)-1) then
                     pidorigin=idx
                     exit
                  endif
               enddo   
            endif
            nev_iq    = nev_(iq) 
            if (mod(iq,10) .ne. 1) call poppr
            ebot = dmin1(ebot,evl(1,jsp))
            i = max(1,nint(qval-qbg)/(3-nspc))
            evtop = max(evtop,evl(i,jsp))
            ecbot = min(ecbot,evl(i+1,jsp))
            if (lmet .eq. 0 .and. iq .eq. 1 .and. jsp .eq. 1) ef0 = evtop
            if(debug) print *,'eeeee44444444444 plbnd=',plbnd
            
            if (plbnd .eq. 0) then
               if (ipr.ge.10 .and. iq.eq.1 .and. ipr.gt.0)
     .           write (stdl,712) (evl(i,jsp),i=1,nev_iq)
 712           format('fp evl',8f8.4)
               if (lwtkb .ne. -1) then ! .and. .not. lwndow) then
                  if (iq .eq. 1 .and. jsp .eq. nsp ) then !
!! Energy-Window reset. dosw=(emin,emax)
                     ef00 = ef0
                     call fixef0(qval-qbg,jsp,1,nev_iq,ndhamx,evl,dosw,ef0)
                     if(master_mpi) then
                        open(newunit=ifid,file='ewindow.'//trim(sname)) 
                        write(ifid,"(3d23.15)") dosw(1:2),ef0
                        close(ifid) !call fclose(ifid)
                     endif
                     if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .                    lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
                        if (master_mpi) call info0(10,1,1,
     .                       ' ... Fermi level reset in second spin'//
     .                       ' channel ... restart band pass')
                        goto99=.true.
                        exit    !this is the case of make co test at ecalj/TestInstall/
                     endif
                  endif

C     Check for cases when nevmx is too small : i=2 => fatal error
                  i = 0
                  if (nevmx.ge.0 .and. lmet .ne. 0) then
                     dum = evl(max(nev_iq,1),jsp)
                     if (.not. ltet .and. ef0+5*bz_w .gt. dum) i=2
                     if (lmet.eq.4 .and. ef0+def+5*bz_w .gt.dum)i=2
                  endif
                  if(i .eq. 2) then
                     write(6,"(a,f13.5,f13.5)")
     &                 'evl(nev='//trim(i2char(nev_iq))//')=',
     &                 evl(max(nev_iq,1),jsp),' but ef0=',ef0
                     call rx('bndfp:... restart with larger efmax or nevmx: bndfp')
                  endif
               endif
            endif
         enddo                  ! end second loop over isp
      enddo                     !end second loop over iq

#ifdef MPIK
c     if (allocated(kpproc)) deallocate(kpproc, stat=ierr)
      call mpi_barrier(MPI_comm_world,ierr) !takao added I think 'safer' by this barrier. jan2015
#endif
      print *,'goto99=',goto99
      if(goto99) goto 99


!!   For now, SO by site, not by class
      if (lwtkb .eq. 1 .and. lso .ne. 0) then
         allocate(ips_iv(nbas))
         do i_spackv=1,nbas
c            call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %spec 
c     .        , 1 , i_spackv + 1 - 1 , ips_iv )
            ips_iv(i_spackv)=ssite ( i_spackv ) %spec 
         enddo
!! orbital moment
         call iorbtm ( sspec , ips_iv , nl , nl , nbas , nsp , orbtm_rv  )
         if (allocated(ips_iv)) deallocate(ips_iv)
      endif
      if (ipr.gt.0) write (stdl,715) nkp,ebot,qval,qbg,bz_w
 715  format('nv nkp',i5,'  ebot',f9.4,'   qval',f10.4,'  qbg',f8.4,'  esmear',f8.4)
      
!! --- Interpolate density to Fermi energy ---
      sev = sumev(1,1)
      if (lmet .eq. 4) then
         call mshn3p ( nbas , ssite , sspec , lmet , lrout , lfrce , qval 
     .     - qbg , ef0 , def , sumqv , sumev , n1 , n2 , n3 , k1 , k2 , 
     .     k3 , srout_zv , sv_p_oqkkl , frc , lrep )
C   ... Store val q & magnetic moment in sumqv(1) and sumqv(2)
         sumqv(2,1) = sumqv(1,1) - sumqv(1,2)
         sumqv(1,1) = sumqv(1,1) + sumqv(1,2)
C   ... Eigenvalue sum including entropy term
         sev = sumev(1,1) + sumev(2,1)
C   ... Remake sev,ef linearly interpolating tabulated sampling DOS
         sev00 = sev
         ef00  = ef0
         if (ldos .ne. 0) then
            call efldos ( qval , nsp , emin , emax , ndos , dos_rv , eferm, sev1 )
            sev   = sev1
            ef0  = eferm
         endif
         if (ipr .gt. 30 .and. ldos .ne. 0)
     .     write(stdo,388) sev00,ef00,sev1,eferm,sev,ef0
 388     format(' ipol:  sev=',f12.6,'   ef=',f12.6:/
     .     ' dos:   sev=',f12.6,'   ef=',f12.6/
     .     ' use:   sev=',f12.6,'   ef=',f12.6)
         if (ipr.gt.0) write (stdl,733) ef00,eferm,ef0,sev00,sev1,sev
 733     format('nf EF:',3f9.5,'    EB:',3f12.5)
c         bz_ndos=ndos
         bz_dosw = dosw
         bz_ef=ef0
         bz_def=def
         if (lrep .eq. 1) then
            ef0 = -1
           write(stdo,*)'Input Fermi energy was too far off, repeat band pass'
#ifdef MPIK
            call mpi_barrier(MPI_comm_world,ierr)  
#endif
            goto 99
         endif
      endif
#if MPIK
      call mpi_barrier(MPI_comm_world,ierr)
#endif
      if (master_mpi) then
         ierr=unlink('MagField') !delete
      endif  
!! --- BZ integration for fermi level, band sum and qp weights ---
      if(lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
c$$$         if (lwndow) then
c$$$            allocate(ww_rv(ndham*nsp*nkp))
c$$$            eferm = min(dosw(1),dosw(2))
c$$$            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
c$$$     .        , ww_rv , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , eferm 
c$$$     .        , 2 , ntet , iv_a_oidtet , sev , dum )
c$$$            eferm = max(dosw(1),dosw(2))
c$$$            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
c$$$     .        , rv_a_owtkb , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , 
c$$$     .        eferm , 2 , ntet , iv_a_oidtet , sev , dum )
c$$$            call daxpy ( ndham * nsp * nkp , - 1d0 , ww_rv , 1 , rv_a_owtkb, 1 )
c$$$            if (allocated(ww_rv)) deallocate(ww_rv)
c$$$         else
         dosrng = 8
         if (bz_n .lt. 0) dosrng = 16
         if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'
         if(bz_fsmommethod == 1) then !takao dec2010
!     ! vnow june22013 !vnow !(in Ry) contains magnetic field
!     ! For eigenvalus, add  -vnow/2 for isp=1, and +vnow/2 for isp=2.
            call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .            , dosrng , rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill,vnow)
         else
!     ! vnow june22013
            call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .           , dosrng , rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill, vnow)
         endif
!     ! june2013 magfield is added
         if(fsmom/=NULLR.and.master_mpi) then
            ifimag=ifile_handle()
            open(ifimag,file='MagField',status='unknown')
            write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow
            close(ifimag)
         endif
!     !         Store val charge & magnetic moment in sumqv(1..2)
         if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
         endif
c$$$         endif
!!
         if (lmet .ne. 4) then
            ef0 = eferm
            bz_ef=ef0
         endif
         if (lmet .gt. 0) then
            if (master_mpi) then
               ifi=ifile_handle()
               open(ifi,file='wkp.'//trim(sname),form='unformatted') 
               i = iobzwt ( 0 , ndhamx , nkp , nspx , eferm , rv_a_owtkb , -ifi )
               close(ifi) !call fclr('wkp',ifi)
            endif
         endif
         goto99=.false.
         if (lwtkb .eq. -1 .and. lrout .gt. 0) then
            call info0(20,0,0,' Start second band pass ...')
            lwtkb = 1
c            if (nspc .eq. 2) lswtk = 1 !commnet out nov2015  right?
            goto99=.true.
         endif
         if (lwtkb .eq. 2 .and. lrout .gt. 0) then
            call info0(20,0,0,' New pass with constrained weights ...')
            goto99=.true.
         endif
#ifdef MPIK
         call mpi_barrier(MPI_comm_world,ierr)  
#endif
c$$$c jan2015 did not work for Testinstall/te
c$$$  call poppr
         if(goto99) goto 99
      endif
c      if(allocated(sfz)) deallocate(sfz)
!! == goto99 ix99 loop over =============================
      
!! ... Save Fermi level, nonmetal or sampling integration
      if (lmet .eq. 0 .or. .not. (lmet .ne. 4 .or. ltet)) then
         if (lmet .eq. 0) then
            ef0 = (evtop+ecbot)/2
            write(6,"(' Highest occ. level = ',f13.5,
     &        ' Lowest unocc. = ',f13.5, ' average =',f13.5)") evtop,ecbot,ef0
         endif
         if (master_mpi) then
            ifi=ifile_handle()
            open(ifi,file='wkp.'//trim(sname),form='unformatted')
            i = iobzwt ( 1 , ndham , nkp , nsp , ef0 , rv_a_owtkb , -ifi )
            close(ifi) !call fclr('wkp',ifi)
         endif
      endif
c$$$ -->!! Here ... Cleanup asa-style moments file, print table of DOS channels

!! ... Average forces so net force on system is zero (APW case)
      if (lfrce .ne. 0 .and. napw .ne. 0) then
         call dpzero(xv,3)
         do  i1 = 1, nbas
            do  i = 1, 3
               xv(i) = xv(i) + frc(i,i1)/nbas
            enddo
         enddo
         do  i1 = 1, nbas
            do  i = 1, 3
               frc(i,i1) = frc(i,i1) - xv(i)
            enddo
         enddo
      endif

!! === Generate DOS on disk ===
c      if(.not. lwndow .and. master_mpi) then
      if( master_mpi) then
         ndos0 = ndos   !it was bz_ndos. I think this is fine
                        !ndos0 is only used as a swith wether DOS or integrated-DOS.
         if ( ldos/=0 ) then
            if(ndos0<0) write(6,"(a)")'... Generating integrated DOS' !, ndos0=,ndos0
            if(ndos0>0) write(6,"(a)")'... Generating total DOS'      !  ndos0=,ndos0
!! Default and external option. 
            if(tdos) then
              dosmax = bz_dosmax
            endif
            if(cmdopt('--tdostetf',10,0,strn)) ltet=.false. !oct2015
            if (ltet) then
!! determine eferm
               call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .           , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw (
     .           2 ) , dos_rv , ndos , eferm , 1 , ntet , iv_a_oidtet , dum , dum )
!! repeat again for given range emax=dosw(2). mar2012
               if(tdos) dosw(2) = dosmax + eferm 
               call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .           , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw (
     .           2 ) , dos_rv , ndos , eferm , 1 , ntet , iv_a_oidtet , dum , dum )!job=1 integrated dos
               if ( ndos0 > 0 ) call xxxdif ( dosw ( 1 ) , dosw ( 2 ) , ndos 
     .           , nspx , 0 , dos_rv ) !derivative 'integrated dos -> DOS'
            else
               if(tdos) dosw(2) = dosmax + eferm  !.or.pdos
c               if ( bz_n .ge. 100) bz_n = mod(bz_n,100)
               if ( ndos0 > 0 ) call makdos ( nkp , ndhamx , ndhamx , nspx
     .           , rv_a_owtkp , evlall , bz_n , bz_w
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
               if ( ndos0 < 0 ) call maknos ( nkp , ndhamx , ndhamx , nspx
     .           , rv_a_owtkp , evlall , bz_n , bz_w 
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
            endif
            if ( nspc .eq. 2 ) call dscal ( ndos , .5d0 , dos_rv , 1 )
            open(newunit=ifi,file='dos.tot.'//trim(sname))
            open(newunit=ifii,file='dosi.tot.'//trim(sname))
            dee=(dosw(2)-dosw(1))/(ndos-1d0)
            do ipts=1,ndos
               eee=dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-eferm
               dosi(1:nsp)= dosi(1:nsp) + dos_rv(ipts,isp)*dee
               write(ifi,"(255(f13.5,x))") eee,  (dos_rv(ipts,isp),isp=1,nsp)
               write(ifii,"(255(f13.5,x))") eee, (dosi(isp),isp=1,nsp)
            enddo
            close(ifi)
            close(ifii)
c            if (allocated(orbtm_rv)) deallocate(orbtm_rv)
c            if (allocated(dos_rv)) deallocate(dos_rv)
         endif
      endif
      if(tdos) call rx0('Done tdos mode:')

!! efermi.lmf 
      if(master_mpi.and.plbnd==0) then
         ifi= ifile_handle()     
         open(unit=ifi,file='efermi.lmf')
         write(ifi,"(d16.8, ' # (Ry) Fermi energy given by lmf')") eferm
         write(ifi,"(d16.8, ' # (Ry) Top of Valence')") evtop
         write(ifi,"(d16.8, ' # (Ry) Bottom of conduction')") ecbot
         write(ifi,"(d16.8, ' # number of electrons qval-qbg=')") qval-qbg
         write(ifi,"(d16.8, ' # mag. mom (nup-ndown)=')") sumqv(2,1) 
         write(ifi,
     &     "('# NOTE: These valuse are overwritten by lmf(-MPIK) in SC loop')")
         write(ifi,"('# NOTE: Not by band plot mode')")
         close(ifi)             
      endif

!! --- Core-level spectroscopy ---
      if (cmdopt0('--cls')) then
         if (master_mpi) then
            call m_clsmode_finalize(bz_ef,ndimh,ndhamx,nspx,nkp,dosw,evlall)
         else
            call fexit(0,0,' ',0)
         endif
      endif
!! --- Assemble output density ---
      do  ib = 1, nbas
        is =  ssite(ib)%spec
        nr =  sspec(is)%nr
        lmxl= sspec(is)%lmxl
        nlml = (lmxl+1)**2
        if (lmxl > -1) then
          allocate( sv_p_orhat1(1,ib)%v(nr*nlml*nsp))
          allocate( sv_p_orhat1(2,ib)%v(nr*nlml*nsp))
          allocate( sv_p_orhat1(3,ib)%v(nr*nsp))
        endif
      enddo
      allocate(qbyl_rv(n0,nsp,nbas))
      allocate(hbyl_rv(n0,nsp,nbas))
      if(debug) write(*,*) 'kino: allocate hbyl_rv,',n0*nbas*nsp
c$$$!!     --window: Put output density into orhat and smrho, and exit
c$$$      if (lwndow) then
c$$$         call mkrout ( ssite , sspec ,  nbas , nsp , ldim  !sham , slat ,
c$$$     .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_a_oorhat , hab_rv , sab_rv 
c$$$     .     , qbyl_rv , hbyl_rv , lrout )
c$$$         call zcopy( k1 * k2 * k3 * nsp , srout_zv , 1 , zv_a_osmrho 
c$$$     .     , 1 )
c$$$         call symrho(ssite , sspec ,  lfrce , zv_a_osmrho , sv_a_oorhat !slat ,
c$$$     .     , qbyl_rv , hbyl_rv , frc )
c$$$         goto 999
c$$$      endif

      call mkrout ( ssite , sspec ,  nbas , nsp , ldim  !sham ,slat , 
     .  , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_orhat1 , hab_rv , sab_rv 
     .  , qbyl_rv , hbyl_rv , lrout )
C
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
         call symrho ( ssite , sspec ,  lfrce , srout_zv , sv_p_orhat1 !slat ,
     .     , qbyl_rv , hbyl_rv , frc )
C   ... New boundary conditions pnu for phi and phidot
C       call pshpr(iprint()-12)
         if (lpnu .gt. 0) then
            call dcopy(size(ham_pmin),ham_pmin,1,xv,1)
            call dcopy(size(ham_pmax),ham_pmax,1,xv(11),1)
            call pnunew ( nbas , nsp , ssite , sspec , xv , xv ( 11 ) , lfrzw 
     .        , hab_rv , sab_rv , qbyl_rv , hbyl_rv, eferm )
         endif
C       call poppr
      endif
!! --- print out qbyl ---
      if(master_mpi) then
         ifqbyl = ifile_handle() 
         open(ifqbyl,file='qbyl.'//trim(sname))
         write(ifqbyl,"(i5)")nbas
         do  ib = 1, nbas
           ispec=ssite(ib)%spec
           write(ifqbyl,"(i4, i3,i2,10f12.6)")ib, ispec, sspec(ispec)%lmxa, (sum(qbyl_rv(il,1:nsp,ib)),il=1,sspec(ispec)%lmxa+1)
         enddo  
         close(ifqbyl)
      endif

!! --- Evaluate Harris energy ---
      call mkehkf(1,sev,sumqv(2,1),sumtv,ehar) !sham,
!! --- Evaluate KS total energy, correction to Harris force ---
      if (lrout .ne. 0) then
!!   ... Correction to harris force
        if(lfrce>0 ) then
          call dfrce (lfrce, ssite , sspec ,  k1 , k2 , k3 , nvl !sctrl , slat ,
     .     , sv_a_oorhat , sv_p_orhat1 , elind , qmom_rv , zv_a_osmrho , 
     .     srout_zv , fh_rv )
        endif
!!   ... Evaluate KS total energy and output magnetic moment
         eks = 0d0
         if (leks .ge. 1) then
            call togpr()
            call mkekin ( nbas , ldim , ssite , sspec ,  sv_p_osig  !sham ,slat , 
     .        , sv_p_otau , sv_p_oppi , lcplxp , sv_p_oqkkl , k1 , k2 , k3
     .        , vconst , zv_a_osmpot , srout_zv , sev , sumtv )
ccccccccccccccccccccccccc
ctakao test mkekin repeat
c         print *,'vvv takao mkekin repeat 000 --- oldpot x oldrho'
c         call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sevxxx,sumtvxxx)
c       stop 'xxxxxxxxxxxxxxxx  111aaa'
ccccccccccccccccccccccccc
ccccccccccccccccccccc
c          print *,'vvv: takao test osmpotbk 000 '
c          call defdc(osmpotbk,k1*k2*k3*2)
c          call dcopy(k1*k2*k3*4,w(osmpot),1,w(osmpotbk),1)
ccccccccccccccccccccc
c          call pshpr(ipr-4)

!! mkpot calls locpot. and locpot calls augmat.
!! augmat calculates sig,tau,ppi.
c          print *,'goto mkpot procid=',procid
            call mkpot (lfrce,
     &          srout_zv , sv_p_orhat1 , qbg , zv_a_osmpot 
     .        , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .        , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .        , vval_rv , fes2_rv , 0 , vorb ) !job=0 means energy term only

           if(cmdopt('--density',9,0,strn)) then
#ifdef MPIK
             call mpi_barrier(MPI_comm_world,ierr)  
#endif
             call rx0('end of --density mode')
           endif  
             
ccccccccccccccccccccccccccccccccccc
ctakao test mkekin repeat
c          print *,'vvv takao mkekin repeat 001x ---newpot x newrho'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),srout_zv,sev,sumtvxxx)
c         print *,'vvv takao mkekin repeat 001 --- newpot x oldrho'
c         call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sev,sumtvxxx)
c          print *,'vvv takao mkekin repeat 002 --- oldpot x oldrho'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpotbk),w(osmrho),sev,sumtvxxx)
c          call rlse(osmpotbk)
c          stop 'xxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccc
            call mkehkf(2,sev,sumqv(2,1),sumtv,eks) !sham,
            call togpr()
         endif
         
!! --- Add together force terms ---
!! fes1_rv: contribution to HF forces from estat + xc potential   This is for input  density  !=3rd term in (B.5) in JPSJ.84.034705
!! fes2_rv: contribution to KS forces from estat + xc potential   This is for output density
!! fh_rv      : 1st term  in (B.5)
!! input frc  : 2nd term in (B.5). Puley term
!! output frc  : totlal
         if(lfrce .gt. 0) call totfrc( nbas , ssite ,  leks !slat ,
     .     , fes1_rv , fes2_rv , fh_rv , frc )
!!   --- Mix input and output densities ---
c         allocate(qm1_rv(nvl))
c         allocate(qm2_rv(nvl))
         call mixrho ( ssite , sspec ,  nsp , iter , sstrnmix 
     .     , qval - qbg , elind , sv_p_orhat1 , sv_a_oorhat , k1 , k2
     .     , k3 , dmxp , srout_zv , zv_a_osmrho )
ccccccccccccccccccccccccc
ctakao test mkekin repeat
cc          print *,'vvv takao mkekin copy test xxx2---'
cc          call dcopy(k1*k2*k3*2,smrhobk,1,w(osmrho),1)
cc          deallocate(smrhobk)
c          print *,'vvv takao mkekin repeat 2 ---'
c          call dcopy(k1*k2*k3*4,w(osmpotbk),1,w(osmpot),1)
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sev,sumtvxxx)
c          stop 'xxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccc
      else
         eks = 0d0
      endif
      ham_ehf=ehar
      ham_ehk=eks
 999  continue
c      if (allocated(qm2_rv)) deallocate(qm2_rv)
c      if (allocated(qm1_rv)) deallocate(qm1_rv)
c      if (allocated(ifbls_iv)) deallocate(ifbls_iv)
      if (allocated(hbyl_rv)) deallocate(hbyl_rv)
      if (allocated(qbyl_rv)) deallocate(qbyl_rv)
      call m_mkpot_deallocate()
      deallocate(evlall)
 9995 continue
      call tcx('bndfp')
      return
!! end of bndfp
      end subroutine bndfp
      end module m_bndfp


      
#if MPIK
      subroutine xmpbnd(kpproc,ndham,nkp,nsp,eb)
C- Collect eb from various processors (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   kpproc
Ci   ndham :leading dimension of eb
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   eb    :energy bands; alias eband
Co Outputs
C ----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer kpproc(0:*),ndham,nkp,nsp
      double precision eb(ndham,*)
      integer i,ista,iend
      integer procid,master
      integer numprocs, ierr
      integer, dimension(:),allocatable :: offset,length
      real(8) ,allocatable :: buf_rv(:)
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0
      do  i = 0, numprocs-1
        ista = kpproc(i)
        iend = kpproc(i+1)-1
        length(i) = (iend - ista + 1)*nsp*ndham
        offset(i+1) = offset(i) + length(i)
      enddo
      ista = kpproc(procid)
      allocate(buf_rv(ndham*nkp*nsp))
#ifdef KINODEBUG
      write(*,*)'kino allocate bug_rv',ndham*nkp*nsp
#endif
      call mpi_allgatherv ( eb ( 1 , 1 + nsp * ( ista - 1 ) ) , length 
     .( procid ) , mpi_double_precision , buf_rv , length , offset 
     ., mpi_double_precision , mpi_comm_world , ierr )
      call dcopy ( ndham * nsp * nkp , buf_rv , 1 , eb , 1 )
      if (allocated(buf_rv)) deallocate(buf_rv)
      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)
      end subroutine xmpbnd
#endif

