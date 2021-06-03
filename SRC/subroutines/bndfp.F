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
!!
      !! emin=dosw(1) emax=dosw(2) dos range
      module m_bndfp
      real(8),protected:: ham_ehf,ham_ehk,bz_dosw(2)
      real(8),protected:: eferm !,bz_def ! revised when called.
      logical,private:: binit=.true.,bbinit=.true.
      
      contains
      subroutine m_bndfp_ef_set(bz_ef00)
      real(8):: bz_ef00
      eferm  = bz_ef00
      binit=.false.
      end
      
      subroutine bndfp( dmxp,iter,maxit, frc, llmfgw) !vorb
      use m_supot,only: zv_a_osmrho,sv_a_oorhat
      
      use m_supot,only: ngabc=>lat_nabc
      use m_suham,only: iv_a_oindxo, nchan=>pot_nlma, nvl=>pot_nlml, ndham=>ham_ndham,ham_ldham
     & ,ndhamx=>ham_ndhamx,nspx=>ham_nspx 
c      use m_hamindex, only: ngrp,symops !,norbmto,ibastab,ltab,ktab,offl, symops_af
      use m_lmfinit, only: n0,nab,nppn,ncutovl,lso,ndos=>bz_ndos,bz_ef0,bz_w,fsmom=>bz_fsmom,
     &     bz_dosmax,lmet=>bz_lmet,bz_fsmommethod,bz_n,
     &     ctrl_nspec,ctrl_pfloat,ldos=>ctrl_ldos,ctrl_zbak,lfrce=>ctrl_lfrce,
     &     pwmode=>ham_pwmode,lrsig=>ham_lsig,epsovl=>ham_oveps,ham_eterms,ham_elind,ham_rsrnge,
     &     ham_scaledsigma,ham_pmin,ham_pmax, !
     &     alat=>lat_alat,stdo,stdl,procid,master, !bz_doswin, 
     &     nkaph,nlmax,nl,nbas,nsp, ham_frzwf, bz_dosmax,
     &     sstrnmix,
     &     lcplxp ,  lekkl,lmaxu,nlibu,lldau,ssite=>v_ssite,sspec=>v_sspec, lpztail,plbnd,leks,lrout
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
      use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
      use m_seneinterp,only: Senebroadcast,Rdsigm2,ndimsig,nspsigm,Sigmainit,Getsenex,sene,Dsene
      use m_seneinterp,only: ham_ndhrs,ham_nqsig,senex
      use m_subzi, only: m_subzi_init, lwtkb,nevmx,numq,lswtk,rv_a_owtkb,rv_a_oswtk !,ldos
      use m_rsibl_ev, only : rsibl_ev ! to plot wavefunction in the fat band mode
      use m_struc_def,only : s_rv1
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq
      use m_lmfgw,only: m_lmfgw_init,sv_p_osigx,sv_p_otaux,sv_p_oppix,spotx,jobgw
      use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate,
     o    zv_a_osmpot, qmom_rv, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o     , vval_rv , fes1_rv
      use m_clsmode,only: m_clsmode_init,m_clsmode_set1,m_clsmode_finalize
      use m_writeham,only: m_writeham_init, m_writeham_write
      use m_qplist,only: m_qplist_init, qplist,nkp,xdatt,labeli,labele,dqsyml,etolc,etolv,
     &     nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml
      use m_shortn3,only: shortn3, shortn3_initialize
      use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_set
      use m_procar,only: m_procar_init,dwgtall,nchanp,m_procar_closeprocar,m_procar_writepdos
      use m_bandcal,only: m_bandcal_init,m_bandcal_2nd,m_bandcal_clean,m_bandcal_allreduce,
     o     evlall,srout_zv,sumqv,sumev,sv_p_oqkkl,sv_p_oeqkkl, frcband, ndimhx_,nev_,dmatu,nevls
      use m_bzintegration,only: m_bzintegration_init
      use m_ldau,only: vorb
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
Co  dmatu :density matrix for LDA+U (changed upon output)
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
C     l         :otherwise 1

      
!!      nspx: number of independent spin channels
!!      nspc is now avoided (memo:nspc=2 for lso==1, nspc=1 for lso/=1 See m_lmfinit)
!!      ndhamx: maximum size of hamiltonian  

C     l   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
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
      integer obuf,oqptr
      integer igetss,lmxh,nlmh,nelt(3),kmax !,nlma
      integer, dimension(:),allocatable :: kpproc
      double precision sttime,entime
c      integer:: pidorigin, ndimhx_iq, nev_iq,idx
c      integer:: work_mpi(10)    ! use only 1 and 2 now
      integer nbx,iter,maxit
      integer:: k1,k2,k3,lpnu=1
      real(8):: dmxp(33),frc(3,1:nbas*numq)
      logical lgors,ltet,cmdopt,cmdopt0,cmdopt2!,lwndow
      character strn*120,strn2*120
      integer:: ipr , iprint , isp , jsp , ispc , isum , nglob 
     .  , i , iq , ismidb  ! ldos  !lmet !, nspc  , nsp
     .  , lrep ,  nev  
     .  , n1 , n2 , n3 !lwtkb !,  lswtk=0 !, nl lwtkb ,
     .  , lgunit , igets , mpsord ,  isw  !stdo , stdl 
     .  , ifi , fxst ,   i1 , i2 !, !plbnd lfrzw ,
     .  , nfilem , iobzwt ,   nll, lnoxc !, lpdiag lwsig  , 
     .  , iv(10) , parg , nk1 , nk2 , nk3 , lshft(3) 
     .  ,  mpipid , ldham(16) , ldim !, ndimh , ndimhx !nfbn(2) , 
      integer:: os , osoptc
      real(8) ,allocatable :: fh_rv(:)
      real(8) ,allocatable :: fes2_rv(:)
      real(8) ,allocatable :: hbyl_rv(:,:,:)
      real(8) ,allocatable :: qbyl_rv(:,:,:)
      complex(8) ,allocatable :: t_zv(:,:)
      integer ,allocatable :: ips_iv(:)
      type(s_rv1) :: sv_p_orhat1(3,nbas)
      double precision sigp(10)
      logical:: readeferm
      double precision esmearm,sev,efermr,
     .  sev00,sev1,sumtv,alfa,qbg,ehar,eks,
     .  dosw(2),dgets !,sumev(2,3),sumqv(3,2)
     .  ,qp(3),shftqp(3),eterms(20),elind,dum,evtop,
     .  ecbot,rsrnge,xv(20),dosrng,dval !,eomin ,vrmt(nbas)
!!    PW basis
c      integer,allocatable:: igv2x(:,:)
      integer ,allocatable :: kv_iv(:)
c      integer pwmode
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
      integer:: nblk,nprow,npcol,nmx
!!     for LDA+U
      real(8) ,allocatable :: ausu_rv(:)
c      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      integer:: i_spackv, i_spacks
      integer ::iwdummy
      logical :: debug=.false.  !debug is only for check write. However, it can cause error
!for 'make check',
!because fp/test/test.fp assumes 'something not inserted in console output'.

      data T /.true./ F /.false./
c      integer,allocatable:: nevls(:,:)
ctakao
      integer :: ifiogw,ifiese,og,ib,ig,iorb,is,nsgrp,uspec_size,nspec
      logical :: nexist,sigmamode,ok !, newsigmasw
      character spid*8
      integer :: ib1,ib2,ii,ik1,ik2,ik3,j1,j2,ikx
      integer:: nhmax , nknknk , jx , nsp_dummy , ndimh_dummy
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),ovlmtoi(:,:),ovliovl(:,:) 
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
      real(8):: qpx(3)          !qpo(3),

      integer:: iqini,iqend
      complex(8),allocatable:: smrhobk(:)
      real(8),allocatable::evl(:,:) ! evlall(:,:,:),
      integer:: idummy
      real(8):: epsnevec
      logical:: onlyonce=.true.
      logical:: isanrg, l_dummy_isanrg,oncewrite !, savez,getz,gotosecondpass,initiqisp
      integer:: getef, i_dummy_getef,ifig!, ndimsig
c      integer:: n_listwf
      integer, allocatable:: listwf(:)
      integer:: ix99 !,nevec unused.
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
c      integer,allocatable:: ndimhx_(:),nev_(:)
      integer:: npri,ibas,ifiqb,ifiwv,lmxa, ifimag
      real(8),parameter::    NULLR =-99999

      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw
      real(8)::ppin(3)

      logical:: fsmode,allband !for --fermisurface for xcrysden.
      integer:: ififm
      real(8):: qshort(3)
! add --fermisurface #1 #2 #3
! only implemented for lmf-MPIK now.
      integer:: ini,iend,nbxx,ival

c mar2014
      integer:: lb1,lb2,onespdummy
      integer:: iprocar,iprocar1,iprocar2,lrsigx
      integer:: l,m,ilm,nlma,iband,hreal,mode,nttabs
      character*256:: sss,s222
      character*1000::ccc

      integer:: ispsigm !nov2014
      logical:: goto99!jan2015
      logical:: fullmesh,PROCARon

      integer:: ntete,init,ient,nevmin,ifip,ifisyml
      real(8):: qx(3)
      real(8),allocatable:: wkd(:), pdosall(:,:,:)
      character*3::charnum3
      character*100::filenm(2)
      real(8):: emin

      character(512):: aaachar
      logical:: massline , writeham=.false.
      logical,save:: siginit=.true. 
      character(10):: i2char
      integer::ipos,ndos_,j
      real(8)::dosmax,posp(3)
      integer:: jspend,jspini,jspp
      integer:: iqread,iqindex,ibread(100)=0,ibindex,nibread, job,ist,ip,ni,ix
   
      integer::npairmx,ikpd,ifqplist
      integer,allocatable:: nqwgt(:,:,:)
      real(8):: tpioa,pi=4d0*datan(1d0)
      complex(8):: aaaa,img=(0d0,1d0)

      real(8):: pnu(n0,2),pnz(n0,2),pmean,symxx(3,3),qvec(3),dee,dosi(2)
      integer:: i_copy_size, ndummy1,ndummy2,ndummy3,ncount,nbandx,ifii
      logical:: phispinsym,ljobgw
      integer:: nr,lmxl,nlml,nelt1,nelt2,nso,ifih,ifile_handle,ispinit,ndimhxx,nnn
      character*8:: charext
      logical:: lfill=.false.,ixx,do2010init(2)
      logical:: spinoffdiag
      real(8):: kkk ,ef0xx,efermxxx      !,eferm !dummy
      real(8),allocatable:: dos_rv(:,:)
!! ----------------------------
      call tcn ('bndfp')
      write(stdo,"(a)") 'bndfp:start'
      if(master_mpi) then ! check Write atomic positions
         write(stdo,"(a)")'   Atomic Positions from rst/ctrl. '
         write(stdo,"(a)")'   Note: rst and ctrl may contain different atomic positions when DYN mode'
         do ib=1,nbas
            is  = ssite(ib)%spec   !species index
            spid= ssite(is)%clabel !atom label
            write(stdo,"('ib clabel pos =',i5,a6,3(x,f13.5))") ib, trim(spid),ssite(ib)%pos
         enddo
      endif
      debug    = cmdopt0('--debugbndfp')
      tdos     = cmdopt0('--tdos')
      fsmode   = cmdopt0('--fermisurface')! Fermi surfece for xcrysden in http://www.xcrysden.org/doc/XSF.html#2l.16
      ipr  = iprint() ! for procid/=master, we set iprint=0 at lmv7.F
      idummy = unlink('RSEQ_ERROR') ! the error flag file caused by rseq. Remove RSEQ_ERROR here because
                                    ! this error is fatal only when it is converged. takaoNov28_2010
      allocate(evl(ndhamx,nspx))
      ldham = ham_ldham
      ldim  = ldham(1)
      if (master_mpi) then
        aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
        write(stdo,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
      endif
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
      if(binit) then ! eferm is from ctrl file when m_bndfp_ef_set in not called (binit=T)
         eferm =bz_ef0
         binit=.false.
      endif
      ltet = ntet>0
      qbg = ctrl_zbak(1) !homogenious background charge
      alfa = 0
!! Set spin-symmetrized pnu. aug2019. See also in pnunew and locpot
      if(cmdopt0('--phispinsym')) call phispinsym_set(ssite) !pnu,pz spin symmetrized      
!! Write out Hamiltonian HamiltonianPMT.* and HamiltonianPMTinfo
      writeham= cmdopt0('--writeham')
      if(writeham) then
         call m_writeham_init()
         open(newunit=ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
      endif      
!! Make one-particle potential without XC part for GWdriver: lmfgw mode 
      if(llmfgw) call m_lmfgw_init() ! Initializer of m_lmfgw. Get sig,tau,ppi,spot without XC(LDA) via mkpot.
!! Generate one-body potential and energy-related quantities. See use m_mkpot:
      call m_mkpot_init(llmfgw) !,vorb) 
!! elind for accelarating convergence
      elind = ham_elind
      if(ham_elind<0d0) elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*ham_elind
      dmxp(33) = elind
!! Setup for BZ integration wtkb
!! NOTE: if (we have wkp.* file).and.lmet==2, wkp file is used for wkpb.
      call m_subzi_init(lrout>0) !,def,bz_n,bz_w,eferm) !get nevmx, and so on
      
!   TK think sspec%ngcut of local orbital is changing when PZ=+10 mode (hankel tail mode)
!   because Hankel's e may be changing. ! lpztail: local orbital of 2nd type(hankel tail) +10 exist or not
      if(lpztail) call sugcut(2, sspec) !sspec%ngcut it was suham2(sspec)
      if(cmdopt0('--cls')) then !clsmode
          call rxx(lso==1,'CLS not implemented in noncoll case') 
         if (lrout .eq. 0) call rx('bndfp: need output density for cls')
         call m_clsmode_init()
      endif
      if(onlyonce) then
         call gen_hamindex(.true., jobgw==0) !index for hamiltonian
         if(writeham) call m_writeham_write() 
         onlyonce=.false.
      endif
!!      
      if(plbnd/=0) then
! Get eferm from efermi.lmf for plbnd mode. Note that efermi.lmf can be different when we use
! different number of k points. For example, you may change NKABC, and run job_pdos.
         open(newunit=ifi,file='efermi.lmf')
         read(ifi,*) eferm
         close(ifi)             
         call mpibc1(eferm,1,4,.false.,'bndfp','eferm')
      endif
      
!! === Get nkp and qplist. fatband mode ---> --mkprocar. plbnd/=0 ---> band plot mode.
!! we can move this outside of iteration cycle.
      if(debug) write(stdo,*) 'bndfp: goto m_qplist_init'
      if(.not.llmfgw) call m_qplist_init(plbnd)
      call mpi_barrier(MPI_COMM_WORLD,ierr) 
!! === \Sigma-Vxc to ham_rv_a_ohrs. Neighbor table ham_iv_a_oiaxs ====
!! ndimsig is the dimension of the self-energy. 
!! ndimsig = ldim means we use only projection onto MTO spaces even when PMT. !aug2012takao
!! --- New sigma inerpolation mode. Simplified.
!! OUTPUT in this sigmainit is used in
!!   call bloch 
!!  where we use ham_iv_a_oiaxs , ham_rv_a_ohrs, ham_ndhrs and so on.
!!  In future, we should make this block to a subroutine which shows clear i/o.
      sigmamode = mod(lrsig,10) .ne. 0
      if( llmfgw.and.jobgw==0 ) then
         continue
      elseif( sigmamode.and.siginit ) then
         call sigmainit()
         siginit=.false. !this is equivalent to broad casting.
      endif              
      if(sigmamode.and.master_mpi)write(stdo,*)' ScaledSigma=',ham_scaledsigma
!! == GW driver and finish ==
      if(llmfgw) then
         call sugw(jobgw) ! NOTE: Sigma-Vxc in real space is give by rv_a_ohrs.
         call tcx('bndfp')
         call rx0('sugw mode')  !exit program here normally.
      endif

!! MPIK k point divider. From iqini to iqend for each processor.    
#if MPIK
      sttime = MPI_WTIME()
      if(master_mpi) write(stdo,"(a)") ' ... Start MPI k-loop ...'
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(nkp,numprocs,1,kpproc(0)) !mpik distribution
      iqini = kpproc(procid)
      iqend = kpproc(procid+1)-1
#else
      iqini=1
      iqend=nkp
#endif
!! Set up Hamiltonian and diagonalization in m_bandcal_init
!! To know outputs, see 'use m_bandcal,only:'. The outputs are evlall, and so on.
!     ! ndimhx_ is protected, but allows writing to other nodes because MPI_BCAST is f77, assumed-size array
      if(.not.(lwtkb==-1.or.lwtkb==0.or.lwtkb==1)) call rx('bndfp: something wrong lwtkb')
! lwtkb=1 accumulate weight
! lwtkb=0 no weight
! lwtkb=-1      
      call m_bandcal_init(iqini,iqend,lrout, eferm,ifih,lwtkb) 
      entime = MPI_WTIME()
      if(master_mpi) write(stdo,"(a,f9.4)") ' ... Done MPI k-loop: elapsed time=',entime-sttime
      sttime = MPI_WTIME()
      if(writeham) close(ifih)
      if(writeham) call rx0('Done --writeham: --fullmesh may be needed. HamiltonianMTO* genereted')
!! Broadcast. Use allreduce to avoid knowing which node to which node.
!! Because mpi is f77, ndimhx_ and so on are not protected.
      nnn = size(ndimhx_)
      call mpibc2_int(ndimhx_,nnn,'bndfp_ndimhx_') !all reduce 
      nnn = size(nev_)
      call mpibc2_int(nev_,nnn,'bndfp_nev_') !all reduce (to avoid which node to which node).
      nnn = size(nevls)
      call mpibc2_int(nevls,nnn,'bndfp_nevls') !all reduce (to avoid which node to which node).
      call xmpbnd(kpproc,ndhamx,nkp,nspx,evlall) !all eigenvalues broadcasted
      evtop=-9999
      ecbot=9999
      do iq=1,nkp
      do jsp=1,nspx          !nspx=1 for SO=1
         if(lso==1) i = max(1,nint(qval-qbg))
         if(lso/=1) i = max(1,nint(qval-qbg)/2)
         evtop = max(evtop,evlall(i,jsp,iq))
         ecbot = min(ecbot,evlall(i+1,jsp,iq))
         if(master_mpi.and.iq==1)write(stdl,"('fp evl',8f8.4)")(evlall(i,jsp,iq),i=1,nev_(iq))
      enddo
      enddo
      
!! pdos mode (--mkprocar and --fullmesh)
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface')
      PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
      if(fullmesh.and.procaron) call m_procar_writepdos(evlall,nev_,eferm,kpproc)
      if(fullmesh.and.procaron) call rx0('Done pdos: --mkprocar & --fullmesh. Check by "grep k-point PROCAR.*.*"')
!! boltztrap data 
      if( cmdopt0('--boltztrap')) call writeboltztrap()
      if( cmdopt0('--boltztrap')) call rx0('Done boltztrap: boltztrap.* are generated')
!! Write bands in bands-plotting case: loop over qp getting evals from array ===
      if(plbnd/=0.and.master_mpi.and.fsmode) call writefs(ndhamx,nspx,nkp,evlall,qplist,eferm)
      if(plbnd/=0.and.master_mpi) then
            write(stdo,*)' Writing bands to bands file ...'
            if(nsyml/=0) call writeband(eferm,evtop,ecbot)
c              call writeband(stdo,nkp,nsyml,nspx,nsp,ndhamx, !,fsmom/=NULLR,vnow
c     i        xdatt,nqp_syml,nqp2n_syml,nevls,evlall,qplist,labeli,labele,alat,eferm,evtop,ecbot,
c     i        sname,nqps_syml,nqpe_syml,nsyml,dqsyml,
c     i        etolv,etolc)
      endif
      if(plbnd/=0.and.fsmom/=NULLR) then
         write(stdo,"(a)")'Caution: fsmom(fixed moment on). In sc cycle, we use additional bias mag. field  '
         write(stdo,"(a)")'Caution: Mag.field is written in MagField. But it is not used for --band mode!'
      endif
      if(plbnd/=0.and.fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
      if(plbnd/=0)            call rx0('plot band mode done')
!! end of band plot mode
      
!! New eferm determined from evlall
      call m_bzintegration_init(evlall,  eferm,sev,sumqv,vnow) !Get new eferm, sev, wtkb...
      if(lmet==0) eferm = (evtop+ecbot)/2d0 !for metal
      if(lmet==0.and. master_mpi)
     &     write(stdo,"('Highest occ., Ef, Lowest unocc.=',3f10.5)")evtop,eferm,ecbot,eferm
      if(master_mpi) then
         open(newunit=ifi,file='efermi.lmf')
         write(ifi,"(d24.16, ' # (Ry) Fermi energy given by lmf')") eferm
         write(ifi,"(d24.16, ' # (Ry) Top of Valence')") evtop
         write(ifi,"(d24.16, ' # (Ry) Bottom of conduction')") ecbot
         write(ifi,"(d24.16,d24.16 ' # number of electrons qval-qbg, qbg=')") qval-qbg,qbg
         write(ifi,"(d24.16, ' # mag. mom (nup-ndown)=')") sumqv(2,1) 
         write(ifi,"(3i10,'    # Used k point to determine Ef')") nkabc
         write(ifi,
     &     "('# NOTE: These valuse are overwritten by lmf(-MPIK) in SC loop')")
         write(ifi,"('# NOTE: Not by band plot mode')")
         close(ifi)             
      endif
!! Generate total DOS 
!!    dos window
      dosw(1)= emin-0.5d0       ! lowest energy limit to plot dos
      dosw(2)= eferm+bz_dosmax  ! lowest energy limit to plot dos
      if(master_mpi.and.(tdos.or.ldos/=0)) then
         write(stdo,"(a,3f9.4)") 'Generating TDOS: efermi, and dos window= ',eferm,dosw
         allocate( dos_rv(ndos,2*nspx))
         if(cmdopt('--tdostetf',10,0,strn)) ltet=.false. !oct2015
         if(ltet) then 
            call bzints(nkabc(1),nkabc(2),nkabc(3), evlall 
     &           , dum , nkp , ndhamx , ndhamx , nspx , dosw(1), dosw(2), dos_rv , ndos , efermxxx ,
     &           1, ntet , iv_a_oidtet , dum , dum ) !job=1 integrated dos dos_rv
            call xxxdif(dosw(1),dosw(2),ndos,nspx,0,dos_rv ) !derivative 'integrated dos -> DOS'
         else
            call makdos (nkp, ndhamx, ndhamx, nspx, rv_a_owtkp, evlall,bz_n, bz_w
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
         endif
         if( lso==1 ) call dscal ( ndos , .5d0 , dos_rv , 1 )
         open(newunit=ifi,file='dos.tot.'//trim(sname))
         open(newunit=ifii,file='dosi.tot.'//trim(sname))
         dee=(dosw(2)-dosw(1))/(ndos-1d0)
         do ipts=1,ndos
            eee=dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-efermxxx !I think efermxxx=eferm
            dosi(1:nsp)= dosi(1:nsp) + dos_rv(ipts,isp)*dee
            write(ifi,"(255(f13.5,x))") eee,  (dos_rv(ipts,isp),isp=1,nsp)
            write(ifii,"(255(f13.5,x))") eee, (dosi(isp),isp=1,nsp)
         enddo
         close(ifi)
         close(ifii)
      endif
      if(tdos) call rx0('Done tdos mode:')
      
!! AllReduce band quantities.  eterms should be cleaned up in future.   
      call m_bandcal_allreduce(lrout,lwtkb)
      eterms=ham_eterms 
      call mpibc2_real(eterms(19),1,'bndfp_rhosig')
      ham_eterms=eterms
      emin=1d9
      do iq = 1, nkp
         qp=qplist(:,iq)
         do isp = 1, nspx 
          jsp = isp
          if(iprint()>20) write(*,'(" bndfp: kpt",i5," of",i5, " k isp=",3f8.4,i2," ndimh nev=",2i5)')
     &      iq,nkp,qp,jsp,ndimhx_(iq),nevls(iq,jsp)
          if(iprint()>=35) write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
          emin= min(minval( evlall(1:nevls(iq,jsp),jsp,iq)),emin)
         enddo
      enddo
!! Accumurate quantities for determined wtkb. Originally called as second pass.
      if(lmet>0 .and. lwtkb==-1 .and. lrout>0) then 
         lwtkb = 1
         call mpi_barrier(MPI_comm_world,ierr)
         call m_bandcal_2nd(iqini,iqend,lrout, eferm) !accumulate srout_zv and so on.
         call m_bandcal_allreduce(lrout,lwtkb)
         eterms=ham_eterms 
         call mpibc2_real(eterms(19),1,'bndfp_rhosig')
         ham_eterms=eterms
      endif
      
!! Core-level spectroscopy --- m_clsmode_set1 is called in m_bandcal
!!    dos window
      dosw(1)= emin-0.5d0       ! lowest energy limit to plot dos
      dosw(2)= eferm+bz_dosmax  ! lowest energy limit to plot dos
      if (cmdopt0('--cls')) then
         if (master_mpi) call m_clsmode_finalize(eferm,ndimh,ndhamx,nspx,nkp,dosw,evlall)
        call rx0('Done cls mode:')
      endif
      
!! contribution to force via band 
      if(lfrce>0) frc=frcband   !part of atomic force from band part
!! write orbital moment      
      if(lwtkb==1 .and. lso /=0) call iorbtm() !write out orbital moment

      
!! Assemble output density, energies and forces ===============================
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
      call mkrout ( ssite , sspec ,  nbas , nsp , ldim  !sham ,slat , 
     .  , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_orhat1 , hab_rv , sab_rv 
     .  , qbyl_rv , hbyl_rv , lrout )
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
         call symrho ( ssite, sspec,  lfrce, srout_zv, sv_p_orhat1, qbyl_rv, hbyl_rv, frc )
C   ... New boundary conditions pnu for phi and phidot
         if (lpnu .gt. 0) then
            call dcopy(size(ham_pmin),ham_pmin,1,xv,1)
            call dcopy(size(ham_pmax),ham_pmax,1,xv(11),1)
            call pnunew ( nbas , nsp , ssite , sspec , xv , xv ( 11 ) , ham_frzwf 
     .        , hab_rv , sab_rv , qbyl_rv , hbyl_rv, eferm )
         endif
      endif
!! --- print out qbyl ---
      if(master_mpi) then
         open(newunit=ifqbyl,file='qbyl.'//trim(sname))
         write(ifqbyl,"(i5)")nbas
         do  ib = 1, nbas
           ispec=ssite(ib)%spec
           write(ifqbyl,"(i4, i3,i2,10f12.6)")ib, ispec, sspec(ispec)%lmxa, (sum(qbyl_rv(il,1:nsp,ib)),il=1,sspec(ispec)%lmxa+1)
         enddo  
         close(ifqbyl)
      endif
      
!! --- Evaluate Harris energy ---
      call mkehkf(1,sev,sumqv(2,1),sumtv,ehar)
!! --- Evaluate KS total energy, correction to Harris force ---
      if(lrout/=0) then
!!   ... Correction to harris force
        if(lfrce>0 ) then
          allocate( fh_rv(3*nbas)  )
          fh_rv= 0d0  
          call dfrce (lfrce, ssite , sspec ,  k1 , k2 , k3 , nvl !sctrl , slat ,
     .     , sv_a_oorhat , sv_p_orhat1 , elind , qmom_rv , zv_a_osmrho , 
     .     srout_zv , fh_rv )
        endif
!!   ... Evaluate KS total energy and output magnetic moment
         eks = 0d0
         if(leks .ge. 1) then
            call togpr()
            call mkekin ( nbas , ldim , ssite , sspec ,  sv_p_osig  !sham ,slat , 
     .        , sv_p_otau , sv_p_oppi , lcplxp , sv_p_oqkkl , k1 , k2 , k3
     .        , vconst , zv_a_osmpot , srout_zv , sev , sumtv )
!! mkpot calls locpot. and locpot calls augmat.
!! augmat calculates sig,tau,ppi.
c          print *,'goto mkpot procid=',procid
            allocate( fes2_rv(3*nbas))
            call mkpot (lfrce,
     &          srout_zv , sv_p_orhat1 , qbg , zv_a_osmpot 
     .        , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .        , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .        , vval_rv , fes2_rv , 0 , vorb ) !job=0 means energy term only

           if(cmdopt('--density',9,0,strn)) then
             call mpi_barrier(MPI_comm_world,ierr)  
             call rx0('end of --density mode')
           endif  
            call mkehkf(2,sev,sumqv(2,1),sumtv,eks) !sham,
            call togpr()
         endif
!! --- Add together force terms ---
!! fes1_rv: contribution to HF forces from estat + xc potential   This is for input  density  !=3rd term in (B.5) in JPSJ.84.034705
!! fes2_rv: contribution to KS forces from estat + xc potential   This is for output density
!! fh_rv      : 1st term  in (B.5)
!! input frc  : 2nd term in (B.5). Puley term
!! output frc  : total
         if(lfrce> 0) then
            call totfrc( nbas , ssite ,  leks, fes1_rv , fes2_rv , fh_rv , frc )
            deallocate(fh_rv)
         endif
!!   --- Mix input and output densities ---
         call mixrho ( ssite , sspec ,  nsp , iter , sstrnmix 
     .     , qval - qbg , elind , sv_p_orhat1 , sv_a_oorhat , k1 , k2
     .     , k3 , dmxp , srout_zv , zv_a_osmrho )
      else
         eks = 0d0
      endif
      
      ham_ehf=ehar
      ham_ehk=eks
      if (allocated(hbyl_rv)) deallocate(hbyl_rv)
      if (allocated(qbyl_rv)) deallocate(qbyl_rv)
      call m_mkpot_deallocate()
      call m_bandcal_clean()
      call tcx('bndfp')
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

      subroutine readindensitymodesetting()
!! If necessary, need to recover this but it should be much simpler.
!! By chaningt wtkb, this will give density for given iq ib.
!!  We have electron density plot by '--density' mode
c$$$!! Density mode. wtkb is the weight for iq,ib. With this setting, 
c$$$!! we can plot |\psi_{iq,ib}(\bfr)|**2 in the manner of electron density
c$$$!     ! lmf --density,iq=12,ib=3,4,5   (here ib=3,4,5 are superposed).
c$$$      if(cmdopt('--density',9,0,strn)) then
c$$$         iqread=0
c$$$         iqindex=index(strn(10:),'iq=')+2
c$$$         if(iqindex==2) then 
c$$$            iqread=-999
c$$$         else
c$$$            read(strn(10+iqindex:),*) iqread
c$$$         endif
c$$$         ibindex=index(strn(10:),'ib=')+2
c$$$         if(ibindex==2) then
c$$$            ibread=-999
c$$$         else
c$$$            do nibread=1,100
c$$$               read(strn(10+ibindex:),*,err=2019,end=2019) ibread(1:nibread)
c$$$c     print *,'xxx ibread=',ibread(1:nibread)
c$$$            enddo
c$$$ 2019       continue
c$$$            nibread=nibread-1
c$$$         endif  
c$$$         if(iqread>0) then
c$$$            if(maxval(ibread(1:nibread))<=0) call rx('--density mode: wrong ib=foobar. Try,e.g. --density,iq=12,ib=5')
c$$$            write(stdo,"('--density bandmode: psi**2 for iq=',i5,' ib=',255i5)") iqread,ibread(1:nibread)
c$$$            rv_a_owtkb(:,:,:)=0d0
c$$$            do ib=1,nibread
c$$$               rv_a_owtkb(ibread(ib),:,iqread)=1d0
c$$$            enddo
c$$$            do ib=1,ndimh
c$$$               if(abs(rv_a_owtkb(ib,jsp,iq))>1d-3) write(stdo,"('ib wtkb=',i5,2f13.6)") ib,rv_a_owtkb(ib,jsp,iq)
c$$$            enddo
c$$$         endif
c$$$      endif   
      end
