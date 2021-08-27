!     ! = Main part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration =
!     !
!     ! From density osmrhod and osmrho, m_mkpot_init gives  osmpot, and sv_p_osig, sv_p_otau, sv_p_oppi (potential)
!     ! Then we do band calcultion m_bandcal_init. It gives osmrho_out and orhoat_out.
!     ! Finally, osmrho and orhoat are returned after mixing procedure (mixrho).
!----------------------------      
!io  smrho :smooth density
!io        :Symmetrized on output
!io  orhoat:vector of offsets containing site density
!io        :Symmetrized on output
!  qbyl  :site- and l-decomposed charges
!  hbyl  :site- and l-decomposed one-electron energies
!o   dmatu:   density matrix in m_bandcal
!o   evalall: eigenvalue in m_bandcal
!o   force: 
!     i  vorb --> use m_ldau in m_mkpot
!!--- key quntities --------------     
!  evlall  : eigenvalues
!  eferm   : Fermi energy
!     eksham   : DFT(Kohn-Sham)  total energy
!     eharris  : Harris foulkner total energy      
      module m_bndfp
      use m_density,only: orhoat,osmrho     !input/output unprotected
      real(8),protected:: ham_ehf, ham_ehk  !output
      real(8),protected:: eferm, qdiff !output
      real(8),allocatable:: force(:,:)
!! other ouput are in modules m_mkpot, m_bandcal       
      logical,private:: binit=.true.,initd=.true.
      contains
      
      subroutine m_bndfp_ef_set(bz_ef00) !called from iors.F. grep 'use m_bndfp'
      real(8):: bz_ef00
      eferm  = bz_ef00
      binit=.false.
      end
      
      subroutine bndfp(iter, llmfgw, plbnd)
! llmfgw=T is for generating eigenfunctions for GW calculations, no iteration.
! plbnd/=0 means band plot mode. no iteration.      
!     ! All read only in bndfp. Data are stored in modules such as m_bandcal, m_mkpot
!     ! For example,, rightafter call m_bandcal_init, we can get evalall, which is used in other modules.
      use m_supot,only: ngabc=>lat_nabc,k1,k2,k3 !for charge mesh
      use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx,nspx=>ham_nspx 
      use m_lmfinit, only: n0,nab,nppn,ncutovl,lso,ndos=>bz_ndos,bz_ef0,bz_w,fsmom=>bz_fsmom,
     &     bz_dosmax,lmet=>bz_lmet,bz_fsmommethod,bz_n,
     &     ctrl_nspec,ctrl_pfloat,ldos=>ctrl_ldos,qbg=>zbak,lfrce=>ctrl_lfrce,
     &     pwmode=>ham_pwmode,lrsig=>ham_lsig,epsovl=>ham_oveps,ham_elind,ham_rsrnge,
     &     ham_scaledsigma,ham_pmin,ham_pmax, !
     &     alat=>lat_alat,stdo,stdl,procid,master, !bz_doswin, 
     &     nkaph,nlmax,nl,nbas,nsp, ham_frzwf, bz_dosmax,
     &     lekkl,lmaxu,nlibu,lldau,sspec=>v_sspec, lpztail,leks,lrout
     &     ,  nchan=>pot_nlma, nvl=>pot_nlml,nspc
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
      use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
      use m_rdsigm2,only: m_rdsigm2_init
      use m_subzi, only: m_subzi_init, lwtkb,rv_a_owtkb,m_subzi_setlwtkb,m_subzi_bzintegration
      use m_rsibl_ev, only : rsibl_ev ! to plot wavefunction in the fat band mode
      use m_struc_def,only : s_rv1
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq
      use m_lmfgw,only: m_lmfgw_init !,jobgw !,sv_p_osigx,sv_p_otaux,sv_p_oppix,spotx
      use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate, m_mkpot_energyterms,m_mkpot_novxc,
     o    osmpot, qmom, vconst, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , sab_rv , qval , qsc , fes1_rv , fes2_rv
      use m_clsmode,only: m_clsmode_init,m_clsmode_set1,m_clsmode_finalize
      use m_writeham,only: m_writeham_init, m_writeham_write
      use m_qplist,only:  qplist,nkp,xdatt,labeli,labele,dqsyml,etolc,etolv,
     &     nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml,
     &  iqini,iqend,ispini,ispend,kpproc    ! MPIK divider. iqini:iqend are node-dependent
      use m_shortn3,only: shortn3, shortn3_initialize
      use m_igv2x,only: napw,ndimh,ndimhx,igv2x
      use m_procar,only: m_procar_init,dwgtall,nchanp,m_procar_closeprocar,m_procar_writepdos
      use m_bandcal,only: m_bandcal_init,m_bandcal_2nd,m_bandcal_clean,m_bandcal_allreduce,
     o     evlall,smrho_out,sv_p_oqkkl,sv_p_oeqkkl, ndimhx_,nev_,nevls,m_bandcal_symsmrho
      use m_mkrout,only: m_mkrout_init,orhoat_out,frcbandsym,hbyl_rv,qbyl_rv
      use m_addrbl,only: m_addrbl_allocate_swtk,swtk
      use m_sugw,only: m_sugw_init
      use m_mkehkf,only: m_mkehkf_etot1,m_mkehkf_etot2
Ci   nbas  : size of basis
Ci   nsp   : number of spins
Ci   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
Ci   lmaxu : lmax for U used to dimension vorb and dmatu. lmaxu=2 if d is, but lmaxu=3 if f is included.
Ci   lldau :lldau(ib)=0 => no U on this site otherwise 
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ssite :struct for site-specific information; see m_struc_def.F
Ci   sspec :struct for species-specific information 
c      
Ci   ndham :dimensioning parameter, at least as large as largest
Ci         :hamiltonian dimension
Ci   leks  :>0 make the Hohnberg-Kohn-Sham energy
Ci         :>1 use the HKS forces
Ci   lrout : =1 generate output density and attendant quantities, =0 do not
Ci   lfrce : 0 suppress generation of forces
Ci   lpnu=1  : =1 means new pnu's
Ci   dmxp  :vector of mixing parameters; see mixrho.f for dmxp(1..25)
Ci         :Additionally:
Ci         :dmxp(33)  is the Lindhard parameter???
Ci   iter  :current iteration number
Ci   maxit :maximum number of iterations
Ci   evlall  : eigenvalues 
Ci  LDA+U inputs and outputs
Co   dmatu :density matrix for LDA+U (changed upon output)
Ci   vorb  :orbital dependent LDA+U potential
!      
Co Outputs
Co   frc   :forces.  Only calculated if lfrce>0.
Co         :If leks<2, forces are HF  forces
Co         :If leks>1, forces are HKS forces
Cl   k1,k2,k3: dimensions smrho,smpot. 
!!      nspx: number of independent spin channels
!!      nspc is now avoided (memo:nspc=2 for lso==1, nspc=1 for lso/=1 See m_lmfinit)
!!      ndhamx: maximum size of hamiltonian  
C     l   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
Cr Remarks ---- now improving...
Cr   Band pass consists of:
Cr   (1) m_mkpot_init make the effective potential,
Cr   (2) m_bandcal_init generate eigenvalues (and eigenvectors if lrout), then m_bzintegration_init
Cr   (3) m_bandcal_2nd  if lrout, assemble the output density by BZ integration
Cr   (4) evaluate hf (and KS, if leks) energy by BZ integration
Cr   (5) mixrho the output density to make a new input density.
Cu Updates before github 2009 was removed. See history github ecalj after 2009.
      implicit none
      
      include "mpif.h"
      integer ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      double precision sttime,entime

      integer:: plbnd
      logical:: llmfgw
      logical:: ltet,cmdopt0,sigmamode,tdos,debug=.false.  
      logical:: fullmesh,PROCARon,writeham=.false.
      logical:: fsmode !for --fermisurface for xcrysden.
      logical,save:: siginit=.true. 
      integer:: iter,i,ifi,ipr,iq,isp,jsp,iprint,ipts
      integer:: idummy, unlink,ifih,ifii,ib
      real(8):: sev,sumtv,eharris,eksham,  dosw(2),elind,dum,evtop,ecbot,qp(3)
      real(8):: fh_rv(3,nbas),vnow,eee,dosi(2),dee,efermxxx,emin,qvalm(2)
      real(8),allocatable:: dosi_rv(:,:),dos_rv(:,:),qmom_in(:)
      real(8),parameter::  NULLR =-99999,pi=4d0*datan(1d0)
!! ----------------------------
      call tcn ('bndfp')
c      if(master_mpi) write(stdo,"(a)") 'bndfp:start'
!  Too complicated to touch it. Probably, it is better to write new mixing routine by yourself
      debug  = cmdopt0('--debugbndfp')
      tdos   = cmdopt0('--tdos')  !total dos mode or not
      fsmode = cmdopt0('--fermisurface')! Fermi surfece for xcrysden in http://www.xcrysden.org/doc/XSF.html#2l.16
      ipr    = iprint() ! for procid/=master, we set iprint=0 at lmv7.F
      ltet = ntet>0 !tetrahedron method or not
      idummy = unlink('RSEQ_ERROR') ! Clean RSEQ_ERROR.
      
!! FermiEnergy eferm is from ctrl.*, rst.* or efermi.lmf.
      if(binit) then 
         eferm = bz_ef0 !! eferm is from ctrl file when m_bndfp_ef_set in not called (binit=T)
         binit=.false.
      endif
!  Get eferm from efermi.lmf for plbnd mode. Note that efermi.lmf can be different from efermi in rst.* file.
!  For example, you may change NKABC, and run job_pdos (lmf-MPIK --quit=band only modify efermi.lmf, without touching rst.*).
      if(plbnd/=0) then
         open(newunit=ifi,file='efermi.lmf')
         read(ifi,*) eferm
         close(ifi)             
         call mpibc1_real(eferm,1,'bndfp_eferm')
      endif
      
!! Set spin-symmetrized pnu. aug2019. See also in pnunew and locpot
      if(cmdopt0('--phispinsym')) call phispinsym_ssite_set() !pnu,pz are spin symmetrized      
!! Write out Hamiltonian HamiltonianPMT.* 
      writeham= cmdopt0('--writeham')
      if(writeham) open(newunit=ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
!! Make one-particle potential without XC part for GWdriver: lmfgw mode 
      if(llmfgw) call m_mkpot_novxc() !Get osigx,otaux oppix spotx (without XC(LDA)) !
               !call m_lmfgw_init() ! Get jobgw. Get sig,tau,ppi,spot without XC(LDA) via mkpot.
!! Generate one-body potential and energy-related quantities. See use m_mkpot:
      call m_mkpot_init(llmfgw) !from (smrho,rhoat), get one-particle potential.
                                !  mkpot->locpot->augmat. augmat calculates sig,tau,ppi.
      if(cmdopt0('--quit=mkpot')) call rx0('--quit=mkpot')
!! Setup for wtkb for BZ integration.
!! NOTE: if (wkp.* exists).and.lmet==2, wkp is used for wkkb.
      call m_subzi_init(lrout>0)
!   Hankel's e of local orbital of PZ>10 (hankel tail mode) is changing.
!   lpztail: if T, local orbital of 2nd type(hankel tail).
      if(lpztail) call sugcut(2, sspec) 
      if(cmdopt0('--cls')) then !clsmode
          call rxx(lso==1,'CLS not implemented in noncoll case') 
         if (lrout == 0) call rx('bndfp: need output density for cls')
         call m_clsmode_init()
      endif
!! === \Sigma-Vxc to ham_rv_a_ohrs. Neighbor table ham_iv_a_oiaxs ====
!! ndimsig is the dimension of the self-energy. 
!! ndimsig = ldim means we use only projection onto MTO spaces even when PMT. !aug2012takao
!! OUTPUTs by sigmainit are stored in m_seneinterp
!! In future, we should make this block to a subroutine which shows clear i/o.
      sigmamode = (lrsig/=0)
c      if( llmfgw.and.jobgw==0 ) then
c         continue
      if( sigmamode.and.siginit ) then
         call m_rdsigm2_init()
         siginit=.false. !this is equivalent to broad casting.
      endif              
      if(sigmamode.and.master_mpi) write(stdo,*)' ScaledSigma=',ham_scaledsigma
!! == GW driver and finish ==
      if(llmfgw) then
         call m_sugw_init() !jobgw) ! NOTE: Sigma-Vxc in real space is give by rv_a_ohrs.
         call tcx('bndfp')
         call rx0('sugw mode')  !exit program here normally.
      endif
!! Set up Hamiltonian and diagonalization in m_bandcal_init
!! To know outputs, see 'use m_bandcal,only:'. The outputs are evlall, and so on.
!     ! ndimhx_ is protected, but allows writing to other nodes because MPI_BCAST is f77, assumed-size array
      if(.not.(lwtkb==-1.or.lwtkb==0.or.lwtkb==1)) call rx('bndfp: something wrong lwtkb')
! lwtkb=1 accumulate weight
! lwtkb=0 no weight
! lwtkb=-1      
      sttime = MPI_WTIME()
      if(nspc==2) call m_addrbl_allocate_swtk(ndham,nsp,nkp)
      call m_bandcal_init(iqini,iqend,ispini,ispend,lrout,eferm,ifih,lwtkb) !All input. get evlall
      entime = MPI_WTIME()
      if(master_mpi) write(stdo,"(a,f9.4)") ' ... Done MPI k-loop: elapsed time=',entime-sttime
      sttime = MPI_WTIME()
      if(writeham) close(ifih)
      if(writeham) call rx0('Done --writeham: --fullmesh may be needed. HamiltonianMTO* genereted')
      
!! Broadcast. Use allreduce to avoid knowing which node to which node.
!! Because mpi is f77, ndimhx_ and so on are not protected.
      call mpibc2_int(ndimhx_,size(ndimhx_),'bndfp_ndimhx_') !all reduce (instead of which node to which node).
      call mpibc2_int(nev_,   size(nev_),   'bndfp_nev_') !all reduce 
      call mpibc2_int(nevls,  size(nevls),  'bndfp_nevls') !all reduce (to avoid which node to which node).
      call xmpbnd2(kpproc,ndhamx,nkp,nspx,evlall) !all eigenvalues broadcasted
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
      if( cmdopt0('--boltztrap')) call writeboltztrap(eferm)
      if( cmdopt0('--boltztrap')) call rx0('Done boltztrap: boltztrap.* are generated')
!! Write bands in bands-plotting case: loop over qp getting evals from array ===
      if(plbnd/=0.and.master_mpi.and.fsmode) call writefs(eferm) !fermi surface
      if(plbnd/=0.and.master_mpi) then
         write(stdo,*)' Writing bands to bands file ...'
         if(nsyml/=0) call writeband(eferm,evtop,ecbot)
      endif
      if(plbnd/=0.and.fsmom/=NULLR) then
       write(stdo,"(a)")'Caution: fsmom(fixed moment). In sc cycle, we use additional bias mag field  '
       write(stdo,"(a)")'Caution: Mag.field is written in MagField. But it is not used for --band mode!'
      endif
      if(plbnd/=0.and.fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
      if(plbnd/=0)            call rx0('plot band mode done')
!! end of band plbnd/=0, that is, band plot mode.
      
!! New eferm and wtkb determined from evlall 
      call m_subzi_bzintegration(evlall,swtk,  eferm,sev,qvalm,vnow) !Get eferm, and wtkb...
      if(lmet==0) eferm = (evtop+ecbot)/2d0 !for metal
      if(lmet==0.and. master_mpi) write(stdo,"(' HOMO; Ef; LUMO =',3f11.6)")evtop,eferm,ecbot
      if(master_mpi) then
         open(newunit=ifi,file='efermi.lmf')
         write(ifi,"(d24.16, ' # (Ry) Fermi energy given by lmf')") eferm
         write(ifi,"(d24.16, ' # (Ry) Top of Valence')") evtop
         write(ifi,"(d24.16, ' # (Ry) Bottom of conduction')") ecbot
         write(ifi,"(d24.16,d24.16 ' #before: number of electrons total at sites:qval, backg:qbg=')") qval,qbg
         write(ifi,"(d24.16, ' #band: charge(nup+nspin), mag.mom(nup-ndown)=')")qvalm(1), qvalm(2)
         write(ifi,"(3i10,'    # Used k point to determine Ef')") nkabc
         write(ifi,"('# NOTE: This file is overwritten by lmf-MPIK SC loop')")
         close(ifi)             
      endif
      
!! Generate total DOS  emin=dosw(1) emax=dosw(2) dos range
      if(master_mpi.and.(tdos.or.ldos/=0)) then
         dosw(1)= emin -0.5d0   ! lowest energy limit to plot dos
         dosw(2)= eferm+bz_dosmax ! lowest energy limit to plot dos
         write(stdo,"(a,3f9.4)") 'Generating TDOS: efermi, and dos window= ',eferm,dosw
         allocate( dosi_rv(ndos,nspx),dos_rv(ndos,nspx)) !for xxxdif
         if(cmdopt0('--tdostetf')) ltet=.false. ! Set tetrahedron=F 
         if(ltet) then 
            call bzints(nkabc(1),nkabc(2),nkabc(3), evlall 
     &           , dum , nkp , ndhamx , ndhamx , nspx , dosw(1), dosw(2), dosi_rv , ndos , efermxxx ,
     &           1, ntet , iv_a_oidtet , dum , dum ) !job=1 give IntegratedDos to dosi_rv
            call xxxdif(dosw(1),dosw(2),ndos,nspx,0,dosi_rv,dos_rv) ! 'integrated DOS dosi_rv -> DOS dos_rv'
         else
            call makdos (nkp, ndhamx, ndhamx, nspx, rv_a_owtkp, evlall,bz_n, bz_w
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
         endif
         if(lso==1) dos_rv=0.5d0*dos_rv !call dscal ( ndos , .5d0 , dos_rv , 1 )
         open(newunit=ifi, file='dos.tot.'//trim(sname) )
         open(newunit=ifii,file='dosi.tot.'//trim(sname))
         dee=(dosw(2)-dosw(1))/(ndos-1d0)
         dosi=0d0
         do ipts=1,ndos
            eee= dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-efermxxx !I think efermxxx=eferm
            dosi(1:nspx)= dosi(1:nspx) + dos_rv(ipts,1:nspx)*dee
            write(ifi,"(255(f13.5,x))") eee,  (dos_rv(ipts,isp),isp=1,nspx) !dos
            write(ifii,"(255(f13.5,x))") eee, (dosi(isp),isp=1,nspx)        !integrated dos
         enddo
         close(ifi)
         close(ifii)
      endif
      if(tdos) call rx0('Done tdos mode:')
      
!! AllReduce band quantities.  
      if(lrout/=0) call m_bandcal_allreduce(lwtkb)
      emin=1d9
      do iq = 1, nkp
         qp=qplist(:,iq)
         do isp = 1, nspx 
          jsp = isp
          if(master_mpi.and.iprint()>=35) then
            write(stdo,'(" bndfp: kpt",i5," of",i5, " k isp=",3f8.4,i2," ndimh nev=",2i5)')
     &            iq,nkp,qp,jsp,ndimhx_(iq),nevls(iq,jsp)
            write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
          endif  
          emin= min(minval( evlall(1:nevls(iq,jsp),jsp,iq)),emin)
         enddo
      enddo
      
!! Accumurate sum in BZ for wtkb
      if(lmet>0 .and. lwtkb==-1 .and. lrout>0) then 
         call m_subzi_setlwtkb(1)
         call mpi_barrier(MPI_comm_world,ierr)
         call m_bandcal_2nd(iqini,iqend,ispini,ispend,lrout, eferm) !accumulate smrho_out and so on.
         call m_bandcal_allreduce(lwtkb)
      endif
      
!! Core-level spectroscopy --- m_clsmode_set1 is called in m_bandcal
!!    dos window
      if(cmdopt0('--cls')) then
         dosw(1)= emin  - 0.5d0    ! lowest energy limit to plot dos  
         dosw(2)= eferm + bz_dosmax ! highest energy limit to plot dos 
         if (master_mpi) call m_clsmode_finalize(eferm,ndimh,ndhamx,nspx,nkp,dosw,evlall)
        call rx0('Done cls mode:')
      endif
      
!! write out orbital moment      
      if(lwtkb==1 .and. lso/=0) call iorbtm() !write out orbital moment
!! Assemble output density, energies and forces ===============================
      if (lrout/=0) call m_bandcal_symsmrho()  !Get smrho_out Symmetrize smooth density
      call m_mkrout_init() !Get force frcbandsym, symmetrized atomic densities orhoat_out, and weights hbyl,qbyl
!!  New boundary conditions pnu for phi and phidot
      if (lrout/=0) call pnunew(eferm) !ssite%pnu ssite%pz are revised.
c$$$  if(master_mpi) call writeqpyl() !if you like to print writeqbyl
      
!! Evaluate Harris-foukner energy 
!     (note: we now use NonMagneticCORE mode as default)
      call m_mkehkf_etot1(sev, eharris)
      
!! Evaluate KS total energy and Force ---
      if(lrout/=0) then
        allocate(qmom_in(nvl))
        qmom_in=qmom !
        elind = ham_elind ! accelarating convergence. Not used so often. Need examination.
        if(ham_elind<0d0) elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*ham_elind
!!   ... Correction to harris force ! qmom is given by m_mkpot_init for input dentisy.
c        if(lfrce>0 ) call dfrce (lfrce,orhoat,orhoat_out,elind,qmom,osmrho,smrho_out,         fh_rv)
!!   ... Evaluate KS total energy and output magnetic moment
        eksham = 0d0
        if(leks>=1) then
           call togpr()
           call mkekin(sv_p_osig,sv_p_otau,sv_p_oppi,sv_p_oqkkl,vconst,osmpot,smrho_out,sev,  sumtv)
           call m_mkpot_energyterms(smrho_out, orhoat_out) !qmom is revised for given orhoat_out
           if(cmdopt0('--density')) then
              call mpi_barrier(MPI_comm_world,ierr)  
              call rx0('end of --density mode')
           endif  
           call m_mkehkf_etot2(sev,sumtv, eksham) 
           call togpr()
        endif
!! --- Add together force terms ---
!! fes1_rv: contribution to HF forces from estat + xc potential
!           This is for input  density !=3rd term in (B.5) in JPSJ.84.034705
!! fes2_rv: contribution to KS forces from estat + xc potential.   This is for output density
!! frcbandsym : 1st term in (B.5)  (puley? need check)
!! fh_rv      : 2nd term in (B.5)  (need check)
!! force : total
        if(lfrce> 0) then
         if(allocated(force)) deallocate(force)
         allocate(force(3,nbas))
         call dfrce (lfrce,orhoat,orhoat_out,elind,qmom_in,osmrho,smrho_out,  fh_rv)
         call totfrc(leks, fes1_rv, fes2_rv, fh_rv, frcbandsym, force)
        endif
        deallocate(qmom_in)
!! Mixing input (osmrho,orhoat) and output (osmrho_out,orhoat_out). Output are orhoat and osmrho.
        call mixrho(iter,qval-qbg,elind, orhoat_out,orhoat, smrho_out,osmrho, qdiff) !mixrho keeps history in it.
      else
        eksham = 0d0
      endif
      ham_ehf= eharris !Harris total energy
      ham_ehk= eksham  !Hohenberg-Kohn-Sham total energy
      call m_mkpot_deallocate()
      call m_bandcal_clean()
      call tcx('bndfp')
      end subroutine bndfp


      subroutine writeqpyl()
      use m_ext,only: sname
      use m_lmfinit,only:nbas,ssite=>v_ssite,sspec=>v_sspec,nsp
      use m_mkrout,only: qbyl_rv
      integer::ib,ispec,ifqbyl,il
      open(newunit=ifqbyl,file='qbyl.'//trim(sname))
      write(ifqbyl,"(i5)") nbas
      do  ib = 1, nbas
         ispec=ssite(ib)%spec
         write(ifqbyl,"(i4, i3,i2,10f12.6)")
     &   ib,ispec,sspec(ispec)%lmxa,( sum(qbyl_rv(il,1:nsp,ib)), il=1,sspec(ispec)%lmxa+1)
      enddo  
      close(ifqbyl)
      end
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
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

      subroutine xxxdif(emin,emax,npts,ndos,mode,dosi,dos)
!! Differentiate Number of States to make DOS
      implicit none
      integer npts,ndos,mode
      real(8):: dosi(npts,ndos), dos(npts,ndos),emin,emax
      integer ip,idos,k
      double precision bin
      bin = 2*(emax-emin)/dble(npts-1)
      do  idos = 1, ndos
        do ip = 2, npts-1
           dos(ip,idos) = (dosi(ip+1,idos) - dosi(ip-1,idos)) / bin
        enddo   
        dos(1,idos)    = dos(2,idos)
        dos(npts,idos) = dos(npts-1,idos)
      enddo  
      end
      
      end module m_bndfp

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine xmpbnd2(kpproc,ndham,nkp,nspx,eb)
      use m_lmfinit, only: stdo,stdl,procid,master
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq
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
      integer:: kpproc(0:*),ndham,nkp,nspx
      double precision eb(ndham,*) 
      integer i,ista,iend, ierr
      integer, dimension(:),allocatable :: offset,length
      real(8) ,allocatable :: buf_rv(:)
      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0
      do  i = 0, numprocs-1
        ista = kpproc(i)
        iend = kpproc(i+1)-1
        length(i) = (iend - ista + 1)*ndham !nsp*ndham
        offset(i+1) = offset(i) + length(i)
      enddo
      ista = kpproc(procid)
      allocate(buf_rv(ndham*nkp*nspx))
      call mpi_allgatherv ( eb(1,ista) , length( procid ), mpi_double_precision, buf_rv , length , offset 
     ., mpi_double_precision , mpi_comm_world , ierr )
      call dcopy ( ndham * nspx * nkp , buf_rv , 1 , eb , 1 )
      if (allocated(buf_rv)) deallocate(buf_rv)
      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)
      end subroutine xmpbnd2

      
