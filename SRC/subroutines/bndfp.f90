!     ! = Main part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration =
!     !
!     ! From density osmrhod and osmrho, m_mkpot_init gives  osmpot, and sv_p_osig, sv_p_otau, sv_p_oppi (potential)
!     ! Then we do band calcultion m_bandcal_init. It gives osmrho_out and orhoat_out.
!     ! Finally, osmrho and orhoat are returned after mixing procedure (mixrho).
!----------------------------
! o  smrho :smooth density
! o        :Symmetrized on output
! o  orhoat:vector of offsets containing site density
! o        :Symmetrized on output
!  qbyl  :site- and l-decomposed charges
!  hbyl  :site- and l-decomposed one-electron energies
!o   dmatu:   density matrix in m_bandcal
!o   evalall: eigenvalue in m_bandcal
!o   force:
!     i  vorb --> use m_ldau in m_mkpot
!!--- key quntities --------------
!  evlall  : eigenvalues
!  eferm   : Fermi energy
!  sev     : sum of band energy
!     eksham   : DFT(Kohn-Sham)  total energy
!     eharris  : Harris foulkner total energy.
!        eharris is not correct for LDA+U case (valv locpot2 do not include U contribution).
module m_bndfp
  use m_ftox
  use m_density,only: orhoat,osmrho !input/output unprotected
  use m_mixrho,only: Mixrho
  real(8),protected:: ham_ehf, ham_ehk, sev  !output
  real(8),protected:: eferm, qdiff  !output
  real(8),protected,allocatable:: force(:,:)
  !! other ouput are in modules m_mkpot, m_bandcal
  logical,private:: binit=.true.,initd=.true.
contains

  subroutine m_bndfp_ef_set(bz_ef00) !called from iors.F. grep 'use m_bndfp'
    real(8):: bz_ef00
    eferm  = bz_ef00
    binit=.false.
  end subroutine m_bndfp_ef_set

  subroutine bndfp(iter, llmfgw, plbnd)
    ! llmfgw=T is for generating eigenfunctions for GW calculations, no iteration.
    ! plbnd/=0 means band plot mode. no iteration.
    !     ! All read only in bndfp. Data are stored in modules such as m_bandcal, m_mkpot
    !     ! For example,, rightafter call m_bandcal_init, we can get evalall, which is used in other modules.
    use m_supot,only: ngabc=>lat_nabc,k1,k2,k3 !for charge mesh
    use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx,nspx=>ham_nspx
    use m_lmfinit, only: n0,nab,nppn,ncutovl,lso,ndos=>bz_ndos,bz_w,fsmom=>bz_fsmom, &
    bz_dosmax,lmet=>bz_lmet,bz_fsmommethod,bz_n, &
         ctrl_nspec,ctrl_pfloat,ldos=>ctrl_ldos,qbg=>zbak,lfrce=>ctrl_lfrce, &
         pwmode=>ham_pwmode,lrsig=>ham_lsig,epsovl=>ham_oveps, &
         ham_scaledsigma, &
         alat=>lat_alat,stdo,stdl,procid,master, &! & bz_doswin,
    nkaph,nlmax,nl,nbas,nsp, ham_frzwf, bz_dosmax, &
         lekkl,lmaxu,nlibu,lldau,lpztail,leks,lrout &
         ,  nchan=>pot_nlma, nvl=>pot_nlml,nspc
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
    use m_rdsigm2,only: M_rdsigm2_init
    use m_subzi, only: M_subzi_init, lwtkb,rv_a_owtkb,M_subzi_setlwtkb,M_subzi_bzintegration
    use m_rsibl, only : Rsibl_ev ! to plot wavefunction in the fat band mode
    use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq,xmpbnd2
    !      use m_lmfgw,only: M_lmfgw_init !,jobgw !,sv_p_osigx,sv_p_otaux,sv_p_oppix,spotx
    use m_mkpot,only: M_mkpot_init,M_mkpot_deallocate, M_mkpot_energyterms,M_mkpot_novxc,& ! & M_mkpot_novxc_dipole,
    osmpot, qmom, vconst, osig,otau,oppi, qval , qsc , fes1_rv , fes2_rv
    use m_clsmode,only: M_clsmode_init,m_clsmode_set1,m_clsmode_finalize
    use m_qplist,only:  qplist,nkp,xdatt,labeli,labele,dqsyml,etolc,etolv, &
         nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml, &
         iqini,iqend,ispini,ispend,kpproc    ! MPIK divider. iqini:iqend are node-dependent
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x
    use m_procar,only: M_procar_init,dwgtall,nchanp,m_procar_closeprocar,m_procar_writepdos
    use m_bandcal,only: M_bandcal_init,M_bandcal_2nd,M_bandcal_clean,M_bandcal_allreduce, &
         evlall,smrho_out,sv_p_oqkkl,sv_p_oeqkkl, ndimhx_,nev_,nevls,m_bandcal_symsmrho
    use m_mkrout,only: M_mkrout_init,orhoat_out,frcbandsym,hbyl_rv,qbyl_rv
    use m_addrbl,only: M_addrbl_allocate_swtk,swtk
    use m_sugw,only: M_sugw_init
    use m_mkehkf,only: M_mkehkf_etot1,M_mkehkf_etot2
    use m_gennlat_sig,only: M_gennlat_init_sig
    use m_dfrce,only: dfrce
    use m_sugcut,only:sugcut
    !i   nbas  : size of basis
    !i   nsp   : number of spins
    !i   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
    !i   lmaxu : lmax for U used to dimension vorb and dmatu. lmaxu=2 if d is, but lmaxu=3 if f is included.
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat in dmats(*,lldau(ib))
    !i   ssite :struct for site-specific information; see m_struc_def.F
    !i   sspec :struct for species-specific information

    !i   ndham :dimensioning parameter, at least as large as largest
    !i         :hamiltonian dimension
    !i   leks  :>0 make the Hohnberg-Kohn-Sham energy
    !i         :>1 use the HKS forces
    !i   lrout : =1 generate output density and attendant quantities, =0 do not
    !i   lfrce : 0 suppress generation of forces
    !i   lpnu=1  : =1 means new pnu's
    !i   dmxp  :vector of mixing parameters; see mixrho.f for dmxp(1..25)
    !i         :Additionally:
    !i         :dmxp(33)  is the Lindhard parameter???
    !i   iter  :current iteration number
    !i   maxit :maximum number of iterations
    !i   evlall  : eigenvalues
    !i  LDA+U inputs and outputs
    !o   dmatu :density matrix for LDA+U (changed upon output)
    !i   vorb  :orbital dependent LDA+U potential

    !o Outputs
    !o   frc   :forces.  Only calculated if lfrce>0.
    !o         :If leks<2, forces are HF  forces
    !o         :If leks>1, forces are HKS forces
    !l   k1,k2,k3: dimensions smrho,smpot.
    !!      nspx: number of independent spin channels
    !!      nspc is now avoided (memo:nspc=2 for lso==1, nspc=1 for lso/=1 See m_lmfinit)
    !!      ndhamx: maximum size of hamiltonian
    !     l   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
    !r Remarks ---- now improving...
    !r   Band pass consists of:
    !r   (1) m_mkpot_init make the effective potential,
    !r   (2) m_bandcal_init generate eigenvalues (and eigenvectors if lrout), then m_bzintegration_init
    !r   (3) m_bandcal_2nd  if lrout, assemble the output density by BZ integration
    !r   (4) evaluate hf (and KS, if leks) energy by BZ integration
    !r   (5) mixrho the output density to make a new input density.
    !u Updates before github 2009 was removed. See history github ecalj after 2009.
    implicit none

    include "mpif.h"
    integer :: ierr, status(MPI_STATUS_SIZE)
    integer :: MAX_PROCS
    parameter (MAX_PROCS = 100)
    integer :: resultlen
    character*(MPI_MAX_PROCESSOR_NAME) name
    character(10) :: shortname(0:MAX_PROCS-1)
    double precision :: sttime,entime

    integer:: plbnd,nk1,nk2,nk3,nx,ny
    logical:: llmfgw,sigx
    logical:: ltet,cmdopt0,sigmamode,tdos,debug=.false.
    logical:: fullmesh,PROCARon,writeham=.false.
    logical:: fsmode !for --fermisurface for xcrysden.
    logical,save:: siginit=.true.
    integer:: iter,i,ifi,ipr,iq,isp,jsp,iprint,ipts
    integer:: idummy, unlink,ifih,ifii,ib
    real(8):: sumtv,eharris,eksham,  dosw(2),dum,evtop,ecbot,qp(3),rydberg,xxx,eeem
    real(8):: fh_rv(3,nbas),vnow,eee,dosi(2),dee,efermxxx,emin,qvalm(2)
    real(8),allocatable:: dosi_rv(:,:),dos_rv(:,:),qmom_in(:)
    real(8),parameter::  NULLR =-99999,pi=4d0*datan(1d0)
    real(8):: elind=0d0
    !! ----------------------------
    call tcn ('bndfp')
    debug  = cmdopt0('--debugbndfp')
    tdos   = cmdopt0('--tdos')  !total dos mode or not
    fsmode = cmdopt0('--fermisurface')!FermiSurfece for xcrysden in http://www.xcrysden.org/doc/XSF.html#2l.16
    ipr    = iprint() ! for procid/=master, we set iprint=0 at lmv7.F
    ltet = ntet>0 !tetrahedron method or not
    !! FermiEnergy eferm is from rst.* or efermi.lmf.
    if(binit) then
       eferm = 0d0 ! if no eferm is read from rst file.
       binit=.false.
    endif
    !  Get eferm from efermi.lmf for plbnd mode.
    !     efermi.lmf can be different from efermi in rst.* file.
    !     For example, you may change NKABC, and run job_pdos (lmf-MPIK --quit=band only modify efermi.lmf, without touching rst.*).
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
    !      if(llmfgw.and.cmdopt0('--dipolematrix')) call m_mkpot_novxc_dipole()
    !! Generate one-body potential and energy-related quantities. See use m_mkpot:
    call m_mkpot_init(llmfgw) !from (smrho,rhoat), get one-particle potential.
    !  mkpot->locpot->augmat. augmat calculates sig,tau,ppi.
    if(cmdopt0('--quit=mkpot')) call rx0('--quit=mkpot')
    !! Setup for wtkb for BZ integration.
    !! NOTE: if (wkp.* exists).and.lmet==2, wkp is used for wkkb.
    call m_subzi_init(lrout>0)
    !   Hankel's e of local orbital of PZ>10 (hankel tail mode) is changing.
    !   lpztail: if T, local orbital of 2nd type(hankel tail).
    if(lpztail) call sugcut(2)
    if(cmdopt0('--cls')) then !clsmode
       call rxx(lso==1,'CLS not implemented in noncoll case')
       if (lrout == 0) call rx('bndfp: need output density for cls')
       call m_clsmode_init()
    endif
    !! === \Sigma-Vxc to hrr
    !! ndimsig is the dimension of the self-energy.
    !  Current version is "ndimsig = ldim", which means we use only projection
    !  onto MTO spaces even when PMT.
    !! OUTPUT of m_rdsigm2_init: ndimsig, nk1,nk2,nk3, and hrr.
    !! These are used in getsenex via bloch2.
    sigmamode = (lrsig/=0)
    if( sigmamode .AND. siginit ) then
       inquire(file='sigm.'//trim(sname),exist=sigx)
       if(sigx) then
          open(newunit=ifi,file='sigm.'//trim(sname),form='unformatted')
          read(ifi) nx,ny,nk1,nk2,nk3
          close(ifi)
          call m_gennlat_init_sig([nk1,nk2,nk3]) !nlat(real-space lattice vectors) for FFT of sigm
       endif
       call m_rdsigm2_init()
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       siginit=.false. !only once
    endif
    if(sigmamode .AND. master_mpi) write(stdo,*)' ScaledSigma=',ham_scaledsigma
    !! == GW driver and finish ==
    if(llmfgw) then
       !         call m_sugw_init(cmdopt0('--dipolematrix'),cmdopt0('--socmatrix'),eferm)
       call m_sugw_init(cmdopt0('--socmatrix'),eferm)
       call tcx('bndfp')
       call rx0('sugw mode')  !exit program here normally.
    endif
    !! Set up Hamiltonian and diagonalization in m_bandcal_init
    !! To know outputs, see 'use m_bandcal,only:'. The outputs are evlall, and so on.
    !! ndimhx_ is protected, but allows writing to other nodes because MPI_BCAST is f77, assumed-size array
    if( .NOT. (lwtkb==-1 .OR. lwtkb==0 .OR. lwtkb==1)) call rx('bndfp: something wrong lwtkb')
    ! lwtkb=1 accumulate weight,
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
    eeem=9999
    do iq=1,nkp
       do jsp=1,nspx          !nspx=1 for SO=1
          if(lso==1) i = max(1,nint(qval-qbg))
          if(lso/=1) i = max(1,nint(qval-qbg)/2)
          evtop = max(evtop,evlall(i,jsp,iq))
          ecbot = min(ecbot,evlall(i+1,jsp,iq))
          eeem = min(eeem,evlall(1,jsp,iq))
          if(master_mpi .AND. iq==1)write(stdl,"('fp evl',8f8.4)")(evlall(i,jsp,iq),i=1,nev_(iq))
       enddo
    enddo
    !! pdos mode (--mkprocar and --fullmesh)
    fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface')
    PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
    if(fullmesh .AND. procaron) call m_procar_writepdos(evlall,nev_,eferm,kpproc)
    if(fullmesh .AND. procaron) call rx0('Done pdos: --mkprocar & --fullmesh. Check by "grep k-point PROCAR.*.*"')
    !! boltztrap data
    if( cmdopt0('--boltztrap')) call writeboltztrap(eferm)
    if( cmdopt0('--boltztrap')) call rx0('Done boltztrap: boltztrap.* are generated')
    !! Write bands in bands-plotting case: loop over qp getting evals from array ===
    if(plbnd/=0 .AND. master_mpi .AND. fsmode) call writefs(eferm) !fermi surface
    if(plbnd/=0 .AND. master_mpi) then
       write(stdo,*)' Writing bands to bands file ...'
       if(nsyml/=0) call writeband(eferm,evtop,ecbot)
    endif
    if(plbnd/=0 .AND. fsmom/=NULLR) then
       write(stdo,"(a)")'Caution: fsmom(fixed moment). In sc cycle, we use additional bias mag field  '
       write(stdo,"(a)")'Caution: Mag.field is written in MagField. But it is not used for --band mode!'
    endif
    if(plbnd/=0 .AND. fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
    if(plbnd/=0) call rx0('plot band mode done') ! end of band plbnd/=0, that is, band plot mode.
    !! New eferm and wtkb determined from evlall
    call m_subzi_bzintegration(evlall,swtk,  eferm,sev,qvalm,vnow) !Get eferm and wtkb in m_subzi
    if(lmet==0) eferm = (evtop+ecbot)/2d0 !for metal
    if(lmet==0 .AND. master_mpi) write(stdo,"(' HOMO; Ef; LUMO =',3f11.6)")evtop,eferm,ecbot
    if(master_mpi) then
       open(newunit=ifi,file='efermi.lmf')
       write(ifi,"(d24.16, ' # (Ry) Fermi energy given by lmf')") eferm
       write(ifi,"(d24.16, ' # (Ry) Top of Valence')") evtop
       write(ifi,"(d24.16, ' # (Ry) Bottom of conduction')") ecbot
       write(ifi,"(2d24.16,' #before: number of electrons total at sites:qval, backg:qbg=')") qval,qbg
       write(ifi,"(2d24.16,' # band: charge(nup+nspin), mag.mom(nup-ndown)')")qvalm(1), qvalm(2)
       write(ifi,"(f24.9 , ' # band: band energy sum (eV)=')") sev*rydberg()
       write(ifi,"(3i10,'    # Used k point to determine Ef')") nkabc
       write(ifi,"(i6,'# iter CAUTION! This file is overwritten by lmf-MPIK SC loop')")iter
       close(ifi)
    endif
    !! Generate total DOS  emin=dosw(1) emax=dosw(2) dos range
    if(master_mpi .AND. (tdos .OR. ldos/=0)) then
       emin = eeem-0.01d0 
       dosw(1)= emin ! lowest energy limit to plot dos
       dosw(2)= eferm+bz_dosmax
       write(stdo,"(a,3f9.4)") 'Generating TDOS: efermi, and dos window= ',eferm,dosw
       allocate( dosi_rv(ndos,nspx),dos_rv(ndos,nspx)) !for xxxdif
       if(cmdopt0('--tdostetf')) ltet= .FALSE. ! Set tetrahedron=F
       if(ltet) then
          call bzints(nkabc(1),nkabc(2),nkabc(3), evlall &
               , dum , nkp , ndhamx , ndhamx , nspx , dosw(1),dosw(2), dosi_rv , ndos ,xxx , &
               1, ntet , iv_a_oidtet , dum , dum ) !job=1 give IntegratedDos to dosi_rv
          dos_rv(2:ndos-1,:)=(dosi_rv(3:ndos,:)-dosi_rv(1:ndos-2,:))/(2d0*(dosw(2)-dosw(1))/(ndos-1))
          dos_rv(1,:)    = dos_rv(2,:)
          dos_rv(ndos,:) = dos_rv(ndos-1,:)
       else
          call makdos (nkp, ndhamx, ndhamx, nspx, rv_a_owtkp, evlall,bz_n, bz_w &
               , -6d0, dosw(1),dosw(2), ndos , dos_rv )
       endif
       if(lso==1) dos_rv=0.5d0*dos_rv !call dscal ( ndos , .5d0 , dos_rv , 1 )
       open(newunit=ifi, file='dos.tot.'//trim(sname) )
       open(newunit=ifii,file='dosi.tot.'//trim(sname))
       dee=(dosw(2)-dosw(1))/(ndos-1d0)
       dosi=0d0
       do ipts=1,ndos
          eee= dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-eferm
          dosi(1:nspx)= dosi(1:nspx) + dos_rv(ipts,1:nspx)*dee
          write(ifi,"(255(f13.5,x))") eee,  (dos_rv(ipts,isp),isp=1,nspx) !dos
          write(ifii,"(255(f13.5,x))") eee, (dosi_rv(ipts,isp),isp=1,nspx)        !integrated dos
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
          if(master_mpi .AND. iprint()>=35) then
             write(stdo,ftox)" bndfp: kpt",iq," of",nkp," k isp=",ftof(qp,4),jsp," ndimh nev=", &
                  ndimhx_(iq),nevls(iq,jsp)
             write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
          endif
          emin= min(minval( evlall(1:nevls(iq,jsp),jsp,iq)),emin)
       enddo
    enddo
    !! Accumurate sum in BZ for wtkb
    if(lmet>0 .AND. lwtkb==-1 .AND. lrout>0) then
       call m_subzi_setlwtkb(1)
       call mpi_barrier(MPI_comm_world,ierr)
       call m_bandcal_2nd(iqini,iqend,ispini,ispend,lrout)!, eferm) !accumulate smrho_out and so on.
       call m_bandcal_allreduce(lwtkb)
    endif
    !! Core-level spectroscopy --- m_clsmode_set1 is called in m_bandcal
    if(cmdopt0('--cls')) then
       dosw(1)= emin  - 0.5d0    ! lowest energy limit to plot dos
       dosw(2)= eferm + bz_dosmax ! highest energy limit to plot dos
       if (master_mpi) call m_clsmode_finalize(eferm,ndimh,ndhamx,nspx,nkp,dosw,evlall)
       call rx0('Done cls mode:')
    endif
    !! write out orbital moment
    if(lwtkb==1 .AND. lso/=0) call iorbtm()  ! write orbital moment
    !! Assemble output density, energies and forces ===============================
    if (lrout/=0) call m_bandcal_symsmrho()  !Get smrho_out Symmetrize smooth density
    call m_mkrout_init() !Get force frcbandsym, symmetrized atomic densities orhoat_out, and weights hbyl,qbyl
    !!  New boundary conditions pnu for phi and phidot
    if (lrout/=0) call pnunew(eferm) !ssite%pnu ssite%pz are revised.
    !  if(master_mpi) call writeqpyl() !if you like to print writeqbyl
    !! Evaluate Harris-foukner energy (note: we now use NonMagneticCORE mode as default)
    call m_mkehkf_etot1(sev, eharris)
    !! Evaluate KS total energy and Force ---
    if(lrout/=0) then
       allocate(qmom_in(nvl))
       qmom_in=qmom !multipole moments.
       eksham = 0d0 !   ... Evaluate KS total energy and output magnetic moment
       if(leks>=1) then
          call mkekin(osig,otau,oppi,sv_p_oqkkl,vconst,osmpot,smrho_out,sev,  sumtv)
          call m_mkpot_energyterms(smrho_out, orhoat_out) !qmom is revised for given orhoat_out
          if(cmdopt0('--density')) then
             call mpi_barrier(MPI_comm_world,ierr)
             call rx0('end of --density mode')
          endif
          call m_mkehkf_etot2(sev,sumtv, eksham)
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
          call dfrce (lfrce,orhoat,orhoat_out,qmom_in,osmrho,smrho_out,  fh_rv)
          call totfrc(leks, fes1_rv, fes2_rv, fh_rv, frcbandsym, force)
       endif
       deallocate(qmom_in)
       ! Mix inputs(osmrho,orhoat) and outputs(osmrho_out,orhoat_out), resulting orhoat and osmrho.
       call mixrho(iter,qval-qbg,orhoat_out,orhoat,smrho_out,osmrho,qdiff)!mixrho keeps history in it.
    else
       eksham = 0d0
    endif
    ham_ehf= eharris !Harris total energy
    ham_ehk= eksham  !Hohenberg-Kohn-Sham total energy
    call m_mkpot_deallocate()
    call m_bandcal_clean()
    call tcx('bndfp')
  end subroutine bndfp
end module m_bndfp

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!      subroutine readindensitymodesetting()
!! If necessary, need to recover this but it should be much simpler.
!! By chaningt wtkb, this will give density for given iq ib.
!!  We have electron density plot by '--density' mode
!$$$!! Density mode. wtkb is the weight for iq,ib. With this setting,
!$$$!! we can plot |\psi_{iq,ib}(\bfr)|**2 in the manner of electron density
!$$$!     ! lmf --density,iq=12,ib=3,4,5   (here ib=3,4,5 are superposed).
!$$$      if(cmdopt('--density',9,0,strn)) then
!$$$         iqread=0
!$$$         iqindex=index(strn(10:),'iq=')+2
!$$$         if(iqindex==2) then
!$$$            iqread=-999
!$$$         else
!$$$            read(strn(10+iqindex:),*) iqread
!$$$         endif
!$$$         ibindex=index(strn(10:),'ib=')+2
!$$$         if(ibindex==2) then
!$$$            ibread=-999
!$$$         else
!$$$            do nibread=1,100
!$$$               read(strn(10+ibindex:),*,err=2019,end=2019) ibread(1:nibread)
!$$$c     print *,'xxx ibread=',ibread(1:nibread)
!$$$            enddo
!$$$ 2019       continue
!$$$            nibread=nibread-1
!$$$         endif
!$$$         if(iqread>0) then
!$$$            if(maxval(ibread(1:nibread))<=0) call rx('--density mode: wrong ib=foobar. Try,e.g. --density,iq=12,ib=5')
!$$$            write(stdo,"('--density bandmode: psi**2 for iq=',i5,' ib=',255i5)") iqread,ibread(1:nibread)
!$$$            rv_a_owtkb(:,:,:)=0d0
!$$$            do ib=1,nibread
!$$$               rv_a_owtkb(ibread(ib),:,iqread)=1d0
!$$$            enddo
!$$$            do ib=1,ndimh
!$$$               if(abs(rv_a_owtkb(ib,jsp,iq))>1d-3) write(stdo,"('ib wtkb=',i5,2f13.6)") ib,rv_a_owtkb(ib,jsp,iq)
!$$$            enddo
!$$$         endif
!$$$      endif
!     end
!end module m_bndfp

