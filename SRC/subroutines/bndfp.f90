!= Main part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration = 2023-jan memo checked
! From density osmrhod and osmrho, m_mkpot_init gives potential,
! osmpot, and sv_p_osig, sv_p_otau, sv_p_oppi (potential site integrals)
! Then we do band calcultion m_bandcal_init. It gives osmrho_out and orhoat_out.
! Finally, osmrho and orhoat are returned after mixing procedure (mixrho).
!-- i/o --------------------------
!o  smrho :smooth density
!o        :Symmetrized on output
!o  orhoat:vector of offsets containing site density
!o        :Symmetrized on output
!o   dmatu:   density matrix in m_bandcal
!o   evalall: eigenvalue in m_bandcal
!o   force:
!i  vorb --> we use m_ldau in m_mkpot
!--- key quntities --------------
!  evlall  : eigenvalues
!  eferm   : Fermi energy
!  sev     : sum of band energy
!     eksham   : DFT(Kohn-Sham)  total energy
!     eharris  : Harris foulkner total energy.
!       NOTE=>eharris is not correct for LDA+U case (valv locpot2 do not include U contribution).
!  qbyl  :site- and l-decomposed charges
!  hbyl  :site- and l-decomposed one-electron energies
!! NOTE: check lmv7->lmfp(iteration loop)->bndfp
module m_bndfp
  use m_ftox
  use m_density,only: orhoat,osmrho !input/output unprotected
  use m_mixrho,only: Mixrho
  real(8),protected,public:: ham_ehf, ham_ehk, sev  !output
  real(8),protected,public:: eferm, qdiff           !output
  real(8),protected,allocatable,public:: force(:,:) !output
  ! NOTE: other ouput of bndfp are stored in modules:
  ! m_mkpot(potential), m_bandcal(band,dmatu), and m_density(density)
  public bndfp, m_bndfp_ef_set
  private
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
         nkaph,nlmax,nl,nbas,nsp, bz_dosmax, &
         lekkl,lmaxu,nlibu,lldau,lpztail,leks,lrout &
         ,  nchan=>pot_nlma, nvl=>pot_nlml,nspc,pnufix
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
    !o         :If leks=1, forces are HF  forces
    !o         :If leks=2, forces are HKS forces
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
    integer:: idummy, unlink,ifih,ifii,ib,ix,ifimag
    real(8):: ekinval,eharris,eksham,  dosw(2),dum,evtop,ecbot,qp(3),rydberg,xxx,eeem
    real(8):: fh_rv(3,nbas),vnow_magfield,eee,dosi(2),dee,efermxxx,emin,qvalm(2)
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
    call m_mkpot_init() !from (smrho,rhoat), get one-particle potential.
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
    do iq=1,nkp
       do jsp=1,nspx     !nspx=1 for SO=1
          do ix=1,ndhamx !2022dec evtop,ecbot for magnetic case
             if(evlall(ix,jsp,iq)<eferm) evtop = max(evtop,evlall(ix,jsp,iq))
             if(evlall(ix,jsp,iq)>eferm) ecbot = min(ecbot,evlall(ix,jsp,iq))
          enddo
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
    if(plbnd/=0 ) then
       if(master_mpi .AND. fsmode) call writefs(eferm) !fermi surface
       if(master_mpi) then
          write(stdo,*)' Writing bands to bands file for gnuplot ...'
          if(nsyml/=0) call writeband(eferm,evtop,ecbot)
       endif
       if(fsmom/=NULLR) then
          write(stdo,"(a)")'Caution: fsmom(fixed moment). In sc cycle, we use additional bias mag field  '
          write(stdo,"(a)")'Caution: Mag.field is written in MagField. But it is not used for --band mode!'
       endif
       if(fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
       call rx0('plot band mode done') ! end of band plbnd/=0, that is, band plot mode.
    endif
    
    !! New eferm and wtkb determined from evlall
    call m_subzi_bzintegration(evlall,swtk, eferm,sev,qvalm,vnow_magfield) !Get the Fermi energy eferm,...
    if(fsmom/=NULLR .AND. master_mpi) then !moved from m_subzi_bzintegration at 2022-dec
       open(newunit=ifimag,file='MagField')
       write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow_magfield 
       close(ifimag)
    endif
    ! -vnow_magfield/2 is added to isp=1, +vnow_magfield/2 to isp=2
    evtop=-9999
    ecbot=9999
    eeem=9999
    do iq=1,nkp
       do jsp=1,nspx     !nspx=1 for SO=1
          do ix=1,ndhamx !2022dec evtop,ecbot for magnetic case
             if(evlall(ix,jsp,iq)<eferm) evtop = max(evtop,evlall(ix,jsp,iq))
             if(evlall(ix,jsp,iq)>eferm) ecbot = min(ecbot,evlall(ix,jsp,iq))
          enddo
       enddo
    enddo
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
       write(stdo,ftox)' bndfp:Generating TDOS: efermi=',ftof(eferm),' dos window emin emax= ',ftof(dosw)
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
          call makdos(nkp,ndhamx,ndhamx,nspx,rv_a_owtkp,evlall,bz_n,bz_w,-6d0,dosw(1),dosw(2),ndos,dos_rv)
       endif
       if(lso==1) dos_rv=0.5d0*dos_rv 
       open(newunit=ifi, file='dos.tot.'//trim(sname) )
       open(newunit=ifii,file='dosi.tot.'//trim(sname))
       dee=(dosw(2)-dosw(1))/(ndos-1d0)
       dosi=0d0
       do ipts=1,ndos
          eee= dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-eferm
          dosi(1:nspx)= dosi(1:nspx) + dos_rv(ipts,1:nspx)*dee
          write(ifi,"(255(f13.5,x))")  eee,(dos_rv(ipts,isp),isp=1,nspx) ! dos
          write(ifii,"(255(f13.5,x))") eee,(dosi_rv(ipts,isp),isp=1,nspx)! integrated dos
       enddo
       close(ifi)
       close(ifii)
    endif
    if(tdos) call rx0('Done tdos mode:')

    !! AllReduce band quantities.
    if(lrout/=0) call m_bandcal_allreduce(lwtkb)
    emin=1d9
    if(master_mpi) write(stdo,ftox)
    do iq = 1,nkp
       qp=qplist(:,iq)
       do isp = 1,nspx
          jsp = isp
          if(master_mpi .AND. iprint()>=35) then
             write(stdo,ftox)" bndfp: kpt",iq," of",nkp," k isp=",ftof(qp,4),jsp," ndimh nev=",&
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
    if(.not.PNUFIX.and.lrout/=0) call pnunew(eferm) !ssite%pnu ssite%pz are revised.
    !  if(master_mpi) call writeqpyl() !if you like to print writeqbyl
    !! Evaluate Harris-foukner energy (note: we now use NonMagneticCORE mode as default)
    call m_mkehkf_etot1(sev, eharris)
    !! Evaluate KS total energy and Force ---
    if(lrout/=0) then
       allocate(qmom_in(nvl))
       qmom_in=qmom !multipole moments.
       eksham = 0d0 !   ... Evaluate KS total energy and output magnetic moment
       if(leks>=1) then
          call mkekin(osig,otau,oppi,sv_p_oqkkl,vconst,osmpot,smrho_out,sev,  ekinval)
          call m_mkpot_energyterms(smrho_out, orhoat_out) !qmom is revised for given orhoat_out
          if(cmdopt0('--density')) then
             call mpi_barrier(MPI_comm_world,ierr)
             call rx0('end of --density mode')
          endif
          call m_mkehkf_etot2(ekinval, eksham)
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
!========================================================  
  subroutine mkekin(sv_p_osig,sv_p_otau,sv_p_oppi,sv_p_oqkkl,vconst,smpot,smrho,sumev, ekinval)
    use m_struc_def
    use m_lmfinit,only:nkaph,nsp,nspc,stdo,nbas,ispec,sspec=>v_sspec,nlmto
    use m_lattic,only: lat_vol
    use m_supot,only: lat_nabc,k1,k2,k3
    use m_orbl,only: Orblib,ktab,ltab,offl,norb,ntab,blks
    !- Evaluate the valence kinetic energy
    !NOTE: When SOC included, ek contains HSO, because ek= Eband- V*n where Eband contains SO contr.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   nlmto :dimension of hamiltonian (MTO)
    !i   lcplxp=1 only now: 0 if ppi is real; 1 if ppi is complex
    !i   osig  :augmentation overlap integrals
    !i   otau  :augmentation kinetic energy integrals
    !i   oppi  :augmentation kinetic + potential integrals
    !i   oqkkl :local density-matrix; see rlocbl.f
    !i   k1..3 :dimensions smpot,smrho
    !i   vconst:constant potential added to hamiltonian
    !i   smpot :smooth input potential on uniform mesh (mkpot.f)
    !i         :smpot = Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
    !i   smrho :smooth output density n0 (no gaussians n0~-n0)
    !i   sumev :sum of eigenvalues
    !o Outputs
    !o   ekinval :kinetic energy
    !l Local variables
    !l   sraugm:sum_ib q * (tau+ppi-tau) : corresponds to valftr in mkpot.f
    !l   smresh:sm rho * sm V ; corresponds to valfsm in mkpot.f
    !r Remarks
    !r   The valence kinetic energy is evaluated in the usual way as
    !r        ekinval = sumev - srhov
    !r   where sumev is the sum-of-eigenvalues and srhov is the integral
    !r   of the output density and input potential.
    !r   Integrals of the density with the xc potential are folded into the
    !r   electrostatic parts, that is:
    !r     V0 = V0es + V0xc  V1 = V1es + V1xc  V2 = V2es + V2xc
    !r   and are not discussed here.
    !r
    !r   mkekin make the electrostatic integral
    !r     int (n0~ Ves~ + n1 Ves1 - n2 Ves2~)                         (40)
    !r   as described in
    !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
    !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
    !r      ed. (Springer-Verlag, Berlin) 2000.
    !r   for an output density defined through (smrho,qkkl) and an input
    !r   potential defined through smpot and the matrix elements (ppi-tau).

    !r   Consider only one atom/cell for simplicity. The output density is
    !r     n = sum_ij Dij F~i F~j with Dij = {sum_n w_n z*_in z_jn}   (32)
    !r   i.e. the contraction of the density matrix with partial densities
    !r     F~i F~j = Fi Fj +
    !r               sum_kLk'L' C^i_kL (P~kLP~k'L' - PkLPk'L') C^j_k'L' (33)
    !r             = n0_ij + n1_ij - n2_ij
    !r   Note that local parts of the partial densities have two `levels'
    !r   of decomposition, namely at the `ij' level as in Eq. 33, or
    !r   at a still finer level in which the (kLk'L') indices are not
    !r   summed over.  Thus
    !r     n{1,2} = sum_ij D_ij n{1,2}_ij
    !r     n{1,2}_ij = sum_kLk'L' C^i_kL n{1,2}_kL,k'L' C^j_k'L'
    !r     n{1,2}_kL,k'L' = PkL Pk'L'
    !r   Note also that the 'k index' may be a sum over polynomials, or when
    !r   function `heads' are dealt with, the function itself, as described
    !r   in augmat.f.  As in making the matrix elements, we have to deal
    !r   with three cases, HH; HP; PP, but this is a inessential detail
    !r   needed only because representing H with sums of polynomials tends
    !r   to be ill-conditioned, and in the description below we ignore it.
    !r
    !r   Densities n0 and n2 have corresponding n0~ and n2~ which include
    !r   the additional multipole terms that guarantee n1 and n2~ have
    !r   the same multipole moments.  Thus:
    !r     n0~ = n0 + sum_M q_M G_M
    !r     n2~ = n2 + sum_M q_M G_M
    !r   where q_M are the difference in multipole moments between n1 and n2
    !r     q_M = int dr Y_M r^m (n1 - n2)
    !r   We can define partial densities for multipole contributions as well
    !r     n2~-n2 = sum_ij D_ij (n2~-n2)_ij
    !r     (n2~-n2)_ij = sum_M Q_ijM G_M
    !r                 = sum_kLk'L'M C^i_kL Q_kkLL'M G_M C^j_k'L'
    !r   with the two forms decomposing q_M into two levels:
    !r     q_M = sum_ij D_ij Q_ijM
    !r     Q_ijM = sum_kLk'L' C^i_kL Q_kkLL'M C^j_k'L'
    !r     Q_kkLL'M = int dr Y_M r^m (P~kL P~k'L' - PkL Pk'L')         (27)
    !r
    !r   Using the identity
    !r     n2~ - n2 = n0~ - n0 = sum_M q_M G_M
    !r   Eq. 40 is evaluated as
    !r     int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + sum_M q_M G_M (Ves0~-Ves2~)
    !r   The first term is evaluated on the mesh and stored in srmesh
    !r   The remaining terms amount to products of the density-matrix
    !r   and the ppi matrix elements.  Thus:
    !r     int n1 Ves1 = sum_ij D_ij int n1_ij Ves1
    !r     int n1_ij Ves1 = sum_kLk'L' C^i_kL int n1_kL,k'L' Ves1 C^j_k'L'
    !r                    = sum_kLk'L' C^i_kL pi1_kk'LL' C^j_k'L'
    !r   where pi1 is the first term of the pi matrix element, Eq. 29:
    !r     pi1_kk'LL' = P~kL V1 P~k'L'
    !r   Similarly for the second term, substituting n2 for n1 and
    !r   Ves2~ for Ves1.
    !r     int n2 Ves2~ = sum_ij D_ij int n2_ij Ves2~
    !r     int n2_ij Ves2~ = sum_kLk'L' C^i_kL int n2_kL,k'L' Ves2~ C^j_k'L'
    !r                    = sum_kLk'L' C^i_kL pi2_kk'LL' C^j_k'L'
    !r     pi2_kk'LL' = P~kL V1 P~k'L'
    !r   The last term just amounts to products of the density-matrix and
    !r   the remaining parts of the ppi matrix element:
    !r     pi_kk'LL'  = pi1_kk'LL' - pi2_kk'LL' + pi3_kk'LL'
    !r     pi3_kk'LL' = sum_M Q_kkLL'M int G_M (Ves0~ - Ves2~)
    !r   Evaluating the last term in the electrostatic integral we have
    !r     rhoV_MP = int sum_M q_M G_M (Ves0~ - Ves2~)
    !r             = int sum_ij D_ij sum_M Q_ijM G_M (Ves0~ - Ves2~)
    !r             = sum_ij D_ij sum_kLk'L'M C^i_kL pi3_kk'LL' C^j_k'L'
    !r   which follows using the relationship between Q_kkLL'M and Q_ijM
    !r   Using the definition of the local density-matrix (see rlocbl.f)
    !r      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
    !r   the electrostatic integral then becomes
    !r     int rhoVes = int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + rhoV_MP
    !r                = int n0 Ves0~
    !r                + sum_ij D_ij sum_kLk'L' C^i_kL pi_kk'LL' C^j_k'L'
    !r                = int n0 Ves0~ + sum_kLk'L' qpp'LL' pi_kk'LL'
    !r
    !u Updates
    !u   29 Jun 05 (MvS) SO hamiltonian not included in d.c. terms when
    !u             evaluating kinetic energy.
    !u   27 Aug 01 Extended to local orbitals.
    !u   18 Jun 00 spin polarized
    !u   20 Jun 00 adapted from nfp get_ekin
    ! ----------------------------------------------------------------------
    implicit none
    type(s_cv1),target :: sv_p_oppi(3,nbas)
    type(s_rv1),target :: sv_p_otau(3,1)
    type(s_rv1),target :: sv_p_osig(3,1)
    type(s_rv1),target :: sv_p_oqkkl(3,nbas)
    real(8),pointer:: OQPP(:), OQHP(:),OQHH(:), &
         OSIGHH(:),OSIGHP(:),OSIGPP(:), &
         OTAUHH(:),OTAUHP(:),OTAUPP(:)
    complex(8),pointer:: OPPIHH(:),OPPIHP(:),OPPIPP(:)
    real(8):: sumev , ekinval , vconst
    complex(8):: smpot(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp)
    integer :: ib,igetss,ipr,is,kmax,lgunit,lmxa,lmxh,n0,n1,n2,n3, &
         ngabc(3),nglob,nkap0,nlma,nlmh
    logical :: lgors
    parameter (n0=10,nkap0=3)
    double precision :: qum1,qum2,sraugm,srhov,srmesh,sum1,sum2,sumh,sumq,sumt,vol,xx
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    include 'mpif.h'
    integer:: procid=0,ier=0,io,iorb
    integer,parameter::master=0
    logical:: iprx
    call mpi_comm_rank(mpi_comm_world,procid,ier)
    iprx=.false.
    if(procid==master) iprx= .TRUE. 
    call tcn('mkekin')
    call getpr(ipr)
    ngabc=lat_nabc
    vol=lat_vol
    ! --- Integral n0(out) (Ves0~ + Vxc0), contribution from mesh ---
    !     Note that it does not include the term (n0~-n0) Ves0~
    srmesh = (dreal(sum(smpot*smrho)) + vconst*dreal(sum(smrho))) *vol/(n1*n2*n3)
    ! --- Integral rhout*veff, part from augmentation ---
    sraugm = 0d0
    do  ib = 1, nbas
       is = ispec(ib) 
       lmxa=sspec(is)%lmxa
       kmax=sspec(is)%kmxt
       lmxh=sspec(is)%lmxb
       if (lmxa == -1) goto 10
       call orblib(ib) ! norb , ltab , ktab , offl ,ntab,blks
       nlma = (lmxa+1)**2
       nlmh = (lmxh+1)**2
       OQHH => sv_p_oqkkl(3,ib)%v !head x head index=3
       OQHP => sv_p_oqkkl(2,ib)%v !head x tail index=2
       OQPP => sv_p_oqkkl(1,ib)%v !tail x tail index=1
       OTAUHH =>sv_p_otau(3,ib)%v
       OTAUHP =>sv_p_otau(2,ib)%v
       OTAUPP =>sv_p_otau(1,ib)%v
       OSIGHH =>sv_p_osig(3,ib)%v
       OSIGHP =>sv_p_osig(2,ib)%v
       OSIGPP =>sv_p_osig(1,ib)%v
       OPPIHH =>sv_p_oppi(3,ib)%cv
       OPPIHP =>sv_p_oppi(2,ib)%cv
       OPPIPP =>sv_p_oppi(1,ib)%cv
       call pvgtkn ( kmax , lmxa , nlma , nkaph , norb , ltab , ktab &
            , blks , lmxh , nlmh , OTAUHH , OSIGHH , OPPIHH &
            , OTAUHP , OSIGHP , OPPIHP &
            , OTAUPP , OSIGPP , OPPIPP  &
            , OQHH , OQHP , OQPP , nsp , nspc , sumt &
            , sumq , sumh )
       !       Add site augmentation contribution to rhout * (ham - ke)
       sraugm = sraugm + sumh - sumt
10     continue
    enddo
    srhov = srmesh + sraugm != n_out*Vin
    ekinval = sumev - srhov   != Eband - nout*Vin (V do not include SO term)
    if (ipr >= 30) write(stdo,"(/a)")' mkekin:'
    if (ipr >= 30) write(stdo,340) srmesh,sraugm,srhov,sumev,ekinval
340 format('   nout*Vin = smpart,onsite,total=:',3f14.6,&
         /'    E_B(band energy sum)=',f12.6,'  E_B-nout*Vin=',f12.6)
    call tcx('mkekin')
  end subroutine mkekin
  subroutine pvgtkn(kmax,lmxa,nlma,nkaph,norb,ltab,ktab,blks,lmxh, &
       nlmh,tauhh,sighh,ppihhz,tauhp,sighp,ppihpz, &
       taupp,sigpp,ppippz,qhh,qhp,qpp,nsp,nspc,&
       sumt,sumq,sumh) 
    !- Local contribution to kinetic energy for one site
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   kmax  :cutoff in PkL expansion
    !i   lmxa  :dimensions sigpp, taupp
    !i   nlma  :L cutoff in PkL expansion
    !i   nkaph :dimensions augmentation matrices
    !i   norb  :number of orbitals for this site
    !i   ltab  :table of l quantum numbers for the orbitals
    !i   ktab  :table of k numbers (orbital type) for the orbitals
    !i   blks  :block size for grouping orbitals into blocks (gtbls1)
    !i   lmxh  :dimensions sighh, sighp, tauhh, tauhp
    !i   nlmh  :dimensions heads ppi and qhh and qhp
    !i   tauhh :head-head kinetic energy integrals (augmat.f)
    !i   sighh :head-head overlap integrals (augmat.f)
    !i   ppihh :head-head kinetic + potential integrals (augmat.f)
    !i   tauhp :head-tail kinetic energy integrals (augmat.f)
    !i   sighp :head-tail overlap integrals (augmat.f)
    !i   ppihp :head-tail kinetic + potential integrals (augmat.f)
    !i   taupp :tail-tail kinetic energy integrals (augmat.f)
    !i   sigpp :tail-tail overlap integrals (augmat.f)
    !i   ppipp :tail-tail potential integrals (augmat.f)
    !i   lcplxp=1 only:0 if ppi is real; 1 if ppi is complex
    !i   qhh   :head-head density matrix for this site
    !i   qhp   :head-tail density matrix for this site
    !i   qpp   :tail-tail density matrix for this site
    !i   nsp   :number of spin channels
    !i   nspc  :2 for coupled spins; otherwise 1
    !o Outputs
    !o   sumt  :site contribution to kinetic energy
    !o   sumq  :site contribution to overlap (charge ?)
    !o   sumh  :site contribution to kinetic energy + potential
    !r Remarks
    !u Updates
    !u    1 Sep 04 Adapted to handle complex ppi
    !u   28 Aug 01 Extended to local orbitals.
    ! ----------------------------------------------------------------------
    implicit none
    integer :: kmax,lmxa,nlma,lmxh,nlmh,nsp,nspc
    integer :: nkaph,norb,ltab(norb),ktab(norb),blks(norb)
    real(8):: &
         tauhh(nkaph,nkaph,0:lmxh,nsp),sighh(nkaph,nkaph,0:lmxh,nsp), &
         tauhp(nkaph,0:kmax,0:lmxh,nsp),sighp(nkaph,0:kmax,0:lmxh,nsp), &
         taupp(0:kmax,0:kmax,0:lmxa,nsp),sigpp(0:kmax,0:kmax,0:lmxa,nsp), &
         qhh(nkaph,nkaph,nlmh,nlmh,nsp), &
         qhp(nkaph,0:kmax,nlmh,nlma,nsp), &
         qpp(0:kmax,0:kmax,nlma,nlma,nsp)
    complex(8)::&
         ppihhz(nkaph,nkaph,nlmh,nlmh,nsp), &
         ppihpz(nkaph,0:kmax,nlmh,nlma,nsp), &
         ppippz(0:kmax,0:kmax,nlma,nlma,nsp)
    integer :: ilm1,ilm2,k1,k2,ll,nlm11,nlm12,nlm21,nlm22,i
    double precision :: sumt,sumq,sumh,xx
    integer :: io1,io2,l1,l2
    ! ... Pkl*Pkl
    sumt = sum([(sum(qpp(:,:,ilm1,ilm1,:)*taupp(:,:,ll(ilm1),:)),ilm1=1,nlma)])
    sumq = sum([(sum(qpp(:,:,ilm1,ilm1,:)*sigpp(:,:,ll(ilm1),:)),ilm1=1,nlma)])
    sumh = sum(qpp(:,:,:,:,:)*ppippz(:,:,:,:,:))
    ! ... Hsm*Hsm
    do  io2 = 1, norb;    if(blks(io2)==0) cycle 
       k2 = ktab(io2)
       nlm21 = ltab(io2)**2+1
       nlm22 = nlm21 + blks(io2)-1
       do  io1 = 1, norb; if(blks(io1)==0) cycle 
          associate(k1=>ktab(io1),nlm11 => ltab(io1)**2+1)
            sumh = sumh+sum([( &
                 sum(qhh(k1,k2,ilm1,nlm21:nlm22,:)*ppihhz(k1,k2,ilm1,nlm21:nlm22,:)) &
                 ,ilm1=nlm11,nlm11+blks(io1)-1) ])
            do  ilm1 = nlm11, nlm11+ blks(io1)-1
               if( nlm21<= ilm1 .and. ilm1<=nlm21+blks(io2)-1) then
                  sumt = sumt + sum(qhh(k1,k2,ilm1,ilm1,:)*tauhh(k1,k2,ll(ilm1),:))
                  sumq = sumq + sum(qhh(k1,k2,ilm1,ilm1,:)*sighh(k1,k2,ll(ilm1),:))
               endif
            enddo
          endassociate
       enddo
    enddo
    ! ... Hsm*Pkl
    do  io1 = 1, norb; if (blks(io1)==0) cycle
       associate(k1=>ktab(io1), nlm11=>ltab(io1)**2+1, nlm11e=>min(ltab(io1)**2+1 +blks(io1)-1,nlma) )
         sumh=sumh+sum([(sum(qhp(k1,:,ilm1,:,:)*ppihpz(k1,:,ilm1,:,:)),ilm1= nlm11,nlm11+blks(io1)-1)])
         sumt=sumt+sum([(sum(qhp(k1,:,ilm1,ilm1,:)*tauhp(k1,:,ll(ilm1),:)),ilm1= nlm11,nlm11e)])
         sumq=sumq+sum([(sum(qhp(k1,:,ilm1,ilm1,:)*sighp(k1,:,ll(ilm1),:)),ilm1= nlm11,nlm11e)])
       endassociate
    enddo
  end subroutine pvgtkn
  subroutine makdos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax, ndos,dos)! Make density of states from bands
    !-----------------------------------------------------------------------
    !i  Input
    !i   nqp   :number of q-points
    !i   nband :number of bands
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   wgts  :band weights
    !i   evl   :band eigenvalues
    !i   n,w   :Methfessel-Paxton order and broadening parameters
    !i   tol   :(tol>0) allowed error in DOS due to truncating the gaussian,
    !i         :        to a finite energy range (number of bins)
    !i         :(tol<0) dimensionless energy window specifying truncation
    !i         :        of gaussian.  Energy window for which gaussian is
    !i         :        taken to be nonzero is set to -tol*w
    !i   emin, emax, ndos: energy range and number of energy mesh points
    !i   nbmx  :leading dimension of evl
    !o  Ouput
    !o    dos: density of states
    !-----------------------------------------------------------------------
    implicit none
    integer :: nqp,nband,nbmx,nsp,n,ndos
    double precision :: wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp),w,emin,emax,tol,wt,emesh
    integer :: i,isp,iband,iq,meshpt,mesh1,mesh2,mrange=999999,iprint
    double precision :: e,x,range,test,step,d,s,xx
    dos=0d0
    step = (emax - emin) / (ndos - 1)
    if ( tol > 0d0 ) then
       do   i = 0, ndos-1
          x = i * step / w
          call delstp(0,x,test,s,xx)
          if ( test < tol ) then
             mrange = i + 1
             goto 3
          endif
       enddo
       if (iprint() > 30) print *,'makdos (warning) : tol too small'
3      continue
       range = 2 * mrange * step
       test = tol
    else
       range = -tol * w
       mrange = range / ( 2 * step )
       call delstp(0,-tol/2,test,s,xx)
    endif
    if (iprint() > 30) write (*,100) range/w,2*mrange,test
    iqloop: do  7  iq = 1, nqp
       wt = abs(wgts(iq)) / nsp
       ibandloop: do  61  iband = 1, nband
          isploop: do   isp = 1, nsp
             e = evl(iband,isp,iq)
             meshpt = (e - emin) / step
             mesh1 = meshpt - mrange
             mesh2 = meshpt + mrange
             if (mesh2 >= ndos) mesh2 = ndos-1
             if (mesh1 < 0) mesh1 = 0
             do   meshpt = mesh1, mesh2
                emesh = emin + meshpt * step
                x = (emesh - e) / w
                call delstp(n,x,d,s,xx)
                dos(meshpt,isp) = dos(meshpt,isp) + wt * d / w
             enddo
          enddo isploop
61     enddo ibandloop
7   enddo iqloop
100 format(/1x,'MAKDOS :  range of gaussians is ',f5.2,'W (',i4,' bins).' &
         /11x,'Error estimate in DOS : ',1pe9.2,' per state.')
  end subroutine makdos
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

