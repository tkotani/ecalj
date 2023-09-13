!= Main part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration = 2023-jan memo checked
! From density osmrho,osmrho, m_mkpot_init gives potential,
! osmpot, and osig, otau, oppi (potential site integrals)
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
  use m_density,only: orhoat,osmrho,eferm  !input/output unprotected  ! NOTE:variable in m_density are not protected
  real(8),protected,public:: ham_ehf, ham_ehk, sev  !output
  real(8),protected,public:: qdiff                  !output
  real(8),protected,allocatable,public:: force(:,:) !output
  !NOTE: other ouputs are in: m_mkpot(potential), m_bandcal(band,dmatu), and m_density(density)
  public bndfp
  private
contains
  subroutine bndfp(iter, llmfgw, plbnd) !Single iteration for given density and dmatu
    ! llmfgw=T is for generating eigenfunctions for GW calculations, no iteration.
    ! plbnd/=0 means band plot mode. no iteration.
    !     ! All read only in bndfp. Data are stored in modules such as m_bandcal, m_mkpot
    !     ! For example,, rightafter call m_bandcal_init, we can get evalall, which is used in other modules.
    use m_ftox
    use m_mixrho,only: mixrho
    use m_bndfp_util,only: mkekin,makdos,phispinsym_ssite_set,iorbtm
    use m_supot,only: n1,n2,n3 !for charge mesh
    use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx,nspx=>ham_nspx
    use m_lmfinit, only: ncutovl,lso,ndos=>bz_ndos,bz_w,fsmom=>bz_fsmom, bz_dosmax,lmet=>bz_lmet,bz_fsmommethod,bz_n, &
         ldos,qbg=>zbak,lfrce,pwmode=>ham_pwmode,lrsig=>ham_lsig,epsovl=>ham_oveps, &
         ham_scaledsigma, alat=>lat_alat,stdo,stdl,procid,master, nkaph,nlmax,nl,nbas,nsp, bz_dosmax, &
         lmaxu,nlibu,lldau,lpztail,leks,lrout,  nchan=>pot_nlma, nvl=>pot_nlml,nspc,pnufix !lmfinit contains fixed input 
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
    use m_rdsigm2,only: m_rdsigm2_init
    use m_subzi, only: m_subzi_init, lwtkb,rv_a_owtkb,m_subzi_setlwtkb,m_subzi_bzintegration
    use m_rsibl, only : Rsibl_ev ! to plot wavefunction in the fat band mode
    use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq,xmpbnd2
    use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate, m_mkpot_energyterms,m_mkpot_novxc,& ! & m_mkpot_novxc_dipole,
         osmpot, qmom, vconst, osig,otau,oppi, qval , qsc , fes1_rv , fes2_rv
    use m_clsmode,only: m_clsmode_init,m_clsmode_set1,m_clsmode_finalize
    use m_qplist,only:  qplist,nkp,xdatt,labeli,labele,dqsyml,etolc,etolv, &
         nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml,iqini,iqend,ispini,ispend,kpproc    ! MPIK divider. iqini:iqend are node-dependent
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x
    use m_procar,only: m_procar_init,dwgtall,nchanp,m_procar_closeprocar,m_procar_writepdos
    use m_bandcal,only: m_bandcal_init,m_bandcal_2nd,m_bandcal_clean,m_bandcal_allreduce, &
         evlall,smrho_out,oqkkl,oeqkkl, ndimhx_,nevls,m_bandcal_symsmrho !,nev_
    use m_mkrout,only: m_mkrout_init,orhoat_out,frcbandsym,hbyl_rv,qbyl_rv
    use m_addrbl,only: m_addrbl_allocate_swtk,swtk
    use m_sugw,only: m_sugw_init
    use m_mkehkf,only: m_mkehkf_etot1,m_mkehkf_etot2
    use m_gennlat_sig,only: m_gennlat_init_sig
    use m_dfrce,only: dfrce
    use m_sugcut,only:sugcut
    ! inputs
    ! main input are osmrho, orhoat (in m_mkpot_init), specifing density, and vorb, which is potential for LDA+U
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
    !l   n1,n2,n3: dimensions smrho,smpot.
    !!      nspx: number of independent spin channels
    !!      nspc is now avoided (memo:nspc=2 for lso==1, nspc=1 for lso/=1 See m_lmfinit)
    !!      ndhamx: maximum size of hamiltonian
    !r How bndfp works?
    !r   (1) m_mkpot_init   make the effective potential,
    !r   (2) m_bandcal_init generate eigenvalues (and eigenvectors if lrout), then m_bzintegration_init
    !r   (3) m_bandcal_2nd  If lrout, assemble the output density by BZ integration
    !r   (4) evaluate hf (and KS, if leks) energy by BZ integration
    !r   (5) mixrho the output density to make a new input density.
    !See history github ecalj after 2009.
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
    integer:: idummy, unlink,ifih,ifii,ib,ix,ifimag,nevmin
    real(8):: ekinval,eharris,eksham,  dosw(2),dum,evtop,ecbot,qp(3),rydberg,xxx,eeem
    real(8):: fh_rv(3,nbas),vnow_magfield,eee,dosi(2),dee,efermxxx,emin,qvalm(2)
    real(8),allocatable:: dosi_rv(:,:),dos_rv(:,:),qmom_in(:)
    real(8),parameter::  NULLR =-99999,pi=4d0*datan(1d0)
    real(8):: elind=0d0
    call tcn ('bndfp')
    debug  = cmdopt0('--debugbndfp')
    tdos   = cmdopt0('--tdos')  !total dos mode or not
    fsmode = cmdopt0('--fermisurface')!FermiSurfece for xcrysden in http://www.xcrysden.org/doc/XSF.html#2l.16
    ipr    = iprint() ! for procid/=master, we set iprint=0 at lmv7.F
    ltet = ntet>0! tetrahedron method or not
    GETefermforplbnd: if(plbnd/=0) then
       open(newunit=ifi,file='efermi.lmf')
       read(ifi,*) eferm  ! efermi.lmf can be different from efermi in rst.* file.
       close(ifi)         !  For example, you may change NKABC, and run job_pdos (lmf-MPIK --quit=band only modify efermi.lmf, without touching rst.*).
       call mpibc1_real(eferm,1,'bndfp_eferm')
    endif GETefermforplbnd
    if(cmdopt0('--phispinsym')) call phispinsym_ssite_set() !pnu,pz are spin symmetrized ! Set spin-symmetrized pnu. aug2019. See also in pnunew and locpot
    writeham= cmdopt0('--writeham') ! Write out Hamiltonian HamiltonianPMT.*
    if(writeham) open(newunit=ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
    if(llmfgw) call m_mkpot_novxc() !Get osigx,otaux oppix spotx (without XC(LDA))  one-particle potential without XC part for GWdriver: lmfgw mode
    !      if(llmfgw.and.cmdopt0('--dipolematrix')) call m_mkpot_novxc_dipole()
    call m_mkpot_init() !from (smrho,rhoat), Get one-particle potential and energy-related quantities. mkpot->locpot->augmat. augmat calculates sig,tau,ppi.
    if(cmdopt0('--quit=mkpot')) call rx0('--quit=mkpot')
    call m_subzi_init(lrout>0) ! Setup for wtkb for BZ integration.  ! NOTE: if (wkp.* exists).and.lmet==2, wkp is used for wkkb.
    if(lpztail) call sugcut(2) !   lpztail: if T, local orbital of 2nd type(hankel tail).
    !    Hankel's e of local orbital of PZ>10 (hankel tail mode) is changing. ==>T.K think current version of PZ>10 might not give so useful advantages.
    if(cmdopt0('--cls')) then !clsmode
       call rxx(lso==1,'CLS not implemented in noncoll case')
       if (lrout == 0) call rx('bndfp: need output density for cls')
       call m_clsmode_init()
    endif
    sigmamode = (lrsig/=0)
    readingsigm: if( sigmamode .AND. siginit ) then !sigm contains \Sigma-Vxc 
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
    endif readingsigm
    !  ndimsig is the dimension of the self-energy. Current version is "ndimsig = ldim", which means we use only projection
    !  onto MTO spaces even when PMT. 
    if(sigmamode .AND. master_mpi) write(stdo,*)' ScaledSigma=',ham_scaledsigma
    ! ==
    GWdriver: if(llmfgw) then        !         call m_sugw_init(cmdopt0('--dipolematrix'),cmdopt0('--socmatrix'),eferm)
       call m_sugw_init(cmdopt0('--socmatrix'),eferm)
       call tcx('bndfp')
       call rx0('sugw mode')  !exit program here normally.
    endif GWdriver
    ! Set up Hamiltonian and diagonalization in m_bandcal_init. To know outputs, see 'use m_bandcal,only:'. The outputs are evlall, and so on.
    if( .NOT. (lwtkb==-1 .OR. lwtkb==0 .OR. lwtkb==1)) call rx('bndfp: something wrong lwtkb')
    sttime = MPI_WTIME()
    if(nspc==2) call m_addrbl_allocate_swtk(ndham,nsp,nkp)
    call m_bandcal_init(iqini,iqend,ispini,ispend,lrout,eferm,ifih,lwtkb) !All input. get evlall ! lwtkb=1 accumulate weight, lwtkb=0 no weight
    entime = MPI_WTIME()
    if(master_mpi) write(stdo,"(a,f9.4)") ' ... Done MPI k-loop: elapsed time=',entime-sttime
    sttime = MPI_WTIME()
    if(writeham) close(ifih)
    if(writeham) call rx0('Done --writeham: --fullmesh may be needed. HamiltonianMTO* genereted')
    call mpibc2_int(ndimhx_,size(ndimhx_),'bndfp_ndimhx_') !all reduce (instead of which node to which node).
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
          if(master_mpi .AND. iq==1)write(stdl,"('fp evl',8f8.4)")(evlall(i,jsp,iq),i=1,nevls(iq,jsp))
       enddo
    enddo
    fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') ! pdos mode (--mkprocar and --fullmesh)
    PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
    nevmin=minval(nevls(1:nkp,1:nspx))
    if(fullmesh .AND. procaron) call m_procar_writepdos(evlall,nevmin,eferm,kpproc) !nev_
    if(fullmesh .AND. procaron) call rx0('Done pdos: --mkprocar & --fullmesh. Check by "grep k-point PROCAR.*.*"')
    if( cmdopt0('--boltztrap')) call writeboltztrap(eferm) ! boltztrap data
    if( cmdopt0('--boltztrap')) call rx0('Done boltztrap: boltztrap.* are generated')
    Writebandmode: if(plbnd/=0 ) then
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
    endif Writebandmode
    if(plbnd/=0 .AND. fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
    if(plbnd/=0) call rx0('plot band mode done') ! end of band plbnd/=0, that is, band plot mode.
    call m_subzi_bzintegration(evlall,swtk, eferm,sev,qvalm,vnow_magfield) !Get the Fermi energy eferm,... New eferm and wtkb determined from evlall
    if(fsmom/=NULLR .AND. master_mpi) then !moved from m_subzi_bzintegration at 2022-dec
       open(newunit=ifimag,file='MagField')
       write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow_magfield 
       close(ifimag)
    endif ! -vnow_magfield/2 is added to isp=1, +vnow_magfield/2 to isp=2
    evtop=-9999
    ecbot=9999
    eeem=9999
    do iq=1,nkp
       do jsp=1,nspx     !nspx=1 for SO=1
          do ix=1,ndhamx !2022dec evtop,ecbot for magnetic case
             if(evlall(ix,jsp,iq)<eferm) evtop = max(evtop,evlall(ix,jsp,iq))
             if(evlall(ix,jsp,iq)>eferm) ecbot = min(ecbot,evlall(ix,jsp,iq))
             eeem = min(eeem,evlall(ix,jsp,iq)) !eeem for emin below
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
    GenerateTotalDOS: if(master_mpi .AND. (tdos .OR. ldos/=0)) then !   emin=dosw(1) emax=dosw(2) dos range
       emin = eeem-0.01d0 
       dosw(1)= emin ! lowest energy limit to plot dos
       dosw(2)= eferm+bz_dosmax !max energy limit to plot dos
       write(stdo,ftox)' bndfp:Generating TDOS: efermi=',ftof(eferm),' dos window emin emax= ',ftof(dosw)
       allocate( dosi_rv(ndos,nspx),dos_rv(ndos,nspx)) !for xxxdif
       if(cmdopt0('--tdostetf')) ltet= .FALSE. ! Set tetrahedron=F
       if(ltet) then
          call bzints(nkabc(1)*nkabc(2)*nkabc(3), evlall, dum ,nkp,nevmin,ndhamx,nspx , dosw(1),dosw(2), dosi_rv, ndos,xxx, & !ndhamx=>nevmin at 2023feb
               1, ntet , iv_a_oidtet , dum , dum ) !job=1 give IntegratedDos to dosi_rv
          dos_rv(2:ndos-1,:)=(dosi_rv(3:ndos,:)-dosi_rv(1:ndos-2,:))/(2d0*(dosw(2)-dosw(1))/(ndos-1))
          dos_rv(1,:)    = dos_rv(2,:)
          dos_rv(ndos,:) = dos_rv(ndos-1,:)
       else
          call makdos(nkp,nevmin,ndhamx,nspx,rv_a_owtkp,evlall,bz_n,bz_w,-6d0,dosw(1),dosw(2),ndos,dos_rv) !ndmahx=>nevmin
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
    endif GenerateTotalDOS
    if(tdos) call rx0('Done tdos mode:')
    if(lrout/=0) call m_bandcal_allreduce(lwtkb) ! AllReduce band quantities.
    emin=1d9
    if(master_mpi) write(stdo,ftox)
    do iq = 1,nkp
       qp=qplist(:,iq)
       do isp = 1,nspx
          jsp = isp
          if(master_mpi .AND. iprint()>=35) then
             write(stdo,ftox)" bndfp: kpt",iq," of",nkp," k isp=",ftof(qp,4),jsp," ndimh nev=",ndimhx_(iq,jsp),nevls(iq,jsp)
             write(stdo,ftox) ftof([(evlall(i,jsp,iq), i=1,nevls(iq,jsp))],4)
          endif
          emin= min(minval( evlall(1:nevls(iq,jsp),jsp,iq)),emin)
       enddo
    enddo
    AccumurateSuminBZforwtkb: if(lmet>0 .AND. lwtkb==-1 .AND. lrout>0) then
       call m_subzi_setlwtkb(1)
       call mpi_barrier(MPI_comm_world,ierr)
       call m_bandcal_2nd(iqini,iqend,ispini,ispend,lrout)!, eferm) !accumulate smrho_out and so on.
       call m_bandcal_allreduce(lwtkb)
    endif AccumurateSuminBZforwtkb
    CorelevelSpectroscopy: if(cmdopt0('--cls')) then !m_clsmode_set1 is called in m_bandcal
       dosw(1)= emin  - 0.5d0    ! lowest energy limit to plot dos
       dosw(2)= eferm + bz_dosmax ! highest energy limit to plot dos
       call m_clsmode_finalize(eferm,ndimh,ndhamx,nspx,nkp,dosw,evlall)
       call rx0('Done cls mode:')
    endif CorelevelSpectroscopy
    if(lwtkb==1 .AND. lso/=0) call iorbtm()  ! write orbital moment
    if (lrout/=0) call m_bandcal_symsmrho()  !Get smrho_out Symmetrize smooth density ! Assemble output density, energies and forces 
    call m_mkrout_init() !Get force frcbandsym, symmetrized atomic densities orhoat_out, and weights hbyl,qbyl
    call pnunew(eferm) !pnuall are revised. !  New boundary conditions pnu for phi and phidot
    !  call writeqpyl() !Set if you like to print writeqbyl
    call m_mkehkf_etot1(sev, eharris) ! Evaluate Harris-foukner energy (note: we now use NonMagneticCORE mode as default)
    EvaluateKStotalEnergyandForce:if(lrout/=0) then
       allocate(qmom_in(nvl))
       qmom_in=qmom !multipole moments.
       eksham = 0d0 !   ... Evaluate KS total energy and output magnetic moment
       if(leks>=1) then
          call mkekin(osig,otau,oppi,oqkkl,vconst,osmpot,smrho_out,sev,  ekinval)
          call m_mkpot_energyterms(smrho_out, orhoat_out) !qmom is revised for given orhoat_out
          if(cmdopt0('--density')) then
             call mpi_barrier(MPI_comm_world,ierr)
             call rx0('end of --density mode')
          endif
          call m_mkehkf_etot2(ekinval, eksham)
       endif
       if(lfrce> 0) then !Add together force terms 
          ! fes1_rv: contribution to HF forces from estat + xc potential.  This is for input  density !=3rd term in (B.5) in JPSJ.84.034705
          ! fes2_rv: contribution to KS forces from estat + xc potential.  This is for output density
          ! frcbandsym : 1st term in (B.5)  (puley? need check)
          ! fh_rv      : 2nd term in (B.5)  (need check)
          if(allocated(force)) deallocate(force)
          allocate(force(3,nbas))
          call dfrce (lfrce,orhoat,orhoat_out,qmom_in,osmrho,smrho_out,  fh_rv)
          call totfrc(leks, fes1_rv, fes2_rv, fh_rv, frcbandsym, force) ! force : total
       endif
       deallocate(qmom_in)
       ! Mix inputs(osmrho,orhoat) and outputs(osmrho_out,orhoat_out), resulting orhoat and osmrho.
       call mixrho(iter,qval-qbg,orhoat_out,orhoat,smrho_out,osmrho,qdiff)!mixrho keeps history in it.
    else
       eksham = 0d0
    endif EvaluateKStotalEnergyandForce
    ham_ehf= eharris !Harris total energy
    ham_ehk= eksham  !Hohenberg-Kohn-Sham total energy
    call m_mkpot_deallocate()
    call m_bandcal_clean()
    call tcx('bndfp')
  end subroutine bndfp
end module m_bndfp
