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
!o   evlall: eigenvalue in m_bandcal
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
!! NOTE: check main_lmf -> lmfp(iteration loop) -> bndfp
module m_bndfp
  use m_lgunit,only:stdo,stdl
  use m_density,only: orhoat,smrho,eferm            !input/output unprotected  ! NOTE:variable in m_density are not protected
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
    !     ! For example, rightafter call m_bandcal_init, we can get evlall, which is used in other modules.
    use m_ftox
    use m_mixrho,only: mixrho
    use m_bndfp_util,only: mkekin,makdos,phispinsym_ssite_set,iorbtm
    use m_supot,only: n1,n2,n3 !for charge mesh
    use m_igv2x,only: nbandmx
    use m_lmfinit,only: ncutovl,lso,ndos=>bz_ndos,bz_w,fsmom=>bz_fsmom, bz_dosmax,lmet=>bz_lmet,bz_fsmommethod,bz_n,nspx
    use m_lmfinit,only: ldos,qbg=>zbak,lfrce,pwmode=>ham_pwmode,lrsig=>ham_lsig,epsovl=>ham_oveps !try to avoid line continuation in fortran
    use m_lmfinit,only: ham_scaledsigma, alat=>lat_alat, nlmax,nbas,nsp, bz_dosmax,nlmxlx,afsym,ispec, phispinsym
    use m_lmfinit,only: lmaxu,nlibu,lldau,lpztail,leks,lrout,  nchan=>pot_nlma, nvl=>pot_nlml,pnufix !lmfinit contains fixed input 
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
    use m_rdsigm2,only: m_rdsigm2_init
    use m_subzi,only: m_subzi_init,m_subzi_bzintegration 
    use m_MPItk,only: master_mpi, strprocid, numprocs=>nsize,xmpbnd2,comm
    use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate, m_mkpot_energyterms,m_mkpot_novxc 
    use m_mkpot,only: osmpot, qmom, vconst, qval , qsc , fes1_rv , fes2_rv, amom
    use m_locpot,only: osig,otau,oppi,vesaverage
    use m_clsmode,only: m_clsmode_init,m_clsmode_set1,m_clsmode_finalize
    use m_qplist,only: qplist,nkp,xdatt,labeli,labele,dqsyml,etolc,etolv
    use m_qplist,only: nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml,kpproc,iqini,iqend    ! MPIK divider. iqini:iqend are node-dependent
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x
    use m_procar,only: dwgtall,nchanp,m_procar_closeprocar,m_procar_writepdos,m_procar_init
    use m_bandcal,only: m_bandcal_init,m_bandcal_2nd,m_bandcal_clean,m_bandcal_allreduce
    use m_bandcal,only: smrho_out,oqkkl,oeqkkl, ndimhx_,nevls,m_bandcal_symsmrho,evlall,spinweightsoc
    use m_mkrout,only: m_mkrout_init,orhoat_out,frcbandsym,hbyl_rv,qbyl_rv
    use m_sugw,only: m_sugw_init
    use m_mkehkf,only: m_mkehkf_etot1,m_mkehkf_etot2
    use m_gennlat_sig,only: m_gennlat_init_sig
    use m_dfrce,only: dfrce
    use m_sugcut,only:sugcut
    use m_bzints,only: bzints
    use m_writeband,only: writeband,writefs,writepdos,writedossawada
    use m_totfrc,only:totfrc
    ! inputs
    ! main input are smrho, orhoat (in m_mkpot_init), specifing density, and vorb, which is potential for LDA+U
    !i   nbas  : size of basis
    !i   nsp   : number of spins
    !i   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
    !i   lmaxu : lmax for U used to dimension vorb and dmatu. lmaxu=2 if d is, but lmaxu=3 if f is included.
    !i   lldau :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat in dmats(*,lldau(ib))
    !i   ndham :dimensioning parameter, at least as large as largest
    !i         :hamiltonian dimension
    !i   leks  :>0 make the Hohnberg-Kohn-Sham energy
    !i         :>1 use the HKS forces
    !i   lrout : =1 generate output density and attendant quantities, =0 do not
    !i   lfrce : =1 generate force (some other options)
    !i   lpnu=1  : =1 means new pnu's
    !i   iter  :current iteration number
    !i   maxit :maximum number of iterations
    !i  LDA+U inputs and outputs
    !o   dmatu :density matrix for LDA+U (changed upon output)
    !i   vorb  :orbital dependent LDA+U potential

    !o Outputs
    !o   frc   :forces.  Only calculated if lfrce>0.
    !o         :If leks=1, forces are HF  forces
    !o         :If leks=2, forces are HKS forces
    !l   n1,n2,n3: dimensions smrho,smpot.
    !!      nspx: number of independent spin channels (nspx=nsp/nspc)
    !       nspc: =1 lso=2, =1 otherwize
    !!      nbandmx: maximum size of hamiltonian: nbandmx can be twiced for lso=1
    !r How bndfp works?
    !r   (1) m_mkpot_init   make the effective potential,
    !r   (2) m_bandcal_init generate eigenvalues (and eigenvectors if lrout), then m_bzintegration_init
    !r   (3) m_bandcal_2nd  If lrout, assemble the output density by BZ integration
    !r   (4) evaluate hf (and KS, if leks) energy by BZ integration
    !r   (5) mixrho the output density to make a new input density.
    
    !See history github ecalj after 2009.
    implicit none
    include "mpif.h"
    character(10) :: head
    integer:: plbnd,nk1,nk2,nk3,nx,ny, iter,i,ifi,ipr,iq,isp,jsp,iprint,ipts,ierr
    integer:: ifih,ifii,ib,ix,ifimag,nevmin,nnn,ikp
    logical:: llmfgw,sigx, ltet,cmdopt0,sigmamode,tdos,debug=.false.
    logical:: fullmesh,PROCARon,writeham=.false.,magexist, fsmode 
    logical,save:: siginit=.true.
    real(8):: sttime,entime,vesav
    real(8):: ekinval,eharris,eksham,  dosw(2),dum(1),evtop,ecbot,qp(3),xxx,dumx
    real(8):: fh_rv(3,nbas),vmag,eee,dosi(2),dee,emin,qvalm(2)
    real(8),allocatable:: dosi_rv(:,:),dos_rv(:,:) 
    real(8),parameter::  pi=4d0*datan(1d0)
    real(8),external:: rydberg
    call tcn ('bndfp')
    debug  = cmdopt0('--debugbndfp')
    tdos   = cmdopt0('--tdos')  !total dos mode or not
    fsmode = cmdopt0('--fermisurface')!FermiSurfece for xcrysden in http://www.xcrysden.org/doc/XSF.html#2l.16
    ipr    = iprint() ! for procid/=master, we set iprint=0 at lmv7.F
    ltet = ntet>0! tetrahedron method or not
    vmag=0d0
    GETefermFORplbndMODE: if(plbnd/=0) then
      open(newunit=ifi,file='efermi.lmf',status='old',err=113)
      read(ifi,*) eferm,vmag ! efermi is consistent with eferm in rst.* file (iors.f90).
      read(ifi,*) 
      read(ifi,*) 
      read(ifi,*) xxx,vesav
      !                        However, efermi.lmf can be modified when we do one-shot dense-mesh calculation as is done in job_band.
      close(ifi)             ! For example, you may change NKABC, and run job_pdos (lmf --quit=band only modify efermi.lmf, without touching rst.*).
      goto 114
113   continue
      call rx('No efermi.lmf: need to repeat sc mode of lmf. --quit=band stops without changing rst file')
114   continue
    endif GETefermFORplbndMODE
    if(phispinsym) call phispinsym_ssite_set() !pnu,pz are spin symmetrized! Set spin-symmetrized pnu.aug2019. See also in pnunew and locpot
    writeham= cmdopt0('--writeham') ! Write out Hamiltonian HamiltonianPMT.*
    if(writeham) open(newunit=ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
    GetPotentialFromDensity: block
      if(llmfgw) call m_mkpot_novxc(smrho,orhoat) !Get osigx,otaux oppix spotx, which are onsite integrals without XC part for GWdriver: lmfgw mode
      call m_mkpot_init(smrho,orhoat)! From smrho and rhoat, get one-particle potential and related quantities. mkpot->locpot->augmat. augmat calculates sig,tau,ppi.
      if(cmdopt0('--quit=mkpot')) call rx0('--quit=mkpot')
    endblock GetPotentialFromDensity
    call m_subzi_init() ! Setup weight wtkb for BZ integration.  ! NOTE: if (wkp.* exists).and.lmet==2, wkp is used for wkkb.
    if(lpztail) call sugcut(2) ! lpztail: if T, local orbital of 2nd type(hankel tail). Hankel's e of local orbital of PZ>10 (hankel tail mode) is changing. ==>T.K think current version of PZ>10 might not give so useful advantages.
    CorelevelSpectroscopyINIT: if(cmdopt0('--cls')) then
      call rxx(lso==1,  'CLS not implemented in noncoll case')
      call rxx(lrout==0,'bndfp: need output density for cls')
      call m_clsmode_init()
    endif CorelevelSpectroscopyINIT
    sigmamode = (lrsig/=0)
    READsigmForQSGW: if( (sigmamode .AND. siginit).and.(.not.cmdopt0('--quitecore')) ) then !sigm contains \Sigma-Vxc 
      inquire(file='sigm.'//trim(sname),exist=sigx)
      if(sigx) then
        open(newunit=ifi,file='sigm.'//trim(sname),form='unformatted')
        read(ifi) nx,ny,nk1,nk2,nk3
        close(ifi)
        call m_gennlat_init_sig([nk1,nk2,nk3]) !nlat(real-space lattice vectors) for FFT of sigm
      endif
      call m_rdsigm2_init()
      call mpi_barrier(comm,ierr)
      siginit=.false. !only once
    endif READsigmForQSGW ! We now set "ndimsig=ldim",which means we use only projection onto MTO spaces even when PMT.
    call rxx(cmdopt0('--wsig_fbz'),'No sigm file from which we do --wsig_fbz')
    if(sigmamode .AND. master_mpi) write(stdo,*)' ScaledSigma=',ham_scaledsigma
    magexist= abs(vmag)>1d-6
    GWdriverExit: if(llmfgw) then
      call m_sugw_init(cmdopt0('--socmatrix'),eferm,vmag,qval)
      call tcx('bndfp')
      call rx0('sugw mode')  !exit program here normally.
    endif GWdriverExit
    ! Set up Hamiltonian and diagonalization in m_bandcal_init. To know outputs, see 'use m_bandcal,only:'. The outputs are evlall, and so on.
    sttime = MPI_WTIME() ! if(nspc==2) call m_addrbl_allocate_swtk(ndham,nsp,nkp)
    if(cmdopt0('--mkprocar')) call m_procar_init()
    GetHamiltonianAndDiagonalize: block
      call m_bandcal_init(lrout,eferm,vmag,ifih) ! Get Hamiltonian and diagonalization resulting evl,evec,evlall.
      entime = MPI_WTIME()                
      if(master_mpi) write(stdo,"(a,f9.4)") ' ... Done MPI k-loop: elapsed time=',entime-sttime
      if(writeham) close(ifih)
      if(writeham) call rx0('Done --writeham: --fullmesh may be needed. HamiltonianMTO* genereted')
      call mpibc2_int(ndimhx_,size(ndimhx_),'bndfp_ndimhx_') 
      call mpibc2_int(nevls,  size(nevls),  'bndfp_nevls')   
      nevmin = minval(nevls(1:nkp,1:nspx))
    endblock GetHamiltonianAndDiagonalize
    BROADCASTevlall:block
      integer:: nspxa
!      if(afsym) then
!        call xmpbnd2(kpproc,nbandmx,nkp,evlall(:,1,:)) !all eigenvalues are distributed nbandmx blocks
!        call xmpbnd2(kpproc,nbandmx,nkp,evlall(:,2,:)) 
!     else
      nspxa=merge(2,nspx,afsym)
      call xmpbnd2(kpproc,nbandmx,nkp*nspxa,evlall)   !all eigenvalues broadcasted !note (iq,isp) order in m_qplist.f90
     !endif
      if(lso==1) call xmpbnd2(kpproc,nbandmx*2,nkp,spinweightsoc)   !all eigenvalues broadcasted
      if(master_mpi) then
        do iq=1,1;do jsp=1,nspx; write(stdl,"('fp evl',8f8.4)")(evlall(i,jsp,iq),i=1,nevls(iq,jsp))
        enddo;        enddo
      endif
    endblock BROADCASTevlall
    BandPLOTmode: block
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') ! pdos mode (--mkprocar and --fullmesh)
      PROCARon = cmdopt0('--mkprocar') 
      if(plbnd/=0.or.(procaron.and.fullmesh).or.cmdopt0('--boltztrap')) then
         ! evlallm is removed. vmag is included in Hamitonian now.  2024-06-14
         evtop=maxval(evlall,mask=evlall<eferm)
         ecbot=minval(evlall,mask=evlall>eferm)
         Procarmode:block
           nevmin = minval(nevls(1:nkp,1:nspx))
           if(fullmesh .AND. procaron) call m_procar_writepdos(evlall,nevmin,eferm,kpproc) 
           if(fullmesh .AND. procaron) call rx0('Done pdos: --mkprocar .and. --fullmesh. Check by "grep k-point PROCAR.*.*"')
         endblock Procarmode
         Boltztrap:if( cmdopt0('--boltztrap')) then
            call writeboltztrap(eferm) ! boltztrap data
            call rx0('Done boltztrap: boltztrap.* are generated')
         endif Boltztrap
         Writebandmode: if(plbnd/=0 ) then 
            if(master_mpi) then
               if(fsmode)  call writefs(evlall,eferm,spinweightsoc)   !fermi surface !bias field vmag added  2023-9-20.
               !Obata bugfix add spinweightsoc at 2024-6-11
               write(stdo,*)' Writing bands to bands file for gnuplot ...'
               if(nsyml/=0)call writeband(evlall,eferm,vesav,evtop,ecbot,spinweightsoc) !bias field added  2023-9-20
            endif
            if(fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
            call rx0('plot band mode done') ! end of band plbnd/=0, that is, band plot mode.
         endif Writebandmode
      endif
    endblock BandPLOTmode
    GetFermiEnergy:block
      call m_subzi_bzintegration(evlall, eferm,sev,qvalm,vmag) !Get new eferm, vmag and wtkb, from evlall ;vmag is given for fsmom(fixedmoment) mode.
      evtop=maxval(evlall,evlall <eferm)
      ecbot=minval(evlall,evlall >eferm)
      emin =minval(evlall) !for plot
      if(lmet==0) eferm = (evtop+ecbot)/2d0 !for metal
    endblock GetFermiEnergy
    WriteEfermiFile: if(master_mpi) then
      if(lmet==0) write(stdo,"(' HOMO; Ef; LUMO =',3f11.6)")evtop,eferm,ecbot
      open(newunit=ifi,file='efermi.lmf')
      write(ifi,"(2d24.16, ' # (Ry) Fermi energy and Bias vmag; -vmag/2 +vmag/2 for each spin,')") eferm,vmag
      write(ifi,"(d24.16, ' # (Ry) Top of Valence')") evtop
      write(ifi,"(d24.16, ' # (Ry) Bottom of conduction')") ecbot
      write(ifi,"(2d24.16,' # (Ry) Efermi-Vesav Vesav')") eferm-vesaverage,vesaverage
      write(ifi,"(2d24.16,' #before: number of electrons total at sites:qval, backg:qbg=')") qval,qbg
      write(ifi,"(2d24.16,' # band: charge(nup+nspin), mag.mom(nup-ndown)')")qvalm(1), qvalm(2)
      write(ifi,"(f24.9 , ' # band: band energy sum (eV)=')") sev*rydberg()
      write(ifi,"(3i10,'    # Used k point to determine Ef')") nkabc
      write(ifi,"(i6,'# iter CAUTION! This file is overwritten by lmf SC loop')")iter
      close(ifi)
    endif WriteEfermiFile
    GenerateTotalDOS: if(master_mpi .AND. (tdos .OR. ldos/=0)) then !   emin=dosw(1) and emax=dosw(2) sets dos range
       dosw(1)= emin-0.01d0                          ! lowest energy limit to plot dos
       dosw(2)= min(maxval(evlall),eferm + 100d0/rydberg())  ! eferm +bz_dosmax !max energy limit to plot dos
       write(stdo,ftox)' bndfp:Generating TDOS: efermi(eV)=',ftof(rydberg()*eferm),&
            ' DOSwindow emin emax(eV)= ',ftof(rydberg()*dosw),'ltet nsp=',ltet,nsp
       allocate( dosi_rv(ndos,nsp),dos_rv(ndos,nsp),source=0d0) !for xxxdif
       if(cmdopt0('--tdostetf')) ltet= .FALSE. ! Set tetrahedron=F
       if(ltet) then
          nnn=nkabc(1)*nkabc(2)*nkabc(3)
          !do iq=1,nkp          !   do jsp=1,nspx
          !      write(stdo,ftox)'iq jsp=',iq,jsp,nevmin,nevls(iq,jsp)
          !      write(stdo,"('fp evl',8f8.4)")(evlall(i,jsp,iq),i=1,nevmin) !ls(iq,jsp))
          !   enddo          !enddo
          call bzints(nnn,evlall,dum,nkp, nevmin,nbandmx,nspx,nsp,dosw(1),dosw(2), dosi_rv,ndos,xxx,1,ntet,iv_a_oidtet,dumx,dumx,&
          !                                                                               job=1 give IntegratedDos to dosi_rv
               spinweightsoc) !2024-5-10           !write(stdo,ftox)'xxx dosi rv=',sum(dosi_rv)
          dos_rv(2:ndos-1,:)=(dosi_rv(3:ndos,:)-dosi_rv(1:ndos-2,:))/(2d0*(dosw(2)-dosw(1))/(ndos-1))
          dos_rv(1,:)    = dos_rv(2,:)
          dos_rv(ndos,:) = dos_rv(ndos-1,:)
          write(stdo,ftox)'iq jsp=',ndos,'dosw=',ftof(dosw),sum(dos_rv),sum(dosi_rv)
       else
          call makdos(nkp,nevmin,nbandmx,nsp,rv_a_owtkp,evlall,bz_n,bz_w,-6d0,dosw(1),dosw(2),ndos,dos_rv) !ndmahx=>nevmin
       endif
       head='.chk'
       if(tdos)head=''
       open(newunit=ifi, file= 'dos.tot.'//trim(sname)//trim(head))
       open(newunit=ifii,file='dosi.tot.'//trim(sname)//trim(head))
       dee=(dosw(2)-dosw(1))/(ndos-1d0)
       dosi=0d0
       do ipts=1,ndos
          eee= dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-eferm
          dosi(1:nsp)= dosi(1:nsp) + dos_rv(ipts,1:nsp)*dee
          write(ifi, "(255(f13.5,x))") eee,(dos_rv(ipts,isp),isp=1,nsp) ! dos
          write(ifii,"(255(f13.5,x))") eee,(dosi_rv(ipts,isp),isp=1,nsp)! integrated dos
       enddo
       close(ifi)
       close(ifii)
    endif GenerateTotalDOS
    if(tdos) call rx0('Done tdos mode:')
    if(master_mpi.and.afsym) write(stdo,ftox)' afsym mode: AF symmetry lets us make bands of isp=2 from isp=1!'
    WRITEeigenvaluesConsole: if(master_mpi .AND. iprint()>=35) then
      if(master_mpi) write(stdo,"('                 ikp isp         q           nev ndimh')")
      do iq  = 1,nkp
        qp  = qplist(:,iq)
        do jsp = 1,nspx
          write(stdo,ftox)' band-efermi(eV):',iq,jsp,ftof(qp,3),'nev=',nevls(iq,jsp),'ndimx=',ndimhx_(iq,jsp)
          write(stdo,"('  ',10f8.3)") rydberg()*(evlall(1:nevls(iq,jsp),jsp,iq)-eferm)
        enddo
      enddo
    endif WRITEeigenvaluesConsole
    GetDensity: if(lrout>0) then 
      call mpi_barrier(comm,ierr)
      call m_bandcal_2nd()  !accumulate products of eigenfunctions to obtian smooth density smrho_out, and so on.
      call m_bandcal_allreduce() 
    endif GetDensity
    CorelevelSpectroscopy2: if(cmdopt0('--cls')) then !m_clsmode_set1 is called in m_bandcal
      dosw(1)= emin  - 0.5d0     ! lowest energy limit to plot dos
      dosw(2)= eferm + bz_dosmax ! highest energy limit to plot dos
      call m_clsmode_finalize(eferm,ndimh,nbandmx,nspx,nkp,dosw,evlall)
      call rx0('Done cls mode:')
    endif CorelevelSpectroscopy2
    if(lso/=0)   call iorbtm() !Write Orbital Moment
    if(lrout/=0) call m_bandcal_symsmrho()  !Get smrho_out. Symmetrize smooth density ! Assemble output density, energies and forces 
    WRITEsmrhoTOxsf: if(cmdopt0('--density') .AND. master_mpi) then ! new density mode
       block
         use m_lmfinit,only:z
         integer::i1,i2,i3,i
         open(newunit=ifi,file='smrho.xsf')
         do isp = 1, nsp
           write(ifi,'("CRYSTAL")');     write(ifi,'("PRIMVEC")')
           write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3),i2=1,3)
           write(ifi,'("PRIMCOORD")');   write(ifi,'(2i5)') nbas,1
           do i = 1, nbas
             write(ifi,'(f8.3,2x,3f10.5)')z(ispec(i)),(pos(i2,i)*alat*0.529177208,i2=1,3)
           enddo
           write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
           write(ifi,'("charge_density_spin_",i1)') isp
           write(ifi,'("BEGIN_DATAGRID_3D_isp_",i1)') isp
           write(ifi,'(3i4)') n1,n2,n3;     write(ifi,'(3f10.5)') 0.,0.,0.
           write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3),i2=1,3)
           write(ifi,'(8e14.6)') (((dble(smrho_out(i1,i2,i3,isp)),i1=1,n1),i2=1,n2),i3=1,n3)
           write(ifi,'("END_DATAGRID_3D_isp_",i1)') isp
           write(ifi,'("END_BLOCK_DATAGRID_3D")')
         enddo
         close(ifi)
       endblock
    endif WRITEsmrhoTOxsf
    call m_mkrout_init() !Get force frcbandsym, symmetrized atomic densities orhoat_out, and weights hbyl,qbyl
    if(.not.pnufix) call pnunew(eferm)   !pnuall are revised. !  New boundary conditions pnu for phi and phidot
    call m_mkehkf_etot1(sev, eharris) !Evaluate_HarrisFoukner_energy (note: we now use NonMagneticCORE mode as default)
    eksham = 0d0
    if(lrout/=0) then
      EvaluateKohnShamTotalEnergyAndForce:block
         real(8):: qmom_in(nlmxlx,nbas)
         qmom_in=qmom !multipole moments.
         eksham = 0d0 !   ... Evaluate KS total energy and output magnetic moment
         if(leks>=1) then
            call mkekin(osig,otau,oppi,oqkkl,vconst,osmpot,smrho_out,sev,  ekinval)
            call m_mkpot_energyterms(smrho_out, orhoat_out) !qmom is revised for given orhoat_out
            if(cmdopt0('--density')) then
               call mpi_barrier(comm,ierr)
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
            call dfrce (lfrce,orhoat,orhoat_out,qmom_in,smrho,smrho_out,  fh_rv)
            call totfrc(leks, fes1_rv, fes2_rv, fh_rv, frcbandsym, force) ! force : total
         endif
         ! Mix inputs(smrho,orhoat) and outputs(smrho_out,orhoat_out), resulting orhoat and smrho.
         call mixrho(iter,qval-qbg,orhoat_out,orhoat,smrho_out,smrho,qdiff)!mixrho keeps history in it.
         call mpi_barrier(comm,ierr) 
         ! write(stdo,ftox)'mixrho: output smrho =',maxval(dreal(smrho)),minval(dreal(smrho)),sum(dreal(smrho))
         ! do ib=1,nbas!write(stdo,ftox)'mixrho: output orhoat=',ib,sum(abs(orhoat(1,ib)%v)),sum(abs(orhoat(2,ib)%v)),sum(abs(orhoat(3,ib)%v))
         ! enddo   
       endblock EvaluateKohnShamTotalEnergyandForce
    endif
    if(magexist) then                   ! Energy by externel mag field = -nup*vmag/2 + ndn*vmag/2 = -(nup-ndn)*vmag/2= -amom*vmag/2
       eharris= eharris + amom*vmag/2d0 ! Subtract this energy 
       eksham = eksham  + amom*vmag/2d0 
    endif
    ham_ehf= eharris  !Harris total energy
    ham_ehk= eksham   !Hohenberg-Kohn-Sham total energy
    call m_mkpot_deallocate()
    call m_bandcal_clean()
    call tcx('bndfp')
  end subroutine bndfp
end module m_bndfp
