      subroutine lmfp(llmfgw)
      use m_lmfinit,only: irs1,irs2,lhf,bz_nevmx,bz_lmet, irs1x10, 
     &     ctrl_mdprm,ctrl_maxit,nbas,nl,nspec,nsp,
     &     qbg=>zbak,ctrl_ldos,ctrl_nitmv,ctrl_nvario,ctrl_quit,ctrl_tol,
     &     ham_seref,ctrl_lfrce,
     &     lat_gam,lat_alat, lat_ldist,
     &     sspec=>v_sspec, ssite=>v_ssite, 
     &     sstrnsymg, symgaf, nat, lldau,nlibu,lmaxu,stdo,lrout,leks,plbnd
      use m_lattic,only: lat_dist,lat_plat,qlat=>lat_qlat,rv_a_opos,rv_a_odlv,rv_a_oqlv,
     &     lat_vol,lat_nkd,lat_awald,lat_nkq
      use m_mksym,only: rv_a_oag,rv_a_osymgr,iv_a_oistab,lat_nsgrp
      use m_ext,only: sname
      use m_MPItk,only: procid,master,mlog,master_mpi
      use m_iors,only: Iors
      use m_supot,only: m_supot_init,k1,k2,k3
      use m_suham,only: m_suham_init
      use m_bndfp,only: Bndfp,   ham_ehf,ham_ehk,qdiff
      use m_mkrout,only: frc_rv=>frc
      use m_ldau,only: m_ldau_init,m_ldau_vorbset
      use m_chkdmu,only: eorb
      use m_bandcal,only: dmatu
      use m_qplist,only: m_qplist_init
      use m_mkpot,only: amom 
      use m_writeham,only: m_writeham_init, m_writeham_write
      use m_lmfgw,only: m_lmfgw_setjobgw
      use m_hamindex,only: m_hamindex_init
c      use m_uspecb,only:uspecb
!! = Main routine of lmf = (following document is roughly checked at May2021)
!! lmfp contains two loops after initialization
!!   1  outer  loop do 2000 is for molecular dynamics (relaxiation).
!!   2. innner loop do 1000 is for electronic structure self-consistency.
!!      Main part of band calculaiton is in bndfp. 
!! (Most of) all data in modules are 'protected'.
!! Thus data in m_lmfinit, m_lattic, m_mksy, m_ext, ... are fixed during iteration.
!! Currently rhoat, smrho, smpot in m_supot are iterated. ham_ehf _ehk changed by iterations.      
!! aug2020. T.kotani removed lshr mode (automatic modification of plat), because
!! We don't maintain the share mode.
!!      Probably, we need to re-design it (maybe outside of fortran code).
!     sspec :struct for species-specific information; see routine uspec
!     ssite :struct for site-specific information; see routine usite
!   sspec and ssite are not clearly protecetd ---> iors can modify it.
!   But, these are fixed before going into do 2000/1000 loop (I think).      
! Memo for LDA+U
!     lmaxu : max l for a U (used for dimensioning)
!     nlibu : total number of U blocks
!     lldau(ib) : U on site ib with dmat beginning at dmats(*,lldau(ib))
!! Unused irlxsh: counter for shear relaxations.  irlxsh=0 => new step
!!      
!! --> history is removed to avoid confusions. See ecalj@github  
      implicit none
      character strn*120, alabl*8, flg*3
      logical :: cmdopt,llmfgw,xyzfrz(3),lbin,cmdopt0 !,lshr=.false.
      integer :: i,ifi,ipr,irs1x, k, nit1,numq, pdim,lsc
      integer :: icom, natrlx , nvrelx , itrlx, irlxsh=0, nitrlx 
      integer:: ibas,unlink,ifipos,ifile_handle, maxit,iter,j,idmatu,iprint
      real(8) :: mdprm(6),gam(4),gam1,bstim,pletot(6,2)
      real(8) :: plat(3,3),xvcart(3),xvfrac(3)
      real(8):: seref,etot(2),qtol,etol
      integer, allocatable :: indrx_iv(:,:) ! for relaxation
      real(8), allocatable :: pos_move(:,:),  ftot_rv(:), wk_rv(:), p_rv(:),w_rv(:)
      real(8), allocatable :: rv_a_omad (:) !   Madelung matrix if necessary
c      complex(8),allocatable::vorbdmat(:),dmatu(:),dmatoa(:) !LDA+U
      integer,parameter:: nm=3
      character(512):: aaachar
      character(10):: i2char
      logical:: writeham
      integer:: jobgw
#if (MPI |MPIK)
      integer:: ierr
      include "mpif.h"
#endif
!! Sep2020 takao: --rs=3 mode is removed. (--rs=3 meand fixed density Harris-foukner MD).
!! --rs=3 => always read from atom file !2020aug.
!! Sep2020 comment " Shorten site positions" removed.
      
      call tcn('lmfp')
      etot(1) = 0 ! Total energy mode --etot ==>moved to m_lmfinit ---
      ipr=iprint() !call getpr(ipr) !ipr is pring verbose control
      maxit = ctrl_maxit
      irs1x = irs1
      if(ipr>=30) call praugm(sspec,0) !Printout properties of species
      
!! atomic position Relaxation setup (MD mode)
      itrlx = 1
      nitrlx= ctrl_nitmv  !num of iteration cycle of relaxiation
      mdprm = ctrl_mdprm  ! c       lshr = nint(mdprm(1)) .gt. 100 !lshr mode switch 
      if( nint(mdprm(1))== 0) nitrlx = 0 !no relaxiation
      if( nint(mdprm(1))>0 .and. nint(mdprm(1))<4 ) call rx('lmf not set up for MD yet')
      if( nitrlx>0 ) then  ! nitrlx >0 is for atomic position relaxiation
        allocate(indrx_iv(2,3*nbas))
        call rlxstp(natrlx, nvrelx, indrx_iv,  xyzfrz, pdim )
        icom = 0
        if (nvrelx .ne. 0) then
          allocate(w_rv(nvrelx*nvrelx))
          allocate(p_rv(pdim))
        endif
        if(master_mpi) then
           open(newunit=ifipos,file='AtomPos.'//trim(sname),form='unformatted',status='new')
           write(ifipos) nbas
           write(ifipos) itrlx,rv_a_opos
        endif   
      endif
      
!!Get jobgw for lmfgw mode
      jobgw=0
      if(master_mpi.and.llmfgw) then
         write(stdo,*)
         write(stdo,*) '=== lmfgw-MPIK: Choose one of following jobs: ==='
         write(stdo,*) '   0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima'
         write(stdo,*) '   1 : GW setup mode'
         write(stdo,*) ' job 0 or 1?'
         read (5,*) jobgw
      endif
      call mpibc1(jobgw,1,2, mlog,'lmfp','jobgw')
      call m_lmfgw_setjobgw(jobgw)
      
!! array allocated in supot rhoat smrho     
      call m_supot_init() ! get G vectors for charge
      call m_suham_init(sspec,ssite)  ! sspec ssite are modified.
      call m_hamindex_init(jobgw==0) !index for hamiltonian gen_hamindex
      
      writeham= cmdopt0('--writeham')
      if(writeham) call m_writeham_init()
      if(writeham) call m_writeham_write()
      
      if( ctrl_quit==8 ) call rx0('quit = ham')
      if (lrout == 0 .and. ctrl_lfrce /= 0) then
        write(stdo,"(a)") 'Error: output density required when forces sought.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call rx('incompatible input. see the end of console output')
      endif
      if (lrout == 0 .and. cmdopt0('--etot')) then
        write(stdo,"(a)") 'Error: output density required with --etot switch.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call rx('incompatible input. see the end of console output')
      endif
      if (lrout == 0 ) maxit = 1
!! LDA+U initialization
      if(nlibu>0) call m_ldau_init()
      if(nlibu>0.and.ctrl_quit==16 ) call rx0('quit = dmat')
      
!!=== Re-entry for lshr (shear mode); but shear mode is currently commented out.
!!   (probably shear mode should be recovered outside of fortran)
      if(master_mpi) write(stdo,*)' lmfp : potential and basis setup ... '
      if(.false.) then 
        allocate(rv_a_omad(nbas**2))
        call madmat (rv_a_omad ) !Monopole Madelung matrix (kept for future).
      endif
!! Initial density Read atm.* or rst.*   rhoat smrho in supot
!!  --- Read restart file or overlap free atom densities ---
!!    irs1x tells what to read and whether to invoke smshft.
C     0    read from atom file  atm
C     1    read from binary     rst 
C     2    read from ascii      rsta
C    +10 -> invoke smshft(1) after file read.
      if (irs1x == 0) then ! read overlap free-atom densities
         call rdovfa()     ! superposition of spherical atomic densities by lmfa.
         nit1 = 0
      elseif( mod(irs1x,4)==1 .or. mod(irs1x,4)==2 ) then 
         lbin = (irs1x==1) ! lbin=T:rst for --rs=1, or lbin=F:rsta for --rs=2,
         k = -1
         if(master_mpi) then
            ifi=ifile_handle() 
            if(lbin)       open(ifi,file='rst.'//trim(sname),form='unformatted')
            if(.not. lbin) open(ifi,file='rsta.'//trim(sname)) 
         endif
         call mpibc1(ifi,1,2,mlog,'lmfp','ifi')
         k = iors(nit1 , lbin , ifi ) ! read rst file. sspec ssite modified  !bz_w in iros removed
         if (master_mpi) close(ifi) 
         call mpibc1(k,1,2,mlog,'lmfp','k')
         do ibas=1,nbas 
            rv_a_opos(:,ibas) = ssite(ibas)%pos
         enddo
         if(k<0) then !initial potential from atm file if rst can not read
            irs1x = 0
            call rdovfa() 
            nit1 = 0
         elseif (irs1x10) then  
            call smshft(1)
         endif
      endif
      
!!  tolerances for convergence check. q charge and e energy
      qtol = ctrl_tol(1)
      etol = ctrl_tol(3)
!!      
      if(cmdopt0('--vbmonly')) then !Get VBM and CBM relative to vaccum (a simple approximaiton. Need fixing if necessary.).
         if(master_mpi) call vbmmode()
         call rx0('--vbmonly mode done')
      endif
      if(cmdopt0('--getq')) then ! Current version is not for spin dependent, with many restrictions.
         if(master_mpi) call getqmode()
         call rx0('--getq mode done')
      endif
      
!! Get q point list at which we do band calculations.
!! we can move this outside of iteration cycle.
      if(.not.llmfgw) call m_qplist_init(plbnd)
      
!! === Main iteration loops ===
      if( nitrlx>=1 ) allocate(pos_move(3,nbas))
      do 2000 itrlx = 1,max(1,nitrlx) ! loop for atomic position relaxiation (molecular dynamics) ===
         if( ipr>=30 ) then ! Write atom positions 
           write(stdo,"(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x,
     &     'pos (multiples of plat)')") 'Basis, after reading restart file'
           do i = 1, nbas
           xvcart = ssite(i)%pos
           xvfrac = matmul(transpose(qlat),xvcart)
           alabl  = sspec(ssite(i)%spec)%name
           write(stdo,"(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)") i, alabl, xvcart,xvfrac
           enddo
         endif
         do 1000 iter = 1,max(1,maxit) ! for electronic structure ===
!     ! Make and diagonalize hamiltonian, make new charge density ---
            if(maxit/=0) then
               if (master_mpi) then
               aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
               write(stdo,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
               endif
               call bndfp(iter, llmfgw) !get total energies ham_ehf and ham_ehk
            endif
!! Check convergence of dmatu (density matrix for LDA+U) and update it. Get vorbdmat (U potential)
            if(nlibu>0.and.lrout>0) call m_ldau_vorbset(ham_ehk,dmatu) !set new vorb from dmat by bndfp-mkpot
!     --density (plot density) is moved to locpot.F(rho1mt and rho2mt) and mkpot.F(smooth part).
!! Write restart file (skip if --quit=band) ---
            if(master_mpi .and. ctrl_quit/=4) then
               if(irs2>0 .and. (lrout>0.or.maxit==0)) then ! .and. irlxsh==0) then
                  lbin = irs2/=2
                  ifi = ifile_handle() 
                  if(lbin)      open(ifi,file='rst.'//trim(sname),form='unformatted')
                  if(.not.lbin) open(ifi,file='rsta.'//trim(sname)) 
                  k = iors ( iter , lbin , - ifi )
                  close(ifi) 
               endif
            endif
!     ! Add to save file; decide on next iteration ---
            if (maxit<=0) goto 9998
            etot(1) = ham_ehf
            etot(2) = ham_ehk + eorb !eorb by LDA+U
            seref   = ham_seref !   ... reference energy
            etot(1) = etot(1) - seref
            if (etot(2)/=0) etot(2) = etot(2) - seref
            flush(6)
            call mpi_barrier(MPI_COMM_WORLD,ierr)
            if(master_mpi) then
               i = 0
!     nwit : kino's memo
!     lsc   :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
!     :1 if not self-consistent, but encountered max. no. iter.
!     :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
!     :3 otherwise
               call nwit( int(ctrl_nvario), iter, maxit, lhf.or.irs1x==0.and.iter==1,
     &              leks+i, etol, qtol, qdiff, 'cxhi', amom, etot, lsc)
c     &              leks+i, etol, qtol, qdiff dmxp(11), 'cxhi', amom, etot, lsc)
            endif
            call mpibc1(lsc,1,2,mlog,'lmfp','lsc')
            if (lsc==2 .and. (.not. lhf) .and. maxit>1) lsc = 3
            if (lsc==1 .and. lrout>0 .or. lsc==3) then
               if (iter .ge. maxit) lsc = 1
               if (iter .lt. maxit) lsc = 3
            endif
            if( ctrl_quit==4 ) call rx0 ( 'lmf : exit (--quit=band)' )
            if( lsc <= 2) exit  !self-consistency exit
 1000    continue               ! ---------------- SCF (iteration) loop end ----
         if(nitrlx==0) exit !no molecular dynamics (=no atomic position relaxation)
!! --- Molecular dynamics. not maintained recently but atomic position relaxation was working
         do ibas=1,nbas
            ssite(ibas)%pos0 = ssite(ibas)%pos
            pos_move(:,ibas) = ssite(ibas)%pos
         enddo
         mdprm = ctrl_mdprm !MD setting
!! Relax atomic positions.   !--> shear mode removed. probably outside of fortran code if necessary
         call relax( ssite , sspec , itrlx , indrx_iv 
     .        , natrlx, nvrelx, frc_rv , p_rv , w_rv , 0 , 0d0 , pos_move, icom )
!!  warn: Updating positions in ssite structure ==> t.kotani think this is confusing because
!!        'positions written in ctrl' and 'positions written in rst' can be different.
         if(master_mpi) write(ifipos) itrlx,pos_move
         do ibas=1,nbas
            ssite(ibas)%pos = pos_move(:,ibas) !warn
         enddo
C     ... Exit when relaxation converged or maximum number of iterations
         if (icom==1) then !takao 
            if (master_mpi) then
               flg = 'C67' !what? 
               call nwitsv(1+2,ctrl_nvario,flg,nsp,amom,etot)
               call tcx('lmfp')
               call fexit(0,111,' LMFP: relaxation converged after %i iteration(s)',itrlx)
            else
               call tcx('lmfp')
               call fexit(0,111,' ',0)
            endif
         endif
C     ... Restore minimum gradient positions if this is last step
         if (itrlx==nitrlx) then
            write(stdo,"(a)")' lmfp: restore positions for minimum g'
            call prelx1(1 , nm , .false. , natrlx , nvrelx , indrx_iv , p_rv, pos_move )
            do ibas=1,nbas !updated positions in site structure
               ssite(ibas)%pos = pos_move(:,ibas)
            enddo
         endif
!! New density after atom shifts 
!! If explicitly told to read from atom files after atom movmment
         call smshft(ctrl_lfrce) !sctrl,sham,slat,
!! Write restart file (to include new positions)
         if (procid .eq. master) then ! .and. .not. lshr) then
            ifi = ifile_handle()
            open(ifi,file='rst.'//trim(sname),form='unformatted') 
            k = iors (iter , .true., - ifi )
            close(ifi)
         endif
         if (procid .eq. master) then
            write(stdo,*)' Delete mixing and band weights files ...'
            ifi = unlink('mixm.'//trim(sname))
            ifi = unlink('wkp.'//trim(sname))
         endif
         call parms0(0,0,0d0,0) !   reset mixing block
!! Exit when maximum number of iterations encountered
         if (itrlx .eq. nitrlx) then
            if (procid .eq. master) then
               call tcx('lmfp')
               call fexit(1,111,' LMFP: relaxation incomplete after %i iteration(s)',nitrlx)
            else
               call tcx('lmfp')
               call fexit(1,111,' ',0)
            endif
         endif
 2000 continue
!! Write positions to file
 9998 continue
      if (cmdopt('--wpos=',7,0,strn) .or. cmdopt('--wpos:mode1:',13,0,strn)) then
         allocate(pos_move(3,nbas))
         do ibas=1,nbas
            pos_move(:,ibas) = ssite(ibas)%pos
         enddo
         gam=lat_gam
         allocate(wk_rv(3*nbas))
         call rdistn (pos_move , wk_rv , nbas , gam (1) , gam(2) , gam ( 3 ) , 1 / gam1 )
         if (allocated(wk_rv)) deallocate(wk_rv)
         if (allocated(pos_move)) deallocate(pos_move)
      endif 
      if(master_mpi.and.nitrlx>0) close(ifipos)
c      if(allocated(ftot_rv)) deallocate(ftot_rv)
      if(allocated(p_rv)) deallocate(p_rv)
      if(allocated(w_rv)) deallocate(w_rv)
      if(allocated(indrx_iv)) deallocate(indrx_iv)
      call tcx('lmfp')
      end subroutine lmfp
