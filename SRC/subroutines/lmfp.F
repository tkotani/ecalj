      subroutine lmfp(llmfgw,irs1x)
      use m_lmfinit,only: irs2,lhf, irs1x10, 
     &     ctrl_mdprm,ctrl_maxit,nbas,nsp,
     &     ctrl_ldos,ctrl_nitmv,ctrl_nvario,ctrl_tol,
     &     ham_seref,ctrl_lfrce,  
     &     sspec=>v_sspec, ssite=>v_ssite, 
     &     nlibu,stdo,lrout,leks,plbnd,lpzex, nitrlx
     &     , indrx_iv,natrlx,xyzfrz,pdim
      use m_lattic,only: qlat=>lat_qlat,rv_a_opos
      use m_chkdmu,only:  eorb
      use m_bandcal,only: dmatu
      use m_mkpot,only:   amom 
      use m_ext,only:     sname
      use m_iors,only:   Iors
      use m_MPItk,only:  mlog,master_mpi
      use m_bndfp,only:  Bndfp,   ham_ehf,ham_ehk,qdiff,force,sev
      use m_ldau,only:   m_ldau_vorbset
      use m_bstrux,only: m_bstrux_init
!! = Main routine of lmf = (following document is roughly checked at May2021)
!! lmfp contains two loops after initialization
!!   1  outer  loop do 2000 is for molecular dynamics (relaxiation).
!!   2. innner loop do 1000 is for electronic structure self-consistency.
!!      Main part of band calculaiton is in bndfp. 
!! (Most of) all data in modules are 'protected'.
!! Thus data in m_lmfinit, m_lattic, m_mksy, m_ext, ... are fixed during iteration.
!! Currently rhoat, smrho, smpot in m_supot are iterated. ham_ehf _ehk changed by iterations.      
!! We currently don't maintain the share (of crystal structure) mode.
!     sspec :struct for species-specific information; see routine uspec
!     ssite :struct for site-specific information; see routine usite
!     sspec and ssite are not protecetd ---> iors can modify it.
!     But, these are fixed before going into do 2000/1000 loop (I think).      
! Memo for LDA+U. See m_ldau module.
!     lmaxu : max l for a U (used for dimensioning)
!     nlibu : total number of U blocks
!     lldau(ib) : U on site ib with dmat beginning at dmats(*,lldau(ib))
!! Unused irlxsh: counter for shear relaxations.  irlxsh=0 => new step
!!      
!! --> history is removed to avoid confusions. See ecalj@github  
!! aug2020. T.kotani removed lshr mode (automatic modification of plat), because
!!      Probably, we need to re-design it (maybe outside of fortran code).
      implicit none
      character strn*120, alabl*8, flg*3
      logical :: cmdopt,llmfgw,lbin,cmdopt0 !,lshr=.false.
      integer :: i,ifi,ipr,irs1x, k, nit1,numq, lsc
      integer :: icom,  nvrelx , itrlx, irlxsh=0 
      
      integer:: ibas,unlink,ifipos,ifile_handle, maxit,iter,j,idmatu,iprint
      real(8) :: gam(4),gam1,bstim,pletot(6,2)
      real(8) :: plat(3,3),xvcart(3),xvfrac(3)
      real(8):: seref,etot(2),qtol,etol
      real(8), allocatable :: pos_move(:,:),  ftot_rv(:), wk_rv(:), p_rv(:),hess(:,:)
      real(8), allocatable :: rv_a_omad (:) !  Madelung matrix if necessary
      integer,parameter:: nm=3
      character(512):: aaachar
      character(10):: i2char
      logical:: cmdopt2
      character(20):: outs=''
      integer:: ierr
      include "mpif.h"
!! --rs=3, always read from atom file, mode is removed. (--rs=3 is for fixed density Harris-foukner MD).
!! Sep2020 comment " Shorten site positions" removed. (we are using shortn3 mainly now)
      call tcn('lmfp')
      ipr = iprint() !call getpr(ipr) !ipr is pring verbose control
!!  tolerances for convergence check. q charge and e energy
      qtol = ctrl_tol(1)
      etol = ctrl_tol(3)
      etot(1) = 0 ! Total energy mode --etot ==>moved to m_lmfinit ---
      maxit = ctrl_maxit
      if(lrout == 0 ) maxit = 1
!! Atomic position Relaxation setup (MD mode)
      if( nitrlx>0 ) then 
        icom = 0
        if(natrlx /= 0) allocate(hess(natrlx,natrlx),p_rv(pdim))
        if(master_mpi) then
           open(newunit=ifipos,file='AtomPos.'//trim(sname),form='unformatted',status='new')
           write(ifipos) nbas
           write(ifipos) 0,rv_a_opos
        endif   
      endif
!! === Main iteration loops ===
      if( nitrlx>=1 ) allocate(pos_move(3,nbas))
      do 2000 itrlx = 1,max(1,nitrlx) ! loop for atomic position relaxiation (molecular dynamics) ===
!! Get all structure constants for nbas and qplist
         if(sum(lpzex)==0) call m_bstrux_init() !this reads ssite%pos %pz and so on.
              ! We can make structure constant (C_akL Eq.(38) in /JPSJ.84.034702) here
              ! if we have no extended local orbital.
              ! When we have extentede local obtail, we run m_bstrux_init after elocp in mkpot.
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
         do 1000 iter = 1,max(1,maxit) ! loop for electronic structure ===
!     ! Make and diagonalize hamiltonian, make new charge density ---
            if(maxit/=0) then
               if (master_mpi) then
                  aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
                  write(stdo,*)
                  write(stdo,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
               endif
               call bndfp(iter,llmfgw,plbnd) !Main part of band calculation. Get total energies ham_ehf and ham_ehk
            endif
!! Check convergence of dmatu (density matrix for LDA+U) and update it. Get vorbdmat (U potential)
            if(nlibu>0.and.lrout>0) call m_ldau_vorbset(ham_ehk,dmatu) !set new vorb from dmat by bndfp-mkpot
!   Things for --density (plot density) are in locpot.F(rho1mt and rho2mt) and mkpot.F(smooth part).
!! Write restart file (skip if --quit=band) ---
            if(master_mpi .and. (.not.cmdopt0('--quit=band'))) then
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
! lsc   :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
!     :1 if not self-consistent, but encountered max. no. iter.
!     :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
!     :3 otherwise
               call nwit( int(ctrl_nvario), iter, maxit, lhf.or.irs1x==0.and.iter==1,
     &              leks+i, etol, qtol, qdiff, 'cxhi', amom, etot, sev,lsc)
            endif
            call mpibc1_int(lsc,1,'lmfp_lsc')
            if (lsc==2 .and. (.not. lhf) .and. maxit>1) lsc = 3
            if (lsc==1 .and. lrout>0 .or. lsc==3) then
               if (iter .ge. maxit) lsc = 1
               if (iter .lt. maxit) lsc = 3
            endif
            if( cmdopt0('--quit=band')) call rx0('lmf-MPIK : exit (--quit=band)')
            if( lsc <= 2) exit  !self-consistency exit
 1000    continue               ! ---------------- SCF (iteration) loop end ----
         if(nitrlx==0) exit !no molecular dynamics (=no atomic position relaxation)
!! --- Molecular dynamics. not maintained recently but atomic position relaxation was working
         do ibas=1,nbas
            ssite(ibas)%pos0 = ssite(ibas)%pos
            pos_move(:,ibas) = ssite(ibas)%pos
         enddo
c         mdprm = ctrl_mdprm !MD setting
!! Relax atomic positions.   !--> shear mode removed. probably outside of fortran code if necessary
         call relax( ssite , sspec , itrlx , indrx_iv 
     .        , natrlx, force , p_rv , hess , 0 , 0d0 , pos_move, icom )
!!  warn: Updating positions in ssite structure ==> t.kotani think this is confusing because
!!        'positions written in ctrl' and 'positions written in rst' can be different.
         if(master_mpi) write(ifipos) itrlx,pos_move
         do ibas=1,nbas
            ssite(ibas)%pos = pos_move(:,ibas) 
         enddo
C     ... Exit when relaxation converged or maximum number of iterations
         if(icom==1) then !takao 
            if(master_mpi) then
               flg = 'C67' !what? 
               call nwitsv(1+2,ctrl_nvario,flg,nsp,amom,etot,sev)
               write(stdo,"(a,i5)")' LMFP: relaxation converged after iteration(s) of ',itrlx
            endif
            exit
         endif
C     ... Restore minimum gradient positions if this is last step
         if (itrlx==nitrlx) then
            write(stdo,"(a)")' lmfp: restore positions for minimum g'
            call prelx1(1 , nm , .false. , natrlx , indrx_iv , p_rv, pos_move )
            do ibas=1,nbas !updated positions in site structure
               ssite(ibas)%pos = pos_move(:,ibas)
            enddo
         endif
!! New density after atom shifts.  
!! If explicitly told to read from atom files after atom movement
         call smshft(ctrl_lfrce) 
!! Write restart file (to include new positions)
         if (master_mpi) then ! .and. .not. lshr) then
            ifi = ifile_handle()
            open(ifi,file='rst.'//trim(sname),form='unformatted') 
            k = iors (iter , .true., - ifi )
            close(ifi)
         endif
         if ( master_mpi) then
            write(stdo,*)' Delete mixing and band weights files ...'
            ifi = unlink('mixm.'//trim(sname))
            ifi = unlink('wkp.'//trim(sname))
         endif
         call parms0(0,0,0d0,0) !   reset mixing block
         if(itrlx==nitrlx.and.master_mpi) write(stdo,"(a)")' LMFP: relaxation incomplete'
 2000 continue
 9998 continue
      if(master_mpi.and.nitrlx>0) close(ifipos)
      if(allocated(p_rv)) deallocate(p_rv)
      if(allocated(hess)) deallocate(hess)
      call tcx('lmfp')
      end subroutine lmfp
