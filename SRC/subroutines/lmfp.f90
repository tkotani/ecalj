subroutine lmfp(llmfgw)
  use m_lmfinit,only: lhf, maxit,nbas,nsp, ctrl_ldos,ctrl_nitmv,ctrl_nvario, &
       ham_seref,ctrl_lfrce,  sspec=>v_sspec, ssite=>v_ssite, &
       nlibu,stdo,lrout,leks,plbnd,lpzex, nitrlx, &
       indrx_iv,natrlx,xyzfrz,pdim,qtol,etol
  use m_lattic,only: qlat=>lat_qlat,rv_a_opos
  use m_bandcal,only: dmatu
  use m_mkpot,only:   amom
  use m_ext,only:     sname
  use m_iors,only:     Iors
  use m_iors_old,only: Iors_old !only for reading vs=1.04 rst file (before 2022-5-14)
  use m_MPItk,only:  mlog,master_mpi
  use m_bndfp,only:  Bndfp,   ham_ehf,ham_ehk,qdiff,force,sev
  use m_ldau,only:   M_ldau_vorbset, eorb
  use m_bstrux,only: M_bstrux_init
  use m_relax,only:  Relax,Prelx1
  use m_mixrho,only: Parms0
  use m_lattic,only: Setopos
  use m_ftox
  !!= Main routine of lmf = (following document is roughly checked at May2021)
  !! lmfp contains two loops after initialization
  !!   1  outer  loop do 2000 is for molecular dynamics (relaxiation).
  !!   2. innner loop do 1000 is for electronic structure self-consistency.
  !!      Main part of band calculaiton is in bndfp.
  !! (Most of) all data in modules are 'protected'.
  !! Thus data in m_lmfinit, m_lattic, m_mksy, m_ext, ... are fixed during iteration.
  !! Currently rhoat, smrho, smpot in m_supot are iterated. ham_ehf _ehk changed by iterations.
  !     sspec :struct for species-specific information. I think fixed during iteration. See m_lmfinit.
  !     ssite :struct for site-specific information
  ! Memo for LDA+U. See m_ldau module.
  !     lmaxu : max l for a U (used for dimensioning)
  !     nlibu : total number of U blocks
  !     lldau(ib) : U on site ib with dmat beginning at dmats(*,lldau(ib))
  !!
  !! We currently don't maintain the share (of crystal structure) mode.
  !! ===> history is removed to avoid confusions. See ecalj@github
  !! aug2020. T.kotani removed lshr mode (automatic modification of plat), because
  !!      Probably, we need to re-design it (maybe outside of fortran code).
  implicit none
  integer,parameter:: nm=3
  character alabl*8, flg*3
  logical :: cmdopt,llmfgw,lbin,cmdopt0 !,lshr=.false.
  integer :: i,ifi,ipr, k, nit1,numq, lsc, icom,  nvrelx , itrlx
  integer:: ibas,unlink,ifipos,ifile_handle,iter,j,idmatu,iprint
  real(8) :: gam(4),gam1,bstim,pletot(6,2), plat(3,3),xvcart(3),xvfrac(3),seref,etot(2),vs,vs1
  real(8),allocatable :: pos_move(:,:),  ftot_rv(:), wk_rv(:), p_rv(:),hess(:,:)
  real(8),allocatable :: rv_a_omad (:) !  Madelung matrix if necessary
  character(512):: aaachar
  character(10):: i2char
  logical:: cmdopt2
  character(20):: outs=''
  integer:: ierr
  integer:: irs(5),irsr(5), irs1,irs2,irs3,irs5
  logical:: irs1x10!,irsrot
  character(256):: strn,strn2
  include "mpif.h"
  call tcn('lmfp')
  ipr = iprint()
  != Reading condition switch [irs1,irs2,irs3,irs5,irs1x10]
  !=    irs1x tells what to read and whether to invoke smshft.
  !=    0    read from atom file  atm;  1  read from binary  rst ;  2 read from ascii  rsta
  !=    +10 -> invoke smshft(1) after file read.
  != --rs=3 mode is removed.(always read from atom file, --rs=3 is for fixed density Harris-foukner MD).
  irs=[1,1,0,0,0]
  if (cmdopt2('--rs=',strn)) then
     strn2=trim(strn)//',999,999,999,999,999'
     read(strn2,*) irsr
     do i=1,5
        if(irsr(i)/=999 ) irs(i)=irsr(i)
     enddo
  endif
  irs1 = mod(irs(1),10) ! irs1 is irrelenvant now . irs1>10 gives irs1x10
  irs2 = irs(2)         ! irs2=0 (no write)
  if(irs1/=1 .AND. irs1/=0 .AND. irs1/=1) call rx('irs1 not \in [0,1,2]')
  if(irs2/=0 .AND. irs2/=1) call rx('irs2/=0 .AND. irs2/=1')
  irs3 = irs(3) ! irs3=1: read site positions from ctrl even when we have rst.
  irs5 = irs(5) ! irs5=1: read pnu from ctrl even when we have rst.
  irs1x10 = (irs(1)/10==1)  ! +10:  smshft after rst/rsta
  if (cmdopt0('--etot')) irs2=0 !not write rst files
  ! xx irsrot  = (irs(1)/100==1) ! +100: rotate local density after file read for shear mode (iors.F)
  ! xx irs4 = irs(4) ! read starting fermi level from ctrl

  !=  Initial density set up. Read atomic and smooth parted os density,
  !=  rhoat smrho in m_density from atm.* or rst.*
  !     Read rst
  k=-1
  vs=2d0
  if(irs1/=0) then
     if(master_mpi) then
        open(newunit=ifi,file='rst.'//trim(sname),form='unformatted')
        read(ifi,end=996,err=996) vs !version id of rst file
        close(ifi)
        write(stdo,ftox)'lmv7: Read rst version ID=',ftof(vs,2)
996     continue
     endif
     call Mpi_barrier(MPI_COMM_WORLD,ierr)
     call Mpibc1_real(vs,1,'lmv7: vs: version id of rst file')
     if(vs==2d0) then       !2020-5-14
        k = iors(nit1,'read',irs3=irs3,irs5=irs5) ! read rst file. sspec ssite maybe modified
     else                   !vs=1.04d0
        k = iors_old(nit1,'read',irs3=irs3,irs5=irs5) ! read rst file. sspec ssite maybe modified
     endif
     call Mpibc1_int(k,1,'lmv7:lmfp_k')
     call Setopos()         ! Set position of atoms read from iors
  endif
  if(k<0) then !irs1==0 .OR. Not reading rst.
     irs1 = 0
     call Rdovfa()  ! Initial potential from atm file if rst can not read
     nit1 = 0
  endif
  if(k>=0 .AND. irs1x10) call Smshft(1)  ! modify denity after reading rst when irs1x10=True

  !! Sep2020 " Shorten site positions" removed.
  etot = 0d0 ! Total energy mode --etot ==>moved to m_lmfinit ---
  if( nitrlx>0 ) then ! Atomic position Relaxation setup (MD mode)
     icom = 0
     if(natrlx /= 0) allocate(hess(natrlx,natrlx),p_rv(pdim))
     if(master_mpi) then
        open(newunit=ifipos,file='AtomPos.'//trim(sname),form='unformatted',status='new')
        write(ifipos) nbas
        write(ifipos) 0,rv_a_opos
     endif
     allocate(pos_move(3,nbas))
  endif

  !==== Main iteration loops ===
  do 2000 itrlx = 1,max(1,nitrlx) ! loop for atomic position relaxiation (molecular dynamics) ===
     !     Get all structure constants for nbas and qplist
     if(sum(lpzex)==0) call M_bstrux_init() !this reads ssite%pos %pz and so on.
     ! We can make structure constant (C_akL Eq.(38) in /JPSJ.84.034702) here
     ! if we have no extended local orbital.
     ! When we have extentede local obtail, we run m_bstrux_init after elocp in mkpot.
     if( ipr>=30 ) then ! Write atom positions
        write(stdo,"(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x, &
             'pos (multiples of plat)')") 'Basis, after reading restart file'
        do i = 1, nbas
           xvcart = ssite(i)%pos !cartesian coordinate
           xvfrac = matmul(transpose(qlat),xvcart) !fractional coodinate
           write(stdo,"(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)")i,trim(sspec(ssite(i)%spec)%name),xvcart,xvfrac
        enddo
     endif
     !===  loop for electronic structure. Atomic force is calculated
     do 1000 iter = 1,max(1,maxit)
        if(maxit/=0) then
           if (master_mpi) then
              aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
              write(stdo,*)
              write(stdo,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
           endif
           call bndfp(iter,llmfgw,plbnd) !Main part of band calculation. Get total energies ham_ehf and ham_ehk
        endif
        ! Check convergence of dmatu (density matrix for LDA+U) and update it. Get vorbdmat (U potential)
        if(nlibu>0 .AND. lrout>0) call m_ldau_vorbset(ham_ehk,dmatu) !set new vorb from dmat by bndfp-mkpot
        !   Things for --density (plot density) are in locpot.F(rho1mt and rho2mt) and mkpot.F(smooth part).
        !! Write restart file (skip if --quit=band) ---
        if(master_mpi .AND. ( .NOT. cmdopt0('--quit=band'))) then
           if(irs2>0 .AND. (lrout>0 .OR. maxit==0)) then
              !                  lbin = .true. !irs2/=2
              ifi = ifile_handle()
              open(ifi,file='rst.'//trim(sname),form='unformatted') !no rsa now
              k = iors ( iter , 'write' ,irs3,irs5)
              !                  k = iors_old ( iter , 'write' ,irs3,irs5)
              close(ifi)
           endif
        endif
        !     Add to save file; decide on next iteration ---
        if (maxit<=0) goto 9998
        etot(1) = ham_ehf        ! etot(1) is not the total energy when LDA+U
        ! because vefv1,vefv2 do not include U contributions
        etot(2) = ham_ehk + eorb !eorb by LDA+U
        seref   = ham_seref !   ... reference energy
        etot(1) = etot(1) - seref
        if (etot(2)/=0) etot(2) = etot(2) - seref
        flush(6)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(master_mpi) then
           i = 0
           !==   convergence check by call nwit
           !     lsc   :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
           !     :1 if not self-consistent, but encountered max. no. iter.
           !     :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
           !     :3 otherwise
           call nwit( int(ctrl_nvario), iter, maxit, (lhf.or.irs1==0).and.iter==1, &
                leks+i, etol, qtol, qdiff, amom, etot, sev,lsc)
        endif
        call mpibc1_int(lsc,1,'lmfp_lsc')
        if (lsc==2 .AND. ( .NOT. lhf) .AND. maxit>1) lsc = 3
        if (lsc==1 .AND. lrout>0 .OR. lsc==3) then
           if (iter >= maxit) lsc = 1
           if (iter < maxit) lsc = 3
        endif
        if( cmdopt0('--quit=band')) call rx0('lmf-MPIK : exit (--quit=band)')
        if( lsc <= 2) exit  !self-consistency exit
1000 enddo               ! ---------------- SCF (iteration) loop end ----
     if(nitrlx==0) exit     !no molecular dynamics (=no atomic position relaxation)
     !==== Molecular dynamics (relaxiation).
     MDblock: Block !Not maintained well recently but atomic position relaxation was working
       do ibas=1,nbas
          ssite(ibas)%pos0 = ssite(ibas)%pos
          pos_move(:,ibas) = ssite(ibas)%pos
       enddo
       !==   Relax atomic positions.   !--> shear mode removed. probably outside of fortran code if necessary
       call Relax(ssite,sspec,itrlx,indrx_iv,natrlx,force,p_rv,hess,0,[0d0],pos_move,icom)
       !     warn: Updating positions in ssite structure ==> t.kotani think this is confusing because
       !     'positions written in ctrl' and 'positions written in rst' can be different.
       if(master_mpi) write(ifipos) itrlx,pos_move
       do ibas=1,nbas
          ssite(ibas)%pos = pos_move(:,ibas)
       enddo
       !==   Exit when relaxation converged or maximum number of iterations
       if(icom==1) then
          if(master_mpi) then
             flg = 'C67' !what?
             call nwitsv(1+2,ctrl_nvario,flg,nsp,amom,etot,sev)
             write(stdo,"(a,i5)")' LMFP: relaxation converged after iteration(s) of ',itrlx
          endif
          exit
       endif
       !==   Restore minimum gradient positions if this is last step
       if (itrlx==nitrlx) then
          write(stdo,"(a)")' lmfp: restore positions for minimum g'
          call Prelx1(1 , nm , .false. , natrlx , indrx_iv , p_rv, pos_move )
          do ibas=1,nbas !updated positions in site structure
             ssite(ibas)%pos = pos_move(:,ibas)
          enddo
       endif
       !==   New density after atom shifts.
       call Smshft(ctrl_lfrce)
       if (master_mpi) then ! .AND. .NOT. lshr) then
          ifi = ifile_handle()
          open(ifi,file='rst.'//trim(sname),form='unformatted')
          k = iors (iter , 'write' ,irs3,irs5)!Write restart file (to include new positions)
          !           k = iors_old (iter , 'write' ,irs3,irs5)! rst version 1.04
          close(ifi)
          write(stdo,*)' Delete mixing and band weights files ...'
          open(newunit=ifi, file='mixm.'//trim(sname)); close(ifi, status='delete')
          open(newunit=ifi, file='wkp.'//trim(sname)); close(ifi, status='delete')
       endif
       call Parms0(0,0,0d0,0) !   reset mixing block
       if(itrlx==nitrlx .AND. master_mpi) write(stdo,"(a)")' LMFP: relaxation incomplete'
     endblock MDblock
2000 enddo
9998 continue
  if(master_mpi .AND. nitrlx>0) close(ifipos)
  if(allocated(p_rv)) deallocate(p_rv)
  if(allocated(hess)) deallocate(hess)
  call tcx('lmfp')
end subroutine lmfp
