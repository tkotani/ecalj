subroutine lmfp(llmfgw)
  use m_lmfinit,only: lhf, maxit,nbas,nsp, ctrl_ldos,ctrl_nitmv,ctrl_nvario, &
       ham_seref,ctrl_lfrce,  sspec=>v_sspec, ssite=>v_ssite, slabl,&
       nlibu,stdo,lrout,leks,plbnd,lpzex, nitrlx, &
       indrx_iv,natrlx,xyzfrz,pdim,qtol,etol,alat
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
  !!   1  outer  MDloop:  do 2000 is for molecular dynamics (relaxiation).
  !!   2. innner Eleloop: do 1000 is for electronic structure self-consistency.
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
  real(8),allocatable :: ftot_rv(:), wk_rv(:), p_rv(:,:),hess(:,:)
  real(8),allocatable :: rv_a_omad (:) !  Madelung matrix if necessary
  character(512):: aaachar
  character(10):: i2char
  logical:: cmdopt2
  character(20):: outs=''
  integer:: ierr,iv
  logical:: irpos,hsign,iatom
  character(256):: strn,strn2
  real(8),allocatable:: poss(:,:),pos0(:,:)
  include "mpif.h"
  
  call tcn('lmfp')
  ipr = iprint()
  allocate(pos0(3,nbas),poss(3,nbas))
  poss = rv_a_opos ! Use atomic positon in m_lattic
  call ReadAtomPos(nbas,poss)! Overwrite pos in AtomPos if it exists. 
     ! Sep2020 " Shorten site positions" here removed.
  etot = 0d0 ! Total energy mode --etot ==>moved to m_lmfinit ---
  if(nitrlx>0 ) then ! Atomic position Relaxation setup (MDloop)
     icom = 0
     if(natrlx /= 0) allocate(hess(natrlx,natrlx),p_rv(natrlx,10))
     if(master_mpi) then
        open(newunit=ifipos,file='AtomPos.'//trim(sname),position='append')
        write(ifipos,ftox) '========'
        write(ifipos,ftox) 0,   '  !itrlx'
        write(ifipos,ftox) nbas,'  !nbas'
        do i=1,nbas
           write(ifipos,ftox) ftof(poss(:,i),16),'   ',i
        enddo
        close(ifipos)
     endif
  endif
  ! switch were here [irs1,irs2,irs3,irs5,irs1x10] ==> removed 2022-6-20

  ! Read atomic- and smooth-part density(rhoat smrho in m_density) from atm.* or rst.*
  vs=2d0
  if(master_mpi) then
     open(newunit=ifi,file='rst.'//trim(sname),form='unformatted')
     read(ifi,end=996,err=996) vs !version id of rst file
     close(ifi)
     write(stdo,ftox)'lmv7: Read rst version ID=',ftof(vs,2)
996  continue
  endif
  call Mpi_barrier(MPI_COMM_WORLD,ierr)
  call Mpibc1_real(vs,1,'lmv7: vs: version id of rst file')
  k=-1 !try to read rst files containing density
  if(vs==2d0) then       !2020-5-14
     k = iors(nit1,'read') ! read rst file. sspec ssite maybe modified
  else                   !vs=1.04d0
     k = iors_old(nit1,'read') ! read rst file. sspec ssite maybe modified
  endif
  call Mpibc1_int(k,1,'lmv7:lmfp_k')
  if(k<0) then 
     call Rdovfa()  ! Initial potential from atm file if rst can not read
     nit1 = 0
     iatom=.true.
  else
     iatom=.false.
  endif
  !smshft is not correct. See T.kotani JPSJ paper for formulation.
  !if(k>=0 .AND. irs1x10) call Smshft(1,poss,poss) ! modify denity after reading rst when irs1x10=True

  !==== Main iteration loops ===
  MDloop: do 2000 itrlx = 1,max(1,nitrlx) !loop for atomic position relaxiation(molecular dynamics)
     call Setopos( poss )  ! Set position of atoms 
     ! We can make structure constant (C_akL Eq.(38) in /JPSJ.84.034702) here
     ! if we have no extended local orbital.
     ! When we have extentede local obtail, we run m_bstrux_init after elocp in mkpot.
     if(sum(lpzex)==0) call M_bstrux_init() ! Get structure constants for nbas and qplist
     if( ipr>=30 ) then ! Write atom positions
        write(stdo,"(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x, &
             'pos (multiples of plat)')") 'Basis, after reading restart file'
        do i = 1, nbas
           xvcart = poss(:,i) !cartesian coordinate
           xvfrac = matmul(transpose(qlat),xvcart) !fractional coodinate
           write(stdo,"(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)")i,trim(slabl(ssite(i)%spec)),xvcart,xvfrac
        enddo
     endif
     Eleloop: do 1000 iter = 1,max(1,maxit)!For electronic structure scf. Atomic forces calculated
        if(maxit/=0) then
           if (master_mpi) then
              aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
              write(stdo,*)
              write(stdo,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
           endif
           call bndfp(iter,llmfgw,plbnd)!Main of band calculation. Get total energies ham_ehf and ham_ehk
        endif
        ! Check conv. of dmatu (density matrix for LDA+U) and update it. Get vorbdmat (U potential)
        if(nlibu>0.AND.lrout>0)call m_ldau_vorbset(ham_ehk,dmatu) !set new vorb from dmat by bndfp-mkpot
        !Things for --density (plot density) are in locpot.F(rho1mt and rho2mt) and mkpot.F(smooth part).
        ! Write restart file (skip if --quit=band) ---
        if(master_mpi .AND. ( .NOT. cmdopt0('--quit=band'))) then
           if(lrout>0 .OR. maxit==0) k=iors(iter, 'write')  ! = iors_old ( iter , 'write' ,irs5)
        endif
        if (maxit<=0) goto 9998  !exit
        etot(1) = ham_ehf        ! etot(1) is not the total energy when LDA+U
        ! because vefv1,vefv2 do not include U contributions
        etot(2) = ham_ehk + eorb !eorb by LDA+U
        seref   = ham_seref !   ... reference energy
        etot(1) = etot(1) - seref
        if (etot(2)/=0) etot(2) = etot(2) - seref
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        if(master_mpi) then
           !==   convergence check by call nwit
           !     lsc   :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
           !     :1 if not self-consistent, but encountered max. no. iter.
           !     :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
           !     :3 otherwise
           hsign= (lhf.or.iatom).and.iter==1
           if(itrlx>1) hsign=.false.
           call nwit(ctrl_nvario,iter,maxit,hsign,leks,etol,qtol,qdiff,amom,etot,sev,lsc)
        endif
        call mpibc1_int(lsc,1,'lmfp_lsc')
        if (lsc==2 .AND. ( .NOT. lhf) .AND. maxit>1) lsc = 3
        if (lsc==1 .AND. lrout>0 .OR. lsc==3) then
           if (iter >= maxit) lsc = 1
           if (iter < maxit) lsc = 3
        endif
        if(master_mpi) flush(stdo)
        if( cmdopt0('--quit=band')) call rx0('lmf-MPIK : exit (--quit=band)')
        if( lsc <= 2) exit  !self-consistency exit
1000 enddo Eleloop              ! ---------------- SCF (iteration) loop end ----
     if(nitrlx==0) exit     !no molecular dynamics (=no atomic position relaxation)
     !==== Molecular dynamics (relaxiation).
     MDblock: Block !Not maintained well recently. Relax and Smshft may need to be corrected.
       pos0 = poss !keep old poss to pos0
       ! Relax atomic positions. Get new poss. Shear mode removed.
       call Relax(ssite,itrlx,indrx_iv,natrlx,force,p_rv,hess,0,[0d0],pos0,poss,icom)
       if (itrlx==nitrlx) then !Set minimum gradient positions if this is last step.
          write(stdo,"(a)")'lmfp: restore positions for minimum g (given by relax)'
          do i = 1, natrlx
             poss(indrx_iv(1,i),indrx_iv(2,i)) = p_rv(i,nm) !nm=3. this is given by relax-gradzr
          enddo !call Prelx1(1 , nm , .false. , natrlx , indrx_iv , p_rv, poss )
       endif
       if(master_mpi) then !new position written to AtomPos.*
          open(newunit=ifipos,file='AtomPos.'//trim(sname),position='append')
          write(ifipos,ftox) '========'
          write(ifipos,ftox) itrlx,'   !itrlx'
          write(ifipos,ftox) nbas, '   !nbas'
          do i=1,nbas
             write(ifipos,ftox) ftof(poss(:,i),16),'   ',i 
          enddo
          close(ifipos)
       endif
       
       ! Smshft is not on wrong idea. See TK paper.
       !   1. nothing to do means nuleus+core shift.
       !   2. simple superposition of atoms rdovfa
       ! call Smshft(ctrl_lfrce,poss,pos0)!New density after atom shifts.
       if (alat*sum((poss-pos0)**2)**.5>0.05d0) then! 0.05 a.u. is intuitively given. 2022-6-22
          write(stdo,*)'RRR: Reset density by atom superposition'
          iv=iprint()
          call Setopos( poss ) ! Set position of atoms before Rdovfa
          call setprint(-1) !supress print out to stdo
          call Rdovfa() !superposition of atom density
          call setprint(iv)
       endif   

       if (master_mpi) then
          write(stdo,*)' Delete mixing and band weights files ...'
          open(newunit=ifi, file='mixm.'//trim(sname)); close(ifi, status='delete')
          open(newunit=ifi, file='wkp.'//trim(sname)); close(ifi, status='delete')
       endif
       if(icom==1) then ! Exit when relaxation converged or maximum number of iterations
          if(master_mpi) then
             call nwitsv(ctrl_nvario,'C',nsp,amom,etot,sev)
             write(stdo,"(a,i5)")' LMFP: relaxation converged after iteration(s) of ',itrlx
          endif
          exit
       endif
       call Parms0(0,0,0d0,0) !   reset mixing block
       if(itrlx==nitrlx .AND. master_mpi) write(stdo,"(a)")' LMFP: relaxation incomplete'
     endblock MDblock
2000 enddo MDloop
9998 continue
  if(allocated(p_rv)) deallocate(p_rv)
  if(allocated(hess)) deallocate(hess)
  call tcx('lmfp')
end subroutine lmfp

subroutine readatompos(nbas,pos)
  use m_ext,only:     sname
  use m_ftox
  real(8):: pos(3,nbas),p(3)
  integer:: ifipos,i,nbas,nbaso
  logical:: irpos
  open(newunit=ifipos,file='AtomPos.'//trim(sname),status='old',err=1010)
  do 
     read(ifipos,*,end=1010)
     read(ifipos,*)
     read(ifipos,*) nbaso
     do i=1,nbaso
        read(ifipos,*) p
        if(i<=nbas) pos(:,i)=p
        !write(6,ftox)i,ftof(p)
     enddo
  enddo
  close(ifipos)
1010 continue
end subroutine readatompos
