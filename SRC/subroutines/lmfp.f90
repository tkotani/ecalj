module m_lmfp !Driver for iteration loop for lmf-MPIK (electronic and MD)
  public lmfp
  private
contains
  subroutine lmfp(llmfgw)
    use m_lmfinit,only: lhf, maxit,nbas,nsp, ham_seref,  sspec=>v_sspec, ispec, slabl,&
         nlibu,stdo,lrout,leks,plbnd,lpzex, nitrlx, indrx_iv,natrlx,xyzfrz,pdim,qtol,etol,alat
    use m_lattic,only: qlat=>lat_qlat,rv_a_opos
    use m_bandcal,only: dmatu
    use m_mkpot,only:   amom
    use m_ext,only:     sname
    use m_iors,only:    iors
    use m_iors_old,only: iors_old !only for reading vs=1.04 rst file (before 2022-5-14)
    use m_MPItk,only:  mlog,master_mpi
    use m_bndfp,only:  Bndfp,   ham_ehf,ham_ehk,qdiff,force,sev
    use m_ldau,only:   m_ldau_vorbset, eorb
    use m_bstrux,only: m_bstrux_init
    use m_relax,only:  relax
    use m_lattic,only: setopos
    use m_ftox
    use m_rdovfa,only:rdovfa
    !!= Main routine of lmf = (following document is roughly checked at May2021)
    !! lmfp contains two loops after initialization
    !!   1  outer  AtomicPositionRelaxiation:  do 2000 is for molecular dynamics (relaxiation).
    !!   2. innner ElectronicStructureSelfConsistencyLoop: do 1000 
    !!      Main part of band calculaiton is in bndfp.
    !! (Most of) all data in modules are 'protected' (m_density contains iterative quantities)
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
    integer :: i,ifi,ipr, k, nit1,numq, lsc, icom,  nvrelx , itrlx,lscx
    integer:: ibas,unlink,ifipos,iter,j,idmatu,iprint
    real(8) :: gam(4),gam1,bstim,pletot(6,2), plat(3,3),xvcart(3),xvfrac(3),seref,etot(2),vs=2d0,vs1
    real(8),allocatable :: rv_a_omad (:) !  Madelung matrix if necessary
    character(10):: i2char
    logical:: cmdopt2
    character(20):: outs=''
    integer:: ierr,iv
    logical:: irpos,hsign,iatom
    character(256):: strn,strn2
    real(8):: hess(natrlx,natrlx),p_rv(natrlx,10),pos0(3,nbas),poss(3,nbas)
    include "mpif.h"
    call tcn('lmfp')
    ipr = iprint()
    poss = rv_a_opos ! Use atomic positon in m_lattic
    call ReadAtomPos(nbas,poss)! Overwrite pos in the file AtomPos.* if it exists.
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    etot = 0d0 ! Total energy mode --etot ==>moved to m_lmfinit ---
    if(nitrlx>0 ) then ! Atomic position Relaxation setup for 
       icom = 0
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
    if(master_mpi) then ! switch were here [irs1,irs2,irs3,irs5,irs1x10] ==> removed 2022-6-20
       open(newunit=ifi,file='rst.'//trim(sname),form='unformatted') ! Read atomic- and smooth-part density(rhoat smrho in m_density) from atm.* or rst.*
       read(ifi,end=996,err=996) vs !version id of rst file
       close(ifi)
       write(stdo,ftox)'lmv7: Read rst version ID=',ftof(vs,2)
996    continue
    endif
    call Mpi_barrier(MPI_COMM_WORLD,ierr)
    call Mpibc1_real(vs,1,'lmv7: vs: version id of rst file')
    k=-1 !try to read rst files containing density
    if(vs==2d0) k = iors(nit1,'read') ! read rst file. sspec ssite maybe modified
    if(vs/=2d0) k = iors_old(nit1,'read') ! read rst file. sspec ssite maybe modified
    call Mpibc1_int(k,1,'lmv7:lmfp_k')
    if(k<0) then 
       call rdovfa()  ! Initial potential from atm file (lmfa) if rst can not read
       nit1 = 0
       iatom=.true.
    else
       iatom=.false.
    endif
    call Mpi_barrier(MPI_COMM_WORLD,ierr)
    AtomicPositionRelaxiation: do 2000 itrlx = 1,max(1,nitrlx) !loop for atomic position relaxiation(molecular dynamics)
       call setopos( poss )  ! Set position of atoms 
       ! We can make structure constant (C_akL Eq.(38) in /JPSJ.84.034702) here when we have no extended local orbital.
       ! NOTE: When we use extenteded local obtail, we run m_bstrux_init after elocp in mkpot.(we may remove extented local orbial in future)
       if(sum(lpzex)==0) call m_bstrux_init() ! Get structure constants for nbas and qplist
       if( ipr>=30) then ! Write atom positions
          write(stdo,"(/1x,'Basis, after reading restart file'/' site spec',8x,'pos (Cartesian coordinates)',9x,&
               'pos (multiples of plat)')")
          do i = 1, nbas
             xvcart = poss(:,i) !cartesian coordinate
             xvfrac = matmul(transpose(qlat),xvcart) !fractional coodinate
             write(stdo,"(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)")i,trim(slabl(ispec(i))),xvcart,xvfrac
          enddo
       endif
       ElectronicStructureSelfConsistencyLoop: do 1000 iter = 1,maxit
          if(master_mpi) write(stdo,*)
          if(master_mpi) write(stdo,"(a)")trim("--- BNDFP:  begin iteration "//trim(i2char(iter))//" of "//trim(i2char(maxit)))
          call bndfp(iter,llmfgw,plbnd)!Main. Band cal. Get total energies ham_ehf and ham_ehk
          if(nlibu>0.AND.lrout>0)call m_ldau_vorbset(ham_ehk, dmatu)!Set new vorb for dmat by given bndfp (mkpot)
          !    --density(plot density) are in locpot.f90(rho1mt and rho2mt) and mkpot.f90(smooth part).
          if(master_mpi.AND.(.NOT.cmdopt0('--quit=band'))) then
             if(lrout>0) k=iors(iter, 'write')  ! = iors_old ( iter , 'write' ,irs5) ! Write restart file (skip if --quit=band) ---
          endif
          if (maxit<=0) goto 9998 
          etot(1) = ham_ehf        
          etot(2) = ham_ehk + eorb !eorb by LDA+U, etot(1) is not the total energy for LDA+U, because vefv1,vefv2 do not include U contributions
          seref   = ham_seref !   ... reference energy supplied from ctrlfile
          etot(1) = etot(1) - seref
          if (etot(2)/=0) etot(2) = etot(2) - seref
          call mpi_barrier(MPI_COMM_WORLD,ierr)
          if(master_mpi) then
             hsign= (lhf.or.iatom).and.iter==1
             if(itrlx>1) hsign=.false.
             call nwit(0,iter,maxit,hsign,leks,etol,qtol,qdiff,amom,etot,sev, lsc) !convergence check. lsc is the self-consistency flag
             !lsc  :0 self-consistency achieved (diffe<=etol, qdiff=dmxp(11)<=qtol)
             !     :1 if not self-consistent, but encountered max. no. iter.
             !     :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
             !     :3 otherwise
          endif
          call mpibc1_int(lsc,1,'lmfp_lsc') !b
          if (lsc==2 .AND. ( .NOT. lhf) .AND. maxit>1) lsc = 3
          if (lsc==1 .AND. lrout>0 .OR. lsc==3) then
             if (iter >= maxit) lsc = 1
             if (iter < maxit) lsc = 3
          endif
          if(master_mpi) flush(stdo)
          if( cmdopt0('--quit=band')) call rx0('lmf-MPIK : exit (--quit=band)')
          if( lsc <= 2) exit  
1000   enddo ElectronicStructureSelfConsistencyLoop
       if(nitrlx==0) exit !no molecular dynamics (=no atomic position relaxation)
       !==== Molecular dynamics (relaxiation). TK think it is better to comine another approach to move/relax atomic positions.
       MDblock: Block !Not maintained well recently. But Testinstall/te test
         pos0 = poss !keep old poss to pos0
         ! Relax atomic positions. Get new poss. Shear mode removed.
         call Relax(itrlx,indrx_iv,natrlx,force,p_rv,hess,pos0,poss,icom)! input:pos0 to output:poss
         if (itrlx==nitrlx) then !Set minimum gradient positions if this is last step.
            write(stdo,"(a)")'lmfp: restore positions for minimum g (given by relax)'
            do i = 1, natrlx
               poss(indrx_iv(1,i),indrx_iv(2,i)) = p_rv(i,nm) !nm=3. this is given by relax-gradzr
            enddo
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
         !   
         ! New density after atom shifts. (for old user. call smshft(ctrl_lfrce,poss,pos0) was on wrong idea,(New user do not need to know it)).
         !   1. nothing to do means nuleus+core shift.
         !   2. simple superposition of atoms rdovfa
         if (alat*sum((poss-pos0)**2)**.5>0.05d0) then! 0.05 a.u. is intuitively given. 2022-6-22
            write(stdo,*)'lmfp-MDMODE: Reset density by atom superposition'
            iv=iprint()
            call Setopos( poss ) ! Set position of atoms before Rdovfa
            call pshpr(-1) !supress print out to stdo
            call Rdovfa() !superposition of atom density
            call poppr()
         endif
         if (master_mpi) then
            write(stdo,*)' Delete mixing and band weights files ...'
            open(newunit=ifi, file='mixm.'//trim(sname)); close(ifi, status='delete')
            open(newunit=ifi, file='wkp.'//trim(sname)); close(ifi, status='delete')
         endif
         if(icom==1) then ! Exit when relaxation converged or maximum number of iterations
            if(master_mpi) then
               call nwit(-99,iter,maxit,hsign,leks,etol,qtol,qdiff,amom,etot,sev,lscx)
               write(stdo,"(a,i5)")' LMFP: relaxation converged after iteration(s) of ',itrlx
            endif
            exit
         endif
         if(itrlx==nitrlx .AND. master_mpi) write(stdo,"(a)")' LMFP: relaxation incomplete'
       endblock MDblock
2000 enddo AtomicPositionRelaxiation
9998 continue
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
       enddo
    enddo
    close(ifipos)
1010 continue
  end subroutine readatompos
  subroutine nwit(nvario,iter,maxit,lhf,lhk,etol,qtol,qdiff,amom,etot,sev, lsc) ! Add to save file, determine whether to continue execution
    use m_lmfinit,only: nsp
    use m_lgunit,only:stdo,stdl
    use m_ext,only:      sname
    use m_ldau,only:  eorb
    use m_ftox
    implicit none 
    intent(in)::  nvario,iter,maxit,lhf,lhk,etol,qtol,qdiff,amom,etot,sev
    intent(out)::                                                          lsc
    !i   nvario:-99 for MD no variables.
    !i   iter  :current iteration number
    !i   maxit :maximum number of iterations
    !i   lhf   :T if result is generated by overlapped free atoms
    !i         :(only affects key character in save file)
    !i   lhk   :1, if HK energy is available
    !i         :Add 10 to use HK energy for convergence, if available
    !i   etol  :energy tolerance for convergence.  If 0, criterion not used
    !i         :Energy change from prior iteration is compared to etol.
    !i   qtol  :charge tolerance for convergence.  If 0, criterion not used
    !i   qdiff :change in charge
    !i   amom  :magnetic moment
    !i   etot  :etot(1) = HF energy; etot(2) = KS energy
    !o Outputs
    !o   lsc   :0 self-consistency achieved (diffe<=etol, qdiff<=qtol)
    !o         :1 if not self-consistent, but encountered max. no. iter.
    !o         :2 Harris energy from overlap of free atoms (iter=1 and lhf=t)
    !o         :3 otherwise
    logical :: lhf,more,letol,lqtol
    integer :: iter,maxit,lhk,lsc,nvario,lbl
    real(8) :: etol,qtol,qdiff,amom,etot(2),diffe,sev,rydberg
    character flg*1
    real(8),save:: ehf1=0d0,ehk1=0d0
    character(360) :: sout,formatc(903:904),lbll
    call tcn('nwit')
    if(nvario==-99) then
       lsc=4
       goto 1010
    endif
    lsc = 3
    more = .true.
    if(iter >= maxit) then
       lsc = 1
       more = .false.
    endif
    if (lhf) lsc = 2
    formatc(903)="(/  ' it',i4,'  of',i4,'    ehf=',f15.6:'   ehk=',f15.6)"
    formatc(904)="(/'   it',i3,'  of',i3,'    ehf=',f15.6:'   ehk=',f15.6)"
    lbl=904 !assign 904 to lbl
    if (maxit > 999 .OR. iter > 999) lbl=903 !assign 903 to lbl
    lbll=formatc(lbl)
    if (mod(lhk,10) == 1) write(stdo,trim(lbll)) iter,maxit,etot(1),etot(2)
    if (mod(lhk,10) == 0) write(stdo,trim(lbll)) iter,maxit,etot(1)
    if (iter > 1 .OR. etol == 0) then
       diffe = etot(1)-ehf1
       if (lhk == 11) diffe = etot(2)-ehk1
       if (eorb/=0d0) diffe = etot(2)-ehk1 !2023-6-21 for LDAU convergence
       letol = dabs(diffe) .le. etol .or. etol .le. 0
       lqtol = dabs(qdiff) .le. qtol .or. qtol .le. 0
       if (letol .AND. lqtol) then
          more = .false.
          lsc = 0
       endif
       if (iter > 1) then
          if (mod(lhk,10) == 1) write(stdo,905) ehf1,ehk1,diffe,qdiff,etol,qtol,more
          if (mod(lhk,10) == 0) write(stdo,906) ehf1,     diffe,qdiff,etol,qtol,more
905       format(' From last iter',4x,'ehf=',f15.6,'   ehk=',f15.6, &
               /' diffe(q)=',f10.6,' (',f8.6,')','    tol=',f9.6,' (',f8.6,')','   more=',l1)
906       format(' From last iter',4x,'ehf=',f15.6,/' diffe(q)=',f10.6,' (',f8.6,')', &
               '    tol=',f9.6,' (',f8.6,')','   more=',l1)
       else
          write(stdo,"(16x,'rms dq=',f10.6,8x,'tol=',f9.6,'   more=',l1)") qdiff,qtol,more
       endif
    endif
    ehk1 = etot(2)
    ehf1 = etot(1)
1010 continue
    block
      integer:: ifi
      real(8):: rydberg
      character(5):: svflg='cxhiC'
      sout = svflg(lsc+1:lsc+1)
      if (nsp == 2)                 sout=trim(sout)//' mmom='//ftof(amom,4)
      if(eorb==0.and.etot(1)/= 0d0) sout=trim(sout)//' ehf(eV)='//ftof(etot(1)*rydberg(),6)
      if(etot(2)/=0d0)              sout=trim(sout)//' ehk(eV)='//ftof(etot(2)*rydberg(),6)
      if(sev/=0)                    sout=trim(sout)//' sev(eV)='//ftof(sev*rydberg(),6)
      open(newunit=ifi,file='save.'//trim(sname),position='append')
      write(ifi, "(a)") trim(sout)
      close(ifi)
      write(stdo,"(a)") trim(sout)
    end block
    call tcx('nwit')
  end subroutine nwit
endmodule m_lmfp
