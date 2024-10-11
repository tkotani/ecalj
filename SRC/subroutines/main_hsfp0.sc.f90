subroutine hsfp0_sc()
  !> Calculates the self-energy \Sigma in GW approximation,  checked 2020jul
  !!  including Off-diagonal components.
  !!  (hsfp0.F is for diagonal part only).
  !! ----------------------------------------
  !!    SEx(q,itp,itpp) = <psi(q,itp) |SEx| psi(q,itpp)>
  !!    SEc(q,itp,itpp) = <psi(q,itp) |SEc| psi(q,itpp)>
  !!    Here SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
  !!
  !! ----------------------------------------
  !! See papers;
  !! [1]T. Kotani and M. van Schilfgaarde, Quasiparticle self-consistent GW method:
  !!     A basis for the independent-particle approximation, Phys. Rev. B, vol. 76, no. 16,
  !!     p. 165106[24pages], Oct. 2007.
  !! [2]T. Kotani, Quasiparticle Self-Consistent GW Method Based on the Augmented Plane-Wave
  !!     and Muffin-Tin Orbital Method, J. Phys. Soc. Jpn., vol. 83, no. 9, p. 094711 [11 Pages], Sep. 2014.
  !!
  !! EIBZ symmetrization;
  !! See [3] C. Friedrich, S. Bl?gel, and A. Schindlmayr,
  !!   Efficient implementation of the GW approximation within the all-electron FLAPW method,
  !!   Physical Review B, vol. 81, no. 12, Mar. 2010.
  !!
  !! Usage: This routine is called from a script for QSGW, ecalj/fpgw/exec/gwsc,
  !! which calls as "echo 2|../exec/hsfp0_sc >lsc" when mode=2.
  !! In the gwsc script, we call hsfp0_sc three times with mode=1,2,3.
  !!   mode= 1: exchange    mode SEx, the exchange part of the self-energy
  !!   mode= 2: correlation mode SEc, the correlated part of the self-energy
  !!   mode= 3: core exchange mode SEXcore
  !! (unused now) mode= 4: plot spectrum function ---See manual ---> this is performed by echo 4|hsfp0
  !!
  !! iSigMode parameter, which determines approximation for self-energy, is given by GWinput file as iSigMode.
  !!     iSigMode==0 SE_nn'(ef)+image integr:delta_nn'(SE_nn(e_n)-SE_nn(ef))
  !!     iSigMode==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
  !!       xxx not support this mode now ... iSigMode==2 SE_nn'((e_n+e_n')/2)
  !! --> iSigMode==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2  <--- this is used as default.
  !!     iSigMode==5 delta_nn' SE_nn(e_n)
  !!     Output file contain hermitean part of SE for energies to be real
  !!    (for example, hermitean conjunction of SE_nn'(e_n) means SE_n'n(e_n')^* )
  !!
  !!     History: We learned so much from GW-LMTO-ASA codeds developed by F.Aryasetiawan.
  !!
  !! Module prograing style
  !!   For writing codes, we recommend our rule of module programing. Minimum are
  !!   1)Classification of variables are in modules
  !!      We define variables in modules. For example, we use natom (number of atoms in the cell)
  !!      as well as pos(3,natom) in the module m_genallcf_v3.
  !!      Anywhere in the prograom, we have to read natom and pos contained in m_genallcf_v3.
  !!      Not make copies of them. Not declear pos(3,natom) except very simple subroutines.
  !!   2) Module variables are 'protected,public' or 'private'. Never use unprotected public variables.
  !!      With 'private' as defaults, it is better to show subrouitnes as public. (default as private).
  !!   3) subroutines which are not in the modules should be very simple mathematical ones.
  !!   4) The intent should be written like this
  !!       subroutine foobar (nctot,ncc,nt0,ntp0,  iclass, phase,
  !!       implicit none
  !!      i                   icore,ncore,nl,nnc,
  !!      o                   zpsi2b)
  !!       intent(in)::       nctot,ncc,nt0,ntp0,  iclass, phase,
  !!      i                   icore,ncore,nl,nnc
  !!       intent(out)::      zpsi2b
  !!           ...
!!!! ----------------------------------------
  use m_readqg,only: READQG0,READNGMX2, ngpmx,ngcmx
  use m_READ_BZDATA,only: READ_BZDATA, nqbz,nqibz,n1,n2,n3,ginv,qbz,wbz,qibz 
  use m_genallcf_v3,only: GENALLCF_V3,Setesmr, natom,nspin,plat,alat,deltaw,esmr_in=>esmr,nctot,ecore,nband, laf
  use m_itq,only: setitq_hsfp0sc,nbandmx, ntq
!  use m_readhbe,only: Readhbe, nband !nprecb,mrecb,mrece,nlmtot,nqbzt,,mrecg
  use m_mpi,only: &
       MPI__Initialize,MPI__root,MPI__Broadcast,MPI__rank,MPI__size,MPI__allreducesum, &
       MPI__consoleout,  MPI__reduceSum,comm
  use m_lgunit,only:m_lgunit_init,stdo
  use m_ftox
  use m_gpu,only: gpu_init
  implicit none
  ! real(8),parameter :: ua  = 1d0 ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
  !
  !\Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
  !  --->   I found COH term method show poor accuracy.
  integer::  ixc, ip, is, nspinmx, i, ix, nq0ix, ngpn1,ngcn1, timevalues(8), irank,isp,nq,ierr
  real(8) :: voltot,valn,efnew,hartree,qreal(3),wgtq0p,quu(3), eftrue,esmref,esmr,ef
  character(128) :: ixcname
  logical:: legas, exonly, iprintx,diagonly=.false.,exchange, hermitianW=.true.
!  integer,allocatable:: irkip(:,:,:,:), nrkip(:,:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), eqt(:), eq(:), eqx(:,:,:),eqx0(:,:,:)
  !complex(8),pointer::zsec(:,:,:)
  complex(8),allocatable:: zsec(:,:,:)
  InitializationBlock:block
    use m_readfreq_r,only:  Readfreq_r
    use m_hamindex,only:    Readhamindex, symgg=>symops,ngrp
    use m_readgwinput,only: ReadGwinputKeys, ebmx_sig,nbmx_sig
    use m_zmel,only: Mptauof_zmel
    use m_readeigen,only: INIT_READEIGEN,INIT_READEIGEN2,LOWESTEVAL
    use m_mem,only:writemem,totalram
    integer:: incwfin
    real(8):: tripl
    logical:: cmdopt2
    character(20):: outs=''
    call MPI__Initialize()
    call gpu_init() 
    call M_lgunit_init()
    call writemem('Start hsfp0: TotalRAM per node='//ftof(totalram(),3)//' GB')
    if(MPI__root) then
       if(cmdopt2('--job=',outs)) then
          read(outs,*) ixc
       else
          write(stdo,*) ' --- Choose modes below ------------'
          write(stdo,*) '  Sx(1) Sc(2) ScoreX(3) '
          write(stdo,*) ' --- Put number above ! ------------'
          read(5,*) ixc
          write(stdo,*) ' ixc=', ixc !computational mode index
       endif
    endif
    call MPI__Broadcast(ixc)
    write(ixcname,"('.mode=',i4.4)")ixc
    call MPI__consoleout('hsfp0_sc'//trim(ixcname)) !Open console output stdout.irank.hsfp0_sc.mode
    call pshpr(60)
    if(ixc==3) then; incwfin= -2 !core exchange mode
    else           ; incwfin= -1 !use 7th colmn for core at the end section of GWIN
    endif
    call GENALLCF_V3(incwfin)    ! readin basic data
    call READ_BZDATA() ! Readin BZ data. See gwsrc/rwbzdata.f ===
!    if(nclass /= natom ) call rx('hsfp0_sc: sanitiy check nclass /= natom ')! CAUTION. ASSUME iclass(iatom)= iatom (because of historical reason)
    !  write(stdo,"(' nqbz nqibz ngrp=',3i12)") nqbz,nqibz,ngrp
    call ReadGwinputKeys()
    call pshpr(30)
    esmref= esmr_in
    if(ixc==1) then
       esmr = esmr_in !read from GWinput
       exchange = .true.
       write(stdo,*) ' --- Exchange mode --- '
    elseif(ixc==2) then
       esmr = esmr_in !read from GWinput
       exchange=.false.
       write(stdo,*) ' --- Correlation mode --- '
    elseif(ixc==3) then
       esmr= 0d0
       exchange = .true.
       write(stdo,*) ' --- CORE Exchange mode --- '
    else
       call rx(' hsfp0_sc: Need input (std input) 1(Sx) 2(Sc) or 3(ScoreX)!')
    endif
    call setesmr(esmr_in=esmr) !set esmr back in genalloc_v3
!    call Readhbe()        ! Read dimensions in m_readhbe
    call Readhamindex()
    call INIT_READEIGEN() ! initialization for readeigen readcphi readgeig.
    call INIT_READEIGEN2()! initialize m_readeigen
    call Mptauof_zmel(symgg,ngrp) ! Put space-group transformation information to m_zmel,call rdpp in this
    call Readngmx2() !return ngpmx and ngcmx in m_readqg
    write(stdo,"(*(g0))")' max number of G for QGpsi and QGcou: ngcmx ngpmx=',ngcmx,ngpmx
    call pshpr(60)
    if(.NOT.exchange) call readfreq_r()  !Readin WV.d and freq_r
    legas=.false. ! if legas=T, homogenius electron gas test case.
    call efsimplef2ax(legas,esmref, valn,ef)!Get num of val electron valn and Fermi energy ef. legas=T give ef for given valn.
    eftrue = ef
    if(ixc==3)ef = LOWESTEVAL() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
    voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3))) ! primitive cell volume
    write(stdo,'(" --- computational conditions --- ")')
    write(stdo,ftox)'  deltaw alat voltot=', ftof([deltaw,alat,voltot])
    write(stdo,ftox)'  ef     esmr   valn=',ftof([ef,esmr,valn])
    nspinmx = nspin
!    call anfcond()
    if(laf) nspinmx=1  ! Antiferro case. Only calculate up spin
    if(mpi__root .AND. mpi__rank/=0) call rx('mpi__root .AND. mpi__rank/=0')
    if(mpi__root) call Setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx) !read or set NTQXX and nbandmx
    call MPI_barrier(comm,ierr) !barrier for writing NTQXX at irank=0
    if(.not.mpi__root) call Setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx) !read NTQXX and nbandmx
  endblock InitializationBlock
  SchedulingSelfEnergyCalculation: Block
    use m_sxcf_count,only: sxcf_scz_count
    if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
    if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
    call sxcf_scz_count(ef,esmr,exchange,ixc,nspinmx) !initializatition quick
  endblock SchedulingSelfEnergyCalculation
  WriteoutInit: block
    use m_readeigen,only: readeval
    real(8):: rydberg
    if(ixc==3) then ! Core-exchange mode. We set Ef just below the valence eigenvalue (to pick up only cores)
       write(stdo,"(a)")'CoreEx mode: We change ef as ef=LOWESTEVAL-1d-3, slightly below the bottom of valence.'
       write(stdo,"(a,f13.5,i5,i5)")' CoreEx mode: ef nspin nctot=',ef,nspin,nctot
       do ix=1,nctot
          write(stdo,"(i4,x,d13.5,x,d13.5)") ix,(ecore(ix,is),is=1,nspin)
       enddo
    endif
    allocate(eqx(ntq,nqibz,nspin),eqt(nband))
    nq  =nqibz
    do is = 1,nspin
       do ip = 1,nqibz
          eqt= READEVAL(qibz(1,ip),is) ! Read eigenvalues given by lmf for writing out files.
          eqx (1:ntq,ip,is) = rydberg()*(eqt(1:ntq)- eftrue)
       enddo
    enddo
    deallocate(eqt)
    hartree=2d0*rydberg()
    call Hswriteinit()    !internal subroutine containd below.  Write initial part of output files
  endblock WriteoutInit
  !! Currently, irkip, involved in eibz (extended BZ) mode is removed.
  !! == irkip:  parallelization is controled by irkip ==
  !! We have to distribute non-zero irkip to all ranks (irkip is dependent on rank).
  !! When irkip(nqibz,ngrp,nq,nspinmx)/=0, we expect grain-size
  !! on each job of (iqibz,igrp,iq,isp) is almost the same.
  !! Our pupose is to calculate the self-energy zsec(itp,itpp,iq).
  ! Remove eibzmode symmetrizer 2023Jan22
  !call Seteibzhs(nspinmx,nq,qibz,iprintx=MPI__root)
  Main4SelfEnergy: Block !time-consuming part Need highly paralellized
    ! use m_sxcf_main,only: sxcf_scz_correlation,sxcf_scz_exchange
    ! use m_sxcf_gemm,only: sxcf_scz_correlation_gemm => sxcf_scz_correlation, sxcf_scz_exchange_gemm => sxcf_scz_exchange
    ! logical:: use_original, cmdopt0
    !    use_original = cmdopt0('--oldsxcf') !2024-7-23
    !    if(use_original) then  
    !write(stdo,ftox)'original old version'
    !if(exchange)      call sxcf_scz_exchange   (ef,esmr,ixc,nspinmx) !main part of job
    !if(.not.exchange) call sxcf_scz_correlation(ef,esmr,ixc,nspinmx) !main part of job
    !    else
    !      write(stdo,ftox) 'gemm version'
    !      if(exchange)      call sxcf_scz_exchange_gemm   (ef,esmr,ixc,nspinmx) !main part of job
    !      if(.not.exchange) call sxcf_scz_correlation_gemm(ef,esmr,ixc,nspinmx) !main part of job
    !    endif
    use m_sxcf_gemm,only: sxcf_scz_correlation, sxcf_scz_exchange
    write(stdo,ftox) 'gemm version'
    if(exchange)      call sxcf_scz_exchange   (ef,esmr,ixc,nspinmx) !main part of job
    if(.not.exchange) call sxcf_scz_correlation(ef,esmr,ixc,nspinmx) !main part of job
  EndBlock Main4SelfEnergy
! Remove eibzmode symmetrizer 2023Jan22 (extended irreducibel BZ mode)
!  SymmetrizeZsec :Block
!    logical:: eibz4sig
!    if(eibz4sig())then
!       do isp=1,nspinmx
!          call zsecsym(zsecall(:,:,:,isp),ntq,nq,nband,nbandmx,nspinmx,nspin,eibzsym,ngrp,tiii,qibz,isp)
!       enddo
!    endif
!  EndBlock SymmetrizeZsec
  Finalizesum: block
!    use m_sxcf_main,only: zsecall
    use m_sxcf_gemm,only: zsecall,reducez
    call reducez(nspinmx)
    if(MPI__root) then
       do is=1,nspinmx
          allocate(zsec,source= cmplx(zsecall(:,:,:,is),kind=8))
          call HsWriteResult() !internal subroutine. write only
          deallocate(zsec)
       enddo
    endif
    call cputid(0)
    if(ixc==1 ) call rx0( ' OK! hsfp0_sc: Exchange mode')
    if(ixc==2 ) call rx0( ' OK! hsfp0_sc: Correlation mode')
    if(ixc==3 ) call rx0( ' OK! hsfp0_sc: Core-exchange mode')
  endblock Finalizesum
  stop
contains 
  subroutine Hswriteinit() !contained in hsfp0_sc. Only write out files, no side effect
    use m_rdpp,only: nbloch !Rdpp ! Generate matrix element for "call get_zmelt".
    implicit none
    integer::ifxc(2)
    write(stdo,*)' ***'
    call READQG0('QGpsi',qibz(1:3,1), quu,ngpn1)
    call READQG0('QGcou',qibz(1:3,1), quu,ngcn1)
    write(stdo,ftox)'nspin nq ntq=',nspin,nq,ntq
    write(stdo,ftox)'spin=',is,'nbloch ngp ngc=',nbloch,ngpn1,ngcn1,'nqbz=',nqbz,'nqibz=',nqibz,'ef=',ftof(ef),'Rydberg'
    write(stdo,ftox)'deltaw(Hartree)=',ftof(deltaw),' alat=',ftof(alat), 'esmr=',ftof(esmr)
    PrintLDAexchangecorrelationXCUXCD: if(ixc==1) then
       allocate(  vxcfp(ntq,nq,nspin) )
       call rsexx(nspin,qibz,ntq,nq, ginv, vxcfp) !add ginv july2011
       MPIroot: if(MPI__root) then
          isploop: do is = 1,nspinmx
             if(is==1) open(newunit=ifxc(1),file='XCU')!//xt(nz))
             if(is==2) open(newunit=ifxc(2),file='XCD')!//xt(nz))
             write (ifxc(is),*) '==================================='
             write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
             write (ifxc(is),*) '==================================='
             call winfo(ifxc(is),nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
             write (ifxc(is),*)' ***'
             write (ifxc(is),"(a)") ' jband   iq ispin                  qibz eigen-Ef (in eV)     LDA XC (in eV)'
             write(stdo,*)
             iploop: do ip = 1,nq
                do i  = 1,ntq
                   write(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                        i,ip,is, qibz(1:3,ip), eqx(i,ip,is), vxcfp(i,ip,is)
                   if(eqx(i,ip,is) <1d20 .AND. vxcfp(i,ip,is)/=0d0) then
                      write(stdo,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f10.4,'  Sxc(LDA)=',f10.4)") &
                           i,ip,is, qibz(1:3,ip), eqx(i,ip,is), vxcfp(i,ip,is)
                   endif
                enddo
             enddo iploop
             close(ifxc(is))
          enddo isploop
       endif MPIroot
       deallocate(vxcfp)
    endif PrintLDAexchangecorrelationXCUXCD
  end subroutine Hswriteinit
  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine HsWriteResult() !contained in hsfp0_sc. Only write out files, no side effect
    use m_rdpp,only: nbloch !Rdpp ! Generate matrix element for "call get_zmelt".
    implicit none
    integer:: iii,ixx,iyy, ifsec(2), ifsex(2), ifsex2(2),ifsec2(2)
    character(1):: keys
    character(4):: kcore
    if(is==1) keys='U'
    if(is==2) keys='D'
    kcore=''
    if(ixc==3) kcore='core'
    if(exchange) then
       open(newunit=ifsex(is), file='SEX'//trim(kcore)//keys) !//xt(nz))
       open(newunit=ifsex2(is),file='SEX'//trim(kcore)//'2'//keys,form='unformatted') !out SEX_nn'
       write(ifsex(is),*) '======================================='
       write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
       write(ifsex(is),*) '======================================='
       call winfo(ifsex(is),nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
       write (ifsex(is),*)' *** '
       write (ifsex(is),"(a)")&
       ' jband   iq ispin                             qibz            eigen-Ef (in eV)           exchange (in eV)'
       write(ifsex2(is)) nspin, nq, ntq,nqbz,nqibz, n1,n2,n3
       write(stdo,*)
       do ip = 1,nq
          do i  = 1,ntq
             write(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                  i,ip,is, qibz(1:3,ip), eqx(i,ip,is), hartree*dreal(zsec(i,i,ip)) 
             if( eqx(i,ip,is)<1d20 .AND. abs(zsec(i,i,ip))/=0d0 ) then 
                write(stdo,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,' eig=',f10.4,'  Sx=',f10.4)") &
                     i,ip,is, qibz(1:3,ip), eqx(i,ip,is), hartree*dreal(zsec(i,i,ip))
             endif
          enddo
          write(ifsex2(is)) is, qibz(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
       enddo
       close(ifsex(is))
       close(ifsex2(is))
    elseif(ixc==2) then
       open(newunit=ifsec(is),file='SEC'//keys) !//xt(nz))
       open(newunit=ifsec2(is),file='SEC2'//keys,form='unformatted') !out SEC_nn'
       write(ifsec2(is)) nspin, nq, ntq ,nqbz,nqibz  ,n1,n2,n3
       write(ifsec(is),*) '=========================================='
       write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)") is
       write(ifsec(is),*) '=========================================='
       call winfo(ifsec(is),nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
       write (ifsec(is),*)' *** '
       write (ifsec(is),"(a)") ' jband   iq ispin                  '// &
            '           qibz            eigen-Ef (in eV)           '// &
            'Re(Sc) 3-points (in eV)                        '// &
            '           In(Sc) 3-points (in eV)                Zfactor(=1)'
       do ip = 1,nq
          do i  = 1,ntq
             if( eqx(i,ip,is)<1d20 .AND. abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                write(stdo,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f8.4,'  Re(Sc) =',f8.4,'  Img(Sc) =',f8.4 )") &
                     i,ip,is, qibz(1:3,ip), eqx(i,ip,is),hartree*dreal(zsec(i,i,ip)),hartree*dimag(zsec(i,i,ip))
             endif
             write(ifsec(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16, 3x,d24.16)") &
                  i,ip,is, qibz(1:3,ip), eqx(i,ip,is),hartree*dreal(zsec(i,i,ip)),hartree*dimag(zsec(i,i,ip))
          end do
          write(ifsec2(is)) is, qibz(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
       end do
       close(ifsec(is))
       close(ifsec2(is))
    endif      
  end subroutine HsWriteResult
end subroutine hsfp0_sc

subroutine rsexx (nspin, q, ntq,nq,ginv, vxco)
  use m_lgunit,only:m_lgunit_init,stdo
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  dimension vxco(ntq,nq,nspin),q(3,nq)!,itq(ntq) !itq is not dependent on q, right?
  real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
  logical ::nocore,lfind
  real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3)
  integer:: ikpx=999999
  write(stdo,*)' OPEN VXCFP '
  open(newunit=ifvxcfp,file='VXCFP',form='unformatted')
  read(ifvxcfp) ldim,nqbz  
  write(stdo,*)' rsexx ldim,nqbz',ldim,nqbz
  allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
  do ikp = 1,nqbz
     read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
     write(stdo,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
  enddo
  close(ifvxcfp)
  do iq=1,nq
     do ikp=1,nqbz
        lfind=.false.
        if(sum( (qqq(1:3,ikp)-q(1:3,iq))**2) <tolq) then
           lfind=.true.
        else
           call rangedq( matmul(ginv,q(1:3,iq)-qqq(:,ikp)), qx)
           if(sum(abs(qx))< tolq) lfind= .TRUE. 
        endif
        if(lfind) then
           ikpx=ikp
           goto 100
        endif
     enddo
     call rx( ' rsexx: not find ikp')
100  continue
     vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(1:ntq,ikpx,1:nspin)
  enddo
end subroutine rsexx
