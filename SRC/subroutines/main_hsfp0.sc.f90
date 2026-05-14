module m_hsfp0_sc
  ! Step WB.3c: state shared across hsfp0_sc_setup / _consume / _writeout phases.
  ! These module variables are populated by _setup and consumed by _consume and
  ! _writeout. The split exists so a streaming caller (main_hrcxq Phase 1-C)
  ! can interleave its own iq-loop W production with the per-kx sxcf consume,
  ! while still reusing _setup for parameter prep and _writeout for SECU/SEX
  ! file emission.
  integer :: hs_ixc
  logical :: hs_exchange
  real(8) :: hs_ef, hs_esmr, hs_eftrue
  integer :: hs_nspinmx, hs_nq, hs_ngpn1, hs_ngcn1
  real(8), allocatable :: hs_eqx(:,:,:)
  public :: hsfp0_sc, hsfp0_sc_setup, hsfp0_sc_consume, hsfp0_sc_writeout
  private :: Hswriteinit, HsWriteResult
contains
  subroutine hsfp0_sc(skip_init, skip_rx0, ixc_in)
    !> Thin wrapper preserving the original entry-point signature. Drives the
    !> three phases sequentially. Streaming callers (main_hrcxq) call
    !> _setup / _consume-equivalent / _writeout themselves so the iq-loop W
    !> production and the kx-loop sxcf consumption can interleave.
    logical, intent(in), optional :: skip_init, skip_rx0
    integer, intent(in), optional :: ixc_in
    call hsfp0_sc_setup(skip_init=skip_init, ixc_in=ixc_in)
    call hsfp0_sc_consume()
    call hsfp0_sc_writeout(skip_rx0=skip_rx0)
  end subroutine hsfp0_sc

  subroutine hsfp0_sc_setup(skip_init, ixc_in)
    !> Phase 1 of hsfp0_sc: read inputs (ixc, GW data, eigenvalues), determine
    !> ef/esmr/nspinmx/eqx, run sxcf_scz_count, emit XCU/XCD if exchange mode.
    !> Outputs flow through module variables hs_*.
    !> Standalone (do_init=.true.): calls MPI__SplitXq(1, mpi__size) to set up
    !> comm_k/comm_b. Streaming callers (hgw) call SplitXq themselves before setup.
    use m_readqg,only: READQG0,READNGMX2, ngpmx,ngcmx
    use m_READ_BZDATA,only: READ_BZDATA, nqbz,nqibz,n1,n2,n3,ginv,qbz,wbz,qibz
    use m_genallcf_v3,only: GENALLCF_V3,Setesmr, natom,nspin,plat,alat,deltaw,esmr_in=>esmr,nctot,ecore,nband, laf
    use m_itq,only: setitq_hsfp0sc,nbandmx, ntq
    use m_mpi,only: &
         MPI__Initialize,MPI__root,MPI__Broadcast,MPI__rank,MPI__size,MPI__allreducesum, &
         MPI__consoleout, MPI__reduceSum, MPI__SplitXq, comm, ipr
    use m_lgunit,only:m_lgunit_init,stdo
    use m_ftox
    use m_gpu,only: gpu_init
    use m_readfreq_r,only:  Readfreq_r
    use m_hamindex,only:    Readhamindex, symgg=>symops,ngrp
    use m_readgwinput,only: ReadGwinputKeys, ebmx_sig,nbmx_sig
    use m_zmel,only: Mptauof_zmel
    use m_readeigen,only: INIT_READEIGEN,INIT_READEIGEN2,LOWESTEVAL,readeval
    use m_mem,only:writemem,totalram
    use m_sxcf_count,only: sxcf_scz_count
    implicit none
    logical, intent(in), optional :: skip_init
    integer, intent(in), optional :: ixc_in
    logical :: do_init
    integer :: incwfin, ip, is, ix, ierr
    real(8) :: voltot, valn, eftrue, esmref, esmr, ef
    real(8), external :: tripl, rydberg
    real(8) :: qreal(3), wgtq0p, quu(3)
    real(8), allocatable :: eqt(:)
    integer :: ixc, nspinmx
    logical :: legas, exchange, cmdopt2
    character(20):: outs=''
    character(3) :: charnum3
    do_init = .true.
    if(present(skip_init)) do_init = .not. skip_init
    InitOnce: if(do_init) then
      call MPI__Initialize()
      call MPI__SplitXq(1, mpi__size)  ! k-parallel; sets comm_k, comm_b for sxcf
      call gpu_init(comm)
      call M_lgunit_init()
      call writemem('Start hsfp0: TotalRAM per node='//ftof(totalram(),3)//' GB')
      if(MPI__root) then
         if(cmdopt2('--job=',outs)) then
            read(outs,*) ixc
         else
            if(ipr) write(stdo,*) ' --- Choose modes below ------------'
            if(ipr) write(stdo,*) '  Sx(1) Sc(2) ScoreX(3) '
            if(ipr) write(stdo,*) ' --- Put number above ! ------------'
            read(5,*) ixc
            if(ipr) write(stdo,*) ' ixc=', ixc !computational mode index
         endif
      endif
      call MPI__Broadcast(ixc)
    else
      ixc = 2
      if(present(ixc_in)) ixc = ixc_in
    endif InitOnce
    call MPI__consoleout('hsfp0_sc.mode'//charnum3(ixc))
    call pshpr(60)
    if(ixc==3) then; incwfin= -2 !core exchange mode
    else           ; incwfin= -1 !use 7th colmn for core at the end section of GWIN
    endif
    InitGenallcf: if(do_init) then
      call GENALLCF_V3(incwfin)
      call READ_BZDATA()
      call ReadGwinputKeys()
    endif InitGenallcf
    call pshpr(30)
    esmref= esmr_in
    if(ixc==1) then
       esmr = esmr_in
       exchange = .true.
       if(ipr) write(stdo,*) ' --- Exchange mode --- '
    elseif(ixc==2) then
       esmr = esmr_in
       exchange=.false.
       if(ipr) write(stdo,*) ' --- Correlation mode --- '
    elseif(ixc==3) then
       esmr= 0d0
       exchange = .true.
       if(ipr) write(stdo,*) ' --- CORE Exchange mode --- '
    else
       call rx(' hsfp0_sc: Need input (std input) 1(Sx) 2(Sc) or 3(ScoreX)!')
    endif
    call setesmr(esmr_in=esmr)
    InitReadEigen: if(do_init) then
      call Readhamindex()
      call INIT_READEIGEN()
      call INIT_READEIGEN2()
    endif InitReadEigen
    call Mptauof_zmel(symgg,ngrp)
    if(do_init) call Readngmx2()
    if(ipr) write(stdo,"(*(g0))")' max number of G for QGpsi and QGcou: ngcmx ngpmx=',ngcmx,ngpmx
    call pshpr(60)
    if(.NOT.exchange) call readfreq_r()
    legas=.false.
    call efsimplef2ax(legas,esmref, valn,ef)
    eftrue = ef
    if(ixc==3) ef = LOWESTEVAL() -1d-3
    voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
    if(ipr) write(stdo,'(" --- computational conditions --- ")')
    if(ipr) write(stdo,ftox)'  deltaw alat voltot=', ftof([deltaw,alat,voltot])
    if(ipr) write(stdo,ftox)'  ef     esmr   valn=',ftof([ef,esmr,valn])
    nspinmx = nspin
    if(laf) nspinmx=1
    if(mpi__root .AND. mpi__rank/=0) call rx('mpi__root .AND. mpi__rank/=0')
    if(mpi__root) call setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx)
    call MPI_barrier(comm,ierr)
    if(.not.mpi__root) call setitq_hsfp0sc(nbmx_sig,ebmx_sig,eftrue,nspinmx)
    SchedulingSelfEnergyCalculation: block
      if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
      if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
      call sxcf_scz_count(ef,esmr,exchange,ixc,nspinmx)
    endblock SchedulingSelfEnergyCalculation
    WriteoutInit: block
      real(8):: rydberg2
      if(ixc==3.and.ipr) then
         write(stdo,"(a)")'CoreEx mode: We change ef as ef=LOWESTEVAL-1d-3, slightly below the bottom of valence.'
         write(stdo,"(a,f13.5,i5,i5)")' CoreEx mode: ef nspin nctot=',ef,nspin,nctot
         do ix=1,nctot
           write(stdo,"(i4,x,d13.5,x,d13.5)") ix,(ecore(ix,is),is=1,nspin)
         enddo
      endif
      if(allocated(hs_eqx)) deallocate(hs_eqx)
      allocate(hs_eqx(ntq,nqibz,nspin), eqt(nband))
      do is = 1,nspin
         do ip = 1,nqibz
            eqt= READEVAL(qibz(1,ip),is)
            hs_eqx(1:ntq,ip,is) = rydberg()*(eqt(1:ntq) - eftrue)
         enddo
      enddo
      deallocate(eqt)
    endblock WriteoutInit
    ! Stash phase outputs into module state for _consume / _writeout.
    hs_ixc      = ixc
    hs_exchange = exchange
    hs_ef       = ef
    hs_esmr     = esmr
    hs_eftrue   = eftrue
    hs_nspinmx  = nspinmx
    hs_nq       = nqibz
    call Hswriteinit()  ! prints summary, emits XCU/XCD for exchange mode.
  end subroutine hsfp0_sc_setup

  subroutine hsfp0_sc_consume()
    !> Phase 2 of hsfp0_sc: invoke sxcf_scz_correlation/_exchange to fill
    !> zsecall in m_sxcf_sc. Streaming callers replace this phase with
    !> their own per-iq production + per-kx step_kx interleaving.
    use m_sxcf_sc,only: sxcf_scz_correlation, sxcf_scz_exchange
    use m_mpi,only: ipr
    use m_lgunit,only:stdo
    use m_ftox
    if(ipr) write(stdo,ftox) 'gemm version'
    if(hs_exchange)      call sxcf_scz_exchange   (hs_ef, hs_esmr, hs_ixc, hs_nspinmx)
    if(.not.hs_exchange) call sxcf_scz_correlation(hs_ef, hs_esmr, hs_ixc, hs_nspinmx)
  end subroutine hsfp0_sc_consume

  subroutine hsfp0_sc_writeout(skip_rx0)
    !> Phase 3 of hsfp0_sc: reduce zsecall to root, write SECU/SEC2U or
    !> SEXU/SEX2U files, optionally rx0 to exit the program.
    use m_sxcf_sc,only: zsecall,reducez
    use m_mpi,only: MPI__root, ipr
    use m_lgunit,only:stdo
    logical, intent(in), optional :: skip_rx0
    logical :: do_rx0
    integer :: is
    complex(8), allocatable :: zsec(:,:,:)
    do_rx0 = .true.
    if(present(skip_rx0)) do_rx0 = .not. skip_rx0
    call reducez(hs_nspinmx)
    if(MPI__root) then
       do is=1,hs_nspinmx
          allocate(zsec, source= cmplx(zsecall(:,:,:,is),kind=8))
          call HsWriteResult(is, zsec)
          deallocate(zsec)
       enddo
    endif
    call cputid(0)
    if(do_rx0) then
       if(hs_ixc==1) call rx0( ' OK! hsfp0_sc: Exchange mode')
       if(hs_ixc==2) call rx0( ' OK! hsfp0_sc: Correlation mode')
       if(hs_ixc==3) call rx0( ' OK! hsfp0_sc: Core-exchange mode')
    else
       if(MPI__root .and. ipr) write(stdo,'(a,i0)') ' OK! hsfp0_sc returning (no rx0), ixc=',hs_ixc
    endif
    if(.not.do_rx0) return
    stop
  end subroutine hsfp0_sc_writeout

  subroutine Hswriteinit()
    !> Print summary line. For exchange mode, also write XCU/XCD (LDA xc).
    use m_readqg,only: READQG0
    use m_READ_BZDATA,only: nqbz, nqibz, qibz, ginv
    use m_genallcf_v3,only: nspin, alat, deltaw
    use m_rdpp,only: nbloch
    use m_itq,only: ntq
    use m_mpi,only: MPI__root, ipr
    use m_lgunit,only:stdo
    use m_ftox
    integer :: is, ip, i
    integer :: ifxc(2)
    real(8) :: quu(3)
    real(8), allocatable :: vxcfp(:,:,:)
    if(ipr) write(stdo,*)' ***'
    call READQG0('QGpsi',qibz(1:3,1), quu, hs_ngpn1)
    call READQG0('QGcou',qibz(1:3,1), quu, hs_ngcn1)
    if(ipr) write(stdo,ftox)'nspin nq ntq=',nspin, hs_nq, ntq
    if(ipr) write(stdo,ftox)'spin=',is,'nbloch ngp ngc=',nbloch,hs_ngpn1,hs_ngcn1, &
         'nqbz=',nqbz,'nqibz=',nqibz,'ef=',ftof(hs_ef),'Rydberg'
    if(ipr) write(stdo,ftox)'deltaw(Hartree)=',ftof(deltaw),' alat=',ftof(alat), 'esmr=',ftof(hs_esmr)
    PrintLDAexchangecorrelationXCUXCD: if(hs_ixc==1) then
       allocate( vxcfp(ntq,hs_nq,nspin) )
       call rsexx(nspin,qibz,ntq,hs_nq, ginv, vxcfp)
       MPIroot: if(MPI__root) then
          isploop: do is = 1,hs_nspinmx
             if(is==1) open(newunit=ifxc(1),file='XCU')
             if(is==2) open(newunit=ifxc(2),file='XCD')
             write (ifxc(is),*) '==================================='
             write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
             write (ifxc(is),*) '==================================='
             call winfo(ifxc(is),nspin,hs_nq,ntq,is,nbloch,hs_ngpn1,hs_ngcn1,nqbz,nqibz,hs_ef,deltaw,alat,hs_esmr)
             write (ifxc(is),*)' ***'
             write (ifxc(is),"(a)") ' jband   iq ispin                  qibz eigen-Ef (in eV)     LDA XC (in eV)'
             if(ipr) write(stdo,*)
             iploop: do ip = 1,hs_nq
                do i  = 1,ntq
                   write(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                        i,ip,is, qibz(1:3,ip), hs_eqx(i,ip,is), vxcfp(i,ip,is)
                   if(hs_eqx(i,ip,is) <1d20 .AND. vxcfp(i,ip,is)/=0d0) then
                      if(ipr) write(stdo,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f10.4,'  Sxc(LDA)=',f10.4)") &
                           i,ip,is, qibz(1:3,ip), hs_eqx(i,ip,is), vxcfp(i,ip,is)
                   endif
                enddo
             enddo iploop
             close(ifxc(is))
          enddo isploop
       endif MPIroot
       deallocate(vxcfp)
    endif PrintLDAexchangecorrelationXCUXCD
  end subroutine Hswriteinit

  subroutine HsWriteResult(is, zsec)
    !> Write SE{X,C}{U,D} text file and SE{X,C}2{U,D} binary file for spin `is`.
    use m_READ_BZDATA,only: nqbz, nqibz, qibz, n1, n2, n3
    use m_genallcf_v3,only: nspin, alat, deltaw
    use m_rdpp,only: nbloch
    use m_itq,only: ntq
    use m_mpi,only: ipr
    use m_lgunit,only:stdo
    integer, intent(in) :: is
    complex(8), intent(in) :: zsec(:,:,:)
    integer:: ip, i, ifsec(2), ifsex(2), ifsex2(2), ifsec2(2)
    real(8):: hartree, rydberg
    character(1):: keys
    character(4):: kcore
    hartree = 2d0*rydberg()
    if(is==1) keys='U'
    if(is==2) keys='D'
    kcore=''
    if(hs_ixc==3) kcore='core'
    if(hs_exchange) then
       open(newunit=ifsex(is), file='SEX'//trim(kcore)//keys)
       open(newunit=ifsex2(is),file='SEX'//trim(kcore)//'2'//keys,form='unformatted')
       write(ifsex(is),*) '======================================='
       write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
       write(ifsex(is),*) '======================================='
       call winfo(ifsex(is),nspin,hs_nq,ntq,is,nbloch,hs_ngpn1,hs_ngcn1,nqbz,nqibz,hs_ef,deltaw,alat,hs_esmr)
       write (ifsex(is),*)' *** '
       write (ifsex(is),"(a)")&
       ' jband   iq ispin                             qibz            eigen-Ef (in eV)           exchange (in eV)'
       write(ifsex2(is)) nspin, hs_nq, ntq, nqbz, nqibz, n1, n2, n3
       if(ipr) write(stdo,*)
       do ip = 1,hs_nq
          do i  = 1,ntq
             write(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                  i,ip,is, qibz(1:3,ip), hs_eqx(i,ip,is), hartree*dreal(zsec(i,i,ip))
             if( hs_eqx(i,ip,is)<1d20 .AND. abs(zsec(i,i,ip))/=0d0 ) then
                if(ipr) write(stdo,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,' eig=',f10.4,'  Sx=',f10.4)") &
                     i,ip,is, qibz(1:3,ip), hs_eqx(i,ip,is), hartree*dreal(zsec(i,i,ip))
             endif
          enddo
          write(ifsex2(is)) is, qibz(1:3,ip), zsec(1:ntq,1:ntq,ip)
       enddo
       close(ifsex(is))
       close(ifsex2(is))
    elseif(hs_ixc==2) then
       open(newunit=ifsec(is),file='SEC'//keys)
       open(newunit=ifsec2(is),file='SEC2'//keys,form='unformatted')
       write(ifsec2(is)) nspin, hs_nq, ntq, nqbz, nqibz, n1, n2, n3
       write(ifsec(is),*) '=========================================='
       write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)") is
       write(ifsec(is),*) '=========================================='
       call winfo(ifsec(is),nspin,hs_nq,ntq,is,nbloch,hs_ngpn1,hs_ngcn1,nqbz,nqibz,hs_ef,deltaw,alat,hs_esmr)
       write (ifsec(is),*)' *** '
       write (ifsec(is),"(a)") ' jband   iq ispin                  '// &
            '           qibz            eigen-Ef (in eV)           '// &
            'Re(Sc) 3-points (in eV)                        '// &
            '           In(Sc) 3-points (in eV)                Zfactor(=1)'
       do ip = 1,hs_nq
          do i  = 1,ntq
             if( hs_eqx(i,ip,is)<1d20 .AND. abs(zsec(i,i,ip))/=0d0 ) then
                if(ipr) write(stdo,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f8.4,'  Re(Sc) =',f8.4,'  Img(Sc) =',f8.4 )") &
                     i,ip,is, qibz(1:3,ip), hs_eqx(i,ip,is), hartree*dreal(zsec(i,i,ip)), hartree*dimag(zsec(i,i,ip))
             endif
             write(ifsec(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16, 3x,d24.16)") &
                  i,ip,is, qibz(1:3,ip), hs_eqx(i,ip,is), hartree*dreal(zsec(i,i,ip)), hartree*dimag(zsec(i,i,ip))
          end do
          write(ifsec2(is)) is, qibz(1:3,ip), zsec(1:ntq,1:ntq,ip)
       end do
       close(ifsec(is))
       close(ifsec2(is))
    endif
  end subroutine HsWriteResult

end module m_hsfp0_sc

subroutine rsexx (nspin, q, ntq,nq,ginv, vxco)
  use m_lgunit,only:m_lgunit_init,stdo
  use m_mpi,only:ipr
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  dimension vxco(ntq,nq,nspin),q(3,nq)!,itq(ntq) !itq is not dependent on q, right?
  real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
  logical ::nocore,lfind
  real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3)
  integer:: ikpx=999999
  if(ipr) write(stdo,*)' OPEN VXCFP '
  open(newunit=ifvxcfp,file='__VXCFP',form='unformatted')
  read(ifvxcfp) ldim,nqbz
  if(ipr) write(stdo,*)' rsexx ldim,nqbz',ldim,nqbz
  allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
  do ikp = 1,nqbz
     read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
     if(ipr) write(stdo,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
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
