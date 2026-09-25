!> Sigma in the MLO representation: read SigRsMLO, Bloch-sum it at any k, and
!> build the PMT matrix element (design: Samples/kBT/sigma_mlo_design.md, stage 4).
!>
!>   Sigma^MLO(k)   = sum_R Sigma^MLO(R) exp(-ikR)          eq (12)
!>   z^MLO(k)       = Hreduction(H^LDA(k), S^PMT(k))        eq (17), rebuilt, never interpolated
!>   A(k)           = S^PMT(k) z^MLO(k),  O^MLO = z^dag A   eq (15)
!>   [Sigma]_mn(k)  = A (O^MLO)^-1 Sigma^MLO(k) (O^MLO)^-1 A^dag   eq (14)
!>
!> SigRsMLO is self-contained (it carries the pair list), so nothing here depends
!> on m_HamPMT's module state.
module m_sigmlo
  use m_lgunit,only: stdo
  use m_ftox
  implicit none
  public :: sigmlo_init, sigmlo_senex, sigmlo_on, read_mloindex, zmlo_frozen, zmlo_store
  public :: zmlo_new_append
  logical, protected :: sigmlo_on = .false.
  private
  integer :: ndimMTO=0, npairmx=0, nspx=0, nbas=0, mlomethod=4, nskip=0
  integer, allocatable :: npair(:,:), nlat(:,:,:,:), nqwgt(:,:,:), ib_tableM(:), ix(:)
  real(8) :: plat(3,3), fff1=2d0, eferm=0d0, ecbot=0d0
  complex(8), allocatable :: sigmlor(:,:,:,:)   !(npairmx, ndimMTO, ndimMTO, nspx)
  logical :: init = .true.
  !--- chi~: TWO SLOTS.  Sigma^MLO_{ab} = <chi~_a|Sigma|chi~_b> is a matrix IN a basis,
  !    so putting it back into the PMT basis is only the same operator if the SAME chi~
  !    is used.  Meanwhile chi~ itself must follow the Hamiltonian, or its energy window
  !    goes stale as QSGW moves the bands.  Both hold at once with two slots:
  !
  !      ZmloSig  the z^MLO the CURRENT SigRsMLO was written in.  getsenex reads Sigma
  !               back through this one.  Loaded here, never appended to.
  !      ZmloNew  the z^MLO sugw's step a' builds from the H of THIS iteration and uses
  !               for c^MLO, so the next Sigma^MLO is written in it.  gwsc promotes it
  !               to ZmloSig right after `mlo` has written the new SigRsMLO.
  !
  !    At a k that is not in ZmloSig (a band plot, or an SCF mesh finer than the Sigma
  !    mesh) there is nothing to read back exactly -- the H of that iteration is gone --
  !    so z is rebuilt from the H in hand and kept in memory for the rest of this run.
  !    ndimh depends on q when pwmode=11, so an entry may only be reused at its own size;
  !    taking the first rows of a zero-padded entry gave a singular O^MLO and NaN.
  !    ECALJ_MLO_NOCACHE=1 disables both slots and rebuilds at every call (the 2026-09-25
  !    failure mode: chi~ then moves inside one SCF loop).
  integer :: nzc = 0
  real(8), allocatable :: qzc(:,:)
  integer, allocatable :: ispzc(:), ndzc(:)   !ndzc: the ndimh each entry was built with.
  complex(8), allocatable :: zcache(:,:,:)
  logical :: nocache = .false., cfirst = .true.
  integer :: nzsig = 0                        !how many entries came from ZmloSig
contains
  !> The MLO index lives in the trailing records of HamRsMLO (not in a __-prefixed
  !> file, which cleargw would delete).  Skip the four data records, then read it.
  subroutine read_mloindex(nd, ld, mm, nsk, ixo, f1, ef, ec)
    use m_cmdopt_registry, only: c0_socmatrix
    integer,intent(out):: nd, ld, mm, nsk
    integer,allocatable,intent(out):: ixo(:)
    real(8),intent(out):: f1, ef, ec
    integer:: ifh, n1,n2,n3
    open(newunit=ifh,file='HamRsMLO',form='unformatted',status='old',action='read')
    read(ifh) n1,n2,n3        !ndimMTO,npairmx,nspx
    read(ifh)                 !hammr
    if(c0_socmatrix) read(ifh) !hammhsor
    read(ifh)                 !ovlmr
    read(ifh)                 !ib_tableM,k_tableM,l_tableM
    read(ifh) nd, ld, mm, nsk
    allocate(ixo(nd)); read(ifh) ixo
    read(ifh) f1, ef, ec
    close(ifh)
  end subroutine read_mloindex

  subroutine sigmlo_init()
    use m_readqplist,only: set_bandedge
    use m_cmdopt_registry,only: c0_mlo
    integer :: ifs, nd2, ld2
    logical :: lex1, lex2
    if(.not.init) return
    init = .false.
    if(.not.c0_mlo) return   !the MLO route is opt-in: leftover SigRsMLO in a directory
                             !must never silently change a conventional run
    inquire(file='SigRsMLO', exist=lex1)
    inquire(file='HamRsMLO', exist=lex2)
    if(.not.(lex1.and.lex2)) return
    call read_mloindex(nd2, ld2, mlomethod, nskip, ix, fff1, eferm, ecbot)
    open(newunit=ifs,file='SigRsMLO',form='unformatted',status='old',action='read')
    read(ifs) ndimMTO, npairmx, nspx, nbas
    if(ndimMTO /= nd2) call rx('m_sigmlo: SigRsMLO and __mloindex disagree on ndimMTO')
    allocate(sigmlor(npairmx,ndimMTO,ndimMTO,nspx))
    read(ifs) sigmlor
    read(ifs) plat
    allocate(npair(nbas,nbas), nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas))
    read(ifs) npair
    read(ifs) nlat
    read(ifs) nqwgt
    allocate(ib_tableM(ndimMTO)); read(ifs) ib_tableM, ix
    close(ifs)
    call set_bandedge(eferm, ecbot)   !Hreduction reads these for the MLO window
    call zmlo_sig_load()              !the chi~ this SigRsMLO was written in
    sigmlo_on = .true.
    write(stdo,ftox)' m_sigmlo: MLO Sigma interpolation ON. ndimMTO nskip=',ndimMTO,nskip, &
         ' |Sigma(R)|=',ftof(sum(abs(sigmlor)))
  end subroutine sigmlo_init

  !> Return the frozen z^MLO(q,isp) if ZmloRef has it; ok=.false. otherwise.  sugw's step
  !> a' must build Sigma^MLO in the SAME chi~ that getsenex will use it in.
  subroutine zmlo_frozen(q0, isp0, nm, nmlo_out, z0, ok)
    real(8),intent(in):: q0(3)
    integer,intent(in):: isp0, nm
    integer,intent(out):: nmlo_out
    complex(8),intent(out):: z0(nm,*)
    logical,intent(out):: ok
    integer:: jc
    ok=.false.; nmlo_out=ndimMTO
    if(nocache) return
    do jc=1,nzc
      if(ispzc(jc)==isp0 .and. ndzc(jc)==nm .and. sum(abs(qzc(:,jc)-q0))<1d-8) then
        z0(1:nm,1:ndimMTO) = zcache(1:nm,1:ndimMTO,jc); ok=.true.; return
      endif
    enddo
  end subroutine zmlo_frozen

  !> Put a z^MLO built elsewhere into the process-local cache (read side only).
  subroutine zmlo_store(q0, isp0, nm, z0)
    real(8),intent(in):: q0(3)
    integer,intent(in):: isp0, nm
    complex(8),intent(in):: z0(nm,ndimMTO)
    if(nocache) return
    call cache_put(q0, isp0, nm, z0)
  end subroutine zmlo_store

  !> ZmloSig.<procid>: the chi~ the current SigRsMLO was written in.  sugw's step a' runs
  !> under MPI and each rank owns a subset of q, so the slot is a SET of per-rank files
  !> (ndimh varies with q when pwmode=11, so a fixed-record direct-access file does not
  !> fit).  The rank count may differ between lmf runs, so scan a range and take what is
  !> there; every rank loads the whole set, which is a few MB.
  subroutine zmlo_sig_load()
    integer,parameter:: maxrank = 4096, gap_stop = 256
    integer:: ifz, nd, nm, isp0, ios, n, ip, miss
    real(8):: q0(3)
    complex(8),allocatable:: z0(:,:)
    logical:: lex
    character(256):: fn
    if(nocache) return
    n=0; miss=0
    do ip = 0, maxrank-1
      write(fn,"('ZmloSig.',i0)") ip
      inquire(file=trim(fn),exist=lex)
      if(.not.lex) then
        miss = miss + 1
        if(miss >= gap_stop) exit
        cycle
      endif
      miss = 0
      open(newunit=ifz,file=trim(fn),form='unformatted',status='old',action='read')
      read(ifz,iostat=ios) nd
      if(ios/=0 .or. nd/=ndimMTO) then
        close(ifz)
        if(ios==0 .and. nd/=ndimMTO) write(stdo,ftox) &
             ' m_sigmlo: ',trim(fn),' has ndimMTO=',nd,' /= ',ndimMTO,' -> ignored'
        cycle
      endif
      do
        read(ifz,iostat=ios) q0, isp0, nm
        if(ios/=0) exit
        allocate(z0(nm,ndimMTO))
        read(ifz,iostat=ios) z0
        if(ios/=0) then; deallocate(z0); exit; endif
        call cache_put(q0, isp0, nm, z0)
        deallocate(z0); n=n+1
      enddo
      close(ifz)
    enddo
    nzsig = n
    if(n>0) write(stdo,ftox)' m_sigmlo: loaded ZmloSig (chi~ of this SigRsMLO), records=',n
  end subroutine zmlo_sig_load

  !> ZmloNew.<procid>: the chi~ of THIS iteration, built by sugw's step a' from the H it
  !> is given.  gwsc renames the set to ZmloSig.* once `mlo` has written the SigRsMLO
  !> that was expressed in it -- that rename is the promotion.
  !> nd is passed in, NOT taken from the module: in the first iteration there is no
  !> SigRsMLO yet, so sigmlo_init returns early and the module's ndimMTO is still 0.
  !> Writing that as the header made the file unreadable (nd /= ndimMTO -> ignored).
  subroutine zmlo_new_append(q0, isp0, nm, nd, z0)
    use m_mpi,only: procid
    real(8),intent(in):: q0(3)
    integer,intent(in):: isp0, nm, nd
    complex(8),intent(in):: z0(nm,nd)
    integer:: ifz
    logical:: lex
    character(256):: fn
    if(nocache) return
    write(fn,"('ZmloNew.',i0)") procid
    inquire(file=trim(fn),exist=lex)
    if(lex) then
      open(newunit=ifz,file=trim(fn),form='unformatted',position='append')
    else
      open(newunit=ifz,file=trim(fn),form='unformatted')
      write(ifz) nd
    endif
    write(ifz) q0, isp0, nm
    write(ifz) z0
    close(ifz)
  end subroutine zmlo_new_append


  subroutine cache_put(q0, isp0, nm, z0)
    real(8),intent(in):: q0(3)
    integer,intent(in):: isp0, nm
    complex(8),intent(in):: z0(nm,ndimMTO)
    complex(8),allocatable:: zt(:,:,:)
    real(8),allocatable:: qt(:,:)
    integer,allocatable:: it(:), nt(:)
    integer:: ld
    if(nzc==0) then
      allocate(zcache(nm,ndimMTO,1), qzc(3,1), ispzc(1), ndzc(1))
      zcache(:,:,1)=z0; qzc(:,1)=q0; ispzc(1)=isp0; ndzc(1)=nm; nzc=1
    else
      ld=max(nm,size(zcache,1))
      allocate(zt(ld,ndimMTO,nzc+1), source=(0d0,0d0))
      allocate(qt(3,nzc+1)); allocate(it(nzc+1)); allocate(nt(nzc+1))
      zt(1:size(zcache,1),:,1:nzc)=zcache; qt(:,1:nzc)=qzc; it(1:nzc)=ispzc; nt(1:nzc)=ndzc
      zt(1:nm,:,nzc+1)=z0; qt(:,nzc+1)=q0; it(nzc+1)=isp0; nt(nzc+1)=nm
      call move_alloc(zt,zcache); call move_alloc(qt,qzc); call move_alloc(it,ispzc)
      call move_alloc(nt,ndzc)
      nzc=nzc+1
    endif
  end subroutine cache_put

  !> senex(ndimh,ndimh) = the PMT matrix element of Sigma, eq (14).
  subroutine sigmlo_senex(qp, isp, ndimh, ovlm, hamm, senex)
    use m_hreduction,only: Hreduction
    real(8),intent(in) :: qp(3)
    integer,intent(in) :: isp, ndimh
    complex(8),intent(in) :: ovlm(ndimh,ndimh), hamm(ndimh,ndimh)
    complex(8),intent(out) :: senex(ndimh,ndimh)
    integer :: i, j, it, ib1, ib2, jsp, nxq
    real(8),parameter :: pi=4d0*atan(1d0)
    complex(8),parameter :: img=(0d0,1d0)
    complex(8) :: sigk(ndimMTO,ndimMTO), omlo(ndimMTO,ndimMTO)
    complex(8) :: hmo(ndimMTO,ndimMTO), omo(ndimMTO,ndimMTO)
    complex(8),allocatable :: zm(:,:), amat(:,:), hl(:,:), ol(:,:), tmp(:,:)
    complex(8) :: ph
    jsp = min(isp, nspx)
    if(.false.) continue
    BlochSum: block  !eq (12); same phase convention as m_mlo_ham::calc_ham_eigen
      sigk = (0d0,0d0)
      do i = 1, ndimMTO
        ib1 = ib_tableM(i)
        do j = 1, ndimMTO
          ib2 = ib_tableM(j)
          do it = 1, npair(ib1,ib2)
            ph = 1d0/dble(nqwgt(it,ib1,ib2)) * exp(-img*2d0*pi*sum(qp*matmul(plat,dble(nlat(:,it,ib1,ib2)))))
            sigk(i,j) = sigk(i,j) + sigmlor(it,i,j,jsp)*ph
          enddo
        enddo
      enddo
      sigk = 0.5d0*(sigk + transpose(dconjg(sigk)))  !kill the residual anti-hermitian part
    endblock BlochSum
    RoundTripCheck: block !ECALJ_SIGMLO_RT=1: at a mesh q, the Bloch sum of Sigma^MLO(R) must
      !reproduce the Sigma^MLO(q) that hqpe_sc wrote (eq 11 -> FFT -> eq 12 is a round trip).
      character(32):: cv
      integer:: st, ifs, nm2,nq2,ns2,n1,n2,n3, ip, isx, iqm
      logical,save:: rt=.false., rfirst=.true.
      real(8),allocatable:: qs(:,:,:)
      complex(8),allocatable:: sq(:,:,:,:)
      logical:: lex
      if(rfirst) then
        rfirst=.false.
        call get_environment_variable('ECALJ_SIGMLO_RT',cv,status=st)
        rt = (st==0 .and. len_trim(cv)>0)
      endif
      if(rt) then
        inquire(file='__SigmMLO.q',exist=lex)
        if(lex) then
          open(newunit=ifs,file='__SigmMLO.q',form='unformatted',status='old',action='read')
          read(ifs) nm2,nq2,ns2,n1,n2,n3
          allocate(qs(3,ns2,nq2), sq(nm2,nm2,nq2,ns2))
          read(ifs) qs; read(ifs) sq; close(ifs)
          iqm=0
          do ip=1,nq2
            if(sum(abs(qs(:,1,ip)-qp))<1d-6) then; iqm=ip; exit; endif
          enddo
          if(iqm>0) then
            write(stdo,"(a,3f9.5,a,f12.4,a,f12.4)")' SIGMLO_RT q=',qp, &
              '  max|Bloch(R)-q-space|[meV]=',maxval(abs(sigk-sq(:,:,iqm,jsp)))*13605.7d0, &
              '   max|q-space|[meV]=',maxval(abs(sq(:,:,iqm,jsp)))*13605.7d0
          endif
          deallocate(qs,sq)
        endif
      endif
    endblock RoundTripCheck
    allocate(zm(ndimh,ndimMTO))
    ChiTildeCache: block
      character(32):: cv
      integer:: st, ic, jc
      if(cfirst) then
        cfirst=.false.
        call get_environment_variable('ECALJ_MLO_NOCACHE',cv,status=st)
        nocache = (st==0 .and. len_trim(cv)>0)
        if(nocache) write(stdo,ftox)' m_sigmlo: chi~ cache DISABLED (ECALJ_MLO_NOCACHE)'
      endif
      ic = 0
      if(.not.nocache) then
        do jc = 1, nzc
          if(ispzc(jc)==isp .and. ndzc(jc)==ndimh .and. sum(abs(qzc(:,jc)-qp))<1d-8) then
            ic = jc; exit
          endif
        enddo
      endif
      if(ic>0) then
        zm = zcache(1:ndimh,1:ndimMTO,ic)          !chi~ frozen: reuse
      else
        block
          complex(8):: hl(ndimh,ndimh), ol(ndimh,ndimh)
          hl = hamm; ol = ovlm                     !Hreduction may modify its arguments
          call Hreduction(mlomethod,.false.,ndimh, hl, ol, ndimMTO, ix, fff1, &
               hmo, omo, qp, nev=nxq, zMLO=zm, nskip_auto=nskip)
        endblock
        if(.not.nocache) call cache_put(qp, isp, ndimh, zm)   !one chi~ for this lmf run
      endif
    endblock ChiTildeCache
    ZmloDump: block !ECALJ_ZMLO_DUMP=1: write z^MLO(k) so it can be compared with the one
      !that sugw's step a' built at the same q.  They must be identical: both are
      !Hreduction(H^LDA(q), S^PMT(q)) with the same frozen index.
      character(32):: cv
      integer:: st, ifz
      logical,save:: dmp=.false., dfirst=.true.
      if(dfirst) then
        dfirst=.false.
        call get_environment_variable('ECALJ_ZMLO_DUMP',cv,status=st)
        dmp = (st==0 .and. len_trim(cv)>0)
        if(dmp) then
          open(newunit=ifz,file='__zmlo_getsenex',form='unformatted')
          close(ifz,status='delete')
        endif
      endif
      if(dmp) then
        open(newunit=ifz,file='__zmlo_getsenex',form='unformatted',position='append')
        write(ifz) qp, isp, ndimh, ndimMTO
        write(ifz) zm(1:ndimh,1:ndimMTO)
        close(ifz)
      endif
    endblock ZmloDump
    allocate(amat(ndimh,ndimMTO))
    amat = matmul(ovlm, zm)                    !A = S^PMT z^MLO           eq (15)
    omlo = matmul(transpose(dconjg(zm)), amat) !O^MLO = z^dag A
    call matcinv(ndimMTO, omlo)
    allocate(tmp(ndimMTO,ndimMTO))
    tmp = matmul(omlo, matmul(sigk, omlo))     !(O^-1) Sigma (O^-1)
    senex = matmul(amat, matmul(tmp, transpose(dconjg(amat))))   !eq (14)
    deallocate(amat,tmp,zm)
  end subroutine sigmlo_senex
end module m_sigmlo
