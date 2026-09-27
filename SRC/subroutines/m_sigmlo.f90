!> Sigma in the MLO representation: read QMLO_SigRs, Bloch-sum it at any k, and
!> build the PMT matrix element (design: Samples/kBT/sigma_mlo_design.md, stage 4).
!>
!>   Sigma^MLO(k)   = sum_R Sigma^MLO(R) exp(-ikR)          eq (12)
!>   z^MLO(k)       = Hreduction(H(k), S^PMT(k))            eq (17), never interpolated
!>   A(k)           = S^PMT(k) z^MLO(k),  O^MLO = z^dag A   eq (15)
!>   [Sigma]_mn(k)  = A (O^MLO)^-1 Sigma^MLO(k) (O^MLO)^-1 A^dag   eq (14)
!>
!> QMLO_SigRs is self-contained (it carries the pair list), so nothing here depends
!> on m_HamPMT's module state.
!>
!> Files of MLO-QSGW.  QMLO_* are kept: results, and what a chain needs to continue.
!> __QMLO_* are rebuilt within one iteration and may be deleted between gwsc calls.
module m_sigmlo
  use m_lgunit,only: stdo
  use m_ftox
  implicit none
  public :: sigmlo_init, sigmlo_senex, sigmlo_on, read_mloindex
  public :: sigmlo_sigq, sigmlo_nmlo
  public :: zmlo_new_open, zmlo_new_write, zmlo_new_close
  character(*),parameter,public:: fn_sigrs  = 'QMLO_SigRs'    !Sigma^MLO(R), written by `mlo` (step 4)
  character(*),parameter,public:: fn_z      = 'QMLO_z'        !z^MLO that QMLO_SigRs was written in
  character(*),parameter,public:: fn_znew   = '__QMLO_zNew'   !z^MLO step a' builds for this iteration
  character(*),parameter,public:: fn_sigq   = '__QMLO_Sig'    !Sigma^MLO(q) of hqpe_sc (step 3), read by `mlo`
  character(*),parameter,public:: fn_mixsig = '__QMLO_mixsig' !Anderson history of Sigma^MLO (hqpe_sc)
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
  !      QMLO_z       the z^MLO the CURRENT QMLO_SigRs was written in.  getsenex reads
  !                   Sigma back through this one.
  !      __QMLO_zNew  the z^MLO sugw's step a' builds from the H of THIS iteration and uses
  !                   for c^MLO, so the next Sigma^MLO is written in it.  `mlo` renames it
  !                   to QMLO_z right after writing the new QMLO_SigRs.
  !
  !    Each slot is ONE file of fixed-length records, one per (q,isp) of the GW q list in
  !    the order of __cmlo.data: record 1+isp+nspx*(iq-1).  Record 1 holds the layout.  A
  !    record is q(3), isp, ndimh and then z(nbandmx,ndimMTO) zero-padded, since ndimh
  !    depends on q when pwmode=11.  Step a' writes it with MPI-IO, every rank its own q.
  !    getsenex reads the index once and afterwards only the records of the k it meets,
  !    with Fortran stream access: a process that owns no k never calls getsenex, so it
  !    could not take part in a collective MPI-IO open.
  !
  !    At a k that is not in QMLO_z (a band plot, or an SCF mesh finer than the Sigma
  !    mesh) there is nothing to read back exactly -- the H of that iteration is gone --
  !    so the run stops unless ECALJ_MLO_ALLOW_REBUILD=1, which rebuilds z from the H in
  !    hand and keeps it in memory for the rest of this run.  ECALJ_MLO_NOCACHE=1 disables
  !    the slot and the memory cache and rebuilds at every call (the 2026-09-25 failure
  !    mode: chi~ then moves inside one SCF loop).
  integer,parameter :: zmagic = 20260926, zhead = 32   !tag of record 1; bytes before z in a record
  logical :: zslot = .false.                !QMLO_z is open and indexed
  integer :: ifzs = -1, nzrec = 0, nzbm = 0
  integer(8) :: lzrec = 0
  real(8), allocatable :: qzs(:,:)
  integer, allocatable :: ispzs(:), ndzs(:) !ndzs=0: a record that was never written
  integer :: ifzn = -1, nznbm = 0, nznd = 0, lznrec = 0   !__QMLO_zNew, written by step a'
  !--- memory cache of the z^MLO this process has met (read from QMLO_z or rebuilt).
  !    ndimh depends on q when pwmode=11, so an entry may only be reused at its own size;
  !    taking the first rows of a zero-padded entry gave a singular O^MLO and NaN.
  integer :: nzc = 0
  real(8), allocatable :: qzc(:,:)
  integer, allocatable :: ispzc(:), ndzc(:)   !ndzc: the ndimh each entry was built with.
  complex(8), allocatable :: zcache(:,:,:)
  logical :: nocache = .false.
  integer :: nmiss = 0                        !k met that were NOT in QMLO_z
  logical :: allowrb = .false., mfirst = .true.
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

  !> withz=.false.: Sigma^MLO(R) only, no z^MLO slot (hqpe_sc needs the Bloch sum alone).
  subroutine sigmlo_init(withz)
    use m_readqplist,only: set_bandedge
    use m_cmdopt_registry,only: c0_mlo
    logical,intent(in),optional :: withz
    integer :: ifs, nd2, ld2, st
    logical :: lex1, lex2, lz
    character(32):: cv
    if(.not.init) return
    init = .false.
    if(.not.c0_mlo) return   !the MLO route is opt-in: leftover QMLO_SigRs in a directory
                             !must never silently change a conventional run
    inquire(file=fn_sigrs, exist=lex1)
    inquire(file='HamRsMLO', exist=lex2)
    if(.not.(lex1.and.lex2)) return
    call read_mloindex(nd2, ld2, mlomethod, nskip, ix, fff1, eferm, ecbot)
    open(newunit=ifs,file=fn_sigrs,form='unformatted',status='old',action='read')
    read(ifs) ndimMTO, npairmx, nspx, nbas
    if(ndimMTO /= nd2) call rx('m_sigmlo: '//fn_sigrs//' and HamRsMLO disagree on ndimMTO')
    allocate(sigmlor(npairmx,ndimMTO,ndimMTO,nspx))
    read(ifs) sigmlor
    read(ifs) plat
    allocate(npair(nbas,nbas), nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas))
    read(ifs) npair
    read(ifs) nlat
    read(ifs) nqwgt
    allocate(ib_tableM(ndimMTO)); read(ifs) ib_tableM, ix
    close(ifs)
    WindowReference: block !Hreduction reads eferm/ecbot of m_readqplist for the MLO window:
      !for a chi~ rebuilt here (diagnostics) and, in lmf --jobgw=1, for the chi~ step a' builds
      !right after the first getsenex.  Take them from the Ef file of the last SCF, as step a'
      !does.  HamRsMLO holds the chain-start (LDA) values; setting those here overrode step
      !a' from iteration 2 on, so its chi~ kept the LDA window (found 2026-09-26).
      use m_readqplist,only: readbandedge
      use m_cmdopt_registry,only: efermi_file
      real(8):: ef
      integer:: ifi, ios
      open(newunit=ifi,file=efermi_file(),status='old',action='read',iostat=ios)
      if(ios==0) then
        read(ifi,*,iostat=ios) ef
        close(ifi)
      endif
      if(ios==0) then
        call set_bandedge(ef, ef)
        call readbandedge()             !ecbot = ef + (ecbot - eferm) of that file
      else
        call set_bandedge(eferm, ecbot)
        write(stdo,ftox)' m_sigmlo: no ',efermi_file(),' -> MLO window from HamRsMLO (chain start)'
      endif
    endblock WindowReference
    call get_environment_variable('ECALJ_MLO_NOCACHE',cv,status=st)
    nocache = (st==0 .and. len_trim(cv)>0)
    if(nocache) write(stdo,ftox)' m_sigmlo: chi~ slot and cache DISABLED (ECALJ_MLO_NOCACHE)'
    lz = .true.
    if(present(withz)) lz = withz
    if(lz) call zslot_open()          !the chi~ this QMLO_SigRs was written in
    sigmlo_on = .true.
    write(stdo,ftox)' m_sigmlo: MLO Sigma interpolation ON. ndimMTO nskip=',ndimMTO,nskip, &
         ' |Sigma(R)|=',ftof(sum(abs(sigmlor)))
  end subroutine sigmlo_init

  integer function sigmlo_nmlo()
    sigmlo_nmlo = ndimMTO
  end function sigmlo_nmlo

  !> Open QMLO_z and read the (q, isp, ndimh) of every record; z itself is read on demand
  !> by zslot_fetch.  Nothing to open is not an error: a band run, or a directory where
  !> the chi~ slot was deliberately left out.
  subroutine zslot_open()
    integer:: magic, nrec, nbm, nd, lrec, idum(3), ios, irec, st
    integer(8):: pos
    logical:: lex
    character(32):: cv
    if(nocache) return
    !ECALJ_MLO_NOSIG=1: ignore the slot and rebuild chi~ at every k from the H in hand.
    !On the Sigma q mesh the slot is the right thing -- it is the basis Sigma was written
    !in.  A band plot at k off that mesh cannot be covered by the slot; draw_mloband.sh
    !(the MLO model of the SCF Hamiltonian) is the way to draw those.  Diagnostic only.
    call get_environment_variable('ECALJ_MLO_NOSIG',cv,status=st)
    if(st==0 .and. len_trim(cv)>0) then
      write(stdo,ftox)' m_sigmlo: ECALJ_MLO_NOSIG -> ',fn_z,' ignored, chi~ rebuilt at every k'
      return
    endif
    inquire(file=fn_z,exist=lex)
    if(.not.lex) return
    open(newunit=ifzs,file=fn_z,access='stream',form='unformatted',status='old',action='read',iostat=ios)
    if(ios/=0) return
    read(ifzs,pos=1,iostat=ios) magic, nrec, nbm, nd, lrec, idum
    if(ios/=0 .or. magic/=zmagic) then
      write(stdo,ftox)' m_sigmlo: ',fn_z,' has no valid header -> ignored'
      close(ifzs); return
    endif
    if(nd/=ndimMTO) then
      write(stdo,ftox)' m_sigmlo: ',fn_z,' has ndimMTO=',nd,' /= ',ndimMTO,' -> ignored'
      close(ifzs); return
    endif
    allocate(qzs(3,nrec), ispzs(nrec), ndzs(nrec))
    do irec = 1, nrec
      pos = int(irec,8)*int(lrec,8) + 1
      read(ifzs,pos=pos,iostat=ios) qzs(:,irec), ispzs(irec), ndzs(irec)
      if(ios/=0) then                      !past the end: never written
        qzs(:,irec) = 0d0; ispzs(irec) = 0; ndzs(irec) = 0
      endif
    enddo
    nzrec = nrec; nzbm = nbm; lzrec = lrec; zslot = .true.
    write(stdo,ftox)' m_sigmlo: ',fn_z,' (chi~ of this ',fn_sigrs,') indexed, records=',count(ndzs>0)
  end subroutine zslot_open

  !> z^MLO(qp,isp) from QMLO_z, if a record with this q, spin and ndimh is there.
  subroutine zslot_fetch(qp, isp, ndimh, zm, found)
    real(8),intent(in):: qp(3)
    integer,intent(in):: isp, ndimh
    complex(8),intent(out):: zm(ndimh,ndimMTO)
    logical,intent(out):: found
    complex(8),allocatable:: zb(:,:)
    integer(8):: pos
    integer:: irec
    found = .false.
    if(.not.zslot) return
    do irec = 1, nzrec
      if(ndzs(irec)==ndimh .and. ispzs(irec)==isp .and. sum(abs(qzs(:,irec)-qp))<1d-8) then
        allocate(zb(nzbm,ndimMTO))
        pos = int(irec,8)*lzrec + 1 + zhead
        read(ifzs,pos=pos) zb
        zm = zb(1:ndimh,1:ndimMTO)
        found = .true.
        return
      endif
    enddo
  end subroutine zslot_fetch

  !> __QMLO_zNew: the chi~ of THIS iteration, built by sugw's step a' from the H it is
  !> given.  `mlo` renames it to QMLO_z once it has written the QMLO_SigRs that was
  !> expressed in it -- that rename is the promotion.  Collective: every rank of step a'.
  !> nd is passed in, NOT taken from the module: in the first iteration there is no
  !> QMLO_SigRs yet, so sigmlo_init returns early and the module's ndimMTO is still 0.
  subroutine zmlo_new_open(nrec, nbm, nd)
    use mpi
    use m_mpi,only: master_mpi, comm
    use m_mpiio,only: openm, writem_buf, mpiio_buf, buf_put
    integer,intent(in):: nrec, nbm, nd
    type(mpiio_buf):: buf
    complex(8),allocatable:: zpad(:)
    integer:: istat, ifi, ierr
    logical:: lex
    if(master_mpi) then  !an MPI-IO open does not truncate: remove a previous file first
      inquire(file=fn_znew,exist=lex)
      if(lex) then
        open(newunit=ifi,file=fn_znew)
        close(ifi,status='delete')
      endif
    endif
    call mpi_barrier(comm,ierr)
    nznbm = nbm; nznd = nd; lznrec = zhead + 16*nbm*nd
    istat = openm(newunit=ifzn, file=fn_znew, recl=lznrec)
    if(master_mpi) then  !record 1: the layout the reader needs before anything else
      allocate(zpad(nbm*nd), source=(0d0,0d0))
      call buf_put(buf, [zmagic, nrec, nbm, nd, lznrec, 0, 0, 0])
      call buf_put(buf, zpad)
      istat = writem_buf(ifzn, 1, buf)
    endif
  end subroutine zmlo_new_open

  !> Record irec (= isp+nspx*(iq-1), as in __cmlo.data) of __QMLO_zNew.
  subroutine zmlo_new_write(irec, q0, isp0, nm, z0)
    use m_mpiio,only: writem_buf, mpiio_buf, buf_put
    integer,intent(in):: irec, isp0, nm
    real(8),intent(in):: q0(3)
    complex(8),intent(in):: z0(nm,nznd)
    type(mpiio_buf):: buf
    complex(8),allocatable:: zp(:,:)
    integer:: istat
    if(nocache) return
    allocate(zp(nznbm,nznd), source=(0d0,0d0))
    zp(1:nm,:) = z0
    call buf_put(buf, q0)
    call buf_put(buf, isp0)
    call buf_put(buf, nm)
    call buf_put(buf, zp)
    istat = writem_buf(ifzn, 1+irec, buf)
  end subroutine zmlo_new_write

  !> Collective, like zmlo_new_open.
  subroutine zmlo_new_close()
    use m_mpiio,only: closem
    integer:: istat
    istat = closem(ifzn)
  end subroutine zmlo_new_close

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

  !> Sigma^MLO(qp) Bloch-summed from QMLO_SigRs, eq (12): the matrix getsenex puts back
  !> into the PMT basis, hermitian part only.
  subroutine sigmlo_sigq(qp, isp, sigk)
    real(8),intent(in) :: qp(3)
    integer,intent(in) :: isp
    complex(8),intent(out) :: sigk(ndimMTO,ndimMTO)
    integer :: i, j, it, ib1, ib2, jsp
    real(8),parameter :: pi=4d0*atan(1d0)
    complex(8),parameter :: img=(0d0,1d0)
    complex(8) :: s
    complex(8),allocatable :: ph(:,:,:)
    jsp = min(isp, nspx)
    allocate(ph(npairmx,nbas,nbas))
    do ib2 = 1, nbas       !the phase depends on the atom pair and R only (it was made for every
      do ib1 = 1, nbas     !(i,j,R): 0.16 s per k for LiTi2O4, 2.5 s of the 4 s of hqpe_sc)
        do it = 1, npair(ib1,ib2)   !same phase convention as m_mlo_ham::calc_ham_eigen
          ph(it,ib1,ib2) = 1d0/dble(nqwgt(it,ib1,ib2)) * exp(-img*2d0*pi*sum(qp*matmul(plat,dble(nlat(:,it,ib1,ib2)))))
        enddo
      enddo
    enddo
    do j = 1, ndimMTO
      ib2 = ib_tableM(j)
      do i = 1, ndimMTO
        ib1 = ib_tableM(i)
        s = (0d0,0d0)
        do it = 1, npair(ib1,ib2)
          s = s + sigmlor(it,i,j,jsp)*ph(it,ib1,ib2)
        enddo
        sigk(i,j) = s
      enddo
    enddo
    sigk = 0.5d0*(sigk + transpose(dconjg(sigk)))  !kill the residual anti-hermitian part
  end subroutine sigmlo_sigq

  !> senex(ndimh,ndimh) = the PMT matrix element of Sigma, eq (14).
  subroutine sigmlo_senex(qp, isp, ndimh, ovlm, hamm, senex)
    use m_hreduction,only: Hreduction
    real(8),intent(in) :: qp(3)
    integer,intent(in) :: isp, ndimh
    complex(8),intent(in) :: ovlm(ndimh,ndimh), hamm(ndimh,ndimh)
    complex(8),intent(out) :: senex(ndimh,ndimh)
    integer :: jsp, nxq
    complex(8) :: sigk(ndimMTO,ndimMTO), omlo(ndimMTO,ndimMTO)
    complex(8) :: hmo(ndimMTO,ndimMTO), omo(ndimMTO,ndimMTO)
    complex(8),allocatable :: zm(:,:), amat(:,:), tmp(:,:)
    jsp = min(isp, nspx)
    call sigmlo_sigq(qp, isp, sigk)
    RoundTripCheck: block !ECALJ_SIGMLO_RT=1: at a mesh q, the Bloch sum of Sigma^MLO(R) must
      !reproduce the Sigma^MLO(q) that hqpe_sc wrote (eq 11 -> FFT -> eq 12 is a round trip).
      character(32):: cv
      integer:: st, ifs, nm2,nq2,ns2,n1,n2,n3, ip, iqm
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
        inquire(file=fn_sigq,exist=lex)
        if(lex) then
          open(newunit=ifs,file=fn_sigq,form='unformatted',status='old',action='read')
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
      logical:: found
      ic = 0
      if(.not.nocache) then
        do jc = 1, nzc
          if(ispzc(jc)==isp .and. ndzc(jc)==ndimh .and. sum(abs(qzc(:,jc)-qp))<1d-8) then
            ic = jc; exit
          endif
        enddo
      endif
      if(ic>0) then
        zm = zcache(1:ndimh,1:ndimMTO,ic)          !met before in this run
      else
        call zslot_fetch(qp, isp, ndimh, zm, found)
        if(.not.found) then
          !A MISS when QMLO_z is there means this k is outside the set that chi~ was
          !built on when the current QMLO_SigRs was written.  Rebuilding here would read
          !Sigma back in a DIFFERENT basis than it was written in -- the 2026-09-25 failure
          !mode, and silent.  Count it and say so.  ECALJ_MLO_ALLOW_REBUILD=1 rebuilds anyway.
          if(zslot) then
            nmiss = nmiss + 1
            if(mfirst) then
              mfirst=.false.
              call get_environment_variable('ECALJ_MLO_ALLOW_REBUILD',cv,status=st)
              allowrb = (st==0 .and. len_trim(cv)>0)
              write(stdo,ftox)' m_sigmlo: WARNING chi~ MISS at q=',ftof(qp),' isp ndimh=',isp,ndimh, &
                   ' -- this k is not in ',fn_z,'; Sigma would be read in a rebuilt basis.', &
                   ' allow_rebuild=',allowrb
            endif
            if(.not.allowrb) call rx('m_sigmlo: k outside '//fn_z//'. '// &
                 'Set ECALJ_MLO_ALLOW_REBUILD=1 to rebuild anyway (basis mismatch).')
          endif
          block
            complex(8):: hl(ndimh,ndimh), ol(ndimh,ndimh)
            hl = hamm; ol = ovlm                     !Hreduction may modify its arguments
            call Hreduction(mlomethod,.false.,ndimh, hl, ol, ndimMTO, ix, fff1, &
                 hmo, omo, qp, nev=nxq, zMLO=zm, nskip_auto=nskip)
          endblock
        endif
        if(.not.nocache) call cache_put(qp, isp, ndimh, zm)   !one chi~ per k for this lmf run
      endif
    endblock ChiTildeCache
    ZmloDump: block !ECALJ_ZMLO_DUMP=1: write z^MLO(k) so it can be compared with the one
      !that sugw's step a' built at the same q (__QMLO_zdump_sugw).  In the SCF lmf of an
      !iteration the two must agree: that z is exactly what step a' wrote.
      character(32):: cv
      integer:: st, ifz
      logical,save:: dmp=.false., dfirst=.true.
      if(dfirst) then
        dfirst=.false.
        call get_environment_variable('ECALJ_ZMLO_DUMP',cv,status=st)
        dmp = (st==0 .and. len_trim(cv)>0)
        if(dmp) then
          open(newunit=ifz,file='__QMLO_zdump_getsenex',form='unformatted')
          close(ifz,status='delete')
        endif
      endif
      if(dmp) then
        open(newunit=ifz,file='__QMLO_zdump_getsenex',form='unformatted',position='append')
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
