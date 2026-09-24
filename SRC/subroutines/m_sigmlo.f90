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
  public :: sigmlo_init, sigmlo_senex, sigmlo_on, read_mloindex
  logical, protected :: sigmlo_on = .false.
  private
  integer :: ndimMTO=0, npairmx=0, nspx=0, nbas=0, mlomethod=4, nskip=0
  integer, allocatable :: npair(:,:), nlat(:,:,:,:), nqwgt(:,:,:), ib_tableM(:), ix(:)
  real(8) :: plat(3,3), fff1=2d0, eferm=0d0, ecbot=0d0
  complex(8), allocatable :: sigmlor(:,:,:,:)   !(npairmx, ndimMTO, ndimMTO, nspx)
  logical :: init = .true.
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
    sigmlo_on = .true.
    write(stdo,ftox)' m_sigmlo: MLO Sigma interpolation ON. ndimMTO nskip=',ndimMTO,nskip, &
         ' |Sigma(R)|=',ftof(sum(abs(sigmlor)))
  end subroutine sigmlo_init

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
    allocate(hl(ndimh,ndimh),ol(ndimh,ndimh),zm(ndimh,ndimMTO))
    hl = hamm; ol = ovlm                       !Hreduction may modify its arguments
    call Hreduction(mlomethod,.false.,ndimh, hl, ol, ndimMTO, ix, fff1, &
         hmo, omo, qp, nev=nxq, zMLO=zm, nskip_auto=nskip)
    deallocate(hl,ol)
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
