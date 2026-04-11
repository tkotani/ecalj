!> read PPOVLGG, PPOVLG, PPOVLI  ngc2, ppx(1:ngc,1:ngc2), ngvecc2(1:3,1:ngc2) are returned.
!! 2026-04-11 Note: refactored. For old code, see git history before this change.
module m_read_ppovl
  use m_lgunit, only: stdo
  use m_nvfortran, only: findloc
  use m_ftox
  implicit none
  integer, public, protected :: nggg, ngcgp, ngcread, nxi, nxe, nyi, nye, nzi, nze, ngc2
  complex(8), public, protected, allocatable :: ppx(:,:), ggg(:), ppovlinv(:,:)
  integer, public, protected, allocatable :: ngvecc2(:,:), nvggg(:,:), nvgcgp2(:,:), ngvecc(:,:), igggi(:,:,:), igcgp2i(:,:,:)
  integer, public, protected :: nnxi, nnxe, nnyi, nnye, nnzi, nnze
  public :: getppx2
  private
  integer :: nqq, ngggmx, nqini, nqnumt
  logical :: ppovlclosed=.true., init=.true.
  logical :: debug=.false.
  real(8), allocatable :: qxtable(:,:)
  logical :: ippovlggooo=.true.
  integer :: ngcmax, ippovl, ippovlg, ippovlp_info
  integer, allocatable :: ngctable(:)
contains
  subroutine getppx2(qi, getngcgp) ! This return nvggg, nvgcgp2, ngvecc,  nggg, ngcgp, ngcread, ggg, ppovlinv
    real(8), intent(in) :: qi(3)
    integer :: iqi, ippovlgg
    integer :: verbose
    integer :: iqi0, igcgp2, iggg
    integer :: istat
    logical, optional :: getngcgp
    if(verbose() >= 100) debug = .TRUE.
    if(present(getngcgp)) then
       open(newunit=ippovlgg, file="__PPOvlpGG", form='unformatted')
       read(ippovlgg) nggg, ngcgp, nqq, nqini, nqnumt
       close(ippovlgg)
       return
    endif
    if(ippovlggooo) then !!  Make igggi inversion table
       open(newunit=ippovlgg, file="__PPOvlpGG", form='unformatted')
       read(ippovlgg) nggg, ngcgp, nqq, nqini, nqnumt
       if(debug) write(stdo, "('Readin getppx2: nggg ngcgp nqq=',3i10)") nggg, ngcgp, nqq
       allocate(nvggg(1:3,1:nggg), ggg(1:nggg), nvgcgp2(1:3,ngcgp))
       read(ippovlgg) nvgcgp2(1:3,1:ngcgp)
       read(ippovlgg) nvggg(1:3,1:nggg)
       read(ippovlgg) ggg(1:nggg)
       close(ippovlgg)
       nxi = minval(nvggg(1,1:nggg))
       nxe = maxval(nvggg(1,1:nggg))
       nyi = minval(nvggg(2,1:nggg))
       nye = maxval(nvggg(2,1:nggg))
       nzi = minval(nvggg(3,1:nggg))
       nze = maxval(nvggg(3,1:nggg))
       allocate(igggi(nxi:nxe, nyi:nye, nzi:nze), source=-100000)
       forall(iggg=1:nggg) igggi(nvggg(1,iggg), nvggg(2,iggg), nvggg(3,iggg)) = iggg
       nnxi = minval(nvgcgp2(1,1:ngcgp))
       nnxe = maxval(nvgcgp2(1,1:ngcgp))
       nnyi = minval(nvgcgp2(2,1:ngcgp))
       nnye = maxval(nvgcgp2(2,1:ngcgp))
       nnzi = minval(nvgcgp2(3,1:ngcgp))
       nnze = maxval(nvgcgp2(3,1:ngcgp))
       allocate(igcgp2i(nnxi:nnxe, nnyi:nnye, nnzi:nnze), source=-100000)
       forall(igcgp2=1:ngcgp) igcgp2i(nvgcgp2(1,igcgp2), nvgcgp2(2,igcgp2), nvgcgp2(3,igcgp2)) = igcgp2 ! inversion table for nvgcgp2
       ippovlggooo = .false.
       allocate(qxtable(3,nqini:nqnumt), ngctable(nqini:nqnumt))
       open(newunit=ippovlp_info, file="__PPOvlp.info", form='unformatted')
       read(ippovlp_info) ngcmax
       do iqi = nqini, nqnumt
         read(ippovlp_info) qxtable(:, iqi), ngctable(iqi)
       enddo
       close(ippovlp_info)
       open(newunit=ippovl,  file="__PPOvlp",  access='direct', recl=16*ngcmax*ngcmax)
       open(newunit=ippovlg, file="__PPOvlpG", access='direct', recl=4*3*ngcmax)
       if(debug) write(stdo, "('init ok!:should be done only once')")
    endif
    ReadPPovlpData: block
      integer :: ngvecc_buf(3,ngcmax)
      iqi = findloc([(sum(abs(qxtable(:,iqi0)-qi))<1d-10, iqi0=nqini,nqnumt)], value=.true., dim=1) + nqini - 1
      if(iqi < nqini) call rx('rppovl.f90: qi is not found. some bug. qi='//ftof(qi))
      ngcread = ngctable(iqi)
      read(ippovlg, rec=iqi-nqini+1) ngvecc_buf
      if(allocated(ngvecc)) deallocate(ngvecc)
      allocate(ngvecc, source=ngvecc_buf(1:3,1:ngcread))
    endblock ReadPPovlpData
  end subroutine getppx2
end module m_read_ppovl
