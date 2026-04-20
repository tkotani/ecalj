module m_mlo_formfactor
  use m_mlo_ham, only: nmlo => ndimMTO, ib_tableM, ib_tableI, nsite
  use m_mlo_scrw, only: nnmlo => nnwf, nnmlo_mask => nnwf_mask
  use m_lgunit, only: stdo
  use m_ftox, only: ftox
  use m_mpi, only: ipr
  implicit none
  public :: get_formfactor_q
  private
contains
  function get_formfactor_q(q, isp, spinflip) result(formfactor)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    complex(8) :: formfactor(nnmlo)
    integer, save :: ifile = -1, isp_prev = -1, nq0i_f = 0
    logical, save :: spinflip_prev = .false.
    real(8), save, allocatable :: q0i_f(:,:)
    complex(8) :: formfactor_f(nmlo, nmlo)
    integer :: iq, ifile_info, nmlo_f, nqbz_f, nspin_f, recl
    logical :: opened
    if(isp /= isp_prev .or. spinflip .neqv. spinflip_prev) then
      open(newunit=ifile_info, file='__MLOFormFactorQ.info', form='unformatted', action='read')
      read(ifile_info) nmlo_f, nqbz_f, nspin_f, nq0i_f
      if(nmlo_f /= nmlo) call rx('get_formfactor_q: nmlo mismatch')
      if(allocated(q0i_f)) deallocate(q0i_f)
      allocate(q0i_f(3, nq0i_f))
      read(ifile_info) q0i_f(:,:)
      close(ifile_info)
      inquire(unit=ifile, opened=opened)
      if(opened) close(ifile)
      recl = nmlo*nmlo*16
      open(newunit=ifile, file=formfactor_q_fname(isp, spinflip), form='unformatted', access='direct', recl=recl, action='read')
      isp_prev = isp
      spinflip_prev = spinflip
    endif
    do iq = 1, nq0i_f
      if(all(abs(q0i_f(:,iq) - q) < 1d-8)) then
        read(ifile, rec=iq) formfactor_f(:,:)
        formfactor(1:nnmlo) = pack(reshape(formfactor_f(:,:), shape=[nmlo*nmlo]), mask = nnmlo_mask)
        return
      endif
    enddo
    call rx('q not found in get_formfactor_q')
  end function get_formfactor_q

  pure function formfactor_q_fname(isp, spinflip) result(fname)
    integer, intent(in) :: isp
    logical, intent(in) :: spinflip
    character(:), allocatable :: fname
    if(isp==1 .and. .not.spinflip) fname = '__MLOFormFactorQ.UP'
    if(isp==2 .and. .not.spinflip) fname = '__MLOFormFactorQ.DN'
    if(isp==1 .and. spinflip)      fname = '__MLOFormFactorQ.UPDN'
    if(isp==2 .and. spinflip)      fname = '__MLOFormFactorQ.DNUP'
  end function formfactor_q_fname
end module m_mlo_formfactor
