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
  function get_formfactor_q(q, ispin1, ispin2) result(formfactor)
    real(8), intent(in) :: q(3)
    integer, intent(in) :: ispin1, ispin2
    complex(8) :: formfactor(nnmlo)
    integer, save :: ifile = -1, nq0i_f = 0
    real(8), save, allocatable :: q0i_f(:,:)
    logical, save :: computed_f(4) = .false.
    complex(8) :: formfactor_f(nmlo, nmlo)
    integer :: iq, ifile_info, nmlo_f, nqbz_f, nspin_f, nbb_f, recl, sidx
    sidx = spin_idx(ispin1, ispin2)
    if(ifile < 0) then
      open(newunit=ifile_info, file='__MLOFormFactorQ.info', form='unformatted', action='read')
      read(ifile_info) computed_f
      read(ifile_info) nmlo_f, nqbz_f, nspin_f, nq0i_f
      if(nmlo_f /= nmlo) call rx('get_formfactor_q: nmlo mismatch')
      if(allocated(q0i_f)) deallocate(q0i_f)
      allocate(q0i_f(3, nq0i_f))
      read(ifile_info) q0i_f(:,:)
      close(ifile_info)
      recl = nmlo*nmlo*16
      open(newunit=ifile, file='__MLOFormFactorQ', form='unformatted', access='direct', recl=recl, action='read')
    endif
    if(.not. computed_f(sidx)) call rx('get_formfactor_q: requested spin pair not computed')
    do iq = 1, nq0i_f
      if(all(abs(q0i_f(:,iq) - q) < 1d-8)) then
        read(ifile, rec=(sidx-1)*nq0i_f+iq) formfactor_f(:,:)
        formfactor(1:nnmlo) = pack(reshape(formfactor_f(:,:), shape=[nmlo*nmlo]), mask = nnmlo_mask)
        return
      endif
    enddo
    call rx('q not found in get_formfactor_q')
  end function get_formfactor_q

  pure function spin_idx(ispin1, ispin2) result(idx)
    integer, intent(in) :: ispin1, ispin2
    integer :: idx
    ! (1,1)→1 UPUP, (2,2)→2 DNDN, (1,2)→3 UPDN, (2,1)→4 DNUP
    idx = 0
    if(ispin1==1 .and. ispin2==1) idx = 1
    if(ispin1==2 .and. ispin2==2) idx = 2
    if(ispin1==1 .and. ispin2==2) idx = 3
    if(ispin1==2 .and. ispin2==1) idx = 4
  end function spin_idx
end module m_mlo_formfactor
