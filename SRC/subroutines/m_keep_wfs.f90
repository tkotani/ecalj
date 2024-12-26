module m_keep_wfs
  use m_lmfinit,only: nspec
  use m_genallcf_v3,only: nband, ndima, nspc
  use m_hamindex,only: ngpmx
  implicit none
  public :: keep_wfs_init, keep_wfs_finalize, set_cphi_from_keep, set_geig_from_keep, update_keep_cphi, update_keep_geig
  private
  integer :: nkeep_wfs = 2 !default value
  logical, parameter :: debug = .false.
  complex(8), allocatable :: keep_cphi(:,:,:), keep_geig(:,:,:)
  integer, allocatable :: qs_cphi(:,:), priority_cphi(:), qs_geig(:,:), priority_geig(:)
  integer :: found_cphi, found_geig, missing_cphi, missing_geig
#ifdef __GPU
  attributes(device) :: keep_cphi, keep_geig
#endif
contains

  subroutine keep_wfs_init()
    use m_keyvalue,only: getkeyvalue
    call getkeyvalue("GWinput", "nkeep_wfs", nkeep_wfs, default=2)
    if(nkeep_wfs < 1) return
    allocate(keep_cphi(ndima*nspc,nband,nkeep_wfs))
    allocate(keep_geig(ngpmx*nspc,nband,nkeep_wfs))
    allocate(qs_cphi(nkeep_wfs,2), qs_geig(nkeep_wfs,2))
    allocate(priority_geig(nkeep_wfs), source = 0)
    allocate(priority_cphi(nkeep_wfs), source = 0)
    found_cphi = 0
    found_geig = 0
    missing_cphi = 0
    missing_geig = 0
  end subroutine keep_wfs_init

  subroutine keep_wfs_finalize()
    use m_lgunit, only:stdo
    use m_ftox
    deallocate(keep_cphi, keep_geig, qs_cphi, qs_geig, priority_cphi, priority_geig)
    write(stdo,ftox) 'keep_wfs: found_cphi=', found_cphi, ' missing_cphi=', missing_cphi, &
                     'ratio=', dble(found_cphi)/dble(found_cphi + missing_cphi)
    write(stdo,ftox) 'keep_wfs: found_geig=', found_geig, ' missing_geig=', missing_geig, &
                     'ratio=', dble(found_geig)/dble(found_geig + missing_geig)
  end subroutine keep_wfs_finalize

  logical function set_cphi_from_keep(iq, is, wfs) result(success)
    integer, intent(in) :: iq, is
    complex(8), intent(out) :: wfs(ndima*nspc, nband)
    integer :: wfs_index
#ifdef __GPU
    attributes(device) :: wfs
#endif
    success = .false.
    if(nkeep_wfs < 1) return
    wfs_index = get_keep_cphi_index(iq, is)
    if(wfs_index > 0) then
      !$acc kernels
      wfs(:,:) = keep_cphi(:, :, wfs_index)
      !$acc end kernels
      priority_cphi(:) = priority_cphi(:) - 1
      priority_cphi(wfs_index) = priority_cphi(wfs_index) + nkeep_wfs
      success = .true.
      found_cphi = found_cphi + 1
    else
      missing_cphi = missing_cphi + 1
    endif
  end function set_cphi_from_keep
  logical function set_geig_from_keep(iq, is, wfs) result(success)
    integer, intent(in) :: iq, is
    complex(8), intent(out) :: wfs(ngpmx*nspc, nband)
    integer :: wfs_index
#ifdef __GPU
    attributes(device) :: wfs
#endif
    success = .false.
    if(nkeep_wfs < 1) return
    wfs_index = get_keep_geig_index(iq, is)
    if(wfs_index > 0) then
      !$acc kernels
      wfs(:,:) = keep_geig(:, :, wfs_index)
      !$acc end kernels
      priority_geig(:) = priority_geig(:) - 1
      priority_geig(wfs_index) = priority_geig(wfs_index) + nkeep_wfs
      success = .true.
      found_geig = found_geig + 1
    else
      missing_geig = missing_geig + 1
    endif
  end function set_geig_from_keep

  subroutine update_keep_cphi(iq, is, wfs)
    integer, intent(in) :: iq, is
    complex(8), intent(in) :: wfs(ndima*nspc, nband)
    integer :: wfs_index, new_index
    integer, allocatable :: min_indices(:)
#ifdef __GPU
    attributes(device) :: wfs
#endif
    if(nkeep_wfs < 1) return
    wfs_index = get_keep_cphi_index(iq, is)
    priority_cphi(:) = priority_cphi(:) - 1
    if(wfs_index > 0) then
      priority_cphi(wfs_index) = priority_cphi(wfs_index) + nkeep_wfs
    else
      min_indices = minloc(priority_cphi)
      new_index = min_indices(1)
      if(new_index < 0 .or. new_index > nkeep_wfs) call rx( 'keep_wfs: bug: out of range of nkeep_wfs ...')
      qs_cphi(new_index,1) = iq
      qs_cphi(new_index,2) = is
      !$acc kernels
      keep_cphi(:,:,new_index) = wfs(:,:)
      !$acc end kernels
      priority_cphi(new_index) = priority_cphi(new_index) + nkeep_wfs
    endif
  end subroutine update_keep_cphi
  subroutine update_keep_geig(iq, is, wfs)
    integer, intent(in) :: iq, is
    complex(8), intent(in) :: wfs(ngpmx*nspc, nband)
    integer :: wfs_index, new_index
    integer, allocatable :: min_indices(:)
#ifdef __GPU
    attributes(device) :: wfs
#endif
    if(nkeep_wfs < 1) return
    wfs_index = get_keep_geig_index(iq, is)
    priority_geig(:) = priority_geig(:) - 1
    if(wfs_index > 0) then
      priority_geig(wfs_index) = priority_geig(wfs_index) + nkeep_wfs
    else
      min_indices = minloc(priority_geig)
      new_index = min_indices(1)
      if(new_index < 0 .or. new_index > nkeep_wfs) call rx( 'keep_wfs: bug: out of range of nkeep_wfs ...')
      qs_geig(new_index,1) = iq
      qs_geig(new_index,2) = is
      !$acc kernels
      keep_geig(:,:,new_index) = wfs(:,:)
      !$acc end kernels
      priority_geig(new_index) = priority_geig(new_index) + nkeep_wfs
    endif
  end subroutine update_keep_geig

  pure integer function get_keep_cphi_index(iq, is) result(wfs_index)
    integer, intent(in) :: iq, is
    integer, allocatable :: indices(:)
    integer :: i
    wfs_index = -1
    indices = pack([(i, i=1, nkeep_wfs)], (qs_cphi(:, 1) == iq) .and. (qs_cphi(:, 2) == is))
    if(size(indices) > 0) wfs_index = indices(1)
  end function get_keep_cphi_index
  pure integer function get_keep_geig_index(iq, is) result(wfs_index)
    integer, intent(in) :: iq, is
    integer, allocatable :: indices(:)
    integer :: i
    wfs_index = -1
    indices = pack([(i, i=1, nkeep_wfs)], (qs_geig(:, 1) == iq) .and. (qs_geig(:, 2) == is))
    if(size(indices) > 0) wfs_index = indices(1)
  end function get_keep_geig_index

end module m_keep_wfs
