module m_sort
  implicit none
  public :: sort_index, lower_bound, upper_bound
contains
  function sort_index(array) result(idx)
    implicit none
    real(8), intent(in) :: array(:)
    integer, allocatable :: idx(:)
    integer :: i, n
    n = size(array)
    if(.not. allocated(idx)) allocate(idx(n))
    idx(1:n) = [(i,i=1,n)]
    call quicksort_recursive(array, idx, 1, n)
  end function
  recursive subroutine quicksort_recursive(array, idx, left, right)
    implicit none
    real(8), intent(in) :: array(:)
    integer, intent(inout) :: idx(:)
    integer, intent(in) :: left, right
    integer :: i, j, pivot, temp
    if(left < right) then
      pivot = idx((left + right) / 2)
      i = left
      j = right
      do
        do while (array(idx(i)) < array(pivot))
          i = i + 1
        end do
        do while (array(idx(j)) > array(pivot))
          j = j - 1
        enddo
        if(i <= j) then
          temp = idx(i)
          idx(i) = idx(j)
          idx(j) = temp
          i = i + 1
          j = j - 1
        end if
        if(i > j) exit
      enddo
      call quicksort_recursive(array, idx, left, j)
      call quicksort_recursive(array, idx, i, right)
    end if
  end subroutine

  integer function lower_bound(array, value, idx)
    implicit none
    real(8), intent(in) :: array(:), value
    integer, intent(in), optional, target :: idx(:)
    integer, pointer :: idx_work(:)
    integer :: n, left, right, mid, i
    n = size(array)
    if (present(idx)) then
      idx_work => idx
    else
      allocate(idx_work(n))
      idx_work = [(i, i=1, n)]
    end if
    left = 1; right = n; lower_bound = n + 1
    do while (left <= right)
      mid = (left + right) / 2
      if (array(idx_work(mid)) >= value) then
        lower_bound = mid
        right = mid - 1
      else
        left = mid + 1
      end if
    end do
    if (.not. present(idx)) deallocate(idx_work)
  end function

  integer function upper_bound(array, value, idx)
    implicit none
    real(8), intent(in) :: array(:)
    real(8), intent(in) :: value
    integer, intent(in), optional, target :: idx(:)
    integer, pointer :: idx_work(:)
    integer :: n, left, right, mid, i
    n = size(array)
    if (present(idx)) then
      idx_work => idx
    else
      allocate(idx_work(n))
      idx_work = [(i, i=1, n)]
    end if
    left = 1; right = n; upper_bound = 0
    do while (left <= right)
      mid = (left + right) / 2
      if (array(idx_work(mid)) <= value) then
        upper_bound = mid
        left = mid + 1
      else
        right = mid - 1
      end if
    end do
    if (.not. present(idx)) deallocate(idx_work)
  end function
end module
