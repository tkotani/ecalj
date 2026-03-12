module m_intg
  !! Numerical integration routines
  implicit none

  interface intg_trapezoidal_nonuniform
     module procedure intg_trapezoidal_nonuniform_d
     module procedure intg_trapezoidal_nonuniform_z
  endinterface
contains

  !> Trapezoidal rule for non-uniform mesh (Real version)
  pure real(8) function intg_trapezoidal_nonuniform_d(x, y) result(res)
    real(8), intent(in) :: x(:), y(:)
    integer :: i, n
    n = min(size(x), size(y))
    res = 0.0d0
    if (n < 2) return
    do i = 1, n - 1
       res = res + 0.5d0 * (y(i) + y(i+1)) * (x(i+1) - x(i))
    enddo
  end function intg_trapezoidal_nonuniform_d

  !> Trapezoidal rule for non-uniform mesh (Complex version)
  pure complex(8) function intg_trapezoidal_nonuniform_z(x, y) result(res)
    real(8),    intent(in) :: x(:)
    complex(8), intent(in) :: y(:)
    integer :: i, n
    n = min(size(x), size(y))
    res = (0.0d0, 0.0d0)
    if (n < 2) return
    do i = 1, n - 1
       res = res + 0.5d0 * (y(i) + y(i+1)) * (x(i+1) - x(i))
    enddo
  end function intg_trapezoidal_nonuniform_z

end module m_intg
