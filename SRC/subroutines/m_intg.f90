module m_intg
  !! Numerical integration routines
  implicit none

  interface intg_trapezoidal_nonuniform
     module procedure intg_trapezoidal_nonuniform_d
     module procedure intg_trapezoidal_nonuniform_z
  endinterface

  interface intg_pade_nonuniform
     module procedure intg_pade_nonuniform_d
     module procedure intg_pade_nonuniform_z
  endinterface

  private :: pade_coeff_3pt_d, pade_coeff_3pt_z, &
             intg_pade_segment_d, intg_pade_segment_z
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

  !> [1/1] Padé-based integration for peaked functions on non-uniform mesh (Real version)
  !! Fits a local rational interpolant f(x)=(Ax+B)/(x+C) through 3 consecutive mesh points
  !! via Thiele's reciprocal differences, then integrates analytically over each sub-interval.
  !! Falls back to trapezoidal if the Padé fit is degenerate or has a pole in the interval.
  real(8) function intg_pade_nonuniform_d(x, y) result(res)
    real(8), intent(in) :: x(:), y(:)
    integer :: i, n, i0, i1, i2
    real(8) :: A, B, C
    logical :: fb
    n = min(size(x), size(y))
    res = 0d0
    if(n < 2) return
    if(n < 3) then
      res = 0.5d0*(y(1)+y(2))*(x(2)-x(1)); return
    endif
    do i = 1, n-1
      if(i < n-1) then
        i0 = i; i1 = i+1; i2 = i+2   ! forward-looking 3-point window
      else
        i0 = i-1; i1 = i; i2 = i+1   ! backward-looking for last interval
      endif
      call pade_coeff_3pt_d(x(i0),y(i0), x(i1),y(i1), x(i2),y(i2), A, B, C, fb)
      if(fb) then
        res = res + 0.5d0*(y(i)+y(i+1))*(x(i+1)-x(i))
      else
        res = res + intg_pade_segment_d(x(i), x(i+1), A, B, C)
      endif
    enddo
  end function intg_pade_nonuniform_d

  !> [1/1] Padé-based integration for peaked functions on non-uniform mesh (Complex version)
  complex(8) function intg_pade_nonuniform_z(x, y) result(res)
    real(8),    intent(in) :: x(:)
    complex(8), intent(in) :: y(:)
    integer :: i, n, i0, i1, i2
    complex(8) :: A, B, C
    logical :: fb
    n = min(size(x), size(y))
    res = (0d0, 0d0)
    if(n < 2) return
    if(n < 3) then
      res = 0.5d0*(y(1)+y(2))*(x(2)-x(1)); return
    endif
    do i = 1, n-1
      if(i < n-1) then
        i0 = i; i1 = i+1; i2 = i+2
      else
        i0 = i-1; i1 = i; i2 = i+1
      endif
      call pade_coeff_3pt_z(x(i0),y(i0), x(i1),y(i1), x(i2),y(i2), A, B, C, fb)
      if(fb) then
        res = res + 0.5d0*(y(i)+y(i+1))*(x(i+1)-x(i))
      else
        res = res + intg_pade_segment_z(x(i), x(i+1), A, B, C)
      endif
    enddo
  end function intg_pade_nonuniform_z

  !> Fit [1/1] Padé f(x)=(A*x+B)/(x+C) to 3 points via Thiele's reciprocal differences (Real)
  !! r01=(x0-x1)/(y0-y1), r12=(x1-x2)/(y1-y2), d2=(x0-x2)/(r01-r12)
  !! A=y0+d2, B=y0*r01*d2-y0*x1-d2*x0, C=r01*d2-x1
  subroutine pade_coeff_3pt_d(x0, y0, x1, y1, x2, y2, A, B, C, fallback)
    real(8), intent(in)  :: x0, y0, x1, y1, x2, y2
    real(8), intent(out) :: A, B, C
    logical, intent(out) :: fallback
    real(8) :: r01, r12, d2, sc
    real(8), parameter :: eps = 1d-10
    fallback = .false.
    sc = max(abs(y0), abs(y1), abs(y2), 1d-100)
    if(abs(y0-y1) < eps*sc .or. abs(y1-y2) < eps*sc) then
      fallback = .true.; return
    endif
    r01 = (x0-x1)/(y0-y1)
    r12 = (x1-x2)/(y1-y2)
    if(abs(r01-r12) < eps*max(abs(r01), abs(r12), 1d-100)) then
      fallback = .true.; return
    endif
    d2 = (x0-x2)/(r01-r12)
    A  = y0 + d2
    B  = y0*r01*d2 - y0*x1 - d2*x0
    C  = r01*d2 - x1
  end subroutine pade_coeff_3pt_d

  !> Fit [1/1] Padé f(x)=(A*x+B)/(x+C) to 3 points via Thiele's reciprocal differences (Complex y, real x)
  subroutine pade_coeff_3pt_z(x0, y0, x1, y1, x2, y2, A, B, C, fallback)
    real(8),    intent(in)  :: x0, x1, x2
    complex(8), intent(in)  :: y0, y1, y2
    complex(8), intent(out) :: A, B, C
    logical,    intent(out) :: fallback
    complex(8) :: r01, r12, d2
    real(8) :: sc
    real(8), parameter :: eps = 1d-10
    fallback = .false.
    sc = max(abs(y0), abs(y1), abs(y2), 1d-100)
    if(abs(y0-y1) < eps*sc .or. abs(y1-y2) < eps*sc) then
      fallback = .true.; return
    endif
    r01 = (x0-x1)/(y0-y1)
    r12 = (x1-x2)/(y1-y2)
    if(abs(r01-r12) < eps*max(abs(r01), abs(r12), 1d-100)) then
      fallback = .true.; return
    endif
    d2 = (x0-x2)/(r01-r12)
    A  = y0 + d2
    B  = y0*r01*d2 - y0*x1 - d2*x0
    C  = r01*d2 - x1
  end subroutine pade_coeff_3pt_z

  !> Analytic integral of (A*x+B)/(x+C) from xa to xb (Real)
  !! = A*(xb-xa) + (B-A*C)*ln|(xb+C)/(xa+C)|
  !! Falls back to trapezoidal if pole (x=-C) is inside the interval.
  real(8) function intg_pade_segment_d(xa, xb, A, B, C) result(res)
    real(8), intent(in) :: xa, xb, A, B, C
    if((xa+C)*(xb+C) < 0d0) then
      ! Pole in interval: trapezoidal fallback
      res = 0.5d0*((A*xa+B)/(xa+C) + (A*xb+B)/(xb+C))*(xb-xa)
    else
      res = A*(xb-xa) + (B-A*C)*log(abs((xb+C)/(xa+C)))
    endif
  end function intg_pade_segment_d

  !> Analytic integral of (A*x+B)/(x+C) from xa to xb (Complex A,B,C; real xa,xb)
  !! Uses complex logarithm: log is well-defined when Im(C) /= 0 (pole off real axis).
  !! Falls back to trapezoidal if C is nearly real and pole is inside the interval.
  complex(8) function intg_pade_segment_z(xa, xb, A, B, C) result(res)
    real(8),    intent(in) :: xa, xb
    complex(8), intent(in) :: A, B, C
    real(8), parameter :: eps = 1d-14
    if(abs(dimag(C)) < eps .and. (xa+dble(C))*(xb+dble(C)) < 0d0) then
      res = 0.5d0*((A*xa+B)/(xa+C) + (A*xb+B)/(xb+C))*(xb-xa)
    else
      res = A*(xb-xa) + (B-A*C)*log((xb+C)/(xa+C))
    endif
  end function intg_pade_segment_z

end module m_intg
