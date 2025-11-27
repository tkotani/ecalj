module m_env
  use iso_c_binding
  implicit none
  private
  public :: setenv, unsetenv
  interface
    function setenv_c(name, value, overwrite) bind(C, name="setenv") result(status)
      use iso_c_binding
      character(kind=c_char), intent(in) :: name(*), value(*)
      integer(c_int), value :: overwrite
      integer(c_int) :: status
    end function setenv_c
    function unsetenv_c(name) bind(C, name="unsetenv") result(status)
      use iso_c_binding
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: status
    end function unsetenv_c
  end interface
contains
  function setenv(name, value, overwrite) result(status)
    character(len=*), intent(in) :: name, value
    logical, intent(in), optional :: overwrite
    integer :: status
    logical :: ow
    ow = .true.
    if (present(overwrite)) ow = overwrite
    status = int(setenv_c(trim(name)//c_null_char, trim(value)//c_null_char, merge(1, 0, ow)))
  end function setenv
  function unsetenv(name) result(status)
    character(len=*), intent(in) :: name
    integer :: status
    status = int(unsetenv_c(trim(name)//c_null_char))
  end function unsetenv
end module m_env
