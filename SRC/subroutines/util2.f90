module util2
  public:: ftoa6
contains
  character(24) function ftoa6(arg)
    real(8):: arg
    write(ftoa6,"(f24.6)") arg
    ftoa6=adjustl(ftoa6)
  end function ftoa6
end module util2

