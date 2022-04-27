module m_util2
  public:: ftoa6
contains
  function ftoa6(arg) result(farg)
    character(:),allocatable:: farg
    character(32):: mmm
    integer::lsize
    real(8):: arg
    write(mmm,"(f32.6)") arg
    mmm  =adjustl(mmm)
    lsize=len(trim(mmm))
    allocate(farg,source=mmm(1:lsize))
  end function ftoa6
end module m_util2
! module m_util2
!   public:: ftoa6
! contains
!   character(24) function ftoa6(arg)
!     real(8):: arg
!     write(ftoa6,"(f24.6)") arg
!     ftoa6=adjustl(ftoa6)
!   end function ftoa6
! end module m_util2

