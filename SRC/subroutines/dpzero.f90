subroutine cinit(array,leng)
  !- Initializes complex array to zero
  integer :: leng
  real :: array(2*leng)
  call ainit(array,leng+leng)
end subroutine cinit
subroutine zinit(array,leng)
  !- Initializes complex array to zero
  integer :: leng
  double precision :: array(2*leng)
  call dpzero(array,leng+leng)
end subroutine zinit
subroutine iinit(array,leng)
  !- Initializes integer array to zero
  integer :: leng,i
  integer :: array(leng)
  do   i=1, leng
     array(i) = 0
  enddo
end subroutine iinit
subroutine ainit(array,leng)
  !- Initializes real array to zero
  integer :: leng,i
  real :: array(leng)
  do   i=1, leng
     array(i) = 0
  enddo
end subroutine ainit
subroutine dpzero(array,leng)
  !- Initializes double precision array to zero
  integer :: leng,i
  double precision :: array(leng)
  do   i=1, leng
     array(i) = 0
  enddo
end subroutine dpzero
real function rval(array,index)
  integer :: index
  !- Returns the real value of ARRAY(INDEX)
  real :: array(index)
  rval = array(index)
end function rval
real(8) function dval(array,index)
  integer :: index
  !- Returns the double precision value of ARRAY(INDEX)
  double precision :: array(index)
  dval = array(index)
END function dval
integer function ival(array,index)
  !- Returns the integer value of ARRAY(INDEX)
  integer :: index
  integer :: array(index)
  ival = array(index)
end function ival
integer function ival2(array,nda,i1,i2)
  !- Returns the integer value of ARRAY(i1,i2)
  !     implicit none
  integer :: nda,i1,i2,array(nda,1)
  ival2 = array(i1,i2)
end function ival2
logical function logval(array,index)
  !- Returns the integer value of ARRAY(INDEX)
  integer :: index
  logical :: array(index)
  logval = array(index)
end function logval
complex function cval(array,index)
  !- Returns the complex value of ARRAY(INDEX)
  integer :: index
  complex :: array(index)
  cval = array(index)
end function cval
subroutine dvset(array,i1,i2,val)
  !- Sets some elements of double precision array to value
  integer :: i1,i2
  double precision :: array(i2),val
  integer :: i
  do   i = i1, i2
     array(i) = val
  enddo
end subroutine dvset
subroutine ivset(array,i1,i2,val)
  !- Sets some elements of integer array to value
  integer :: i1,i2,array(1),val,i
  do   i = i1, i2
     array(i) = val
  enddo
end subroutine ivset
subroutine lvset(array,i1,i2,val)
  !- Sets some elements of logical array to value
  integer :: i1,i2,i
  logical :: array(1),val
  do  i = i1, i2
     array(i) = val
  enddo
end subroutine lvset
!$$$      subroutine redfrr(oname,leng)
!$$$C- Release to pointer oname, reallocate double oname of length leng
!$$$C     implicit none
!$$$      integer oname,leng
!$$$      call rlse(oname)
!$$$      call defrr(oname,leng)
!$$$      end
!$$$      subroutine redfi(oname,leng)
!$$$C- Release to pointer oname, reallocate double oname of length leng
!$$$C     implicit none
!$$$      integer oname,leng
!$$$      call rlse(oname)
!$$$      call defi(oname,leng)
!$$$      end

