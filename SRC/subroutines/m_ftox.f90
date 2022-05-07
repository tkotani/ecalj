! General module for formatted write'
! ftox: general format. usage write(6,ftox) foobar,xxx
! ftom(drealx) truncate right-hand-side zeros. drealx can be real(8) arrays.
! ftof(drealx,m) show m digits below period. m is option (six when no m)
!==============
! program test_ftox
!   use m_ftox
!   real(8):: w1,w2,aaa(5)
!   character(12):: fff
!   n=-3132922
!   w1=12245.807807807d0
!   w2=-243.901901901d0
!   aaa=[1.24,4.56,117.1234567,4566666.9,90808.]
!   write(6,ftox)'nnn1',n,'aaa',ftof(w1),'bbb',ftof(w1,5),'bbb',ftod(w2,12)
!   write(6,ftox)'nnn2',ftof(aaa(1:3))
!   write(6,ftox)'nnn3',ftod(aaa,2)
! end program test_ftoX
!=========================
!We have
!>nnn1 -3132922 aaa 12245.807808 bbb 12245.80781 bbb -0.243901901901D+03
!>nnn2 1.240000 4.560000 117.123459
!>nnn3 0.12D+01 0.46D+01 0.12D+03 0.46D+07 0.91D+05
module m_FtoX
  public:: ftof,ftod,ftom !ftom is removing right-end zeros below decimal point.
  character(11),public:: ftox='(*(g0,x))'
  interface ftof
     module procedure ftof,ftofv
  endinterface ftof
  interface ftod
     module procedure ftod,ftodv
  endinterface ftod
  interface ftom
     module procedure ftom,ftomv
  endinterface ftom
  private
contains
  function ftodv(argv,ixx) result(farg)
    intent(in):: argv,ixx
    real(8):: argv(:)
    integer::i,ix
    character(:),allocatable:: farg
    integer,optional:: ixx
    character(1000):: mmm
    ix=6
    if(present(ixx)) ix=ixx
    write(mmm,"(*(g0,x))") (ftod(argv(i),ix),i=1,size(argv))
    mmm=adjustl(mmm)
    if(allocated(farg)) deallocate(farg)
    allocate(farg,source=mmm(1:len(trim(mmm))))
  end function ftodv

  function ftofv(argv,ixx) result(farg)
    intent(in):: argv,ixx
    real(8):: argv(:)
    integer::i,ix
    character(:),allocatable:: farg
    integer,optional:: ixx
    character(1000):: mmm
    ix=6
    if(present(ixx)) ix=ixx
    write(mmm,"(*(g0,x))") (ftof(argv(i),ix),i=1,size(argv))
    mmm=adjustl(mmm)
    if(allocated(farg)) deallocate(farg)
    allocate(farg,source=mmm(1:len(trim(mmm))))
  end function ftofv

  function ftomv(argv,ixx) result(farg)
    intent(in):: argv,ixx
    real(8):: argv(:)
    integer::i,ix
    character(:),allocatable:: farg
    integer,optional:: ixx
    character(1000):: mmm
    ix=6
    if(present(ixx)) ix=ixx
    write(mmm,"(*(g0,x))") (ftom(argv(i),ix),i=1,size(argv))
    mmm=adjustl(mmm)
    if(allocated(farg)) deallocate(farg)
    allocate(farg,source=mmm(1:len(trim(mmm))))
  end function ftomv

  function ftof(arg,ixx) result(farg)
    intent(in)::arg,ixx
    real(8):: arg
    character(:),allocatable:: farg
    character(32):: mmm,fmt
    integer::lsize,ix
    integer,optional:: ixx
    ix=6
    if(present(ixx)) ix=ixx
    write(mmm,"("//"f32."//charnum3(ix)//")") arg
    mmm=adjustl(mmm)
    if(allocated(farg)) deallocate(farg)
    allocate(farg,source=mmm(1:len(trim(mmm))))
  end function ftof
  !
  function ftod(arg,ixx) result(farg)
    intent(in)::arg,ixx
    real(8):: arg
    character(:),allocatable:: farg
    character(32):: mmm
    integer::lsize,ix
    integer,optional:: ixx
    ix=6
    if(present(ixx)) ix=ixx
    write(mmm,"(d32."//charnum3(ix)//")") arg
    mmm  =adjustl(mmm)
    if(allocated(farg)) deallocate(farg)
    allocate(farg,source=mmm(1:len(trim(mmm))))
  end function ftod
  !  
  function ftom(arg,ixx) result(farg) !arg =3.45600000 is '3.45', trucates to rightside zeros'
    intent(in)::arg,ixx
    real(8):: arg
    character(:),allocatable:: farg
    character(32):: mmm,fmt
    integer::lsize,ix,i,j
    integer,optional:: ixx
    write(mmm,"(*(g0,x))")ftof(arg,12)
    j=len(trim(mmm))
    do i=len(mmm),1,-1
       if(mmm(i:i)==' ') cycle
       if(mmm(i:i)=='.') then
          j=i-1
          exit
       elseif(mmm(i:i)/='0') then
          j=i
          exit
       endif
    enddo
    if(allocated(farg)) deallocate(farg)
    allocate(farg,source=mmm(1:j))
  end function ftom

  
  character(3) function charnum3(num)
    integer(4) ::num
    charnum3=char(48+mod(num/100,10))//char(48+mod(num/10,10))//char(48+mod(num,10))
  end function charnum3
end module m_FtoX

! program test_ftox
!   use m_ftox
!   real(8):: w1,w2,aaa(5)
!   character(12):: fff
!   n=-3132922
!   w1=12245.807807807d0
!   w2=-243.901901901d0
!   aaa=[1.24,4.56,117.1234567,4566666.9,90808.]
!   write(6,ftox)'nnn1',n,'aaa',ftof(w1),'bbb',ftof(w1,5),'bbb',ftod(w2,12)
!   write(6,ftox)'nnn2',ftof(aaa(1:3))
!   write(6,ftox)'nnn3',ftod(aaa,2)
! end program test_ftoX
