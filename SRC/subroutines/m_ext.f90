!> Get arglist. Get extension of ctrl.foobar
module m_args
  character(120),public,protected,allocatable::arglist(:)
  character(1024),public,protected:: argall
  public:: m_setargs,m_setargsc
  integer:: narg
  logical:: init=.true.
  integer,private:: nx=64
contains
  subroutine m_setargs()
    integer:: iarg,iargc
    character(120) :: strn
    if(.not.init) return
    if(allocated(arglist)) return !If arglist is already allocated, do nothing
    narg = iargc()
    allocate(arglist(narg))
    argall=''
    do iarg=1,iargc()
       call getarg(iarg,strn)
       arglist(iarg)=trim(strn)
       argall=trim(argall)//' '//trim(strn)
    enddo
    init=.false.
  endsubroutine m_setargs
  subroutine m_setargsc(cname,nsize) bind(C) !Pass narg and arglist from python instead of m_setargs
    implicit none
    integer:: nsize,i,n
    character(1):: cname(nsize)
    character(512):: string=''
    character(120),allocatable::arglist2(:)
    forall(i=1:nsize) string(i:i)=cname(i)
    string=adjustl(string)
    argall=string
    if(allocated(arglist2)) deallocate(arglist2)
    allocate(arglist2(nx))
    narg=0
    do
       n=index(string,' ')
       if(n==1) exit
       narg=narg+1
       arglist2(narg)=trim(string(:n))
       string=adjustl(string(n+1:))
    enddo
    if(allocated(arglist)) deallocate(arglist)
    allocate(arglist(1:narg))
    call move_alloc(from=arglist2,to=arglist)
    do i=1,narg
       write(*,*)'m_setargsc=',i, trim(arglist(i))
    enddo   
  end subroutine m_setargsc
end module m_args

module m_ext
  use m_args,only: m_setargs,arglist,narg
  character(512),public,protected::sname='temp'
  public:: m_ext_init
contains
  subroutine m_ext_init()
    logical :: master
    integer:: ifi,ipos,i,na
    character*256:: sss,s222,argv
    call m_setargs()
    do i = 1, narg
      !write(*,*)'m_ext_init=',i, trim(arglist(i))
       if(arglist(i)(1:1)/='-') then
          sname=trim(arglist(i))
          goto 999
       endif
    enddo
    call rx0('no args: Need to set foobar for ctrl.foobar')
999 continue
  end subroutine m_ext_init
end module m_ext
logical function cmdopt0(argstr)! Check a command-line argument exist. 
  use m_args,only: m_setargs,arglist,narg
  !i Inputs  argstr: command-line string to search; search to strln chars
  !o Outputs cmdopt: T if argument found, else F
  implicit none
  character(*):: argstr
  integer ::     nargs,strln !dummy
  logical :: lsequ
  integer :: iarg,nargf,idum,nxarg,strlnx
  character(120) :: strn
  cmdopt0 = .false.
  call m_setargs()
  do iarg=1,narg
     strlnx = len_trim(argstr) !override input strln
     if(arglist(iarg)(1:strlnx)==trim(argstr)) then
        cmdopt0 = .true.
        return
     endif
  enddo
end function cmdopt0
logical function cmdopt2(argstr,outstr)  ! return it in outstr
  use m_args,only: m_setargs,arglist,narg
  !i Inputs argstr: command-line string to search; search to strln chars
  !o Outputs cmdopt: T if argument found, else F
  !o   outstr: output string
  implicit none
  character(*):: argstr,outstr
  integer ::     nargs,strln !dummy
  logical :: lsequ
  integer :: iarg,nargf,idum,nxarg,strlnx
  character(120) :: strn
  cmdopt2 = .false.
    call m_setargs()
  do iarg=1,narg
     strlnx = len_trim(argstr) !override input strln
     if(arglist(iarg)(1:strlnx)==trim(argstr)) then
        cmdopt2 = .true.
        outstr = arglist(iarg)(strlnx+1:)
        return
     endif
  enddo
end function cmdopt2
