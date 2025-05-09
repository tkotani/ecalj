!> Some utilities. Get arglist. 
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
  subroutine m_setargsc(cname,prt) bind(C) !Pass narg and arglist from python instead of m_setargs
    implicit none
    integer:: i,n
    logical:: prt
    character(1):: cname(*)
    character(1024)::convcchar,string
    character(120),allocatable::arglist2(:)
    argall= trim(convcchar(cname))
    string=argall
    argall=' '//argall
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
    if(prt) then
    do i=1,narg
       write(*,*)'m_setargsc=',i, trim(arglist(i))
    enddo
    endif
  end subroutine m_setargsc
end module m_args

module m_ext
  use m_args,only: m_setargs,arglist,narg
  character(512),public,protected::sname='tempext',dirname
  public:: m_ext_init
contains
  subroutine m_ext_init() bind(C)
    logical :: master
    integer:: ifi,ipos,i,na,getcwd
    character*256:: sss,s222,argv
    i = getcwd(dirname)
    do i = 1, narg
       !write(*,*)'m_ext_init=',i, trim(arglist(i))
       if(arglist(i)(1:5)=='ctrl.') then
          sname=trim(arglist(i)(6:))
          goto 999
       endif
       if(arglist(i)(1:1)/='-') then
          sname=trim(arglist(i))
          goto 999
       endif
    enddo
    write(6,"( &
         /'Usage: lmf,lmfa,lmchk [--OPTION] [-vfoobar] [extension]'&
         /' Some options:'&
         /'  --help',      t17,'Show this document'&
         /'  --pr=#1',     t17,'Set the verbosity (stack) to values #1' &
         /'  -vfoobar=expr',  t17,'Define numerical variable foobar' &
         /'  --time=#1,#2',t17,'Print timing info to # levels (#1=summary; #2=on-the-fly,e.g. --time=5,5)' &
         /'  --jobgw=0 or --jobgw=1  lmf-MPIK works as the GW driver (previous lmfgw-MPIK)' &
         /'  --quit=band, --quit=mkpot or --quit=dmat: Stop points. Surpress writing rst' &
         /'  NOTE: For description of ctrl file, see ecalj/Document/help_lmf.org and so on!' &
         /'  NOTE: lmf read rst.* prior to atm.* file (Removed --rs options at 2022-6-20)' &
         /'  NOTE: Other command-line-options => Search call cmdopt in SRC/*/*.f90'  )")
    call rx0('no args: Need to set foobar for ctrl.foobar ')
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

! module m_prgnam
!    character(32):: prgnamx = ''
! contains
!    subroutine set_prgnam(prgnam)
!       character(*):: prgnam
!       prgnamx = prgnam
!    end subroutine set_prgnam
!    subroutine set_prgnamc(prgnamc) bind(C)
!       character(1024):: convcchar
!       character(1):: prgnamc(*)
!       prgnamx = trim(convcchar(prgnamc))
!       write (*, *) 'prgnamx=', trim(prgnamx)
!     end subroutine set_prgnamc
! end module m_prgnam

function convcchar(instr) result(outstr) !convert char(1) to char(1024)
   use iso_c_binding
   integer:: i,nend
   character(1):: instr(1024)
   character(1024):: outstr,instr2
   forall(i=1:1024) instr2(i:i) = instr(i)
   nend = index(instr2, c_null_char)-1
   outstr=''
   outstr(1:nend)= instr2(1:nend)
end function convcchar
