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
    use m_cmdopt_registry, only: load_cmdopt2_registry
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
    ! Parse cmdopt2 entries (--jobgw=, --nb=, ...) into typed module
    ! variables in m_cmdopt_registry. Idempotent; recursive cmdopt2
    ! calls re-enter m_setargs and bail out at the `init` check above.
    call load_cmdopt2_registry()
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
    ! Parse cmdopt2 entries into typed module variables (same as m_setargs).
    block
      use m_cmdopt_registry, only: load_cmdopt2_registry
      call load_cmdopt2_registry()
    end block
  end subroutine m_setargsc
end module m_args

module m_ext
  use m_args,only: m_setargs,arglist,narg
  character(512),public,protected::sname='tempext',dirname
  public:: m_ext_init, print_usage_and_quit
contains
  subroutine m_ext_init() bind(C)
    logical :: master
    integer:: ifi,ipos,i,na,getcwd,ios,nmatch,dotpos,getpid
    character*256:: sss,s222,argv,line,fname,candidate,tmpfile
    character*32 :: pidstr
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
    ! No positional arg supplied. GW-side programs (qg4gw, heftet, hbasfp0,
    ! hvccfp0, hx0fp0, hwmatK_MPI, ...) don't take sname on the command line.
    ! Auto-detect: if cwd has exactly one ctrlg.<x>.toml, use that <x>.
    ! Per-rank tmp file (PID-suffixed) to avoid races across MPI ranks.
    write(pidstr,'(i0)') getpid()
    tmpfile = '.ext_glob_'//trim(pidstr)
    call execute_command_line('ls -1 ctrlg.*.toml 2>/dev/null > '//trim(tmpfile), wait=.true.)
    open(newunit=ifi, file=trim(tmpfile), status='old', action='read', iostat=ios)
    if (ios == 0) then
       nmatch = 0
       candidate = ''
       do
          read(ifi, '(a)', iostat=ios) line
          if (ios /= 0) exit
          if (len_trim(line) == 0) cycle
          nmatch = nmatch + 1
          candidate = trim(line)
       enddo
       close(ifi, status='delete')
       if (nmatch == 1) then
          fname = trim(candidate)
          dotpos = index(fname, '.toml', back=.true.)
          if (fname(1:6) == 'ctrlg.' .and. dotpos > 7) then
             sname = fname(7:dotpos-1)   ! chars between "ctrlg." and ".toml"
             goto 999
          endif
       endif
    endif
    ! Couldn't auto-detect; leave sname at default ('tempext').
    ! Callers that need a real sname must check and abort themselves.
    return
999 continue
  end subroutine m_ext_init

  !> Print a short usage banner and exit. Called from main_lmf / main_lmfa /
  !  main_lmchk when --help is passed; the full option catalogue lives in
  !  ecaljdoc (manual/cmdopts) because the source has ~135 cmdopt sites
  !  scattered across SRC/, and keeping a copy in sync here is hopeless.
  subroutine print_usage_and_quit(prgnam)
    character(*), intent(in) :: prgnam
    write(6,'(a)') ''
    write(6,'(a)') 'Usage: '//trim(prgnam)//' <sname> [--option ...] [--ctrlg:<path>=<value> ...]'
    write(6,'(a)') ''
    write(6,'(a)') '  <sname> is the extension of the control file (ctrlg.<sname>.toml).'
    write(6,'(a)') '  GW-side utilities (qg4gw, heftet, hbasfp0, hvccfp0, hx0fp0, ...) take no'
    write(6,'(a)') '  positional <sname>; they auto-detect from the unique ctrlg.*.toml in cwd.'
    write(6,'(a)') ''
    write(6,'(a)') 'Inputs (lmf / lmfa / lmchk, only files actually read):'
    write(6,'(a)') '  ctrlg.<sname>.toml          main control file (TOML schema)'
    write(6,'(a)') '  PB.<sname>.toml             product-basis table (GW path only)'
    write(6,'(a)') '  syml.<sname>                k-line for --band'
    write(6,'(a)') '  atmpnu.{1,2,3}.<sname>      atomic radial wfns (from lmfa)'
    write(6,'(a)') '  rst.<sname>                 density restart (carries from previous SCF)'
    write(6,'(a)') '  sigm.<sname>                self-energy (QSGW)'
    write(6,'(a)') ''
    write(6,'(a)') 'Run-time TOML override:'
    write(6,'(a)') '  --ctrlg:<dotted.path>=<value>   override a key in ctrlg.<sname>.toml.'
    write(6,'(a)') '    Examples: --ctrlg:verbose=50            --ctrlg:time=[5,5]'
    write(6,'(a)') '              --ctrlg:bz.nkabc=[8,8,8]      --ctrlg:ham.scaledsigma=0.8'
    write(6,'(a)') '              --ctrlg:ham.so=1              --ctrlg:ham.nspin=2'
    write(6,'(a)') '              --ctrlg:ham.phispinsym=true   --ctrlg:spec.1.r=2.5'
    write(6,'(a)') '    Values are TOML-typed: bool lowercase (true/false), strings quoted,'
    write(6,'(a)') '    arrays in [...]. Applied in memory; the file on disk is untouched.'
    write(6,'(a)') '    Each applied override is logged on rank 0.'
    write(6,'(a)') ''
    write(6,'(a)') '    Retired (now abort): -v..., --[<path>]=..., --<a.b>=...,'
    write(6,'(a)') '                         --pr=N, --time=..., --phispinsym.'
    write(6,'(a)') '                         Use --ctrlg:<path>=<value>.'
    write(6,'(a)') ''
    write(6,'(a)') 'Full documentation:'
    write(6,'(a)') '  manual:           https://ecalj.github.io/ecaljdoc/manual/lmf'
    write(6,'(a)') '  --foo catalogue:  https://ecalj.github.io/ecaljdoc/manual/cmdopts'
    write(6,'(a)') '  TOML override:    https://ecalj.github.io/ecaljdoc/manual/toml_migration'
    write(6,'(a)') ''
    write(6,'(a)') 'Most-asked subset:'
    write(6,'(a)') '  --help                this banner'
    write(6,'(a)') '  --ctrlg:verbose=N      console verbosity (e.g. =50 traces, =70 debug)'
    write(6,'(a)') '  --ctrlg:time=[N,M]     CPU timing log: depth, on-the-fly (e.g. [5,5])'
    write(6,'(a)') '  --band                band plot along syml.<sname>'
    write(6,'(a)') '  --jobgw={0,1}         run as GW driver (replaces lmfgw-MPIK)'
    write(6,'(a)') '  --quit={show,ham,mkpot,dmat,band}   staged stop points'
    write(6,'(a)') '  --writeham            emit HamiltonianPMT.* for downstream tools'
    write(6,'(a)') '  --mkprocar --fullmesh PROCAR / Fermi-surface mesh'
    write(6,'(a)') '  --gpu                 GPU path (grabs /tmp/gpu.lock)'
    write(6,'(a)') ''
    flush(6)
    call exit(0)   ! bypass rx0 because MPI may not be initialized yet
  end subroutine print_usage_and_quit
end module m_ext
logical function cmdopt0(argstr)! Check a command-line argument exist. 
  use m_args,only: m_setargs,arglist,narg
  !i Inputs  argstr: command-line string to search; search to strln chars
  !o Outputs cmdopt: T if argument found, else F
  implicit none
  character(*):: argstr
  integer :: iarg
  cmdopt0 = .false.
  call m_setargs()
  do iarg=1,narg
     if(trim(arglist(iarg)) == trim(argstr)) then
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
