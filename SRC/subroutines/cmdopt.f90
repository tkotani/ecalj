logical function cmdopt0(argstr)
  !! new version of cmdopt
  !! Check a command-line argument exist. Then return it in outstr
  ! ----------------------------------------------------------------
  !i Inputs
  !i   argstr: command-line string to search; search to strln chars
  !o Outputs
  !o   cmdopt: T if argument found, else F
  !o   outstr: output string
  ! ----------------------------------------------------------------
  !     implicit none
  character(*):: argstr
  integer ::     nargs,strln !dummy
  ! Local parameters
  logical :: lsequ
  integer :: iarg,nargf,idum,nxarg,strlnx
  character(120) :: strn
  cmdopt0 = .false.
  do iarg=1,iargc()
     call getarg(iarg,strn)
     strlnx = len_trim(argstr) !override input strln
     if(strn(1:strlnx)==trim(argstr)) then
        cmdopt0 = .true.
        return
     endif
  enddo
end function cmdopt0

logical function cmdopt2(argstr,outstr)
  !! new version of cmdopt
  !! Check a command-line argument exist. Then return it in outstr
  ! ----------------------------------------------------------------
  !i Inputs
  !i   argstr: command-line string to search; search to strln chars
  !o Outputs
  !o   cmdopt: T if argument found, else F
  !o   outstr: output string
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  character(*):: argstr
  character(*):: outstr
  integer ::      nargs,strln !dummy
  ! Local parameters
  logical :: lsequ
  integer :: iarg,nargf,idum,nxarg,strlnx
  character(120) :: strn
  cmdopt2 = .false.
  do iarg=1,iargc()
     call getarg(iarg,strn)
     strlnx = len_trim(argstr) !override input strln
     if(strn(1:strlnx)==trim(argstr)) then
        cmdopt2 = .true.
        outstr = strn(strlnx+1:)
        return
     endif
  enddo
end function cmdopt2

logical function cmdopt(argstr,strln,nargs,outstr)
  !!                                   dummy,dummy
  !! Check a command-line argument exist. Then return it in outstr
  ! ----------------------------------------------------------------
  !i Inputs
  !i   argstr: command-line string to search; search to strln chars
  !o Outputs
  !o   cmdopt: T if argument found, else F
  !o   outstr: output string
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  character(*):: argstr
  character(*):: outstr
  integer ::     nargs,strln !dummy
  ! Local parameters
  logical :: lsequ
  integer :: iarg,nargf,idum,nxarg,strlnx
  character(120) :: strn
  cmdopt = .false.
  do iarg=1,iargc()
     call getarg(iarg,strn)
     strlnx = len_trim(argstr) !override input strln
     if(trim(strn(1:strlnx))==trim(argstr)) then
        cmdopt = .true.
        outstr = strn
        return
     endif
  enddo
end function cmdopt
