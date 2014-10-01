      subroutine getarf(iarg,strn)
C- Returns a command-line argument
C ----------------------------------------------------------------------
Ci Inputs
Ci   iarg   :index to command-line argument (0 for executable)
Co Outputs
Co   strn   :string containing iargth command-line argument
Cr Remarks
Cr   This routine returns the number of command-line arguments
Cr   (a fortran-callable version of getarg).  Fortran has no
Cr   standard for this function, so it is machine-dependent. 
Cr   This routine and function nargf below are interfaces that
Cr   confine the machine dependence to this file.
Cr
Cr   It may be that the compiler supplies this function in a system call
Cr   (e.g. getarg).  If this approach works, it is the preferred method.
Cr
Cr   A possible alternative is available if the linker passes the
Cr   information through to the entry point (main or similarly named
Cr   function; see fmain.c), it can be extracted by calling
Cr   gtargc(iarg,strn) defined in fmain.c.
Cr
Cr   Even if fortran-supplied system calls are available, it can happen
Cr   that the call fails unless the main entry point is a fortran
Cr   routine The SGI and SUN system call work properly only when the
Cr   main program is fortran.  In that case gtargc() works.
Cr
Cr   Compilers running under linux don't start with C main;
Cr   however the system call often works and it can be used instead.
Cu Updates
Cu   17 Jul 01 Written to accomodate problems with HP f90 compiler
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer iarg
      character*(*) strn
#if HASGTARGC
      call gtargc(iarg,strn)
#elif HASGETARG
      call getarg(iarg,strn)
#else HASGETCOMMANDARGUMENT
c -------  fortran2003std
      call get_command_argument(iarg,strn)
C#else
C#error "need a subroutine to get command line argument"
#endif

      end
      integer function nargf()
C- Returns number of command-line argument, including command name
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Co   iargc  :index of the last command line argument.
Cr Remarks
Cr   This routine returns the number of command-line arguments,
Cr   and is a companion to getarf above.
Cr
Cr   This routine returns the number of command-line arguments
Cr   (a fortran-callable version of getarg).  Fortran has no
Cr   standard for this function, so it is machine-dependent. 
Cr   This routine and function getarf above are interfaces that
Cr   confine the machine dependence to this file.
Cr
Cr   It may be that the compiler supplies this function in a system
Cr   call.  If this approach works, it is the preferred method.
Cr
Cr   A possible alternative is available if the linker passes the
Cr   information through to the entry point (main or similarly named
Cr   function; see fmain.c), it can be extracted by calling nargc()
Cr   defined in fmain.c.
Cr
Cr   Even if fortran-supplied system calls are available, it can happen
Cr   that the call fails unless the main entry point is a fortran
Cr   routine The SGI and SUN system call work properly only when the
Cr   main program is fortran.  In that case nargc() in fmain.c works.
Cu Updates
Cu    3 Aug 04 First written
C ----------------------------------------------------------------------
C     implicit none

C ... This uses nargc defined in fmain.c, and requires that the linker
C     passes command-line arguments to the C main entry point.
#if HASNARGC
      integer nargc
      nargf = nargc()
C
C ... This works if a fortran-callable substitute for C iargc available
C     (e.g. DEC fort)
#elif HASIARGC
      integer iargc
c      external iargc
      nargf = iargc() + 1

C ... Another possible fortran-callable substitute for C iargc 
C     (e.g. Intel ifort)
#elif HASNARGS
      integer nargs
c      external nargs
      nargf = nargs()

#else HASCOMMANDARGUMENTCOUNT
c ---------  fortran2003std
      integer:: command_argument_count
      nargf = command_argument_count()+1
C#else
C#error "need a subroutine to get command line argument count."
#endif

      end

      subroutine ftime(datim)
C- Returns date and time
C ----------------------------------------------------------------------
Ci Inputs
Co Outputs
Co   datim  : string containing date and time.  
Co          : On unix systems, this is a 26-character string.
Cr Remarks
Cr   Fortran has no standard for this function, so it is
Cr   machine-dependent.  It may be that the compiler supplies this
Cr   function in a system call.  If this approach works, it is the
Cr   preferred method.
Cr   
Cr   A possible alternative is to invoke 'cftime' in fsubs.c which can
Cr   be invoked from fortran.  If neither of these works, user may wish
Cr   to find an alternative.
Cu Updates
Cu    5 Aug 04 First written
C ----------------------------------------------------------------------
C     implicit none
      character datim*(*)

C ... This uses FCTIME defined in fsubs.c
#if FCTIME
      call fctime(datim)
C
C ... fortran-callable date and time
#elif FDATE
      call fdate(datim)
      datim=datim(1:24) !takao. If this is not, write(6,*) gives CR at the ene of datim*26.
#elif DOL_FDATA
      call $fdate(datim)
#else 
      datim = ' '
#endif

      end
#if TESTF
      program test
      integer i
      character strn*50
      i = nargf()
      print *, 'no command line arguments + 1 = ', i
      call getarf(1,strn)
      print *, '1st command line argument = ', strn
      end
#endif
#if TESTC
      subroutine fmain
      integer i
      character strn*50
      i = nargf()
      print *, 'no command line arguments + 1 = ', i
      call getarf(1,strn)
      print *, '1st command line argument = ', strn
      end
#endif

