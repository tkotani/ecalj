      logical function pars1v(strn,size,nam,term,cast,val)
C- Parses a string for a specified variable
C ----------------------------------------------------------------
Ci Inputs
Ci   recrd(0:*): string recrd is parsed from i to size-1
Ci   iopt:   0: var declaration of existing variable is supressed
Ci           1: var declaration of existing variable supersedes
Ci   mxdecl: parsing ends if number of declarations exceeds mxdecl
Co Outputs
Co   i:      last character parsed
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer size,cast
      character*(*) strn,nam,term*1
      double precision val
C Local parameters
      integer i,j
      logical parstr,a2bin

      pars1v = .false.
      i = 0
      if (.not. parstr(strn,nam,size,len(nam),term,i,j)) return
      if (a2bin(strn,val,cast,0,' ',j,size)) pars1v = .true.
      end

