      integer function garg(tok,iarg,cast,sep,itrm,narg,nread,it,res)
C- Returns vector numerical vals for a specified command-line argument
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iarg: index to command-line argument
Ci         Use <=0 to call nxtarg
Ci   cast: 0=logical, 1=char, 2=int, 3=real, 4=double
Ci   sep:  string of chars that separate and/or terminate expressions.
Ci   itrm: index to which of the separators in 'sep' to be treated as
Ci         terminators.  garg stops looking for more arguments whenever
Ci         any sep(i>=itrm) is encountered.
Ci   narg: number of args to convert.  If NO args to convert,
Ci         garg still distinguishes appropriate terminator.
Co Outputs:
Co   res:  Vector of numbers that were converted
Co   it:   vector of indices, one for each entry in res, labelelling
Co         which char in 'sep' terminated the expr. for that number
Co   garg,nread: 0 if tok does not match with command-line arg.
Co               1 if tok does match, AND narg=0 
Co               n if token match and converted sans error narg numbers
Co              -n if error on conversion of argument n
Co            -999 if iarg missing from command line
Cr Remarks
Cr   Intended to parse command-line arguments like -x=1,3 or 1:3#4.
Cr   The command-line argument is compared to tok (unless tok is blank).
Cr   There is further parsing for the following numerical arguments 
Cr   when the command-line argument starts with tok.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer iarg,ip,cast,narg,itrm,it(1),nread
      character*(*) tok,sep
      double precision res(narg)
C Local Variables
      integer a2vec,nsep,nxargs
      character strn*80
      logical nxtarg,cmdstr,lsequ,ltmp

      nsep = len(sep)
C --- Read argument, possible error if missing ---
      strn = ' '
      nread = -999
      if (iarg .le. 0) then
        if (iarg .lt. 0) ltmp = nxargs(nxargs(-1)-1) .eq. 0
        if (.not. nxtarg(strn)) goto 99
      else
        if (.not. cmdstr(iarg,strn)) goto 99
      endif

      ip = 0
      nread = 0
      ltmp = lsequ(strn,tok,len(tok),' ',ip) .or. tok .eq. ' '
      if (.not. ltmp) then
        goto 99
      endif
      ip = ip-1
      if (tok .eq. ' ') ip = 0

      if (narg .eq. 0) then
        nread = -1
        call chrps2(strn,sep,nsep,len(strn),ip,it)
        if (it(1) .ne. 0) nread = 1
        goto 99
      endif
      nread = a2vec(strn,len(strn),ip,cast,sep,nsep,itrm,narg,it,res)

   99 garg = nread
      end
      integer function parg(tok,cast,strn,ip,lstr,sep,itrm,narg,it,res)
C- Returns vector of binary values from a string
C ----------------------------------------------------------------------
Ci Inputs
Ci   tok:  token marking input
Ci  cast:  0=logical, 2=int, 3=real, 4=double
Ci  strn(ip:lstr):string to parse, from (ip=0 for first char)
Ci  lstr:  length of strn
Ci  sep:   string of characters, each of which separates arguments
Ci  itrm:  characters sep(itrm:*) signal the last argument
Ci  narg:  number of values to parse.
Cio Inputs/Outputs
Co   ip:   on input, position in strn where to start parsing
Co         on ouput, position in strn on exit.
Co Outputs
Co   res:  Vector of numbers that were converted
Co   it:   Vector of indices, one for each entry in res, labeling
Co         which char in 'sep' terminated the expr. for that entry.
Co parg:   0 if token is not matched in strn.
Co         n if token match and converted sans error narg numbers
Co           (for narg=0, returns 1 if token matched)
Co        -n if error on conversion of argument n
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer lstr,ip,cast,narg,itrm,it(1)
      character*(*) tok,sep,strn
      double precision res(narg)
C Local Variables
      logical ldum,parstr
      character term*1
      integer jp,np,nsep,lentok,a2vec


      nsep = len(sep)
      lentok = len(tok)
      term = tok(lentok:lentok)

C --- Find end of string ---
      jp = ip
      it(1) = 0
      if (itrm .le. nsep)
     .call chrps2(strn,sep(itrm:nsep),nsep-itrm+1,lstr,jp,it)
      if (it(1) .ne. 0) then
        np = jp
      else
        np = lstr
      endif

C     print *, 'np,lstr=',np,lstr,ip,jp

C --- Parse for tok within string strn, returning 0 if missing  ---
C      print *, tok
C      print *, strn

      if (tok .ne. ' ') then
        if (narg .eq. 0 .and. np .eq. lentok) then
          ip = np+1
          parg = 0
          if (strn(1:np) .eq. tok(1:np)) parg = 1
          return
        elseif (.not. parstr(strn,tok,np-lentok,lentok,term,ip,jp)) then
          parg = 0
          ip = np+1
          return
        endif
      else
        jp = ip
      endif

C --- Parse for vector of binary values to convert
      if (narg .eq. 0) then
        ip = jp
        parg = 1
        return
      endif

      ip = jp
      parg = a2vec(strn,np,ip,cast,sep,nsep,itrm,narg,it,res)

      end


