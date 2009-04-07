      subroutine macset(args,err)
C- Routines for macro definitions
C ----------------------------------------------------------------------
Ci Inputs
Ci   args   :ascii string defining macro, of form:
Ci          :macro_name(arg1,arg2,..) expand_string
Ci   err
Co Outputs
Cl Local variables
Cl   nmaca  : nmaca(i) = number of argments for macro i; limited to nmxarg
Cl   macarag:macarag(j,i) is the variable names for argument j to macro i
Cl          :First argument is macro name
Cr Remarks
Cr
Cu Updates
Cu   08 Jan 07  Increase the max number of macros
Cu   19 Dec 02  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nmac,nmacx,nmxarg,err
      character*(*) args,strout
      
C ... Local parameters
      parameter (nmacx=10,nmxarg=5)
C     logical savvar(nmacx)
C     double precision locval(nmacx)
      integer nmaca(nmacx),nchar(nmxarg,2)
      character*16 macarg(nmxarg,nmacx),argsub(nmacx)
      character*256 macexp(nmacx),strn,strn2
      integer j,j1,j2,i1,i2,iarg,imac,i,ls,errin
      save nmac,nmaca,macarg,macexp
      data nmac /0/

C ... too many macros ... exit err=-1
      err = -1
      if (nmac .ge. nmacx) goto 99

C ... No macro defined ... exit err=-2
      call word(args,1,j1,j2)
      err = -2
      if (j1 .gt. j2) goto 99
C ... macro name doesn't start with a letter ... exit err=-3
      err = -3
      call wordg(args(j1:),10,'A-Za-z_0-9',1,i1,i2)
      if (i2 .lt. i1) goto 99
C ... Assign macro name
      macarg(1,nmac+1) = args(j1+i1-1:j1+i2-1)
C ... First letter after macro name isn't '(' ... exit err=-4
      err = -4
      i1 = i2+j1
      call nword(args,1,i1,i2)
      if (args(i1:i1) .ne. '(') goto 99

C --- Assign macro arguments ---
      i1 = i1+1
      iarg = 1
   10 continue

C ... Error if macro has too many arguments
      iarg = iarg+1
      err = -5
      if (iarg .gt. nmxarg) goto 99

C ... Error if no terminator to macro
      call nwordg(args,1,'),',1,i1,i2)
      err = -6
      if (i2 .ge. len(args)) goto 99

      macarg(iarg,nmac+1) = args(i1:i2)
      i1 = i2+2
      if (args(i2+1:i2+1) .ne. ')') goto 10

C ... get macro string ... error if none
      err = -7
      call nword(args,1,i1,i2)
      if (i2 .lt. i1) goto 99
      macexp(nmac+1) = args(i1:i2)

C ... Macro arguments complete ... increment number of macros by 1
      nmac = nmac+1
      nmaca(nmac) = iarg
      err = 0
      goto 99

C --- macro evaluation ---
      entry macevl(args,strout,err)
C Input err:   0 -> just return index to macro in err if it exists
C Input err: <>0 -> return index to macro, and strout = expansion

      errin = err

C ... No macro defined ... exit err=-2
      call word(args,1,j1,j2)
      err = -2
      if (j1 .gt. j2) goto 99

C ... No match to existing macro names ... exit err=-1
      err = -1
      call wordg(args(j1:),10,'A-Za-z_0-9',1,i1,i2)
      if (i2 .lt. i1) goto 99
      do  imac = 1, nmac
        if (macarg(1,imac) .eq. args(j1+i1-1:j1+i2-1)) goto 15
      enddo
      goto 99

C ... macro found
   15 continue
      if (errin .eq. 0) then
        err = imac
        return
      endif

C ... First letter after macro name isn't '(' ... exit err=-4
      err = -4
      i1 = i2+j1
      call nword(args,1,i1,i2)
      if (args(i1:i1) .ne. '(') goto 99

C --- Assign macro arguments ---
      i1 = i1+1
      iarg = 1
   20 continue

C ... Error if macro has too many arguments
      iarg = iarg+1
      err = -5
      if (iarg .gt. nmaca(imac)) goto 99

C ... Error if no terminator to macro
      call nwordg(args,1,'),',1,i1,i2)
      err = -6
      if (i2 .ge. len(args)) goto 99

C ... For this variable, evaluate argument and load var table with res
      argsub(iarg) = args(i1:i2)
      i1 = i2+2
      if (args(i2+1:i2+1) .ne. ')') goto 20

C ... count number of characters in each macro string
      do  i = 2, nmaca(imac)
        call word(macarg(i,imac),1,j1,j2)
        nchar(i,1) = j2
        call word(argsub(i),1,j1,j2)
        nchar(i,2) = j2
      enddo

C     i = which string holds current value of macro
C     j = pointer to current index in string
      i = 1
      j = 1
      strn = macexp(imac)
      ls = len(strn)
      call skpblb(strn,len(strn),i)
      ls = i+1

C --- macro substitution ---
   30 continue
      do  i = 2, nmaca(imac)
        if (strn(j:j+nchar(i,1)-1) .eq. macarg(i,imac)) then
C         Also next character cann be part of a word
          call wordg(strn(j+nchar(i,1)-1:),10,'A-Za-z_0-9',1,i1,i2)
          if (i2 .le. i1) then
          strn2(1:ls) = strn
          strn(j:j+nchar(i,2)-1) = argsub(i)
          ls = ls + nchar(i,2)-nchar(i,1)
          strn(j+nchar(i,2):ls) = strn2(j+nchar(i,1):)
          if (nchar(i,2)-nchar(i,1) .lt. 0) then
            strn(ls+1:ls+1) = ' '
          endif
          j = j + nchar(i,2)-1
          goto 35
          endif
        endif
      enddo
   35 continue
C     Move past current word 
      call nwordg(strn,10,'A-Za-z_0-9',1,j,j2)
      if (j2 .gt. j) then
        j = j2+1
      else
        j = j+1
      endif
      if (j .le. ls) goto 30
      err = imac
      strout = strn(1:ls)

C --- Error exit ---
   99 continue

      end

