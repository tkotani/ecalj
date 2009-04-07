      logical function a2bina(instr,res,cast,count,term,j,jmaxi)
C- Convert ASCII to logical, integer, real, double, with possible assignment
C ----------------------------------------------------------------
Ci Inputs
Ci   instr: contains string to be converted
Ci   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
Ci   j:     offset to first character in string to read
Ci   term:  character that terminates the expression
Ci   jmaxi: j is not to exceed jmax
Ci          Also, jmaxi should not exceed (length-of-instr)-1
Co Outputs
Co   count'th element of array res is set to converted value
Co   j is put past last character read
Cr Remarks
Cr   a2bina is identical to the standard expression parser, a2bin, with
Cr   the addition that a valid expression may include assignment
Cr   operators.  The following adds 2 to variable x, and returns value 2
Cr       x+=2
Cr  Sequences of expressions separated by commas, are permitted, viz:
Cr       y=6,x=2,y-=3,y/x
Cr  assigns y to 6, x to 2, subtracts 3 from y, and returns 1.5
Cr
Cr  Allowed assignment operators  'op'  for syntax "var op expr"  are
Cr   var op expr  variable assignment
Cr     =  assign expr to variable
Cr    *=  multiply expr into variable
Cr    /=  divide expr into variable
Cr    +=  add expr into variable
Cr    -=  subtract expr from variable
Cu Updates
Cu   27 Sep 04 Accomodate macro calls
C --------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*1 instr(0:*),term
      integer cast,count,j,jmaxi
      double precision res(0:1)
C ... Local parameters
      integer tj,jmax,nextop,i0,i,k
      character*1 aops(9), aterm(2), ctmp*1, a*80, lterm
      logical lassgn,a2bin,nxtexp,namchr,lmac
      double precision dum,dum2,v2dbl
C ... Allowed characters in a name
      namchr(ctmp) =ctmp .ge. 'a' .and. ctmp .le. 'z' .or. ctmp .eq. '_'
     .         .or. ctmp .ge. 'A' .and. ctmp .le. 'Z'
     .         .or. ctmp .ge. '0' .and. ctmp .le. '9'
      data aops/' ','=','*','/','+','-','^',' ',','/
      data aterm/' ',','/

      aops(1) = term
      aterm(1) = term
      a2bina = .false.
      tj = j
      jmax = jmaxi
      if (jmax .lt. 0) jmax = 999999

C ... re-entry point for multiple expressions
   10 continue
      lmac = .false.
      call skipbl(instr,jmax,tj)

C ... Check for assignment statement.  1st char must be allowed name char
      if (tj .gt. jmax) return
      j = tj
      lassgn =
     . namchr(instr(j)) .and. (instr(j) .lt. '0' .or. instr(j) .gt. '9')
C     Look for first non-word character
      if (lassgn) then
      do  12  tj = j, jmax
      a(tj-j+1:tj-j+2) = instr(tj) // ' '
   12 if (.not. namchr(instr(tj))) goto 14
C     Reached last character; no non-word character found
      lassgn = .false.
   14 continue

C ... If 1st word is a macro, statement is not an assignment
      if (lassgn) then
        k = 0
        call macevl(a,' ',k)
        if (k .gt. 0) lassgn = .false.
        if (k .gt. 0) lmac = .true.
      endif

C ... Next check for assignment statement. Inspect first non-word char
      if (lassgn) then
      call skipbl(instr,jmax,tj)
      call chrps2(instr,aops,9,tj,tj,nextop)
C     Set lassgn=.true. if some assignment operator was found
      lassgn = .false.
      if (tj .lt. jmax .and. nextop .gt. 1 .and. nextop .lt. 8) then
        i0 = tj
        lassgn = nextop .gt. 2 .and. instr(tj+1) .eq. '='
        if (lassgn) i0 = tj+1
        lassgn = lassgn .or. nextop .eq. 2 .and. instr(tj+1) .ne. '='
        if (lassgn) then
          a = ' '
          call strncp(a,instr(j),1,1,tj-j)
          j = i0+1
        endif
      endif
      endif
      endif
C     End of checking for assignment statement. If true, lassgn=.true.
C     j is updated to first character in expression

C ... See if there is an expression following this one
      nxtexp = .false.
      lterm = term
      tj = j
      if (.not. lmac) then
        call chrps2(instr,aterm,2,jmax,tj,i)
        if (i .eq. 2)  then
          lterm = aterm(2)
          nxtexp = .true.
        endif
      endif
      a2bina = a2bin(instr,res,cast,count,lterm,j,jmaxi)
      if (.not. a2bina) return

      if (lassgn) then
        dum = v2dbl(res,res,res,res,cast,count+1)
        dum2 = 0
        call getsyv(a,dum2,i)
        if (i .eq. 0) call addsyv(a,0d0,i)
        if (nextop-1 .eq. 2) dum = dum2*dum
        if (nextop-1 .eq. 3) dum = dum2/dum
        if (nextop-1 .eq. 4) dum = dum2+dum
        if (nextop-1 .eq. 5) dum = dum2-dum
        if (nextop-1 .eq. 6) dum = dum2**dum
        call chsyv(a,dum,i)
C       call shosyv(0,0,0,6)
      endif

      tj = j
      if (nxtexp) goto 10
      end
      double precision function v2dbl(resL,resI,resR,resd,cast,n)
C- Converts a number of different casts into a double
Ci   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
C     implicit none
      integer cast,n
      logical resL(n)
      integer resI(n)
      real resR(n)
      double precision resD(n)

      if (cast .eq. 0) then
        v2dbl = 0
        if (resL(n)) v2dbl = 1
      elseif (cast .eq. 2) then
        v2dbl = resI(n)
      elseif (cast .eq. 3) then
        v2dbl = dble(resR(n))
      elseif (cast .eq. 4) then
        v2dbl = resD(n)
      else
        call rx('v2dbl: bad cast')
      endif
      end

