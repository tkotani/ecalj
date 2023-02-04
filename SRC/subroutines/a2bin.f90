! no hope to improve this routine ! Replace this by python or something else!
! a2bin, a2vec are used at rdfiln,symvar, m_gtv
integer function a2vec(str,lstr,ip,cast,sep,nsep,itrm,nvec,ix,res)
  !- Parses ascii string for a vector of binary values
  !example:  parts(1:lparts)='(1/2-.749777179)*.727140463'
  !          ixx = a2vec(parts(1:lparts)//' ',lparts+1,ip,4,' ',1,1,1,ix,res)
  !          returns res=-0.18162309358489384
  ! ----------------------------------------------------------------------
  ! Inputs:
  !   str(ip:lstr-1):string to parse, from (ip=0 for first char)
  !   cast:        0=logical, 2=int, 3=real, 4=double
  !   sep(1:nsep): class of chars that separate arguments
  !   itrm,nvec:   index to sep which indicates last arg: parse terminates
  !                when sep(i>=itrm) encountered or nvec arguments parsed.
  !                nvec<0: sign used to allow nam=expr as an expression
  !                which has the additional effect of loading the nam
  !                with expr into the variables table.
  !                itrm<0:  sign used to skip over white space before
  !                parsing for next token.
  ! Outputs:
  !   res:         binary result
  !   ip:          position in str on exit
  !   ix:          indices to which characters terminated expression
  !   a2vec:       number of values actually converted (- if error)
  ! Remarks
  !u Updates
  !u   26 May 07   When last expr ends at lstr w/out term, return w/out
  !u               error only if a2bin successfully parses expr
  !u   02 Feb 01   return w/out error when last expr ends at lstr w/out term
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed Parameters
  integer :: ip,lstr,cast,nsep,nvec,itrm,ix(1)
  character str*(*), sep(*)*1
  character c_sep
  double precision :: res(1)
  ! Local Variables
  integer :: iv,ip0,iitrm
  logical :: a2bin

  iitrm = iabs(itrm)

  ! --- For each element do ---
  a2vec = 0
  do  10  iv = 1, iabs(nvec)

     !   --- Determine appropriate terminator, possibly exit ---
     !   ... ip0 points to start of next token
12   if (itrm < 0) call skipbl(str,lstr,ip)
     if (ip >= lstr) return
     ip0 = ip
     call chrps2(str,sep,nsep,lstr,ip0,ix(iv))
     if (ip0 >= lstr) ix(iv) = 0
     !       No terminator was found: possibly exit with error
     if (ix(iv) == 0) then
        c_sep=' ' ! any value is OK.
     else
        c_sep=sep(ix(iv))
     endif
     if (ix(iv) == 0) then
        !         Still ok if expr runs to last char in string with no term
        ip0 = ip
        ! i          if (.not.a2bin(str,res,cast,iv-1,sep(ix(iv)),ip0,lstr-1)) then
        if ( .NOT. a2bin(str,res,cast,iv-1,c_sep,ip0,lstr-1)) then
           a2vec = -iv
           return
        endif
     endif
     ! i        if (.not. a2bin(str,res,cast,iv-1,sep(ix(iv)),ip,lstr-1)) then
     if ( .NOT. a2bin(str,res,cast,iv-1,c_sep,ip,lstr-1)) then
        ! ...     Try for nam=expr
        if (nvec < 0) then
           call skipbl(str,lstr,ip)
           ip0 = ip
           call parsyv(str,lstr,1,0,ip)
           if (ip == ip0) then
              a2vec = -iv
              return
           endif
           !           call shosyv(0,0,0,6)
           if (ix(iv) >= iitrm) return
           goto 12
        endif
        ! ...     Error exit
        a2vec = -iv
        return
     endif
     a2vec = iv
     if (ix(iv) >= iitrm) return
10 enddo
end function a2vec

logical function a2bin(instr,res,cast,count,term,j,jmaxi)
  !- Convert ASCII to logical, integer, real, double
  ! ----------------------------------------------------------------
  !i Inputs
  !i   cast:  0=logical, 1=char, 2=int, 3=real, 4=double
  !i   j:     offset to first character in string to read
  !i   term:  character that terminates the expression
  !i   jmaxi: j is not to exceed jmax
  !i          Also, jmaxi should not exceed (length-of-instr)-1
  !o Outputs
  !o   count'th element of array res is set to converted value
  !o   j is put past last character read
  !r Remarks
  !r   An ASCII string is converted to a double which is then cast to
  !r   the appropriate type and assigned to the count'th element of res.
  !r   A2BIN conforms to the 1977 ANSI FORTRAN standard.  Without using
  !r   recursion an ASCII string is interpreted with all normal operator
  !r   precidence and associativity through the use of a stack for
  !r   numbers and a stack for deferred operations.  The array OPRULE is
  !r   used to say when operations are deferred and when they are
  !r   executed.  The rules of precidence and associativity for all
  !r   operators are derived solely from the values in OPRULE.  Function
  !r   calls (sin, log, sqrt, ...) are parsed as unary operators.
  !r   OPRULE is initialized in a block data program and its contents
  !r   are this:
  !r        Current operator       Current operator is the most recently
  !r                             parsed operator from the ASCII string.
  !r T    ~ ^ * + < & | ? : ( )  Top operator is the operator on the top
  !r o  ~ F T T T T T T T T F T  of the operator stack.  If OPRULE is
  !r p  ^ F F T T T T T T T F T  .true. for a situation then the top
  !r    * F F T T T T T T T F T  operator is popped and executed.  This
  !r o  + F F F T T T T T T F T  continues until there are no more
  !r p  < F F F F T T T T T F T  operators or until OPRULE yields
  !r e  & F F F F F T T T T F T  .false., then the current operator is
  !r r  | F F F F F F T T T F T  pushed onto the operator stack.
  !r a  ? F F F F F F F F F F T  Parentheses are treated as operators.
  !r t  : F F F F F F F T T F T    A2BIN essentially does an on-the-fly
  !r o  ( F F F F F F F F F F F  conversion from infix to postfix and
  !r r  ) T T T T T T T T T T T  valuates the postfix expression as it is
  !r                             converted.
  !r
  !r   The operators are in order of precedence with associativity:
  !r 1> - (arithmetic negative), ~ (logical negative, .NOT.), and
  !r     functions abs(), exp(), log(), sin(), asin(), sinh(), cos(),
  !r      acos(), cosh(), tan(), atan(), tanh(), flor(), ceil(), erfc(),
  !r      and sqrt() (flor() rounds to the next lowest integer;
  !r      ceil() rounds up, erfc() is the the error function complement).
  !r 2> ^ (exponentiation)
  !r 3< * (times), / (divide), % (modulus)
  !r 4< + (add), - (subtract)
  !r 5< < (.lt.); > (.gt.); = (.eq.); <> (.ne.); <= (.le.);  >= (.ge.)
  !r 6< & (.and.)
  !r 7< | (.or.)
  !r 8&9  ?: conditional operators
  !r 10&11 parentheses
  !r
  !r   There is a granularity to the comparison operators (<, >, =, ...)
  !r   that is given by the variable MACHEP.  E.g., A=B if A and B are
  !r   within MACHEP of each other.  MACHEP is a regular variable; if
  !r   not set, it takes a default value of 1E-10.
  !r
  !r   The conditional operators work as follows, the expression
  !r   <logical>?<exp1>:<exp2> has the value <exp1> if <logical> is
  !r   .true.  or <exp2> if <logical> is .false.  If you put conditional
  !r   expressions inside any of the three expressions then A2BIN should
  !r   parse it in the common sense manner but parentheses should be
  !r   used just to be on the safe side.
  !r
  !r   The variables PTKTYP and PREVOP store the most recent token type
  !r   and operator. This is used by GETTOK to distinguish unary
  !r   arithmetic negation from binary subtraction, since they both use
  !r   the same operator.
  !r
  !r   PARENC is used to determine where the expression ends when the
  !r   terminator character is also an operator.
  !u Updates
  !u   19 Dec 02 Added macro handling
  ! --------------------------------------------------------------------
  !     implicit none
  character(1) :: instr(0:*),term
  integer :: cast,count,j,jmaxi
  double precision :: res(0:*)
  ! ... Local parameters
  integer :: namlen,parenc
  integer :: opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok,maxsiz
  double precision :: numtok,numstk(0:32),machep
  parameter (maxsiz=72,namlen=40)
  double precision :: dum
  character*(namlen) vartok,macstr*256,strn*256
  logical   oprule(0:10,0:10)
  common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp, &
       prevop,opstk,optok,parenc,oprule,vartok
  logical ::   gettok,fstpas,lerr
  integer ::   tj,ival,jmax,tm,k,t0,jmaxm
  save fstpas
  data fstpas /.true./

  ! --- Initialization: stacks empty, '-' must be unary negation ---
  a2bin= .false.
  tm = -1
  tj = j
  jmax = jmaxi
  if (jmax < 0) jmax = 999999
  numnum = 0
  opnum = 0
  parenc = 0
  machep = -1D0
  call numsyv(ival)
  if (ival == 0 .OR. fstpas) then
     call getsyv('t',numtok,ival)
     if (ival == 0) call addsyv('t',1D0,ival)
     call getsyv('f',numtok,ival)
     if (ival == 0) call addsyv('f',0D0,ival)
     call getsyv('pi',numtok,ival)
     if (ival == 0) call addsyv('pi',4*datan(1d0),ival)
     fstpas = .false.
  endif
  toktyp = 3
  optok = 9
  call skipbl(instr,jmax,tj)
  ! --- Get next token (variable, number, or operator) from ASCII string
10 continue
  if (tm < 0) then
     t0 = tj
     if ( .NOT. gettok(instr,tj,term,jmax)) goto 20
  else
     t0 = tm
     if ( .NOT. gettok(macstr,tm,' ',jmaxm)) goto 20
     !       end of macro expansion
     if (tm >= jmaxm) tm = -1
  endif
  goto (30,40,50,50,60),toktyp
  return
  ! ... handle variable token: get value and treat as number token
30 call getsyv(vartok,numtok,ival)
  write(6,*)'vvvtok  ',vartok,numtok,ival
  if (ival == 0) return
  ! ... handle number token: put it on the number stack
40 numnum = numnum+1
  if (numnum > 32) stop 'A2BIN: expression too complex.'
  numstk(numnum) = numtok
  goto 10
  ! ... handle operator token: pop a few ops and push current op
50 continue
  ! ... handle vector token by assigning ival to # stack
  if (toktyp == 4) then
     call getsvv(vartok,ival,0,1,1,dum)
     if (ival == 0) return
     numnum = numnum+1
     if (numnum > 32) stop 'A2BIN: expression too complex.'
     numstk(numnum) = ival
  endif
  if ((opnum == 0) .OR. ( .NOT. oprule(mod(opstk(opnum),16), &
       mod(optok,16)))) goto 56
  opnum = opnum-1
  ! ... doop() expects the old op to already be popped, hence opnum+1
  call doop(opstk(opnum+1),lerr)
  if (lerr) return
  goto 50
56 opnum = opnum+1
  if (opnum > 32) stop 'A2BIN: expression too complex.'
  if (optok == 9) parenc = parenc+1
  if (optok == 10) parenc = parenc-1
  opstk(opnum) = optok
  goto 10
  ! ... handle macro token
60 continue
  do  k = t0, jmax
     strn(k-t0+1:k-t0+1) = instr(k)
     !       print *, k, strn(1:40)
     if (instr(k) == ')') then
        tj = k+1
        strn(k-t0+2:k-t0+2) = ' '
        goto 61
     endif
  enddo
61 continue
  k = 1
  call macevl(strn,macstr,k)
  if (k < 0) call rx('a2bin: could not parse macro')
  call word(macstr,1,k,jmaxm)
  tm = 0
  goto 10

  ! --- end expression, pop remaining operators ---
20 if (opnum == 0) goto 25
  opnum = opnum-1
  call doop(opstk(opnum+1),lerr)
  if (lerr) return
  goto 20
  ! --- If expression passes the last few syntax checks, return it
25 continue
  !     if (numnum .eq. 0) return
  if (numnum /= 1) return
  ! ... j,res are not modified unless there is a valid expression
  a2bin= .true.
  j = tj
  call setans(res,res,res,res,cast,count,numstk(1))
end function a2bin

subroutine doop(op,lerr)
  !- Perform the passed operation on the top 1, 2, or 3 numbers of the
  !- number stack.
  ! -------------------------------------------------------------------
  !i Inputs
  !i   op:   operator code
  !r Remarks
  !r   lowest 4 bits of op are the operator class, next 28 bits are
  !r   the operator number.  The '?' and '(' operators do nothing
  !r   and should not be passed to DOOP in a syntactically correct
  !r   expression.  All the ')' operator does is pop the associated
  !r   '(' off the operator stack.
  ! -------------------------------------------------------------------
  !     implicit none
  integer ::   op
  integer ::   namlen,parenc
  integer ::   opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok
  double precision :: numtok,numstk(0:32),machep
  logical   oprule(0:10,0:10),lerr
  parameter (namlen=40)
  character*(namlen) vartok
  common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp, &
       prevop,opstk,optok,parenc,oprule,vartok
  integer ::   t1,t2,ival
  double precision ::  v0=1d99,v1=1d99,log2d,derfc
  real :: ran1
! #if HANSMR
!   integer :: ixx,l
!   double precision :: e,rsm,xi(0:10),phi(0:10)
! #endif

  ! --- Get operator class and number ---
  lerr = .false.
  t1 = op/16
  t2 = mod(op,16)
  ! ... Make sure that there enough numbers on the number stack to do
  ! ... the operation
  if ((t2 == 0) .AND. (numnum == 0)) goto  1090
  if ((t2 == 8) .AND. (numnum < 3)) goto  1090
  if ((t2 > 0) .AND. (t2 < 7) .AND. (numnum < 2)) goto  1090
  ! ... Load the top 2 numbers into v0 and v1 to speed things up
  v0 = numstk(numnum)
  if (numnum /= 1) v1 = numstk(numnum-1)
  ! ... Except for unary ops and ')', pop v0 off the number stack
  if ((t2 /= 0) .AND. (op /= 10)) numnum=numnum-1
  ! ... Array element is really a binop, so decrement numnum
  if ((t2 == 0) .AND. (t1 == 19)) numnum=numnum-1

  ! --- Go to appropriate operator class code
  goto ( 1000, 1010, 1020, 1030, 1040, 1050, 1060, 1070, &
       1080, 1090, 1100),t2+1

  ! --- Go to the appropriate unary op/function code
1000 goto (100,101,102,103,104,105,106,107,108,109,110,111,112,113, &
       114,115,116,117,118,119),t1+1
  ! ... Arithmetic negation
100 numstk(numnum) = -v0
  return
  ! ... Logical not
101 numstk(numnum) = log2d(v0.eq.0D0)
  return
  ! ... Functions
102 numstk(numnum) = dabs(v0)
  return
103 numstk(numnum) = dexp(v0)
  return
104 numstk(numnum) = dlog(v0)
  return
105 numstk(numnum) = dsin(v0)
  return
106 numstk(numnum) = dasin(v0)
  return
107 numstk(numnum) = dsinh(v0)
  return
108 numstk(numnum) = dcos(v0)
  return
109 numstk(numnum) = dacos(v0)
  return
110 numstk(numnum) = dcosh(v0)
  return
111 numstk(numnum) = dtan(v0)
  return
112 numstk(numnum) = datan(v0)
  return
113 numstk(numnum) = dtanh(v0)
  return
  ! ... flor()
114 if (v0 /= int(v0)) then
     if (v0 < 0D0) then
        numstk(numnum) = int(v0)-1D0
     else
        numstk(numnum) = int(v0)
     endif
  endif
  return
  ! ... ceil()
115 if (v0 /= int(v0)) then
     if (v0 > 0D0) then
        numstk(numnum) = int(v0)+1D0
     else
        numstk(numnum) = int(v0)
     endif
  endif
  return
  ! ... erfc()
116 numstk(numnum) = derfc(v0)
  return
117 continue
  if (v0 < 0D0) call rx('A2BIN: out of domain in SQRT')
  numstk(numnum) = sqrt(v0)
  return
  ! ... ran()  (to set seed, call ran(#), with # nonzero integer)
118 continue
  if (v0 /= 0d0) call ran1in(nint(v0))
  numstk(numnum) = ran1()
  return
  ! ... Return vec[v1], element v0
119 continue
  !      ii1 = v1+dble(1)/2
  !      ii2 = v0+dble(1)/2
  !      call getsvv(' ',ii1,1,ii2,ii2,numstk(numnum))
  if (numnum < 1) then
     lerr = .true.
     return
  endif
  call getsvv(' ',nint(v1),1,nint(v0),nint(v0),numstk(numnum))
  return
  ! --- Binary operators ---
  ! ... exponentiation
1010 numstk(numnum) = v1**v0
  return
  ! ... *, /, and %
1020 goto (200,201,202),t1+1
200 numstk(numnum) = v1*v0
  return
201 numstk(numnum) = v1/v0
  return
202 numstk(numnum) = mod(v1,v0)
  return
  ! ... + and - (binary subtraction)
1030 continue
  if (t1 == 0) then
     numstk(numnum) = v1+v0
  else
     numstk(numnum) = v1-v0
  endif
  return
  ! ... Conditional operators
1040 if (machep < 0D0) then
     call getsyv('machep',machep,ival)
     if (ival == 0) machep = 1d-10
  endif
  goto (400,401,402,403,404,405),t1+1
400 numstk(numnum) = log2d(v1.lt.v0)
  if (numstk(numnum) /= 0D0) goto 405
  return
401 numstk(numnum) = log2d(v1.gt.v0)
  if (numstk(numnum) /= 0D0) goto 405
  return
402 numstk(numnum) = log2d(abs(v1-v0).le.machep)
  return
403 numstk(numnum) = log2d(v1.lt.v0)
  if (numstk(numnum) == 0D0) goto 402
  return
404 numstk(numnum) = log2d(v1.gt.v0)
  if (numstk(numnum) == 0D0) goto 402
  return
405 numstk(numnum) = log2d(abs(v1-v0).gt.machep)
  return
  ! ... Logical and
1050 numstk(numnum) = v1*v0
  return
  ! ... Logical or
1060 numstk(numnum) = abs(v0)+abs(v1)
  return
  ! ... '?' or '(' mean syntax error in expression
1070 continue
1090 continue
  stop 'A2BIN: expression syntax.'
  ! ... ')', pop the associated '(' and do nothing else
1100 continue
  lerr = ((opnum .eq. 0) .or. (opstk(opnum) .ne. 9))
  opnum = opnum-1
  return
  ! ... Conditional expr: pop v1 and load either v0 or v1 onto the stack
1080 numnum = numnum-1
  if (numstk(numnum) == 0D0) then
     numstk(numnum) = v0
  else
     numstk(numnum) = v1
  endif
  ! ... Make sure that there is a '?' and pop it
  lerr = ((opnum.eq.0) .or. (opstk(opnum).ne.7))
  opnum = opnum-1

end subroutine doop

subroutine setans(resL,resI,resR,resD,cast,count,val)
  !- cast val appropriately and copy to count'th element of res*
  ! -------------------------------------------------------------------
  !i Inputs
  !i   cast:   0=logical, 1=char, 2=int, 3=real, 4=double
  !i   count:  element to copy to
  !i   val:    value to copy in double form
  !o Outputs
  !o   resL:   logical array, count'th element is copied over
  !o   resI:   integer array, count'th element is copied over
  !o   resR:   real array, count'th element is copied over
  !o   resD:   double array, count'th element is copied over
  !r Remarks
  !r   What's the purpose of casting to char? I don't know, so you
  !r   get an error.
  ! -------------------------------------------------------------------
  !     implicit none
  logical   resL(0:*)
  integer ::   resI(0:*),cast,count
  real ::      resR(0:*)
  double precision ::  resD(0:*),val

  goto (10,20,30,40,50),cast+1
20 stop 'A2BIN: can''t cast result.'
10 resL(count) = val .ne. 0D0
  return
30 resI(count) = int(val)
  return
40 resR(count) = val
  return
50 resD(count) = val
  return
end subroutine setans

real(8) function log2d(l)
  !- Convert logical to double
  ! -------------------------------------------------------------------
  !i Inputs
  !i   l:   logical to be converted
  !o Outputs
  !o   Returns 1D0 for .true. and 0D0 for .false.
  ! -------------------------------------------------------------------
  !     implicit none
  logical ::   l
  log2d = 0D0
  if (l) log2d = 1D0
  return
END function log2d
logical function a2d(strn,lskp,j,jmax,term,res)
  !- Convert ASCII string to number
  ! -------------------------------------------------------------------
  !i Inputs
  !i   strn:  string containing number
  !i   j:     location in string to start parsing number
  !i   jmax:  maximum value of j
  !i   lskp:  true, when a2d should skip blanks in number
  !i   term:  exit if a2d encounters character term
  !o Outputs
  !o   j:     location just after number in strn
  !o   res:   double precision number
  !o   a2d:   true if parse successful, false if not
  !r Remarks
  !r   a2d has a generic parser; however it is faster and more accurate
  !r   to use the Fortran internal read.  An unformatted read is
  !r   preferable, but does not conform to the ANSII 77 standard.  Do
  !r   Do ccomp -dUNFORMATTED_READ a2bin.f for this version
  !r   Some machines, such as the Cray, will correctly read an arbitrary
  !r   floating-point number using En.0 format; for this version do
  !r   comp -dFORMATTED_READ a2bin.f
  !r   Use the generic parser if your machine accepts neither of these.
  !r   ISFRAC and E1 are used to parse to the left of the decimal point.
  !r   OK is used to make sure that there is at least one decimal digit
  !r   in the mantissa and the exponent (if there is one).
  ! -------------------------------------------------------------------
  !     implicit none
  integer ::   i,j0,j,jmax,iprint
  logical ::   lskp
  character(1) ::       t,strn(0:*),term
  double precision :: res,sgn
  ! if FORMATTED_READ
  ! else
  !      logical isfrac,neg,ok
  !      integer e,e1
  ! endif

  ! if FORMATTED_READ | UNFORMATTED_READ
  integer :: maxs,strs
  parameter (maxs=72)
  character*(maxs) strn2

  ! --- Find end of string; exit to 20 w/ strn(j-1) = last char of num ---
  a2d = .false.
  j0 = j
  if (lskp) call skipbl(strn,1000,j0)
  j = j0-1
10 j = j+1
  if (jmax > 0 .AND. j > jmax) goto 20
  if (j > j0+maxs) goto 99
  t = strn(j)
  if (t == term) goto 20
  if (t >= '0' .AND. t <= '9') goto 10
  if (t == '.') goto 10
  if (t == 'd' .OR. t == 'D' .OR. t == 'e' .OR. t == 'E') then
     if (lskp .AND. strn(j) == ' ') call skipbl(strn,1000,j)
     if ((strn(j+1) == '+') .OR. (strn(j+1) == '-')) j = j+1
     goto 10
  endif
  if (j == j0 .AND. (t == '+' .OR. t == '-')) goto 10
  if (lskp .AND. (t == ' ' .OR. t == '        ')) goto 10
  ! --- Copy strn to strn2 for FORTRAN read ---
20 continue
  if (lskp .AND. j > j0 .AND. &
       (strn(j-1) == ' ' .OR. strn(j-1) == '        ')) then
     j = j-1
     goto 20
  endif
  strs = j-j0
  strn2 =  ' '
  call strcop(strn2(maxs+1-strs:maxs),strn(j0),strs,'z',i)
  ! if FORMATTED_READ
  !      read(strn2,'(E72.0)',err=99) res
  ! else
  read(strn2,*,err=99) res
  ! endif
  a2d = .true.
  if (iprint() >= 130) print 333, strn2(maxs+1-strs:maxs), res
333 format(' a2d: converted ', a,' to',g24.16)
  return
  ! --- Error handling ---
99 continue
  if (iprint() >= 60) print 332, (strn(i), i=j0, j-1)
332 format(' a2d: parse error for string ',72a1)
  !$$$#else
  !$$$c #else  (! FORMATTED_READ | UNFORMATTED_READ) <---takao think this is wrong.
  !$$$C --- Generic parser: initialization ---
  !$$$      sgn = 1
  !$$$      res = 0
  !$$$      e = 0
  !$$$      e1 = 0
  !$$$      isfrac = .false.
  !$$$      ok = .false.
  !$$$      j0 = j
  !$$$      if (lskp) call skipbl(strn,1000,j0)
  !$$$      j = j0
  !$$$      a2d = .false.
  !$$$C --- Parse mantissa ---
  !$$$   10 t = strn(j)
  !$$$      if (lskp .and. t .eq. ' ') goto 20
  !$$$      if (t .eq. '.') then
  !$$$        if (isfrac) return
  !$$$        isfrac = .true.
  !$$$        goto 20
  !$$$      elseif (j.eq.j0 .and. t.eq.'-') then
  !$$$        sgn = -1
  !$$$        goto 20
  !$$$      elseif (j.eq.j0 .and. t.eq.'+') then
  !$$$        goto 20
  !$$$      endif
  !$$$      if ((t.lt.'0') .or. (t.gt.'9')) goto 30
  !$$$      ok = .true.
  !$$$      res = res*10D0 + ichar(t)-ichar('0')
  !$$$      if (isfrac) e1 = e1-1
  !$$$   20 j = j+1
  !$$$      goto 10
  !$$$C --- Parse exponent ---
  !$$$   30 continue
  !$$$      if (.not. ok) return
  !$$$      ok = .false.
  !$$$      if ((t.eq.'d').or.(t.eq.'D').or.(t.eq.'e').or.(t.eq.'E')) then
  !$$$        j = j+1
  !$$$        if (lskp) call skipbl(strn,1000,j0)
  !$$$        if ((strn(j).eq.'+') .or. (strn(j).eq.'-')) then
  !$$$          neg = strn(j) .eq. '-'
  !$$$          j = j+1
  !$$$          if (lskp) call skipbl(strn,1000,j0)
  !$$$        else
  !$$$          neg = .false.
  !$$$        endif
  !$$$   35   t = strn(j)
  !$$$        if ((t.lt.'0') .or. (t.gt.'9')) goto 38
  !$$$        ok = .true.
  !$$$        e = e*10 + ichar(t)-ichar('0')
  !$$$        j = j+1
  !$$$        if (lskp) call skipbl(strn,1000,j0)
  !$$$        goto 35
  !$$$   38   continue
  !$$$        if (.not. ok) return
  !$$$        if (neg) e= -e
  !$$$      endif
  !$$$      res = sgn*res*(10D0**(e+e1))
  !$$$      a2d = .true.
  !$$$
  !$$$      print *, 'a2d: converted ',(strn(i),i=j0,j-1), ' to ',res
  !$$$#endif
end function a2d

logical function gettok(str,j,term,jmax)
  !- Get next variable, number, or operator token
  ! -------------------------------------------------------------------
  !i Inputs
  !i   str:   string containing token
  !i   j,jmax:location in str to start parsing token, not to exceed jmax
  !i   term:  expression terminator character
  !o Outputs
  !o   j:     location just after token
  !o   Returns whether a token was found
  !o   toktyp: 1, a variable
  !o           2, a number
  !o           3, an operator or a function
  !o           4, a vector
  !o           5, a macro
  !o           0, the 'assignment' operator (causes a2bin to return F)
  !r Remarks
  !r   If something that isn't some kind of token is encountered then
  !r   GETTOK returns .false.
  ! -------------------------------------------------------------------
  !     implicit none
  character(1) :: term,str(0:*)
  integer ::   j,jmax,namlen,maxsiz,parenc,nunop,k
  integer ::   opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok
  double precision :: numtok,numstk(0:32),machep,res(0:1)
  parameter (namlen=40,maxsiz=72,nunop=17)
  character*(namlen) vartok
  character*(maxsiz) strn
  logical   oprule(0:10,0:10)
  common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp, &
       prevop,opstk,optok,parenc,oprule,vartok
  character(5) :: unop(0:nunop-1), unops(0:10),ctmp
  integer ::   i,s,t
  logical :: a2d,namchr
  ! ... Allowed characters in a name
  namchr(ctmp) =ctmp .ge. 'a' .and. ctmp .le. 'z' .or. ctmp .eq. '_' &
       .or. ctmp .ge. 'A' .and. ctmp .le. 'Z' &
       .or. ctmp .ge. '0' .and. ctmp .le. '9'
  data unops/'A~','^','*/%','+-','<>=','&','|','?',':','(',')'/
  data unop/'abs','exp','log','sin','asin','sinh','cos','acos', &
       'cosh','tan','atan','tanh','flor','ceil','erfc','sqrt', &
       'ran'/

  ! --- Check for terminator or end of expression
  gettok = .false.
  if (j > jmax) return
  if ((str(j) == term) .AND. (parenc == 0)) then
     j = j+1
     return
  endif
  ! --- Save previous token type and operator code so that GETTOK can
  !     distinguish unary negation from binary subtraction
  ptktyp = toktyp
  prevop = optok
  gettok = .true.
  ! --- Check for variable token first, variables are a letter followed
  !     by a string of letters and digits
  if (((str(j) >= 'A') .AND. (str(j) <= 'Z')) .OR. &
       ((str(j) >= 'a') .AND. (str(j) <= 'z'))) then
     s = j
10   j = j+1
     if (j > jmax) then
        if (s > jmax) j = j-1
        goto 11
     endif
     if (namchr(str(j))) goto 10
     !        if (((str(j).ge.'A') .and. (str(j).le.'Z')) .or.
     !     .      ((str(j).ge.'a') .and. (str(j).le.'z')) .or.
     !     .      ((str(j).ge.'0') .and. (str(j).le.'9'))) goto 10
11   strn = ' '
     if (j-s > namlen) call rx('A2BIN: variable name too long.')
     ! --- Get variable (isn't FORTRAN string handling great?)
     call strcop(strn(1:j-s),str(s),j-s,' ',i)
     vartok = strn(1:j-s)
     !   --- Check if variable is actually a macro
     if (str(j) == '(') then
        k = 0
        call macevl(strn,' ',k)
        if (k > 0) then
           toktyp = 5
           return
        endif
     endif
     ! --- Check if variable is actually a function or a vector
     call tokmat(vartok,unop,nunop,5,' ',i,s,.false.)
     if (i >= 0) then
        toktyp = 3
        optok = 16*(i+2)
     else
        call getsvv(vartok,i,0,1,1,res)
        if (i > 0) then
           toktyp = 4
           optok = 16*(17+2)
        else
           toktyp = 1
        endif
     endif
     return
     ! --- If it isn't a variable then check if it's a number
  else if ((str(j) == '.') .OR. &
       ((str(j) >= '0')  .AND. (str(j) <= '9'))) then
     if ( .NOT. a2d(str, .FALSE. ,j,jmax,term,numtok)) goto 99
     toktyp = 2
     return
     ! --- Otherwise it's an operator or a mistake
  else
     toktyp = 3
     do  30  i = 0, 10
        t = 0
        ! --- Check each operator class i for the operator
        call chrpos(unops(i),str(j),3,t)
        if (t /= 3) then
           ! --- Operator found, calculate the operator code
           optok = t*16 + i
           j = j+1
           !           Valid special case unop is ')' and j=jmax+1
           if (j <= jmax+1 .AND. i == 10) return
           if (j > jmax) goto 99
           ! --- See if it's unary negation
           if ((optok == 19) .AND. (ptktyp == 3) .AND. (prevop /= 10)) &
                optok = 0
           ! --- Distinguish '=' from '==': '=' only causes a2bin=F
           if (optok == 36) then
              if (str(j) /= '=') then
                 toktyp = 0
                 gettok = .true.
                 return
              endif
              j = j+1
              ! --- Distinguish '<=' from '<'
           else if ((optok == 4) .AND. (str(j) == '=')) then
              optok = 52
              j = j+1
              if (j > jmax) goto 99
              ! --- Distinguish '>=' from '>'
           else if ((optok == 20) .AND. (str(j) == '=')) then
              optok = 68
              j = j+1
              ! --- Distinguish '<>' from '<'
           else if ((optok == 4) .AND. (str(j) == '>')) then
              optok = 84
              j = j+1
           endif
           if (j > jmax) goto 99
           return
        endif
30   enddo
  endif
  ! --- If we fell through the op, check loop then an illegal character
  !     is in the ASCII string.  Cause gettok to return .false.
99 numnum = 0
  opnum = 0
  gettok = .false.
end function gettok

subroutine sstyle(style)
  !- Sets style of expression parsing
  ! -------------------------------------------------------------------
  !i Inputs
  !i   style:   0=old style, 1=new style
  !r Remarks
  !r   Old style parsing means expressions are evaluated from left to
  !r   right with no operator precidence, but parenthesis work as
  !r   expected.  New style parsing means expressions are parsed with
  !r   algebraic precidence and associativity.  For old style parsing
  !r   the rule array OPRULE is changed to:
  !r
  !r        Current operator
  !r
  !r T    ~ ^ * + < & | ? : ( )
  !r o  ~ F T T T T T T T T F T
  !r p  ^ F T T T T T T T T F T
  !r    * F T T T T T T T T F T
  !r o  + F T T T T T T T T F T
  !r p  < F T T T T T T T T F T
  !r e  & F T T T T T T T T F T
  !r r  | F T T T T T T T T F T
  !r a  ? F F F F F F F F F F T
  !r t  : F T T T T T T T T F T
  !r o  ( F F F F F F F F F F F
  !r r  ) T T T T T T T T T T T
  ! -------------------------------------------------------------------
  !     implicit none
  integer :: style,i,j
  logical oprule(0:10,0:10),orule(0:10,0:10),nrule(0:10,0:10)
  integer :: opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok,namlen
  integer :: parenc
  double precision :: numtok,numstk(0:32),machep
  parameter (namlen=40)
  character*(namlen) vartok
  common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp, &
       prevop,opstk,optok,parenc,oprule,vartok
  data nrule/10*.false.,2*.true.,9*.false.,4*.true.,7*.false., &
       5*.true.,6*.false.,6*.true.,5*.false.,7*.true., &
       4*.false.,8*.true.,3*.false.,8*.true.,.false.,.true., &
       .false.,8*.true.,.false.,.true.,.false.,.true., &
       10*.false.,10*.true.,.false.,.true./
  data orule/10*.false.,8*.true.,.false.,.true.,.false., &
       8*.true.,.false.,.true.,.false.,8*.true.,.false.,.true., &
       .false.,8*.true.,.false.,.true.,.false.,8*.true.,.false., &
       .true.,.false.,8*.true.,.false.,.true.,.false.,8*.true., &
       .false.,.true.,.false.,8*.true.,.false.,.true.,.false., &
       .true.,10*.false.,10*.true.,.false.,.true./

  if (style == 0) then
     do    i = 0, 10
        do    j = 0, 10
           oprule(i,j) = orule(i,j)
        enddo
     enddo
  else
     do    i = 0, 10
        do    j = 0, 10
           oprule(i,j) = nrule(i,j)
        enddo
     enddo
  endif
  return
end subroutine sstyle

block data da2bin
   !- Block data for A2BIN
   ! -------------------------------------------------------------------
   !p Purpose
   !p   Initialize the array OPRULE once.
   ! -------------------------------------------------------------------
   integer :: opnum,numnum,toktyp,ptktyp,prevop,opstk(0:32),optok,namlen
   integer :: parenc
   double precision :: numtok,numstk(0:32),machep
   parameter (namlen=40)
   character*(namlen) vartok
   logical   oprule(0:10,0:10)
   common /a2binc/ machep,numtok,numstk,opnum,numnum,toktyp,ptktyp, &
        prevop,opstk,optok,parenc,oprule,vartok
   data oprule/10*.false.,2*.true.,9*.false.,4*.true.,7*.false., &
        5*.true.,6*.false.,6*.true.,5*.false.,7*.true., &
        4*.false.,8*.true.,3*.false.,8*.true.,.false.,.true., &
        .false.,8*.true.,.false.,.true.,.false.,.true., &
        10*.false.,10*.true.,.false.,.true./
END block data
function ran1()
  !- Return a random deviate between 0.0 and 1.0.
  ! ----------------------------------------------------------------
  !i Inputs
  !o Outputs
  !o   ran1
  !r Remarks
  !r   Algorithm from Knuth; adapted here from Numerical Recipes, chapter 7.
  !r   Uses three linear
  !r   congruential generators (two for high and low order, a third
  !r   to shuffle).  Use ran1in to initialize.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Local
  real :: ran1
  integer :: j
  integer :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
  parameter (m1=259200,ia1=7141,ic1=54773)
  parameter (m2=134456,ia2=8121,ic2=28411)
  parameter (m3=243000,ia3=4561,ic3=51349)
  ! Static
  integer :: ix1,ix2,ix3
  real :: r(97),rm1,rm2
  common /dran1/ ix1,ix2,ix3,rm1,rm2,r

  ! Generate next number for each sequence;
  ! use third to generate random integer between 1 and 97
  ix1 = mod(ia1*ix1+ic1,m1)
  ix2 = mod(ia2*ix2+ic2,m2)
  ix3 = mod(ia3*ix3+ic3,m3)
  j = 1 + (97*ix3)/m3
! #if TEST
!   if (j > 97 .OR. j < 1) pause
! #endif
  ! Return the table entry ...
  ran1 = r(j)
  ! And refill it.
  r(j) = (float(ix1) + float(ix2)*rm2)*rm1
  return
end function ran1
subroutine ran1in(iseed)
  !- A simple one-parameter initializer for ran1
  ! ----------------------------------------------------------------
  !i Inputs
  !i   iseed
  !o Outputs
  !o   ran1 is set up
  !r Remarks
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed
  integer :: iseed
  ! Local
  integer :: j
  integer :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
  parameter (m1=259200,ia1=7141,ic1=54773)
  parameter (m2=134456,ia2=8121,ic2=28411)
  parameter (m3=243000,ia3=4561,ic3=51349)
  ! To preserve
  integer :: ix1,ix2,ix3
  real :: r(97),rm1,rm2
  common /dran1/ ix1,ix2,ix3,rm1,rm2,r
  rm1 = 1./m1
  rm2 = 1./m2
  ! Seed the first, second and third sequences
  ix1 = mod(ic1-iseed,m1)
  ix1 = mod(ia1*ix1+ic1,m1)
  ix2 = mod(ix1,m2)
  ix1 = mod(ia1*ix1+ic1,m1)
  ix3 = mod(ix1,m3)
  ! Fill table with sequential uniform deviates generated by first two
  do  11  j = 1, 97
     ix1 = mod(ia1*ix1+ic1,m1)
     ix2 = mod(ia2*ix2+ic2,m2)
     r(j) = (float(ix1) + float(ix2)*rm2)*rm1
11 enddo
  return
end subroutine ran1in
