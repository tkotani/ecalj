subroutine cpstr(strn,lstr,opt,delim,is,io,sout)
  !- Copy a string, excluding delimiters
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   strn    string
  !i   lstr    string length
  !i   opt     1s digit
  !i           1: skip initial blanks
  !i         100s digit
  !i           1: compress repeated blanks into a single blank
  !i        1000s digit
  !i           1: purge delimiter from sout (last char copied unles eos)
  !i   delim   string delimiters marking end of string (any character
  !i           in delim counts as a delimiter
  !i   is      first character to copy (string starts at character 1)
  !o Outputs
  !o   sout    strn copied to sout; see Remarks
  !o   is      smaller of index in strn of last char copied and lstr+1
  !o   io      position of terminator, or 1+ last character copied if
  !o           limit of string reached beforehand
  !r Remarks
  !r   delimiters between " " or ' ' are excluded.  pairs " " and ' '
  !r   are excised in the output string sout
  !u Updates
  !u   01 Nov 01 delimiter may contain more than one character
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: lstr,is,io,opt
  character *(*) strn, delim*(*), sout*(*)
  character ch*60
  integer :: i2,opt0,opt2,opt3,it,ia,nch,ldelim

  sout = ' '
  opt0 = mod(opt,10)
  opt2 = mod(opt/100,10)
  opt3 = mod(opt/1000,10)
  ldelim = len(delim)
  ch = delim // '"'' '
  nch = ldelim + 2
  if (opt2 /= 0) nch = nch+1
  is = is-1
  ia = 1
  ! ... Skip past leading blanks
  if (opt0 /= 0) then
     call skipbl(strn,lstr,is)
  endif

  ! ... Find i2 : points to past last char of the argument
  i2 = is
12 is = i2
  call chrps2(strn,ch,nch,lstr,i2,it)
  if (it <= ldelim .AND. i2 < lstr) i2 = i2+1
  if (i2 > is) sout(ia:ia+i2-is-1) = strn(is+1:i2)
  if (it <= 0 .AND. i2 == lstr) i2 = i2+1
  ia = ia+i2-is
  ! ... A blank encountered ... skip blanks and continue
  if (it == ldelim+3 .AND. i2 < lstr) then
     call skipbl(strn,lstr,i2)
     sout(ia:ia) = ' '
     ia = ia+1
     goto 12
     ! ... A quote encountered ... find match and continue
  elseif (it > ldelim .AND. i2 < lstr) then
     i2 = i2+1
     is = i2
     call chrpos(strn,ch(it:it),lstr,i2)
     call strncp(sout,strn,ia,is+1,i2-is)
     ia = ia+i2-is
     i2 = i2+1
     goto 12
  endif

  io = ia-1
  if (is >= lstr) io = io+1
  !      if (is .ge. lstr) print *, 'hi',io
  if (opt3 == 1) then
     if (i2 <= lstr) sout(io:io) = ' '
  endif
  is = min(i2,lstr+1)

end subroutine cpstr

subroutine eostr(strn,lstr,opt,delim,is)
  !- Mark end of a string
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   strn,is string, and first character (1st character starts at 1)
  !i   delim   string delimiter marking end of string
  !i   opt     1s  digit 1: exclude delimiters between " " or ' '
  !i           10s digit 1: skip initial blanks
  !i                     2: use strn(is) as delimiter;  delim is not used
  !o Outputs
  !i   is      smaller of position in strn of delimiter and lstr+1
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: lstr,is,opt
  character *(*) strn, delim*1
  character ch*3
  integer :: i2,opt0,opt1,it,k
  data ch /' "'''/


  opt0 = mod(opt,10)
  opt1 = mod(opt/10,10)
  ch(1:1) = delim
  if (opt1 == 2) then
     ch(1:1) = strn(is:is)
  else
     is = is-1
  endif

  ! ... Skip past leading blanks
  if (opt1 == 1) then
     call skipbl(strn,lstr,is)
  endif

  ! ... Find i2 : points to past last char of the argument
  k = 1
  if (opt0 /= 0) k = 3
  i2 = is
12 is = i2

  call chrps2(strn,ch,k,lstr,i2,it)
  if (it == 1 .AND. i2 < lstr) i2 = i2+1
  ! ... A quote encountered ... find match and continue
  if (it > 1 .AND. i2 < lstr) then
     i2 = i2+1
     is = i2
     call chrpos(strn,ch(it:it),lstr,i2)
     i2 = i2+1
     goto 12
  endif

  is = min(i2,lstr+1)

end subroutine eostr
!      subroutine fmain
!      implicit none
!      integer is,io
!      character*80 strn,sout
!      strn = '  test" this string "now" and again"'
!      strn = '  test" this string "now" and again ",   and     '//
!     .  '" -until the very- "   very end'

!      print *, ' This is the entire string which will be copied:'
!      print '(a,a)', '1234567890123456789012345678901234567890'//
!     .  '1234567890123456789012345678901234567890',
!     .  ' (a ruler)'
!      print '(a)', strn

!C      is = 1
!C      call cpstr(strn,len(strn),1,' ',is,io,sout)
! c     call eostr(strn,len(strn),10,' ',is)
!C      print 333, '80 1',is,io,sout(1:io)

!      print *, ' '
!      print '(''lstr opt del      is  io sout...'')'
!      print '(''--------------------------------'')'

!      is = 1
!      call cpstr(strn,11,1,'s',is,io,sout)
!      print 333, '11    1  s    ',is,io,sout(1:io)
!      print 334, strn(1:min(is,11))

!      is = 1
!      call cpstr(strn,11,1001,'s',is,io,sout)
!      print 333, '11 1001  s    ',is,io,sout(1:io)
!      print 334, strn(1:min(is,11))

!      is = 1
!      call cpstr(strn,5,1,'z',is,io,sout)
!      print 333, '5     1  z    ',is,io,sout(1:io)
!      print 334, strn(1:min(is,5))


!      is = 1
!      call cpstr(strn,11,1,'z',is,io,sout)
!      print 333, '11    1  z    ',is,io,sout(1:io)
!      print 334, strn(1:min(is,11))

!      is = 1
!      call cpstr(strn,11,1,' ',is,io,sout)
!      print 333, '11    1  blank',is,io,sout(1:io)
!      print 334, strn(1:min(is,11))

!      is = 1
!      call cpstr(strn,80,1,' ',is,io,sout)
!      print 333, '80    1  blank',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))

!      is = 1
!      call cpstr(strn,80,11,'wz',is,io,sout)
!      print 333, '80    1  zw   ',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))

!      is = 1
!      call cpstr(strn,80,11,'zv',is,io,sout)
!      print 333, '80   11  zv   ',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))

!      is = 1
!      call cpstr(strn,80,1011,'zv',is,io,sout)
!      print 333, '80 1011  zv   ',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))

!      is = 1
!      call cpstr(strn,80,11,'v,',is,io,sout)
!      print 333, '80   11  v,   ',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))

!      is = 1
!      call cpstr(strn,80,11,'v,o',is,io,sout)
!      print 333, '80   11  v,o  ',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))

!      is = 1
!      call cpstr(strn,80,1011,'v,o',is,io,sout)
!      print 333, '80 1011  v,o  ',is,io,sout(1:io)
!      print 334, strn(1:min(is,80))


!      print *, '-------- end of tests of cpstr -------------'
!      print *, ' '
!      print *, ' '
!      print *, '-------- test eostr -------------'

!      strn = '  test" eostr for delim / this string "now" and again"'//
!     .  '  ? '
!      print *, ' This is the string we will test'
!      print '(a,a)', '1234567890123456789012345678901234567890'//
!     .  '123456789012345678901234567890',
!     .  ' (a ruler)'
!      print '(a)', strn
!      print *, ' '
!      is = 1
!      call eostr(strn,80,11,'/',is)
!      print 335, is, ' ''', strn(1:min(is,80)), ''''
!      is = 1
!      call eostr(strn,80,11,'?',is)
!      print 335, is, ' ''', strn(1:min(is,80)), ''''
!      is = 1
!      call eostr(strn,80,11,'"',is)
!      print 335, is, ' ''', strn(1:min(is,80)), ''''

!      strn = '"  test eostr for delim / this string "now" and again"'//
!     .  '  ? '
!      is = 1
!      call eostr(strn,80,21,'"',is)
!      print 335, is, ' ''', strn(1:min(is,80)), ''''

!  333 format(a,2x,2i4,1x,'|',a)
!  334 format(25x,'|',a)
!  335 format(i5,a,a,a)

!      end

