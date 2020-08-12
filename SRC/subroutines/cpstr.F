      subroutine cpstr(strn,lstr,opt,delim,is,io,sout)
C- Copy a string, excluding delimiters
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn    string
Ci   lstr    string length
Ci   opt     1s digit
Ci           1: skip initial blanks
Ci         100s digit
Ci           1: compress repeated blanks into a single blank
Ci        1000s digit
Ci           1: purge delimiter from sout (last char copied unles eos)
Ci   delim   string delimiters marking end of string (any character
Ci           in delim counts as a delimiter
Ci   is      first character to copy (string starts at character 1)
Co Outputs
Co   sout    strn copied to sout; see Remarks
Co   is      smaller of index in strn of last char copied and lstr+1
Co   io      position of terminator, or 1+ last character copied if
Co           limit of string reached beforehand
Cr Remarks
Cr   delimiters between " " or ' ' are excluded.  pairs " " and ' '
Cr   are excised in the output string sout
Cu Updates
Cu   01 Nov 01 delimiter may contain more than one character
C ----------------------------------------------------------------------
C     implicit none
      integer lstr,is,io,opt
      character *(*) strn, delim*(*), sout*(*)
      character ch*60
      integer i2,opt0,opt2,opt3,it,ia,nch,ldelim

      sout = ' '
      opt0 = mod(opt,10)
      opt2 = mod(opt/100,10)
      opt3 = mod(opt/1000,10)
      ldelim = len(delim)
      ch = delim // '"'' '
      nch = ldelim + 2
      if (opt2 .ne. 0) nch = nch+1
      is = is-1
      ia = 1
C ... Skip past leading blanks
      if (opt0 .ne. 0) then
        call skipbl(strn,lstr,is)
      endif

C ... Find i2 : points to past last char of the argument
      i2 = is
   12 is = i2
      call chrps2(strn,ch,nch,lstr,i2,it)
      if (it .le. ldelim .and. i2 .lt. lstr) i2 = i2+1
      if (i2 .gt. is) sout(ia:ia+i2-is-1) = strn(is+1:i2)
      if (it .le. 0 .and. i2 .eq. lstr) i2 = i2+1
      ia = ia+i2-is
C ... A blank encountered ... skip blanks and continue
      if (it .eq. ldelim+3 .and. i2 .lt. lstr) then
        call skipbl(strn,lstr,i2)
        sout(ia:ia) = ' '
        ia = ia+1
        goto 12
C ... A quote encountered ... find match and continue
      elseif (it .gt. ldelim .and. i2 .lt. lstr) then
        i2 = i2+1
        is = i2
        call chrpos(strn,ch(it:it),lstr,i2)
        call strncp(sout,strn,ia,is+1,i2-is)
        ia = ia+i2-is
        i2 = i2+1
        goto 12
      endif

      io = ia-1
      if (is .ge. lstr) io = io+1
C      if (is .ge. lstr) print *, 'hi',io
      if (opt3 .eq. 1) then
        if (i2 .le. lstr) sout(io:io) = ' '
      endif
      is = min(i2,lstr+1)

      end

      subroutine eostr(strn,lstr,opt,delim,is)
C- Mark end of a string
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn,is string, and first character (1st character starts at 1)
Ci   delim   string delimiter marking end of string
Ci   opt     1s  digit 1: exclude delimiters between " " or ' '
Ci           10s digit 1: skip initial blanks
Ci                     2: use strn(is) as delimiter;  delim is not used
Co Outputs
Ci   is      smaller of position in strn of delimiter and lstr+1
C ----------------------------------------------------------------------
C     implicit none
      integer lstr,is,opt
      character *(*) strn, delim*1
      character ch*3
      integer i2,opt0,opt1,it,k
      data ch /' "'''/


      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      ch(1:1) = delim
      if (opt1 .eq. 2) then
        ch(1:1) = strn(is:is)
      else
        is = is-1
      endif

C ... Skip past leading blanks
      if (opt1 .eq. 1) then
        call skipbl(strn,lstr,is)
      endif

C ... Find i2 : points to past last char of the argument
      k = 1
      if (opt0 .ne. 0) k = 3
      i2 = is
   12 is = i2

      call chrps2(strn,ch,k,lstr,i2,it)
      if (it .eq. 1 .and. i2 .lt. lstr) i2 = i2+1
C ... A quote encountered ... find match and continue
      if (it .gt. 1 .and. i2 .lt. lstr) then
        i2 = i2+1
        is = i2
        call chrpos(strn,ch(it:it),lstr,i2)
        i2 = i2+1
        goto 12
      endif

      is = min(i2,lstr+1)

      end
C      subroutine fmain
C      implicit none
C      integer is,io
C      character*80 strn,sout
C      strn = '  test" this string "now" and again"'
C      strn = '  test" this string "now" and again ",   and     '//
C     .  '" -until the very- "   very end'
C
C      print *, ' This is the entire string which will be copied:'
C      print '(a,a)', '1234567890123456789012345678901234567890'//
C     .  '1234567890123456789012345678901234567890',
C     .  ' (a ruler)'
C      print '(a)', strn
C
CC      is = 1
CC      call cpstr(strn,len(strn),1,' ',is,io,sout)
CCc     call eostr(strn,len(strn),10,' ',is)
CC      print 333, '80 1',is,io,sout(1:io)
C
C      print *, ' '
C      print '(''lstr opt del      is  io sout...'')'
C      print '(''--------------------------------'')'
C
C      is = 1
C      call cpstr(strn,11,1,'s',is,io,sout)
C      print 333, '11    1  s    ',is,io,sout(1:io)
C      print 334, strn(1:min(is,11))
C
C      is = 1
C      call cpstr(strn,11,1001,'s',is,io,sout)
C      print 333, '11 1001  s    ',is,io,sout(1:io)
C      print 334, strn(1:min(is,11))
C
C      is = 1
C      call cpstr(strn,5,1,'z',is,io,sout)
C      print 333, '5     1  z    ',is,io,sout(1:io)
C      print 334, strn(1:min(is,5))
C
C
C      is = 1
C      call cpstr(strn,11,1,'z',is,io,sout)
C      print 333, '11    1  z    ',is,io,sout(1:io)
C      print 334, strn(1:min(is,11))
C
C      is = 1
C      call cpstr(strn,11,1,' ',is,io,sout)
C      print 333, '11    1  blank',is,io,sout(1:io)
C      print 334, strn(1:min(is,11))
C
C      is = 1
C      call cpstr(strn,80,1,' ',is,io,sout)
C      print 333, '80    1  blank',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C      is = 1
C      call cpstr(strn,80,11,'wz',is,io,sout)
C      print 333, '80    1  zw   ',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C      is = 1
C      call cpstr(strn,80,11,'zv',is,io,sout)
C      print 333, '80   11  zv   ',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C      is = 1
C      call cpstr(strn,80,1011,'zv',is,io,sout)
C      print 333, '80 1011  zv   ',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C      is = 1
C      call cpstr(strn,80,11,'v,',is,io,sout)
C      print 333, '80   11  v,   ',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C      is = 1
C      call cpstr(strn,80,11,'v,o',is,io,sout)
C      print 333, '80   11  v,o  ',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C      is = 1
C      call cpstr(strn,80,1011,'v,o',is,io,sout)
C      print 333, '80 1011  v,o  ',is,io,sout(1:io)
C      print 334, strn(1:min(is,80))
C
C
C      print *, '-------- end of tests of cpstr -------------'
C      print *, ' '
C      print *, ' '
C      print *, '-------- test eostr -------------'
C
C      strn = '  test" eostr for delim / this string "now" and again"'//
C     .  '  ? '
C      print *, ' This is the string we will test'
C      print '(a,a)', '1234567890123456789012345678901234567890'//
C     .  '123456789012345678901234567890',
C     .  ' (a ruler)'
C      print '(a)', strn
C      print *, ' '
C      is = 1
C      call eostr(strn,80,11,'/',is)
C      print 335, is, ' ''', strn(1:min(is,80)), ''''
C      is = 1
C      call eostr(strn,80,11,'?',is)
C      print 335, is, ' ''', strn(1:min(is,80)), ''''
C      is = 1
C      call eostr(strn,80,11,'"',is)
C      print 335, is, ' ''', strn(1:min(is,80)), ''''
C
C      strn = '"  test eostr for delim / this string "now" and again"'//
C     .  '  ? '
C      is = 1
C      call eostr(strn,80,21,'"',is)
C      print 335, is, ' ''', strn(1:min(is,80)), ''''
C
C  333 format(a,2x,2i4,1x,'|',a)
C  334 format(25x,'|',a)
C  335 format(i5,a,a,a)
C
C      end

