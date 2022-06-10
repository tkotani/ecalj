integer function awrite(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
  !- Formatted output, with ascii conversion of binary numbers
  !i ifi: <>0, local output string written to abs(ifi);
  !i      <=0, sout copied to local output string initially;
  !i           local output string copied back to sout on exit
  !i       >0, sout unaltered on exit
  !i mxln: abs(mxln) = maximum number of characters to copy
  !i mxln < 0: suppress trailing blanks when writing to logical unit.
  !o sout:     output string (see ifi, above)
  !r Characters are copied from fmt to the output string sout, which is
  !r   then (optionally) written to logical unit ifi.  Pointer ip keeps
  !r   track of the current position for writing to sout.  Copy
  !r   is literal except when a control char % is encountered.
  !r   % characters do one of several functions:
  !r   %l writes to sout an ascii representation of logical argument a_j
  !r      (NB: j is 1 for first conversion, 2 for second, etc).
  !r   %i writes an integer argument a_j
  !r   %d, %e, %g, %G, %D and %F write ascii rep of double precision a_j
  !r     'd' writes in decimal notation
  !r     'e' writes in exponential notation
  !r     'g' and 'G' take the minimum size of 'd' and 'e'; 'g' is to
  !r                 specify relative precision, 'G' absolute precision.
  !r     All of the above generate 'pretty' ascii representations
  !r     (see pretty.f)
  !r     'D' mimics the fortran 'f' format and is intended for output in
  !r     fixed columns.
  !r     'F' puts the number in a specified space, using whatever form
  !r     produces the most decimal places of precision.
  !r   %% quotes a "%" literally
  !r   %a shifts ip past last nonblank character
  !r   %f shifts ip forward
  !r   %b shifts ip backward
  !r   %p sets   ip to a fixed value
  !r   %t is obsolete
  !r   %x blanks the output string
  !r   %z can suppress leading the leading zero in a decimal fraction
  !r   %W shifts ip forward until a whitespace is encountered
  !r   %w shifts ip forward until a non-whitespace is encountered
  !r   %c closes up whitespace around ip
  !r   %o opens  up whitespace around ip
  !r   %u if numerical argument = NULLI output 'null' instead of number
  !r      Optional argument n1:
  !r      0 turn off null option, for this and future calls
  !r      1 (or default) set null option, for this and future calls
  !r     >1 set null option for this call only
  !r     <0 turn off null option for this call only
  !r   %? conditionally parses one of two strings (see below)
  !r   %j increments argument jumps over call arguments
  !r   %N is turned into a newline, calling nlchar to get newline
  !r Most control characters have optional arguments.
  !r For d,e,g,G,D,F,l,i the general syntax is:
  !r   %[n1][:n2][,n3][;n4][#n5]x, with x=d,e,g,G or F
  !r Here n1..n5 are integer expressions:
  !r   n1 number of values to convert (a_j is regarded as a vector)
  !r   n2 number of blank spaces preceding first character
  !r      n2<0 => subtract one space if argument is negative
  !r   n3 minimum number of digits to display (after '.' for 'd'
  !r      and 'G' and total number for 'e' and 'g')
  !r   n4 round to n4 decimal places (absolute for 'd' and 'G',
  !r      and relative for 'e' and 'g')
  !r   n5 if conversion uses less than n5 characters, append trailing
  !r      blanks to fill (used for lining data in columns)
  !r For D the meanings of n2..n4 are different:
  !r   n2 is not used
  !r   n3 number of digits after decimal
  !r   n4 field width
  !r For F:
  !r   n2 is not used
  !r   n3 is not used
  !r   n4 is the field width
  !r For l and i:
  !r   n3 is the field width
  !r For z, j, p, a, f, o, b, and and the general syntax is:
  !r   %[n1]x, with x=z, p, a, f, b, u
  !r   n1 repeats (f, b)
  !r   n1 1=>suppresses leading 0, 0=>ensures its presence (z)
  !r   For u, see above
  !r NB: there is an option to substitute for any of n1..n4 one of the
  !r arguments a_j.  This is done by using the character 'n' instead
  !r of some integer expression (eg %n:n,5d).  awrite uses the
  !r next argument a_j is used for n, and increments j.  Thus, %n:n,5d
  !r consumes the next three a_j, the first describing the number
  !r of elements to write, the second the number of spaces between
  !r arguments.
  !r For ? the general syntax is
  !r   %?QexprQstr1Qstr2Q
  !r   str1 is parsed if "expr" evaluates to nonzero, otherwise str2 is.
  !r   Q is some character, eg ';'.  It should NOT be some character
  !r   that may be confused as part of "expr", like '?' or '/'.
  !r   The next argument argument a_j is temporarily set to an integer
  !r   value and temporarily named `n', which may be used in 'expr'.
  !r   Also the current number of characters in the string is temporarily
  !r   assigned to `p'.  Finally, as a special case for a expression
  !r   involving strings, the following is permitted:
  !r     %c==X
  !r   where X is some character.  This expression evaluates to nonzero
  !r   if the character output string at the current position is equal
  !r   to X; otherwise it evaluates to zero.
  !r   Example:
  !r     call awrit2('three plus one is %?;n==1;%i;four;, no?',mxlen,
  !r                 s,mxlen,-i1mach(2),m,4)
  !r     prints out "three plus one is 4, no?" if m equals 1; otherwise
  !r     prints out "three plus one is four, no?"
  !u Updates
  !u   13 Oct 07 Modified lnull, for permanent option
  !u   02 Aug 07 Added %u: outputs null string when arg=NULLI (bin2av)
  !u   27 Feb 02 Added %c==X type of conditional expression
  !u    8 May 01 addition of n5 modifier described above
  ! ----------------------------------------------------------------
  implicit none
  integer :: ifi,mxln
  character*(*) fmt,sout
  double precision :: a1(*),a2(*),a3(*),a4(*),a5(*),a6(*),a7(*),a8(*)
  integer :: i,lfmt,ip,jp,cast,ia,iff,i2,iterm,ndec,iv(5),nblk,nx,j, &
       mxl,ls,icond,ivawrt,lens,nsyv,ires,fw
  equivalence (iv(1),i2),(iv(2),nblk),(iv(3),nx),(iv(4),ndec), (iv(5),fw)
  logical :: a2bin,ltmp,lnull,lnulls
  double precision :: holdn,holdp,xx
  parameter (lens=1024)
  character*(lens) s,ss,fchr*29,fm*20,cc*1,ccond*1
  save lnulls
  data fchr /' :;,#irdegGltapfbzDwWcoxu?jNF'/
  data lnulls /.false./
  logical:: l_exec
  ! ... ia,iff,ip: indices to current argument, pos in fmt, pos in s
  mxl = iabs(mxln)
  s = ' '
  if (ifi <= 0) s = sout 
  ip = 0  ! ... ip holds the current position in the output string
  iff = 0 ! ... iff holds the current position in the format string
  ia = 1 ! ... index to current argument in the argument list
  icond = 0 ! ... icond nonzero when in the middle of a conditional expression
  ccond = ' ' ! ... ccond is the terminating character for conditional expression
  call numsyv(nsyv) ! ... hold on to n,p in vars table; we need them as local variables
  call getsyv('p',holdp,j)
  call getsyv('n',holdn,j)
  lfmt = len(fmt)
  ls = len(s)
  lnull = lnulls
 ! --- Parse next character in fmt ---
19 ia = ia-1
20 iff = iff+1
  !  ...  End of fmt
  if (iff > lfmt) goto 10
  !  ...  Character terminating conditional string
  if (icond > 0 .AND. fmt(iff:iff) == ccond) then
     if (icond == 2) then
        call chrpos(fmt,ccond,lfmt,iff)
        iff = iff+1
     endif
     icond = 0
     goto 20
  endif
  !  ...  Any non-% character
  ! ino      if (fmt(iff:iff) .ne. '%' .or. fmt(iff:iff+1) .eq. '%%') then
  l_exec=.false.
  if (fmt(iff:iff) /= '%') l_exec= .TRUE. 
  if (iff+1<=len(fmt)) then
     if (fmt(iff:iff+1) == '%%')  l_exec= .TRUE. 
  endif
  if (l_exec) then
     !         print *, 'parsing non-%:', iff, fmt(1:iff)
     ip = ip+1
     if (ip <= min(ls,mxl)) s(ip:ip) = fmt(iff:iff)
     if (iff+1<=len(fmt)) then
        if (fmt(iff:iff+1) == '%%') iff = iff+1
     endif
     goto 20
     !   --- Parse % ---
  else
     !     ... Default values for %command  !  print *, 'now parsing %:', iff, fmt(iff:min(lfmt,iff+10))
     nblk = 0
     fw = 0
     nx = 99
     ndec = 0
     i2 = 1
     ia = ia+1
     iff = iff+1
     !     ... iterm flags whether cc is ':;,#', to use later
     j = 0
     call chrps2(fmt(iff:iff),fchr,len(fchr),0,j,iterm)
     !     ... Re-entry if to parse another argument
25   if (iterm >= 2 .AND. iterm <= 5) iff = iff+1
     cc = fmt(iff:iff)
     ires = 1 ! ires is integer argument; set to 1 for default value
     if (cc == '(') then ! Check for an integer expression () preceding command:
        j = 1
        do  35  i = iff+1, lfmt
           if (fmt(i:i) == '(') j=j+1
           if (fmt(i:i) == ')') j=j-1
           if (j == 0) then
              xx = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
              call lodsyv('n',1,xx,j)
              call lodsyv('p',1,dble(ip),j)
              j = iff-1
              ltmp = .not. a2bin(fmt,ires,2,0,' ',j,i-1)
              call rxx(ltmp,'awrite: failed to parse () in format')
              iff = i+1
              goto 36
           endif
35      enddo
        call rx('awrite: missing matching () in format')
36      continue
     else if (cc == 'n') then ! ... Prior integer argument if next char 'n':
        ires = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
        ia = ia+1
        iff = iff+1
        !     ... Prior integer argument if char an integer:
     else if (cc >= '0' .AND. cc <= '9' .OR. cc == '-') then
        do  22  i = iff, lfmt-1
           j = i+1
           if (fmt(j:j) >= '0' .AND. fmt(j:j) <= '9') goto 22
           j = iff-1
           !             call pshpr(130)
           ltmp = .not. a2bin(fmt,ires,2,0,fmt(i+1:i+1),j,i)
           call rxx(ltmp,'awrite: failed to parse format')
           iff = i+1
           goto 23
22      enddo
23      continue
     endif
     !          print 335, ires,iff,fmt(1:iff)
     ! 335     format('*now ires=',i2,' parsed to ',i3,': ',a)
     !     ... If this was an argument to one of ':;,#'
     if (iterm >= 2 .AND. iterm <= 5) then
        iv(iterm) = ires
        !     ... Otherwise ires is an argument to command
     else
        iv(1) = ires
     endif
     cc = fmt(iff:iff)!     ... Next character is the terminator
     j = 0
     call chrps2(cc,fchr,len(fchr),0,j,iterm)
     !     ... If an argument, run through parse again
     if (iterm >= 2 .AND. iterm <= 5) goto 25
     !     ... Otherwise a command:
     cast = 99
     cc = fmt(iff:iff)
     if (cc == 'l') cast=0
     if (cc == 'i') cast=2
     if (cc == 'r') cast=3
     if (cc == 'd' .OR. cc == 'e' .OR. cc == 'D' .OR. &
          cc == 'g' .OR. cc == 'G' .OR. cc == 'F') &
          cast=4
     if (cc == 't') then
        call rx('awrite: use p, not t')
     endif
     if (cc == 'z') then
        call bin2a0(i2)
        goto 19
     elseif (cc == 'j') then
        ia = ia+i2
        goto 19
     elseif (cc == 'a') then
        call skpblb(s,ls,ip)
        ip = ip+i2
        goto 19
     elseif (cc == 'p') then
        ip = i2
        goto 19
     elseif (cc == 'f') then
        ip = ip+i2
        goto 19
     elseif (cc == 'b') then
        ip = ip-i2
        goto 19
     elseif (cc == 'N') then
        call nlchar(1,s(ip+1:ip+1))
        ip = ip+1
        goto 19
        ! ---     Entry point for conditional expression ---
     elseif (cc == '?') then
        if (icond /= 0) call rx('awrite encountered nested "%?"')
        icond = 1
        call lodsyv('p',1,dble(ip),j)
        xx = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
        call lodsyv('n',1,xx,j)
        ia = ia+1
        iff = iff+1
        !     ...   ccond is character terminating conditional string
        ccond = fmt(iff:iff)
        !     ...   If next char is '%', expression is of the string type:
        if (fmt(iff+1:iff+1) == '%') then
           if (fmt(iff+1:iff+4) == '%c==') then
              ltmp = fmt(iff+5:iff+5) .eq. s(ip:ip)
              iff = iff+6
           else
              call rxs('awrite: failed to parse : ',fmt(iff:))
           endif
           !     ...   Parse expression
        elseif ( .NOT. a2bin(fmt,ltmp,0,0,ccond,iff,lfmt)) then
           call rx('awrite: failed to parse conditional expr')
        endif
        !     ...   Use first string, or skip to second string
        if (ltmp) then
           icond = 2
        else
           icond = 3
           call chrpos(fmt,ccond,lfmt,iff)
           iff = iff+1
        endif
        goto 19
        !     ... clear string
     elseif (cc == 'x') then
        s = ' '
        goto 19
        !     ... toggle on 'null' option
     elseif (cc == 'u') then
        if (i2 < 0) then
           lnull = .false.
        elseif (i2 == 0) then
           lnulls = .false.
           lnull = .false.
        elseif (i2 == 1) then
           lnulls = .true.
           lnull = .true.
        else
           lnull = .true.
        endif
        goto 19
        ! ...     pad whitespace around ip
     elseif (cc == 'o') then
        if (ip == 0) ip = 1
        ss = s(ip:ls)
        s(ip:ip+i2-1) = ' '
        s(ip+i2:ls) = ss
        ip = ip+i2-1
        goto 19
        ! ...     close up whitespace around ip
     elseif (cc == 'c') then
        do  13  j = ip, 1, -1
           if (s(j:j) /= ' ' .AND. s(j:j) /= '        ') goto 14
           ip = j
13      enddo
14      continue
        jp = ip-1
        do  15  j = jp+1, mxl
           if (s(j:j) /= ' ' .AND. s(j:j) /= '        ') goto 16
           jp = j
15      enddo
16      continue
        if (jp-ip+1 > 0) then
           ss = s(jp+1:ls)
           s(ip:ls) = ss
        endif
        goto 19
        ! ...     skip to next nw
     elseif (cc == 'w') then
        do  17  j = ip, mxl
           if (s(j:j) /= ' ' .AND. s(j:j) /= '        ') goto 19
           ip = j
17      enddo
        ! ...     skip to next whitespace
     elseif (cc == 'W') then
        do  18  j = ip, mxl
           if (s(j:j) == ' ' .OR. s(j:j) == '        ') goto 19
           ip = j
18      enddo
     endif
  endif
  if (cast == 99) call rx('awrite: unknown control: ' // cc)
  ! ---   Generate format for bin2a ---
  fm = ' '
  if (cast == 4) then
     fm = cc
     if (cc == 'G') fm = 'g'
     j = 1
     if (nx /= 99) call bin2a(' ',0,0,nx,2,0,20,fm,j)
     if (cc == 'G') call bin2a(':20',0,0,0,1,0,20,fm,j)
  endif
  ! ---   Convert binary numbers ---
  i2 = i2-1
  if (ia == 1) call bin2av(fm,fw,nblk,ndec,a1,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 2) call bin2av(fm,fw,nblk,ndec,a2,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 3) call bin2av(fm,fw,nblk,ndec,a3,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 4) call bin2av(fm,fw,nblk,ndec,a4,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 5) call bin2av(fm,fw,nblk,ndec,a5,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 6) call bin2av(fm,fw,nblk,ndec,a6,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 7) call bin2av(fm,fw,nblk,ndec,a7,cast,0,i2,' ',mxl,lnull,s,ip)
  if (ia == 8) call bin2av(fm,fw,nblk,ndec,a8,cast,0,i2,' ',mxl,lnull,s,ip)
  goto 20
10 continue
  ip = min(ip,mxl)
  if (mxln < 0) then
     call skpblb(s,ip,ip)
     ip = ip+1
  endif
  ia = iabs(ifi)
  if (ifi /= 0 .AND. ip > 0) write(ia,333) s(1:ip)
333 format(a)
  if (ifi <= 0 .AND. ip > 0) sout = s(1:ip)
  awrite = ip
  call lodsyv('p',1,holdp,j) !! --- Restore or undo symbolic variables p,n ---
  call lodsyv('n',1,holdn,j)
  call clrsyv(nsyv)
  return
end function awrite

subroutine bin2av(fmt,w,nblk,ndec,res,cast,i1,i2,sep,mxln,lnull, &
     outs,ip)
  !- Write out a vector of of numbers using bin2a
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   fmt   :format passed to bin2a
  !i    w    :unused if zero.  If >0,
  !i         :w = minimum spacing between successive numbers
  !i   nblk  :number of blanks preceding each value, or if nblk < 0,
  !i         :|nblk| spaces are prepended for positive numbers
  !i         :|nblk|-1 spaces are prepended for negative numbers
  !i   ndec  :retain a mininimum ndec digits after decimal (see bin2a)
  !i   res   :vector of binaries to convert to ascii string
  !i   cast  :0=logical, 1=char, 2=int, 3=real, 4=double
  !i   i1    :convert numbers res(i1..i2)
  !i   i2    :convert numbers res(i1..i2)
  !i   sep   :separator between numbers
  !i   mxln  :maximum allowed value of ip
  !i   ip    :string position pointer
  !i   lnull :if T, numbers equal to NULLI are turned into NULL
  !o Outputs
  !o   outs  :string containing ascii rep'sn of binary numbers
  !i   ip    :string position pointer updated to end of string
  !r Remarks
  !u Updates
  !u   01 Aug 07 new lnull
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) fmt,outs,sep*1
  integer :: nblk,cast,i1,i2,ip,mxln,ndec,w
  double precision :: res(0:1)
  logical :: lnull
  ! ... Local parameters
  logical :: lneg,llnull
  integer :: nblk2,ival,i,k,ip0,NULLI
  real :: rval
  double precision :: dval
  parameter (NULLI=-99999)

  if (mxln <= 0) return
  do  i = i1, i2
     nblk2 = nblk
     if (nblk < 0) then
        nblk2 = -nblk
        lneg = .false.
        if (cast == 2) lneg = ival(res,i+1) < 0
        if (cast == 3) lneg = rval(res,i+1) < 0
        if (cast == 4) lneg = dval(res,i+1) < 0
        if (lneg) nblk2 = nblk2-1
     endif
     !       Set flag llnul if lnull is ON and argument matches NULLI
     llnull = .false.
     if (lnull) then
        if (cast == 2) llnull = ival(res,i+1) == NULLI
        if (cast == 3) llnull = rval(res,i+1) == NULLI
        if (cast == 4) llnull = dval(res,i+1) == dble(NULLI)
        if (llnull) then
           call skpblb(fmt,len(fmt),ip0)
           fmt(2+ip0:) = ':n'
        endif
     endif
     ip0 = ip
     call bin2a(fmt,nblk2,ndec,res,cast,i,mxln,outs,ip)
     if (llnull) then
        !         If fixed width, leave position of null as is
        if (fmt(1:1) == 'D' .OR. fmt(1:1) == 'F' .OR. &
             (cast == 2 .AND. ndec > 0)) then
           !         Skip if not sufficient space for leading blanks + null
        else if (ip-3 <= 1+ip0+iabs(nblk)) then
           !         Otherwise rewrite null starting at 1+ip0+iabs(nblk)
        else
           outs(1+ip0+iabs(nblk):ip) = 'NULL'
           ip = 4+ip0+iabs(nblk)
        endif
        !         print *, outs(1:ip)
     endif
     if (sep /= ' ' .AND. i < i2) then
        ip = ip+1
        outs(ip:ip) = sep
     endif
     if (w /= 0) then
        do  k = ip+1, ip0+w
           outs(k:k) = sep
           ip = ip+1
        enddo
     endif
  enddo
end subroutine bin2av

subroutine awrit8(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
  !- Subroutine versions of integer function awrite
  !     implicit none
  double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
  character*(*) sout,fmt
  integer :: ifi,mxln,ip,jp,awrite
  save ip
  entry awrit7(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7)
  entry awrit6(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6)
  entry awrit5(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5)
  entry awrit4(fmt,sout,mxln,ifi,a1,a2,a3,a4)
  entry awrit3(fmt,sout,mxln,ifi,a1,a2,a3)
  entry awrit2(fmt,sout,mxln,ifi,a1,a2)
  entry awrit1(fmt,sout,mxln,ifi,a1)
  entry awrit0(fmt,sout,mxln,ifi)
  ip = awrite(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
  return
  entry awrip(jp)
  jp = ip
end subroutine awrit8

subroutine vwrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8,cast,ires,res)
  !- Writes either integer or double into ires or res, depending on cast
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ia    :indicates which of arrays a1..a8 to extract element from
  !i   n     :which entry in array a_ia
  !i   a1..a8:element is extracted from one of these arrays
  !i   cast  :array cast
  !o Outputs
  !o   ires  :if cast is integer, result poked into ires
  !o   res   :if cast is double, result poked into res
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: ia,n,cast,ivawrt,ires
  double precision :: dvawrt,res
  double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
  if (cast == 2) then
     ires = ivawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
  elseif (cast == 4) then
     res = dvawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
  else
     call rxi('vwrt: cannot handle cast',cast)
  endif
end subroutine vwrt

integer function ivawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
  !     implicit none
  integer :: ia,n,ival
  double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
  ivawrt=99999
  if (ia == 1) ivawrt = ival(a1,n)
  if (ia == 2) ivawrt = ival(a2,n)
  if (ia == 3) ivawrt = ival(a3,n)
  if (ia == 4) ivawrt = ival(a4,n)
  if (ia == 5) ivawrt = ival(a5,n)
  if (ia == 6) ivawrt = ival(a6,n)
  if (ia == 7) ivawrt = ival(a7,n)
  if (ia == 8) ivawrt = ival(a8,n)
end function ivawrt

real(8) function dvawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
  !     implicit none
  integer :: ia,n
  double precision :: a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
  dvawrt=1d99
  if (ia == 1) dvawrt = a1(n)
  if (ia == 2) dvawrt = a2(n)
  if (ia == 3) dvawrt = a3(n)
  if (ia == 4) dvawrt = a4(n)
  if (ia == 5) dvawrt = a5(n)
  if (ia == 6) dvawrt = a6(n)
  if (ia == 7) dvawrt = a7(n)
  if (ia == 8) dvawrt = a8(n)
END function dvawrt

subroutine bin2a(fmt,nblk,ndec,res,cast,count,mxlen,outstr,ip)
  !- Converts number to ascii format, stripping leading blanks, trailing 0
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   fmt: cast=1:   holds the string to be appended to outstr)
  !i        cast=0,2: not used
  !i        cast=3,4: syntax X[#][:sw] where X is one of
  !i                  'd' to write in decimal representation
  !i                  'e' to write in exponential format
  !i                  'g' to use the smaller of 'd' and 'e'
  !i                  '(n.m)' fixed format, mimics fortran fmt (fn.m)
  !i                  'D' also mimics fortran fmt (Fn.m)
  !i                      D# => supplies n=#; arg ndec supplies m
  !i                  'F' fixed field, picking between d and e that
  !i                      F# => # is field width
  !i                      generates the most significant digits
  !i                  See Remarks for further description
  !i   nblk:  strip leading blanks, leaving a maximum of nblk
  !i   ndec:  (cast=3,4 only) retain a mininimum ndec digits after decimal
  !i          point, i.e. do not suppress trailing zeros to ndec.
  !i          ndec=0 does nothing.  ndec should not exceed precsn.
  !i          (cast=2 only): ndec specifies a field width
  !i   res:   binary value to be converted into ascii string
  !i   cast:  cast of res: 0=logical, 1=char, 2=int, 3=real, 4=double
  !i   count: res(count) is to be converted.  NB: count=0 for first entry
  !i   mxlen: maximum length of outstr
  ! o Inputs/Outputs
  ! o  ip:    on input, starting position in outstr for write
  ! o         NB: ip=0 points to first character in string
  ! o  ip:    on output, position of final character written to outstr
  !o  Outputs
  !o   outstr:binary res(count) written in ascii form to outstr(ip:..)
  !r Remarks
  !r  *The string representation of floating point numbers is generated
  !r   by a "prettified" modification of the fortran write statement
  !r   (pretty.f), which includes suppression of trailing zeros and the
  !r   option to include or suppress the leading zero in decimal
  !r   fractions less than 1.  Floating-point formats include:
  !r     'd[n][:sw]' for decimal representation,
  !r     'e[n][:sw]' for exponential representation,
  !r     'g[n][:sw]' uses the minimum length of 'd' and 'e'
  !r     'D[n][:sw]' simulates the standard fortran format fn.m
  !r                 Here n follows D, ndec the role of m.  Or:
  !r     'Fn'        fixed field, picking between d and e that generates
  !r                 the most significant digits
  !r      (n.m)      also simulates the standard fortran format.
  !r
  !r  *Optional modifier 'n' is a number specifying how many decimals of
  !r   precision (n=6 if not specified). By default, n means:
  !r      for 'd' format, the absolute precision: i.e.
  !r        number of digits after the decimal point
  !r     for 'e' format, the relative precision , i.e.
  !r        number of digits printed
  !r     for 'D' format, it is the field width n in fortran format fn.m
  !r  *Optional modifier sw is a compound of the 1's and 10's digits.
  !r       1's digit of sw can overwrite the default meaning of 'n' above.
  !r                 sw=0 => n corresponds to absolute precision
  !r                 sw=1 => n corresponds to relative precision
  !r       10's digit nonzero suppresses leading blanks.
  !r  *Entry bin2a0 allows the user to set the default of sw.
  !r  *Examples:
  !r     call bin2a('d2',1,3,1.234951d0,...)    => 1.23
  !r     call bin2a('d4',1,4,1.234951d0,...)    => 1.2350
  !r     call bin2a('d3:11',1,0,1.2349501d-6,4) => .00000123
  !r     call bin2a('e2',1,3,1.2499d7,...)      => 1.2e7
  !r     call bin2a('e5',1,5,1.2349510d7,...)   => 1.2350e7
  !r     call bin2a('e5:0',1,4,1.2349501d5,...) => 1.234950100e5
  !r     call bin2a('g',1,0,1.23d-5,...)        => 1.23e-5
  !r     call bin2a('g3:10',1,3,1.24996d-5,...) => .000
  !r     call bin2a('g4:10',1,4,1.24996d-5,...) => 1e-5
  !r     call bin2a('f4:10',1,4,1.24996d-5,...) => 1e-5
  !u Updates
  !u   02 Aug 07 Added :n outputs null string when res=NULLI
  ! ----------------------------------------------------------------------
  implicit none
  ! Passed Parameters

  !c kino's correctio for ifort was
  !c     character(mxlen):: outstr ! ?---> !character(*) can not check size of outstr.
  !c However, because of a bug in grortran4.3.4, this is not allowed. Thus I now use character(*).
  !c Sep2010
  character(*):: outstr
  !c
  character(*):: fmt
  double precision :: res(0:*)
  integer :: nblk,cast,count,ip,mxlen,ndec,is
  ! Local Variables
  logical :: lD,lS,lF,lnull,llnull,parstr
  integer :: i,j,k,iprint,lsmx,n1,n2,np,precsn,fw, &
       ix(4),iv(4),a2vec,p,isw,isw0,getdig,m,ndig,ndige
  parameter (lsmx=80)
  character(20) :: lfmt, strn*(lsmx), strn2*(lsmx), ss*(lsmx)
  character:: fm*20
  real :: rval
  double precision :: xx
  integer :: NULLI
  parameter (NULLI=-99999)
  save isw0
  data isw0 /0/

  !     write(*,"('enter bin2a: cast,fmt=',i4,1x,a$)") cast,fmt

  ! --- Convert binary to ascii representation (log, int, char) ---
  lnull = .false.
  llnull = .false.
  goto (10,11,12,20,20), cast+1
  call rx('bin2a: bad cast')
10 continue
  call bin2al('(L8)',res,count,strn2)
  goto 15
11 continue
  strn2 = fmt
  goto 15
12 continue
  call bin2ai('(I16)',res,count,strn2,lnull)
  goto 15
  ! --- copy strn2 to strn with appropriate number of blanks ---
15 continue
  i = 0
  call skipbl(strn2,lsmx,i)
  strn = ' '
  !     If a field width specified, overwrite spillover with '*'
  if (ndec /= 0) then
     call skpblb(strn2,lsmx,j)
     j = j-ndec+1
     if (j > i) then
        strn2(j+1:j+ndec) = '****************'
     endif
     i  = j
  endif
  strn(1+nblk:lsmx) = strn2(i+1:lsmx)
  call skpblb(strn,lsmx,n1)
  n1 = n1+1
  if (lnull .AND. fmt /= ' ') then
     i = 0
     if (parstr(fmt,':n',len(fmt)-1,2,'n',i,j)) then
        llnull = .true.
     endif
  endif
  goto 50

  ! --- Entry for setting up or determinining defaults ---
  entry bin2a0(is)
  if (is >= 0) isw0 = is
  if (is < 0) is = isw0
  return

  ! --- Binary->ascii representation, floating-point ---
20 continue
  if (cast == 3) xx = rval(res,count+1)
  if (cast == 4) xx = res(count)
  lnull = xx .eq. dble(NULLI)

  ! ... Determine appropriate format
  lfmt = fmt
  i = 0
  call skipbl(fmt,len(fmt),i)
  if (i >= len(fmt)) then
     lfmt = 'g'
  else
     lfmt = fmt(i+1:len(fmt))
  endif
  i = 0
  if (parstr(lfmt,':n',len(lfmt)-1,2,'n',i,j)) then
     lfmt(i+1:) = ' '
     llnull = .true.
  endif
  ! --- Do the conversion, floating point ---
  if (lfmt(1:1) == '(') then
     write(ss,lfmt) xx
     call pretty(ss,nblk,ndec,20,isw0,strn,n1)
  else
     strn  = ' '
     strn2 = ' '
     lD = .false.
     lF = .false.
     j = 0
     !   ... i=1 => 'd'  i=2 =>  'e'  i=3 => 'g'
     call chrps2(lfmt,'degDF',5,len(lfmt),j,i)
     if (i <= 0) call rx('bin2a: bad format: '//lfmt)
     if (i == 5) then
        i = 3
        lF = .true.
     elseif (i == 4) then
        i = 1
        lD = .true.
     endif
     !   ... Get precsn (or field width for D or F), in iv(1), sw in iv(2)
     j = j+1
     np = a2vec(lfmt,len(lfmt),j,2,': ',2,2,2,ix,iv)
     isw = 1 + isw0
     if (i == 1) isw = 0 + isw0
     !   ... Simulated fortran format: precsn dictated by ndec
     if (lF) then
        if (np <= 0) call rx('bin2a: bad format: '//lfmt)
        fw = iv(1)
     elseif (lD) then
        precsn = ndec
        fw = -1
        if (np >= 1) fw = iv(1)
        !   ... if precsn explicit, use it
     elseif (np >= 1) then
        precsn = iv(1)
        !   ... This is the default
     else
        precsn = 6
     endif
     if (np >= 2) isw = iv(2)
     if (isw >= 20) isw = mod(isw,10) + isw0
     !  21   continue
     !   ... p is the exponent
     p = 0
     if (xx /= 0) then
        p = int(dlog10(dabs(xx)))
        if (dabs(xx) < 1) p = p-1
     endif
     !   ... fortran 'f' format
     if (i == 1 .OR. i == 3) then
        !     ... Estimate total width of format statement for fortran write
        if (lF) then
           !       ... m is the space consumed by a '-' sign
           m = (1-int(dsign(1d0,xx)))/2
           !       ... precsn = # rhs dec = field width - '-' - '.' - (p+1)
           precsn = fw - m - 1 - max(p+1,1)
           !       ... Only works on some compilers
           !            if (mod(isw,10) .ne. 0)
           !     .      precsn = fw - m - 1 - max(p+1,0)
           !       ... ndig = how many nonzero decimals printed
           ndig = max(precsn+p+1,0)
           !       ... Exclude 'f' if it doesn't fit
           if (precsn < 0) then
              ndig = -1
              !       ... Exclude 'e' if it does, and number isn't small
           else if (p > -2) then
              i = 1
           endif
           !       ... Determine how many digits we get from 'e' format
           if (i /= 1) then
              write(ss,'(1pe20.0)') xx
              !             print *, ss
              !         ... We want at least 1 more digit than f format
              call pretty(ss,0,max(ndig+1,1),max(ndig+1,1),1,strn,j)
              !         ... Tack on trailing 'e0' if pretty discarded it
              k = 0
              call chrpos(strn,'e',j,k)
              if (k >= j) j = j+2
              !         ... How many decimals for 'e' format
              ndige = max(ndig+1,1) + fw - j
              !         ... If pretty suppresses '.', add it back if ndige>1
              if (ndige > 1) then
                 k = 0
                 call chrpos(strn,'.',j,k)
                 if (k >= j) ndige=ndige-1
              endif
              !             print *, strn
           else
              ndige = ndig-1
           endif
           !       ... Generate string for F format here.
           if (ndig < 0 .AND. ndige < 0) then
              strn = ' '
              strn(nblk+1:nblk+fw) = '********************************'
              n1 = fw+nblk
              goto 50
           else if (ndig >= ndige) then
              i = 1
           else
              i = 2
              precsn = ndige
              goto 35
           endif
        elseif ( .NOT. lD .OR. (lD .AND. fw == -1)) then
           fw = max(p+3,5) + precsn
           if (getdig(isw,0,10) == 1) &
                fw = max(p+3,3) + max(precsn-p-1,0)
           fw = max(fw,10)
           if (fw > min(lsmx-2,99)) then
              strn = ' '
              strn(nblk+1:nblk+1) = '*'
              n1 = 1+nblk
              goto 35
           endif
        endif
        j = fw
        !     ... Insert leading blanks
        !         if (lF) then
        if (lF .OR. lD) then
           j = j+nblk
        endif
        if (j >= 10) write(fm,'(''(f'',i2,''.'')') j
        if (j < 10) write(fm,'(''( f'',i1,''.'')') j
        k = j
        !     ... Number of decimals for fortran write
        j = precsn
        if ( .NOT. (lD .OR. lF)) then
           if (getdig(isw,0,10) == 1) j = precsn-p-1
           j = max(j,0)
        endif
        !         decimals can't exceed field width - 1
        j = max(min(k-1,j),0)
        if (j >= 10) write(fm(6:8),'(i2,'')'')') j
        if (j < 10) write(fm(6:7),'(i1,'')'')') j
        write(ss,fm) xx
        if (lD .OR. lF) then
           if (nblk <= 0) then
           elseif (ss(1:nblk) /= ' ') then
              ss(1:k) = '*****************************************************'
              ss(1:nblk) = ' '
           endif
           strn = ss
           call skpblb(strn,lsmx,n1)
           n1 = n1+1

           k = 0
           call chrps2(strn,'-.0123456789',12,n1,k,j)
           j = j-1
           lS = j .eq. 0
           if (lS) call chrps2(strn,'.0123456789',11,n1,k,j)
           !     ...   Case fraction should have a leading '0'
           if (j == 1 .AND. getdig(isw,1,10) == 0) then
              if (lS .AND. k > 1) strn(k-1:k) = '-0'
              if ( .NOT. lS .AND. k > 0) strn(k:k) = '0'
              !     ...   Case fraction should have no leading '0'
           elseif (j == 2 .AND. getdig(isw,1,10) /= 0) then
              if (lS) strn(k:k+1) = ' -'
              if ( .NOT. lS)strn(k+1:k+1) = ' '
           endif
        else
           !           print *, 'before pretty ...', fm, ss
           call pretty(ss,nblk,ndec,precsn,isw,strn,n1)
        endif
35      continue
     endif
     !    .. fortran 'e' format
     if (i == 2 .OR. i == 3) then
        j = p + precsn
        if (getdig(isw,0,10) == 1) j = precsn-1
        if (j > 22) then
           strn2 = ' '
           strn2(nblk+1:nblk+1) = '*'
           n2 = 1+nblk
           goto 45
        endif
        j = min(max(j,0),99)
        if (j >= 10) write(fm,'(''(1pe30.'',i2,'')'')') j
        if (j < 10) write(fm,'(''(1pe30.'',i1,'')'')') j
        write(ss,fm) xx
        !         print *, 'before pretty ...', fm, ss
        j = ndec
        if (lF) j = precsn
        call pretty(ss,nblk,j,precsn,isw,strn2,n2)
        !     ... Tack on trailing 'e0' if pretty discarded it
        j = 0
        call chrpos(strn2,'e',n2,j)
        if (j >= n2 .AND. i == 2) then
           strn2(n2+1:n2+2) = 'e0'
           n2 = n2+2
        endif
        !     ... Sometimes the '.' is suppressed; make fw right
        if (lF .AND. n2 < fw+nblk) n2 = fw+nblk
45      continue
     endif
     if (i == 2 .OR. i == 3 .AND. &
          (n2 < n1 .OR. strn(nblk+1:nblk+1) == '*')) then
        strn = strn2
        n1 = n2
     endif
  endif
  ! --- Copy to outstr ---
50 continue
  n1 = max(n1,0)
  n2 = max(min(n1,mxlen-ip),0)
  !     Handle null number: replace ascii string with 'NULL'
  if (lnull .AND. llnull) then
     strn(1:n2) = ' '
     i = max(n2-3,1+nblk)
     strn(i:n2) = 'NULL'
  endif
  if (n2 > 0) outstr(ip+1:ip+n2) = strn(1:n2)
  ip = ip+n2
  if (ip == mxlen .AND. n2 < n1) outstr(ip:ip) = '|'
  if (iprint() > 120) print '(1x,a,a)', 'bin2a:',outstr(1:ip)
end subroutine bin2a

subroutine bin2al(fmt,res,count,strn)
  character*(*) fmt, strn
  integer :: count
  logical res(0:*)
  write(strn,fmt) res(count)
end subroutine bin2al

subroutine bin2ai(fmt,res,count,strn,lnull)
  !- Conversion of integer to ascii string
  ! ----------------------------------------------------------------------
  !i Inputs
  !l         :
  !i   fmt   : fortran format
  !i   res   : res(count) is converted
  !i   count : index to res: res(count) is converted
  !o Outputs
  !o   strn  : ascii representation of integer
  !i   lnull : true if res(count)=NULLI, otherwise false
  ! ----------------------------------------------------------------------
  !     implicit none
  character*(*) fmt, strn
  integer :: count
  integer :: res(0:count)
  logical :: lnull
  integer :: NULLI
  parameter (NULLI=-99999)
  write(strn,fmt) res(count)
  lnull = res(count) .eq. NULLI
end subroutine bin2ai

subroutine nlchar(ich,ps)
  integer :: ich
  character(*) ps
  if (ich==1) then
     ps(1:1) =  char(10)
  endif
end subroutine nlchar

integer function getdig(n,i,base)
  !- Extracts one digit from an integer
  ! ----------------------------------------------------------------
  !i Inputs
  !i   n,i,base
  !o Outputs
  !o   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
  ! ----------------------------------------------------------------
  implicit none
  integer :: n,i,base
  getdig = mod(n/base**i,base)
end function getdig

subroutine pretty(sin,nblk,ndec,precsn,sw,sout,nout)
  !- Prettifies an ascii representation of a floating-point number
  !i  nblk: number of leading blanks in output string
  !i  ndec: minimum number of decimals to display (see sw)
  !i  precsn: truncate to this many decimals of precision (see sw)
  !i  sw: 1's digit: 0 for ndec and precsn absolute (number of places
  !i                   to the right of '.').
  !i                 1 for ndec and precsn number of significant digits.
  !i     10's digit: 1 to suppress possible leading zero
  !o  sout: reformatted ascii representation of number
  !o  nout: length of sout
  !r  ndec and precsn must be >= 0
  !r  If ndec>precsn, ndec truncated to precsn (precsn>0)

  !     implicit none
  !      character*(*) sin,sout*40
  character*(*) sin,sout
  integer :: ndec,precsn,nout,sw,nblk
  integer :: ls,is,iex,ixp,i,j,k,m,i0,i1,ip,nd,im,id,lsin,getdig
  integer :: ndig,idec
  logical :: a2bin,ltmp,isdig
  character ss*40, sexp*5, ch*1, sss*40, fmt*8
  double precision :: xx,dround


  ! --- Copy input string to local ---
  !     print *, 'sin=',sin
  ss = ' '
  ss(2:39) = sin
  sout = ' '
  nout = 1+nblk
  i0 = nout
  lsin = len(sin)

  ! --- Setup, normalize variations of input string ---
  ! ... i1=pos of first char, ip = pos of decimal point
1 continue
  i1 = 0
  ls = len(ss)
  call skipbl(ss,ls,i1)
  call skpblb(ss,ls,is)
  i1 = i1+1
  ! ... Remove leading '+'
  if (ss(i1:i1) == '+') i1 = i1+1
  if (ss(i1:i1) == '-') then
     sout(nout:nout) = '-'
     i1 = i1+1
     nout = nout+1
     i0 = nout
  endif
  ! ... Check for exceptions
  if (ss(i1:i1) == '*') then
     sout(nout:nout) = '*'
     return
  endif
  ! ... Check for 'NaN or Infinity'
  if (ss(i1:ls) == 'nan' .OR. ss(i1:ls) == 'NaN' .OR. &
       ss(i1:ls) == 'NAN' .OR. ss(i1:ls) == 'QNAN') then
     sout(nout:nout+2) = 'NaN'
     nout = nout+2
     return
  endif
  if (ss(i1:ls) == 'inf' .OR. ss(i1:ls) == 'Inf' .OR. &
       ss(i1:ls) == 'INF' .OR. ss(i1:ls) == 'Infinity') then
     sout(nout:nout+2) = 'Inf'
     nout = nout+2
     return
  endif

  is = is+1
  ! ... Location of decimal point (ip); prepend 0 if missing
  ip = 0
  call chrpos(ss,'.',is,ip)
  ip = ip+1
  if (ip > is) then
     is = is+1
     ss(is:is) = '.'
  endif
  if (ip == i1) then
     i1 = i1-1
     ss(i1:i1) = '0'
  endif

  ! --- Find end of mantissa; get exponent (binary in ixp) ---
  iex = 0
  im = ip-1
  ixp = 0
  call chrps2(ss,'dDeE+-',6,is,im,j)
  if (j /= 0) then
     iex = im+2
     if (j >= 5) iex = iex-1
     if (ss(iex:iex) == '+' .OR. ss(iex:iex) == '-') iex = iex+1
     do  10  j = iex, is
        k = j
        if (ss(j:j) /= '0') goto 12
10   enddo
12   continue
     if (ss(iex-1:iex-1) == '-') then
        k = k-1
        ss(k:k) = '-'
     endif
     iex = k
     k = iex-1
     if ( .NOT. a2bin(ss,ixp,2,0,' ',k,-1)) &
          call rx('pretty: parse error')
     !        print *, iex,is, ':', ss(1:iex-1), '|', ss(iex:is), '|', ixp
  endif

  ! --- Truncate string at precsn decimal places, possibly rounding ---
  if (precsn > 0) then
     k = ip+precsn+ixp
     m = i1
     ! ...   If relative precision, k offset from first significant digit
     if (getdig(sw,0,10) == 1) then
        do  14  i = i1, im
           m = i
           if (ss(i:i) >= '1' .AND. ss(i:i) <= '9') goto 15
14      enddo
15      continue
        k = m-1+precsn
        if (m <= ip .AND. k >= ip) k = k+1
     endif
     k = max(min(k,im),m)
     !        print *, k, ' ##|', ss(i1:k), '|'
     ! ...   Round upwards, append exponent and start over
     if (ss(k+1:k+1) >= '5' .AND. ss(k+1:k+1) <= '9') then
        im = min(k,im)
        sss = ss(i1:max(im+1,ip))
        read(sss,'(e30.20)') xx
        fmt = '(f30.'
        if (im-ip >= 10) write(fmt(6:8),'(i2,'')'')') im-ip
        if (im-ip < 10) write(fmt(6:7),'(i1,'')'')') max(im-ip,0)
        ! ...     Next line makes some recalcitrant machines round properly
        if (im > ip) xx = xx + 4d0/10d0**(1+im-ip)
        if (im > ip) write(sss,fmt) xx
        if (im <= ip) write(sss,fmt) dround(xx,k+1-i1)
        !          print *, sss, im-ip
        ! ...     Append exponent
        if (j > 0) then
           im = 0
           call chrps2(sin,'0123456789',10,lsin,im,j)
           call chrps2(sin,'dDeE+-',6,lsin,im,j)
           call skpblb(sss,40,k)
           call strcop(sss(k+2:40),sin(im+1:lsin),lsin-im,' ',j)
           !            print *, sss, im-ip
        endif
        ! ...     Start over
        ss = sss
        goto 1
     endif
     if (k < ip) then
        do   j = k+1, ip-1
           ss(j:j) = '0'
        enddo
     endif
     !        print *, precsn, '(precision):', ss(1:im), '|', k,im
     im = max(min(k,im),ip)
  endif

  ! --- Patch 0.mmmEnnn so that leading digit is nonzero ---
  if (iex /= 0 .AND. ixp /= 0) then
     if (ss(i1:i1+1) == '0.') then
        ixp = ixp-1
        i1 = i1+1
        ss(i1:i1) = ss(i1+1:i1+1)
        ss(i1+1:i1+1) = '.'
        ip = i1+1
        !          print *, 'patch', ss(i1:is)
     endif
     sexp = ' '
     write(sexp,'(i5)') ixp
     !        print *, 'exponent', sexp
  endif

  ! --- Copy mantissa to sout, counting decimals ---
  idec = 0
  ndig = 0
  ltmp = getdig(sw,1,10) .ne. 0
  do  20  j = i1, im
     isdig = ss(j:j) .ge. '0' .and. ss(j:j) .le. '9'
     if (j == i1 .AND. ltmp .AND. ss(j:j) == '0') goto 20
     if (j == ip .OR. isdig) then
        sout(nout:nout) = ss(j:j)
        nout = nout+1
        if (isdig .AND. j > i1) ndig = ndig+1
        if (isdig .AND. j > ip) idec = idec+1
     else
        goto 22
     endif
20 enddo
22 continue
  nout = nout-1
  idec = idec-ixp
  !      print *, 'string before truncation:', sout(1:nout),
  !     .  '...', idec, ' decimals', ndig, ' digits; ndec=',ndec

  ! --- Append or strip trailing zeros, preserving ndec decimals ---
  ip = 0
  call chrpos(sout,'.',nout,ip)
  ip = ip+1
  if (ip > nout) then
     ss = sin
     i1 = 0
     ls = len(ss)
     call skipbl(ss,ls,i1)
     call skpblb(ss,ls,is)
     call rx('pretty:  missing ''.'' in string'//ss(i1:is+1))
  endif
  ! ... Find i1 which points to 1st significant digit
  i1 = 0
  call chrps2(sout,'123456789',9,nout,i1,j)
  i1 = i1+1
  ! ... Case number is zero
  if (i1 > nout) then
     i1 = 0
     call chrpos(sout,'0',nout,i1)
     i1 = i1+1
     if (i1 > nout) call rx('pretty:  missing digit')
  endif
  id = 0
  nd = ndec
  if (precsn > 0) nd = min(nd,precsn)
  if (ndec > 0) then
     id = ip+nd+ixp
     if (getdig(sw,0,10) == 1) id = i1+nd-1
     if (getdig(sw,0,10) == 1 .AND. i1 < ip) id = id+1
  endif
  ! ... Append
  if (nout < id) then
     do   k = nout+1, id
        sout(k:k) = '0'
     enddo
     nout = id
  endif
  ! ... Strip
  if (0<id .AND. id<=len(sout))then
     if (sout(id:id) == '.') id = id-1
  endif
  do  41  k = nout, i0, -1
     j = k
     !        print *, '!!!',sout(1:j),' ',j,id
     ch = sout(j:j)
     if (j <= id) goto 42
     if (ch == '0' .OR. ch == '.') sout(j:j) = ' '
     if (ch /= '0' .OR. ch == '.') goto 42
41 enddo
42 continue
  ! ... Fix up a completely missing mantissa
  nout = j
  if (nout == i0 .AND. sout(i0:i0) == ' ') then
     sout(i0:i0) = '0'
     return
  endif
  if (sout(nout:nout) == ' ') nout = nout-1

  ! --- Append exponent if there is one ---
  if (ixp /= 0) then
     j = 0
     call skipbl(sexp,5,j)
     !        print *, j,':',sexp(j+1:5)
     sout(nout+1:nout+7-j) = 'e' // sexp(j+1:5)
     nout = nout+6-j
  endif

  !      print *, 'final string:', sout(1:nout)
  !      print *, '--------------'

end subroutine pretty

double precision function dround(x,n)
  !- Rounds double precision x after n digits
  !     implicit none
  integer :: n,is,i
  double precision :: x,s
! #if QUAD
!   double precision :: xnint
! #endif
! #if QUAD
!   s = qlog10(dble(abs(x)))
! #else
  dround=1d99
  s = dlog10(dabs(x))
!#endif
  is = s
  if (is > 0) then
     is = -int(is) + n-1
  else
     is = int(-is) + n
  endif
  s = 1
  do    i = 1, iabs(is)
     s = 10*s
  enddo
!#if QUAD
!  if (is > 0) dround = xnint(s*dble(x))/s
!  if (is < 0) dround = xnint(dble(x)/s)*s
!#else
  if (is > 0) dround = dnint(s*x)/s
  if (is < 0) dround = dnint(x/s)*s
!#endif
END function dround
