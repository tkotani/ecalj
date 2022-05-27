      integer function awrite(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
C- Formatted output, with ascii conversion of binary numbers
Ci ifi: <>0, local output string written to abs(ifi);
Ci      <=0, sout copied to local output string initially;
Ci           local output string copied back to sout on exit
Ci       >0, sout unaltered on exit
Ci mxln: abs(mxln) = maximum number of characters to copy
Ci mxln < 0: suppress trailing blanks when writing to logical unit.
Co sout:     output string (see ifi, above)
Cr Characters are copied from fmt to the output string sout, which is
Cr   then (optionally) written to logical unit ifi.  Pointer ip keeps
Cr   track of the current position for writing to sout.  Copy
Cr   is literal except when a control char % is encountered.
Cr   % characters do one of several functions:
Cr   %l writes to sout an ascii representation of logical argument a_j
Cr      (NB: j is 1 for first conversion, 2 for second, etc).
Cr   %i writes an integer argument a_j
Cr   %d, %e, %g, %G, %D and %F write ascii rep of double precision a_j
Cr     'd' writes in decimal notation
Cr     'e' writes in exponential notation
Cr     'g' and 'G' take the minimum size of 'd' and 'e'; 'g' is to
Cr                 specify relative precision, 'G' absolute precision.
Cr     All of the above generate 'pretty' ascii representations
Cr     (see pretty.f)
Cr     'D' mimics the fortran 'f' format and is intended for output in
Cr     fixed columns.
Cr     'F' puts the number in a specified space, using whatever form
Cr     produces the most decimal places of precision.
Cr   %% quotes a "%" literally
Cr   %a shifts ip past last nonblank character
Cr   %f shifts ip forward
Cr   %b shifts ip backward
Cr   %p sets   ip to a fixed value
Cr   %t is obsolete
Cr   %x blanks the output string
Cr   %z can suppress leading the leading zero in a decimal fraction
Cr   %W shifts ip forward until a whitespace is encountered
Cr   %w shifts ip forward until a non-whitespace is encountered
Cr   %c closes up whitespace around ip
Cr   %o opens  up whitespace around ip
Cr   %u if numerical argument = NULLI output 'null' instead of number
Cr      Optional argument n1:
Cr      0 turn off null option, for this and future calls
Cr      1 (or default) set null option, for this and future calls
Cr     >1 set null option for this call only
Cr     <0 turn off null option for this call only
Cr   %? conditionally parses one of two strings (see below)
Cr   %j increments argument jumps over call arguments
Cr   %N is turned into a newline, calling nlchar to get newline
Cr Most control characters have optional arguments.
Cr For d,e,g,G,D,F,l,i the general syntax is:
Cr   %[n1][:n2][,n3][;n4][#n5]x, with x=d,e,g,G or F
Cr Here n1..n5 are integer expressions:
Cr   n1 number of values to convert (a_j is regarded as a vector)
Cr   n2 number of blank spaces preceding first character
Cr      n2<0 => subtract one space if argument is negative
Cr   n3 minimum number of digits to display (after '.' for 'd'
Cr      and 'G' and total number for 'e' and 'g')
Cr   n4 round to n4 decimal places (absolute for 'd' and 'G',
Cr      and relative for 'e' and 'g')
Cr   n5 if conversion uses less than n5 characters, append trailing
Cr      blanks to fill (used for lining data in columns)
Cr For D the meanings of n2..n4 are different:
Cr   n2 is not used
Cr   n3 number of digits after decimal
Cr   n4 field width
Cr For F:
Cr   n2 is not used
Cr   n3 is not used
Cr   n4 is the field width
Cr For l and i:
Cr   n3 is the field width
Cr For z, j, p, a, f, o, b, and and the general syntax is:
Cr   %[n1]x, with x=z, p, a, f, b, u
Cr   n1 repeats (f, b)
Cr   n1 1=>suppresses leading 0, 0=>ensures its presence (z)
Cr   For u, see above
Cr NB: there is an option to substitute for any of n1..n4 one of the
Cr arguments a_j.  This is done by using the character 'n' instead
Cr of some integer expression (eg %n:n,5d).  awrite uses the
Cr next argument a_j is used for n, and increments j.  Thus, %n:n,5d
Cr consumes the next three a_j, the first describing the number
Cr of elements to write, the second the number of spaces between
Cr arguments.
Cr For ? the general syntax is
Cr   %?QexprQstr1Qstr2Q
Cr   str1 is parsed if "expr" evaluates to nonzero, otherwise str2 is.
Cr   Q is some character, eg ';'.  It should NOT be some character
Cr   that may be confused as part of "expr", like '?' or '/'.
Cr   The next argument argument a_j is temporarily set to an integer
Cr   value and temporarily named `n', which may be used in 'expr'. 
Cr   Also the current number of characters in the string is temporarily
Cr   assigned to `p'.  Finally, as a special case for a expression
Cr   involving strings, the following is permitted:
Cr     %c==X
Cr   where X is some character.  This expression evaluates to nonzero
Cr   if the character output string at the current position is equal
Cr   to X; otherwise it evaluates to zero.
Cr   Example:
Cr     call awrit2('three plus one is %?;n==1;%i;four;, no?',mxlen,
Cr                 s,mxlen,-i1mach(2),m,4)
Cr     prints out "three plus one is 4, no?" if m equals 1; otherwise
Cr     prints out "three plus one is four, no?"
Cu Updates
Cu   13 Oct 07 Modified lnull, for permanent option
Cu   02 Aug 07 Added %u: outputs null string when arg=NULLI (bin2av)
Cu   27 Feb 02 Added %c==X type of conditional expression
Cu    8 May 01 addition of n5 modifier described above
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,mxln
      character*(*) fmt,sout
      double precision
     .a1(*),a2(*),a3(*),a4(*),a5(*),a6(*),a7(*),a8(*)
C ... Local parameters
      integer i,lfmt,ip,jp,cast,ia,iff,i2,iterm,ndec,iv(5),nblk,nx,j,
     .mxl,ls,icond,ivawrt,lens,nsyv,ires,fw
      equivalence (iv(1),i2),(iv(2),nblk),(iv(3),nx),(iv(4),ndec),
     .(iv(5),fw)
      logical a2bin,ltmp,lnull,lnulls
      double precision holdn,holdp,xx
      parameter (lens=1024)
      character*(lens) s,ss,fchr*29,fm*20,cc*1,ccond*1
      save lnulls
      data fchr /' :;,#irdegGltapfbzDwWcoxu?jNF'/
      data lnulls /.false./

ckino
      logical:: l_exec

C --- Setup ---
C ... ia,iff,ip: indices to current argument, pos in fmt, pos in s
      mxl = iabs(mxln)
      s = ' '
      if (ifi .le. 0) s = sout
C ... ip holds the current position in the output string
      ip = 0
C ... iff holds the current position in the format string
      iff = 0
C ... index to current argument in the argument list
      ia = 1
C ... icond nonzero when in the middle of a conditional expression
      icond = 0
C ... ccond is the terminating character for conditional expression
      ccond = ' '
C ... hold on to n,p in vars table; we need them as local variables
      call numsyv(nsyv)
      call getsyv('p',holdp,j)
      call getsyv('n',holdn,j)
      lfmt = len(fmt)
      ls = len(s)
      lnull = lnulls
#if DEBUG
      print '(''entering awrite, fmt='',a)', fmt
#endif

C --- Parse next character in fmt ---
   19 ia = ia-1
   20 iff = iff+1
C  ...  End of fmt
      if (iff .gt. lfmt) goto 10
#if DEBUG
      print *, 'now parsing #', iff, ' char=',fmt(iff:iff)
#endif
C  ...  Character terminating conditional string
      if (icond .gt. 0 .and. fmt(iff:iff) .eq. ccond) then
        if (icond .eq. 2) then
          call chrpos(fmt,ccond,lfmt,iff)
          iff = iff+1
        endif
        icond = 0
        goto 20
      endif
C  ...  Any non-% character
ckino      if (fmt(iff:iff) .ne. '%' .or. fmt(iff:iff+1) .eq. '%%') then
        l_exec=.false.
        if (fmt(iff:iff) .ne. '%') l_exec=.true.
        if (iff+1<=len(fmt)) then
           if (fmt(iff:iff+1) .eq. '%%')  l_exec=.true.
        endif
      if (l_exec) then
*         print *, 'parsing non-%:', iff, fmt(1:iff)
        ip = ip+1
        if (ip .le. min(ls,mxl)) s(ip:ip) = fmt(iff:iff)
        if (iff+1<=len(fmt)) then
        if (fmt(iff:iff+1) .eq. '%%') iff = iff+1
        endif
        goto 20
C   --- Parse % ---
      else
C     ... Default values for %command
*         print *, 'now parsing %:', iff, fmt(iff:min(lfmt,iff+10))
        nblk = 0
        fw = 0
        nx = 99
        ndec = 0
        i2 = 1
        ia = ia+1
        iff = iff+1
C     ... iterm flags whether cc is ':;,#', to use later
        j = 0
        call chrps2(fmt(iff:iff),fchr,len(fchr),0,j,iterm)
C     ... Re-entry if to parse another argument
   25   if (iterm .ge. 2 .and. iterm .le. 5) iff = iff+1
        cc = fmt(iff:iff)
C     ... ires is integer argument; set to 1 for default value
        ires = 1
C     ... Check for an integer expression () preceding command:
        if (cc .eq. '(') then
          j = 1
          do  35  i = iff+1, lfmt
            if (fmt(i:i) .eq. '(') j=j+1
            if (fmt(i:i) .eq. ')') j=j-1
            if (j .eq. 0) then
              xx = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
              call lodsyv('n',1,xx,j)
              call lodsyv('p',1,dble(ip),j)
              j = iff-1
              ltmp = .not. a2bin(fmt,ires,2,0,' ',j,i-1)
              call rxx(ltmp,'awrite: failed to parse () in format')
              iff = i+1
              goto 36
            endif
   35     continue
          call rx('awrite: missing matching () in format')
   36     continue
C     ... Prior integer argument if next char 'n':
        else if (cc .eq. 'n') then
          ires = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
          ia = ia+1
          iff = iff+1
C     ... Prior integer argument if char an integer:
        else if (cc .ge. '0' .and. cc .le. '9' .or. cc .eq. '-') then
          do  22  i = iff, lfmt-1
            j = i+1
            if (fmt(j:j) .ge. '0' .and. fmt(j:j) .le. '9') goto 22
            j = iff-1
*             call pshpr(130)
            ltmp = .not. a2bin(fmt,ires,2,0,fmt(i+1:i+1),j,i)
            call rxx(ltmp,'awrite: failed to parse format')
            iff = i+1
            goto 23
   22     continue
   23     continue
        endif
*          print 335, ires,iff,fmt(1:iff)
* 335     format('*now ires=',i2,' parsed to ',i3,': ',a)
C     ... If this was an argument to one of ':;,#'
        if (iterm .ge. 2 .and. iterm .le. 5) then
          iv(iterm) = ires
C     ... Otherwise ires is an argument to command
        else
          iv(1) = ires
        endif
C     ... Next character is the terminator
        cc = fmt(iff:iff)
        j = 0
        call chrps2(cc,fchr,len(fchr),0,j,iterm)
C     ... If an argument, run through parse again
        if (iterm .ge. 2 .and. iterm .le. 5) goto 25
C     ... Otherwise a command:
        cast = 99
        cc = fmt(iff:iff)
        if (cc .eq. 'l') cast=0
        if (cc .eq. 'i') cast=2
        if (cc .eq. 'r') cast=3
        if (cc .eq. 'd' .or. cc .eq. 'e' .or. cc .eq. 'D' .or.
     .  cc .eq. 'g' .or. cc .eq. 'G' .or. cc .eq. 'F')
     .  cast=4
        if (cc .eq. 't') then
          call rx('awrite: use p, not t')
        endif
        if (cc .eq. 'z') then
          call bin2a0(i2)
          goto 19
        elseif (cc .eq. 'j') then
          ia = ia+i2
          goto 19
        elseif (cc .eq. 'a') then
          call skpblb(s,ls,ip)
          ip = ip+i2
          goto 19
        elseif (cc .eq. 'p') then
          ip = i2
          goto 19
        elseif (cc .eq. 'f') then
          ip = ip+i2
          goto 19
        elseif (cc .eq. 'b') then
          ip = ip-i2
          goto 19
        elseif (cc .eq. 'N') then
          call nlchar(1,s(ip+1:ip+1))
          ip = ip+1
          goto 19
C ---     Entry point for conditional expression ---
        elseif (cc .eq. '?') then
          if (icond .ne. 0) call rx('awrite encountered nested "%?"')
          icond = 1
          call lodsyv('p',1,dble(ip),j)
          xx = ivawrt(ia,1,a1,a2,a3,a4,a5,a6,a7,a8)
          call lodsyv('n',1,xx,j)
          ia = ia+1
          iff = iff+1
C     ...   ccond is character terminating conditional string
          ccond = fmt(iff:iff)
C     ...   If next char is '%', expression is of the string type:
          if (fmt(iff+1:iff+1) .eq. '%') then
            if (fmt(iff+1:iff+4) .eq. '%c==') then              
              ltmp = fmt(iff+5:iff+5) .eq. s(ip:ip)
              iff = iff+6
            else
              call rxs('awrite: failed to parse : ',fmt(iff:))
            endif
C     ...   Parse expression
          elseif (.not. a2bin(fmt,ltmp,0,0,ccond,iff,lfmt)) then
            call rx('awrite: failed to parse conditional expr')
          endif
C     ...   Use first string, or skip to second string
          if (ltmp) then
            icond = 2
          else
            icond = 3
            call chrpos(fmt,ccond,lfmt,iff)
            iff = iff+1
          endif
          goto 19
C     ... clear string
        elseif (cc .eq. 'x') then
          s = ' '
          goto 19
C     ... toggle on 'null' option
        elseif (cc .eq. 'u') then
          if (i2 .lt. 0) then
            lnull = .false.
          elseif (i2 .eq. 0) then
            lnulls = .false.
            lnull = .false.
          elseif (i2 .eq. 1) then
            lnulls = .true.
            lnull = .true.
          else
            lnull = .true.
          endif
          goto 19
C ...     pad whitespace around ip
        elseif (cc .eq. 'o') then
          if (ip .eq. 0) ip = 1
          ss = s(ip:ls)
          s(ip:ip+i2-1) = ' '
          s(ip+i2:ls) = ss
          ip = ip+i2-1
          goto 19
C ...     close up whitespace around ip
        elseif (cc .eq. 'c') then
          do  13  j = ip, 1, -1
            if (s(j:j) .ne. ' ' .and. s(j:j) .ne. '        ') goto 14
            ip = j
   13     continue
   14     continue
          jp = ip-1
          do  15  j = jp+1, mxl
            if (s(j:j) .ne. ' ' .and. s(j:j) .ne. '        ') goto 16
            jp = j
   15     continue
   16     continue
          if (jp-ip+1 .gt. 0) then
            ss = s(jp+1:ls)
            s(ip:ls) = ss
          endif
          goto 19
C ...     skip to next nw
        elseif (cc .eq. 'w') then
          do  17  j = ip, mxl
            if (s(j:j) .ne. ' ' .and. s(j:j) .ne. '        ') goto 19
            ip = j
   17     continue
C ...     skip to next whitespace
        elseif (cc .eq. 'W') then
          do  18  j = ip, mxl
            if (s(j:j) .eq. ' ' .or. s(j:j) .eq. '        ') goto 19
            ip = j
   18     continue
        endif
      endif
      if (cast .eq. 99) call rx('awrite: unknown control: ' // cc)

C ---   Generate format for bin2a ---
      fm = ' '
      if (cast .eq. 4) then
        fm = cc
        if (cc .eq. 'G') fm = 'g'
        j = 1
        if (nx .ne. 99) call bin2a(' ',0,0,nx,2,0,20,fm,j)
        if (cc .eq. 'G') call bin2a(':20',0,0,0,1,0,20,fm,j)
      endif

C ---   Convert binary numbers ---
      i2 = i2-1
C        fw = 0
C        if (nblk .lt. 0) then
C          fw = -nblk
C          nblk = -1
C        endif
      if (ia.eq.1)
     .call bin2av(fm,fw,nblk,ndec,a1,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.2)
     .call bin2av(fm,fw,nblk,ndec,a2,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.3)
     .call bin2av(fm,fw,nblk,ndec,a3,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.4)
     .call bin2av(fm,fw,nblk,ndec,a4,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.5)
     .call bin2av(fm,fw,nblk,ndec,a5,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.6)
     .call bin2av(fm,fw,nblk,ndec,a6,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.7)
     .call bin2av(fm,fw,nblk,ndec,a7,cast,0,i2,' ',mxl,lnull,s,ip)
      if (ia.eq.8)
     .call bin2av(fm,fw,nblk,ndec,a8,cast,0,i2,' ',mxl,lnull,s,ip)

      goto 20
   10 continue


C --- Finish up and exit ---
      ip = min(ip,mxl)
      if (mxln .lt. 0) then
        call skpblb(s,ip,ip)
        ip = ip+1
      endif
      ia = iabs(ifi)
      if (ifi .ne. 0 .and. ip .gt. 0) write(ia,333) s(1:ip)
  333 format(a)
      if (ifi .le. 0 .and. ip .gt. 0) sout = s(1:ip)
      awrite = ip

C --- Restore or undo symbolic variables p,n ---
*     call shosyv(0,0,0,6)
      call lodsyv('p',1,holdp,j)
      call lodsyv('n',1,holdn,j)
      call clrsyv(nsyv)
*     call shosyv(0,0,0,6)
      return

      end
      subroutine bin2av(fmt,w,nblk,ndec,res,cast,i1,i2,sep,mxln,lnull,
     .outs,ip)
C- Write out a vector of of numbers using bin2a
C ----------------------------------------------------------------------
Ci Inputs
Ci   fmt   :format passed to bin2a
Ci    w    :unused if zero.  If >0,
Ci         :w = minimum spacing between successive numbers
Ci   nblk  :number of blanks preceding each value, or if nblk < 0,
Ci         :|nblk| spaces are prepended for positive numbers
Ci         :|nblk|-1 spaces are prepended for negative numbers
Ci   ndec  :retain a mininimum ndec digits after decimal (see bin2a)
Ci   res   :vector of binaries to convert to ascii string
Ci   cast  :0=logical, 1=char, 2=int, 3=real, 4=double
Ci   i1    :convert numbers res(i1..i2)
Ci   i2    :convert numbers res(i1..i2)
Ci   sep   :separator between numbers
Ci   mxln  :maximum allowed value of ip
Ci   ip    :string position pointer
Ci   lnull :if T, numbers equal to NULLI are turned into NULL
Co Outputs
Co   outs  :string containing ascii rep'sn of binary numbers
Ci   ip    :string position pointer updated to end of string
Cr Remarks
Cu Updates
Cu   01 Aug 07 new lnull
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) fmt,outs,sep*1
      integer nblk,cast,i1,i2,ip,mxln,ndec,w
      double precision res(0:1)
      logical lnull
C ... Local parameters
      logical lneg,llnull
      integer nblk2,ival,i,k,ip0,NULLI
      real rval
      double precision dval
      parameter (NULLI=-99999)

      if (mxln .le. 0) return
      do  i = i1, i2
        nblk2 = nblk
        if (nblk .lt. 0) then
          nblk2 = -nblk
          lneg = .false.
          if (cast .eq. 2) lneg = ival(res,i+1) .lt. 0
          if (cast .eq. 3) lneg = rval(res,i+1) .lt. 0
          if (cast .eq. 4) lneg = dval(res,i+1) .lt. 0
          if (lneg) nblk2 = nblk2-1
        endif
C       Set flag llnul if lnull is ON and argument matches NULLI
        llnull = .false.
        if (lnull) then
          if (cast .eq. 2) llnull = ival(res,i+1) .eq. NULLI
          if (cast .eq. 3) llnull = rval(res,i+1) .eq. NULLI
          if (cast .eq. 4) llnull = dval(res,i+1) .eq. dble(NULLI)
          if (llnull) then
            call skpblb(fmt,len(fmt),ip0)
            fmt(2+ip0:) = ':n'
          endif
        endif
        ip0 = ip
        call bin2a(fmt,nblk2,ndec,res,cast,i,mxln,outs,ip)
        if (llnull) then
C         If fixed width, leave position of null as is
          if (fmt(1:1) .eq. 'D' .or. fmt(1:1) .eq. 'F' .or.
     .    (cast .eq. 2 .and. ndec .gt. 0)) then
C         Skip if not sufficient space for leading blanks + null
          else if (ip-3 .le. 1+ip0+iabs(nblk)) then
C         Otherwise rewrite null starting at 1+ip0+iabs(nblk)
          else
            outs(1+ip0+iabs(nblk):ip) = 'NULL'
            ip = 4+ip0+iabs(nblk)
          endif
C         print *, outs(1:ip) 
        endif
        if (sep .ne. ' ' .and. i .lt. i2) then
          ip = ip+1
          outs(ip:ip) = sep
        endif
        if (w .ne. 0) then
          do  k = ip+1, ip0+w
            outs(k:k) = sep
            ip = ip+1
          enddo
        endif
      enddo
      end

      subroutine awrit8(fmt,sout,mxln,ifi,a1,a2,a3,a4,a5,a6,a7,a8)
C- Subroutine versions of integer function awrite
C     implicit none
      double precision a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
      character*(*) sout,fmt
      integer ifi,mxln,ip,jp,awrite
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

      end
      subroutine vwrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8,cast,ires,res)
C- Writes either integer or double into ires or res, depending on cast
C ----------------------------------------------------------------------
Ci Inputs
Ci   ia    :indicates which of arrays a1..a8 to extract element from
Ci   n     :which entry in array a_ia
Ci   a1..a8:element is extracted from one of these arrays
Ci   cast  :array cast
Co Outputs
Co   ires  :if cast is integer, result poked into ires
Co   res   :if cast is double, result poked into res
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer ia,n,cast,ivawrt,ires
      double precision dvawrt,res
      double precision a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)

      if (cast .eq. 2) then
        ires = ivawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
      elseif (cast .eq. 4) then
        res = dvawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
      else
        call rxi('vwrt: cannot handle cast',cast)
      endif
      end

      integer function ivawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
C     implicit none
      integer ia,n,ival
      double precision a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)

      if (ia .eq. 1) ivawrt = ival(a1,n)
      if (ia .eq. 2) ivawrt = ival(a2,n)
      if (ia .eq. 3) ivawrt = ival(a3,n)
      if (ia .eq. 4) ivawrt = ival(a4,n)
      if (ia .eq. 5) ivawrt = ival(a5,n)
      if (ia .eq. 6) ivawrt = ival(a6,n)
      if (ia .eq. 7) ivawrt = ival(a7,n)
      if (ia .eq. 8) ivawrt = ival(a8,n)
      end
      double precision function dvawrt(ia,n,a1,a2,a3,a4,a5,a6,a7,a8)
C     implicit none
      integer ia,n
      double precision a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)

      if (ia .eq. 1) dvawrt = a1(n)
      if (ia .eq. 2) dvawrt = a2(n)
      if (ia .eq. 3) dvawrt = a3(n)
      if (ia .eq. 4) dvawrt = a4(n)
      if (ia .eq. 5) dvawrt = a5(n)
      if (ia .eq. 6) dvawrt = a6(n)
      if (ia .eq. 7) dvawrt = a7(n)
      if (ia .eq. 8) dvawrt = a8(n)
      end
C      subroutine awsnot(icond)
C      integer icond
Cc     print *, 'hi'
C      end

      subroutine bin2a(fmt,nblk,ndec,res,cast,count,mxlen,outstr,ip)
C- Converts number to ascii format, stripping leading blanks, trailing 0
C ----------------------------------------------------------------------
Ci Inputs:
Ci   fmt: cast=1:   holds the string to be appended to outstr)
Ci        cast=0,2: not used
Ci        cast=3,4: syntax X[#][:sw] where X is one of
Ci                  'd' to write in decimal representation
Ci                  'e' to write in exponential format
Ci                  'g' to use the smaller of 'd' and 'e'
Ci                  '(n.m)' fixed format, mimics fortran fmt (fn.m)
Ci                  'D' also mimics fortran fmt (Fn.m)
Ci                      D# => supplies n=#; arg ndec supplies m
Ci                  'F' fixed field, picking between d and e that
Ci                      F# => # is field width
Ci                      generates the most significant digits
Ci                  See Remarks for further description
Ci   nblk:  strip leading blanks, leaving a maximum of nblk
Ci   ndec:  (cast=3,4 only) retain a mininimum ndec digits after decimal
Ci          point, i.e. do not suppress trailing zeros to ndec.
Ci          ndec=0 does nothing.  ndec should not exceed precsn.
Ci          (cast=2 only): ndec specifies a field width
Ci   res:   binary value to be converted into ascii string
Ci   cast:  cast of res: 0=logical, 1=char, 2=int, 3=real, 4=double
Ci   count: res(count) is to be converted.  NB: count=0 for first entry 
Ci   mxlen: maximum length of outstr
Cio Inputs/Outputs
Cio  ip:    on input, starting position in outstr for write
Cio         NB: ip=0 points to first character in string
Cio  ip:    on output, position of final character written to outstr
Co  Outputs
Co   outstr:binary res(count) written in ascii form to outstr(ip:..)
Cr Remarks
Cr  *The string representation of floating point numbers is generated
Cr   by a "prettified" modification of the fortran write statement
Cr   (pretty.f), which includes suppression of trailing zeros and the
Cr   option to include or suppress the leading zero in decimal
Cr   fractions less than 1.  Floating-point formats include:
Cr     'd[n][:sw]' for decimal representation,
Cr     'e[n][:sw]' for exponential representation,
Cr     'g[n][:sw]' uses the minimum length of 'd' and 'e'
Cr     'D[n][:sw]' simulates the standard fortran format fn.m
Cr                 Here n follows D, ndec the role of m.  Or:
Cr     'Fn'        fixed field, picking between d and e that generates
Cr                 the most significant digits
Cr      (n.m)      also simulates the standard fortran format.
Cr
Cr  *Optional modifier 'n' is a number specifying how many decimals of
Cr   precision (n=6 if not specified). By default, n means:
Cr      for 'd' format, the absolute precision: i.e. 
Cr        number of digits after the decimal point
Cr     for 'e' format, the relative precision , i.e.
Cr        number of digits printed
Cr     for 'D' format, it is the field width n in fortran format fn.m 
Cr  *Optional modifier sw is a compound of the 1's and 10's digits.
Cr       1's digit of sw can overwrite the default meaning of 'n' above.
Cr                 sw=0 => n corresponds to absolute precision
Cr                 sw=1 => n corresponds to relative precision
Cr       10's digit nonzero suppresses leading blanks.
Cr  *Entry bin2a0 allows the user to set the default of sw.
Cr  *Examples:
Cr     call bin2a('d2',1,3,1.234951d0,...)    => 1.23
Cr     call bin2a('d4',1,4,1.234951d0,...)    => 1.2350
Cr     call bin2a('d3:11',1,0,1.2349501d-6,4) => .00000123
Cr     call bin2a('e2',1,3,1.2499d7,...)      => 1.2e7
Cr     call bin2a('e5',1,5,1.2349510d7,...)   => 1.2350e7
Cr     call bin2a('e5:0',1,4,1.2349501d5,...) => 1.234950100e5
Cr     call bin2a('g',1,0,1.23d-5,...)        => 1.23e-5
Cr     call bin2a('g3:10',1,3,1.24996d-5,...) => .000
Cr     call bin2a('g4:10',1,4,1.24996d-5,...) => 1e-5
Cr     call bin2a('f4:10',1,4,1.24996d-5,...) => 1e-5
Cu Updates
Cu   02 Aug 07 Added :n outputs null string when res=NULLI
C ----------------------------------------------------------------------
      implicit none
C Passed Parameters

cc kino's correctio for ifort was
cc     character(mxlen):: outstr ! ?---> !character(*) can not check size of outstr.  
cc However, because of a bug in grortran4.3.4, this is not allowed. Thus I now use character(*).
cc Sep2010
      character(*):: outstr 
cc
      character(*):: fmt
      double precision res(0:*)
      integer nblk,cast,count,ip,mxlen,ndec,is
C Local Variables
      logical lD,lS,lF,lnull,llnull,parstr
      integer i,j,k,iprint,lsmx,n1,n2,np,precsn,fw,
     .ix(4),iv(4),a2vec,p,isw,isw0,getdig,m,ndig,ndige
      parameter (lsmx=80)
      character*20 lfmt*20, fm*20, strn*(lsmx), strn2*(lsmx), ss*(lsmx)
      real rval
      double precision xx
      integer NULLI
      parameter (NULLI=-99999)
      save isw0
      data isw0 /0/

C     write(*,"('enter bin2a: cast,fmt=',i4,1x,a$)") cast,fmt

C --- Convert binary to ascii representation (log, int, char) ---
      lnull = .false.
      llnull = .false.
      goto (10,11,12,20,20), cast+1
      call rx('bin2a: bad cast')
 10   continue
      call bin2al('(L8)',res,count,strn2)
      goto 15
 11   continue
      strn2 = fmt
      goto 15
 12   continue
      call bin2ai('(I16)',res,count,strn2,lnull)
      goto 15
C --- copy strn2 to strn with appropriate number of blanks ---
   15 continue
      i = 0
      call skipbl(strn2,lsmx,i)
      strn = ' '
C     If a field width specified, overwrite spillover with '*'
      if (ndec .ne. 0) then
        call skpblb(strn2,lsmx,j)
        j = j-ndec+1
        if (j .gt. i) then
          strn2(j+1:j+ndec) = '****************'
        endif
        i  = j
      endif
      strn(1+nblk:lsmx) = strn2(i+1:lsmx)
      call skpblb(strn,lsmx,n1)
      n1 = n1+1
      if (lnull .and. fmt .ne. ' ') then
        i = 0
        if (parstr(fmt,':n',len(fmt)-1,2,'n',i,j)) then
          llnull = .true.
        endif
      endif
      goto 50

C --- Entry for setting up or determinining defaults ---
      entry bin2a0(is)
      if (is .ge. 0) isw0 = is
      if (is .lt. 0) is = isw0
      return

C --- Binary->ascii representation, floating-point ---
   20 continue
      if (cast .eq. 3) xx = rval(res,count+1)
      if (cast .eq. 4) xx = res(count)
      lnull = xx .eq. dble(NULLI)

C ... Determine appropriate format
      lfmt = fmt
      i = 0
      call skipbl(fmt,len(fmt),i)
      if (i .ge. len(fmt)) then
        lfmt = 'g'
      else
        lfmt = fmt(i+1:len(fmt))
      endif
      i = 0
      if (parstr(lfmt,':n',len(lfmt)-1,2,'n',i,j)) then
        lfmt(i+1:) = ' '
        llnull = .true.
      endif
C --- Do the conversion, floating point ---
      if (lfmt(1:1) .eq. '(') then
        write(ss,lfmt) xx
        call pretty(ss,nblk,ndec,20,isw0,strn,n1)
      else
        strn  = ' '
        strn2 = ' '
        lD = .false.
        lF = .false.
        j = 0
C   ... i=1 => 'd'  i=2 =>  'e'  i=3 => 'g'
        call chrps2(lfmt,'degDF',5,len(lfmt),j,i)
        if (i .le. 0) call rx('bin2a: bad format: '//lfmt)
        if (i .eq. 5) then
          i = 3
          lF = .true.
        elseif (i .eq. 4) then
          i = 1
          lD = .true.
        endif
C   ... Get precsn (or field width for D or F), in iv(1), sw in iv(2)
        j = j+1
        np = a2vec(lfmt,len(lfmt),j,2,': ',2,2,2,ix,iv)
        isw = 1 + isw0
        if (i .eq. 1) isw = 0 + isw0
C   ... Simulated fortran format: precsn dictated by ndec
        if (lF) then
          if (np .le. 0) call rx('bin2a: bad format: '//lfmt)
          fw = iv(1)
        elseif (lD) then
          precsn = ndec
          fw = -1
          if (np .ge. 1) fw = iv(1)
C   ... if precsn explicit, use it
        elseif (np .ge. 1) then
          precsn = iv(1)
C   ... This is the default
        else
          precsn = 6
        endif
        if (np .ge. 2) isw = iv(2)
        if (isw .ge. 20) isw = mod(isw,10) + isw0
C  21   continue
C   ... p is the exponent
        p = 0
        if (xx .ne. 0) then
          p = int(dlog10(dabs(xx)))
          if (dabs(xx) .lt. 1) p = p-1
        endif
C   ... fortran 'f' format
        if (i .eq. 1 .or. i .eq. 3) then
C     ... Estimate total width of format statement for fortran write
          if (lF) then
C       ... m is the space consumed by a '-' sign
            m = (1-int(dsign(1d0,xx)))/2
C       ... precsn = # rhs dec = field width - '-' - '.' - (p+1)
            precsn = fw - m - 1 - max(p+1,1)
C       ... Only works on some compilers
C            if (mod(isw,10) .ne. 0)
C     .      precsn = fw - m - 1 - max(p+1,0)
C       ... ndig = how many nonzero decimals printed
            ndig = max(precsn+p+1,0)
C       ... Exclude 'f' if it doesn't fit
            if (precsn .lt. 0) then
              ndig = -1
C       ... Exclude 'e' if it does, and number isn't small
            else if (p .gt. -2) then
              i = 1
            endif
C       ... Determine how many digits we get from 'e' format
            if (i .ne. 1) then
              write(ss,'(1pe20.0)') xx
C             print *, ss
C         ... We want at least 1 more digit than f format
              call pretty(ss,0,max(ndig+1,1),max(ndig+1,1),1,strn,j)
C         ... Tack on trailing 'e0' if pretty discarded it
              k = 0
              call chrpos(strn,'e',j,k)
              if (k .ge. j) j = j+2
C         ... How many decimals for 'e' format
              ndige = max(ndig+1,1) + fw - j
C         ... If pretty suppresses '.', add it back if ndige>1
              if (ndige .gt. 1) then
                k = 0
                call chrpos(strn,'.',j,k)
                if (k .ge. j) ndige=ndige-1
              endif
C             print *, strn
            else
              ndige = ndig-1
            endif
C       ... Generate string for F format here.
            if (ndig .lt. 0 .and. ndige .lt. 0) then
              strn = ' '
              strn(nblk+1:nblk+fw) = '********************************'
              n1 = fw+nblk
              goto 50
            else if (ndig .ge. ndige) then
              i = 1
            else
              i = 2
              precsn = ndige
              goto 35
            endif
          elseif (.not. lD .or. (lD .and. fw .eq. -1)) then
            fw = max(p+3,5) + precsn
            if (getdig(isw,0,10) .eq. 1)
     .      fw = max(p+3,3) + max(precsn-p-1,0)
            fw = max(fw,10)
            if (fw .gt. min(lsmx-2,99)) then
              strn = ' '
              strn(nblk+1:nblk+1) = '*'
              n1 = 1+nblk
              goto 35
            endif
          endif
          j = fw
C     ... Insert leading blanks
C         if (lF) then
          if (lF .or. lD) then
            j = j+nblk
          endif
          if (j .ge. 10) write(fm,'(''(f'',i2,''.'')') j
          if (j .lt. 10) write(fm,'(''( f'',i1,''.'')') j
          k = j
C     ... Number of decimals for fortran write
          j = precsn
          if (.not. (lD .or. lF)) then
            if (getdig(isw,0,10) .eq. 1) j = precsn-p-1
            j = max(j,0)
          endif
C         decimals can't exceed field width - 1
          j = max(min(k-1,j),0)
          if (j .ge. 10) write(fm(6:8),'(i2,'')'')') j
          if (j .lt. 10) write(fm(6:7),'(i1,'')'')') j
          write(ss,fm) xx
          if (lD .or. lF) then
            if (nblk .le. 0) then
            elseif (ss(1:nblk) .ne. ' ') then
              ss(1:k) =
     .        '*****************************************************'
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
C     ...   Case fraction should have a leading '0'
            if (j .eq. 1 .and. getdig(isw,1,10) .eq. 0) then
              if (lS .and. k .gt. 1) strn(k-1:k) = '-0'
              if (.not. lS .and.k .gt. 0) strn(k:k) = '0'
C     ...   Case fraction should have no leading '0'
            elseif (j .eq. 2 .and. getdig(isw,1,10) .ne. 0) then
              if (lS) strn(k:k+1) = ' -'
              if (.not. lS)strn(k+1:k+1) = ' '
            endif
          else
*           print *, 'before pretty ...', fm, ss
            call pretty(ss,nblk,ndec,precsn,isw,strn,n1)
          endif
   35     continue
        endif
C    .. fortran 'e' format
        if (i .eq. 2 .or. i .eq. 3) then
          j = p + precsn
          if (getdig(isw,0,10) .eq. 1) j = precsn-1
          if (j .gt. 22) then
            strn2 = ' '
            strn2(nblk+1:nblk+1) = '*'
            n2 = 1+nblk
            goto 45
          endif
          j = min(max(j,0),99)
          if (j .ge. 10) write(fm,'(''(1pe30.'',i2,'')'')') j
          if (j .lt. 10) write(fm,'(''(1pe30.'',i1,'')'')') j
          write(ss,fm) xx
*         print *, 'before pretty ...', fm, ss
          j = ndec
          if (lF) j = precsn
          call pretty(ss,nblk,j,precsn,isw,strn2,n2)
C     ... Tack on trailing 'e0' if pretty discarded it
          j = 0
          call chrpos(strn2,'e',n2,j)
          if (j .ge. n2 .and. i .eq. 2) then
            strn2(n2+1:n2+2) = 'e0'
            n2 = n2+2
          endif
C     ... Sometimes the '.' is suppressed; make fw right
          if (lF .and. n2 .lt. fw+nblk) n2 = fw+nblk
   45     continue
        endif
*        if (i .eq. 3)
*     .    print *, n1,n2,'compare |', strn(1:n1), '|', strn2(1:n2), '|'
        if (i .eq. 2 .or. i .eq. 3 .and.
     .  (n2 .lt. n1 .or. strn(nblk+1:nblk+1) .eq. '*')) then
          strn = strn2
          n1 = n2
        endif
      endif

C --- Copy to outstr ---
   50 continue
      n1 = max(n1,0)
      n2 = max(min(n1,mxlen-ip),0)
C     Handle null number: replace ascii string with 'NULL'
      if (lnull .and. llnull) then
        strn(1:n2) = ' '
        i = max(n2-3,1+nblk)
        strn(i:n2) = 'NULL'
      endif
      if (n2 .gt. 0) outstr(ip+1:ip+n2) = strn(1:n2)
      ip = ip+n2
      if (ip .eq. mxlen .and. n2 .lt. n1) outstr(ip:ip) = '|'
      if (iprint() .gt. 120) print '(1x,a,a)', 'bin2a:',outstr(1:ip)

      end
      subroutine bin2al(fmt,res,count,strn)
      character*(*) fmt, strn
      integer count
      logical res(0:*)
      write(strn,fmt) res(count)
      end
      subroutine bin2ai(fmt,res,count,strn,lnull)
C- Conversion of integer to ascii string
C ----------------------------------------------------------------------
Ci Inputs
Cl         :
Ci   fmt   : fortran format
Ci   res   : res(count) is converted
Ci   count : index to res: res(count) is converted
Co Outputs
Co   strn  : ascii representation of integer
Ci   lnull : true if res(count)=NULLI, otherwise false
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) fmt, strn
      integer count
      integer res(0:count)
      logical lnull
C ... Local parameters
      integer NULLI
      parameter (NULLI=-99999)
      write(strn,fmt) res(count)
      lnull = res(count) .eq. NULLI
      end

