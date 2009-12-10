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
C ... Patch to handle strings of length longer than 267 on the CRAY.
#if CRAY
      integer nstrn,ind,nmxstr
      parameter (nmxstr=267)
#endif
      character*(*) fmt,sout
      double precision
     .  a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
C ... Local parameters
      integer i,lfmt,ip,jp,cast,ia,iff,i2,iterm,ndec,iv(5),nblk,nx,j,
     .  mxl,ls,icond,ivawrt,lens,nsyv,ires,fw
      equivalence (iv(1),i2),(iv(2),nblk),(iv(3),nx),(iv(4),ndec),
     .  (iv(5),fw)
      logical a2bin,ltmp,lnull,lnulls
      double precision holdn,holdp,xx
      parameter (lens=1024)
      character*(lens) s,ss,fchr*29,fm*20,cc*1,ccond*1
      save lnulls
      data fchr /' :;,#irdegGltapfbzDwWcoxu?jNF'/
      data lnulls /.false./

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
        if (fmt(iff:iff) .ne. '%' .or. fmt(iff:iff+1) .eq. '%%') then
*         print *, 'parsing non-%:', iff, fmt(1:iff)
          ip = ip+1
          if (ip .le. min(ls,mxl)) s(ip:ip) = fmt(iff:iff)
          if (fmt(iff:iff+1) .eq. '%%') iff = iff+1
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
   25     if (iterm .ge. 2 .and. iterm .le. 5) iff = iff+1
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
   35       continue
            call rx('awrite: missing matching () in format')
   36       continue
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
   22       continue
   23       continue
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
     .        cc .eq. 'g' .or. cc .eq. 'G' .or. cc .eq. 'F')
     .      cast=4
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
              if (s(j:j) .ne. ' ' .and. s(j:j) .ne. '	') goto 14
              ip = j
   13       continue
   14       continue
            jp = ip-1
            do  15  j = jp+1, mxl
              if (s(j:j) .ne. ' ' .and. s(j:j) .ne. '	') goto 16
              jp = j
   15       continue
   16       continue
            if (jp-ip+1 .gt. 0) then
              ss = s(jp+1:ls)
              s(ip:ls) = ss
            endif
            goto 19
C ...     skip to next nw
          elseif (cc .eq. 'w') then
            do  17  j = ip, mxl
              if (s(j:j) .ne. ' ' .and. s(j:j) .ne. '	') goto 19
              ip = j
   17       continue
C ...     skip to next whitespace
          elseif (cc .eq. 'W') then
            do  18  j = ip, mxl
              if (s(j:j) .eq. ' ' .or. s(j:j) .eq. '	') goto 19
              ip = j
   18       continue
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
     .    call bin2av(fm,fw,nblk,ndec,a1,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.2)
     .    call bin2av(fm,fw,nblk,ndec,a2,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.3)
     .    call bin2av(fm,fw,nblk,ndec,a3,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.4)
     .    call bin2av(fm,fw,nblk,ndec,a4,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.5)
     .    call bin2av(fm,fw,nblk,ndec,a5,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.6)
     .    call bin2av(fm,fw,nblk,ndec,a6,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.7)
     .    call bin2av(fm,fw,nblk,ndec,a7,cast,0,i2,' ',mxl,lnull,s,ip)
        if (ia.eq.8)
     .    call bin2av(fm,fw,nblk,ndec,a8,cast,0,i2,' ',mxl,lnull,s,ip)

        goto 20
   10 continue


C --- Finish up and exit ---
      ip = min(ip,mxl)
      if (mxln .lt. 0) then
        call skpblb(s,ip,ip)
        ip = ip+1
      endif
      ia = iabs(ifi)
#if CRAY
      if (ifi .ne. 0) then
        nstrn = ip / nmxstr
        do  30  i = 1, nstrn
          ind = (i - 1)*nmxstr + 1
          write(ia,333) s(ind:(ind+nmxstr-1))
   30   continue
        if (nstrn*nmxstr .lt. ip) write(ia,333) s((nstrn*nmxstr+1):ip)
      endif
#else
      if (ifi .ne. 0 .and. ip .gt. 0) write(ia,333) s(1:ip)
#endif
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
     .  outs,ip)
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
     .       (cast .eq. 2 .and. ndec .gt. 0)) then
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

