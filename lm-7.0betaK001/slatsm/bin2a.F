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
C     implicit none
C Passed Parameters
      character*(*) outstr
      character*(*) fmt
      double precision res(0:1)
      integer nblk,cast,count,ip,mxlen,ndec,is
C Local Variables
      logical lD,lS,lF,lnull,llnull,parstr
      integer i,j,k,iprint,lsmx,n1,n2,np,precsn,fw,
     .  ix(4),iv(4),a2vec,p,isw,isw0,getdig,m,ndig,ndige
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
   10 call bin2al('(L8)',res,count,strn2)
      goto 15
   11 strn2 = fmt
      goto 15
   12 call bin2ai('(I16)',res,count,strn2,lnull)
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
     .        fw = max(p+3,3) + max(precsn-p-1,0)
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
     .          '*****************************************************'
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
     .    (n2 .lt. n1 .or. strn(nblk+1:nblk+1) .eq. '*')) then
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
      logical res(0:1)
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
      integer res(0:1)
      logical lnull
C ... Local parameters
      integer NULLI
      parameter (NULLI=-99999)
      write(strn,fmt) res(count)
      lnull = res(count) .eq. NULLI
      end

