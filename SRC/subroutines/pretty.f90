      subroutine pretty(sin,nblk,ndec,precsn,sw,sout,nout)
C- Prettifies an ascii representation of a floating-point number
Ci  nblk: number of leading blanks in output string
Ci  ndec: minimum number of decimals to display (see sw)
Ci  precsn: truncate to this many decimals of precision (see sw)
Ci  sw: 1's digit: 0 for ndec and precsn absolute (number of places
Ci                   to the right of '.').
Ci                 1 for ndec and precsn number of significant digits.
Ci     10's digit: 1 to suppress possible leading zero
Co  sout: reformatted ascii representation of number
Co  nout: length of sout
Cr  ndec and precsn must be >= 0
Cr  If ndec>precsn, ndec truncated to precsn (precsn>0)

C     implicit none
C      character*(*) sin,sout*40
      character*(*) sin,sout
      integer ndec,precsn,nout,sw,nblk
      integer ls,is,iex,ixp,i,j,k,m,i0,i1,ip,nd,im,id,lsin,getdig
      integer ndig,idec
      logical a2bin,ltmp,isdig
      character ss*40, sexp*5, ch*1, sss*40, fmt*8
      double precision xx,dround


C --- Copy input string to local ---
*     print *, 'sin=',sin
      ss = ' '
      ss(2:39) = sin
      sout = ' '
      nout = 1+nblk
      i0 = nout
      lsin = len(sin)

C --- Setup, normalize variations of input string ---
C ... i1=pos of first char, ip = pos of decimal point
    1 continue
      i1 = 0
      ls = len(ss)
      call skipbl(ss,ls,i1)
      call skpblb(ss,ls,is)
      i1 = i1+1
C ... Remove leading '+'
      if (ss(i1:i1) .eq. '+') i1 = i1+1
      if (ss(i1:i1) .eq. '-') then
        sout(nout:nout) = '-'
        i1 = i1+1
        nout = nout+1
        i0 = nout
      endif
C ... Check for exceptions
      if (ss(i1:i1) .eq. '*') then
        sout(nout:nout) = '*'
        return
      endif
C ... Check for 'NaN or Infinity'
      if (ss(i1:ls) .eq. 'nan' .or. ss(i1:ls) .eq. 'NaN' .or.
     .ss(i1:ls) .eq. 'NAN' .or. ss(i1:ls) .eq. 'QNAN') then
        sout(nout:nout+2) = 'NaN'
        nout = nout+2
        return
      endif
      if (ss(i1:ls) .eq. 'inf' .or. ss(i1:ls) .eq. 'Inf' .or.
     .ss(i1:ls) .eq. 'INF' .or. ss(i1:ls) .eq. 'Infinity') then
        sout(nout:nout+2) = 'Inf'
        nout = nout+2
        return
      endif

      is = is+1
C ... Location of decimal point (ip); prepend 0 if missing
      ip = 0
      call chrpos(ss,'.',is,ip)
      ip = ip+1
      if (ip .gt. is) then
        is = is+1
        ss(is:is) = '.'
      endif
      if (ip .eq. i1) then
        i1 = i1-1
        ss(i1:i1) = '0'
      endif

C --- Find end of mantissa; get exponent (binary in ixp) ---
      iex = 0
      im = ip-1
      ixp = 0
      call chrps2(ss,'dDeE+-',6,is,im,j)
      if (j .ne. 0) then
        iex = im+2
        if (j .ge. 5) iex = iex-1
        if (ss(iex:iex) .eq. '+' .or. ss(iex:iex) .eq. '-') iex = iex+1
        do  10  j = iex, is 
          k = j
          if (ss(j:j) .ne. '0') goto 12
 10     continue
   12   continue
        if (ss(iex-1:iex-1) .eq. '-') then
          k = k-1
          ss(k:k) = '-'
        endif
        iex = k
        k = iex-1
        if (.not. a2bin(ss,ixp,2,0,' ',k,-1))
     .  call rx('pretty: parse error')
*        print *, iex,is, ':', ss(1:iex-1), '|', ss(iex:is), '|', ixp
      endif

C --- Truncate string at precsn decimal places, possibly rounding ---
      if (precsn .gt. 0) then
        k = ip+precsn+ixp
        m = i1
C ...   If relative precision, k offset from first significant digit
        if (getdig(sw,0,10) .eq. 1) then
          do  14  i = i1, im
            m = i
            if (ss(i:i) .ge. '1' .and. ss(i:i) .le. '9') goto 15
 14       continue
   15     continue
          k = m-1+precsn
          if (m .le. ip .and. k .ge. ip) k = k+1
        endif
        k = max(min(k,im),m)
*        print *, k, ' ##|', ss(i1:k), '|'
C ...   Round upwards, append exponent and start over
        if (ss(k+1:k+1) .ge. '5' .and. ss(k+1:k+1) .le. '9') then
          im = min(k,im)
          sss = ss(i1:max(im+1,ip))
          read(sss,'(e30.20)') xx
          fmt = '(f30.'
          if (im-ip .ge. 10) write(fmt(6:8),'(i2,'')'')') im-ip
          if (im-ip .lt. 10) write(fmt(6:7),'(i1,'')'')') max(im-ip,0)
C ...     Next line makes some recalcitrant machines round properly
          if (im .gt. ip) xx = xx + 4d0/10d0**(1+im-ip)
          if (im .gt. ip) write(sss,fmt) xx
          if (im .le. ip) write(sss,fmt) dround(xx,k+1-i1)
*          print *, sss, im-ip
C ...     Append exponent
          if (j .gt. 0) then
            im = 0
            call chrps2(sin,'0123456789',10,lsin,im,j)
            call chrps2(sin,'dDeE+-',6,lsin,im,j)
            call skpblb(sss,40,k)
            call strcop(sss(k+2:40),sin(im+1:lsin),lsin-im,' ',j)
*            print *, sss, im-ip          
          endif
C ...     Start over
          ss = sss
          goto 1
        endif
        if (k .lt. ip) then
          do   j = k+1, ip-1
             ss(j:j) = '0'
          enddo   
        endif
*        print *, precsn, '(precision):', ss(1:im), '|', k,im
        im = max(min(k,im),ip)
      endif

C --- Patch 0.mmmEnnn so that leading digit is nonzero ---
      if (iex .ne. 0 .and. ixp .ne. 0) then
        if (ss(i1:i1+1) .eq. '0.') then
          ixp = ixp-1
          i1 = i1+1
          ss(i1:i1) = ss(i1+1:i1+1)
          ss(i1+1:i1+1) = '.'
          ip = i1+1
*          print *, 'patch', ss(i1:is)
        endif
        sexp = ' '
        write(sexp,'(i5)') ixp
*        print *, 'exponent', sexp
      endif

C --- Copy mantissa to sout, counting decimals ---
      idec = 0
      ndig = 0
      ltmp = getdig(sw,1,10) .ne. 0
      do  20  j = i1, im
        isdig = ss(j:j) .ge. '0' .and. ss(j:j) .le. '9'
        if (j .eq. i1 .and. ltmp .and. ss(j:j) .eq. '0') goto 20
        if (j .eq. ip .or. isdig) then
          sout(nout:nout) = ss(j:j)
          nout = nout+1
          if (isdig .and. j .gt. i1) ndig = ndig+1
          if (isdig .and. j .gt. ip) idec = idec+1
        else
          goto 22
        endif
   20 continue
   22 continue
      nout = nout-1
      idec = idec-ixp
*      print *, 'string before truncation:', sout(1:nout),
*     .  '...', idec, ' decimals', ndig, ' digits; ndec=',ndec

C --- Append or strip trailing zeros, preserving ndec decimals ---
      ip = 0
      call chrpos(sout,'.',nout,ip)
      ip = ip+1
      if (ip .gt. nout) then
        ss = sin
        i1 = 0
        ls = len(ss)
        call skipbl(ss,ls,i1)
        call skpblb(ss,ls,is)
        call rx('pretty:  missing ''.'' in string'//ss(i1:is+1))
      endif
C ... Find i1 which points to 1st significant digit
      i1 = 0
      call chrps2(sout,'123456789',9,nout,i1,j)
      i1 = i1+1
C ... Case number is zero
      if (i1 .gt. nout) then
        i1 = 0
        call chrpos(sout,'0',nout,i1)
        i1 = i1+1
        if (i1 .gt. nout) call rx('pretty:  missing digit')
      endif
      id = 0
      nd = ndec
      if (precsn .gt. 0) nd = min(nd,precsn)
      if (ndec .gt. 0) then
        id = ip+nd+ixp
        if (getdig(sw,0,10) .eq. 1) id = i1+nd-1
        if (getdig(sw,0,10) .eq. 1 .and. i1 .lt. ip) id = id+1
      endif
C ... Append 
      if (nout .lt. id) then
        do   k = nout+1, id
           sout(k:k) = '0'
        enddo   
        nout = id
      endif
C ... Strip
      if (0<id .and. id<=len(sout))then
      if (sout(id:id) .eq. '.') id = id-1
      endif
      do  41  k = nout, i0, -1
        j = k
*        print *, '!!!',sout(1:j),' ',j,id
        ch = sout(j:j)
        if (j .le. id) goto 42
        if (ch .eq. '0' .or. ch .eq. '.') sout(j:j) = ' '
        if (ch .ne. '0' .or. ch .eq. '.') goto 42
   41 continue
   42 continue
C ... Fix up a completely missing mantissa
      nout = j
      if (nout .eq. i0 .and. sout(i0:i0) .eq. ' ') then
        sout(i0:i0) = '0'
        return
      endif
      if (sout(nout:nout) .eq. ' ') nout = nout-1

C --- Append exponent if there is one ---
      if (ixp .ne. 0) then
        j = 0
        call skipbl(sexp,5,j)
*        print *, j,':',sexp(j+1:5)
        sout(nout+1:nout+7-j) = 'e' // sexp(j+1:5)
        nout = nout+6-j
      endif
      
*      print *, 'final string:', sout(1:nout)
*      print *, '--------------'

      end

      double precision function dround(x,n)
C- Rounds double precision x after n digits
C     implicit none
      integer n,is,i
      double precision x,s
#if QUAD
      double precision xnint
#endif
#if QUAD
      s = qlog10(dble(abs(x)))
#else
      s = dlog10(dabs(x))
#endif
      is = s
      if (is .gt. 0) then
        is = -int(is) + n-1
      else
        is = int(-is) + n
      endif
      s = 1
      do    i = 1, iabs(is)
         s = 10*s
      enddo   
#if QUAD
      if (is .gt. 0) dround = xnint(s*dble(x))/s
      if (is .lt. 0) dround = xnint(dble(x)/s)*s
#else
      if (is .gt. 0) dround = dnint(s*x)/s
      if (is .lt. 0) dround = dnint(x/s)*s
#endif
      end
