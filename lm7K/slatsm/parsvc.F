      logical function parsvc(iopt,t,ip,v)
C- Parses string, converting to a vector, or vice-versa
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iopt   -1 to convert string t to vector,
Ci           1 to convert vector v to string t
Ci           2 same as 1, but use shorthand 'D', 'X', 'Y', 'Z'.
Ci   ip      position in t first char (0 for first char)
Ci Inputs/Outputs:
Cio  t       string representation of vector (see Remarks)
Cio  v       vector
Co Outputs:
Co   ip      position after last char: t(ip-1) is last char
Co   parsvc  T if parse was successful; F if not.
Cr Remarks
Cr   The general string represntation has the form '(a,b,c)'
Cr   where a, b and c are expressions.
Cr   The following shorthand notations are allowed:
Cr   'D' for  (1,1,1)
Cr   'X', 'Y', 'Z' for (1,0,0), (0,1,0), (0,0,1).
Cb Bugs
Cb   no check is made on the length of t
C ----------------------------------------------------------------------
C     implicit none
      integer ip
      double precision v(3)
      character*(1) t(0:1)
C Local variables
      double precision tiny,x,y,z,d
      character rchr*9, sout*50
      integer itrm,a2vec,awrite,ix(3),ich,iopt,m,i,iz,id
      logical lveq0(3),lveq1(3),a2bin
      parameter (tiny=1d-4)

      data rchr /'(XxYyZzDd'/

C --- Convert t to vec ---
      if (iopt .eq. -1) then
        call dpzero(v,3)
        ich = 0
        call chrps2(t(ip),rchr,len(rchr),ich,ich,itrm)
C   ... First char not in rchr: not a recognizable vector
        parsvc = .false.
        if (itrm .eq. 0) return
C   ... '('
        if (itrm .eq. 1) then
          ip = ip+1
c this doesn't handle complicated cases
c         parsvc = a2vec(t,ip+100,ip,4,',)',2,3,3,ix,v) .eq. 3
          if (.not. a2bin(t,v,4,0,',',ip,-1)) return
          if (.not. a2bin(t,v,4,1,',',ip,-1)) return
          if (.not. a2bin(t,v,4,2,')',ip,-1)) return
          parsvc = .true.
          return
        else
C         ... 'd'
          if (itrm .ge. 8) then
            v(1) = 1
            v(2) = 1
            v(3) = 1
C         ... 'x' 'y' or 'z'
          else
            v(itrm/2) = 1
          endif
          parsvc = .true.
          ip = ip+1
        endif
        return

C --- Convert vec to t ---
      elseif (mod(iopt,4) .eq. 1 .or. mod(iopt,4) .eq. 2) then
        parsvc = .true.
        do  10  i = 1, 3
          lveq0(i) = dabs(v(i)) .lt. tiny
   10   lveq1(i) = dabs(v(i)-1) .lt. tiny
C   ... Exclude shorthand
        if (mod(iopt,4) .eq. 1) then
          lveq0(3) = .false.
          lveq1(3) = .false.
        endif
        if     (lveq1(1).and.lveq1(2) .and. lveq1(3)) then
          t(ip) = 'd'
        elseif (lveq1(1).and.lveq0(2) .and. lveq0(3)) then
          t(ip) = 'x'
        elseif (lveq0(1).and.lveq1(2) .and. lveq0(3)) then
          t(ip) = 'y'
        elseif (lveq0(1).and.lveq0(2) .and. lveq1(3)) then
          t(ip) = 'z'
        else
          t(ip) = '('
          do  12  i = 1, 3
            m = awrite('%x%;7d%?#n<>3#,#)#',
     .      sout,len(sout),0,v(i),i,i,i,i,i,i,i)
            x = v(i)
            y = v(i)*dsqrt(3d0)*4
            if (abs(x) .gt. tiny .and. dabs(dabs(x)-1) .gt. tiny
     .      .and. iopt .le. 4) then
              if (abs(1/x-nint(1/x)) .lt. tiny) then
                if (x .gt. 0) m = awrite('%x1/%d%?#n<>3#,#)#',
     .          sout,len(sout),0,1/x,i,i,i,i,i,i,i)
                if (x .lt. 0) m = awrite('%x-1/%d%?#n<>3#,#)#',
     .          sout,len(sout),0,-1/x,i,i,i,i,i,i,i)
              elseif (abs(y-nint(y)) .lt. tiny .and.
     .        abs(y) .gt. .5d0) then
                d = 12
                iz = nint(abs(y))
                id = 12
                if (iz .eq. 2) id = 6
                if (iz .eq. 4) id = 3
                if (iz .eq. 6) id = 2
                if (iz .eq. 12) id = 1
                y = y/(12/id)
                m = awrite('%x%?#n==1#%j#%d*#sqrt(3)/%i%?#n<>3#,#)#',
     .          sout,len(sout),0,nint(y),y,id,i,i,i,i,i)
              endif
            endif
            call strncp(t(ip+1),sout,1,1,m)
            ip = ip+m
   12     continue
        endif
        ip = ip+1
      endif
      end
C      subroutine fmain
C      implicit none
C      double precision v(3)
C      integer ip
C      logical parsvc,lsw
C      character *50 t
C
C
C      call dpzero(v,3)
C      v(1) = 1
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      call dpzero(v,3)
C      v(2) = 1
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      call dpzero(v,3)
C      v(3) = 1
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      call dpzero(v,3)
C      v(1) = 1
C      v(2) = 1
C      v(3) = 1
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      call dpzero(v,3)
C      v(1) = 1d0/3
C      v(2) = .25d0
C      v(3) = .3d0
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t//'%a '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      call dpzero(v,3)
C      v(1) = 1d0/3
C      v(2) = .25d0
C      v(3) = .3d0
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(5,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1;7d '//t//'%a '//t(ip:ip+1),
C     .  ' ',80,6,lsw,ip,v)
C
C      call dpzero(v,3)
C      v(1) = 5*sqrt(3d0)/12
C      v(2) = -sqrt(3d0)/2
C      v(3) = sqrt(3d0)/6
C      ip = 1
C      t = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
C      lsw = parsvc(1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1;7d '//t//'%a '//t(ip:ip+1),
C     .  ' ',95,6,lsw,ip,v)
C
C      t = ' xa'
C      ip = 1
C      lsw = parsvc(-1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      t = ' yb'
C      ip = 1
C      lsw = parsvc(-1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      t = ' zc'
C      ip = 1
C      lsw = parsvc(-1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      t = ' Dd'
C      ip = 1
C      lsw = parsvc(-1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      t = ' (.1,2-1,pi/4)e'
C      ip = 1
C      lsw = parsvc(-1,t,ip,v)
C      call awrit3(' lsw=%l  ip=%i  v=%3:1d '//t(ip:ip+1),' ',
C     .  80,6,lsw,ip,v)
C
C      end

