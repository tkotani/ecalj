      subroutine wrdsg(str,mode,sep,nw)
C- Counts words in str, where word is any char not in 'sep'
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   mode  : 1s digit nonzero => skip blanks before start of new word
Ci         :10s digit nonzero => word is any char in sep
Ci   sep   :string containing list of word delimiters
Co Outputs
Co   nw    :number of sep-delimited words in str
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) str,sep
      integer mode,nw
C ... Local parameters
      integer ifblk,ifinc,i0,i,lstr,lsep,iterm,ip
      character*1 cl,ch
C ... External calls
Ckino 2011.05 del external,       external chrps2,skipbl

      nw = 0
      i0 = 1
      lstr = len(str)
      lsep = len(sep)
      ifblk = mod(mode,10)
      ifinc = mod(mode/10,10)
   99 i = i0-1
      if (ifblk .ne. 0) call skipbl(str,lstr-1,i)
      if (i .ge. lstr) return
      if (ifinc .eq. 0) then
        call chrps2(str,sep,lsep,lstr-1,i,iterm)
        i0 = i+2
      else
   10   ip = 0
        i = i+1
   12   ip = ip+1
C       ip2 = ip
        if (ip .le. lsep) then
          cl = sep(ip:ip)
          ch = cl
          if (ip+2 .le. lsep) then
            if (sep(ip+1:ip+1) .eq. '-') ch = sep(ip+2:ip+2)
            ip = ip+2
          elseif (ip+1 .eq. lsep) then
            if (sep(ip+1:ip+1) .eq. '-') then
              ch = '-'
              ip = ip+1
            endif
          endif
C     ... Case this character belongs to a word
          if (cl .le. str(i:i) .and. ch .ge. str(i:i)) then
            goto 10
          endif
*          print *, i,str(i:i),cl,ch,ip
          goto 12
        endif
C   ... This character not in lsep; end of word encountered
        i0 = i+1
      endif
      nw = nw+1
      goto 99
      end
      subroutine wordg(str,mode,sep,iw,i1,i2)
C- Returns i1,i2 so that str(i1:i2) is the iw-th word from beginning,
C- where word is any sequence of characters not in 'sep'
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   mode  :  1s digit nonzero => skip blanks before start of new word
Ci         : 10s digit nonzero => delimiter is any char not in sep
Ci         :100s if iw'th word is not found, return i2=-1
Ci         :     (not checked if 10s digit of mode also set)
Ci   sep   :string containing list of word delimiters
Ci   iw    :find iw-th word
Co Outputs
Co   i1    :str(i1:i2) is iw-th word
Co   i2    :-//-
Cl Local variables
Cr Remarks
Cb Bugs
Cb    mode=100 not checked when used in conjunction with mode=10
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) str,sep
      integer mode,iw,i1,i2
      i1 = 1
      call nwordg(str,mode,sep,iw,i1,i2)
      end

      subroutine nwordg(str,mode,sep,iw,i1,i2)
C- Returns i1,i2 so that str(i1:i2) is the iw-th word from current pos,
C- where word is any sequence of characters not in 'sep'
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   mode  :  1s digit nonzero => skip blanks before start of new word
Ci         : 10s digit nonzero => delimiter is any char not in sep
Ci         :100s if iw'th word is not found, return i2=-1
Ci   sep   :string containing list of word delimiters
Ci   iw    :find iw-th word
Ci   i1    :start search from str(i1:)
Co Outputs
Co   i1    :str(i1:i2) is iw-th word
Co   i2    :-//-
Cl Local variables
Cr Remarks
Cb Bugs
Cb    mode=100 not checked when used in conjunction with mode=10
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) str,sep
      integer mode,iw,i1,i2
C ... Local parameters
      integer ifblk,ifinc,nw,i0,i,lstr,lsep,iterm,ip
      character*1 cl,ch
C ... External calls
Ckino 2011.05 del external,       external chrps2,skipbl

      nw = 0
      i0 = i1
      i2 = 0
      lstr = len(str)
      lsep = len(sep)
      ifblk = mod(mode,10)
      ifinc = mod(mode/10,10)
C ... Start of next search for word.  i is offset for str(0:)
   99 i = i0-1
      if (ifblk .ne. 0) call skipbl(str,lstr-1,i)
      if (i .ge. lstr) return
      if (nw+1 .eq. iw) i1 = i+1
C ... Case word is any char not in sep
      if (ifinc .eq. 0) then
C   ... Find the next separator
        call chrps2(str,sep,lsep,lstr-1,i,iterm)
C   ... If at the end, go past the end
        if (iterm .eq. 0 .and. mod(mode/100,10) .ne. 0) then
          i2 = -1
          return
        endif
        if (i .eq. lstr-1 .and. iterm .eq. 0) i=i+1
        i0 = i+2
C       print '(a)', str(1:i+1)
C ... Case word is any char in sep
      else
   10   ip = 0
        i = i+1
   12   ip = ip+1
C       ip2 = ip
        if (ip .le. lsep) then
          cl = sep(ip:ip)
          ch = cl
          if (ip+2 .le. lsep) then
            if (sep(ip+1:ip+1) .eq. '-') then
              ch = sep(ip+2:ip+2)
              ip = ip+2
            endif
C          elseif (ip+1 .eq. lsep) then
C            if (sep(ip+1:ip+1) .eq. '-') then
C              ch = '-'
C              ip = ip+1
C            endif
          endif
C     ... Case this character belongs to a word
          if (cl .le. str(i:i) .and. ch .ge. str(i:i)) then
            goto 10
          endif
*         print *, i,str(i:i),cl,ch,ip
          goto 12
        endif
C   ... This character not in lsep; end of word encountered
        i0 = i+1
*       print *, 'last char',i,'nw=',nw+1,' ',str(1:i)
      endif
      nw = nw+1
      if (nw .eq. iw) i2 = i0-2
      if (nw .lt. iw) goto 99

      end

