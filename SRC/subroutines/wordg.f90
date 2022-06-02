subroutine wrdsg(str,mode,sep,nw)
  !- Counts words in str, where word is any char not in 'sep'
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   str   :string
  !i   mode  : 1s digit nonzero => skip blanks before start of new word
  !i         :10s digit nonzero => word is any char in sep
  !i   sep   :string containing list of word delimiters
  !o Outputs
  !o   nw    :number of sep-delimited words in str
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) str,sep
  integer :: mode,nw
  ! ... Local parameters
  integer :: ifblk,ifinc,i0,i,lstr,lsep,iterm,ip
  character(1) :: cl,ch
  ! ... External calls
  external chrps2,skipbl

  nw = 0
  i0 = 1
  lstr = len(str)
  lsep = len(sep)
  ifblk = mod(mode,10)
  ifinc = mod(mode/10,10)
99 i = i0-1
  if (ifblk /= 0) call skipbl(str,lstr-1,i)
  if (i >= lstr) return
  if (ifinc == 0) then
     call chrps2(str,sep,lsep,lstr-1,i,iterm)
     i0 = i+2
  else
10   ip = 0
     i = i+1
12   ip = ip+1
     !       ip2 = ip
     if (ip <= lsep) then
        cl = sep(ip:ip)
        ch = cl
        if (ip+2 <= lsep) then
           if (sep(ip+1:ip+1) == '-') ch = sep(ip+2:ip+2)
           ip = ip+2
        elseif (ip+1 == lsep) then
           if (sep(ip+1:ip+1) == '-') then
              ch = '-'
              ip = ip+1
           endif
        endif
        !     ... Case this character belongs to a word
        if (cl <= str(i:i) .AND. ch >= str(i:i)) then
           goto 10
        endif
        !          print *, i,str(i:i),cl,ch,ip
        goto 12
     endif
     !   ... This character not in lsep; end of word encountered
     i0 = i+1
  endif
  nw = nw+1
  goto 99
end subroutine wrdsg
subroutine wordg(str,mode,sep,iw,i1,i2)
  !- Returns i1,i2 so that str(i1:i2) is the iw-th word from beginning,
  !- where word is any sequence of characters not in 'sep'
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   str   :string
  !i   mode  :  1s digit nonzero => skip blanks before start of new word
  !i         : 10s digit nonzero => delimiter is any char not in sep
  !i         :100s if iw'th word is not found, return i2=-1
  !i         :     (not checked if 10s digit of mode also set)
  !i   sep   :string containing list of word delimiters
  !i   iw    :find iw-th word
  !o Outputs
  !o   i1    :str(i1:i2) is iw-th word
  !o   i2    :-//-
  !l Local variables
  !r Remarks
  !b Bugs
  !b    mode=100 not checked when used in conjunction with mode=10
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) str,sep
  integer :: mode,iw,i1,i2
  i1 = 1
  call nwordg(str,mode,sep,iw,i1,i2)
end subroutine wordg

subroutine nwordg(str,mode,sep,iw,i1,i2)
  !- Returns i1,i2 so that str(i1:i2) is the iw-th word from current pos,
  !- where word is any sequence of characters not in 'sep'
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   str   :string
  !i   mode  :  1s digit nonzero => skip blanks before start of new word
  !i         : 10s digit nonzero => delimiter is any char not in sep
  !i         :100s if iw'th word is not found, return i2=-1
  !i   sep   :string containing list of word delimiters
  !i   iw    :find iw-th word
  !i   i1    :start search from str(i1:)
  !o Outputs
  !o   i1    :str(i1:i2) is iw-th word
  !o   i2    :-//-
  !l Local variables
  !r Remarks
  !b Bugs
  !b    mode=100 not checked when used in conjunction with mode=10
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) str,sep
  integer :: mode,iw,i1,i2
  ! ... Local parameters
  integer :: ifblk,ifinc,nw,i0,i,lstr,lsep,iterm,ip
  character(1) :: cl,ch
  ! ... External calls
  external chrps2,skipbl

  nw = 0
  i0 = i1
  i2 = 0
  lstr = len(str)
  lsep = len(sep)
  ifblk = mod(mode,10)
  ifinc = mod(mode/10,10)
  ! ... Start of next search for word.  i is offset for str(0:)
99 i = i0-1
  if (ifblk /= 0) call skipbl(str,lstr-1,i)
  if (i >= lstr) return
  if (nw+1 == iw) i1 = i+1
  ! ... Case word is any char not in sep
  if (ifinc == 0) then
     !   ... Find the next separator
     call chrps2(str,sep,lsep,lstr-1,i,iterm)
     !   ... If at the end, go past the end
     if (iterm == 0 .AND. mod(mode/100,10) /= 0) then
        i2 = -1
        return
     endif
     if (i == lstr-1 .AND. iterm == 0) i=i+1
     i0 = i+2
     !       print '(a)', str(1:i+1)
     ! ... Case word is any char in sep
  else
10   ip = 0
     i = i+1
12   ip = ip+1
     !       ip2 = ip
     if (ip <= lsep) then
        cl = sep(ip:ip)
        ch = cl
        if (ip+2 <= lsep) then
           if (sep(ip+1:ip+1) == '-') then
              ch = sep(ip+2:ip+2)
              ip = ip+2
           endif
           !          elseif (ip+1 .eq. lsep) then
           !            if (sep(ip+1:ip+1) .eq. '-') then
           !              ch = '-'
           !              ip = ip+1
           !            endif
        endif
        !     ... Case this character belongs to a word
        if (cl <= str(i:i) .AND. ch >= str(i:i)) then
           goto 10
        endif
        !         print *, i,str(i:i),cl,ch,ip
        goto 12
     endif
     !   ... This character not in lsep; end of word encountered
     i0 = i+1
     !       print *, 'last char',i,'nw=',nw+1,' ',str(1:i)
  endif
  nw = nw+1
  if (nw == iw) i2 = i0-2
  if (nw < iw) goto 99

end subroutine nwordg

