subroutine words(str,nw) !- Count blank-delimited words in str
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   str   :string
  !o Outputs
  !o   nw    :number of blank-delimited words in str
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) str
  integer :: nw
  ! ... Local parameters
  integer :: i1,i2,i0,i

  nw = 0
  i1 = 0
  i2 = 0
  i0 = 1
99 do  10  i = i0, len(str)
     if(str(i:i) /= ' ') then
        i1 = i
        goto 90
     endif
10 enddo
  return
90 nw = nw+1
  do  20  i = i1,len(str)
     if(str(i:i) == ' ') then
        i2 = i
        goto 91
     endif
20 enddo
  return
91 i0 = i2
  goto 99
end subroutine words
subroutine word(str,iw,j1,j2)
  !- Returns j1,j2 so that str(j1:j2) is the iw-th word from beginning
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   str   :string
  !i   iw    :find iw-th word
  !o Outputs
  !o   j1    :str(j1:j2) is iw-th word
  !o   j2    :-//-
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character*(*) str
  integer :: iw,j1,j2
  ! ... External calls
  external nword
  j1 = 1
  call nword(str,iw,j1,j2)
end subroutine word

subroutine nword(str,iw,j1,j2)
  !- Returns j1,j2 so that str(j1:j2) is the iw-th word from current pos
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   str   :string
  !i   iw    :find iw-th word
  !i   j1    :start search from str(j1:)
  !o Outputs
  !o   j1    :str(j1:j2) is iw-th word
  !o   j2    :-//-
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: iw,j1,j2
  character*(*) str
  ! ... Local parameters
  integer :: nw,i1,i2,i0,i
  nw = 0
  i1 = 0
  i2 = 0
  i0 = j1
  j2 = -1
99 do  10  i = i0, len(str)
     !   ... skip until nonblank char
     if(str(i:i) /= ' ') then
        i1 = i
        goto 90
     endif
10 enddo
  return
  !   ... skip until a blank char
90 nw = nw+1
  if (nw == iw) j1 = i1
  do  20  i = i1, len(str)
     if(str(i:i) == ' ') then
        i2 = i
        goto 91
     endif
20 enddo
  ! ... We have reached the end of the string
  if (nw == iw) j2 = len(str)
  return
  ! ... cleanup: exit if word sought, else try again
91 i0 = i2
  if (nw == iw) then
     j2 = i2-1
     return
  endif
  goto 99
end subroutine nword

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
