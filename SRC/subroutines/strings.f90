subroutine strcop(dest,source,len,term,i)
  !- copy one string to another
  ! ----------------------------------------------------------------
  !i Inputs/Outputs
  !i   dest,source: source and destination strings
  !i   len:   maximum number of characters to copy
  !i   term:  terminator
  !o   i:     number of characters copied (including terminator)
  !r Remarks
  !r   string copy continues until term encountered or
  !r   len characters are checked.
  ! ----------------------------------------------------------------
  implicit none
  integer :: len
  character(1) :: dest(len),source(len),term
  integer :: i
  i = 0
  if (len == 0) return
  do
     i = i+1
     dest(i) = source(i)
     if (dest(i) /= term .AND. i < len) cycle
     exit
  enddo
end subroutine strcop
subroutine strncp(dest,src,idest,isrc,len)  !- Copy a string of given length from src to dest.
  integer :: idest,isrc,len
  integer :: i
  character*(1) dest(*),src(*)
  do    i = 0, len-1
     dest(idest+i) = src(isrc+i)
  enddo
end subroutine strncp
subroutine chrps2(s,ch,nch,maxch,ich,iterm)  !- Finds position of any of a set of characters in string
  ! ----------------------------------------------------------------
  !i Inputs
  !i   s:   string (declared as s(0:*), following C conventions)
  !i   ch,nch:  set of characters sought, and number in set
  !i   ich: start search at s(ich)
  !i   maxch: see ich
  !o Outputs
  !o   ich: position of character ch, not to exceed maxch
  !o   iterm: index to char in ch that terminated search (0 if none)
  !r Remarks
  !r    seeks match at string(i0), string(i0+1) until ch is found or until
  !r    ich = maxch.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: ich,maxch,iterm,nch
  character(1) :: ch(nch),s(0:*)
  integer :: i
!  write(6,*)'mmmmixrho ccccc',s(ich:maxch)
  iterm = 0
10 do  20  i = 1, nch
     if (s(ich) /= ch(i)) goto 20
     iterm = i
     return
20 enddo
  if (ich >= maxch) return
  ich = ich+1
  goto 10
end subroutine chrps2
subroutine skipbl(t,nt,i)  !- Parses string for first nonblank character
  !     implicit none
  integer :: nt,i
  character(1) :: t(0:nt)
  if (i >= nt) return
99 if (t(i) /= ' ') return
  i = i+1
  if (i >= nt) return
  goto 99
end subroutine skipbl
