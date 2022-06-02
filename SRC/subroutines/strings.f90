logical function lsequ(s1,s2,len,term,i)
  !- Determine whether two strings are equal or not
  ! ----------------------------------------------------------------
  !i Inputs
  !i   s1,s2: strings to compare
  !i   len:   maximum length of string
  !i   term:  terminator
  !o Outputs
  !o   lsequ: returned true or false
  !o   i:     number of characters tokened (including terminator)
  !r Remarks
  !r   string comparison continues until terminator encountered or
  !r   len characters are checked.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: i,len
  character(1) :: s1(len),s2(len),term
  !     lsequ = .true.
  !     do  10  i = 1, len
  !       lsequ = (lsequ .and. s1(i) .eq. s2(i))
  !       if (s1(i) .eq. term) return
  !  10 continue
  lsequ = .false.
  do  10  i = 1, len
     if (s1(i) /= s2(i)) return
     if (s1(i) == term) goto 15
10 enddo
15 lsequ = .true.
  return
end function lsequ
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
  !     implicit none
  integer :: len
  character(1) :: dest(len),source(len),term
  integer :: i
  i = 0
  if (len == 0) return
10 i = i+1
  dest(i) = source(i)
  if (dest(i) /= term .AND. i < len) goto 10
end subroutine strcop
subroutine strcat(s1,len1,term1,s2,len2,term2,i)
  !- concatenate one string to another
  ! ----------------------------------------------------------------
  !i Inputs/Outputs
  !i   s1,s2: source string and string to concatenate
  !i   len1:  maximum length of s1
  !i   term1: terminator for s1
  !i   len2:  maximum length of s2
  !i   term2: terminator for s2
  !o Outputs
  !o   i:     number of characters in s1 (including terminator)
  !r Remarks
  !r   concatenation continues until term encountered or
  !r   len2 characters are concatenated.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: len1,len2
  character(1) :: s1(len1),s2(len2),term1,term2
  integer :: i,j

  i = 0
10 i = i+1
  if (s1(i) /= term1 .AND. i < len1) goto 10
  j = 0
  if (s1(i) == term1) i = i-1
20 j = j+1
  i = i+1
  s1(i) = s2(j)
  if (s2(j) /= term2 .AND. j < len2) goto 20
end subroutine strcat
subroutine strncp(dest,src,idest,isrc,len)
  !- Copy a string of given length from src to dest.
  !  No check is made on lengths of strings.
  !     implicit none
  integer :: idest,isrc,len
  integer :: i
  character*(1) dest(*),src(*)
  do    i = 0, len-1
     dest(idest+i) = src(isrc+i)
  enddo
end subroutine strncp
logical function parstr(s1,s2,reclen,len,term,i,j)
  !- Find a substring within a given string
  ! ----------------------------------------------------------------
  !i Inputs
  !i   s1: string in which to seek match
  !i   s2: match string
  !i   reclen: (maximum size of s1) - len; see Bugs, Remarks
  !i   len: size of match string
  !i   i: character where search within string should begin
  !i   term: character that terminates the match string
  !o Outputs
  !o   i: index to first position of token.
  !o      NB: this routine follows C conventions; first char at i=0
  !o   j: index to first position after token
  !o      NB: this routine follows C conventions; first char at j=0
  !b Bugs
  !b   reclen is inappropriately named, since the size of the
  !b   match string is not taken into account.
  !r Remarks
  !r    seeks match at s1(i), s1(i+1), ... s1(reclen-1) until
  !r    match is found.  returns false if no match found.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: reclen,len,i,j
  character(1) :: s1(0:*),s2(0:*),term
  logical :: lsequ

  parstr = .false.
  i = i-1
10 i = i+1
  if (i >= reclen) return
  if ( .NOT. lsequ(s1(i),s2,len,term,j)) goto 10
  parstr = .true.
  j = j + i
end function parstr
subroutine chrpos(s,ch,maxch,ich)
  !- Finds position of character in string
  ! ----------------------------------------------------------------
  !i Inputs
  !i   s:   string (declared as s(0:*), following C conventions)
  !i   ch:  character sought
  !i   ich: start search at s(ich)
  !i   maxch: see ich
  !o Outputs
  !o   ich: position of character ch, not to exceed maxch
  !r Remarks
  !r    seeks match at string(i0), string(i0+1) until ch is found or until
  !r    ich = maxch.
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: ich,maxch
  character(1) :: ch,s(0:*)

10 if (ich == maxch  .OR.  s(ich) == ch) return
  ich = ich+1
  goto 10
end subroutine chrpos
subroutine chrps2(s,ch,nch,maxch,ich,iterm)
  !- Finds position of any of a set of characters in string
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
subroutine skipbl(t,nt,i)
  !- Parses string for first nonblank character
  !     implicit none
  integer :: nt,i
  character(1) :: t(0:nt)
  if (i >= nt) return
99 if (t(i) /= ' ') return
  i = i+1
  if (i >= nt) return
  goto 99
end subroutine skipbl
subroutine skpblb(t,nt,i)
  !- Parses string for first nonblank character, right to left
  !     implicit none
  integer :: nt,i
  character(1) :: t(0:*)
  i = nt
99 i = i-1
  if (i < 0) return
  if (t(i) /= ' ') return
  goto 99
end subroutine skpblb
subroutine skp2bl(t,nt,i)
  !- Parses string for first blank character
  !     implicit none
  integer :: nt,i
  character(1) :: t(0:nt)
  if (i >= nt) return
99 if (t(i) == ' ') return
  i = i+1
  if (i < nt) goto 99
end subroutine skp2bl
subroutine tokmat(string,token,n,len,term,itoken,ltoken,lopt)
  !- compare a string to a list of strings
  ! ----------------------------------------------------------------
  !i Inputs
  !i   string: test string
  !i   token: vector of strings to compare
  !i   n,len: number and length of strings in token
  !i   term : string terminator
  !i   lopt:  if true, tokmat stops with error message when no match found
  !o Outputs
  !o   itoken: index to string matched, -1 if none
  !o   ltoken: length of token matched
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: n,len,itoken,ltoken,i
  logical :: lsequ,lopt
  character(1) :: string(*), term
  character*(*) token(0:*)
  do  10  itoken = 0, n-1
     if (lsequ(string,token(itoken),len,term,ltoken)) return
10 enddo
  itoken = -1
  if (lopt) then
     print *, 'TOKMAT: unmatched ', (string(i), i=1,len)
     call cexit(-1,1)
  endif
end subroutine tokmat
!      subroutine fmain
!      character*42 s
!      integer i,j
!      i=0
!      s = 'hello there world'
!      call skipbl(s,len(s),i)
!      print *, i
!      call skp2bl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skp2bl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skipbl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skipbl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skp2bl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skipbl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skp2bl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skipbl(s,len(s),i)
!      print *, i, s(1:i), '|'
!      call skp2bl(s,len(s),i)
!      print *, i, s(1:i), '|'

!      end

