subroutine words(str,nw)
  !- Count blank-delimited words in str
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

