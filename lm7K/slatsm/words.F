      subroutine words(str,nw)
C- Count blank-delimited words in str
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Co Outputs
Co   nw    :number of blank-delimited words in str
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) str
      integer nw
C ... Local parameters
      integer i1,i2,i0,i

      nw = 0
      i1 = 0
      i2 = 0
      i0 = 1
   99 do  10  i = i0, len(str)
        if(str(i:i) .ne. ' ') then
          i1 = i
          goto 90
        endif
   10 continue
      return
   90 nw = nw+1
      do  20  i = i1,len(str)
        if(str(i:i) .eq. ' ') then
          i2 = i
          goto 91
        endif
   20 continue
      return
   91 i0 = i2
      goto 99
      end
      subroutine word(str,iw,j1,j2)
C- Returns j1,j2 so that str(j1:j2) is the iw-th word from beginning
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   iw    :find iw-th word
Co Outputs
Co   j1    :str(j1:j2) is iw-th word
Co   j2    :-//-
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) str
      integer iw,j1,j2
C ... External calls
Ckino 2011.05 del external,       external nword
      j1 = 1
      call nword(str,iw,j1,j2)
      end

      subroutine nword(str,iw,j1,j2)
C- Returns j1,j2 so that str(j1:j2) is the iw-th word from current pos
C ----------------------------------------------------------------------
Ci Inputs
Ci   str   :string
Ci   iw    :find iw-th word
Ci   j1    :start search from str(j1:)
Co Outputs
Co   j1    :str(j1:j2) is iw-th word
Co   j2    :-//-
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer iw,j1,j2
      character*(*) str
C ... Local parameters
      integer nw,i1,i2,i0,i
      nw = 0
      i1 = 0
      i2 = 0
      i0 = j1
      j2 = -1
  99  do  10  i = i0, len(str)
C   ... skip until nonblank char
        if(str(i:i) .ne. ' ') then
          i1 = i
          goto 90
        endif
  10  continue
      return
C   ... skip until a blank char
  90  nw = nw+1
      if (nw .eq. iw) j1 = i1
      do  20  i = i1, len(str)
        if(str(i:i) .eq. ' ') then
          i2 = i
          goto 91
        endif
   20 continue
C ... We have reached the end of the string
      if (nw .eq. iw) j2 = len(str)
      return
C ... cleanup: exit if word sought, else try again
   91 i0 = i2
      if (nw .eq. iw) then
        j2 = i2-1
        return
      endif
      goto 99
      end

