      logical function lsequ(s1,s2,len,term,i)
C- Determine whether two strings are equal or not
C ----------------------------------------------------------------
Ci Inputs
Ci   s1,s2: strings to compare
Ci   len:   maximum length of string
Ci   term:  terminator
Co Outputs
Co   lsequ: returned true or false
Co   i:     number of characters tokened (including terminator)
Cr Remarks
Cr   string comparison continues until terminator encountered or
Cr   len characters are checked.
C ----------------------------------------------------------------
C     implicit none
      integer i,len
      character*1 s1(len),s2(len),term
c     lsequ = .true.
c     do  10  i = 1, len
c       lsequ = (lsequ .and. s1(i) .eq. s2(i))
c       if (s1(i) .eq. term) return
c  10 continue
      lsequ = .false.
      do  10  i = 1, len
        if (s1(i) .ne. s2(i)) return
        if (s1(i) .eq. term) goto 15
   10 continue
   15 lsequ = .true.
      return
      end
      subroutine strcop(dest,source,len,term,i)
C- copy one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   dest,source: source and destination strings
Ci   len:   maximum number of characters to copy
Ci   term:  terminator
Co   i:     number of characters copied (including terminator)
Cr Remarks
Cr   string copy continues until term encountered or
Cr   len characters are checked.
C ----------------------------------------------------------------
C     implicit none
      integer len
      character*1 dest(len),source(len),term
      integer i
      i = 0
      if (len .eq. 0) return
   10 i = i+1
      dest(i) = source(i)
      if (dest(i) .ne. term .and. i .lt. len) goto 10
      end
      subroutine strcat(s1,len1,term1,s2,len2,term2,i)
C- concatenate one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   s1,s2: source string and string to concatenate
Ci   len1:  maximum length of s1
Ci   term1: terminator for s1
Ci   len2:  maximum length of s2
Ci   term2: terminator for s2
Co Outputs
Co   i:     number of characters in s1 (including terminator)
Cr Remarks
Cr   concatenation continues until term encountered or
Cr   len2 characters are concatenated.
C ----------------------------------------------------------------
C     implicit none
      integer len1,len2
      character*1 s1(len1),s2(len2),term1,term2
      integer i,j

      i = 0
   10 i = i+1
      if (s1(i) .ne. term1 .and. i .lt. len1) goto 10
      j = 0
      if (s1(i) .eq. term1) i = i-1
   20 j = j+1
      i = i+1
      s1(i) = s2(j)
      if (s2(j) .ne. term2 .and. j .lt. len2) goto 20
      end
      subroutine strncp(dest,src,idest,isrc,len)
C- Copy a string of given length from src to dest.
C  No check is made on lengths of strings.
C     implicit none
      integer idest,isrc,len
      integer i
      character*(1) dest(*),src(*)
      do  10  i = 0, len-1
   10 dest(idest+i) = src(isrc+i)
      end
      logical function parstr(s1,s2,reclen,len,term,i,j)
C- Find a substring within a given string
C ----------------------------------------------------------------
Ci Inputs
Ci   s1: string in which to seek match
Ci   s2: match string
Ci   reclen: (maximum size of s1) - len; see Bugs, Remarks
Ci   len: size of match string
Ci   i: character where search within string should begin
Ci   term: character that terminates the match string
Co Outputs
Co   i: index to first position of token. 
Co      NB: this routine follows C conventions; first char at i=0
Co   j: index to first position after token
Co      NB: this routine follows C conventions; first char at j=0
Cb Bugs
Cb   reclen is inappropriately named, since the size of the
Cb   match string is not taken into account.
Cr Remarks
Cr    seeks match at s1(i), s1(i+1), ... s1(reclen-1) until
Cr    match is found.  returns false if no match found.
C ----------------------------------------------------------------
C     implicit none
      integer reclen,len,i,j
      character*1 s1(0:*),s2(0:*),term
      logical lsequ

      parstr = .false.
      i = i-1
   10 i = i+1
      if (i .ge. reclen) return
      if (.not. lsequ(s1(i),s2,len,term,j)) goto 10
      parstr = .true.
      j = j + i
      end
      subroutine chrpos(s,ch,maxch,ich)
C- Finds position of character in string
C ----------------------------------------------------------------
Ci Inputs
Ci   s:   string (declared as s(0:*), following C conventions)
Ci   ch:  character sought
Ci   ich: start search at s(ich)
ci   maxch: see ich
Co Outputs
Co   ich: position of character ch, not to exceed maxch
Cr Remarks
Cr    seeks match at string(i0), string(i0+1) until ch is found or until
Cr    ich = maxch.
C ----------------------------------------------------------------
C     implicit none
      integer ich,maxch
      character*1 ch,s(0:*)

   10 if (ich .eq. maxch  .or.  s(ich) .eq. ch) return
      ich = ich+1
      goto 10
      end
      subroutine chrps2(s,ch,nch,maxch,ich,iterm)
C- Finds position of any of a set of characters in string
C ----------------------------------------------------------------
Ci Inputs
Ci   s:   string (declared as s(0:*), following C conventions)
Ci   ch,nch:  set of characters sought, and number in set
Ci   ich: start search at s(ich)
ci   maxch: see ich
Co Outputs
Co   ich: position of character ch, not to exceed maxch
Co   iterm: index to char in ch that terminated search (0 if none)
Cr Remarks
Cr    seeks match at string(i0), string(i0+1) until ch is found or until
Cr    ich = maxch.
C ----------------------------------------------------------------
C     implicit none
      integer ich,maxch,iterm,nch
      character*1 ch(nch),s(0:*)
      integer i

      iterm = 0
   10 do  20  i = 1, nch
        if (s(ich) .ne. ch(i)) goto 20
        iterm = i
        return
   20 continue
      if (ich .ge. maxch) return
      ich = ich+1
      goto 10
      end
      subroutine skipbl(t,nt,i)
C- Parses string for first nonblank character
C     implicit none
      integer nt,i
      character*1 t(0:nt)
      if (i .ge. nt) return
   99 if (t(i) .ne. ' ') return
      i = i+1
      if (i .ge. nt) return
      goto 99
      end
      subroutine skpblb(t,nt,i)
C- Parses string for first nonblank character, right to left
C     implicit none
      integer nt,i
      character*1 t(0:*)
      i = nt
   99 i = i-1
      if (i .lt. 0) return
      if (t(i) .ne. ' ') return
      goto 99
      end
      subroutine skp2bl(t,nt,i)
C- Parses string for first blank character
C     implicit none
      integer nt,i
      character*1 t(0:nt)
      if (i .ge. nt) return
   99 if (t(i) .eq. ' ') return
      i = i+1
      if (i .lt. nt) goto 99
      end
      subroutine tokmat(string,token,n,len,term,itoken,ltoken,lopt)
C- compare a string to a list of strings
C ----------------------------------------------------------------
Ci Inputs
Ci   string: test string
Ci   token: vector of strings to compare
Ci   n,len: number and length of strings in token
Ci   term : string terminator
Ci   lopt:  if true, tokmat stops with error message when no match found
Co Outputs
Co   itoken: index to string matched, -1 if none
Co   ltoken: length of token matched
C ----------------------------------------------------------------
C     implicit none
      integer n,len,itoken,ltoken,i
      logical lsequ,lopt
      character*1 string(*), term
      character*(*) token(0:*)
      do  10  itoken = 0, n-1
        if (lsequ(string,token(itoken),len,term,ltoken)) return
   10 continue
      itoken = -1
      if (lopt) then
        print *, 'TOKMAT: unmatched ', (string(i), i=1,len)
        call cexit(-1,1)
      endif
      end
C      subroutine fmain
C      character*42 s
C      integer i,j
C      i=0
C      s = 'hello there world'
C      call skipbl(s,len(s),i)
C      print *, i
C      call skp2bl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skp2bl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skipbl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skipbl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skp2bl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skipbl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skp2bl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skipbl(s,len(s),i)
C      print *, i, s(1:i), '|'
C      call skp2bl(s,len(s),i)
C      print *, i, s(1:i), '|'
C
C      end

