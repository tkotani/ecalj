      subroutine mkilst(strn,nlist,list)
C- Resolve list (ascii string) into a vector of integers
C ----------------------------------------------------------------------
Ci Inputs
Ci   strn  :string holding list of integers
Co Outputs
Co   nlist :number of integers in list
Co         :nlist<0 => mkilst failed to parse list
Co   list  :list of integers
Cr Remarks
Cr   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
Cr   low:high:step
Cr   low, high, and step are integer expressions specifying the sequence
Cr     low, low+step, low+2*step, ... high.
Cr   If :step is missing, the step size defaults to 1.  If also :high
Cr   is missing,  the sequence reduces to a single integer. Thus,
Cr     '5+1'       becomes a single number, 6.
Cr     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
Cr     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
Cr   Sequences may be strung together separated by commas, eg
Cr     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
Cu Updates
Cu   02 Feb 01 strn is now a character string
C ----------------------------------------------------------------------
C     implicit none
      integer list(*),nlist
      character*(*) strn
      integer it(512),iv(512),a2vec,ip,i,j,k

      ip = 0
      nlist = -1
      call skipbl(strn,len(strn),ip)
      k = a2vec(strn,len(strn),ip,2,',: ',3,3,100,it,iv)
      if (k .lt. 1) return
      if (k .ge. 99) call rx('mkilst: increase size of iv')
      it(k+1) = 0
      iv(k+1) = iv(k)
C ... loop over all iv
      nlist = 0
      i = 0
   14 i = i+1
C ... Case iv => a single number
      if (it(i) .ne. 2) then
        nlist = nlist+1
        list(nlist) = iv(i)
C ... Case iv => n1:n2:n3
      elseif (it(i+1) .eq. 2) then
        do  12  j = iv(i), iv(i+1), iv(i+2)
          nlist = nlist+1
   12   list(nlist) = j
        i = i+2
C ... Case iv => n1:n2
      else
        do  17  j = iv(i), iv(i+1)
          nlist = nlist+1
   17   list(nlist) = j
        i = i+1
      endif
      if (i .lt. k) goto 14
      end

      subroutine mkilss(iopt,slst,nlist,list)
C- Resolve a list into a sorted vector of integers, paring duplicates
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt  :1s digit
Ci         : 0, do not sort list or pare duplicates
Ci         : 1, sort list, paring duplicates
Ci         :10s digit
Ci         : 0 return nlist<0 if fail to parse list
Ci         : 1 abort with error message if fail to parse list
Ci   slst  :string holding list of integers
Co Outputs
Co   nlist :number of integers in list
Co         :nlist<0 => mkilss failed to parse list
Co   list  :list of integers
Cr Remarks
Cr   Syntax: Na,Nb,... where each of the Na, Nb, etc ... has a syntax
Cr   low:high:step
Cr   low, high, and step are integer expressions specifying the sequence
Cr     low, low+step, low+2*step, ... high.
Cr   If :step is missing, the step size defaults to 1.  If also :high
Cr   is missing,  the sequence reduces to a single integer. Thus,
Cr     '5+1'       becomes a single number, 6.
Cr     '5+1:8+2'   becomes a sequence of numbers, 6 7 8 9 10
Cr     '5+1:8+2:2' becomes a sequence of numbers, 6 8 10.
Cr   Sequences may be strung together separated by commas, eg
Cr     '11,2,5+1:8+2:2' becomes a list 11 2 6 8 10.
Cu Updates
Cu   20 Oct 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer iopt,list(*),nlist
      character*(*) slst
C ... Local parameters
      integer i,j

      call mkilst(slst,nlist,list)
      if (nlist .lt. 0 .and. (iopt/10 .gt. 0))
     .call rxs('mkilss: failed to parse site list ',slst)
      if (nlist .lt. 0) return
      if (mod(iopt,10) .eq. 0) return

C     Sort list and pare duplicates
      call ishell(nlist,list)
      j = 1
      do  16  i = 2, nlist
        if (list(i) .gt. list(j)) then
          list(j+1) = list(i)
          j = j+1
        endif
   16 continue
      nlist = j
      end

      subroutine mkils0(strn,nlist,ip)
C- Like mkilst, but returns size of list and last char parsed in strn
C     implicit none
      integer nlist
      character*(*) strn
      integer it(100),iv(100),a2vec,ip,i,j,k

      ip = 0
      nlist = -1
      call skipbl(strn,len(strn),ip)
      if (ip .ge. len(strn)) return
      k = a2vec(strn,len(strn),ip,2,',: ',3,3,100,it,iv)
      if (k .lt. 1) return
      if (k .ge. 99) call rx('mkilst: increase size of iv')
      it(k+1) = 0
      iv(k+1) = iv(k)
C ... loop over all iv
      nlist = 0
      i = 0
   14 i = i+1
C ... Case iv => a single number
      if (it(i) .ne. 2) then
        nlist = nlist+1
C ... Case iv => n1:n2:n3
      elseif (it(i+1) .eq. 2) then
        do  12  j = iv(i), iv(i+1), iv(i+2)
   12   nlist = nlist+1
        i = i+2
C ... Case iv => n1:n2
      else
        do  17  j = iv(i), iv(i+1)
   17   nlist = nlist+1
        i = i+1
      endif
      if (i .lt. k) goto 14

      end
      subroutine ilst2a(list,nlist,strn)
C- Make a an ascii list from a vector of integers
C ----------------------------------------------------------------------
Ci Inputs
Co   list  :list of integers
Co   nlist :number of integers in list
Co Outputs
Ci   strn  :string holding list of integers
Cr Remarks
Cr   ilst2a performs the inverse operation to mkilst,
Cr   making a compact ascii representation of a list
Cu Updates
Cu   13 Sep 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlist,list(nlist)
      character*(*) strn
C ... Local variables
      integer i,i1,ls,awrite,ip,j
C     integer iilst(9)
C     data iilst /1,3,5,2,4,6,8,10,99/

      strn = ' '
      ip = 1
      if (nlist .eq. 0) return
      ls = len(strn)
      i1 = 1
      i = 0
   10 continue
      i = i+1
      if (i .gt. nlist) then
        return
      endif
      if (ip .gt. 1) then
        strn(ip:ip) = ','
        ip = ip+1
      endif

C   ... Near the end of list --- no sequences to check
      if (i .gt. nlist-2) then
        ip = ip+awrite('%i',strn(ip:),ls,0,list(i),0,0,0,0,0,0,0)
        goto 10
      endif

C   ... Case linear sequence of at least three
      if (list(i+1)-list(i) .eq. list(i+2)-list(i+1)) then
        i1 = i
        j = list(i+1)-list(i)
   20   continue
        if (i .le. nlist-2) then
          i = i+1
          if (j .eq. list(i+2)-list(i+1)) goto 20
        endif
        if (j .eq. 1) then
          ip = ip+awrite('%i:%i',strn(ip:),ls,0,
     .    list(i1),list(i+1),j,0,0,0,0,0)
        else
          ip = ip+awrite('%i:%i:%i',strn(ip:),ls,0,
     .    list(i1),list(i+1),j,0,0,0,0,0)
        endif
        i = i+1
        goto 10
      endif

C   ... Case no linear sequence ... just append single element
      ip = ip+awrite('%i',strn(ip:),ls,0,list(i),0,0,0,0,0,0,0)
      goto 10
      end

#if TEST
      subroutine fmain
      implicit none
      character*20 strn
      integer nlist,list(20),i

      strn = '                 2,1'
      call mkilst(strn,nlist,list)
      print *, nlist, (list(i), i=1,nlist)
      strn = '                2,1 '
      call mkilst(strn,nlist,list)
      print *, nlist, (list(i), i=1,nlist)
      strn = '             22:33:3'
      call mkilst(strn,nlist,list)
      print *, nlist, (list(i), i=1,nlist)

      end
#endif

