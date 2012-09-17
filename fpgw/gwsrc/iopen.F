c-- For digital alpha routines -----------------------
      integer function iopen (nam,iform,istat,mrecl)

c opens a file with name >>nam<<  with unit file iopen
c the integer function ifile (nam) gives  the file unit of file nam
c if the file is not opened, ifile = 0
c the integer function iclose(nam) closes the file unit of file nam
c if the file is not opened, it does not matter

c nam   = file name, up to 32 characters
c iform = 0 ==> unformatted, otherwise formatted
c istat = 0 ==> old,
c       = 1 ==> new,
c       = 2 ==> scratch, otherwise unknown
c mrecl = maximum record length for direct access (in bytes)
c         if = 0 ==> sequential access

c iopen = file number, starting from unit 10

      implicit real*8(a-h,o-z)
      parameter     (mxfil=1000) !ifi=10 to 10+mxfil-1

      character*(*) nam
      character*32  fnamn
      character*64  namt
      character*11  fmt,stat
      dimension     fnamn(mxfil),iunit(mxfil)
      save          fnamn,iunit,namn

      logical :: namnx
      data          namn /0/
      integer(4):: namnmax=0,verbose
      logical :: debug=.false.
c check that file is not already opened oct2005
      if(debug) then
        write(6,*) '------iopenstart ---------'
        do   i = 1,namnmax
          write (6,"(i3,' ',a)") i+9,fnamn(i)
        enddo
      endif

c skip leading blanks
c      namt       = nam
c      i          = 0
c      call rmvbl (namt,32,i)
c      namt       = namt(i+1:32)
       namt= trim(adjustl(nam))

c check that file is not already opened
      do       i = 1,namnmax
        if(namt .eq. fnamn(i)) then
          write (*,6000) fnamn(i)
          stop 'iopen: file already exists'
        endif
      end do

!takao Aug2005
      naminit= 1
      do ix= naminit, naminit+ mxfil
        namn=ix
        inquire(unit=namn+9, opened=namnx)
c        print *,' namnxxx=',namn,namnx,fnamn(namn)
        if(.not.namnx) goto 1012
      enddo
      stop 'iopen: enlarge mxfil...'
 1012 continue
      if(verbose()>91) write(6,"(' opening: ',i5,1x,a)") namn+9,namt
!      namn       = namn + 1
!      if(namn .gt. mxfil) stop 'iopen: too many files'

      fnamn(namn)= namt
      iunit(namn)= 9 + namn
      if(namn>namnmax) namnmax=namn
cccccccccccc
c      print *,'namn,namnmax=',namn,namnmax
ccccccccccccccc
ctakao
      if(debug) then
        do ix=9+1,namnmax
          write(6,"(' iopen: namn iunit fnamn =',i5,i5,' ',a32)")
     & ix,iunit(ix),fnamn(ix)
        enddo
      endif

c format
      fmt        = 'formatted'
      if(iform .eq. 0) fmt = 'unformatted'

c status
      stat       = 'unknown'
      if(istat .eq. 0) stat = 'old'
      if(istat .eq. 1) stat = 'new'
      if(istat .eq. 2) stat = 'scratch'

c open file
c> sequential
      iopen      = iunit(namn)
      if (mrecl .eq. 0) then
        open(iopen,file=namt,form=fmt,status=stat)
      endif

c> direct access
      if (mrecl .gt. 0) then
        print *,' mrelc.ne.0 file=',namt,mrecl
        open(iopen,file=namt,form=fmt,status=stat,access='direct',
     .     recl=mrecl)
      endif

c formats
 6000 format(a)
      return

c----------------------------------------------------------------------
c identify the file unit for the given file name >>nam<<
      entry ifile (nam)

c skip leading blanks
c      namt       = nam
c      i          = 0
c      call rmvbl (namt,32,i)
c      namt       = namt(i+1:32)
      namt=trim(adjustl(namt))

c identify file unit
      ifile      = 0
      do       i = 1,namn
        if(namt .eq. fnamn(i))then
          ifile      = iunit(i)
          return
        endif
      end do
c     if(ifile .eq. 0)stop 'ifile: cannot find the file'

c----------------------------------------------------------------------
c close file unit for the given file name >>nam<<
      entry iclose(nam)

c check that file is not already opened oct2005
      if(debug) then
        write(6,*) '------iclosestart ---------'
        print *,'namn,namnmax=',namn,namnmax
        write (6,"(' ',a)") nam
        do   i = 1,namnmax
          write (6,"(2i3,' ',a)") i+9,iunit(i),fnamn(i)
        enddo
      endif

c skip leading blanks
c      namt       = nam
c      i          = 0
c      call rmvbl (namt,32,i)
c      namt       = namt(i+1:32)
      namt = trim(adjustl(nam))

c identify file unit
      iclose     = 0
      do       i = 1,namnmax !takao fix sep2012. it was
        if(namt .eq. fnamn(i))then
ctakao
c      write(6,1033) i,iunit(i),fnamn(i)
c 1033 format(' iclose: namn fnamn iunit=',i5,i5,' ',a32)
          iclose     = iunit(i)
          fnamn(i)   = 'close_djfoaafai' !takao Aug2005
c
          close (iclose)
          if(verbose()>91) write(6,"(' closed : ',i5,1x,a)") iclose,namt
c          write(6,"(' closed : ',i5,1x,a)") iclose,namt
          return
        endif
      end do
c     if(iclose .eq. 0)stop 'iclose: cannot find the file'

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cctaken from system.f
c      subroutine rmvbl(t,nt,i)
cC- Parses string T(I) for blanks
c      integer nt,i
c      character*1 t(0:nt)
c   99 if (t(i) .ne. ' ') return
c      i = i + 1
c      if (i .ge. nt) return
c      goto 99
c      end
