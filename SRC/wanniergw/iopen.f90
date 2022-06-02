!-- For digital alpha routines -----------------------
integer function iopen (nam,iform,istat,mrecl)

  ! opens a file with name >>nam<<  with unit file iopen
  ! the integer function ifile (nam) gives  the file unit of file nam
  ! if the file is not opened, ifile = 0
  ! the integer function iclose(nam) closes the file unit of file nam
  ! if the file is not opened, it does not matter

  ! nam   = file name, up to 32 characters
  ! iform = 0 ==> unformatted, otherwise formatted
  ! istat = 0 ==> old,
  !       = 1 ==> new,
  !       = 2 ==> scratch, otherwise unknown
  ! mrecl = maximum record length for direct access (in bytes)
  !         if = 0 ==> sequential access

  ! iopen = file number, starting from unit 10

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  parameter     (mxfil=1000) !ifi=10 to 10+mxfil-1

  character*(*) nam
  character(32) ::  fnamn
  character(64) ::  namt
  character(11) ::  fmt,stat
  dimension     fnamn(mxfil),iunit(mxfil)
  save          fnamn,iunit,namn

  logical :: namnx
  data          namn /0/
  integer(4):: namnmax=0,verbose
  logical :: debug=.false.
  !$$$c check that file is not already opened oct2005
  !$$$      if(debug) then
  !$$$        write(6,*) '------iopenstart ---------'
  !$$$        do   i = 1,namnmax
  !$$$          write (6,"(i3,' ',a)") i+9,fnamn(i)
  !$$$        enddo
  !$$$      endif

  ! skip leading blanks
  !      namt       = nam
  !      i          = 0
  !      call rmvbl (namt,32,i)
  !      namt       = namt(i+1:32)
  namt= trim(adjustl(nam))

  ! check that file is not already opened
  do       i = 1,namnmax
     if(namt == fnamn(i)) then
        write (*,6000) fnamn(i)
        ! top2rx 2013.08.09 kino          stop 'iopen: file already exists'
        call rx( 'iopen: file already exists')
     endif
  end do

  ! akao Aug2005
  naminit= 1
  do ix= naminit, naminit+ mxfil
     namn=ix
     inquire(unit=namn+9, opened=namnx)
     !        write(6,*)' namnxxx=',namn,namnx,fnamn(namn)
     if( .NOT. namnx) goto 1012
  enddo
  ! top2rx 2013.08.09 kino      stop 'iopen: enlarge mxfil...'
  call rx( 'iopen: enlarge mxfil...')
1012 continue
  if(verbose()>91) write(6,"(' opening: ',i5,1x,a)") namn+9,namt
  !      namn       = namn + 1
  !      if(namn .gt. mxfil) stop 'iopen: too many files'

  fnamn(namn)= namt
  iunit(namn)= 9 + namn
  if(namn>namnmax) namnmax=namn
  ! cccccccccc
  !      write(6,*)'namn,namnmax=',namn,namnmax
  ! ccccccccccccc
  ! akao
  if(debug) then
     do ix=9+1,namnmax
        write(6,"(' iopen: namn iunit fnamn =',i5,i5,' ',a32)") &
             ix,iunit(ix),fnamn(ix)
     enddo
  endif

  ! format
  fmt        = 'formatted'
  if(iform == 0) fmt = 'unformatted'

  ! status
  stat       = 'unknown'
  if(istat == 0) stat = 'old'
  if(istat == 1) stat = 'new'
  if(istat == 2) stat = 'scratch'

  ! open file
  !> sequential
  iopen      = iunit(namn)
  if (mrecl == 0) then
     open(iopen,file=namt,form=fmt,status=stat)
  endif

  !> direct access
  if (mrecl > 0) then
     write(6,*)' mrelc.ne.0 file=',namt,mrecl
     open(iopen,file=namt,form=fmt,status=stat,access='direct', &
          recl=mrecl)
  endif

  ! formats
6000 format(a)
  return

  !$$$c----------------------------------------------------------------------
  !$$$c identify the file unit for the given file name >>nam<<
  !$$$      entry ifile (namt) !comment out Sep201(related to 'KPNT' and 'PRODUCT').
  !$$$
  !$$$c skip leading blanks
  !$$$c      namt       = nam
  !$$$c      i          = 0
  !$$$c      call rmvbl (namt,32,i)
  !$$$c      namt       = namt(i+1:32)
  !$$$      namt=trim(adjustl(namt))
  !$$$
  !$$$c identify file unit
  !$$$      ifile      = 0
  !$$$      do       i = 1,namn
  !$$$        if(namt .eq. fnamn(i))then
  !$$$          ifile      = iunit(i)
  !$$$          return
  !$$$        endif
  !$$$      end do
  !$$$c     if(ifile .eq. 0)stop 'ifile: cannot find the file'

  !----------------------------------------------------------------------
  ! close file unit for the given file name >>nam<<
  entry iclose(nam)

  ! check that file is not already opened oct2005
  if(debug) then
     write(6,*) '------iclosestart ---------'
     write(6,*)'namn,namnmax=',namn,namnmax
     write (6,"(' ',a)") nam
     do   i = 1,namnmax
        write (6,"(2i3,' ',a)") i+9,iunit(i),fnamn(i)
     enddo
  endif

  ! skip leading blanks
  !      namt       = nam
  !      i          = 0
  !      call rmvbl (namt,32,i)
  !      namt       = namt(i+1:32)
  namt = trim(adjustl(nam))

  ! identify file unit
  iclose     = 0
  do       i = 1,namnmax !takao fix sep2012. it was
     if(namt == fnamn(i))then
        ! akao
        !      write(6,1033) i,iunit(i),fnamn(i)
        ! 1033 format(' iclose: namn fnamn iunit=',i5,i5,' ',a32)
        iclose     = iunit(i)
        fnamn(i)   = 'close_djfoaafai' !takao Aug2005

        close (iclose)
        if(verbose()>91) write(6,"(' closed : ',i5,1x,a)") iclose,namt
        !          write(6,"(' closed : ',i5,1x,a)") iclose,namt
        return
     endif
  end do
  !     if(iclose .eq. 0)stop 'iclose: cannot find the file'

end function iopen

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! taken from system.f
!      subroutine rmvbl(t,nt,i)
!C- Parses string T(I) for blanks
!      integer nt,i
!      character*1 t(0:nt)
!   99 if (t(i) .ne. ' ') return
!      i = i + 1
!      if (i .ge. nt) return
!      goto 99
!      end
