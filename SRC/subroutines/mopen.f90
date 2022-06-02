!= Taken from NFP code =============================
! --- mopen
subroutine mopen(ifi,ext,type)
  integer:: ifi,jpr,nid
  character(120) :: id,fname
  character*(*) ext
  character(1) :: type
  common /fdefs/ jpr,nid,id
  write (fname,100) id(1:nid),ext
100 format(a,'.',a)
  if (type == 'f' .OR. type == 'F') then
     if(jpr >= 2) write(6,300) ifi,fname
300  format(' == open  file',i3,'   ',a16,'  (formatted)')
     open (ifi, file=fname, form='FORMATTED',status='UNKNOWN')
  else if (type == 'u' .OR. type == 'U') then
     if(jpr >= 2) write(6,301) ifi,fname
301  format(' == open  file',i3,'   ',a16,'  (unformatted)')
     open (ifi, file=fname, form='UNFORMATTED',status='UNKNOWN')
  else
     !         call rx('mopen: type is not f or u') ! or a-->takao')
     ! top2rx 2013.08.09 kino        stop 'mopen: type is not f or u'
     call rx( 'mopen: type is not f or u')
  endif
end subroutine mopen

! --- mfdef
subroutine mfdef(ipr)
  !  Sets file identifier from env-variable LMJOBID.
  !  ipr controls verbosity for subsequent opens and closes.
  character(120) :: id,id1
  common /fdefs/ jpr,nid,id
  integer:: ipr,jpr,nw,nid,i1,i2
  jpr=ipr
  call getenv('LMJOB',id1)
  call words(id1,nw)
  !      if(nw.eq.0)
  !     .   call rx('mfdef: environment variable LMJOB is not set')
  if(nw == 0) id1 = 'dat'
  call word(id1,1,i1,i2)
  id=id1(i1:i2)
  nid=i2-i1+1
  if (jpr >= 1) write(6,300) id(1:nid)
300 format(' I/O file identifier: ',a)
end subroutine mfdef



