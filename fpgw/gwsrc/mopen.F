C= Taken from NFP code =============================
c --- mopen
      subroutine mopen(ifi,ext,type)
      character*120 id,fname
      character*(*) ext
      character*1 type
      common /fdefs/ jpr,nid,id
      write (fname,100) id(1:nid),ext
  100 format(a,'.',a)
      if (type.eq.'f' .or. type.eq.'F') then
        if(jpr.ge.2) write(6,300) ifi,fname
  300   format(' == open  file',i3,'   ',a16,'  (formatted)')
        open (ifi, file=fname, form='FORMATTED',status='UNKNOWN')
      else if (type.eq.'u' .or. type.eq.'U') then
        if(jpr.ge.2) write(6,301) ifi,fname
  301   format(' == open  file',i3,'   ',a16,'  (unformatted)')
        open (ifi, file=fname, form='UNFORMATTED',status='UNKNOWN')
      else
c         call rx('mopen: type is not f or u') ! or a-->takao')
Cstop2rx 2013.08.09 kino        stop 'mopen: type is not f or u'
        call rx( 'mopen: type is not f or u')
      endif
      end
c --- strip
      subroutine strip(str,i1,i2)
C!  strip(str,i1,i2) -- Returns i1,i2 so that str(i1:i2) stripped.
      character*(*) str
      i1=0
      do 1 i=1,len(str)
        if(str(i:i).ne.' ') then
          i1=i
          goto 2
        endif
  1   continue
      i1=1
      i2=0
      return
  2   i2=len(str)+1
      do 3 i=len(str),1,-1
        if(str(i:i).ne.' ') then
          i2=i
          goto 4
        endif
  3   continue
  4   continue
      end


c ------ words ------
      subroutine words(str,nw)
C  counts blank-delimted words in str
      character*(*) str
      nw=0
      i1=0
      i2=0
      i0=1
  99  do i=i0,len(str)
        if(str(i:i).ne.' ') then
          i1=i
          goto 90
        endif
  10  enddo
      return
  90  nw=nw+1
      do i=i1,len(str)
        if(str(i:i).eq.' ') then
          i2=i
          goto 91
        endif
      enddo
      return
  91  i0=i2
      goto 99
      end
c --- mfdef
      subroutine mfdef(ipr)
c  Sets file identifier from env-variable LMJOBID.
c  ipr controls verbosity for subsequent opens and closes.
      character*120 id,id1
      common /fdefs/ jpr,nid,id
      jpr=ipr
      call getenv('LMJOB',id1)
      call words(id1,nw)
C      if(nw.eq.0)
C     .   call rx('mfdef: environment variable LMJOB is not set')
      if(nw.eq.0) id1 = 'dat'
      call word(id1,1,i1,i2)
      id=id1(i1:i2)
      nid=i2-i1+1
      if (jpr.ge.1) write(6,300) id(1:nid)
  300 format(' I/O file identifier: ',a)
      end
c ------ word ------
      subroutine word(str,iw,j1,j2)
C  returns i1,i2 so that str(j1:j2) is the iw-th word
      character*(*) str
      nw=0
      i1=0
      i2=0
      i0=1
  99  continue
      do i=i0,len(str)
        if(str(i:i).ne.' ') then
          i1=i
          goto 90
        endif
      enddo
      j1=0
      j2=0
      return
  90  nw=nw+1
      if (nw.eq.iw) j1=i1
      do i=i1,len(str)
        if(str(i:i).eq.' ') then
          i2=i
          goto 91
        endif
      enddo
      i2=len(str)+1
  91  i0=i2
      if (nw.eq.iw) then
        j2=i2-1
        return
      endif
      goto 99
      end






