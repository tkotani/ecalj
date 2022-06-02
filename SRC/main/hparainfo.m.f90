      program parainfo
c$$$C - Split QPNT into QPNT and get kinit and kend
c$$$      use keyvalue
c$$$      implicit real*8(a-h,o-z)
c$$$      character*8 charext
c$$$      character*100 qhead(8)
c$$$      integer(4),allocatable:: kinit(:),kend(:)
c$$$      real(8),allocatable:: q(:,:)
c$$$
c$$$      integer(4):: ret
c$$$c      logical:: readgwinput
c$$$
c$$$      write(6,*) ' number of machine? '
c$$$      read(5,*)   nmachine0
c$$$      call headver('hparainfo',machine0)
c$$$      ifix=ifile_hanlde()
c$$$
c$$$c --- goto X0KDIV section ---
c$$$      write(6,*) ' === Generate X0KDIV section === '
c$$$      open(ifix,file='NQIBZ')
c$$$      read(101,*) nqibz,nq0i
c$$$      close(ifix)
c$$$
c$$$      nqtot = nqibz +nq0i-1
c$$$      nmachine = nmachine0
c$$$      if(nqtot < nmachine ) nmachine = nqtot !correct number of nmachine when nmachine> number of k
c$$$      kadd  = nqtot/ nmachine
c$$$      nnn   = nqtot - kadd *nmachine
c$$$
c$$$      allocate(kinit(nmachine+1),kend(nmachine+1))
c$$$      kinit (1) = 2
c$$$      do i=1,nmachine
c$$$        kend(i) = kinit(i) + kadd -1
c$$$        if( i<=nnn ) kend(i) = kend(i) + 1
c$$$        kinit(i+1) = kend(i) +1
c$$$      enddo
c$$$      open(ifix,file='X0KDIV')
c$$$      write(ifix,*) nmachine
c$$$      write(ifix,"(100i5)") (kinit(i),i=1,nmachine)
c$$$      write(ifix,"(100i5)") (kend (i),i=1,nmachine)
c$$$      close(ifix)
c$$$      deallocate(kinit,kend)
c$$$
c$$$
c$$$c --- goto QPNT section ---
c$$$      write(6,*) ' === Generate QPNT.{number} section === '
c$$$c      ifqpnt    = iopen('QPNT',1,0,0)
c$$$c      if(readgwinput()) then
c$$$      call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
c$$$c      else
c$$$c        ifqpnt    = iopen('QPNT',1,0,0)
c$$$c      endif
c$$$
c$$$      do i=1,8
c$$$        read (ifqpnt,"(a)") qhead(i)
c$$$c        write (6,"(a)") qhead(i)
c$$$      enddo
c$$$      rewind(ifqpnt)
c$$$      call readx (ifqpnt,1000 )
c$$$      call readx (ifqpnt,100)
c$$$      call readx (ifqpnt,100)
c$$$      read (ifqpnt,*) nqtot
c$$$      allocate(q(3,nqtot))
c$$$      do       k = 1, nqtot
c$$$        read (ifqpnt,      *)  i,q(1,k),q(2,k),q(3,k)
c$$$c        write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
c$$$      enddo
c$$$c
c$$$      nmachine = nmachine0
c$$$      if(nqtot < nmachine ) nmachine = nqtot !correct number of nmachine when nmachine> number of k
c$$$      kadd  = nqtot/ nmachine
c$$$      nnn   = nqtot - kadd *nmachine
c$$$      allocate(kinit(nmachine+1),kend(nmachine+1))
c$$$      kinit (1) = 1
c$$$      do i = 1,nmachine
c$$$        open (ifix,file ='QPNT.'//charext(i))
c$$$        kend(i) = kinit(i) + kadd -1
c$$$        if( i<=nnn ) kend(i) = kend(i) + 1
c$$$        do ix=1,8
c$$$          write (ifix,"(a)") qhead(ix)
c$$$        enddo
c$$$        write (ifix,*) kend(i) - kinit(i) + 1
c$$$        do ix = kinit(i), kend(i)
c$$$          write(ifix,'(i3,3f23.16)') ix,q(1:3,ix)
c$$$        enddo
c$$$        close(ifix)
c$$$        kinit(i+1) = kend(i) +1
c$$$      enddo
      call rx0( ' OK! parainfo QPNT.{number} and X0KDIV generated')
      end
