      program parainfo
C - Split QPNT into QPNT and get kinit and kend
      use keyvalue
      implicit real*8(a-h,o-z)
      character*8 charext
      character*100 qhead(8)
      integer(4),allocatable:: kinit(:),kend(:)
      real(8),allocatable:: q(:,:)

      integer(4):: ret
c      logical:: readgwinput

      write(6,*) ' number of machine? '
      read(5,*)   nmachine0
      call headver('hparainfo',machine0)
      ifix=101

c --- goto X0KDIV section ---
      write(6,*) ' === Generate X0KDIV section === '
      open(ifix,file='NQIBZ')
      read(101,*) nqibz,nq0i
      close(ifix)

      nqtot = nqibz +nq0i-1
      nmachine = nmachine0
      if(nqtot < nmachine ) nmachine = nqtot !correct number of nmachine when nmachine> number of k
      kadd  = nqtot/ nmachine
      nnn   = nqtot - kadd *nmachine

      allocate(kinit(nmachine+1),kend(nmachine+1))
      kinit (1) = 2
      do i=1,nmachine
        kend(i) = kinit(i) + kadd -1
        if( i<=nnn ) kend(i) = kend(i) + 1
        kinit(i+1) = kend(i) +1
      enddo
      open(ifix,file='X0KDIV')
      write(ifix,*) nmachine
      write(ifix,"(100i5)") (kinit(i),i=1,nmachine)
      write(ifix,"(100i5)") (kend (i),i=1,nmachine)
      close(ifix)
      deallocate(kinit,kend)


c --- goto QPNT section ---
      write(6,*) ' === Generate QPNT.{number} section === '
c      ifqpnt    = iopen('QPNT',1,0,0)
c      if(readgwinput()) then
      call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
c      else
c        ifqpnt    = iopen('QPNT',1,0,0)
c      endif

      do i=1,8
        read (ifqpnt,"(a)") qhead(i)
c        write (6,"(a)") qhead(i)
      enddo
      rewind(ifqpnt)
      call readx (ifqpnt,1000 )
      call readx (ifqpnt,100)
      call readx (ifqpnt,100)
      read (ifqpnt,*) nqtot
      allocate(q(3,nqtot))
      do       k = 1, nqtot
        read (ifqpnt,      *)  i,q(1,k),q(2,k),q(3,k)
c        write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
      enddo
c
      nmachine = nmachine0
      if(nqtot < nmachine ) nmachine = nqtot !correct number of nmachine when nmachine> number of k
      kadd  = nqtot/ nmachine
      nnn   = nqtot - kadd *nmachine
      allocate(kinit(nmachine+1),kend(nmachine+1))
      kinit (1) = 1
      do i = 1,nmachine
        open (ifix,file ='QPNT.'//charext(i))
        kend(i) = kinit(i) + kadd -1
        if( i<=nnn ) kend(i) = kend(i) + 1
        do ix=1,8
          write (ifix,"(a)") qhead(ix)
        enddo
        write (ifix,*) kend(i) - kinit(i) + 1
        do ix = kinit(i), kend(i)
          write(ifix,'(i3,3f23.16)') ix,q(1:3,ix)
        enddo
        close(ifix)
        kinit(i+1) = kend(i) +1
      enddo
Cstop2rx 2013.08.09 kino      stop ' OK! parainfo QPNT.{number} and X0KDIV generated'
      call rx0( ' OK! parainfo QPNT.{number} and X0KDIV generated')
      end
