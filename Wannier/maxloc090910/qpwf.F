      program hqpwf
c-------------------------------------------------------------
c output 'Q0Pwf' file
c
c------------------------------------------------------------
      use keyvalue
      implicit real*8(a-h,o-z)

      character(len=150):: recrdxxx
      character(len=10) :: keyw1='unit_2pioa',keyw2
      logical ::unit2=.false. !  readgwinput,
      integer(4)::nnn(3),ret

      real*8 :: q0wf(3),dq0wf(3),qpwf(3)
c------------------------------------


c QPNT data
      call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret) 
      print *,' ifqpnt ret=',ifqpnt,ret
c
      call readx   (ifqpnt,10)
      read (ifqpnt,*) iqall,iaf
      call readx   (ifqpnt,100)
      read (ifqpnt,*) nwf,q0wf(1:3),dq0wf(1:3)
      qpwf = q0wf + dq0wf

      close(ifqpnt)


c      ifi = 6661
c      open (ifi, file='Q0Pwf')
c      write(ifi,"(i5)") 1
c      write(ifi,"(d24.16,3x, 3d24.16)" )
c     &     (1.0d0,dq0wf(1:3))
c      close(ifi)

      ifi = 6661
      open (ifi, file='QPwf')
      write(ifi,"(i5)") 1
      write(ifi,"(d24.16,3x, 6d24.16)" )
     &     1.0d0,qpwf(1:3),dq0wf(1:3)
      close(ifi)

      call cputid(0)
      stop ' OK! q0pwf'
      end

c-----------------------------------------------------------------------
