      subroutine rangedq(qin, qout)
      real(8):: qin(3),qout(3),qx(3),tolq2=1d-4
      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
      qout  = qx - idint(qx)
      qout  = qout                            - tolq2
      if(qout(1) <=0d0) qout(1) = 1d-16
      if(qout(2) <=0d0) qout(2) = 1d-16
      if(qout(3) <=0d0) qout(3) = 1d-16
c      print *,'zzzzzzzzzzzzzzzzzzzzzzzzzz'
c      print *, qx
c      print *, qout
      end

      subroutine rangedq_xxx(qin, qout)
      real(8):: qin,qout,qx,tolq2=1d-4
      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
      qout  = qx - idint(qx)
      qout  = qout                            - tolq2
      if(qout<=0d0) qout=1d-16
c      print *,'zzzzzzzzzzzzzzzzzzzzzzzzzz'
c      print *,'qx=  ',qx
c      print *,'qout=',qout
      end

      subroutine rangedq3(qin, qout)
      real(8):: qin(3),qout(3),qx(3),tolq2=1d-4
      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
      qout  = qx - idint(qx)
      qout  = qout                            - tolq2
      if(qout(1) <=0d0) qout(1) = 1d-10
      if(qout(2) <=0d0) qout(2) = 1d-10
      if(qout(3) <=0d0) qout(3) = 1d-10
      end

