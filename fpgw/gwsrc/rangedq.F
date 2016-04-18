      subroutine rangedq(qin, qout)
      real(8), intent(in) :: qin(3)
      real(8), intent(out) :: qout(3)
      real(8):: qx(3),tolq2=1d-4
      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
      qout  = qx - idint(qx)
      qout  = qout                            - tolq2
      if(qout(1) <=0d0) qout(1) = 1d-16
      if(qout(2) <=0d0) qout(2) = 1d-16
      if(qout(3) <=0d0) qout(3) = 1d-16
      end

      subroutine rangedq_xxx(qin, qout)
! this removes integer part. qout is just the fractional part of qin.
! probably we will replace this with nint.
      real(8):: qin,qout,qx,tolq2=1d-4
      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
      qout  = qx - idint(qx)
      qout  = qout                            - tolq2
      if(qout<=0d0) qout=1d-16
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

