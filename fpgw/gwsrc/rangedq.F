      subroutine rangedq(qin, qout)
      real(8):: qin(3),qout(3),tol=1d-12 !,qx(3),tolq2=1d-6 ! !,tolq2=1d-4
      qout= qin-nint(qin)
      do ix=1,3
        if(qout(ix)>0.5d0-tol) qout(ix)=-0.5d0 !this is needed to distinguish 0.5d0 and -0.5d0.
      enddo 
c$$$      real(8), intent(in) :: qin(3)
c$$$      real(8), intent(out) :: qout(3)
c$$$      real(8):: qx(3),tolq2=1d-4
c$$$      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
c$$$      qout  = qx - idint(qx)
c$$$      qout  = qout                            - tolq2
c$$$      if(qout(1) <=0d0) qout(1) = 1d-16
c$$$      if(qout(2) <=0d0) qout(2) = 1d-16
c$$$      if(qout(3) <=0d0) qout(3) = 1d-16
      end


c$$$      subroutine rangedq3(qin, qout)
c$$$      real(8):: qin(3),qout(3)
c$$$      qout= qin-nint(qin)
c$$$      do ix=1,3
c$$$        if(qout(ix)>0.5d0-tol) qout(ix)=-0.5d0
c$$$      enddo
c$$$c$$$      real(8):: qin(3),qout(3),qx(3),tolq2=1d-4
c$$$c$$$      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
c$$$c$$$      qout  = qx - idint(qx)
c$$$c$$$      qout  = qout                            - tolq2
c$$$c$$$      if(qout(1) <=0d0) qout(1) = 1d-10
c$$$c$$$      if(qout(2) <=0d0) qout(2) = 1d-10
c$$$c$$$      if(qout(3) <=0d0) qout(3) = 1d-10
c$$$      end
c
c
c$$$      subroutine rangedq_xxx(qin, qout)
c$$$! this removes integer part. qout is just the fractional part of qin.
c$$$! probably we will replace this with nint.
c$$$      real(8):: qin,qout,qx,tolq2=1d-4
c$$$      qx    = qin + dble(int(abs(qin))) + 1d0 + tolq2
c$$$      qout  = qx - idint(qx)
c$$$      qout  = qout                            - tolq2
c$$$      if(qout<=0d0) qout=1d-16
c$$$      end

