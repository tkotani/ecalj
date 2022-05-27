      subroutine freq01 (nx,ua,   freqx,freqw,wx,expa)
c generates a gaussian point x between (0,1) and w = (1-x)/x
c and the weights in x
c also generates expa = exp(-ua^2 w^2)
c nx    = no. gaussian points
c  ua   = s.o.
c freqx = gaussian points
c freqw = (1-x)/x
c wx    = weights of gaussian points x
c expa  = s.o.
c originally 92.02.27 Ferdi.Aryasetiawan
      implicit real*8 (a-h,o-z)
      integer:: nx,ix
      dimension freqx(nx),freqw(nx),wx(nx),expa(nx)
c generate gaussian points
      call gauss   (nx,0.d0,1.d0,freqx,wx)
c calculate w = 1/(1+x)
      ua2        = ua*ua
      do      ix = 1,nx
        freqw(ix)  = (1.d0 - freqx(ix)) / freqx(ix)
        expa(ix)   = dexp(-ua2*freqw(ix)*freqw(ix))
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine freq01x (nx,    freqx,freqw,wx) !,expa)
c 92.02.27
c generates a gaussian point x between (0,1) and w = (1-x)/x
c and the weights in x
c also generates expa = exp(-ua^2 w^2)
c nx    = no. gaussian points
c  ua   = s.o.
c freqx = gaussian points
c freqw = (1-x)/x
c wx    = weights of gaussian points x
c expa  = s.o.
      implicit real*8 (a-h,o-z)
      dimension freqx(nx),freqw(nx),wx(nx),expa(nx)
      integer:: nx,ix
c generate gaussian points
      call gauss   (nx,0.d0,1.d0,freqx,wx)
c calculate w = 1/(1+x)
c      ua2        = ua*ua
      write(6,"(' --- freq01x:  ix    x    freqw(a.u.)---')")
      do      ix = 1,nx
        freqw(ix)  = (1.d0 - freqx(ix)) / freqx(ix)
        write(6,"('            ',i4,2f9.4)") ix,freqx(ix),freqw(ix)
      end do
      return
      end
