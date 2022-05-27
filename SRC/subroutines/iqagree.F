      real(8) function polinta2(xx, x,y,n, xmin, xmax, d0min, d0max)
c Interpolation considering the zero derivative at xmin (and/or) at xmax.
c  d0min=T means that the derivative of the interpolated function at xmin is zero.
c  d0max=T means that the derivative of the interpolated function at xmax is zero.
      implicit none
      integer(4):: n ,imin,imax, id1,id2, nnn,ni,nx,i
      real(8) :: xx, x(n), y(n),xin(3*n),yin(3*n),xmin,xmax,
     &           polinta,eps=1d-13
      logical :: d0min, d0max
      if(d0min) then
        if(abs(x(1)-xmin)<eps) then
          ni=n-1
        else
          ni=n
        endif
        do i=1,ni
          xin(i) = xmin -(x(n+1-i)-xmin)
          yin(i) = y(n+1-i)
        enddo
      else
        ni = 0
      endif
      xin(ni+1:ni+n)= x(1:n)
      yin(ni+1:ni+n)= y(1:n)
      if(d0max) then
        if(abs(x(n)-xmax)<eps) then
          nx=n-1
        else
          nx=n
        endif
        do i=1,nx
          xin(i+ni+n) = xmax + (xmax - x(1+nx-i))
          yin(i+ni+n) = y(1+nx-i)
        enddo
      else
        nx = 0
      endif
      nnn = ni+ n + nx
      polinta2 = polinta(xx,xin,yin,nnn)
      end
