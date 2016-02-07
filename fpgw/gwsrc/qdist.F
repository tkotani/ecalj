      subroutine qdist(q,q1,ux,uy,uz,gam)
c  new version of 03.10.89. u=(ux,uy,uz) gives direction,
c  gam is multiplier in real space along u, volume is conserved.
c  if gam<0 then do not maintain volume, just scale along u.
      implicit real*8 (a-h,p-z)
      dimension q(3),q1(3)
      g=1.d0/gam
      u2=ux*ux+uy*uy+uz*uz
      a=(ux*q(1)+uy*q(2)+uz*q(3))/u2
      xxx=1d0/dsqrt(dabs(g))
      if (g.lt.0d0) xxx=1d0
      yyy=a*(g-xxx)
      q1(1)=xxx*q(1)+yyy*ux
      q1(2)=xxx*q(2)+yyy*uy
      q1(3)=xxx*q(3)+yyy*uz
      return
      end

      subroutine rdist(v,v1,ux,uy,uz,g)
      implicit real*8 (a-h,p-z)
      dimension v(3),v1(3)
      u2=ux*ux+uy*uy+uz*uz
      a=(ux*v(1)+uy*v(2)+uz*v(3))/u2
      xxx=1d0/dsqrt(dabs(g))
      if (g.lt.0d0) xxx=1d0
      yyy=a*(g-xxx)
      v1(1)=xxx*v(1)+yyy*ux
      v1(2)=xxx*v(2)+yyy*uy
      v1(3)=xxx*v(3)+yyy*uz
      return
      end
