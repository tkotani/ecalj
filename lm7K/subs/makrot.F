      subroutine makrot(rot,r)
C- Make rotation vector and angle into deformation gradient tensor
C     implicit none
      double precision rot(4),r(3,3)
      double precision n(3),cost,sint,rnorm,dcos,dsin,dsqrt
      integer i,j,k,e,mod
      cost = dcos(rot(4))
      sint = dsin(rot(4))
      rnorm = dsqrt(rot(1)**2 + rot(2)**2 + rot(3)**2)
      do  1  k = 1, 3
        n(k) = rot(k) / rnorm
    1 continue
      do  2  i = 1, 3
        do  2  j = 1, 3
          r(i,j) = (1 - cost)*n(i)*n(j)
          if (i .eq. j) then
            r(i,j) = r(i,j) + cost
          else
            k = 3 - mod(i+j,3)
            if (i .lt. k) then
              e = 1
              if (k-i .eq. 2) e = -1
            else
              e = -1
              if (i-k .eq. 2) e = 1
            endif
            r(i,j) = r(i,j) + sint*n(k)*e
          endif
    2 continue
      end



