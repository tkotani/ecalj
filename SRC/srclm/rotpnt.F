      subroutine rotpnt(p,q,g)
C     implicit none
      double precision p(3),q(3),g(3,3),h(3)
      integer i,j
      do 1 i=1,3
        h(i)=0.d0
      do 1 j=1,3
        h(i)=h(i)+g(i,j)*p(j)
    1 continue
      do 2 i=1,3
        q(i)=h(i)
    2 continue
      end

