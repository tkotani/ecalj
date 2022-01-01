      SUBROUTINE GRPOP(V,V1,G,I)
C     implicit none
      double precision G(3,3,*),V(3),V1(3)
      integer i
      V1(1) = G(1,1,I)*V(1) + G(1,2,I)*V(2) + G(1,3,I)*V(3)
      V1(2) = G(2,1,I)*V(1) + G(2,2,I)*V(2) + G(2,3,I)*V(3)
      V1(3) = G(3,1,I)*V(1) + G(3,2,I)*V(2) + G(3,3,I)*V(3)
      RETURN
      END
