SUBROUTINE GTBVEC(K,B,SHIFT,V)
  !     implicit none
  integer :: K(3)
  double precision :: V(3),B(3,3),SHIFT(3)
  V(1) = SHIFT(1) + K(1)*B(1,1) + K(2)*B(1,2) + K(3)*B(1,3)
  V(2) = SHIFT(2) + K(1)*B(2,1) + K(2)*B(2,2) + K(3)*B(2,3)
  V(3) = SHIFT(3) + K(1)*B(3,1) + K(2)*B(3,2) + K(3)*B(3,3)
END SUBROUTINE GTBVEC

