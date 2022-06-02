! takao modified this based on delw_2.2.
subroutine dsumdf(n,scal,a1,ofa1,l1,a2,ofa2,l2)
  !- Returns scaled sum and difference of two vectors
  ! ----------------------------------------------------------------
  !i Inputs
  !i   n    :number elements to scale and combine
  !i   scal :scale sum and difference by scal; see Outputs
  !i   a1   :first vector
  !i   ofa1 :offset to first entry in a1
  !i   l1   :skip length in a1
  !i   a2   :second vector
  !i   ofa2 :offset to first entry in a2
  !i   l2   :skip length in a2
  !o Outputs
  !o   a1   :a1 <- scal*(a1+a2)
  !o   a2   :a1 <- scal*(a1-a2)
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n,l1,l2,ofa1,ofa2
  double precision :: scal, a1(1), a2(1)
  ! Local parameters
  real(8) ,allocatable :: a_rv(:)
  ! --- a1-a2-> temp;  a1+a2 -> a1;  temp -> a2 ---
  allocate(a_rv(n))
  call dcopy ( n , a1 ( 1 + ofa1 ) , l1 , a_rv , 1 )
  call daxpy ( n , - 1d0 , a2 ( 1 + ofa2 ) , l2 , a_rv , 1 )
  call daxpy (n,1d0,a2(1+ofa2),l2,a1(1+ofa1),l1)
  call dcopy ( n , a_rv , 1 , a2 ( 1 + ofa2 ) , l2 )
  deallocate(a_rv)
  if (scal == 1) return
  call dscal(n,scal,a1(1+ofa1),l1)
  call dscal(n,scal,a2(1+ofa2),l1)
end subroutine dsumdf

