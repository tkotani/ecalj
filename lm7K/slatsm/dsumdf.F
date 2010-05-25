c takao modified this based on delw_2.2.
      subroutine dsumdf(n,scal,a1,ofa1,l1,a2,ofa2,l2)
C- Returns scaled sum and difference of two vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   n    :number elements to scale and combine
Ci   scal :scale sum and difference by scal; see Outputs
Ci   a1   :first vector
Ci   ofa1 :offset to first entry in a1
Ci   l1   :skip length in a1
Ci   a2   :second vector
Ci   ofa2 :offset to first entry in a2
Ci   l2   :skip length in a2
Co Outputs
Co   a1   :a1 <- scal*(a1+a2)
Co   a2   :a1 <- scal*(a1-a2)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n,l1,l2,ofa1,ofa2
      double precision scal, a1(1), a2(1)
C Local parameters
      real(8) ,allocatable :: rv_w_oa(:)
C --- a1-a2-> temp;  a1+a2 -> a1;  temp -> a2 ---
      allocate(rv_w_oa(n))
      call dcopy ( n , a1 ( 1 + ofa1 ) , l1 , rv_w_oa , 1 ) 
      call daxpy ( n , - 1d0 , a2 ( 1 + ofa2 ) , l2 , rv_w_oa , 1 ) 
      call daxpy (n,1d0,a2(1+ofa2),l2,a1(1+ofa1),l1)
      call dcopy ( n , rv_w_oa , 1 , a2 ( 1 + ofa2 ) , l2 ) 
      deallocate(rv_w_oa)
      if (scal .eq. 1) return
      call dscal(n,scal,a1(1+ofa1),l1)
      call dscal(n,scal,a2(1+ofa2),l1)
      end

