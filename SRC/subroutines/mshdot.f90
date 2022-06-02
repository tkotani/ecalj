subroutine mshdot(vol,n,n1,n2,n3,k1,k2,k3,c,d,sum1,sum2)
  use m_lgunit,only:stdo
  !- Scalar product of two functions tabulated on the real-space mesh
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   vol   :cell volume
  !i n1,n2,n3:size of mesh
  !i k1,k2,k3:dimensions arrays c,d
  !i    c    :first function
  !i    d    :second function
  !o Outputs
  !o   sum1  :real part of integral
  !o   sum2  :imaginary part of integral
  !r Remarks
  !u Updates
  !u   21 Apr 00 Adapted from nfp meshdot
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n,n1,n2,n3,k1,k2,k3
  double precision :: vol,sum1,sum2
  double complex c(k1,k2,k3,n),d(k1,k2,k3,n)
  ! ... Local parameters
  integer :: i1,i2,i3,iprint,i
  double complex csum

  !      stdo = lgunit(1)
  csum = 0d0
  do  i = 1, n
     do  i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              csum = csum + c(i1,i2,i3,i)*d(i1,i2,i3,i)
           enddo
        enddo
     enddo
  enddo

  csum = csum*vol/(n1*n2*n3)
  sum1 = dble(csum)
  sum2 = dimag(csum)

  if (iprint() >= 50) write(stdo,724) sum1,sum2
724 format(/' meshdot: integral',2f12.6)

end subroutine mshdot

