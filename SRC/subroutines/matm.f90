!> interface for matrix multiplication  c=a*b -------------------------
      subroutine matm(a,b,c,n1,n2,n3)
      integer(4), intent(in) :: n1,n2,n3
      complex(8), intent(in) :: a(n1,n2), b(n2,n3)
      complex(8), intent(out) :: c(n1,n3)
      c=matmul(a,b)
      end
      
      subroutine matmaw(a,b,c,n1,n2,n3,ww)
      integer(4), intent(in) :: n1,n2,n3
      complex(8), intent(in) :: a(n1,n2), b(n2,n3)
      real(8), intent(in) :: ww(n1)
      complex(8), intent(out) :: c(n1,n3)
      complex(8):: aa(n1,n2)
      integer:: ix
      do ix=1,n2
         aa(:,ix) = ww(:)*a(:,ix)
      enddo   
      c= c+ matmul(aa,b)
      end
      subroutine matma(a,b,c,n1,n2,n3)
      integer(4), intent(in) :: n1,n2,n3
      complex(8), intent(in) :: a(n1,n2), b(n2,n3)
      complex(8), intent(out) :: c(n1,n3)
      c= c+ matmul(a,b)
      end

      subroutine matcinv(n,a)
      implicit none
      integer(4) :: n, info, ipiv(n)
      complex(8):: a(n,n)
      complex(8),allocatable:: work(:)
      call zgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matcinv: zegtrf info=',info
        call rx( ' matcinv: zegtrf ')
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matcinv: zegtri info=',info
        call rx( 'matcinv: zegtri ')
      endif
      end

      subroutine matinv(n,a)
      implicit none
      integer(4) :: n, info, ipiv(n)
      real(8):: a(n,n)
      real(8),allocatable:: work(:)
      call dgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matinv: degtrf info=',info
        call rx( ' matinv: degtrf ')
      endif
      allocate(work(n*n))
      call dgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matinv: degtri info=',info
        call rx( 'matinv: degtri ')
      endif
      end

      subroutine matinv2(n,a,info)
      implicit none
      integer :: n, info, ipiv(n)
      real(8):: a(n,n)
      real(8),allocatable:: work(:)
      call dgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matinv: degtrf info=',info
        return
      endif
      allocate(work(n*n))
      call dgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matinv: degtri info=',info
        return
      endif
      end

