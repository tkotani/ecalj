      subroutine matinv(n,a)
      implicit none
      integer(4) :: n, info, ipiv(n)
      real(8):: a(n,n)
      real(8),allocatable:: work(:)
      call dgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matinv: degtrf info=',info
        stop    ' matinv: degtrf '
      endif
      allocate(work(n*n))
      call dgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matinv: degtri info=',info
        stop    'matinv: degtri '
      endif
      end

      subroutine matinv2(n,a,info)
      implicit none
      integer(4) :: n, info, ipiv(n)
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


      subroutine matcinv(n,a)
C --- a inverse is returned.
      implicit none
      integer(4) :: n, info, ipiv(n)
      complex(8):: a(n,n)
      complex(8),allocatable:: work(:)
      call zgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matcinv: zegtrf info=',info
        stop    ' matcinv: zegtrf '
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matcinv: zegtri info=',info
        stop    'matcinv: zegtri '
      endif
      end
