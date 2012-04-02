      subroutine matm(a,b,c,n1,n2,n3)
c-- interface for matrix multiplication  c=a*b -------------------------
c assumed size-array
c
c!!! takao 2000 Sep.  If you use this routine with cxml of zgemm(GOTO),
c you might have to allocate a(n1,n2+1) to avoid the stack segmanation fault due to the complilar bug.
c
      integer(4)::  n1,n2,n3
      complex(8) :: a(n1,n2),b(n2,n3), c(n1,n3)
c
c      print *,"sumcheck=",sum(a),sum(b),n1,n2,n3
      call ZGEMM ( "N", "N", n1, n3, n2, (1d0,0d0), 
     &             a, n1, 
     &             b, n2,
     &             (0d0,0d0), c, n1 )
      end

      subroutine matcinv(n,a)
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

      subroutine matinv(n,a)
      implicit none
      integer(4) :: n, info, ipiv(n)
      real(8):: a(n,n)
      real(8),allocatable:: work(:)
      call dgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        print *,' matinv: zegtrf info=',info
        stop    ' matinv: zegtrf '
      endif
      allocate(work(n*n))
      call dgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        print *,'matinv: zegtri info=',info
        stop    'matinv: zegtri '
      endif
      end

c--------------------------------------------------------------------------
c      subroutine mml(a,b,c) This would not work why?
c      Maybe becase the size information can not be passed to the lower level.
c-- interface for matrix multiplication  c=a*b--------------------------------
c assumed shape-arrayc
c     integer(4) :: n1,n2,n3
c      complex(8) :: a(:,:),b(:,:), c(:,:)
c     n1=size(a,dim=1)
c      n2=size(a,dim=2)
c      n3=size(b,dim=2)
c     print *, n1,n2,n3
c      if(n2/=size(b,dim=1)) stop "mml: 1"
c      if(n1/=size(c,dim=1)) stop "mml: 2"
c      if(n3/=size(c,dim=2)) stop "mml: 3"
c      call ZGEMM ( "N", "N", n1, n3, n2, (1d0,0d0),
c     &             a, n1,
c     &             b, n2,
c     &             (0d0,0d0), c, n1 )
c      end
