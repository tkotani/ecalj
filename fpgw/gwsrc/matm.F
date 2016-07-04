      subroutine matmmsparse(a,b,c,zout,n,etol,iele)
      implicit none
      integer::  n,nmax,i1x,i2x,iele,i,iy
      complex(8) :: a(n,n),b(n,n), c(n,n),zout(n,n)
      complex(8),allocatable:: zz(:),aa(:,:),cc(:,:)
      real(8):: etol
      integer,allocatable:: i1(:),i2(:)
      nmax=n**2 !this is safe but loss of memory.
      allocate(i1(nmax),i2(nmax),zz(nmax))
      iele=0
      do i1x=1,n
      do i2x=1,n
        if(abs(b(i1x,i2x))>etol) then
         iele=iele+1
         if(iele>nmax) call rx( ' error in matmsparce: enlarge nmax')
         i1(iele)=i1x
         i2(iele)=i2x
         zz(iele)=b(i1x,i2x)
        endif
      enddo
      enddo
      allocate(aa(n,iele),cc(n,iele))
      do i=1,iele
      aa(:,i)=a(:,i1(i))
      cc(:,i)=c(i2(i),:)
      enddo
c      print *,'matmmsparse n iele=',n,iele
      zout=0d0
      do i=1,iele
      do iy=1,n
cc        zout(:,iy)= zout(:,iy) + a(:,i1(i))*zz(i)*c(i2(i),iy)
        zout(:,iy)= zout(:,iy) + aa(:,i)*zz(i)*cc(iy,i)
c        call zaxpy(n,zz(i)*cc(iy,i),aa(1,i),1,zout(1,iy),1) 
      enddo
      enddo
      deallocate(i1,i2,zz)
      end


      subroutine matmsparse(a,b,c,n,isparse,etol)
! not test yet.
      integer::  n,nmax,i1x,i2x,iele,i
      complex(8) :: a(n,n),b(n,n), c(n,n)
      complex(8),allocatable:: zz(:)
      integer,allocatable:: i1(:),i2(:)
      nmax=10*n !this is a test case.
      allocate(i1(nmax),i2(nmax),zz(nmax))
      if(isparse==1) then
         iele=0
         do i1x=1,n
         do i2x=1,n
            if(abs(a(i1x,i2x))>etol) then
               iele=iele+1
               if(iele>nmax) call rx( ' error in matmsparce: enlarge nmax')
               i1(iele)=i1x
               i2(iele)=i2x
               zz(iele)=a(i1x,i2x)
            endif
         enddo
         enddo
         c=0d0
         do i=1,iele
           c(i1(i),:)= zz(i)*b(i2(i),:)
         enddo
      elseif(isparse==2) then
         iele=0
         do i1x=1,n
         do i2x=1,n
            if(abs(b(i1x,i2x))>etol) then
               iele=iele+1
               if(iele>nmax) call rx( ' error in matmsparce: enlarge nmax')
               i1(iele)=i1x
               i2(iele)=i2x
               zz(iele)=b(i1x,i2x)
            endif
         enddo
         enddo
         c=0d0
         do i=1,iele
           c(:,i2(i))= a(:,i1(i))*zz(i)
         enddo
      endif
      deallocate(i1,i2,zz)
      end


!> interface for matrix multiplication  c=a*b -------------------------
      subroutine matm(a,b,c,n1,n2,n3)
      integer(4), intent(in) :: n1,n2,n3
      complex(8), intent(in) :: a(n1,n2), b(n2,n3)
      complex(8), intent(out) :: c(n1,n3)
      call ZGEMM ( "N", "N", n1, n3, n2, dcmplx(1d0,0d0), 
     &             a, n1, 
     &             b, n2,
     &             dcmplx(0d0,0d0), c, n1 )
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
