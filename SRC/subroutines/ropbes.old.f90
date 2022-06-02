subroutine ropbes(r,e,lmax,y,h,xi,n)
  !- Radial bessel functions divided by r**l (vectorizes)
  ! ----------------------------------------------------------------
  !i Inputs
  !i   r    list of points
  !i   e    energy
  !i   y,h  work vectors of length n each
  !o Outputs
  !o   xi   J(r,l)/r**l, according to standard definition
  !r Remarks
  !r   J(r,lmax,lmax-1)  are calculated by a power series.
  !r   J for lower l are calculated by downward recursion
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: lmax,n
  double precision :: e,r(n),xi(n,0:lmax),h(n),y(n)
  integer :: i,l
  double precision :: xx

  if (lmax < 0) return
  do  8  i = 1, n
     y(i) = r(i)*r(i)*e
8 enddo
  ! --- Power series expansion for lmax, lmax-1 ---
  call ropbs1(y,lmax,xi(1,lmax),h,n)
  if (lmax >= 1) call ropbs1(y,lmax-1,xi(1,lmax-1),h,n)
  ! --- Downward recursion ---
  do  10  l = lmax-2, 0, -1
     xx = 2*l+3
     do  1  i = 1, n
        xi(i,l) = xx*xi(i,l+1) - y(i)*xi(i,l+2)
1    enddo
10 enddo
end subroutine ropbes
subroutine ropbs1(y,l,xi,h,n)
  !- Evaluates bessel function for one l using power series expansion
  !     implicit none
  integer :: n,l
  double precision :: xi(n),h(n),y(n)
  integer :: i,k
  double precision :: df,tol,top,xx

  tol = 1d-12
  df = 1d0
  do  1  k = 3, 2*l+1,2
     df = df*k
1 enddo
  do  2  i = 1, n
     xi(i) = 1d0
     h(i) = 1d0
2 enddo
  do  10  k = 1, 500
     xx = -1d0/( 2*k*(2*k+2*l+1) )
     do  3  i = 1, n
        h(i) = xx*h(i)*y(i)
        xi(i) = xi(i)+h(i)
3    enddo
     top = 0d0
     do  4  i = 1, n
        top = dmax1(top,dabs(h(i)))
4    enddo
     if (top <= tol) goto 11
10 enddo
  call rx('ropbes: power series failed to converge')
11 xx = 1d0/df
  do  5  i = 1, n
     xi(i) = xi(i)*xx
5 enddo
end subroutine ropbs1
! testing ...
!      subroutine fmain
!      implicit none
!      integer nr,lmax
!      parameter (nr=8,lmax=4)
!      double precision xi(nr,0:lmax),ri(nr),e,y(nr),h(nr)
!      double precision phi(0:lmax+1),psi(0:lmax+1)
!      integer i,l
!      do  10  i = 1, nr
!   10 ri(i) = 2*dble(i)/nr
!      e = -1.5
!      call ropbes(ri,e,lmax,y,h,xi,nr)

!      do  20  i = 1, nr
!        call bessl2(e*ri(i)**2,0,lmax,phi,psi)
!        print 333, i, ri(i),(xi(i,l), l=0,lmax)
!        print 333, i, ri(i),(phi(l), l=0,lmax)
!        print 334, ((xi(i,l)-phi(l))*1e12, l=0,lmax)
!  333   format(' i',i3,7f12.6)
!  334   format(5x,12x,7f12.3)
!   20 continue

!      end

