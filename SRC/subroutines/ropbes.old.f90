      subroutine ropbes(r,e,lmax,y,h,xi,n)
C- Radial bessel functions divided by r**l (vectorizes)
C ----------------------------------------------------------------
Ci Inputs
Ci   r    list of points
Ci   e    energy
Ci   y,h  work vectors of length n each
Co Outputs
Co   xi   J(r,l)/r**l, according to standard definition
Cr Remarks
Cr   J(r,lmax,lmax-1)  are calculated by a power series.
Cr   J for lower l are calculated by downward recursion
C ----------------------------------------------------------------
C     implicit none
      integer lmax,n
      double precision e,r(n),xi(n,0:lmax),h(n),y(n)
      integer i,l
      double precision xx

      if (lmax .lt. 0) return
      do  8  i = 1, n
        y(i) = r(i)*r(i)*e
    8 continue
C --- Power series expansion for lmax, lmax-1 ---
      call ropbs1(y,lmax,xi(1,lmax),h,n)
      if (lmax .ge. 1) call ropbs1(y,lmax-1,xi(1,lmax-1),h,n)
C --- Downward recursion ---
      do  10  l = lmax-2, 0, -1
        xx = 2*l+3
        do  1  i = 1, n
          xi(i,l) = xx*xi(i,l+1) - y(i)*xi(i,l+2)
    1   continue
   10 continue
      end
      subroutine ropbs1(y,l,xi,h,n)
C- Evaluates bessel function for one l using power series expansion
C     implicit none
      integer n,l
      double precision xi(n),h(n),y(n)
      integer i,k
      double precision df,tol,top,xx

      tol = 1d-12
      df = 1d0
      do  1  k = 3, 2*l+1,2
        df = df*k
    1 continue
      do  2  i = 1, n
        xi(i) = 1d0
        h(i) = 1d0
    2 continue
      do  10  k = 1, 500
        xx = -1d0/( 2*k*(2*k+2*l+1) )
        do  3  i = 1, n
          h(i) = xx*h(i)*y(i)
          xi(i) = xi(i)+h(i)
    3   continue
        top = 0d0
        do  4  i = 1, n
          top = dmax1(top,dabs(h(i)))
    4   continue
        if (top .le. tol) goto 11
   10 continue
      call rx('ropbes: power series failed to converge')
   11 xx = 1d0/df
      do  5  i = 1, n
        xi(i) = xi(i)*xx
    5 continue
      end
C testing ...
C      subroutine fmain
C      implicit none
C      integer nr,lmax
C      parameter (nr=8,lmax=4)
C      double precision xi(nr,0:lmax),ri(nr),e,y(nr),h(nr)
C      double precision phi(0:lmax+1),psi(0:lmax+1)
C      integer i,l
C      do  10  i = 1, nr
C   10 ri(i) = 2*dble(i)/nr
C      e = -1.5
C      call ropbes(ri,e,lmax,y,h,xi,nr)
C
C      do  20  i = 1, nr
C        call bessl2(e*ri(i)**2,0,lmax,phi,psi)
C        print 333, i, ri(i),(xi(i,l), l=0,lmax)
C        print 333, i, ri(i),(phi(l), l=0,lmax)
C        print 334, ((xi(i,l)-phi(l))*1e12, l=0,lmax)
C  333   format(' i',i3,7f12.6)
C  334   format(5x,12x,7f12.3)
C   20 continue
C
C      end

