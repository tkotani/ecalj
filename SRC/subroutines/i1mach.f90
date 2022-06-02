!> taken from https://community.intel.com/t5/Intel-Fortran-Compiler/Weird-Fortran/td-p/1185072?
!> for f90      
      function i1mach(i) result(s)
      implicit none
      integer :: i,s
      if(i==2.or.i==4) then
         s=6
      elseif(i==9) then
         s=huge(0)
      else
         call rx('i1mach not defined')
      endif   
      end function i1mach
c$$$Cr Remarks
c$$$C   dmach(1-3) are as returned by the BLAS subroutine dmach and are
c$$$C   defined as follows.
c$$$C        b = base of arithmetic
c$$$C        t = number of base b digits
c$$$C        l = smallest possible exponent
c$$$C        u = largest possible exponent
c$$$C   dmach(1): eps = b**(1-t)
c$$$C   dmach(2): tiny = 100.0*b**(-l+t)
c$$$C   dmach(3): huge = 0.01*b**(u-t)
c$$$C
c$$$C   d1mach(1-5) are as returned by the BLAS subroutine d1mach and are
c$$$C   defined as follows.
c$$$C   d1mach(1) = b**(l-1), the smallest positive magnitude.
c$$$C   d1mach(2) = b**(u*(1 - b**(-t))), the largest magnitude.
c$$$C   d1mach(3) = b**(-t), the smallest relative spacing.
c$$$C   d1mach(4) = b**(1-t), the largest relative spacing.
c$$$C   d1mach(5) = log10(b)
c$$$C   d1mach and dmach call the C segment mkcon found in fsubs.c
c$$$C ----------------------------------------------------------------------
      function d1mach(i) result(s)
      implicit none
      integer:: i
      double precision s,dm(5)
      logical :: beg = .true.
      save dm
      if(i.lt.1.or.i.gt.5)stop 'D1MACH(arg < 1 or arg > 5)'
      if(beg)then
         beg=.false.
         dm(1) = tiny(0d0)
         dm(2) = huge(0d0)
         dm(3) = epsilon(0d0)/radix(0d0)
         dm(4) = epsilon(0d0)
         dm(5) = log10(2d0)
      end if
      s = dm(i)
      return
      end function d1mach
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      function dmach(i) result(s)
      implicit none
      integer:: i
      double precision s,dm(5)
      logical :: beg = .true.
      real(8):: t,b,l,u,dlamch,eps
      save dm
      if(i.lt.1.or.i.gt.3)stop 'DMACH(arg < 1 or arg > 3)'
      if(beg)then
         beg=.false.
c     machine constant from lapack routine
         b = dlamch('b')        !base
         eps = dlamch('p')      !eps*base
         l = dlamch('m')        !emin
         u = dlamch('l')        !emax
         t = int(1d0-(log(eps)/log(b)))
         dm(1) = b**(1d0-t)
         dm(2) = 100d0*b**(l+t)
         dm(3) = (b**(u-t))/100d0
      endif
      s = dm(i)
      return
      end function dmach
