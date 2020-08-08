#define RIEEE 1
      double precision function d1mach(i)
C- Machine constants, C-version
C ----------------------------------------------------------------------
Ci Inputs: i
Ci   
Co Outputs:
Co       Returns machine constant values, called as d1mach or dmach
Cr Remarks
C   dmach(1-3) are as returned by the BLAS subroutine dmach and are
C   defined as follows.
C        b = base of arithmetic
C        t = number of base b digits
C        l = smallest possible exponent
C        u = largest possible exponent
C   dmach(1): eps = b**(1-t)
C   dmach(2): tiny = 100.0*b**(-l+t)
C   dmach(3): huge = 0.01*b**(u-t)
C
C   d1mach(1-5) are as returned by the BLAS subroutine d1mach and are
C   defined as follows.
C   d1mach(1) = b**(l-1), the smallest positive magnitude.
C   d1mach(2) = b**(u*(1 - b**(-t))), the largest magnitude.
C   d1mach(3) = b**(-t), the smallest relative spacing.
C   d1mach(4) = b**(1-t), the largest relative spacing.
C   d1mach(5) = log10(b)
C   d1mach and dmach call the C segment mkcon found in fsubs.c
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer i
C Local variables
      double precision dmach,dm(3),d1m(5)

      if (i .lt. 1  .or.  i .gt. 5) then
        call rx ('d1mach i out of bounds')
      endif
      call mkdcon(dm,d1m)
      d1mach = d1m(i)
      return

      entry dmach(i)
      call mkdcon(dm,d1m)
      if (i .lt. 1  .or.  i .gt. 3) then
        call rx ('dmach i out of bounds')
      endif
      dmach = dm(i)
      end
C      subroutine fmain
C      implicit none
C      integer i
C      double precision d1mach,dmach
C      print *, 'testing d1mach ... '
C      do  1  i = 1, 5
C        print *, 'i, d1mach(i)=', i, d1mach(i)
C    1 continue 
C      print *, 'testing dmach ... '
C      do  2  i = 1, 3
C        print *, 'i, dmach(i)=', i, dmach(i)
C    2 continue 
C      end

