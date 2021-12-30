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
