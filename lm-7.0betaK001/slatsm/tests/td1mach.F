      subroutine fmain
c     program test
      implicit none
      integer i
      double precision d1mach,dmach
      character*3  yn
      print *, 'testing d1mach ...  numbers should be:'
      print *, '1,2: smallest positive and largest magnitude;'
      print *, '3,4: the smallest and largest relative spacing;'
      print *, '5: log(2)/log(10).'
      print *, ' '
      do  10  i = 1, 5
   10 print *, 'i, d1mach(i)=', i, d1mach(i)

      print *, ' '
      print *, 'testing dmach ...  numbers should be:'
      print *, '1 machine epsilon,  2 tiny number,  3 huge number'

      do  20  i = 1, 3
   20 print *, 'i, dmach(i)=', i, dmach(i)

      print *, ' '
      yn = 'no'
      if (d1mach(1) .le. dmach(2)) yn = 'yes'
      print *, 'td1mach: d1mach(1) .le. dmach(2) ? ', yn
      yn = 'no'
      if (d1mach(2) .ge. dmach(3)) yn = 'yes'
      print *, 'td1mach: d1mach(2) .ge. dmach(3) ? ', yn
      yn = 'no'
      if (d1mach(3) .le. dmach(1)) yn = 'yes'
      print *, 'td1mach: d1mach(3) .le. dmach(1) ? ', yn

      if (d1mach(1).gt.dmach(2) .or. d1mach(1).lt.0) call cexit(-1,1)
      if (d1mach(2).lt.dmach(3) .or. d1mach(2).lt.0) call cexit(-1,1)
      if (d1mach(3).gt.dmach(1) .or. d1mach(3).lt.0) call cexit(-1,1)


      end
