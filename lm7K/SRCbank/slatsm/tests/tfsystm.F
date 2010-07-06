      subroutine fmain

      print *, '... testing: call fsystm(''echo hi there'')'
      call fsystm('echo hi there',j)
      print *, 'fsystm returned j=',j
      print *, ' '
      print *, '... testing: call fsystm(''ls *.f'')'
      call fsystm('ls *.f',j)
      print *, 'fsystm returned j=',j
      end

