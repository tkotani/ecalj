      subroutine fmain
      implicit none
      integer ls,n,i
      logical cmdopt
      character strn*72,ch*1
      ls = len('''hi there, johnny!'' and now, what?')
      call acmdop('''hi there, johnny!'' and now, what?',ls,0)
      ls = len('here are some "new args" ')
      call acmdop('here are some "new args" ',ls,2)
      call spchar(1,ch)
      ls = len('funny ''  ~'//ch//'{W'//ch//'}''   ')
      call acmdop(' funny ''  ~'//ch//'{W'//ch//'}''   ',ls,2)
      call ncmdop(n)
      do  10  i = 1, n+1
        call gcmdop(i,strn)
        call awrit1(' %i '//strn,' ',80,6,i)
   10 continue
      print *, ' '

      if (cmdopt('abc',3,0,strn)) then      
        print *, 'we found "abc"'
      else
        print *, 'we did not find "abc"'
      endif

      if (cmdopt('and',3,0,strn)) then      
        print *, 'we found "and"'
      else
        print *, 'we did not find "and"'
      endif

      if (cmdopt('what!',5,0,strn)) then      
        print *, 'we found "what!"'
      else
        print *, 'we did not find "what!"'
      endif

      if (cmdopt('what',4,0,strn)) then      
        print *, 'we found "what"'
      else
        print *, 'we did not find "what"'
      endif


      end
