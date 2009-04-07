      subroutine fmain
C test s8tor8
      implicit none
      double precision x
      character strn5*5,strn6*6,strn7*7,strn8*8,strn9*9

      strn9 = '123456789'

      print *, 'storing strn9 ', strn9
      print *, 'recover strn5..9 from x:'

      call s8tor8(strn9,x)
      call r8tos8(x,strn5)
      print *, 'strn5 = ', '|',strn5,'|'
      call r8tos8(x,strn6)
      print *, 'strn6 = ', '|',strn6,'|'
      call r8tos8(x,strn7)
      print *, 'strn7 = ', '|',strn7,'|'
      call r8tos8(x,strn8)
      print *, 'strn8 = ', '|',strn8,'|'
      call r8tos8(x,strn9)
      print *, 'strn9 = ', '|',strn9,'|'

      print *
      print *, 'storing strn5 ', strn5
      call s8tor8(strn5,x)
      print *, 'recover strn5..9 from x:'
      call r8tos8(x,strn5)
      print *, 'strn5 = ', '|',strn5,'|'
      call r8tos8(x,strn6)
      print *, 'strn6 = ', '|',strn6,'|'
      call r8tos8(x,strn7)
      print *, 'strn7 = ', '|',strn7,'|'
      call r8tos8(x,strn8)
      print *, 'strn8 = ', '|',strn8,'|'
      call r8tos8(x,strn9)
      print *, 'strn9 = ', '|',strn9,'|'


      end

