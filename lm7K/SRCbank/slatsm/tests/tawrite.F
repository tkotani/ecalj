      subroutine fmain
      implicit none
      integer i1mach,awrite
      double precision a1(10),b(3),big,small,ax(10)
      logical a2bin
      integer a2,j,ix(2)
      integer mxlen
      parameter (mxlen=120)
      character*(mxlen) s,ch
      logical lx(2)



C      ch = ' '
C      call nlchar(1,ch)
C      print *, ch

      j = 80
      s = ' '
      a1(1) = 1
      a1(2) = 2
      a1(3) = 3
      lx(1) = .true.
      lx(2) = .false.
      a2 = -7
      do  30  j = 1, 10
   30 ax(j) = (-1d0)**(j+1)/dble(j)

C --- Test %?expr;strn1;strn2; ---
      call awrit2('three plus one is %?;n==1;%i;four;, no?',
     .  s,mxlen,-i1mach(2),1,4)
      call awrit2('three plus one is %?;n==1;%i;four;, no?',
     .  s,mxlen,-i1mach(2),0,4)
      call awrit2('This string is longer than 28%?#p>28#, I see#, no?#',
     .  s,mxlen,-i1mach(2),0,4)
      call awrit2('This string is longer than 29%?#p>29#, I see#, no?#',
     .  s,mxlen,-i1mach(2),0,4)
      call awrit2('This string is longer than 28%?#p>28#%N so newline!#'
     .  ,s,mxlen,-i1mach(2),0,4)

      s = ' '
      call awrit2('check conditional %%c==X: '//
     .  'cursor now at X?%?#%c==X# ... Yes# ... No!#'
     .  ,s,mxlen,-i1mach(2),0,4)
      s = ' '
      call awrit2('check conditional %%c==X: '//
     .  'cursor now at X?%b%?#%c==X# ... Yes# ... No!#'
     .  ,s,mxlen,-i1mach(2),0,4)

C --- Test %n, %j, %(expr) ---
      call awrit3('test %%n: %n:n,5d',s,mxlen,-i1mach(2),3,4,a1)
      call awrit5('test %%j: %2j%n:n,5d',s,mxlen,-i1mach(2),7,8,3,4,a1)
      call awrit5('test %%j: %(3-1)j%n:n,5d',s,mxlen,-i1mach(2),7,8,3,4
     .  ,a1)
      call awrit4('%xtabs2%((n-p%n)%n)f%jtabs5%((n-p%n)%n)f%j'//
     .  'tabs4%((n-p%n)%n)f%jtabs3%((n-p%n)%n)f%(-1)j'//
     .  'tabs2%((2-p%2)%2)ftabs4%((4-p%4)%4)ftabs3%((3-p%3)%3)ftabs!',
     .  s,mxlen,-i1mach(2),2,5,4,3)

C --- Test %:-arg ---
      ix(1) = -2
      ix(2) = 3
      b(1) =  1.23d0
      b(2) = -4.56d0
      b(3) =  7.89d0
      call awrit2('test '':-nblk'' :%2:-3,4;5d %2:-3,4i',
     .  s,mxlen,-i1mach(2),b,ix)
      call awrit2('test '':-nblk'' :%2:-3,4;5d %2:,4i',
     .  s,mxlen,-i1mach(2),b,ix)
      ix(2) = -2
      ix(1) = 3
      call awrit2('test '':-nblk'' :%2:-3,4;5d %2:-3,4i',
     .  s,mxlen,-i1mach(2),b,ix)
      call awrit2('test '':-nblk'' :%2:-3,4;5d %2:,4i',
     .  s,mxlen,-i1mach(2),b,ix)

      call awrit2('test ''5:-1,3;5#11d'' :%5:-1,3;5#11d',
     .  s,mxlen,-i1mach(2),ax,ix)
      call awrit2('test ''5:-1;5#11d''   :%5:-1;5#11d',
     .  s,mxlen,-i1mach(2),ax,ix)
      call awrit2('test ''5:2,3;5#11d''  :%5:2,3;5#11d',
     .  s,mxlen,-i1mach(2),ax,ix)
      call awrit2('test ''5:-1,2;5#11g'' :%5:-1,2;5#11g',
     .  s,mxlen,-i1mach(2),ax,ix)

C      call awrit1('%:-1,5;5g',s,mxlen,-i1mach(2),4.2265d-01)
C      call awrit1('%:-1,5;5g',s,mxlen,-i1mach(2),-1.1325d-01)
C      call awrit1('%:-1,5;5g',s,mxlen,-i1mach(2),3.0345d-04)
C      call awrit1('%:-1,5;5g',s,mxlen,-i1mach(2),-3.0345d-04)
C      stop

C --- Test %p, %a, %f, %W, %w, %c, %b and %x ---
      s = 'abcdefghijklmnopqrstuvwxyz'
      call awrit1('try to %aappend this',s,mxlen,-i1mach(2),a1)
      call awrit1('%12p:insert:%3aat the end',s,mxlen,-i1mach(2),a1)
      call awrit1('%12p%3fX%a%fZ',s,mxlen,-i1mach(2),a1)
      call awrit1('%12p%3f%b%i%a',s,mxlen,-i1mach(2),a2)
      call awrit1('%10p%W something else! %W1%W2',s,mxlen,-i1mach(2),a2)
      call awrit0('let    us      try %%w ....',s,mxlen,-i1mach(2))
      call awrit0('let    us      try %%w ....%4p%wme%a now!',
     .  s,mxlen,-i1mach(2))
C     call awrit0('%5p%c%f%c%W%f%wc%a',s,mxlen,-i1mach(2))
      call awrit0('%5p%c%f%c%W%f%w%c%1oc%a',s,mxlen,-i1mach(2))
      s = 'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'
      call awrit0('is it clear yet?%a',s,mxlen,-i1mach(2))
      call awrit0('%xhow about now?%a',s,mxlen,-i1mach(2))

C --- Test logicals ---
      call awrit1('here come two (%2:1l) logicals',s,mxlen,i1mach(2),lx)

C --- Test %d ---
      a2=11
      big = -1234.987654321d0
      big = -1234567.8987654321d0
      big = -123498.76543212345d0
      small = 3e-4
      call awrit3(
     .  'big number %;2d  small%:2d  to 3 decimals%:2,3;3d',
     .  s,mxlen,i1mach(2),big,small,small)
      call awrit1('hello there%3:2d dolly',s,mxlen,i1mach(2),a1)
      call awrit2('hello there%3d dolly %%%i--',s,mxlen,i1mach(2),a1,a2)

      print *, '... try awrite to string s:, subsequent write to stdout'
      call awrit1('hello there%3:2d dolly',s,mxlen,0,a1)
      call awrit0(s,' ',-mxlen,i1mach(2))
      call awrit4('try %%,5i format: |%,5i| |%,5i| |%,5i| |%,5i|',
     .  s,mxlen,i1mach(2),1,12,123,123456)

C --- Test d format ---
      call awrit2('d (big) %1,0;2d %1,2;2d',s,mxlen,0,1.23d7,1.23d7)
      print 333, s(1:72)
      call awrit4('d (middle) %d %,2;2d %,0;4d %,4;4d',
     .  s,mxlen,0,1.2349d0,1.2349d0,1.234951d0,1.234951d0)
      print 333, s(1:72)
      call awrit4('d (small) %d %,2;2d %,0;10d %,14;15d',
     .  s,mxlen,0,1.2349501d-6,1.2349501d-6,1.2349501d-6,1.2349501d-6)
      print 333, s(1:72)
      call awrit4('d (-small) %d %,2;2d %,0;10d %,14;15d',s,mxlen,0,
     .  -1.2349501d-6,-1.2349501d-6,-1.2349501d-6,-1.2349501d-6)
      print 333, s(1:72)

C --- Test D format ---
      s = ' '
      call awrit2('D (big) %,0;12D %,2;12D',s,mxlen,0,1.23d7,1.23d7)
      print 333, s(1:72)
      call awrit4('D (middle) %D %,2;10D %,0;6D %,4;6D',
     .  s,mxlen,0,1.2349d0,1.2349d0,1.234951d0,1.234951d0)
      print 333, s(1:72)
      call awrit4('D (small) %D %,2;10D %,0;10D %,14;15D',
     .  s,mxlen,0,1.2349501d-6,1.2349501d-6,1.2349501d-6,1.2349501d-6)
      print 333, s(1:72)


C --- Test e format ---
      call awrit6('e (big) %e %,4e %;2e %;3e %,3;4e %,5;5e',
     .  s,mxlen,0,1.23d7,1.23d7,1.23d7,1.23d7,1.23d7,1.234951d7)
      call awrit2('%a %;2e %;2e',s,mxlen,0,1.2499d7,1.2501d7)
      print 333, s(1:72)
      call awrit6('e (middle) %e %,4e %;2e %;3e %,3;4e %,5;5e',
     .  s,mxlen,0,1.23d1,1.23d1,1.23d1,1.23d1,1.23d1,1.234951d1)
      print 333, s(1:72)
      j = awrite('e (small) %e %,4e %;2e %;3e %,3;4e %,5;5e',
     .  s,mxlen,0,1.23d-5,1.23d-5,1.23d-5,1.23d-5,1.23d-5,1.234951d-5,
     .  0,0)
      print 333, s(1:j)

C --- Test F format ---
      s = ' '
      call awrit5('F (12) %;12F %;12F %;12F %;12F %;12F',
     .  s,mxlen,0,1.23456789d11,1.23456789d10,1.23456789d2,
     .  1.23456789d-3,1.23456789d-4)
      print 333, s(1:80)
      call awrit7('F (7)  %;7F %;7F %;7F %;7F %;7F %;7F %;7F',
     .  s,mxlen,0,1.23456789d6,1.23456789d5,1.23456789d2,1.23456789d1,
     .  1.23456789d-1,1.23456789d-3,1.23456789d-4)
      print 333, s(1:72)
      s = ' '
      call awrit6('F (s)  %;6F %;5F %;4F %;3F %;2F %;6F',s,mxlen,0,
     .  -1.4926251754d-3,-1.4926251754d-3,-1.4926251754d-3,
     .  -1.4926251754d-3,-1.4926251754d-3,-1.9971712021d-4)
      print 333, s(1:72)

C --- Test g format ---
      s = ' '
      call awrit5('g (very big) %g %,4g %;1g %;2g %;2g',
     .  s,mxlen,0,1.23d197,1.23d97,1.2499d97,1.2499d97,1.2501d97)
      call skpblb(s,72,j)
      print 333, s(1:j)
  333 format(a)
      s = ' '
      call awrit5('g (big) %g %,4g %;1g %,2;2g %,2;2g',
     .  s,mxlen,0,1.23d7,1.23d7,1.2499d7,1.2499d7,1.2501d7)
      print 333, s(1:72)
      call awrit5('g (middle) %g %,4g %;1g %,2;2g %,2;2g,',
     .  s,mxlen,0,1.23d0,1.23d0,1.2499d0,1.2499d0,1.2501d0)
      call awrit2('%a %,0;5g %,5;5g',s,mxlen,0,1.23d0,1.23d0)
      print 333, s(1:72)
      call awrit5('g (small) %10z%g %,4g %;1g %,2;2g %,2;2g',
     .  s,mxlen,0,1.23d-5,1.23d-4,1.2499d-5,1.2499d-5,1.2501d-5)
      call awrit4('%a G %;2G %,3;3G %,4;4G %;4G',s,mxlen,0,
     .  1.23d-5,1.24996d-5,1.24996d-5,1.24996d-5)
      print 333, s(1:72)
      call awrit5('g (zero) %g %,2;2g %,2;2G %0z%,2;2g %,2;2e',
     .  s,mxlen,0,0d0,0d0,0d0,0d0,0d0)
      print 333, s(1:72)

      call cexit(0,1)
      end
