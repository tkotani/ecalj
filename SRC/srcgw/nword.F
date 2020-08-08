      integer(4) function nword()
      real(8):: r
      integer(4):: len
c      nword = 4 ! in alpha
c      nword =1 ! in AIX
c      nword=NWORD_RECORDSIZE
      inquire(iolength=len) r
      nword = 8/len
c      write(6,*)' nword=',nword
      end
