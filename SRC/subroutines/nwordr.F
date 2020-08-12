      integer(4) function nwordr()
      real(8):: r
      integer(4):: len
c      nword = 4 ! in alpha
c      nword =1 ! in AIX
c      nword=NWORD_RECORDSIZE
      inquire(iolength=len) r
      nwordr = 8/len
c      write(6,*)' nword=',nword
      end
