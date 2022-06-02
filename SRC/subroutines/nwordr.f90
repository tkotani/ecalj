integer(4) function nwordr()
  real(8):: r
  integer(4):: len
  !      nword = 4 ! in alpha
  !      nword =1 ! in AIX
  !      nword=NWORD_RECORDSIZE
  inquire(iolength=len) r
  nwordr = 8/len
  !      write(6,*)' nword=',nword
END function nwordr
