      integer function nlme(nl)
C- Returns no. of Slater-Koster matrix elements, given no. of l's nl
      implicit none
C Passed parameters
      integer nl
C Local parameters
      nlme = (nl+2)*(nl+1)*nl/6
      end
