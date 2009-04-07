      subroutine zmscop(mode,nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .  src,dest)
C- Copy a subblock block of a complex matrix to a destination
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode : 0 copies src into dest
Ci          1 adds src into dest
Ci   nlma : row subblock size
Ci   nlmb : col subblock size
Ci   ndas : row dimension of source matrix
Ci   ndad : row dimension of destination matrix
Ci   ofas : offset to first row of source matrix
Ci   ofbs : offset to first col of source matrix
Ci   ofad : offset to first row of destination matrix
Ci   ofbd : offset to first col of destination matrix
Ci   src  : source matrix
Co Outputs
Co   dest : destination matrix
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd
      double precision src(2,ndas,1),dest(2,ndad,1)
C ... Local parameters
      double precision alpha(2),beta(2)

      alpha(1) = 1
      alpha(2) = 0
      beta(1) = 0
      beta(2) = 0
      if (mode .eq. 1) beta(1) = 1

      call zmsadd(nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .  alpha,beta,src,dest)

      end

