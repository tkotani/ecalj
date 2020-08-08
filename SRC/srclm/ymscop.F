      subroutine ymscop(mode,nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .src,offsi,dest,offdi)
C- Copy a subblock block of a complex matrix to a destination
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode : 1s digit
Ci          0 copies src into dest
Ci          1 adds src into dest
Ci          2 copies -src into dest
Ci          3 adds -src into dest
Ci      add 4 if src and dest matrix are real
Ci   nlma : row subblock size
Ci   nlmb : col subblock size
Ci   ndas : effective leading dimension of source matrix
Ci        : = number of elements separating columns in source matrix
Ci        : Thus, ndas:
Ci        : = leading dimension of src if real, imag fully separated
Ci        : = 2* leading dimension of src if real, imag sep by cols
Ci   ndad : effective leading dimension of destination matrix
Ci        : = number of elements separating columns in destination matrix
Ci        : Thus, ndad:
Ci        : = leading dimension of dest if real, imag fully separated
Ci        : = 2* leading dimension of dest if real, imag sep by cols
Ci   ofas : offset to first row of source matrix
Ci   ofbs : offset to first col of source matrix
Ci   ofad : offset to first row of destination matrix
Ci   ofbd : offset to first col of destination matrix
Ci   src  : source matrix
Ci   offsi: separation between real,imaginary parts of src
Ci   offdi: separation between real,imaginary parts of dest
Co Outputs
Co   dest : destination matrix
Cr Remarks
Cu Updates
Cu   28 Jan 01 Added '2' to mode for copying -1*matrix
Cu   02 Jan 01 Added '4' to mode for copying real matrices
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,offsi,offdi
      double precision src(ndas,1),dest(ndad,1)
C ... Local parameters
      integer mode4
      double precision alpha(2),beta(2)

      alpha(1) = 1
      alpha(2) = 0
      beta(1) = 0
      mode4 = mod(mode,4)
      if (mode4 .ge. 2) alpha(1) = -1
      if (mode4 .eq. 1 .or. mode4 .eq. 3) beta(1) = 1
      beta(2) = 0

      if (mode .eq. 4) then
        call dmscop(dest,ndad,src,ndas,ofas+1,ofas+nlma,ofbs+1,ofbs+
     .  nlmb,ofad+1,ofbd+1,alpha)
      else if (mode .eq. 5) then
        call dmsadd(dest,ndad,src,ndas,ofas+1,ofas+nlma,ofbs+1,ofbs+
     .  nlmb,ofad+1,ofbd+1,alpha)
      else if (mode .gt. 5) then
        call rx('ymscop: mode not implemented')
      else
        call ymsadd(nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .  alpha,beta,src,offsi,dest,offdi)
      endif

      end

