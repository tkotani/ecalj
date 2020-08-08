      subroutine zmsadd(nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .alpha,beta,src,dest)
C- Add a subblock block of a double complex matrix to a destination
C ----------------------------------------------------------------------
Ci Inputs
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
Co   dest : destination matrix: dest is overwritten by
Co        : dest <- alpha * src + beta * dest
Cr Remarks
Cb Bugs
Cb   This routine has not been thoroughly checked ...
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd
      double precision src(2,ndas,1),dest(2,ndad,1),alpha(2),beta(2)
C ... Local parameters
      integer ia,ib

      if (beta(1) .eq. 0 .and. beta(2) .eq. 0) then
        if (alpha(1) .eq. 1 .and. alpha(2) .eq. 0) then
          do  10  ib = 1, nlmb
            do  10  ia = 1, nlma
              dest(1,ia+ofad,ib+ofbd) = src(1,ia+ofas,ib+ofbs)
              dest(2,ia+ofad,ib+ofbd) = src(2,ia+ofas,ib+ofbs)
   10     continue
        else
          do  12  ib = 1, nlmb
            do  12  ia = 1, nlma
              dest(1,ia+ofad,ib+ofbd) =
     .        alpha(1)*src(1,ia+ofas,ib+ofbs) -
     .        alpha(2)*src(2,ia+ofas,ib+ofbs)
              dest(2,ia+ofad,ib+ofbd) =
     .        alpha(1)*src(2,ia+ofas,ib+ofbs) +
     .        alpha(2)*src(1,ia+ofas,ib+ofbs)
   12     continue
        endif
      elseif (beta(1) .eq. 1 .and. beta(2) .eq. 0) then
        if (alpha(1) .eq. 1 .and. alpha(2) .eq. 0) then
          do  20  ib = 1, nlmb
            do  20  ia = 1, nlma
              dest(1,ia+ofad,ib+ofbd) = dest(1,ia+ofad,ib+ofbd) +
     .        src(1,ia+ofas,ib+ofbs)
              dest(2,ia+ofad,ib+ofbd) = dest(2,ia+ofad,ib+ofbd) +
     .        src(2,ia+ofas,ib+ofbs)
   20     continue
        else
          do  22  ib = 1, nlmb
            do  22  ia = 1, nlma
              dest(1,ia+ofad,ib+ofbd) = dest(1,ia+ofad,ib+ofbd) +
     .        alpha(1)*src(1,ia+ofas,ib+ofbs) -
     .        alpha(2)*src(2,ia+ofas,ib+ofbs)
              dest(2,ia+ofad,ib+ofbd) = dest(2,ia+ofad,ib+ofbd) +
     .        alpha(1)*src(2,ia+ofas,ib+ofbs) +
     .        alpha(2)*src(1,ia+ofas,ib+ofbs)
   22     continue
        endif
      else
        do  30  ib = 1, nlmb
          do  30  ia = 1, nlma
            dest(1,ia+ofad,ib+ofbd) =
     .      beta(1)*dest(1,ia+ofad,ib+ofbd) -
     .      beta(2)*dest(2,ia+ofad,ib+ofbd) +
     .      alpha(1)*src(1,ia+ofas,ib+ofbs) -
     .      alpha(2)*src(2,ia+ofas,ib+ofbs)
            dest(2,ia+ofad,ib+ofbd) =
     .      beta(1)*dest(2,ia+ofad,ib+ofbd) +
     .      beta(2)*dest(1,ia+ofad,ib+ofbd) +
     .      alpha(1)*src(2,ia+ofas,ib+ofbs) +
     .      alpha(2)*src(1,ia+ofas,ib+ofbs)
   30   continue
      endif
      end

