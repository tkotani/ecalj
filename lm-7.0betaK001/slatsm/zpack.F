      subroutine zpack(mode,modes,moded,nr1,nr2,nc1,nc2,kl,ku,lds,ldd,
     .  src,dest)
C- Copy and pack/unpack a double complex matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode : 0 copies src into dest
Ci          1 adds src into dest
Ci   modes :defines how src is packed
Ci          0 matrix is not packed
Ci          1 matrix is packed in LAPACK banded form; see Remarks
Ci            src(kl+ku+1+i-j,j) points to entry src(i,j)
Ci            for max(1,j-ku)<=i<=min(m,j+kl)
Ci   moded :defines how dest is packed; see modes for definitions
Ci  nr1,nr2:copy array rows nr1..nr2
Ci          (a subset of this if src or dest is in banded form)
Ci  nc1,nc2:copy array columns nc1..nc2
Ci   kl    :parameter kl for banded storage (see Remarks)
Ci          (not used if neither src nor dest is in banded form)
Ci   ku    :parameter ku for banded storage (see Remarks)
Ci          (not used if neither src nor dest is in banded form)
Ci   lds   :leading dimension of src.  If modes=1, lds >= 2*kl+ku+1
Ci   ldd   :leading dimension of dest. If moded=1, ldd >= 2*kl+ku+1
Ci   src   :source matrix
Co Outputs
Co   dest  :destination matrix
Cr Remarks
Cr   mode=1: LAPACK banded storage of matrix a_ij, stored in matrix A
Cr     Let kl = the number of subdiagonals within the band of a_ij
Cr     Let ku = the number of superdiagonals within the band of a_ij
Cr     a_ij is stored in rows kl+1 to 2*kl+ku+1 of A
Cr     rows 1 to kl of A are not used in storage.
Cr
Cr     The j-th column of a_ij is stored in the j-th column of A as:
Cr     A(kl+ku+1+i-j,j) = a_ij for max(1,j-ku) <= i <= min(m,j+kl)
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,modes,moded,lds,nr1,nr2,nc1,nc2,ldd,kl,ku
      double complex src(lds,nc2),dest(ldd,nc2)
C ... Local parameters
      integer i,j,ofc,nr,nc
      double precision fac

      fac = 0
      if (mode .eq. 1) fac = 1
      if (modes .eq. 0 .and. moded .eq. 0) then
        nr = nr2-nr1+1
        nc = nc2-nc1+1
        call zmscop(mode,nr,nc,lds,ldd,0,0,0,0,src(nr1,nc1),
     .              dest(nr1,nc1))
      elseif (modes .eq. 0 .and. moded .eq. 1) then
        if (ldd .lt. 2*kl+ku+1) call xerbla('ZPACK ',11)
        ofc = kl+ku+1
        do  10  j = nc1, nc2
        do  10  i = max(nr1,j-ku), min(nr2,j+kl)
   10   dest(i+ofc-j,j) = fac*dest(i+ofc-j,j) + src(i,j)
      elseif (modes .eq. 1 .and. moded .eq. 0) then
        if (lds .lt. 2*kl+ku+1) call xerbla('ZPACK ',10)
        ofc = kl+ku+1
        do  20  j = nc1, nc2
        do  20  i = max(nr1,j-ku), min(nr2,j+kl)
   20   dest(i,j) = fac*dest(i,j) + src(i+ofc-j,j)
      elseif (modes .eq. 1 .and. moded .eq. 1) then
        if (lds .lt. 2*kl+ku+1) call xerbla('ZPACK ',11)
        if (ldd .lt. 2*kl+ku+1) call xerbla('ZPACK ',10)
        ofc = kl+ku+1
        do  30  j = nc1, nc2
        do  30  i = max(nr1,j-ku), min(nr2,j+kl)
   30   dest(i+ofc-j,j) = fac*dest(i+ofc-j,j) + src(i+ofc-j,j)
      else
        call fexit2(-1,111,' Exit -1 ZPACK: '//
     .    'modes=%i and moded=%i not implemented',modes,moded)
      endif
      end

