      subroutine dmsadd(dest,ldd,src,lds,n1rs,n2rs,n1cs,n2cs,
     .  n1rd,n1cd,fac)
C- Shift and add a subblock of one matrix to another
C  Inputs: 
Ci   dest,ldd  destination matrix and leading dimension
Ci   src, lds  source matrix and leading dimension
Ci   n1rs,n2rs rows src matrix to copy
Ci   n1cs,n2cs cols src matrix to copy
Ci   n1rd,n1cd starting row, col to copy to dest.
Ci   fac       add fac*src into dest
Co Outputs
Co    fac * src(n1rs:n2rs,n1cs:n2cs) is added to 
Co    dest(n1rd:n1rd+n2rs-n1rs-1,n1cd:n1cd+n2cs-n1cs-1)
C ----------------------------------------------------------------
C     implicit none
      integer ldd,lds,n1rs,n2rs,n1cs,n2cs,
     .  n1rd,n1cd
      double precision dest(ldd,1),src(lds,1),fac
C Local variables
      integer i,j,offr,offc

      offr = n1rd-n1rs
      offc = n1cd-n1cs
      if (fac .eq. 1d0) then
        do  10  j = n1cs, n2cs
        do  10  i = n1rs, n2rs
   10   dest(i+offr,j+offc) = dest(i+offr,j+offc) + src(i,j)
      else
        do  20  j = n1cs, n2cs
        do  20  i = n1rs, n2rs
   20   dest(i+offr,j+offc) = dest(i+offr,j+offc) + fac*src(i,j)
      endif
      end

