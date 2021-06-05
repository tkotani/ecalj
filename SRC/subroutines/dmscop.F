      subroutine dmscop(dest,ldd,src,lds,n1rs,n2rs,n1cs,n2cs,
     .n1rd,n1cd,fac)
C- Shift and copy a matrix
C  Inputs: 
Ci   dest,ldd  destination matrix and leading dimension
Ci   src, lds  source matrix and leading dimension
Ci   n1rs,n2rs rows src matrix to copy
Ci   n1cs,n2cs cols src matrix to copy
Ci   n1rd,n1cd starting row, col to copy to dest.
Ci   fac       scale dest by factor fac
Co Outputs
Co    dest(n1rd:n1rd+n2rs-n1rs-1,n1cd:n1cd+n2cs-n1cs-1)
Co             is copied from fac * src(n1rs:n2rs,n1cs:n2cs)
C ----------------------------------------------------------------
C     implicit none
      integer ldd,lds,n1rs,n2rs,n1cs,n2cs,
     .n1rd,n1cd
      double precision dest(ldd,1),src(lds,1),fac
C Local variables
      integer i,j,offr,offc

      offr = n1rd-n1rs
      offc = n1cd-n1cs
      if (fac .eq. 1d0) then
        do    j = n1cs, n2cs
          do  i = n1rs, n2rs
             dest(i+offr,j+offc) = src(i,j)
          enddo
       enddo
      else
         do    j = n1cs, n2cs
            do  i = n1rs, n2rs
               dest(i+offr,j+offc) = fac*src(i,j)
            enddo
         enddo
      endif
      end
C      subroutine fmain
C      implicit none
C      integer ldd,lds,n1rs,n2rs,n1cs,n2cs,n1rd,n1cd,i,j
C      parameter (ldd=10,lds=8)
C      double precision dest(ldd,ldd),src(lds,lds),fac
C
C      n1rs = 2
C      n2rs = min(ldd,lds)-1
C      n1cs = 3
C      n2cs = min(ldd,lds)-2
C      n1rd= 4
C      n1cd= 5
C
C      do  10  i = 1, ldd
C      do  10  j = 1, ldd
C   10 dest(i,j) = -1*0
C
C      do  20  i = 1, lds
C      do  20  j = 1, lds
C   20 src(i,j) = i+10*j ! dble(i)/dble(j)
C
C      fac = 10d0
C      call dmscop(dest,ldd,src,lds,n1rs,n2rs,n1cs,n2cs,
C     .  n1rd,n1cd,fac)
C
C      call ywrm(0,' ',1,6,'(10f13.8)',dest,ldd*ldd,ldd,ldd,ldd)
C      end

