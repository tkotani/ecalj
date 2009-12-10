      subroutine lgand(sname,struc,lval,mask)
C- Returns a logical AND of a vector of
C ----------------------------------------------------------------
Ci Inputs
Ci   lval is a logical T or F
Ci   mask should be an integer multiple of 2.
Ci        Only the lowest bit of mask is used.
Co Outputs
Co  struc element corresponding to label 'name' is modified.
Co        The mask bit of that entry is set to lval.
C ----------------------------------------------------------------
C     implicit none
      logical lval
      integer mask
      character*(*) sname
      integer ix1(10),maski,nmask,t2n
      double precision struc(1),x2(10),x3(10),x4(10),x5(10)

      if (mask .le. 0) call rxi('lsets: mask must be >0 but given',mask)
      call spack(0,sname,struc,ix1,x2,x3,x4,x5)
      maski = mask
      nmask = ix1(1)
      t2n = 1
C ... Find lowest nonzero bit of mask, corresponding ix1
   10 continue
      if (mod(maski,2) .eq. 0) then
        t2n = t2n*2
        maski = maski/2
        ix1(1) = ix1(1)/2
        goto 10
      endif
      if (lval) then
        nmask = nmask + t2n*(1-mod(ix1(1),2))
      else
        nmask = nmask + t2n*(0-mod(ix1(1),2))
      endif
      call spack(1,sname,struc,nmask,x2,x3,x4,x5)

      end

