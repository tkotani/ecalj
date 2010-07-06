C$$$      subroutine lgand(sname,struc,lval,mask)
C$$$C- Returns a logical AND of a vector of
C$$$C ----------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   lval is a logical T or F
C$$$Ci   mask should be an integer multiple of 2.
C$$$Ci        Only the lowest bit of mask is used.
C$$$Co Outputs
C$$$Co  struc element corresponding to label 'name' is modified.
C$$$Co        The mask bit of that entry is set to lval.
C$$$C ----------------------------------------------------------------
C$$$C     implicit none
C$$$      logical lval
C$$$      integer mask
C$$$      character*(*) sname
C$$$      integer ix1(10),maski,nmask,t2n
C$$$      double precision struc(1),x2(10),x3(10),x4(10),x5(10)

C$$$      if (mask .le. 0) call rxi('lsets: mask must be >0 but given',mask)
C$$$      call spack(0,sname,struc,ix1,x2,x3,x4,x5)
C$$$      maski = mask
C$$$      nmask = ix1(1)
C$$$      t2n = 1
C$$$C ... Find lowest nonzero bit of mask, corresponding ix1
C$$$   10 continue
C$$$      if (mod(maski,2) .eq. 0) then
C$$$        t2n = t2n*2
C$$$        maski = maski/2
C$$$        ix1(1) = ix1(1)/2
C$$$        goto 10
C$$$      endif
C$$$      if (lval) then
C$$$        nmask = nmask + t2n*(1-mod(ix1(1),2))
C$$$      else
C$$$        nmask = nmask + t2n*(0-mod(ix1(1),2))
C$$$      endif
C$$$      call spack(1,sname,struc,nmask,x2,x3,x4,x5)

C$$$      end

