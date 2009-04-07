      subroutine sp2cls(sname,sspec,ics,nelt,is1,is2,oarr)
C- Allocate and extract one element of a structure for a range of species
C ----------------------------------------------------------------------
Ci Inputs
Ci   sname :a species-based structure and the element to extract, e.g.
Ci         :'spec idmod'
Ci   sspec :array holding structure sname
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   nelt  :number of row entries this element has in structure
Ci   is1   :extract elements starting at species is1
Ci   is2   :extract elements ending at species is2
Co Outputs
Co   oarr  :Array oarr is allocated and loaded,
Co         :corresponding to arr(1:nelt,is1:is2)
Cr Remarks
Cr   
Cu Updates
Cu   26 Nov 02 
C ----------------------------------------------------------------------
C     implicit none
      integer is2,oarr,ics(is2)
      double precision sspec(1)
C ... Heap
      integer w(1)
      common /w/ w
C ... Local variables
      character*(*) sname
      integer cast,nelt,i,ic,is,offs,ssize,is1,offi,offa,nelti
      call spack(2,sname,sspec,1,offi,cast,nelti,i)
      ssize = int(sspec(1))
      oarr = 0
      if (cast .eq. 2) call defi(oarr,(is2-is1+1)*nelt)
      if (cast .eq. 4) call defrr(oarr,(is2-is1+1)*nelt)
      call rxx(oarr.eq.0,'sp2cls: bad cast')
      do  10  ic = is1, is2
        if (ics(1) .eq. 0) then
          is = ic
        else
          is = ics(ic)
        endif
        offs = offi + ssize*(is-1)
        offa = 1 + (ic-is1)*nelt
        call peekx(nelt,offa,offs,1,cast,w(oarr),i,i,i,i,sspec)
   10 continue
C      i = (is2-is1+1)*nelt
C      if (cast.eq.2) call awrit2('%n:1i',' ',100,6,i,w(oarr))
C      if (cast.eq.4) call awrit2('%n:1d',' ',100,6,i,w(oarr))
      end

