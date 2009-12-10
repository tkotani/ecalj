      subroutine slinz(volwgt,ec,emin,emax,dosi,nr)
C- Adds to number-of-states for one tetrahedron
C ----------------------------------------------------------------
Ci Inputs
Ci   volwgt:weight on tetrahedron
Ci   ec    :energies at corners of tetrahedron
Ci   emin, emax, energy window
Ci   nr    :number of bins + 1
Co Outputs
Co   dosi(k), integrated density in kth bin from tethdn.
Co   ec:   sorted on output
Cr Remarks
Cr
Cu Updates
Cu   19 Jun 04 Put in guard against integer overflow
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nr
      double precision ec(4),dosi(nr),volwgt,emin,emax
C Local parameters
      integer i,i01,i02,i03,i04,i1,i2,i3,i4,j
      double precision c0,c1,c2,c3,cc,de,e,e1,e2,e3,e4,x
      double precision fuzz,mxmin,xx
      parameter (fuzz=1d-8)
C     Guard against overflow, underflow. Make -1 < i0[1..4] < nr+1
      mxmin(xx) = min(max(xx,-2*de),nr*de)/de+1.9999999d0

C --- Sort the ec ---
      do  3  i = 1, 3
        do  3  j = 1, 4-i
        if (ec(j) .gt. ec(j+1)) then
          e = ec(j)
          ec(j) = ec(j+1)
          ec(j+1) = e
        endif
    3 continue
      e1 = ec(1)
      e2 = ec(2)
      e3 = ec(3)
      e4 = ec(4)
      if (e4 .lt. emin+0*fuzz) then
        i4 = 1
        go to 26
      endif
      de = (emax-emin)/(nr-1)
C      i01 = (e1-emin)/de + 1.9999999d0
C      i02 = (e2-emin)/de + 1.9999999d0
C      i03 = (e3-emin)/de + 1.9999999d0
C      i04 = (e4-emin)/de + 1.9999999d0
C
C      print *,'slinz',max(i01,1), int(mxmin(e1-emin)),
C     .                max(i01,1)-int(mxmin(e1-emin))
C      print *,'slinz',max(i02,1), int(mxmin(e2-emin)),
C     .                max(i02,1)-int(mxmin(e2-emin))
C      print *,'slinz',max(i03,1), int(mxmin(e3-emin)),
C     .                max(i03,1)-int(mxmin(e3-emin))
C      print *,'slinz',max(i04,1), int(mxmin(e4-emin)),
C     .                max(i04,1)-int(mxmin(e4-emin))
      i01 = mxmin(e1-emin)
      i02 = mxmin(e2-emin)
      i03 = mxmin(e3-emin)
      i04 = mxmin(e4-emin)

C --------------------------------
      i1 = max0(i01,1)
      i2 = min0(i02-1,nr)
      if (i1 .le. i2) then
        cc = volwgt/((e2-e1)*(e3-e1)*(e4-e1))
        do  20  i = i1, i2
          x = emin - e1 + (i-1)*de
          dosi(i) = dosi(i) + cc*x**3
   20   continue
      endif
      i2 = max0(i02,1)
      i3 = min0(i03-1,nr)
      if (i2 .le. i3) then
        c3 = volwgt*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
        c2 = volwgt*3d0/((e3-e1)*(e4-e1))
        c1 = c2*(e2-e1)
        c0 = c1*(e2-e1)/3d0
        do  21  i = i2, i3
          x = emin - e2 + (i-1)*de
          dosi(i) = dosi(i) + c0 + x*(c1 + x*(c2 + x*c3))
   21   continue
      endif
      i3 = max0(i03,1)
      i4 = min0(i04-1,nr)
      if (i3 .le. i4) then
        cc = volwgt/((e3-e4)*(e2-e4)*(e1-e4))
        do  22  i = i3, i4
          x = emin - e4 + (i-1)*de
          dosi(i) = dosi(i) + volwgt - cc*x**3
   22   continue
      endif
      i4 = max0(i04,1)
   26 continue
      do  25  i = i4, nr
   25 dosi(i) = dosi(i) + volwgt
      end

