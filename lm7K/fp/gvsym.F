      subroutine gvsym(ng,gv,ips0,bgv,c,csym)
C- Symmetrize a function c, given in the form of a list
C ----------------------------------------------------------------------
Ci Inputs
Ci   gv,ng   Lattice vectors, and number
Ci   ips0    pointer to first vector in star of this vector; see sgvsym
Ci   bgv     phase factor sum; see sgvsym
Ci   c       unsymmetrized function
Co Outputs:
Co   csym    symmetrized function
Cr Remarks:
Cu    7 Sep 98 Adapted from nfp gvsym.f
C ----------------------------------------------------------------------
C     implicit none
      integer ng,ips0(ng)
      double precision gv(ng,3)
      double complex bgv(ng),c(ng),csym(ng)
      integer i,j,i0,kstar,ipr,iprint,nstar

C ... Sum up coefficients for first vector in each star
      do  10  i = 1, ng
        csym(i) = (0d0,0d0)
   10 continue
      do  12  i = 1, ng
        j = ips0(i)
        csym(j) = csym(j) + bgv(i)*c(i)
   12 continue

C ... Normalize
      do  20  i0 = 1, ng
        if (ips0(i0) .eq. i0) then
          kstar = 0
          do  22  i = i0, ng
            if (ips0(i) .eq. i0) kstar = kstar+1
   22     continue
          csym(i0) = csym(i0)/kstar
        endif
   20 continue

C ... Make all the coefficients
      do  30  i = 1, ng
        j = ips0(i)
        csym(i) = csym(j)*dconjg(bgv(i))
   30 continue

C ... Printout
      ipr = iprint()
      if (ipr .le. 55) return
      print 255
      nstar = 0
      do  40  i0 = 1, ng
        if (ips0(i0) .eq. i0) then
          nstar = nstar+1
          if (ipr .ge. 60) print *, ' '
          do  44  i = i0, ng
            if (ips0(i) .eq. i0) then
              if (i .eq. i0) then
                print 251, nstar,i,gv(i,1),gv(i,2),gv(i,3),
     .             c(i),csym(i)
              else
                if (ipr .ge. 60) 
     .             print 252, i,gv(i,1),gv(i,2),gv(i,3),
     .             c(i),csym(i)
              endif
  251         format(i4,i5,3f6.1,2f12.8,1x,2f12.8)
  252         format(4x,i5,3f6.1,2f12.8,1x,2f12.8)
  255         format(/' star  ig',8x,'recip',17x,'c_in',20x,'c_sym')
            endif
   44     continue
        endif
   40 continue

      end

