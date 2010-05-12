      subroutine rsort(nr,rsq,rc1,rc2,lsort,n1,n2,idx,wk)
C- Partition a mesh of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsq,nr  : vector of points r**2, and number of points.
Ci   rc1,rc2 : boundaries of bins in terms of r**2
Ci   lsort   : if .true. then rsq must be ordered as rsq(i-1).le.rsq(i)
Co Outputs
Co   n1,n2   : indecies corresponding to boundaries of the bins
Co   idx     : permutation indecies; not used if lsort = .true.
Co   wk      : an arry containing permutted rsq; not used if
Co           : lsort = .true.
Co   lsort   : if .false. but rsq is nevertheless ordered then lsort
Co           : is set to .true.
Cr Remarks
Cr   Points are partitioned into three bins:
Cr     r < rc1       <=> rsq(1..n1-1)
Cr     rc1 < r < rc2 <=> rsq(n1..n2-1)
Cr     rc2 < r       <=> rsq(n2..nr)
Cr   If rsq are not ordered (lsort = .false.) then wk plays the role
Cr   of rsq, whereas the info about the correspondence between points
Cr   in rsq and wk is held in idx: rsq(i) = wk(idx(i)).
Cr   idx2 is a work array
C ----------------------------------------------------------------------
C     implicit none
C Passed variables
      integer nr,n1,n2,idx(nr)
      double precision rsq(nr),wk(nr),rc1,rc2
      logical lsort
C Local variables
      integer ir,j,k,n0
C Work auto-arrays
      integer idx2(nr)


      n0 = 0
      n1 = 0
      n2 = nr+1
C ... Case points already sorted.  Find n1,n2.
      if (lsort) then

c ... verify that rsq is indeed ordered
        do  ir = 2, nr
          if (rsq(ir) .lt. rsq(ir-1)) 
     .    call rx(' rsort: rsq is not ordered. Check input switches.')
        enddo

        n1 = 1
        if (nr .eq. 1) then
          if (rsq(1) .lt. rc1) n1 = 2
          if (rsq(1) .gt. rc2) n2 = 1
        else
          if (rsq(1) .ge. rc1) then
            n1 = 1
          else
            call huntx(rsq,nr,rc1,0,n1)
            n1 = n1+1
          endif
          if (rsq(nr) .le. rc2) then
            n2 = nr+1
          else
            n2 = nr
            call huntx(rsq,nr,rc2,0,n2)
            n2 = n2+1
          endif
        endif

C ... Case points not sorted (iwk, wk required now)
      else
C     On output, lsort is true if points already sorted.
        lsort = .true.
        do  12  ir = 1, nr
C     n1 is offset to block rc1<r<rc2,  n2 offset to block r>rc2
C     idx is a map of original list, separating into the three groups
C     wk is a table of r**2 for permuted list of points
          if (rsq(ir) .lt. rc2) then
            if (rsq(ir) .lt. rc1) then
              n0 = n0+1
              wk(n0) = rsq(ir)
              idx(ir) = n0
            else
              n1 = n1+1
              idx(ir) = n1
              idx2(n1) = ir
            endif
          else
            n2 = n2-1
            wk(n2) = rsq(ir)
            idx(ir) = n2
          endif
          if (ir .eq. 1) goto 12
          if (rsq(ir) .lt. rsq(ir-1)) lsort = .false.
   12   continue
C ... Now we can poke wk for the n1 intermediate points
        if (.not. lsort) then
          do  14  j = 1, n1
            k = idx2(j)
            idx(k) = idx(k)+n0
            wk(n0+j) = rsq(k)
   14     continue
        endif
        n1 = n0+1
      endif

      end

