      subroutine gtbsl2(l1,lmxh,eh,rsmh,l2)
C- Returns the highest l with rsm,e common to those of a given l
C ----------------------------------------------------------------------
Ci Inputs
Ci   l1    :current l
Ci   lmxh :basis l-cutoff
Ci   eh    :energy of smoothed Hankel
Ci   rsmh  :smoothing radius of smoothed hankel
Co Outputs
Co   l2    :large l for which eh and rsmh l1..l2 are in common
Cr Remarks
Cr   Routine used group functions, strux in blocks.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer l1,l2,lmxh
      double precision rsmh(0:lmxh),eh(0:lmxh)
      double precision e,rsm

      e = eh(l1)
      rsm = rsmh(l1)
      l2 = l1
   10 if (l2 .ge. lmxh) return
      if (rsmh(l2+1).ne.rsm .or. eh(l2+1).ne.e) return
      l2 = l2+1
      goto 10

      end

