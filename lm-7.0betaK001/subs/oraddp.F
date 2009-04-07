      subroutine oraddp(mode,iax,is1,is2,ib1,ib2,nl2,iwk,iprm)
C- Alters iprm to include orbitals connected to sites ib1..ib2.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode   :1, initially exclude all orbitals i associated with
Ci           ib1..ib2 by setting iprm(i) to negative definite.
Ci   iax    :neighbor table containing pair information (pairc.f)
Ci   is1,is2:range of orbital orbital pairs to consider
Ci   ib1,ib2:range of allowed site indices
Ci   nl2    :spacing between starting offsets in iprm table
Ci   iwk    :work array of size at least as large as max ib in iax.
Ci   iprm   :permutation indices ordering orbitals
Ci           iprm(i) => this orbital is excluded from basis
Co Outputs
Ci   iprm   :some set of iprm(i) which are connected to sites
Ci           ib1..ib2 set positive definite.
Cr Remarks
Cr   Orbital pairs is=is1..is2 connect a set of sites {ib}, which is
Cr   the collection of all distinct iax(2,is).  oraddp sets the sign
Cr   of iprm for those orbitals belonging the set of {ib} sites.
Cr   The sign of iprm(i) for orbitals i belonging to any site ib
Cr   for which   ib1<=ib<=ib2   is set positive.  Thus, the sign
Cr   of iprm flags which orbitals are connected to the range sites
Cr   ib1..ib2.  To exclude orbitals not in this list, set initially
Cr   all iprm(i)<0.  Obtaining orbitals for groups of ranges
Cr   ib1..ib2,ib1'..ib2',... can be obtained by successive calls.
Cu Updates
Cu   13 Dec 01 iwk need not be dimensioned larger than ib2.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax,mode,ib1,ib2,is1,is2,nl2,iprm(*),iwk(ib2)
      parameter (niax=10)
      integer iax(niax,is2)
C ... Local parameters
      integer ib,offi,lm,ibmin,ibmax,is

      if (mode .eq. 1) then
        do  5  ib = ib1, ib2
        offi = nl2*(ib-1)
        do  5  lm = offi+1, offi+nl2
          iprm(lm) = -iabs(iprm(lm))
    5   continue
      endif

      if (is2 .lt. is1) return

      ibmin = iax(2,is1)
      ibmax = iax(2,is1)
      do  10  is = is1, is2
        ib = iax(2,is)
        ibmin = min(ibmin,iax(2,is))
        ibmax = max(ibmax,iax(2,is))
   10 continue
      ibmax = min(ibmax,ib2)
      call ivset(iwk,ibmin,ibmax,0)
      do  11  is = is1, is2
        ib = iax(2,is)
        if (ib .le. ibmax) iwk(ib) = 1
   11 continue

      do  20  ib = ib1, ib2
        if (iwk(ib) .ne. 0) then
          offi = nl2*(ib-1)
          do  22  lm = offi+1, offi+nl2
            iprm(lm) = iabs(iprm(lm))
   22     continue
        endif
   20 continue

      end

