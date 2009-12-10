      subroutine suldau(nbas,sspec,ssite,nlibu,lmaxu,lldau)
C- Finds lda+U sites and counts number of blocks
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Co Outputs
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Co   nlibu :number of LDA+U blocks
Co   lmaxu :highest l for which a U block found, used as 
Co         :dimensioning parameter for U matrix
Cr Remarks
Cr
Cu Updates
Cu   27 Apr 05 (Lambrecht) first created
C------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nlibu,lmaxu,lldau(nbas),igetss,is,ib,l,lmxa,idu(4)
      double precision sspec(1), ssite(1)
C ... Local parameters
C     integer stdo,nglob
C     stdo = nglob('stdo')

      nlibu = 0
      lmaxu = 0
      do  ib = 1, nbas
        lldau(ib) = 0
        is = igetss('site spec',ib,ssite)
        call upack('spec lmxa idu',sspec,is,lmxa,idu,0,0)
        do  l = 0, min(lmxa,3)
          if (idu(l+1) .ne. 0) then
            if (lldau(ib) .eq. 0) lldau(ib) = nlibu+1
            nlibu = nlibu+1
            lmaxu = max(lmaxu,l)
          endif
        enddo
      enddo

      if (nlibu .ne. 0) then
        call info2(10,1,0,
     .    ' suldau:  %i U block(s)  lmaxu = %i',nlibu,lmaxu)
      endif

      end

