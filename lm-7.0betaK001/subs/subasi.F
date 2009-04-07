      subroutine subasi(sctrl,sspec,sham)
C- Read some parameters defining hamiltonian indpendent of potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct containing parameters governing program flow
Ci   sspec :struct containing species-specific information
Ci   sham  :struct containing information about the hamiltonian
Co Outputs
Co   Parameters defining basis are set in sspec
Co   Global parameters nkaph and mxorb are set
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr  *It generates and packs a table of hamiltonian offsets offH,
Cr   orbital permutation indices oindxo.
Cr
Cr  *For the ASA 2nd generation LMTO:
Cr   Extract order of potential function from gfopts
Cr   Transform pp's to alpha representation
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sctrl(1),sspec(1),sham(1)
C ... Local parameters
      double precision dglob,xx
      integer fopna,ifi,igetss,is,lmxax,nbasp,nkaph,nkmax,nlmax,
     .  nspec,lfp

C     call upack('ctrl lfp nbasp nspec',sctrl,lfp,nbasp,nspec,0,0)

C --- Get the maximum L-cutoff ---
C      Moved to rdccat.f
C      lmxax = -1
C      do  is = 1, nspec
C        lmxax = max(lmxax,igetss('spec lmxa',is,sspec))
C      enddo
C      nlmax = (lmxax+1)**2
C      xx = dglob('nkaph',1d0,1)
C      xx = dglob('mxorb',dble(nlmax),1)

      end

