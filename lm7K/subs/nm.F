      integer function nm(nl,nl2,l)
C- returns number of m quantum numbers for looping
C ----------------------------------------------------------------
Ci Inputs
Ci   nl,nl2,l
Co Outputs
Co   m=0 if nl=nl2; m=l if nl=nl2**2; stops with error otherwise
Cr Remarks
Cr   Routine is used to determine number of m quantum numbers for a
Cr   specified l; provides facility to for a routine to refer
Cr   to arrays of dimension (nl,*) or of (nlm,*) without routine
Cr   having to specify.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nl2,l
C Local parameters
      nm = 0
      if (nl .eq. nl2) return
      if (nl**2 .ne. nl2)
     .  call fexit(-1,119, 'NM: mismatch between nl and nl2',0)
      nm = l
      end

