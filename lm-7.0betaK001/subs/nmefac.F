      subroutine nmefac(nmto,kaph,f)
C- NMTO energy factors
C ----------------------------------------------------------------------
Ci Inputs
Ci   nmto  :number of NMTO energies
Ci   kaph  :NMTO kinetic energy
Co Outputs
Co   f     :f^N_i = prod_{i'=0, \neq i}^{N} (kaph_i-kaph_i')^-1
Co         :See, OKA, Mt. St. Odile book, Eq. 91
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nmto
      double precision kaph(nmto),f(nmto)
C ... Local parameters
      integer ik,jk

C --- Factors f^N_i ---
      do  ik = 1, nmto
        f(ik) = 1
        do  jk = 1, nmto
          if (ik .ne. jk) f(ik) = f(ik)/(kaph(ik)-kaph(jk))
        enddo
      enddo

      end


