      subroutine madmat(nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg,dmad)
      use m_lmfinit,only: stdo


C- Coefficients to Madelung matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg
Co Outputs
Co   dmad
Cr Remarks
Cr   The potential from a basis of more than one atom is obtained by
Cr   superposing the potential from each sublattice.  Matrix element
Cr   dmad(i,j) is the (1/2)potential at position tau(i) from unit
Cr   charges on sublattice j, compensated by a uniform background of
Cr   density 1/vol.  (if tau(i)=0, the potential from the singularity
Cr   0 is not included.)
Cr
Cr   Call lattc to generate inputs awald,vol,dlat,nkd,glat,nkg.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nkd,nkg
      double precision awald,alat,vol
      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg)
      double precision dmad(nbas,nbas)
C local variables
      integer ibas,jbas,iprint,i
      double precision tau(3)

C --- Generate Madelung matrix ---
      do  10  ibas = 1, nbas
      do  10  jbas = ibas, nbas
        do  15  i = 1, 3
          tau(i) = bas(i,jbas)-bas(i,ibas)
   15   continue
        call shortn(tau,tau,dlat,nkd)
        call strx00(tau,awald,alat,vol,glat,nkg,dlat,nkd,
     .    dmad(ibas,jbas))
        dmad(jbas,ibas) = dmad(ibas,jbas)
   10 continue

C --- Printout ---
      if (iprint() .lt. 90) return
Changenglob      stdo = nglob('stdo')
c      stdo = globalvariables%stdo
      write(stdo,100) nbas
  100 format(/' Madelung matrix: nbas=',i5)
      do  20  ibas = 1, nbas
        print 110, (dmad(ibas,jbas),jbas=1,nbas)
  110   format(7f11.7)
   20 continue
      end

