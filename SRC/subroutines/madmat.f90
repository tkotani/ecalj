c      subroutine madmat(nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg,dmad)
c        call madmat ( nbas , rv_a_opos , lat_awald , lat_alat , lat_vol , rv_a_odlv
c     .      , lat_nkd , rv_a_oqlv , lat_nkq , rv_a_omad )
      subroutine madmat(dmad)
      use m_lmfinit,only: stdo,nbas,alat=>lat_alat
      use m_lattic,only: awald=>lat_awald,vol=>lat_vol,dlat=>rv_a_odlv
     &  ,nkd=>lat_nkd,glat=>rv_a_oqlv, nkg=>lat_nkq,bas=>rv_a_opos
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
c      integer nbas,nkd,nkg
c      double precision awald,alat,vol
c      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg)
      double precision dmad(nbas,nbas)
C local variables
      integer ibas,jbas,iprint,i
      double precision tau(3)
      dmad=0.0d0
C --- Generate Madelung matrix ---
      do  101  ibas = 1, nbas
      do  10  jbas = ibas, nbas
        do  15  i = 1, 3
          tau(i) = bas(i,jbas)-bas(i,ibas)
   15   continue
        call shortn(tau,tau,dlat,nkd)
        call strx00(tau,awald,alat,vol,glat,nkg,dlat,nkd,
     .    dmad(ibas,jbas))
        dmad(jbas,ibas) = dmad(ibas,jbas)
   10 continue
 101  continue
C --- Printout ---
      if (iprint() .lt. 90) return
      write(stdo,100) nbas
  100 format(/' Madelung matrix: nbas=',i5)
      do  20  ibas = 1, nbas
        print 110, (dmad(ibas,jbas),jbas=1,nbas)
  110   format(7f11.7)
   20 continue
      end

