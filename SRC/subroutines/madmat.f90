!      subroutine madmat(nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg,dmad)
!        call madmat ( nbas , rv_a_opos , lat_awald , lat_alat , lat_vol , rv_a_odlv
!     .      , lat_nkd , rv_a_oqlv , lat_nkq , rv_a_omad )
subroutine madmat(dmad)
  use m_lmfinit,only: stdo,nbas,alat=>lat_alat
  use m_lattic,only: awald=>lat_awald,vol=>lat_vol,dlat=>rv_a_odlv &
       ,nkd=>lat_nkd,glat=>rv_a_oqlv, nkg=>lat_nkq,bas=>rv_a_opos
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
  !- Coefficients to Madelung matrix
  ! ----------------------------------------------------------------
  !i Inputs
  !i   nbas,bas,awald,alat,vol,dlat,nkd,glat,nkg
  !o Outputs
  !o   dmad
  !r Remarks
  !r   The potential from a basis of more than one atom is obtained by
  !r   superposing the potential from each sublattice.  Matrix element
  !r   dmad(i,j) is the (1/2)potential at position tau(i) from unit
  !r   charges on sublattice j, compensated by a uniform background of
  !r   density 1/vol.  (if tau(i)=0, the potential from the singularity
  !r   0 is not included.)
  !r
  !r   Call lattc to generate inputs awald,vol,dlat,nkd,glat,nkg.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  !      integer nbas,nkd,nkg
  !      double precision awald,alat,vol
  !      double precision bas(3,nbas),dlat(3,nkd),glat(3,nkg)
  double precision :: dmad(nbas,nbas)
  ! local variables
  integer :: ibas,jbas,iprint,i
  double precision :: tau(3)
  real(8):: pp(3)
  dmad=0.0d0
  ! --- Generate Madelung matrix ---
  do  101  ibas = 1, nbas
     do  10  jbas = ibas, nbas
        do  15  i = 1, 3
           tau(i) = bas(i,jbas)-bas(i,ibas)
15      enddo
        !call shortn(tau,tau,dlat,nkd)
        pp = matmul(transpose(qlat),tau)
        call shortn3_plat(pp)
        tau = matmul(plat,pp+nlatout(:,1))
        call strx00(tau,awald,alat,vol,glat,nkg,dlat,nkd, &
             dmad(ibas,jbas))
        dmad(jbas,ibas) = dmad(ibas,jbas)
10   enddo
101 enddo
  ! --- Printout ---
  if (iprint() < 90) return
  write(stdo,100) nbas
100 format(/' Madelung matrix: nbas=',i5)
  do  20  ibas = 1, nbas
     print 110, (dmad(ibas,jbas),jbas=1,nbas)
110  format(7f11.7)
20 enddo
end subroutine madmat

