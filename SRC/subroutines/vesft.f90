subroutine vesft(ng,gv,kv,cv,smrho,smpot,ssum)
  !Make electrostatic potential of density given in recip space
  use m_lmfinit,only:alat=>lat_alat
  use m_lattic,only: vol=>lat_vol
  use m_supot,only: k1,k2,k3
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   cv    :work array holding smrho and smpot in glist form
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
  !i   smrho :FT of smooth density on uniform mesh
  !o Outputs
  !o   smpot :FT of smooth electrostatic potential
  !o   sum   :integral pot*density
  implicit none
  integer :: ng,i, kv(ng,3)
  real(8):: gv(ng,3) , ssum, tpiba,g2
  complex(8):: smrho(k1,k2,k3),smpot(k1,k2,k3),cv(ng),ccc(ng)
  real(8),parameter:: pi = 4d0*datan(1d0), pi8  = 8d0*pi
  call tcn('vesft')
  tpiba=2d0*pi/alat
  call gvgetf(ng,1,kv,k1,k2,k3,smrho,cv) ! ... Gather density coefficients
  ! ... smpot(G) = 8 pi /G**2 smrho(G)
  ccc(1)=0d0
  ccc(2:ng) = [((pi8/(tpiba**2*sum(gv(i,:)**2)))*cv(i),i=2,ng)]
  ssum = vol*sum(ccc(2:ng)*cv(2:ng))
  call gvputf(ng,1,kv,k1,k2,k3,ccc,smpot)! ... Scatter smooth potential into array
  call tcx('vesft')
end subroutine vesft


