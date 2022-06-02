subroutine vesft(ng,gv,kv,cv,k1,k2,k3,smrho,smpot,sum) !slat,
  use m_struc_def           !Cgetarg
  use m_lmfinit,only:lat_alat
  use m_lattic,only: lat_vol


  !- Make electrostatic potential of density given in recip space
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   slat  :struct containing information about the lattice
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   cv    :work array holding smrho and smpot in glist form
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
  !i   smrho :FT of smooth density on uniform mesh
  !o Outputs
  !o   smpot :FT of smooth electrostatic potential
  !o   sum   :integral pot*density
  !r Remarks
  !u Updates
  !u   22 Apr 00 Adapted from nfp ves_ft.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng
  integer :: kv(ng,3),k1,k2,k3
  real(8):: gv(ng,3) , sum
  !      type(s_lat)::slat
  double complex smrho(k1,k2,k3),smpot(k1,k2,k3),cv(ng)
  ! ... Local parameters
  integer :: i
  double precision :: pi,pi8,alat,vol,tpiba,g2
  double complex ccc
  call tcn('vesft')
  alat=lat_alat
  vol=lat_vol
  pi   = 4d0*datan(1d0)
  pi8  = 8*pi
  tpiba=2*pi/alat

  ! ... Gather density coefficients
  call gvgetf(ng,1,kv,k1,k2,k3,smrho,cv)

  ! ... smpot(G) = 8 pi /G**2 smrho(G)
  !     call info2(30,1,0,' vesft:  smooth density coeff to (0,0,0) '//
  !    .  'is %;12,6D',cv(1),0)
  ! cccccccccccc
  !      print *,'ng=',ng
  !      stop 'xxxxxxxxxx vesft'
  ! cccccccccccc
  sum = 0d0
  cv(1) = 0
  do  i = 2, ng
     g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
     ccc = (pi8/g2)*cv(i)
     sum = sum + dconjg(cv(i))*ccc
     cv(i) = ccc
  enddo
  sum = vol*sum

  ! ... Scatter smooth potential into array
  call gvputf(ng,1,kv,k1,k2,k3,cv,smpot)
  call tcx('vesft')
end subroutine vesft


