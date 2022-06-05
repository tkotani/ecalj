subroutine suylg(ltop,alat,ng,gv,g,g2,yl)
  use m_ropyln,only: ropyln
  !- Set up vectors g, g2, yl from list of vectors gv
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ltop  :l-cutoff for YL
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !o Outputs
  !o   g     :gv scaled by (2 pi / alat)
  !o   g2    :square of g
  !o   yl    :YL(g)
  !u Updates
  !u   30 May 00 adapted from nfp su_ylg
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ltop,ng
  double precision :: alat,gv(ng,3),g(ng,3),yl(ng,1),g2(ng)
  ! ... Local parameters
  integer :: i
  double precision :: pi,tpiba

  ! ... Make (2*pi/alat)*gv in g
  pi = 4d0*datan(1d0)
  tpiba = 2d0*pi/alat
  do  i = 1, ng
     g(i,1) = tpiba*gv(i,1)
     g(i,2) = tpiba*gv(i,2)
     g(i,3) = tpiba*gv(i,3)
  enddo

  ! ... Make the yl's and g2
  call ropyln(ng,g(1,1),g(1,2),g(1,3),ltop,ng,yl,g2)

end subroutine suylg

