subroutine fpiint(nx,np,nxp,x,w)
  ! akao. Maybe it will be better to replace this with spherical design
  ! s http://www2.research.att.com/~njas/sphdesigns/

  !- Points and weights for integration on a sphere surface
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nx    :number of points in polar angle (Legendre integration)
  !i          Use nx<0 for special points; see Remarks.
  !i   np    :number of points in phi (uniform points)
  !i         :np=0 => np depends on nx making dphi approx constant.
  !i          for nx<0, np is not used.
  !o Outputs
  !o   nxp   :total number of number of points in quadrature
  !o   x     :cartesian coordinates of points on unit sphere
  !o   w     :weights corresponding to x
  !r Remarks
  !r   fpiint generates a mesh of points on a unit sphere for angular
  !r   integration, either using a set of special generated from the
  !r   Platonic solids, or by integrating the polar angle with Legendre
  !r   gaussian quadrature and the azimuthal angle with a set of evenly
  !r   spaced points on a circle.
  !r   For special points, invoke fpiint with one of the following:
  !     nx= -4 integrates any ylm<=2 exactly (tetrahedron)
  !     nx= -6 integrates any ylm<=3 exactly (faces of cube)
  !     nx= -8 integrates any ylm<=3 exactly (cube)
  !        -12 integrates any ylm<=5 exactly (icosahedron)
  !        -20 integrates any ylm<=5 exactly (faces of icosahedron)
  !        -30 integrates any ylm<=5 exactly (sides of icosahedron)
  !        -60 integrates any ylm<=5 exactly (buckeyball)
  !        -32 integrates any ylm<=9 exactly  (combination of 12,20)
  !        -62 integrates any ylm<=11 exactly (combination of 12,20,30)
  !        -92 integrates any ylm<=11 exactly (combination of 12,20,60)
  !       -122 integrates any ylm<=15 exactly (combination of 12,20,30,60)
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nx,np,nxp
  double precision :: x(3,*),w(*),fpi
  integer :: i,iprint

  if (nx >= 0) then
     call nintsp(nx,np,nxp,x,w)
  else
     fpi = 16*datan(1d0)
     nxp = -nx
     if (nx == -32) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        do   i = 1, 12
           w(i) = 5d0*fpi/(14*12)
        enddo
        do    i = 13, 32
           w(i) = 9d0*fpi/(14*20)
        enddo
     elseif (nx == -62) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        call platsl(x(1,33),30)
        do  i = 1, 12
           w(i) = 125d0*fpi/(14*12*33)
        enddo
        do   i = 13, 32
           w(i) = 81d0*fpi/(14*20*33)
        enddo
        do   i = 33, 62
           w(i) = 256d0*fpi/(14*30*33)
        enddo
     elseif (nx == -92) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        call platsl(x(1,33),60)
        do    i = 1, 12
           w(i) = 1/12.34817490904537d0*fpi/12
        enddo
        do   i = 13, 32
           w(i) = 2.986997567806883d0/12.34817490904537d0*fpi/20
        enddo
        do    i = 33, 92
           w(i) = 8.361177341238484d0/12.34817490904537d0*fpi/60
        enddo
     elseif (nx == -122) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        call platsl(x(1,33),30)
        call platsl(x(1,63),60)
        do   i = 1, 12
           w(i) = (0.0939463041645901d0)*fpi/12
        enddo
        do    i = 13, 32
           w(i) = (0.2373458837681504d0)*fpi/20
        enddo
        do    i = 33, 92
           w(i) = (0.0378880378880377d0)*fpi/30
        enddo
        do    i = 63, 122
           w(i) = (0.6308197741792218d0)*fpi/60
        enddo
     else
        call platsl(x,nxp)
        do   i = 1, nxp
           w(i) = fpi/nxp
        enddo
     endif
  endif
  ! --- Printout ---
  if (iprint() < 80) return
  print '(/'' fpiint:'',i5,'' points generated:'')', nxp
  do  i = 1, nxp
     print 333, i, x(1,i), x(2,i), x(3,i), w(i)
  enddo
333 format(i3,4f20.15)
end subroutine fpiint

