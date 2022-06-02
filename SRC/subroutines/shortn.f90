subroutine shortn(p,p1,dlat,nkd)
  !- Get p1 = shortest vector such that p1-p is a lattice vector.
  ! ----------------------------------------------------------------
  !i Inputs
  !i   p     :vector to shorten
  !i   dlat  :lattice vectors, sorted by increasing length
  !i   nkd   :number of dlat
  !o Outputs
  !o   p1    :shortened vector
  !r Remarks
  !r   A slightly skewed norm is used to make result unique.
  !r   The first vector in the list must be the zero vector.
  ! ----------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nkd
  double precision :: p(3),p1(3),dlat(3,nkd)
  ! ... Local parameters
  integer :: irep,k0,k
  double precision :: anrm2,x,y,z,p2,critk0,dd,crit
  anrm2(x,y,z) = x*x*1.00001d0 + y*y*1.00002d0 + z*z*1.00003d0 &
       -x*0.000004d0 - y*0.000003d0 - z*0.000002d0

  p1(1) = p(1)
  p1(2) = p(2)
  p1(3) = p(3)
  do  88  irep = 1, 20
     p2 = anrm2(p1(1),p1(2),p1(3))
     k0 = 1
     critk0 = p2
     do  52  k = 1, nkd
        dd = dlat(1,k)**2 + dlat(2,k)**2 + dlat(3,k)**2
        !         Only valid if dlat is sorted:
        if (dd > p2*4d0) goto 53
        crit = anrm2(p1(1)+dlat(1,k),p1(2)+dlat(2,k),p1(3)+dlat(3,k))
        if (crit < critk0) then
           k0 = k
           critk0 = crit
        endif
52   enddo
53   if (k0 == 1) return
     p1(1) = p1(1) + dlat(1,k0)
     p1(2) = p1(2) + dlat(2,k0)
     p1(3) = p1(3) + dlat(3,k0)
88 enddo
  call rx('shortn: shortest vector not found')
end subroutine shortn

