subroutine projql(qlat)
  use m_lmfinit,only: stdo


  !- Project 3D reciprocal lattice vecs on surface (2D recip lattice vecs)
  ! ----------------------------------------------------------------------
  !i Input qlat(3,3)  3D reciprocal lattice vectors
  !o Output qlat(3,3) 2D reciprocal lattice vectors
  !r if bi i=1,3 are 3D vectors
  !r    Bi i=1,2 are 2D vectors
  !r    Bi(u)=bi(u)-(Bi.B3) B3(u)/(B3.B3) for u=x,y,z
  ! ----------------------------------------------------------------------
  implicit none
  double precision :: qlat(3,3),b3b3,bib3
  integer :: i,j,iprint
  b3b3 = qlat(1,3)**2 + qlat(2,3)**2 + qlat(3,3)**2
  do  20  i = 1, 2
     bib3 = qlat(1,i)*qlat(1,3) + qlat(2,i)*qlat(2,3) &
          + qlat(3,i)*qlat(3,3)
     do  21  j = 1, 3
        qlat(j,i) = qlat(j,i) - qlat(j,3)*bib3/b3b3
21   enddo
20 enddo
  ! ... Ignore qlat(3,.) by making it zero
  !      do  22  j = 1, 3
  !   22 qlat(j,3)=0d0
  if (iprint() >= 20) then
     ! hangenglob        stdo = nglob('stdo')
     !        stdo = globalvariables%stdo
     write(stdo,*) ' 2D-reciprocal vectors:'
     do  23  i = 1, 2
        write(stdo,222) (qlat(j,i),j=1,3)
23   enddo
  endif
222 format(3f12.6)
end subroutine projql

