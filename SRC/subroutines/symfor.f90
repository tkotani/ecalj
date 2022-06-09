subroutine symfor(nbas,mode,g,ng,istab,fwk,f)
  !- Symmetrize forces
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas:  number of atoms in basis
  !i   mode:  1 symmetrize by f(  ib ) = sum_g(ig) f(R(ib))
  !i          2 symmetrize by f(R(ib)) = sum_g(ig) f(  ib)
  !i   g,ng:  symmetry operations, and number
  !i   istab: site into which g,ag transforms site i
  !i   fwk:   work array of the same dimensions as f
  !o Outputs
  !o   f:  forces are symmetrized
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nbas,ng, istab(nbas,1),i,j,ib,ig,jb,mode
  double precision :: g(3,3,1),fwk(3,nbas),f(3,nbas)
  fwk=f/ng !!call dpcopy(f,fwk,1,3*nbas,1d0/ng)
  f=0d0
  if (mode == 1) then
     do  101  ig = 1, ng
        do  10  ib = 1, nbas
           jb = istab(ib,ig)
           do    i = 1, 3
              do    j = 1, 3
                 f(i,ib) = f(i,ib) + g(i,j,ig)*fwk(j,jb)
              enddo
           enddo
10      enddo
101  enddo
  else
     do  201  ig = 1, ng
        do  20  ib = 1, nbas
           jb = istab(ib,ig)
           do    i = 1, 3
              do    j = 1, 3
                 f(i,jb) = f(i,jb) + g(i,j,ig)*fwk(j,ib)
              enddo
           enddo
20      enddo
201  enddo
  endif
end subroutine symfor

