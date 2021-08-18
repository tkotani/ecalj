      subroutine symfor(nbas,mode,g,ng,istab,fwk,f)
C- Symmetrize forces
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas:  number of atoms in basis
Ci   mode:  1 symmetrize by f(  ib ) = sum_g(ig) f(R(ib))
Ci          2 symmetrize by f(R(ib)) = sum_g(ig) f(  ib)
Ci   g,ng:  symmetry operations, and number
Ci   istab: site into which g,ag transforms site i
Ci   fwk:   work array of the same dimensions as f
Co Outputs
Co   f:  forces are symmetrized
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,ng
      integer istab(nbas,1)
      double precision g(3,3,1),fwk(3,nbas),f(3,nbas)
C Local variables
      integer i,j,ib,ig,jb,mode

C     call prmx('f-in',f,3,3,nbas)
      call dpcopy(f,fwk,1,3*nbas,1d0/ng)
      call dpzero(f,3*nbas)

      if (mode .eq. 1) then
        do  101  ig = 1, ng
        do  10  ib = 1, nbas
          jb = istab(ib,ig)
          do    i = 1, 3
          do    j = 1, 3
             f(i,ib) = f(i,ib) + g(i,j,ig)*fwk(j,jb)
          enddo
          enddo
   10   continue
 101    continue
      else
        do  201  ig = 1, ng
        do  20  ib = 1, nbas
          jb = istab(ib,ig)
          do    i = 1, 3
          do    j = 1, 3
             f(i,jb) = f(i,jb) + g(i,j,ig)*fwk(j,ib)
          enddo
          enddo
   20   continue
 201    continue
      endif

C     call prmx('f-sym',f,3,3,nbas)

      end

