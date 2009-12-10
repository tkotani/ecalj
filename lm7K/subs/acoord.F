      subroutine acoord(ib,it,alat,plat,pos,iax,ntab,coord) 
C- Find absolute coordinates of atom it in cluster ib 
C ----------------------------------------------------------------
Ci Inputs:
Ci   ib    :head of the cluster 
Ci   it    :atom index in cluster ib as specified in the iax table  
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat 
Ci   pos   :site positions 
Co   iax   :neighbor table containing information about each pair (pairc.f)
Co   iax(1)   :site index to basis atoms ib=source;
Co             all pairs with common ib are contiguous
Co   iax(2)   :site index to jb=field of each pair
Co   iax(3..5):multiples of plat added the difference in site positions
Co             that connect the pair.
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: no matching pairs outside (ib1..ib2) can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(8)   :not used
Co   iax(9)   :not used 
Co   iax(10)  :effective site index; see siteid.f
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Co Outputs:
Co   coord :Cartesian coordinates of atom it in cluster ib
Cr Remarks
Cu Updates
Cu   15 May 08 (S. Lozovoi) First written
C ----------------------------------------------------------------
C     implicit none
      integer ib,it,niax 
      parameter (niax=10)
      integer iax(niax,*),ntab(*)
      double precision alat,plat(3,3),pos(3,*),coord(3)
C Local variables
      integer offR,it1,ix
      double precision d,drr2,dr(3)

      offR = ntab(ib)
      it1 = offR + it 
      if (iax(1,it1) .ne. ib) 
     .  call rxi('acoord: mismatch in iax table, record # ',it1)

      d = drr2(plat,pos(1,iax(1,it1)),pos(1,iax(2,it1)),
     . iax(3,it1),iax(4,it1),iax(5,it1),dr)

      do ix = 1, 3
        coord(ix) = alat*(pos(ix,ib) + dr(ix))
      enddo
      
      end

