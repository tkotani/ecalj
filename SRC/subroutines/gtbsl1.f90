!> Marks blocks of contiguous l for which rsm and e are unchanged
subroutine gtbsl8(norb,ltab,ktab,rsmh,eh, ntab,blks)
  implicit none
  integer :: norb,ltab(norb),ktab(norb), ntab(norb),blks(norb)
  integer,parameter::nkap0=3,n0=10
  double precision :: rsmh(n0,nkap0),eh(n0,nkap0)
  integer :: iorb,jorb,ki,kj,li,lj,lk,kk,init,io
  ntab= [(iorb,iorb=1,norb)]
  blks =[(sum( [ (2*ltab(io)+1, io=iorb,ntab(iorb)) ] ),iorb=1,norb)]
  where( [(rsmh(ltab(iorb)+1,ktab(iorb)) <= 0, iorb=1,norb)] ) blks = 0 
end subroutine gtbsl8
!> Marks blocks of contiguous l for which rsm and e are unchanged
subroutine gtbsl1(mode,norb,ltab,ktab,rsmh,eh, ntab,blks)
  !  requires l be consecutive and kappa index be constant  in contiguous block
  !  orbital types with rsmh<=0 are excluded
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :=1 rsm must match for each orbital type in contiguous block
  !i             of eh  must match for each orbital type in contiguous block
  !i          =0 not check the matches
  !i   norb  :number of orbital types, i.e. kinds of radial functions 
  !i   ltab  :table of l-quantum numbers for each type (orbl.f)
  !i   ktab  :table of energy index for each type (orbl.f)
  !i   rsmh  :table of smoothing radii in basis (uspecb.f)
  !i   eh    :energy of Hankels in basis 
  !o Outputs
  !o   ntab  :table of upper ranges for each orbital block
  !o   blks  :blks(iorb) = size of contiguous block of orbitals for
  !o         :orbital type iorb. It can be \sum_i (2ltab(iorb)+1) where i=iorb,iorb+1,iorb+2...
  !o         :This enables consecutive orbital types to be grouped into
  !o         :a single block when appropriate, depending on mode.
  !o         :blks(iorb)=0 if the block have been subsumed into prior block 
  !r Remarks
  !r   Routine groups envelope functions, strux into larger blocks.
  !r   for efficient generation of strux and matrix elements.
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nkap0,n0,mode,norb
  integer :: ltab(norb),ktab(norb), ntab(norb),blks(norb)
  parameter (nkap0=3,n0=10)
  double precision :: rsmh(n0,nkap0),eh(n0,nkap0)
  integer :: iorb,jorb,ki,kj,li,lj,lk,kk,init,io
  double precision :: x1,x2
  logical :: neq,deq,lreqe,lreqr,lreql,lreqz,lnever
  real(8),parameter:: eps=1d-8
  lreqr = mode==1
  ntab=-1 
  do  iorb = 1, norb
     if (ntab(iorb)==0) cycle !if iorb block is already sumsumed.
     ki = ktab(iorb) ! Radial funciton index
     li = ltab(iorb) ! l index
     ntab(iorb) = iorb
     if (rsmh(li+1,ki)<= 0) cycle 
     lk = li
     kk = ki
     do  jorb = iorb+1, norb
        kj = ktab(jorb)
        lj = ltab(jorb)
        if (lreqr.and.(abs(eh(lj+1,kj)-eh(li+1,ki))>eps.or.abs(rsmh(lj+1,kj)-rsmh(li+1,ki))>eps)) exit
        if (.not.(lj==lk+1.and.kj==kk)) exit
        ntab(iorb) = jorb !   Increment upper limit
        ntab(jorb) = 0
        lk = lj
        kk = kj
     enddo
  enddo
  blks =[(sum( [ (2*ltab(io)+1, io=iorb,ntab(iorb)) ] ),iorb=1,norb)]
  where( [(rsmh(ltab(iorb)+1,ktab(iorb)) <= 0, iorb=1,norb)] ) blks = 0 
end subroutine gtbsl1
