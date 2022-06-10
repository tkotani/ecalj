subroutine gtbsl1(mode,norb,ltab,ktab,rsmh,eh,ntab,blks)
  !- Marks blocks of contiguous l for which rsm and e are unchanged
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1  rsm must match for each orbital type in contiguous block
  !i         :2  eh  must match for each orbital type in contiguous block
  !i         :4  requires l be consecutive and kappa index be constant
  !i         :   in contiguous block
  !i         :8  never make contiguous blocks
  !i         :16 orbital types with rsmh<=0 are excluded
  !i         :   any combination is allowed;
  !i             however, if 8 is set, 1,2,4 have no effect.
  !i   norb  :number of orbital types, i.e. kinds of radial functions 
  !i   ltab  :table of l-quantum numbers for each type (orbl.f)
  !i   ktab  :table of energy index for each type (orbl.f)
  !i   rsmh  :table of smoothing radii in basis (uspecb.f)
  !i         :unused if 1's digit of mode is not set
  !i   eh    :energy of Hankels in basis (uspecb.f)
  !i         :unused if 2's bit of mode is not set
  !i         :       or if 8's bit of mode is set
  !o Outputs
  !o   ntab  :table of upper ranges for each orbital block
  !o   blks  :blks(iorb) = size of contiguous block of orbitals for
  !o         :orbital type iorb. It can be \sum_i (2ltab(iorb)+1) where i=iorb,iorb+1,iorb+2...
  !o         :This enables consecutive orbital types to be grouped into
  !o         :a single block when appropriate, depending on mode.
  !o         :blks flags orbital types jorb that have been subsumed into
  !o         :prior block by assigning blks(jorb)=0.
  !r Remarks
  !r   Routine groups envelope functions, strux into larger blocks.
  !r   for efficient generation of strux and matrix elements.
  !u Updates
  !u   14 Aug 04 Added mode=8,16
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nkap0,n0,mode,norb
  integer :: ltab(norb),ktab(norb),ntab(norb),blks(norb)
  parameter (nkap0=3,n0=10)
  double precision :: rsmh(n0,nkap0),eh(n0,nkap0)
  integer :: iorb,jorb,ki,kj,li,lj,lk,kk
  double precision :: x1,x2
  logical :: neq,deq,lreqe,lreqr,lreql,lreqz,lnever,bittst,match
  logical:: l_dummy_isanrg,isanrg
  neq(x1,x2) = dabs(x1-x2) > 1d-8
  lreqr = mod(mode,2)==1   ! +1
  lreqe = mod(mode/2,2)==1 ! +2
  lreql = mod(mode/4,2)==1 ! +4
  lnever= mod(mode/8,2)==1 ! +8
  lreqz = mod(mode/16,2)==1! +16
  blks=-1 
  do  iorb = 1, norb
     if (blks(iorb)==0) cycle !if iorb block is already sumsumed.
     ki = ktab(iorb) ! Radial funciton index
     li = ltab(iorb) ! l index
     if (lreqz.and. rsmh(li+1,ki)<= 0) then
        blks(iorb) = 0
        ntab(iorb) = iorb
        cycle 
     endif
     lk = li
     kk = ki
     blks(iorb) = 2*li+1 !  original block size for iorb = 2*l+1 
     ntab(iorb) = iorb
     do  jorb = iorb+1, norb
        kj = ktab(jorb)
        lj = ltab(jorb)
        if (lnever) exit 
        if (lreqe.and.neq(eh(lj+1,kj),eh(li+1,ki))) exit
        if (lreqr.and.neq(rsmh(lj+1,kj),rsmh(li+1,ki))) exit
        if (lreql.and..not.(lj==lk+1  .and. kj==kk)) exit
        ntab(iorb) = jorb !           Increment upper limit
        blks(iorb) = blks(iorb) + 2*lj+1 !  Increment block size(iorb) by size jorb
        blks(jorb) = 0 ! Flag that jorb is subsumed
        ntab(jorb) = 0
        lk = lj
        kk = kj
     enddo
  enddo
end subroutine gtbsl1
