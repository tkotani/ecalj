      subroutine gtbsl1(mode,norb,ltab,ktab,rsmh,eh,ntab,blks)
C- Marks blocks of contiguous l for which rsm and e are unchanged
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1  rsm must match for each orbital type in contiguous block
Ci         :2  eh  must match for each orbital type in contiguous block
Ci         :4  requires l be consecutive and kappa index be constant
Ci         :   in contiguous block
Ci         :8  never make contiguous blocks
Ci         :16 orbital types with rsmh<=0 are excluded
Ci         :   any combination is allowed;
Ci             however, if 8 is set, 1,2,4 have no effect.
Ci   norb  :number of orbital types, i.e. kinds of radial functions
Ci         :(orbl.f)
Ci   ltab  :table of l-quantum numbers for each type (orbl.f)
Ci   ktab  :table of energy index for each type (orbl.f)
Ci   rsmh  :table of smoothing radii in basis (uspecb.f)
Ci         :unused if 1's digit of mode is not set
Ci   eh    :energy of Hankels in basis (uspecb.f)
Ci         :unused if 2's bit of mode is not set
Ci         :       or if 8's bit of mode is set
Co Outputs
Co   ntab  :table of upper ranges for each orbital block
Co   blks  :blks(iorb) = size of contiguous block of orbitals for
Co         :orbital type iorb and possibly iorb+1,iorb+2...
Co         :This enables consecutive orbital types to be grouped into
Co         :a single block when appropriate, depending on mode.
Co         :blks flags orbital types jorb that have been subsumed into
Co         :prior block by assigning blks(jorb)=0.
Cr Remarks
Cr   Routine groups envelope functions, strux into larger blocks.
Cr   for efficient generation of strux and matrix elements.
Cu Updates
Cu   14 Aug 04 Added mode=8,16
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkap0,n0,mode,norb
      integer ltab(norb),ktab(norb),ntab(norb),blks(norb)
      parameter (nkap0=3,n0=10)
      double precision rsmh(n0,nkap0),eh(n0,nkap0)
C ... Local parameters
      integer iorb,jorb,ki,kj,li,lj,lk,kk
      double precision x1,x2
      logical deq,lreqe,lreqr,lreql,lreqz,lnever,bittst,match
      logical:: l_dummy_isanrg,isanrg
      deq(x1,x2) = dabs(x1-x2) .lt. 1d-8

Ckino isanrg is logical function,       call isanrg(mode,1,31,'gtbsl1:','mode',.true.)
      l_dummy_isanrg=isanrg(mode,1,31,'gtbsl1:','mode',.true.)
      lreqr = bittst(mode,1)
      lreqe = bittst(mode,2)
      lreql = bittst(mode,4)
      lnever= bittst(mode,8)
      lreqz = bittst(mode,16)

C     Set all block sizes to -1 to flag they haven't been touched
      call ivset(blks,1,norb,-1)

      do  iorb = 1, norb
C       If iorb not subsumed into prior orbitals, then mark
        if (blks(iorb) .ne. 0) then
          ki = ktab(iorb)
          li = ltab(iorb)
          if (lreqz) then
            if (rsmh(li+1,ki) .le. 0) then
              blks(iorb) = 0
              ntab(iorb) = iorb
              goto 10
            endif
          endif

          lk = li
          kk = ki
C         Initial block size = size of this l
          blks(iorb) = 2*li+1
          ntab(iorb) = iorb
          do  jorb = iorb+1, norb
            kj = ktab(jorb)
            lj = ltab(jorb)
            match = .not. lnever
            if (.not. match) goto 10
            if (lreqe) then
              match = match .and. deq(eh(lj+1,kj),eh(li+1,ki))
            endif
            if (lreqr) then
              match = match .and. deq(rsmh(lj+1,kj),rsmh(li+1,ki))
            endif
            if (lreql) then
              match = match .and. lj .eq. lk+1
              match = match .and. kj .eq. kk
            endif
            if (match) then
C             Increment upper limit
              ntab(iorb) = jorb
C             Increment block size(iorb) by size jorb
              blks(iorb) = blks(iorb) + 2*lj+1
C             Flag that jorb is subsumed
              blks(jorb) = 0
              ntab(jorb) = 0
              lk = lj
              kk = kj
            else
              goto 10
            endif
          enddo
   10     continue
        endif
      enddo

      end

