      integer function nnrl(mode,ib1,ib2,iprmb,ndim)
      use m_lmfinit,only: nkaph,nl


C- Returns number of RL channels between sites ib1, ib2
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0 number of RL channels
Ci          1 number of Rl channels (contracted over m)
Ci          2 number of sites containing a basis (contracted over l and m)
Ci         10s digit
Ci          0 count number of channels in ib1..ib2
Ci          1 find maximum number of channels for sites ib1..ib2
Ci   ib1   :starting site index
Ci   ib2   :ending site index
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ndim  :cutoff for iprmb: orbitals for which iprmb<ndim are included
Co Outputs
Co   nnrl  :Number of RL or Rl channels, depending on mode.
Cr Remarks
Cu Updates
Cu   10 Jan 09 Handles 1s digit mode=2
Cu   29 Jun 00 Handles multiple-kappa case; also added 10s digit mode
C ----------------------------------------------------------------------
C     implicit none
      integer mode,ib1,ib2,iprmb(1),ndim
      integer nglob,lmri,ib,nnrli,li,ik,mode0,nnrlx
      logical lmark
      mode0 = mod(mode,10)

      nnrli = 0
      nnrlx = 0
      do  3  ib = ib1, ib2
        lmark = .false.
        lmri = nl*nl*nkaph*(ib-1)
        do  2  ik = 1, nkaph
        do  2  li = 0, nl-1
          lmri = lmri + 2*li+1
          if (iprmb(lmri) .gt. ndim) goto 2
          if (mode0 .eq. 0) then
            nnrli = nnrli + 2*li+1
          elseif (mode0 .eq. 1) then
            nnrli = nnrli + 1
          else
            if (.not. lmark) then
              nnrli = nnrli + 1
              lmark = .true.
            endif
          endif
    2   continue
        if (mode .ge. 10) then
          nnrlx = max(nnrlx,nnrli)
          nnrli = 0
        endif
    3 continue

      if (mode .ge. 10) then
        nnrl = nnrlx
      else
        nnrl = nnrli
      endif

      end

