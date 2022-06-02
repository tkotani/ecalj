integer function nnrl(mode,ib1,ib2,iprmb,ndim)
  use m_lmfinit,only: nkaph,nl


  !- Returns number of RL channels between sites ib1, ib2
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i          0 number of RL channels
  !i          1 number of Rl channels (contracted over m)
  !i          2 number of sites containing a basis (contracted over l and m)
  !i         10s digit
  !i          0 count number of channels in ib1..ib2
  !i          1 find maximum number of channels for sites ib1..ib2
  !i   ib1   :starting site index
  !i   ib2   :ending site index
  !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
  !i   ndim  :cutoff for iprmb: orbitals for which iprmb<ndim are included
  !o Outputs
  !o   nnrl  :Number of RL or Rl channels, depending on mode.
  !r Remarks
  !u Updates
  !u   10 Jan 09 Handles 1s digit mode=2
  !u   29 Jun 00 Handles multiple-kappa case; also added 10s digit mode
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: mode,ib1,ib2,iprmb(1),ndim
  integer :: nglob,lmri,ib,nnrli,li,ik,mode0,nnrlx
  logical :: lmark
  mode0 = mod(mode,10)

  nnrli = 0
  nnrlx = 0
  do  3  ib = ib1, ib2
     lmark = .false.
     lmri = nl*nl*nkaph*(ib-1)
     do  21  ik = 1, nkaph
        do  2  li = 0, nl-1
           lmri = lmri + 2*li+1
           if (iprmb(lmri) > ndim) goto 2
           if (mode0 == 0) then
              nnrli = nnrli + 2*li+1
           elseif (mode0 == 1) then
              nnrli = nnrli + 1
           else
              if ( .NOT. lmark) then
                 nnrli = nnrli + 1
                 lmark = .true.
              endif
           endif
2       enddo
21   enddo
     if (mode >= 10) then
        nnrlx = max(nnrlx,nnrli)
        nnrli = 0
     endif
3 enddo

  if (mode >= 10) then
     nnrl = nnrlx
  else
     nnrl = nnrli
  endif

end function nnrl

