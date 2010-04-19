      subroutine tbdiag(ltb,nbas,nl,nsp,nspc,nsp1,ipc,nsites,npr,initc,
     .                  qnu,nelts,del,delta,hrs,h0,ov)
C- Adds diagonal part of hamiltonian and overlap
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas,nl,nspc (coupled spins), nsp (uncoupled spins), ipc
Ci   nsites, total number of neighbors in all clusters
Ci   npr(1,i): offset in accumulated list of all neighbors to the
Ci             ith cluster associated with the ith atom
Ci   qnu: holds \eps_0, form of emom for compatibility with LM programs
Ci   nelts: first dimension of del
Ci   del: holds electrostatic increments to the diagonal matrix
Ci          elements read in from the SITE category
Ci   delta: increments for L >= 0, if UL is set
Co Outputs
Co   hrs: diagonal part of hamiltonian set up
Co   ov: 1 added to diagonal part of overlap matrix if lov=T
Cr Remarks
Cr   Spin splitting is achieved through the start parameters, qnu.
Cr   For nspc=2 there are spin-spin blocks
Cr   For TBU the spins are uncoupled
Cr   Overlap is not spin-dependent in TB+U
Cr   h0 is H^in. i.e., H without the e'static increments (TB-L and TB+U)
Cr   if employing the old MRS theory, U1=T, then h0 should point to hrs
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ltb,nbas,nl,nsp,nspc,nsp1,nsites,nelts
      integer ipc(1),npr(0:1,nbas),initc(1)
      double precision qnu(3,nl,nsp1,2),del(nelts,1),
     .  hrs(nl**2,nl**2,1),h0(nl**2,nl**2,1),ov(nl**2,nl**2,1),
     .  delta(nl**2,nl**2,nbas,nsp1)
C Local parameters
      integer i,j,l,lp,ll,isp,ind,offset
      logical lov,bittst,UL,TBU,nonSC

      lov = bittst(ltb,1)
      UL = bittst(ltb,2**15)
      TBU = bittst(ltb,2**13)
      nonSC = .not. UL .and. .not. TBU

C --- copy spin down into spin up ---
      if (nsp .eq. 2) then
        call dcopy(nl**4*nsites,hrs,1,hrs(1,1,nsites+1),1)
      endif

C --- keep a copy of H_0 for making the first order band energy ---
      if (.not. nonSC) then
        call dcopy(nl**4*nsites*nspc**2*nsp,hrs,1,h0,1)
      endif

C --- offset between spin up and down block in real space hamiltonian
      if (nspc .eq. 2) then
        offset = 3
      else
        offset = 1
      endif

      do  isp = 1, nsp1
        do  i = 1, nbas
          if (mod(initc(ipc(i)),2) .ne. 1)
     .      call rxi('tbdiag: missing diagonal H for class no.',ipc(i))
          j = npr(1,i) + (isp-1)*offset*nsites + 1
          do  l = 1, nl**2
            ind = nl*(isp-1) + ll(l) + 1
            if (nonSC) then
              hrs(l,l,j) = hrs(l,l,j) + qnu(2,ll(l)+1,isp,ipc(i))
     .          + del(ind,i)
            else
              do  lp = 1, nl**2
                if (l .eq. lp) then
                  h0(l,lp,j) = h0(l,lp,j) + qnu(2,ll(l)+1,isp,ipc(i))
                  hrs(l,lp,j) = hrs(l,lp,j) + delta(l,lp,i,isp)
     .                        + qnu(2,ll(l)+1,isp,ipc(i))
     .                        + del(ind,i)
                else
                  hrs(l,lp,j) = hrs(l,lp,j) + delta(l,lp,i,isp)
                endif
              enddo
            endif
            if (lov) then
              if ((isp .eq. 1) .or. (isp .eq. 2 .and. nspc .eq. 2)) then
                ov(l,l,j) = ov(l,l,j) + 1d0
              endif
            endif
          enddo
        enddo
      enddo

      end
