c$$$      subroutine makidx(nl,nkaph,ib1,ib2,iprma,sspec,ips,offH,iprmb,
c$$$     .ldham)
c$$$      use m_struc_def  !Cgetarg
c$$$      use m_lmfinit,only: mxorb
c$$$
c$$$
c$$$!! takao think this is too complicated...
c$$$!! I am not sure whether this is OK for PMT mode or not.
c$$$!!
c$$$C- Makes hamiltonian offsets and permutation indices
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   nl    :(global maximum l) + 1
c$$$Ci   nkaph :number of types of one l-quantum number in the basis
c$$$Ci   ib1   :Make offsets for range ib1..ib2
c$$$Ci   ib2   :Make offsets for range ib1..ib2
c$$$Ci   iprma :if supplied, iprma holds a permutation table of site
c$$$Ci         :indices.  Orbital of the hamiltonian are ordered by
c$$$Ci         :sites iprma(1),iprma(2),iprma(3),...
c$$$Ci         :iprma(1)=0=> iprma is not supplied.  Then the orbitals
c$$$Ci         :are ordered by sites 1,2,3,...
c$$$Ci   sspec :struct for species-specific information; see routine uspec
c$$$Ci     Elts read: idxdn lmxb lmxa
c$$$Ci     Stored:    idxdn
c$$$Ci   ips   :species table: site ib belongs to species ips(ib)
c$$$Ci         :Used here to point to idxdn appropriate for site ib
c$$$Co Outputs
c$$$Co   offH  :Tables of hamiltonian offsets in the lower  intermediate
c$$$Co         :and higher blocks.
c$$$Co         : offH(m,ik,ib) contains the total number of orbitals in
c$$$Co         : downfolding block envelope function ik, site ib.
c$$$Co         :See Remarks for further description
c$$$Co         :NB: offH is not set unless offH(1) is >0 on input.
c$$$Co         :NB: if permutation array iprma is supplied, it is the
c$$$Co         :caller's responsibility to initialized offH to zero
c$$$Co
c$$$Co   iprmb :contains a permutation table of orbitals; see Remarks.
c$$$Co
c$$$Co   ldham(i) holds size of block i, plus all the preceding blocks.
c$$$Co
c$$$Cl Local variables
c$$$Cr Remarks
c$$$Cr   This is the main routine that sets up how orbitals are ordered
c$$$Cr   within the hamiltonian.  the arrays it returns also supply
c$$$Cr   information about hamiltonian offsets.
c$$$Cr
c$$$Cr   Each basis function is labelled by an l quantum number and a species
c$$$Cr   index (a function of this type is centered at every site corresponding
c$$$Cr   to the species index).  Also, there may be more than one kind of basis
c$$$Cr   function per site and l, but there can be at most nkaph of such kinds.
c$$$Cr   Thus the maximum possible number of function types associated with a
c$$$Cr   particular site is nkaph*lmxa.
c$$$Cr
c$$$Cr  *idxdn(1..lmxb+1,1:nkaph) keeps track of how functions are used
c$$$Cr   in the construction of the basis.  For a particular 0<=l<=lmxb and
c$$$Cr   1<=ik<nkaph, idxdn takes one of the following:
c$$$Cr    value   Signifies
c$$$Cr      0     Role of this orbital has not yet been determined
c$$$Cr      1     Orbital is included as "active" orbitals, which means
c$$$Cr            they are included in the hamiltonian and diagonalized
c$$$Cr      2     Orbital is treated as an "intermediate" orbital, which
c$$$Cr            means it is downfolded and included in a perturbative way
c$$$Cr      3     Orbital is treated as an "high" orbital, which means
c$$$Cr            means it is included in augmentation for tails of other
c$$$Cr            orbitals, but is otherwise not part of the basis.
c$$$Cr      4     Orbital is completely neglected (no orbital is defined)
c$$$Cr     10     Orbital is a local orbital whose value and slope are
c$$$Cr            constructed to be zero at the MT boundary.
c$$$Cr            It is included in the basis.
c$$$Cr     11     Orbital is a local orbital with a smooth Hankel tail
c$$$Cr            and it is included in the basis.
c$$$Cr     12     Orbital is a local orbital with a smooth Hankel tail
c$$$Cr            It is incorporated perturbatively to the hamiltonian
c$$$Cr            and is not assembled as part of the hamiltonian matrix
c$$$Cr
c$$$Cr  *offH and iprmb both contain information about the nature
c$$$Cr   and ordering of orbitals making up the hamiltonian matrix.
c$$$Cr
c$$$Cr  *iprmb is a permutation table grouping orbitals according into
c$$$Cr   lower, intermediate, higher and neglected blocks.  In the
c$$$Cr   following code segment, iprmb(lmr)  points which the row in
c$$$Cr   the hamiltonian orbital lm corresponding to site ib is found:
c$$$Cr
c$$$Cr     mxorb = nglob('mxorb')
c$$$Cr     lmr = mxorb*(ib1-1)
c$$$Cr     do  ib = ib1, ib2
c$$$Cr     do  lm = 1, mxorb
c$$$Cr       lmr = lmr+1
c$$$Cr       orbital-row = iprmb(lmr)
c$$$Cr       if (orbital-row is in the bounds of block you want) then
c$$$Cr         ...
c$$$Cr       endif
c$$$Cr
c$$$Cr   Alternatively, you can think of iprmb as a two-dimensional array
c$$$Cr   iprmb(mxorb,nbas).  The leading dimension `mxorb' is a global
c$$$Cr   variable, set to nkaph * nl**2.
c$$$Cr
c$$$Cr   Whether an orbital is present in or missing from a block can be
c$$$Cr   inferred from whether `orbital-row' is within the bounds of
c$$$Cr   the desired block (defined by offH or ldham).
c$$$Cr
c$$$Cr   Some routines needs a list of which l-quantum numbers are
c$$$Cr   specified in a block; use routine `orbl' to extract it, eg:
c$$$Cr     ldim = ldham(1)
c$$$Cr     do  ib = 1, nbas
c$$$Cr       call upack('site spec pos',ssite,ib,is,p,0,0)
c$$$Cr       call orbl(ib,0,ldim,iprmb,norb,ltab,ktab,xx,offl,xx)
c$$$Cr       do  iorb = 1, norb
c$$$Cr         l   = ltab(iorb)  <- l quantum number
c$$$Cr         ik  = ktab(iorb)  <- not needed for single-kappa hamiltonians
c$$$Cr         off = offl(iorb)  <- hamiltonian offset
c$$$Cr         ...
c$$$Cr
c$$$Cr  *norb contains the number of (l,kappa) combinations belonging to a
c$$$Cr   site (ie number of orbitals supressing m degeneracy).
c$$$Cr
c$$$Cr  *offH(m,ik,ib) contains the total number of orbitals in block
c$$$Cr   m for all sites i=1..ib-1.  m is a block index:
c$$$Cr     m = 1  => lower block
c$$$Crxxx     m = 2  => intermediate block
c$$$Cr     m = 3  => higher block
c$$$Crxxx     m = 4  => lower+intermediate block
c$$$Crxxx     m = 5  => lower+intermediate+higher block
c$$$Cr
c$$$Cr   offH(m,ik,ib) is set for ib=ib1..ib2+1.  Thus the size of block m
c$$$Cr   for site ib = offH(m,1,ib+1) - offH(m,1,ib).  The size of a block
c$$$Cr   of a particular kappa is offH(m,ik+1,ib),offH(m,ik,ib)
c$$$Cr   NB: if a site permutation array is supplied, block sizes must
c$$$Cr   be computed from offH(m,1,iprma(ib+1)) - offH(m,1,iprma(ib))
c$$$Cr
c$$$Cb Bugs
c$$$Cb   if a permutation table iprma is supplied, and ib1 and ib2 do
c$$$Cb   not span the whole range of sites, offH is not initialized
c$$$Cu Updates
c$$$Cu   24 Jun 04 Adapted to new definition for idxdn (see Remarks)
c$$$Cu   10 Apr 02 Redimensioned idxdn to accomodate larger lmax
c$$$Cu   28 Jun 00 Automatic switching to high orbitals with l>lmxb
c$$$Cu    8 Jun 00 Extended to multiple-kappa case.  New argument list.
c$$$Cu   19 Oct 99 Routine completely rewritten
c$$$Cu    4 May 98 Bug fix.  Now requires priorities to be numbered 1..3
c$$$C ----------------------------------------------------------------
c$$$      implicit none
c$$$C Passed parameters
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$      integer nl,nkaph,ib1,ib2,nkap0,n0,n0H,ips(ib2),iprma(ib2),i_copy_size
c$$$      parameter (nkap0=3,n0=10,n0H=5)
c$$$      integer offh(n0H,nkap0,ib2+1),iprmb(*),ldham(4) !offh(n0H,nkap0,4) bugfix sep2010
c$$$      ! offh(n0H,nkap0,5) : bug fix sep2010
c$$$C Local parameters
c$$$      integer indx,iposn,l,m,iprior,j,lgunit,idxdn(n0,nkap0),is,ib,i,
c$$$     .ltab(n0),nd,ndim,ipr,maxp,ik,ndimi(n0H,nkap0),isum,jb,
c$$$     .nglob,stdo,stdl,lmxb,lmxa,ibp1,idxpr
c$$$      logical lx,lprmib,lz
c$$$      character*80 outs,strn*5
c$$$      parameter (maxp=3)
c$$$
c$$$      logical:: l_dummy_isanrg, isanrg
c$$$
c$$$C --- Setup ---
c$$$c      if (nkaph .gt. nkap0) call rx('makidx: increase nkap0')
c$$$c      lprmib = iprma(1) .gt. 0
c$$$
c$$$      
c$$$c$$$      if (offH(1,1,1) .eq. -1) goto 100
c$$$c$$$      if (lprmib) then
c$$$c$$$      else
c$$$c$$$        call iinit(offH(1,1,ib1),n0H*nkap0*(ib2-ib1+2))
c$$$c$$$      endif
c$$$c$$$      nd = dlog(dble(ib2)) / dlog(10d0)
c$$$c$$$      if (ipr .ge. 50) write(stdo,333)
c$$$c$$$  333 format(/' Makidx:  basis arranged in downfolding order:'
c$$$c$$$     ./'  ib',5x,'low',6x,'intermed',7x,'high',8x,'.. offH ..')
c$$$c$$$
c$$$c$$$C --- For each site ib = ib1 ... ib2, make offH(ib) ---
c$$$c$$$      do  jb = ib1, ib2
c$$$c$$$        if (lprmib) then
c$$$c$$$          ib = iprma(jb)
c$$$c$$$          ibp1 = iprma(jb+1)
c$$$c$$$        else
c$$$c$$$          ib = jb
c$$$c$$$          ibp1 = jb+1
c$$$c$$$        endif
c$$$c$$$        is = ips(ib)
c$$$c$$$
c$$$c$$$        i_copy_size=size(sspec(is)%idxdn)
c$$$c$$$        call icopy(i_copy_size,sspec(is)%idxdn,1,idxdn,1)
c$$$c$$$        lmxb=sspec(is)%lmxb
c$$$c$$$        lmxa=sspec(is)%lmxa
c$$$c$$$
c$$$c$$$C       orbs > lmxb always high; orbs > lmxa always neglected
c$$$c$$$        do  ik = 1, nkaph
c$$$c$$$          do  l = lmxb+1, nl-1
c$$$c$$$            if (idxpr(idxdn(l+1,ik)) .lt. 3) idxdn(l+1,ik) = 3
c$$$c$$$            if (idxpr(idxdn(l+1,ik)) .lt. 4 .and. l .gt. lmxa)
c$$$c$$$     .      idxdn(l+1,ik) = 4
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$c$$$        i_copy_size=size(sspec(is)%idxdn)
c$$$c$$$        call icopy(i_copy_size,idxdn,1,sspec(is)%idxdn,1)
c$$$c$$$ccccccccccccccccccccc
c$$$c$$$c        cycle
c$$$c$$$cccccccccccccccccccc        
c$$$c$$$c        sspec(is)%idxdn = idxdn
c$$$c$$$        ndimi=0 !call iinit(ndimi,n0H*nkap0)
c$$$c$$$        if (nl .gt. n0) call rx('makidx need more space for idxdn')
c$$$c$$$
c$$$c$$$C   ... Assemble offH(1..maxp,1..nkaph) for this site
c$$$c$$$        do  ik = 1, nkaph
c$$$c$$$          call awrit2('%x%,ni+',outs,80,0,nd+4,ib)
c$$$c$$$          if (ik .ne. 1) call awrit2('%x%npk%i+',outs,80,0,nd+2,ik)
c$$$c$$$          call mkidx2(nl,maxp,idxdn(1,ik),ltab,ndimi(1,ik),outs)
c$$$c$$$          do  iprior = 1, maxp
c$$$c$$$            ndim = offH(iprior,ik,ib) + ndimi(iprior,ik)
c$$$c$$$            if (ik .lt. nkaph) offH(iprior,ik+1,ib) = ndim
c$$$c$$$            if (ik .eq. nkaph) offH(iprior,1,ibp1) = ndim
c$$$c$$$          enddo
c$$$c$$$          if (isum(maxp,ndimi(1,ik),1) .ne. 0 .and. ipr .ge. 50)
c$$$c$$$     .    call awrit4('%(n*13-11)p%j%,5i%,5i%,5i',
c$$$c$$$     .    outs,80,-stdo,maxp+1,offH(1,ik,ib),offH(2,ik,ib),
c$$$c$$$     .    offH(3,ik,ib))
c$$$c$$$C
c$$$c$$$C     .      call awrit4('%(n*13-11)p %,5i%,5i%,5i',
c$$$c$$$C     .      outs,80,-stdo,maxp+1,offH(1,ik,ib),offH(2,ik,ib),
c$$$c$$$C     .      offH(3,ik,ib))
c$$$c$$$        enddo
c$$$c$$$
c$$$c$$$C   ... Fill out offH for kappa's between nkaph and nkap0
c$$$c$$$        do  ik = nkaph+1, nkap0
c$$$c$$$          do  iprior = 1, maxp
c$$$c$$$            offH(iprior,ik,ib) = offH(iprior,1,ibp1)
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$c$$$
c$$$c$$$C   ... Concatenate l+i and l+i+h blocks
c$$$c$$$        do  ik = 1, nkaph
c$$$c$$$          offH(maxp+1,ik,ibp1) = offH(maxp+1,1,ib)
c$$$c$$$     .    + isum(nkaph,ndimi(1,1),n0H)
c$$$c$$$     .    + isum(nkaph,ndimi(2,1),n0H)
c$$$c$$$     .    + isum(ik-1,ndimi(1,1),1)
c$$$c$$$     .    + isum(ik-1,ndimi(2,1),1)
c$$$c$$$          offH(maxp+2,ik,ibp1) = offH(maxp+2,1,ib)
c$$$c$$$     .    + isum(nkaph,ndimi(1,1),n0H)
c$$$c$$$     .    + isum(nkaph,ndimi(2,1),n0H)
c$$$c$$$     .    + isum(nkaph,ndimi(3,1),n0H)
c$$$c$$$     .    + isum(ik-1,ndimi(1,1),1)
c$$$c$$$     .    + isum(ik-1,ndimi(2,1),1)
c$$$c$$$     .    + isum(ik-1,ndimi(3,1),1)
c$$$c$$$        enddo
c$$$c$$$        do  ik = nkaph+1, nkap0
c$$$c$$$          do  iprior = maxp+1, maxp+2
c$$$c$$$            offH(maxp+1,ik,ib) = offH(maxp+1,nkaph,ib)
c$$$c$$$     .      + ndimi(1,nkaph) + ndimi(2,nkaph)
c$$$c$$$            offH(maxp+2,ik,ib) = offH(maxp+2,nkaph,ib)
c$$$c$$$     .      + ndimi(1,nkaph) + ndimi(2,nkaph) + ndimi(3,nkaph)
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$c$$$      enddo
c$$$c$$$      
c$$$c$$$C     This should never be needed
c$$$c$$$C     do  ik = 1, nkaph-1
c$$$c$$$C       do  iprior = 1, maxp
c$$$c$$$C         offH(iprior,ik+1,ib2+1) = offH(iprior,1,ib2+1)
c$$$c$$$C       enddo
c$$$c$$$C     enddo
c$$$
c$$$C --- Make iprmb ---
c$$$c  100 continue
c$$$C     For now, require that morb = nl*nkaph
c$$$Changenglob      mxorb = nglob('mxorb')
c$$$c      mxorb = globalvariables%mxorb
c$$$c      i = nl**2*nkaph
c$$$cCkino isanrg is logical function,       call isanrg(mxorb,i,i,'makidx:','mxorb', .true.)
c$$$c      l_dummy_isanrg=isanrg(mxorb,i,i,'makidx:','mxorb', .true.)
c$$$c     call iinit(ndimi,n0H*nkap0)
c$$$c      ndimi=0
c$$$      call getpr(ipr)
c$$$      stdo = lgunit(1)
c$$$      stdl = lgunit(2)
c$$$      indx = 0
c$$$      do  112  iprior = 1, maxp+1
c$$$        do  110  jb = ib1, ib2
c$$$          ib = jb
c$$$          ibp1 = jb+1
c$$$          is = ips(ib)
c$$$          iposn = mxorb*(ib-1)
c$$$          idxdn=sspec(is)%idxdn
c$$$          do  1121 ik = 1, nkaph
c$$$          do  111  l = 0, nl-1
c$$$             lz =.false.
c$$$             if( (idxdn(l+1,ik)==10.or.idxdn(l+1,ik)==11) .and. iprior==1) lz=.true.
c$$$             lx = lz .or. (idxdn(l+1,ik)==iprior)
c$$$            do  120  m = -l, l
c$$$              iposn = iposn+1
c$$$              if (lx) then
c$$$                indx = indx+1
c$$$                iprmb(iposn) = indx
c$$$c                ndimi(iprior,ik) = ndimi(iprior,ik)+1
c$$$              endif
c$$$  120       continue
c$$$  111     continue
c$$$ 1121     continue
c$$$  110   continue
c$$$        ldham(iprior) = indx
c$$$  112 continue
c$$$
c$$$c$$$C --- Printout ---
c$$$c$$$      if (ipr .lt. 10) return
c$$$c$$$      if (nkaph .eq. 1) then
c$$$c$$$        call awrit4('dm  L %i  I %i  H %i  N %i ',' ',80,stdl,ldham(1),
c$$$c$$$     .  ldham(2)-ldham(1),ldham(3)-ldham(2),ldham(4)-ldham(3))
c$$$c$$$      else
c$$$c$$$        call icopy(nkaph,ndimi(1,1),n0H,ltab,1)
c$$$c$$$        call awrit3('%xdm  L %i (%n:1i )',outs,80,0,ldham(1),nkaph,
c$$$c$$$     .  ltab)
c$$$c$$$        call icopy(nkaph,ndimi(2,1),n0H,ltab,1)
c$$$c$$$        call awrit3('%a  I %i (%n:1i )',outs,80,0,ldham(2)-ldham(1),
c$$$c$$$     .  nkaph,ltab)
c$$$c$$$        call icopy(nkaph,ndimi(3,1),n0H,ltab,1)
c$$$c$$$        call awrit4('%a  H %i (%n:1i )  N %i',outs,80,-stdl,
c$$$c$$$     .  ldham(3)-ldham(2),nkaph,ltab,ldham(4)-ldham(3))
c$$$c$$$      endif
c$$$c$$$      if (ipr .lt. 20) return
c$$$c$$$      call awrit5('%?;(n>=30);%N;; Makidx:  hamiltonian dimensions'//
c$$$c$$$     .' Low, Int, High, Negl: %i %i %i %i',' ',80,stdo,ipr,ldham(1),
c$$$c$$$     .ldham(2)-ldham(1),ldham(3)-ldham(2),ldham(4)-ldham(3))
c$$$c$$$      if (ipr .ge. 40 .or. (ipr .gt. 30 .and. nkaph .gt. 1)) then
c$$$c$$$        write(stdo,345)
c$$$c$$$  345   format(' kappa',3x,
c$$$c$$$     .  'Low',3x,'Int',3x,'High',2x,'L+I',2x,'L+I+H',2x,'Neglected')
c$$$c$$$        if (nkaph .gt. 1) then
c$$$c$$$          do  ik = 1, nkaph
c$$$c$$$            i = isum(3,ndimi(1,ik),1)
c$$$c$$$            write(stdo,346) ik, (ndimi(j,ik), j=1,3),
c$$$c$$$     .      ndimi(1,ik)+ndimi(2,ik),i,nl**2*(ib2-ib1+1)-i
c$$$c$$$  346       format(i4,2x,5i6,i8)
c$$$c$$$  347       format(1x,a5,5i6,i8)
c$$$c$$$          enddo
c$$$c$$$          strn = ' all '
c$$$c$$$        else
c$$$c$$$          strn = '  -  '
c$$$c$$$        endif
c$$$c$$$        write(stdo,347) strn, ldham(1),ldham(2)-
c$$$c$$$     .  ldham(1),ldham(3)-ldham(2),ldham(2),ldham(3),ldham(4)-
c$$$c$$$     .  ldham(3)
c$$$c$$$      endif
c$$$c$$$ccccccccccccccccccccccccccccccccccc
c$$$c$$$      do  ib = ib1, ib2+1
c$$$c$$$         print '(i5,5i5,2x,5i5,2x,5i5)', ib,(offH(j,1,ib), j=1,5)
c$$$c$$$      enddo
c$$$c      i = 0
c$$$c      do  210  ib = ib1, ib2
c$$$c      do  210  ik = 1,  nkaph
c$$$c      do  210  m = 1,  nl**2
c$$$c        i = i+1
c$$$c        write(6,"('iperbbb ',6i5)") ib, ik, m, i, iprmb(i)
c$$$c  210 continue
c$$$c      call rx('done')
c$$$cccccccccccccccccccccccccccccccccc
c$$$      end subroutine makidx
c$$$
c$$$c$$$      subroutine mkidx2(nl,maxp,idxdn,ltab,ndimi,outs)
c$$$c$$$
c$$$c$$$C ----------------------------------------------------------------------
c$$$c$$$Ci Inputs
c$$$c$$$Ci   nl    :(global maximum l) + 1
c$$$c$$$Ci   maxp  :maximum number
c$$$c$$$Ci   idxdn :controls how orbital is to be included in the basis
c$$$c$$$Co Outputs
c$$$c$$$Co   ltab  :list of l's in permuted order
c$$$c$$$Co   ndimi :sum (ltab(i)+1)**2 for each priority
c$$$c$$$Co   outs
c$$$c$$$Cr Remarks
c$$$c$$$Cu Updates
c$$$c$$$C ----------------------------------------------------------------------
c$$$c$$$C     implicit none
c$$$c$$$C ... Passed parameters
c$$$c$$$      integer n0,n0H
c$$$c$$$      parameter (n0=10,n0H=5)
c$$$c$$$      integer nl,maxp,ltab(n0),idxdn(n0),ndimi(n0H)
c$$$c$$$C ... Local parameters
c$$$c$$$      integer iprior,i2,i1,norb,iwk(n0),ndim,ipr,i,idxdn2(n0),idxpr
c$$$c$$$      character*80 outs
c$$$c$$$      character onam(9)*1,ostr*10
c$$$c$$$      data onam /'s','p','d','f','g','h','i','j','l'/
c$$$c$$$
c$$$c$$$C     call icopy(nl,idxdn,1,iwk2,1)
c$$$c$$$      call getpr(ipr)
c$$$c$$$      ipr = 51
c$$$c$$$      do  10  i1 = 1, nl
c$$$c$$$        idxdn2(i1) = idxpr(idxdn(i1))
c$$$c$$$   10 continue
c$$$c$$$      call ivheap(1,nl,idxdn2,iwk,101)
c$$$c$$$      call iinit(ndimi,maxp)
c$$$c$$$      iprior = 0
c$$$c$$$      i2 = 0
c$$$c$$$   11 i1 = i2+1
c$$$c$$$      iprior = iprior+1
c$$$c$$$C ... We have exhausted this group of orbitals
c$$$c$$$      if (i1 .gt. nl .or. iprior .gt. maxp) return
c$$$c$$$
c$$$c$$$      i2 = i1-1
c$$$c$$$      norb = 0
c$$$c$$$      ndim = 0
c$$$c$$$C ... Keep looping here until find all orbitals of this priority
c$$$c$$$   12 i2 = i2+1
c$$$c$$$      if (i2 .le. nl) then
c$$$c$$$        if (idxdn2(iwk(i2)) .eq. iprior) then
c$$$c$$$          norb = norb+1
c$$$c$$$          ndim = ndim + 2*(iwk(i2)-1) + 1
c$$$c$$$          ltab(norb) = iwk(i2)-1
c$$$c$$$          goto 12
c$$$c$$$        endif
c$$$c$$$      endif
c$$$c$$$      ndimi(iprior) = ndim
c$$$c$$$C     offH(iprior,1,ibp1) = offH(iprior,1,ibp1) + ndim
c$$$c$$$      if (ipr .ge. 50) then
c$$$c$$$        write(ostr,'(10a1)')(' ', i=norb+1,5),(onam(ltab(i)+1),i=1,norb)
c$$$c$$$        call awrit2('%a%b '//ostr(1:5)//' (%i)%(n*13+2)p%j+',
c$$$c$$$     .  outs,80,0,ndim,iprior)
c$$$c$$$      endif
c$$$c$$$      i2 = i2-1
c$$$c$$$      goto 11
c$$$c$$$
c$$$c$$$      end subroutine mkidx2
c$$$
c$$$      integer function idxpr(idxdn)
c$$$
c$$$C- Given a value idxdn, assign priority to orbital (kernel called by makidx)
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   idxdn :downfolding indices
c$$$Co Outputs
c$$$Co   idxpr :priority in the hamiltonian corresponding to idxdn
c$$$Cl Local variables
c$$$Cl         :
c$$$Cr Remarks
c$$$Cr   idxdn indicates how a basis function is to be used in the
c$$$Cr   construction of the basis;
c$$$Cr   idxpr uses idxdn to assign a priority according to what level
c$$$Cr   of treatment it has:
c$$$Cr   idxdn idxpr  Signifies
c$$$Cr     0    0     Role of this orbital has not yet been determined
c$$$Cr     1    1     Orbital is included as "active" orbitals, which means
c$$$Cr                they are included in the hamiltonian and diagonalized
c$$$Cr     2    2     Orbital is treated as an "intermediate" orbital, which
c$$$Cr                means it is downfolded and included in a perturbative way
c$$$Cr     3    3     Orbital is treated as an "high" orbital, which means
c$$$Cr                means it is included in augmentation for tails of other
c$$$Cr                orbitals, but is otherwise not part of the basis.
c$$$Crxxx     4    4     Orbital is completely neglected (no orbital is defined)
c$$$Cr    10    1     Orbital is a local orbital whose value and slope are
c$$$Cr                constructed to be zero at the MT boundary.
c$$$Cr                It is included in the basis.
c$$$Cr    11    1     Orbital is a local orbital with a smooth Hankel tail
c$$$Cr                and it is included in the basis.
c$$$Crxxx    12    2     Orbital is a local orbital with a smooth Hankel tail
c$$$Crxxx                It is incorporated perturbatively to the hamiltonian
c$$$Crxxx                and is not assembled as part of the hamiltonian matrix
c$$$Cr
c$$$Cu Updates
c$$$Cu   24 Jun 04 First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer idxdn
c$$$C ... Local parameters
c$$$
c$$$      if (idxdn .le. 4) then
c$$$        idxpr = idxdn
c$$$      elseif (idxdn .eq. 10 .or. idxdn .eq. 11) then
c$$$        idxpr = 1
c$$$      elseif (idxdn .eq. 12) then
c$$$        idxpr = 2
c$$$      else
c$$$        call rxi('makidx: bad idxdn',idxdn)
c$$$      endif
c$$$
c$$$      end function idxpr
c$$$
c$$$
