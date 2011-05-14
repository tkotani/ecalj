      subroutine makidx(nl,nkaph,ib1,ib2,iprma,sspec,ips,offH,iprmb,
     .ldham)

      use m_struc_def  !Cgetarg
      use m_globalvariables

C- Makes hamiltonian offsets and permutation indices
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   nkaph :number of types of one l-quantum number in the basis
Ci   ib1   :Make offsets for range ib1..ib2
Ci   ib2   :Make offsets for range ib1..ib2
Ci   iprma :if supplied, iprma holds a permutation table of site
Ci         :indices.  Orbital of the hamiltonian are ordered by
Ci         :sites iprma(1),iprma(2),iprma(3),...
Ci         :iprma(1)=0=> iprma is not supplied.  Then the orbitals
Ci         :are ordered by sites 1,2,3,...
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: idxdn lmxb lmxa
Ci     Stored:    idxdn
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci         :Used here to point to idxdn appropriate for site ib
Co Outputs
Co   offH  :Tables of hamiltonian offsets in the lower, intermediate
Co         :and higher blocks.
Co         : offH(m,ik,ib) contains the total number of orbitals in
Co         : downfolding block envelope function ik, site ib.
Co         :See Remarks for further description
Co         :NB: offH is not set unless offH(1) is >0 on input.
Co         :NB: if permutation array iprma is supplied, it is the
Co         :caller's responsibility to initialized offH to zero
Co
Co   iprmb :contains a permutation table of orbitals; see Remarks.
Co
Co   ldham(i) holds size of block i, plus all the preceding blocks.
Co
Cl Local variables
Cr Remarks
Cr   This is the main routine that sets up how orbitals are ordered
Cr   within the hamiltonian.  the arrays it returns also supply
Cr   information about hamiltonian offsets.
Cr
Cr   Each basis function is labelled by an l quantum number and a species
Cr   index (a function of this type is centered at every site corresponding
Cr   to the species index).  Also, there may be more than one kind of basis
Cr   function per site and l, but there can be at most nkaph of such kinds.
Cr   Thus the maximum possible number of function types associated with a
Cr   particular site is nkaph*lmxa.
Cr
Cr  *idxdn(1..lmxb+1,1:nkaph) keeps track of how functions are used
Cr   in the construction of the basis.  For a particular 0<=l<=lmxb and
Cr   1<=ik<nkaph, idxdn takes one of the following:
Cr    value   Signifies
Cr      0     Role of this orbital has not yet been determined
Cr      1     Orbital is included as "active" orbitals, which means
Cr            they are included in the hamiltonian and diagonalized
Cr      2     Orbital is treated as an "intermediate" orbital, which
Cr            means it is downfolded and included in a perturbative way
Cr      3     Orbital is treated as an "high" orbital, which means
Cr            means it is included in augmentation for tails of other
Cr            orbitals, but is otherwise not part of the basis.
Cr      4     Orbital is completely neglected (no orbital is defined)
Cr     10     Orbital is a local orbital whose value and slope are
Cr            constructed to be zero at the MT boundary.
Cr            It is included in the basis.
Cr     11     Orbital is a local orbital with a smooth Hankel tail
Cr            and it is included in the basis.
Cr     12     Orbital is a local orbital with a smooth Hankel tail
Cr            It is incorporated perturbatively to the hamiltonian
Cr            and is not assembled as part of the hamiltonian matrix
Cr
Cr  *offH and iprmb both contain information about the nature
Cr   and ordering of orbitals making up the hamiltonian matrix.
Cr
Cr  *iprmb is a permutation table grouping orbitals according into
Cr   lower, intermediate, higher and neglected blocks.  In the
Cr   following code segment, iprmb(lmr)  points which the row in
Cr   the hamiltonian orbital lm corresponding to site ib is found:
Cr
Cr     mxorb = nglob('mxorb')
Cr     lmr = mxorb*(ib1-1)
Cr     do  ib = ib1, ib2
Cr     do  lm = 1, mxorb
Cr       lmr = lmr+1
Cr       orbital-row = iprmb(lmr)
Cr       if (orbital-row is in the bounds of block you want) then
Cr         ...
Cr       endif
Cr
Cr   Alternatively, you can think of iprmb as a two-dimensional array
Cr   iprmb(mxorb,nbas).  The leading dimension `mxorb' is a global
Cr   variable, set to nkaph * nl**2.
Cr
Cr   Whether an orbital is present in or missing from a block can be
Cr   inferred from whether `orbital-row' is within the bounds of
Cr   the desired block (defined by offH or ldham).
Cr
Cr   Some routines needs a list of which l-quantum numbers are
Cr   specified in a block; use routine `orbl' to extract it, eg:
Cr     ldim = ldham(1)
Cr     do  ib = 1, nbas
Cr       call upack('site spec pos',ssite,ib,is,p,0,0)
Cr       call orbl(ib,0,ldim,iprmb,norb,ltab,ktab,xx,offl,xx)
Cr       do  iorb = 1, norb
Cr         l   = ltab(iorb)  <- l quantum number
Cr         ik  = ktab(iorb)  <- not needed for single-kappa hamiltonians
Cr         off = offl(iorb)  <- hamiltonian offset
Cr         ...
Cr
Cr  *norb contains the number of (l,kappa) combinations belonging to a
Cr   site (ie number of orbitals supressing m degeneracy).
Cr
Cr  *offH(m,ik,ib) contains the total number of orbitals in block
Cr   m for all sites i=1..ib-1.  m is a block index:
Cr     m = 1  => lower block
Cr     m = 2  => intermediate block
Cr     m = 3  => higher block
Cr     m = 4  => lower+intermediate block
Cr     m = 5  => lower+intermediate+higher block
Cr
Cr   offH(m,ik,ib) is set for ib=ib1..ib2+1.  Thus the size of block m
Cr   for site ib = offH(m,1,ib+1) - offH(m,1,ib).  The size of a block
Cr   of a particular kappa is offH(m,ik+1,ib),offH(m,ik,ib)
Cr   NB: if a site permutation array is supplied, block sizes must
Cr   be computed from offH(m,1,iprma(ib+1)) - offH(m,1,iprma(ib))
Cr
Cb Bugs
Cb   if a permutation table iprma is supplied, and ib1 and ib2 do
Cb   not span the whole range of sites, offH is not initialized
Cu Updates
Cu   24 Jun 04 Adapted to new definition for idxdn (see Remarks)
Cu   10 Apr 02 Redimensioned idxdn to accomodate larger lmax
Cu   28 Jun 00 Automatic switching to high orbitals with l>lmxb
Cu    8 Jun 00 Extended to multiple-kappa case.  New argument list.
Cu   19 Oct 99 Routine completely rewritten
Cu    4 May 98 Bug fix.  Now requires priorities to be numbered 1..3
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      type(s_spec)::sspec(*)

      integer nl,nkaph,ib1,ib2,nkap0,n0,n0H,ips(ib2),iprma(ib2)
      parameter (nkap0=3,n0=10,n0H=5)
      integer offh(n0H,nkap0,ib2+1),iprmb(*),ldham(4) !offh(n0H,nkap0,4) bugfix sep2010
            ! offh(n0H,nkap0,5) : bug fix sep2010
C Local parameters
      integer indx,iposn,l,m,iprior,j,lgunit,idxdn(n0,nkap0),is,ib,i,
     .ltab(n0),nd,ndim,ipr,maxp,ik,ndimi(n0H,nkap0),isum,mxorb,jb,
     .nglob,stdo,stdl,lmxb,lmxa,ibp1,idxpr
      logical lx,lprmib
      character*80 outs,strn*5
      parameter (maxp=3)

C --- Setup ---
      if (nkaph .gt. nkap0) call rx('makidx: increase nkap0')
      stdo = lgunit(1)
      stdl = lgunit(2)
      lprmib = iprma(1) .gt. 0

      call getpr(ipr)
      if (offH(1,1,1) .eq. -1) goto 100
      if (lprmib) then
C        do  jb = ib1, ib2
C          ib = iprma(jb)
C          if (ib .le. 0) call rx('makidx: bad permutation array')
C          call iinit(offH(1,1,ib),n0H*nkap0)
C        enddo
C        ib = iprma(ib2+1)
C        call iinit(offH(1,1,ib),n0H*nkap0)
      else
        call iinit(offH(1,1,ib1),n0H*nkap0*(ib2-ib1+2))
      endif
      nd = dlog(dble(ib2)) / dlog(10d0)
      if (ipr .ge. 50) write(stdo,333)
  333 format(/' Makidx:  basis arranged in downfolding order:'
     ./'  ib',5x,'low',6x,'intermed',7x,'high',8x,'.. offH ..')

C --- For each site ib = ib1 ... ib2, make offH(ib) ---
      do  jb = ib1, ib2
        if (lprmib) then
          ib = iprma(jb)
          ibp1 = iprma(jb+1)
        else
          ib = jb
          ibp1 = jb+1
        endif
        is = ips(ib)

        i_copy_size=size(sspec(is)%idxdn)
        call i8icopy(i_copy_size,sspec(is)%idxdn,1,idxdn,1)
        lmxb=sspec(is)%lmxb
        lmxa=sspec(is)%lmxa

C       orbs > lmxb always high; orbs > lmxa always neglected
        do  ik = 1, nkaph
          do  l = lmxb+1, nl-1
            if (idxpr(idxdn(l+1,ik)) .lt. 3) idxdn(l+1,ik) = 3
            if (idxpr(idxdn(l+1,ik)) .lt. 4 .and. l .gt. lmxa)
     .      idxdn(l+1,ik) = 4
          enddo
        enddo

        i_copy_size=size(sspec(is)%idxdn)
        call ii8copy(i_copy_size,idxdn,1,sspec(is)%idxdn,1)

        call iinit(ndimi,n0H*nkap0)
        if (nl .gt. n0) call rx('makidx need more space for idxdn')

C   ... Assemble offH(1..maxp,1..nkaph) for this site
        do  ik = 1, nkaph
          call awrit2('%x%,ni+',outs,80,0,nd+4,ib)
          if (ik .ne. 1) call awrit2('%x%npk%i+',outs,80,0,nd+2,ik)
          call mkidx2(nl,maxp,idxdn(1,ik),ltab,ndimi(1,ik),outs)
          do  iprior = 1, maxp
            ndim = offH(iprior,ik,ib) + ndimi(iprior,ik)
            if (ik .lt. nkaph) offH(iprior,ik+1,ib) = ndim
            if (ik .eq. nkaph) offH(iprior,1,ibp1) = ndim
          enddo
          if (isum(maxp,ndimi(1,ik),1) .ne. 0 .and. ipr .ge. 50)
     .    call awrit4('%(n*13-11)p%j%,5i%,5i%,5i',
     .    outs,80,-stdo,maxp+1,offH(1,ik,ib),offH(2,ik,ib),
     .    offH(3,ik,ib))
C
C     .      call awrit4('%(n*13-11)p %,5i%,5i%,5i',
C     .      outs,80,-stdo,maxp+1,offH(1,ik,ib),offH(2,ik,ib),
C     .      offH(3,ik,ib))
        enddo

C   ... Fill out offH for kappa's between nkaph and nkap0
        do  ik = nkaph+1, nkap0
          do  iprior = 1, maxp
            offH(iprior,ik,ib) = offH(iprior,1,ibp1)
          enddo
        enddo

C   ... Concatenate l+i and l+i+h blocks
        do  ik = 1, nkaph
          offH(maxp+1,ik,ibp1) = offH(maxp+1,1,ib)
     .    + isum(nkaph,ndimi(1,1),n0H)
     .    + isum(nkaph,ndimi(2,1),n0H)
     .    + isum(ik-1,ndimi(1,1),1)
     .    + isum(ik-1,ndimi(2,1),1)
          offH(maxp+2,ik,ibp1) = offH(maxp+2,1,ib)
     .    + isum(nkaph,ndimi(1,1),n0H)
     .    + isum(nkaph,ndimi(2,1),n0H)
     .    + isum(nkaph,ndimi(3,1),n0H)
     .    + isum(ik-1,ndimi(1,1),1)
     .    + isum(ik-1,ndimi(2,1),1)
     .    + isum(ik-1,ndimi(3,1),1)
        enddo
        do  ik = nkaph+1, nkap0
          do  iprior = maxp+1, maxp+2
            offH(maxp+1,ik,ib) = offH(maxp+1,nkaph,ib)
     .      + ndimi(1,nkaph) + ndimi(2,nkaph)
            offH(maxp+2,ik,ib) = offH(maxp+2,nkaph,ib)
     .      + ndimi(1,nkaph) + ndimi(2,nkaph) + ndimi(3,nkaph)
          enddo
        enddo
      enddo
C     This should never be needed
C     do  ik = 1, nkaph-1
C       do  iprior = 1, maxp
C         offH(iprior,ik+1,ib2+1) = offH(iprior,1,ib2+1)
C       enddo
C     enddo

C --- Make iprmb ---
  100 continue
      indx = 0
C     For now, require that morb = nl*nkaph
Changenglob      mxorb = nglob('mxorb')
      mxorb = globalvariables%mxorb
      i = nl**2*nkaph
Ckino isanrg is logical function,       call isanrg(mxorb,i,i,'makidx:','mxorb', .true.)
      l_dummy_isanrg=isanrg(mxorb,i,i,'makidx:','mxorb', .true.)
      call iinit(ndimi,n0H*nkap0)
      do  112  iprior = 1, maxp+1
C       iposn = 0
        do  110  jb = ib1, ib2
          if (lprmib) then
            ib = iprma(jb)
            ibp1 = iprma(jb+1)
          else
            ib = jb
            ibp1 = jb+1
          endif
          is = ips(ib)
          iposn = mxorb*(ib-1)

          i_copy_size=size(sspec(is)%idxdn)
          call i8icopy(i_copy_size,sspec(is)%idxdn,1,idxdn,1)

          do  111 ik = 1, nkaph
            do  111  l = 0, nl-1
              lx = iabs(idxpr(idxdn(l+1,ik))) .eq. iprior .or.
     .        iabs(idxpr(idxdn(l+1,ik))) .gt. iprior .and.
     .        iprior .gt. maxp
              do  120  m = -l, l
                iposn = iposn+1
                if (lx) then
                  indx = indx+1
                  iprmb(iposn) = indx
                  ndimi(iprior,ik) = ndimi(iprior,ik)+1
                endif
  120         continue
  111     continue
  110   continue
        ldham(iprior) = indx
  112 continue

C --- Printout ---
      if (ipr .lt. 10) return
      if (nkaph .eq. 1) then
        call awrit4('dm  L %i  I %i  H %i  N %i ',' ',80,stdl,ldham(1),
     .  ldham(2)-ldham(1),ldham(3)-ldham(2),ldham(4)-ldham(3))
      else
        call icopy(nkaph,ndimi(1,1),n0H,ltab,1)
        call awrit3('%xdm  L %i (%n:1i )',outs,80,0,ldham(1),nkaph,
     .  ltab)
        call icopy(nkaph,ndimi(2,1),n0H,ltab,1)
        call awrit3('%a  I %i (%n:1i )',outs,80,0,ldham(2)-ldham(1),
     .  nkaph,ltab)
        call icopy(nkaph,ndimi(3,1),n0H,ltab,1)
        call awrit4('%a  H %i (%n:1i )  N %i',outs,80,-stdl,
     .  ldham(3)-ldham(2),nkaph,ltab,ldham(4)-ldham(3))
      endif
      if (ipr .lt. 20) return
      call awrit5('%?;(n>=30);%N;; Makidx:  hamiltonian dimensions'//
     .' Low, Int, High, Negl: %i %i %i %i',' ',80,stdo,ipr,ldham(1),
     .ldham(2)-ldham(1),ldham(3)-ldham(2),ldham(4)-ldham(3))
      if (ipr .ge. 40 .or. (ipr .gt. 30 .and. nkaph .gt. 1)) then
        write(stdo,345)
  345   format(' kappa',3x,
     .  'Low',3x,'Int',3x,'High',2x,'L+I',2x,'L+I+H',2x,'Neglected')
        if (nkaph .gt. 1) then
          do  ik = 1, nkaph
            i = isum(3,ndimi(1,ik),1)
            write(stdo,346) ik, (ndimi(j,ik), j=1,3),
     .      ndimi(1,ik)+ndimi(2,ik),i,nl**2*(ib2-ib1+1)-i
  346       format(i4,2x,5i6,i8)
  347       format(1x,a5,5i6,i8)
          enddo
          strn = ' all '
        else
          strn = '  -  '
        endif
        write(stdo,347) strn, ldham(1),ldham(2)-
     .  ldham(1),ldham(3)-ldham(2),ldham(2),ldham(3),ldham(4)-
     .  ldham(3)
      endif

C      do  200  ib = ib1, ib2+1
C  200 print '(i5,5i5,2x,5i5,2x,5i5)', ib,(offH(j,1,ib), j=1,15)
C      i = 0
C      do  210  ib = ib1, ib2
C      do  210  ik = 1,  nkaph
C      do  210  m = 1,  nl**2
C        i = i+1
C        print '(6i5)', ib, ik, m, i, iprmb(i)
C  210 continue
C      call rx('done')

      end subroutine makidx

      subroutine mkidx2(nl,maxp,idxdn,ltab,ndimi,outs)

C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   maxp  :maximum number
Ci   idxdn :controls how orbital is to be included in the basis
Co Outputs
Co   ltab  :list of l's in permuted order
Co   ndimi :sum (ltab(i)+1)**2 for each priority
Co   outs
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n0,n0H
      parameter (n0=10,n0H=5)
      integer nl,maxp,ltab(n0),idxdn(n0),ndimi(n0H)
C ... Local parameters
      integer iprior,i2,i1,norb,iwk(n0),ndim,ipr,i,idxdn2(n0),idxpr
      character*80 outs
      character onam(9)*1,ostr*10
      data onam /'s','p','d','f','g','h','i','j','l'/

C     call icopy(nl,idxdn,1,iwk2,1)
      call getpr(ipr)
      ipr = 51
      do  10  i1 = 1, nl
        idxdn2(i1) = idxpr(idxdn(i1))
   10 continue
      call ivheap(1,nl,idxdn2,iwk,101)
      call iinit(ndimi,maxp)
      iprior = 0
      i2 = 0
   11 i1 = i2+1
      iprior = iprior+1
C ... We have exhausted this group of orbitals
      if (i1 .gt. nl .or. iprior .gt. maxp) return

      i2 = i1-1
      norb = 0
      ndim = 0
C ... Keep looping here until find all orbitals of this priority
   12 i2 = i2+1
      if (i2 .le. nl) then
        if (idxdn2(iwk(i2)) .eq. iprior) then
          norb = norb+1
          ndim = ndim + 2*(iwk(i2)-1) + 1
          ltab(norb) = iwk(i2)-1
          goto 12
        endif
      endif
      ndimi(iprior) = ndim
C     offH(iprior,1,ibp1) = offH(iprior,1,ibp1) + ndim
      if (ipr .ge. 50) then
        write(ostr,'(10a1)')(' ', i=norb+1,5),(onam(ltab(i)+1),i=1,norb)
        call awrit2('%a%b '//ostr(1:5)//' (%i)%(n*13+2)p%j+',
     .  outs,80,0,ndim,iprior)
      endif
      i2 = i2-1
      goto 11

      end subroutine mkidx2

      integer function idxpr(idxdn)

C- Given a value idxdn, assign priority to orbital (kernel called by makidx)
C ----------------------------------------------------------------------
Ci Inputs
Ci   idxdn :downfolding indices
Co Outputs
Co   idxpr :priority in the hamiltonian corresponding to idxdn
Cl Local variables
Cl         :
Cr Remarks
Cr   idxdn indicates how a basis function is to be used in the
Cr   construction of the basis;
Cr   idxpr uses idxdn to assign a priority according to what level
Cr   of treatment it has:
Cr   idxdn idxpr  Signifies
Cr     0    0     Role of this orbital has not yet been determined
Cr     1    1     Orbital is included as "active" orbitals, which means
Cr                they are included in the hamiltonian and diagonalized
Cr     2    2     Orbital is treated as an "intermediate" orbital, which
Cr                means it is downfolded and included in a perturbative way
Cr     3    3     Orbital is treated as an "high" orbital, which means
Cr                means it is included in augmentation for tails of other
Cr                orbitals, but is otherwise not part of the basis.
Cr     4    4     Orbital is completely neglected (no orbital is defined)
Cr    10    1     Orbital is a local orbital whose value and slope are
Cr                constructed to be zero at the MT boundary.
Cr                It is included in the basis.
Cr    11    1     Orbital is a local orbital with a smooth Hankel tail
Cr                and it is included in the basis.
Cr    12    2     Orbital is a local orbital with a smooth Hankel tail
Cr                It is incorporated perturbatively to the hamiltonian
Cr                and is not assembled as part of the hamiltonian matrix
Cr
Cu Updates
Cu   24 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer idxdn
C ... Local parameters

      if (idxdn .le. 4) then
        idxpr = idxdn
      elseif (idxdn .eq. 10 .or. idxdn .eq. 11) then
        idxpr = 1
      elseif (idxdn .eq. 12) then
        idxpr = 2
      else
        call rxi('makidx: bad idxdn',idxdn)
      endif

      end function idxpr


