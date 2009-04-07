      subroutine strscr(loka,nbas,npadl,npadr,alat,plat,bas,rwats,nl,
     .  kap,nkap,lmaxw,siinv,cy,cg,indxcg,jcg,lpgf,ldot,leqv,ntab,iax,
     .  alpha,adot,tral,trad,s,sdot)
C- Calculates 2nd gen or NMTO structure constants for all sites
C ----------------------------------------------------------------
Ci Inputs
Ci   loka:       conventions for Hankels and Bessel functions
Ci   nbas:       number of atoms in the basis
Ci   nl:         maximum l quantum number+1
Ci   alat,plat:  primitive lattice translation vectors, and scale
Ci   bas:        basis vectors
Ci   alpha:      screening parameters
Ci   adot:       energy derivative of screening parameters
Ci   tral:       transformation matrix for head and tail functions
Ci   trad:       (kappa*avw)^2-derivative of tral
Ci   cg,indxcg,jcg: Clebsch Gordon coefficients
Ci   cy          Normalization of YL
Ci   ntab,iax:   neighbor table information (see Remarks)
Ci   kap         kappa's for which to calc strux, scaled by avw
ci   nkap        number of kappa's
Ci   ldot        T, make sdot
Ci   lmaxw:      l-cutoff for Watson sphere.
Ci   rwats       Watson sphere radii
Ci   siinv:      parameters for iterative inversion of strux
Ci               1: number of sites for cutoff
Ci               2: number of iterations
Ci               3: allowed tolerance
Ci   lpgf:       T, strx for planar Green's function (see Remarks)
Ci   ldot:       T, make sdot
Ci   nbasp:      Number of elements in the padded basis,
Ci               usually nbas except when lpgf=T
Ci   leqv:  T,   look for equivalent strx
Cio  iax(*,i):   information concerning pair i connecting R,R'
Ci   iax(1,*)    site index for 'head' site R'
Ci   iax(2,*)    site index for 'augmentation' site R
Ci   iax(3..5)   number of primitive lattice translation vectors
Ci               separating the unit cells corresponding to R and R'.
Ci   iax(6,*)    index to conjugate (jb,ib) pair matching (ib,jb).
Ci   iax(7,*)    permutation table ordering cluster by increasing
Ci               effective site index (see pairc.f; not used here)
Co   iax(8,*)    in the case strux are equivalent, points to the
Co               site it is equivalent to.
Ci   iax(9,i)    holds (lmx+1)**2 for basis atom iax(2,i) (mkiaxd)
Co Outputs
Co   s           screened strux in (nl**2,nl*2) blocks for all pairs
Co               1..ntab(nbas+1)
Co   sdot        avw^-2 times energy derivative of s
Cr Remarks
Cr strscr generates screened structure constants and sdot (nkap=1) for
Cr a specified neighbor table (nbas,ntab,iax).  Structure constants are
Cr strung together in one long vector s of nl**2 x nl**2 blocks for each
Cr pair in the iax table. For one pair s is stored as S(R'L',RL), the
Cr augmentation being at R'.
Cr
Cr ntab(ib) contains the total number of pairs in the table before site
Cr ib; ntab(ib+1)-ntab(ib) the number of pairs in ib.  iax contains the
Cr geometrical information for the cluster around each ib.  R and R' are
Cr separated by some primitive translation vectors, plus some difference
Cr of basis vectors within the unit cell:
Cr
Cr Convention for s:  let ib,jb be iax(1,isite) and iax(2,isite).
Cr Then s(L',L,isite) are the strux for a one-center expansion
Cr around site jb of a screened screened Hankel centered at ib:
Cr   H_L(r-r_b) = sum_La s(La,Lb,isite) J_La'(r-r_a)
Cr Here 'a' the the augmentation index, 'b' is the basis.
Cr
Cr Mixed Cholesky decomposition/ Dyson expansion:
Cr structure constants are obtained as follows.  Definining the vector
Cr beta as the solution to
Cr   (alpha^-1 - S^0)_uv beta_v = delta_u0, where
Cr the central atom is column zero (u,v composite RL indices), this
Cr is done for an inner block by Cholesky decomposition.  The outer
Cr block and inner block are then converged in a Dyson expansion by
Cr the pair of equations
Cr   beta_u =  \sum _v (alpha S^0)_uv beta_v
Cr for u belonging to outer block, and v summed over both blocks, and
Cr   beta_u = (alpha_uu \delta_uv - S^0)_uu^-1 (1 + S^0_uv beta_v)
Cr for u the inner block and v the outer block.
Cr the latter solved again using Cholesky decomposed matrix.
Cr These equations are iterated, each iteration adding a correction to
Cr the approximate inverse.
Cr
Cr Structure constants for layer geometry:
Cr   alpha, adot, bas must be doubly padded by the addition of 
Cr   left and right bulk layers.
Cr
Cr Watson sphere:
Cr   See Andersen and Woolley, Mol. Phys. 26, 905 (1973).
Cr   A Watson sphere with radius rmx+deltr centered around the central
Cr   atom is introduced.  The region outside the Watson sphere is
Cr   called Omega.  Orbitals are introduced, which are regular outside
Cr   Omega, i.e., inside the Watson sphere; they are Bessel functions
Cr   with a l up to lmaxw.  These Bessel functions are expanded into
Cr   Bessel functions inside the spheres centered around all atoms in
Cr   the cluster.  Likewise the tails of the Neumann functions are
Cr   expanded into Neumann functions inside Omega.  Now screened
Cr   orbitals can be introduced which fulfill the same expansion rules
Cr   as the non-screened orbitals but with tails, which are now linear
Cr   combinations of Bessel and Hankel functions.  The coefficients of
Cr   mixing alpha and alphom are chosen so that all tails vanish at a
Cr   screening core radius sigma*wsr and on the Watson sphere.
Cr
Cr   Formally this is done by increasing the size of S^0 which now has
Cr   the structure:
Cr
Cr      ( S^0_{R',R}        S^0_{R',Omega}    )
Cr      (                                     )
Cr      ( S^0_{Omega,R}     S^0_{Omega,Omega} )
Cr
Cr   where S^0{R',R}         is the structure matrix Neumann-> Neumann
Cr         S^0_{R',Omega}    is the expansion matrix Bessel->Bessel
Cr         S^0_{Omega,R}     is the expansion matrix Neumann->Neumann
Cr         S^0_{Omega,Omega} is zero
Cr
Cu Updates
Cu   06 Aug 06 Redesigned so routines compatible with 2-kappa strux
Cu   19 May 04 make strux work with Methfessel conventions for H,J
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lpgf,leqv,ldot
      integer nbas,npadl,npadr,nl,nkap,ntab(nbas+1),niax,lmaxw,loka
      parameter (niax=10)
      integer iax(niax,1)
      double precision kap(nkap),s(nl*nl,nl*nl,1),sdot(nl*nl,nl*nl,1)
      double precision alpha(nl**2,*),adot(nl**2,*),
     .  tral(4,nl**2,*),trad(4,nl**2,*),
     .  alat,plat(9),bas(3,nbas),rwats(nbas),siinv(5)
C For Clebsch Gordon coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      logical cmdopt
      integer i,ik,lqinv,iat,ndimL,ndimW,nl2,ns,nsk,nitab,nitbck,
     .  isum,offik,nbasp,nbaspp,nsite0,nttab,nclus,iclus,streqv,stdo,
     .  ipr,nglob,nmto,nenv
      integer osadot,os0a,os0d,osalph,osflg,osflgd,oalphv,oadotv
      double precision xx
      character*80 outs
C heap:
      integer w(1)
      common /w/ w

C     call wkprnt(1)
      stdo = nglob('stdo')
      call getpr(ipr)
      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas
      nttab = ntab(nbasp+1)
      call rxx(.not.lpgf .and. nbas.ne.nbasp,'strscr: nbas ne nbasp')
      nl2 = nl**2
      nmto = 0
      if (nkap .gt. 1) nmto = 1

C --- Temporarily scale the lattice and basis vectors ---
      call dscal(9,alat,plat,1)
      call dscal(3*nbaspp,alat,bas,1)

C --- Determine order in which strux are to be calculated
C      if (.true.) then
C        call defi(oiwk,nttab)
C        call defdr(opos,3*nttab)
C        call siteid(iax,nbas,ntab,plat,bas,w(opos),w(oiwk),i)
C        print *, i
C        call rlse(opos)
C        call defi(osid,nttab)
C        call ordsit(iax,nbas,ntab,w(oiwk),w(osid))
C        stop
C      endif

C --- For each kappa, do ---
      nsk = 0
      do  20  ik = 1, nkap
C ... nitab is number of strux actually in the table
      nitab  = 0
      nitbck= 0
      nsite0 = 0

      if (ipr .gt. 40) write(stdo,'(1x)')

C --- Make the structure constants for all sites, this kappa ---
      offik = 1+(ik-1)*nbasp
      do  30  iat = 1, nbasp

        call info5(41,0,0,' strscr:  start strux for site %i'//
     .    '%?#(n>1)#, E*avw**2=%d##',
     .    iat,nkap,kap(ik),0,0)

C       Number of atoms in this cluster
        nclus = ntab(iat+1)-ntab(iat)
C       Offset to current cluster
        iclus = ntab(iat)+1
C       The true dimension of screened S-alpha
        ndimL = isum(nclus,iax(9,iclus),niax)
C       The dimension of strux matrix, including the Watson sphere
        ndimW = ndimL + (max(lmaxw,-1)+1)**2

C ...   Look for an already-calculated equivalent strux
        if (leqv) then
          if (nkap .gt. 1) call rx('check equiv strux for kap>1')
          i = streqv(ntab,ldot,nl2,alpha,adot,iat,0,plat,bas,iax,nitab)
          if (ipr .gt. 40 .and. i .gt. 0)
     .      call awrit2(' strscr:  map strux for site%,4i  from site'//
     .      '  %i',' ',80,stdo,iat,i)
          if (i .gt. 0) goto 30
        else
          do  32  i = iclus, ntab(iat+1)
          nitab = nitab+1
   32     iax(8,i) = nitab
        endif

C ...   Work arrays ...
        call defdr(osalph,nl2*ndimW)
        call defdr(oalphv, ndimW)
        if (ldot) then
          call defdr(oadotv, ndimW)
          call defdr(osadot,-nl2*ndimW)
        else
          oadotv = 1
        endif

        lqinv = 0
        if (cmdopt('--qinv',6,0,outs)) lqinv = 1
C       if (cmdopt('--qinv',6,0,outs) .or. siinv(1) .ne. 0) lqinv = 1
        call defdr(os0a,ndimW*ndimW)
        call defdr(os0d,ndimW*ndimW)

C ...   Make vector of alpha
        call mkalpv(alpha(1,offik),adot(1,offik),
     .    iax(1,iclus),0,nl2,nbasp,1,nclus,ndimW,w(oalphv))
        if (ldot) call mkalpv(alpha(1,offik),adot(1,offik),
     .    iax(1,iclus),2,nl2,nbasp,1,nclus,ndimW,w(oadotv))
        if (lmaxw .ge. 0) call alwats(ndimL,kap(ik),2,lmaxw,
     .    rwats(iat),w(oalphv),w(oadotv))

C ...   Calculate salpha and/or salpha-dot
C        call salph1(loka,lqinv,ldot,nbaspp,nl,siinv,ndimW,lmaxW,nclus,
C     .    plat,bas,w(oalphv),w(oadotv),tral(1,1,offik),trad(1,1,offik),
C     .    iax(1,iclus),cy,cg,indxcg,jcg,kap(ik),w(os0a),w(os0d),
C     .    w(osalph),w(osadot))
        nenv = iax(9,iclus)
        call salph1(nbaspp,nenv,lqinv,ldot,nmto*10+loka,siinv,ndimW,
     .    lmaxW,nclus,plat,bas,w(oalphv),w(oadotv),iax(1,iclus),cy,cg,
     .    indxcg,jcg,kap(ik),w(os0a),w(os0d),w(osalph),w(osadot))

C   ... NMTO: scale salpha,sdot according to tral matrix
        if (nmto .eq. 1) call scalsr(iax(1,iclus),ldot,ndimW,nl2,nclus,
     .    trad(1,1,offik),tral(1,1,offik),w(osalph),w(osadot))
C       call prmx('balpha after scaling',w(osalph),ndimW,ndimW,nenv)
C       call prmx('bdot after scaling',w(osadot),ndimW,ndimW,nenv)

C ...   Append salpha and/or salpha-dot to structure constant s
        ns = nitab-nclus
        call addtos(nl2,ndimW,iax(1,iclus),nl2,1,nclus,w(osalph),ns+nsk,
     .    s)
        if (ldot) then
          ns = nitab-nclus
          call addtos(nl2,ndimW,iax(1,iclus),nl2,1,nclus,w(osadot),ns+
     .      nsk,sdot)
        endif

        if (iat .le. nbas) nsite0 = nsite0 + nclus
        nitbck = nitbck + nclus
        call rlse(osalph)
C       call wksav(.false.)

C ... Reentry for case equivalent strx
   30 continue

C --- Sanity check ---
      if (leqv .and. nitbck .ne. nitab) call rx('bug in strscr')

C --- Symmetrize structure constants ---
      call defi(osflg,-nitab)
      call defi(osflgd,-nitab)
C      if (lpgf) then
C        call sympad(nl2,nbas,npadl,npadr,nsite0,nitab,iax,1,
C     .    w(osflg),s(1,1,nsk))
C        if (ldot) call sympad(nl2,nbas,npadl,npadr,nsite0,nitab,iax,1,
C     .    w(osflgd),sdot)
C      endif
      call symstr(0,nl2,nttab,iax,1,1,w(osflg),s(1,1,1+nsk),s,xx)
      if (ldot)
     .call symstr(0,nl2,nttab,iax,1,1,w(osflgd),sdot(1,1,1+nsk),sdot,xx)
      nsk = nsk+ns+nclus
   20 continue

C --- Undo scaling of the lattice and basis vectors ---
      call dscal(9,1/alat,plat,1)
      call dscal(3*nbaspp,1/alat,bas,1)

C --- Info printout ---
      if (ipr .ge. 30) call awrit2(' strscr:  generated %i '//
     .  'inequivalent strux from %i total',' ',80,stdo,nitab,nttab)

      end

