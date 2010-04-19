      subroutine strrs(nbas,npadl,npadr,alat,plat,bas,rwats,
     .  nl,el,nel,lmaxW,siinv,cy,cg,indxcg,jcg,lpgf,ldot,leqv,loka,
     .  ntab,iax,alpha,adot,s,sdot)
C- Structure matrix for screened hankels, one site
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas:       number of atoms in the basis
Ci   nl:         maximum l quantum number+1
Ci   alat,plat:  primitive lattice translation vectors, and scale
Ci   bas:        basis vectors
Ci   alpha:      screening parameters
Ci   adot:       energy derivative of screening parameters
Ci   cg,indxcg,jcg: Clebsch Gordon coefficients
Ci   cy          Normalization of YL
Ci   ntab,iax:   neighbor table information (see Remarks)
ci   nel         number of LMTO energies
Ci   el          LMTO energies for which to calc strux, in a.u.
Ci   ldot        T, make sdot
Ci   lmaxW:      l-cutoff for Watson sphere.
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
Ci   loka:       conventions for Hankels and Bessel functions
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
Ci   iax(9,*)    holds (lmx+1)**2 for basis atom of pair i
Cl Local variables
Cl   nclus       size of current cluster
Cl   nenv        number of envelope function on head sphere in current
Cl               cluster.  For now nenv=nl**2
Cb Bugs
Cb   This routine should have with the following distinct dimensions
Cb   ndimL:total number of envelope functions in cluster
Cb   ndimW:ndimC + number of functions in Watson sphere
Cb   ndimS:number of envelope functions in cluster to be determined
Cb        :by direct inversion
Cb   nenv :number of envelope functions on head sphere.
Cb
Cb   At the moment, makalp and makalv, salph2 assume orbitals 1..nl**2
Cb   are included for each site.       
Cb                                     
Cb   ldot = T is implemented (but not tested) only for nel=1 
Co Outputs
Co   s           screened strux in (nl**2,nl*2) blocks for all pairs
Co               1..ntab(nbas+1)
Co   sdot        avw^-2 times energy derivative of s
Cr Remarks
Cr strrs generates screened structure constants and optionally sdot for
Cr a single site.  Structure constants are stored as a vector s 
Cr of blocks dimensioned nds**2nel**2, one block for each
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
Cr   with a l up to lmaxW.  These Bessel functions are expanded into
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
Cu   23 Jul 08 (S Lozovoi) Updated to work with species-dependent lmax
Cu   03 Aug 06 Adapted from strscr
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lpgf,leqv,ldot
      integer nbas,npadl,npadr,nl,nel,ntab(nbas+1),niax,lmaxW,loka
      parameter (niax=10)
      integer iax(niax,1)
      double precision el(nel),s(nl*nl,nl*nl,1),sdot(nl*nl,nl*nl,1)
      double precision alpha(nl**2,nbas,nel,nel),adot(nl**2,*),
     .  alat,plat(9),bas(3,nbas),rwats(nbas),siinv(5)
C For Clebsch Gordon coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      logical cmdopt
      integer i,lqinv,iat,ndimL,ndimW,nl2,ns,nitab,nitbck,
     .  isum,nbasp,nbaspp,nsite0,nttab,nclus,iclus,streqv,stdo,
     .  ipr,nglob,nds,nenv
      integer osadot,os0a,os0d,osalph,osflg,osflgd,oalphv,oadotv,os01,
     .  os02,okpvt,owk
      double precision xx
      character*80 outs
C heap:
      integer w(1)
      common /w/ w

C     call wkprnt(1)
      call tcn('strrs')
      stdo = nglob('stdo')
      call getpr(ipr)
      if (nel .eq. 1) then
        if (ldot) call info0(10,0,0,
     .             ' strrs: Warning! sdot were not tested for SSSW')
      elseif (nel .eq. 2) then
        if (ldot)  
     .    call rx(' strrs: sdot are not implemented for 2-kappa SSSW')
      else
        call rxi(' strrs: nel should be either 1 or 2. nel = %i',nel) 
      endif
      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas
      nttab = ntab(nbasp+1)
      call rxx(.not.lpgf .and. nbas.ne.nbasp,'strscr: nbas ne nbasp')
      nl2 = nl**2
C     Dimension of s
      nds = nl2

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
C ... nitab is number of strux actually in the table
      nitab  = 0
      nitbck= 0
      nsite0 = 0

C --- Make the structure constants for all sites, this kappa ---
      do  iat = 1, nbasp

C       Number of atoms in this cluster
        nclus = ntab(iat+1)-ntab(iat)
C       Offset to current cluster
        iclus = ntab(iat)+1
C       The true dimension of screened S-alpha
        ndimL = isum(nclus,iax(9,iclus),niax)
C       The dimension of strux matrix, including the Watson sphere
        ndimW = ndimL + (max(lmaxW,-1)+1)**2
C       Number of head functions in this sphere
        nenv = iax(9,iclus)

C ...   Look for an already-calculated equivalent strux
        if (leqv) then
          call rx('leqv should be moved out of strrs ... stop for now')
C          if (nel .gt. 1) call rx('check equiv strux for el>1')
C          i = streqv(ntab,ldot,nl2,alpha,adot,iat,0,plat,bas,iax,nitab)
C          if (ipr .gt. 40 .and. i .gt. 0)
C     .      call awrit2(' strscr:  map strux for site%,4i  from site'//
C     .      '  %i',' ',80,stdo,iat,i)
C          if (i .gt. 0) goto 30
        else
          do  i = iclus, ntab(iat+1)
            nitab = nitab+1
            iax(8,i) = nitab
          enddo
        endif

C ...   Work arrays ...
        call defdr(osalph,-nl2*ndimW*nel**2)
        call defdr(oalphv, ndimW*nel**2)
        if (ldot) then
          call defdr(oadotv, ndimW*nel**2)
          call defdr(osadot,-nl2*ndimW*nel**2)
        else
          oadotv = 1
          osadot = 1
        endif

        lqinv = 0
        if (cmdopt('--qinv',6,0,outs)) lqinv = 1
C       if (cmdopt('--qinv',6,0,outs) .or. siinv(1) .ne. 0) lqinv = 1
        call defdr(os0a,ndimW*ndimW)
        call defdr(os0d,ndimW*ndimW)

C ...   Make vector of alpha^-1
        if (nel .eq. 1) then
          call mkalpv(alpha,adot,iax(1,iclus),0,nds,nbas,nel,nclus,
     .      ndimW,w(oalphv))
          if (ldot) call mkalpv(alpha,adot,iax(1,iclus),2,nds,nbas,nel,
     .      nclus,ndimW,w(oadotv))
          if (lmaxw .ge. 0) call alwats(ndimL,el,2,lmaxw,
     .      rwats(iat),w(oalphv),w(oadotv))
        else
          call mkalpv(alpha,adot,iax(1,iclus),1,nds,nbas,nel,nclus,
     .      ndimW,w(oalphv))
          if (lmaxw .ge. 0)
     .      call rx('watson sphere with 2-kap not implemented')
        endif

C ...   Calculate b and/or b-dot
        if (nel .eq. 1) then
          call salph1(nbaspp,nenv,lqinv,ldot,loka,siinv,ndimW,lmaxW,
     .      nclus,plat,bas,w(oalphv),w(oadotv),iax(1,iclus),cy,cg,
     .      indxcg,jcg,el,w(os0a),w(os0d),w(osalph),w(osadot))
        else
C         These next two for Dyson expansion in L (not implemented)
C         call defdr(os01,-ndimL*(ndimL-ndimS))
C         call defdr(os02,(ndimL-ndimS)*nl2)
          call defdr(os01,1)
          call defdr(os02,1)
          call defi(okpvt,ndimW*nel)
          call defdr(owk,ndimW*2 + ndimW*4*nl2)
          call salph2(nbaspp,nenv,lqinv,ldot,loka,siinv,ndimW,lmaxW,
     .      nl2,nclus,plat,bas,w(oalphv),w(oadotv),iax(1,iclus),cy,cg,
     .      indxcg,jcg,el,w(os0a),w(os0d),w(os01),w(os02),w(owk),
     .      w(okpvt),w(osalph),w(osadot))
        endif
        call rlse(os0a)

C ...   Append salpha and/or salpha-dot to structure constant s
        ns = nitab-nclus
        call addtos(nds,ndimW,iax(1,iclus),nl2,nel,nclus,w(osalph),ns,s)
        if (ldot) then
          ns = nitab-nclus
          call addtos(nds,ndimW,iax(1,iclus),nl2,nel,nclus,w(osadot),
     .      ns,sdot)
        endif

        if (iat .le. nbas) nsite0 = nsite0 + nclus
        nitbck = nitbck + nclus
        call rlse(osalph)
C       call wksav(.false.)

C   ... Reentry for case equivalent strx
C   30   continue
      enddo

C --- Sanity check ---
C     if (leqv .and. nitbck .ne. nitab) call rx('bug in strrs')

      if (.not. cmdopt('--nosyx',7,0,outs) .and. nel .eq. 1) then
C ---   Symmetrize structure constants ---
        call defi(osflg,-nitab)
        call defi(osflgd,-nitab)
C        if (lpgf) then
C          call sympad(nl2,nbas,npadl,npadr,nsite0,nitab,iax,1,
C     .      w(osflg),s)
C          if (ldot) call sympad(nl2,nbas,npadl,npadr,nsite0,nitab,iax,1,
C     .      w(osflgd),sdot)
C        endif
        xx = 0
        call symstr(0,nl2,nttab,iax,1,nel,w(osflg),s(1,1,1),s,xx)
          call info2(30,0,0,' strrs: strux are symmetrized;'//
     .                      ' max error %g',xx,0)
        if (ldot) then 
          call symstr(0,nl2,nttab,iax,1,nel,w(osflgd),sdot(1,1,1),
     .      sdot,xx)
          call info2(30,0,0,' strrs: s-dot are symmetrized;'//
     .                      ' max error %g',xx,0)
        endif
      endif

C --- Undo scaling of the lattice and basis vectors ---
      call dscal(9,1/alat,plat,1)
      call dscal(3*nbaspp,1/alat,bas,1)

C --- Info printout ---
      if (ipr .ge. 30) call awrit2(' strrs:  generated %i '//
     .  'inequivalent strux from %i total',' ',80,stdo,nitab,nttab)

      call tcx('strrs')

      end

