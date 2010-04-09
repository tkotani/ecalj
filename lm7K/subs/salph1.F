      subroutine salph1(nbas,nenv,lqinv,ldot,loka,siinv,ndimW,lmaxW,npr,
     .  plat,bas,alphv,adotv,iax,cy,cg,indxcg,jcg,el,s0a,s0d,
     .  balph,bdot)

C- Calculate one-center expansion coefficients for one cluster and energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nenv  :number of envelope function on head sphere; see Bugs
Ci   lqinv :switch governing how matrix inversion is accomplished
Ci          0 uses lapack dsytrf
Ci          1 uses dqinvb
Ci   ldot  :T, make energy derivatives of bdot also
Ci   loka: :1s digit old Andersen conventions for Hankel and Bessels
Ci         :10s digit => NMTO; return (alpha^-1 + S0)^-1 unscaled
Ci         :see scalsr.f for NMTO-scaling according to Tank conventions
Ci   siinv: parameters for iterative inversion of strux
Ci          1: number of sites for cutoff
Ci          2: number of iterations
Ci          3: allowed tolerance
Ci   ndimW :dimensions of strux, including space for Watson sphere
Ci   lmaxw :maximum l for Watson-sphere
Ci   npr   :number of pairs in cluster about this site
Ci   plat  :primitive lattice vectors, in atomic units
Ci   bas   :basis vectors, in atomic units
Ci   alphv :vector of alphas, of length ndimW
Ci   adotv :vector of adot/alpha**2, of length ndimW
Ci         :(ldot=F) not used
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   el    :energy of Hankels, atomic units
Ci         :NB: if loka=1, el should be scaled by avw**2
Ci   s0a   :a work array, holding unscreened strux B0, dim. ndimW*ndimW
Ci   s0d   :(ldot=T) a work array holding energy derivative of s0a
Ci         :(ldot=F) not used
Cl Local variables
Cl   ndimS :size of inner block when using iterative inversion.
Cl         :ndimS = ndimW when iterative inversion is not used.
Co Outputs
Co   balph :structure constants in a two-dimensional array ndimW*nenv
Co   bdot  :energy derivative of salpha * avw**-2
Co   alphv: is returned inverted
Cr Remarks
Cr   Returns screened strux B^alpha for a single energy
Cr     B^alpha = alp^-1 -  alp^-1 * (alp^-1 - s0)^-1 * alp^-1
Cr   and optionally its energy derivative.
Cr
Cr   Conventions and symbols in the LMTO literature have been various:
Cr      2nd generation notation refers to B^alpha as S^alpha
Cr      Tanks' implementation of Andersen's NMTO does the same.
Cr      (Phys. Stat. Sol. B217, 89 (2000))
Cr   However, the appropriate symbol is Balpha, since this matrix
Cr   represents the coefficients to the 1-center expansion of
Cr   screened functions around neighboring sites.
Cr
Cr   10s digit of loka is included for compatibility with Tank's
Cr   implentation of NMTO.  His definition requires an auxilliary
Cr   matrix called the 'tral' matrix; see mktra2.f
Cr   For compatibility with Tank's approach, this routine exits
Cr   after constructing (alp^-1 - s0)^-1.  The user then
Cr   scales this matrix with a call to scalsr:
Cr     call scalsr(iax,ldot,ndimW,nenv,npr,trad,tral,balph,bdot)
Cu Updates
Cu   06 Aug 06 Closer compatibility to 2-kappa case
Cu   07 Aug 06 Eliminated NMTO-specific scaling, for SSSW compatibility
Cu   19 May 04 make strux work with Methfessel conventions for H,J
Cu    9 Sep 00 iterative branch made to work for lqinv=f and ldot=t
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical ldot
      integer lqinv,nbas,ndimW,lmaxw,npr,niax,loka
      parameter (niax=10)
      integer iax(niax,npr)
      double precision el,bas(3,nbas),plat(3,3),alphv(ndimW),
     .  adotv(ndimW),s0a(ndimW,ndimW),s0d(ndimW,ndimW),siinv(5),
     .  balph(ndimW,*),bdot(ndimW,*)
C For Clebsch Gordon coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      integer j,i,ipr,nenv,ndwk,owk,owk2,owk3,os0i,orhs,lgunit,
     .  ndimS,lm2
      integer nblk,offs,offd,oipm,nito,ncut,niter
      double precision dsmx1,dsmx2,stol,stoll
      parameter (stoll=1d-8,nito=1)
#if DEBUG
      integer ndef,fopna
#endif
C heap:
      integer w(1)
      common /w/ w

C ... Setup and memory allocation
      ncut = siinv(1)
      niter = siinv(2)
      stol = siinv(3)
      ndwk = max(nenv,(lmaxw+1)**2)
      call defdr(owk,max(ndwk*ndwk,ndimW*nenv))

      call getpr(ipr)

      call tcn('salph1')

C     call prmx('alphv',alphv,ndimW,ndimW,1)

C ... make alpha^-1 + S0
      do  10  i = 1, ndimW
   10 alphv(i) = 1/alphv(i)

      i = 11
      if (ldot) i = 30
      if (mod(loka,10) .eq. 1) i = i+1200
      call defdr(owk,max(ndwk*ndwk,ndimW*nenv))
      call maks0(npr,plat,bas,indxcg,jcg,cg,cy,
     .  alphv,adotv,iax,el,i,lmaxw,s0a,s0d,ndimW)

C     call prmx('alpha^-1+S0',s0a,ndimW,ndimW,ndimW)
      call defi(offs,2*npr+2)
      call defi(offd,2*npr+2)
      call defi(oipm,2*npr+2)
      call pvsal0(ndimW,ncut,npr,iax,lmaxw,w(offs),w(offd),w(oipm),nblk,
     .  ndimS)

      if (ndimS .lt. ndimW) then
        call defrr(owk2,ndimW**2)
        call pmblk(nblk,w(oipm),w(offs),w(offd),ndimW,s0a,8100,0d0,0,
     .    ndimW,w(owk2),j)
        if (j .ne. ndimW) call rx('bug in pvsal0')
C        i = fopna('out',-1,4)
C        rewind i
C        print *, 's0 after permutation (binary)'
C        call ywrm(1,' ',1,i,' ',w(owk2),ndimw*ndimw,ndimw,ndimw,ndimw)
C        call fclose(i)
C        pause

        call dcopy(ndimW**2,w(owk2),1,s0a,1)
C       call prmx('maks0: permuted s0',s0,ndimW,ndimW,ndimW)
        if (ldot) then
          call pmblk(nblk,w(oipm),w(offs),w(offd),ndimW,s0d,8100,0d0,0,
     .      ndimW,w(owk2),j)
          if (j .ne. ndimW) call rx('bug in pvsal0')
          call dcopy(ndimW**2,w(owk2),1,s0d,1)
        endif
        call rlse(owk2)
      endif

      if (lqinv .eq. 0) then
        call defdr(owk2,ndimW)
      else
        call defrr(owk2, ndimS*(ndimS+1))
      endif
      if (ndimS .lt. ndimW .or. lqinv .ne. 0) then
        call defrr(os0i, ndimS**2)
        call defrr(orhs, ndimW*nenv)
      endif

C ... Add term delta_RL,R'L' to rhs
      call dpzero(balph,ndimW*nenv)
      do  15  j = 1, nenv
   15 balph(j,j) = 1

#if DEBUG
      i = fopna('out',-1,4)
      rewind i
c      print *, 's0 before inversion (binary)'
      call ywrm(1,' ',1,i,' ',s0a,ndimw*ndimw,ndimw,ndimw,ndimw)
      call fclose(i)
      pause
#endif

      call tcn('Invert S0')
C --- alpha^-1 - S0 for inner block ---
      if (lqinv .eq. 0 .and. ndimS .eq. ndimW) then
C       call dpack(s0a,ndimW)
C       call dspfa(s0a,ndimW,w(owk2),ndef)
C       do  17  j = 1, nenv
C  17   call dspsl(s0a,ndimW,w(owk2),balph(1,j))
C       call dsifa(s0a,ndimW,ndimW,w(owk2),i)
C       do  17  j = 1, nenv
C  17   call dsisl(s0a,ndimW,ndimW,w(owk2),balph(1,j))
        j = ndimW*64
        call defrr(owk3,j)
        call dsytrf('L',ndimW,s0a,ndimW,w(owk2),w(owk3),j,i)
        if (i .ne. 0) call rx('SALPH1: matrix singular')
        call dsytrs('L',ndimW,nenv,s0a,ndimW,w(owk2),balph,ndimW,i)
        call rlse(owk3)
      else
C       call prmx('s0',s0a,ndimW,ndimW,ndimW)
        call dmcpy(s0a,ndimW,1,w(os0i),ndimS,1,ndimS,ndimS)
        if (lqinv .eq. 0) then
          j = ndimW*64
          call defrr(owk3,j)
          call dsytrf('L',ndimS,w(os0i),ndimS,w(owk2),w(owk3),j,i)
          if (i .ne. 0) call rx('SALPH1: matrix singular')
          call dsytrs('L',ndimS,nenv,w(os0i),ndimS,w(owk2),balph,ndimW,
     .      i)
          call rlse(owk3)
        elseif (lqinv .eq. 1) then
          call dqinvb('s',w(os0i),ndimS,2,ndimS,nenv,w(owk2),ndimS,
     .      w(owk2),balph,ndimW,i)
        else
          call dqinv('s',w(os0i),ndimS,2,ndimS,w(owk2),ndimS,i)
          call dmcpy(w(os0i),ndimS,1,balph,ndimW,1,ndimS,nenv)
        endif
C       call prmx('s0 after inverse',w(os0i),ndimS,ndimS,ndimS)
      endif
      if (i .ne. 0) call rx('SALPH1: matrix singular')
      call tcx('Invert S0')

C --- Iterative solution for alpha^-1 - S0 ---
      if (ndimS .lt. ndimW) then
        call tcn('iterative Balpha')
C   ... No. Dyson iterative refinements before updating inner block
        call dpzero(w(orhs),nenv*nenv)
        call dcopy(nenv,1d0,0,w(orhs),nenv+1)
        do  20  i = 1, niter
        call pvsal1(ndimW,ndimS,nenv,s0a,nito,balph,w(owk),w(orhs),
     .    nenv,dsmx1)
        call tcn('dlnref')
        call dlnref(lqinv,0,ndimS,ndimW,nenv,w(os0i),ndimS,s0a,ndimW,
     .    w(orhs),nenv,'sb',w(owk2),balph,w(owk),dsmx2)
        call tcx('dlnref')
        call info5(50,0,0,' salph1 iter %i: ndimS=%i  ndimW=%i  '//
     .    'ds(S)=%,3;3g  ds(W)=%,3;3g',i,ndimS,ndimW,dsmx1,dsmx2)
        if (ipr .gt. 50) then
C        call awrit5(' salph1 iter %i: ndimS=%i  ndimW=%i  ds(S)=%,3;3g'
C     .    //'  ds(W)=%,3;3g',' ',80,lgunit(1)
        endif
        if (max(dsmx1,dsmx2) .lt. stoll) goto 21
   20   continue
   21   continue
        call info2(40,0,0,' salph1 converged strux to tol='//
     .    '%,3;3g in %i iter',max(dsmx1,dsmx2),i-1)
C        if (ipr .gt. 40) call awrit2(' salph1 converged strux to tol='//
C     .    '%,3;3g in %i iter',' ',80,lgunit(1),max(dsmx1,dsmx2),i-1)
        if (max(dsmx1,dsmx2) .gt. stol) then
          call fexit2(-1,111,' Exit -1 SALPH1: sought tol=%,3;3g but'
     .      //' reached only %,3;3g',stol,max(dsmx1,dsmx2))
        endif
        call tcx('iterative Balpha')
      endif
c     call prmx('balph',balph,ndimW,ndimW,nenv)

C --- Make screened s-dot ---
      if (ldot) then
        call tcn('make Bdot')
C   ... Make (-adot*a**-2 - Bdot^0) * (a^-1 - S^0)^-1
C       call prmx('s0d',s0d,ndimW,ndimW,ndimW)
        call dgemm('N','N',ndimW,nenv,ndimW,1d0,s0d,ndimW,balph,ndimW,
     .             0d0,bdot,ndimW)
C   ... - d/dE s =
C   ... (a^-1 - S^0)^-1 * (-adot a^-2 - Bdot^0) * (a^-1 - S^0)^-1
        if (lqinv .eq. 0 .and. ndimS .eq. ndimW) then
C         do  18  j = 1, nenv
C  18     call dspsl(s0a,ndimW,w(owk2),bdot(1,j))
C         do  18  j = 1, nenv
C  18     call dsisl(s0a,ndimW,ndimW,w(owk2),bdot(1,j))
          call dsytrs('L',ndimW,nenv,s0a,ndimW,w(owk2),bdot,ndimW,i)
        else
          if (ndimS .lt. ndimW) call dcopy(ndimW*nenv,bdot,1,w(orhs),1)
          if (lqinv .eq. 0) then
          call dsytrs('L',ndimS,nenv,w(os0i),ndimS,w(owk2),bdot,ndimW,i)
          else
          call dqinvb('sb',w(os0i),ndimS,2,ndimS,nenv,w(owk2),ndimS,
     .      w(owk2),bdot,ndimW,i)
          endif
        endif
C ...   Iterative improvement of bdot
        if (ndimS .lt. ndimW) then
C   ... No. Dyson iterative refinements before updating inner block
          do  32  j = 1, nenv
          do  32  i = ndimS+1, ndimW
            bdot(i,j) = 0
   32     continue
          do  30  i = 1, niter
          call pvsal1(ndimW,ndimS,nenv,s0a,nito,bdot,w(owk),w(orhs),
     .        ndimW,dsmx1)
          call tcn('dlnref')
          call dlnref(lqinv,0,ndimS,ndimW,nenv,w(os0i),ndimS,s0a,ndimW,
     .      w(orhs),ndimW,'sb',w(owk2),bdot,w(owk),dsmx2)
          call tcx('dlnref')
          if (ipr .gt. 50) then
          call awrit5(' salph1 iter %i: ndimS=%i  ndimW=%i'
     .        //'  ds(S)=%,3;3g  ds(W)=%,3;3g',' ',80,lgunit(1),i,ndimS,
     .        ndimW,dsmx1,dsmx2)
          endif
          if (max(dsmx1,dsmx2) .lt. stoll) goto 31
   30     continue
   31     continue
          call info2(41,0,0,' salph1 converged bdot to tol=%,3;3g'//
     .      ' in %i iter',max(dsmx1,dsmx2),i-1)
C          if (ipr .gt. 40)
C     .     call awrit2(' salph1 converged bdot to tol=%,3;3g in %i iter'
C     .      ,' ',80,lgunit(1),max(dsmx1,dsmx2),i-1)
          if (max(dsmx1,dsmx2) .gt. stol) then
            call fexit2(-1,111,' Exit -1 SALPH1: sought tol=%,3;3g but'
     .        //' reached only %,3;3g',stol,max(dsmx1,dsmx2))
          endif
        endif
C       call prmx('bdot',bdot,ndimW,ndimW,nenv)

C   ... Finally, bdot = d/dE s
        call dscal(ndimW*nenv,-1d0,bdot,1)

        call tcx('make Bdot')
      endif
      call rlse(owk)

C --- Undo permutation in balph, bdot ---
      if (ndimS .lt. ndimW) then
        call defrr(owk2,ndimW*nenv)
C       call prmx('balph (P)',balph,ndimW,ndimW,nenv)
        call pmblk(nblk,w(oipm),w(offd),w(offs),ndimW,balph,4200,0d0,
     .    nenv,ndimW,w(owk2),j)
        if (j .ne. ndimW) call rx('bug in pvsal0')
        call dcopy(ndimW*nenv,w(owk2),1,balph,1)
        if (ldot) then
          call pmblk(nblk,w(oipm),w(offd),w(offs),ndimW,bdot,4200,0d0,
     .    nenv,ndimW,w(owk2),j)
          call dcopy(ndimW*nenv,w(owk2),1,bdot,1)
        endif
        call rlse(owk2)
      endif

C      call prmx('balpha before scaling',balph,ndimW,ndimW,nenv)
C     call prmx('bdot before scaling',bdot,ndimW,ndimW,nenv)

C --- NMTO : scaling is done externally; see scalsr.f ---
      if (loka/10 .eq. 1) return

C --- Balpha = alp^-1 * (alp^-1 - s0)^-1 * alp^-1  and ---
C     Bdot <=  alp^-1 * [bdot] * alp^-1 -
C              adot * alp^-1 * Balpha  +  Balpha * alp^-1 * adot ---
      do  lm2 = 1, nenv
        do  i = 1, ndimW
          balph(i,lm2) = alphv(i)*balph(i,lm2)*alphv(lm2)
        enddo
        if (ldot) then
        do  i = 1, ndimW
          bdot(i,lm2)  = alphv(i)*bdot(i,lm2)*alphv(lm2)
          bdot(i,lm2)  = bdot(i,lm2) + balph(i,lm2)*
     .                   (adotv(i)/alphv(i) + adotv(lm2)/alphv(lm2))
        enddo
        endif
      enddo

C --- Balpha = alp^-1 - alp^-1 * (alp^-1 - s0)^-1 * alp^-1
C     Bdot <=  adot * alp^-2 - [bdot]
      call dscal(ndimW*nenv,-1d0,balph,1)
      if (ldot) call dscal(ndimW*nenv,-1d0,bdot,1)
      do  lm2 = 1, nenv
        balph(lm2,lm2) = balph(lm2,lm2) + alphv(lm2)
      enddo
      if (ldot) then
      do  lm2 = 1, nenv
        bdot(lm2,lm2) = bdot(lm2,lm2) + adotv(lm2)
      enddo
      endif

C      call prmx('balpha after scaling',balph,ndimW,ndimW,nenv)
C      call prmx('bdot after scaling',bdot,ndimW,ndimW,nenv)

C --- Andersen 2nd gen conventions are -B ---
      if (mod(loka,10) .eq. 1) then
        call dscal(ndimW*nenv,-1d0,balph,1)
      endif

      call tcx('salph1')
Cgetarg       end
       end subroutine salph1 

      subroutine pvsal0(ndimW,ncut,npr,iax,lmaxw,
     .  offs,offd,ipm,nblk,ndimS)

C- Set up subblock lists and permutation arrays to permute strux.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimW :
Ci   ncut  :number of sites inside sphere where full information
Ci         :is kept.
Ci   npr   :total number of sites
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci         :iax(9,i) = number of orbitals associated w/ site i
Ci   lmaxw :maximum l for Watson-sphere
Co Outputs
Co   ... See Remarks for explanation
Co       Here,  n_i = # orbitals at site i
Co   offs  :offs(1) = 0
Co         :offs(2) = total number of orbitals of sites inside ncut
Co         :        = sum_i<ncut n_i
Co         :offs(2i+1) = offs(2) + sum_j<i n_(i+ncut) + 1
Co         :offs(2i+2) = offs(2) + sum_j<=i n_(i+ncut)
Co   offd  :offd(1) = 0
Co         :offd(2) = offs(2)
Co         :offd(i+2) = offd(2) + sum_j<i (1)
Co   ipm   :ipm(1) = 1
Co         :ipm(2) = 2
Co         :ipm(i) = 2*i+2
Co   ndimS :source dimension
Co   nblk  :Number of permutation blocks
Cl Local variables
Cl         :
Cr Remarks
Cr   This routine assembles permutation pointers to reorder the
Cr   orbitals in the strux array into all orbitals within the inner
Cr   inner block, followed by just the s orbitals belonging to sites
Cr   in the outer block, followed by Watson sphere orbitals, followed
Cr   finally by the remaining orbitals in the outer block.  These
Cr   last orbitals are not included in the direct inversion of s;
Cr   their contribution to the inverse of s is obtained iteratively.
Cr
Cr   Thus, offs contains a sequence of offsets that partition orbitals
Cr   at each site in the outer block into two blocks: the block for the
Cr   s orbital, followed by a block for the remaining orbitals.
Cr   The first entries in offs encompass the entire inner block.
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ndimW,lmaxw,ncut,npr,niax,ndimS,nblk
      parameter (niax=10)
      integer iax(niax,npr),offs(2*npr+2),offd(2*npr+2),ipm(2*npr+2)
C ... Local parameters
      integer i,nout,n,noutw,iprint,lgunit

C ... Setup: count number of permutation blocks
      ndimS = ndimW
      if (ncut .eq. 0) return
      nout = npr-ncut
      noutw = nout
      if (lmaxw .ge. 0) noutw = noutw+1

C ... Make the inner block offsets
      offs(1) = 0
      offs(2) = 0
      do  10  i = 1, ncut
   10 offs(2) = offs(2) + iax(9,i)
      offd(1) = 0
      offd(2) = offs(2)
      ipm(1) = 1
      ipm(2) = 2
      n = 2
      nblk = 2
      do  12  i = 1, nout
        offs(n+1) = offs(n) + min(1,iax(9,ncut+i))
        offs(n+2) = offs(n) + iax(9,ncut+i)
        offd(nblk+1) = offd(nblk) + min(1,iax(9,ncut+i))
        n = n+2
        nblk = nblk+1
        ipm(nblk) = n
   12 continue
      n = n-1
      nblk = nblk-1
      if (lmaxw .ge. 0) then
        n = n+1
        nblk = nblk+1
        offs(n+1) = offs(n) + (lmaxw+1)**2
        offd(nblk+1) = offd(nblk) + (lmaxw+1)**2
        ipm(nblk) = n
      endif
      ndimS = offd(nblk+1)
      n = 1
      do  16  i = 1, nout
        nblk = nblk+1
        n = n+2
        offd(nblk+1) = offd(nblk) + max(0,iax(9,ncut+i)-1)
        ipm(nblk) = n
   16 continue

C      call awrit3('offs %i: %n:1i',' ',250,6,nblk,nblk+1,offs)
C      call awrit3('offd %i: %n:1i',' ',250,6,nblk,nblk+1,offd)
C      call awrit3('ipm  %i: %n:1i',' ',250,6,nblk,nblk,ipm)

      if (iprint() .ge. 50) call awrit3(' pvsal0: ncut = %i  ndimW = %i'
     .  //'  ndimS = %i',' ',80,lgunit(1),ncut,ndimW,ndimS)

Cgetarg       end
       end subroutine pvsal0 

      subroutine pvsal1(ndimW,ndimS,nenv,s0,niter,balph,dsalp,b,nb,dsmx)

C- Gauss-Seidel iteration towards iterative inversion
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimW :dimensions of strux, including space for Watson sphere
Ci   ndimS :size of inner block for which salpha is not to be updated
Ci   nenv  :column dimension of balph
Ci   s0    : alpha^-1 - s0
Ci   balph : estimate for salpha = (alpha^-1 - s0)^-1 * b (nenv columns)
C    b,nb  :right-hand side for (alpha^-1 - s0) * salpha = b, and dimension
Ci   niter :number of iterations
Co Outputs
Co   balph : updated estimate in subblock (ndimS+1..ndimW,1..nenv)
Co   dsalp : work array holding changes in salpha
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer ndimW,ndimS,nenv,nb
      double precision s0(ndimW,ndimW),balph(ndimW,nenv),
     .  dsalp(ndimW,nenv),b(nb,nenv),dsmx,dsmxi
      integer i1,i2,ii,i,j,nrow,mrow,iter,niter
      parameter (nrow=54)

C --- Estimate ndimS+1,ndimW block by matrix Gauss-Seidel ---
      do  10  iter = 1, niter
      dsmx = 0d0
      do  12  i1 = ndimS+1, ndimW, nrow
        mrow = min(ndimW-i1+1,nrow)
        i2 = i1+mrow-1

C   ... (i1,i2) subblock of (a^-1 - s0) (trial alpha) - b
        call dgemm('N','N',i2-i1+1,nenv,ndimW,1d0,s0(i1,1),
     .    ndimW,balph,ndimW,0d0,dsalp(i1,1),ndimW)
        if (i1 .le. nb) then
          ii = min(nb,i2)
          do  20  j = 1, nenv
          do  20  i = i1, ii
            dsalp(i,j) = dsalp(i,j) - b(i,j)
  20     continue
        endif
        dsmxi = 0d0
        do  30  j = 1, nenv
        do  30  i = i1, i2
          dsalp(i,j) = dsalp(i,j) / s0(i,i)
          dsmxi = max(dsmxi,abs(dsalp(i,j)))
C         For overrelaxation
C         q = max(q,abs(dsalp(i,j)/dsalp(i,j)))
   30   continue

c       call prmx('dsalph',dsalp,ndimW,ndimW,nenv)

        do  40  j = 1, nenv
        do  40  i = i1, i2
          balph(i,j) = balph(i,j) - dsalp(i,j)
  40   continue
        dsmx = max(dsmx,dsmxi)

   12 continue

C      call awrit2(' pvsal1: iter %i  max ds=%,3;3g',
C     .  ' ',80,6,iter,dsmx)

C      Overrelaxation w
C      if (q .lt. 1 .and. l .eq. 2) then
C        w = 2/(1 + sqrt(1 - ((q+w-1)/w)**2/q))
C      endif

   10 continue
C     call prmx('balph after niter',balph,ndimW,ndimW,nenv)
Cgetarg       end
       end subroutine pvsal1 


