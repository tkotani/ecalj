      subroutine rothrm(opt,ldmpa,ldima,iprmb,rotm,ib1,ib2,ldU,U,Uz)
C- Generate rotation matrix for hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :  0 use real harmonics; output U is real
Ci         :  1 use true spherical harmonics; output U is complex
Ci         :  2 use real harmonics, but return U as complex
Ci   ldmpa :offset to first orbital in downfolding block of h to rotate
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   rotm  :3x3 cartesian rotation matrix
Ci   ib1   :U  generated for (ib1..ib2) site subblock of h
Ci   ib2   : -//-
Ci   ldU  :leading dimension of U
Co Outputs
Co   U,Uz :Rotation matrix for (ib1..ib2) subblock of h
Co        :U is the real form (opt=0); Uz is the complex form.
Co        :Only of of (U,Uz) is returned.
Cl Local variables
Cl   rYL  :rotation matrix that rotates Y_lm
Cr Remarks
Cr   This routine generates rotation matrix U for the hamiltonian matrix
Cr   which consists of blocks which transform in the row and colum indices
Cr   as Y_lm.
Cr
Cr   The unrotated and rotated hamiltonians h and hr are related by :
Cr      hr = U h U+
Cr
Cr   See also rothph, below
Cb Bugs
Cb   This routine has only been checked for ib1=1 and ib2=nbas
Cu Updates
Cu   17 Dec 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer opt,ldU,iprmb(1),ldmpa,ldima,ib1,ib2
      double precision rotm(3,3),U(ldU,ldU)
      double complex Uz(ldU,ldU)
C ... Local parameters
      integer n0,nkap0,nlx,nl2,i,ibas,iorb,j,li,
     .  nlm,nlmi,norb,offri,offi
      parameter (n0=10,nkap0=3,nlx=9,nl2=nlx*nlx)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      double precision xx,rYL(nl2,nl2)
      double complex rYLz(nl2,nl2)

      if (opt .ne. 0 .and. opt .ne. 2)
     .  call rx('rothrm not ready for this opt')

C ... Rotation matrix for spherical harmonics
      call ylmrtg(nl2,rotm,rYL)
C     call prmx('rothrm: rYL',rYL,nl2,nl2,nl2)

      if (opt .eq. 1) then
        call ylmrtg(nl2,rotm,rYL)
        call s2sph(1,nlx,nlx,rYL,nl2,nl2,nl2,nl2,rYLz)
C       call yprm('rYL',3,rYL,nl2*nl2,nl2,nl2,nl2)
      endif

C ... Initialize U
      if (opt .eq. 0) then
        call dpzero(U,ldU*ldU)
      else
        call dpzero(Uz,ldU*ldU*2)
      endif

C --- For each (ibas,l) block, poke rotation matrix rotating Ylm(block)
      do  ibas = ib1, ib2

C   ... Row, column offsets to h
        call orbl(ibas,ldmpa,ldima,iprmb,norb,ltab,ktab,xx,offl,nlm)

C   --- For each (ibas,l) block, poke rYL(i,j) into U(i,j) ---
C   ... loop over pairs of orbital blocks
        do  iorb = 1, norb
          li = ltab(iorb)
          nlmi = 2*li + 1
          offi = offl(iorb) - ldmpa
          offri = li**2
          if (li+1 .gt. nlx)
     .      call rxi('rothrm: increase nlx, need nlx=',(li+1)**2)

C     ... loop over pairs of m channels within this block
          if (opt .eq. 0) then
            do  j = 1, nlmi
            do  i = 1, nlmi
              U(i+offi,j+offi) = rYL(i+offri,j+offri)
            enddo
            enddo
          elseif (opt .eq. 1) then
            do  j = 1, nlmi
            do  i = 1, nlmi
              Uz(i+offi,j+offi) = rYLz(i+offri,j+offri)
            enddo
            enddo
          else
            do  j = 1, nlmi
            do  i = 1, nlmi
              Uz(i+offi,j+offi) = rYL(i+offri,j+offri)
            enddo
            enddo
          endif

        enddo

      enddo

C      if (opt .eq. 0) then
C        call prmx('U',U,ldU,ldU,ldU)
C      else
C        call zprm('Uz',2,Uz,ldU,ldU,ldU)
C      endif
      end

      subroutine rothph(opt,q,tbas,ldima,iprmb,ib1,ib2,ldH,H)
C- Extra scaling of hamiltonian by a phase
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0 scale on left by phase
Ci         :1 scale on right by conjugate of phase
Ci         :2 scale both left and right
Ci         :10s digit
Ci         :0 phase = exp(2*pi*tbas*q)
Ci         :1 phase = exp(-2*pi*tbas*q)
Ci   q     :k-point for current h
Ci         :Note: q and tbas must be in reciprocal units
Ci   tbas  :set of basis vector shifts; see Remarks
Ci         :Note: q and tbas must be in reciprocal units
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ib1   :H  generated for (ib1..ib2) site subblock of h
Ci   ib2   : -//-
Ci   ldh  :leading dimension of H
Co Outputs
Co   H    :H is scaled by a q-dependent phase; see Remarks
Cl Local variables
Cr Remarks
Cr   when the hamiltonian is constructed by Ewald sums, and a basis
Cr   vector is shortened by dpos, the hamiltonin is scaled by
Cr   a q-dependent phase.
Cr   hamiltonian rows belonging to that site are scaled by
Cr      phase  = exp(2*pi*i * q * dpos)
Cr   hamiltonian columns belonging to that site are scaled by
Cr      phase  = exp(-2*pi*i * q * dpos)
Cb Bugs
Cb   This routine has only been checked for ib1=1 and ib2=nbas
Cb   This routine assumes hamiltonian is contiguous block betw/ (ib1..ib2)
Cb   This routine uses off from orbl and assigns ldimp=0.
Cb     These assume only lower block is being rotated
Cu Updates
Cu   17 Dec 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer opt,ldH,iprmb(1),ldima,ib1,ib2
      double precision q(3),tbas(3,ib2)
      double complex H(ldH,ldH)
C ... Local parameters
      logical nphase
      integer n0,nkap0,i,ib,j,j1,j2,ndimi,norb,off,ipr,PRDBG
      double precision twopi,sp
      double complex phase
      parameter (n0=10,nkap0=3,PRDBG=80)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)

      if (q(1) .eq. 0 .and. q(2) .eq. 0 .and. q(3) .eq. 0) return

      twopi = 8d0*datan(1d0)
      nphase = .false.
      call getpr(ipr)
      if (mod(opt/10,10) .ne. 0) nphase = .true.

C ... Offsets to starting and ending regions of h(ib1..ib2)
      call orbl(ib1,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
      j1 = off+1
      call orbl(ib2,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
      j2 = off + ndimi

C --- For each site ibas, construct poke rotation matrix rotating Ylm(ibas)
      do  ib = ib1, ib2

C   ... Row, column offsets to h
        call orbl(ib,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)

C   ... Phase
        sp = twopi*(tbas(1,ib)*q(1) + tbas(2,ib)*q(2) + tbas(3,ib)*q(3))
        if (nphase) sp = -sp
        phase = dcmplx(dcos(sp),dsin(sp))
        if (ipr .ge. PRDBG) then
          print 333, ib,sp*360/twopi,q
  333     format(' rothph: ib=',i4,'  phase=',f7.1,' deg.  q=',3f7.1)
        endif

C   --- For all orbitals belonging to site ib, h->phase*h and/or h->h*phase+
        if (mod(opt,10) .eq. 0 .or. mod(opt,10) .eq. 2) then
          do  j = j1, j2
            do  i = off+1, off+ndimi
              h(i,j) = phase*h(i,j)
            enddo
          enddo
        endif

        if (mod(opt,10) .eq. 1 .or. mod(opt,10) .eq. 2) then
          phase = dconjg(phase)
          do  i = off+1, off+ndimi
            do  j = j1, j2
              h(j,i) = h(j,i)*phase
            enddo
          enddo
        endif

      enddo

      end

