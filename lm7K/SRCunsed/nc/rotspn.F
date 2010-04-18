      subroutine rotspn(mode,nbas,nl,indxsh,eula,neul,theta,d1,d2,ldd,
     .  ldima,ldimb,lds,lds2,s,sr)
C- Returns a structure matrix or set of vectors rotated by spinor
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit pertains to storage of hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (not implemented)
Ci          10s digit distinguishes how complex arithmetic is handled.
Ci           0: s,sr have real, imaginary separated, eg
Ci              s = s0(lds,lds2,2), with s0(*,*,1..2) = real..imag
Ci           1: s,sr are in complex*16 format, i.e. for real s
Ci              s = s1(2,lds,lds2), with s1(1..2,*,*) = real..imag
Ci           2: s,sr  have real, imaginary separated by columns, eg
Ci              s = s2(lds,2,lds2), with s2(*,1..2,*) = real..imag
Ci        100s digit
Ci           1 if to reverse the sense of rotation, i.e. U -> U+
Ci       1000s digit (unused when sr is rotated as a set of vectors)
Ci           1 Add d1 to diagonal (diagonal in RL and spin)
Ci           2 Scale sr_RLs,R'L's' by d2(RLs) d2(R'L's')
Ci             or when rotating vectors, scale vectors sr_RLs,* by d2(RLs)
Ci             before rotation
Ci           3 combination of 1+2
Ci           4 Add d1 to diagonal (diagonal in RL and spin),
Ci             plus rel. off diag, i.e. d_RL,RL+1 for spin 12 block
Ci                                 and  d_RL+1,RL for spin 21 block
Ci      10000s digit (options are mutually exclusive)
Ci           1 if to add to sr (by default sr is initialized to zero)
Ci           2 Input s is a spin spiral held in s(k-q/2) and s(k+q/2)
Ci             Use this mode to transforms ll block of s
Ci           3 Input is already spinor stored in sr; s is not used.
Ci           4 Input is a set of vectors stored in sr, and sr is
Ci             overwritten by as U sr.  s,d1,d2 are not used.
Ci           5 (SS) Make spin spiral il block from s(k-q/2) and s(k+q/2)
Ci             Only Sss is made; no subsequent rotations or scalings
Ci     100000s digit (options are mutually exclusive)
Ci             (not applicable when sr is rotated as a set of vectors)
Ci           0 s is not diagonal
Ci           1 s is diagonal in R but not L (not implemented)
Ci           2 s is diagonal in R and L; then lds2 should be 1
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   eula  :Euler angles for specifying spinor
Ci   neul  :1  if Euler angles are l-independent
Ci          nl if Euler angles are l-dependent
Ci          nl*nl if Euler angles are lm-dependent
Ci   theta :(ss only): spin spiral rotation angle
Ci   ...   The next 5 arguments apply to rotations of structure
Ci         matrices only (10000s digit mode < 4)
Ci   d1    :complex diagonal to be (optionally) added to sr
Ci   d2    :double precision diagonal to b (optionally) multiplied on sr
Ci   ...   The next 3 arguments apply when 10000s digit mode < 4:
Ci   ldd   :leading dimension of d1
Ci   ldima :cutoff of lower block of hamiltonian (leading, or
Ci         :augmentation dimension); see also ldimb.  Orbital pairs
Ci         :(i,j) with indxsh(i)>ldima or indxsh(j)>ldimb are discarded
Ci   ldimb :cutof of lower block of hamiltonian (second,
Ci         :or basis dimension); see ldima.
Ci   ...   The next 3 arguments apply to rotations of vectors only
Ci         (10000s digit mode = 4)
Ci   ldd   :number of vectors to rotate
Ci   ldima :Bounds which elements are to be included in rotation.
Ci         :Elements i with ldimb<indxsh(i)<ldima are ignored.
Ci   ldimb :see ldima.
Ci
Ci   lds   :leading dimension of s,sr,d1,d2
Ci   lds2  :second dimension of s,sr.
Ci          Set lds2=1 if s and sr are diagonal in R and L
Ci   s     :unrotated structure matrix, of dimension (lds,lds2)
Ci          See 10s digit of mode for how the complex part is placed.
Ci          In the case of a spin spiral, s consists of two such
Ci          matrices, thus s = s(lds,lds2,2).  The index to the
Ci          last dimension is 1 for s(k-q/2) and 2 for s(k+q/2).
Co Outputs
Co   sr    :rotated matrix, of dimension (lds,2,lds2,2) i.e.
Co         :sr <-  U+ s  U       if 10000s digit mode = 0 or 2
Co         :   <-  U+ s  U  + sr if 10000s digit mode = 1
Co         :   <-  U+ sr U       if 10000s digit mode = 3
Co         :   <-  U  sr         if 10000s digit mode = 4
Co         :Optionally, sr gets further scaled:
Co         : sr -> d2(RL) sr d2(R'L')  if 1000s digit mode >= 2
Co         : sr ->  sr + d1(RL)        if 1000s digit mode has 1s bit
Cl Local variables
Cl   w(ou) :work array holding spinor rotation matrices
Cr Remarks
Cu Updates
Cu   04 Jun 08 Generate Sil block for SS (10000s digit mode = 5)
Cu   08 Jul 03 Diagonal matrix added to S can be of fully relativistic
Cu             form (1000s digit = 4)
Cu   16 Jan 01 Add ability to rotate a set of vectors.
Cu   27 Apr 00 Separate augmentation and basis cutoffs
Cu             d1 is now a complex matrix
Cu            *argument list changed
C ----------------------------------------------------------------------
C     implicit none
      integer mode,ldd,lds,lds2,ldima,ldimb,nbas,neul,nl,indxsh(*)
      double precision theta,d1(2,ldd,2),d2(lds,2,2),eula(nbas,neul,3),
     .  s(lds,lds2,2),sr(lds,2,lds2,2,2)
C Local variables
      integer mode2,mode4,lmode,modei,ou
C heap:
      integer w(1)
      common /w/ w

      call tcn('rotspn')
      mode2 = mod(mode/100,10)
      mode4 = mod(mode/10000,10)
      if (mode4.eq.0 .or. mode4.eq.2 .or. mode4.eq.5)
     .  call dpzero(sr,lds*2*lds2*2*2)

C ... Set up spinor rotation matrices for input Euler angles
      call defdc(ou,4*nl*nl*nbas)
      call rotspu(mode2,nbas,nl,eula,neul,w(ou))

C     modei : strip mode of mode2 and mode4
      modei = mode - 10000*mode4 - 100*mode2

C ... Rotate spin spiral s to sr
      if (mode4 .eq. 2 .or. mode4 .eq. 5) then
        lmode = modei
        if (mode4 .eq. 5) lmode = modei + 100
        call rotss(lmode,nbas,nl,indxsh,ldima,ldimb,lds,lds2,theta,
     .    s,s,s,sr,sr,sr)
C       Case il block: we are finished
        if (mode4 .eq. 5) goto 99
C       Flag s as being already rotated; subsequent operations apply to sr
        lmode = modei + 10000
      elseif (mode4 .eq. 3) then
C       Flag s as being already rotated; subsequent operations apply to sr
        lmode = modei + 10000
      else
        lmode = modei
      endif

C ... Replace sr with U sr
      if (mode4 .eq. 4) then
        call rotsp2(lmode,nbas,nl,indxsh,w(ou),d2,ldd,ldimb,ldima,
     .    lds,lds2,sr,sr,sr)
C ... Replace sr with U s U+ or U sr U+
      else
        call rotsp1(lmode,nbas,nl,indxsh,w(ou),d1,d2,ldd,ldima,ldimb,
     .    lds,lds2,s,s,s,sr,sr,sr)
      endif

   99 continue
      call rlse(ou)
      call tcx('rotspn')
      end
      subroutine rotss(mode,nbas,nl,indxsh,ldima,ldimb,lds,lds2,theta,
     .  s0,s1,s2,sr0,sr1,sr2)
C- Rotates a structure matrix by a spin spiral
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit pertains to storage of hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (not implemented)
Ci          10s digit distinguishes how complex arithmetic is handled.
Ci              Only one pair of the set (s0,sr0) (s1,sr1), (s2,sr2)
Ci              is used, depending on this digit.
Ci           0: s,sr have real, imaginary separated, eg
Ci              s = s0(lds,lds2,2), with s0(*,*,1..2) = real..imag
Ci           1: s,sr are in complex*16 format, i.e. for real s
Ci              s = s1(2,lds,lds2), with s1(1..2,*,*) = real..imag
Ci           2: s,sr  have real, imaginary separated by columns, eg
Ci              s = s2(lds,2,lds2), with s2(*,1..2,*) = real..imag
Ci        100s digit:
Ci           0: S is ll block
Ci           1: S is il block; subtract ldimb from ldima in permutation
Ci       1000s digit not used
Ci      10000s digit not used
Ci     100000s digit (options are mutually exclusive)
Ci           0 s is not diagonal in R or L
Ci           1 s is diagonal in R but not L (not implemented)
Ci           2 s is diagonal in R and L; then lds2 should be 1
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   ldima :dimension of lower block of hamiltonian (leading, or
Ci         :augmentation dimension); see also ldimb.  Orbital pairs
Ci         :(i,j) with indxsh(i)>ldima or indxsh(j)>ldimb are discarded
Ci   ldimb :dimension of lower block of hamiltonian (second,
Ci         :or basis dimension); see ldima.
Ci   lds   :leading dimensions of s,sr
Ci   lds2  :second dimensions of s,sr
Ci   s*    :unrotated structure matrix; one of the three representations
Ci          is used, depending on 10s digit of mode.  It is dimensioned
Ci          s(lds,lds2,2).  The index to the last dimension
Ci          is 1 for s(k-q/2) and 2 for s(k+q/2).
Co Outputs
Co   sr*   :rotated structure matrix.; one of the three representations
Co         :is used, depending on 10s digit of mode
Cl Local variables
Cr Remarks
Cr   The SS rotation is
Cr
Cr   S(k,q,theta) = S(k-q/2) u1  +  S(k+q/2) u2
Cr
Cr     with
Cr
Cr        ( cos^2(theta/2) -sin(theta)/2 )      ( sin^2(theta/2)  sin(theta)/2 )
Cr   u1 = (                              ) u2 = (                              )
Cr        (-sin(theta)/2  sin^2(theta)/2 )      ( sin(theta)/2  cos^2(theta)/2 )
Cr
Cu Updates
Cu   04 Jun 08 Enable this routine to generate Sil block for SS (100s digit mode)
C ----------------------------------------------------------------------
C     implicit none
      integer mode,ldima,ldimb,lds,lds2,nbas,nl,indxsh(*)
      double precision theta
      double precision s0(lds,lds2,2,2),sr0(lds,2,lds2,2,2)
      double precision s2(lds,2,lds2,2),sr2(2,lds,2,lds2,2)
      double complex   s1(lds,lds2,2),sr1(lds,2,lds2,2)
C Local variables
      logical ldiag
      integer iprint
      integer lmi,ib,i,il,im,lmj,jb,j,jl,jm,i1,j1,ilm,jlm,kcplx,ldimp,
     .  mode2
      double precision s,cc,ss,u1(2,2),u2(2,2),s1r,s1i,s2r,s2i

      mode2 = mod(mode/100,10)
      call isanrg(mode2,0,1,'rotss','100s digit mode',.true.)
      kcplx = mod(mode/10,10)
      i = mod(mode/100000,10)
      if (i .eq. 1) call rxi('rotspn: mode not implemented',mode)
      ldiag = i .eq. 2
      ldimp = 0
      if (mode2 .eq. 1) ldimp = ldimb

C      if (mode2) then
C        call yprm('rotspn: S1',kcplx+2,s0,lds*lds2,lds,lds,lds2)
C        call yprm('rotspn: S2',kcplx+2,s0(1,1,1,2),lds*lds2,lds,lds,
C     .    lds2)
C      endif


C --- Rotation matrix for SS ---
      cc = dcos(theta/2)**2
      ss = dsin(theta/2)**2
      s  = dsin(theta)/2
      u1(1,1) = cc
      u1(1,2) = -s
      u1(2,1) = -s
      u1(2,2) = ss
      u2(1,1) = ss
      u2(1,2) =  s
      u2(2,1) =  s
      u2(2,2) = cc

      lmi = 0
      do  20  ib = 1, nbas
        ilm = 0
        do  25  il = 0, nl-1
        do  26  im = -il, il
        lmi = lmi+1
        ilm = ilm+1
C       Modify i to handle sil block
        i = indxsh(lmi)
        if (i .le. ldimp .or. i .gt. ldima) goto 26
        i = i - mode2*ldimb
        lmj = 0
        do  30  jb = 1, nbas
          jlm = 0
          do  35  jl = 0, nl-1
          do  36  jm = -jl, jl
            lmj = lmj+1
            jlm = jlm+1
            if (ldiag) then
              if (lmj .ne. 1) goto 30
              j = 1
            else
              j = indxsh(lmj)
              if (j .le. 0 .or. j .gt. ldimb) goto 36
            endif

C      ...  Rotate s into sr
            do  40  i1 = 1, 2
            do  40  j1 = 1, 2

              if (kcplx .eq. 0) then
                s1r = s0(i,j,1,1)
                s1i = s0(i,j,2,1)
                s2r = s0(i,j,1,2)
                s2i = s0(i,j,2,2)
              elseif (kcplx .eq. 1) then
                call rx('rotss not ready for this kcplx')
              else
                call rx('rotss not ready for this kcplx')
              endif

              sr0(i,i1,j,j1,1) = s1r*u1(i1,j1) + s2r*u2(i1,j1)
              sr0(i,i1,j,j1,2) = s1i*u1(i1,j1) + s2i*u2(i1,j1)

   40       continue
   36     continue
   35     continue
   30   continue
   26 continue
   25 continue
   20 continue

C     if (iprint() .ge. 110 .or. mode2) then
      if (iprint() .ge. 110) then
        call yprm('rotspn: s(ss)',kcplx+2,sr0,lds*lds2*4,lds*2,
     .    (ldima-ldimp)*2,ldimb*2)
      endif

      end
      subroutine rotsp1(mode,nbas,nl,indxsh,u,d1,d2,ldd,ldima,ldimb,lds,
     .  lds2,s0,s1,s2,sr0,sr1,sr2)
C- Rotates a structure matrix by spinor matrices.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit pertains to storage of hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (not implemented)
Ci          10s digit distinguishes how complex arithmetic is handled.
Ci              Only one pair of the set (s0,sr0) (s1,sr1), (s2,sr2)
Ci              is used, depending on this digit.
Ci           0: s,sr have real, imaginary separated, eg
Ci              s = s0(lds,lds2,2), with s0(*,*,1..2) = real..imag
Ci           1: s,sr are in complex*16 format, i.e. for real s
Ci              s = s1(2,lds,lds2), with s1(1..2,*,*) = real..imag
Ci           2: s,sr  have real, imaginary separated by columns, eg
Ci              s = s2(lds,2,lds2), with s2(*,1..2,*) = real..imag
Ci        100s digit
Ci           0: S is ll block
Ci           1: S is il block; subtract ldimb from ldima in permutation
Ci       1000s digit
Ci           1 Add d1 to diagonal
Ci           2 Scale sr_RLs,R'L's' by d2(RLs) d2(R'L's')
Ci           3 combination of 1+2
Ci      10000s digit
Ci           1 Input s is already spinor sr; here, s is not used.
Ci     100000s digit (options are mutually exclusive)
Ci           0 s is not diagonal
Ci           1 s is diagonal in R but not L (not implemented)
Ci           2 s is diagonal in R and L; then lds2 should be 1
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   u     :spinor rotation matrices
Ci   d1    :diagonal, in permuted order, to be added to diagonal of sr.
Ci          d1 is only used for specified values of mode.
Ci   d2    :diagonal, in permuted order, to scale sr.
Ci          d2 is only used for specified values of mode.
Ci   ldd   :leading dimension of d1
Ci   ldima :dimension of lower block of hamiltonian (leading, or
Ci         :augmentation dimension); see also ldimb.  Orbital pairs
Ci         :(i,j) with indxsh(i)>ldima or indxsh(j)>ldimb are discarded
Ci   ldimb :dimension of lower block of hamiltonian (second,
Ci         :or basis dimension); see ldima.
Ci   lds   :leading dimensions of s,sr,d1,d2
Ci   lds2  :second dimensions of s,sr
Ci   s*    :unrotated structure matrix; only one of the three forms of
Ci         :complex storage is used, depending on 10s digit of mode.
Co Outputs
Co   sr*   :rotated structure matrix, rotated by U+ s U.
Co          Only one of the three forms of complex storage is used,
Co          depending on 10s digit of mode
Cl Local variables
Cr Remarks
Cu Updates
Cu   04 Jun 08 Enable this routine to rotate Sil block (100s digit mode)
Cu   08 Jul 03 Diagonal matrix added to S can be of fully relativistic
Cu             form (1000s digit = 4)
C ----------------------------------------------------------------------
C     implicit none
      integer mode,ldima,ldimb,lds,lds2,ldd,nbas,nl,indxsh(*)
      double precision s0(lds,lds2,2),sr0(lds,2,lds2,2,2),d1(2,ldd,2,2)
      double precision s2(lds,2,lds2),sr2(lds,2,2,lds2,2),d2(lds,2)
      double complex   s1(lds,lds2),sr1(lds,2,lds2,2),u(2,2,nl*nl,nbas)
C Local variables
      double complex s11,s12,s21,s22,s,su(2,2),usu(2,2)
      logical l2,ls,ldiag
      integer lmi,ib,i,ii,il,im,lmj,jb,j,jl,jm,i1,j1,ilm,jlm,kcplx,
     .  mode0,mode2,mode3,mode4,iprint,k,ldimp

      mode0 = mod(mode,10)
      mode2 = mod(mode/100,10)
      kcplx = mod(mode/10,10)
      mode3 = mod(mode/1000,10)
      mode4 = mod(mode/10000,10)
      i = mod(mode/100000,10)
      l2 = mod(mode3,4) .ge. 2
      ls = mode4 .eq. 0
      ldiag = i .eq. 2
      call isanrg(mode0,0,0,'rotsp1','1s digit mode',.true.)
C     ? looks as though kcplx=2 is ok.  Need check
      call isanrg(kcplx,0,1,'rotsp1','kcplx',.true.)
      call isanrg(mode2,0,1,'rotsp1','100s digit mode',.true.)
      if (i .eq. 1) call rxi('rotspn: mode not implemented',mode)
      ldimp = 0
      if (mode2 .eq. 1) ldimp = ldimb
      if (mode2 .eq. 1 .and. l2)
     .  call rxi('rotsp1: this mode not implemented:',mode)

C --- sr <- U+ s U or d2 U+ s U d2 or U+ sr U or d2 U+ sr U d2 ---
      lmi = 0
      do  20  ib = 1, nbas
        ilm = 0
        do  25  il = 0, nl-1
        do  26  im = -il, il
        lmi = lmi+1
        ilm = ilm+1
        ii = indxsh(lmi)
        if (ii .le. ldimp .or. ii .gt. ldima) goto 26
        i = ii - mode2*ldimb
        lmj = 0
#ifndef DEBUG
        do  30  jb = 1, nbas
          jlm = 0
          do  35  jl = 0, nl-1
          do  36  jm = -jl, jl
            lmj = lmj+1
            jlm = jlm+1

            j = indxsh(lmj)
            if (j .le. 0 .or. j .gt. ldimb) goto 36

C           k=j except in diagonal case; then k=1
            if (ldiag) then
              if (ii .ne. j) goto 36
              k = 1
            else
              k = j
            endif

C      ...  Rotate s into sr = U+ s U
            if (ls) then
              if (kcplx .eq. 0) then
                s = dcmplx(s0(i,k,1),s0(i,k,2))
              elseif (kcplx .eq. 1) then
                s = s1(i,k)
              else
                s = dcmplx(s2(i,1,k),s2(i,2,k))
              endif
              do  112  i1 = 1, 2
              do  112  j1 = 1, 2
                usu(i1,j1) = (dconjg(u(1,i1,ilm,ib))*u(1,j1,jlm,jb) +
     .            dconjg(u(2,i1,ilm,ib))*u(2,j1,jlm,jb))*s
                if (l2) usu(i1,j1) = usu(i1,j1)*(d2(ii,i1)*d2(j,j1))
                if (kcplx .eq. 0) then
                  sr0(i,i1,k,j1,1) = sr0(i,i1,k,j1,1) + dble(usu(i1,j1))
                  sr0(i,i1,k,j1,2) = sr0(i,i1,k,j1,2) +dimag(usu(i1,j1))
                elseif (kcplx .eq. 1) then
                  sr1(i,i1,k,j1) = sr1(i,i1,k,j1) + usu(i1,j1)
                else
                  sr2(i,i1,1,k,j1) = sr2(i,i1,1,k,j1) + dble(usu(i1,j1))
                  sr2(i,i1,2,k,j1) = sr2(i,i1,2,k,j1) +dimag(usu(i1,j1))
                endif
  112         continue
C      ...  Overwrite sr with its rotation
            else
              if (kcplx .eq. 0) then
                s11 = dcmplx(sr0(i,1,k,1,1),sr0(i,1,k,1,2))
                s12 = dcmplx(sr0(i,1,k,2,1),sr0(i,1,k,2,2))
                s21 = dcmplx(sr0(i,2,k,1,1),sr0(i,2,k,1,2))
                s22 = dcmplx(sr0(i,2,k,2,1),sr0(i,2,k,2,2))
              elseif (kcplx .eq. 1) then
                s11 = sr1(i,1,k,1)
                s12 = sr1(i,1,k,2)
                s21 = sr1(i,2,k,1)
                s22 = sr1(i,2,k,2)
              else
                s11 = dcmplx(sr2(i,1,1,k,1),sr2(i,1,2,k,1))
                s12 = dcmplx(sr2(i,1,1,k,2),sr2(i,1,2,k,2))
                s21 = dcmplx(sr2(i,2,1,k,1),sr2(i,2,2,k,1))
                s22 = dcmplx(sr2(i,2,1,k,2),sr2(i,2,2,k,2))
              endif
              su(1,1) = s11*u(1,1,jlm,jb) + s12*u(2,1,jlm,jb)
              su(2,1) = s21*u(1,1,jlm,jb) + s22*u(2,1,jlm,jb)
              su(1,2) = s11*u(1,2,jlm,jb) + s12*u(2,2,jlm,jb)
              su(2,2) = s21*u(1,2,jlm,jb) + s22*u(2,2,jlm,jb)
              do  212  i1 = 1, 2
              do  212  j1 = 1, 2
                usu(i1,j1) = dconjg(u(1,i1,ilm,ib))*su(1,j1) +
     .                       dconjg(u(2,i1,ilm,ib))*su(2,j1)
                if (l2) usu(i1,j1) = usu(i1,j1)*(d2(ii,i1)*d2(j,j1))
                if (kcplx .eq. 0) then
                  sr0(i,i1,k,j1,1) = dble(usu(i1,j1))
                  sr0(i,i1,k,j1,2) = dimag(usu(i1,j1))
                elseif (kcplx .eq. 1) then
                  sr1(i,i1,k,j1) = usu(i1,j1)
                else
                  sr2(i,i1,1,k,j1) = dble(usu(i1,j1))
                  sr2(i,i1,2,k,j1) = dimag(usu(i1,j1))
                endif
  212         continue
            endif

   36     continue
   35     continue
   30   continue
C ... debugging ... overwrite sr with u
#else
      do  i1 = 1, 2
      do  j1 = 1, 2
        sr0(i,i1,i,j1,1) = dble(u(i1,j1,ilm,ib))
        sr0(i,i1,i,j1,2) = dimag(u(i1,j1,ilm,ib))
      enddo
      enddo
#endif
   26 continue
   25 continue
   20 continue

#if DEBUG
C      call yprm('rotspn: U',kcplx+2,sr0,lds*lds2*4,lds*2,
C     .  (ldima-ldimp)*2,ldimb*2)
#endif


C --- Add diagonal function d1 ---
      if (mod(mode3,2) .eq. 1) then
        do  50  i1 = 1, 2
          if (kcplx .eq. 0) then
            do  40  i = 1, min(ldima,ldimb)
            sr0(i,i1,i,i1,1) = sr0(i,i1,i,i1,1) + d1(1,i,i1,1)
   40       sr0(i,i1,i,i1,2) = sr0(i,i1,i,i1,2) + d1(2,i,i1,1)
          elseif (kcplx .eq. 1) then
            do  41  i = 1, min(ldima,ldimb)
   41       sr1(i,i1,i,i1) = sr1(i,i1,i,i1) +
     .                       dcmplx(d1(1,i,i1,1),d1(2,i,i1,1))
          else
            do  42  i = 1, min(ldima,ldimb)
            sr2(i,i1,1,i,i1) = sr2(i,i1,1,i,i1) + d1(1,i,i1,1)
   42       sr2(i,i1,2,i,i1) = sr2(i,i1,2,i,i1) + d1(2,i,i1,1)
          endif
   50 continue

      else if (mod(mode3/4,2) .eq. 1) then
        do  150  i1 = 1, 2
          if (kcplx .eq. 0) then
            do  140  i = 1, min(ldima,ldimb)
            sr0(i,i1,i,i1,1) = sr0(i,i1,i,i1,1) + d1(1,i,i1,i1)
  140       sr0(i,i1,i,i1,2) = sr0(i,i1,i,i1,2) + d1(2,i,i1,i1)
          elseif (kcplx .eq. 1) then
            do  141  i = 1, min(ldima,ldimb)
  141       sr1(i,i1,i,i1) = sr1(i,i1,i,i1) +
     .                       dcmplx(d1(1,i,i1,i1),d1(2,i,i1,i1))
          else
            do  142  i = 1, min(ldima,ldimb)
            sr2(i,i1,1,i,i1) = sr2(i,i1,1,i,i1) + d1(1,i,i1,i1)
  142       sr2(i,i1,2,i,i1) = sr2(i,i1,2,i,i1) + d1(2,i,i1,i1)
          endif
  150   continue
       do  240  i = 1+1, min(ldima,ldimb)
         if (kcplx .eq. 0) then
           sr0(i-1,1,i,2,1) = sr0(i-1,1,i,2,1) + d1(1,i-1,1,2)
           sr0(i-1,1,i,2,2) = sr0(i-1,1,i,2,2) + d1(2,i-1,1,2)
           sr0(i,2,i-1,1,1) = sr0(i,2,i-1,1,1) + d1(1,i,2,1)
           sr0(i,2,i-1,1,2) = sr0(i,2,i-1,1,2) + d1(2,i,2,1)
         elseif (kcplx .eq. 1) then         
           sr1(i-1,1,i,2) = sr1(i-1,1,i,2) +
     .       dcmplx(d1(1,i-1,1,2),d1(2,i-1,1,2))
           sr1(i,2,i-1,1) = sr1(i,2,i-1,1) +
     .       dcmplx(d1(1,i,2,1),d1(2,i,2,1))
         else
           sr2(i-1,1,1,i,2) = sr2(i-1,1,1,i,2) + d1(1,i-1,1,2)
           sr2(i-1,1,2,i,2) = sr2(i-1,1,2,i,2) + d1(2,i-1,1,2)
           sr2(i,2,1,i-1,1) = sr2(i,2,1,i-1,1) + d1(1,i,2,1)
           sr2(i,2,2,i-1,1) = sr2(i,2,2,i-1,1) + d1(2,i,2,1)
         endif
  240 continue
      endif

C     if (iprint() .ge. 110 .or. mode2) then
      if (iprint() .ge. 110) then
        call yprm('rotspn: sr',kcplx+2,sr0,lds*lds2*4,lds*2,
     .    (ldima-ldimp)*2,ldimb*2)
      endif

      end
      subroutine rotsp2(mode,nbas,nl,indxsh,u,d2,nvec,ldimp,ldima,lds,
     .  lds2,sr0,sr1,sr2)
C- Rotates a list of vectors by spinor matrices.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit pertains to storage of hamiltonian (not used)
Ci          10s digit distinguishes how complex arithmetic is handled.
Ci              Only one of the set sr = (sr0,sr1,sr2)
Ci              is used, depending on this digit.
Ci           0: sr has real, imaginary separated, i.e.
Ci              sr = s0(lds,lds2,2), with s0(*,*,1..2) = real..imag
Ci           1: sr are in complex*16 format, i.e.
Ci              sr = s1(2,lds,lds2), with s1(1..2,*,*) = real..imag
Ci           2: sr  has real, imaginary separated by columns, i.e.
Ci              sr = s2(lds,2,lds2), with s2(*,1..2,*) = real..imag
Ci        100s digit not used
Ci       1000s digit
Ci           1 Scale sr_RLs,j by d1(RLs) after rotation (unimplemented)
Ci           2 Scale sr_RLs,j by d2(RLs) before rotation
Ci      10000s digit not used
Ci     100000s digit not used
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutation indices ordering orbitals in downfolding order
Ci   u     :spinor rotation matrices
Ci   d2    :diagonal, in permuted order, to scale sr.
Ci         :d2 is only used for specified values of mode.
Ci   nvec  :number of vectors to rotate
Ci   ldimp :hamiltonian block consists of orbitals betw. ldimp and ldima
Ci   ldima :offset to last orbital in block to be rotated:
Ci         :Only elements i with ldimp<indxsh(i)<=ldima are rotated.
Ci   lds   :leading dimension of sr and d2
Ci   lds2  :second dimension of sr
Cio Inputs/Outputs
Cio   sr*  :On input, unrotated vectors sr(*,1..nvec)
Cio        :On output, overwritten by rotated vector U sr(*,1..nvec)
Cio        :Only one of the three forms of complex storage is used,
Cio        :(sr0,sr1,sr2) depending on 10s digit of mode.
Cl Local variables
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer mode,ldima,ldimp,lds,lds2,nvec,nbas,nl,indxsh(*)
      double precision sr0(lds,2,lds2,2),sr2(lds,2,2,lds2),d2(lds,2)
      double complex   sr1(lds,2,lds2),u(2,2,nl*nl,nbas)
C Local variables
      double complex s1,s2,us(2)
      integer lmi,ib,i,il,im,j,ilm,mode0,kcplx,mode3,iprint

      mode0 = mod(mode,10)
      kcplx = mod(mode/10,10)
      mode3 = mod(mode/1000,10)
c     mode4 = mod(mode/10000,10)
      call isanrg(mode0,0,0,'rotsp2','1s digit mode',.true.)
C     ? looks as though kcplx=2 is ok.  Need check
      call isanrg(kcplx,0,1,'rotsp2','kcplx',.true.)

C --- sr <- U sr ---
      lmi = 0
      do  20  ib = 1, nbas
        ilm = 0
        do  25  il = 0, nl-1
        do  26  im = -il, il
        lmi = lmi+1
        ilm = ilm+1
        i = indxsh(lmi)
        if (i .le. ldimp .or. i .gt. ldima) goto 26
        do  30  j = 1, nvec

C          if (ib .eq. 2 .and. j .eq. 7) then
C            print *, 'ha'
C            print *, dcmplx(sr0(i,1,j,1),sr0(i,1,j,2))
C            print *, dcmplx(sr0(i,2,j,1),sr0(i,2,j,2))
C            print *, u(1,1,ilm,ib)
C            print *, u(1,2,ilm,ib)
C            print *, u(2,1,ilm,ib)
C            print *, u(2,2,ilm,ib)
C          endif

C    ...  Overwrite sr with U sr
          if (kcplx .eq. 0) then
            s1 = dcmplx(sr0(i,1,j,1),sr0(i,1,j,2))
            s2 = dcmplx(sr0(i,2,j,1),sr0(i,2,j,2))
            if (mode3 .eq. 2) then
              s1 = d2(i,1) * s1
              s2 = d2(i,2) * s2
            endif
            us(1) = u(1,1,ilm,ib)*s1 + u(1,2,ilm,ib)*s2
            us(2) = u(2,1,ilm,ib)*s1 + u(2,2,ilm,ib)*s2
            sr0(i,1,j,1) = dble(us(1))
            sr0(i,1,j,2) = dimag(us(1))
            sr0(i,2,j,1) = dble(us(2))
            sr0(i,2,j,2) = dimag(us(2))
          elseif (kcplx .eq. 1) then
            s1 = sr1(i,1,j)
            s2 = sr1(i,2,j)
            us(1) = u(1,1,ilm,ib)*s1 + u(1,2,ilm,ib)*s2
            us(2) = u(2,1,ilm,ib)*s1 + u(2,2,ilm,ib)*s2
            sr1(i,1,j) = us(1)
            sr1(i,2,j) = us(2)
          else
            s1 = dcmplx(sr2(i,1,1,j),sr2(i,1,2,j))
            s2 = dcmplx(sr2(i,1,1,j),sr2(i,1,2,j))
            us(1) = u(1,1,ilm,ib)*s1 + u(1,2,ilm,ib)*s2
            us(2) = u(2,1,ilm,ib)*s1 + u(2,2,ilm,ib)*s2
            sr2(i,1,1,j) = dble(us(1))
            sr2(i,1,2,j) = dimag(us(1))
            sr2(i,2,1,j) = dble(us(2))
            sr2(i,2,2,j) = dimag(us(2))
          endif
   30   continue
   26 continue
   25 continue
   20 continue

      if (iprint() .ge. 110) then
       call yprm('rotspn: Us',kcplx+2,sr0,lds*lds2*2,lds*2,ldima*2,nvec)
      endif

      end

