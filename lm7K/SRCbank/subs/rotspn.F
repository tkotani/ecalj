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

      call rx('rotspn not implemented')

Cgetarg       end
       end subroutine rotspn 


