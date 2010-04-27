Cgetarg...info...           structure ['sblham', 'ham', 'sham']
Cgetarg...info...           structure ['sblham', 'pot', 'spot']
Cgetarg...info...           use_to_add ['sblham', 'm_struc_def']
      subroutine sblham(sctrl,sham,spot,fnam,ib1,ib2,nbasx,lbloch,plat,
     .  ldim,idim,ldl,ldi,ldl2,klu,zp,qp,sll,sil,sii)
       
       use m_struc_def  !Cgetarg

C- Generate subblock of Bloch summed hamiltonian from file of strux
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct containing parameters governing program flow
Ci   sham  :struct containing information about the hamiltonian
Ci   spot  :struct containing information about the potential
Ci   fnam  :file name
Ci   ib1   :the neighbor list of pairs for which to accumulate
Ci          sum starts at beginning of neighbor list for ib1
Ci   ib2   :the neighbor list of pairs for which to accumulate
Ci          sum ends at end of neighbor list for ib2
Ci   nbasx :At least the largest site which will be touched in this
Ci          subblock.  It should not exceed the number of sites
Ci          for which indexing information is available
Ci   lbloch:passed through to bloch.
Ci          1s digit concerns storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci         10s digit distinguishes how complex arithmetic is handled
Ci           1: s is returned complex*16 format:
Ci              s = s(2,ldl,ldl2), with s(1,*) = real, s(2,*) = imag
Ci           2: s has real, imaginary separated by columns
Ci              s = s(ldl,2,ldl2), with s(*,1..2,*) = real..imag
Ci           3: s has real, imaginary separated
Ci              s = s(ldl,ldl2,2), with s(*,*,1..2) = real..imag
Ci        100s digit
Ci           1 if to add to s (s not initialized to zero)
Ci           2 subtract from s
Ci           3 combination of 1+2
Ci       1000s digit 1 if to convert s to spherical harmonics
Ci      10000s digit 1 ?? Suppress accumulation into Bloch sum pairs
Ci                       whose site index is out of (ib1,ib2) range
Ci   ldim  :dimension of basis in the lower set.  Orbitals for which
Ci         :0<iprm(i)<=ldim are accumulated into sll.
Ci   idim  :dimension of downfolded sii
Ci   ldl   :leading dimension of sll, for dimensioning purposes
Ci          Must have
Ci   ldi   :leading and second dimension of sii
Ci   ldl2  :second dimension of sll and sil
Ci   klu   :size of sub- and super-diagonal, if sll is in banded format
Ci   zp    :complex energy (not used now)
Ci   qp    :k-point
Co Outputs
Co   sll   :lower-lower bloch of Bloch summed matrix
Co   sil   :lower-intermediate bloch of Bloch summed matrix
Co   sii   :intermediate-intermediate bloch of Bloch summed matrix
Cr Remarks
Cr   This routine generates a subblock of the Bloch summed hamiltonian,
Cr   by including only orbitals associated with sites ib1..ib2.
Cr
Cr   The band form follows LAPACK band storage conventions:
Cr   s(i,j) is stored in location (kl+ku+1+i-j,j)
Cr   with kl,ku = size of sub- and super-diagonal.
Cr   Here we take kl=ku=klu.
Cu Updates
C    23 Aug 05 (A. Chantis) fixed for spherical harmonics.
Cu   10 Jan 02 slbham reads pf from pot->opalp instead of pot->opf
C ----------------------------------------------------------------------
C     implicit none
      character*(*) fnam
Cgetarg       double precision sctrl(1),sham(1),spot(1)
       type(s_ctrl)::sctrl
       type(s_ham)::sham
       type(s_pot)::spot

      double precision plat(3,3),qp(3),zp(2)
      integer ib1,ib2,nbasx,klu,lbloch,ldi,ldl,ldl2
C     real + imaginary storage mode
      double precision sll(ldl,ldl2,2),sil(ldi,ldl2,2),sii(ldi,ldi,2)
C     complex*16 storage mode
C     double precision sll(2,ldl,ldl2),sil(2,ldi,ldl2),sii(2,ldi,ldi)
C     real + imaginary in columns storage mode
C     double precision sll(ldl,2,ldl2),sil(ldi,2,ldl2),sii(ldi,2,ldi)
C     For potential
      integer opf
C     For iostr
      logical iostr,ltmp, bittst
      integer itral,lmaxw,nitab,nkap,oalph,oiax,ontab
      double precision ekap2(20),ckbas
C     For bloch
      integer idim,is1,is2,ldim,nbasp,nl,oidxsh,oiprm,oiwk,os,lham
C     Heap
      integer w(1)
      common /w/ w

      call tcn('sblham')

C --- Read strux from disk file ---
      ltmp = iostr(2+8,fnam,nl,nbasp,nkap,ekap2,itral,ckbas,lmaxw,nitab,
     .  oalph,oiax,ontab,os)

C --- Bloch sum strux ---
      is1 = w(ontab+ib1-1)+1
      is2 = w(ontab+ib2)
C ... Restrict pairs outside of block
Cgetarg...info...           integer :: oidxsh , integer(8) :: sham%oindxo 1
Cgetarg...info...           integer :: lham , integer(8) :: sham%lham 1
Cgetarg       call upack('ham oindxo lham',sham,oidxsh,lham,0,0,0)
       
       oidxsh=sham%oindxo
       lham=sham%lham

      call defi(oiwk,nbasx)
      call defi(oiprm,nbasx)
      call icopy(nbasx*nl*nl,w(oidxsh),1,w(oiprm),1)
C     First call excludes all orbitals
      call oraddp(1,w(oiax),1,0,1,nbasx,nl*nl,w(oiwk),w(oiprm))
C     Second call includes orbitals in range ib1..ib2
      call oraddp(0,w(oiax),is1,is2,ib1,ib2,nl*nl,w(oiwk),w(oiprm))
C     A third call could exclude those we didn't want in the first place
C     (not implemented)
      call rlse(oiwk)
      if (bittst(lham,256)) then
         lbloch = lbloch + 1000
      endif   
      call bloch(lbloch,qp,nl,plat,nl**2,w(oiprm),is1,is2,w(oiax),w(os),
     .  nl**2,1,1,ldim,ldim,idim,ldl,ldi,ldl2,klu,sll,sil,sii)

C --- ASA : Add -P to sll with (S-P) ---
Cgetarg...info...           integer :: opf , integer(8) :: spot%opalp 1
Cgetarg       call upack1('pot opalp',spot,opf)
       
       opf=spot%opalp


      call sblhm1(lbloch,nl*nl,ib1,ib2,w(oiprm),w(opf),0,ldim,idim,ldl,
     .  ldi,ldl2,klu,sll,sil,sii)

      call rlse(oalph)

      call tcx('sblham')
Cgetarg       end
       end subroutine sblham 


      subroutine sblhm1(lbloch,nl2,ib1,ib2,iprm,pfun,offp,ldim,idim,ldl,
     .  ldi,ldl2,klu,sll,sil,sii)

C- Add or subtract diagonal P from hamiltonian subblock
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:bits are compatible with those of routine bloch.
Ci          1s digit concerns storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci         10s digit distinguishes how complex arithmetic is handled
Ci           1: s is returned complex*16 format:
Ci              s = s(2,ldl,ldl2), with s(1,*) = real, s(2,*) = imag
Ci           2: s has real, imaginary separated by columns
Ci              s = s(ldl,2,ldl2), with s(*,1..2,*) = real..imag
Ci           3: s has real, imaginary separated
Ci              s = s(ldl,ldl2,2), with s(*,*,1..2) = real..imag
Ci        100s digit
Ci           0 subtract P from s
Ci           2 add P to s
Ci   nl2   :spacing between offsets to pfun in successive sites:
Ci          offset to pfun for site ib is nl2*(ib-1).
Ci   ib1   :the neighbor list of pairs for which to accumulate
Ci          sum starts at beginning of neighbor list for ib1
Ci   ib2   :the neighbor list of pairs for which to accumulate
Ci          sum ends at end of neighbor list for ib2
Ci   iprm  :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   pfun  :vector of potential functions, in iprm order
Ci   offp  :offset to potential function, ib=1
Ci   ldim  :dimension of lower set of orbitals. See iprm, above.
Ci   idim  :dimension of intermediate set. See iprm, above.
Ci   ldl   :leading dimension of sll
Ci   ldi   :leading and second dimension of sii
Ci   ldl2  :second dimension of sll and sil
Ci   klu   :size of sub- and super-diagonal, if s stored banded form
Co Outputs
Co   sll   :lower-lower block of hamiltonian matrix
Co   sil   :lower-intermediate block of hamiltonian matrix
Co   sii   :intermediate-intermediate block of hamiltonian matrix
Cr Remarks
Cr   Complex diagonal pfun is added to sll,sil,sii
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lbloch,nl2,ib1,ib2,ldim,idim,ldl,ldi,ldl2,klu,offp,iprm(*)
      double precision pfun(2,nl2)
C     sll,sil,sii assuming real + imaginary storage mode
      double precision sll(ldl,ldl2,2),sil(ldi,ldl2,2),sii(ldi,ldi,2)
C ... Local parameters
      integer nlmx
      parameter (nlmx=81)
      double precision pwk(nlmx,2,nlmx),fac
      integer ia,ib,lblchi,kcplx,ld11,ld12,ld13,ld21,ld22,ld23,lidim,
     .  offa,offb,offpi,lm,oi,ip

      if (nl2 .gt. nlmx) call rxi('sblhm1: need nlmx at least',nl2)
      lidim = ldim+idim

C     Pick up true dimensions of sll,sil,sii from formal ones
      kcplx = mod(lbloch/10,10)
      fac = -1
      if (mod(lbloch/100,10) .ge. 2) fac = 1
      call cplxdm(kcplx,ldl,ldl2,ld11,ld21,oi,oi)
      call cplxdm(kcplx,ldi,ldl2,ld12,ld22,oi,oi)
      call cplxdm(kcplx,ldi,ldi,ld13,ld23,oi,oi)
C     Tell pblch1 that P is diagonal and complex
      lblchi = mod(lbloch,1000) + 11000

C --- For each site in the subblock, do ---
      do  10  ib = ib1, ib2

C  ...  Copy pfun(iprm order) to the diagonal of matrix pwk (L order)
        ia = ib
        offpi = nl2*(ia-1)
        do  12  lm = 1, nl2
          ip = iprm(offpi+lm)
          if (ip .le. lidim) then
            pwk(lm,1,lm) = pfun(1,ip+offp)
            pwk(lm,2,lm) = pfun(2,ip+offp)
          endif
   12   continue

C   --- Lower-lower block ---
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        call pblch1(lblchi,nl2,offa,offb,ld11,ld21,klu,iprm,0,ldim,
     .    0,ldim,pwk,pwk,nlmx,fac,0d0,sll)

C   --- Intermediate-lower block ---
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        if (idim .eq. 0) goto 10
        call pblch1(lblchi,nl2,offa,offb,ld12,ld22,klu,iprm,ldim,lidim,
     .    0,ldim,pwk,pwk,nlmx,fac,0d0,sil)

C   --- Intermediate-intermediate block ---
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        if (idim .eq. 0) goto 10
        call pblch1(lblchi,nl2,offa,offb,ld13,ld23,klu,iprm,ldim,lidim,
     .    ldim,lidim,pwk,pwk,nlmx,fac,0d0,sii)

   10 continue
Cgetarg       end
       end subroutine sblhm1 


