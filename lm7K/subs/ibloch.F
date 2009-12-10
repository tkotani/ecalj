      subroutine ibloch(lbloch,qp,wtqp,plat,mxorb,iprmb,is1,is2,iax,
     .  s,sc,sz,nds,isp,nsp,ldima,ldimb,ldl,ldl2,sll)
C- Accumulate unsymmetrized inverse Bloch transform from one q
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:1s digit pertains to r.s. transform:
Ci           0 if r.s. transform is real
Ci           1 if r.s. transform is complex
Ci
Ci          10s digit distinguishes how complex arithmetic is handled
Ci              (parameter kcplx below)
Ci           0: sll or s has real, imaginary separated
Ci              sll = sll(ldl,ldl2,2), with sll(*,*,1..2) = real..imag
Ci           1: sll or s are in complex*16 format:
Ci              sll = sll(2,ldl,ldl2), with sll(1..2,*,*) = real..imag
Ci           2: sll or s have real, imaginary separated by columns
Ci              sll = sll(ldl,2,ldl2), with sll(*,1..2,*) = real..imag
Ci
Ci        100s digit:
Ci           2 Accumulate inverse transform of -sll
Ci
Ci       1000s digit not used
Ci      10000s digit not used
Ci
Ci     100000s digit pertains to whether s has permuted orbital order
Ci           0 if s has normal order
Ci           1 if s has a permuted orbital order
Ci
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci
Ci   qp    :k-point
Ci   wtqp  :BZ weight for qp
Ci
Ci   plat  :primitive lattice vectors, in units of alat
Ci
Ci   mxorb :leading dimension of iprmb
Ci
Ci   iprmb :permutation indices ordering orbitals in downfolding order.
Ci          Unpermuted orbitals run from 1 ... n, with indices to
Ci          orbitals corresponding to site ib starting at 1+mxorb(ib-1).
Ci          Orbital indices i for which 0<iprmb(i)<=ldim are accumulated
Ci          into the lower set; those for which ldim<iprmb(i)<=ldim+idim
Ci          are accumulated the intermediate set.  Indices for which
Ci          iprmb(i) lie outside this range are not accumulated.
Ci
Ci   is1,is2:Bloch sum contribution from pairs is1..is2
Ci
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci          For each pair i, the following portion is used by bloch:
Ci          iax(1,i): basis atom for source (column) index
Ci                    If <= 0, bloch excludes this pair from the sum
Ci          iax(2,i): basis atom for augmentation (row) index
Ci                    If <= 0, bloch excludes this pair from the sum
Ci          iax(3..5,i): lattice vectors separating the two sites
Ci                        as multiples of plat
Ci          iax(8,i): points to an equivalent pair, if nonzero
Ci
Ci   sll   :Bloch matrix for this qp, from which inverse transform is
Ci         :accumulated.  sll and zsk hold the same information; only one
Ci         :is used depending on how complex arithmetic is specified
Ci         :(see lbloch).  zsk is in complex*16 format; sll has real
Ci         :and imaginary parts separated
Ci
Ci   nds   :leading dimension of s
Ci
Ci   ldima :cutoff for lower set of orbitals in the augmentation
Ci          (row) dimension.  See iprmb, above.
Ci
Ci   ldimb :cutoff for lower set of orbitals, in the source
Ci          (column) dimension.  See iprmb, above.
Ci
Ci   ldl   :leading dimension of sll
Ci
Ci   ldl2  :second dimension of sll and sil
Ci
Co Outputs
Co   s,zs  :The inverse Bloch transform is accumulated in s.  s and zs
Co         :hold the same information; only one is used depending on how
Co         :complex arithmetic is specified (see lbloch).  zs is in
Co         :complex*16 format; sll has real and imaginary parts separated
Cr Remarks:
Cr   ibloch produces an unsymmetrized transform, the inverse of the one
Cr   defined in subroutine bloch:
Cr     s(r1,l1,T+r2,l2) = sum_k w_k s(k;r1,l1,r2,l2) exp(-i k . T)
Cr   It is correct only if the k-sum is taken over the full Brillouin
Cr   zone.  If the irreducible BZ is used, the unsymmetrized transform
Cr   must be symmetrized; see rsmsym.f.
Cr
Cr   ibloch has a structure similar to routine bloch (forward bloch
Cr   transform).  Some differences with bloch:
Cr
Cr     * No banded format for sll
Cr     * No intermediate or higher blocks
Cr     * No rotation to spherical harmonics
Cr     * No onsite-only mode
Cr
Cr   Note: if the point group operations include those that come from
Cr   time-reversal symmetry, but not from the space group, weights
Cr   from those points which include in their star points which are
Cr   only related by time-reversal symmetry should be flagged with
Cr   a (-) sign.  See bzmesh.f
Cu Updates
Cu   01 Mar 03 Extensively revised
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lbloch,nds,is1,is2,ldima,ldimb,ldl,ldl2,niax,mxorb,isp,nsp
      parameter (niax=10)
      integer iprmb(*),iax(niax,is2)
      double precision qp(3),wtqp,plat(3,3)
      double precision sll(ldl,ldl2,2)
C     These three arrays can all point
C     Arrays s,sc,sz hold the same data but assume different complex
C     formats (see 10s digit of lbloch above).  Only one form is
C     used here; thus all can point to the same address space.
C     s is used when kcplx=0 OR when output s is REAL
      double precision s(nds,nds,nsp,is2)
C     sc is used when kcplx=2 AND when s is COMPLEX
      double precision sc(nds,2,nds,nsp,is2)
C     sz is used when kcplx=1 AND when s is COMPLEX
      double precision sz(2,nds,nds,nsp,is2)
C ... Local parameters
      integer ia,ib,isite,j,k,kcplx,ksite,ksite2,lblchi,lblchp,ld11,
     .  ld21,offa,offb,oi,scplx
      double precision twopi,TdotK,cosT,sinT,wc,ws
C     double precision wc0,ws0

C --- Setup ---
      call tcn('ibloch')
      twopi = 8*datan(1d0)
      lblchi = mod(lbloch,10000)
C     onsite = mod(lblchi/1000,10) .ge. 2
      j = mod(lbloch/10000,10)

C     Pick up true dimensions of sll,sil,sii from formal ones
      kcplx = mod(lblchi/10,10)
      scplx = mod(lblchi,10)
      call isanrg(scplx,0,1,'ibloch:','1s digit lbloch',.true.)

      call cplxdm(kcplx,ldl,ldl2,ld11,ld21,oi,oi)
C     call cplxdm(kcplx,ldi,ldl2,ld12,ld22,oi,oi)
C     call cplxdm(kcplx,ldi,ldi,ld13,ld23,oi,oi)
C     lblchp passed to pblch1. s cplx -> shift flag from 40 to 1000
      lblchp = lblchi
      if (scplx .eq. 1) then
        lblchp = lblchp + 1000 - 1
      endif
C     s has permuted orbital order
      j = mod(lbloch/100000,10)
      lblchp = lblchp + j*100000

C --- Add contribution to inverse Bloch sum to each R,R' pair ---
      do  100  isite = is1, is2

        ia = iax(2,isite)
        ib = iax(1,isite)
C       Any nonpositive site indices are excluded from sum
        if (ia .le. 0 .or. ib .le. 0) goto 100
C        if (onsite) then
C          if (ia .ne. ib .or. iax(3,isite) .ne. 0 .or.
C     .      iax(4,isite) .ne. 0 .or. iax(5,isite) .ne. 0) goto 100
C        endif

C   ... (wc,ws) = wtqp*exp(-i k . T)
        TdotK = 0
        do  30  j = 1, 3
        do  30  k = 1, 3
   30   TdotK = TdotK + twopi*qp(j)*plat(j,k)*iax(2+k,isite)
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
        if (mod(lblchi/100,10) .ge. 2) then
          cosT = -cosT
          sinT = -sinT
        endif
        wc  =  abs(wtqp)/2*cosT
        ws  = -abs(wtqp)/2*sinT
        if (wtqp .lt. 0) then
          wc = wc/2
          ws = ws/2
        endif

C   ... Use equivalent of isite to some other site, if it exists
        ksite = isite
        if (iax(8,isite) .ne. 0) ksite = iax(8,isite)

C   --- Rotate s to spherical harmonics ---
C        if (mod(mod(lblchi/1000,10),2) .ne. 0)
C     .    call s2sph(2,nl,nl,s(1,1,isp,ksite),nds,nds,ndss,ndss,sph)

        offb = mxorb*(ib-1)
        offa = mxorb*(ia-1)

C       scaling if s is complex arithmetic
        ksite2 = (scplx+1)*(ksite-1) + 1
C       if (isite .eq. 3) then
C         wc0 = wc
C         ws0 = ws
C       endif
        call iblch1(lblchp,mxorb,offa,offb,ld11,ld21,iprmb,0,ldima,0,
     .    ldimb,s(1,1,isp,ksite2),sc(1,1,1,isp,ksite),
     .    sz(1,1,1,isp,ksite),nds,wc,ws,sll)
        if (wtqp .lt. 0) then
          offb = mxorb*(ib-1)
          offa = mxorb*(ia-1)
          call iblch1(lblchp+100,mxorb,offa,offb,ld11,ld21,iprmb,0,
     .      ldima,0,ldimb,s(1,1,isp,ksite2),sc(1,1,1,isp,ksite),
     .      sz(1,1,1,isp,ksite),nds,wc,-ws,sll)
        endif

  100 continue

C      isite = 3
C      print 333, ' ibloch: isite,qp',isite,qp,' wc,ws=',wc0,ws0
C  333 format(a,i4,3f8.4,a,2f11.6)
C      call pvtro9(scplx,1,1,nds,s(1,1,isp,isite),sz(1,1,1,isp,isite))

      call tcx('ibloch')
      end

      subroutine iblch1(lbloch,mxorb,offa,offb,lds,ld2,iprmb,
     .  hdpa,hdna,hdpb,hdnb,s,sc,sz,nds,wc,ws,sk)
C- Contribution of one pair to Bloch sum of strux
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch: 1s digit concerns storage of Bloch summed sk
Ci           0: sk is stored in normal, unpacked form
Ci           1: sk is stored in banded form (see Remarks)
Ci              The band form follows LAPACK band storage conventions:
Ci              sk(i,j) is stored in location (kl+ku+1+i-j,j)
Ci              with kl,ku = size of sub- and super-diagonal.
Ci              Here we take kl=ku=klu.
Ci          10s digit distinguishes how complex arithmetic is handled
Ci           0: sk has real, imaginary separated
Ci              sk = sk(ldl,ldl2,2), with sk(*,*,1..2) = real..imag
Ci           1: sk is returned complex*16 format:
Ci              sk = sk(2,ldl,ldl2), with sk(1..2,*,*) = real..imag
Ci           2: sk has real, imaginary separated by columns
Ci              sk = sk(ldl,2,ldl2), with sk(*,1..2,*) = real..imag
Ci        100s digit
Ci           0  use sk
Ci           1  use sk+
Ci       1000s digit 1 if real-space s is complex
Ci             In this case, sc is used in place of s.
Ci      10000s digit 1 if s(or sc) is diagonal.
Ci     100000s digit 1 if s(or sc) has a permuted orbital order
Ci   mxorb :number of orbital channels for this pair, including
Ci         :lower, intermediate and high blocks.  Only orbitals in the
Ci         :appropriate subblock (defined by set hdpa,hdna,hdpb,hdnb)
Ci         :are added to sk.
Ci   offa  :offset to iprmb array for start of this block, 1st dimension
Ci         :offa is DESTROYED on output
Ci   offb  :offset to iprmb array for start of this block, 2nd dimension
Ci         :offb is DESTROYED on output
Ci   lds   :leading dimension of sk. Its value depends on the complex
Ci          storage format.  For specified 10s digit lbloch, use:
Ci          0  lds = leading dimension of sk
Ci          1  lds = 2
Ci          2  lds = leading dimension of sk
Ci   ld2   :second dimension of sk. Its value depends on the complex
Ci          storage format.  For specified 10s digit lbloch, use:
Ci          0  lds = formal second dimension of sk
Ci          1  lds = formal leading dimension of sk
Ci          2  lds = 2
Ci   klu   :size of sub- and super-diagonal when matrix is stored in
Ci          banded form.
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   hdpa  :last orbital in prior downfolding subblock, first dimension
Ci         :hdpa is zero for lowest subblock; see Remarks
Ci   hdna  :last orbital in this downfolding subblock, first dimension
Ci   hdpb  :last orbital in prior downfolding subblock, second dimension
Ci         :hdpb is zero for lowest subblock; see Remarks
Ci   hdnb  :last orbital in this downfolding subblock,second dimension
Ci   sk    :contribution sk * exp(-i TdotK) is added into s
Ci   nds   :leading dimension of s or sc
Ci   wc    :cos of phase factor -T . q, with T = lattice vector,
Ci         :multiplied by k-point weight
Ci   ws    :sin of phase factor -T . q, with T = lattice vector
Ci         :multiplied by k-point weight
Co Outputs
Co   s,sc,sz :real-space matrix for which to accumulated Bloch sum
Co         :Orbitals are in RL order.  Arrays s,sc,sz hold the same
Co         :data but in different complex formats (see lcmplx above).
Co         :All should point to the same address space.
Cl Local variables
Cl  rcmplx :T when input s is complex
Cl  ldiag  :T when input s is diagonal
Cl  lbnd   :T when banded storage scheme is sought
Cl  kcplx  :10s = complex storage mode type (10s digit of lbloch)
Cl  ipa    :offset to hamiltonian subblock for this pair, augmentation
Cl  ipb    :offset to hamiltonian subblock for this pair, basis
Cl  lma1   :lma loop over lma1..lma2
Cl  lma2   :lma loop over lma1..lma2
Cl  lma    :loops over augmentation orbitals
Cl  lmb    :loops over basis orbitals
Cl  ndima  :augmentation dimension of this downfolding subblock
Cl  ndimb  :basis dimension of this downfolding subblock
Cl  ofbnd  :the additional offset for banded storage
Cr Remarks
Cr    hdpa,hdna,hdpb,hdnb define the range of the downfolding
Cr    subblocks in the first and second dimensions.  Orbitals
Cr    outside this range are not accumulated in the Bloch sum.
Cu Updates
Cu   18 Jul 02 Various changes to accomodate fp input.  New argument list
Cu   23 Jun 02 Various changes to accomodate fp input.  New argument list
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lbloch,mxorb,lds,ld2,offa,offb,nds,iprmb(*),hdpa,
     .  hdna,hdpb,hdnb
      double precision wc,ws,sk(lds,ld2,*)
      double precision s(nds,nds,2),sc(nds,2,nds),sz(2,nds,nds)
C ... Local parameters
      logical rcmplx,lbnd,ldiag,lprmr,ldag
      integer ipa,ipb,kcplx,lma1,lma2,lma,lmb,ndima,ndimb,ofbnd,offa0
      integer isa,isb,offra,offrb,klu

      klu = 0
      lbnd   = mod(lbloch,10) .ne. 0
      kcplx  = mod(lbloch/10,10)
      ldag   = mod(lbloch/100,10) .ne. 0

      rcmplx = mod(mod(lbloch/1000,10),2) .ne. 0
      ldiag  = mod(lbloch/10000,10) .ne. 0
      lprmr  = mod(lbloch/100000,10) .ne. 0
      if (ldiag .and. lprmr)
     .  call rx('iblch1 not ready for ldiag and lprmr')
      ndimb = hdnb - hdpb
      ndima = hdna - hdpa
      offa0 = offa
      ofbnd = 0
C     Range of lma when s is not diagonal
      lma1 = 1
      lma2 = mxorb

C     Case permute orbital order in s: offr[ab] = offset to orbitals
      if (lprmr) then
        do  3  lmb = 1, mxorb
          offrb = iprmb(offb+lmb)
          if (offrb.gt.0 .and. offrb.le.ndimb) goto 4
    3   continue
C       No columns, nothing to copy
        return
    4   continue
        offrb = offrb-1
        do  5  lma = 1, mxorb
          offra = iprmb(offa+lma)
          if (offra.gt.0 .and. offra.le.ndima) goto 6
    5   continue
C       No rows, nothing to copy
        return
    6   continue
        offra = offra-1
      endif

      if (kcplx .eq. 0) then

      if (ldag) call rx('ibloch: not implemented this mode in dagger')

C --- For each basis orbital, do (complex storage mode 0) ---
      do  10  lmb = 1, mxorb
        offb = offb+1
        ipb = iprmb(offb) - hdpb
        if (ipb .le. 0 .or. ipb .gt. ndimb) goto 10
        isb = lmb
        if (lprmr) isb = iprmb(offb) - offrb
        if (lbnd) ofbnd = 2*klu+1-ipb
        if (ldiag) then
          lma1 = lmb
          lma2 = lmb
        endif
        offa = offa0 + lma1-1
C   ... Loop over augmentation orbitals
        do  12  lma = lma1, lma2
          offa = offa+1
          ipa = iprmb(offa) - hdpa
          if (ipa .le. 0 .or. ipa .gt. ndima) goto 12
          isa = lma
          if (lprmr) isa = iprmb(offa) - offra
          s(isa,isb,1) = s(isa,isb,1) +
     .    sk(ipa+ofbnd,ipb,1)*wc - sk(ipa+ofbnd,ipb,2)*ws
          if (rcmplx)
     .    s(isa,isb,2) = s(isa,isb,2) +
     .    sk(ipa+ofbnd,ipb,2)*wc + sk(ipa+ofbnd,ipb,1)*ws
   12   continue
   10 continue

      elseif (kcplx .eq. 1) then

C --- For each basis orbital, do (complex storage mode 1) ---
      do  110  lmb = 1, mxorb
        offb = offb+1
        ipb = iprmb(offb) - hdpb
        if (ipb .le. 0 .or. ipb .gt. ndimb) goto 110
        isb = lmb
        if (lprmr) isb = iprmb(offb) - offrb
        if (lbnd) ofbnd = 2*klu+1-ipb
        if (ldiag) then
          lma1 = lmb
          lma2 = lmb
        endif
        offa = offa0 + lma1-1
C   ... Loop over augmentation orbitals, case real-space s is complex
        if (rcmplx .and. .not. ldag) then
          do  112  lma = lma1, lma2
            offa = offa+1
            ipa = iprmb(offa) - hdpa
            if (ipa .le. 0 .or. ipa .gt. ndima) goto 112
            isa = lma
            if (lprmr) isa = iprmb(offa) - offra
            sz(1,isa,isb) = sz(1,isa,isb) +
     .        sk(1,ipa+ofbnd,ipb)*wc - sk(2,ipa+ofbnd,ipb)*ws
            sz(2,isa,isb) = sz(2,isa,isb) +
     .        sk(2,ipa+ofbnd,ipb)*wc + sk(1,ipa+ofbnd,ipb)*ws
  112     continue
C   ... Loop over augmentation orbitals, case real-space s is real
        elseif (.not. ldag) then
          do  114  lma = lma1, lma2
            offa = offa+1
            ipa = iprmb(offa) - hdpa
            if (ipa .le. 0 .or. ipa .gt. ndima) goto 114
            isa = lma
            if (lprmr) isa = iprmb(offa) - offra
            s(isa,isb,1) = s(isa,isb,1) +
     .        sk(1,ipa+ofbnd,ipb)*wc - sk(2,ipa+ofbnd,ipb)*ws
  114     continue
C   ... Loop over dagger of augmentation orbitals, r.s. s is complex
        elseif (rcmplx) then
          do  116  lma = lma1, lma2
            offa = offa+1
            ipa = iprmb(offa) - hdpa
            if (ipa .le. 0 .or. ipa .gt. ndima) goto 116
            isa = lma
            if (lprmr) isa = iprmb(offa) - offra
            sz(1,isa,isb) = sz(1,isa,isb) +
     .        sk(1,ipb,ipa+ofbnd)*wc - sk(2,ipb,ipa+ofbnd)*ws
            sz(2,isa,isb) = sz(2,isa,isb) +
     .        sk(2,ipb,ipa+ofbnd)*wc + sk(1,ipb,ipa+ofbnd)*ws
  116     continue
C   ... Loop over dagger of augmentation orbitals, r.s. s is real
        else
          do  118  lma = lma1, lma2
            offa = offa+1
            ipa = iprmb(offa) - hdpa
            if (ipa .le. 0 .or. ipa .gt. ndima) goto 118
            isa = lma
            if (lprmr) isa = iprmb(offa) - offra
            s(isa,isb,1) = s(isa,isb,1) +
     .        sk(1,ipb,ipa+ofbnd)*wc - sk(2,ipb,ipa+ofbnd)*ws
  118     continue
        endif

  110 continue

      elseif (kcplx .eq. 2) then

      if (ldag) call rx('ibloch: not implemented this mode in dagger')

C --- For each basis orbital, do (complex storage mode 2) ---
      do  210  lmb = 1, mxorb
        offb = offb+1
        ipb = iprmb(offb) - hdpb
        if (ipb .le. 0 .or. ipb .gt. ndimb) goto 210
        isb = lmb
        if (lprmr) isb = iprmb(offb) - offrb
        if (lbnd) ofbnd = 2*klu+1-ipb
        if (ldiag) then
          lma1 = lmb
          lma2 = lmb
        endif
        offa = offa0 + lma1-1
C   ... Loop over augmentation orbitals, case real-space s is complex
        if (rcmplx) then
          do  212  lma = lma1, lma2
            offa = offa+1
            ipa = iprmb(offa) - hdpa
            if (ipa .le. 0 .or. ipa .gt. ndima) goto 212
            isa = lma
            if (lprmr) isa = iprmb(offa) - offra
            sc(isa,1,isb) = sc(isa,1,isb) +
     .        sk(ipa+ofbnd,1,ipb)*wc - sk(ipa+ofbnd,2,ipb)*ws
            sc(isa,2,isb) = sc(isa,2,isb) +
     .        sk(ipa+ofbnd,2,ipb)*wc + sk(ipa+ofbnd,1,ipb)*ws
  212     continue
C   ... Loop over augmentation orbitals, case real-space s is real
        else
          do  214  lma = lma1, lma2
            offa = offa+1
            ipa = iprmb(offa) - hdpa
            if (ipa .le. 0 .or. ipa .gt. ndima) goto 214
            isa = lma
            if (lprmr) isa = iprmb(offa) - offra
            s(isa,isb,1) = s(isa,isb,1) +
     .        sk(ipa+ofbnd,1,ipb)*wc - sk(ipa+ofbnd,2,ipb)*ws
  214     continue
        endif
  210 continue

      endif

      end

