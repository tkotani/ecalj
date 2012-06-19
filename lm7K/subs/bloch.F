      subroutine bloch(lbloch,qp,nl,plat,mxorb,iprmb,is1,is2,iax,s,
     .nds,isp,nsp,ldima,ldimb,idim,ldl,ldi,ldl2,klu,sll,sil,sii)

C- Bloch transform of real-space matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:1s digit pertains to storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci
Ci          10s digit distinguishes how complex arithmetic is handled
Ci           0: sll has real, imaginary separated
Ci              sll = sll(ldl,ldl2,2), with sll(*,*,1..2) = real..imag
Ci           1: sll is returned complex*16 format:
Ci              sll = sll(2,ldl,ldl2), with sll(1..2,*,*) = real..imag
Ci           2: sll has real, imaginary separated by columns
Ci              sll = sll(ldl,2,ldl2), with sll(*,1..2,*) = real..imag
Ci              By default, input s is real
Ci           4: Input s is in complex*16 format:
Ci              This bit can be taken in combination with bits 1,2
Ci              NB: these conventions apply to sll, sil, sii
Ci
Ci        100s digit:
Ci           0 copy to sll, sil or sii (ie initialize array to zero)
Ci           1 add to sll, sil or sii (ie do not initialize to zero)
Ci           2 Make Bloch transform of (-s)
Ci           4 scale Bloch phase factor (i k . T) by -1
Ci             Any combination of the above is allowed
Ci
Ci       1000s digit:
Ci           1 if to convert s to spherical harmonics
Ci           2 to restrict s to ib=jb, no translation vector
Ci           4 if to use transpose of s(r1,l1,T+r2,l2) in place of s
Ci             Any combination of the above is allowed
Ci
Ci      10000s digit pertains to which of sll,sil,sii are generated
Ci           0 generate sll, sil, sii
Ci           1 suppress generation of sll
Ci           2 suppress generation of sil
Ci           4 suppress generation of sii
Ci             switches 1,2,4 may be combined
Ci
Ci     100000s digit pertains to whether s has permuted orbital order
Ci           0 if s has normal order
Ci           1 if s has a permuted orbital order.
Ci             In this case, the first row (column) in each R' (R)
Ci             block of s corresponds to the first (permuted) orbital
Ci             associated with site R' (R).
Ci
Ci   qp    :k-point
Ci
Ci   nl    :(global maximum l) + 1
Ci         :only used converting s to s. harmonics; see Bugs
Ci
Ci   plat  :primitive lattice vectors, in units of alat
Ci
Ci   mxorb :leading dimension of iprmb
Ci
Ci   iprmb :permutation indices ordering orbitals in sll, sil, sii
Ci         :in downfolding order.
Ci          Unpermuted orbitals run from 1 ... n, with indices to
Ci          orbitals corresponding to site ib starting at 1+nl*nl*(ib-1).
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
Ci   s     :real-space matrix to be Bloch summed
Ci
Ci   nds   :leading dimension of s
Ci
Ci   ldima :cutoff for lower set of orbitals in the augmentation
Ci          (row) dimension.  See iprmb, above.
Ci
Ci   ldimb :cutoff for lower set of orbitals, in the source
Ci          (column) dimension.  See iprmb, above.
Ci
Ci   idim  :dimension of intermediate set. See iprmb, above.
Ci
Ci   ldl   :leading dimension of sll
Ci
Ci   ldi   :leading and second dimension of sii
Ci
Ci   ldl2  :second dimension of sll and sil
Ci
Ci   klu   :size of sub- and super-diagonal, if s stored banded form
Ci
Co Outputs
Co   sll   :lower-lower block of Bloch summed matrix
Co
Co   sil   :lower-intermediate block of Bloch summed matrix
Co
Co   sii   :intermediate-intermediate block of Bloch summed matrix
Co
Cr Remarks
Cr  *This routine assembles a bloch sum of a real-space matrix, viz
Cr     s(k;r1,l1,r2,l2) = sum_T s(r1,l1,T+r2,l2) * exp(i k . T)
Cr   where r1 and r2 are basis vectors and T = t2-t1 is the difference
Cr   in primitive lattice translation vectors.
Cr
Cr   For pair i,  T is obtained from iax(3..5,i).
Cr
Cr   Contribution from pair i in the iax table may be suppressed
Cr   by setting iax(1,i) or iax(2,i) <= 0
Cr
Cl Local variables
Cl   isite  :index to current pair
Cl   lblchi :a local copy of lbloch with digits>1000 stripped
Cl   lblchp :a local copy of lbloch suitable for pblch1
Cl   lsph   :T rotate strux to spherical harmonics (add 1000 to lbloch)
Cl   onsite :T sum only diagonal parts of s (add 2000 to lbloch)
Cl   ltrans :T Bloch sum of transpose of s (add 4000 to lbloc)
Cl   ndss   :dimension of sc.
Cl          :NB: with automatic arrays, ndss is always just nds.
Cl   scplx  :Input R.S. s is in complex*16 format
Cl
Cb Bugs
Cb   conversion to spherical harmonics assumes simple ordering
Cb   s,p,d,.. of s.  This should be changed, argument nl eliminated.
Cu Updates
Cu   10 Oct 03 Rotation to s-harm can be for complex s
Cu   09 May 03 Bloch transform s_transpose (lbloch 4000)
Cu   30 Mar 03 Switch for Bloch phase = exp(-i k.T)
Cu   11 Jan 03 Bug fix complex s, spin polarized case
Cu   18 Jul 02 Additional changes to accomodate fp input
Cu   23 Jun 02 Various changes to accomodate fp input
Cu             and input s to be complex*16 format.
Cu             New argument list.
Cu   20 Jul 99 Routine was revised with a changed argument list.
Cu             The original bloch was renamed to blcho.
Cu   17 Dec 99 Added on-site restriction (1000's digit)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lbloch,nds,nl,is1,is2,ldima,ldimb,idim,ldl,ldi,ldl2,niax,
     .klu,mxorb,iprmb(*),isp,nsp
      parameter (niax=10)
      integer iax(niax,is2)
      double precision qp(3),plat(3,3)
      double precision s(nds,nds,nsp,is2)
C     real + imaginary storage mode
      double precision sll(ldl,ldl2,2),sil(ldi,ldl2,2),sii(ldi,ldi,2)
C     complex*16 storage mode
C     double precision sll(2,ldl,ldl2),sil(2,ldi,ldl2),sii(2,ldi,ldi)
C     real + imaginary in columns storage mode
C     double precision sll(ldl,2,ldl2),sil(ldi,2,ldl2),sii(ldi,2,ldi)
C ... Local parameters
      logical bittst,maskii,maskil,maskll,onsite,ltrans,lsph
      integer ia,ib,iprint,isite,j,k,kcplx,ksite,lblchi,lblchp,ld11,
     .ld12,ld13,ld21,ld22,ld23,lidim,nds1,ndss,offa,offb,oi,scplx
      double precision twopi,TdotK,cosT,sinT
      double precision sc(nds,2,nds),swk(nds,2,nds)
C     parameter (ndss=49)
C     double precision sc(ndss,ndss,2),swk(ndss,2,ndss)
C ... External calls
      external cplxdm,dpzero,pblch1,rx,s2sph,tcn,tcx,yprm,zmscop,ztoy
C ... Heap
C     integer w(1)
C     common /w/ w

C --- Setup ---
      call tcn('bloch')
      ndss = nds
      twopi = 8*datan(1d0)
      lblchi = mod(lbloch,10000)
      lsph   = mod(mod(lblchi/1000,10),2) .ne. 0
      onsite = mod(mod(lblchi/1000,10),4) .ge. 2
      ltrans = mod(lblchi/1000,10) .ge. 4
      j = mod(lbloch/10000,10)
      maskll = bittst(j,1)
      maskil = bittst(j,2)
      maskii = bittst(j,4)
      lidim = ldima+idim
      if (mod(mod(lblchi/100,10),2) .eq. 0) then
        if (.not. maskll) call dpzero(sll,2*ldl*ldl2)
        if (.not. maskil) call dpzero(sil,2*ldi*ldl2)
        if (.not. maskii) call dpzero(sii,2*ldi**2)
      endif

C     Pick up true dimensions of sll,sil,sii from formal ones
      kcplx = mod(mod(lblchi/10,10),4)
      scplx = mod(lblchi/10,10)/4
      call cplxdm(kcplx,ldl,ldl2,ld11,ld21,oi,oi)
      call cplxdm(kcplx,ldi,ldl2,ld12,ld22,oi,oi)
      call cplxdm(kcplx,ldi,ldi,ld13,ld23,oi,oi)
C     lblchp passed to pblch1. s cplx -> shift lblchi flag 40->1000
      lblchp = mod(lblchi,1000)
      if (scplx .eq. 1) lblchp = lblchp - 40
      if (scplx .eq. 1 .or. ltrans .or. lsph) lblchp = lblchp + 1000
C     s has permuted orbital order
      j = mod(lbloch/100000,10)
      lblchp = lblchp + j*100000

C     Leading and second dimension of s as used by pblch1
      if (mod(mod(lblchp/1000,10),2) .ne. 0) then
C       If pblch1 uses sc, it's ndss
        nds1 = ndss
      else
C       If pblch1 uses s, it's nds
        nds1 = nds
      endif

C --- For each RR' pair, add contribution to Bloch sum ---
      do  100  isite = is1, is2

        ia = iax(2,isite)
        ib = iax(1,isite)
C       Any nonpositive site indices are excluded from sum
        if (ia .le. 0 .or. ib .le. 0) goto 100
        if (onsite) then
          if (ia .ne. ib .or. iax(3,isite) .ne. 0 .or.
     .    iax(4,isite) .ne. 0 .or. iax(5,isite) .ne. 0) goto 100
        endif

        TdotK = 0
        do  30  j = 1, 3
        do  30  k = 1, 3
          TdotK = TdotK + twopi*qp(j)*plat(j,k)*iax(2+k,isite)
   30   continue
        if (mod(lblchi/100,10) .ge. 4) then
          TdotK = -TdotK
        endif
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
        if (mod(mod(lblchi/100,10),4) .ge. 2) then
          cosT = -cosT
          sinT = -sinT
        endif

C   ... Use equivalent of isite to some other site, if it exists
        ksite = isite
        if (iax(8,isite) .ne. 0) ksite = iax(8,isite)

C   --- Handle all cases when s is to be copied to sc ---
C   ... Input s is complex
        if (scplx .eq. 1) then
          if (mod(lblchp/1000,10) .ne. 1) call rx('oops!')
          j = 2*ksite-1
C         Rotate complex s to spherical harmonics
          if (lsph) then
            call zmscop(0,nds,nds,nds,ndss,0,0,0,0,s(1,1,2*isp-1,j),swk)
            call ztoyy(swk,ndss,ndss,ndss,ndss,1,0)
            if (nl**2.ne.ndss) call rx('bloch not ready for sph. harm')
            call s2sph(12,nl,nl,swk,nds,nds,ndss,ndss,sc)
C         Direct copy
          else
            call zmscop(0,nds,nds,nds,ndss,0,0,0,0,s(1,1,2*isp-1,j),sc)
C            if (isite .eq. 13) then
C            call zprm('s',2,sc,ndss,nds,nds)
C            endif
C           call zprm('s',2,sc,ndss,nds,nds)
            call ztoy(sc,ndss,nds,nds,0)
C            if (isite .eq. 13) then
C            print *, 'ksite=',ksite
C            call yprm('s',4,sc,ndss,ndss,nds,nds)
C            endif
          endif
C   ... Rotate s to spherical harmonics
        else if (lsph) then
          if (mod(lblchp/1000,10) .ne. 1) call rx('oops!')
          call s2sph(2,nl,nl,s(1,1,isp,ksite),nds,nds,ndss,ndss,sc)

C   ... If otherwise to use transpose of s
        else if (ltrans) then
          if (mod(lblchp/1000,10) .ne. 1) call rx('oops!')
          call rx('bloch not ready for scplx=1 and ltrans')
        endif
C   ... Take transpose of sc
        if (ltrans) then
C         print *, isite
C         call yprm('s before transpose',4,sc,ndss,ndss,nds,nds)
          call ymscop(0,nds,nds,2*ndss,2*ndss,0,0,0,0,sc,ndss,swk,ndss)
          call ymtrns(1,swk,2*ndss,1,ndss,sc,2*ndss,1,ndss,1,nds,1,nds)
C         call yprm('copy of s',4,swk,ndss,ndss,nds,nds)
C         call yprm('s after transpose',4,sc,ndss,ndss,nds,nds)
        endif

C      if (iax(1,isite) .eq. 1 .and. iax(2,isite) .eq. 1) then
C        print 338, isite, sc(1,1,1)
C  338   format(i4,f12.6)
C      endif

C   --- Lower-lower block ---
        if (ldima .ne. 0 .and. .not. maskll) then
          offb = mxorb*(ib-1)
          offa = mxorb*(ia-1)
          if (ltrans) then
            offa = mxorb*(ib-1)
            offb = mxorb*(ia-1)
          endif
          call pblch1(lblchp,mxorb,offa,offb,ld11,ld21,klu,iprmb,0,
     .    ldima,0,ldimb,s(1,1,isp,ksite),sc,nds1,cosT,sinT,sll)
        endif

C   --- Intermediate-lower block ---
        if (idim .ne. 0 .and. .not. maskil) then
          offb = mxorb*(ib-1)
          offa = mxorb*(ia-1)
          call pblch1(lblchp,mxorb,offa,offb,ld12,ld22,klu,iprmb,ldima,
     .    lidim,0,ldimb,s(1,1,isp,ksite),sc,nds1,cosT,sinT,sil)
        endif

C   --- Intermediate-intermediate block ---
        if (idim .ne. 0 .and. .not. maskii) then
          offb = mxorb*(ib-1)
          offa = mxorb*(ia-1)
          call pblch1(lblchp,mxorb,offa,offb,ld13,ld23,klu,iprmb,ldima,
     .    lidim,ldimb,lidim,s(1,1,isp,ksite),sc,nds1,cosT,sinT,sii)
        endif

C       call snott(isite,sll,ld21,cost)
  100 continue

      if (iprint() .ge. 110) then
        k = mod(mod(lblchi/10,10),4)
        if (.not. maskll)
     .  call yprm('bloch: Sll',12+k,sll,ldima*ldimb,ldl,ldima,ldimb)

        if (idim .gt. 0) then
          if (.not. maskil)
     .    call yprm('bloch: Sil',2+k,sil,idim*ldimb,ldi,idim,ldimb)
          if (.not. maskii)
     .    call yprm('bloch: Sii',12+k,sii,idim*idim,ldi,idim,idim)
        endif
      endif

      call tcx('bloch')
      end subroutine bloch


      subroutine pblch1(lblchp,mxorb,offa,offb,lds,ld2,klu,iprmb,
     .hdpa,hdna,hdpb,hdnb,s,sc,nds,cosT,sinT,sk)

C- Contribution of one pair to Bloch sum of strux
C ----------------------------------------------------------------------
Ci Inputs
Ci   lblchp: 1s digit concerns storage of Bloch summed sk
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
Ci        100s digit is not used
Ci       1000s digit 1 if real-space s is complex
Ci             In this case, sc is used in place of s.
Ci      10000s digit 1 if s(or sc) is diagonal.
Ci     100000s digit 1 if s(or sc) has a permuted orbital order
Ci   mxorb :number of orbital channels for this pair, including
Ci         :lower, intermediate and high blocks.  Only orbitals in the
Ci         :appropriate subblock (defined by set hdpa,hdna,hdpb,hdnb)
Ci         :are added to sk.
Ci   offa  :offset to iprmb array for start of this block, 1st dimension
Ci   offb  :offset to iprmb array for start of this block, 2nd dimension
Ci   lds   :leading dimension of sk. Its value depends on the complex
Ci          storage format.  For specified 10s digit lblchp, use:
Ci          0  lds = leading dimension of sk
Ci          1  lds = 2
Ci          2  lds = leading dimension of sk
Ci   ld2   :second dimension of sk. Its value depends on the complex
Ci          storage format.  For specified 10s digit lblchp, use:
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
Ci   s     :real-space matrix for which to accumulated Bloch sum
Ci          orbitals in s are in RL order
Ci   sc    :same as s, but in complex form of the kcplx=2 type.
Ci          orbitals in sc are in RL order
Ci   nds   :leading dimension of s or sc
Ci   cosT  :cos of phase factor T . q, with T = lattice vector
Ci   sinT  :sin of phase factor T . q, with T = lattice vector
Co Outputs
Co   sk    :contribution from this pair is added into sk(offa,offb)
Co          s * exp(i TdotK) is added into sk for this matrix subblock
Co          orbitals in sk  are in downfolding order
Cl Local variables
Cl  rcmplx :T when input s is complex
Cl  ldiag  :T when input s is diagonal
Cl  lbnd   :T when banded storage scheme is sought
Cl  kcplx  :10s = complex storage mode type (10s digit of lblchp)
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
      integer lblchp,mxorb,lds,ld2,offa,offb,nds,klu,iprmb(*),hdpa,
     .hdna,hdpb,hdnb
      double precision cosT,sinT,sk(lds,ld2,2),s(nds,1),sc(nds,2,1)
C ... Local parameters
      logical rcmplx,lbnd,ldiag,lprmr
      integer ipa,ipb,kcplx,lma1,lma2,lma,lmb,ndima,ndimb,ofbnd,offa0
      integer isa,isb,offra,offrb

      lbnd   = mod(lblchp,10) .ne. 0
      kcplx  = mod(lblchp/10,10)
      rcmplx = mod(mod(lblchp/1000,10),2) .ne. 0
      ldiag  = mod(lblchp/10000,10) .ne. 0
      lprmr  = mod(lblchp/100000,10) .ne. 0
      if (ldiag .and. lprmr)
     .call rx('pblch1 not ready for ldiag and lprmr')
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
C   ... Loop over augmentation orbitals, case real-space s is complex
          if (rcmplx) then
            do  12  lma = lma1, lma2
              offa = offa+1
              ipa = iprmb(offa) - hdpa
              if (ipa .le. 0 .or. ipa .gt. ndima) goto 12
              isa = lma
              if (lprmr) isa = iprmb(offa) - offra
              sk(ipa+ofbnd,ipb,1) = sk(ipa+ofbnd,ipb,1) +
     .        sc(isa,1,isb)*cosT - sc(isa,2,isb)*sinT
              sk(ipa+ofbnd,ipb,2) = sk(ipa+ofbnd,ipb,2) +
     .        sc(isa,1,isb)*sinT + sc(isa,2,isb)*cosT
   12       continue
C   ... Loop over augmentation orbitals, case real-space s is real
          else
            do  14  lma = lma1, lma2
              offa = offa+1
              ipa = iprmb(offa) - hdpa
              if (ipa .le. 0 .or. ipa .gt. ndima) goto 14
              isa = lma
              if (lprmr) isa = iprmb(offa) - offra
              sk(ipa+ofbnd,ipb,1) = sk(ipa+ofbnd,ipb,1) + s(isa,isb)*cosT
              sk(ipa+ofbnd,ipb,2) = sk(ipa+ofbnd,ipb,2) + s(isa,isb)*sinT
   14       continue
          endif
   10   continue

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
          if (rcmplx) then
            do  112  lma = lma1, lma2
              offa = offa+1
              ipa = iprmb(offa) - hdpa
              if (ipa .le. 0 .or. ipa .gt. ndima) goto 112
              isa = lma
              if (lprmr) isa = iprmb(offa) - offra
              sk(1,ipa+ofbnd,ipb) = sk(1,ipa+ofbnd,ipb) +
     .        sc(isa,1,isb)*cosT - sc(isa,2,isb)*sinT
              sk(2,ipa+ofbnd,ipb) = sk(2,ipa+ofbnd,ipb) +
     .        sc(isa,1,isb)*sinT + sc(isa,2,isb)*cosT
  112       continue
C   ... Loop over augmentation orbitals, case real-space s is real
          else
            do  114  lma = lma1, lma2
              offa = offa+1
              ipa = iprmb(offa) - hdpa
              if (ipa .le. 0 .or. ipa .gt. ndima) goto 114
              isa = lma
              if (lprmr) isa = iprmb(offa) - offra
              sk(1,ipa+ofbnd,ipb) = sk(1,ipa+ofbnd,ipb) + s(isa,isb)*cosT
              sk(2,ipa+ofbnd,ipb) = sk(2,ipa+ofbnd,ipb) + s(isa,isb)*sinT
  114       continue
          endif
  110   continue

      elseif (kcplx .eq. 2) then

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
              sk(ipa+ofbnd,1,ipb) = sk(ipa+ofbnd,1,ipb) +
     .        sc(isa,1,isb)*cosT - sc(isa,2,isb)*sinT
              sk(ipa+ofbnd,2,ipb) = sk(ipa+ofbnd,2,ipb) +
     .        sc(isa,1,isb)*sinT + sc(isa,2,isb)*cosT
  212       continue
C   ... Loop over augmentation orbitals, case real-space s is real
          else
            do  214  lma = lma1, lma2
              offa = offa+1
              ipa = iprmb(offa) - hdpa
              if (ipa .le. 0 .or. ipa .gt. ndima) goto 214
              isa = lma
              if (lprmr) isa = iprmb(offa) - offra
              sk(ipa+ofbnd,1,ipb) = sk(ipa+ofbnd,1,ipb) + s(isa,isb)*cosT
              sk(ipa+ofbnd,2,ipb) = sk(ipa+ofbnd,2,ipb) + s(isa,isb)*sinT
  214       continue
          endif
  210   continue

      endif

      end subroutine pblch1


