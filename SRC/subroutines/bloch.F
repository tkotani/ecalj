      subroutine bloch(qp,nl,plat,mxorb,iprmb,is2,iax,s,
     .nds,isp,nsp, ldl, sll)!,sil,sii)
C- Bloch transform of real-space matrix
C ----------------------------------------------------------------------
Ci Inputs
c      lbloch = 100000 + 4000 + 40*1 + 10   is fixed.
c      lblchi = 4000 + 40*1 + 10
Ci   lbloch:1s digit pertains to storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Cixxx           1: s is stored in banded form (see Remarks)
Ci
Ci          10s digit distinguishes how complex arithmetic is handled
Cixxx           0: sll has real, imaginary separated
Cixxx              sll = sll(ldl,ldl2,2), with sll(*,*,1..2) = real..imag
Ci           1: sll is returned complex*16 format:
Ci              sll = sll(2,ldl,ldl2), with sll(1..2,*,*) = real..imag
Cixxx           2: sll has real, imaginary separated by columns
Cixxx              sll = sll(ldl,2,ldl2), with sll(*,1..2,*) = real..imag
Cixxx              By default, input s is real
Cixxx           4: Input s is in complex*16 format:
Cixxx              This bit can be taken in combination with bits 1,2
Cixxx              NB: these conventions apply to sll, sil, sii
Ci
Ci        100s digit:
Ci           0 copy to sll, sil or sii (ie initialize array to zero)
Cixx           1 add to sll, sil or sii (ie do not initialize to zero)
Cixx           2 Make Bloch transform of (-s)
Cixx           4 scale Bloch phase factor (i k . T) by -1
Cixx             Any combination of the above is allowed
Ci
Ci       1000s digit:
Cixx           1 if to convert s to spherical harmonics
Cixx           2 to restrict s to ib=jb, no translation vector
Ci           4 if to use transpose of s(r1,l1,T+r2,l2) in place of s
Ci             Any combination of the above is allowed
Ci
Ci      10000s digit pertains to which of sll,sil,sii are generated
Ci           0 generate sll, sil, sii
Cixx           1 suppress generation of sll
Cixx           2 suppress generation of sil
Cixx           4 suppress generation of sii
Ci             switches 1,2,4 may be combined
Ci
Ci     100000s digit pertains to whether s has permuted orbital order
Cixxx           0 if s has normal order
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
      double precision sll(ldl,ldl,2)!,sil(ldi,ldl2,2),sii(ldi,ldi,2)
C ... Local parameters
      integer ia,ib,iprint,isite,j,k,kcplx,ksite,ld11,
     .ld12,ld13,ld21,ld22,ld23,lidim,nds1,ndss,offa,offb !,oi,scplx
      double precision TdotK,cosT,sinT
      double precision sc(nds,2,nds),swk(nds,2,nds)
      logical:: debug=.false.
      real(8),parameter:: twopi = 8*datan(1d0)
C --- Setup ---
      call tcn('bloch')
      ndss = nds
      nds1 = ndss
c      lidim = ldima+idim
      sll=0d0
c      sil=0d0
c      sii=0d0
c      ld11 = 2
c      ld21 = ldl
c      ld12 = 2
c      ld22 = ldi
c      ld13 = 2
c      ld23 = ldi
C --- For each RR' pair, add contribution to Bloch sum ---
      do  100  isite = 1, is2
        ia = iax(2,isite)
        ib = iax(1,isite)
        if (ia <0 .or. ib<0) goto 100
        TdotK = twopi*sum(qp*matmul(plat,iax(3:5,isite)))
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
C   ... Use equivalent of isite to some other site, if it exists
        ksite = isite
        if (iax(8,isite) .ne. 0) ksite = iax(8,isite)
C-- Handle all cases when s is to be copied to sc --- ... Input s is complex
        j = 2*ksite-1
        call zmscop(0,nds,nds,nds,ndss,0,0,0,0,s(1,1,2*isp-1,j),sc)
        call ztoy(sc,ndss,nds,nds,0)
        sc(:,1,:)=transpose(sc(:,1,:))
        sc(:,2,:)=transpose(sc(:,2,:))
C   --- Lower-lower block  Sll ---
        offa = mxorb*(ib-1)     !ltrans mode
        offb = mxorb*(ia-1)
        call pblch1(mxorb,offa,offb,2,ldl,0,iprmb,
     .       ldl,sc,nds1,cosT,sinT,sll)
  100 continue
      call tcx('bloch')
      end subroutine bloch

      subroutine pblch1(mxorb,offa,offb,lds,ld2,klu,iprmb,!lblchp,
     .hdn,sc,nds,cosT,sinT,sk)
C- Contribution of one pair to Bloch sum of strux
C ----------------------------------------------------------------------
Ci Inputs 101010
Ci   lblchp: 1s digit concerns storage of Bloch summed sk
Cix           0: sk is stored in normal, unpacked form
Cix           1: sk is stored in banded form (see Remarks)
Cix              The band form follows LAPACK band storage conventions:
Cix              sk(i,j) is stored in location (kl+ku+1+i-j,j)
Cix              with kl,ku = size of sub- and super-diagonal.
Cix              Here we take kl=ku=klu.
Ci          10s digit distinguishes how complex arithmetic is handled
Cixxx           0: sk has real, imaginary separated
Cixxx              sk = sk(ldl,ldl2,2), with sk(*,*,1..2) = real..imag
Ci on          1: sk is returned complex*16 format:
Ci               sk = sk(2,ldl,ldl2), with sk(1..2,*,*) = real..imag
Cixxx            2: sk has real, imaginary separated by columns
Cixxx                sk = sk(ldl,2,ldl2), with sk(*,1..2,*) = real..imag
Cixxx            100s digit is not used
Ci on  1000s digit 1 if real-space s is complex
Ci              In this case, sc is used in place of s.
Cixxx          10000s digit 1 if s(or sc) is diagonal.
Ci on 100000s digit 1 if s(or sc) has a permuted orbital order
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
      implicit none
      integer lblchp,mxorb,lds,ld2,offa,offb,nds,klu,iprmb(*),hdpa,
     .hdna,hdpb,hdn
      double precision cosT,sinT,sk(lds,ld2,2),sc(nds,2,1)
      logical rcmplx,lbnd,ldiag,lprmr
      integer ipa,ipb,kcplx,lma1,lma2,lma,lmb,ndima,ndimb,ofbnd,offa0
      integer isa,isb,offra,offrb
c      lbnd   = .false. !mod(lblchp,10) .ne. 0
c      kcplx  =  1 !mod(lblchp/10,10)
c      rcmplx = .true. !mod(mod(lblchp/1000,10),2) .ne. 0
c      ldiag  = .false. !mod(lblchp/10000,10) .ne. 0
c      lprmr  = .true. !mod(lblchp/100000,10) .ne. 0
c      if (ldiag .and. lprmr) call rx('pblch1 not ready for ldiag and lprmr')
      ndimb = hdn !- hdpb
      ndima = hdn !- hdpa
      offa0 = offa
c      ofbnd = 0
c      lma1 = 1
c      lma2 = mxorb
C     Case permute orbital order in s: offr[ab] = offset to orbitals
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
        return
    6   continue
        offra = offra-1
c        print *,'oooooffr=',offra,offrb
C --- For each basis orbital, do (complex storage mode 1) ---
        do  110  lmb = 1, mxorb
          offb = offb+1
          ipb = iprmb(offb) !- hdpb
          isb = iprmb(offb) - offrb
          if (ipb .le. 0 .or. ipb .gt. ndimb) cycle
          offa = offa0 !+ lma1-1
C   ... Loop over augmentation orbitals, case real-space s is complex
            do  112  lma = 1, mxorb
              offa = offa+1
              ipa = iprmb(offa) !- hdpa
              isa = iprmb(offa) - offra
              if (ipa .le. 0 .or. ipa .gt. ndima) cycle
              sk(1,ipa,ipb) = sk(1,ipa,ipb) +  sc(isa,1,isb)*cosT - sc(isa,2,isb)*sinT
              sk(2,ipa,ipb) = sk(2,ipa,ipb) +  sc(isa,1,isb)*sinT + sc(isa,2,isb)*cosT
  112       continue
  110   continue
      end subroutine pblch1
