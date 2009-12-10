      subroutine blochi(lbloch,qp,nl,plat,iprmb,is1,is2,iax,nkaps,ndoti,
     .  kaps,nkapk,kapk,s,sd,nds,nitab,ldmpa,ldima,ldmpb,ldimb,ldl,ldl2,
     .  klu,sk)
C- Energy-interpolate and bloch transform real-space structure matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:1s digit pertains to storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci
Ci          10s digit distinguishes how complex arithmetic is handled
Ci           0: sk has real, imaginary separated
Ci              sk = sk(ldl,ldl2,2), with sk(*,*,1..2) = real..imag
Ci           1: sk is returned complex*16 format:
Ci              sk = sk(2,ldl,ldl2), with sk(1..2,*,*) = real..imag
Ci           2: sk has real, imaginary separated by columns
Ci              sk = sk(ldl,2,ldl2), with sk(*,1..2,*) = real..imag
Ci
Ci        100s digit:
Ci           1 if to add to s (s not initialized to zero)
Ci           2 subtract from s
Ci           3 combination of 1+2
Ci
Ci       1000s digit:
Ci           1 if to convert s to spherical harmonics
Ci           2 to restrict s to ib=jb, no translation vector
Ci
Ci      10000s digit nonzero: generate energy derivative in sk
Ci
Ci   qp    :k-point
Ci
Ci   nl    :(global maximum l) + 1
Ci
Ci   plat  :primitive lattice vectors, in units of alat
Ci
Ci   iprmb :permutation indices ordering orbitals in downfolding order.
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
Ci   ndoti :number of derivatives sd to include in energy interpolation
Ci
Ci   nkaps :number of strux s to include in energy interpolation
Ci
Ci   kaps  :energies of strux
Ci
Ci   nkapk :number of energies at which to evaluate Bloch sum
Ci
Ci   kapk  :energies of Bloch-summed strux
Ci
Ci   s     :real-space matrix to be Bloch summed
Ci
Ci   sd    :energy derivative of s
Ci
Ci   nds   :leading dimension of s
Ci
Ci   nitab :third dimension of s
Ci
Ci   ldmpa :offset to first orbital in the augmentation (row) dimension
Ci          for the matrix subblock to be Bloch summed, or equivalently
Ci          the last orbital in the prior subblock (0 for lowest block)
Ci
Ci   ldima :last orbital in the augmentation (row) dimension for the
Ci           matrix subblock to be Bloch summed.  See iprmb, above.
Ci
Ci   ldmpb :offset to first orbital in the source (column) dimension
Ci          for the matrix subblock to be Bloch summed, or equivalently
Ci          the last orbital in the prior subblock (0 for lowest block)
Ci
Ci   ldimb :last orbital in the source (column) dimension for the
Ci           matrix subblock to be Bloch summed.  See iprmb, above.
Ci
Ci   ldl   :leading dimension of sk
Ci
Ci   ldl2  :second dimension of sk
Ci
Ci   klu   :size of sub- and super-diagonal, if s stored banded form
Ci
Co Outputs
Co   sk   :specified block of Bloch summed matrix
Co
Cr Remarks
Cr  *This routine assembles a bloch sum of a real-space matrix, viz
Cr      s(k;r1,l1,r2,l2) = sum_T s(r1,l1,T+r2,l2) * exp(i k . T)
Cr   where r1 and r2 are basis vectors and T = t2-t1 is the difference
Cr   in primitive lattice translation vectors.
Cr
Cr   For pair i,  T is obtained from iax(3..5,i).
Cr
Cr   Contribution from pair i in the iax table may be suppressed
Cr   by setting iax(1,i) or iax(2,i) <= 0
Cr
Cl Local variables
Ci  isite  :current pair over
Cl  lblchi :a local copy of lbloch with digits>1000 stripped
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lbloch,nds,nitab,nl,is1,is2,ldmpa,ldima,ldmpb,ldimb,ldl,
     .  ldl2,niax,klu,iprmb(*),nkaps,nkapk,ndoti
      parameter (niax=10)
      integer iax(niax,is2)
      double precision qp(3),plat(3,3),s(nds,nds,nitab,nkaps),
     .  sd(nds,nds,nitab,nkaps),kaps(nkaps),kapk(nkapk)
C     real + imaginary storage mode
      double precision sk(ldl,ldl2,2,nkapk)
C     double complex storage mode
C     double precision sk(2,ldl,ldl2,nkapk)
C     real + imaginary in columns storage mode
C     double precision sk(ldl,2,ldl2,nkapk)
C Local parameters
      logical onsite,ldot
      integer ia,ib,iprint,isite,j,k,kcplx,ksite,lblchi,ld11,
     .  ld21,ndss,nl2,offa,offb,oi,ik,jk,ip,n0,isw
      character*40 outs
      parameter (ndss=49,n0=10)
      integer iprmk(n0,n0)
      double precision twopi,TdotK,cosT,sinT
      double precision sek(ndss,ndss),sph(ndss,ndss,2)
      double precision Pm(n0,n0),Pn(n0,n0),Pmp(n0,n0),Pnp(n0,n0),
     .                 a(n0,n0),b,c,xx(n0),y,yp

C --- Setup ---
      call tcn('blochi')
      twopi = 8*datan(1d0)
      lblchi = mod(lbloch,10000)
      ldot   = mod(lbloch/10000,10) .ne. 0
      onsite = mod(lblchi/1000,10) .ge. 2
      nl2 = nl*nl
      if (ldima-ldmpa .le. 0 .or. ldimb-ldmpb .le. 0) return
      if (mod(mod(lblchi/100,10),2) .eq. 0) then
        call dpzero(sk,2*ldl*ldl2*nkapk)
      endif
C     Pick up true dimensions of sk,sil,sii from formal ones
      kcplx = mod(lblchi/10,10)
      call cplxdm(kcplx,ldl,ldl2,ld11,ld21,oi,oi)

C --- Setup coefficients for energy interpolation ---
      do  jk = 1, nkapk
C       iprmk orders kaps by increasing distance from kapk(jk)
        do  ik = 1, nkaps
          xx(ik) = kaps(ik) - kapk(jk)
        enddo
        call dvheap(1,nkaps,xx,iprmk(1,jk),0d0,11)
        call hrmins(kaps,iprmk(1,jk),nkaps,min(nkaps,ndoti),kapk(jk),
     .    Pm(1,jk),Pn(1,jk),Pmp(1,jk),Pnp(1,jk),a(1,jk))
      enddo

C --- For each RR' pair, add contribution to Bloch sum ---
      do  100  isite = is1, is2

        ia = iax(2,isite)
        ib = iax(1,isite)
C       Any nonpositive site indices are excluded from sum
        if (ia .le. 0 .or. ib .le. 0) goto 100
        if (onsite) then
          if (ia .ne. ib .or. iax(3,isite) .ne. 0 .or.
     .      iax(4,isite) .ne. 0 .or. iax(5,isite) .ne. 0) goto 100
        endif

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

C   ... Use equivalent of isite to some other site, if it exists
        ksite = isite
        if (iax(8,isite) .ne. 0) ksite = iax(8,isite)

C   --- For each of the kapk, do ---
        do  jk = 1, nkapk

C     --- Energy-interpolate s(ksite) ---
          do  j = 1, nds
          do  k = 1, nds
            y = 0
            yp = 0
            do  ik = 1, min(ndoti,nkaps)
              ip = iprmk(ik,jk)
              b = a(ik,jk)*s(j,k,ksite,ip) - sd(j,k,ksite,ip)
              c = (kaps(ip) - kapk(jk))*b + s(j,k,ksite,ip)
              y = y + Pm(ik,jk)*Pn(ik,jk)*c
              yp = yp - Pm(ik,jk)*Pn(ik,jk)*b +
     .            (Pmp(ik,jk)*Pn(ik,jk)+Pm(ik,jk)*Pnp(ik,jk))*c
            enddo
            do  ik = min(ndoti,nkaps)+1, nkaps
              ip = iprmk(ik,jk)
              y = y + Pm(ik,jk) * Pn(ik,jk) * s(j,k,ksite,ip)
              yp = yp + s(j,k,ksite,ip) *
     .             (Pmp(ik,jk)*Pn(ik,jk) + Pm(ik,jk)*Pnp(ik,jk))
            enddo
            sek(j,k) = y
            if (ldot) sek(j,k) = yp
          enddo
          enddo

C    --- Rotate s to spherical harmonics ---
          if (mod(mod(lblchi/1000,10),2) .ne. 0)
     .      call s2sph(2,nl,nl,sek,ndss,ndss,ndss,ndss,sph)

          offb = nl2*(ib-1)
          offa = nl2*(ia-1)
          call pblch1(lblchi,nl2,offa,offb,ld11,ld21,klu,iprmb,ldmpa,
     .      ldima,ldmpa,ldimb,sek,sph,ndss,cosT,sinT,sk(1,1,1,jk))

        enddo
  100 continue

      if (iprint() .ge. 110) then
        k = 12 + mod(lblchi/10,10)
        do  jk = 1, nkapk
          call awrit2('%xSk%?#n#-dot##, e=%d',outs,len(outs),0,
     .      isw(ldot),kapk(jk))
          call yprm(outs,k,sk(1,1,1,jk),ldima*ldimb,ldl,ldima,ldimb)
        enddo
      endif

      call tcx('blochi')
      end

