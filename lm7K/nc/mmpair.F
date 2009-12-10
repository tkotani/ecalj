      subroutine mmpair(nbas,nttab,iax,ebil,ibil,aamom,e,
     .  etot,amag,pcorr,f)
C- Pairwise forces on all spins, micromagnetics hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ebil  :bilinear coefficients
Ci   ibil  :index to pair correlation function for this pair
Ci   e     :direction vector, cartesian coordinates
Ci   aamom :magnitude of local moments
Co Outputs
Co   amag  :system magnetization along the three direction vectors
Co   pcorr :pair correlation function
Co   f     :bilinear contribution to forces are added to f
Co   etot  :bilinear contribution to total energy added to etot
Cr Remarks
Cr
Cu Updates
Cu   26 Nov 02 amag is scaled by aamom
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax
      parameter (niax=10)
      integer nbas,nttab,iax(niax,1),ibil(1)
      double precision f(nbas,3),ebil(nttab),e(nbas,3),pcorr(0:1),
     .  aamom(nbas)
C ... Local parameters
      integer i,i1,i2
      double precision xx1,xx2,xx3,etot,amag(3),etpair

      etpair = 0

C --- Bilinear term ---
      do  10  i = 1, nttab
        i1 = iax(1,i)
        i2 = iax(2,i)
c        if (i2 .ne. 3) goto 10
        xx1 = ebil(i)*e(i2,1)
        xx2 = ebil(i)*e(i2,2)
        xx3 = ebil(i)*e(i2,3)
        f(i1,1) = f(i1,1) + xx1
        f(i1,2) = f(i1,2) + xx2
        f(i1,3) = f(i1,3) + xx3
        etpair = etpair + xx1*e(i1,1) + xx2*e(i1,2) + xx3*e(i1,3)
        pcorr(ibil(i)) = pcorr(ibil(i)) +
     .    e(i1,1)*e(i2,1)+ e(i1,2)*e(i2,2) + e(i1,3)*e(i2,3)
   10 continue
C ... Add into etot, remembering double-counting
      etot = etot + etpair/2

C --- Average magnetization ---
      call dpzero(amag,3)
      do  20  i = 1, nbas
        amag(1) = amag(1) + aamom(i)*e(i,1)
        amag(2) = amag(2) + aamom(i)*e(i,2)
        amag(3) = amag(3) + aamom(i)*e(i,3)
   20 continue
      amag(1) = amag(1)/nbas
      amag(2) = amag(2)/nbas
      amag(3) = amag(3)/nbas

      end

