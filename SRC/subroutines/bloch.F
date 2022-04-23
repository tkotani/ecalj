      subroutine bloch2(qp,isp, sll)
      use m_rdsigm2,only: hrr,ndimsig,sfz,nk1,nk2,nk3
      use m_mkqp,only: bz_nabc
      use m_hamindex, only: ib_table
      use m_gennlat,only: nlatS,nlatE,nlat
      use m_lattic,only: plat=>lat_plat
      use m_mpitk,only: master_mpi,procid
!     !  Bloch transform of real-space matrix
!! sll(j,i) =\sum_Ts hhrs(i,j,isp,T)*exp(ikT)
      implicit none
      integer:: isp,iset
      complex(8):: sll(ndimsig,ndimsig)
      real(8),parameter:: twopi = 8*datan(1d0)
      complex(8):: phase,img=(0d0,1d0)
      real(8):: qp(3)
      integer:: ik1,ik2,ik3,i,j,ib1,ib2,ix,nS,nE,nnn
      call tcn('bloch2')
      write(6,*) 'bbbbb bloch2',procid
      sll=0d0
cccccccccccccccccccccccccccc
      do i=1,ndimsig
         ib1 = ib_table(i)      !atomic-site index in the primitive cell
c         write(6,*)'bbbbbbbbbb',i,ib1
      enddo
      nnn=nk1*nk2*nk3
      do ik1=1,nk1
      do ik2=1,nk2
      do ik3=1,nk3
        do i=1,ndimsig !ndimMTO
        do j=1,ndimsig !ndimMTO
         ib1 = ib_table(i) !atomic-site index in the primitive cell
         ib2 = ib_table(j)
         nS=nlatS(ik1-1,ik2-1,ik3-1,ib1,ib2)
         nE=nlatE(ik1-1,ik2-1,ik3-1,ib1,ib2)
         phase=0d0
         do ix=nS,nE
         phase=phase+exp(-img*twopi*sum(qp*matmul(plat,nlat(:,ix,ib1,ib2))))
         enddo
         phase= phase/(nE-nS+1) 
c         if(master_mpi) print *,'ns ne=',ns,ne,phase,
         sll(i,j) = sll(i,j) + hrr(ik1,ik2,ik3,i,j,isp)/nnn*phase
      enddo
      enddo
      enddo
      enddo
      enddo
      call tcx('bloch2')
      end subroutine bloch2
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine bloch(qp,nl,plat,mxorb,iprmb,is2,iax,hhrs,
     .nds,isp,nsp, ldl, sll)!,sil,sii)
!!  Bloch transform of real-space matrix
!! sll(j,i) =\sum_Ts hhrs(i,j,isp,T)*exp(ikT)
      implicit none
      integer lbloch,nds,nl,is1,is2,ldima,ldimb,idim,ldl,ldi,ldl2,niax,klu,mxorb,iprmb(*),isp,nsp
      parameter (niax=10)
      integer iax(niax,is2)
      double precision qp(3),plat(3,3)
      complex(8):: hhrs(nds,nds,nsp,*)
      complex(8):: sll(ldl,ldl)!,sil(ldi,ldl2,2),sii(ldi,ldi,2)
      integer ia,ib,iprint,isite,j,k,kcplx,ksite,ld11,
     .ld12,ld13,ld21,ld22,ld23,lidim,nds1,ndss,offa,offb !,oi,scplx
      real(8),parameter:: twopi = 8*datan(1d0)
      complex(8):: phase,img=(0d0,1d0)
      integer:: ipa,ipb,lma,lmb,offa0, isa,isb,offra,offrb
      call tcn('bloch')
      sll=0d0
      do  100  isite = 1, is2   !For each RR' pair, add contribution to Bloch sum ---
         ia = iax(2,isite) !sll(Rb,Ra) colmn
         ib = iax(1,isite) !row
         if (ia <0 .or. ib<0) goto 100
         phase= exp (img*twopi*sum(qp*matmul(plat,iax(3:5,isite))) )
         ksite = isite          !Use equivalent of isite to some other site, if it exists
         if (iax(8,isite) .ne. 0) ksite = iax(8,isite)
         offb = mxorb*(ib-1)  !row !offa = mxorb*(ib-1)  
         offa = mxorb*(ia-1)  !column   offb = mxorb*(ia-1)
ccccccccccccccccccccccccccccccccccccccccccccccc
         offa0 = offa
         do  3  lmb = 1, mxorb
            offrb = iprmb(offb+lmb)
            if (offrb.gt.0 .and. offrb.le.ldl) goto 4
 3       continue
C     No columns, nothing to copy
         return
 4       continue
         offrb = offrb-1
         
         do  5  lma = 1, mxorb
            offra = iprmb(offa+lma)
            if (offra.gt.0 .and. offra.le.ldl) goto 6
 5       continue
         return
 6       continue
         offra = offra-1
         
c     print *,'oooooffr=',offra,offrb
         do  110  lmb = 1, mxorb
            offb = offb+1
            ipb = iprmb(offb)   !- hdpb
            isb = iprmb(offb) - offrb
            if (ipb .le. 0 .or. ipb .gt. ldl) cycle
            offa = offa0        !+ lma1-1
            do  112  lma = 1, mxorb
               offa = offa+1
               ipa = iprmb(offa) !- hdpa
               isa = iprmb(offa) - offra
               if (ipa .le. 0 .or. ipa .gt. ldl) cycle
               sll(ipb,ipa) = sll(ipb,ipa) +  hhrs(isb,isa,isp,ksite)*phase
 112        continue
 110     continue
ccccccccccccccccccccccccccccccccccccccccccccccc
 100  continue
      call tcx('bloch')
      end subroutine bloch

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
