      subroutine intdrv(sctrl,soptic,sbz,eband,nbmax,nsp,nspc,efermi,
     .  idtet,vol,nfilm,nempm,optmt)
C- Driver for SHG Integrals
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: loptc
Ci     Stored:   *
Ci     Passed to:*
Ci   soptic:struct containing parameters for optical ME
Ci     Elts read: mode ne window ocrng unrng esciss nchi2 axes
Ci     Stored:   *
Ci     Passed to:*
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: nkabc nkp ntet
Ci     Stored:   *
Ci     Passed to:*
Ci   eband :energy bands for irreducible part of BZ
Ci   nbmax :leading dimension of eband
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   efermi:Fermi level
Ci   idtet :(0,i) no. of tetrahedra of the i'th kind
Ci         :(1-4,i) identifies the i'th tetrahedron in terms of
Ci         :four irreducible k-points
Ci   vol   :volume
Ci   nfilm,nempm: dimensions optmt
Ci   optmt :<i|grad|j> connecting occ i with unocc j
Co Outputs:
Co   Im(eps) or joint density of states written to file 'optdf'
Cr Remarks
Cr   Adapted from bzints to make joint density of states or Im(eps)
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cu Updates
Cu   02 Mar 01 Added scissors operator
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbmax,nsp,nspc,idtet(0:4,*),npts,nfilm,nempm
      double precision sctrl(1),sbz(1),soptic(1),efermi,vol
      double precision optmt(3,nfilm,nempm,nsp/nspc,*),
     .                 eband(nbmax,nsp,*)
    
C Local parameters
      logical bittst
      integer nfilo,nfiup,nemlo,nemup,nkp,ntet
      integer ifio,fopn,loptic,mode
      integer ocrng(2),unrng(2)
      integer nkabc(3),n1,n2,n3
      double precision optrng(2),emin,emax
      equivalence (emin,optrng(1)),(emax,optrng(2))
      double precision esciss
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      integer nchi2, iabc(3,6)

      call rx('INTDRV not installed; optics library required')


      end
      

