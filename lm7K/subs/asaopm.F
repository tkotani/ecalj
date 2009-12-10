      subroutine asaopm(sctrl,sarray,spot,slat,soptic,nl,isp,nsp,nspc,
     .  nclass,nbas,ldimx,nev,qp,ikp,nkp,eband,nbmax,accob,nfilm,nempm,
     .  optmt,velmt)
C- Makes electric dipole and second harmonic optical matrix elements
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas lham loptc
Ci     Stored:   
Ci   sarray:structure containing offsets to various arrays
Ci           sarray->ogrrme offset to matrix elements of 
Ci                          grad phi,phidot; see rgrme.f
Ci     Elts read: oipc
Ci     Stored:   
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: ogrrme
Ci     Stored:   
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nsgrp oistab oag osymgr
Ci     Stored:   
Ci   soptic:struct for parameters for optical ME; see routine uoptic
Ci     Elts read: mode ocrng unrng nchi2 axes
Ci     Stored:   
Ci   nl    :(global maximum l) + 1
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nclass:number of inequivalent classes
Ci   nbas  :size of basis
Ci   ldimx :hamiltonian lower block dimension, (2x dim, noncoll case)
Ci   nev   :actual number of eigenvectors generated
Ci   qp    :k-point
Ci   ikp   :k-point index, used to specify where to poke optmt,velmt
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   eband :energy bands; alias eb (sec*.f)
Ci   nbmax :maximum number of bands
Co Outputs
Co   accob        is transformed to spherical harmonics, and cmplx*16
Co                contains decomposition of norm (aka dnpp)
Co   optmt(*,ikp) ASA <i|grad|j> connecting occ i with unocc j
Co   velmt(*,ikp) ASA <i|grad|i>
Cr Remarks
Cr   Optics package adapted from Sergey Rashkeev with Walter Lambrecht,
Cr   which was adapted from an earlier version by V. Antropov.
Cu Updates
Cu   22 Feb 04 Redefined flags in loptc switch
Cu   20 Nov 02 (jek) extended to metals
Cu   20 Dec 00 (wrl) extended to noncollinear case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nl,nsp,nclass,ldimx,nbmax,ikp,nkp,nfilm,nempm,isp,
     .  nbas,nspc,nev
      double precision sctrl(1),sarray(1),spot(1),slat(1),soptic(1),
     .  qp(3,*),eband(nbmax)
      double precision accob(0:nl*nl-1,nbas,nspc,2,2,nev)
      double precision optmt(3,nfilm,nempm,nsp/nspc,nkp),
     .                 velmt(3,nfilm,nsp/nspc,nkp)

      call rx('ASAOPM not installed; optics library required')

      end

