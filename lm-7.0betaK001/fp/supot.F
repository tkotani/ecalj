#define FP 1
      subroutine supot(mode,sctrl,slat,spot)
C- Initialization for potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make Madelung matrix for monopoles (ASA)
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: nbas nbasp nspin lcd,4
Ci   slat  :struct containing information about the lattice
Ci     Elts read: alat vol opos awald nkd nkq odlv oqlv plat nabc gmax
Ci                nsgrp osymgr oag
Ci     Stored:    ng ogv okv gmax nabc oips0 obgv
Ci   spot  :struct containing information about the potential
Ci     Stored:    omad osmrho osmpot oorhat
Co Outputs
Cr Remarks
Cr   Permanent arrays are created to handle the potential, depending
Cr   on the mode
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode
      double precision sctrl(1),slat(1),spot(1)
C ... Local parameters
      integer nbas,nbasp,nbaspp,nsp,omad,opos,nkd,nkq,odlv,oqlv,igets,
     .  ngabc(3),n1,n2,n3,k1,k2,k3,ngmx,ogv,okv,ng,osmrho,osmpot,oorhat,
     .  oips0,obgv,ngrp,og,oag,iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
      integer w(1)
      common /w/ w

      call info0(50,1,0,' supot : potential setup ... ')

      call upack('ctrl nbas nbasp nspin',sctrl,nbas,nbasp,nsp,0,0)
      call upack('lat alat vol opos',slat,alat,vol,opos,0,0)
      call upack('lat awald nkd nkq odlv oqlv',slat,awald,nkd,nkq,odlv,
     .  oqlv)

      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
        call defdr(omad,nbaspp**2)
        call madmat(nbaspp,w(opos),awald,alat,vol,w(odlv),nkd,
     .    w(oqlv),nkq,w(omad))
        call pack1('pot omad',spot,omad)
      endif

#if FP
C --- Setup for FT charge density, potential representation ---
      if (igets('ctrl lcd,4',sctrl) .ne. 0) then
        call upack('lat alat plat nabc gmax',slat,alat,plat,ngabc,gmax,
     .    0)

C   ... Generate energy cutoff gmax or n1..n3
        call pshpr(iprint()-0)
        call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
        call poppr
        call fftz30(n1,n2,n3,k1,k2,k3)

C   ... Make list of lattice vectors within cutoff
        call defrr(ogv, ngmx*3)
        call defi (okv, ngmx*3)
        call gvlist(alat,plat,w,n1,n2,n3,gmax,8,ngmx,ng,w(okv),
     .    w(ogv),xx,xx)
        if (ng .ne. ngmx) call rx('supot: bug in gvlist')

        call pack5('lat ng ogv okv gmax nabc',slat,ng,ogv,okv,gmax,
     .    ngabc)
C       call prmx('gv',w(ogv),ngmx,ng,3)

        call defcc(osmrho, -k1*k2*k3*nsp)
        call defcc(osmpot, -k1*k2*k3*nsp)
        call defi(oorhat,3*nbaspp)
        call pack5('pot osmrho osmpot oorhat',spot,osmrho,osmpot,oorhat,
     .    0,0)

C   ... Setup for symmetrization
        call defi (oips0,   ng)
        call defcc(obgv,    ng)
        call upack('lat nsgrp osymgr oag',slat,ngrp,og,oag,0,0)
        call sgvsym(ngrp,w(og),w(oag),ng,w(ogv),w(oips0),w(obgv))
        call pack2('lat oips0 obgv',slat,oips0,obgv)

      endif
#endif
      end

