      subroutine asaddq(sctrl,slat,sham,spot,sarray,ikp,nkp,
     .  isp,nev,nbmax,eband,zll,
     .  nlibu,lmaxu,idu,vorb,
     .  nfilem,lwtkb,wtkb,wtkp,zval,metal,ldnpp,
     .  nqpp,ldos,lrhos,nsite,lsite,ldmat,qnu,qpp,dnpp,rhos,orbtm)
C- Accumulate contributions from this qp byproducts of ASA eigenvectors
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: lncol lham loptc nbas nclass nl nspec nspin
Ci     Stored:    *
Ci
Ci   slat  :struct containing information about the lattice
Ci     Elts read: avw alat plat ocg ojcg oidxcg nsgrp oistab oag osymgr
Ci                vol
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: ldham oindxo oeula neula qss
Ci     Stored:    *
Ci     Passed to: *
Ci
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opp
Ci     Stored:    *
Ci     Passed to: *
Ci
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: oipc onrc
Ci     Stored:    *
Ci     Passed to: *
Ci   ikp   :index to current k-point
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   isp   :current spin channel (1 or 2)
Ci   nev   :actual number of eigenvectors generated
Ci   nbmax :maximum number of bands
Ci   eband :energy bands; alias eb (sec*.f)
Ci   zll   :ll block of eigenvectors
Ci   nfilem:logical unit for moments file.
Ci         :If nfilem=0, nothing is written to file
Ci   lwtkb :>0 accumulate moments (and rhos) directly using weights wtkb;
Ci         :   accmom, accsm are not written to disk.
Ci   wtkb  :integration weights, used if lwtkb is set (bzwts.f)
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   zval : valence charge
Ci   metal :T => metal, F => nonmetal
Ci   ldnpp :T, return phi,phid decomposition of the norm
Ci   nqpp  :>0, accumulate qpp; then nqpp is leading dimension of qpp
Ci   ldos  :a combination of the following integers
Ci         : 1 make dos
Ci         : 2 generate weights for partial dos
Ci         : 4 generate weights for m-decompos'n of pdos
Ci   lrhos :T, accumulate accsm for site diagonal spin-density matrix
Ci   lsite :list of sites for which to accumulate partial dos.  In
Ci         :this mode, dos weights are stored by site, not by class
Ci   nsite :number of sites in lsite. :nsite=0 => no list and
Ci         :asaddq uses default mode for storing for storing doswts.
Ci   ldmat :>0, make density-matrix (not implemented)
Co Outputs
Co   qnu   :energy-weighted moments of the sphere charges
Co   qpp   :wave function products, needed for nonspherical density.
Co   dnpp  :phi,phid decomposition of the norm
Co   rhos  :spin density-matrix
Cl Local variables
Cl  w(ozil):il block of eigenvectors (secmat.f)
Cl   w(oc) :small h * eigenvectors (secmat.f)
Cl w(odoswt):partial DOS weights
Cr Remarks
Cr   asaddq was designed for the ASA, and is oriented towards quantities
Cr   that depend on R and l but not m.  But the m-resolved weights may
Cr   needed, e.g. for ldos=4.  These were patched in;
Cr   at present, asaddq does not properly symmetrize the dos weights
Cr   and the caller must avoid symmetry operations when making m-resolved
Cr   weights.
Cr
Cr   The moments file containing the partial weights can either contain
Cr   partial weights resolved by Rl (the default) or Rlm; it can also
Cr   contain just the weights, or energy moments of the weights.  The
Cr   latter is needed if doswts are later used to make the moments qnu.
Cr   Only if partial dos resolved by m are needed (ldos=4) will the
Cr   doswts by resolved by m.  Additional `spin' channels may be required
Cr   in the  noncollinear case.  Subroutine iomomn decides on how many
Cr   channels need be defined in a particular context, and documents what
Cr   the channels are for.  Subroutine makwts creates the weights;
Cr   subroutine iomoms does the I/O for them.
Cr
Cu Updates
Cu   04 Jun 08 Output density in SS with downfolding case
Cu   27 May 08 generate orbtm for case lwtkb>0
Cu   08 Nov 07 (J. Xu) LDA+U implementation; qpp is complex
Cu   03 Nov 04 No longer impose lorth condition for two-C hamiltonian
Cu   21 Apr 04 Additions for an m-dependent spin-density matrix
Cu   14 Feb 03 bug fix: correctly dimension dosw2
Cu    9 Jan 03 nfilem=0 -> nothing written to moments file
Cu   18 Jan 02 Redesigned to handle accumulation of m-resolved weights
Cu   17 Apr 01 Added switch lkkrm
Cu   03 Feb 01 logical dos changed into ldos
Cu   19 Jan 01 extended to handle downfolded orbitals, noncollinear case
Cu   20 Dec 00 dnpp may be accumulated for the noncollinear case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical metal,lrhos,ldnpp
      integer isp,nev,ikp,nkp,nfilem,nbmax,nqpp,ldmat,ldos,lwtkb,
     .  nsite,lsite(*)
      double precision sctrl(1),sham(1),spot(1),slat(1),sarray(1)
      double precision eband(nbmax,*),zll(*),wtkb(nbmax),dnpp(*),
     .  wtkp(nkp),zval,qnu(1),rhos(2,*),orbtm(*)
      double complex qpp(nqpp,4,*)
      integer nlibu,lmaxu,idu(4,*)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
C ... Local parameters
      logical bittst,lorth,lkkrm,lgors,llwtkb
      integer o,odosws,odosw2,odoswt,oc,oeula,oidxsh,oipc,onrc,opp,opph,
     .  owk,owk2,owk3,ozil,oqnu
      integer i,idim,l2,ldham(16),ldim,ldimx,lham,li,lidim,
     .  lihdim,lihdmx,lncol,loptic,nbas,nclass,neul,nl,nlo,nsp,nspc,
     .  nmom,nspec,iw,nfstg,j,iomoms,nschan,ng,nchan,nchan2,nrhos
      double precision qss(4),xx
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (ldimx,ldham(5)),(lihdmx,ldham(7)),(nspc,ldham(4))
C ... Heap
      integer w(1)
      common /w/ w

C     No eigenvectors; nothing to accumulate
      if (nev .le. 0) return

C --- Setup ---
      call tcn('asaddq')
      call defi(o,1)
      call upack('ctrl lncol lham loptc',sctrl,lncol,lham,loptic,0,0)
      lkkrm = .not. lgors('ctrl lasa,256',sctrl)
      call upack('ctrl nbas nclass nl nspec nspin',sctrl,nbas,nclass,nl,
     .  nspec,nsp)
      call upack('ham ldham oindxo',sham,ldham,oidxsh,0,0,0)
      call upack('ham oeula neula qss',sham,oeula,neul,qss,0,0)
      call upack('array oipc onrc',sarray,oipc,onrc,0,0,0)
      call upack('pot opp nrhos',spot,opp,nrhos,0,0,0)
C     Some dimensioning parameters
      idim = lidim - ldim
      li = ldim * idim * nspc**2
      l2 = ldimx**2
C     lorth not necessary for twoc if use 3C secmt2 to assemble H
C     lorth = bitand(lham,3) .ne. 0 .or. bittst(lham,32)
      lorth = bittst(lham,32)
      nlo = nl
      if (ldos .ge. 4) nlo = nl*nl
      llwtkb = lwtkb .gt. 0
      nmom = nspc
      if (lwtkb .eq. 0 .and. metal) nmom = nspc+2
      call iomomn(metal,ldos,lrhos,nev,nspc,nmom,lwtkb,nfstg)
      nschan = mod(nfstg/10,10)
      call upack('lat nsgrp',slat,ng,0,0,0,0)
      if (nsite .gt. 0) then
        nchan = nlo*nbas
        nchan2 = nrhos*nbas
        if (lrhos) call rx('ASADDQ: --pdos not implemented with lrhos')
      else
        nchan = nlo*nclass
        nchan2 = nrhos*nclass
      endif
      call rxx(ldos.gt.1.and..not.metal,
     .  ' Switch metal=F incompatible with partial dos')

C ... Vector of potential parameters
      call defdr(opph,5*lihdim*nsp)
      call makpph(nl,nsp,nbas,lihdim,w(oipc),w(oidxsh),w(opp),w(opph))
C     call prmx('pph before adding V(LDA+U)',w(opph),5,5,2*lihdim)
      call u2pph(0,nbas,lmaxu,nsp,nlibu,idu,vorb,w(oidxsh),ldim,
     .  lihdim,1d0,w(opph))
C     call prmx('pph after adding V(LDA+U)',w(opph),5,5,2*lihdim)

      if (lncol .ne. 0) call saiwav(ldim,lidim,5,lihdim,w(opph))

C --- Make zil and h * zll ---
      if (.not. lorth) then
        call defdc(ozil,max(li,1))
        call defdc(oc,l2)
        call asadq1(nbas,nl,w(oidxsh),ldim,idim,lihdim,lncol,w(oeula),
     .    neul,qss,isp,nev,w(opph),zll,w(ozil),w(oc))
      else
        call defrr(ozil,1)
        call defrr(oc,1)
      endif

C ... Copy pph(spin 2) to pph(spin 1)
      if (isp .eq. 2)
     .  call dpscop(w(opph),w(opph),5*lihdim,1+5*lihdim,1,1d0)

C --- Decompose evecs into sites; accumulate data spec'd by caller ---
      call defdc(owk, lihdmx)
      call defdc(owk2,lihdmx)
      call defdc(owk3,lihdmx)
C     iw = offset in wtkb to element wtkb(ikp,isp)
      iw = 1 + ldim*((isp-1)+nsp*(ikp-1))

C ... Allocate memory for doswt
      call defdr(odoswt,-nchan*nev*max(nschan,nspc**2))
      i = 1
      if (lrhos) i = nchan2*nev*nspc*3
      call defdc(odosw2,-i)

C ... doswts for specified list: make for all sites, then copy those sought
      if (nsite .ne. 0) then
        call defrr(oqnu,-3*nl*nsp*nbas)
        call makwts(metal,llwtkb,lrhos,ldos,lkkrm,lorth,nl,nlo,nrhos,
     .    isp,nsp,nspc,nbas,ikp,nkp,nev,wtkp,zval,w(opph),ldim,idim,
     .    lihdim,w(oidxsh),ldnpp,dnpp,nqpp,qpp,zll,w(ozil),w(oc),0,nbas,
     .    w,w(owk),w(owk2),w(owk3),eband,wtkb(iw),w(oqnu),w(odoswt),
     .    w(odosw2),w)

        nchan = nlo*nsite
        call defdr(odosws,-nchan*nev*nschan)
        call asadq4(nl,nlo,isp,nsp,nbas,nev,nschan,nsite,w(oipc),
     .    w(onrc),lsite,w(odoswt),w(oqnu),w(odosws),qnu)
        i = odosws
        odosws = odoswt
        odoswt = i

C ... doswts for all classes
      else
        call makwts(metal,llwtkb,lrhos,ldos,lkkrm,lorth,nl,nlo,nrhos,
     .    isp,nsp,nspc,nbas,ikp,nkp,nev,wtkp,zval,w(opph),ldim,idim,
     .    lihdim,w(oidxsh),ldnpp,dnpp,nqpp,qpp,zll,w(ozil),w(oc),
     .    w(oipc),nclass,w(onrc),w(owk),w(owk2),w(owk3),eband,wtkb(iw),
     .    qnu,w(odoswt),w(odosw2),rhos)

        if (llwtkb .and. nspc .eq. 2) then
          xx = 0
          call moment(2,nl,nlo,nrhos,isp,nsp,nspc,nkp,ldim,nev,nev,ikp,
     .      wtkb,nclass,w(odoswt),w(odosw2),qnu,rhos,orbtm,xx)
        endif

      endif

C ... Save energy bands and possibly doswt to disk
      if (nfilem .ne. 0) then
        j = iomoms(-nfilem,nl,nsp,nspc,2,ldim,nfstg,nschan,1,1,nev,nev,
     .    nchan,nchan2,nev,eband,w(odoswt),w(odoswt),w(odosw2),xx,xx)
      endif

C --- Accumulate density matrix  ---
      if (ldmat .ne. 0) then
        call rx('asaddq not ready to make density matrix')
      endif

      call rlse(o)
      call tcx('asaddq')
      end

      subroutine asadq4(nl,nlo,isp,nsp,nbas,nev,nschan,nsite,ipc,nrclas,
     .  lsite,doswt,qnub,doswts,qnu)
C- Copies dos weights ordered by classes to one ordered by a list of sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl    :number of l channels per site
Ci   nlo   :number of dos (l- or lm-) channels per site
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis; doswt has been made for all sites
Ci   nev   :number of eigenvectors generated
Ci   nschan:number of 'spin channels' in doswt
Ci   nsite :number of sites in subset
Ci   nrclas:nrclas(i) = number of atoms in the ith class
Ci   lsite :list of sites in subset
Ci   doswt :dos weights, for all sites
Ci   qnub  :energy-weighted moments of the sphere charges,
Ci         :for all sites
Co Outputs
Co   doswts:dos weights for list of sites lsite
Co   qnu   :moments by class, averaged by sites
Cu Updates
Cu   17 Jan 02
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nl,nlo,isp,nsp,nbas,nev,nschan,nsite,lsite(nsite),ipc(1),
     .  nrclas(*)
      double precision
     .  doswt(nlo,nbas,nev,nschan),qnub(3,nl,nsp,nbas),
     .  doswts(nlo,nsite,nev,nschan),qnu(3,nl,nsp,*)
C ... Local parameters
      integer i,ib,ic,j
      double precision fac

C ... Copy doswts for list of sites
      do  i = 1, nsite
        if (lsite(1) .le. 0) then
          ib = i
        else
          ib = lsite(i)
        endif

        do  j  = 1, nev*nschan
          call dcopy(nlo,doswt(1,ib,j,1),1,doswts(1,i,j,1),1)
        enddo
      enddo

C ... Accumulate qnu by class
      do  ib = 1, nbas
        ic = ipc(ib)
        fac = 1d0/nrclas(ic)
        call daxpy(3*nl,fac,qnub(1,1,isp,ib),1,qnu(1,1,isp,ic),1)
      enddo

      end

      subroutine asadq1(nbas,nl,indxsh,ldim,idim,lihdim,lncol,eula,
     .  neul,qss,isp,nev,pph,zll,zil,c)
C- Kernel called by asaddq to make zil and c = h * zll
C     implicit none
      integer nbas,nl,indxsh,ldim,idim,lihdim,lncol,nev,neul
      double precision zll(ldim,*),zil(idim,2,nev),pph(5,lihdim,*),
     .  c(ldim,*),eula(nbas,neul,3),qss(4)
C Local
      logical bittst,lss
      integer fopnT,nsss,idimx,iprint,isp,l2,ldimx,li,ncsw,nfilet,
     .  nspc,od1,osil,osilss,osll,od2,ozwk,nevw,lidim,mode
C     Used to dimension workspace, size ldim*nevw*nspc**2*2
      parameter (nevw=64)
      double precision dum
C ... Heap
      integer w(1)
      common /w/ w

      nspc = 1
      if (lncol .ne. 0) nspc = 2
      li = ldim*idim*nspc**2
      lss = bittst(lncol,2)
      ldimx = ldim*nspc
      idimx = idim*nspc
      l2 = ldimx**2
C     nfilet = fopn('TMP')
      nfilet = fopnT('TMP',-1,4,0)
      rewind nfilet
      lidim = ldim+idim

C --- Make zil= H_il*z_ll; sil must be in file tmp ---
      nsss = 1
      if (lss) nsss = 2
      if (idim .ne. 0) then
C       if (lss) call rx('ASADDQ: not ready for downfolding with SS')
        call tcn('zil')
        call defdc(osil,li)
C   ... Read in sil
        call dpdump(w(osil),idim*ldim*nsss*2,nfilet)
C       call yprm('Read Sil',02,w(osil),idim*ldim,idim,idim,ldim)
        if (lncol .ne. 0) then
C         Skip past iwaves contribution to h
          call dpdump(dum,1,nfilet)
          call defrr(od2,ldimx+idimx)
          call defdc(ozwk,ldimx*nevw)
          
          if (lss) then
            call defdc(osilss,-li)
            call rotspn(50000,nbas,nl,indxsh,eula,neul,qss(4),
     .        w,w,1,lidim,ldim,idim,ldim,w(osil),w(osilss))
C           call yprm('Silss',2,w(osilss),idimx*ldimx,idimx,idimx,ldimx)
            mode = 11
          else
            mode = 1
          endif
          call asadq3(mode,ldim,idim,lihdim,nev,nevw,nbas,nl,indxsh,
     .      eula,neul,qss,pph,w(osil),w(osilss),zll,w(ozwk),w(od2),zil)
          call rlse(od2)
        else
          call makrcz(ldim,idim,nev,pph(1,1,isp),w(osil),zll,zil)
        endif
        if (iprint() .ge. 110) then
          call yprm('zil',2,zil,idimx*ldimx,idimx,idimx,nev)
        endif
        call rlse(osil)
        call tcx('zil')
      endif

C --- Make c = small h * eigenvectors ---
      call tcn('h*z')
C ... Read in sll
      call defdc(osll,l2)
      call dpdump(w(osll),ldim*ldim*nsss*2,nfilet)
      if (lncol .ne. 0) then
#if NC
        if (lss) then
          ncsw = 3000 + 20000
          call dcopy(l2,w(osll),1,c,1)
        else
          ncsw = 3000
          call dcopy(ldim**2*2,w(osll),1,c,1)
        endif
C       Brute-force branch.  Makes h directly; then h * zll
C       Works for both SS and otherwise
C       Note: with 4 Jun08 revisions, probably not needed.
        if (lss) then
          call defcc(od1,-ldimx)
          call defdr(od2,-ldimx)
C         d2 <- srdel
          call asadq2(ldim,lihdim,3,1,1d0,pph,w(od2))
C         d1 <- C-enu
          call asadq2(ldim,lihdim,2,2,1d0,pph,w(od1))
          call asadq2(ldim,lihdim,1,2,-1d0,pph,w(od1))
C         Make h -> w(osll).  Note c contains sll
          call rotspn(ncsw,nbas,nl,indxsh,eula,neul,qss(4),
     .      w(od1),w(od2),ldim,ldim,ldim,ldim,ldim,c,w(osll))
C         call yprm('h',2,w(osll),l2,ldimx,ldimx,ldimx)
          call zmpy(w(osll),ldimx,1,l2,zll,ldimx,1,l2,c,ldimx,1,l2,
     .      ldimx,nev,ldimx)
          call rlse(od1)
C       Fast branch.  Not implemented for SS.
        else
          call defrr(od2,ldimx)
          call defdc(ozwk,ldimx*nevw)
          call asadq3(0,ldim,ldim,lihdim,nev,nevw,nbas,nl,indxsh,eula,
     .      neul,qss,pph,w(osll),w,zll,w(ozwk),w(od2),c)
          call rlse(od2)
        endif
#endif
      else
        call makdsd(1,ldim,ldim,ldim,ldim,0,0,pph(1,1,isp),w(osll),
     .    w(osll))
        call zmpy(w(osll),ldimx,1,l2,zll,ldimx,1,l2,c,ldimx,1,l2,
     .    ldimx,nev,ldimx)
      endif
      if (iprint() .ge. 110) then
        call yprm('h*z',2,c,l2,ldimx,ldimx,nev)
      endif
      call rlse(osll)
      call tcx('h*z')

      end
      subroutine asadq2(ldim,lihdim,ipp,icplx,fac,pph,ppv)
C- Adds fac * (ldim part of pph) into vector ppv
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimension of pph
Ci   ipp   :which element in pph to copy
Ci   icplx :1 if ppv is real, 2 if ppv is complex
Ci   fac   :add fac*pph into ppv
Ci   pph   :potential parameters in downfolding order (makpph.f)
Co Outputs
Co   ppv   :fac * one column of pp's added into vector ppv
C ----------------------------------------------------------------------
C     implicit none
      integer ipp,ldim,lihdim,icplx
      double precision fac,pph(5,lihdim,2),ppv(icplx,ldim,2)

      call daxpy(ldim,fac,pph(ipp,1,1),5,ppv,icplx)
      call daxpy(ldim,fac,pph(ipp,1,2),5,ppv(1,1,2),icplx)
      end

      subroutine asadq3(mode,ldim,idim,lihdim,nev,nevw,nbas,nl,indxsh,
     .  eula,neul,qss,pph,sil,silss,zll,zwk,d2,zil)
C- Generate decorated strux * zll for noncollinear eigenvectors zll
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 strux are sil block, and output is zil; see Remarks
Ci         :1 strux are sll block, and output is h*zll; see Remarks
ci         :10 same as mode 0, but for SS.
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (mode=0)
Ci         :same as ldim (mode=1)
Ci   lihdim:number of lower+intermediate+higher orbitals; dimensions pph
Ci   nev   :number of eigenvectors for which to generated zil
Ci   nevw  :size of work array zwk; see zwk
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   qss   :spin spiral parameters (not now implemented)
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :iwaves are OVERWRITTEN
Ci   sil   :unrotated structure structure constants for:
Ci         :intermediate-lower block (mode=0)
Ci         :lower-lower block, aka sll (mode=1)
Ci   silss :unrotated structure structure constants for SS (lss=T)
Ci         :intermediate-lower block (mode=0)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   zll   :eigenvectors 1..nev.
Ci         :zll=zll(ldim,2,ldim,2,2)
Ci   zwk   :complex work array holding some portion of modified zll
Ci         :zwk = zwk(ldim,2,nevw)
Ci         :nevw may be chosen any size >=1; recommended nevw=64.
Ci         :Larger nevw demands more memory; small nevw runs slower
Ci   d2
Co Outputs
Co   zil   :eigenvectors for downfolded orbitals (mode=0)
Co         :c = h*zll (mode=1)
Cl Local variables
Cl   ozl   :offset of imaginary part to zll
Cl   ozi   :offset of imaginary part to zi
Cl   ozwi  :offset of imaginary part to zwk
Cl   offi  :offset of i block (ldim,0 for mode=0,1)
Cl   offh  :last element in i block (ldim+idim,ldim for mode=0,1)
Cr Remarks
Cr   The noncollinear case is analogous to the collinear case, except
Cr   that the structure constants sil are rotated by spinor U.  Thus:
Cr   (mode=0) The intermediate eigenvectors are geenerated from zll as
Cr     zil = di Ui+ Sil Ul dl zll                 (see makrcz.f)
Cr   where di and dl are the same as in the collinear case; see makrcz.
Cr   (mode=1) The ll product c = h*z is
Cr       c = (dl Ul+ Sll Ul dl + C-enu) zll       (see makwts.f)
Cr            dl Ul+ Sll Ul dl zll + (C-enu) zll
Cr   The most efficient implementation is to evaluate matrix products
Cr   right-to-left, eg   zil = di (Ui+ (Sil (Ul (dl zll))))
Cr   because it takes advantage of the spin-diagonal character of sil.
Cr   This routine does the operations in that order for both modes.
Cr   As a consequence, it is not suitable for SS.
Cu Updates
Cu   04 Jun 08 Make Zil for SS case
Cu   19 Jan 01 first implementation, used to make nocollinear zil
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ldim,idim,lihdim,nevw,nbas,nl,indxsh(*),neul,nev
      double precision pph(5,lihdim,2),d2(ldim,2),
     .  sil(idim,ldim,2),zll(ldim,2,nev),zwk(ldim,2,nevw,2),
     .  silss(idim,2,ldim,2,2),eula(nbas,neul,3),qss(4),zil(idim,2,nev)
C ... Local parameters
      logical lss
      integer ldimx,nspc,j,idimx,nevi,iev,ozi,ozwi,offi,offh,isp,i,
     .  ozl,mode0
      double precision xx

      mode0 = mod(mode,10)
      lss = mode .ge. 10
      nspc = 2
      idimx = nspc*idim
      ldimx = nspc*ldim
c     ozl is offset of imaginary part to zll
      ozl = ldimx*ldimx
c     ozi is offset of imaginary part to zi
      ozi = idimx*ldimx
c     ozwi is offset of imaginary part to zwk
      ozwi = ldimx*nevw
c     offi is offset of i block relative to first block
c     offh is last element in i block
      if (mode0 .eq. 0) then
        offi = 0
        offh = ldim
      else
        offi = ldim
        offh = ldim+idim
        do  j = ldim+1, ldim+idim
          pph(3,j,1) = pph(3,j,1) / ( pph(2,j,1) - pph(1,j,1) )
          pph(3,j,2) = pph(3,j,2) / ( pph(2,j,2) - pph(1,j,2) )
        enddo
      endif
      call dpzero(d2,ldimx)
      call asadq2(ldim,lihdim,3,1,1d0,pph,d2)

C ... For each block of eigenvectors sized nevi, do
      do  iev = 1, nev, nevw
      nevi = min(nevw,nev-iev+1)

C     Copy piece of eigenvector to zwk
      call ymscop(0,ldimx,nevi,ldimx,ldimx,0,iev-1,0,0,zll,ozl,zwk,ozwi)
C     call yprm('zll(i)',2,zwk,ozwi,ldimx,ldimx,nevi)
C ... zwk <- U dl zll
      call rotspn(42000,nbas,nl,indxsh,eula,neul,qss(4),xx,d2,
     .  nevi,ldim,0,ldim,nevw,xx,zwk)
C     call yprm('U d zll',2,zwk,ozwi,ldimx,ldimx,nevi)
C ... zil <- sil U dl zll
      if (lss) then
C       call yprm('silss',2,silss,idimx*ldimx,idimx,idimx,ldimx)
        call yygemm('N','N',idim*2,nevi,ldim*2,1d0,silss,
     .    silss(1,1,1,1,2),idim*2,zwk,zwk(1,1,1,2),ldimx,0d0,
     .    zil(1,1,iev),zil(1+ozi,1,iev),idimx)
      else
C       call yprm('sil',2,sil,idim*ldim,idim,idim,ldim)
        call yygemm('N','N',idim,nevi,ldim,1d0,sil,sil(1,1,2),idim,zwk,
     .    zwk(1,1,1,2),ldimx,0d0,zil(1,1,iev),zil(1+ozi,1,iev),idimx)
        call yygemm('N','N',idim,nevi,ldim,1d0,sil,sil(1,1,2),idim,
     .    zwk(1,2,1,1),zwk(1,2,1,2),ldimx,0d0,
     .    zil(1,2,iev),zil(1+ozi,2,iev),idimx)
      endif
C     call yprm('sil U dl zll',2,zil(1,1,iev),ozi,idimx,idimx,nevi)
C ... zil <- Ui+ sil U dl zll
      call rotspn(40100,nbas,nl,indxsh,eula,neul,qss(4),xx,xx,
     .  nevi,offh,offi,idim,ldimx,xx,zil(1-offi,1,iev))
C     call yprm('U+ sil U dl zll',2,zil(1,1,iev),ozi,idimx,idimx,nevi)
C ... zil overwritten by true zil = di U+ sil U dl zll for il case
C                           or by   dl U+ sll U dl zll for ll case
      call makdsd(20,idimx,ldimx,idim,nevi,0,offi,pph,zil(1,1,iev),
     .  zil(1,1,iev))
      call makdsd(20,idimx,ldimx,idim,nevi,0,offi,pph(1,1,2),
     .  zil(1,2,iev),zil(1,2,iev))
C     call yprm('zil(iev)',2,zil(1,1,iev),ozi,idimx,idimx,nevi)

C ... ll case: add (C-enu) zll to
      if (mode0 .eq. 0) then
        do  isp = 1, 2
          do  i = 1, ldim
            xx = pph(2,i,isp) - pph(1,i,isp)
            do  j = iev, iev+nevi-1
              zil(i,isp,j) = zil(i,isp,j) + xx*zll(i,isp,j)
              zil(i+ozi,isp,j) = zil(i+ozi,isp,j) + xx*zll(i+ozl,isp,j)
            enddo
          enddo
        enddo

      endif

      enddo
C     call yprm('zil',2,zil,ozi,idimx,idimx,nev)
      end

C      subroutine asadq4(ssite,sspec,slat,nbas,nclass,nsp,nl,nev,nspc,
C     .  qnub,doswtb,qnu,doswt,dosw2)
CC- Symmetrize the dos weights
CC ----------------------------------------------------------------------
CCi Inputs
CCi   ssite :struct containing site-specific information
CCi   sspec :struct containing species-specific information
CCi   slat  :struct containing information about the lattice
CCi   nbas  :size of basis
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCi Inputs/Outputs
CCio  orhoat:vector of offsets containing site density
CCio        :Symmetrized on output
CCio  qbyl  :site- and l-decomposed charges
CCio        :Symmetrized on output
CCio  hbyl  :site- and l-decomposed one-electron energies
CCio        :Symmetrized on output
CCio  f     :forces
CCio        :Symmetrized on output
CCr Remarks
CCu Updates
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nbas,nsp,nl,nev,nspc,nclass
C      double precision ssite(1),sspec(1),slat(1),
CC    .  f(3,nbas),qbyl(n0,nsp,nbas),hbyl(n0,nsp,nbas),
C     .  qnub(3,nl,nsp,nbas),qnu(3,nl,nsp,*),
C     .  doswtb(nl,nl,nbas,nev,nspc,nspc),
C     .  doswt(nl,nl,nclass,nev,nspc,nspc),dosw2(2,nl,nl,nclass,nev,3)
CC ... Local parameters
C      integer oag,ocg,ocy,og,oidxcg,oipa,oipc,oips,ojcg,opos,opos0,orho,
C     .  osym,owk
C      integer stdo,ib0,ic,ipr,iprint,is,lgunit,
C     .  ngrp,nlml,nlmx,nr,nrclas,igetss,mxint,oistab,ival,n1,n2
C      double precision plat(9),qlat(9)
CC ... Heap
C      integer w(1)
C      common /w/ w
C
C      call tcn('asadq4')
C      stdo = lgunit(1)
C      ipr = iprint()
C
C      call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
C      call upack('lat plat qlat oistab',slat,plat,qlat,oistab,0,0)
C      call upack('lat nsgrp osymgr oag',slat,ngrp,og,oag,0,0)
C
C      call defi (oips, nbas)
C      call defi (oipc, nbas)
C      call defrr(opos0,3*nbas)
C
C      call spackv(10,'site class',ssite,1,nbas,w(oipc))
C      call spackv(10,'site pos',ssite,1,nbas,w(opos0))
C      if (nclass .ne. mxint(nbas,w(oipc)))
C     .  call rx('asaddq : class mismatch')
C
CC --- Start loop over classes ---
C      call defi  (oipa,    nbas)
C      call defrr (opos,    3*nbas)
C
C      print 222, (doswtb(3,3,ic,16,1,1),ic=1,3)
C  222 format('a',3f14.7,$)
C
C      do  ic = 1, nclass
C        call psymr0(-2,ic,nbas,w(oipc),w(opos0),w(opos),w(oipa),nrclas)
C        ib0 = ival(w(oipa),1)
C        is = igetss('site spec',ib0,ssite)
CC       call upack('spec lmxa nr',sspec,is,lmxa,nr,0,0)
C        nlml = nl*nl
C        if (ipr .ge. 40) write(stdo,800) ic,nrclas,nlml
C  800   format(/' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)
C
CC   ... Make the projectors; make to at least to l=1 for forces
C        nlmx = nlml
C        call defrr (osym,   nlmx*nlmx*nrclas)
C        call symprj(nrclas,nlmx,ngrp,nbas,w(oistab),w(og),w(oag),plat,
C     .    qlat,w(opos),w(osym))
C
CC   ... Apply the projectors to doswtb; add to doswt
C        n1 = 1
C        n2 = nev*nspc*nspc
C        call defrr(owk,  nlml*n2)
C        call symr1(nrclas,w(oipa),n1,nlml,n2,nbas,nlmx,nbas,w(osym),
C     .    doswtb,w(owk))
C
C        call symr2(w(oipa),ic,1,nlml,n2,nbas,nclass,doswtb,doswt)
C        call rlse(owk)
C
CC   ... Apply the projectors to qnub, (rotate as nlml=1); add to qnu
C        n1 = 3*nl*nsp
C        n2 = 1
C        call defrr(owk,  n1*1)
C        call symr1(nrclas,w(oipa),n1,1,n2,1,nlmx,nbas,w(osym),qnub,
C     .    w(owk))
C        call rlse(owk)
C        call symr2(w(oipa),ic,n1,1,1,nbas,nclass,qnub,qnu)
C
C        call rlse (osym)
C      enddo
C      call rlse (oips)
C
C      call tcx('asadq4')
C
C      end
C
C      subroutine symr1(nrclas,ipa,n1,nlml,n2,ns2,nlmx,nbas,sym,a,wk)
CC- Symmetrize array for one class of sites
CC  array has form a(n1,nlml,nbas,n2)
C      implicit none
CC ... Passed parameters
C      integer nrclas,n2,ns2,nbas
C      integer ipa(nrclas),nlmx,n1,nlml
C      double precision sym(nlmx,nlmx,nrclas),wk(n1,nlml,n2)
C      double precision a(n1,nlml,nbas,n2)
CC ... Local parameters
C      integer stdo,ia,ib,iprint,lgunit,nn,i2
C      double precision wgt
C
CC     stdo = lgunit(1)
C
CC ... Accumulate symmetrized true density on first site
C      call dpzero(wk, n1*nlml*n2)
C      do  ia = 1, nrclas
C        ib = ipa(ia)
C        call pxsmr1(1d0,n1,nlml,n2,ns2,sym(1,1,ia),a(1,1,ib,1),wk,nn)
C      enddo
C
CC ... Copy to all sites in class
C      wgt = nrclas
C      do  ia = 1, nrclas
C      ib = ipa(ia)
C      do  i2 = 1, n2
C        call dpzero(a(1,1,ib,i2), n1*nlml)
C      enddo
C      call pysmr1(wgt,n1,nlml,n2,ns2,sym(1,1,ia),wk,a(1,1,ib,1),nn)
C      enddo
C
CC      if (iprint() .ge. 40) write(stdo,100) nn,nlml*nlml
CC  100 format(' psymr: did',i5,'  of',i5)
C
C      end
C
C      subroutine symr2(ib,ic,n1,nlml,n2,nbas,nclass,as,ac)
CC- Add array ordered by site for one class to array ordered by class
CC  array has form a(n1,nlml,nbas,n2)
C      implicit none
CC ... Passed parameters
C      integer ib,ic,n1,nlml,n2,nbas,nclass
C      double precision as(n1,nlml,nbas,n2),ac(n1,nlml,nclass,n2)
CC ... Local parameters
C      integer i2
C
C      do  i2 = 1, n2
C        call dpadd(ac(1,1,ic,i2),as(1,1,ib,i2),1,n1*nlml,1d0)
C      enddo
C
C      end

