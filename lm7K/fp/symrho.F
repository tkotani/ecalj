      subroutine symrho(ssite,sspec,slat,lf,smrho,orhoat,qbyl,hbyl,f)
C- Symmetrize charge density and related quantities
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    class pos
Ci     Passed to: spackv
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl lmxa nr
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy plat qlat oistab nsgrp osymgr oag
Ci     Stored:    *
Ci     Passed to: *
Ci   lf    :>0 symmetrize forces
Ci Inputs/Outputs
Cio  smrho :smooth density
Cio        :Symmetrized on output
Cio  orhoat:vector of offsets containing site density
Cio        :Symmetrized on output
Cio  qbyl  :site- and l-decomposed charges
Cio        :Symmetrized on output
Cio  hbyl  :site- and l-decomposed one-electron energies
Cio        :Symmetrized on output
Cio  f     :forces
Cio        :Symmetrized on output
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   19 Jun 00 Packaged from nfp symrat.f and symsmr.f
C ----------------------------------------------------------------------
C     implicit none
      integer lf,orhoat(3,*)
      double precision ssite(1),sspec(1),slat(1),f(*),qbyl(*),hbyl(*)
      double complex smrho(*)
      integer nsp,nbas,ngabc(3),n1,n2,n3,k1,k2,k3,nglob,iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))

      call tcn('symrho')
      call info(30,1,0,' Symmetrize density..',0,0)
      nbas = nglob('nbas')
      nsp = nglob('nsp')
      call upack1('lat nabc',slat,ngabc)
      call fftz30(n1,n2,n3,k1,k2,k3)

      call symrat(ssite,sspec,slat,nbas,nsp,lf,orhoat,qbyl,hbyl,f)
      if (iprint() .gt. 50) call prrhat(nbas,ssite,sspec,orhoat)
      call symsmr(slat,nsp,k1,k2,k3,smrho)
      call tcx('symrho')
      end

      subroutine symrat(ssite,sspec,slat,nbas,nsp,lf,orhoat,qbyl,hbyl,f)
C- Symmetrize the atomic charge densities and the forces.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    class pos
Ci     Passed to: spackv
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl lmxa nr
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy plat qlat oistab nsgrp osymgr oag
Ci     Stored:    *
Ci     Passed to: *
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lf    :>0 symmetrize forces
Ci Inputs/Outputs
Cio  orhoat:vector of offsets containing site density
Cio        :Symmetrized on output
Cio  qbyl  :site- and l-decomposed charges
Cio        :Symmetrized on output
Cio  hbyl  :site- and l-decomposed one-electron energies
Cio        :Symmetrized on output
Cio  f     :forces
Cio        :Symmetrized on output
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,orhoat(3,nbas),n0,lf
      parameter (n0=10)
      double precision ssite(1),sspec(1),slat(1),
     .  f(3,nbas),qbyl(n0,nsp,nbas),hbyl(n0,nsp,nbas)
C ... Local parameters
      integer stdo,ib0,ic,ipr,iprint,is,lgunit,lmxa,lmxl,nclass,
     .  ngrp,nlml,nlmx,nr,nrclas,oag,ocg,ocy,og,oidxcg,oipa,oipc,oips,
     .  ojcg,opos,opos0,orho,osym,igetss,mxint,oistab,ival
      double precision plat(9),qlat(9)
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('symrat')
      stdo = lgunit(1)
      ipr = iprint()

      call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
      call upack('lat plat qlat oistab',slat,plat,qlat,oistab,0,0)
      call upack('lat nsgrp osymgr oag',slat,ngrp,og,oag,0,0)

      call defi (oips, nbas)
      call defi (oipc, nbas)
      call defrr(opos0,3*nbas)

C ... Separate sites into symmetry classes
C      call symcls(nbas,ssite,sspec,slat,ngrp,w(og),w(oag),
C     .  w(oips),w(opos0),nclass,w(oipc))
C      if (ipr .ge. 40) write(stdo,300) nclass
C  300 format(/' symrat: number of symmetry classes is',i3)

      call spackv(10,'site class',ssite,1,nbas,w(oipc))
      call spackv(10,'site pos',ssite,1,nbas,w(opos0))
      nclass = mxint(nbas,w(oipc))

C --- Start loop over classes ---
      call defi  (oipa,    nbas)
      call defrr (opos,    3*nbas)

      do  ic = 1, nclass
        call psymr0(-2,ic,nbas,w(oipc),w(opos0),w(opos),w(oipa),nrclas)
        if (nrclas .gt. 0) then
        ib0 = ival(w(oipa),1)
        is = igetss('site spec',ib0,ssite)
        call upack('spec lmxl lmxa nr',sspec,is,lmxl,lmxa,nr,0)
        nlml = (lmxl+1)**2
        if (ipr .ge. 40) write(stdo,800) ic,nrclas,nlml
  800   format(/' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)

C   ... Make the projectors; make to at least to l=1 for forces
        nlmx = max0(nlml,4)
        call defrr (osym,   nlmx*nlmx*nrclas)
        call symprj(nrclas,nlmx,ngrp,nbas,w(oistab),w(og),w(oag),plat,
     .    qlat,w(opos),w(osym))

C   ... Apply the projectors to rhoat
        if (lmxl .gt. -1) then
        call defrr(orho,  nr*nlml*nsp)
        call psymr1(nrclas,w(oipa),nr,nlml,nsp,nlmx,w(osym),w(orho),
     .    orhoat,1)
        call psymr1(nrclas,w(oipa),nr,nlml,nsp,nlmx,w(osym),w(orho),
     .    orhoat,2)

C   ... Symmetrize site charges and eigval sum
        call psymrq(nrclas,nsp,w(oipa),lmxa,qbyl,hbyl)
        endif

C   ... Symmetrize the forces
        if (lf .ne. 0) call psymrf(nrclas,w(oipa),nlmx,w(osym),f)

        call rlse (osym)
        endif
      enddo

      call rlse (oips)

      call tcx('symrat')

      end

      subroutine psymrf(nrclas,ipa,nlmx,s,f)
C- Symmetrize forces
C     implicit none
C ... Passed parameters
      integer nrclas,nlmx,ipa(nrclas)
      double precision s(nlmx,nlmx,nrclas),f(3,1)
C ... Local parameters
      integer ia,ib
      double precision x(3)

      x(1) = 0d0
      x(2) = 0d0
      x(3) = 0d0
      do  ia = 1, nrclas
        ib = ipa(ia)
        x(1)= x(1)+s(4,4,ia)*f(1,ib)+s(4,2,ia)*f(2,ib)+s(4,3,ia)*f(3,ib)
        x(2)= x(2)+s(2,4,ia)*f(1,ib)+s(2,2,ia)*f(2,ib)+s(2,3,ia)*f(3,ib)
        x(3)= x(3)+s(3,4,ia)*f(1,ib)+s(3,2,ia)*f(2,ib)+s(3,3,ia)*f(3,ib)
      enddo
      do  ia = 1, nrclas
        ib = ipa(ia)
        f(1,ib) = (s(4,4,ia)*x(1)+s(2,4,ia)*x(2)+s(3,4,ia)*x(3))*nrclas
        f(2,ib) = (s(4,2,ia)*x(1)+s(2,2,ia)*x(2)+s(3,2,ia)*x(3))*nrclas
        f(3,ib) = (s(4,3,ia)*x(1)+s(2,3,ia)*x(2)+s(3,3,ia)*x(3))*nrclas
      enddo
      end

      subroutine psymrq(nrclas,nsp,ipa,lmxa,qbyl,hbyl)
C- Symmetrize l-decomposed site charges and eval sums
C     implicit none
C ... Passed parameters
      integer nrclas,nsp,lmxa,ipa(nrclas),n0
      parameter (n0=10)
      double precision qbyl(n0,nsp,1),hbyl(n0,nsp,1)
C ... Local parameters
      integer stdo,ia,ib,iprint,l,lgunit,isp
      double precision qsum(n0,2),hsum(n0,2),fac

      stdo = lgunit(1)
      call dpzero(qsum,2*n0)
      call dpzero(hsum,2*n0)
      fac = 1d0/nrclas
      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
        do  l = 0, lmxa
          qsum(l+1,isp) = qsum(l+1,isp) + qbyl(l+1,isp,ib)*fac
          hsum(l+1,isp) = hsum(l+1,isp) + hbyl(l+1,isp,ib)*fac
        enddo
        enddo
      enddo

      if (iprint() .ge. 40) then
        write(stdo,770) (ipa(ia),ia = 1,nrclas)
  770   format(' symmetrized qbyl,hbyl for class containing ib=',20i3)
        if (nsp .eq. 1) write(stdo,780)
     .    (l,qsum(l+1,1),hsum(l+1,1), l=0,lmxa)
        if (nsp .eq. 2) write(stdo,781)
     .    (l,(qsum(l+1,isp),hsum(l+1,isp),isp=1,nsp), l=0,lmxa)
  780   format(i5,2f12.6)
  781   format(i5,2f12.6,'   spin 2',2f12.6)
      endif

      do  ia = 1, nrclas
        ib = ipa(ia)
        do  isp = 1, nsp
        do  l = 0, lmxa
          qbyl(l+1,isp,ib) = qsum(l+1,isp)
          hbyl(l+1,isp,ib) = hsum(l+1,isp)
        enddo
        enddo
      enddo

      end

      subroutine psymr1(nrclas,ipa,nr,nlml,nsp,nlmx,sym,rho,orhoat,icmp)
C- Symmetrize density for one class of atoms
C     implicit none
C ... Passed parameters
      integer nrclas,nsp
      integer ipa(nrclas),orhoat(3,1),nlmx,nr,nlml,icmp
      double precision sym(nlmx,nlmx,nrclas),rho(nr,nlml,nsp)
C ... Local parameters
      integer stdo,ia,ib,iprint,lgunit,nn
      double precision wgt
C ... Heap
      integer w(1)
      common /w/ w

C ... Accumulate symmetrized true density on first site
      stdo = lgunit(1)
      call dpzero(rho, nr*nlml*nsp)
      do  ia = 1, nrclas
        ib = ipa(ia)
        call pxsmr1(1d0,nr,nlml,nsp,sym(1,1,ia),w(orhoat(icmp,ib)),rho,
     .    nn)
      enddo

C ... Copy to all sites in class
      wgt = nrclas
      do  ia = 1, nrclas
      ib = ipa(ia)
      call dpzero(w(orhoat(icmp,ib)), nr*nlml*nsp)
      call pysmr1(wgt,nr,nlml,nsp,sym(1,1,ia),rho,w(orhoat(icmp,ib)),nn)
      enddo

      if (iprint() .ge. 40) write(stdo,100) nn,nlml*nlml
  100 format(' psymr: did',i5,'  of',i5)

      end

      subroutine symsmr(slat,nsp,k1,k2,k3,smrho)
C- Symmetrize the smooth charge density
C     implicit none
C ... Passed parameters
      integer nsp,k1,k2,k3
      double precision slat(1)
      double complex smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer n1,n2,n3,ng,ngrp,obgv,ocsym,ocv,ogv,oips0,okv,ngabc(3),
     .  isp
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C ... Heap
      integer w(1)
      common /w/ w

      call tcn('symsmr')
      call upack('lat nabc ng ogv okv',slat,ngabc,ng,ogv,okv,0)
      call upack('lat nsgrp oips0 obgv',slat,ngrp,oips0,obgv,0,0)
      if (ngrp .gt. 1) then
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3)
        call defcc(ocv,     ng)
        call defcc(ocsym,   ng)
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
C       call zprm3('smrho before poke',smrho,k1,k2,k3)
        do  isp = 1, nsp
          call gvgetf(ng,1,w(okv),k1,k2,k3,smrho(1,1,1,isp),w(ocv))
          call gvsym(ng,w(ogv),w(oips0),w(obgv),w(ocv),w(ocsym))
          call dpadd(w(ocsym),w(ocv),1,ng*2,-1d0)
          call gvaddf(ng,w(okv),k1,k2,k3,w(ocsym),smrho(1,1,1,isp))
C         call zprm3('smrho after poke',smrho,k1,k2,k3)
        enddo
        call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
        call rlse(ocv)

C ... Force density to be real and positive
C       call rhopos(smrho,k1,k2,k3,n1,n2,n3)
C        do  10  i23 = 1, k2*k3
C        do  10  i1  = 1, k1
C   10   smrho(i1,i23,1) = dble(smrho(i1,i23,1))

      else
        call info(30,1,1,' Smooth density not symmetrized (ngrp=1)',0,0)
      endif

      call tcx('symsmr')

      end

      subroutine rhopos(smrho,k1,k2,k3,n1,n2,n3)
C- Make smrho real and positive
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3
      double complex smrho(k1,k2,k3)
C ... Local parameters
      integer stdo,lgunit,i1,i2,i3,nneg
      double precision rmin,xx

      stdo = lgunit(1)
      nneg = 0
      rmin = 999
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1
        xx = dble(smrho(i1,i2,i3))
        rmin = min(rmin,xx)
        if (xx .lt. 0) then
          nneg = nneg+1
          xx = 1d-8
        endif
        smrho(i1,i2,i3) = xx
   10 continue

      if (nneg .gt. 0) write(stdo,333) nneg,rmin
  333 format(' rhopos (warning): mesh density negative at',i6,
     .  ' points.  min=',f13.8)

      end

      subroutine rhoqm(smrho,k1,k2,k3,n1,n2,n3,nsp,vol,qsum)
C- Return charge, magnetic moment of smooth density
C ----------------------------------------------------------------------
Ci Inputs
Ci   smrho :smooth density on uniform mesh 
Ci   k1..k3:
Ci   n1..n3:
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   vol   :cell volume
Co Outputs
Co   qsum  :qsum(1) = smrho(+) + smrho(-)
Co         :qsum(2) = smrho(+) - smrho(-) (for nsp=2 only)
Cl Local variables
Cl         :
Cr Remarks
Cr   Input smrho is assumed to be (rho1, rho2) 
Cr   If instead smrho=(rho1+rho2,rho1-rho2) => qsum(1,2) = q+amom, q-amom
Cu Updates
Cu   13 Dec 08 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k1,k2,k3,n1,n2,n3,nsp
      double complex smrho(k1,k2,k3,nsp)
      double precision vol,qsum(2) 
C ... Local parameters
      integer i,i1,i2,i3
      double precision sumi,q1,fac

      qsum(1) = 0
      qsum(2) = 0
      fac = vol/(n1*n2*n3)
      q1 = 0
      do  i = 1, nsp
        sumi = 0
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          sumi = sumi + dble(smrho(i1,i2,i3,i))
        enddo
        enddo
        enddo
        if (i .eq. 2) qsum(2) = qsum(2) + q1-sumi
        q1 = sumi
        qsum(1) = qsum(1) + sumi
      enddo
      qsum(1) = fac*qsum(1)
      qsum(2) = fac*qsum(2)
C     write(*,333) qsum
C 333 format(' rhoqm : istl charge, moment = ',2f13.7)

      end

