      subroutine smshft(ssite,sspec,slat,sctrl,sham,orhoat,smrho)
C- Estimate the smooth density for a shift in atomic positions.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos pos0
Ci     Duplicate: spec spec
Ci     Stored:    pos
Ci     Passed to: pvsms1 rhgcmp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z p pz lmxa lmxl a nr rmt nxi exi chfa rsmfa
Ci     Duplicate: z p pz lmxa
Ci     Stored:    *
Ci     Passed to: pvsms1 gtpcor rhgcmp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat qlat nabc ng ogv okv vol
Ci     Duplicate: alat vol
Ci     Stored:    *
Ci     Passed to: pvsms1 symsmr rhgcmp
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: lfrce
Ci     Stored:    *
Ci     Passed to: *
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: elind
Ci     Stored:    *
Ci     Passed to: *
Ci   orhoat:vector of offsets containing site density
Co Outputs
Co   smrho :a perturbation is added to smrho, depending on job
Cr Remarks
Cr   job describes which ansatz for charge shift is used for correction
Cr     <=0  do not calculate correction to force
Cr       1  shift in free-atom density
Cr       2  shift in core+nuclear density
Cr     +10  to screen the rigid shift by the Lindhard function
Cr   (job taken from ctrl->lfrce)
Cu Updates
Cu   17 Sep 01 Adapted for local orbitals.  Altered argument list
Cu    3 Jul 00 Adapted from nfp smshft.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer orhoat(3,1)
      double precision ssite(1),sspec(1),slat(1),sctrl(1),sham(1)
      double complex smrho(*)
C ... Local parameters
      integer i,ib,igets,igetss,iprint,is,k1,k2,k3,kcor,lcor,lmxa,n0,n1,
     .  n2,n3,nbas,ng,ngabc(3),nglob,nsp
      integer ocgr,ocgs,ocwk,ogv,okv
      integer kmax,job
      parameter (n0=10)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,dgets,elind,pi,qc,qsc,qcor(2),plat(3,3),
     .  qlat(3,3),qv,qval,tpiba,vol,z,pnu(n0,2),pnz(n0,2)
C ... Heap
      integer w(1)
      common /w/ w

C     stdo = lgunit(1)
      job = igets('ctrl lfrce',sctrl)
      if (job .le. 0) return

      call tcn('smshft')

C --- Setup and printout ---
      nsp  = nglob('nsp')
      nbas = nglob('nbas')
      call upack('lat alat plat qlat',slat,alat,plat,qlat,0,0)
      call upack('lat nabc ng ogv okv vol',slat,ngabc,ng,ogv,okv,vol)
      call fftz30(n1,n2,n3,k1,k2,k3)
C     call zprm3('input smrho',0,smrho,k1,k2,k3)

C ... Hold on to original smrho (cgr)
      call defcc (ocgr,  ng*nsp)
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
      call gvgetf(ng,nsp,w(okv),k1,k2,k3,smrho,w(ocgr))

C --- Shift in unscreened density at the two positions ---
      call defcc (ocgs, ng*nsp)
      call defcc (ocwk, ng)
      kmax = 0
      call pvsms1(ssite,sspec,slat,nbas,nsp,kmax,
     .  ng,w(ogv),orhoat,w(ocwk),w(ocgs),job)
      call rlse(ocwk)

C ... Debugging: print unscreened shift in pseudo core density
C      call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocgs),smrho)
C      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C      call zprm3('unscreened local density',0,smrho,k1,k2,k3)

C --- Screened shift ---
      if (job .gt. 10) then
C       Compute elind if not given
        qval = 0d0
        do  12  ib = 1, nbas
          is = igetss('site spec',ib,ssite)
          call upack('spec z p pz lmxa',sspec,is,z,pnu,pnz,lmxa)
          if (lmxa .eq. -1) goto 12
          call gtpcor(sspec,is,kcor,lcor,qcor)
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
          qval = qval+qv
   12   continue
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
        elind = dgets('ham elind',sham)
        if (elind .lt. 0d0) elind = -(3*pi**2*qval/vol)**.66666d0*elind
        if (nsp.eq.2) call dsumdf(ng*2,1d0,w(ocgs),0,1,w(ocgs),ng*2,1)
C        call lindxx(122,n1,n2,n3,k1,k2,k3,ng,w(okv),w(ocgs),w(ogv),
C     .    tpiba,elind,w,w,w,w)
        call lindsc(2,ng,w(ogv),tpiba,elind,w(ocgs))
C        call dscal(2*ng,.001d0,w(ocgs),1)

        if (nsp.eq.2) call dsumdf(ng*2,.5d0,w(ocgs),0,1,w(ocgs),ng*2,1)
      endif

C ... Debugging: show delta smrho
      call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocgs),smrho)
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C      call zprm3('screened delta smrho',0,smrho,k1,k2,k3)
C      print *, 'returning with delta smrho'
C      return

C --- Add shift to smrho, ensuring no shift in <rho> ---
      do  i = 1, nsp
        call dvset(w(ocgs),1+ng*(i-1),1+ng*(i-1),0d0)
      enddo
      call dpadd(w(ocgr),w(ocgs),1,ng*2*nsp,1d0)
      call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocgr),smrho)
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)

C --- Symmetrize the shifted density ---
      call symsmr(slat,nsp,k1,k2,k3,smrho)

      if (iprint() .gt. 100)
     .  call zprm3('shifted smrho',0,smrho,k1,k2,k3*nsp)

      call rlse(ocgr)
      call tcx('smshft')

      end

      subroutine pvsms1(ssite,sspec,slat,nbas,nsp,kmax,ng,
     .  gv,orhoat,cwk,cg,job)
C- Shift in smoothed density according to job.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos pos0
Ci     Stored:    pos
Ci     Passed to: rhgcmp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl z p pz lmxa a nr rmt nxi exi chfa rsmfa
Ci     Stored:    name
Ci     Passed to: spacks gtpcor rhgcmp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat vol
Ci     Stored:    *
Ci     Passed to: rhgcmp
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   cy    :Normalization constants for spherical harmonics
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   orhoat:vector of offsets containing site density
Ci   job :describes which ansatz for charge shift is used for correction
Ci         <=0  do not calculate correction to force
Ci           1  shift in free-atom density
Ci           2  shift in core+nuclear density
Ci         +10  to screen the rigid shift by the Lindhard function
Co Outputs
Co  cg   coefficients to FT of shifted density
Cl Local variables
Cl  qloc   :difference in true and smoothed local charge.
Cl         :If the charge is assembled from overlapping atom-centered
Cl         :densities, qloc is the difference between the smoothed
Cl         :and head densities.
Cr Remarks
Cr   Shift of the "free atom densities."
Cr   The table below shows the densities and corresponding charges, and
Cr   parameters that hold their representations (true and smoothed
Cr   approximate forms):
Cr      density     charge    reps'n     smooth -reps'n
Cr      rho(smH)     qfat    cofh,ceh    already smooth
Cr      rho1-rho2    qloc     rhoat      qkl
Cr      rhoc         qc
Cr   This routine constructs the following difference
Cr   rhat(final) - rhat(initial)  positions, in the smooth reps'n, where
Cr      rhat = rho(smH) + qg * g(r)
Cr   where
Cr       qg = qval+qsc-qfat-qloc
Cr
Cr   In the special case rho is assembled from a superposition of
Cr   free-atom densities, and rho1-rho2 = rhoval(free-atm)-rho(smH)
Cr   (see ovlcor.f).  Thus in this case:
Cr      rho(free-atom) = rho(smH) + rho1-rho2 + rhoc
Cr   with the corresponding integrated charges
Cr       qval=z-qc     = qfat     + qloc      - qc
Cr   Thus in this special case qg=0: the only shift comes from rho(smH).
Cr   Because the local density (which contains the remaining part of
Cr   the free-atom density) will automatically be shifted, it follows
Cr   that the shifted smooth density will correspond to the
Cr   smooth sum-of-FA densities constructed at the shifted positions.
Cr
Cr   In the general case, qg is not zero.  By shifting the a gaussian
Cr   along with the sm-Hankels, the integrated total density of charge
Cr   shifted (local density + mesh density) is neutral.
Cr
Cr   Improvements: if the tail density were also shifted inside each
Cr   augmentation sphere, the total density would correspond exactly
Cr   to the sum-of-FA densities at the shifted positions, when the
Cr   starting density is also a sum-of-FA densities.
Cr
Cr   Shift of the "core + valence densities."
Cr
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,ng,job,orhoat(3,1)
      double precision ssite(1),sspec(1),slat(1),gv(ng,3)
      double complex cg(ng,nsp),cwk(ng)
C Local variables
      integer ib,i,is,iv0,kmax,lmxl,igetss,lmxa,nr,nrmx,n0,
     .  nlml,nxi,ie,ixi,ig,ipr,iprint,kcor,lcor,stdo,lgunit
      parameter (nrmx=1501, n0=10)
      double precision a,aa,alat,df(0:20),e,
     .  exi(n0),gam,hfc(n0,2),pi,pnew(3),pnu(n0,2),pnz(n0,2),pold(3),
     .  pp,qall,qc,qcor(2),qsc,qfat,qg,qloc,qval,rmt,
     .  rsmfa,rwgt(nrmx),scalp,sum,tpiba,v(3),v2,vol,volsp,y0,z
      character*35 strn,spid*8
      double complex phase
C ... Heap
      integer w(1)
      common /w/ w

C ... Setup
      call tcn('pvsms1')
      call stdfac(20,df)
      call upack('lat alat vol',slat,alat,vol,0,0,0)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      tpiba = 2*pi/alat
      ipr = iprint()
      stdo = lgunit(1)
      volsp = vol*nsp

C --- For each site, accumulate shift in density ---
      call dpzero(cg,2*ng*nsp)
      iv0 = 0
      strn = 'free atom densities'
      if (job .eq. 11) strn = 'screened free atom densities'
      if (job .eq. 12) strn = 'screened core+multipole densities'
      if (ipr .ge. 30) write(stdo,339) strn
  339 format(/' smshft:  add shifted ',a/'   site',16x,'old pos',22x,
     .  'new pos',14x,'shift')

      do  10  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        call spacks(0,'spec name',sspec,spid,is,is)
        lmxl = igetss('spec lmxl',is,sspec)
        nlml = (lmxl+1)**2
        if (lmxl .eq. -1) goto 10

        call upack('site spec pos pos0',ssite,ib,is,pnew,pold,0)
        pp = alat*dsqrt((pnew(1)-pold(1))**2 + (pnew(2)-pold(2))**2
     .                + (pnew(3)-pold(3))**2)
        if (ipr .ge. 30) write(stdo,340) ib,spid,pold,pnew,pp/alat
  340   format(i4,':',a,f8.5,2f9.5,2x,3f9.5,2x,f9.6)
C       Skip this site if shift is negligible
        if (pp .le. 1d-6) goto 18

C   --- Shift in mesh density, job 1 ---
        if (mod(job,10) .eq. 1) then
          call upack('spec z p pz lmxa',sspec,is,z,pnu,pnz,lmxa)
          call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
          call upack('spec nxi exi chfa rsmfa',sspec,is,nxi,exi,hfc,
     .      rsmfa)
          call gtpcor(sspec,is,kcor,lcor,qcor)
          if (nr .gt. nrmx) call rx('dfrce: nr gt nrmx')
          call radwgt(rmt,a,nr,rwgt)
          call radsum(nr,nr,nlml,nsp,rwgt,w(orhoat(1,ib)),qloc)
          call radsum(nr,nr,nlml,nsp,rwgt,w(orhoat(2,ib)),sum)
          qloc = (qloc-sum)/y0
          qfat = 0d0
          do  12  i  = 1, nsp
          do  12  ie = 1, nxi
            gam  = 0.25d0*rsmfa**2
            qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
            qfat = qfat + hfc(ie,i)*qall
   12     continue
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
C         Excess sphere charge.  See Remarks above.
          qg = qval+qsc-qfat-qloc

C     ... Shift in smoothed free atom density
          do  14    i = 1, nsp
          do  14  ixi = 1, nxi
          e = exi(ixi)
          do  15 ig = 1, ng
            v(1) = gv(ig,1)*tpiba
            v(2) = gv(ig,2)*tpiba
            v(3) = gv(ig,3)*tpiba
            v2 = v(1)**2+v(2)**2+v(3)**2
            aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
            scalp = -alat*(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
            phase = dcmplx(dcos(scalp),dsin(scalp))
            scalp = -alat*(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
            phase = phase - dcmplx(dcos(scalp),dsin(scalp))
            cg(ig,i) = cg(ig,i) + hfc(ixi,i)*aa*phase*y0/vol
   15     continue
   14     continue

C     ... Add gaussian to conserve local charge; see Remarks
          do  16   i = 1, nsp
          do  16  ig = 1, ng
            v(1) = gv(ig,1)*tpiba
            v(2) = gv(ig,2)*tpiba
            v(3) = gv(ig,3)*tpiba
            v2 = v(1)**2+v(2)**2+v(3)**2
            scalp = -alat*(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
            phase = dcmplx(dcos(scalp),dsin(scalp))
            scalp = -alat*(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
            phase = phase - dcmplx(dcos(scalp),dsin(scalp))
            cg(ig,i) = cg(ig,i) + qg*phase*dexp(-gam*v2)/volsp
   16     continue

C   --- Shift in mesh density, job 12 ---
        elseif (job .eq. 12) then

C     ... Core + valence at old position
          call dpzero(cwk,ng*2)
          call pack2('site pos',ssite,ib,pold)
          call rhgcmp(131,ib,ib,ssite,sspec,slat,orhoat,kmax,ng,cwk)
          call dscal(ng*2,-1d0,cwk,1)
C     ... Core + valence at new position
          call pack2('site pos',ssite,ib,pnew)
          call rhgcmp(131,ib,ib,ssite,sspec,slat,orhoat,kmax,ng,cwk)

C     ... Add to cg
          do  i = 1, nsp
            call daxpy(ng*2,1d0/nsp,cwk,1,cg(1,i),1)
          enddo

        else
          call rxi('smshft: bad job:',job)
        endif

   18   iv0 = iv0+nlml
   10 continue

      call tcx('pvsms1')
      end

      subroutine pvsms2(ssite,sspec,rotm,nbas,nsp,orhoat)
C- Rotate local densities by specified rotation
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: nr lmxl
Ci     Stored:    name
Ci     Passed to: spacks
Ci   rotm  :3x3 cartesian rotation matrix
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   orhoat:vector of offsets containing site density
Co Outputs
Co   orhoat:On output the different m-channels of rhoat(1) and rhoat(2)
Co         :are mixed by the rotation
Cl Local variables
Cr Remarks
Cr   For a rotation matrix R, The density is stored in the 1-center form
Cr      rho_l(r) YL(rhat)
Cr   Given a rotation matrix R, this it transforms as
Cr      rho_l(r) YL(R rhat) = rho_l(r) rYL(rhat)
Cr   where rYL is made by ylmrtg
Cr
Cb Bugs
Cb   No ability is supplied when the Yl are true instead of real
Cb   spherical harmonics
Cu Updates
Cu   21 Dec 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,orhoat(3,1)
      double precision ssite(1),sspec(1)
      double precision rotm(3,3)
C Local variables
      integer ib,i,j,is,lmxl,igetss,nr,nlml,ipr,stdo,nlx,nl2,nglob
      parameter (nlx=9, nl2=nlx*nlx)
      double precision rYL(nl2,nl2)
      character spid*8
C ... Heap
      integer w(1)
      common /w/ w

C ... Setup
      call getpr(ipr)
      stdo = nglob('stdo')

C ... Rotation matrix for real spherical harmonics
C     call prmx('pvsms2 rotm',rotm,3,3,3)
      call ylmrtg(nl2,rotm,rYL)
C     call prmx('rYL',rYL,nl2,nl2,nl2)

C --- For each site and l, rotate the m-components ---
      if (ipr .ge. 20) then
        call info0(20,0,0,' Rotate local densities using R=')
        write (stdo,350) ((rotm(i,j),j=1,3),i=1,3)
  350   format(3f11.6)
      endif

      do  10  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        call spacks(0,'spec name',sspec,spid,is,is)
        call upack('spec nr lmxl',sspec,is,nr,lmxl,0,0)
        if (lmxl .eq. -1) goto 10
        nlml = (lmxl+1)**2
        if (nlml .gt. nl2) call rx('increase nl2 in pvsms2')

        call pvsms3(nr,nr,nlml,nsp,rYL,nl2,w(orhoat(1,ib)))
        call pvsms3(nr,nr,nlml,nsp,rYL,nl2,w(orhoat(2,ib)))

   10 continue

      end

      subroutine pvsms3(nrx,nr,nlml,nsp,rYL,nl2,rho)
C- Rotation of an l-dependent density
C ----------------------------------------------------------------------
Ci Inputs
Ci   nrx   :leading dimension of rho
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nsp   :2 for spin-polarized case, otherwise 1
Cl   rYL   :rotation matrix that rotates Y_lm
Ci   nl2   :leading dimension of rYL
Co Outputs
Co   rho   :On output the different m-channels of rho are
Co         :mixed by rYL
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   21 Dec 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nrx,nr,nlml,nsp,nl2
      double precision rho(nrx,nlml,nsp)
      double precision rYL(nl2,nl2)
C ... Local parameters
      integer isp
C     integer l,lmax,ll,nlmi,offri
      double precision rwk(nrx,nlml)

C     lmax = ll(nlml)

C     call prmx('starting rho',rho,nrx,nr,nlml*nsp)

      if (nlml .eq. 0) return
      do  isp = 1, nsp

        call dgemm('N','T',nr,nlml,nlml,1d0,rho(1,1,isp),nrx,
     .    rYL,nl2,0d0,rwk,nrx)
        call dcopy(nrx*nlml,rwk,1,rho(1,1,isp),1)

C        faster if done l-by-l
C        do  l = 0, lmax
C
C          nlmi = 2*l + 1
C          offri = l**2
C          print *, l, nlmi,offri
C
C          call dgemm
C        enddo
      enddo

C     call prmx('ending rho',rho,nrx,nr,nlml*nsp)

      end

