#define BLAS3 1
      subroutine secm2c(sctrl,ckbas,plat,nl,nsp,nbas,nclass,ipc,indxsh,
     .  qss,eula,neul,ikp,nkp,ldim,lidim,lihdim,ndim,qp,pp,sop,vmtz,wsr,
     .  avw,isp,pti,pph,nfilet,nevmx,efmax,nev,addsx,z,eb)
C- Set up and diagonalize ASA Hamiltonian, 2-center approximation.
C ----------------------------------------------------------------
Ci Inputs:
Ci   plat,nl,nsp,nbas,lmx,pp,vmtz,indxsh
Ci   eula: Euler angles for spin rotations
Ci   ldim: dimension of l - wave block of hamiltonian matrix
Ci   nevmx:-2 Do not diagonalize, but copy h into z
Co Outputs:
Co   eigenvalues and eigenvectors are returned in eb, z
Co   pph : vector of pot. par's. in alpha rep'n; pp are returned
Co   in alpha representation.
Cr Remarks
Cr   Folding down not allowed here.
Cr   z must be available as a work space of size ldim**2*2, even if
Cr   evecs are not calculated.
C ----------------------------------------------------------------
C     implicit none
      integer nl,nsp,nbas,nclass,neul,isp,nfilet,ldim,lidim,lihdim,ndim,
     .  nevmx,nev,ipc(nbas),indxsh(ndim),ikp,nkp,addsx
      double precision ckbas,plat(3,3),eb(ldim*2),pp(6,nl,nsp,0:*),vmtz,
     .  avw,qp(3),z(ldim,ldim*2),pti(ndim,nsp),eula(nbas,neul,3),wsr(*),
     .  pph(5,lihdim,nsp),efmax,qss(4),sop(0:nl-1,nsp,nsp,9,*),sctrl(1)
C heap:
      integer w(1)
      common /w/ w
C local variables
      integer i,j,ipr,ii,lgunit,nsite,i2,li,lncol,lham,ld2,ld22,n,
     .  bit,nlspc,nl2,i1mach,linv,ldimx,oiax,ontab,os,owk,owk2,owk3,
     .  odiawk,oo,occd,oalph,oadot,oa,osd,oh,osii,osil,oxsidt,ogma,idim,
     .  lbloch,isw
      logical ndnfld,ccor,lx,lss,lnc,lso,bittst,lgors,ltmp,iostr
      double precision kap2(20),qpq(3),ddot,ercc,elin
C     double precision Ze(2),Gij(2,ldim,ldim)
c     integer oxsi,ode,onewa,osz
      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) .eq. bit)

C --- Setup ---
      call upack('ctrl lham lncol elin',sctrl,lham,lncol,elin,0,0)
      call getpr(ipr)
      nl2 = nl**2
      idim = lidim - ldim
      if (idim .ne. 0) call fexit(-1,111,'  Exit -1 SECM2C: '//
     .  '%i intermediate waves sought but none allowed',idim)
      nlspc = nl * nsp * nclass

C     Noncollinear switches
      lnc = lncol .ne. 0
      lss = bittst(lncol,2)
      lso = bittst(lncol,4)

C     Possibly rotate to spherical harmonics when making Bloch sum
      lbloch = 0
      if (bittst(lham,256)) lbloch = 1000

C     Combined correction required if downfolding
      ndnfld = idim .eq. 0
      ccor = .not. ndnfld .or. lgors('ctrl lasa,4',sctrl)
C     if (.not. ccor  .and.  ipr  .ge. 30  .and.  ikp .eq. 1)
C    . print *, 'SECM2C : Combined Correction switched off'
C     if (.not. ndnfld .and. bittst(lham,128))
C    .  call rx('secmat: no downfolding in gamma rep')

C     Diagonalize by inverse iteration, or not
      linv = 0
      if (nevmx .gt. 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1

C ... Sanity checks
      if (.not. ndnfld .and. bittst(lham,128))
     .  call rx('SECMAT: no downfolding in gamma rep')
      if (lnc .and. .not. ndnfld)
     .  call rx('noncollinear magnetism not implemented with dnfolding')
      call isanrg(addsx,0,0,'secm2c:','addsx',.true.)
      call isanrg(isw(lso),0,0,'secm2c:','spin-orbit',.true.)


C     Some dimensioning parameters and memory allocation
      ldimx = ldim
      ld2 = ldim**2
      if (lnc) ldimx = 2*ldim
      ld22 = 2*ld2
      if (ndnfld) li = 1
      if (ndnfld) i2 = 1
      call zinit(z,ldimx**2)
      call defdr(owk,2*lihdim)
      call defdc(osii,i2)
      call defdc(osil,li)
      call defdc(oa,li)
      call defdc(oh,-ldimx**2)
      call dcopy(3,qp,1,qpq,1)

C --- Get screened sdot from disc and Bloch-transform ---
      osd = 1
      if (ccor) then
        if (lss) call defdc(osd,ld22)
        if (lnc .and. .not. lss) call defdc(osd,ld2)
        ltmp = iostr(8,'SDOT',nl,nbas,1,kap2,0,ckbas,-1,nsite,
     .    oadot,oiax,ontab,os)
        nsite = w(ontab+nbas)
        if (lss) then
          qpq(1) = qp(1) + qss(1)/2
          qpq(2) = qp(2) + qss(2)/2
          qpq(3) = qp(3) + qss(3)/2
          call bloch(0,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),w(os),
     .      nl2,1,1,ldim,ldim,0,ldim,0,ldim,0,w(osd),w,w)
          if (ipr .ge. 110)
     .    call yprm('Sd(2)',12,w(osd),ldim*ldim,ldim,ldim,ldim)
          call pvsec2(ld22,w(osd))
          qpq(1) = qp(1) - qss(1)/2
          qpq(2) = qp(2) - qss(2)/2
          qpq(3) = qp(3) - qss(3)/2
        endif
        if (lss .or. lnc) then
          call bloch(0,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),w(os),
     .      nl2,1,1,ldim,ldim,0,ldim,0,ldim,0,w(osd),w,w)
          if (ipr .ge. 110)
     .    call yprm('Sd(1)',12,w(osd),ldim,ldim,ldim,ldim)
        else
          call bloch(0,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),w(os),
     .      nl2,1,1,ldim,ldim,0,ldim,0,ldim,0,z,w,w)
          call dscal(ld22,-avw**2,z,1)
        endif
        call rlse(oiax)
        call defdr(occd,3*lihdim)
      else
        occd = 1
      endif

C --- Get screened strux from disc and Bloch-transform them ---
      ltmp = iostr(8,'STR',nl,nbas,1,kap2,0,ckbas,-1,nsite,oalph,
     .  oiax,ontab,os)
      nsite = w(ontab+nbas)
      if (lss) then
        qpq(1) = qp(1) + qss(1)/2
        qpq(2) = qp(2) + qss(2)/2
        qpq(3) = qp(3) + qss(3)/2
        call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),
     .    w(os),nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z(1+ld22,1),w,w)
        qpq(1) = qp(1) - qss(1)/2
        qpq(2) = qp(2) - qss(2)/2
        qpq(3) = qp(3) - qss(3)/2
        if (ipr .ge. 110)
     .    call yprm('Sll(2)',12,z(1+ld22,1),ldim*ldim,ldim,ldim,ldim)
      endif
      call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),w(os),
     .  nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,w(oh),w,w)
      call rlse(oiax)

C --- Transform to gamma representation ---
      if (bittst(lham,128)) then
C   ... s,sdot are correctly rotated but ccor should be done in pert th
        if (ccor) call rx('secm2c: ccor not working for Gamma repsn')
        call defdr(ogma,ldim)
C   ... Put pp in gamma; make vector of gamma's
        call defdr(oo,nlspc)
        i = -3
        if (isp .eq. 2) i = -4
        if (lnc) i = -5
        call pptrns(i,nl,ipc,nclass,nsp,w(ogma),nbas,pp,w(oo))
        call rlse(oo)
        call makpph(nl,nsp,nbas,lihdim,ipc,indxsh,pp,pph)
C       call prm('gamma',0,w(ogma),lihdim,lihdim,1)
C       call prm('alpha',0,w(oalph),lihdim,lihdim,1)
C   ... put gamma-alpha in gam, gamma in alpha
        call daxpy(ldim,-1d0,w(oalph),1,w(ogma),1)
        call daxpy(ldim, 1d0,w(ogma),1,w(oalph),1)
C   ... Transform sll = w(oh) and sdot into gamma
        call defdr(owk3,ld22)
        call defdr(oxsidt,-ldim)
C   ... Make adot so that 3C contribution to cc is exactly zero
*       call prmx('gamma-alpha',w(ogma),ldim,ldim,1)
*       call prmx('adot',w(oadot),ldim,ldim,1)
        if (lss) then
C     ... Use copies of gamma-alpha, xsidot since mksbet overwrites them
          call dpcopy(w(ogma),w(owk),1,ldim,1d0)
          call defdr(owk2,ldim)
          call dpcopy(w(oxsidt),w(owk2),1,ldim,1d0)
C     ... swap 1st and 2nd B.T. sd, to avoid more workspace
          if (ccor) call pvsec2(ld22,w(osd))
          call mksbet(ccor,ldim,w(owk),w(owk2),z(1+ld22,1),w(osd),
     .      w(owk3))
          if (ipr .ge. 110) then
            if (ccor) call yprm('Sd(2) (gam)',12,w(osd),ldim*ldim,ldim,
     .        ldim,ldim)
            call yprm('Sll(2)(gam)',12,z(1+ld22,1),ldim*ldim,ldim,ldim,
     .        ldim)
          endif
          if (ccor) call pvsec2(ld22,w(osd))
        endif
        call mksbet(ccor,ldim,w(ogma),w(oxsidt),w(oh),w(osd),w(owk3))
        if (ipr .ge. 110) then
          if (ccor) call yprm('Sd (gam)',12,w(osd),ldim*ldim,ldim,ldim,
     .      ldim)
          call yprm('Sll (gamma)',12,w(oh),ldim*ldim,ldim,ldim,ldim)
        endif
        call rlse(owk3)
C --- Or else transform potential parameters to alpha rep'n ---
      else
        call defdr(oo,nlspc)
        call pptrns(0,nl,ipc,nclass,nsp,w(oalph),1,pp,w(oo))
        call rlse(oo)
      endif

C ... Copy h into z for noncoll case
      if (lnc .or. lss) then
        call dpcopy(w(oh),z,1,ld22,1d0)
      endif

C ... Scale sdot.  All scalings should be unified
      if (ccor .and. osd .ne. 1)
     .  call dscal((ldimx/ldim)*ld22,-avw**2,w(osd),1)

C --- Diagonal matrices for combined correction---
      if (ccor) then
        call makdia(nl,nbas,lihdim,indxsh,lidim,ipc,wsr,avw,w(oalph),
     .    w(oadot),w(occd))
        call dscal(3*lihdim,avw**2,w(occd),1)
        ercc = (ddot(3*ldim,w(occd),1,w(occd),1) -
     .    ddot(2*ldim,w(occd),1,w(occd),1))/ldim
        if (ercc.gt.1d-10 .and. ipr.ge.10 .and. ikp.eq.1)
     .    call awrit1(' SECM2C (warning): 3-C CC should be zero '//
     .    'but has rms %1,5;5d',' ',80,lgunit(1),dsqrt(ercc))
      endif

C --- Non-collinear two-center hamiltonian ---
      if (nsp .eq. 2 .and. (lss .or. lnc)) then
        call defdr(owk2,ldim*2)
        call hml2nc(nbas,nl,indxsh,qss,eula,neul,pph,
     .    ccor,lss,lnc,w(occd),w(owk2),vmtz,elin,ldim,z,w(osd),w(oh))
        call rlse(owk2)
C --- Collinear two-center hamiltonian ---
      else
C    ... ASA 2-center hamiltonian + ccor + terms of order o*enu ---
        if (.true. .or. ccor) then
          call hmlt2c(ccor,w(occd),vmtz,elin,ldim,lihdim,pph(1,1,isp),
     .      w(oh),z,w(owk))
        else
C    ... Simple ASA 2-center hamiltonian; uses no work space
         call makdsd(0,ldim,ldim,ldim,ldim,0,0,pph(1,1,isp),w(oh),w(oh))
         call daxpy(ldim,1d0,pph(2,1,isp),5,w(oh),ldim+1)
        endif
        call rlse(oalph)
      endif
      if (ipr .ge. 110)
     .  call yprm('H',12,w(oh),ldimx*ldimx,ldimx,ldimx,ldimx)

C --- Return with hamiltonian in z if nevmx is -2 ---
      if (nevmx .eq. -2) then
        call dcopy(ldimx**2*2,w(oh),1,z,1)
        return
      endif

C --- Eigenvalues and eigenvectors of 2C Hamiltonian ---
#if BLAS3
      lx = .true.
#else
      lx = .false.
#endif
      if (linv .ne. 0) then
        call defdr(odiawk,ldimx*11)
      else
C ...   Use 5*ldim for parallel implementations ...
        call defdr(odiawk,ldimx*5)
      endif
C      call prmx('h',w(oh),ldim,ldim,ldim)
      call diagno(ldimx,w(oh),w(1),w(odiawk),lx,
     .  .false.,linv,nevmx,efmax,nev,z,eb)

C     call prmx('evl',eb,ldim,ldim,1)

C  999 print *, 'zp= ?'
C      read (*,*) ze
C      call ev2cgf(isp,ldim,ldim,ldim,nev,eb,z,ze,Gij)
C      goto 999

      call rlse(odiawk)

      if (ipr .ge. 30) then
        j = min(9,ldimx)
        if (ipr .ge. 35) j = ldimx
#if LINUX_PGI
        do  18  ii = 1, 1
#else
        do  18  ii = 1, 2
#endif
        call awrit3(' SECM2C:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,qp)
   18   write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
        if (ipr .ge. 36 .and. nev .gt. 0) call awrit5(
     .    ' nev, nevmx, ldim=  %i  %i  %i  ev(nev) = %1;5d  efmax '//
     .    '= %1;5d',' ',80,i1mach(2),nev,nevmx,ldimx,eb(nev),efmax)
        call ftflsh(lgunit(1))
      endif
      if (ipr .ge. 110)
     .  call yprm('evec',2,z,ldimx*ldimx,ldimx,ldimx,ldimx)

C      if (switch(25)) then
C        if (lnc) call rx('pert cor not implemented')
C        if (nevmx .le. 0)
C     .    call rx('SECM2C: pert correction not allowed with nevmx=0')
C         call defrr(ode,-lihdim)
C         ltmp = iostr(8,'STR',nl,nbas,1,kap2,0,ckbas,-1,nsite,
C     .     oalph,oiax,ontab,os)
C        nsite = w(ontab+nbas)
C        call blcho(0,1,1,qp,nbas,plat,offH,1,nsite,w(oiax),
C     .    w(os),nl2,ldim,0,ldim,w(oh),w(osil),w(osii))
C        call rlse(oiax)
C        call defdc(osz,ld2)
C        call hmlt2p(ccor,w(occd),vmtz,elin,
C     .    ldim,lihdim,pph(1,1,isp),eb,nev,w(oh),z,w(osz),w(owk),w(ode))
C        call rlse(osz)
C        if (ccor) then
C          ltmp = iostr(8,'SDOT',nl,nbas,1,kap2,0,ckbas,-1,nsite,
C     .      oadot,oiax,ontab,os)
C          nsite = w(ontab+nbas)
C          call blcho(0,1,1,qp,nbas,plat,offH,1,nsite,w(oiax),
C     .      w(os),nl2,ldim,0,ldim,w(oh),w(1),w(1))
C          call dscal(ld22,-avw**2,w(oh),1)
C          call rlse(oiax)
C          call defdc(osz,ld2)
C          call hmltdp(elin,ldim,pph(1,1,isp),eb,nev,w(oh),z,w(osz),w(ode))
C        endif
C        call daxpy(nev,1d0,w(ode),1,eb,1)
C        if (ipr .ge. 1) then
C          do  28  ii = 1, 2
C            write(lgunit(ii),345) (eb(i),i=1,j)
C  345       format(' Incl. perturbation correction'/(255(9f8.4:/)))
C   28     continue
C        endif
CC       call rlse(ode)
C      endif
      call rlse(owk)
      if (ipr .ge. 110) call query('V<110 to skip matrix printing',-1,0)
      end
C      subroutine ev2cgf(isp,ni,nj,ldim,nbmx,eval,evc,ze,Gij)
CC- Make GF from eigenvectors
C      implicit none
C      integer isp,ldim,nbmx,ni,nj
C      double precision evc(ldim,ldim,2),eval(nbmx,isp)
C      double complex Gij(ni,nj),Ze
C
C      integer i,j,nu
C      double complex wk
C
C      do  20  i = 1, ni
C      do  20  j = 1, nj
C        wk = 0
C        do  22  nu = 1, ldim
C          wk = wk + (dcmplx(evc(i,nu,1),evc(i,nu,2)) *
C     .               dcmplx(evc(j,nu,1),-evc(j,nu,2))) /
C     .              (Ze-eval(nu,isp))
C   22   continue
C        Gij(i,j) = wk
C
C   20 continue
C
C      call zprm('gf',2,Gij,ldim,ldim,ldim)
C
C      end
C

