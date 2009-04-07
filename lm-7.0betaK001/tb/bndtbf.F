      subroutine bndtbf(sctrl,sarray,slat,nterm,nlmesp,nset,nclass,nl,
     .  nsp,nspc,nkp,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,iam,npm,
     .  ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,dcocf,
     .  tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,eband,epar,
     .  dpar)
C- Get eigenvalues and their parameter derivatives for each k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   nterm,nlmesp,nset,nclass,nl,nsp,nkp,ldim,mxnbr,memode
Ci   nelts: first dimension of delta
Ci   nfit: number of eigenvalues to fit out of ldim*nkp total bands
Ci   nvar: number of parameters to vary out of all parameters
Ci   iam,npm,lmx,indxsh
Ci   ip1(nterm+1,nlmesp,nset,4): pointer to locations in full list of
Ci     variables for parameter types 1-4 (Hamiltonian and overlap
Ci     parameters).  A first index of nterm+1 corresponds to the decay
Ci     parameter (memode=2 or 3 only).
Ci   ip2(nl,nsp,nclass,5:6): pointer to locations in full list of
Ci     variables for parameter types 5 and 6 (diagonal Hamiltonian and
Ci     spin-orbit parameters).  For the spin-orbit parameters a first
Ci     index of 1 corresponds to p, 2 to d, etc.  A first index of nl
Ci     is not used (there is no spin-orbit coupling for s-states).
Ci   ifit(1,i),ifit(2,i): range of bands to fit for ith k-point
Ci   ivar(1,i): points to the position in the full list of ith variable
Ci   ivar(2,i): paramter type of ith variable, with types 1 to 6:
Ci     1: Hamiltonian parameter
Ci     2: Hamiltonian crystal field parameter
Ci     3: overlap parameter
Ci     4: overlap crystal field parameter
Ci     5: diagonal Hamiltonian parameter
Ci     6: spin-orbit parameter
Ci   rmaxh,efmax,qp,hso
Ci   decay,deccf,decov,dcocf,tabme,tabcf,tabov,tbocf,qnu,vso: parameters
Ci   delta: holds electrostatic increments to the diagonal matrix
Ci   eband: work array to hold bands
Co Outputs
Co   epar: eigenvalues for current parameter set for all k-points
Co   dpar: derivatives of eigenvalues wrt TB parameters for all k-points
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision sctrl(1),sarray(1),slat(1)
      integer nterm,nlmesp,nset,nclass,nl,nsp,nkp,nelts,ldim,nfit,mxnbr,
     .  memode,nvar
      integer ip1(nterm+1,nlmesp,nset,4),iam(3,1),npm(2,nclass),
     .  ip2(nl,nsp,nclass,5:6),lmx(nclass),ifit(2,nkp),indxsh(1),
     .  ivar(2,nvar)
      logical rl
      double precision rmaxh,efmax
      double precision decay(nlmesp,nset),deccf(nlmesp,nset),
     .  decov(nlmesp,nset),dcocf(nlmesp,nset),tabme(nterm,nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),tabov(nterm,nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),qnu(3,nl,nsp,nclass),delta(nelts,1),
     .  qp(3,nkp),vso(nl,nclass),hso(nl**2,nl**2,4,2),
     .  eband(ldim*nsp,nkp),epar(nfit),dpar(nvar,nfit),pot0(*)
C Heap
      integer w(1)
      common /w/ w
C Local parameters
      integer i,j,iprint,mxcsiz,nttab,nsite,i1mach,ipr,nspc,ldimx,iev,
     .  ikp,iprx,nev,ltb,lncol,lpgf(2),nbas,npadl,npadr,nclasp,nbasp,
     .  fopn,isp
      integer ontab,oiax,onpr,oh,oh0,odh,oov,owk,os,ozll,oips,oipc,
     .  oinitc,opos
      double precision ckbas,cksumf,plat(3,3),kap2(20)
      logical lov,cryf,ocryf,addsll,lso,fitpar,bittst,ldum1,ldum2,iostr,
     .  cmdopt
      character*120 outs

C --- Setup ---
      call upack4('ctrl ltb lncol lpgf nbas',sctrl,ltb,lncol,lpgf,nbas)
      call upack('array npadl npadr nclasp oips oipc',sarray,npadl,
     .  npadr,nclasp,oips,oipc)
      nbasp = nbas + npadl + npadr
      call upack1('array ohave',sarray,oinitc)
      call upack2('lat plat opos',slat,plat,opos)
      ckbas = cksumf(w(opos),3*nbas)

C --- Initialization ---
      lov    = bittst(ltb,1)
      cryf   = bittst(ltb,2)
      ocryf  = bittst(ltb,4)
      addsll = bittst(ltb,8) .and. lov
      lso    = bittst(lncol,4)

      call rxx(lso .and. nsp .ne. 2,'BNDTBF: must set NSPIN=2 for SO=T')

C --- Get neighbor table iax for each atom in the cluster ---
      if (lpgf(1) .ne. 0) then
        i = 2
        j = 1
      else
        i = 3
        j = -1
      endif
C ... Make nttab,ontab,oiax
      mxcsiz = mxnbr
      call pairs(nbas,nbasp,1d0,plat,rmaxh/2,w(opos),
     .  -1,i,j,w,nttab,ontab,oiax,mxcsiz)
C ... Patch iax(6) for the padded basis
      if (nbasp .gt. nbas)
     .  call pairp6(nbas,npadl,npadr,w(oiax),w(ontab))
      nsite = nttab
      call mkiaxd(nsite,lmx,w(oips),w(oiax))

C ... tb programs not converted to ntab from npr yet ...
      call defi(onpr,2*nbas)
      call npr2tb(1,nbasp,w(onpr),w(ontab))

C --- Slater-Koster Hamiltonian ---
      fitpar = .true.
      oov = 1
      if (nspc .eq. 2) call rx('update tbham with nsp=2, nspc=2')
      call tbham(nsp,nspc,sctrl,sarray,slat,fitpar,nvar,ip1,ivar,nlmesp,
     .  memode,decay,deccf,decov,dcocf,iam,npm,nterm,nset,
     .  tabme,tabcf,tabov,tbocf,nsite,w(onpr),oiax,oh,oh0,odh,oov,1,1,1)

C ... diagonal elements
      ldum1 = bittst(ltb,2**15)
      if (ldum1) ltb = ltb - 2**15
      if (nspc .eq. 2) call rx('update tbdiag; eliminate nsp')
      call tbdiag(ltb,nbasp,nl,nsp,nspc,nsp,w(oipc),nsite,w(onpr),
     .  w(oinitc),qnu,nelts,delta,w(1),w(oh),w(oh0),w(oov))
      if (ldum1) ltb = ltb + 2**15

C --- Accumulate parameter derivatives for diagonal hamiltonian MEs ---
      call tbfdia(nbas,nl,nsp,nclass,nsite,w(onpr),ip2,w(oipc),w(odh))

C ... What about fit for spin-orbit parameters?

C --- Add to Hamiltonian MEs: hLL' -> hLL' + ebarLL' * sLL' ---
      if (addsll) then
C ... Fit does not work here yet, should it?
        ldum1 = bittst(ltb,16)
        ldum2 = bittst(ltb,128)
        if (ldum1) ltb = ltb - 16
        if (ldum2) ltb = ltb - 128
        call defdr(owk,nl*nbas)
        call addes(nsite,nl,nsp,nbasp,ltb,w(onpr),w(oiax),w(oov),
     .    0,w(owk),w(oh),0)
        call rlse(owk)
        if (ldum1) ltb = ltb + 16
        if (ldum2) ltb = ltb + 128
      endif

C --- Write Hamiltonian to STR file ---
      if (lpgf(1) .ne. 0 .or. cmdopt('-dumph',6,0,outs)) then
        call npr2tb(0,nbasp,w(onpr),w(onpr))
        call defdr(owk,-nbasp*nl**2)
        ldum1 = iostr(1,'STR',nl,nbasp,1,kap2,0,ckbas,-1,nsite,owk,oiax,
     .    onpr,oh)
        nsite = w(onpr+nbasp)
        call fclose(fopn('STR'))
        call npr2tb(1,nbasp,w(onpr),w(onpr))
        call rlse(owk)
      endif

C --- Printout ---
      if (iprint() .ge. 50) then
        call defdr(owk,-nbasp*nl**2)
C   ... Printout Hamiltonian
        print *
        print *,'Real-space TB hamiltonian :'
        do  10  j = 1, nsp**2
          if (nsp .eq. 2) print *
          if (j .eq. 1 .and. nsp .eq. 2)
     .      print *,'Spin: Up-Up, w/o spin-orbit'
          if (j .eq. 2) print *,'Spin: Up-Down, w/o spin-orbit'
          if (j .eq. 3) print *,'Spin: Down-Up, w/o spin-orbit'
          if (j .eq. 4) print *,'Spin: Down-Down, w/o spin-orbit'
          os = oh + (j-1)*nl**4*nsite*i1mach(18)
          call shostr(nl**2,nsite,nbasp,plat,w(opos),0,w(owk),w(oiax),
     .      w(ontab),w(os),1,1,0d0,1,1d0)
C     ... Printout Overlap
          if (lov) then
            print *,'overlap matrix'
            os = oov + (j-1)*nl**4*nsite*i1mach(18)
            call shostr(nl**2,nsite,nbasp,plat,w(opos),0,w(owk),w(oiax),
     .        w(ontab),w(os),1,1,0d0,1,1d0)
          endif
   10   continue
C   ... Printout parameter derivatives
        do  30  i = 1, nvar
          print *
          if (ivar(2,i) .eq. 1 .or. ivar(2,i) .eq. 5
     .                        .or. (ivar(2,i) .eq. 6 .and. lso)) then
            write(*,500) ivar(1,i)
  500       format(' Derivatives for Hamiltonian parameter: ',i4)
          elseif (ivar(2,i) .eq. 3 .and. lov) then
            write(*,510) ivar(1,i)
  510       format(' Derivatives for overlap parameter: ',i4)
          elseif (ivar(2,i) .eq. 2 .and. cryf) then
            write(*,520) ivar(1,i)
  520       format(' Derivatives for crystal field parameter: ',i4)
          elseif (ivar(2,i) .eq. 4 .and. ocryf) then
            write(*,530) ivar(1,i)
  530       format(' Derivatives for overlap crystal field parameter: ',
     .        i4)
          endif
          do  20  j = 1, nsp**2
            if (nsp .eq. 2) print *
            if (j .eq. 1 .and. nsp .eq. 2) print *,'Spin: Up-Up'
            if (j .eq. 2) print *,'Spin: Up-Down'
            if (j .eq. 3) print *,'Spin: Down-Up'
            if (j .eq. 4) print *,'Spin: Down-Down'
            os = odh
     .        + (nsp**2*(i-1)+(j-1))*nl**4*nsite*i1mach(18)
            call shostr(nl**2,nsite,nbasp,plat,w(opos),0,w(owk),w(oiax),
     .        w(ontab),w(os),1,1,0d0,1,1d0)
   20     continue
   30   continue
        call rlse(owk)
      endif

C --- Determine verbosity for secmtb ---
      ipr = 0
      if (iprint() .ge. 40) ipr = 1
      if (iprint() .ge. 45) ipr = 2
      if (iprint() .ge. 50) ipr = 3

C --- Loop over k-points ---
      ldimx = ldim*nspc
      call defdc(ozll,ldimx*ldimx)
      call dpzero(dpar,nvar*nfit)
      iev = 0
      do  40  ikp = 1, nkp
        if (ifit(1,ikp) .gt. ifit(2,ikp)) goto 40
        iprx = ipr
        if (ikp .ge. 1) iprx = max(ipr-1,0)
        if (iprint() .ge. 40 .and. mod(ikp,100) .eq. 0) iprx=max(ipr,1)

C --- Make bands and eigenvectors ---
        if (nspc .eq. 2 .or. nsp .eq. 2) call rx('need isp here')
        isp = 1
        call secmtb(sctrl,plat,nbas,nl,nspc,nsp,isp,lmx,w(oipc),indxsh,
     .    ldimx,ifit(2,ikp),efmax,ikp,nkp,qp(1,ikp),nsite,w(oiax),
     .    w(onpr),w(oh),vso,hso,w(oov),pot0,rl,iprx,.true.,nev,
     .    w(ozll),eband(1,ikp),w,w)

C --- Make derivatives of eigenvalues wrt TB parameters ---
        call defdc(owk,ldimx*ldimx)
        call derivp(nkp,nbas,nclass,nsite,nl,nspc,ldimx,nvar,nfit,iprx,
     .    ikp,nev,ifit,w(oipc),ivar,lmx,indxsh,w(oiax),w(onpr),
     .    qp(1,ikp),plat,w(odh),w(ozll),eband(1,ikp),lso,lov,cryf,ocryf,
     .    w(owk),iev,dpar)
        call rlse(owk)
   40 continue
      call rlse(onpr)

C --- Copy bands to be fit into a single array ---
      call xcopb(nkp,ldimx,nfit,ifit,eband,epar)

      end
