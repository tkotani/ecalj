      subroutine fitmrq(sctrl,sarray,slat,stb,nterm,nlmesp,nset,nl,nsp,
     .  nspc,nclass,npair,nkp,ldim,nelts,memode,mxnbr,idim,dclabl,itab,
     .  idec,itbcf,idcf,itbov,idov,itocf,idocf,iq1,ivso,idxdn,iam,npm,
     .  lmx,indxsh,rmaxh,eband,delta,qp,hso,v0,pnu,ifit,tabme,decay,
     .  tabcf,deccf,tabov,decov,tbocf,dcocf,qnu,vso,pot0,rl)
C- Perform Levenberg-Marquardt nonlinear least-squares fit of TB params
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: ltb lncol maxit tol
Ci       maxit: maximum number of Levenberg-Marquardt iterations
Ci       tol:   tolerance for chi-squared of fit
Ci     Stored:   
Ci     Passed to: xgetp mrqmin
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read:
Ci     Stored:   
Ci     Passed to: mrqmin
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat
Ci     Stored:   
Ci     Passed to: mrqmin
Ci   stb
Ci     Elts read: fmode nbfit ebfit
Ci       fmode: 0 fit parameters used in tbham, 1 not implemented yet
Ci       nbfit: lower and upper band indices to be fit for all k-points
Ci     Stored:   
Ci     Passed to: mrqmin
Ci   nterm: number of parameters for each matrix element;
Ci          also dimensions parameter tables.
Ci   nlmesp: (number of matrix elts/spin channel) * (# spin channels)
Ci   nset  : number of rules specifying matrix elements
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :Note: this needs to be fixed for nspc=2
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nclass:number of inequivalent classes
Ci   npair: number of matrix element pairs
Ci   nkp   :number of k-points (bzmesh.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nelts: first dimension of delta
Ci   memode:Defines functional form for matrix elements (rdtbh.f)
Ci   mxnbr :upper limit to total number of pairs in neighbor table
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   dclabl:class name, packed as a real number
Ci   itab,idec,itbcf,idcf,itbov,idov,itocf,idocf,iq1,ivso:
Ci         : switches to freeze parameters: =0 do not freeze, =1 freeze,
Ci         :=2 stick to symmetric value
Ci   idxdn :controls how orbital is to be included in the basis
Ci   iam   :iam(1,kk) and iam(2,kk) are the two classes for the kkth ME pair
Ci   npm
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   rmaxh :radius cutoff for neighbor table
Ci   eband: eigenvalues to be fit for each band and k-point
Ci   delta: holds electrostatic increments to the diagonal matrix
Ci   qp    :list of k-points
Ci   hso   :spin-orbit hamiltonian
Ci   v0    :pair potential (for printout only)
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   ifit  :work array for range of bands to fit at each k-point
Ci   tabme :a table of parameters, associated with a particular functional
Ci         :form (specified by memode) for TB matrix elments (rdtbh.f)
Ci   decay,deccf,decov,dcocf: exponential or power decay parameters
Ci         :matrix element [memode = 2, v_ij d^(-b); 3, v_ij exp(-c d)]
Co   tabov :table of overlap matrix elements, structured as tabme
Co   tbocf :table of crystal field MEs for overlap, structured as tabme
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   vso   :table of spin-orbit parameters
Co   pot0  :monopole potential at site R (eq. 7.81, Finnis)
Ci   rl    :if true, estat calculated for L>0
Ci   pnu: pnu's (for printout)
Co Outputs
Co   tabme,decay,tabcf,deccf,tabov,decov,tbocf,dcocf,qnu,vso:
Co     fitted TB parameters (initial guesses on input)
Cl Local variables
Cl   cnvg: tolerance for chi-squared of fit (unpacked from ctrl->tol)
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      double precision sctrl(1),sarray(1),slat(1),stb(1)
      integer nterm,nlmesp,nset,nl,nsp,nspc,nclass,npair,nkp,ldim,nelts,
     .  memode,mxnbr,idim
      integer itab(nterm,nlmesp,nset),idec(nlmesp,nset),
     .  itbcf(nterm,nlmesp,nset),idcf(nlmesp,nset),
     .  itbov(nterm,nlmesp,nset),idov(nlmesp,nset),
     .  itocf(nterm,nlmesp,nset),idocf(nlmesp,nset),
     .  iq1(nl,nsp,nclass),ivso(nl,nclass),idxdn(nl,nclass),
     .  iam(3,npair),npm(2,nclass),lmx(nclass),indxsh(1),ifit(2,nkp)
      double precision rmaxh
      double precision dclabl(nclass),eband(ldim*nsp,nkp),
     .  delta(nelts,1),qp(3,nkp),hso(nl**2,nl**2,4,2),v0(9,nset),
     .  pnu(nl,nsp,nclass),tabme(nterm,nlmesp,nset),decay(nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),deccf(nlmesp,nset),
     .  tabov(nterm,nlmesp,nset),decov(nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),dcocf(nlmesp,nset),
     .  qnu(3,nl,nsp,nclass),vso(nl,nclass),pot0(*)
      logical rl
C Heap
      integer w(1)
      common /w/ w
C Local parameters
      integer ltb,lncol,fmode,nbfit(2),maxit
      double precision ebfit(2),alat,cnvg,tol(3)
      integer nfac,npme,npq1,npdec,npso,npar,nvar,nvar2,ldimx,nbnd,
     .  nfit,iq,ilow,ihi,ib,nbmin,nbmax,i1mach,ifin,ierr,n,last,iprint,
     .  i,ifi,lgunit,iset,np,ip,il,jl,ic,it,isp,jc
      integer oivar,oivar2,oip1,oip2,opar,oiwk1,oefit,oiwk2,oiwk3,
     .  osig,oepar,odpar,obeta,otryp,odelp,oalp,ocov,oiwk4
      double precision efmax,d1mach,avgnb,alam,ochi,savc,chi0,rms0,chi,
     .  rmsdel,e1,e2,dval,dif
      double precision eps,eps2
      parameter(eps = 1d-4, eps2 = 1d-8)
      character*80 outs,strn1,strn2,clabli*8
      character*1 spd(5)
      character*3 sk(10)
      logical lov,cryf,ocryf,ldec,lso,low,cmdopt,bittst

      data spd /'s','p','d','f','g'/
      data sk /'ss0','sp0','pp0','pp1','sd0','pd0','pd1','dd0','dd1',
     .  'dd2'/

      call upack4('ctrl ltb lncol maxit tol',sctrl,ltb,lncol,maxit,tol)
      cnvg = tol(1)
      call upack3('tb fmode nbfit ebfit',stb,fmode,nbfit,ebfit)
      call upack1('lat alat',slat,alat)

      call rxx(fmode .ne. 0,'FITMRQ: only FMODE=0 for now')

      lov   = bittst(ltb,1)
      cryf  = bittst(ltb,2)
      ocryf = bittst(ltb,4)
      ldec  = memode .eq. 2 .or. memode .eq. 3
      lso   = bittst(lncol,4)

C --- Determine total number of parameters ---
      nfac = 1
      if (cryf) nfac = nfac + 1
      if (lov) then
        nfac = nfac + 1
        if (ocryf) nfac = nfac + 1
      endif
      npme = nterm*nlmesp*nset
      npq1 = nl*nsp*nclass
      npdec = 0
      npso = 0
      if (ldec) npdec = nlmesp*nset
      if (lso) npso = (nl-1)*nclass
      npar = (npme + npdec)*nfac + npq1 + npso

C --- Determine parameters to be varied, put parameters in array ---
      call defi(oivar,  -2*npar)
      call defi(oivar2, -2*npar)
      call defi(oip1, -4*(nterm+1)*nlmesp*nset)
      call defi(oip2, -2*nl*nsp*nclass)
      call defdr(opar, npar)
      call defi(oiwk1,   3*npar)
      call xgetp(sctrl,0,nterm,nlmesp,nset,nl,nsp,nclass,npair,npar,
     .  memode,itab,idec,itbcf,idcf,itbov,idov,itocf,idocf,iq1,ivso,iam,
     .  npm,idxdn,tabme,decay,tabcf,deccf,tabov,decov,tbocf,dcocf,qnu,
     .  vso,w(oiwk1),nvar,nvar2,w(oivar),w(oivar2),w(oip1),w(oip2),
     .  w(opar))
      call rlse(oiwk1)

C --- Determine eigenvalues to fit and set up array of band limits ---
      efmax = 999d0
      ldimx = ldim*nspc
      if (nbfit(1) .le. 0) nbfit(1) = 1
      if (nbfit(2) .lt. nbfit(1) .or. nbfit(2) .gt. ldimx)
     .  nbfit(2) = ldimx
      call icopy(nkp,nbfit(1),0,ifit,2)
      call icopy(nkp,nbfit(2),0,ifit(2,1),2)
      nbnd = nbfit(2) - nbfit(1) + 1
      nfit = nbnd*nkp
      if (ebfit(2)-ebfit(1) .gt. d1mach(3)) then
        nfit = 0
        nbmin = ldimx
        nbmax = 0
        do  30  iq = 1, nkp
          low = .true.
          ilow = ldimx + 1
          ihi = ldimx
          do  10  ib = 1, ldimx
            if (low .and. eband(ib,iq) .ge. ebfit(1)) then
              low = .false.
              ilow = ib
            endif
            if (eband(ib,iq) .gt. ebfit(2)) then
              ihi = ib - 1
              goto 20
            endif
   10     continue
   20     if (ilow .gt. nbfit(1)) ifit(1,iq) = ilow
          if (ihi  .lt. nbfit(2)) ifit(2,iq) = ihi
          if (ifit(1,iq) .gt. ifit(2,iq)) then
            ifit(1,iq) = 1
            ifit(2,iq) = 0
          endif
          nbnd = ifit(2,iq) - ifit(1,iq) + 1
          nfit = nfit + nbnd
          nbmin = min0(nbmin,nbnd)
          nbmax = max0(nbmax,nbnd)
   30   continue
        call rxx(nfit .eq. 0,'FITMRQ: no bands in energy range')
        avgnb = dble(nfit) / dble(nkp)
      endif

C --- Copy bands to be fit into a single array ---
      call defdr(oefit, nfit)
      call xcopb(nkp,ldimx,nfit,ifit,eband,w(oefit))

C --- Repeat Levenberg-Marquardt iterations until converged ---
      print *
      call awrit6(' FITMRQ: fmode=%i  nvar=%i  nvar2=%i  npar=%i  '
     .  //'nfit=%i  maxit=%i',outs,-80,i1mach(2),fmode,nvar,nvar2,npar,
     .  nfit,maxit)
      if (ebfit(2)-ebfit(1) .gt. d1mach(3)) then
        call awrit5('         nbmin=%i  nbmax=%i  nbavg=%d  ebmin=%d  '
     .    //'ebmax=%d',outs,-80,i1mach(2),nbmin,nbmax,avgnb,ebfit(1),
     .    ebfit(2))
      else
        call awrit3('         ibmin=%i  ibmax=%i  nbnd=%i',outs,-80,
     .    i1mach(2),nbfit(1),nbfit(2),nbnd)
      endif
      call defi(oiwk1,  npar)
      call defi(oiwk2,  npar)
      call defi(oiwk3,  npar)
      call defdr(osig,  nfit)
      call defdr(oepar, nfit)
      call defdr(odpar, nvar*nfit)
      call defdr(obeta, npar)
      call defdr(otryp, npar)
      call defdr(odelp, npar)
      call defdr(oalp,  npar*npar)
      call defdr(ocov,  npar*npar)
      call dvset(w(osig),1,nfit,1d0)
      alam = -1d0
      ochi = 1d99
      savc = 1d99
      ifin = 0
      ierr = 0
      do  40  n = 1, maxit+1
        call mrqmin(sctrl,sarray,slat,stb,npar,nterm,nlmesp,nset,
     .    nclass,nl,nsp,nspc,nkp,nelts,nfit,ldim,nvar,nvar2,mxnbr,
     .    memode,idim,w(oivar),w(oivar2),w(oip1),iam,npm,w(oip2),lmx,
     .    ifit,indxsh,alam,rmaxh,efmax,delta,qp,hso,w(oefit),w(osig),
     .    w(oiwk1),w(oiwk2),w(oiwk3),eband,w(oepar),w(odpar),w(obeta),
     .    w(otryp),w(odelp),decay,deccf,decov,dcocf,tabme,tabcf,tabov,
     .    tbocf,qnu,vso,pot0,rl,w(opar),chi0,rms0,chi,rmsdel,savc,
     .    w(oalp),w(ocov))
        if (iprint() .ge. 20) then
          print *
          call awrit4(' FITMRQ: it=%i  nit=%i  chi^2=%g  rms=%g',
     .      outs,-80,i1mach(2),n,maxit,chi,rmsdel)
        endif
        if ((chi .lt. ochi .and. dabs(ochi/chi - 1d0) .lt. eps)
     .    .or. (n .gt. maxit/2 .and. chi .lt. eps2)
     .    .or. rmsdel .lt. cnvg) then
          if (ifin .eq. 1) goto 50
          ifin = 1
        endif
        ochi = chi
   40 continue

C --- Fit did not converge ---
      ierr = 99
      write(*,500) n-1,chi,rmsdel,chi0,rms0
  500 format(/' FITMRQ (warning): Max iterations exceeded'//
     .        ' After ',i5,' iterations, chi^2 = ',g16.10,
     .        ', rms = ',g16.10/
     .        ' Prior to 1st iteration, chi^2 = ',g16.10,
     .        ', rms = ',g16.10)

C --- Final call to set covariance and curvature (alpha) matrices ---
   50 last = min0(n,maxit)
C      alam = 0
C      call mrqmin(sctrl,sarray,slat,stb,npar,nterm,nlmesp,nset,
C     .  nclass,nl,nsp,nkp,nelts,nfit,ldim,nvar,nvar2,mxnbr,memode,idim,
C     .  w(oivar),w(oivar2),w(oip1),iam,npm,w(oip2),lmx,ifit,indxsh,alam,
C     .  rmaxh,efmax,delta,qp,hso,w(oefit),w(osig),w(oiwk1),w(oiwk2),
C     .  w(oiwk3),eband,w(oepar),w(odpar),w(obeta),w(otryp),w(odelp),
C     .  decay,deccf,decov,dcocf,tabme,tabcf,tabov,tbocf,qnu,vso,w(opar),
C     .  chi0,rms0,chi,rmsdel,savc,w(oalp),w(ocov))

C --- Printout ---
      if (ierr .eq. 0) write(*,510) last,chi,rmsdel,chi0,rms0
  510 format(/' After ',i5,' iterations, chi^2 = ',g16.10,
     .        ', rms = ',g16.10/
     .        ' Prior to 1st iteration, chi^2 = ',g16.10,
     .        ', rms = ',g16.10)

C --- Print parameters ---
      if (iprint() .gt. 30) then
        if (ierr .ne. 0) write(*,520)
  520   format(/' The numbers below are for the best fit so far:')
        call shotbm(tabme,decay,v0,nset,npair,nterm,nl,nsp,nlmesp,
     .    memode,nclass,dclabl,iam,0,'Hamiltonian')
        if (cryf .and. memode .ne. 1) call shotbm(tabcf,deccf,v0,nset,
     .    npair,nterm,nl,nsp,nlmesp,memode,nclass,dclabl,iam,0,
     .    'Hamiltonian crystal field')
        if (lov .and. memode .ne. 1) call shotbm(tabov,decov,v0,nset,
     .    npair,nterm,nl,nsp,nlmesp,memode,nclass,dclabl,iam,0,
     .    'overlap')
        if (ocryf .and. memode .ne. 1) call shotbm(tbocf,dcocf,v0,nset,
     .    npair,nterm,nl,nsp,nlmesp,memode,nclass,dclabl,iam,0,
     .    'overlap crystal field')
      endif

C --- Print eigenvalues ---
      if (iprint() .ge. 40) then
        write(*,530)
        do  60  i = 1, nfit
          e1 = dval(w(oefit),i)
          e2 = dval(w(oepar),i)
          dif = e2 - e1
          write(*,540) i,e1,e2,dif
   60   continue
  530   format(/4x,'I',5x,'EIGENVALUE',10x,'FIT',12x,'DIFFERENCE')
  540   format(i5,2f16.10,g20.10)
      endif

C --- Parameters to LOG file (new ME category) ---
      ifi = lgunit(2)
      write(ifi,550) memode
  550 format(/'Fitted TB parameters:'//'ME ',i1)

C ... Determine classes for each parameter set
      if (cmdopt('-akm',4,0,outs)) then
        call defi(oiwk4,  2*nclass)
        call icopy(2*nclass,npm,1,w(oiwk4),1)
      endif
      do  90  iset = 1, nset
        np = 0
        call iinit(npm,2*nclass)
        do  70  ip = 1, npair
          if (iam(3,ip) .eq. iset) then
            np = np + 1
            npm(1,iam(1,ip)) = 1
            npm(2,iam(2,ip)) = 1
          endif
   70   continue
        call rxx(np .eq. 0,'FITMRQ: bad ME pointers')
        strn1 = ' '
        strn2 = ' '
        il = 0
        jl = 0
        do  80  ic = 1, nclass
          if (npm(1,ic) .eq. 1) then
            if (il .eq. 0) then
              il = 1
              call awrit1('%a%i',strn1,80,0,ic)
            else
              call awrit1('%a,%i',strn1,80,0,ic)
            endif
          endif
          if (npm(2,ic) .eq. 1) then
            if (jl .eq. 0) then
              jl = 1
              call awrit1('%a%i',strn2,80,0,ic)
            else
              call awrit1('%a,%i',strn2,80,0,ic)
            endif
          endif
   80   continue
        call skpblb(strn1,80,il)
        call skpblb(strn2,80,jl)

C ... Hamiltonian parmeters in current set
        write(ifi,560) strn1(1:il+1)//' '//strn2(1:jl+1)//' |'
        write(ifi,570) ((tabme(it,il,iset),it=1,nterm),il=1,nlmesp)
        if (ldec) write(ifi,580) (decay(il,iset),il=1,nlmesp)
        if (nvar2 .eq. 0) then
          write(ifi,590) ((itab(it,il,iset).ne.0,it=1,nterm),
     .      il=1,nlmesp)
          if (ldec) write(ifi,600) (idec(il,iset).ne.0,il=1,nlmesp)
        else
          write(ifi,595) ((itab(it,il,iset),it=1,nterm),il=1,nlmesp)
          if (ldec) write(ifi,605) (idec(il,iset),il=1,nlmesp)
        endif

C ... Crystal field Hamiltonian parmeters
        if (cryf) then
          write(ifi,560) '  &'
          write(ifi,570) ((tabcf(it,il,iset),it=1,nterm),il=1,nlmesp)
          if (ldec) write(ifi,580) (deccf(il,iset),il=1,nlmesp)
          write(ifi,590) ((itbcf(it,il,iset).ne.0,it=1,nterm),
     .      il=1,nlmesp)
          if (ldec) write(ifi,600) (idcf(il,iset).ne.0,il=1,nlmesp)
        endif

C ... Overlap parmeters
        if (lov) then
          write(ifi,560) '  @'
          write(ifi,570) ((tabov(it,il,iset),it=1,nterm),il=1,nlmesp)
          if (ldec) write(ifi,580) (decov(il,iset),il=1,nlmesp)
          if (nvar2 .eq. 0) then
            write(ifi,590) ((itbov(it,il,iset).ne.0,it=1,nterm),
     .        il=1,nlmesp)
            if (ldec) write(ifi,600) (idov(il,iset).ne.0,il=1,nlmesp)
          else
            write(ifi,595) ((itbov(it,il,iset),it=1,nterm),il=1,nlmesp)
            if (ldec) write(ifi,605) (idov(il,iset),il=1,nlmesp)
          endif

C ... Crystal field overlap parmeters
          if (ocryf) then
            write(ifi,560) '  %'
            write(ifi,570) ((tbocf(it,il,iset),it=1,nterm),il=1,nlmesp)
            if (ldec) write(ifi,580) (dcocf(il,iset),il=1,nlmesp)
            write(ifi,590) ((itocf(it,il,iset).ne.0,it=1,nterm),
     .        il=1,nlmesp)
            if (ldec) write(ifi,600) (idocf(il,iset).ne.0,il=1,nlmesp)
          endif
        endif

C ... Pair potential parmeters
        write(ifi,560) '  !'
        write(ifi,570) (v0(it,iset),it=1,9)

  560   format(3x,a)
  570   format(6x,4g16.8)
  580   format(7x,'DECAY='/(6x,4g16.8))
  590   format(7x,'FRZ=   ',28l2:/(14x,28l2))
  595   format(7x,'FRZ=   ',28i2:/(14x,28i2))
  600   format(7x,'FRZDEC=',28l2:/(14x,28l2))
  605   format(7x,'FRZDEC=',28i2:/(14x,28i2))
   90 continue

C --- Parameters to LOG file (new START category) ---
      call awrit2('START   NIT=%i CNVG=%g CNTROL=T',outs,-80,ifi,maxit,
     .  cnvg)
      do  100  ic = 1, nclass
        call r8tos8(dclabl(ic),clabli)
        write(ifi,610) clabli,(pnu(il,1,ic),il=1,nl)
        if (nsp .eq. 2) write(ifi,620) (pnu(il,2,ic),il=1,nl)
        write(ifi,630) (((qnu(it,il,1,ic),it=1,3),il=1,nl),isp=1,nsp)
  100 continue
  610 format(8x,'ATOM=',a8/10x,'P=',6f9.4:/(12x,6f9.4))
  620 format(12x,6f9.4)
  630 format(10x,'Q=',3g16.8:/(12x,3g16.8))

C --- Parameters to LOG file using AKM format ---
      if (cmdopt('-akm',4,0,outs)) then
C ... Site energies
        write(ifi,640) alat,alat*rmaxh,rmaxh,nl*nclass
        do  120  ic = 1, nclass
          do  110  il = 1, nl
            call awrit1('%x    %ie'//spd(il)//' ',outs,80,0,ic)
            write(ifi,650) outs(1:8),qnu(2,il,1,ic)
  110     continue
  120   continue
  640   format(/'      alat      rmax     rmxoa                idnt'/
     .    3f10.6,6x,'          '/'Site energies',4x,i3)
  650   format(a8,f14.9)

C ... Hamiltonian parameters
        call icopy(2*nclass,w(oiwk4),1,npm,1)
        write(ifi,660) 'Hamiltonian',nlmesp
        do  150  ic = 1, nclass
          do  140  jc = ic, nclass
            call meptr(ic,jc,iam,npm,iset)
            if (iset .eq. 0) goto 140
            do  130  il = 1, nlmesp
              call awrit2('%x   %i%i'//sk(il)//'  ',outs,80,0,ic,jc)
              write(ifi,670) outs(1:10),(tabme(it,il,iset), it=1,6)
  130       continue
  140     continue
  150   continue
  660   format(a11,6x,i3)
  670   format(a10,2(e17.9,2f8.4))

C ... Overlap parmeters
        if (lov) then
          write(ifi,660) 'Overlap    ',nlmesp
          do  180  ic = 1, nclass
            do  170  jc = ic, nclass
              call meptr(ic,jc,iam,npm,iset)
              if (iset .eq. 0) goto 170
              do  160  il = 1, nlmesp
                call awrit2('%x   %i%i'//sk(il)//'  ',outs,80,0,ic,jc)
                write(ifi,670) outs(1:10),(tabov(it,il,iset), it=1,6)
  160         continue
  170       continue
  180     continue
        endif

      endif

C --- Exit ---
      call rlse(oivar)

      end

      subroutine xgetp(sctrl,mod,nterm,nlmesp,nset,nl,nsp,nclass,npair,
     .  npar,memode,itab,idec,itbcf,idcf,itbov,idov,itocf,idocf,iq1,
     .  ivso,iam,npm,idxdn,tabme,decay,tabcf,deccf,tabov,decov,tbocf,
     .  dcocf,qnu,vso,iwrk,nvar,nvar2,ivar,ivar2,ip1,ip2,par)
C- Determine parameters to be varied, put or get parameters from array
C ----------------------------------------------------------------------
Ci Inputs
Ci   mod: =0 put parameters into a single array, =1 get them from array
Ci   nterm,nlmesp,nset,nl,nsp,nclass,npair,memode,iam,npm,idxdn
Ci   npar: total number of parameters
Ci   itab,idec,itbcf,idcf,itbov,idov,itocf,idocf,iq1,ivso:
Ci     switches to freeze parameters: =0 do not freeze, =1 freeze,
Ci     =2 stick to symmetric value
Ci   tabme,decay,tabcf,deccf,tabov,decov,tbocf,dcocf,qnu,vso: parameters
Ci   iwrk: work array for sticky parameter indices
Co Outputs
Co   nvar: number of parameters to vary out of npar total parameters
Co   nvar2: number of sticky parameters
Co   ivar(1,i): points to the position in the full list of ith variable
Co   ivar(2,i): paramter type of ith variable, with types 1 to 6:
Co     1: Hamiltonian parameter
Co     2: Hamiltonian crystal field parameter
Co     3: overlap parameter
Co     4: overlap crystal field parameter
Co     5: diagonal Hamiltonian parameter
Co     6: spin-orbit parameter
Co   ivar2(1,i): points to position in full list of ith sticky parameter
Co   ivar2(2,i): points to position in full list of parameter to stick to
Co   ip1(nterm+1,nlmesp,nset,4): pointer to locations in full list of
Co     variables for parameter types 1-4 (Hamiltonian and overlap
Co     parameters).  A first index of nterm+1 corresponds to the decay
Co     parameter (memode=2 or 3 only).
Co   ip2(nl,nsp,nclass,5:6): pointer to locations in full list of
Co     variables for parameter types 5 and 6 (diagonal Hamiltonian and
Co     spin-orbit parameters).  For the spin-orbit parameters a first
Co     index of 1 corresponds to p, 2 to d, etc.  A first index of nl
Co     is not used (there is no spin-orbit coupling for s-states).
Co   par: all of the parameters in a single array
Cr Remarks
Cr   The ordering of the parameters in the full list is:
Cr     1. Hamiltonian MEs by (term,L,spin)
Cr     If memode=2 or 3:
Cr       2. Hamiltonian decay parameters by (L,spin)
Cr     If crystal field terms:
Cr       3. Crystal field MEs by (term,L,spin)
Cr       If memode=2 or 3:
Cr         4. Crystal field decay parameters by (L,spin)
Cr     If nonorthogonal TB:
Cr       5. Overlap MEs by (term,L,spin)
Cr       If memode=2 or 3:
Cr         6. Overlap decay parameters by (L,spin)
Cr       If overlap crystal field terms :
Cr         7. Overlap crystal field MEs by (term,L,spin)
Cr         If memode=2 or 3:
Cr           8. Overlap crystal field decay parameters by (L,spin)
Cr     1 to 8 are repeated for each ME set
Cr     9. Diagonal Hamiltonian MEs by (L,spin,class)
Cr     If spin-orbit and nl>1:
Cr       10. Spin-orbit MEs by (L,nclass) excluding L=0
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision sctrl(1)
      integer mod,nterm,nlmesp,nset,nl,nsp,nclass,npair,npar,memode,
     .  nvar,nvar2
      integer itab(nterm,nlmesp,nset),idec(nlmesp,nset),
     .  itbcf(nterm,nlmesp,nset),idcf(nlmesp,nset),
     .  itbov(nterm,nlmesp,nset),idov(nlmesp,nset),
     .  itocf(nterm,nlmesp,nset),idocf(nlmesp,nset),
     .  iq1(nl,nsp,nclass),ivso(nl,nclass),iam(3,npair),npm(2,nclass),
     .  idxdn(nl,nclass),iwrk(3,npar),ivar(2,npar),ivar2(2,npar),
     .  ip1(nterm+1,nlmesp,nset,4),ip2(nl,nsp,nclass,5:6)
      double precision tabme(nterm,nlmesp,nset),decay(nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),deccf(nlmesp,nset),
     .  tabov(nterm,nlmesp,nset),decov(nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),dcocf(nlmesp,nset),
     .  qnu(3,nl,nsp,nclass),vso(nl,nclass),par(npar)
C Local parameters
      integer ltb,lncol,nsave,ip,iset,il,it,ic1,ic2,iset2,iv,ic,isp
      logical cryf,lov,ocryf,ldec,lso,bittst

      call upack2('ctrl ltb lncol',sctrl,ltb,lncol)

      lov   = bittst(ltb,1)
      cryf  = bittst(ltb,2)
      ocryf = bittst(ltb,4)
      ldec  = memode .eq. 2 .or. memode .eq. 3
      lso   = bittst(lncol,4)

      nvar = 0
      nvar2 = 0
      nsave = 0
      ip = 1

      do  30  iset = 1, nset
        do  20  il = 1, nlmesp
C --- Hamiltonian MEs ---
          do  10  it = 1, nterm
            if (mod .eq. 0) then
              par(ip) = tabme(it,il,iset)
              if (itab(it,il,iset) .eq. 0) then
                nvar = nvar + 1
                ivar(1,nvar) = ip
                ivar(2,nvar) = 1
                ip1(it,il,iset,1) = nvar
              elseif (itab(it,il,iset) .eq. 2) then
                call xgetc(iset,npair,iam,ic1,ic2)
                call meptr(ic2,ic1,iam,npm,iset2)
                if (iset .ne. iset2 .and. itab(it,il,iset2) .eq. 0) then
                  nvar2 = nvar2 + 1
                  ivar2(1,nvar2) = ip
                  iwrk(1,nvar2) = it
                  iwrk(2,nvar2) = il
                  iwrk(3,nvar2) = iset2
                endif
              endif
            else
              tabme(it,il,iset) = par(ip)
            endif
            ip = ip + 1
   10     continue

C --- Hamiltonian decay parameters ---
          if (ldec) then
            if (mod .eq. 0) then
              par(ip) = decay(il,iset)
              if (idec(il,iset) .eq. 0) then
                nvar = nvar + 1
                ivar(1,nvar) = ip
                ivar(2,nvar) = 1
                ip1(nterm+1,il,iset,1) = nvar
              elseif (idec(il,iset) .eq. 2) then
                call xgetc(iset,npair,iam,ic1,ic2)
                call meptr(ic2,ic1,iam,npm,iset2)
                if (iset .ne. iset2 .and. idec(il,iset2) .eq. 0) then
                  nvar2 = nvar2 + 1
                  ivar2(1,nvar2) = ip
                  iwrk(1,nvar2) = nterm + 1
                  iwrk(2,nvar2) = il
                  iwrk(3,nvar2) = iset2
                endif
              endif
            else
              decay(il,iset) = par(ip)
            endif
            ip = ip + 1
          endif
   20   continue
   30 continue

C --- Sticky Hamiltonian parameters ---
      if (mod .eq. 0) then
        do  35  iv = nsave+1, nvar2
          it = iwrk(1,iv)
          il = iwrk(2,iv)
          iset2 = iwrk(3,iv)
          ivar2(2,iv) = ivar(1,ip1(it,il,iset2,1))
   35   continue
        nsave = nvar2
      endif

C --- Crystal field MEs ---
      if (cryf) then
        do  60  iset = 1, nset
          do  50  il = 1, nlmesp
            do  40  it = 1, nterm
              if (mod .eq. 0) then
                par(ip) = tabcf(it,il,iset)
                if (itbcf(it,il,iset) .eq. 0) then
                  nvar = nvar + 1
                  ivar(1,nvar) = ip
                  ivar(2,nvar) = 2
                  ip1(it,il,iset,2) = nvar
                endif
              else
                tabcf(it,il,iset) = par(ip)
              endif
              ip = ip + 1
   40       continue

C --- Crystal field decay parameters ---
            if (ldec) then
              if (mod .eq. 0) then
                par(ip) = deccf(il,iset)
                if (idcf(il,iset) .eq. 0) then
                  nvar = nvar + 1
                  ivar(1,nvar) = ip
                  ivar(2,nvar) = 2
                  ip1(nterm+1,il,iset,2) = nvar
                endif
              else
                deccf(il,iset) = par(ip)
              endif
              ip = ip + 1
            endif
   50     continue
   60   continue
      endif

C --- Overlap MEs ---
      if (lov) then
        do  90  iset = 1, nset
          do  80  il = 1, nlmesp
            do  70  it = 1, nterm
              if (mod .eq. 0) then
                par(ip) = tabov(it,il,iset)
                if (itbov(it,il,iset) .eq. 0) then
                  nvar = nvar + 1
                  ivar(1,nvar) = ip
                  ivar(2,nvar) = 3
                  ip1(it,il,iset,3) = nvar
                elseif (itbov(it,il,iset) .eq. 2) then
                  call xgetc(iset,npair,iam,ic1,ic2)
                  call meptr(ic2,ic1,iam,npm,iset2)
                  if (iset .ne. iset2
     .                .and. itbov(it,il,iset2) .eq. 0) then
                    nvar2 = nvar2 + 1
                    ivar2(1,nvar2) = ip
                    iwrk(1,nvar2) = it
                    iwrk(2,nvar2) = il
                    iwrk(3,nvar2) = iset2
                  endif
                endif
              else
                tabov(it,il,iset) = par(ip)
              endif
              ip = ip + 1
   70       continue

C --- Overlap decay parameters ---
            if (ldec) then
              if (mod .eq. 0) then
                par(ip) = decov(il,iset)
                if (idov(il,iset) .eq. 0) then
                  nvar = nvar + 1
                  ivar(1,nvar) = ip
                  ivar(2,nvar) = 3
                  ip1(nterm+1,il,iset,3) = nvar
                elseif (idov(il,iset) .eq. 2) then
                  call xgetc(iset,npair,iam,ic1,ic2)
                  call meptr(ic2,ic1,iam,npm,iset2)
                  if (iset .ne. iset2 .and. idov(il,iset2) .eq. 0) then
                    nvar2 = nvar2 + 1
                    ivar2(1,nvar2) = ip
                    iwrk(1,nvar2) = nterm + 1
                    iwrk(2,nvar2) = il
                    iwrk(3,nvar2) = iset2
                  endif
                endif
              else
                decov(il,iset) = par(ip)
              endif
              ip = ip + 1
            endif
   80     continue
   90   continue

C --- Sticky overlap parameters ---
        if (mod .eq. 0) then
          do  95  iv = nsave+1, nvar2
            it = iwrk(1,iv)
            il = iwrk(2,iv)
            iset2 = iwrk(3,iv)
            ivar2(2,iv) = ivar(1,ip1(it,il,iset2,3))
   95     continue
          nsave = nvar2
        endif

C --- Overlap crystal field MEs ---
        if (ocryf) then
          do  120  iset = 1, nset
            do  110  il = 1, nlmesp
              do  100  it = 1, nterm
                if (mod .eq. 0) then
                  par(ip) = tbocf(it,il,iset)
                  if (itocf(it,il,iset) .eq. 0) then
                    nvar = nvar + 1
                    ivar(1,nvar) = ip
                    ivar(2,nvar) = 4
                    ip1(it,il,iset,4) = nvar
                  endif
                else
                  tbocf(it,il,iset) = par(ip)
                endif
                ip = ip + 1
  100         continue

C --- Overlap crystal field decay parameters ---
              if (ldec) then
                if (mod .eq. 0) then
                  par(ip) = dcocf(il,iset)
                  if (idocf(il,iset) .eq. 0) then
                    nvar = nvar + 1
                    ivar(1,nvar) = ip
                    ivar(2,nvar) = 4
                    ip1(nterm+1,il,iset,4) = nvar
                  endif
                else
                  dcocf(il,iset) = par(ip)
                endif
                ip = ip + 1
              endif
  110       continue
  120     continue
        endif
      endif

C --- Diagonal Hamiltonian MEs ---
      if (mod .eq. 0) then
        call dcopy(nl*nsp*nclass,qnu(2,1,1,1),3,par(ip),1)
        do  150  ic = 1, nclass
          do  140  isp = 1, nsp
            do  130  il = 1, nl
              if (iq1(il,isp,ic) .eq. 0 .and. idxdn(il,ic) .eq. 1) then
                nvar = nvar + 1
                ivar(1,nvar) = ip
                ivar(2,nvar) = 5
                ip2(il,isp,ic,5) = nvar
              endif
              ip = ip + 1
  130       continue
  140     continue
  150   continue
      else
        call dcopy(nl*nsp*nclass,par(ip),1,qnu(2,1,1,1),3)
        ip = ip + nl*nsp*nclass
      endif

C --- Spin-orbit MEs ---
      if (lso .and. nl .ge. 2) then
        do  170  ic = 1, nclass
          do  160  il = 1, nl-1
            if (mod .eq. 0) then
              par(ip) = vso(il,ic)
              if (ivso(il,ic) .eq. 0 .and. idxdn(il+1,ic) .eq. 1) then
                nvar = nvar + 1
                ivar(1,nvar) = ip
                ivar(2,nvar) = 6
                ip2(il,1,ic,6) = nvar
              endif
            else
              vso(il,ic) = par(ip)
            endif
            ip = ip + 1
  160     continue
  170   continue
      endif

      ip = ip - 1
      call rxx(ip .ne. npar,'XGETP: ip ne npar')
      call rxx(mod .eq. 0 .and. nvar .gt. npar,'XGETP: nvar gt npar')

      end

      subroutine xgetc(k,npair,iam,ic1,ic2)
C- Get two classes associated with the kth matrix element index
C ----------------------------------------------------------------------
Ci Inputs
Ci   k: matrix element index
Ci   npair: number of matrix element pairs
Ci   iam: see remarks
Co Outputs
Co   ic1,ic2: two classes associated with kth matrix element index
Cr Remarks
Cr   iam(1,i) and iam(2,i) are the two classes for the ith ME pair
Cr   k = iam(3,i) is a pointer to tabme(*,k) which holds the
Cr   matrix elements <iam(1,i) | H | iam(2,i)>
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer k,npair,ic1,ic2
      integer iam(3,npair)
C Local parameters
      integer i

      do  10  i = 1, npair
        if (iam(3,i) .eq. k) then
          ic1 = iam(1,i)
          ic2 = iam(2,i)
          return
        endif
   10 continue
      call rx('XGETC: no classes found')

      end

      subroutine xcopb(nkp,ldim,nfit,ifit,eband,efit)
C- Copy bands to be fit into a single array
C ----------------------------------------------------------------------
Ci Inputs
Ci   nkp,ldim
Ci   nfit: number of eigenvalues to fit out of ldim*nkp total bands
Ci   ifit(1,i),ifit(2,i): range of bands to fit for ith k-point
Ci   eband: eigenvalues for each band and k-point
Co Outputs
Co   efit: eigenvalues to be fit in a single array
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nkp,ldim,nfit
      integer ifit(2,nkp)
      double precision eband(ldim,nkp),efit(nfit)
C Local parameters
      integer ip,iq,ib

      ip = 0
      do  20  iq = 1, nkp
        do  10  ib = ifit(1,iq), ifit(2,iq)
          ip = ip + 1
          efit(ip) = eband(ib,iq)
   10   continue
   20 continue
      call rxx(ip .ne. nfit,'XCOPB: ip ne nfit')

      end
