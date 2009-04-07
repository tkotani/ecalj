      subroutine mrqmin(sctrl,sarray,slat,stb,npar,nterm,nlmesp,nset,
     .  nclass,nl,nsp,nspc,nkp,nelts,nfit,ldim,nvar,nvar2,mxnbr,memode,
     .  idim,ivar,ivar2,ip1,iam,npm,ip2,lmx,ifit,indxsh,alam,rmaxh,
     .  efmax,delta,qp,hso,efit,sig,iwk1,iwk2,iwk3,eband,epar,dpar,beta,
     .  tryp,delp,decay,deccf,decov,dcocf,tabme,tabcf,tabov,tbocf,qnu,
     .  vso,pot0,rl,par,chi0,rms0,chi,rmsdel,savc,alp,cov)
C- Performs one Levenberg-Marquardt iteration
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read:
Ci     Stored:   
Ci     Passed to: mrqcof xgetp
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read:
Ci     Stored:   
Ci     Passed to: mrqcof
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read:
Ci     Stored:   
Ci     Passed to: mrqcof
Ci   stb
Ci     Elts read: alam alsc
Ci     Stored:   
Ci     Passed to:
Ci   npar  :total number of TB parameters
Ci   nterm :number of parameters for each matrix element;
Ci         :also dimensions parameter tables.
Ci   nlmesp: (number of matrix elts/spin channel) * (# spin channels)
Ci   nset  : number of rules specifying matrix elements
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci         :Note: this needs to be fixed for nspc=2
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nkp   :number of k-points (bzmesh.f)
Ci   nelts :first dimension of delta
Ci   nfit  :number of eigenvalues to fit out of ldim*nkp total bands
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nvar  :number of parameters to vary out of npar total parameters
Ci   nvar2 :number of sticky parameters
Ci   mxnbr :upper limit to total number of pairs in neighbor table
Ci   memode:Defines the functional form for TB matrix elements (rdtbh.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   ivar(1,i): points to the position in the full list of ith variable
Ci   ivar(2,i): paramter type of ith variable, with types 1 to 6:
Ci   ivar2(1,i): points to position in full list of ith sticky parameter
Ci   ivar2(2,i): points to position in full list of parameter to stick to
Ci   ip1: pointer to locations in full list of variables
Ci   ip2: pointer to locations in full list of variables
Ci   ifit(1,i),ifit(2,i): range of bands to fit for ith k-point
Ci   iam   :iam(1,kk) and iam(2,kk) are the two classes for the kkth ME pair
Ci   npm
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   alam  :lambda variable, see remarks and Numerical Recipes
Ci   alam0 :initial value of alam, see Numerical Recipes
Ci   alsc  :scale factor for alam, see Numerical Recipes
Ci   rmaxh :cutoff radius for pair neighbor table
Ci   efmax :largest eigenvalue for which to find eigenvectors
Ci   delta :holds electrostatic increments to the diagonal matrix
Ci   qp    :list of k-points
Ci   hso   :spin-orbit hamiltonian
Ci   efit: eigenvalues to be fit in a single array (see routine xcopb)
Ci   sig: standard deviations for each eigenvalue to be fit
Ci   iwk1,iwk2,iwk3: integer work arrays for call to gaussj
Ci   eband: work array to hold bands
Ci   epar: work array for eigenvalues corresponding to current parameters
Ci   dpar: work array for derivatives of eigenvalues wrt parameters
Ci   beta: work array for beta vector, see Numerical Recipes
Ci   tryp: work array for trial parameters
Ci   delp: work array for change in parameters
Ci   pot0  :monopole potential at site R (eq. 7.81, Finnis)
Ci   rl    :if true, estat calculated for L>0
Co Outputs
Co   decay,deccf,decov,dcocf: exponential or power decay parameters
Co         :matrix element [memode = 2, v_ij d^(-b); 3, v_ij exp(-c d)]
Co   tabme :a table of parameters, associated with a particular functional
Co          form (specified by memode) for TB matrix elments (rdtbh.f)
Co          updated here
Co   qnu   :energy-weighted moments of the sphere charges
Co          updated here?
Co   vso   :table of spin-orbit parameters
Co          updated here
Co   tabcf: table of crystal field MEs 
Co          updated here
Co   tabov: table of overlap matrix elements, structured as tabme
Co          updated here
Co   tbocf: table of crystal field MEs for overlap, structured as tabme
Co          updated here
Co   par: updated TB parameters in a single array
Co   ... All of the above are updated in this routine
Co   chi0: chi-squared value for fit with initial parameters
Co   rms0: RMS for fit with initial parameters
Co   chi: chi-squared value for fit with current parameters
Co   rmsdel: RMS for fit with current parameters
Co   savc: chi-squared value from previous iteration
Co   alp: curvature (alpha) matrix, see Numerical Recipes
Co   cov: covariance matrix, see Numerical Recipes
Co Outputs
Co   decay,deccf,decov,dcocf,tabme,tabcf,tabov,tbocf,qnu,vso:
Co     updated TB parameters
Cr Remarks
Cr   See pp 521-528 in Numerical Recipes.
Cr   Set alam < 0 for initialization.  Set alam = 0d0 for final call to
Cr   set covariance and curvature (alpha) matrices.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      double precision sctrl(1),sarray(1),slat(1),stb(1)
      integer npar,nterm,nlmesp,nset,nclass,nl,nsp,nspc,nkp,nelts,nfit,
     .  ldim,nvar,nvar2,mxnbr,memode,idim
      integer ivar(2,npar),ivar2(2,npar),ip1(nterm+1,nlmesp,nset,4),
     .  iam(3,1),npm(2,nclass),ip2(nl,nsp,nclass,5:6),lmx(nclass),
     .  ifit(2,nkp),indxsh(1),iwk1(npar),iwk2(npar),iwk3(npar)
      double precision alam,rmaxh,efmax,chi0,rms0,chi,rmsdel,savc
      double precision delta(nelts,1),qp(3,nkp),hso(nl**2,nl**2,4,2),
     .  efit(nfit),sig(nfit),eband(ldim*nsp,nkp),epar(nfit),
     .  dpar(nvar,nfit),beta(npar),tryp(npar),delp(npar),
     .  decay(nlmesp,nset),deccf(nlmesp,nset),decov(nlmesp,nset),
     .  dcocf(nlmesp,nset),tabme(nterm,nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),tabov(nterm,nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),qnu(3,nl,nsp,nclass),vso(nl,nclass),
     .  par(npar),alp(npar,npar),cov(npar,npar),pot0(*)
      logical rl
C Local parameters
      integer kk,j,ihit,k,iprint,i1mach
      double precision alam0,alsc,dif
      character*80 outs
      character*16 ptype(6)
      data ptype /'Hamiltonian','Crystal Field','Overlap',
     .  'Overlap Crys. F.','Ham. diagonal','Ham. spin-orbit'/

      call upack2('tb alam alsc',stb,alam0,alsc)

C --- Initialization ---
      if (alam .lt. 0d0) then
        kk = nvar + 1

C --- Does ivar contain a proper permutation of the parameters ---
        do  20  j = 1, npar
          ihit = 0
          do  10  k = 1, nvar
            if (ivar(1,k) .eq. j) ihit = ihit + 1
   10     continue
          if (ihit .eq. 0) then
            ivar(1,kk) = j
            kk = kk + 1
          elseif (ihit .gt. 1) then
            call rx('MRQMIN: duplicate parameter in ivar')
          endif
   20   continue
        call rxx(kk .ne. (npar+1),
     .    'MRQMIN: Improper permutation in ivar')
        alam = alam0
        call mrqcof(sctrl,sarray,slat,nterm,nlmesp,nset,nclass,nl,nsp,
     .    nspc,nkp,npar,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,iam,
     .    npm,ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,
     .    dcocf,tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,
     .    efit,sig,eband,epar,dpar,chi,alp,beta)
        savc = chi
        chi0 = chi
        rmsdel = dsqrt(chi / nfit)
        rms0 = rmsdel
        call dcopy(npar,par,1,tryp,1)
C   ... Printout
        if (iprint() .ge. 30) then
          call awrit3('%N MRQMIN: Initial values: '
     .      //'chi^2=%g  rms=%g  alam=%g',outs,-80,i1mach(2),
     .      chi,rmsdel,alam)
          if (iprint() .ge. 50) then
            write(*,500)
            do  30  j = 1, nfit
              dif = epar(j) - efit(j)
              write(*,510) j,efit(j),epar(j),dif
   30       continue
  500       format(/4x,'I',5x,'EIGENVALUE',10x,'FIT',12x,'DIFFERENCE')
  510       format(i5,2f16.10,g20.10)
          endif
        endif
      endif

C --- Alter linearized fitting matrix by augmenting diagonal elements ---
      do  50  j = 1, nvar
        do  40  k = 1, nvar
          cov(j,k) = alp(j,k)
   40   continue
        cov(j,j) = alp(j,j)*(1d0 + alam)
        delp(j) = beta(j)
   50 continue

C --- Matrix solution ---
      call gaussj(nvar,npar,1,1,iwk1,iwk2,iwk3,cov,delp)

C --- Once converged evaluate covariance matrix with alam=0 ---
      if (alam .eq. 0) then
        call covsrt(nvar,npar,ivar,cov)
        return
      endif

C --- Trial parameters ---
      do  60  j = 1, nvar
        tryp(ivar(1,j)) = par(ivar(1,j)) + delp(j)
   60 continue
      do  70  j = 1, nvar2
        tryp(ivar2(1,j)) = tryp(ivar2(2,j))
   70 continue
C ... Printout
      if (iprint() .ge. 40) then
        write(*,520)
        do  80  j = 1, nvar
          write(*,530) ivar(1,j),par(ivar(1,j)),tryp(ivar(1,j)),delp(j),
     .      ptype(ivar(2,j))
   80   continue
  520   format(/4x,'I',4x,'OLD PARAM',8x,'TRIAL PARAM',9x,'DIFFERENCE',
     .    9x,'TYPE')
  530   format(i5,3g18.10,3x,a)
      endif
C ... Put trial parameters tryp into tables
      call xgetp(sctrl,1,nterm,nlmesp,nset,nl,nsp,nclass,0,npar,memode,
     .  0,0,0,0,0,0,0,0,0,0,0,0,0,tabme,decay,tabcf,deccf,tabov,decov,
     .  tbocf,dcocf,qnu,vso,0,j,k,0,0,0,0,tryp)
      call mrqcof(sctrl,sarray,slat,nterm,nlmesp,nset,nclass,nl,nsp,
     .  nspc,nkp,npar,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,iam,
     .  npm,ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,
     .  dcocf,tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,
     .  efit,sig,eband,epar,dpar,chi,cov,delp)

      if (iprint() .ge. 30) call awrit2('%x%N MRQMIN: '
     .  //'old chi^2=%g  new chi^2=%g',outs,80,0,savc,chi)
      if (chi .lt. savc*0.9999d0) then
C --- Success, accept new solution and decrease alam ---
        alam = alam / alsc
        savc = chi
        do  100  j = 1, nvar
          do  90  k = 1, nvar
            alp(j,k) = cov(j,k)
   90     continue
          beta(j) = delp(j)
          par(ivar(1,j)) = tryp(ivar(1,j))
  100   continue
        do  110  j = 1, nvar2
          par(ivar2(1,j)) = tryp(ivar2(2,j))
  110   continue

      else
C --- Failure, increase alam ---
        alam = alsc*alam
        chi = savc
      endif
      rmsdel = dsqrt(chi / nfit)
C ... Printout
      if (iprint() .ge. 30) then
        call awrit1('%a  alam=%g',outs,-80,-i1mach(2),alam)
        if (iprint() .ge. 50) then
          write(*,500)
          do  120  j = 1, nfit
            dif = epar(j) - efit(j)
            write(*,510) j,efit(j),epar(j),dif
  120     continue
        endif
      endif

C --- Put parameters par into tables ---
      call xgetp(sctrl,1,nterm,nlmesp,nset,nl,nsp,nclass,0,npar,memode,
     .  0,0,0,0,0,0,0,0,0,0,0,0,0,tabme,decay,tabcf,deccf,tabov,decov,
     .  tbocf,dcocf,qnu,vso,0,j,k,0,0,0,0,par)

      end
