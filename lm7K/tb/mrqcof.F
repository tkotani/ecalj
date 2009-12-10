      subroutine mrqcof(sctrl,sarray,slat,nterm,nlmesp,nset,nclass,nl,
     .  nsp,nspc,nkp,npar,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,
     .  iam,npm,ip2,lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,
     .  dcocf,tabme,tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,efit,
     .  sig,eband,epar,dpar,chi,alp,beta)
C- Evaluate linearized fitting matrix alpha and vector beta
C ----------------------------------------------------------------------
Ci Inputs
Ci   nterm,nlmesp,nset,nclass,nl,nsp,nkp,ldim,mxnbr,memode,idim
Ci   npar: total number of TB parameters
Ci   nelts: first dimension of delta
Ci   nfit: number of eigenvalues to fit out of ldim*nkp total bands
Ci   nvar: number of parameters to vary out of npar total parameters
Ci   iam,npm,lmx,indxsh
Ci   ip1: pointer to locations in full list of variables
Ci   ip2: pointer to locations in full list of variables
Ci   ifit(1,i),ifit(2,i): range of bands to fit for ith k-point
Ci   ivar(1,i): points to the position in the full list of ith variable
Ci   ivar(2,i): paramter type of ith variable, with types 1 to 6:
Ci   rmaxh,efmax,qp,hso
Ci   decay,deccf,decov,dcocf,tabme,tabcf,tabov,tbocf,qnu,vso: parameters
Ci   delta: holds electrostatic increments to the diagonal matrix
Ci   efit: eigenvalues to be fit in a single array
Ci   sig: standard deviations for each eigenvalue to be fit
Ci   eband: work array to hold bands
Ci   epar: work array for eigenvalues corresponding to current parameters
Ci   dpar: work array for derivatives of eigenvalues wrt parameters
Co Outputs
Co   chi: chi-squared value for fit with current parameters
Co   alp: alpha matrix, see Numerical Recipes
Co   beta: beta vector, see Numerical Recipes
Cr Remarks
Cr   Called by mrqmin.  See Eq. 14.4.8, p. 522 in Numerical Recipes.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      double precision sctrl(1),sarray(1),slat(1)
      integer nterm,nlmesp,nset,nclass,nl,nsp,nspc,nkp,npar,nelts,ldim,
     .  nfit,mxnbr,memode,nvar
      integer ip1(nterm+1,nlmesp,nset,4),iam(3,1),npm(2,nclass),
     .  ip2(nl,nsp,nclass,5:6),lmx(nclass),ifit(2,nkp),indxsh(1),
     .  ivar(2,npar)
      double precision rmaxh,efmax,chi
      double precision decay(nlmesp,nset),deccf(nlmesp,nset),
     .  decov(nlmesp,nset),dcocf(nlmesp,nset),tabme(nterm,nlmesp,nset),
     .  tabcf(nterm,nlmesp,nset),tabov(nterm,nlmesp,nset),
     .  tbocf(nterm,nlmesp,nset),qnu(3,nl,nsp,nclass),delta(nelts,1),
     .  qp(3,nkp),vso(nl,nclass),hso(nl**2,nl**2,4,2),efit(nfit),
     .  sig(nfit),eband(ldim*nsp,nkp),epar(nfit),dpar(nvar,nfit),
     .  alp(npar,npar),beta(npar),pot0(*)
      logical rl
C Local parameters
      integer i,j,k
      double precision del,sig2i,wt

C --- Get eigenvalues and derivatives for current parameter set ---
      call bndtbf(sctrl,sarray,slat,nterm,nlmesp,nset,nclass,nl,nsp,
     .  nspc,nkp,nelts,ldim,nfit,mxnbr,memode,nvar,ip1,iam,npm,ip2,
     .  lmx,ifit,indxsh,ivar,rmaxh,efmax,decay,deccf,decov,dcocf,tabme,
     .  tabcf,tabov,tbocf,qnu,delta,qp,vso,hso,pot0,rl,eband,epar,dpar)

C --- Initialize beta and (symmetric) alpha ---
      do  20  j = 1, nvar
        do  10  k = 1, j
          alp(j,k) = 0d0
   10   continue
        beta(j) = 0d0
   20 continue

C --- Loop over data and find chi-squared ---
      chi = 0d0
      do  50  i = 1, nfit
        del = efit(i) - epar(i)
        sig2i = 1d0 / (sig(i)*sig(i))
        do  40  j = 1, nvar
          wt = dpar(j,i)*sig2i
          do  30  k = 1, j
            alp(j,k) = alp(j,k) + wt*dpar(k,i)
   30     continue
          beta(j) = beta(j) + del*wt
   40   continue
        chi = chi + del*del*sig2i
   50 continue

C --- Fill in symmetric side of alpha ---
      do  70  j = 2, nvar
        do  60  k = 1, j - 1
          alp(k,j) = alp(j,k)
   60   continue
   70 continue

      end
