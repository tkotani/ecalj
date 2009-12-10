      subroutine strg(ivl,nbas,nl,kmx,rsma,alat,plat,pos,ips,
     .  rsm,rmt,ehvl,cy,cg,indxcg,jcg,ntab,iax,ntabg,sg)
C- Structure matrix for value-Laplacian basis, all sites
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies the functions used to built the strux
Ci          ivl = 0 use G0 + G1
Ci          ivl = 1 use G0 + Hs  
Ci          ivl = 2 use G0 + Hs-dot  
Ci   nbas  :number of atoms in the basis 
Ci   nl    :maximum l quantum number + 1 
Ci   kmx   :kmax for decomposition of Gaussians and sm. Hankels into 
Ci         :polynomials P_kL for each species                     
Ci   rsma  :polynomial smoothing radii for each species 
Ci   plat,alat :primitive lattice translation vectors and scale  
Ci   pos   :basis vectors
Ci   ips   :index to which species each site belongs
Ci   rsm   :smoothing radii for each species and angular momentum 
Ci   rmt   :augmentation radii for each species and angular momentum 
Ci   ehvl  :fixed energies for Hankels or their energy derivatives 
Ci         :for each species and angular momentum (see ivl)
Ci         :not used if ivl = 0                               
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordon coefficients
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntabg :ntabg(ib) no. of neighbors in cluster centered at ib (pairg.f)
Cl Local variables
Cl   nclus :size of current cluster
Co Outputs
Co   sg    :screened strux in (nl**2,nl**2) blocks for all pairs
Co          (1..ntabg(ib),ib=1..nbas), see Remarks
Co Bugs
Co   For the moment all strux are treated as inequivalent
Co   Pass here lmxcg to check if the size of the CG arrays is appropriate
Cr Remarks
Cr   sg has actually 5 indecies which are not introduced here explicitly 
Cr   (see addtog.f)
Cr   
Cr   sg(L',L,i',i,R") defines the screening transformation from Gaussians
Cr   G_iRL to the value-Laplacian set of functions U_iRL as:
Cr     U_iRL = \sum_i'R'L' sg(L',L,i',i,R") * G_i'R'L'             (*)
Cr   where R" = R' - R
Cr
Cr   If ivl = 1, then smoothed Hankels Hs_RL are used instead  of G_1RL in (*)
Cr   If ivl = 2, then smoothed Hs-dot_RL are used in place of G_1RL in (*)
Cr
Cr   U_iRL have the property that the average value (i=0) or Laplacian (i=1)
Cr   of U_iRL*Y_L' over sphere at R is \delta_LL', and zero for all other 
Cr   spheres (see salphg.f for extended comments). 
Cr
Cr   sg is stored only for the head node of each of nbas clusters, ie 
Cr   for R = pos(ib)
Cu Updates
Cu   05 May 08 (S. Lozovoi) Hs and Hs-dot options added
Cu   27 Feb 08 (S. Lozovoi) Adapted from strrs.f 
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ivl,nbas,nl,niax
      parameter (niax=10)
      integer iax(niax,*),ntab(nbas+1),ntabg(nbas)
      integer ips(nbas),kmx(*)
      double precision rmt(nl,*),rsm(nl,*),ehvl(nl,*),rsma(*)
      double precision sg(nl*nl,nl*nl,*)
      double precision alat,plat(9),pos(3,nbas)
      integer indxcg(*),jcg(*)
      double precision cy(*),cg(*)
C Local parameters
      integer iat,ndimL,ndimL2,nl2,nitab,
     .  isum,nttab,nclus,iclus,stdo,
     .  ipr,nglob,nds
      integer osgi,osga

C heap:
      integer w(1)
      common /w/ w

      call tcn('strg')
      stdo = nglob('stdo')
      call getpr(ipr)
      nttab = ntab(nbas+1)
      nl2 = nl**2
C     Leading dimensions of sg
      nds = nl2

C --- Temporarily scale the lattice and basis vectors ---
c     call dscal(9,alat,plat,1)
c     call dscal(3*nbaspp,alat,bas,1)

C ... nitab is number of strux actually in sg array
      nitab  = 0

C --- Make the structure constants for all sites ---
      do  iat = 1, nbas

C       Number of atoms in this cluster
        nclus = ntabg(iat)
C       Offset to current cluster
        iclus = ntab(iat)+1
C       The true dimension of subblocks of sg 
        ndimL = isum(nclus,iax(9,iclus),niax)
C       Number of head functions in this sphere
C       nenv = iax(9,iclus)

C ...   Work arrays ...
        ndimL2 = 4*ndimL*ndimL
        call defdr(osgi, ndimL2)
        call defdr(osga, ndimL2)
         
        call salphg(ivl,nclus,nl,ndimL,kmx,rsma,plat,alat,pos,ips,
     .    rsm,rmt,ehvl,cy,cg,indxcg,jcg,iax(1,iclus),w(osga),w(osgi))

C ...   Append sgi to structure constant sg
        call addtog(nds,ndimL,iax(1,iclus),nclus,w(osgi),nitab,sg)

        nitab = nitab + nclus
        call rlse(osgi)
      enddo

C --- Symmetrize structure constants ---
c     call defi(osflg,-nitab)
c     call defi(osflgd,-nitab)
c     call symstr(0,nl2,nttab,iax,1,nel,w(osflg),s(1,1,1),s,xx)
c     if (ldot) call symstr(0,nl2,nttab,iax,1,nel,w(osflgd),
c    .  sdot(1,1,1),sdot,xx)

C --- Undo scaling of the lattice and basis vectors ---
c     call dscal(9,1/alat,plat,1)
c     call dscal(3*nbaspp,1/alat,bas,1)

C --- Info printout ---
      if (ipr .ge. 30) call awrit2(' strg:  generated %i '//
     .  'inequivalent strux from %i total',' ',80,stdo,nitab,nttab)

      call tcx('strg')

      end

