      subroutine strg1c(ivl,nds,nlmy,nl,iat,ixi,ips,rsm,ehvl,rmt,kmx,
     .  rsma,alat,plat,pos,iax,ntab,ntabg,cy,cg,indxcg,jcg,sg,slj,slk)
C- 1c decomposition of screened basis functions and their Laplacians 
C  via the polynomial expansion (for testing)                                     
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies the functions used in the screened basis
Ci          ivl = 0 Gaussians and their Laplacians, G0 & G1
Ci          ivl = 1 Gaussians and sm. Hankels, G0 & Hs  
Ci          ivl = 2 Gaussians and energy derivatives of 
Ci                  sm. Hankels, G0 & Hs-dot  
Ci   nds   :leading dimensions of sg
Ci   nlmy  :leading dimensions of slj and slk
Ci   nl    :upper bound for l+1, leading dimension of rsm and rmt 
Ci   iat   :site R at which the screened basis functions are defined
Ci   ixi   :Indicates site R' for one-center expansion.  ixi is index
Ci         :to pair table relative to start of cluster at iat
Ci         :Thus ixi=1 points to head; ixi=2 points to 1st neighbor, etc
Ci   ips   :index to which species each site belongs
Ci   rsm   :Gaussian/Hankel smoothing radii for each species and 
Ci          angular momentum 
Ci   ehvl  :energies for sm. Hankels (ivl=1) or their energy derivatives 
Ci         :(ivl=2) for each species and angular momentum (see ivl)
Ci         :Not referenced if ivl = 0                               
Ci   rmt   :augmentation radii for each species and angular momentum 
Ci   kmx   :kmx(is) is kmax for decomposition of Gaussians or sm. Hankels  
Ci         :into polynomials P_kL for species 'is'                    
Ci   rsma  :P_kL smoothing radii for each species  
Ci   alat  :lattice constant
Ci   plat  :primitive lattice translation vectors in units of alat
Ci   pos   :basis vectors
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Ci   ntabg :ntabg(ib) no. of neighbors in cluster centered at ib (pairg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordon coefficients
Ci   sg    :structure constant matrix defining the screened basis
Co Outputs
Co   slj   :slj(i,L',L) is L'-component of 1c expansion of U_iRL around site ixi
Co          i = 0 and i = 1 correspond to value and Laplacian parts
Co          of the value-Laplacian unit basis derived for site iat 
Co   slk   :slk(i,L',L) = \lap slj(i,L',L) 
Cb Bugs
Cb   need to pass lmxcg to check if CG arrays are large enough
Cr Remarks
Cr   Screened basis U_iRL is the unit value-Laplacian set of functions having
Cr   the property that their values/laplacians at sphere R' are zero except
Cr   when indecies match:
Cr     \lap^i' U_iRL(R'L') = \delta_ii' * \delta_LL' * \delta_RR'
Cr   
Cr   U_iRL are constructed from a double set of functions F_0, F_1 as
Cr      U_iRL = \sum_i"R"L" sg_i"R"L",iRL * F_i"R"L"                  (*)
Cr   where the sum over i" runs over 0,1, R" and L" are sites and L channels
Cr   in a given cluster, and R is alwayas the center of the cluster. 
Cr   
Cr   F_0 are always the Gaussians G0, whereas F_1 are either G1, Hs, or Hs-dot     
Cr   depending on switch ivl (see above).
Cr   
Cr   Program evaluates U_iRL(R'L') as array slj, and \lap U_iRL(R'L') as 
Cr   array slk with R = iat and R' = ixi using (*) and 1c decomposition of F_i.
Cu Updates
Cu   05 May 08 (S. Lozovoi) Hs and Hs-dot added
Cu   27 Feb 08 (S. Lozovoi) First written 
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ivl,iat,nds,nlmy,nl,ixi
      integer niax,nkap
      parameter (niax=10, nkap=2)
      integer iax(niax,*),ntab(*),ntabg(*)
      integer ips(*),kmx(*)
      double precision alat,plat(3,3),pos(3,*)
      double precision rmt(nl,*),rsm(nl,*),ehvl(nl,*),rsma(*)
      double precision sg(nds,nds,nkap,nkap,*)
      double precision slj(nlmy,nlmy,nkap),slk(nlmy,nlmy,nkap)
      integer indxcg(*),jcg(*)
      double precision cy(*),cg(*)
C Local parameters
      integer nclus,nlmh,nlmg,ll,lh,lg,kmax 
      integer offR,isg,it,isi,isj,ib,jb,itt,ig,ih,il
      integer nlmf,if,ikap                          
      integer n0,nclusm,ndim0,pmx                       
      parameter (n0=10, nclusm=200, ndim0=(n0+1)**2, pmx=1) 
      double precision pclus(3,nclusm)  
      double precision rsmg(0:n0),rhc(0:n0),eg(0:n0)  
      double precision gvl(ndim0,ndim0,0:pmx),hvl(ndim0,ndim0,0:pmx)
      double precision sg1,sg2

      nclus = ntabg(iat)
C ... offR and isg are the offsets to iax table and sg matrix, respectively      
      offR = ntab(iat)
      isg = 0
      if (iat .gt. 1) then
        do it = 1, iat-1
          isg = isg + ntabg(iat)
        enddo
      endif

C --- Checks ---
      if (nds .gt. n0**2)
     .  call rxi('strg1c: nds is bigger than n0**2. nds = ',nds)  
      if (nclus .le. 0)
     .  call rxi('strg1c: empty cluster encountered for iat = ',iat)  
      if (nclus .gt. ntab(iat+1)-ntab(iat))
     .  call rxi('strg1c: size of the cluster is too big.'//  
     .  ' nclus = ',nclus)  
      if (nclus .gt. nclusm)
     .  call rxi('strg1c: size of the cluster exceeds the maximum.'//  
     .  ' nclus = ',nclus)  

C --- ixi-independent setup ---
C ... Head of the cluster
C     ic = iax(1,offR+1)      
      nlmf = iax(9,offR+1)      

c ... find coordinates of all atoms in the cluster
      do it = 1, nclus
        call acoord(iat,it,alat,plat,pos,iax,ntab,pclus(1,it))
      enddo

c ... Expansion site      
      ib = iax(2,offR+ixi)      
      isi = ips(ib)
      nlmh = iax(9,offR+ixi)      
      lh = ll(nlmh) 
      do il = 0, lh
        rhc(il) = rmt(il+1,isi)
      enddo
      kmax = kmx(isi)

      call dpzero(slj,nkap*nlmy*nlmy)
      call dpzero(slk,nkap*nlmy*nlmy)

C --- begin loop over all atoms in the cluster
      do it = 1, nclus
        itt = it+offR
        jb = iax(2,itt)      
        isj = ips(jb)

        nlmg = iax(9,itt)      
        lg = ll(nlmg) 
        do il = 0, lg
          rsmg(il) = rsm(il+1,isj)
        enddo
        if (ivl .ne. 0) then
          do il = 0, lg
            eg(il) = ehvl(il+1,isj)
          enddo
        endif 

c ... 1c-decomposition of F_0 and F_1 around ixi 
        call gh1c(ivl,pclus(1,it),rsmg,nlmg,eg,
     .    pclus(1,ixi),rsma(isi),nlmh,kmax,rhc,
     .    cy,cg,indxcg,jcg,ndim0,gvl,hvl)

c ... combine with sg to make slj and slk 
c index-site correspondence: iat --> if, ixi --> ih, it --> ig

        do ig = 1, nlmg
          do if = 1, nlmf
            do ikap = 1, nkap
              sg1 = sg(ig,if,1,ikap,isg+it)
              sg2 = sg(ig,if,2,ikap,isg+it)
              do ih = 1, nlmh
                slj(ih,if,ikap) = slj(ih,if,ikap) +
     .            gvl(ih,ig,0)*sg1 + hvl(ih,ig,0)*sg2 
                slk(ih,if,ikap) = slk(ih,if,ikap) +
     .            gvl(ih,ig,1)*sg1 + hvl(ih,ig,1)*sg2 
              enddo
            enddo
          enddo
        enddo

c ... end loop over cluster atoms
      enddo

      end

