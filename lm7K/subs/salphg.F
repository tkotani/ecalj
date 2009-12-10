      subroutine salphg(ivl,nclus,nl,lsga,kmx,rsma,
     .  plat,alat,pos,ips,rsm,rhc,ehvl,cy,cg,indxcg,jcg,iax,sga,sgi)
C- Value-Laplacian transformation of a double set of functions:
C  Gaussians G_0L and either G_1L, Hsm_L, or Hsm-dot_L for one cluster 
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies functions used to build the value-Laplacian set
Ci          0 Gaussians and their Laplacians, G0 & G1
Ci          1 Gaussians and sm. Hankels, G0 & Hsm
Ci          2 Gaussians and energy derivatives of sm. Hankels, G0 & Hsm-dot
Ci   nclus :number of atoms in the current cluster
Ci   nl    :maximum angular momentum + 1, leading dimension of 
Ci          rsm, rhc, and ehvl
Ci   lsga  :total number of channels in the cluster
Ci          2*lsga is the leading dimension of sga and sgi
Ci   kmx   :kmax for the expansion of basis functions into polynomials P_kL 
Ci   rsma  :polynomial smoothing radii, depends on species but not on l 
Ci   plat  :primitive lattice vectors, in units of alat
Ci   alat  :lattice constant 
Ci   pos   :basis vectors, in units of alat
Ci   ips   :index to which species each site belongs
Ci   rsm   :l- and species-dependent smoothing radii for both Gaussians 
Ci          and Hankels
Ci   rhc   :l- and species-dependent augmentation radii 
Ci   ehvl  :l- and species-dependent energies for sm. Hankels (ivl = 1) 
Ci          or their energy derivatives (ivl = 2)
Ci          Not referenced if ivl = 0
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordon coefficients
Ci   iax   :neighbor table containing pair information (pairc.f)
Co Outputs
Co   sgi   :structure constants to construct the value-Laplacian set
Co          out of the given set of functions (see Remarks)
Co   sga   :sga is the inverse of sgi until the call to dgesv, then
Co          used as a work array (see Remarks)
Cl Local variables
Cl   lqinv :switch governing how matrix inversion is accomplished 
Cl          0 uses lapack dgesv
Cl          1 uses dqinvb
Cl          lqinv = 0 unless switch --qinv is passed through the
Cl          command line
Cl   nenv  :number of envelope functions on the head sphere of the cluster
Cl          For now nenv = (lmax+1)**2
Cb Bugs
Cb   lqinv = 1 is not implemented 
Cb   No check as to whether the CG arrays are large enough
Cr Remarks
Cr   1. The value-Laplacian transformation of a double set of functions
Cr   {F_iRL,i=0,1, where R and L run over all sites and L channels}  
Cr   brings F_iRL into the "value-Laplacian" set of functions U_jRL as
Cr           U_jRL = \sum_iR'L' sg(L',L,i,j,R') * F_iR'L'          (1)
Cr   where the summation runs over all sites and L channels of a given
Cr   cluster.
Cr
Cr   2. U_jRL has the property that its L-resolved value/Laplacian over  
Cr   the surface of a sphere centered at R' is zero unless the sites and 
Cr   indecies match:                                                    '
Cr      \lap^i U_jRL|    = \delta_ij \delta_LL' \delta_RR'         (2)
Cr                  |R'L'  
Cr
Cr   3. The significant speed up is achieved by recognizing that we do
Cr   not need property (2) to hold for ALL sites R in the cluster, but 
Cr   only for the center of the cluster R0. Thus we only have to build 
Cr   the U_jL (= U_jRL with R=R0) subset so that  
Cr      \lap^i U_jL|    = \delta_ij \delta_LL' \delta_R0,R'        (3)
Cr                 |R'L'  
Cr   Hence only index R' is kept in sg(...R') in Eq.(1) 
Cr
Cr   4. Transformation matrix sg(L',L,i,j,R') is stored here as a 
Cr   rectangular array sgi(L'iR',Lj). Notation LiR means a linear string 
Cr   with index ordering (((ilm=1,Lmax),ikap=1,2),ir=1,nclus)  
Cr   sgi is converted to the "standard" form (1) later on, when strux for
Cr   all clusters are glued together (addtog.f)
Cr
Cr   5. In terms of sgi(L'iR',Lj), Eq.(1) reads 
Cr           U_jL = \sum_i"R"L" sgi(L"i"R",Lj) * F_i"R"L"          (4)
Cr   Averaging (4) over sphere surfaces as in (3), we have   
Cr      \lap^i U_jL|    = \sum_i"R"L" sga(L'iR',L"i"R") * sgi(L"i"R",Lj)  
Cr                 |R'L'     
Cr                      =  \delta_ij \delta_LL' \delta_R0,R'       (5) 
Cr   where the square plan matrix sga is defined as 
Cr          sga(L'iR',L"i"R") := \lap^i F_i"R"L"|                  (6)  
Cr                                              |R'L'
Cr   6. The way to find sgi in practice is then 
Cr      (1) build sga 
Cr      (2) invert it, and
Cr      (3) keep only first 2*nenv columns provided that the cluster 
Cr   center always comes first in the list of sites (is in array iax).  
Cr   
Cr   7. Set ncheck = 1 to verify that sga*sgi is the unit matrix, ie
Cr   no accuracy loss occurs during inversion    
Cu Updates
Cu   15 May 08 (S. Lozovoi) Adapted from salph1.f
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ivl,nclus,nl,lsga,kmx(*),ips(*)                    
      integer niax
      parameter (niax=10)
      integer iax(niax,*)
      double precision rhc(nl,*),rsm(nl,*),ehvl(nl,*),rsma(*)
      double precision alat,pos(3,*),plat(3,3)
      double precision sga(2*lsga,2*lsga),sgi(2*lsga,2*lsga)
C For Clebsch Gordon coefficients:
      double precision cy(*),cg(*)
      integer jcg(*),indxcg(*)
C Local parameters
      logical cmdopt
      integer lqinv
      integer n0,ndim0,nclusm
      parameter (n0=10, ndim0=(n0+1)**2, nclusm=200)
      integer lclus,ipr,nenv,lsga2
      integer opivot,osg2,ounit
      integer nlmh,nlmg
C     integer pmax
      integer it,jt,ib,jb,ii,il
      integer iih,iig,iihm,iigm,ilmh,ilmg
      integer ll,ierr,ncheck                             
      integer isclus(nclusm),isj                   
      double precision ph(3),pg(3)
      double precision pclus(3,nclusm),
     .  rsmc(0:n0,nclusm),ehvlc(0:n0,nclusm)
      double precision gvl(ndim0,ndim0,0:1),hvl(ndim0,ndim0,0:1)
      double precision rr,drr2,distm,xx
      double precision tol
      character*80 outs
      parameter (tol=1d-12)

C Heap:
      integer w(1)
      common /w/ w

      call tcn('salphg')
C     stdo = nglob('stdo')
      call getpr(ipr)

#if DEBUG
C ncheck = 1 checks if the matrix inversion is OK
       ncheck = 1
#else
       ncheck = 0
#endif

C --- Checks ---
      if (nclus .lt. 1)
     .  call rxi('salphg: empty cluster encountered. nclus = ',nclus)
      if (nclus .gt. nclusm) call
     .  rxi('salphg: cluster is too big. Increase nclusm up to ',nclus)
      if (nl .gt. n0-1)
     .  call rxi('salphg: lmax is bigger then n0. lmax = ',nl-1)
c add here the check for CG lmax ...

C --- Preliminary setup ---
C      if (ivl .eq. 0) then
C        pmax = 2
C      else
C        pmax = 1
C      endif

      ib = iax(1,1)
      nenv = iax(9,1)
      do ii = 1, 3
        ph(ii) = pos(ii,ib)
      enddo

      lclus = 0
      do it = 1, nclus
        jb = iax(2,it)
        isj = ips(jb)
c ... coordinates of all atoms in the cluster
        rr = drr2(plat,pos(1,ib),pos(1,jb),
     .    iax(3,it),iax(4,it),iax(5,it),pg)
        do ii = 1, 3
          pclus(ii,it) = (pg(ii)+ph(ii))*alat
        enddo
c ... species index, rsm, and ehvl for all atoms in the cluster
        isclus(it) = isj
        nlmh = iax(9,it)
        do il = 0, ll(nlmh)
          rsmc(il,it) = rsm(il+1,isj)          
        enddo
        if (ivl .ne. 0) then
          do il = 0, ll(nlmh)
            ehvlc(il,it) = ehvl(il+1,isj)          
          enddo
        endif
        lclus = lclus + nlmh
      enddo
      if (lclus .ne. lsga)
     .  call rx('salphg: matrix size mismatch')

C --- Make matrix sga_ji = value/laplacian of F_Ri at sphere Rj, 
C     see Eq.(6) in Remarks ---
      iigm = 0
      do jt = 1, nclus
        isj = isclus(jt)
        nlmg = iax(9,jt)
C       lmaxg = ll(nlmg)
         
        iihm = 0
        do it = 1, nclus
C         isi = isclus(it)
          nlmh = iax(9,it)
C         lmaxh = ll(nlmh)

          call gh1c(ivl,pclus(1,it),rsmc(0,it),nlmh,ehvlc(0,it),
     .    pclus(1,jt),rsma(isj),nlmg,kmx(isj),rhc(1,isj),
     .    cy,cg,indxcg,jcg,ndim0,gvl,hvl)

c ... Append to sga
          iig = iigm
          do ilmg = 1, nlmg
            iig = iig + 1

            iih = iihm
            do ilmh = 1, nlmh
              iih = iih + 1
              sga(iig,iih)           = gvl(ilmg,ilmh,0)
              sga(iig,iih+nlmh)      = hvl(ilmg,ilmh,0)
              sga(iig+nlmg,iih)      = gvl(ilmg,ilmh,1)
              sga(iig+nlmg,iih+nlmh) = hvl(ilmg,ilmh,1)
            enddo
          enddo
          iihm = iihm + 2*nlmh
        enddo
        iigm = iigm + 2*nlmg
      enddo

C ... Check if size of sga is correct
      if (iihm .ne. iigm .or. iigm .ne. 2*lsga)
     .  call rx('salphg: matrix sga is not filled properly') 

C --- Invert sga ---
      lqinv = 0
      if (cmdopt('--qinv',6,0,outs)) lqinv = 1 

      if(lqinv .eq. 0) then
c ... initialize sgi as unit matrix        
        call dpzero(sgi, 4*lsga*nenv)
        do it = 1, 2*nenv 
          sgi(it,it) = 1d0
        enddo
c Save sga in w(osg2) and unit matrix in w(ounit) for future checks      
        if (ncheck .ne. 0) then
          lsga2 = (2*lsga)**2
          call defdr(ounit, lsga2)
          call defdr(osg2, lsga2)
          call dcopy(lsga2,sga,1,w(osg2),1)
          call dcopy(lsga2,sgi,1,w(ounit),1)
        endif
        ierr = 0
        call defi(opivot, 2*lsga)
        call dgesv(2*lsga,2*nenv,sga,2*lsga,w(opivot),sgi,2*lsga,ierr)
        call rlse(opivot)
      else
        call rx('salphg: direct inversion is not implemented.'//
     .    ' Use the default')
      endif
      if (ierr .ne. 0) 
     .  call rxi('salphg: failed to invert sga. ierr = ',ierr)
      
      if (ncheck .ne. 0) then
        xx = distm(2*lsga,w(osg2),2*lsga,sgi,2*lsga,w(ounit),
     .             2*nenv,2*lsga,2*nenv) 
        if (xx .gt. tol) call rx1(
     .    'salphg: loss of accuracy in sga inversion. |sga*sgi-I| =%;4g'
     .    ,xx)
        call rlse(ounit)
      endif

      call tcx('salphg')
      end

      double precision function distm(nla,a,nlb,b,nlc,c,n,m,k)
C- Find to what precision A(nxm)*B(mxk)=C(nxk) is satisfied elementwise
C ----------------------------------------------------------------
Ci Inputs
Ci   nla   :leading dimension of a
Ci   nlb   :leading dimension of b
Ci   nlc   :leading dimension of c
Ci   a,b,c :matrices to check
Ci   n,m,k :actual dimensions of a(n x m), b(m x k), and c(n x k)
Co Outputs
Co   distm :max(abs(A*B-C)) 
Cl Local variables
Cb Bugs
Cr Remarks
Cu Updates
Cu   26 Jun 08 First written 
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nla,nlb,nlc,n,m,k
      double precision a(nla,*),b(nlb,*),c(nlc,*)
C Local parameters
      integer it,jt
      double precision xx

      call dgemm('n','n',n,k,m,1d0,a,nla,b,nlb,-1d0,c,nlc)

      distm = 0d0
      do it = 1, n
        do jt = 1, k
          xx = abs(c(it,jt))
          if (xx .gt. distm) distm = xx
        enddo
      enddo

      end 

