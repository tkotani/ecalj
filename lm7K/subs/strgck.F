      subroutine strgck(ivl,nds,nl,iat,ips,rsm,ehvl,alat,plat,pos,
     .  iax,ntab,ntabg,cy,xmrp,sg,slj,slk)
C- Evaluate value-Laplacian functions at point xmrp without using
C  one-center expansion (for testing)
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies the functions used to built the value-Laplacian 
Ci          set U_iRL (see Remarks)
Ci          ivl = 0: G0 & G1
Ci          ivl = 1: G0 & Hs
Ci          ivl = 2: G0 & Hs-dot
Ci   nds   :leading dimensions of sg 
Ci   nl    :maximum l quantum number + 1, leading dimension of rsm,ehvl
Ci   iat   :center R of the cluster for which the value-Laplacian functions 
Ci          are constructed
Ci   ips   :index to which species each site belongs
Ci   rsm   :smoothing radii for each species and angular momentum
Ci   ehvl  :Hankel energies for each species and angular momentum
Ci          not used if ivl = 0
Ci   plat,alat :primitive lattice translation vectors and scale
Ci   pos   :basis vectors
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairc.f)
Ci   ntabg :ntabg(ib) no. of neighbors in cluster centered at ib (pairg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   xmrp  :Cartesian coordinates of the point at which to evaluate the 
Ci          value-Laplacian functions
Ci   sg    :screened strux in (nl**2,nl**2) blocks for all pairs
Ci          (1..ntabg(ib),ib=1..nbas)
Co Outputs
Co   slj(i,L) :approximate U_iRL, kappa- and L-resolved value-Laplacian 
Co            functions at xmrp   
Co   slk(i,L) :slk(i,L) = \lap slj(i,L)
Cr Remarks
Cr   slj are calculated by applying the value-Laplacian strux           
Cr   to a set of actial functions evaluated at point r (= xmrp) rather 
Cr   then to their expansion into polynomials
Cr           slj(i,L) = \sum_jR'L' sg(L',L,j,i,R') * F_jR'L'(r)    (1)
Cr   where the summation runs over all sites and L channels of a given
Cr   cluster, and F_jR'L'(r) = F_jL'(r-R') are Gaussians if j=0 and
Cr   G1/Hs/Hs-dot (depending on ivl) if j=1. 
Cr
Cr   Since sg were built for functions expanded into polynomials,   
Cr   slj and slk are not exactly U_iRL and \lap U_iRL, but tend 
Cr   to them as the polynomial expansion becomes more and more accurate.
Cu Updates
Cu   15 May 08 Adapted from strck.f
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ivl,iat,nds,nl
      integer niax,nkap
      parameter (niax=10, nkap=2)
      integer iax(niax,*),ntab(*),ntabg(*),ips(*)
      double precision alat,plat(3,3),pos(3,*),rsm(nl,*),ehvl(nl,*)
      double precision cy(*),xmrp(3),slj(nkap,nds),slk(nkap,nds),
     .  sg(nds,nds,nkap,nkap,*)
C Local parameters
      integer nclus,nlmh,nlmg,ll,lmaxg,ikap,pmax
      integer offR,isg,it,isj,ib,jb,ii,itt,ilmg,ilmh,il,jhd 
      integer ip,ill
      integer n0,pmaxx,nlm0
      parameter (n0=10, pmaxx=2, nlm0=(n0+1)**2)
      double precision ddot,xmrpp(3),xmrp0(3),pclus(3),yl(nlm0)
      double precision rsmc(0:n0),ehvlc(0:n0),xcheck,tol,xm1,xmm1,xm2
      double precision rsx,ehx,rsm0,eh0,cc,fac
      double precision gex(0:pmaxx,nlm0),ggrad(3,nlm0)
      double precision hkl(0:pmaxx,0:n0),hex(0:pmaxx,nlm0)
      double precision hs(0:n0),dhs(0:n0),ddhs(0:n0)
      double precision hsp(0:n0),dhsp(0:n0),ddhsp(0:n0)
      data xmrp0/1d-5,2d-5,3d-5/,tol/1d-12/

         
C --- Offset to iax table for cluster connected to R
      offR  = ntab(iat)
      nclus = ntabg(iat)
C --- Offset to strux for this cluster 
      isg = 0
      if (iat .gt. 1) then
        do it = 1, iat-1
          isg = isg + ntabg(it)
        enddo
      endif  

      if (ivl .eq. 0) then
        pmax = 2
      else
        pmax = 1
        jhd = 10*(ivl-1)+2
      endif

C --- Checks ---
      if (nds .gt. n0*n0)
     .  call rxi('strgck: nds is bigger than n0^2. nds = ',nds)
      if (nclus .lt. 1)
     .  call rxi('strgck: empty cluster encountered. nclus = ',nclus)
      if (nclus .gt. ntab(iat+1)-ntab(iat))
     .  call rxi('strgck: cluster is too big. nclus = ',nclus)

C ... sort, lmax, and coordinates/alat of the head of the cluster 
      ib = iax(1,offR+1)
C     isi = ips(ib)
      nlmh = iax(9,offR+1)
           
      call dpzero(slj,nkap*nds)
      call dpzero(slk,nkap*nds)

c ... begin loop over all atoms in the cluster
      do it = 1, nclus
        itt = it + offR      
        jb = iax(2,itt)
        isj = ips(jb)
        nlmg = iax(9,itt)
        call acoord(ib,it,alat,plat,pos,iax,ntab,pclus)
          
        lmaxg = ll(nlmg) 
        do il = 0, lmaxg
          rsmc(il) = rsm(il+1,isj)
        enddo
        if (ivl .ne. 0) then
          do il = 0, lmaxg
            ehvlc(il) = ehvl(il+1,isj)
          enddo
        endif

c ... coordinates xmrp relative to atom it: xmrpp = x - R'  
        do ii = 1, 3
          xmrpp(ii) = xmrp(ii) - pclus(ii)
        enddo
      
c ... if xmrp accidentally hits a regular site, make a small offset
c     and print a warning
        xcheck = ddot(3,xmrpp,1,xmrpp,1)
        if (xcheck .le. tol*tol) then
          do ii = 1, 3
            xmrpp(ii) = xmrp0(ii) 
          enddo
          call info5(10,1,0,
     .      ' strgck: Warning! |xmrpp(R,R'')|=%g for R=%i and R''=%i',
     .      dsqrt(xcheck),iat,it,0,0) 
          call info2(10,0,1,
     .      '%11f Point is offset by vector xmrp0 = %3:1;3g',
     .      xmrp0,0)
        endif

c ... make solid Gaussians at xmrpp
        call solgsg(xmrpp,rsmc,lmaxg,pmax,cy,pmaxx,gex,ggrad)
        
        if (ivl .eq. 0) then
          do ip = 0, pmax-1
            do ilmg = 1, nlmg
              hex(ip,ilmg) = gex(ip+1,ilmg)
            enddo
          enddo
        else
c ... make solid Hsm or Hsm-dot at xmrpp
          call sylm(xmrpp,yl,lmaxg,xm2)
          xm1 = dsqrt(xm2)

          rsx = -1d2
          ehx = 1d2
c ... radial part
          do il = lmaxg, 0, -1
            rsm0 = rsmc(il)
            eh0 = ehvlc(il)
            if (dabs(rsm0-rsx) + dabs(eh0-ehx) .gt.
     .        tol) then
              call hanszd(jhd,xm1,eh0,rsm0,il,
     .          hs,dhs,ddhs,hsp,dhsp,ddhsp)
              if (ivl .eq. 1) then
                do ill = 0, il
                  hkl(0,ill) = hs(ill)
                  hkl(1,ill) = ddhs(ill)
                enddo
              else
                do ill = 0, il
                  hkl(0,ill) = hsp(ill)
                  hkl(1,ill) = ddhsp(ill)
                enddo
              endif
            endif
            rsx = rsm0
            ehx = eh0
          enddo
c ... scale by r^{-l}
          if (lmaxg .ge. 1) then
            xmm1 = 1d0/xm1
            fac = 1d0
            do il = 1, lmaxg
              fac = fac*xmm1 
              do ip = 0, pmax
                hkl(ip,il) = hkl(ip,il)*fac
              enddo
            enddo
          endif

c ... multiply by YL 
          do  ilmg = 1, nlmg
            il = ll(ilmg)
            cc = yl(ilmg)*cy(ilmg)
            do ip = 0, pmax
              hex(ip,ilmg) = hkl(ip,il)*cc
            enddo
          enddo
        endif

c ... combine with sg to evaluate the value-Laplacian set slj
c     and their Laplacians slk
        do ilmg = 1, nlmg
          do ikap = 1, nkap
            do ilmh = 1, nlmh
              slj(ikap,ilmh) = slj(ikap,ilmh) +
     .                       gex(0,ilmg)*sg(ilmg,ilmh,1,ikap,isg+it) +
     .                       hex(0,ilmg)*sg(ilmg,ilmh,2,ikap,isg+it) 
              slk(ikap,ilmh) = slk(ikap,ilmh) +
     .                       gex(1,ilmg)*sg(ilmg,ilmh,1,ikap,isg+it) +
     .                       hex(1,ilmg)*sg(ilmg,ilmh,2,ikap,isg+it) 
            enddo
          enddo
        enddo
        

c ... end loop over cluster atoms
      enddo
      
      end

