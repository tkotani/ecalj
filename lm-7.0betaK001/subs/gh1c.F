      subroutine gh1c(ivl,ph,rsmh,nlmh,ehvl,pg,rsmg,nlmg,kmax,rhc,
     .  cy,cg,indxcg,jcg,ndim,gvl,hvl)
C- L-resolved value and Laplacian of Gaussians/Hankels centered 
C  at ph over a sphere centered at pg using 1c decomposition 
C ----------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies the functions used in the screened basis
Ci          ivl = 0 Gaussians and their Laplacians, G0 & G1
Ci          ivl = 1 Gaussians and sm. Hankels, G0 & Hs  
Ci          ivl = 2 Gaussians and energy derivatives of 
Ci                  sm. Hankels, G0 & Hs-dot  
Ci   ph    :head site 
Ci   rsmh  :l-dependent smoothing radii for Gaussians and Hankels 
Ci          centered at the head site 
Ci   nlmh  :L-cutoff for Gaussians / sm. Hankels at the head site 
Ci   ehvl  :energies for sm. Hankels (ivl=1) or their energy derivatives 
Ci         :(ivl=2) for each angular momentum.
Ci         :Not referenced if ivl = 0                               
Ci   pg    :expansion site 
Ci   rsmg  :P_kL smoothing radius at the expansion site 
Ci   nlmg  :L-cutoff for P_kL at the expansion site
Ci   nds   :leading dimensions of sg
Ci   kmax  :Gaussians or sm. Hankels are decomposed into 
Ci         :polynomials P_kL up to k = kmax
Ci   rhc   :l-dependent augmentation radii at the expansion site
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg,indxcg,jcg: Clebsch Gordon coefficients
Ci   ndim  :Leading dimensions of gvl and hvl
Co Outputs
Co   gvl   :L-decomposed value and Laplacians of generalised Gaussians G_pL
Co   hvl   :L-decomposed value and Laplacian of sm. Hankel (ivl=1) or
Co          the energy derivative of sm. Hankel (ivl=2) 
Ci          Not referenced if ivl = 0                               
Cl Local variables
Cl  spg    :Expansion coefficients of G_kL in terms of P_kL at ph.
Cl  shr    :Expansion coefficients of H_kL in terms of P_kL at ph, 
Cl         :real part
Cl  shi    :Expansion coefficients of H_kL in terms of P_kL at ph, 
Cl         :imaginary part
Cl  sdr    :Expansion coefficients of dot-H_kL in terms of P_kL at ph
Cb Bugs
Cb   need to pass lmxcg to check if CG arrays are sufficiently large
Cr Remarks
Cr   gvl(L',L,p) :L'th component of G_pL centered at ph 
Cr   hvl(L',L,p) :the meaning depends on ivl 
Cr                ivl = 0 L'th component of G_p+1L centered at ph, 
Cr                        ie hvl(L',L,p) = \lap gvl(L',L,p) 
Cr                ivl = 1 L'th component of Hsm_pL centered at ph 
Cr                ivl = 2 L'th component of Hsm-dot_pL centered at ph 
Cr    
Cr   Irrespective of a function type (Gaussians, sm. Hankels, or Hsm-dot), 
Cr   gvl and hvl are obtained through the polynomial expansion of function 
Cr   in question at a neighboring site. If ph and pg sites coincide, then
Cr   the actual function rather than its polynomial decomposition is used.
Cu Updates
Cu   05 May 08 Hs and Hs-dot added
Cu   27 Feb 08 First written 
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ivl,nlmh,nlmg,kmax,ndim
      double precision ph(3),pg(3)
      double precision rsmh(0:*),rsmg,ehvl(0:*),rhc(0:*)
      double precision gvl(ndim,ndim,0:1),hvl(ndim,ndim,0:1)
      integer indxcg(*),jcg(*)
      double precision cy(*),cg(*)
C Local parameters
      integer n0,pmaxx
      parameter (n0=10, pmaxx=2)
      integer ndimx,kmaxx                       
c     parameter (ndimx=(n0+1)**2, kmaxx=n0)    
      integer lmaxh,lmaxg,ilmg
      integer pmax,ll,jvsl,jhd
      integer ii,ip,im,il,ill 
      double precision spg(0:kmax,ndim,0:pmaxx,ndim),
     . shr(0:kmax,ndim,0:pmaxx-1,ndim),
c    . shi(0:kmax,ndim,0:n0,ndim), 
     . shdr(0:kmax,ndim,0:pmaxx-1,ndim)
      double precision pkl(0:kmax,0:n0),gpkl(0:kmax,0:n0)
      double precision gkl(0:pmaxx,0:n0),hkl(0:1,0:n0),
     .  rsml(0:n0),rhcl(0:n0),dums(ndim,ndim)
      double precision hs(0:n0),dhs(0:n0),ddhs(0:n0)
      double precision hsp(0:n0),dhsp(0:n0),ddhsp(0:n0)
      double precision rsx,dhx,dhm0,rsm0,fac,sss,ehvl0,ehx
      double precision g0,g1,h0,h1
      double precision tol
      parameter (tol=1d-15)

C Heap:
      integer w(1)
      common /w/ w

      ndimx = ndim
      kmaxx = kmax
      lmaxg = ll(nlmg)
      lmaxh = ll(nlmh)

C --- Checks ---
      if (max(nlmh,nlmg) .gt. ndimx)
     .  call rxi('gh1c: nlm exceeds ndimx. nlm = ',max(nlmh,nlmg))
      if (max(lmaxh,lmaxg) .gt. n0)
     .  call rxi('gh1c: lmax exceeds n0. lmax = ',max(lmaxh,lmaxg))
c     if (ndim .gt. ndimx)
c    .  call rxi('gh1c: ndim is bigger than ndimx. ndim = ',ndim)  
c     if (nlmh .gt. ndimx)
c    .  call rxi('gh1c: nlmh is bigger than ndimx. nlmh = ',nlmh)  
c     if (nlmg .gt. ndimx)
c    .  call rxi('gh1c: nlmg is bigger than ndimx. nlmg = ',nlmh)  
c     if (kmax .gt. kmaxx)
c    .  call rxi('gh1c: kmax is bigger than n0. kmax = ',kmax)  

C --- setup ---
      if (ivl .eq. 0) then
        pmax = 2
      else
        pmax = 1
        jhd = 10*(ivl-1)+1
      endif 
c ... stop if there are any positive energies 
      if (ivl .ne. 0) then
        do il = 0, lmaxh
          if (ehvl(il) .gt. 0) 
     .    call rxi('gh1c: positive energy for l = ',il)  
        enddo
      endif

C ... Check if sites coinside
      sss = 0d0 
      do ii = 1, 3
        sss = sss + (ph(ii)-pg(ii))**2
      enddo
c ... ph \= pg => proceed with 1c expansion 
      if (sss .gt. tol) then
        
c ... make 1c coefficients
        call gxpml(ph,pg,rsmh,-rsmg,pmax,kmax,nlmh,nlmg,kmaxx,ndimx,
     .  pmaxx,cg,indxcg,jcg,cy,spg)
        if (ivl .ge. 1) then
          call hdxpml(jhd,ph,pg,ehvl,rsmh,-rsmg,pmax,kmax,nlmh,nlmg,
     .      kmaxx,ndimx,pmaxx-1,cg,indxcg,jcg,cy,shr,shr,shdr)
        endif 

c ... find values and Laplacians
        jvsl = 11
        call vsl(jvsl,rhc,rsmg,kmax,lmaxg,1,lmaxh,kmaxx,ndimx,
     .   pmaxx,spg,gvl,dums,gvl(1,1,1),pkl,gpkl)
        jvsl = 10
        if (ivl .eq. 0) then
          call vsl(jvsl,rhc,rsmg,kmax,lmaxg,pmax,lmaxh,kmaxx,ndimx,
     .     pmaxx,spg(0,1,1,1),hvl,dums,hvl(1,1,1),pkl,gpkl)
        elseif (ivl .eq. 1) then
          call vsl(jvsl,rhc,rsmg,kmax,lmaxg,pmax,lmaxh,kmaxx,ndimx,
     .     pmaxx-1,shr,hvl,dums,hvl(1,1,1),pkl,gpkl)
        elseif (ivl .eq. 2) then
          call vsl(jvsl,rhc,rsmg,kmax,lmaxg,pmax,lmaxh,kmaxx,ndimx,
     .     pmaxx-1,shdr,hvl,dums,hvl(1,1,1),pkl,gpkl)
        else
          call rxi('gh1c: ivl must be between 0 and 2. ivl = ',ivl)  
        endif

c ... ph = pg => do explicitly (only diagonal terms needed) 
      else
        
C ... Handle negative smoothing radii and rmt
        if (rsmh(0) .lt. 0d0) then
          call dvset(rsml(0),1,lmaxh+1,-rsmh(0))
        else
          call dcopy(lmaxh+1,rsmh(0),1,rsml(0),1)
        endif
        if (rhc(0) .lt. 0d0) then
          call dvset(rhcl(0),1,lmaxh+1,-rhc(0))
        else
          call dcopy(lmaxh+1,rhc(0),1,rhcl(0),1)
        endif
 
C make on-site Gaussians for each l
        rsx = -1d2
        dhx = -1d2
        do il = lmaxh, 0, -1
          dhm0 = rhcl(il)
          rsm0 = rsml(il)
          if (dabs(rsm0-rsx) + dabs(dhm0-dhx) .gt. tol) then
            call radgkl(dhm0,rsm0,pmax,il,pmaxx,gkl) 
C scale Gaussians by r^l
            if (il .ge. 1) then
              fac = 1.d0
              do ill = 1, il
                fac = fac*dhm0
                do ip = 0, pmax
                  gkl(ip,ill) = gkl(ip,ill)*fac 
                enddo
              enddo
            endif
          endif
          rsx = rsm0
          dhx = dhm0
        enddo
C if ivl > 0, also make on-site Hankels or H-dots 
        if (ivl. gt. 0) then
          rsx = -1d2
          dhx = -1d2
          ehx = 1d2
          do il = lmaxh, 0, -1
            dhm0 = rhcl(il)
            rsm0 = rsml(il)
            ehvl0 = ehvl(il)
            if (dabs(rsm0-rsx) + dabs(dhm0-dhx) + dabs(ehvl0-ehx) .gt. 
     .        tol) then
              call hanszd(jhd+1,dhm0,ehvl0,rsm0,il,
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
            dhx = dhm0
            ehx = ehvl0
          enddo
        endif

c fill in gvl and hvl
        call dpzero(gvl,ndim*ndim*2) 
        call dpzero(hvl,ndim*ndim*2) 
        ilmg = 0
        do il = 0, lmaxh
          g0 = gkl(0,il)
          g1 = gkl(1,il)
          if (ivl .eq. 0) then
            h0 = g1
            h1 = gkl(2,il)
          else
            h0 = hkl(0,il)
            h1 = hkl(1,il)
          endif
          do im = -il, il
            ilmg = ilmg + 1
            gvl(ilmg,ilmg,0) = g0 
            gvl(ilmg,ilmg,1) = g1 
            hvl(ilmg,ilmg,0) = h0 
            hvl(ilmg,ilmg,1) = h1 
          enddo 
        enddo 

c end the if-loop over sites 
      endif

      end

