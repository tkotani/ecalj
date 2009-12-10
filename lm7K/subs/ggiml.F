      subroutine ggiml(dp,rsmh,rsmg,nlmh,nlmg,pmax,kmax,
     .  ktop0,ndimg,ptop0,cg,indxcg,jcg,cy,s)
C- Integrals between energy-independent gaussians \int G_kL(r-Rg)*G_pM(r-Rh) d^3 r 
C ----------------------------------------------------------------------
Ci Inputs
Ci   dp    :Connecting vector between the sites at which the gaussians are centered       
Ci         :dp = Rg - Rh (see Remarks)  
Ci   rsmh  :vector of l-dependent smoothing radii of the gaussians to be expanded 
Ci         :rsmh must be specified for 0..ll(nlmh)
Ci   rsmg  :vector of l-dependent smoothing radii of the gaussians associated with 
Ci         :polynomials P_kL centred at the origin of the expansion 
Ci         :rsmg must be specified for 0..ll(nlmg)
Ci   nlmh  :L-cutoff for gaussians 
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   pmax  :max power of the laplace operator for gaussians to be expanded 
Ci   kmax  :polynomial cutoff
Ci   ktop0,ndimg,ptop0:leading dimensions of s
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   s     :integrals of gaussian; see Remarks
Cl Local variables
Cl   yl    :unnormalised sph. harmonic polynomials
Cr Remarks
Cr   s(k,L,p,M) = \int G_kL(e=0,rsg;r-Rg) * G_pM(e=0,rsh;r-Rh) d^3 r 
Cr              = (-)^l \sum_N C_LMN G_qN(e=0,rs;Rg-Rh), 
Cr   where q = k + p + (l+m-n)/2; rs = sqrt(rsg^2 + rsh^2), 
Cr   Rg = pg (centre of expansion) and 
Cr   Rh = ph (centre of the Gaussian to be expanded)  
Cr   See MvS's notes on smooth hankels 
Cu Updates
Cu   30 Aug 06  adapted from ghibl.f (S. Lozovoi)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlmg,nlmh,pmax,kmax,ktop0,ndimg,ptop0
      integer jcg(1),indxcg(1)
      double precision rsmg(0:1),rsmh(0:1)
      double precision dp(3),cg(1),cy(1)
      double precision s(0:ktop0,ndimg,0:ptop0,nlmh)
C ... Local parameters
      integer nlm0,ktopm,lmx0
      integer icg,icg1,icg2,ii,ilg,ilh,ilm,indx,mh,mg,pkmax,
     . k,ip,ik,lg,lh,ll,lm,lmaxg,lmaxh,lmxx,ione,ione0
      parameter (lmx0=12, nlm0=(lmx0+1)**2, ktopm=31)
      double precision cc,rmxx,rmx2,dp1,dp2
      double precision rmh,rmh2,rmg,rmg2
      double precision gkl(0:ktopm,0:lmx0),yl(nlm0)
      double precision tol
      data tol/1.d-16/

      if (nlmh .eq. 0 .or. nlmg .eq. 0) return

C ... rsm-independent setup
      lmaxh = ll(nlmh)
      lmaxg = ll(nlmg)
      pkmax = pmax + kmax 
      lmxx = lmaxg + lmaxh
      if (lmxx .gt. lmx0) call rxi('ggiml: increase lmx0 to',lmxx)
      if (max0(lmaxg,lmaxh)+pkmax .gt. ktopm) 
     . call rxi('ggiml: increase ktopm to',max0(lmaxg,lmaxh)+pkmax)

C ... make sph. harmonics for connecting vector dp
      call sylm(dp,yl,lmxx,dp2)
      dp1 = dsqrt(dp2)
      do ilm = 1, (lmxx+1)**2
        yl(ilm) = yl(ilm)*cy(ilm)
      enddo
 
C ... Initialise s 
      call dpzero(s(0,1,0,1),(ktop0+1)*ndimg*(ptop0+1)*nlmh)

C --- Loops over lh and lg  ---
      ione0 = (-1)**lmaxg
      rmxx = -1.d2
      do lh = lmaxh, 0, -1
        rmh = rsmh(lh)
        rmh2 = rmh*rmh
        ione = ione0
        do lg =  lmaxg, 0, -1
          rmg = rsmg(lg)
          rmg2 = rmg*rmg
          rmx2 = rmg2 + rmh2
c ... if total rsm changed, update radial gaussians
          if (dabs(rmx2-rmxx) .gt. tol) then
            rmxx = rmx2
            call radgkl(dp1,dsqrt(rmx2),max0(lh,lg)+pkmax,
     .                  lh+lg,ktopm,gkl)
          endif
C ... Combine with Clebsch-Gordan coefficients
          ilh = lh*lh
          do mh = -lh, lh
            ilh = ilh + 1
            ilg = lg*lg
            do mg = -lg, lg
              ilg = ilg + 1
              ii = max0(ilg,ilh)
              indx = (ii*(ii-1))/2 + min0(ilg,ilh)
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1)-1
              do  icg = icg1, icg2
                ilm = jcg(icg)
                lm = ll(ilm)
                k = (lg+lh-lm)/2
                cc = ione*cg(icg)*yl(ilm)                                      
                do  ip = 0, pmax
                  do  ik = 0, kmax
                    s(ik,ilg,ip,ilh) = s(ik,ilg,ip,ilh) + 
     .                 gkl(k+ik+ip,lm)*cc
                  enddo 
                enddo 
              enddo
            enddo
          enddo
          ione = -ione
        enddo
      enddo

      end


