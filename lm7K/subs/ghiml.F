      subroutine ghiml(job,dp,elh,rsmh,rsmg,nlmh,nlmg,pmax,kmax,
     .  ktop0,ndimg,ptop0,cg,indxcg,jcg,cy,sr,si)
C- Integrals between energy-independent gaussians and smoothed hankels 
C  centered at Rg and Rh, respectively:  \int G_kL(r-Rg)*H_pM(e,r-Rh) d^3 r 
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :if job=0, returns integrals of gaussians * smoothed hankels 
Ci         :if job=1, returns only the real part of the above integrals
Ci         :if job=2, returns integrals of gaussians * difference between 
Ci         :smoothed and bare hankels 
Ci   dp    :Connecting vector between the sites: dp = Rg - Rh (see Remarks)  
Ci   elh,rsmh:vectors of l-dependent energies and smoothing radii of 
Ci         :smoothed Hankels; must be specified for 0..ll(nlmh)
Ci   rsmg  :vector of l-dependent smoothing radii of the gaussians associated  
Ci         :with polynomials P_kL centred at the origin of the expansion 
Ci         :rsmg must be specified for 0..ll(nlmg)
Ci   nlmh  :L-cutoff for the hankels 
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   pmax  :max power of the laplace operator for hankels to be expanded 
Ci   kmax  :polynomial cutoff
Ci   ktop0,ndimg,ptop0:leading dimensions of sr and si
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   sr,si :real and imaginary part of s, s = sr + i*si (see Remarks)
Co         :for negative energies si set to 0  
Co         :if job>1, si is not referenced 
Cl Local variables
Cl   yl    :unnormalised sph. harmonic polynomials
Cr Remarks
Cr   1. Integrals are defined as:
Cr      s(k,L,p,M) = \int G_kL(e=0,rsg;r-Rg) * H_pM(eh,rsh;r-Rh) d^3 r 
Cr              = (-)^l exp(-rsg^2*eh/4) \sum_N C_LMN H_qN(eh,rs;Rg-Rh),  (*)
Cr   where q = k + p + (l+m-n)/2; rs = sqrt(rsg^2 + rsh^2), 
Cr   Rg = pg (centre of expansion) and 
Cr   Rh = ph (centre of the sm. Hankel to be expanded)  
Cr   (See MvS's notes on smooth hankels) 
Cr   2. For the difference between smoothed and bare hankels (job=2), 
Cr   H_qN(eh,rs;Rg-Rh) in (*) are replaced with: 
Cr   [H_qN(eh,rs;Rg-Rh) - H_qN(eh,rsg;Rg-Rh)] 
Cb Bugs
Cu Updates
Cu   16 Jul 07 (S. Lozovoi) positive energies, laplacians, and
Cu             (smooth - bare) hankel options included 
Cu   23 Jan 07 (S. Lozovoi) Adapted from ghibl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlmg,nlmh,pmax,kmax,ktop0,ndimg,ptop0
      integer job,jcg(1),indxcg(1)
      double precision rsmg(0:1),rsmh(0:1),elh(0:1)
      double precision dp(3),cg(1),cy(1)
      double precision sr(0:ktop0,ndimg,0:ptop0,nlmh)
      double precision si(0:ktop0,ndimg,0:ptop0,nlmh)
C ... Local parameters
      integer nlm0,ktopm,lmx0
      integer icg,icg1,icg2,ii,ilg,ilh,ilm,indx,mh,mg,
     . k,ip,ik,ktop,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,ione
      parameter (lmx0=12, nlm0=(lmx0+1)**2, ktopm=31)
      logical lpos
      double precision cc,chg
      double precision rmxx,rmx1,rmx2,dp1,dp2,eh,ehx,beta
      double precision rmg,rmg2,rmgx,rsmv(0:lmx0)
      double precision gkl(0:ktopm,0:lmx0),yl(nlm0)
      double precision xi(0:lmx0,0:ktopm),xi0(0:lmx0,0:ktopm)
      double precision fi(0:lmx0,0:ktopm)
      double precision fpi,tol
      data tol/1.d-16/

      if (nlmh .eq. 0 .or. nlmg .eq. 0) return

C ... rsm-independent setup
      fpi = 16d0*datan(1d0)
      lmaxh = ll(nlmh)
      lmaxg = ll(nlmg)
      lmaxx = lmaxg+lmaxh
      ktop = max0(lmaxg,lmaxh) + pmax + kmax
      if (lmaxx .gt. lmx0) call rxi('ghiml: increase lmx0 to',lmaxx)
      if (ktop .gt. ktopm) call rxi('ghiml: increase ktopm to',ktop)

C ... Make sph. harmonics for connecting vector dp

      call sylm(dp,yl,lmaxx,dp2)
      dp1 = dsqrt(dp2)
      do ilm = 1, (lmaxx+1)**2
        yl(ilm) = yl(ilm)*cy(ilm)
      enddo

 
C ... Initialise s 

      call dpzero(sr(0,1,0,1),(ktop0+1)*ndimg*(ptop0+1)*nlmh)
      if(job .eq. 0) 
     .  call dpzero(si(0,1,0,1),(ktop0+1)*ndimg*(ptop0+1)*nlmh)

C --- Outer loops over lh and lg  ---
      rmgx = -1.d2
      rmxx = -1.d2
      ione = -1
      do lg =  0, lmaxg
        ione = -ione
        rmg = rsmg(lg)
        rmg2 = rmg*rmg
        ehx = 1.d20
        do lh = 0, lmaxh
          eh = elh(lh)
          lpos = (eh .gt. 0d0) 
          rmx2 = rmg2 + rsmh(lh)**2
c     ... if total rsm or energy changed, update smoothed hankels
          if (dabs(eh-ehx)+dabs(rmx2-rmxx) .gt. tol) then
            ehx = eh
            rmxx = rmx2
            beta = fpi*dexp(0.25d0*eh*rmx2)
            chg = dexp(-0.25d0*eh*rmg2) 
            rmx1 = dsqrt(rmx2)
            rsmv(0) = -rmx1
            call hansmz(dp1,eh,rsmv,xi(0,0),fi(0,0),lmaxx) 
c     ...   make laplacians H_k+1,L = -e*H_kL - 4*pi*G_kL 
            if(ktop .gt. 0) then
              call radgkl(dp1,rmx1,ktop-1,lmaxx,ktopm,gkl)
              do ik = 1, ktop
                do ii = 0, lmaxx
                  xi(ii,ik) = -eh*xi(ii,ik-1) - beta*gkl(ik-1,ii)
                enddo
              enddo
              if (job .eq. 0 .and. lpos) then
                do ik = 1, ktop
                  do ii = 0, lmaxx
                    fi(ii,ik) = -eh*fi(ii,ik-1) 
                  enddo
                enddo
              endif
            endif
c     ...   if job=2, update bare hankels if necessary and replace xi --> (xi-xi0)
            if (job .eq. 2) then 
              if (dabs(rmg-rmgx) .gt. tol) then 
                rmgx = rmg
                rsmv(0) = -rmg
                call hansmz(dp1,eh,rsmv,xi0,fi,lmaxx) 
                if(ktop .gt. 0) then
                  call radgkl(dp1,rmg,ktop-1,lmaxx,ktopm,gkl)
                  do ik = 1, ktop
                    do ii = 0, lmaxx
                      xi0(ii,ik) = -eh*xi0(ii,ik-1) - beta*gkl(ik-1,ii)
                    enddo
                  enddo
                endif
              endif
              call daxpy((lmaxx+1)*ktop,-1d0,xi0,1,xi,1)
            endif 
          endif
C     ... Combine with Clebsch-Gordan coefficients
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
                cc = ione*chg*cg(icg)*yl(ilm)                                      
                do  ip = 0, pmax
                  do  ik = 0, kmax
                    sr(ik,ilg,ip,ilh) = sr(ik,ilg,ip,ilh) + 
     .                                  xi(lm,k+ik+ip)*cc
                  enddo 
                enddo 
                if (job .eq. 0 .and. lpos) then
                  do  ip = 0, pmax
                    do  ik = 0, kmax
                      si(ik,ilg,ip,ilh) = si(ik,ilg,ip,ilh) + 
     .                                    fi(lm,k+ik+ip)*cc
                    enddo 
                  enddo 
                endif 
              enddo
            enddo
          enddo
        enddo
      enddo

      end


