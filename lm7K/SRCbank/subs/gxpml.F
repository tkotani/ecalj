      subroutine gxpml(ph,pg,rsmh,rsmg,pmax,kmax,nlmh,nlmg,k0,ndim,p0,
     .  cg,indxcg,jcg,cy,s)
C- Coefficients to expand generalised energy-independent Gaussians G_pL
C  centred at ph into a sum of polynomials P_kL centred at pg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ph    :function to expand is centred at ph; see Remarks
Ci   pg    :origin of the expansion; see Remarks
Ci   rsmh  :vector of l-dependent smoothing radii of Gaussians to expand
Ci         :EITHER must be specified for 0..ll(nlmh)
Ci         :OR     rsmh(0)<0.  Implies rsmh(l) = -const for all l
Ci   rsmg  :vector of l-dependent smoothing radii for polynomials
Ci         :same convention as for rsmh
Ci   kmax  :polynomial cutoff
Ci   pmax  :max power of the Laplace operator (lap)^p G_L := G_pL
Ci   nlmh  :L-cutoff for Gaussians being expanded
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   k0    :leading dimension of coefficient array s
Ci   ndim  :second dimension of coefficient array s
Ci   p0    :third dimension of coefficient array s
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalisation constants for spherical harmonics
Co Outputs
Co   s     :coefficients s(k,M,p,L); see Remarks
Co         :here k=0..kmax, p=0..pmax, M=1..nlmg, L=1..nlmh
Cr Remarks
Cr   Expansion is:  G_pL(r-ph) = sum(kM) s(k,M,p,L) * P_kM(r-pg)
Cr   See Sec.XII in J. Math. Phys. {\bf 39},3393 (1998) and MvS's notes
Cr   on smooth Hankels
Cr
Cb Bugs
Cb   Does not handle the case pg = ph
Cu Updates
Cu   30 Aug 06 adapted from fp/hxpbl.f (S. Lozovoi)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer p0,k0,pmax,kmax,ndim,nlmg,nlmh,jcg(*),indxcg(*)
      double precision rsmh(0:*),rsmg(0:*)
      double precision ph(3),pg(3),cg(*),cy(*)
      double precision s(0:k0,ndim,0:p0,*)
C ... Local parameters
      integer n0
      integer ilmg,ilmh,k,p,il,ll,lmaxg,lmaxh,m,nm,i
      parameter (n0=10)
      double precision dfact,fac,factk,fpi,rs,rs2,facrs
      double precision rsmhl(0:n0), rsmgl(0:n0),dp(3)

      if (nlmg .eq. 0 .or. nlmh .eq. 0) return
      fpi = 16d0*datan(1d0)

C ... Handle negative smoothing radii
      lmaxh = ll(nlmh)
      if (rsmh(0) .lt. 0d0) then
        call dvset(rsmhl(0),1,lmaxh+1,-rsmh(0))
      else
        call dcopy(lmaxh+1,rsmh(0),1,rsmhl(0),1)
      endif
      lmaxg = ll(nlmg)
      if (rsmg(0) .lt. 0d0) then
        call dvset(rsmgl(0),1,lmaxg+1,-rsmg(0))
      else
        call dcopy(lmaxg+1,rsmg(0),1,rsmgl(0),1)
      endif

C ... Only the connecting vector matters
      do i = 1, 3
        dp(i) = pg(i) - ph(i)
      enddo

C ... Integrals of two gaussians
      call ggiml(dp,rsmhl,rsmgl,nlmh,nlmg,pmax,kmax,
     .  k0,ndim,p0,cg,indxcg,jcg,cy,s)


C ... Scale to get coefficients of the P_kL
      ilmg = 0
      dfact = 1d0
      do il = 0, lmaxg
        nm = 2*il+1
        dfact = dfact*nm
        rs = rsmgl(il)
        rs2 = rs*rs
        do m = 1, nm
          ilmg = ilmg+1
          factk = 1d0
          facrs = rs**il
          do k = 0, kmax
            fac = fpi*facrs / (factk*dfact)
            do ilmh = 1, nlmh
              do p = 0, pmax
                s(k,ilmg,p,ilmh) = s(k,ilmg,p,ilmh)*fac
              enddo
            enddo
            factk = 4*(k+1)*factk
            facrs = rs2*facrs
          enddo
        enddo
      enddo

      end

