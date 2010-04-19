      subroutine hxpml(job,ph,pg,eh,rsmh,rsmg,pmax,kmax,nlmh,nlmg,
     .  k0,ndim,p0,cg,indxcg,jcg,cy,sr,si)
C- Coefficients to expand smoothed hankels H_pL centered at ph
C  into a sum of polynomials P_kL centered at pg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :if job=0, returns real and imaginary parts of the expansion 
Ci         :if job=1, returns only the real part of the expansion 
Ci         :if job=2, returns the expansion of the difference between
Ci         :smoothed and bare hankels (always real)
Ci   ph    :Function is centered at ph; see Remarks
Ci   pg    :origin of the expansion; see Remarks
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 0..ll(nlmh)
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed Hankel
Ci         :EITHER must be specified for 0..ll(nlmh)
Ci         :OR rsmh(0) = const <0. Implies rsmh(l) = -const for all l
Ci   rsmg  :vector of l-dependent smoothing radii for polynomials
Ci         :same convention as for rsmh
Ci   pmax  :max power of the Laplace operator (lap)^p H_L := H_pL
Ci   kmax  :polynomial cutoff
Ci   nlmh  :L-cutoff for smoothed Hankel functions being expanded
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   k0,ndim,p0:leading dimensions of coefficient arrays sr and si
Ci   cg    :Clebsch Gordan coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordan coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Co Outputs
Co   sr,si :real and imaginary parts of cofficients s(k,M,p,L); see Remarks
Co         :here k=0..kmax, M=1..nlmg, p=0..pmax, L=1..nlmh
Co         :if job=0, si = 0 for negative energies 
Co         :if job>0, si is not referenced 
Cr Remarks
Cr   Expansion is:  H_pL(r-ph) = sum(kM) s(k,M,p,L) * P_kM(r-pg)
Cr   for p=0 case, see J. Math. Phys. {\bf 39}, 3393 (1998), Eq. 12.17
Cr
Cr   As rsmg -> 0, expansion turns into a Taylor series of H_L.
Cr   As rsmg increases the error is spread out over a progressively
Cr   larger radius, thus reducing the error for larger r while
Cr   reducing accuracy at small r.
Cb Bugs
Cb   Does not handle the case pg = ph
Cb   Does not handle zero smoothing case rsmh(l) or rsmg(l) = 0 
Cu Updates
Cu   23 Jan 07 Adapted from hxpbl.f
Cu   17 Jul 07 Updated to accomodate changes in ghiml.f (S. Lozovoi)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,p0,k0,pmax,kmax,ndim,nlmg,nlmh,jcg(1),indxcg(1)
      double precision eh(0:1),rsmh(0:1),rsmg(0:1)
      double precision ph(3),pg(3),cg(1),cy(1)
      double precision sr(0:k0,ndim,0:p0,1),si(0:k0,ndim,0:p0,1)
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

C ... Integrals of gaussians * sm. hankels
      call ghiml(job,dp,eh,rsmhl,rsmgl,nlmh,nlmg,pmax,kmax,
     .  k0,ndim,p0,cg,indxcg,jcg,cy,sr,si)

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
                sr(k,ilmg,p,ilmh) = sr(k,ilmg,p,ilmh)*fac
              enddo
            enddo
C ...       Imaginary part is required only if job=0 and e positive 
            if (job .eq. 0) then
              do ilmh = 1, nlmh
                if (eh(ll(ilmh)) .gt. 0) then
                  do p = 0, pmax
                    si(k,ilmg,p,ilmh) = si(k,ilmg,p,ilmh)*fac
                  enddo
                endif
              enddo
            endif
            factk = 4*(k+1)*factk
            facrs = rs2*facrs
          enddo
        enddo
      enddo

      end

