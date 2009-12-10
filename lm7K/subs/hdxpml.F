      subroutine hdxpml(job,ph,pg,eh,rsmh,rsmg,pmax,kmax,nlmh,nlmg,
     .  k0,ndim,p0,cg,indxcg,jcg,cy,sr,si,sdr)
C- Coefficients to expand smoothed hankels H_pL centered at ph
C  into a sum of polynomials P_kL centered at pg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :if 1s digit = 0, returns real and imaginary parts of the expansion 
Ci         :if 1s digit = 1, returns only the real part of the expansion 
Ci         :if 1s digit = 2, returns the expansion of the difference between
Ci         :smoothed and bare hankels (always real)
Ci         :if 10s digit > 0, same as above but in addition expands energy
Ci                            derivative of sm. Hankels
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
Co         :if job 1s digit = 0, si = 0 for negative energies 
Co         :if job 1s digit > 0, si is not referenced 
Co   sdr   :real part of s-dot(k,M,p,L), ie the coefficients of 
Co         :H-dot expansion into P_kL
Cr Remarks
Cr   Expansion is:  H_pL(r-ph) = sum(kM) s(k,M,p,L) * P_kM(r-pg)
Cr   for p=0 case, see J. Math. Phys. {\bf 39}, 3393 (1998), Eq. 12.17
Cr
Cr   Same form for H-dot:
Cr          H-dot_pL(r-ph) = sum(kM) s-dot(k,M,p,L) * P_kM(r-pg)
Cr
Cr   As rsmg -> 0, the expansion turns into a Taylor series of H_L.
Cr   As rsmg increases the error is spread out over a progressively
Cr   larger radius, thus reducing the error for larger r while
Cr   reducing accuracy at small r.
Cb Bugs
Cb   Does not handle the case pg = ph
Cb   Does not handle zero smoothing case rsmh(l) or rsmg(l) = 0 
Cb   H-dot option is currently implemented only for job = 11
Cu Updates
Cu   23 Jan 07 (S. Lozovoi) Adapted from hxpbl.f
Cu   17 Jul 07 (S. Lozovoi) Updated to accomodate changes in ghiml.f 
Cu   15 May 08 (S. Lozovoi) H-dot added
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,p0,k0,pmax,kmax,ndim,nlmg,nlmh,jcg(*),indxcg(*)
      double precision eh(0:*),rsmh(0:*),rsmg(0:*)
      double precision ph(3),pg(3),cg(*),cy(*)
      double precision sr(0:k0,ndim,0:p0,*),si(0:k0,ndim,0:p0,*)
      double precision sdr(0:k0,ndim,0:p0,*)
C ... Local parameters
      logical ldot
      integer n0
      parameter (n0=10)
      integer mode0,ilmg,ilmh,k,p,il,ll,lmaxg,lmaxh,m,nm,i
      double precision dfact,fac,factk,fpi,rs,rs2,facrs
      double precision rsmhl(0:n0), rsmgl(0:n0),dp(3)

      mode0 = mod(job,10)
      ldot = (mod(job/10,10) .gt. 0)

      if (ldot .and. (mode0 .ne. 1)) call rx('hdxpml: H-dot option'//
     .  ' is currently implemented only with mod(job,10) = 1')

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
      call ghdiml(job,dp,eh,rsmhl,rsmgl,nlmh,nlmg,pmax,kmax,
     .  k0,ndim,p0,cg,indxcg,jcg,cy,sr,si,sdr)

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
c               if (ldot) sdr(k,ilmg,p,ilmh) = sdr(k,ilmg,p,ilmh)*fac
              enddo
            enddo
c ...       make H-dot if required
            if (ldot) then
              do ilmh = 1, nlmh
                do p = 0, pmax
                  sdr(k,ilmg,p,ilmh) = sdr(k,ilmg,p,ilmh)*fac
                enddo
              enddo
            endif
C ...       Imaginary part is required only if job=0 and e positive 
            if (mode0 .eq. 0) then
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

