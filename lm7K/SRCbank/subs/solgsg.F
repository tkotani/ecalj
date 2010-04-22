      subroutine solgsg(px,rsm,lmx,pmax,cy,p0,gex,ggrad)
C- Solid generalized Gaussians  G_pL = (lap)^p G_L and the gradients \nabla G_L
C  at a given point X
C ----------------------------------------------------------------------
Ci Inputs
Ci   px    :coordinates of X relative to the centre of the Gaussian
Ci   rsm   :vector of l-dependent smoothing radii of the Gaussian
Ci         :must be specified for 0..ll(nlm)
Ci   lmax  :max angular momentum
Ci   pmax  :max power of the Laplace operator
Ci   cy    :Normalisation constants for spherical harmonics
Ci   p0    :leading dimension of gex
Co Outputs
Co   gex   : G_pL(X), p = 0,...,pmax, L = 1,...,(lmax+1)^2
Co   ggrad : \nabla G_L(X), L = 1,...,(lmax+1)^2
Cr Remarks
Cb Bugs
Cu Updates
Cu   23 Jul 08 handles px=0 (S. Lozovoi)
Cu   22 Jun 07 gradents added (S. Lozovoi)
Cu   31 Aug 06 first written (S. Lozovoi)
C ----------------------------------------------------------------------
C     implicit none
C Input parameters
      double precision px(3),rsm(0:*),cy(*)
      integer lmx,pmax,p0
C Output parameters
      double precision gex(0:p0,*),ggrad(3,*)
C Local variables
      integer n0
      parameter (n0 = 10)
      integer ll,ilm,ilr,ip,il,nlmp,lmax
      double precision rsml(0:n0),rsm0,rsx
      double precision gkl(0:n0,0:n0)
      double precision gexx(0:p0,(n0+1)**2)
      double precision yl((n0+1)**2)
      double precision r1
      integer kz,kx1,kx2,ky1,ky2
      double precision cz,cx1,cx2,cy1,cy2
      double precision tol
      data tol/1.d-14/

      if(lmx. ge. n0 .or. pmax .gt. n0)
     .  call rxi(' solgsg: Increase n0 up to at least ',max0(lmx,pmax))
      if(pmax. lt. 1) call
     .  rxi(' solgsg: pmax = %i should be at least 1 for gradients',pmax)

      lmax = lmx
      call sylm(px,yl,lmax+1,r1)
      r1 = dsqrt(r1)

c ... if r1 = 0, only l=0 term survives in G_pL and
c                     l=1 terms survive in \grad G_L
      if(dabs(r1) .lt. tol .and. lmx .gt. 1) then
        lmax = 1
        call dpzero(gex, (p0+1)*(lmx+1)**2)
        call dpzero(ggrad, 3*(lmx+1)**2)
      endif

c ... need lmax+1 for gradients
      nlmp = (lmax+2)**2
      do  ilm = 1, nlmp
        yl(ilm) = yl(ilm)*cy(ilm)
      enddo

C ... Handle negative smoothing radii
      if (rsm(0) .lt. 0d0) then
        call dvset(rsml(0),1,lmax+1,-rsm(0))
      else
        call dcopy(lmax+1,rsm(0),1,rsml(0),1)
      endif

c ... start big loop over smoothing radii

      rsx = -1d2
      do ilr = lmax, 0, -1
        rsm0 = rsml(ilr)
        if (dabs(rsm0-rsx) .gt. tol) then

          call radgkl(r1,rsm0,pmax,lmax+1,n0,gkl)
          rsx = rsm0

c ... make solid Gaussians and their laplacians for given rsm0 up to lmax+1
          do  ilm = 1, nlmp
            il = ll(ilm)
            do ip = 0, pmax
              gexx(ip,ilm) = gkl(ip,il)*yl(ilm)
            enddo
          enddo

c ... make gradients \nabla G_L up to current lmax
          call dpzero(ggrad, 3*(ilr+1)**2)

          do  ilm = 1, (ilr+1)**2
            call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
            ggrad(1,ilm) = ggrad(1,ilm) - cx1*gexx(0,kx1)
     .                                  - cx2*gexx(0,kx2)
            ggrad(2,ilm) = ggrad(2,ilm) - cy1*gexx(0,ky1)
     .                                  - cy2*gexx(0,ky2)
            ggrad(3,ilm) = ggrad(3,ilm) - cz*gexx(0,kz)

            if (ilm .le. ilr*ilr) then
              ggrad(1,kx1) = ggrad(1,kx1) - cx1*gexx(1,ilm)
              ggrad(1,kx2) = ggrad(1,kx2) - cx2*gexx(1,ilm)
              ggrad(2,ky1) = ggrad(2,ky1) - cy1*gexx(1,ilm)
              ggrad(2,ky2) = ggrad(2,ky2) - cy2*gexx(1,ilm)
              ggrad(3,kz)  = ggrad(3,kz)  - cz *gexx(1,ilm)
            endif
          enddo

c ... save G_pL up to current lmax only
          call dcopy((p0+1)*(ilr+1)**2,gexx(0,1),1,gex(0,1),1)

        endif
      enddo

      end

