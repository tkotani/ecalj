subroutine hklgbl(p,rsm,e,q,kmax,nlm,k0,nlm0,cy,hkl,ghkl)!slat,
  !- Bloch-sums of k,L-dependent smooth hankel functions and gradients
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for hkl
  !i   k0    :leading dimension of hkl
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   hkl   :Bloch-summed smoothed Hankels
  !o   ghkl  :gradients of hkl
  !r Remarks
  !r   H_kL = laplace^k H_L
  !r   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
  !r   H_kL are made to kmax+1, lmax+1 in order to assemble grads.
  !u Updates
  !u   25 May 00 Adapted from nfp hklg_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nlm,nlm0
  double precision :: q(3),p(3),cy(1),e,rsm!slat(1),
  double complex hkl(0:k0,nlm0),ghkl(0:k0,nlm0,3)
  ! ... Local parameters
  integer :: ilm,k,kx1,kx2,ky1,ky2,kz,ll,lmax,m,nlm1
  double precision :: cx1,cx2,cy1,cy2,cz

  if (nlm == 0) return
  lmax = ll(nlm)
  nlm1 = (lmax+2)**2
  if (nlm1 > nlm0) call rxi('hklgbl: need nlm0 ge ',nlm1)
  if (kmax+1 > k0) call rxi('hklgbl: need k0 ge ',kmax+1)
  do  m = 1, 3
     do  ilm = 1, nlm
        do  k = 0, kmax
           ghkl(k,ilm,m) = (0d0,0d0)
        enddo
     enddo
  enddo

  ! ... Make Hkl's up to one higher in l and k
  call hklbl(p,rsm,e,q,kmax+1,nlm1,k0,cy,hkl)!,slat

  ! ... Assemble gradients using Clebsh-Gordans for p functions
  do  ilm = 1, nlm
     call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
     do  k = 0, kmax
        ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl(k,kx1)-cx2*hkl(k,kx2)
        ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl(k,ky1)-cy2*hkl(k,ky2)
        ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl(k,kz)
        if (ilm <= lmax*lmax) then
           ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl(k+1,ilm)
           ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl(k+1,ilm)
           ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl(k+1,ilm)
           ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl(k+1,ilm)
           ghkl(k,kz,3) =  ghkl(k,kz,3)  - cz *hkl(k+1,ilm)
        endif
     enddo
  enddo

end subroutine hklgbl

