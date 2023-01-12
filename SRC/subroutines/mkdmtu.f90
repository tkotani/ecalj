subroutine mkdmtu(isp,iq,qp,nev,evec,dmatu)
  use m_lmfinit,only: ispec,sspec=>v_sspec,nbas,nlmax,nsp,nspc,nl,n0,nppn,nlibu,lmaxu,nlibu,lldau,idu
  use m_mkpot,only: ppnl=>ppnl_rv
  use m_subzi, only: wtkb=>rv_a_owtkb
  use m_igv2x,only: ndimh
  use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx
  use m_makusq,only: makusq
  use m_locpot,only: rotp
  !- Calculate density matrix for LDA+U channels
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
  !i     Stored: *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxa idu
  !i     Stored: *
  !i     Passed to: *
  !i   wtkb  :eigenvalue weights for BZ integration of occupied states
  !i   isp   :current spin channel (1 or 2)
  !i   iq    :qp index, used only to address element in wtkb
  !i         :NB: aus is stored only for current qp
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   ndham :dimensions wtkb,aus
  !i   nlmax :1st dimension of aus (maximum nlma over all sites)
  !i   nbas  :size of basis
  !i   nev   :actual number of eigenvectors generated
  !i   ppnl  :nmto-like pot pars
  !i   aus   :coefficients to phi and phidot made previously by makusqldau
  !i  lldau  :lldau(ib)=0 => no U on this site otherwise
  !i         :U on site ib with dmat beginning at dmats(*,lldau(ib))
  !o Outputs
  !o   dmatu :density matrix for specified LDA+U channels
  !b Bugs
  !b   Never checked for noncollinear case
  !u Updates
  !u   09 Nov 05 Convert dmat to complex form
  !u   28 Jun 05 bug fix for nspc=2
  !u   09 Jun 05 (MvS) extended to local orbitals
  !u   30 Apr 05 (WRL) first created
  ! ----------------------------------------------------------------------
  implicit none
  integer :: isp,iq,nev
  double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
  double complex add,au,as,az,ap1,ap2
  double precision :: dlphi,rmt,dlphip,phi,phip,dphi,dphip,r(2,2),det,phz,dphz
  integer :: lmxa,ilm1,ilm2,l,iv,m1,m2,ib,is,igetss,iblu,ispc, ksp
  complex(8),allocatable ::aus(:,:,:,:,:)
  double complex evec(ndimh,nsp,nev)
  real(8)::qp(3)
  complex(8):: auas(2)
  allocate(aus(nlmax,ndham*nspc,3,nsp,nbas))! aus(2*nlmax*ndhamx*3*nsp*nbas))
  aus=0d0
  call makusq(nbas,[0] , nev,  isp, 1, qp, evec, aus )
  iblu = 0
  do  ib = 1, nbas
     if (lldau(ib) == 0) goto 10
     is = ispec(ib) !ssite(ib)%spec
     lmxa=sspec(is)%lmxa
     !idu= sspec(is)%idu
     rmt= sspec(is)%rmt
     do  l = 0, min(lmxa,3)
        if (idu(l+1,is) /= 0) then
           iblu = iblu+1
           !           In noncollinear case, isp=1 always => need internal ispc=1..2
           !           ksp is the current spin index in both cases:
           !           ksp = isp  in the collinear case
           !               = ispc in the noncollinear case
           !           ispc=1 for independent spins, and spin index when nspc=2
           do  ispc = 1, nspc
              ksp = max(ispc,isp)
              !             For rotation (u,s) to (phi,phidot)
!              dlphi  = ppnl(3,l+1,ksp,ib)/rmt
!              dlphip = ppnl(4,l+1,ksp,ib)/rmt
!              phi    = ppnl(5,l+1,ksp,ib)
!              phip   = ppnl(6,l+1,ksp,ib)
!              dphi   = phi*dlphi/rmt
!              dphip  = dlphip/rmt*phip
!              det    = phi*dphip - dphi*phip
!              r(1,1) = dphip/det
!              r(1,2) = -dphi/det
!              r(2,1) = -phip/det
!              r(2,2) = phi/det
              !             For projection from loc. orb. onto (u,s)
              phz    = ppnl(11,l+1,ksp,ib)
              dphz   = ppnl(12,l+1,ksp,ib)
              ilm1 = l*l
              do  m1 = -l, l
                 ilm1 = ilm1+1
                 ilm2 = l*l
                 do  m2 = -l, l
                    ilm2 = ilm2+1
                    add = (0d0,0d0)
         !  Since (au,as,az) are coefficients to (u,s,gz), (gz is local orbital with val=slo=0 at MT)
         !  Local orbital contribution adds to u,s
         !  deltau = -phi(rmax) * az   deltas = -dphi(rmax) * az
                    do  iv = 1, nev
                       az = aus(ilm1,iv,3,ksp,ib)
                       au = aus(ilm1,iv,1,ksp,ib) - phz*az
                       as = aus(ilm1,iv,2,ksp,ib) - dphz*az !u,s components
                       auas= matmul([au,as],rotp(l,ksp,:,:,ib))
                       ap1 = auas(1) !au*r(1,1) + as*r(2,1) !projection to phi components.
                       az = aus(ilm2,iv,3,ksp,ib)
                       au = aus(ilm2,iv,1,ksp,ib) - phz*az
                       as = aus(ilm2,iv,2,ksp,ib) - dphz*az
                       auas= matmul([au,as],rotp(l,ksp,:,:,ib))
                       ap2 = auas(1) !au*r(1,1) + as*r(2,1)
                       add = add + ap1*dconjg(ap2)*wtkb(iv,isp,iq)
                    enddo
                    dmatu(m1,m2,ksp,iblu) = dmatu(m1,m2,ksp,iblu) + add
                 enddo
              enddo
           enddo
        endif
     enddo
10   continue
  enddo
  deallocate(aus)
end subroutine mkdmtu
