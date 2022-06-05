subroutine rstr0(nxi,lxi,exi,nlmx,np,x,y,z,lmxa,iop,hl,hd)
  use m_ropyln,only: ropyln
  !- Reduced strux for a vector of energies, standard definition IV-43.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nxi   :number of energies for each connecting vector
  !i   lxi   :lmax for each energy and point
  !i   exi   :list of nxi energies
  !i   nlmx  :leading dimension of hl,hd
  !i   np    :number of connecting vectors
  !i   x,y,z :connecting vectors
  !i   lmxa  :Make hl for lxi + lmxa
  !i   iop   :1: calculate h only
  !i          any other number: calculate both h and hdot
  !o Outputs
  !o   hl    :reduced strux for each energy, point, to lxi + lmxa
  !o   hd    :reduced energy derivative
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: np,nxi,lxi(nxi,np),lmxa,nlmx,iop
  double precision :: hl(nlmx,nxi,np),hd(nlmx,nxi,np), &
       exi(nxi),x(np),y(np),z(np)
  ! ... Local parameters
  integer:: ie , ip , lmax , ilm , l , m , lmxx , lmxy
  real(8) ,allocatable :: ylm_rv(:)
  real(8) ,allocatable :: r2_rv(:)

  parameter (lmxx=20)
  double precision :: psi(-1:lmxx),phi(-1:lmxx),yl((lmxx+1)**2),r2, &
       rfac,xxx,fac2l(0:lmxx),psid,psil
  ! ... Heap

  ! --- Make ylm for all points, up to lmxy, which find ---
  lmxy = -1
  do   ip = 1, np
     do   ie = 1, nxi
        lmxy = max(lmxy,lxi(ie,ip))
     enddo
  enddo
  lmxy = lmxy+lmxa
  fac2l(0) = 1
  do  2  l = 1, lmxy
     fac2l(l) = fac2l(l-1) * (2*l-1)
2 enddo
  if (lmxy > lmxx) call rxi('rstr0: lmax exceeds lmxx:',lmxy)
  allocate(ylm_rv((lmxy+1)**2*np))
  allocate(r2_rv(np))
  call ropyln ( np , x , y , z , lmxy , np , ylm_rv , r2_rv )
  ! --- For each point, do ---
  do  5  ip = 1, np
     call pvstr0 ( lmxy , ip , np , ylm_rv , r2_rv , yl , r2 )
     if (r2 < 1d-10) then
        call dpzero(hl(1,1,ip),nlmx*nxi)
        if (iop /= 1) call dpzero(hd(1,1,ip),nlmx*nxi)
        goto 5
     endif
     ! --- Reduced strx hl, or hl and hd for all energies, this point ---
     if (iop == 1) then
        do  10  ie = 1, nxi
           lmax = lmxa+lxi(ie,ip)
           !     ... Not worth vectorizing (I think)
           !         call bessl2(exi(ie)*r2,0,lmax,phi(0),psi(0))
           call besslr(exi(ie)*r2,0,-1,lmax,phi,psi)
           ilm = 0
           rfac = dsqrt(r2)
           xxx = 1d0/r2
           do  201  l = 0, lmax
              rfac = rfac*xxx
              !     ... undo fac2l scaling, to recover standard MSM IV-43.
              !         psil = rfac * psi(l) * fac2l(l)
              psil = rfac * psi(l)
              do  20  m = -l, l
                 ilm = ilm+1
                 hl(ilm,ie,ip) = psil*yl(ilm)
20            enddo
201        enddo
10      enddo
     else
        do  110  ie = 1, nxi
           lmax = lmxa+lxi(ie,ip)
           !         call bessl2(exi(ie)*r2,-1,lmax,phi(-1),psi(-1))
           call besslr(exi(ie)*r2,0,-1,lmax,phi(-1),psi(-1))
           ilm = 0
           rfac = dsqrt(r2)
           xxx = 1d0/r2
           do  1201  l = 0, lmax
              rfac = rfac*xxx
              !         psil = rfac * psi(l) * fac2l(l)
              !         psid = rfac * psi(l-1)*r2/(4*l-2) * fac2l(l)
              psil = rfac * psi(l)
              psid = rfac * psi(l-1)*r2/2
              do  120  m = -l, l
                 ilm = ilm+1
                 hl(ilm,ie,ip) = psil*yl(ilm)
                 hd(ilm,ie,ip) = psid*yl(ilm)
120           enddo
1201       enddo
110     enddo
     endif
5 enddo
  if (allocated(r2_rv)) deallocate(r2_rv)
  if (allocated(ylm_rv)) deallocate(ylm_rv)

end subroutine rstr0
subroutine pvstr0(lmax,ip,nd,ylm,r2,yl,rsq)
  !     implicit none
  integer :: lmax,ip,nd,nlm,ilm
  double precision :: ylm(nd,1), yl(1), r2(ip), rsq

  nlm = (lmax+1)**2
  rsq = r2(ip)
  do  10  ilm = 1, nlm
     yl(ilm) = ylm(ip,ilm)
10 enddo
end subroutine pvstr0

