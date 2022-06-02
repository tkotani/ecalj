real(8) function nocctotg2(ispin, ef,esmr,qbz,wbz, &
     nband,nqbz)
  use m_readeigen, only: readeval
  ! Count the total number of electrons under Ef.
  ! use readeval
  ! ispin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  implicit none
  integer ::it,is,k,ispin,nqbz,nband
  real(8):: wbz(nqbz),qbz(3,nqbz),ekt(nband),esmr,ef,wgt,wiocc
  nocctotg2 = 0d0
  wgt       = 0d0
  do is = 1,ispin
     do k   = 1,nqbz
        ekt= readeval(qbz(:,k),is)
        wiocc = 0d0
        do it = 1,nband
           if(    ekt(it)  + 0.5d0*esmr < ef  ) then
              wiocc = wiocc  + 1d0
           elseif(ekt(it) - 0.5d0*esmr < ef  ) then
              wiocc  = wiocc + (ef- (ekt(it)-0.5d0*esmr))/esmr
           endif
        enddo
        nocctotg2 = nocctotg2 + wbz(k)* wiocc
        wgt       = wgt       + wbz(k)
     enddo
  enddo
  if(ispin==1) nocctotg2 = nocctotg2*2
  write(6,*)' Ef=',ef
  write(6,*)' wgt nocc=',wgt,nocctotg2
END function nocctotg2

!------------------------------------------------------------------

real(8) function nocctotg(ifev,ispin, ef,esmr,wbz, &
     nband,nqbz)
  ! Count the total number of electrons under Ef.
  ! ifev(2) = direct access unit file for eigenvalues
  ! ispin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  implicit real*8 (a-h,o-z)
  real(8):: wbz(nqbz)
  integer::   ifev(2),ispin,is,k,nqbz,nband,it
  real(8):: ekt(nband,nqbz)
  nocctotg  = 0d0
  wgt      = 0d0
  do is = 1,ispin
     call rwdd   (ifev(is),  nband,nqbz, &  !  read eigenvalues
          ekt)
     do k  = 1,nqbz
        wiocc = 0d0
        do it = 1,nband

           if(    ekt(it,k)  + 0.5d0*esmr < ef  ) then
              wiocc = wiocc  + 1d0
           elseif(ekt(it,k) - 0.5d0*esmr < ef  ) then
              wiocc  = wiocc + (ef- (ekt(it,k)-0.5d0*esmr))/esmr
           endif

        enddo
        nocctotg = nocctotg + wbz(k)* wiocc
        wgt      = wgt      + wbz(k)
     end do
  end do
  if(ispin==1) nocctotg = nocctotg*2
  write(6,*)' Ef=',ef
  write(6,*)' wgt nocc=',wgt,nocctotg
END function nocctotg
!------------------------------------------------------------------
real(8) function nocctotf(ifev,ispin, efm,efp,wfacef,wbz, &
     nband,nqbz)
  ! Count the total number of electrons under Ef.
  ! ifev(2) = direct access unit file for eigenvalues
  ! ispin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  implicit real*8 (a-h,o-z)
  real(8):: wbz(nqbz)
  integer::   ifev(2),ispin,nqbz,nband,is,it,k
  real(8):: ekt(nband,nqbz)
  nocctotf  = 0d0
  wgt      = 0d0
  do is = 1,ispin
     call rwdd   (ifev(is),  & ! & read eigenvalues
     nband,nqbz, &
          ekt)
     do k  = 1,nqbz
        wiocc = 0d0
        do it = 1,nband
           if( ekt(it,k) < efm  ) then
              wiocc = wiocc +1d0
           elseif(efm < ekt(it,k) .AND. ekt(it,k)< efp) then
              wiocc = wiocc + wfacef
           endif
           if(efm ==ekt(it,k) .OR. ekt(it,k)==efp ) &
                                ! top2rx 2013.08.09 kino     &      stop ' nocctotf: efm or efp coincides with ekt'
                call rx( ' nocctotf: efm or efp coincides with ekt')
        enddo
        nocctotf = nocctotf + wbz(k)* wiocc
        wgt      = wgt      + wbz(k)
     end do
  end do
  if(ispin==1) nocctotf = nocctotf*2d0
  !     write(6,*)' Ef=',ef
  !     write(6,*)' wgt nocc=',wgt,nocctotf
END function nocctotf
!---------------------------------------
! inverse of
subroutine invkibzx(irk,nqibz,ngrp,nqbz, &
     ibzx)
  ! find k in IBZ for given k in FBZ.
  integer :: irk(nqibz,ngrp),ibzx(nqbz),nqbz,ngrp,iqx,iqi,ig,nqibz
  do iqx  = 1,nqbz
     do iqi= 1,nqibz
        do ig = 1,ngrp
           if(irk(iqi,ig)==iqx) then
              ibzx(iqx)=iqi
              goto 999
           endif
        enddo
     enddo
     ! top2rx 2013.08.09 kino        stop ' invkibzx: can not find ibzx'
     call rx( ' invkibzx: can not find ibzx')
999  continue
     !        write(6,*)' ibzx ',iqx, ibzx(iqx), irk(ibzx(iqx),ig)
  enddo
  !      stop ' invkibzx:'
end subroutine invkibzx






