      real(8) function nocctotg2(ispin, ef,esmr,qbz,wbz,
     d                         nband,nqbz)
      use m_readeigen, only: readeval
c Count the total number of electrons under Ef.
c use readeval
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level
c nband   = no. states
c nqbz    = no. k-points
      implicit none
      integer(4) ::it,is,k,ispin,nqbz,nband
      real(8):: wbz(nqbz),qbz(3,nqbz),ekt(nband),esmr,ef,wgt,wiocc
      nocctotg2 = 0d0
      wgt       = 0d0
      do is = 1,ispin
        do k   = 1,nqbz
          call readeval(qbz(:,k),is,  ekt)
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
      print *,' Ef=',ef
      print *,' wgt nocc=',wgt,nocctotg2
      end

c------------------------------------------------------------------

      real(8) function nocctotg(ifev,ispin, ef,esmr,wbz,
     d                         nband,nqbz)
c Count the total number of electrons under Ef.
c ifev(2) = direct access unit file for eigenvalues
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level
c nband   = no. states
c nqbz    = no. k-points
      implicit real*8 (a-h,o-z)
      real(8):: wbz(nqbz)
      integer(4)::   ifev(2)
      real(8):: ekt(nband,nqbz)
      nocctotg  = 0d0
      wgt      = 0d0
      do is = 1,ispin
        call rwdd   (ifev(is),   !read eigenvalues
     d              nband,nqbz, 
     o              ekt)
        do k  = 1,nqbz
          wiocc = 0d0
          do it = 1,nband
c
            if(    ekt(it,k)  + 0.5d0*esmr < ef  ) then
              wiocc = wiocc  + 1d0
            elseif(ekt(it,k) - 0.5d0*esmr < ef  ) then
              wiocc  = wiocc + (ef- (ekt(it,k)-0.5d0*esmr))/esmr
            endif
c
          enddo
          nocctotg = nocctotg + wbz(k)* wiocc
          wgt      = wgt      + wbz(k)
        end do
      end do
      if(ispin==1) nocctotg = nocctotg*2
      print *,' Ef=',ef
      print *,' wgt nocc=',wgt,nocctotg
      end
c------------------------------------------------------------------
      real(8) function nocctotf(ifev,ispin, efm,efp,wfacef,wbz,
     d                         nband,nqbz)
c Count the total number of electrons under Ef.
c ifev(2) = direct access unit file for eigenvalues
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level
c nband   = no. states
c nqbz    = no. k-points
      implicit real*8 (a-h,o-z)
      real(8):: wbz(nqbz)
      integer(4)::   ifev(2)
      real(8):: ekt(nband,nqbz)
      nocctotf  = 0d0
      wgt      = 0d0
      do is = 1,ispin
        call rwdd   (ifev(is),   !read eigenvalues
     d              nband,nqbz, 
     o              ekt)
        do k  = 1,nqbz
          wiocc = 0d0
          do it = 1,nband
            if( ekt(it,k) < efm  ) then
              wiocc = wiocc +1d0
            elseif(efm < ekt(it,k) .and. ekt(it,k)< efp) then
              wiocc = wiocc + wfacef
            endif
            if(efm ==ekt(it,k) .or. ekt(it,k)==efp )
Cstop2rx 2013.08.09 kino     &      stop ' nocctotf: efm or efp coincides with ekt'
     &      call rx( ' nocctotf: efm or efp coincides with ekt')
          enddo
          nocctotf = nocctotf + wbz(k)* wiocc
          wgt      = wgt      + wbz(k)
        end do
      end do
      if(ispin==1) nocctotf = nocctotf*2d0
c     print *,' Ef=',ef
c     print *,' wgt nocc=',wgt,nocctotf
      end
c---------------------------------------
c inverse of
      subroutine invkibzx(irk,nqibz,ngrp,nqbz,
     o  ibzx)
c find k in IBZ for given k in FBZ.
      integer(4) :: irk(nqibz,ngrp),ibzx(nqbz)
      do iqx  = 1,nqbz
        do iqi= 1,nqibz
          do ig = 1,ngrp
            if(irk(iqi,ig)==iqx) then
              ibzx(iqx)=iqi
              goto 999
            endif
          enddo
        enddo
Cstop2rx 2013.08.09 kino        stop ' invkibzx: can not find ibzx'
        call rx( ' invkibzx: can not find ibzx')
  999   continue
c        write(6,*)' ibzx ',iqx, ibzx(iqx), irk(ibzx(iqx),ig)
      enddo
c      stop ' invkibzx:'
      end






