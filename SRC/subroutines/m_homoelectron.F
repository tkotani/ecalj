      module m_homoelectron
!!      use m_genallcf_v3,only: genallcf_v3,
!     !    & alat, plat
      use m_keyvalue,only: getkeyvalue
      implicit none
      public:: Read_qgband,Efermi_egas

      private
      logical:: init=.true.
      real(8):: xmx2(3),rlatp(3,3)
      contains
      
!!! generate |q+G| bands
      !!! calculate |q+G|^2
      subroutine read_qgband(alat,plat,q,ngv,ngvec,isp,gsq,qshort)
      use m_shortn3,only: shortn3_initialize,shortn3
      intent(in)::           alat,plat,q,ngv,ngvec,isp
      intent(out)::                                    gsq,qshort
      integer :: ngv,isp, ngvec(3,*)
      real(8) :: q(3),plat(3,3),alat
      real(8) :: gsq(:),qshort(3)
      real(8),allocatable :: qgw(:)
      real(8) :: tpioa,pi,spene !spin polarization
      integer:: igv
      logical::gamma
      real(8):: qqin(3),qlat(3,3)
      integer:: nlatout(3,48),noutmx,nout,iout

      pi=4d0*datan(1d0)
      tpioa=2d0*pi/alat

            ! get spin polarization : Delta
      call getkeyvalue("GWinput","ene_sppola",spene, default=0d0 )

      !!!q+G
      allocate(qgw(1:3))
      call minv33tp(plat,qlat)

!     ! Get shortest q
      if(init) then
         call shortn3_initialize(qlat)
         init=.false.
      endif   
      noutmx=48
      qqin = matmul(transpose(plat),q)
      call shortn3(qqin,noutmx, nout,nlatout)
      iout=1
      qshort(1:3)=  matmul(qlat(:,:), qqin+nlatout(:,iout))
      
!!! For check
      gsq=0; qgw=0
      do igv=1,ngv
         qgw=tpioa*(qshort+matmul(qlat,ngvec(1:3,igv)))
         gsq(igv)=sum(qgw**2)

         !!! spin polarization
         if (isp==1) then
            gsq(igv)=gsq(igv)-spene
         else
            gsq(igv)=gsq(igv)+spene
         endif
            
      enddo
      end subroutine

      subroutine efermi_egas(ntot,alat,plat,efz)
      intent(in)::           ntot,alat,plat
      intent(out)::                         efz
!     integer(4),intent(in):: ntot
      real(8):: ntot,alat,plat(3,3)
      real(8):: efz
      real(8):: voltot,pi,alpha,det33 !,tripl
      real(8):: qfermi,spene,rydberg,rs

      pi=4d0*datan(1d0)
      voltot = abs(alat**3*det33(plat))
      alpha=(9*pi/4d0)**(1d0/3d0)
      call getkeyvalue("GWinput","ene_sppola",spene, default=0d0 )
      
      efz=(ntot*3*pi**2/voltot)**(2d0/3d0) - spene**2*3/8d0  !ef is calculated from ntot.   
      qfermi= dsqrt(efz)
      rs    = alpha/qfermi

      write(6,"('efermi_egas:: efermi[ryd],efermi[eV], qfermi, rs',5f9.5)")
     &     efz,efz*rydberg(),qfermi,rs
      end subroutine efermi_egas

c--------------------------------
      real(8) function det33(am)
      implicit none
      real(8),intent(in) :: am(3,3)
      det33= am(1,1)*am(2,2)*am(3,3)
     &       -am(1,1)*am(3,2)*am(2,3)
     &       -am(2,1)*am(1,2)*am(3,3)
     &       +am(2,1)*am(3,2)*am(1,3)
     &       +am(3,1)*am(1,2)*am(2,3)
     &       -am(3,1)*am(2,2)*am(1,3)
      end function det33
      
      end module m_homoelectron
