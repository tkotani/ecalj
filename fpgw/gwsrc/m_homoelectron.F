      module m_homoelectron

!!      use m_genallcf_v3,only: genallcf_v3,
!!    & alat, plat
      implicit none

      logical,private:: init=.true.
      real(8),private:: xmx2(3),rlapt(3,3)


      contains
!!! generate |q+G| bands
      !!! calculate |q+G|^2
      subroutine read_qgband(alat,plat,q,ngv,ngvec,isp,spene,gsq,qshort)
      
      integer(4), intent(in):: ngv,isp
      integer(4),intent(in) :: ngvec(3,*)
      real(8), intent(in) :: q(3),plat(3,3),alat,spene !spin polarization
      real(8), intent(out) :: gsq(:),qshort(3)
      real(8),allocatable :: qgw(:)
      real(8) :: tpioa,qlat(3,3),pi
      integer:: igv
      logical::gamma

      real(8):: qqin(3),rlatp(3,3)
      integer:: nlatout(3,48),noutmx,nout,iout

      pi=4d0*datan(1d0)
      tpioa=2d0*pi/alat
      
      !!!q+G
      allocate(qgw(1:3))
      call minv33tp(plat,qlat)

!     ! Get shortest q
      if(init) then
         call shortn3_initialize(qlat,rlatp,xmx2)
         init=.false.
      endif   
      noutmx=48
      qqin = matmul(transpose(plat),q)
      call shortn3(rlatp,xmx2, qqin,noutmx, nout,nlatout)
      iout=1
      qshort(1:3)=  matmul(qlat(:,:), qqin+nlatout(:,iout))
      
!!! For check
      gsq=0; qgw=0
      do igv=1,ngv
         qgw=tpioa*(qshort+matmul(qlat,ngvec(1:3,igv)))
         print *, "0822, |QGw|",sqrt(sum(qgw**2))

c$$$         if (gamma) then
c$$$            gsq(igv)=sqrt(sum((tpioa*matmul(qlat,ngvec(1:3,igv)))**2))

         gsq(igv)=sum(qgw**2)

         !!! spin polarization
         if (isp==1) then
            gsq(igv)=gsq(igv)-spene
         else
            gsq(igv)=gsq(igv)+spene
         endif
            
         print *, "0822, |G|^2",gsq(igv)
      enddo
      end subroutine

c$$$      subroutine efermi_egas
c$$$
c$$$      real(8)::ntot,bzvol
c$$$      
c$$$c$$$ あと,BZの体積をvoltotとして,
c$$$c$$$ バレンス電子数にntotに対してefermi,qfermiを与えるには,
c$$$c$$$  efz=(ntot*3*pi**2/voltot)**(2d0/3d0)
c$$$!     ef is calculated from ntot.
c$$$c$$$         qfermi= dsqrt(efz)
c$$$c$$$         alpha = (9*pi/4d0)**(1d0/3d0)
c$$$c$$$         rs    = alpha/qfermi
c$$$c$$$ ここでrsは電子ガスのパラメーターで電子間の距離を表します.
      
      end module m_homoelectron
