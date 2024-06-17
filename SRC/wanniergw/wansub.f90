!-----------------------------------------------------
subroutine calc_npw(nfac, npw)
  use m_QG,only: ngvecp,qqqa,nqnum,ngp
  use m_genallcf_v3,only: alat,plat
  implicit none
  ! input
  integer :: nfac
  ! output
  integer :: npw(3)
  ! local
  integer :: iq,ig,id,itmp(3),ntmp(3)
  double precision :: pi,gtmp(3),gcutmax,gcuttmp,at(3,3),g(3,3)
  logical :: debug=.false.
  write(*,"(a)") '--- calc_npw ---'
  call mytranspose(plat,At,3,3)
  call myinv3(At,G)
  pi=4.0d0*atan(1.0d0)
  ntmp(1:3)=0
  do iq=1,nqnum
     gcutmax=-1.0d0
     do ig=1,ngp(iq)
        call mymatvec(G,dble(ngvecp(1:3,ig,iq)),gtmp,3,3)
        gtmp(1:3)=gtmp(1:3)+qqqa(1:3,iq)
        gtmp(1:3)=gtmp(1:3)*2.0d0*pi/alat
        gcuttmp=sqrt(sum(gtmp(1:3)**2))
        if (gcutmax < gcuttmp) gcutmax=gcuttmp
        do id=1,3
           itmp(id)=abs(ngvecp(id,ig,iq))
           if (ntmp(id) < itmp(id)) ntmp(id)=itmp(id)
        enddo
     enddo
     if(debug) write(*,"(a,2i5,f10.5)") '# iq ngp gcutmax= ',iq,ngp(iq),gcutmax
  enddo
  !      npw(1:3)=2*ntmp(1:3)+2
  npw(1:3)=nfac*ntmp(1:3)+2
  write(*,"(a,3i6)") '# npw(1:3)=',npw(1:3)
end subroutine calc_npw

!-----------------------------------------------------
! Linear interpolation of gx/r
double precision function calc_gxr(r,l,n,ic,isp)
  !      use m_LMTO
  use m_genallcf_v3,only: alat,plat
  use m_lmf2gw,only: bb,nr,aa,gx=>gx_d
  implicit none
  ! input
  double precision :: r
  integer :: l,n,ic,isp
  ! local
  double precision :: r1,r2
  integer :: ir
  ir=1+int(log(r/bb(ic)+1.0d0)/aa(ic))
  if (ir < 1) stop 'ir < 1'
  if (ir > nr(ic)-1) stop 'ir > nr(ic)-1'
  r1=bb(ic)*(exp((ir-1)*aa(ic))-1d0)
  r2=bb(ic)*(exp((ir  )*aa(ic))-1d0)
  if (r1 > r) stop 'r1 > r'
  if (r2 <= r) stop 'r2 <= r'
  calc_gxr=(r-r2)/(r1-r2)*gx(ir,l,n,ic,isp) &
       + (r-r1)/(r2-r1)*gx(ir+1,l,n,ic,isp)
  calc_gxr=calc_gxr/(r+1.0d-20)
end function calc_gxr

!-----------------------------------------------------
subroutine b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,npw, &
     phi,wan)
  !! Make Wannier functions from Bloch functions in real space representation.
  use m_lmf2gw,only: bb,nr,aa
  use m_genallcf_v3,only: alat,nsp=>nspin,plat
  implicit none
  integer :: nq_wfn,nband_wfn,npw(3),bindx_wfn(nband_wfn),tlat(3)
  double precision :: q_wfn(3,nq_wfn),tvec(3),phase,pi,rtmp(3)
  double complex :: &
       phi(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nq_wfn,nsp), &
       wan(npw(1)+1,npw(2)+1,npw(3)+1,nband_wfn,nsp) &
       ,ephase
  integer :: iq,isp
  ! debug:
  !      wan(:,:,:,:,1) = phi(:,:,:,:,2,1)
  !      return
  pi = 4.0d0*atan(1.d0)
  rtmp(:) = dble(tlat(:))
  call mymatvec(plat,rtmp,tvec,3,3)
  tvec(1:3)=alat*tvec(1:3)
  wan = (0.0d0,0.0d0)
  do isp = 1,nsp
     do iq = 1,nq_wfn
        phase=2.0d0*pi/alat*sum(q_wfn(1:3,iq)*tvec(1:3))
        ephase=dcmplx(cos(phase),-sin(phase))
        wan(:,:,:,:,isp) = wan(:,:,:,:,isp) + phi(:,:,:,:,iq,isp)*ephase
     enddo ! iq
  enddo ! isp
  wan = wan / dble(nq_wfn)
end subroutine b2w

!--------------------------------------------------------
subroutine chkinv33(a,b)
  implicit none
  integer :: i,j
  double precision :: a(3,3),b(3,3),c(3,3),r,eps
  eps = 1.d-6
  !      c = matmul(a,b)
  do i=1,3
     do j=1,3
        c(i,j) = sum(a(:,i)*b(:,j))
     enddo
  enddo
  do i=1,3
     c(i,i) = c(i,i)-1.0d0
  enddo
  do i=1,3
     do j=1,3
        r = abs(c(i,j))
        if (r > eps) stop 'chkinv33 error'
     enddo
  enddo
end subroutine chkinv33
!--------------------------------------------------------
subroutine calc_rho_2(alat_ang,nq_wfn,nband_wfn,mesh,rini,rfin, &
     phipw,phiaug,phitot)
  !      use m_LMTO
  use m_lmf2gw,only: bb,nr,aa
  use m_genallcf_v3,only: alat,nsp=>nspin,plat
  implicit none
  ! input
  integer :: nq_wfn,nband_wfn,mesh(3)
  double precision :: alat_ang,rini(3),rfin(3)
  double complex :: &
       phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
       phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
       phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

  double complex,allocatable :: rhopw(:,:,:), &
       rhoaug(:,:,:),rhotot(:,:,:)


  integer :: isp,iq,ib,i1,i2,i3
  double precision :: nel

  double precision :: qdum(3),vol
  integer :: bindxdum
  write(*,*) '--- calc_rho ---'
  call mydet3(plat,VOL)
  VOl=abs(VOL)*alat**3

  ! Allocate rho
  allocate(rhopw(mesh(1)+1,mesh(2)+1,mesh(3)+1))
  allocate(rhoaug(mesh(1)+1,mesh(2)+1,mesh(3)+1))
  allocate(rhotot(mesh(1)+1,mesh(2)+1,mesh(3)+1))

  rhopw(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0
  rhoaug(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0
  rhotot(1:mesh(1)+1,1:mesh(2)+1,1:mesh(3)+1)=0.0d0

  do isp=1,nsp
     do iq=1,nq_wfn
        do ib=1,nband_wfn
           do i3=1,mesh(3)+1
              do i2=1,mesh(2)+1
                 do i1=1,mesh(1)+1
                    rhopw(i1,i2,i3)=rhopw(i1,i2,i3)+ &
                         abs(phipw(i1,i2,i3,ib,iq,isp))**2
                    rhoaug(i1,i2,i3)=rhoaug(i1,i2,i3)+ &
                         abs(phiaug(i1,i2,i3,ib,iq,isp))**2
                    rhotot(i1,i2,i3)=rhotot(i1,i2,i3)+ &
                         abs(phitot(i1,i2,i3,ib,iq,isp))**2

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  rhopw(:,:,:)= &
       rhopw(:,:,:)/dble(nq_wfn)
  rhoaug(:,:,:)= &
       rhoaug(:,:,:)/dble(nq_wfn)
  rhotot(:,:,:)= &
       rhotot(:,:,:)/dble(nq_wfn)


  nel=0.0d0
  do i3=1,mesh(3)
     do i2=1,mesh(2)
        do i1=1,mesh(1)
           nel=nel+rhotot(i1,i2,i3)
        enddo
     enddo
  enddo
  nel=nel*dble(3-nsp)*VOL/dble(mesh(1)*mesh(2)*mesh(3))
  write(*,*) 'nel = ',nel

  qdum(1:3)=0.0d0
  bindxdum=0
  call wfn2dx_2(alat_ang,plat,1,1,1,qdum,bindxdum, &
       mesh,rini,rfin,rhopw,rhoaug,rhotot)

  deallocate(rhopw)
  deallocate(rhoaug)
  deallocate(rhotot)
end subroutine calc_rho_2

