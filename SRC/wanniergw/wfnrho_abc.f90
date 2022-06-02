!      module m_wfrho_abc
!      contains
subroutine calc_phiall_abc2(nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
     npw,mesh,nsp,nband,ldim2,ngpmx, &
     geig,cphi,nwf, &
     phipw,phiaug,phitot)
  !      use m_readeigen,only: readcphif,readgeig
  ! cccccccccccccccccccccccccccccccccccc
  use m_QG,only:ngp
  ! cccccccccccccccccccccccccccccccccccc
  !      use m_LMTO
  !      use m_FFT3D
  implicit none
  ! inputs
  integer :: nq_wfn,nband_wfn,bindx_wfn(nband_wfn),nsp
  double precision :: q_wfn(3,nq_wfn)
  integer :: npw(3),mesh(3),ldim2,nband,ngpmx
  ! outputs
  double complex :: &
       phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
       phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
       phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)
  integer :: isp,iq_wfn,ib,i1,i2,i3
  integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
  double complex :: phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
       phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)
  double complex :: eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
       eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
  real(8):: q(3),quu(3)
  logical :: debug=.false.
  integer:: nwf
  complex(8):: geig(ngpmx,nwf,nq_wfn,nsp),cphi(ldim2,nwf,nq_wfn,nsp)

  !------------------------------------------------------------------
  if(debug) write(*,"(a)") '--- calc_phiall_abc2 ---'

  !      print *,'eeee dim =', mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp

  !      call fft_init(npw,'B')
  call calc_augregion_abc(mesh(1),mesh(2),mesh(3),augregion)
  !      allocate(geig2(ngpmx,nband))
  !      allocate(cphi2(ldim2,nband))
  ! omp parallel do private(iq, eikr,eikT, phipwtmp,phiaugtmp )
  do iq_wfn=1,nq_wfn
     q = q_wfn(1:3,iq_wfn)
     !         write(6,"('ffff iq q=',i5,3f10.4)") iq_wfn,q
     call calc_eikreikT_abc(q,mesh,augregion,eikr,eikT)
     do isp=1,nsp
        !            call readgeig(q,ngpmx,isp, quu, geig2)
        !            if(sum(abs(q-quu))>1d-6) stop 'mmlf111eeeee'
        !            call readcphi(q,ldim2,isp, quu, cphi2)
        !            if(sum(abs(q-quu))>1d-6) stop 'mmlf222eeeee'
        do ib=1,nband_wfn
           ! omp critical
           !     write(*,"(a,i2,2i5,3f10.4,i5)")
           !     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq,qtt(1:3,iq),ib
           !               write(*,"(a,i2,2i5,3f10.4,i5)")
           !     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq_wfn,q_wfn(1:3,iq_wfn),ib
           ! ccccccccccccccccccccccccc
           !      if(iq_wfn>=63) write(6,*)'bbbbb222     ',ib,iq_wfn,isp,ngp(iq_wfn)
           !      if(iq_wfn>=63) write(6,*)'bbbbb222 sum ',sum(geig(1:ngp(iq_wfn),ib,iq_wfn,isp))
           ! ccccccccccccccccccccccccc

           ! omp end critical
           call calc_phi_abc2(geig(1,ib,iq_wfn,isp),ngpmx,cphi(1,ib,iq_wfn,isp),ldim2, &
                iq_wfn,npw,mesh,augregion, &
                phipwtmp,phiaugtmp)
           !               write(6,*)'sumccc=',sum(abs(phipwtmp)),sum(abs(phiaugtmp))
           do i3=1,mesh(3)+1
              do i2=1,mesh(2)+1
                 do i1=1,mesh(1)+1 !   bloch function
                    ! cccccccccccccccccccccccc
                    !         write(6,"('eeeee iq =',100i5)") iq_wfn,nq_wfn,i1,i2,i3,ib,isp
                    ! cccccccccccccccccccccccc
                    phipw(i1,i2,i3,ib,iq_wfn,isp) =  eikr(i1,i2,i3) *phipwtmp(i1,i2,i3)
                    phiaug(i1,i2,i3,ib,iq_wfn,isp)=  eikT(i1,i2,i3) *phiaugtmp(i1,i2,i3)
                    phitot(i1,i2,i3,ib,iq_wfn,isp)= &
                         phipw(i1,i2,i3,ib,iq_wfn,isp)+ &
                         phiaug(i1,i2,i3,ib,iq_wfn,isp)
                    ! cccccccccccccccccccccccc
                    !         write(6,"('ggggg iq =',100i5)") iq_wfn,nq_wfn,i1,i2,i3,ib,isp
                    ! cccccccccccccccccccccccc
                 enddo
              enddo
           enddo
           !       write(*,'(6f10.5)') phitot(:,:,:,ib,iq_wfn,isp)
        enddo               !ib
     enddo                  !isp

  enddo                     !iq
  !      deallocate(geig2,cphi2)
  !      print *,'eeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
  ! cccccccccccc
end subroutine calc_phiall_abc2
! ccccccccccccccccccccccccccccccccccccccccccccc
subroutine calc_augregion_abc(n1,n2,n3,augregion)
  use m_lmf2gw,only: bb,nr,aa,alat,iclass,nclass,bas,nbas,plat
  !      use m_genallcf_v3,only: nbas=>natom, bas=>pos,plat
  implicit none
  ! input
  integer :: n1,n2,n3
  ! output
  integer :: augregion(4,n1+1,n2+1,n3+1)
  ! local
  integer :: nshell
  parameter (nshell=4)
  integer :: i1,i2,i3,j1,j2,j3,ibas,ic
  double precision :: rmax,ratom(3),r(3),rtmp(3),dr
  logical:: debug=.false.

  write(*,*) '--- calc_augregion ---',nclass
  augregion(:,:,:,:)=0

  do ibas=1,nbas
     ic=iclass(ibas)
     rmax = bb(ic)*(exp((nr(ic)-1)*aa(ic))-1d0)
     write(6,*)'ibas, rmax=',ibas,ic,rmax
     do j1=-nshell,nshell
        do j2=-nshell,nshell
           do j3=-nshell,nshell
              rtmp(1)=j1
              rtmp(2)=j2
              rtmp(3)=j3
              call mymatvec(plat,rtmp,ratom,3,3)
              ratom(1:3)=alat*(ratom(1:3)+bas(1:3,ibas))

              do i3=1,n3+1
                 do i2=1,n2+1
                    do i1=1,n1+1

                       rtmp(1)=(i1-1)/dble(n1)
                       rtmp(2)=(i2-1)/dble(n2)
                       rtmp(3)=(i3-1)/dble(n3)
                       !            call mymatvec(plat,rtmp,r,3,3)
                       !            r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
                       r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
                       r(1:3)=alat*r(1:3)
                       dr=sqrt(sum((r(1:3)-ratom(1:3))**2))
                       if (dr < rmax) then
                          if (augregion(4,i1,i2,i3) /= 0) then
                             stop 'calc_augregion_abc: Overlap in augmented region!'
                          endif
                          augregion(1,i1,i2,i3)=j1
                          augregion(2,i1,i2,i3)=j2
                          augregion(3,i1,i2,i3)=j3
                          augregion(4,i1,i2,i3)=ibas
                       endif
                    enddo !i1
                 enddo !i2
              enddo !i3
           enddo !j3
        enddo !j2
     enddo !j1
  enddo !ibas
end subroutine calc_augregion_abc
! ccccccccccccccccccccccccccccccccccc
subroutine calc_eikreikT_abc(kvec,mesh, &
     augregion,eikr,eikT)
  !      use m_LMTO
  use m_lmf2gw,only: bb,nr,aa,alat,nsp,plat
  implicit none
  ! input
  double precision :: kvec(3)
  integer :: mesh(3),augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
  ! output
  double complex ::  eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
       eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
  ! local
  integer :: i1,i2,i3
  double precision :: rtmp(3),r(3),tvec(3)
  double precision :: phase,pi
  pi=4.0d0*atan(1.0d0)
  if( .FALSE. ) write(*,*) 'kvec=',kvec (1:3)
  ! Calculate e^{ikr}
  do i3=1,mesh(3)+1
     do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
           rtmp(1)=(i1-1)/dble(mesh(1))
           rtmp(2)=(i2-1)/dble(mesh(2))
           rtmp(3)=(i3-1)/dble(mesh(3))
           !        call mymatvec(plat,rtmp,r,3,3)
           !        r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
           r(:)= plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
           r(1:3)=alat*r(1:3)
           phase=2.0d0*pi/alat*sum(kvec(1:3)*r(1:3))
           eikr(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
        enddo
     enddo
  enddo

  ! Calculate e^{ikT}
  do i3=1,mesh(3)+1
     do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1

           if (augregion(4,i1,i2,i3) /= 0) then
              rtmp(1:3)=augregion(1:3,i1,i2,i3)
              !             tvec(i) =plat(i,j)*rtmp(j)
              call mymatvec(plat,rtmp,tvec,3,3)
              tvec(1:3)=alat*tvec(1:3)
              !  2 pi  k(i)*tvec(i)
              phase=2.0d0*pi/alat*sum(kvec(1:3)*tvec(1:3))
              eikT(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
           else
              eikT(i1,i2,i3)=0.0d0
           endif
        enddo
     enddo
  enddo

end subroutine calc_eikreikT_abc
! ccccccccccccccccccccccccccccccccccc

subroutine calc_phi_abc2(geig,ngpmx,cphi,ldim2,iq, npw,mesh, &
     !!-- Plane wave expansion of an eigenfunciton (geig,cphi).
     augregion, &
     phipwtmp,phiaugtmp)
  use m_QG,only:ngvecp,ngp
  use m_lmf2gw,only:mnla,iclass,nbas,bas,plat,alat
  !      use m_genallcf_v3,only: nbas=>natom, bas=>pos,plat,alat
  implicit none
  integer :: isp,iq,npw(3),mesh(3),ngpmx,ldim2 !,iband
  integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
  double precision :: qlat(3,3)
  double complex :: eigr,ci = (0.0d0,1.0d0)
  double complex :: &
       phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1), &
       phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)
  integer :: itmp(3),ig,id,i1,i2,i3,j1,j2,j3,ii
  double precision :: rtmp(3),r(3),r0(3) !points to plot
  double precision :: ratom(3) ! atomic points
  double precision :: dr(3)
  complex(8):: geig(ngpmx),cphi(ldim2)
  integer,parameter :: lmax=6
  double complex :: Y(2*lmax+1,lmax+1)
  double precision :: Yreal(2*lmax+1,lmax+1)
  double precision :: calc_gxr
  double precision :: drlength,theta,pphi,sintheta
  integer :: idim,il,mtmp,ntmp,ltmp, ibas
  double complex :: phia
  double precision :: pi=4.0d0*atan(1.0d0),tpi = 8d0*atan(1.0d0)
  logical :: debug=.false.
  if(debug) write(6,*)' --- calc_phi_abc2 ---'
  ! ccccccccccccccccccccccccc
  !      if(iq>=63) write(6,*)'aaaaaa222     ',ngp(iq)
  !      if(iq>=63) write(6,*)'aaaaaa222 sum ',sum(geig(1:ngp(iq)))
  ! ccccccccccccccccccccccccc

  call minv33tp(plat,qlat)
  !      call chkinv33(plat,qlat)
  phipwtmp = 0d0
  !      write(6,*)'aaaaaa1',mesh
  do i3=1,mesh(3)+1
     do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
           !         write(6,*)'aaaaaa',i1,i2,i3
           rtmp(1)=(i1-1)/dble(mesh(1))
           rtmp(2)=(i2-1)/dble(mesh(2))
           rtmp(3)=(i3-1)/dble(mesh(3))
           !         r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
           r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
           !         r0(:) = matmul(qlat,r)
           do ii=1,3
              r0(ii) = sum(qlat(:,ii)*r(:))
           enddo ! ii
           !   r0(i)=G0(j,i)*r(j)*
           !   G(i)= G0(j,i)*nG(i)
           !   exp (i 2 pi  G(i)*r(i) )
           !         if(iq>=63) write(6,*)'aaaaaa222     ',i1,i2,i3,iq,ngp(iq)
           !         if(iq>=63) write(6,*)'aaaaaa222 sum ',sum(geig(1:ngp(iq)))

           !, size(ngvecp(:,:,iq)), size(ngvecp(1,:,iq))
           do ig=1,ngp(iq)
              !           if(iq>=63) write(6,*)'ig ngpmx',ig
              eigr=exp(ci*tpi*sum(r0(:)*dble(ngvecp(:,ig,iq))))
              !           if(iq>=63) write(6,*)'eigr',eigr
              !           if(iq>=63) write(6,*)'geig',geig(ig)
              phipwtmp(i1,i2,i3) = phipwtmp(i1,i2,i3) &
                   + eigr*geig(ig) !,iband) !,iq,isp)
              !           if(iq>=63) write(6,*)'end phip'
              !        phipwtmp(i1,i2,i3)=out_fft(mod(i1-1,npw(1))+1,
              !     &       mod(i2-1,npw(2))+1,mod(i3-1,npw(3))+1)
           enddo ! ig
           !         if(iq>=63) write(6,*)'aaaaaa3333',i1,i2,i3,ngp(iq)
        enddo ! i1
     enddo ! i2
  enddo ! i3
  ! Augmented part
  if(debug) write(6,*)' ----- goto augmented part ------------------'
  phiaugtmp(:,:,:)=0.0d0
  do i3=1,mesh(3)+1
     do i2=1,mesh(2)+1
        do i1=1,mesh(1)+1
           if (augregion(4,i1,i2,i3) /= 0) then
              !          write(6,*)i1,i2,i3,mesh(1)+1,mesh(2)+1,mesh(3)+1
              ! set plane-wave part to zero
              !          phiaugtmp(i1,i2,i3)=0.0d0
              rtmp(1)=(i1-1)/dble(mesh(1))
              rtmp(2)=(i2-1)/dble(mesh(2))
              rtmp(3)=(i3-1)/dble(mesh(3))
              !          call mymatvec(plat,rtmp,r,3,3)
              !          r(1:3)=alat*r(1:3)
              !          r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
              r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
              r(1:3)=alat*r(1:3)
              rtmp(1:3)=augregion(1:3,i1,i2,i3)
              call mymatvec(plat,rtmp,ratom,3,3)
              ratom(1:3)=alat*(ratom(1:3)+bas(1:3,augregion(4,i1,i2,i3)))
              dr(1:3)=r(1:3)-ratom(1:3)
              drlength=sqrt(sum(dr(1:3)**2))
              !---
              !c          call calc_phiaug(dr,augregion(4,i1,i2,i3),
              !c     &         phiaugtmp(i1,i2,i3),isp,iq,iband)
              ! x=r*sin(theta)*cos(pphi)
              ! y=r*sin(theta)*sin(pphi)
              ! z=r*cos(theta)
              theta    = acos(dr(3)/(drlength+1.0d-15))
              sintheta = sqrt(1.0d0-cos(theta)**2)
              pphi     = acos(dr(1)/(drlength*sintheta+1.0d-15))
              if (dr(2) < 0.0d0) pphi=2*pi-pphi
              do il=0,lmax
                 call calc_Ylm(il,theta,pphi, &
                      Y(1:2*il+1,il+1), &
                      Yreal(1:2*il+1,il+1))
              enddo
              !          phia=0.0d0
              do idim=1,ldim2
                 if (mnla(4,idim) == ibas) then
                    mtmp=mnla(1,idim)
                    ntmp=mnla(2,idim)
                    ltmp=mnla(3,idim)
                    if (ltmp > lmax) then
                       stop 'ltmp.gt.lmax!'
                    endif
                    phiaugtmp(i1,i2,i3)=phiaugtmp(i1,i2,i3) +cphi(idim) &
                    *calc_gxr(drlength,ltmp,ntmp,iclass(ibas),isp) &
                         *Yreal(mtmp+ltmp+1,ltmp+1)
                 endif
              enddo
           endif
        enddo
     enddo
  enddo
  if(debug) write(6,*)'---- end of calc_phi_abc2 ----------'
end subroutine calc_phi_abc2


! cccccccccccccccccccccccccccccccccc
subroutine calc_rho_abc(alat_ang,nq_wfn,nband_wfn,mesh, &
     phipw,phiaug,phitot)
  !      use m_LMTO
  use m_lmf2gw,only: bb,nr,aa,alat,nsp,plat
  implicit none
  ! input
  integer :: nq_wfn,nband_wfn,mesh(3)
  double precision :: alat_ang
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
  call wfn2dx_abc(alat_ang,plat,1,1,1,qdum,bindxdum, &
       mesh,rhopw,rhoaug,rhotot)

  deallocate(rhopw)
  deallocate(rhoaug)
  deallocate(rhotot)
end subroutine calc_rho_abc
! cccccccccccccccccccccccccccc
!      end module m_wfrho_abc
