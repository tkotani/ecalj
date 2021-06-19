c      module m_wfrho_abc 
c      contains
      subroutine calc_phiall_abc2(nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     i    npw,mesh,nsp,nband,ldim2,ngpmx,  
     i    geig,cphi,nwf,  
     o    phipw,phiaug,phitot)
c      use m_readeigen,only: readcphif,readgeig
cccccccccccccccccccccccccccccccccccccc
      use m_QG,only:ngp
cccccccccccccccccccccccccccccccccccccc
c      use m_LMTO
c      use m_FFT3D
      implicit none
c inputs
      integer :: nq_wfn,nband_wfn,bindx_wfn(nband_wfn),nsp
      double precision :: q_wfn(3,nq_wfn)
      integer :: npw(3),mesh(3),ldim2,nband,ngpmx
c outputs
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)
c local
      integer :: isp,iq_wfn,ib,i1,i2,i3
      integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
      double complex :: ! temp. array
     &     phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)
      double complex :: ! e^{ikr} and e^{ikT} (T is lattice vector)
     &     eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
      real(8):: q(3),quu(3)
c      double complex,allocatable :: geig2(:,:),cphi2(:,:)
      logical :: debug=.false.
      integer:: nwf
      complex(8):: geig(ngpmx,nwf,nq_wfn,nsp),cphi(ldim2,nwf,nq_wfn,nsp)

!------------------------------------------------------------------
      if(debug) write(*,"(a)") '--- calc_phiall_abc2 ---'

c      print *,'eeee dim =', mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp

!      call fft_init(npw,'B')
      call calc_augregion_abc(mesh(1),mesh(2),mesh(3),augregion)
c      allocate(geig2(ngpmx,nband))
c      allocate(cphi2(ldim2,nband))
!$omp parallel do private(iq, eikr,eikT, phipwtmp,phiaugtmp )
      do iq_wfn=1,nq_wfn
         q = q_wfn(1:3,iq_wfn) 
c         write(6,"('ffff iq q=',i5,3f10.4)") iq_wfn,q
         call calc_eikreikT_abc(q,mesh,augregion,eikr,eikT)
         do isp=1,nsp
c            call readgeig(q,ngpmx,isp, quu, geig2)
c            if(sum(abs(q-quu))>1d-6) stop 'mmlf111eeeee'
c            call readcphi(q,ldim2,isp, quu, cphi2)
c            if(sum(abs(q-quu))>1d-6) stop 'mmlf222eeeee'
            do ib=1,nband_wfn
!$omp critical
c     write(*,"(a,i2,2i5,3f10.4,i5)")
c     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq,qtt(1:3,iq),ib
c               write(*,"(a,i2,2i5,3f10.4,i5)")
c     &         '# isp,iq_wfn,iq,q,ib=',isp,iq_wfn,iq_wfn,q_wfn(1:3,iq_wfn),ib
ccccccccccccccccccccccccccc
c      if(iq_wfn>=63) write(6,*)'bbbbb222     ',ib,iq_wfn,isp,ngp(iq_wfn)
c      if(iq_wfn>=63) write(6,*)'bbbbb222 sum ',sum(geig(1:ngp(iq_wfn),ib,iq_wfn,isp))
ccccccccccccccccccccccccccc

!$omp end critical
               call calc_phi_abc2(geig(1,ib,iq_wfn,isp),ngpmx,cphi(1,ib,iq_wfn,isp),ldim2,
     i              iq_wfn,npw,mesh,augregion,
     o              phipwtmp,phiaugtmp)
c               write(6,*)'sumccc=',sum(abs(phipwtmp)),sum(abs(phiaugtmp))
               do i3=1,mesh(3)+1
               do i2=1,mesh(2)+1
               do i1=1,mesh(1)+1 !   bloch function
cccccccccccccccccccccccccc
c         write(6,"('eeeee iq =',100i5)") iq_wfn,nq_wfn,i1,i2,i3,ib,isp
cccccccccccccccccccccccccc
                  phipw(i1,i2,i3,ib,iq_wfn,isp) =  eikr(i1,i2,i3) *phipwtmp(i1,i2,i3)
                  phiaug(i1,i2,i3,ib,iq_wfn,isp)=  eikT(i1,i2,i3) *phiaugtmp(i1,i2,i3)
                  phitot(i1,i2,i3,ib,iq_wfn,isp)=
     &                 phipw(i1,i2,i3,ib,iq_wfn,isp)+
     &                 phiaug(i1,i2,i3,ib,iq_wfn,isp)
cccccccccccccccccccccccccc
c         write(6,"('ggggg iq =',100i5)") iq_wfn,nq_wfn,i1,i2,i3,ib,isp
cccccccccccccccccccccccccc
               enddo
               enddo
               enddo
c       write(*,'(6f10.5)') phitot(:,:,:,ib,iq_wfn,isp)
            enddo               !ib
         enddo                  !isp

      enddo                     !iq
c      deallocate(geig2,cphi2)
c      print *,'eeeeeeeeeeeeeeeeeeeeeeeeeeeeee'
cccccccccccccc
      end subroutine calc_phiall_abc2
ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_augregion_abc(n1,n2,n3,augregion)
      use m_DATA4GW,only: bb,nr,aa,alat,iclass,nclass,bas,nbas,plat
c      use m_genallcf_v3,only: nbas=>natom, bas=>pos,plat
      implicit none
c input
      integer :: n1,n2,n3
c output
      integer :: augregion(4,n1+1,n2+1,n3+1)
c local
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
c            r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
            r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
            r(1:3)=alat*r(1:3)
            dr=sqrt(sum((r(1:3)-ratom(1:3))**2))
            if (dr.lt.rmax) then
               if (augregion(4,i1,i2,i3).ne.0) then
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
ccccccccccccccccccccccccccccccccccccc
      subroutine calc_eikreikT_abc(kvec,mesh,
     &                           augregion,eikr,eikT)
c      use m_LMTO
      use m_DATA4GW,only: bb,nr,aa,alat,nsp,plat
      implicit none
c input
      double precision :: kvec(3)
      integer :: mesh(3),augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
c output
      double complex :: ! e^{ikr} and e^{ikT} (T is lattice vector)
     &     eikr(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     eikT(mesh(1)+1,mesh(2)+1,mesh(3)+1)
c local
      integer :: i1,i2,i3
      double precision :: rtmp(3),r(3),tvec(3)
      double precision :: phase,pi
      pi=4.0d0*atan(1.0d0)
      if(.false.) write(*,*) 'kvec=',kvec (1:3) 
c Calculate e^{ikr}
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1
        rtmp(1)=(i1-1)/dble(mesh(1))
        rtmp(2)=(i2-1)/dble(mesh(2))
        rtmp(3)=(i3-1)/dble(mesh(3))
!        call mymatvec(plat,rtmp,r,3,3)
c        r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
        r(:)= plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
        r(1:3)=alat*r(1:3)
        phase=2.0d0*pi/alat*sum(kvec(1:3)*r(1:3))
        eikr(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
      enddo
      enddo      
      enddo

c Calculate e^{ikT}
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1

        if (augregion(4,i1,i2,i3).ne.0) then
          rtmp(1:3)=augregion(1:3,i1,i2,i3)
c             tvec(i) =plat(i,j)*rtmp(j)
          call mymatvec(plat,rtmp,tvec,3,3)          
          tvec(1:3)=alat*tvec(1:3)
c  2 pi  k(i)*tvec(i) 
          phase=2.0d0*pi/alat*sum(kvec(1:3)*tvec(1:3))
          eikT(i1,i2,i3)=dcmplx(cos(phase),sin(phase))
        else
          eikT(i1,i2,i3)=0.0d0
        endif
      enddo
      enddo      
      enddo

      end subroutine calc_eikreikT_abc
ccccccccccccccccccccccccccccccccccccc

      subroutine calc_phi_abc2(geig,ngpmx,cphi,ldim2,iq, npw,mesh,
!!-- Plane wave expansion of an eigenfunciton (geig,cphi).
     &     augregion,
     o     phipwtmp,phiaugtmp)
      use m_QG,only:ngvecp,ngp
      use m_DATA4GW,only:mnla,iclass,nbas,bas,plat,alat
c      use m_genallcf_v3,only: nbas=>natom, bas=>pos,plat,alat
      implicit none
      integer :: isp,iq,npw(3),mesh(3),ngpmx,ldim2 !,iband
      integer :: augregion(4,mesh(1)+1,mesh(2)+1,mesh(3)+1)
      double precision :: qlat(3,3)
      double complex :: eigr,ci = (0.0d0,1.0d0)
      double complex :: 
     &     phipwtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1),
     &     phiaugtmp(mesh(1)+1,mesh(2)+1,mesh(3)+1)
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
ccccccccccccccccccccccccccc
c      if(iq>=63) write(6,*)'aaaaaa222     ',ngp(iq)
c      if(iq>=63) write(6,*)'aaaaaa222 sum ',sum(geig(1:ngp(iq)))
ccccccccccccccccccccccccccc

      call minv33tp(plat,qlat)
c      call chkinv33(plat,qlat)
      phipwtmp = 0d0
c      write(6,*)'aaaaaa1',mesh
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1
c         write(6,*)'aaaaaa',i1,i2,i3
         rtmp(1)=(i1-1)/dble(mesh(1))
         rtmp(2)=(i2-1)/dble(mesh(2))
         rtmp(3)=(i3-1)/dble(mesh(3))
c         r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
          r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
!         r0(:) = matmul(qlat,r)
         do ii=1,3
           r0(ii) = sum(qlat(:,ii)*r(:))
         enddo ! ii
c   r0(i)=G0(j,i)*r(j)*
c   G(i)= G0(j,i)*nG(i) 
c   exp (i 2 pi  G(i)*r(i) ) 
c         if(iq>=63) write(6,*)'aaaaaa222     ',i1,i2,i3,iq,ngp(iq)
c         if(iq>=63) write(6,*)'aaaaaa222 sum ',sum(geig(1:ngp(iq)))

!, size(ngvecp(:,:,iq)), size(ngvecp(1,:,iq))
         do ig=1,ngp(iq)
c           if(iq>=63) write(6,*)'ig ngpmx',ig
           eigr=exp(ci*tpi*sum(r0(:)*dble(ngvecp(:,ig,iq))))
c           if(iq>=63) write(6,*)'eigr',eigr
c           if(iq>=63) write(6,*)'geig',geig(ig)
           phipwtmp(i1,i2,i3) = phipwtmp(i1,i2,i3)
     &        + eigr*geig(ig) !,iband) !,iq,isp)
c           if(iq>=63) write(6,*)'end phip'
!        phipwtmp(i1,i2,i3)=out_fft(mod(i1-1,npw(1))+1,
!     &       mod(i2-1,npw(2))+1,mod(i3-1,npw(3))+1)
         enddo ! ig
c         if(iq>=63) write(6,*)'aaaaaa3333',i1,i2,i3,ngp(iq)
      enddo ! i1
      enddo ! i2
      enddo ! i3
c Augmented part
      if(debug) write(6,*)' ----- goto augmented part ------------------'
      phiaugtmp(:,:,:)=0.0d0
      do i3=1,mesh(3)+1
      do i2=1,mesh(2)+1
      do i1=1,mesh(1)+1
        if (augregion(4,i1,i2,i3).ne.0) then
c          write(6,*)i1,i2,i3,mesh(1)+1,mesh(2)+1,mesh(3)+1
c set plane-wave part to zero
c          phiaugtmp(i1,i2,i3)=0.0d0
          rtmp(1)=(i1-1)/dble(mesh(1))
          rtmp(2)=(i2-1)/dble(mesh(2))
          rtmp(3)=(i3-1)/dble(mesh(3))
!          call mymatvec(plat,rtmp,r,3,3)
!          r(1:3)=alat*r(1:3)
c          r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
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
c x=r*sin(theta)*cos(pphi)
c y=r*sin(theta)*sin(pphi)
c z=r*cos(theta)
          theta    = acos(dr(3)/(drlength+1.0d-15))
          sintheta = sqrt(1.0d0-cos(theta)**2)
          pphi     = acos(dr(1)/(drlength*sintheta+1.0d-15))
          if (dr(2).lt.0.0d0) pphi=2*pi-pphi
          do il=0,lmax
             call calc_Ylm(il,theta,pphi,
     &            Y(1:2*il+1,il+1),
     &            Yreal(1:2*il+1,il+1))
          enddo
c          phia=0.0d0
          do idim=1,ldim2
             if (mnla(4,idim).eq.ibas) then
                mtmp=mnla(1,idim)
                ntmp=mnla(2,idim)
                ltmp=mnla(3,idim)
                if (ltmp.gt.lmax) then
                   stop 'ltmp.gt.lmax!'
                endif
                phiaugtmp(i1,i2,i3)=phiaugtmp(i1,i2,i3) +cphi(idim) !,iband) !,iq,isp)
     &               *calc_gxr(drlength,ltmp,ntmp,iclass(ibas),isp)
     &               *Yreal(mtmp+ltmp+1,ltmp+1)
             endif
          enddo
        endif
      enddo
      enddo
      enddo
      if(debug) write(6,*)'---- end of calc_phi_abc2 ----------'
      end subroutine calc_phi_abc2


cccccccccccccccccccccccccccccccccccc
      subroutine calc_rho_abc(alat_ang,nq_wfn,nband_wfn,mesh,
     &   phipw,phiaug,phitot)
c      use m_LMTO
      use m_DATA4GW,only: bb,nr,aa,alat,nsp,plat
      implicit none
c input
      integer :: nq_wfn,nband_wfn,mesh(3)
      double precision :: alat_ang
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

      double complex,allocatable :: rhopw(:,:,:),
     &     rhoaug(:,:,:),rhotot(:,:,:)
c
      integer :: isp,iq,ib,i1,i2,i3
      double precision :: nel
      double precision :: qdum(3),vol
      integer :: bindxdum

      write(*,*) '--- calc_rho ---'
      call mydet3(plat,VOL)
      VOl=abs(VOL)*alat**3
c Allocate rho
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
          rhopw(i1,i2,i3)=rhopw(i1,i2,i3)+
     &         abs(phipw(i1,i2,i3,ib,iq,isp))**2
          rhoaug(i1,i2,i3)=rhoaug(i1,i2,i3)+
     &         abs(phiaug(i1,i2,i3,ib,iq,isp))**2
          rhotot(i1,i2,i3)=rhotot(i1,i2,i3)+
     &         abs(phitot(i1,i2,i3,ib,iq,isp))**2

        enddo
        enddo
        enddo
      enddo
      enddo
      enddo
      rhopw(:,:,:)=
     &     rhopw(:,:,:)/dble(nq_wfn)
      rhoaug(:,:,:)=
     &     rhoaug(:,:,:)/dble(nq_wfn)
      rhotot(:,:,:)=
     &     rhotot(:,:,:)/dble(nq_wfn)


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
      call wfn2dx_abc(alat_ang,plat,1,1,1,qdum,bindxdum,
     &     mesh,rhopw,rhoaug,rhotot)

      deallocate(rhopw)
      deallocate(rhoaug)
      deallocate(rhotot)
      end subroutine calc_rho_abc
cccccccccccccccccccccccccccccc
c      end module m_wfrho_abc 
