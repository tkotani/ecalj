c contains all information of FPLMTO-LDA.
      module m_DATA4GW
      implicit none
c scalars
      integer :: nsp,   nbas,   nclass, nrmx,   ncoremx,
     &     lmxamx,ngpmx,  nband,  ldim2,   nqtt,
     &     nphi, nphimx
      double precision :: alat,efermi
      logical :: nocore

c arrays
      integer,allocatable :: iclass(:),lmxa(:),
     &     nr(:),konf(:,:),ncore(:)

      integer,allocatable :: nindx(:),lindx(:),ibasindx(:)
      double precision,allocatable :: zz(:),aa(:),bb(:),
     &     bas(:,:),plat(:,:),qtt(:,:),ec(:,:,:),evl(:,:,:),
     &     vxclda(:,:,:),gx(:,:,:,:,:),gcore(:,:,:,:)

      complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:)

c 
      integer,allocatable :: mnla(:,:)
      contains
ccccccccccccccccccccccccccccccccccccc
      subroutine read_DATA4GW()
      implicit none
      integer :: ifi
      write(6,*) '# --- read_DATA4GW ---'
      ifi = 1000
      open(ifi,file='DATA4GW_V2',form='unformatted',status='old',
     &     action='read')
c from gwinput0_v2
      read(ifi)
     &     nsp,   nbas,   nclass, nrmx,   ncoremx,
     &     lmxamx,ngpmx,  nband,  ldim2,   nqtt,
     &     nphi, nphimx

      write(6,"(a,5i5)") 'nsp,nbas,nclass,nrmx,ncoremx=',
     &     nsp,   nbas,   nclass, nrmx,   ncoremx

      write(6,"(a,7i5)")
     &     'lmxamx,ngpmx,nband,ldim2,nqtt,nphi,nphimx=',
     &     lmxamx,ngpmx,  nband,  ldim2,   nqtt,
     &     nphi, nphimx

c allocate arrays
      allocate(iclass(nbas),lmxa(nclass),nr(nclass),
     &     konf(0:lmxamx,nclass),ncore(nclass))
      allocate(zz(nclass),aa(nclass),bb(nclass),bas(3,nbas),
     &     plat(3,3),qtt(3,nqtt),ec(ncoremx, nclass, nsp),
     &     evl(nband, nqtt, nsp),vxclda (nband, nqtt, nsp),
     &     gx (nrmx, 0:lmxamx, nphimx, nclass,nsp),
     &     gcore(nrmx, ncoremx, nclass,nsp))
      allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
c from gwinput0_v2x
      read(ifi)
     &              iclass,lmxa,nr,konf,ncore,
     &              zz,aa,bb,bas,alat,plat,
     &              qtt,efermi,ec,evl,vxclda,gx,gcore,nocore !omit cphi and geig
cv2
     &              ,nindx,lindx,ibasindx
      close(ifi)
c      write(*,*) 'ibasindx=',ibasindx
c      write(*,*) 'efermi=',efermi
      end subroutine read_DATA4GW
ccccccccccccccccccccccccccccccccccccc
c$$$      subroutine read_CphiGeig()
c$$$      implicit none
c$$$      integer :: ificg,ikp,isp,ngp
c$$$
c$$$      write(6,*) '#--- read_CphiGeig() ---'
c$$$      ificg=1000
c$$$
c$$$      allocate(
c$$$     &     geig(ngpmx,  nband, nqtt,nsp),
c$$$     &     cphi(ldim2, nband, nqtt,nsp) 
c$$$     &     )
c$$$      open(ificg,file='CphiGeig',form='unformatted',
c$$$     &     status='old',action='read')      
c$$$      do ikp = 1,nqtt
c$$$        do isp =1,nsp
c$$$
c$$$          read(ificg) cphi(1:ldim2,1:nband,ikp,isp) 
c$$$          geig(1:ngpmx,1:nband,ikp,isp)=0d0
c$$$          read(ificg) ngp
c$$$          write(6,"(a,i4,i5,3f18.10)") '  ikp ngp,q=',
c$$$     &         ikp,ngp,qtt(1:3,ikp)
c$$$ccccccccccccccccc
c$$$          read(ificg) geig(1:ngp,1:nband,ikp,isp)
c$$$        enddo
c$$$      enddo
c$$$      close(ificg)
c$$$      end subroutine read_CphiGeig
ccccccccccccccccccccccccccccccccccccc
      subroutine set_mnla()
      implicit none
      integer :: ix
      integer :: ibas,lx,nx,mx,ic
      integer,allocatable :: nvmax(:,:)
c
c
      write(*,*) '--- set_mnla ---'
      allocate(nvmax(0:lmxamx,nclass))
      nvmax(:,:)=0
      do ix =1,ldim2
        nx = nindx(ix)  
        lx   = lindx(ix)
        ibas =ibasindx(ix)
        ic = iclass(ibas)
        if( nx> nvmax(lx,ic) ) nvmax(lx,ic) = nx
      enddo

      write(6,"(5a5)") 'm','n','l','atom','ix'
      allocate(mnla(4,ldim2))
      ix=0

c      do nx = 1,nvmax(lx,ic)
      do nx = 1,3
        do ibas = 1,nbas
          ic = iclass(ibas)
          do lx = 0,lmxa(ic)
            if (nx.gt.nvmax(lx,ic)) then
              cycle
            endif

            do mx = -lx,lx
              ix=ix+1
              mnla(1,ix)=mx
              mnla(2,ix)=nx
              mnla(3,ix)=lx
              mnla(4,ix)=ibas
              write(6,"(5i5)") mnla(1:4,ix),ix
            enddo
          enddo
        enddo
      enddo

      if (ix.ne.ldim2) then
        write(6,*) 'Error in set_mnla: ix!=ldim2'
        write(6,*) 'ix2,ldim2=',ix,ldim2
        stop 'Error in set_mnla: ix!=ldim2'
      endif
      deallocate(nvmax)
      end subroutine set_mnla
ccccccccccccccccccccccccccccccccccccc
      end module m_DATA4GW
