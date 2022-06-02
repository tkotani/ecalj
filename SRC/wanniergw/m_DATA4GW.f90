! contains all information of FPLMTO-LDA.
module m_DATA4GW
  implicit none
  ! scalars
  integer :: nsp,   nbas,   nclass, nrmx,   ncoremx, &
       lmxamx,ngpmx,  nband,  ldim2,   nqtt, &
       nphi, nphimx
  double precision :: alat,efermi
  logical :: nocore
  integer,allocatable :: iclass(:),lmxa(:), &
       nr(:),konf(:,:),ncore(:)
  integer,allocatable :: nindx(:),lindx(:),ibasindx(:)
  double precision,allocatable :: zz(:),aa(:),bb(:), &
       bas(:,:),plat(:,:),qtt(:,:),ec(:,:,:),evl(:,:,:), &
       vxclda(:,:,:),gx(:,:,:,:,:),gcore(:,:,:,:)
  complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:)
  integer,allocatable :: mnla(:,:)
contains
  ! ccccccccccccccccccccccccccccccccccc
  subroutine read_DATA4GW()
    implicit none
    integer :: ifi,ifile_handle
    write(6,*) '# --- read_DATA4GW ---'
    ifi = ifile_handle()
    open(ifi,file='DATA4GW_V2',form='unformatted',status='old', &
         action='read')
    ! from gwinput0_v2
    read(ifi) &
         nsp,   nbas,   nclass, nrmx,   ncoremx, &
         lmxamx,ngpmx,  nband,  ldim2,   nqtt, &
         nphi, nphimx

    write(6,"(a,5i5)") 'nsp,nbas,nclass,nrmx,ncoremx=', &
         nsp,   nbas,   nclass, nrmx,   ncoremx

    write(6,"(a,7i5)") &
         'lmxamx,ngpmx,nband,ldim2,nqtt,nphi,nphimx=', &
         lmxamx,ngpmx,  nband,  ldim2,   nqtt, &
         nphi, nphimx

    ! allocate arrays
    allocate(iclass(nbas),lmxa(nclass),nr(nclass), &
         konf(0:lmxamx,nclass),ncore(nclass))
    allocate(zz(nclass),aa(nclass),bb(nclass),bas(3,nbas), &
         plat(3,3),qtt(3,nqtt),ec(ncoremx, nclass, nsp), &
         evl(nband, nqtt, nsp),vxclda (nband, nqtt, nsp), &
         gx (nrmx, 0:lmxamx, nphimx, nclass,nsp), &
         gcore(nrmx, ncoremx, nclass,nsp))
    allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
    ! from gwinput0_v2x
    read(ifi) iclass,lmxa,nr,konf,ncore, &
         zz,aa,bb,bas,alat,plat, &
         qtt,efermi,ec,evl,vxclda,gx,gcore,nocore,nindx,lindx,ibasindx
    close(ifi)
    !      write(*,*) 'ibasindx=',ibasindx
    !      write(*,*) 'efermi=',efermi
  end subroutine read_DATA4GW
  ! ccccccccccccccccccccccccccccccccccc
  !$$$      subroutine read_CphiGeig()
  !$$$      implicit none
  !$$$      integer :: ificg,ikp,isp,ngp
  !$$$
  !$$$      write(6,*) '#--- read_CphiGeig() ---'
  !$$$      ificg=1000
  !$$$
  !$$$      allocate(
  !$$$     &     geig(ngpmx,  nband, nqtt,nsp),
  !$$$     &     cphi(ldim2, nband, nqtt,nsp)
  !$$$     &     )
  !$$$      open(ificg,file='CphiGeig',form='unformatted',
  !$$$     &     status='old',action='read')
  !$$$      do ikp = 1,nqtt
  !$$$        do isp =1,nsp
  !$$$
  !$$$          read(ificg) cphi(1:ldim2,1:nband,ikp,isp)
  !$$$          geig(1:ngpmx,1:nband,ikp,isp)=0d0
  !$$$          read(ificg) ngp
  !$$$          write(6,"(a,i4,i5,3f18.10)") '  ikp ngp,q=',
  !$$$     &         ikp,ngp,qtt(1:3,ikp)
  !$$$ccccccccccccccccc
  !$$$          read(ificg) geig(1:ngp,1:nband,ikp,isp)
  !$$$        enddo
  !$$$      enddo
  !$$$      close(ificg)
  !$$$      end subroutine read_CphiGeig
  ! ccccccccccccccccccccccccccccccccccc
  subroutine set_mnla()
    implicit none
    integer :: ix
    integer :: ibas,lx,nx,mx,ic
    integer,allocatable :: nvmax(:,:)


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

    !      do nx = 1,nvmax(lx,ic)
    do nx = 1,3
       do ibas = 1,nbas
          ic = iclass(ibas)
          do lx = 0,lmxa(ic)
             if (nx > nvmax(lx,ic)) then
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

    if (ix /= ldim2) then
       write(6,*) 'Error in set_mnla: ix!=ldim2'
       write(6,*) 'ix2,ldim2=',ix,ldim2
       stop 'Error in set_mnla: ix!=ldim2'
    endif
    deallocate(nvmax)
  end subroutine set_mnla
  ! ccccccccccccccccccccccccccccccccccc
end module m_DATA4GW
