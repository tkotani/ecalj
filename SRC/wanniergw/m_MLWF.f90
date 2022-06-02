module m_MLWF
  !      use m_DATA4GW
  !      use m_QG
  implicit none
  integer :: nqbz,nwf,iko_ix,iko_fx
  double precision :: a
  double complex,allocatable :: dnk(:,:,:,:)
  real(8),allocatable:: qreg(:,:)
  !      complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:)
contains
  !$$$ccccccccccccccccccccccccccccccccccccccc
  !$$$      subroutine setup_MLWF()
  !$$$      use m_readeigen,only: readcphi,readgeig,init_readeigen,init_readeigen2
  !$$$      use m_readQG
  !$$$      implicit none
  !$$$      integer :: isp,iqbz,ikp,ib,iwf,nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg,ifhbed,iopen
  !$$$      double complex,allocatable :: geig2(:,:),cphi2(:,:)
  !$$$      real(8)::quu(3),ginv(3,3),det,qlat(3,3)
  !$$$      write(6,*) '--- setup_MLWF ---'
  !$$$      call read_MLWF()
  !$$$      call dinv33x (plat,qlat)  !it was dinv33(plat,1,qlat) by Ferdi
  !$$$      call dinv33(qlat,0,ginv,det)
  !$$$      ifhbed     = iopen('hbe.d',1,0,0)
  !$$$      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
  !$$$      call init_readeigen(ginv,nsp,nband,mrece) !initialization of readEigen
  !$$$      write(6,*)'1111111rrrr'
  !$$$      call init_readeigen2(mrecb,ldim2,mrecg) !initialize m_readeigen
  !$$$      write(6,*)'2222222rrrr'
  !$$$c replace geig and cphi with those for WF
  !$$$c geig2 = cphi2 = 0 for ik > nqbz
  !$$$      call readngmx('QGpsi',ngpmx)
  !$$$      write(6,*)'333333333',ngpmx,nband,ldim2
  !$$$      allocate(geig2(ngpmx,nband))
  !$$$      allocate(cphi2(ldim2,nband))
  !$$$      write(6,*)'4444444444333333333'
  !$$$c
  !$$$c      geig2 = geig
  !$$$c      cphi2 = cphi
  !$$$c      deallocate(geig,cphi)
  !$$$      allocate(geig(ngpmx,nwf,nqbz,nsp))
  !$$$      allocate(cphi(ldim2,nwf,nqbz,nsp))
  !$$$      geig = 0d0
  !$$$      cphi = 0d0
  !$$$      do ikp = 1,nqbz
  !$$$      do isp = 1,nsp
  !$$$         write(6,"('setup_MLWF: ikp isp=',2i5,3f12.5)") ikp,isp,qreg(:,ikp)
  !$$$         call readcphi(qreg(:,ikp),ldim2,isp, quu, cphi2)
  !$$$         if(sum(abs(qreg(:,ikp)-quu))>1d-6) stop 'mmlf111eeeee'
  !$$$         call readgeig(qreg(:,ikp),ngpmx,isp, quu, geig2)
  !$$$c         write(6,*)'22222 ikp=',ikp,isp
  !$$$         if(sum(abs(qreg(:,ikp)-quu))>1d-6) stop 'mmlf222eeeee'
  !$$$         do iwf = 1,nwf
  !$$$         do ib = iko_ix,iko_fx
  !$$$            geig(:,iwf,ikp,isp) = geig(:,iwf,ikp,isp) +
  !$$$     &           geig2(:,ib)*dnk(ib,iwf,ikp,isp)
  !$$$            cphi(:,iwf,ikp,isp) = cphi(:,iwf,ikp,isp) +
  !$$$c     &           cphi2(:,ib,ikp,isp)*dnk(ib,iwf,ikp,isp)
  !$$$     &           cphi2(:,ib)*dnk(ib,iwf,ikp,isp)
  !$$$         enddo ! ib
  !$$$         enddo ! iwf
  !$$$      enddo ! isp
  !$$$      enddo ! ikp
  !$$$      deallocate(geig2,cphi2,dnk)
  !$$$      nband = nwf
  !$$$c      stop 'xxxxxxxxxxxxxxxxx end fo read_MLWF xxxxxxxxx'
  !$$$      end subroutine setup_MLWF
  !$$$ccccccccccccccccccccccccccccccccccccccc
  !$$$      subroutine read_MLWF()
  !$$$      implicit none
  !$$$      double precision :: q(3),eps
  !$$$      parameter (eps=1d-4)
  !$$$      integer :: ifi
  !$$$      integer :: isp,iqbz
  !$$$      integer :: iqbz2
  !$$$
  !$$$      write(6,*) '--- read_MLWF ---'
  !$$$
  !$$$      do isp=1,nsp
  !$$$c file open
  !$$$         ifi = 1000
  !$$$         if (isp.eq.1) then
  !$$$           open(ifi,file='MLWU',form='unformatted',status='old',
  !$$$     &          action='read')
  !$$$         else
  !$$$           open(ifi,file='MLWD',form='unformatted',status='old',
  !$$$     &          action='read')
  !$$$         endif
  !$$$c nqbz mesh-points
  !$$$         read(ifi)nqbz,nwf,iko_ix,iko_fx
  !$$$         if (isp.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),qreg(3,nqbz))
  !$$$         do iqbz = 1,nqbz
  !$$$            read(ifi)iqbz2,q(1:3)
  !$$$            write(6,"(i5,3f13.5)") iqbz,q(:)
  !$$$            qreg(:,iqbz)=q
  !$$$c            q(:) = q(:) - qtt(:,iqbz)
  !$$$c            if (sum(q(:)**2).gt.eps) stop 'MLWU/D: qbz error'
  !$$$c            if (iqbz2.ne.iqbz) stop 'MLWU/D: iqbz error'
  !$$$            read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
  !$$$         enddo
  !$$$
  !$$$c fileclose
  !$$$         close(ifi)
  !$$$      enddo ! isp
  !$$$      end subroutine read_MLWF
  !$$$cccccccccccccccccccccccccccccccccccccccc
end module m_MLWF
