      module m_MLWF
c      use m_DATA4GW
c      use m_QG
      implicit none
      integer :: nqbz,nwf,iko_ix,iko_fx
      double precision :: a
      double complex,allocatable :: dnk(:,:,:,:)
      real(8),allocatable:: qreg(:,:)
c      complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:)
      contains
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$      subroutine setup_MLWF()
c$$$      use m_readeigen,only: readcphi,readgeig,init_readeigen,init_readeigen2
c$$$      use m_readQG
c$$$      implicit none
c$$$      integer :: isp,iqbz,ikp,ib,iwf,nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg,ifhbed,iopen
c$$$      double complex,allocatable :: geig2(:,:),cphi2(:,:)
c$$$      real(8)::quu(3),ginv(3,3),det,qlat(3,3)
c$$$      write(6,*) '--- setup_MLWF ---'
c$$$      call read_MLWF()
c$$$      call dinv33x (plat,qlat)  !it was dinv33(plat,1,qlat) by Ferdi
c$$$      call dinv33(qlat,0,ginv,det)
c$$$      ifhbed     = iopen('hbe.d',1,0,0)
c$$$      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
c$$$      call init_readeigen(ginv,nsp,nband,mrece) !initialization of readEigen
c$$$      print *,'1111111rrrr'
c$$$      call init_readeigen2(mrecb,ldim2,mrecg) !initialize m_readeigen
c$$$      print *,'2222222rrrr'
c$$$c replace geig and cphi with those for WF
c$$$c geig2 = cphi2 = 0 for ik > nqbz
c$$$      call readngmx('QGpsi',ngpmx)
c$$$      print *,'333333333',ngpmx,nband,ldim2
c$$$      allocate(geig2(ngpmx,nband))
c$$$      allocate(cphi2(ldim2,nband))
c$$$      print *,'4444444444333333333'
c$$$c
c$$$c      geig2 = geig
c$$$c      cphi2 = cphi
c$$$c      deallocate(geig,cphi)
c$$$      allocate(geig(ngpmx,nwf,nqbz,nsp))
c$$$      allocate(cphi(ldim2,nwf,nqbz,nsp))
c$$$      geig = 0d0
c$$$      cphi = 0d0
c$$$      do ikp = 1,nqbz
c$$$      do isp = 1,nsp
c$$$         write(6,"('setup_MLWF: ikp isp=',2i5,3f12.5)") ikp,isp,qreg(:,ikp)
c$$$         call readcphi(qreg(:,ikp),ldim2,isp, quu, cphi2)
c$$$         if(sum(abs(qreg(:,ikp)-quu))>1d-6) stop 'mmlf111eeeee'
c$$$         call readgeig(qreg(:,ikp),ngpmx,isp, quu, geig2)
c$$$c         print *,'22222 ikp=',ikp,isp
c$$$         if(sum(abs(qreg(:,ikp)-quu))>1d-6) stop 'mmlf222eeeee'
c$$$         do iwf = 1,nwf
c$$$         do ib = iko_ix,iko_fx
c$$$            geig(:,iwf,ikp,isp) = geig(:,iwf,ikp,isp) + 
c$$$     &           geig2(:,ib)*dnk(ib,iwf,ikp,isp)
c$$$            cphi(:,iwf,ikp,isp) = cphi(:,iwf,ikp,isp) + 
c$$$c     &           cphi2(:,ib,ikp,isp)*dnk(ib,iwf,ikp,isp)
c$$$     &           cphi2(:,ib)*dnk(ib,iwf,ikp,isp)
c$$$         enddo ! ib
c$$$         enddo ! iwf
c$$$      enddo ! isp
c$$$      enddo ! ikp
c$$$      deallocate(geig2,cphi2,dnk)
c$$$      nband = nwf
c$$$c      stop 'xxxxxxxxxxxxxxxxx end fo read_MLWF xxxxxxxxx'
c$$$      end subroutine setup_MLWF
c$$$ccccccccccccccccccccccccccccccccccccccc
c$$$      subroutine read_MLWF()
c$$$      implicit none
c$$$      double precision :: q(3),eps
c$$$      parameter (eps=1d-4)
c$$$      integer :: ifi
c$$$      integer :: isp,iqbz
c$$$      integer :: iqbz2
c$$$
c$$$      write(6,*) '--- read_MLWF ---'
c$$$
c$$$      do isp=1,nsp
c$$$c file open
c$$$         ifi = 1000
c$$$         if (isp.eq.1) then
c$$$           open(ifi,file='MLWU',form='unformatted',status='old',
c$$$     &          action='read')
c$$$         else
c$$$           open(ifi,file='MLWD',form='unformatted',status='old',
c$$$     &          action='read')
c$$$         endif
c$$$c nqbz mesh-points
c$$$         read(ifi)nqbz,nwf,iko_ix,iko_fx
c$$$         if (isp.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),qreg(3,nqbz))
c$$$         do iqbz = 1,nqbz
c$$$            read(ifi)iqbz2,q(1:3)
c$$$            write(6,"(i5,3f13.5)") iqbz,q(:)
c$$$            qreg(:,iqbz)=q
c$$$c            q(:) = q(:) - qtt(:,iqbz)
c$$$c            if (sum(q(:)**2).gt.eps) stop 'MLWU/D: qbz error'
c$$$c            if (iqbz2.ne.iqbz) stop 'MLWU/D: iqbz error'
c$$$            read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
c$$$         enddo
c$$$
c$$$c fileclose
c$$$         close(ifi)
c$$$      enddo ! isp
c$$$      end subroutine read_MLWF
c$$$cccccccccccccccccccccccccccccccccccccccc
      end module m_MLWF
