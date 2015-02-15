      program read_hmlwf
      implicit none

      integer(4) :: nspin,nwf,natom,nqbz,nrws,is
      real(8) :: alat,plat(3,3)
      integer(4),allocatable :: irws(:)
      real(8),allocatable :: pos(:,:),qbz(:,:),wbz(:),rws(:,:)
      complex(8),allocatable :: hk(:,:,:,:)

      integer(4) :: ifi,ir,iq
      real(8) :: pi,rk,q(3)
      complex(8) :: ci,ceikr
      complex(8),allocatable :: hr(:,:,:,:),hq(:,:,:)

      ifi = 99
      open(ifi,file='HMLWF')
      read(ifi,*)nspin,natom,nwf,nqbz,nrws
! nspin = number of spins
! natom = number of atoms contained in the unit cell
! nwf   = number of Wannnier functions
! nqbz  = number of k points
! nrws  = number of Wigner-Seitz (WS) supercell (equal or larger than nqbz)
      allocate(pos(3,natom),qbz(3,nqbz),wbz(nqbz),rws(3,nrws),
     &         irws(nrws),hk(nwf,nwf,nqbz,nspin))
      read(ifi,*)alat ! lattice constant
      read(ifi,*)plat ! unit vector in units of alat
      read(ifi,*)pos  ! atomic coordinates in units of alat
      read(ifi,*)qbz  ! k mesh in units of alat/2pi
      read(ifi,*)wbz  ! weight of k points (= 1/nqbz)
      read(ifi,*)rws  ! positon of the WS supercell
      read(ifi,*)irws ! number of equivalent WS supercells
      do is = 1,nspin
         read(ifi,*)hk(:,:,:,is) ! Hamiltonian in the Wannier gauge
      enddo
      close(99)

!      stop 'ok: read data'

! to get the real space representation,
      pi = 4d0*atan(1d0)
      ci = (0d0,1d0)
      allocate(hr(nwf,nwf,nrws,nspin))
      hr = (0d0,0d0)
      ir = 0
      do is = 1,nspin
      do ir = 1,nrws
      do iq = 1,nqbz
         rk = 2d0*pi*sum(rws(:,ir)*qbz(:,iq))
         ceikr = exp(-ci*rk)
         hr(:,:,ir,is) = hr(:,:,ir,is) + ceikr*hk(:,:,iq,is)/dble(nqbz)
      enddo ! iq 
      enddo ! ir 
      enddo ! is

! to get the Hamiltinian at an arbitary k point,
      allocate(hq(nwf,nwf,nspin))
      hq = (0d0,0d0)
      q = 0
!     q = ???  ! input q(1:3)
      do is = 1,nspin
      do ir = 1,nrws
         rk = 2d0*pi*sum(rws(:,ir)*q(:))
         ceikr = exp(ci*rk) / dble(irws(ir))
         hq(:,:,is) = hq(:,:,is) + ceikr * hr(:,:,ir,is)
      enddo ! ir
      enddo ! is

      end
