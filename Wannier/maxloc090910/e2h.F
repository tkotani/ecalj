c eval(ib,ik,is): GW eigenvalue including the imaginary part,
c   measured from the Fermi level
c   Im part is positive (negative) when the real part is negative (positive)
c evec(ii,ib,ik,is) = < psiw(ii,ik,is) | psih(ib,ik,is) >
c   psiw(ii,ik,is) is the Bloch function in the Wannier gauge
c   psih(ib,ik,is) is the eivenfunction of the (reduced) Hamiltonian
c ib: band index
c ik: k point in the full BZ
c is: spin
      implicit real*8 (a-h,o-z)
      complex(8),allocatable :: eval(:,:,:),evec(:,:,:,:),
     &   zw1(:,:),zw2(:,:),evec2(:,:),zwc(:,:,:,:)
      real(8),allocatable :: qbz(:,:),zfac(:,:,:),eval2(:)

      open(99,file='EGW')
      read(99,*)nspin,nqbz,nwf
      allocate(eval(nwf,nqbz,nspin),evec(nwf,nwf,nqbz,nspin),
     &         qbz(3,nqbz),zfac(nwf,nqbz,nspin))
      read(99,*)qbz
      read(99,*)eval
      read(99,*)evec
      read(99,*)zfac
      close(99)

      allocate(zw1(nwf,nwf),zw2(nwf,nwf),evec2(nwf,nwf),eval2(nwf),
     &         zwc(nwf,nwf,nqbz,nspin))
      zwc = 0.0d0
      do is = 1,nspin
      do iq = 1,nqbz
         zw1 = 0.0d0
         do in = 1,nwf ! Wannier gauge
         do im = 1,nwf ! Wannier gauge
         do ii = 1,nwf ! Ham. gauge
           zw1(im,in) = zw1(im,in)
     &    + evec(im,ii,iq,is)*dreal(eval(ii,iq,is))
     &    *dconjg(evec(in,ii,iq,is))
         enddo ! ii  
         enddo ! im
         enddo ! in
c 
         call hermite(zw1,nwf)
         call diag_hm(zw1,nwf,eval2,evec2)

         do in = 1,nwf ! Wannier gauge
         do im = 1,nwf ! Wannier gauge
         do ii = 1,nwf ! Ham. gauge
            zwc(im,in,iq,is) = zwc(im,in,iq,is) + 
     &    + evec2(im,ii)*eval(ii,iq,is)*dconjg(evec2(in,ii))
         enddo ! ii  
         enddo ! im
         enddo ! in
      enddo ! iq
      enddo ! is

      open(98,file='HGWC_hermite')
      write(98,*)nspin,nqbz,nwf
      write(98,*)qbz
      write(98,*)zwc
      write(98,*)zfac
      close(98)

      stop 'read_egw ok'

      end
c-----------------------------------------------------------------------
      subroutine hermite(zw,n)
      implicit real*8(a-h,o-z)
      complex(8) :: zw(n,n),zw2(n,n)

      zw2 = zw
      do i=1,n
      do j=i,n
         zw(i,j) = zw2(i,j) + dconjg(zw2(j,i))
         zw(j,i) = dconjg(zw(i,j))
      enddo   
      enddo   

      return
      end
c-----------------------------------------------------------------------
      subroutine diag_hm(zmat,ndim,eval,evecc)
      implicit real*8(a-h,o-z)

      complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
      complex(8):: zmat(ndim,ndim),evecc(ndim,ndim)
      real(8):: eval(ndim),wk(ndim,11)
      integer iwk(ndim)

      allocate(zmat2(ndim,ndim),ovlpc(ndim,ndim))

      nev  = ndim
      nmx  = ndim

      zmat2 = zmat

      ovlpc = (0d0,0d0)
      do i=1,ndim
         ovlpc(i,i) = (1d0,0d0)
      enddo

      evecc = (0d0,0d0)
      eval = 0d0

!      call diagno(ndim,zmat2,ovlpc,wk,iwk,evecc,eval)
      call diagcv(ovlpc,zmat2, evecc, ndim, eval, nmx, 1d99, nev)

      deallocate(zmat2,ovlpc)

      return
      end
c-----------------------------------------------------------------------
      subroutine wan2ham(zin,evecc,ndim,zout)
      implicit real*8(a-h,o-z)

      integer ii,jj
      complex(8) :: zin(ndim,ndim),evecc(ndim,ndim),zout(ndim,ndim),
     &              zw(ndim,ndim)

      zw = matmul(zin,evecc)
      do jj=1,ndim
      do ii=1,ndim
         zout(ii,jj)=sum(dconjg(evecc(:,ii))*zw(:,jj))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
