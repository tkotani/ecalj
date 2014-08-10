      program read_mmn
      implicit real*8(a-h,o-z)
      real(8),allocatable :: bb(:,:),wbb(:),qbz(:,:),wbz(:)
     &                     , r2n(:)
      complex(8),allocatable :: mmn(:,:,:,:)

      ifi = 99
      open(ifi,file='mmn.up')

      read(ifi,*)
      read(ifi,*)nwf,nqbz,nbb
      allocate( mmn(nwf,nwf,nbb,nqbz),bb(3,nbb),wbb(nbb),
     &          qbz(3,nqbz),wbz(nqbz) )
      read(ifi,*)
      read(ifi,*)mmn
      read(ifi,*)
      read(ifi,*)bb
      read(ifi,*)
      read(ifi,*)wbb
      read(ifi,*)
      read(ifi,*)qbz
      read(ifi,*)
      read(ifi,*)wbz

c calculate <n0| r^2 | n0>
      allocate( r2n(nwf) )
      r2n = 0d0
      do in = 1,nwf
      do iq  = 1,nqbz
      do ibb = 1,nbb
         r2n(in) = r2n(in) + 
     &            (1d0 - abs(mmn(in,in,ibb,iq))**2
     &          + (dimag(log(mmn(in,in,ibb,iq))))**2 )
     &          * wbb(ibb) * wbz(iq)
      enddo ! ib
      enddo ! iq
      enddo ! in

      close(ifi)

      end
