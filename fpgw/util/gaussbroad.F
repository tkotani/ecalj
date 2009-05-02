      integer(4),parameter:: nx=10000,nd=1
      real(8):: omg(nx), dat(1:nd,nx),dat2(1:nd),sigm
      character*100 :: a
      read(5,"(a)") a
      do ix=1,nx
         read(5,*,end=889) omg(ix),dat(1:nd,ix)
      enddo
 889  continue
      sigm=0.6
      write(6,"(a)") a
      nnx=ix-1
      do ix =1,nnx
        dat2=0d0
        do ixx=1,nnx
          dat2(:) =dat2(:)
     &     + exp(-(omg(ix)-omg(ixx))**2/2d0/sigm**2) *dat(:,ixx) 
        enddo
        write(6,"(100f13.6)") omg(ix),dat2(:)
      enddo
      end
