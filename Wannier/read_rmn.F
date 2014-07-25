      program read_rmn
      implicit real*8(a-h,o-z)
      real(8),allocatable :: rt(:,:)
      complex(8),allocatable :: rmn(:,:,:,:)

      ifi = 99
      open(ifi,file='rmn.up')
      read(ifi,*)
      read(ifi,*)nwf,nr
      allocate(rt(3,nr),rmn(3,nwf,nwf,nr))
c rt(x,cell) = R(x,cell)
      read(ifi,*)
      read(ifi,*)rt
c rmn(x,m,n,cell) = <0m | r_x | nR>
      read(ifi,*)
      read(ifi,*)rmn

      write(*,*)rmn

      end
