      program aaa
      real(8):: bb(3,3,3),cc(3,3)
      do i=1,3
      do j=1,3
      do k=1,3
        bb(i,j,k)=100*i+10*j+k
      enddo
      enddo
      enddo
      cc=sum(bb(:,:,:),dim=3)
      do i=1,3
      do j=1,3
         print *,i,j,cc(i,j),cc(i,j)-sum(bb(i,j,1:3))
      enddo
      enddo
      end
