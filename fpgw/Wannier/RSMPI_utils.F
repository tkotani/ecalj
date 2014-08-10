      integer function get_non_blank_rsmpi(bufsize,buf)
      integer :: bufsize
      character :: buf(bufsize)

      integer :: i

      do i=bufsize,1,-1
         if (buf(i).ne.'') then
            get_non_blank_rsmpi = i
            return
         endif
      enddo
c 
      get_non_blank_rsmpi = 0
      return
      end

      subroutine set_index_rsmpi(N,m,N_local,N_local_index)
      implicit none

      integer :: N, m
      integer :: N_local(m), N_local_index(m,N/m + 1)

c local
      integer :: i,j,k

      N_local(:) = N/m
      

      do i=1,mod(N,m)
        N_local(i) = N_local(i) + 1
      enddo

      N_local_index(:,:) = 0
      k=0

      do i=1,m
        do j=1,N_local(i)
          k = k + 1
          N_local_index(i,j) = k

        enddo
      enddo

      if (k .ne. N) then
        stop "Error in set_index(): k .ne. N"
      endif
      end subroutine set_index_rsmpi

c--------------------------------------
c RS: set pointer of VCCFP
c
      subroutine set_vcoul_rsmpi(ifvcfpout,iq_in)
      implicit none
      integer,intent(in) :: ifvcfpout,iq_in
c local
      integer :: ndummy1,ndummy2,iq,nn
      complex(8),allocatable :: vcoul(:,:)


      rewind ifvcfpout
      read(ifvcfpout) ndummy1, ndummy2

      do iq=1,iq_in - 1
         read(ifvcfpout) nn 
         allocate(vcoul(nn,nn))
         read(ifvcfpout) vcoul(1:nn,1:nn)
         deallocate(vcoul)
      enddo
      end subroutine set_vcoul_rsmpi

