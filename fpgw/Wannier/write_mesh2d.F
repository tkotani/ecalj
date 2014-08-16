      subroutine write_mesh2d(nq_wfn,nband_wfn,nsp,
     i   rini, fini,
     i   mesh, phipw,phiaug,phitot)
      implicit none
      integer:: nq_wfn,nband_wfn,nsp
      integer:: mesh(3)
      complex(8),dimension(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn, nsp) ::
     .   phipw,phiaug,phitot 
      real(8):: rini(3),fini(3)

      integer:: i,ifile=100,j,k,iq
      real(8):: pos(3)

      k=2
      do iq=1,nq_wfn
      do j=1,mesh(2)+1
      do i=1,mesh(1)+1
      pos(1)= rini(1)+(real(i-1)/mesh(1))*(fini(1)-rini(1))
      pos(2)= rini(2)+(real(j-1)/mesh(2))*(fini(2)-rini(2))
      pos(3)= rini(3)+(real(k-1)/mesh(3))*(fini(3)-rini(3))
c      write(*,*) i,j,k, pos
      write(ifile+iq,'(20f10.5)') pos(1:3),
     .   phipw(i,j,k,1,iq,1), phiaug(i,j,k,1,iq,1),phitot(i,j,k,1,iq,1)
      enddo 
       write(ifile+iq,*) ' '
      enddo
      enddo

      end subroutine write_mesh2d

