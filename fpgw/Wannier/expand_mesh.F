      module m_expand_mesh 
      contains
      subroutine expand_mesh(plat,
     i    nq_wfn,nband_wfn,q_wfn,nsp,
     i    rini,rfin,
     i    mesh0,phipw0,phiaug0,phitot0,
     o    mesh, phipw,phiaug,phitot )
      
      implicit none
      integer,intent(in):: nq_wfn,nband_wfn,nsp
      real(8),intent(in)::  q_wfn(3,nq_wfn)
      integer:: rini(3), rfin(3)
      integer:: mesh0(3),mesh(3)
      double complex:: phipw0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,
     &         nband_wfn,nq_wfn,nsp)
      double complex:: phiaug0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,
     &         nband_wfn,nq_wfn,nsp)
      double complex:: phitot0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,
     &         nband_wfn,nq_wfn,nsp)

      double complex::phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,
     &         nband_wfn,nq_wfn,nsp)
      double complex::phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,
     &         nband_wfn,nq_wfn,nsp)
      double complex::phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,
     &         nband_wfn,nq_wfn,nsp)

      real(8):: pi ,cell_l(3),phase
      integer:: i,j,k,iq
      complex(8):: ci=(0.0d0,1.0d0), expikt
      real(8):: sss(3) ,plat(3,3)!
      pi=2.0d0*asin(1.0d0) 

c  phase factor, pi*sum_i(q_wfn(i)*cell_lbound(i))

      iqloop: do iq=1,nq_wfn

      do i=rini(1),rfin(1)-1
      do j=rini(2),rfin(2)-1
      do k=rini(3),rfin(3)-1
c        cell_l=(/i,j,k/)
        sss=plat(:,1)*i+plat(:,2)*j+plat(:,3)*k !
c similar to exp (ikT) 
        phase=sum(q_wfn(1:3,iq)*sss) !bugfix jul2014 cell_l(1:3) )
        expikt=exp(ci*2.0d0*pi*phase )  ! phase of the cell 
!bugfix july2014 mesh0* 1+(k-rini(1))*mesh0(1) ---> 1+(k-rini(3))*mesh0(3) 
        phipw(1+(i-rini(1))*mesh0(1):1+(1+i-rini(1))*mesh0(1), 
     .        1+(j-rini(2))*mesh0(2):1+(1+j-rini(2))*mesh0(2),
     .        1+(k-rini(3))*mesh0(3):1+(1+k-rini(3))*mesh0(3),:,iq,:) =
     .    expikt*  phipw0(1:1+mesh0(1),
     .                    1:1+mesh0(2),
     .                    1:1+mesh0(3),:,iq,:) 
        phiaug(1+(i-rini(1))*mesh0(1):1+(1+i-rini(1))*mesh0(1), 
     .        1+(j-rini(2))*mesh0(2):1+(1+j-rini(2))*mesh0(2),
     .        1+(k-rini(3))*mesh0(3):1+(1+k-rini(3))*mesh0(3),:,iq,:) =
     .    expikt*  phiaug0(1:1+mesh0(1),
     .                    1:1+mesh0(2),
     .                    1:1+mesh0(3),:,iq,:) 
        phitot(1+(i-rini(1))*mesh0(1):1+(1+i-rini(1))*mesh0(1), 
     .        1+(j-rini(2))*mesh0(2):1+(1+j-rini(2))*mesh0(2),
     .        1+(k-rini(3))*mesh0(3):1+(1+k-rini(3))*mesh0(3),:,iq,:) =
     .    expikt*  phitot0(1:1+mesh0(1),
     .                    1:1+mesh0(2),
     .                    1:1+mesh0(3),:,iq,:) 
      enddo
      enddo
      enddo

      enddo iqloop 

      end subroutine expand_mesh
      end module m_expand_mesh 
