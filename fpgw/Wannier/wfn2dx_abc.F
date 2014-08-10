      subroutine wfn2dx_abc(alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     mesh,phipw,phiaug,phitot)
      implicit none
c input
      double precision :: alat,plat(3,3)
      integer :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
      double precision :: q_wfn(3,nq_wfn)
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

c local
      integer :: isp,iq,ib,iband,ifile,i1,i2,i3
      double precision :: rtmp(3),r(3)
      character*(23) :: fn

      ifile=100
      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        iband=bindx_wfn(ib)
        write(fn,"(a4,i1,a1,i4.4,a1,i4.4,a8)")
     &       'phis',isp,'q',iq,'b',iband,'.general'

        write(*,*) 'open ',fn
        open(ifile,file=fn)

        write(ifile,"(a,i2)") '#isp=',isp
        write(ifile,"(a,i5,3f10.4)") '#iq_wfn,q=',iq,q_wfn(1:3,iq)
        write(ifile,"(a,2i5)") '#ib,bindx=',ib,iband
        write(ifile,"(a)") 'header = marker "DATA\n"'
        write(ifile,"(a,i4,a,i4,a,i4)")
     &       'grid = ', mesh(1)+1,' x ',mesh(2)+1,' x ',mesh(3)+1
        write(ifile,"(a)") 'format = ascii'      
        write(ifile,"(a)") 'interleaving = field'
        write(ifile,"(a)") 'majority = row'
        write(ifile,"(a,a,a)") 'field = locations, ',
     &       'Re_phipw, Im_phipw, Re_phiaug, ',
     &       'Im_phiaug, Re_phitot, Im_phitot, |phitot|'
        write(ifile,"(a,a)") 'structure = 3-vector,',
     &       ' scalar, scalar, scalar, scalar, scalar, scalar, scalar'
        write(ifile,"(a,a)") 'type = float,',
     &       'float, float, float, float, float, float, float'
        write(ifile,"(a,a,a)") 'dependency = positions,',
     &       ' positions, positions,',
     &       ' positions, positions, positions, positions, positions'
        write(ifile,"(a)") 'end'
        write(ifile,"(a)") 'DATA'
        do i1=1,mesh(1)+1
          do i2=1,mesh(2)+1
            do i3=1,mesh(3)+1
              rtmp(1)=(i1-1)/dble(mesh(1))
              rtmp(2)=(i2-1)/dble(mesh(2))
              rtmp(3)=(i3-1)/dble(mesh(3))
!              call mymatvec(plat,rtmp,r,3,3)
c              r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
          r(:) = plat(:,1)*rtmp(1)+plat(:,2)*rtmp(2)+plat(:,3)*rtmp(3)
              r(1:3)=alat*(r(1:3))

              write(ifile,"(10f20.14)") r(1:3),
     &             phipw(i1,i2,i3,ib,iq,isp),
     &             phiaug(i1,i2,i3,ib,iq,isp),
     &             phitot(i1,i2,i3,ib,iq,isp),
     &             abs(phitot(i1,i2,i3,ib,iq,isp))
            enddo
          enddo
        enddo
        close(ifile)
      enddo
      enddo
      enddo
      end subroutine wfn2dx_abc
ccccccccccccccccccccccccccccccccccc
