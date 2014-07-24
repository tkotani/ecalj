      module m_cubeformat
      private::  wrt_pos_xyz  
      contains
      subroutine wrt_cube(
     i     basename,
     i     alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     i     mesh,rini,rfin,phipw,phiaug,phitot,
     i     natom,apos,nclass,iclass,zz )
      implicit none
c input
      character(*),intent(in):: basename
      double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3)
      integer,intent(in) :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
      double precision,intent(in) :: q_wfn(3,nq_wfn)
      double complex,intent(in) :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

      integer,intent(in) :: natom,nclass,iclass(natom)
      double precision,intent(in) :: apos(3,natom),zz(nclass)

      character(200):: filename
      integer:: ifile=150

      integer :: isp,iq,ib,iband,i1,i2,i3,natomall
      double precision :: rtmp(3),r(3)
      double precision,parameter:: zero=0.0d0 

      double precision :: Z(natom)
      integer:: ic,ia,nrange(3),idim, iimg


      do ia=1,natom
        ic=iclass(ia)
        Z(ia)=zz(ic)
      enddo

      do ia=1,3
      nrange(ia)= max( -floor(rini(ia)), ceiling(rfin(ia)) ) 
      enddo 
      write(*,*) 'nrange=',nrange(1:3) 

      call wrt_pos_xyz('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
     o  natomall )
c      call wrt_pos_2('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
c     o  natomall )


      do iimg=1,2
      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        iband=bindx_wfn(ib)
        if (iimg.eq.1) then
        write(filename,"(a,a1,i1,a1,i4.4,a1,i4.4,a6)")
     &       basename(:len_trim(basename)), 's',
     &       isp,'q',iq,'b',iband,'r.cube'
        else
        write(filename,"(a,a1,i1,a1,i4.4,a1,i4.4,a6)")
     &       basename(:len_trim(basename)), 's',
     &       isp,'q',iq,'b',iband,'i.cube'
        endif
        write(*,*) 'open ',filename
        open(ifile,file=filename,status='unknown')
        write(ifile,'(a)') 'wavefunction'
        write(ifile,'(a,4I5)') 'isp,iq,ib,iband=',isp,iq,ib,iband

        r(:)= rini(:) *alat
        write(ifile,100) natomall, r(1:3) 

        r = (rfin-rini)*alat
        idim=1
        write(ifile,100) mesh(idim)+1,r(idim)/mesh(idim),zero,zero
        idim=2
        write(ifile,100) mesh(idim)+1,zero,r(idim)/mesh(idim),zero
        idim=3
        write(ifile,100) mesh(idim)+1,zero,zero,r(idim)/mesh(idim)

      call wrt_pos_xyz('write',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
     o  natomall )
c      call wrt_pos_2('write',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
c     o  natomall )

        do i1=1,mesh(1)+1
          do i2=1,mesh(2)+1
             if (iimg.eq.1) then
              write(ifile,200) 
     &           (  real(phitot(i1,i2,i3,ib,iq,isp)), i3=1,mesh(3)+1)
             else
              write(ifile,200) 
     &           (  imag(phitot(i1,i2,i3,ib,iq,isp)), i3=1,mesh(3)+1)
             endif
          enddo
        enddo

        close(ifile)

      enddo ! ib
      enddo ! iq
      enddo ! isp 

      enddo ! iimg
                
 100    format(i6,4f20.10)
 200    format(6f20.10)
      end  subroutine wrt_cube

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#if 0
      subroutine wrt_pos_2(
     i    job,  ! = 'write' or 'query'
     i    ifile,alat,plat, rini,rfin,natom,apos, Z, nrange ,
     o    natomall )
      implicit none
      integer:: ifile
      character(*):: job
      double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3),
     &     apos(3,natom),Z(natom)
      integer,intent(in) :: natom,nrange(3)

      integer,intent(out) :: natomall

      integer:: natomx
      integer :: i,i1,i2,i3
      double precision :: v1(3),v2(3),aini(3),afin(3),eps

      natomall=natom*(2*nrange(1)+1)*(2*nrange(2)+1)*(2*nrange(3)+1)

      if (job.eq.'write' .or. job.eq.'output') then

      do i=1,natom
        do i1=-nrange(1),nrange(1)
        do i2=-nrange(2),nrange(2)
        do i3=-nrange(3),nrange(3)
          v1(1)=dble(i1)
          v1(2)=dble(i2)
          v1(3)=dble(i3)
          call mymatvec(plat,v1,v2,3,3)
          v2(1:3)=alat*(v2(1:3)+apos(1:3,i))
          write(ifile,'(i6,4E20.12)') int(z(i)), z(i), v2(1:3)
        enddo
        enddo
        enddo
      enddo

      endif

      end 
#endif

      subroutine wrt_pos_xyz(
     i    job,  ! = 'write' or 'query'
     i    ifile,alat,plat, rini,rfin,natom,apos, Z, nrange ,
     o    natomall )
      implicit none
      integer:: ifile
      character(*):: job 
      double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3),
     &     apos(3,natom),Z(natom)
      integer,intent(in) :: natom,nrange(3)

      integer,intent(out) :: natomall

      integer:: natomx
      integer :: i,i1,i2,i3
      double precision :: v1(3),v2(3),aini(3),afin(3),eps
      double precision,allocatable :: rall(:,:),zall(:)

      eps = 0.05d0
      aini = alat*(rini-eps)
      afin = alat*(rfin+eps)

      natomx=natom*(2*nrange(1)+1)*(2*nrange(2)+1)*(2*nrange(3)+1)
      allocate(rall(3,natomx),zall(natomx))
      natomall = 0
      do i=1,natom
        do i1=-nrange(1),nrange(1)
        do i2=-nrange(2),nrange(2)
        do i3=-nrange(3),nrange(3)
          v1(1)=dble(i1)
          v1(2)=dble(i2)
          v1(3)=dble(i3)
          call mymatvec(plat,v1,v2,3,3)
          v2(1:3)=alat*(v2(1:3)+apos(1:3,i))
          if ( (v2(1).ge.aini(1).and.v2(1).le.afin(1))
     &    .and.(v2(2).ge.aini(2).and.v2(2).le.afin(2))
     &    .and.(v2(3).ge.aini(3).and.v2(3).le.afin(3)) ) then
             natomall = natomall + 1
             rall(1:3,natomall) = v2(1:3)
             zall(natomall) = Z(i)
          endif   
        enddo
        enddo
        enddo
      enddo

      if (job.eq.'write' .or. job.eq.'output') then
      do i=1,natomall
         write(ifile,"(i6,4E20.12)") int(zall(i)),zall(i),rall(1:3,i)
      enddo
      endif 

      deallocate(rall,zall)

      end  subroutinewrt_pos_xyz

      end module m_cubeformat

