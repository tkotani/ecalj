      module m_xsfformat 
      private:: wrt_pos_xyz
      contains 
      subroutine wrt_xsf(
     i     basename,vis_unit,
     i     alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     i     mesh,rini,rfin,phipw,phiaug,phitot,
     i     natom,apos,nclass,iclass,zz )
      implicit none
c input
      character(*),intent(in):: basename,vis_unit
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
      print *,'--- goto wrt_xsf basename=',basename,natom
      do ia=1,natom
        ic=iclass(ia)
        Z(ia)=zz(ic)
      enddo
      do ia=1,3
         print *,'rini rfin=',ia,rini(ia),rfin(ia)
         nrange(ia)= max( -floor(rini(ia)), ceiling(rfin(ia)) ) 
      enddo 
      write(*,*) 'nrange=',nrange(1:3) 
      call wrt_pos_xyz('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
     o     natomall )      
      write(filename,'(a,a)')  basename(:len_trim(basename)), '.xsf'
      write(*,*) 'open ',filename
      open(ifile,file=filename,status='unknown')
      write(ifile,'(a)') '# wavefunction'
      write(ifile,'(a)') 'PRIMVEC'
      do i1=1,3 
         write(ifile,'(3f20.10)') plat(:,i1)*alat
      enddo 
      call wrt_pos_xyz('write',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange,
     o  natomall )
      write(ifile,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
      write(ifile,'(a)') basename 
      do iimg=1,2
      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        iband=bindx_wfn(ib)
c        write(ifile,'(a,i1,a,i3.3,a,i3.3,a,i3.3,a,i1)') 
c     .    'isp',isp,'_iq',iq,'_ib',ib,'_',iband,'_ri',iimg
        write(ifile,'(a,i1,a,i3.3,a,i3.3,a,i3.3,a,i1)') 
     .    'BEGIN_DATAGRID_3D_isp',isp,'_iq',iq,'_ib',ib,'_',iband,'_ri',iimg
        write(ifile,'(3i5)') mesh(1:3)+1
        if (trim(vis_unit).eq.'alat') then
           r(:)= rini(:) *alat
        else if (trim(vis_unit).eq.'abc') then
           r(:) = plat(:,1)*rini(1)+ plat(:,2)*rini(2)+ plat(:,3)*rini(3)
           r= r*alat 
        endif 
        write(ifile,'(3f20.5)')  r(1:3) 
        if (vis_unit.eq.'alat') then
           r = (rfin-rini)*alat
        else if (vis_unit.eq.'abc') then
           r = plat(:,1)*(rfin(1)-rini(1))+ plat(:,2)*(rfin(2)-rini(2))+
     .          plat(:,3)*(rfin(3)-rini(3))
           r=r*alat 
        endif
        idim=1
        write(ifile,'(3f20.5)') alat*plat(:,1)*(rfin(1)-rini(1)) !r(idim),zero,zero
        idim=2
        write(ifile,'(3f20.5)') alat*plat(:,2)*(rfin(2)-rini(2)) !zero,r(idim),zero
        idim=3
        write(ifile,'(3f20.5)') alat*plat(:,3)*(rfin(3)-rini(3)) !zero,zero,r(idim)
        do i3=1,mesh(3)+1
          do i2=1,mesh(2)+1
             if (iimg.eq.1) then
              write(ifile,200) 
     &           (  real(phitot(i1,i2,i3,ib,iq,isp)), i1=1,mesh(1)+1)
             else
              write(ifile,200) 
     &           (  imag(phitot(i1,i2,i3,ib,iq,isp)), i1=1,mesh(1)+1)
             endif
          enddo
        enddo
        write(ifile,'(a)') 'END_DATAGRID_3D'
      enddo ! ib
      enddo ! iq
      enddo ! isp 
      enddo ! iimg
      write(ifile,'(a)') 'END_BLOCK_DATAGRID_3D'
      close(ifile)
 100  format(i6,4f20.10)
 200  format(6E20.10)
      end  subroutine wrt_xsf
!--------------------------------------------------------------------
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
      write(ifile,'(a)') 'ATOMS' 
      do i=1,natomall
         write(ifile,"(i3,4F20.5)") int(zall(i)),rall(1:3,i)
      enddo
      endif 
      deallocate(rall,zall)
      end subroutine wrt_pos_xyz
      end module m_xsfformat 


