module m_xsfformat
  private:: wrt_pos_xyz
contains
  subroutine wrt_xsf( &
       basename,vis_unit, &
       alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
       mesh,rini,rfin,phipw,phiaug,phitot, &
       natom,apos,nclass,iclass,zz )
    implicit none
    character(*),intent(in):: basename,vis_unit
    double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3)
    integer,intent(in) :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
    double precision,intent(in) :: q_wfn(3,nq_wfn)
    double complex,intent(in) :: &
         phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
         phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

    integer,intent(in) :: natom,nclass,iclass(natom)
    double precision,intent(in) :: apos(3,natom),zz(nclass)
    character(200):: filename
    integer:: ifile,ifile_handle
    integer :: isp,iq,ib,iband,i1,i2,i3,natomall
    double precision :: rtmp(3),r(3)
    double precision,parameter:: zero=0.0d0
    double precision :: Z(natom)
    integer:: ic,ia,nrange(3),idim, iimg,ix
    character(3)::ril
    write(6,*)'wrt_xsf: basename='//trim(basename)
    do ia=1,natom
       ic=iclass(ia)
       Z(ia)=zz(ic)
    enddo
    do ix=1,3
       write(6,*)'rini rfin=',ix,rini(ix),rfin(ix)
       nrange(ix)= max( -floor(rini(ix)), ceiling(rfin(ix)) )
    enddo
    write(*,*) 'natom nrange=',natom,nrange(1:3)
    !      do ia=1,natom
    !         write(6,*)'bbb=',ia,apos(:,ia)
    !      enddo
    call wrt_pos_xyz('query natom',ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, &
         natomall )
    write(filename,'(a,a)')  basename(:len_trim(basename)), '.xsf'
    write(*,*) 'open ',filename
    ifile=ifile_handle()
    open(ifile,file=filename,status='unknown')
    write(ifile,'(a)') '# wavefunction'
    write(ifile,'(a)') 'PRIMVEC'
    !      close(ifile)
    do i1=1,3
       write(ifile,'(3f20.10)') plat(:,i1)*alat
    enddo
    call wrt_pos_xyz('write',ifile,alat,plat, rini,rfin,natom,apos,Z,nrange, &
         natomall )
    write(ifile,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
    write(ifile,'(a)') basename
    do iimg=1,2
       if(iimg==1) ril='_Re'
       if(iimg==2) ril='_Im'
       do isp=1,nsp
          do iq=1,nq_wfn
             do ib=1,nband_wfn
                iband=bindx_wfn(ib)
                !        write(ifile,'(a,i1,a,i3.3,a,i3.3,a,i3.3,a,i1)')
                !     .    'isp',isp,'_iq',iq,'_ib',ib,'_',iband,'_ri',iimg
                write(ifile,'(a,i1,a,i3.3,a,i3.3,a,i3.3,a)') &
                     'BEGIN_DATAGRID_3D_isp',isp,'_iq',iq,'_ib',ib,'_',iband,ril
                write(ifile,'(3i5)') mesh(1:3)+1
                if (trim(vis_unit) == 'alat') then
                   r(:)= rini(:) *alat
                   write(ifile,'(3f20.5)') r(1:3) !origin
                   write(ifile,'(3f20.5)') alat*(rfin(1)-rini(1))
                   write(ifile,'(3f20.5)') alat*(rfin(2)-rini(2))
                   write(ifile,'(3f20.5)') alat*(rfin(3)-rini(3))
                elseif (trim(vis_unit) == 'abc') then
                   r(:) = plat(:,1)*rini(1)+ plat(:,2)*rini(2)+ plat(:,3)*rini(3)
                   r= r*alat
                   write(ifile,'(3f20.5)') r(1:3) !origin
                   write(ifile,'(3f20.5)') alat*plat(:,1)*(rfin(1)-rini(1))
                   write(ifile,'(3f20.5)') alat*plat(:,2)*(rfin(2)-rini(2))
                   write(ifile,'(3f20.5)') alat*plat(:,3)*(rfin(3)-rini(3))
                endif
                do i3=1,mesh(3)+1
                   do i2=1,mesh(2)+1
                      if (iimg == 1) then
                         write(ifile,200) &
                              (  real(phitot(i1,i2,i3,ib,iq,isp)), i1=1,mesh(1)+1)
                      else
                         write(ifile,200) &
                              (  imag(phitot(i1,i2,i3,ib,iq,isp)), i1=1,mesh(1)+1)
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
100 format(i6,4f20.10)
200 format(6E20.10)
  end subroutine wrt_xsf

  !--------------------------------------------------------------------
  subroutine wrt_pos_xyz( job, ifile,alat,plat, rini,rfin,natom,apos, Z, nrange, natomall )
    !! atoms in range rini:rfin is given. These are in fractional coordinate.
    !! nrange is allowed range in the unit of plat
    implicit none
    integer:: ifile
    character(*):: job
    double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3), apos(3,natom),Z(natom)
    integer,intent(in) :: natom,nrange(3)
    integer,intent(out) :: natomall
    integer:: natomx
    integer :: i,i1,i2,i3
    double precision :: v1(3),v2(3),aini(3),afin(3),eps
    double precision ,allocatable :: rall(:,:),zall(:)

    integer::ix
    real(8)::qlat(3,3)

    call minv33tp(plat,qlat)
    eps = 0.05d0
    !      aini = alat*(rini-eps)
    !      afin = alat*(rfin+eps)

    write(6,*)' xxxxxxxxx nrange=',nrange
    write(6,*)' xxxxxxxxx rini=',rini
    write(6,*)' xxxxxxxxx rfin=',rfin

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
                !          call mymatvec(plat,v1,v2,3,3)
                !          v2(1:3)=alat*(v2(1:3)+apos(1:3,i)) !
                do ix=1,3
                   v2(ix)= v1(ix) + sum( apos(1:3,i)*qlat(1:3,ix) ) !atomic position in fractional coordinate
                enddo
                if(  rini(1)-eps<v2(1) .AND. v2(1)<rfin(1)+eps .AND. &
                     rini(2)-eps<v2(2) .AND. v2(2)<rfin(2)+eps .AND. &
                     rini(3)-eps<v2(3) .AND. v2(3)<rfin(3)+eps ) then
                   !          if ( (v2(1).ge.aini(1).and.v2(1).le.afin(1))
                   !     &    .and.(v2(2).ge.aini(2).and.v2(2).le.afin(2))
                   !     &    .and.(v2(3).ge.aini(3).and.v2(3).le.afin(3)) ) then
                   natomall = natomall + 1
                   call mymatvec(plat,v1,v2,3,3)
                   rall(1:3,natomall) = alat*(v2(1:3)+apos(1:3,i))
                   zall(natomall) = Z(i)
                endif
             enddo
          enddo
       enddo
    enddo

    !      if (job.eq.'write' .or. job.eq.'output') then
    !      do i=1,natomall
    !         write(ifile,"(i6,4E20.12)") int(zall(i)),zall(i),rall(1:3,i)
    !      enddo
    !      endif
    if (job == 'write' .OR. job == 'output') then
       write(ifile,'(a)') 'ATOMS'
       do i=1,natomall
          write(ifile,"(i3,4F20.5)") int(zall(i)),rall(1:3,i)
       enddo
    endif
    deallocate(rall,zall)
  end  subroutine wrt_pos_xyz

end module m_xsfformat


