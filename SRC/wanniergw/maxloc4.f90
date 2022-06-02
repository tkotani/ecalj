
!  write *rws and hrotr to ifh

subroutine write_hrotr(ifh, hrotr, &
     rws,irws,drws, &
     nwf,nrws )
  implicit none
  complex(8),intent(in) :: hrotr(nwf,nwf,nrws)
  integer,intent(in):: ifh, nwf,nrws
  real(8),intent(in) :: rws(3,nrws),drws(nrws)
  integer,intent(in) :: irws(nrws)

  integer:: i,j,k
  integer:: ir,im,in
  real(8):: rtmp,heps2

  write(ifh,*) 'nwf ',nwf
  write(ifh,*) 'nrws ',nrws

  write(ifh,*) '<rws>'
  do i=1,nrws
     write(ifh,*) i, rws(:,i), drws(i) , irws(i)
  enddo
  write(ifh,*)'</rws>'

  write(ifh,*) '<hrotr>'
  do i=1,nwf
     do j=1,nwf
        write(ifh,*)  i,j, 'i,j , the next line is hrotr(i,j,:)'
        write(ifh,'(10E20.10)')   hrotr(i,j,:)
     enddo
  enddo
  write(ifh,*) '</hrotr>'

  write(ifh,*) '<hrotr.abs>'
  do i=1,nwf
     do j=1,nwf
        write(ifh,*)  i,j, 'i,j , the next line is hrotr(i,j,:)'
        write(ifh,'(10E20.10)')  ( abs(hrotr(i,j,k)),k=1,nrws)
     enddo
  enddo
  write(ifh,*) '</hrotr.abs>'


end subroutine write_hrotr


subroutine read_hrotr(filename,nwf,nrws, &
     hrotr)
  use m_keyvalue,only: getkeyvalue

  implicit none
  character(*),intent(in):: filename
  integer,intent(in):: nwf,nrws
  complex(8):: hrotr(nwf,nwf,nrws)

  integer:: nwf_, nrws_
  character(10):: thisfunc='read_hrotr'
  integer:: ierror, ifh
  integer:: i,j,i_,j_
  character(120):: str

  write(*,*) 'reading ',filename
  call getkeyvalue(filename,'nwf',nwf_)
  call getkeyvalue(filename,'nrws',nrws_)

  ierror=0
  if ( nwf /= nwf_) then
     write(*,*) thisfunc,': data inconsistent nwf=', nwf, ' nwf(file)=',nwf_
     ierror=ierror+1
  endif

  if ( nrws /= nrws_ ) then
     write(*,*) thisfunc,': data inconsistent nrws=', nrws, ' nrws(file)=',nrws_
     ierror=ierror+1
  endif

  if (ierror /= 0) then
     goto 999
  endif

  call getkeyvalue(filename,'<hrotr>',unit=ifh,status=ierror,errstop='on')
  write(*,*) 'ifh,ierror=',ifh,ierror
  if (ierror == 0) then
     write(*,*) thisfunc,': failed to read <hrotr>'
     goto 999
  endif

  do i=1,nwf
     do j=1,nwf
        read(ifh,'(a120)')  str
        write(*,*) 'str=',str(:len_trim(str))
        read(str,*)i_,j_
        write(*,*) '1)',i_,j_
        read(ifh,'(10E20.10)')   hrotr(i,j,1:nrws)
        !          read(ifh,*)   hrotr(i,j,1:nrws)
        write(*,*) '2)',i_,j_
     enddo
  enddo

  close(ifh)

  return
999 write(*,*) 'abnormal exit'
  stop 'in read_hrotr'

end subroutine read_hrotr


subroutine make_hrotrcut( hrotr, &
     rws,irws,drws, &
     rcut,heps, &
     nwf,nrws, &
     hrotrcut )
  implicit none
  complex(8):: hrotr(nwf,nwf,nrws)
  real(8):: rws(3,nrws),drws(nrws)
  integer:: irws(nrws)
  real(8):: rcut, heps
  integer:: nwf,nrws
  complex(8),intent(out):: hrotrcut(nwf,nwf,nrws)
  integer:: ir,im,in
  real(8):: heps2,rtmp
  heps2=heps*heps
  hrotrcut=hrotr
  do ir = 1,nrws
     rtmp = dsqrt(sum(rws(:,ir)**2))  ! unit of alat
     write(*,"('cut:',i5,2f10.3)") ir,rtmp,rcut
     if (rtmp>rcut) then
        hrotrcut(:,:,ir)=0.0d0
     endif
     do im = 1,nwf
        do in = 1,nwf
           if ( im == in ) continue
           write(*,"('  ',2i5,'(',d13.5,',',d13.5,')',d13.5,d13.5)") &
                im,in,hrotr(im,in,ir) ,dble(hrotr(im,in,ir))**2 + dimag(hrotr(im,in,ir))**2, heps2
           if ( dble(hrotr(im,in,ir))**2 + dimag(hrotr(im,in,ir))**2 < heps2 ) then
              hrotrcut(im,in,ir)=0.0d0
           endif

        enddo
     enddo
  enddo
end subroutine make_hrotrcut





subroutine cart_to_fract (cart, fract_coord, qlat)
  implicit none
  real(8) :: qlat(3,3)
  real(8) :: cart(3), fract_coord(3)
  integer :: i
  do i=1,3
     fract_coord(i) = sum(cart(1:3)*qlat(1:3,i))
  enddo

  return

end subroutine cart_to_fract


subroutine write_hopping_output(is, hrotr, &
     rws,irws,alat,plat,qlat,pos,natom, &
     ibasiwf, nwf,nrws, spid, m_indx, l_indx, &
     nphix, iphi, ldim2)
  implicit none
  integer:: natom, is
  real(8) :: alat,plat(3,3),pos(3,natom),qlat(3,3)
  real(8),allocatable :: cart_coord(:,:), fract_coord(:,:)
  integer:: i,j,k, ldim2
  complex(8),intent(in) :: hrotr(nwf,nwf,nrws)
  integer,intent(in):: nwf,nrws
  real(8),intent(in) :: rws(3,nrws)
  integer,intent(in) :: irws(nrws)

  integer(4) :: ibasiwf(nwf), nphix, iphi(nphix,nwf)
  integer :: quantum_l, quantum_sym
  character(4) :: quantum_n, spin
  character(9), dimension(0:3,-5:5) :: orbital_sym
  character(8) :: spid(natom)
  integer(4) ::  m_indx(ldim2), l_indx(ldim2)
  integer:: ir, iwf, iwf1, iwf2, ifh1,ifh

  open(newunit=ifh, file='Hopping.dat')
  if(is==1)open(newunit=ifh1,file='Hopping.up')
  if(is==2)open(newunit=ifh1,file='Hopping.dn')

  write(ifh,"('Name : Need to be modified')")
  write(ifh1,"('Name : Need to be modified')")
  allocate (cart_coord(3,3))
  do i=1,3
     cart_coord(1:3,i) = alat*plat(1:3,i)*0.529177249
     write(ifh,"(3F16.9)") cart_coord(1:3,i)
     write(ifh1,"(3F16.9)") cart_coord(1:3,i)
  enddo
  write(ifh,"(2I6,I12)") nwf,nrws,nwf*nwf*nrws
  write(ifh,"(1A,I5)") "spin",is
  write(ifh1,"(2I6,I12)") nwf,nrws,nwf*nwf*nrws
  write(ifh1,"(1A,I5)") "spin",is
  deallocate (cart_coord)

  write(ifh,"(9(1A,7X))") &
       "leg","atom","n","l","sym","spin","x","y","z"
  write(ifh1,"(9(1A,7X))") &
       "leg","atom","n","l","sym","  spin","x  ","y  ","z  "


  orbital_sym(0,0)="s"
  orbital_sym(1,-1)="y"
  orbital_sym(1,0)="z"
  orbital_sym(1,1)="x"
  orbital_sym(2,-2)="xy"
  orbital_sym(2,-1)="yz"
  orbital_sym(2,0)="z2"
  orbital_sym(2,1)="xz"
  orbital_sym(2,2)="x2y2"
  orbital_sym(3,-3)="y(3x2-y2)" !see wikipedia!
  orbital_sym(3,-2)="xyz"
  orbital_sym(3,-1)="y(5z2-1)"
  orbital_sym(3,0)="z(5z2-1)"
  orbital_sym(3,1)="x(5z2-1)"
  orbital_sym(3,2)="z(x2-y2)"
  orbital_sym(3,3)="x(x2-3y2)"

  quantum_n = "--"

  allocate (fract_coord(3,nwf))
  do iwf=1,nwf
     if (is == 1) spin = "up"
     if (is == 2) spin = "down"

     call cart_to_fract(pos(1:3,ibasiwf(iwf)), fract_coord(1:3,iwf), &
          qlat)

     write(ifh,"(I3,1A16,A6,I6,2A12,3F10.5)") &
          iwf, spid(ibasiwf(iwf)), quantum_n, l_indx(iphi(1,iwf)), &
          orbital_sym(l_indx(iphi(1,iwf)),m_indx(iphi(1,iwf))), &
          spin, fract_coord(1:3,iwf)

     write(ifh1,"(I3,1A16,A6,I6,2A12,3F10.5)") &
          iwf, spid(ibasiwf(iwf)), quantum_n, l_indx(iphi(1,iwf)), &
          orbital_sym(l_indx(iphi(1,iwf)),m_indx(iphi(1,iwf))), &
          spin, fract_coord(1:3,iwf)

  end do
  deallocate (fract_coord)

  allocate (fract_coord(3,nrws))
  do ir=1,nrws
     call cart_to_fract(rws(:,ir), fract_coord(1:3,ir), &
          qlat)

     do iwf1=1,nwf
        do iwf2=1,nwf
           write(ifh,"(1x, 3i6, 3f12.6, 2x,2i4,3f12.6)") &
                nint(fract_coord(:,ir)), rws(:,ir), iwf1, iwf2, &
                hrotr(iwf1,iwf2,ir)*13.605698066
           write(ifh1,"(1x, 3i6, 3f12.6, 2x,2i4,3f12.6)") &
                nint(fract_coord(:,ir)), rws(:,ir), iwf1, iwf2, &
                hrotr(iwf1,iwf2,ir)*13.605698066

        enddo
     enddo
     !        do iwf1=1,nwf
     !          do iwf2=1,nwf
     !            write(ifh,"(1x, 6f12.6, 2i4,3f12.6)")
     !     &             rws(:,ir), fract_coord(:,ir), iwf1, iwf2,  hrotr(iwf1,iwf2,ir)
     !          enddo
     !        enddo
  enddo
  deallocate (fract_coord)
  close(ifh1)
  close(ifh)
end subroutine write_hopping_output


