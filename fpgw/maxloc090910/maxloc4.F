c
c  write *rws and hrotr to ifh
c
      subroutine write_hrotr(ifh, hrotr,
     i                  rws,irws,drws,
     d                  nwf,nrws )
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
 

      subroutine read_hrotr(filename,nwf,nrws,
     o          hrotr)
       use keyvalue
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
      if ( nwf.ne.nwf_) then
         write(*,*) thisfunc,': data inconsistent nwf=', nwf, ' nwf(file)=',nwf_
         ierror=ierror+1
      endif

      if ( nrws .ne.nrws_ ) then
         write(*,*) thisfunc,': data inconsistent nrws=', nrws, ' nrws(file)=',nrws_
         ierror=ierror+1
      endif

      if (ierror.ne.0) then
         goto 999
      endif

      call getkeyvalue(filename,'<hrotr>',unit=ifh,status=ierror,errstop='on') 
      write(*,*) 'ifh,ierror=',ifh,ierror
      if (ierror.eq.0) then
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
c          read(ifh,*)   hrotr(i,j,1:nrws)
          write(*,*) '2)',i_,j_
        enddo
      enddo

      close(ifh)
      
      return 
 999     write(*,*) 'abnormal exit'
         stop 'in read_hrotr'

      end subroutine read_hrotr


      subroutine make_hrotrcut( hrotr,
     i                  rws,irws,drws,
     i                  rcut,heps,
     d                  nwf,nrws,
     o                  hrotrcut )
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
        write(*,*) 'cut:',ir,rtmp,rcut 
        if (rtmp>rcut) then
          hrotrcut(:,:,ir)=0.0d0
        endif

        do im = 1,nwf
          do in = 1,nwf
        
            if ( im.eq.in ) continue 
             write(*,*) im,in,hrotr(im,in,ir) ,dble(hrotr(im,in,ir))**2 + dimag(hrotr(im,in,ir))**2, heps2
            if ( dble(hrotr(im,in,ir))**2 + dimag(hrotr(im,in,ir))**2 < heps2 ) then
                hrotrcut(im,in,ir)=0.0d0
            endif

          enddo
        enddo

      enddo 

      end subroutine make_hrotrcut
 
