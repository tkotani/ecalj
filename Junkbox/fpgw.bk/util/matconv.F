      implicit none
      integer(4):: i,irow,icol
      character(10000) ::iline
      real(8),allocatable:: dmat(:,:)
      read(5,"(a)") iline
      i = index(iline,"%")   ;  if(i/=1) stop 'No % line at first'
      i = index(iline,"cols");  read(iline(i+4:),*) icol
      i = index(iline,"rows");  read(iline(i+4:),*) irow
c      print *,' irow icol=',irow,icol
      allocate(dmat(icol,irow))
      write(6,"('% rows ',i10, '   cols ',i10)") irow,icol
      read(5,*) dmat
      do i=1,irow
        write(6,"(100(f12.4,x))") dmat(:,i)
      enddo
      end
