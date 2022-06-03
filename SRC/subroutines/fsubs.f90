subroutine cexit(pv,ps)
  implicit none
  include 'mpif.h'
  integer:: pv,ps,i
  integer:: status,ierr
  if (ps /= 0) then
     if (pv == 0) then
        call MPI_finalized(status,ierr)
        if (status == 0) then
           call MPI_finalize(ierr)
        endif
     endif
     call exit(pv)
  endif
end subroutine cexit

subroutine locase(ps)
  character(*) ps
  integer::i,n,shift
  n=len_trim(ps)
  shift=-ichar('A')+ichar('a')
  do i=1,n
     if ( ichar(ps(i:i)) >= ichar('A') &
          .AND. ichar(ps(i:i)) <= ichar('Z') ) then
        ps(i:i) = char( ichar(ps(i:i))+ shift )
     endif
  enddo
end subroutine locase

