      subroutine cexit(pv,ps)
      implicit none
#if MPI|MPIK
      include 'mpif.h'
#endif
      integer:: pv,ps,i
      integer:: status,ierr
      if (ps.ne.0) then
#if MPI|MPIK
        if (pv.eq.0) then
          call MPI_finalized(status,ierr)
          if (status.eq.0) then
            call MPI_finalize(ierr)
          endif
        endif
#endif
        call exit(pv)
      endif
      end subroutine cexit

      subroutine locase(ps)
      character(*) ps
      integer::i,n,shift
      n=len_trim(ps)
      shift=-ichar('A')+ichar('a') 
      do i=1,n
        if ( ichar(ps(i:i)) >= ichar('A')
     .  .and. ichar(ps(i:i)) <= ichar('Z') ) then
          ps(i:i) = char( ichar(ps(i:i))+ shift )
        endif
      enddo
      end subroutine locase

