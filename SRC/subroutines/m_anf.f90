! !> Antiferro conditions set by call anfcond
! module m_anf
!   implicit none
!   logical,protected:: laf !! - laf: antiferro switch
!   integer,allocatable,protected:: ibasf(:) !! - ibasf(ibas) specify AF pair atom.
! contains
!   subroutine anfcond()
!     implicit none
!     integer:: ifi, natom
!     open(newunit=ifi,file='MLOindex',form='unformatted')
!     read(ifi) natom
!     read(ifi)
!     read(ifi)
!     allocate(ibasf(natom))
!     read(ifi) laf,ibasf
!     close(ifi)
!     if(laf) write(6,"(a,100i4)") ' Antiferromode=',ibasf
!   end subroutine anfcond
! end module m_anf

