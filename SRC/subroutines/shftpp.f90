subroutine shftpp(nc,nlsp,pp,vold,vnew,oshft,nshft)
  !- Shift or undo shift of pp's by constant potential
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   nc,pp,ves:
  !r Remarks:   shift pp's as follows:
  !r     oshft\nshft   F                 T
  !r           ---------------------------------------
  !r       F  |   do nothing       shift by vnew
  !r       T  | shift by -vold   shift by vnew-vold
  !r           ---------------------------------------
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed Parameters
  logical :: oshft,nshft
  integer :: nc,nlsp
  double precision :: pp(6,nlsp,nc),vold(nc),vnew(nc)
  ! Local variables
  integer :: ic

  ! --- Shift enu and c for each class ---
  do  10  ic = 1, nc
     if (oshft) then
        call daxpy(nlsp,-1d0,vold(ic),0,pp(1,1,ic),6)
        call daxpy(nlsp,-1d0,vold(ic),0,pp(2,1,ic),6)
     endif

     if (nshft) then
        call daxpy(nlsp,1d0,vnew(ic),0,pp(1,1,ic),6)
        call daxpy(nlsp,1d0,vnew(ic),0,pp(2,1,ic),6)
     endif
10 enddo

end subroutine shftpp





