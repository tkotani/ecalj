integer function iprint()!print verbose setting
!  use m_MPItk,only: master_mpi
  implicit none
  integer :: mpipid,procid,ierr,comm
  integer:: verbose_in,setprint,ix,set0,setprint0,vb
  integer,save:: verbose0=30,verbose=30
  include "mpif.h"
  comm=MPI_COMM_WORLD
  call MPI_COMM_RANK(comm, procid, ierr )
  iprint = verbose
  if(procid/=0) iprint=0 !write only at master node
  return
  entry setprint0(verbose_in)!base
  verbose0=verbose_in
  verbose=verbose0
  return
  entry set0()
  verbose=verbose0
  return
  entry setprint(ix)
  verbose=ix
  return
end function iprint
subroutine pshpr(vb) !temporary set iprint
  integer:: setprint,vb,i
  i=setprint(vb)
end subroutine pshpr
subroutine poppr()   !pop to iprint by setpr0
  integer:: set0,i
  i=set0()
end subroutine poppr
subroutine getpr(ix) 
  integer:: iprint,ix
  ix=iprint()
end subroutine getpr
subroutine setpr0(ix)
  integer:: setprint0,ix,i
  i=setprint0(ix)
endsubroutine
