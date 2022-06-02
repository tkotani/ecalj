subroutine info(jpr,l1,l2,string,a1,a2)
  use m_lgunit,only:stdo
  !- Printout when ipr>=jpr, with some blank lines
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   jpr   :print with verbosity is at least abs(jpr)
  !i         :jpr<0 -> print only if master node (MPI)
  !i   l1    :number of blank lines to print out before string
  !i   l2    :number of blank lines to print out after string
  !i         :l2<0 => write without carriage return
  !i   string:string to print
  !o Outputs
  !o   string is printed on stdo
  !r Remarks
  !u Updates
  !u   24 Nov 02 New info2, info5, info8
  !u   26 Jun 02 Add arguments a1,a2 so info can write data
  !u   02 Nov 01 Use awrite instead of write statement
  !u             Special handling of l2=-1
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: jpr,l1,l2,a1,a2,a3,a4,a5,a6,a7,a8
  character*(*) string
  ! ... Local parameters
  integer :: i,iprint,recl,awrite,mpipid
  parameter (recl=1024)
  character*(recl) lstr
  entry info0(jpr,l1,l2,string)
  entry info2(jpr,l1,l2,string,a1,a2)
  entry info5(jpr,l1,l2,string,a1,a2,a3,a4,a5)
  entry info8(jpr,l1,l2,string,a1,a2,a3,a4,a5,a6,a7,a8)
  if (iprint() < iabs(jpr)) return
  if (jpr < 0) then
     i = mpipid(1)
     if (i /= 0) return
  endif
  do  10  i = 1, l1
     write(stdo,100)
10 enddo
  if (l2 >= 0) then
     call awrit8(string,' ',recl,stdo,a1,a2,a3,a4,a5,a6,a7,a8)
  else
     i = awrite(string,lstr,recl,0,a1,a2,a3,a4,a5,a6,a7,a8)
     call cwrite(lstr,0,i-1,0)
  endif
  do  20  i = 1, l2
     write(stdo,100)
20 enddo
100 format(1x)
  return
end subroutine info

subroutine cwrite(ps,i1,i2,newln)
  character(*) ps
  integer:: i1,i2,newln,i
  integer,save:: iii
  if (newln /= 0) then
     write(*,'(a,$)') ps(i1+1:i2+1)
  else
     write(*,'(a)')   ps(i1+1:i2+1)
  endif
end subroutine cwrite
