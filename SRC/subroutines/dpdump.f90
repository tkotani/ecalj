subroutine dpdump(array,length,ifile)
  !     - Binary I/O of an array
  integer:: length,ifile
  double precision :: array(length)
  if (ifile > 0) read(ifile) array
  if (ifile < 0) write(-ifile) array
end subroutine dpdump
subroutine dpsdmp(array,n1,n2,ifile)
  !- Binary I/O of an array segment
  integer :: n1,n2,ifile,length
  double precision :: array(n2)
  length = n2-n1+1
  if (length > 0) call dpdump(array(n1),length,ifile)
end subroutine dpsdmp
logical function lddump(array,length,ifile)
  !- Binary I/O of an array, returning T if I/O without error or EOF
  !     implicit none
  integer :: length,ifile
  double precision :: array(length),xx,yy
  lddump = .true.
  if (ifile > 0) then
     yy = array(length)
     !       (some random number)
     xx = -1.9283746d0*datan(1d0)
     array(length) = xx
     read(ifile,end=90,err=91) array
     if (xx /= array(length)) return
     array(length) = yy
     goto 90
90   continue
91   continue
     lddump = .false.
  else
     write(-ifile) array
  endif
end function lddump
