module m_ctrl2ctrlp
contains
subroutine ConvertCtrl2CtrlpByPython() bind(C)
  use m_args,only: argall
  use m_ext,only :sname        ! sname contains extension. foobar of ctrl.foobar
  use m_cmdpath,only:cmdpath
!  use m_lgunit,only: stdo,stdl
  implicit none
  character(512):: cmdl
  logical:: fileexist
  integer::ifi
  inquire(file='ctrl.'//trim(sname),exist=fileexist)
  if(.NOT.fileexist) call rx("No ctrl file found!! ctrl."//trim(sname))
  cmdl=trim(cmdpath)//'ctrl2ctrlp.py '//trim(argall)//'<ctrl.'//trim(sname)//' >ctrlp.'//trim(sname)
!  write(stdo,"(a)")'cmdl for python='//trim(cmdl)
  call system(cmdl) !See  results ctrlp.* given by ctrl2ctrl.py 
end subroutine ConvertCtrl2CtrlpByPython
endmodule
