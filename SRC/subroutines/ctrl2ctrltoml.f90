module m_ctrl2ctrltoml
  ! Replaces the old m_ctrl2ctrlp module. Spawns ctrl2ctrltoml.py to convert
  ! ctrl.<sname> -> ctrl.<sname>.toml. The Fortran consumer (m_lmfinit) then
  ! loads the TOML via toml-f and reformats to recrd(:) for rval2.
contains
  subroutine ConvertCtrl2CtrltomlByPython() bind(C)
    use m_args,only: argall
    use m_ext,only :sname
    use m_cmdpath,only:cmdpath
    implicit none
    character(512):: cmdl
    logical:: fileexist
    inquire(file='ctrl.'//trim(sname),exist=fileexist)
    if(.NOT.fileexist) call rx("No ctrl file found!! ctrl."//trim(sname))
    cmdl = trim(cmdpath)//'ctrl2ctrltoml.py '//trim(argall)// &
         ' <ctrl.'//trim(sname)//' >ctrl.'//trim(sname)//'.toml'
    call system(cmdl)
  end subroutine ConvertCtrl2CtrltomlByPython
end module m_ctrl2ctrltoml
