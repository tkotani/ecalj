module m_ctrl2ctrltoml
  ! Ensures ctrl.<sname>.toml exists before m_lmfinit's ReadCtrlp.
  !
  ! Dispatch:
  !   1. ctrl.<sname> exists  -> ALWAYS regenerate ctrl.<sname>.toml from
  !      ctrl.<sname> via ctrl2ctrltoml.py (LEGACY path is master; the
  !      regeneration must happen on every invocation because the user
  !      may pass different -v overrides between successive calls).
  !   2. ctrl.<sname>.toml only -> use it as-is (NEW TOML path; the .toml
  !      itself is master, e.g. authored by ctrlgenToml.py or hand-edited).
  !   3. Neither -> error.
  !
  ! Both decisions are logged so the active path is visible.
contains
  subroutine ConvertCtrl2CtrltomlByPython() bind(C)
    use m_args,    only: argall
    use m_ext,     only: sname
    use m_cmdpath, only: cmdpath
    use m_lgunit,  only: stdo
    use m_MPItk,   only: master_mpi
    implicit none
    character(512):: cmdl
    logical:: ctrlexist, tomlexist
    inquire(file='ctrl.'//trim(sname),         exist=ctrlexist)
    inquire(file='ctrl.'//trim(sname)//'.toml', exist=tomlexist)
    if (ctrlexist) then
       if (master_mpi) write(stdo,'(a)') &
            ' m_ctrl2ctrltoml: ctrl.'//trim(sname)// &
            ' present; regenerating ctrl.'//trim(sname)// &
            '.toml via ctrl2ctrltoml.py (LEGACY path)'
       cmdl = trim(cmdpath)//'ctrl2ctrltoml.py '//trim(argall)// &
            ' <ctrl.'//trim(sname)//' >ctrl.'//trim(sname)//'.toml'
       call system(cmdl)
       return
    endif
    if (tomlexist) then
       if (master_mpi) write(stdo,'(a)') &
            ' m_ctrl2ctrltoml: using existing ctrl.'//trim(sname)// &
            '.toml (NEW TOML path; no ctrl.'//trim(sname)//' present)'
       return
    endif
    call rx('m_ctrl2ctrltoml: neither ctrl.'//trim(sname)// &
         '.toml nor ctrl.'//trim(sname)//' was found')
  end subroutine ConvertCtrl2CtrltomlByPython
end module m_ctrl2ctrltoml
