program main
  use m_args, only: m_setargs
  use m_ext,  only: m_ext_init
  use m_gwinit, only: gwinit_v2
  call m_setargs()
  call m_ext_init()    ! pick up sname so PB.<sname>.toml.tmp gets the right name
  call gwinit_v2()
end program main
