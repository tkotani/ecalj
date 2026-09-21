module m_cmdopt_registry   ! stub for the standalone contour test
  logical, public, protected, save :: c0_WVR2ptRaxis = .false.
end module m_cmdopt_registry
subroutine rx(msg)
  character(*) :: msg
  write(*,*) 'rx: ', msg; stop 1
end subroutine rx
