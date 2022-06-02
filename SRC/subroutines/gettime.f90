subroutine gettime(datim)
  !     implicit none
  character datim*(*)
  datim = ' '
  call ftime(datim)
end subroutine gettime

