subroutine sopen(output) bind(c)
  character(1):: output(*)
  character(1024):: convcchar
  open(6,file=trim(convcchar(output)))
end subroutine sopen
subroutine sclose() bind(c)
  close(6)
end subroutine sclose
