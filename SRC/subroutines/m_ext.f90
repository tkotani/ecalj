!> get extension of ctrl.foobar
module m_ext 
  character(512),public,protected::sname='temp'
  public:: m_ext_init
contains
  subroutine m_ext_init()
    logical :: master
    integer:: ifi,ipos,i,na
    character*256:: sss,s222,argv
    na= iargc()
    do i = 1, na
       call getarg( i, argv )
       if(argv(1:1)/='-') then
          sname=trim(argv)
          exit
       endif
    enddo
  end subroutine m_ext_init
end module m_ext
