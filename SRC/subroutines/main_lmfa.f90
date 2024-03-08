subroutine lmfa() bind(C)
  use m_args,only: m_setargs
  use m_ext,only:      m_Ext_init,sname
  use m_MPItk,only:    m_MPItk_init,nsize,master_mpi
  use m_lgunit,only:   m_lgunit_init, stdo,stdl
  use m_cmdpath,only: setcmdpath
  use m_lmfinit,only:m_Lmfinit_init
  use m_freeat,only:   Freeat
  use m_prgnam,only:   set_prgnam
  implicit none
  include "mpif.h" 
  logical:: cmdopt0
  character:: aaa*512
  character(8) :: prgnam='LMFA', charext
  integer:: ierr
  call m_setargs()
  call set_prgnam(prgnam)
  call m_MPItk_init() 
  call m_ext_init()  ! Get sname, e.g. trim(sname)=si of ctrl.si
  call m_lgunit_init()
  if(nsize/=1) call rx('Current lmfa is only for single core')
  aaa= '=== START LFMA ==='
  if(master_mpi) call show_programinfo(stdo)
  if(master_mpi) write(stdo,"(a)") trim(aaa)
  if(master_mpi) write(stdl,"(a)") trim(aaa)
  if(master_mpi) write(stdo,*) 'mpisize=',nsize
  if(master_mpi) write(stdl,*) 'mpisize=',nsize
  if(master_mpi) call setcmdpath()
  if(cmdopt0('--help')) then  !help and quit
     call m_lmfinit_init(prgnam) ! show help and quit for --input
     call rx0('end of help mode')
  endif
  call ConvertCtrl2CtrlpByPython(prgnam)
  if(cmdopt0('--quit=ctrlp')) call rx0('--quit=ctrlp')
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)
  call m_lmfinit_init(prgnam) ! Computational settings.
  call Freeat()  !Spherical atom calculation
  write(stdo,*)"OK! end of "//trim(prgnam)//" ======================"
end subroutine lmfa
