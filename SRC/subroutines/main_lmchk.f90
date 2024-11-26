module m_lmchk
  contains
subroutine lmchk(commin) bind(C)
  use m_args,only: argall,m_setargs
  use m_ext,only:     m_Ext_init,sname
  use m_MPItk,only:   m_MPItk_init,nsize,master_mpi
  use m_lgunit,only:  m_lgunit_init, stdo,stdl
  use m_cmdpath,only: setcmdpath
  use m_lmfinit,only: m_Lmfinit_init,nlibu,plbnd
  use m_lattic,only:  m_Lattic_init
  use m_mksym,only:   m_mksym_init
  use m_lmaux,only:   lmaux
  implicit none
  include "mpif.h" 
  logical:: cmdopt0,cmdopt2
  character:: outs*20,aaa*512,sss*128
  integer,optional::commin
  integer:: iarg,k,ierr,comm
  character(32):: prgnam='LMCHK'
  comm=MPI_COMM_WORLD
  if(present(commin)) comm= commin
  call m_MPItk_init(comm) 
  call m_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
  call m_lgunit_init()
  if(nsize/=1) call rx('Current lmchk is only for single core')
  aaa=argall
  aaa= '=== START '//trim(prgnam)//'  '//trim(aaa)//' ==='
  if(master_mpi) write(stdo,"(a)") trim(aaa)
  if(master_mpi) write(stdl,"(a)") trim(aaa)
  if(master_mpi) write(stdo,*) 'mpisize=',nsize
  if(master_mpi) write(stdl,*) 'mpisize=',nsize
  if(master_mpi) call setcmdpath() !set self-command path
  if(cmdopt0('--help')) then  !help and quit
     call m_lmfinit_init(prgnam) ! show help and quit for --input
     call rx0('end of help mode')
  endif
  if(master_mpi) call ConvertCtrl2CtrlpByPython()
  call MPI_BARRIER( comm, ierr)
  call m_lmfinit_init(prgnam,comm) ! Computational settings.
  !call m_mksym_init()  !symmetry go into m_lattic and m_mksym
  if(cmdopt2('--pr=',outs)) then
     read(outs,*) k
     call Setprint(k)
  endif
  if( .NOT. master_mpi) call setprint(-100) !iprint() is negative except master
  call m_lattic_init() !lattice setup (for ewald sum)
  call m_mksym_init()  !symmetry go into m_lattic and m_mksym
  call lmaux()              !Main part. check crystal structure
  call rx0("OK! end of "//trim(prgnam)//" ======================")
end subroutine lmchk
end module m_lmchk
