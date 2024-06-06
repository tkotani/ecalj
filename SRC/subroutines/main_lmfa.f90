module m_lmfa
contains
  subroutine lmfa(commin) bind(C)
    use m_args,only: argall
    use m_ext,only:      m_Ext_init,sname
    use m_MPItk,only:    m_MPItk_init,nsize,master_mpi
    use m_lgunit,only:   m_lgunit_init, stdo,stdl
    use m_lmfinit,only:m_Lmfinit_init
    use m_freeat,only:   Freeat
    implicit none
    include "mpif.h" 
    logical:: cmdopt0
    character:: aaa*512
    character(8) :: prgnam='LMFA', charext
    integer:: ierr,info,nsizex
    integer,optional:: commin
    integer:: comm ,ifi
    comm=MPI_COMM_WORLD
    if(present(commin)) comm= commin
    call m_MPItk_init(comm) 
    call m_ext_init()  ! Get sname, e.g. trim(sname)=si of ctrl.si
    call m_lgunit_init()
    if(nsize/=1) call rx('Current lmfa is only for single core')
    if(master_mpi) call show_programinfo(stdo)
    aaa= '=== START LFMA ==='
    if(master_mpi) write(stdo,"(a)") trim(aaa)
    if(master_mpi) write(stdl,"(a)") trim(aaa)
    if(master_mpi) write(stdo,*) 'mpisize=',nsize
    if(master_mpi) write(stdl,*) 'mpisize=',nsize
    if(cmdopt0('--help')) then  !help and quit
       call m_lmfinit_init(prgnam) ! show help and quit for --input
       call rx0('end of help mode')
    endif
    open(newunit=ifi,file='save.'//trim(sname),position='append')
    write(ifi,"(a)")'Start '//trim(prgnam)//trim(argall)
    close(ifi)
    if(master_mpi) call ConvertCtrl2CtrlpByPython()
    if(cmdopt0('--quit=ctrlp')) call rx0('--quit=ctrlp')
    call MPI_BARRIER( comm, ierr)
    call m_lmfinit_init(prgnam,comm) ! Computational settings.
    call Freeat()  !Spherical atom calculation
    write(stdo,"(a)")"OK! end of "//trim(prgnam)//" ======================"
  end subroutine lmfa
end module m_lmfa
