program lmfa
  use m_lmfinit,only:m_Lmfinit_init
  use m_ext,only:      m_Ext_init,sname
  use m_MPItk,only:    m_MPItk_init,m_MPItk_finalize,nsize,master_mpi
  use m_lgunit,only:   m_lgunit_init, stdo,stdl
  use m_freeat,only:   Freeat
  use m_cmdpath,only: setcmdpath
  implicit none
  integer:: k
  logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2
  logical:: writeham,lbin
  character:: outs*20,aaa*512,sss*128,uuu*128
  character prgnam*8
  character(8) :: charext
  integer:: iarg,jobgw,iprint,nit1,ifi,irs1x
  data prgnam /'LMFA'/
  call m_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
  call m_MPItk_init(prgnam)
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
     call Rx0('end of help mode')
  endif
  call m_lmfinit_init(prgnam) ! Computational settings.
  call Freeat()  !Spherical atom calculation
  call m_MPItk_finalize()
  if(master_mpi) write(6,"(a)") "OK! end of "//trim(prgnam)//" ======================"
  call Exit(0)
end program Lmfa

include "show_programinfo.fpp"
