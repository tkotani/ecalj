program lmchk
  use m_lmfinit,only:m_Lmfinit_init,nlibu,plbnd
  use m_ext,only:      m_Ext_init,sname
  use m_lattic,only:   m_Lattic_init
  use m_MPItk,only:    m_MPItk_init,m_MPItk_finalize,nsize,master_mpi
  use m_lmaux,only: lmaux
  use m_mksym,only: m_mksym_init
  use m_lgunit,only:   m_lgunit_init, stdo,stdl
  use m_cmdpath,only: setcmdpath
  implicit none
  logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2,writeham,lbin
  character:: outs*20,aaa*512,sss*128,uuu*128, prgnam*8, charext*8
  integer:: k,iarg,jobgw,iprint,nit1,ifi,ifile_handle
  data prgnam /'LMCHK'/
  call m_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
  call m_MPItk_init(prgnam)
  call m_lgunit_init()
  if(nsize/=1) call rx('Current lmchk is only for single core')
  aaa=''
  do iarg=1,iargc()
     call getarg(iarg,sss)
     aaa=aaa//' '//trim(sss)
  enddo
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
  call m_lmfinit_init(prgnam) ! Computational settings go into m_lmfinit
  if(cmdopt2('--pr=',outs)) then
     read(outs,*) k
     call Setprint(k)
  endif
  if( .NOT. master_mpi) call setprint(-100) !iprint() is negative except master
  call m_lattic_init() ! lattice setup (for ewald sum)
  call m_mksym_init(prgnam) ! symmetry go into m_lattic and m_mksym
  call lmaux()      !check crystal structure
  call m_MPItk_finalize()
  if(master_mpi) write(6,"(a)") "OK! end of "//trim(prgnam)//" ======================"
  call exit(0)
end program lmchk
