      program lmchk
      use m_lmfinit,only:m_Lmfinit_init,nlibu,plbnd
      use m_ext,only:      m_Ext_init,sname
      use m_lattic,only:   m_Lattic_init
      use m_rdfiln,only:   m_Rdfiln_init
      use m_MPItk,only:    m_MPItk_init,m_MPItk_finalize,nsize,master_mpi
      use m_lmaux,only: lmaux
      use m_mksym,only: m_mksym_init
      implicit none
      integer:: k,lgunit,stdo,stdl
      logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2
      logical:: writeham,lbin
      character:: outs*20,aaa*512,sss*128,uuu*128
      character prgnam*8
      character*8 charext
      integer:: iarg,jobgw,iprint,nit1,ifi,ifile_handle
      data prgnam /'LMCHK'/
      call M_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
      call M_MPItk_init(prgnam)
      if(nsize/=1) call rx('Current lmchk is only for single core')
      stdo = lgunit(1)          !standard output stdo=6 usually
      stdl = lgunit(2)
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
      call M_rdfiln_init()
      call M_lmfinit_init(prgnam) ! Computational settings go into m_lmfinit
      if(cmdopt2('--pr=',outs)) then 
        read(outs,*) k
        call Setprint(k)
      endif
      if(.not.master_mpi) call setprint(-100) !iprint() is negative except master
      call m_lattic_init() ! lattice setup (for ewald sum)
      call m_mksym_init(prgnam) ! symmetry go into m_lattic and m_mksym
      call Lmaux()      !check crystal structure
      call m_MPItk_finalize()
      if(master_mpi) write(6,"(a)") "OK! end of "//trim(prgnam)//" ======================"
      call Exit(0)
      end program Lmchk
