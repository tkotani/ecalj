      program lmfa
      use m_lmfinit,only:m_Lmfinit_init
      use m_ext,only:      m_Ext_init,sname
      use m_rdfiln,only:   m_Rdfiln_init
      use m_MPItk,only:    m_MPItk_init,m_MPItk_finalize,nsize,master_mpi
      implicit none
      integer:: k,lgunit,stdo,stdl
      logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2
      logical:: writeham,lbin
      character:: outs*20,aaa*512,sss*128,uuu*128
      character prgnam*8
      character*8 charext
      integer:: iarg,jobgw,iprint,nit1,ifi,irs1x,ifile_handle
      data prgnam /'LMFA'/
      call M_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
      call M_MPItk_init(prgnam)
      if(nsize/=1) call rx('Current lmfa is only for single core')
      stdo = lgunit(1)          !standard output stdo=6 usually
      stdl = lgunit(2)
      aaa= '=== START LFMA ==='
      if(master_mpi) call show_programinfo(stdo)
      if(master_mpi) write(stdo,"(a)") trim(aaa)
      if(master_mpi) write(stdl,"(a)") trim(aaa)
      if(master_mpi) write(stdo,*) 'mpisize=',nsize
      if(master_mpi) write(stdl,*) 'mpisize=',nsize
      if(cmdopt0('--help')) then  !help and quit
         call M_lmfinit_init(prgnam) ! show help and quit for --input
         call Rx0('end of help mode')
      endif
      call M_rdfiln_init()        ! read ctrl file
      call M_lmfinit_init(prgnam) ! Computational settings. Set data from rdfiln.
      call Freeat()  !Spherical atom calculation
      inquire(file='RSEQ_ERROR',exist=fileexist) ! error file exist or not
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
      call m_MPItk_finalize()
      if(master_mpi) write(6,"(a)") "OK! end of "//trim(prgnam)//" ======================"
      call Exit(0)
      end program Lmfa

      include "show_programinfo.fpp"
