!! == Main program for lmf part. Now we use MPI for all program, but lmfa is only for single core.
!  lmfa,lmchk,lmf-MPIK,lmfgw-MPIK are from this code
!     We use module-based programing.
!     In principle, all the data are generared and stored in some modules with 'protection'.
!     We can not modify data in a module by other modules.
!      (but this module-based programming is still in progress). 2021may t.kotani
      program lmf
      use m_lmfinit,only:m_Lmfinit_init
      use m_ext,only:    m_Ext_init,sname
      use m_lattic,only: m_Lattic_init
      use m_mkqp,only:   m_Mkqp_init
      use m_rdfiln,only: m_Rdfiln_init
      use m_MPItk,only:  m_MPItk_init,m_MPItk_finalize,nsize,master_mpi
      implicit none
      integer:: k,lgunit,stdo,stdl
      logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2
      character:: outs*20,aaa*512,sss*128,uuu*128
      character prgnam*8
      character*8 charext
      integer:: fext,iarg
#if LMFGWD
      data prgnam /'LMFGWD'/
#elif LMFA
      data prgnam /'LMFA'/
#elif LMCHK
      data prgnam /'LMCHK'/
#else
      data prgnam /'LMF'/
#endif
      call m_ext_init()        ! Get sname, e.g. trim(sname)=si of ctrl.si
      call m_MPItk_init(prgnam)
      stdo = lgunit(1)          !standard output stdo=6 usually
      stdl = lgunit(2)
      aaa=''
      do iarg=1,iargc()
         call getarg(iarg,sss)
         aaa=aaa//' '//trim(sss)
      enddo
      uuu=' ====================================='
      aaa= '===START '//trim(prgnam)//'  '//trim(aaa)//trim(uuu)
      if(master_mpi) write(stdo,"(a)") trim(aaa)
      if(master_mpi) write(stdl,"(a)") trim(aaa)
      if(master_mpi) write(stdo,*) 'mpisize=',nsize
      if(master_mpi) write(stdl,*) 'mpisize=',nsize
      
!! Set iprint(); the verbosity from cmd line      
      if(cmdopt2('--pr',outs).or.cmdopt2('-pr',outs)) then 
        read(outs(5:),*) k
        call Setpr(k) 
      endif
      if(.not.master_mpi) call setpr(-100) !iprint() is negative except master
!! help and quit      
      if (cmdopt0('--help')) then 
         call m_lmfinit_init(prgnam) ! show help and quit for --input
         call Rx0('end of help mode')
      endif   
!! New pdos mode (use --mkprocar and --fullmesh together). and exit ===
!! We use all k points (--fullmesh), instead of using crystal symmetry. See job_pdos
      if( Cmdopt0('--writepdos') ) then
         write(6,*) '... Doing writepdos mode. Wait a while ...'
         write(6,*) '... See job_pdos to know how to call --writepdos mode'
         call Writepdos(trim(sname))
         call Rx0('done: end of --writepdos mode.')
      endif
!! Sawada's simple mode and exit
      if( Cmdopt0('--wdsawada') ) then
         write(6,*) '... write Dos from tetraf.dat and eigenf.dat. '
         write(6,*) '...  eigenf.dat is for qplistf.dat '
         call Writedossawada()
         call Rx0('done: end of --wdsawada mode.')
      endif
      
!! Read ctrl.* into m_rdfiln
!! m_rdfiln is too complicated to maintain. We will use simple reader.
!! 'math operation and -v substitution' may be convenient, but probably we will separate them away.
      call m_rdfiln_init()
      
!! Read recrd of ctrl file into m_lmfinit
!!    Set all initial condition to the module m_lmfinit.
!!    it allow expression like ALAT=7.88*1.1
      call m_lmfinit_init(prgnam) !--> all initianl data go into m_lmfinit
      
!! Lattice and Symmetry
      call Latsym(prgnam) !---> all data go into m_lattic and m_mksym
!! q points info     
      if(trim(prgnam) == 'LMF') then
         call m_mkqp_init() !all data of BZ go into m_mkqp
      endif
      
!! Main part of lmf ------------------------------------
!! Computations are started from variables in m_lmfinit, m_lattic, m_mksym, m_mkqp
      if(prgnam=='LMCHK') then
        if(nsize/=1) call rx('Current lmchk is only for single core')
        call Lmaux()  !check crystal structure
      elseif(prgnam=='LMFA') then
        if(nsize/=1) call rx('Current lmfa is only for single core')
        call Freeat() !spherical atom calculation
      else              
         if(prgnam =='LMFGWD') llmfgw=.true.
         call Lmfp(llmfgw)   !lmfp: self-consistent calculation. lmfgw for GW driver
      endif
      
!! Finalization ----      
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
      call m_MPItk_finalize()
      if(master_mpi) write(6,"(a)") "OK! end of "//trim(prgnam)//" ======================"
      call exit(0)
      end program Lmf
