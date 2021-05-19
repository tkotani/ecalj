!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_lmfinit,only: m_lmfinit_init
      use m_ext,only:    m_ext_init,sname
      use m_lattic,only: m_Lattic_init
      use m_mkqp,only:   m_mkqp_init
      use m_rdfiln,only: m_rdfiln_init
      use m_MPItk,only:  m_MPItk_init,m_MPItk_finalize, procid,master,nsize
      implicit none
      integer:: k,lgunit,stdo 
      logical:: fileexist,cmdopt,llmfgw=.false.,cmdopt0,cmdopt2
      character:: outs*20
      character prgnam*8
      character*8 charext
      integer fext
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
      stdo = lgunit(1) !standard output stdo=6 usually
      if(procid==master) write(stdo,"(a)") ' ---------  START '//trim(prgnam)//'  -------------'
      call Show_programinfo(stdo)
      if(cmdopt2('--pr',outs).or.cmdopt2('-pr',outs)) then 
        read(outs(5:),*) k
        call Setpr(k) !Set iprint(); the verbosity from cmd line
      endif
      if(procid/=master) call setpr(-100) !iprint() is negative except master
!! help and quit      
      if (cmdopt0('--input')) then 
         call m_lmfinit_init(prgnam) ! show help and quit for --input
         call Rx0('end of --input(--help) mode')
      endif   
!! New pdos mode (use --mkprocar and --fullmesh together). and exit ===
!! We use all k points (--fullmesh), instead of using crystal symmetry.
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
!! Read ctrl.//sname into m_rdfiln
!! m_rdfiln is too complicated to maintain. We will use simple reader.
!!   'math operation and -v substitution' may be convenient, but probably we will separate away them
      call m_rdfiln_init()
!! Read recrd of ctrl file into m_lmfinit
!!    rdctrl2 sets all initial condition to the module m_lmfinit.
!!    rdctrl2 allow expression like ALAT=7.88*1.1
      call m_lmfinit_init(prgnam) !--> all initianl data go into m_lmfinit
!! Lattice and Symmetry
      call Latsym(prgnam) !---> all data go into m_lattic m_mksym
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
      call Rx0s(prgnam) !for single core exit
      end program Lmf
