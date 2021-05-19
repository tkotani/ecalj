!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_lmfinit,only: Rdctrl2,v_sspec,bz_lmet,ctrl_ldos
      use m_ext,only:    Getext,sname
      use m_lattic,only: Lattic
      use m_mkqp,only:   Mkqp
      use m_rdfiln,only: Readctrl
      use m_MPItk,only: MPItk_init,MPItk_finalize, procid,master,nsize
      implicit none
      integer:: ifr,ix,
     &     i,j,k,lgunit,i1mach,auxmod,lc,stdo, 
     &     nfilin,nrecs,ic,ig,ib,ibas,nspec,inumaf,iprint,ioff,irecrd
      character(5) :: charnum5
      logical:: fileexist,ipr10,lshow,lshowp,swtmp,cmdopt,ltet,llmfgw=.false.,cmdopt0,cmdopt2
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
      call Getext()        ! Return sname, e.g. trim(sname)=si of ctrl.si
      call MPItk_init(prgnam)
      call Finits() 
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
         call Rdctrl2(prgnam) !recrd,recln,nrecs,prgnam) !prgnam dependence
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
      call Readctrl()
!! Read recrd of ctrl file into m_lmfinit
!!    rdctrl2 sets all initial condition to the module m_lmfinit.
!!    rdctrl2 allow expression like ALAT=7.88*1.1
      call Rdctrl2(prgnam) !--> all initianl data go into m_lmfinit
!! Lattice and Symmetry
      call Latsym(prgnam) !---> all data go into m_lattic m_mksym
!! q points info     
      if(trim(prgnam) == 'LMF') then
         ltet = bz_lmet/=0 .or. ctrl_ldos/=0
         call Mkqp(ltet,.False.,1,-2) !all data go into m_mkqp
      endif
!! Main part of lmf ------------------------------------
!! Computations are started from variables in m_lmfinit, m_lattic, m_mksym, m_mkqp
      if(prgnam=='LMCHK') then
        if(nsize/=1) call rx('Current lmchk is only for single core')
        call Lmaux()  !check crystal structure
      elseif(prgnam=='LMFA') then
        if(nsize/=1) call rx('Current lmfa is only for single core')
        call Freeat(v_sspec) !spherical atom calculation
      else              
         if(prgnam =='LMFGWD') llmfgw=.true.
         call Lmfp(llmfgw)   !lmfp: self-consistent calculation. lmfgw for GW driver
      endif
!! Finalization ----      
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
#if MPIK
      call MPItk_finalize()
#else
      call Rx0s(prgnam) !for single core exit
#endif
      end program Lmf
