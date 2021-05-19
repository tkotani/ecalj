!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_lmfinit,only:  Rdctrl2,nbas, v_sspec,v_ssite,bz_lmet,ctrl_ldos
      use m_ext,only: Getext,sname
      use m_lattic,only: Lattic,rv_a_opos
      use m_mkqp,only: Mkqp
      use m_rdfiln,only: Readctrl
      implicit none
      integer:: procid=0, master=0, mpipid, nproc,ifr,ix,
     &     i,j,k,lgunit,i1mach,auxmod,lc,stdo, 
     &     nfilin,nrecs,ic,ig,ib,ibas,nspec,inumaf,iprint,ioff,irecrd
      real(8)::pwmode,pwemin,pwemax, bz_w,dval
      character(5) :: charnum5
      logical:: fileexist,ipr10,lshow,lshowp,swtmp,cmdopt,ltet,llmfgw=.false.,cmdopt0,cmdopt2
      character:: outs*20
      character prgnam*8, ext*100
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
      
#if MPIK
      include "mpif.h"
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10,allocatable:: shortname(:) 
      character*26 datim
      double precision starttime, endtime
      logical mlog
      integer  nsize, id
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      if(id==0) write(*,*) 'mpi size=',nsize
      allocate(shortname(0:nsize-1))
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      shortname(procid) = trim(name)
      procid = mpipid(1) !
      nproc  = mpipid(0) ! num of processors
      call Gettime(datim)
      call Finits() 
      call Getext()        ! Return sname, e.g. trim(sname)=si of ctrl.si
      mlog = cmdopt0('--mlog') !! set log for --mlog (not maintained well)
      if(mlog) write(lgunit(3),"(a)")' lmf '//datim//' Process '
     &  //trim(charnum5(procid))//' of '//trim(charnum5(nproc-1))//' on '//trim(shortname(procid)) 
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if(procid == master) ext=''
      if(procid /= master) ext='_'//trim(charext(procid))
      ierr = fext('.'//trim(sname)//trim(ext)) !set ext for fopen
#else
      call Finits() 
      call Getext()        ! Return sname, e.g. trim(sname)=si of ctrl.si
#endif
      stdo = lgunit(1) !standard output stdo=6 usually
      if(procid == master) then
        write(stdo,"(a)") ' --------------  START '//trim(prgnam)//'  ------------------'
      endif
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
      call Rdctrl2(prgnam) !recrd,recln,nrecs,prgnam) !prgnam dependence
!! Lattice and Symmetry
      call Latsym(prgnam)
      if(trim(prgnam) == 'LMF') then
         ltet = bz_lmet/=0 .or. ctrl_ldos/=0
         call Mkqp(ltet,.False.,1,-2)
      endif
!! We prepared all input for main part until now.
!! These are in module variables in m_lmfinit, m_lattic, m_mksym, m_mkqp
      
!! Main part of lmf ------------------------------------
      if(prgnam=='LMCHK') then
        if(nsize/=1) call rx('Current lmchk is only for single core')
        call Lmaux() 
      elseif(prgnam=='LMFA') then
        if(nsize/=1) call rx('Current lmfa is only for single core')
        call Freeat(v_sspec) 
      else  
         if(prgnam =='LMFGWD') llmfgw=.true.
         call Lmfp(llmfgw)      !lmfp start from data in m_lmfinit
      endif
      
!! Finalization ----      
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
#if MPIK
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if( procid == master ) call rx0(prgnam//' on '//trim(shortname(procid)))
      call Mpi_finalize(ierr)
#else
      call Rx0s(prgnam) !for single core exit
#endif
      end program Lmf
