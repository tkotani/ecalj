!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_lmfinit,only:  mxrecs,recln0,rdctrl2,  v_sspec !,v_ssite,slabl_
      use m_ext,only: getext,sname
      implicit none
      integer:: procid, master, mpipid, nproc,ifr,ix,
     &     i,j,k,lgunit,i1mach,auxmod,lc,stdo, 
     &     nfilin,nrecs,ic,ig,ib,ibas,nspec,inumaf,iprint,nbas,ioff,irecrd
      real(8)::pwmode,pwemin,pwemax, bz_w,dval
      character*(recln0):: recrd(mxrecs)
      character(8) :: alabl
      character(5) :: charnum5
      logical:: fileexist,ipr10,lshow,lshowp,swtmp,cmdopt,ltet,llmfgw=.false.,cmdopt0,cmdopt2
      character:: strn*1000,outs*20
      character prgnam*8, ext*100
      character*8 charext
      integer fext
      
c#if LMFRCD 
c      data prgnam /'LMFRCD'/ #preprocessor of ctrl file to recrd
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
#endif
      call show_programinfo(6)
      call getext()        ! Return sname, e.g. trim(sname)=si of ctrl.si
      ext = adjustl(sname) ! Get extension of ctrl file to sname \in m_ext
      stdo   = lgunit(1) !standard output stdo=6 usually
      master = 0
      procid = mpipid(1) !
      nproc  = mpipid(0) ! num of processors
      call finits() 
      if(cmdopt2('--pr',outs).or.cmdopt2('-pr',outs)) then! Set iprint(); the verbosity from cmd line
        read(outs(5:),*) k
        call setpr(k)
      endif
      if(procid/=master) call setpr(-100) !iprint() is negative except master
     
!! Show input info
      if (cmdopt0('--input')) then !('--help')) then
         call rdctrl2(recrd,recln0,nrecs,prgnam) !prgnam dependence
         call rx0('end of --input(--help) mode')
      endif   

!! New pdos mode (use --mkprocar and --fullmesh together). and exit ===
!! We use all k points (--fullmesh), instead of using crystal symmetry.
      if( cmdopt0('--writepdos') ) then
         write(6,*) '... Doing writepdos mode. Wait a while ...'
         write(6,*) '... See job_pdos to know how to call --writepdos mode'
         call writepdos(trim(sname))
         call rx0('done: end of --writepdos mode.')
      endif
      
!! Sawada's simple mode and exit
      if( cmdopt0('--wdsawada') ) then
         write(6,*) '... write Dos from tetraf.dat and eigenf.dat. '
         write(6,*) '...  eigenf.dat is for qplistf.dat '
         call writedossawada()
         call rx0('done: end of --wdsawada mode.')
      endif
      
!! Read ctrl file into recrd
!!  rdfile is too complicated to maintain. We will use simple reader.
!!   'math operation and -v substitution' may be convenient, but probably we will separate away them
      nrecs = 0
      if (procid==master) then
         inquire(file='ctrl.'//trim(sname),exist=fileexist)
         if(.not.fileexist) call rx("No ctrl file found! ctrl."//trim(sname)) 
         open(newunit=nfilin,file='ctrl.'//trim(sname)) 
         call findctrlstart(nfilin) ! if a tag 'ctrlstart' in ctrl, ctrl is read from the tag.
         alabl = '#{}% ct '
         call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs) !read ctrl into recrd
         close(nfilin)
      endif
      
!  recrd already replace foobar with foobar2 if we have an command-line option
!    -vfoobar=foobar2 option (replacement of %const defined in ctrl file).
! Following ctrlrecrd is just for check recrd contains simplified ctrl with -vfoobar replacement.
c      if(procid==master.and.prgnam=='LMFRCD') then
c         open(newunit=ifr,file='ctrlrcd.'//trim(sname))
c         ilen=0
c         do i=1,nrecs
c           if(len_trim(recrd(i))/=0) ilen=ilen+1
c         enddo
c         write(ifr,"(a)") ilen
c         do i=1,nrecs
c           if(len_trim(recrd(i))/=0) write(ifr,"(a)") trim(recrd(i))
c         enddo
c         close(ifr)
c      else   
c         open(newunit=ifr,file='ctrlrcd.'//trim(sname))
c         do i=1,nrecs
c           read(ifr,"(a)",end=1019) recrd(i)
c        enddo
c 1019   continue
c        close(ifr)
c     endif
      
!! Show or not show readin ctrl file.      
      lshowp = cmdopt0('--showp')
      lshow  = cmdopt0('--show')
      if(procid==master.and.(lshow.or.lshowp) ) then
         write(stdo,"('---- preprocessed ctrl file -------')")
         do i = 1, nrecs
            if(trim(recrd(i))=='') cycle
            write(stdo,"(a)")'%% '//trim(recrd(i))//' %%'
         enddo
      endif
      if(procid==master) write(stdo,"(a,' preprocesses ctrl file. nrecs=', i5)") trim(prgnam),nrecs
      if(lshowp) call cexit(0,1)

!! mpi setup      
#if MPIK
      allocate(shortname(0:nsize-1))
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      shortname(procid) = trim(name)
      if(procid == master) then
        write(stdo,"(a)") ' --------------  START '//trim(prgnam)//'  ------------------'
        ext=''
      else   
        ext='_'//trim(charext(procid))
      endif   
      ierr = fext('.'//trim(sname)//trim(ext)) !set ext for fopen
      call gettime(datim)
      mlog = cmdopt0('--mlog') !! set log for --mlog (not maintained well)
      if(mlog) write(lgunit(3),"(a)")' lmf '//datim//' Process '
     &  //trim(charnum5(procid))//' of '//trim(charnum5(nproc-1))//' on '//trim(shortname(procid)) 
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      call mpibc1( nrecs,1,2,mlog,'main','nrecs')
      call MPI_BCAST( recrd,recln0*(nrecs+1),MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
#endif
      
!! Read recrd of ctrl file into m_lmfinit
!!    rdctrl2 sets all initial condition to the module m_lmfinit.
!!    rdctrl2 allow expression like ALAT=7.88*1.1
      call rdctrl2(recrd,recln0,nrecs,prgnam) !prgnam dependence

!! Main part of lmf ------------------------------------
      if(prgnam=='LMCHK') then
        if(nsize/=1) call rx('Current lmchk is only for single core')
        call lmaux() 
      elseif(prgnam=='LMFA') then
        if(nsize/=1) call rx('Current lmfa is only for single core')
        call freeat(v_sspec)
      else  
         if(prgnam =='LMFGWD') llmfgw=.true.
         call lmfp(llmfgw)      !lmfp start from data in m_lmfinit
      endif
      
!! Finalization ----      
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
#if MPIK
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if( procid == master ) call rx0(prgnam//' on '//trim(shortname(procid)))
      call mpi_finalize(ierr)
#else
      call rx0s(prgnam) !for single core exit
#endif
      end program lmf
