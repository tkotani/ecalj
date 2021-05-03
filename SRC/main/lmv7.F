!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_lmfinit, only:  mxrecs,recln0,rdctrl2
      use m_ext,only: getext,sname
      implicit none
      integer:: procid, master, mpipid, nproc,
     &     i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg, igets,a2vec,
     &     nfilin,nrecs,fopna,ic,ig,ib,ibas,nspec,inumaf,iprint,nbas,ioff
      real(8)::pwmode,pwemin,pwemax, bz_w,dval
      character*(recln0):: recrd(mxrecs)
      character(8) :: alabl
      character(5) :: charnum5
      logical:: fileexist,ipr10,lshow,lshowp,swtmp,cmdopt,ltet,llmfgw=.false.
      character:: strn*1000,outs*20
      character prgnam*8, ext*100 
#if LMFGWD
      data prgnam /'LMFGWD'/
#else
      data prgnam /'LMF'/
#endif
#if MPIK
      include "mpif.h"
      integer fext
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
#endif
      call show_programinfo(6)
      if(cmdopt('--h',3,0,outs)) then
        call lmhelp(prgnam) 
      endif
      stdo   = lgunit(1)
      master = 0
      procid = mpipid(1)
      nproc  = mpipid(0)
      call finits()!2,0,0,i)
      call getext() !return sname, e.g. trim(sname)=si of ctrl.si
      ext = adjustl(sname)      ! Get extension of ctrl file to sname \in m_ext
!!
!! === New pdos mode (use --mkprocar and --fullmesh together). ===
!! We use all k points (--fullmesh), instead of using crystal symmetry.
      if( cmdopt('--writepdos',11,0,strn) ) then
         write(6,*) '... Doing writepdos mode. Wait a while ...'
         write(6,*) '... See job_pdos about how to call --writepdos mode'
         call writepdos(trim(sname))
         call rx0('done: end of --writepdos mode.')
      endif
!! sawada's simple mode      
      if( cmdopt('--wdsawada',10,0,strn) ) then
         write(6,*) '... write Dos from tetraf.dat and eigenf.dat. '
         write(6,*) '...  eigenf.dat is for qplistf.dat '
         call writedossawada()
         call rx0('done: end of --wdsawada mode.')
       endif
!!      
      if (cmdopt('--input',6,0,strn)) goto 1111
      nrecs = 0
      if (procid==master) then
         inquire(file='ctrl.'//trim(sname),exist=fileexist)
         if(.not.fileexist) call rx("No ctrl file found!") 
!! Read ctrl file into recrd
!   rdfile is too complicated to maintain. We will use simple reader.
!   We will give up math operation in ctrl.   
         open(newunit=nfilin,file='ctrl.'//trim(sname)) ! = fopna('CTRL',-1,1)
         call findctrlstart(nfilin) ! ctrl is read from a tag 'ctrlstart' if it exist in ctrl file.
         alabl = '#{}% ct '
         call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs) !read control file into recrd
         close(nfilin)
      endif
#if MPIK
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      shortname(procid)=trim(name)
      mlog = cmdopt('--mlog',6,0,strn) !! set log file for --mlog (not maintained well)
      if (procid .eq. master) then
         call headl2(prgnam,0,stdo)
         i = fextg(ext)
      endif
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      if (procid .eq. master) then
        call gettime(datim)
        if (mlog) i = fopn('MLOG')
        call poseof(fopn('LOG'))
        if (mlog) then
          write(lgunit(3),"(a)")' lmf '//datim//' Process '
     &    //trim(charnum5(procid))//' of '//trim(charnum5(nproc))//' on '
     .    //trim(shortname(procid)) //' is master'
        endif
      else !for procid !=0
        call strcat(ext,20,' ','_',1,' ',i) !mpi file extention
        call bin2a(' ',0,0,procid,2,0,20,ext,i)
        ierr = fext(ext(1:i+1))
        if (mlog) ierr = fopn('MLOG')
        ierr = fextg(ext)
        call gettime(datim)
        if (mlog) then
          write(lgunit(3),"(a)")' lmf '//datim//' Process '
     &    //trim(charnum5(procid))//' of '//trim(charnum5(nproc))//' on '
     .    //trim(shortname(procid))
     .    //' file extension is '//ext(2:i+1)
        endif
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      call mpibc1(nrecs,1,2,mlog,'main','nrecs')
      call MPI_BCAST(recrd,recln0*(nrecs+1),MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
#endif
!! Set iprint(); the top-level verbosity if specified from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        print *,'outs',outs
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k) ! --pr=xxx is read as k=xxx
        if (i .eq. 1) call setpr(k)
      endif
      if(procid/=master) call setpr(-100) !iprint() is negative except master
!! Show readin ctrl file.      
      lshowp = cmdopt('--showp',7,0,strn)
      lshow  = cmdopt('--show', 6,0,strn)
      if(procid==master.and.(lshow.or.lshowp) ) then
         write(stdo,"('---- preprocessed ctrl file -------')")
         do i = 1, nrecs
            if(trim(recrd(i))=='') cycle
            write(stdo,"(a)")'%% '//trim(recrd(i))//' %%'
         enddo
      endif
      if(procid==master) write(stdo,"(a,' preprocesses ctrl file. nrecs=', i5)") trim(prgnam),nrecs
      if(lshowp) call cexit(0,1)
 1111 continue
!! Set starting data.  rdctrl2 sets all initial condition to m_lmfinit module.
      call rdctrl2(recrd,recln0,nrecs,prgnam)
!! Main part of lmf ------------------------------------
      if(prgnam =='LMFGWD') llmfgw=.true.
      call lmfp(llmfgw) !lmfp start from data in m_lmfinit
!! Finalization ----      
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
#if MPIK
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if ( procid .eq. master ) then
        call rx0(prgnam//' on '//trim(shortname(procid)))
      else
        call fexit(0,0,' ',0)
      endif
      call mpi_finalize(ierr)
#else
      call rx0s(prgnam) !for single core exit
#endif
      end program lmf
