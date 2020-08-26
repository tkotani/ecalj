!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_struc_def,only: s_lat,s_pot,s_mix,s_array,s_spec,s_site !,s_move,s_str,s_bz,s_ctrl ,s_ham
      use m_lmfinit, only:
     &     NULLI,NULLR,T,F,mxrecs,recln0,
     &     slabl_,v_spot,v_slat,v_slat2,v_smix,v_sarry,v_sarry2,v_sspec,v_ssite,
     &     v_ssite2,sstrnmix,sstrnsymg,
     &     symgaf !v_sbz,v_sctrl,v_sctrl2,v_sham,
      use m_hamindex,only: 
     &     iclasstaf,symops_af,ag_af,ngrpaf !we give them in this routine and used in sugw.F.
      use m_ext,only: getext,sname
      use m_globalvariables,only: globalvariables
      use m_work,only: rv_a_opos
      implicit none

      integer:: procid, master, mpipid, nproc,
     &     i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg, igets,a2vec,
     &     nfilin,nrecs,fopna,ic,ig,ib,ibas,nspec,inumaf,iprint,nbas,ioff
      real(8)::pwmode,pwemin,pwemax, bz_w,dval
      character*(recln0):: recrd(mxrecs)
      character(8) :: alabl
      character(5) :: charnum5
      logical:: fileexist,ipr10,lshow,lshowp,swtmp,cmdopt,ltet
      character:: strn*1000,outs*20,strn2*1000
      character prgnam*8, ext*100 
#if LMFGWD
      data prgnam /'LMFGWD'/
#else
      data prgnam /'LMF'/
#endif
#if MPI | MPIK
      include "mpif.h"
      integer fext
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
c      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog
      integer  nsize, id
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      if(id==0) write(*,*) 'mpi size=',nsize
#endif
      call show_programinfo(6)
      if(cmdopt('--h',3,0,outs)) then
        call lmhelp(prgnam) 
      endif
      stdo   = lgunit(1)
      master = 0
      procid = mpipid(1)
      nproc  = mpipid(0)
      globalvariables%stdo = dble(stdo) 
      i = dble(stdo)
      call finits(2,0,0,i)
!! Get extension of ctrl file to sname \in m_ext 
      call getext()
      ext = adjustl(sname)
      if (cmdopt('--sfill',7,0,strn)) call rx('no --sfill option in lm7K')
!! set log file for --mlog (not maintained well)
#if MPI | MPIK
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      shortname(procid)=trim(name)
      mlog = cmdopt('--mlog',6,0,strn)
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
#endif
!! Set iprint(); the top-level verbosity if specified from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k) ! --pr=xxx is read as k=xxx
        if (i .eq. 1) call setpr(k)
      endif
      if(procid/=master) call setpr(-100) !iprint() is negative except master
!! Abort if ctrl file is missing
      if (cmdopt('--input',6,0,strn)) goto 1111
      if (procid==master) then
         inquire(file='ctrl.'//trim(sname),exist=fileexist)
         if(.not.fileexist) call rx("No ctrl file found!") 
      endif
!! Read ctrl file into recrd
!   rdfile is too complicated to maintain. We will use simple reader.
!   We will give up math operation in ctrl.   
      nrecs = 0
      if(procid == master) then
         open(newunit=nfilin,file='ctrl.'//trim(sname)) ! = fopna('CTRL',-1,1)
         call findctrlstart(nfilin) ! ctrl is read from a tag 'ctrlstart' if it exist in ctrl file.
         alabl = '#{}% ct '
         call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs) 
         close(nfilin)
      endif
#if MPI | MPIK
      call mpibc1(nrecs,1,2,mlog,'main','nrecs')
      call MPI_BCAST(recrd,recln0*(nrecs+1),MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
#endif
      
!! Show readin ctrl file.      
      lshowp = cmdopt('--showp',7,0,strn)
      lshow  = cmdopt('--show',6,0,strn)
      if(procid==master.and.(lshow.or.lshowp) ) then
         write(stdo,"('---- preprocessed ctrl file -------')")
         do i = 1, nrecs
            if(trim(recrd(i))=='') cycle
            write(stdo,"(a)")'%% '//trim(recrd(i))//' %%'
         enddo
      endif
      write(stdo,"(a,' preprocesses ctrl file. nrecs=', i5)") trim(prgnam),nrecs
      if(lshowp) call cexit(0,1)
 1111 continue
      
!! rdctrl2 sets all initial condition to m_lmfinit module.
      call rdctrl2(recrd,recln0,nrecs,prgnam)
!! main part of lmf ------------------------------------
      call lmfp(prgnam,v_spot,v_slat, !v_sctrl,,v_sham
     .  v_smix,v_sspec,v_ssite,trim(sstrnmix)) !,v_smove,v_sstr,,v_sarry v_sbz,
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
      
!! finalization ----      
#if MPI | MPIK
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
