!! == Main program for lmf and lmfgw (single core or -dMPIK works).
      program lmf
      use m_struc_def,only: s_bz,s_ctrl,s_lat,s_ham,s_pot,s_mix,s_array,s_spec,s_site !,s_move,s_str,
      use m_lmfinit, only: 
     &     slabl_,v_sbz,v_sctrl,v_sctrl2,v_sham,v_spot,v_slat,v_slat2,v_smix,v_sarry,v_sarry2,v_sspec,v_ssite,
     &     v_ssite2,sstrnmix,sstrnsymg
      use m_lmfinit,only: symgaf !read in call rdctrl2
      use m_hamindex,only: 
     &     iclasstaf,symops_af,ag_af,ngrpaf !we give them in this routine and used in sugw.F.
c     &  symops,invgx,ag,
c     &  shtvg,ngrp,miat,tiat,iclasst,
c     use m_globalvariables
      use m_ext,only: getext,sname
      use m_globalvariables,only: globalvariables
      use m_work,only: rv_a_opos
      implicit none
c      type(s_bz):: v_sbz
c      type(s_ctrl):: v_sctrl,v_sctrl2
c      type(s_lat):: v_slat,v_slat2
c      type(s_ham):: v_sham
c      type(s_pot):: v_spot
c      type(s_mix):: v_smix
c      type(s_array):: v_sarry,v_sarry2
c      type(s_move):: v_smove
c      type(s_str):: v_sstr
c      type(s_spec),allocatable:: v_sspec(:)
c      type(s_site),allocatable:: v_ssite(:),v_ssite2(:)
c      integer,parameter:: lstrn=10000
c      character(lstrn):: sstrnmix,sstrnsymg
      integer:: procid, master, mpipid, nproc,
     &     i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg, igets,a2vec,
     &     nfilin,nrecs,fopna,ic,ig,ib,ibas,nspec,inumaf,iprint,nbas,ioff
      logical,parameter:: T=.true., F=.false.
      integer,parameter::  NULLI=-99999,mxrecs=10000,recln0=511
      real(8),parameter::  NULLR =-99999
      real(8)::pwmode,pwemin,pwemax, bz_w,dval
      character*(recln0):: recrd(mxrecs)
      character(8) :: alabl
      character(5) :: charnum5
      character prgnam*8, ext*100 
      character:: strn*1000,outs*20,strn2*1000
      logical:: fileexist,ipr10,lshow,lshowp,swtmp,cmdopt,ltet
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
      if (procid .eq. master) call headl2(prgnam,0,stdo)
c      call finits(2,0,0,i)
      if (procid .eq. master) then
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
!   We will give up mathematical operation in ctrl.   
      nrecs = 0
      if(procid == master) then
         open(newunit=nfilin,file='ctrl.'//trim(sname)) ! = fopna('CTRL',-1,1)
         call findctrlstart(nfilin) ! ctrl is read from 'ctrlstart' if it exist.
         alabl = '#{}% ct '
         call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs) 
         close(nfilin)
      endif
#if MPI | MPIK
      call mpibc1(nrecs,1,2,mlog,'main','nrecs')
      call MPI_BCAST(recrd,recln0*(nrecs+1),MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
#endif
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
      call rdctrl2(recrd,recln0,nrecs,prgnam) 
!!     . slabl_,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite, !v_sstr,
!!     . v_sarry,sstrnmix,sstrnsymg  
#if MPI | MPIK
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif
      print *,'end of rdctrl2'
!!  -- End of readin ctrl file. ----------------------------
!!--- Lattice setup:deformation --- 
      call lattic() !this gives v_slat
      do ibas=1,v_sctrl%nbas
        v_ssite(ibas)%pos = rv_a_opos(3*(ibas-1)+1:3*(ibas-1)+3) !atom positions in the primitive cell
      enddo
!! --- We store all settings in module variables of m_rdctrl.
      print *,'end of lattic'
     
      
      
!!--- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
      if (cmdopt('--nosym',7,0,outs).or.cmdopt('--pdos',6,0,outs)) strn = ' '
!! when lmfgw, 1st digit of lc is zero--> no inversion added in mksym.
      lc = 10
#ifndef LMFGWD 
!! iand(2**n,i) means 2**n ditit of i.
      if( v_sctrl%noinv==0 ) lc = lc+2 ! inversion for lmf
#endif
      
c$$$!! Read AF case 2015, SPGGRAF category and SITE_ATOM_AF
c$$$      ipr10= iprint()>10 !this is only for master
c$$$      if(len_trim(symgaf)>0) then
c$$$        nbas = v_sctrl%nbas
c$$$        if(ipr10) then
c$$$          write(6,*)
c$$$          write(6,"(a)")       ' AF: ======================================== '
c$$$          write(6,"(a)")       ' AF: Antiferro mode: SPGGRAF='//trim(symgaf)
c$$$          write(6,"(a)")       ' AF:  (neglct waring in GENSYM) '
c$$$          write(6,"(a)")       ' AF:  Caution: I use v_sarry2=v_arry. Is this allowed in fortran?'
c$$$          do j=1,nbas
c$$$            write(6,"(a,2i3)") ' AF:  ibas,AF=',j,v_ssite(j)%iantiferro
c$$$          enddo
c$$$        endif
c$$$c        v_sarry2=v_sarry
c$$$        allocate(v_ssite2(nbas))
c$$$        v_ssite2=v_ssite
c$$$        inumaf=0
c$$$        do j=1,nbas
c$$$          do k=j,nbas
c$$$            if( v_ssite2(j)%iantiferro+v_ssite2(k)%iantiferro==0) then
c$$$              v_sarry2%iv_a_oips(k) =v_sarry2%iv_a_oips(j)
c$$$              inumaf=inumaf+1
c$$$              exit
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$        v_sctrl2 = v_sctrl
c$$$        v_sctrl2%nspec= v_sctrl%nspec - inumaf
c$$$        v_slat2=v_slat
c$$$        strn2=trim(strn)//' '//trim(symgaf)
c$$$        if(procid==master) call pshpr(60)
c$$$        call mksym(lc,slabl_,strn2 ,v_sctrl2,v_slat2,v_ssite2,v_sarry2)
c$$$        if(procid==master) call poppr()
c$$$        if(ipr10) write(6,"(a)")' AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
c$$$        ngrpaf=v_slat2%nsgrp
c$$$        allocate(iclasstaf(nbas),symops_af(3,3,ngrpaf),ag_af(3,ngrpaf))
c$$$        call dcopy ( ngrpaf * 9 , v_slat2%rv_a_osymgr , 1 , symops_af , 1 )
c$$$        call dcopy ( ngrpaf * 3 , v_slat2%rv_a_oag ,    1 , ag_af , 1 )
c$$$        if(ipr10) write(6,"(a,i3)") ' AF: ngrpaf=',ngrpaf
c$$$        do ib=1,nbas
c$$$          iclasstaf(ib)=v_ssite2(ib)%class
c$$$        enddo
c$$$!! We use  symops_af, ngrpaf and so on in suham.F to make a table in sugw.F-gen_hamindex.F.
c$$$!!   It is passed through m_hamindex
c$$$!!   The table is for rotation of the self-energy file sigm.*, used in rdsigm2.F-hamfb3k-rotsig.
c$$$!!
c$$$c$$$        call mptauof ( symops_af , ngrpaf , v_slat2%plat , nbas , v_slat2%rv_a_opos , iclasstaf
c$$$c$$$     .   , miat , tiat , invgx , shtvg )
c$$$c$$$        do ig = 1,ngrpaf
c$$$c$$$          write(*,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
c$$$c$$$          do i=1,3
c$$$c$$$            write(*,"(3e24.16)") symops_af(i,1:3,ig)
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$          if(ipr10) write(6,"(a)")       ' AF: ===== end of AF section================= '
c$$$          if(ipr10) write(6,"(a)")       
c$$$       endif
       
!! ---------------------------
      if(procid==master) call pshpr(60)
      call mksym(lc,slabl_,strn,v_sctrl,v_slat,v_ssite)!,v_sarry)
      if(procid==master) call poppr
c$$$!! test1 ------------------
c$$$      ngrp=v_slat%nsgrp
c$$$      nbas=v_sctrl%nbas
c$$$      allocate(iclasst(nbas),invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp),symops(3,3,ngrp),ag(3,ngrp))
c$$$      call dcopy ( ngrp * 9 , v_slat%rv_a_osymgr , 1 , symops , 1 )
c$$$      call dcopy ( ngrp * 3 , v_slat%rv_a_oag , 1 , ag , 1 )
c$$$      write(*,*) 'ngrp=',ngrp
c$$$      do ib=1,nbas
c$$$        iclasst(ib)=v_ssite(ib)%class
c$$$      enddo
c$$$      call mptauof ( symops , ngrp , v_slat%plat , nbas , v_slat%rv_a_opos , iclasst 
c$$$     . , miat , tiat , invgx , shtvg )
c$$$      do ig = 1,ngrp
c$$$        write(*,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
c$$$        do i=1,3
c$$$          write(*,"(3e24.16)") symops(i,1:3,ig)
c$$$        enddo
c$$$      enddo

!!  Set default values for species data
      call defspc(v_sspec)
#ifndef LMFGWD
      ltet = iand(1,int(v_sctrl%lmet)) .ne.0 .or. iand(4+2+1,int(v_sctrl%ldos)) .ne.0
!! iand(2**n,i) means 2**n ditit of i.
!! iand(4+1,i) = 1,0,0,4,5,4,5,0, 1,0,0,4,5,4,5,0, ... for i=1,2,3,4,5,6 (I think).
      call mkqp(v_sctrl,v_sbz,v_slat,ltet,F,1,-2)
#endif

      
!! show and quit
      if ( int(v_sctrl%quit) == 1 ) then
        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif
!! main part of lmf ------------------------------------
      call lmfp(prgnam,v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .  v_smix,v_sspec,v_ssite,trim(sstrnmix)) !,v_smove,v_sstr,,v_sarry
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
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
