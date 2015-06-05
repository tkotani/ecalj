      program lmf
      use m_struc_def,only: s_bz,s_ctrl,s_lat,s_ham,s_pot,s_mix,s_array,s_move,s_str,s_spec,s_site !,s_tb
      use m_rdctrl2_func, only: rdctrl2

      use m_rdctrl,only: symgaf !read in call rdctrl2
      use m_hamindex,only: 
     & iclasstaf,symops_af,ag_af,ngrpaf !we give them in this routine and used in sugw.F.
c     &  symops,invgx,ag,
c     &  shtvg,ngrp,miat,tiat,iclasst,

c      use m_globalvariables
#if MPI | MPIK
#if MPE
      use m_events !Event numbers: include "events.ins"
#endif
#endif
      implicit none
!! == Main program for lmf and so on. ==
!! For lmf, lmfgw. Single core, MPI, and MPIK modes.
!! ============================
      integer procid, master, mpipid, nproc
#if MPI | MPIK
      include "mpif.h"
#if MPE
      include "mpef.h"
#endif
#if MPI | MPIK
      integer fext
#endif
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog
c#if MPE
cC Event numbers:
c      include "fp/events.ins"
c#endif MPE
#endif

C ... Controls for IO
c      integer lstrn
c      parameter (lstrn=20000)

      character prgnam*8, vrsion(2)*6, ext*100
C ... Structure for strings
c      character*(lstrn) sstrn

C ... For structures
      logical lgors
      integer mxspec
c      type(s_sstrn):: v_sstrn
      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl,v_sctrl2
      type(s_lat):: v_slat,v_slat2
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_mix):: v_smix
      type(s_array):: v_sarry,v_sarry2
      type(s_move):: v_smove
      type(s_str):: v_sstr
c      type(s_tb):: v_stb
      type(s_spec),allocatable:: v_sspec(:)
      type(s_site),allocatable:: v_ssite(:),v_ssite2(:)
      parameter (mxspec=256)
      character*8 slabl(mxspec) !,slablaf(mxspec)
C ... miscellaneous local variables
      character strn*1000,outs*20,strn2*1000
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     .igets,a2vec,NULLI
      double precision vn(2),vn2(2),dglob
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)
C ... Program-dependent name and help
#if LMFGWD
      data prgnam /'LMFGWD'/
#elif LMF
      data prgnam /'LMF'/
#elif  LMDOS
      data prgnam /'LMDOS'/
#endif
      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=120)
      character*8 alabl
      character*(mxrecs*recln0):: recrd
      real(8),parameter::    NULLR =-99999
      real(8)::pwmode,pwemin,pwemax
      integer::ic
      real(8):: bz_w,dval!test
      logical:: fileexist,ipr10

      integer,parameter:: lstrn=10000
      character(lstrn):: sstrnmix,sstrnsymg !,jobid !header,

      integer:: ig,ib,ibas,nspec,inumaf,iprint,nbas

!! =================================================================
#if MPI|MPIK
      integer  nsize, id
      call mpi_init(ierr)
cstart--- check mpi. You can delete this.
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      if(id==0) write(*,*) 'mpi size=',nsize
cend-----
#endif

C --- Version ---
      call show_programinfo(6)
      vrsion(1) = 'LM'
      vrsion(2) = ' '
      vn(1) = 7.00d0
      vn2(1) = 0
      vn(2) = 0.00d0
      vn2(2) = 0
#if LMF | LMFGWD
      vrsion(2) = 'FP'
      vn(2) = 7.00d0
      vn2(2) = 0
#endif
      vn(1) = vn(1) + vn2(1)/1d5
      vn(2) = vn(2) + vn2(2)/1d5
      stdo = lgunit(1)
      i = dglob('stdo',dble(stdo),1)
      master = 0
      procid = mpipid(1)
      nproc  = mpipid(0)

C --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs))
     .call lmhelp(prgnam,vn(2),vn2(2),0) !no common/w
      if (cmdopt('--version',9,0,outs)) then
        write(*,'(f4.2)') vn(2)
        call rx0('')
      endif

C --- Dynamic memory allocation and other initialization ---
#if MPI | MPIK
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      mlog = cmdopt('--mlog',6,0,strn)
#if MPE
      ierr = MPE_INIT_LOG()
      EVENT_START_RDCTRL = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_RDCTRL   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_UGCOMP = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_UGCOMP   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_FSMBL  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_FSMBL    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_PZHEV  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_PZHEV    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_KLOOP  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_KLOOP    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_MIXRHO = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_MIXRHO   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_SMHSBL = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_SMHSBL   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_AUGMBL = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_AUGMBL   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_HSIBL  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_HSIBL    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_RSIBL  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_RSIBL    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_RLOCBL = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_RLOCBL   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_DFRCE  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_DFRCE    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_BCAST  = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_BCAST    = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_ALLRED = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_ALLRED   = MPE_LOG_GET_EVENT_NUMBER()
      EVENT_START_BARRIER= MPE_LOG_GET_EVENT_NUMBER()
      EVENT_END_BARRIER  = MPE_LOG_GET_EVENT_NUMBER()
#endif
C --- Dynamic memory allocation and other initialization ---
      if (procid .eq. master) call headl2(prgnam,0,stdo)
      call finits(2,0,0,i)
      if (procid .eq. master) then
        i = fextg(ext)
      endif
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      if (procid .eq. master) then
#if MPE
        ierr = MPE_DESCRIBE_STATE(EVENT_START_RDCTRL,EVENT_END_RDCTRL,"rdctrl","pink")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_UGCOMP,EVENT_END_UGCOMP,"ugcomp","maroon")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_FSMBL,EVENT_END_FSMBL,  "fsmbl","aquamarine")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_PZHEV,EVENT_END_PZHEV,  "pzhev","brown")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_KLOOP,EVENT_END_KLOOP,  "k-loop","brown")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_MIXRHO,EVENT_END_MIXRHO,"mixrho","orange")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_SMHSBL,EVENT_END_SMHSBL,"smhsbl","blue")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_AUGMBL,EVENT_END_AUGMBL,"augmbl","cyan")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_HSIBL,EVENT_END_HSIBL,  "hsibl","gray")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_RSIBL,EVENT_END_RSIBL,  "rsibl","red")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_RLOCBL,EVENT_END_RLOCBL, "rlocbl","green")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_DFRCE,EVENT_END_DFRCE,  "dfrce","magenta")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_BCAST,EVENT_END_BCAST,  "broadcast","coral")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_ALLRED,EVENT_END_ALLRED, "allreduce","purple")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_BARRIER,EVENT_END_BARRIER, "barrier","yellow")
#endif
        call gettime(datim)
        if (mlog) i = fopn('MLOG')
        call poseof(fopn('LOG'))
        if (mlog) then
          call awrit2(' lmf '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//' is master',' ',
     .    256,lgunit(3),procid,nproc)
        endif
      else !for procid !=0
        call strcat(ext,20,' ','_',1,' ',i) !mpi file extention
        call bin2a(' ',0,0,procid,2,0,20,ext,i)
        ierr = fext(ext(1:i+1))
c        print *,' xxxx file extension: ',procid,ext(1:i+1)
        if (mlog) ierr = fopn('MLOG')
        ierr = fextg(ext)
        call gettime(datim)
        if (mlog) then
          call awrit2(' lmf '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' file extension is '//ext(2:i+1),' ',
     .    256,lgunit(3),procid,nproc)
        endif
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#else
C     Add version to variables table?
      call headl2(prgnam,0,stdo)
      call finits(2,0,0,i)
      call poseof(fopn('LOG'))
      i = fextg(ext)
      call word(ext,1,i,j)
      if (ext(i:i) .eq. '.') i=i+1
      if (ext(j:j) .eq. '.') j=j-1
      if (j .ge. i) call ptenv('EXT='//ext(i:j))
#endif

C ... Abort with error message if ctrl file is missing (swtmp = .true.)
      swtmp = .false.
      if (procid==master) then
      if (.not.cmdopt('--input',6,0,strn)) then
        if (fxst('CTRL') .ne. 1) swtmp=.true.
      endif
      endif
      call mpibc1(swtmp,1,1,.false.,'lmf','error')
      if(swtmp) call rx("No ctrl file found!") 

C ... Set special file directory for temporary files
C     User may which to customize the directory
C     Default is to use the standard directory
C     sttmpd is located at the bottom of this file.
      call sttmpd

C ... File logical units
      i = fadd('TMP',-1,4)
      i = fadd('BAND',-1,4)

!! Set the top-level verbosity if specified from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif
!!   iprint() is negative except master
      if(procid/=master) call setpr(-100) 

C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
      if (procid .eq. master) then
        if (.not.cmdopt('--input',7,0,strn)) then
          nfilin = fopna('CTRL',-1,1)
          call findctrlstart(nfilin) !takao
          alabl = '#{}% ct '
          if (cmdopt('--show',6,0,strn)) alabl = '#{}% ctp'
          call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
          i = 60
          if (cmdopt('--show',6,0,strn)) then
            i = 1
            call info0(i,0,0,' ------------------------ End of input file ----------------------')
          endif
          call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     .    ' input file',nrecs,0)
          if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
        endif
      endif
#if MPI | MPIK
      call mpibc1(nrecs,1,2,mlog,'main','nrecs')
      call MPI_BCAST(recrd,recln0*(nrecs+1),MPI_CHARACTER,
     .master,MPI_COMM_WORLD,ierr)
#endif
c      if( cmdopt('--pdos',6,0,strn)) then
c        write(*,*) ' --pdos: Enforce SYMGRP e. Read sigm_fbz for RDSIG. Set --quit=band.'
c      endif

C      write(stdo,*)' === TEST out recrd === '
C      j = 0
C      do i = 1, nrecs
C        write(stdo,"(i2,1x,256a1)") procid,recrd(j:j+recln0-1)
C        j = j + recln0
C      enddo
C      write(stdo,*)' --- end of rdfile --- '
C      call rx('done')

!! --- Read recrd ---
      call rdctrl2(recrd,recln0,nrecs,prgnam,vrsion,vn,vn2,F,
     . slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     . v_sarry,v_smove, sstrnmix,sstrnsymg)  !,v_stb

#if MPI | MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_RDCTRL,procid,"rdctrl")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
#endif
#endif
C --- Lattice setup ---
c      print *,'goto setcg'
      call setcg(v_slat,8,12)
      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
C --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
!! mar2013 --pdos Also in rdctrl2.F
      if (cmdopt('--nosym',7,0,outs).or.cmdopt('--pdos',6,0,outs)) strn = ' '
      lc = 10

#ifndef LMFGWD
      if ( .not. iand(1,int(v_sctrl%lqp)) .ne.0 ) lc = lc+2
#endif
      if (cmdopt('--sfill',7,0,strn)) then
        call rx('no --sfill option in lm7K')
      endif

!! read AF case 
      ipr10= iprint()>10 !this is only for master
      if(len_trim(symgaf)>0) then
        nbas = v_sctrl%nbas
        if(ipr10) then
          write(6,*)
          write(6,"(a)")       ' AF: ======================================== '
          write(6,"(a)")       ' AF: Antiferro mode: SPGGRAF='//trim(symgaf)
          write(6,"(a)")       ' AF:  (neglct waring in GENSYM) '
          write(6,"(a)")       ' AF:  Caution: I use v_sarry2=v_arry. Is this allowed in fortran?'
          do j=1,nbas
            write(6,"(a,2i3)") ' AF:  ibas,AF=',j,v_ssite(j)%iantiferro
          enddo
        endif
        v_sarry2=v_sarry
        allocate(v_ssite2(nbas))
        v_ssite2=v_ssite
        inumaf=0
        do j=1,nbas
          do k=j,nbas
            if( v_ssite2(j)%iantiferro+v_ssite2(k)%iantiferro==0) then
              v_sarry2%iv_a_oips(k) =v_sarry2%iv_a_oips(j)
              inumaf=inumaf+1
              exit
            endif
          enddo
        enddo
        v_sctrl2=v_sctrl
        v_sctrl2%nspec= v_sctrl%nspec-inumaf
        v_slat2=v_slat
        strn2=trim(strn)//' '//trim(symgaf)
        if(procid==master) call pshpr(60)
        call mksym(lc,slabl,strn2 ,v_sctrl2,v_slat2,v_ssite2,v_sarry2)
        if(procid==master) call poppr()
        if(ipr10) write(6,"(a)")' AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
        ngrpaf=v_slat2%nsgrp
        allocate(iclasstaf(nbas),symops_af(3,3,ngrpaf),ag_af(3,ngrpaf))
c    ,invgx(ngrpaf),miat(nbas,ngrpaf),tiat(3,nbas,ngrpaf), shtvg(3,ngrpaf),
        call dcopy ( ngrpaf * 9 , v_slat2%rv_a_osymgr , 1 , symops_af , 1 )
        call dcopy ( ngrpaf * 3 , v_slat2%rv_a_oag ,    1 , ag_af , 1 )
        if(ipr10) write(6,"(a,i3)") ' AF: ngrpaf=',ngrpaf
        do ib=1,nbas
          iclasstaf(ib)=v_ssite2(ib)%class
        enddo
!! We use  symops_af, ngrpaf and so on in suham.F to make a table in sugw.F-gen_hamindex.F.
!!   It is passed through m_hamindex
!!   The table is for rotation of the self-energy file sigm.*, used in rdsigm2.F-hamfb3k-rotsig.
!!
c$$$        call mptauof ( symops_af , ngrpaf , v_slat2%plat , nbas , v_slat2%rv_a_opos , iclasstaf
c$$$     .   , miat , tiat , invgx , shtvg )
c$$$        do ig = 1,ngrpaf
c$$$          write(*,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
c$$$          do i=1,3
c$$$            write(*,"(3e24.16)") symops_af(i,1:3,ig)
c$$$          enddo
c$$$        enddo
          if(ipr10) write(6,"(a)")       ' AF: ===== end of AF section================= '
          if(ipr10) write(6,"(a)")       
      endif


!! ---------------------------
      if(procid==master) call pshpr(60)
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
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


C ... Set default values for species data
      call defspc(v_sspec)

C ... Patch for now ... maybe replace
#if LMF
#ifndef LMFGWD
      ltet = iand(1,int(v_sctrl%lmet)) .ne.0 .or. iand(4+2+1,int(v_sctrl%ldos)) .ne.0
!! iand(4+1,i) = 1,0,0,4,5,5,4,0,  1,0,0,4,5,5,4,...
      call mkqp(v_sctrl,v_sbz,v_slat,ltet,F,1,-2)
#endif
#endif

C ... quit after SHOW
      if ( int(v_sctrl%quit) == 1 ) then
        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif

#if LMF
      call lmfp(prgnam,v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .v_smix,v_sspec,v_ssite,v_sstr,v_smove,v_sarry,trim(sstrnmix))
      inquire(file='RSEQ_ERROR',exist=fileexist)
      if(fileexist) call rx('RSEQ_ERROR: maybe in phidx.F See takao comment in it.')
#endif
!! -------------- End of program -------------
 1000 continue
#if MPI | MPIK
#if MPE
      if (procid .eq. master) i = fextg(ext)
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      i = 0
      call skp2bl(ext,20,i)
      ierr = MPE_FINISH_LOG(ext(2:i))
#endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if ( procid .eq. master ) then
        call rx0(prgnam//' on '//shortname(procid)(1:namelen(procid)))
      else
        call fexit(0,0,' ',0)
      endif
#else
      call rx0(prgnam)
#endif
#if MPI|MPIK
      call mpi_finalize(ierr)
#endif
      end program lmf


      subroutine lmhelp(prgnam,vn,vn2,wksize)
C- Help printout
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   vn    :program main version
Ci   vn2   :program minor version
Ci   wksize:worksize
Co Outputs
Co   message written to stdout
Cr Remarks
Cu Updates
Cu   11 Apr 03
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      character prgnam*8
      double precision vn,vn2
      integer wksize
C ... Local parameters
      integer i1,i2
      character ch*1,outs*1000

      call locase(prgnam)
      call info0(0,0,0,' usage:  '//prgnam//
     .'%a [--OPTION] [-var-assign] [ext]')

      print 343
      print 344
  343 format(/' --h'/' --help',t17,'Print this message, and quit'
     ./' --input',t17,
     .'List categories, tokens, and data program expects, and quit'
     ./' --show',t17,
     .'Print control file after parsing by preprocessor,'/t17,
     .'and echo input data as read from the control file'
     ./' --showp',t17,
     .'Same as --show, but quit after input parsed'
     ./' --iactiv',t17,'(--no-iactiv) ',
     .'Turn on (off) interactive mode'/
     .t17,'This switch overrides input file setting',
     ./' --pr#1[,#2...]',t17,
     .'Set the verbosity (stack) to values #1,#2, ...'
     ./' --time=#1[,#2]',t17,
     .'Print timing info to # levels (#1=summary; #2=on-the-fly)'/
     ./' -vnam=expr',t17,
     .'Define numerical variable "nam"; set to result of ''expr'''
     ./' -cnam=strn',t17,
     .'Define character variable "nam"; set to ''strn'''
     .)


  344 format(
     ./' --rpos=filnam',t17,
     .'After reading input file, read site positions from "filnam"'/
     .' --fixpos[:tol=#]',
     .' Adjust positions slightly, rendering them as'/
     .t17,'exactly as possible consistent with the symmetry group')

      if (.true.) then
        outs = '%N '//prgnam//'%a-specific options:'
        call strip(outs,i1,i2)
        call info0(0,0,0,outs(1:i2))

        if (prgnam .eq. 'lmfgwd') then
          call info0(0,0,0,
     .    '%N%1f ')
        endif

        if (prgnam .eq. 'lmf') then
          call info0(0,0,0,
     .    '%N%1f--rs=#1,#2,#3,#4,#5'//
     .    '%N%6f#1=0 start from atm file; 1 from rst file;'//
     .    ' 2 from rsta file'//
     .    '%N%11fadd 10 to shift sm-rho 1st iter'//
     .    '%N%11fadd 100 to rotate local rho 1st iter'//
     .    '%N%6f#2=1 save rst file'//
     .    '%N%6f(#3,#4,#5)=0 read (pos,E_f,pnu) from rst file'//
     .    '%N%N --band --wden --rdbasp --optbas --etot'//
     .    ' (see documentation)')
        endif

      endif

      if (nint(vn2) .ne. 0) then
        ch = char(nint(vn2)+ichar('a')-1)
      else
        ch = ' '
      endif
      call info2(0,0,0,'%N version %,2d'//ch//'   worksize %i K',
     .vn-vn2/1d5,wksize/1000)
      call fexit(0,0,' ',0)
      end subroutine lmhelp


      subroutine sttmpd
C- Creates special tmpdir for saving temporary files
C  User may wish to customize this routine.
      character tmpdir*100
      integer fopnT,ifi
C ... for henry, lm-MPIK
C      integer i1,i2,nw
C      character*40 strn

C     return
C ... Set customization of temporary directory here, if desired
C     This is usual default (current working directory)
      tmpdir = ' '
C     call gtenv('HOME',tmpdir)
C     call gtenv('TMPDIR',tmpdir)
C     call getenv('HOME',tmpdir)
C     call getenv('TMPDIR',tmpdir)
C     tmpdir = '/home/tmp/'

C ... for henry, lm-MPIK specific
C      call getenv('HOME',strn)
C      call strip(strn,i1,i2)
C      call wrdsg(strn(i1:i2),0,'/',nw)
C      call wordg(strn,0,'/',nw,i1,i2)
C      strn = strn(i1:)
C      if (strn .eq. 'markv') strn = 'svan2'
C      call word(strn,1,i1,i2)
C      tmpdir = '/home/' // strn(i1:i2) // '/tmp'

C ... Set the directory
      ifi = fopnT(tmpdir,0,0,11)

C     debugging check
C      ifi = fopnT('tmp' ,-1,0,0)
C      call fshow
C      print *, ifi
C      write(ifi,*) 'test'
C      call rx('done')
Cgetarg       end
      end subroutine sttmpd



      subroutine findctrlstart(nfilin)
! if we find 'ctrlstart', locate reading at the next line of ctrlstart.
! this is useful if you like to use script, GWinput, ctrl in a file.
      character(len=9):: strn, ccc
      do
        read(nfilin,"(a)",err=1010,end=1010) strn
        if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
      end subroutine findctrlstart

