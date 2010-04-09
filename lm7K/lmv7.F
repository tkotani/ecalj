Cgetarg...info...           structure ['fmain', 'lat', 'v_slat']
Cgetarg...info...           structure ['fmain', 'ctrl', 'v_sctrl']
Cgetarg...info...           structure ['fmain', 'ctrl', 'v_sbz']
      program main
      implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer ierr, nsize, id 
      character(14):: testbuf
      call mpi_init(ierr)
cstart--- check mpi. You can delete this. 
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      write(*,*) 'mpi: size and id=',nsize,id
cend-----
#endif      
      call fmain()
#if MPI|MPIK
      call mpi_finalize(ierr)
#endif
      end program main

C -------------- Main program for LM and related programs ----------
C This file contains the main program for most of the collection
C of main programs in the ASA package and its extensions.
C
C The various main programs are created from this one using the ccomp
C utility.
C
C It is documented in doc/lmto.html.
C Updates
C   13 May 04 New main lmmc.f
C ----------------------------------------------------
C ... Operating system permitting, The true main is in slatsm.a (fmain.c)
#define unix 1
#if unix
      subroutine fmain

#endif
      use m_struc_def
      use m_rdctrl2_func, only: rdctrl2
      implicit none

      integer procid, master, mpipid, nproc
#if MPI | MPIK
      include "mpif.h"
#if MPE
      include "mpef.h"
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
#if MPE
C Event numbers:
      include "fp/events.ins"
#endif MPE
#endif

C Heap allocation
      integer wksize
#if MPI | MPIK
      parameter(wksize= 300 000 000)
#else
      parameter(wksize= 300 000 000)
#endif
      integer w(wksize)
C     Next two lines guarantee w is aligned along a d.p. boundary
      double precision ws
      equivalence (ws,w(1))
      common /w/ w

C ... Controls for IO
      integer lstrn
      parameter (lstrn=2000)

      character prgnam*8, vrsion(2)*6, ext*100
C ... Structure for strings
      character*(lstrn) sstrn

C ... For structures
      logical lgors
      integer mxspec
CKi      integer v_sbz,v_sctrl,v_slat,v_sham,v_spot,v_smix,v_sspec,v_ssite,v_sarry,
CKi     .  v_smove,v_sstr,v_stb
      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl
      type(s_lat):: v_slat
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_mix):: v_smix
      type(s_array):: v_sarry
      type(s_move):: v_smove
      type(s_str):: v_sstr
      type(s_tb):: v_stb

      type(s_spec),pointer:: v_sspec(:)
      type(s_site),pointer:: v_ssite(:)

      parameter (mxspec=256)
      character*8 slabl(mxspec)

C ... miscellaneous local variables
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     .  igets,a2vec,NULLI
      double precision vn(2),vn2(2),dglob
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)
#if LMCTL
      integer fopnn
#endif
#if MPI | MPIK
      integer fext
#endif
#if TBE
      integer ltbe
      double precision mdprm(6)
#endif

C ... Program-dependent name and help
#if LMPG
#elif LM
#elif LMFA
      data prgnam /'LMFA'/
#elif LMFGWD
      data prgnam /'LMFGWD'/
#elif LMF
      data prgnam /'LMF'/
#elif LMAQU
      data prgnam /'LMAQU'/
#elif  LMDOS
      data prgnam /'LMDOS'/
#elif  LMCTL
      data prgnam /'LMCTL'/
#elif  LMIMP
      data prgnam /'LMIMP'/
#elif  LMCOR
      data prgnam /'LMCOR'/
#elif  LMFIT
      data prgnam /'LMFIT'/
#elif  LMPLAN
      integer nbas,nbasp,nbaspp,nl,opnu,oqnu,oqc,oqt,oves
      double precision emad,trumad,vmtz(2)
      data prgnam /'LMPLAN'/
#elif  LMAVGM
      data prgnam /'LMAVGM'/
#elif  LMMIX
      data prgnam /'LMMIX'/
#elif  LMSHF
      data prgnam /'LMSHF'/
      data auxmod /8/
#elif  TBPG
      data prgnam /'TBPG'/
#elif  TBE
      data prgnam /'TBE'/
#elif  TBBND
      data prgnam /'TBBND'/
#elif  TBFIT
      data prgnam /'TBFIT'/
#elif  TBDOS
      data prgnam /'TBDOS'/
#elif  MMAG
      data prgnam /'MMAG'/
#endif

C ... Program-dependent cagetories
#if LMPG
      data prgnam /'LMPG'/
#elif LMGPOL
      data prgnam /'LMGPOL'/
      data auxmod /32/
#elif BLM
      data prgnam /'BLM'/
#elif LM
      data prgnam /'LM'/
#elif LMMC
      data prgnam /'LMMC'/
#elif LMAQU | LMMIX
      call rx('missing setup for lmaqu, lmmix')  ! Why???TK
#elif LMGF
      data prgnam /'LMGF'/
#elif LMSTR
      data prgnam /'LMSTR'/
#elif LMCHK
      data prgnam /'LMCHK'/ auxmod /1/
#elif LMSCELL
      logical rdstrn
      integer nkd,odlv,nbx,ix(9)
      double precision plat(3,3),plx(3,3),ddot,dmach
      data prgnam /'LMSCELL'/
#elif LMXBS
      data prgnam /'LMXBS'/
      data auxmod /4/
#endif

      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=120)
      character*8 alabl
c      character*1 recrd(0:1)
c      equivalence( recrd(0),w(1))
CKi      character*1,allocatable:: recrd(:)
      character:: recrd_(0:mxrecs*recln0-1)
      character*(1000):: recrd
      equivalence ( recrd, recrd_ )


      real(8),parameter::    NULLR =-99999
      real(8)::pwmode,pwemin,pwemax


      real(8):: bz_w !test

C -------------- First executable statement ---------------
C     Link in routines that cause potential library conflicts
      call nada

c toksw was here. MOved to rdccat

C --- Version ---
      vrsion(1) = 'LM'
      vrsion(2) = ' '
      vn(1) = 7.00d0
      vn2(1) = 0
      vn(2) = 0.00d0
      vn2(2) = 0
#if LM | LMPG | LMGF | LMAVGM | LMIMP | LMCTL | LMSTR
      vrsion(2) = 'ASA'
      vn(2) = 7.00d0
      vn2(2) = 0
#elif LMFA | LMF | LMFGWD
      vrsion(2) = 'FP'
      vn(2) = 7.00d0
      vn2(2) = 0
#elif LMMC
      vrsion(2) = 'MOL'
      vn(2) = 3.00d0
      vn2(2) = 0
#elif TBE
      vrsion(2) = 'TB'
      vn(2) = 9.00d0
      vn2(2) = 0
#elif MMAG
      vrsion(2) = 'MM'
      vn(2) = 2.00d0
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
     .  call lmhelp(prgnam,vn(2),vn2(2),wksize)
      if (cmdopt('--version',9,0,outs)) then
        write(*,'(f4.2)') vn(2)
        call cexit(0,1)
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
      if (procid .eq. master) call headl2(prgnam,wksize,stdo)
      call finits(2,0,0,i)
      call pshpr(0)
      call wkinit(wksize)
      call wkfast(T)
      call poppr
      if (procid .eq. master) then
        i = fextg(ext)
      endif
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      if (procid .eq. master) then
#if MPE
        ierr = MPE_DESCRIBE_STATE(EVENT_START_RDCTRL,EVENT_END_RDCTRL,
     .                            "rdctrl","pink")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_UGCOMP,EVENT_END_UGCOMP,
     .                            "ugcomp","maroon")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_FSMBL,EVENT_END_FSMBL,
     .                            "fsmbl","aquamarine")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_PZHEV,EVENT_END_PZHEV,
     .                            "pzhev","brown")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_KLOOP,EVENT_END_KLOOP,
     .                            "k-loop","brown")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_MIXRHO,EVENT_END_MIXRHO,
     .                            "mixrho","orange")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_SMHSBL,EVENT_END_SMHSBL,
     .                            "smhsbl","blue")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_AUGMBL,EVENT_END_AUGMBL,
     .                            "augmbl","cyan")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_HSIBL,EVENT_END_HSIBL,
     .                            "hsibl","gray")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_RSIBL,EVENT_END_RSIBL,
     .                            "rsibl","red")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_RLOCBL,EVENT_END_RLOCBL,
     .                            "rlocbl","green")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_DFRCE,EVENT_END_DFRCE,
     .                            "dfrce","magenta")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_BCAST,EVENT_END_BCAST,
     .                            "broadcast","coral")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_ALLRED,EVENT_END_ALLRED,
     .                            "allreduce","purple")
        ierr = MPE_DESCRIBE_STATE(EVENT_START_BARRIER,EVENT_END_BARRIER,
     .                            "barrier","yellow")
#endif
        call gettime(datim)
        if (mlog) i = fopn('MLOG')
        call poseof(fopn('LOG'))
        if (mlog) then
          call awrit2(' lmf '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//' is master',' ',
     .      256,lgunit(3),procid,nproc)
        endif
      else
        call strcat(ext,20,' ','_',1,' ',i)
        call bin2a(' ',0,0,procid,2,0,20,ext,i)
        ierr = fext(ext(1:i+1))
        if (mlog) ierr = fopn('MLOG')
        ierr = fextg(ext)
        call gettime(datim)
        if (mlog) then
          call awrit2(' lmf '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' file extension is '//ext(2:i+1),' ',
     .      256,lgunit(3),procid,nproc)
        endif
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#else
C     Add version to variables table?
C     j = 0
C     swtmp = a2bin('T',ltet,1,0,' ',j,-1)
C     call addsyv('version',vn,i)
      call headl2(prgnam,wksize,stdo)
      call finits(2,0,0,i)
      call pshpr(0)
      call wkinit(wksize)
      call wkfast(T)
      call poppr
      call poseof(fopn('LOG'))
      i = fextg(ext)
      call word(ext,1,i,j)
      if (ext(i:i) .eq. '.') i=i+1
      if (ext(j:j) .eq. '.') j=j-1
      if (j .ge. i) call ptenv('EXT='//ext(i:j))
#endif

#if BLM
      call suctrl(prgnam)
      call fexit(0,0,' ',0)
#else

C ... Abort with error message if ctrl file is missing (swtmp = .true.)
      swtmp = .false.
#if MPI | MPIK
      if (procid .eq. master) then
#endif
      if (cmdopt('--input',6,0,strn)) then
        if (nproc .gt. 0) call rx('--input not allowed with MPI')
      else
        if (fxst('CTRL') .ne. 1) then
          call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     .      i1mach(2))
          swtmp = .true.
        endif
      endif
#if MPI | MPIK
      endif
#endif
      call mpibc1(swtmp,1,1,.false.,'lmf','error')
      if (swtmp) call cexit(-1,1)

C ... Set special file directory for temporary files
C     User may which to customize the directory
C     Default is to use the standard directory
C     sttmpd is located at the bottom of this file.
      call sttmpd

C ... File logical units
      i = fadd('TMP',-1,4)
      i = fadd('BAND',-1,4)
#if TBE
      i = fadd('STRT',-1,4)
      i = fadd('QMOM',-1,4)
#endif
#if LM | LMDOS | LMGF
      i = fadd('MOMS',-1,4)
      i = fadd('CLS',-1,4)
#endif

C --- Set the top-level verbosity if spec'd from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif

C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
Cki      allocate( recrd( 0:mxrecs*recln0-1 ) )
      if (procid .eq. master) then
      if (.not.cmdopt('--input',7,0,strn)) then
        nfilin = fopna('CTRL',-1,1)
ctakao
        call findctrlstart(nfilin)

        alabl = '#{}% ct '
        if (cmdopt('--show',6,0,strn)) alabl = '#{}% ctp'
        call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
c To show variables.
c        call shosyv(0,-1,0,6)

        i = 60
        if (cmdopt('--show',6,0,strn)) then
          i = 1
          call info0(i,0,0,' ------------------------ End '//
     .                     'of input file ----------------------')
        endif
        call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     .    ' input file',nrecs,0)
        if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
      endif
      endif
#if MPI | MPIK
      call mpibc1(nrecs,1,2,mlog,'main','nrecs')
      call MPI_BCAST(recrd,recln0*(nrecs+1),MPI_CHARACTER,
     .  master,MPI_COMM_WORLD,ierr)
#endif

C      write(stdo,*)' === TEST out recrd === '
C      j = 0
C      do i = 1, nrecs
C        write(stdo,"(i2,1x,256a1)") procid,recrd(j:j+recln0-1)
C        j = j + recln0
C      enddo
C      write(stdo,*)' --- end of rdfile --- '
C      call rx('done')

C --- Read recrd ---
      call rdctrl2(recrd,recln0,nrecs,prgnam,vrsion,vn,vn2,F,
     .  slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     .  v_sarry,v_smove,v_stb,sstrn)

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

#if LMSCELL
      call setcg(v_slat,8,12)
Cgetarg...info...           real(8) :: dmach , real(8) :: v_slat%tol 1
Cgetarg       call pack5('lat as tol nkdmx nkqmx',v_slat,1d0,dmach(2),150000,
Cgetarg      .  150000,0)
       
       v_slat%as=1d0 
       v_slat%tol=dmach(2) 
       v_slat%nkdmx=150000 
       v_slat%nkqmx=150000 

      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
Cgetarg...info...           integer :: nkd , integer(8) :: v_slat%nkd 1
Cgetarg...info...           integer :: odlv , integer(8) :: v_slat%odlv 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: v_slat%plat 9
Cgetarg...info...           real(8) :: plx(3,3) , real(8) :: v_slat%slat 9
Cgetarg       call upack('lat nkd odlv plat slat',v_slat,nkd,odlv,plat,plx,0)
       
       nkd=v_slat%nkd
       odlv=v_slat%odlv
       i_copy_size=size(v_slat%plat) 
       call dcopy(i_copy_size,v_slat%plat,1,plat,1) 
       i_copy_size=size(v_slat%slat) 
       call dcopy(i_copy_size,v_slat%slat,1,plx,1) 


C ... Set default values for species data
      call defspc(v_sspec)

      if (plx(1,1) .eq. NULLR .or. ddot(9,plx,1,plx,1) .eq. 0) then
   88   continue
        write(stdo,'(/a)') ' Enter lattice vectors of supercell:'
        strn = ' '
        i = i1mach(1)
        swtmp = rdstrn(i,strn,len(strn),F)
        call word(strn,1,i,j)
        if (i .gt. j .or. strn(i:i) .eq. '/') then
          call dpzero(plx,9)
        else
          i = 0
          if (a2vec(strn,len(strn),i,4,', ',2,-3,9,ix,plx) .ne. 9) then
            print *, 'failed to read 9 values from string ... redo'
            goto 88
          endif
        endif
        if (ddot(9,plx,1,plx,1) .eq. 0) call dcopy(9,plat,1,plx,1)
      endif

      call supcel(0,v_sctrl,v_slat,v_ssite,slabl,v_sspec,plx,
     .  nkd,w(odlv),nbx)
      goto 1000
#endif

C ... Input to scale MT radii (define RESIZE)
#if LM | LMCHK | LMPOL | LMPLAN | LMXBS | LMCOR | LMSTR | LMDOS
#define RESIZE 1
#endif
#if LMGF | LMPG | LMF | LMFGWD | LMFA | LMMC
#define RESIZE 1
#endif

C --- Lattice setup ---
#if ! (LMCTL | LMMC)
ctakao will add |LMFA)
      call setcg(v_slat,8,12)
      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
#endif

C --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      call upacks('strn symg',i,j)
      if (j .ge. i) strn = sstrn(i:j+1)
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      lc = 20
#if TBE
Cgetarg...info...           real(8) :: mdprm(6) , real(8) :: v_sctrl%mdprm 6
Cgetarg       call upack('ctrl mdprm',v_sctrl,mdprm,0,0,0,0)
       
       i_copy_size=size(v_sctrl%mdprm) 
       call dcopy(i_copy_size,v_sctrl%mdprm,1,mdprm,1) 

      if (nint(mdprm(1)) .eq. 1 .or. nint(mdprm(1)) .eq. 2) then
        strn = ' '
      endif
      lc = 40
#endif
#if MMAG
      lc = 0
#endif
#if LMFA
C      lc = 0
#endif
#if LMF
      lc = 10
#endif
#ifndef LMFGWD
Cgetarg       if (.not. lgors('ctrl lqp,1',v_sctrl)) lc = lc+2
       if ( .not. iand(1,int(v_sctrl%lqp)) .ne.0 ) lc = lc+2 

#endif
#if ! LMMC
ctakao will add |LMFA)
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
#else
      lc = 0
#endif

C --- Allocate permanent class arrays, maps and other initialization ---
      if (lc .ge. 20)
     .  call clsprm(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .  v_slat,v_sbz,v_sstr,sstrn)

C --- Read available class parameters from file ---
#if LMCTL | LMCHK
      call aiocls(.false.,0,v_sarry,v_sctrl,v_sham,v_spot,
     .  v_sspec,v_slat,1,0)
#elif LM | TBE | LMGF | LMSHF | LMPLAN
#ifndef TBE
      call aiocls(.false.,12,v_sarry,v_sctrl,v_sham,v_spot,
     .  v_sspec,v_slat,1,0)
#endif
c$$$      call rdctrl2(recrd,recln0,nrecs,prgnam,vrsion,vn,vn2,T,
c$$$     .  slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
c$$$     .  v_sarry,v_smove,v_stb,sstrn)
#if MPI | MPIK
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif

#ifndef TBE
      call aiocls(.false.,17,v_sarry,v_sctrl,v_sham,v_spot,
     .  v_sspec,v_slat,1,0)
      if (lc .ge. 20)
     .  call clsprp(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .  v_slat,v_sbz,v_sstr,sstrn)
#endif
#endif
C     Discard contents of input file
Cki      deallocate(recrd)

C --- Optionally resize spheres ---
#if RESIZE
      if (cmdopt('--sfill',7,0,strn)) then
        print *, '*** --sfill DEPRECATED .. see lmto.html'
        call sfill(strn(8:),slabl,v_sctrl,v_sarry,v_slat,
     .    v_sspec,v_ssite)
      else
        call sfill('auto',slabl,v_sctrl,v_sarry,v_slat,v_sspec,
     .    v_ssite)
      endif
#endif

C ... Set default values for species data
      call defspc(v_sspec)

C ... Patch for now ... maybe replace
#if LM | TBE | LMGF | LMF
#ifndef LMFGWD
Cgetarg       ltet = lgors('ctrl lmet,1',v_sctrl) .or.
Cgetarg      .       lgors('ctrl ldos,4+2+1',v_sbz)
       ltet = iand(1,int(v_sctrl%lmet)) .ne.0 .or. iand(4+2+1,int(v_sctrl%ldos)) .ne.0 

      call mkqp(v_sctrl,v_sbz,v_slat,ltet,F,1,-2)
#endif
#endif

C ... quit after SHOW
Cgetarg       if (igets('ctrl quit',v_sctrl) .eq. 1) then
       if ( int(v_sctrl%quit) .eq. 1 ) then 

        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif

#if LMFA
      call freeat(v_sctrl,v_sspec,v_sham,v_spot)
#endif

#if LMCHK & FP
      if (cmdopt('--fp',4,0,strn)) then
        call fpchk(v_sspec,v_ssite)
        call cexit(0,1)
      endif
#endif
#if LMCHK | LMXBS | LMGPOL | LMSHF
      if (cmdopt('--findes',8,0,strn)) auxmod = 128
      call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
     .  v_ssite,v_sstr,v_sarry,v_sbz,slabl,sstrn,auxmod,wksize)
#endif

#if LM
      call lmasa(prgnam,v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .  v_smix,v_sspec,v_ssite,v_sstr,v_sarry,sstrn)
#endif

#if LMF
      call lmfp(prgnam,v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .  v_smix,v_sspec,v_ssite,v_sstr,v_smove,v_sarry,sstrn)
#endif

#if LMMC
      if (cmdopt('--fit',5,0,strn)) then
        call lmcfit('MCFIT',v_sctrl,v_sham,v_spot,v_slat,
     .    v_sbz,v_smix,v_sspec,v_ssite,v_sstr,v_smove,
     .    v_sarry,sstrn)
      elseif (cmdopt('--atom',6,0,strn)) then
        call lmca('MCA',v_sctrl,v_sham,v_spot,v_slat,v_sbz,
     .    v_smix,v_sspec,v_ssite,v_sstr,v_smove,v_sarry,
     .    sstrn)
      elseif (cmdopt('--xbs',5,0,strn)) then
        call lmxbs('MCXBS',v_sctrl,v_sham,v_spot,v_slat,v_sbz,
     .    v_smix,v_sspec,v_ssite,v_sstr,v_smove,v_sarry,
     .    sstrn)
      elseif (cmdopt('--rho',5,0,strn)) then
        call lmrho('MCXBS',v_sctrl,v_sham,v_spot,v_slat,v_sbz,
     .    v_smix,v_sspec,v_ssite,v_sstr,v_smove,v_sarry,
     .    sstrn)
      else
        call lmce('MCE',v_sctrl,v_sham,v_spot,v_slat,v_sbz,
     .    v_smix,v_sspec,v_ssite,v_sstr,v_smove,v_sarry,
     .    sstrn)
      endif
#endif

#if LMSTR
      call asastr(prgnam,v_sctrl,v_slat,v_sspec,v_sstr,
     .  v_sarry,v_sham)
#endif

#if LMCTL
      i = fopnn('LOG')
      j = 0
      if (cmdopt('-mad',4,0,strn)) j = 1
      if (cmdopt('-spin1',6,0,strn)) j = j+10
      if (cmdopt('-spinf',6,0,strn)) j = j+20
      if (cmdopt('-enu',4,0,strn)) j = j+100
      call shoctl(v_sctrl,v_sspec,v_spot,v_sarry,j,i)
      call fexit(0,0,' ',0)
#endif

#if LMGF
Cgetarg       i = igets('ctrl lcgf',v_sctrl)
       i = int(v_sctrl%lcgf) 

C ... Special branch for exchange coupling
      if (i .ge. 10 .and. i .lt. 30) then
      call exasa(v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .  v_sspec,v_ssite,v_sarry,sstrn,slabl)
      else
C ... This lmasa must have gfasa uncommented
      call lmasa(prgnam,v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .  v_smix,v_sspec,v_ssite,v_sstr,v_sarry,sstrn)
      endif
#endif

#if LMDOS
      call asados(prgnam,v_sbz,v_sctrl,v_sham,v_slat,
     .  v_sspec,v_ssite,v_sarry)
#endif

#if TBE
      ltbe = 1
      if (cmdopt('--band',6,0,outs)) ltbe = 2
      if (cmdopt('--fit',5,0,outs)) ltbe = 4
      call tbzint(prgnam,v_sbz,v_sctrl,v_sham,v_spot,v_slat,
     .  v_smix,v_sspec,v_ssite,v_sstr,v_sarry,v_stb,sstrn,
     .  ltbe)
#endif

#if LMIMP
      call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
     .  v_ssite,v_sstr,v_sarry,v_sbz,slabl,sstrn,2**6,wksize)
#endif

#if MMAG
C     Uncomment these lines for Monte-Carlo atomistic simulations
C      call mcasim(v_sctrl,v_sarry,v_sspec,v_slat,v_smove,
C     .  sstrn)
      call mmag(v_sctrl,v_sham,v_sarry,v_sspec,v_slat,
     .  v_smove,sstrn)
#endif

#if LMPLAN
      call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
     .  v_ssite,v_sstr,v_sarry,v_sbz,slabl,sstrn,2**1,wksize)
#endif
#endif

C -------------- End of program -------------
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

Cgetarg       end
       end subroutine fmain 

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
     .  '%a [--OPTION] [-var-assign] [ext]')

      print 343
      print 344
  343 format(/' --h'/' --help',t17,'Print this message, and quit'
     .  /' --input',t17,
     .  'List categories, tokens, and data program expects, and quit'
     .  /' --show',t17,
     .  'Print control file after parsing by preprocessor,'/t17,
     .  'and echo input data as read from the control file'
     .  /' --showp',t17,
     .  'Same as --show, but quit after input parsed'
     .  /' --iactiv',t17,'(--no-iactiv) ',
     .  'Turn on (off) interactive mode'/
     .  t17,'This switch overrides input file setting',
     .  /' --pr#1[,#2...]',t17,
     .  'Set the verbosity (stack) to values #1,#2, ...'
     .  /' --time=#1[,#2]',t17,
     .  'Print timing info to # levels (#1=summary; #2=on-the-fly)'/
     .  /' -vnam=expr',t17,
     .  'Define numerical variable "nam"; set to result of ''expr'''
     .  /' -cnam=strn',t17,
     .  'Define character variable "nam"; set to ''strn'''
     .  )


  344 format(
     .  /' --rpos=filnam',t17,
     .  'After reading input file, read site positions from "filnam"'/
     .  ' --fixpos[:tol=#]',
     .  ' Adjust positions slightly, rendering them as'/
     .  t17,'exactly as possible consistent with the symmetry group')

      if (.true.) then
        outs = '%N '//prgnam//'%a-specific options:'
        call strip(outs,i1,i2)
        call info0(0,0,0,outs(1:i2))

        if (prgnam .eq. 'lmfa') then
          call info0(0,0,0,
     .      '%N%1f ')
        endif

        if (prgnam .eq. 'lmfgwd') then
          call info0(0,0,0,
     .      '%N%1f ')
        endif

        if (prgnam .eq. 'lmf') then
          call info0(0,0,0,
     .      '%N%1f--rs=#1,#2,#3,#4,#5'//
     .      '%N%6f#1=0 start from atm file; 1 from rst file;'//
     .      ' 2 from rsta file'//
     .      '%N%11fadd 10 to shift sm-rho 1st iter'//
     .      '%N%11fadd 100 to rotate local rho 1st iter'//
     .      '%N%6f#2=1 save rst file'//
     .      '%N%6f(#3,#4,#5)=0 read (pos,E_f,pnu) from rst file'//
     .      '%N%N --band --wden --rdbasp --optbas --etot'//
     .      ' (see documentation)')
        endif

        if (prgnam .eq. 'lmmc') then
          call info0(0,0,0,
     .      '%N%1f--atom --fit --rs --st'//
     .       '%N --atom invokes the free atom program'//
     .       '%N --fit  invokes the two-center fit'
     .      )
        endif

        if (prgnam .eq. 'lm') then
          call info0(0,0,0,
     .      '%N%1f--rs=#1,#2 -mix=#1[,#2] -onesp -map[0]')
        endif

        if (prgnam .eq. 'lmgf') then
          call info0(0,0,0,
     . '%N%1f-ef=#%N Exchange-mode-specific options:'//
     . '%N  --sites[:pair]:site-list'//
     . '%N Further switches for exchange printout mode:'//
     . '%N --wrsj[:j00][:amom][:sscl][:[g]scl=#][:tol=#] --wmfj -rcut=#'
     .      )
        endif

        if (prgnam .eq. 'lmdos') then
          call info0(0,0,0,
     .      '%N%1f--dos:options  modifies number and kinds of dos '//
     .      'generated;%N%16fsee documentation (doc/lmto.html)')
        endif

        if (prgnam .eq. 'lmchk') then
          call info0(0,0,0,
     .  '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
     .  '[:pairs:pair-list]...'//
     .  '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
     .  '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
     .  '%N --findes'//
     .  '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
        endif

        if (prgnam .eq. 'lmscell') then
          call info0(0,0,0,
     .  '%N%1f--wsite=name'
     .  //'%N%1f--wsitex=name'
     .  //'%N --wrsj[:fn=name][:scl=#]'
     .  //'%N --shorten'
     .  //'%N --first'
     .  //'%N --pl:expr'
     .  //'%N --sites:site-list'
     .  //'%N --sort:expr [expr] [expr]'
     .  //'%N --disp:fname:site-list')
C     . disp#rcut#expr-center-list#expr-neighbor-list
        endif

        if (prgnam .eq. 'lmstr') then
          call info0(0,0,0,
     .      '%N%1f--chk%11fcompares file str.ext with str1.ext'//
     .       '%N --plot[:con|:line[,v1x..z,v2x..z]|onec] '//
     .       'plots envelope'//
     .       '%N --pltg[:con|:line[,v1x..z,v2x..z]|onec] '//
     .       'plots envelope, val-lap')
        endif

        if (prgnam .eq. 'lmxbs') then
          call info0(0,0,0,
     .      '%N%1f-shift=x1,x1,x3 -spec '//
     .      '-dup=d1,d2,d3[,expr] -bs=val -ss=val')
        endif

        if (prgnam .eq. 'lmctl') then
          call info0(0,0,0,
     .      '%N%1f -spin1 -spinf -mad -enu')
        endif

        if (prgnam .eq. 'lmpg') then
          call info0(0,0,0,
     .      '%N%1f -map -onesp')
        endif

        if (prgnam .eq. 'lmplan') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'lmshf') then
          call info0(0,0,0,
     .      '%N%1f-enu=expr   linearize pot. pars around ''expr''')
        endif

        if (prgnam .eq. 'lmimp') then
          call info0(0,0,0,
     .      '%N%1f -rs -4 -5 -3s -4s -47u -5s (1 is required)')
        endif

        if (prgnam .eq. 'lmcor') then
          call info0(0,0,0,
     .      '%N%1f -findr')
        endif

        if (prgnam .eq. 'lmavgm') then
          call info0(0,0,0,
     .      '%N%1f -spin1')
        endif

        if (prgnam .eq. 'lmmix') then
          call info0(0,0,0,
     .      '%N%1f-fn=mix-file-name -bin2a or -a2bin')
        endif

        if (prgnam .eq. 'tbe') then
          call info0(0,0,0,
     .      '%N%1f--band --wpos=fnam -cont -dumph'//
     .      ' --st --md=# --mv=# --xyz=#'
     .    //'%N%3fThe last 4 switches apply to'//
     .      ' molecular dynamics simulations')
        endif

        if (prgnam .eq. 'mmag') then
          call info0(0,0,0,'%N%1f--cont --wrsj[:fn=name]')
        endif

      endif

      if (nint(vn2) .ne. 0) then
        ch = char(nint(vn2)+ichar('a')-1)
      else
        ch = ' '
      endif
      call info2(0,0,0,'%N version %,2d'//ch//'   worksize %i K',
     .  vn-vn2/1d5,wksize/1000)
      call fexit(0,0,' ',0)
Cgetarg       end
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

      character(len=9):: strn, ccc
      do
        read(nfilin,"(a)",err=1010,end=1010) strn
        if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
Cgetarg       end
       end subroutine findctrlstart 

