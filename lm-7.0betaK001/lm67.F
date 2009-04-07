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

C#define LM
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
C#define unix
C#ifdef unix
      subroutine fmain
C#endif
      implicit none

C#ifdefC MPI | MPIK
C      include "mpif.h"
C#ifdefC MPE
C      include "mpef.h"
C#endifC
C      integer numprocs, ierr, status(MPI_STATUS_SIZE)
C      integer procid, master
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C      double precision starttime, endtime
C      logical mlog
C#ifdefC MPE
CC Event numbers:
C      include "fp/events.ins"
C#endifC MPE
C#endif MPI | MPIK

C Heap allocation
      integer wksize
C#ifdefC MPI | MPIK
C      parameter(wksize= 120 000 000)
C#else
      parameter(wksize= 80 000 000)
C#endif
      integer w(wksize)
C     Next two lines guarantee w is aligned along a d.p. boundary
      double precision ws
      equivalence (ws,w(1))
      common /w/ w

C ... Controls for IO
      integer lstrn
      parameter (lstrn=2000)

      character prgnam*8, vrsion(2)*6, ext*20
C ... Structure for strings
      character*(lstrn) sstrn

C ... For structures
      logical lgors
      integer mxspec
      integer osbz,osctrl,oslat,osham,ospot,osmix,osspec,ossite,osarry,
     .  osmove,osstr,ostb
      parameter (mxspec=256)
      character*8 slabl(mxspec)

C ... miscellaneous local variables
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     .  igets,a2vec
      double precision vn(2),vn2(2),dglob
      logical T,F,swtmp,cmdopt,ltet,a2bin
      parameter (T=.true., F=.false.)
C#ifdefC LMCTL
C      integer fopnn
C#endif
C#ifdefC MPI | MPIK
C      integer fext
C#endif
C#ifdefC TBE
C      double precision mdprm(6)
C#endif

C ... Program-dependent name and help
C#ifdefC LMPG
C#elseif LM
C#elseifC LMFA
C      data prgnam /'LMFA'/
C#elseifC LMFGWD
C      data prgnam /'LMFGWD'/
C#elseifC LMF
C      data prgnam /'LMF'/
C#elseifC LMAQU
C      data prgnam /'LMAQU'/
C#elseifC  LMDOS
C      data prgnam /'LMDOS'/
C#elseifC  LMCTL
C      data prgnam /'LMCTL'/
C#elseifC  LMIMP
C      data prgnam /'LMIMP'/
C#elseifC  LMCOR
C      data prgnam /'LMCOR'/
C#elseifC  LMFIT
C      data prgnam /'LMFIT'/
C#elseifC  LMPLAN
C      integer nbas,nbasp,nbaspp,nl,opnu,oqnu,oqc,oqt,oves
C      double precision emad,trumad,vmtz(2)
C      data prgnam /'LMPLAN'/
C#elseifC  LMAVGM
C      data prgnam /'LMAVGM'/
C#elseifC  LMMIX
C      data prgnam /'LMMIX'/
C#elseifC  LMSHF
C      data prgnam /'LMSHF'/
C      data auxmod /8/
C#elseifC  TBPG
C      data prgnam /'TBPG'/
C#elseifC  TBE
C      data prgnam /'TBE'/
C#elseifC  TBBND
C      data prgnam /'TBBND'/
C#elseifC  TBFIT
C      data prgnam /'TBFIT'/
C#elseifC  TBDOS
C      data prgnam /'TBDOS'/
C#elseifC  MMAG
C      data prgnam /'MMAG'/
C#endif

C ... Program-dependent cagetories
C#ifdefC LMPG
C      data prgnam /'LMPG'/
C#elseifC LMGPOL
C      data prgnam /'LMGPOL'/
C      data auxmod /32/
C#elseifC BLM
C      data prgnam /'BLM'/
C#elseif LM
      data prgnam /'LM'/
C#elseifC LMMC
C      data prgnam /'LMMC'/
C#elseifC LMAQU | LMMIX
C      call rx('missing setup for lmaqu, lmmix')  ! Why???TK
C#elseifC LMGF
C      data prgnam /'LMGF'/
C#elseifC LMSTR
C      data prgnam /'LMSTR'/
C#elseifC LMCHK
C      data prgnam /'LMCHK'/ auxmod /1/
C#elseifC LMSCELL
C      integer nkd,odlv,nbx
C      double precision plat(3,3),plx(3,3),ddot,dmach
C      data prgnam /'LMSCELL'/
C#elseifC LMXBS
C      data prgnam /'LMXBS'/
C      data auxmod /4/
C#endif

      integer:: nfilin,offset,mxrecs,recln0,nrecs
     &  ,recoff,iprint,fopna
      parameter (mxrecs=10000,recln0=120) !this is also in rdctrl2
      character*8 alabl
c      character*1 recrd(0:1)
c      equivalence( recrd(0),w(1))
      character*1,allocatable:: recrd(:)


C -------------- First executable statement ---------------
C     Link in routines that cause potential library conflicts
      call nada

c toksw was here. MOved to rdccat

C --- Version ---
      vrsion(1) = 'LM'
      vn(1) = 7.00d0
      vn2(1) = 0
      vn(2) = 0.00d0
      vn2(2) = 0
C#ifdef LM | LMPG | LMGF | LMAVGM | LMIMP | LMCTL | LMSTR
      vrsion(2) = 'ASA'
      vn(2) = 7.00d0
      vn2(2) = 0
C#elseifC LMFA | LMF | LMFGWD
C      vrsion(2) = 'FP'
C      vn(2) = 7.00d0
C      vn2(2) = 0
C#elseifC LMMC
C      vrsion(2) = 'MOL'
C      vn(2) = 3.00d0
C      vn2(2) = 0
C#elseifC TBE
C      vrsion(2) = 'TB'
C      vn(2) = 9.00d0
C      vn2(2) = 0
C#elseifC MMAG
C      vrsion(2) = 'MM'
C      vn(2) = 2.00d0
C      vn2(2) = 0
C#endif
      vn(1) = vn(1) + vn2(1)/1d5
      vn(2) = vn(2) + vn2(2)/1d5
      stdo = lgunit(1)
      i = dglob('stdo',dble(stdo),1)

C --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs))
     .  call lmhelp(prgnam,vn(2),vn2(2),wksize)
      if (cmdopt('--version',9,0,outs)) then
        write(*,'(f4.2)') vn(2)
        call cexit(0,1)
      endif

C --- Dynamic memory allocation and other initialization ---
C#ifdefC MPI | MPIK
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      mlog = cmdopt('--mlog',6,0,strn)
C#ifdefC MPE
C      ierr = MPE_INIT_LOG()
C      EVENT_START_RDCTRL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_RDCTRL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_UGCOMP = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_UGCOMP   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_FSMBL  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_FSMBL    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_PZHEV  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_PZHEV    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_KLOOP  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_KLOOP    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_MIXRHO = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_MIXRHO   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_SMHSBL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_SMHSBL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_AUGMBL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_AUGMBL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_HSIBL  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_HSIBL    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_RSIBL  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_RSIBL    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_RLOCBL = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_RLOCBL   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_DFRCE  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_DFRCE    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_BCAST  = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_BCAST    = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_ALLRED = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_ALLRED   = MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_START_BARRIER= MPE_LOG_GET_EVENT_NUMBER()
C      EVENT_END_BARRIER  = MPE_LOG_GET_EVENT_NUMBER()
C#endifC
CC --- Dynamic memory allocation and other initialization ---
C      if (procid .eq. master) call headl2(prgnam,wksize,stdo)
C      call finits(2,0,0,i)
C      call pshpr(0)
C      call wkinit(wksize)
C      call wkfast(T)
C      call poppr
C      if (procid .eq. master) then
C        i = fextg(ext)
C      endif
C      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
C      if (procid .eq. master) then
C#ifdefC MPE
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_RDCTRL,EVENT_END_RDCTRL,
C     .                            "rdctrl","pink")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_UGCOMP,EVENT_END_UGCOMP,
C     .                            "ugcomp","maroon")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_FSMBL,EVENT_END_FSMBL,
C     .                            "fsmbl","aquamarine")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_PZHEV,EVENT_END_PZHEV,
C     .                            "pzhev","brown")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_KLOOP,EVENT_END_KLOOP,
C     .                            "k-loop","brown")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_MIXRHO,EVENT_END_MIXRHO,
C     .                            "mixrho","orange")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_SMHSBL,EVENT_END_SMHSBL,
C     .                            "smhsbl","blue")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_AUGMBL,EVENT_END_AUGMBL,
C     .                            "augmbl","cyan")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_HSIBL,EVENT_END_HSIBL,
C     .                            "hsibl","gray")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_RSIBL,EVENT_END_RSIBL,
C     .                            "rsibl","red")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_RLOCBL,EVENT_END_RLOCBL,
C     .                            "rlocbl","green")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_DFRCE,EVENT_END_DFRCE,
C     .                            "dfrce","magenta")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_BCAST,EVENT_END_BCAST,
C     .                            "broadcast","coral")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_ALLRED,EVENT_END_ALLRED,
C     .                            "allreduce","purple")
C        ierr = MPE_DESCRIBE_STATE(EVENT_START_BARRIER,EVENT_END_BARRIER,
C     .                            "barrier","yellow")
C#endifC
C        call gettime(datim)
C        if (mlog) i = fopn('MLOG')
C        call poseof(fopn('LOG'))
C        if (mlog) then
C          call awrit2(' lmf '//datim//' Process %i of %i on '
C     .      //shortname(procid)(1:namelen(procid))//
C     .      ' is master',' ',256,
C     .      lgunit(3),procid,numprocs)
C        endif
C      else
C        call strcat(ext,20,' ','_',1,' ',i)
C        call bin2a(' ',0,0,procid,2,0,20,ext,i)
C        ierr = fext(ext(1:i+1))
C        if (mlog) ierr = fopn('MLOG')
C        ierr = fextg(ext)
C        call gettime(datim)
C        if (mlog) then
C          call awrit2(' lmf '//datim//' Process %i of %i on '
C     .      //shortname(procid)(1:namelen(procid))//
C     .      ' file extension is '//ext(2:i+1),' ',
C     .      256,lgunit(3),procid,numprocs)
C        endif
C      endif
C      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#else
C      procid = 0
C      master = 0

      j = 0
      swtmp = a2bin('T',ltet,1,0,' ',j,-1)
      call addsyv('version',vn,i)
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
C#endif

C#ifdefC BLM
C      call suctrl(prgnam)
C      call fexit(0,0,' ',0)
C#else

C ... Abort with error message if ctrl file is missing
C#ifdefC MPI | MPIK
C      if (procid .eq. master) then
C#endif

      if (cmdopt('--input',6,0,strn)) then
c        if (ncount(1) .eq. 0) ncount(1) = -1
      else
        if (fxst('CTRL') .ne. 1) then
          call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     .      i1mach(2))
          call cexit(-1,1)
        endif
      endif
C#ifdefC MPI | MPIK
C      endif
C#endif

C ... Set special file directory for temporary files
C     User may which to customize the directory
C     Default is to use the standard directory
C     sttmpd is located at the bottom of this file.
      call sttmpd

C ... File logical units
      i = fadd('TMP',-1,4)
      i = fadd('BAND',-1,4)
C#ifdefC TBE
C      i = fadd('STRT',-1,4)
C      i = fadd('QMOM',-1,4)
C#endif
C#ifdef LM | LMDOS | LMGF
      i = fadd('MOMS',-1,4)
      i = fadd('CLS',-1,4)
C#endif

C --- Set the top-level verbosity if spec'd from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif

C --- Converter CTRL+console into recrd  ---
C   recrd, nrecs are obtained.
      nrecs  = 0
      allocate( recrd( 0:mxrecs*recln0-1 ) )
      if (.not.cmdopt('--input',7,0,strn)) then
        nfilin = fopna('CTRL',-1,1)
ctakao
        call findctrlstart(nfilin)

        alabl = '#{}% ct '
        if (cmdopt('--show',6,0,strn)) alabl = '#{}% ctp'
        call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
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
c      write(6,*)' === TEST out recrd === '
c      ioff=0
c      do i=1,nrecs
c        write(6,"(256a1)") recrd(ioff:ioff+recln0-1)
c        ioff= ioff + recln0
c      enddo
c      write(6,*)' --- end of rdfile --- '


C --- Read recrd ---
      call rdctrlchk(recrd,recln0,nrecs,
     .  prgnam,vrsion,vn,vn2,F,
     .  slabl,osbz,osctrl,osham,ospot,oslat,osmix,osspec,ossite,osstr,
     .  osarry,osmove,ostb,sstrn)
      deallocate(recrd)
      stop 'at lmft'

c      call rdctrl3(recrd,recln0,nrecs,
c     .  prgnam,vrsion,vn,vn2)


C#ifdefC MPI | MPIK
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_RDCTRL,procid,"rdctrl")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endifC
C      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C#endifC
C      if (procid .ne. master) then
C        call pshpr(0)
C        call pshpr(0)
C        do  i = 1, 4
C          call sprt(i,0)
C        enddo
C      endif
C#endif

C#ifdefC LMSCELL
C      call setcg(w(oslat),8,12)
C      call pack5('lat as tol nkdmx nkqmx',w(oslat),1d0,dmach(2),150000,
C     .  150000,0)
C      call lattic(w(oslat),w(osctrl),w(ossite),w(osarry))
C      call upack('lat nkd odlv plat slat',w(oslat),nkd,odlv,plat,plx,0)
C
CC ... Set default values for species data
C      call defspc(w(osspec))
C
C      if (ddot(9,plx,1,plx,1) .eq. 0) then
C        write(stdo,'(/a)') ' Enter basis vectors of supercell:'
C        read(5,*)  plx
C        if (ddot(9,plx,1,plx,1) .eq. 0) call dcopy(9,plat,1,plx,1)
C      endif
C
C      call supcel(0,w(osctrl),w(oslat),w(ossite),slabl,w(osspec),plx,
C     .  nkd,w(odlv),nbx)
C      goto 1000
C#endif

C --- Lattice setup ---
C#ifndef LMCTL | LMMC
      call setcg(w(oslat),8,12)
      call lattic(w(oslat),w(osctrl),w(ossite),w(osarry))
C#endif

C --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      call upacks('strn symg',i,j)
      if (j .ge. i) strn = sstrn(i:j+1)
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      lc = 20
C#ifdefC TBE
C      call upack('ctrl mdprm',w(osctrl),mdprm,0,0,0,0)
C      if (nint(mdprm(1)) .eq. 1 .or. nint(mdprm(1)) .eq. 2) then
C        strn = ' '
C      endif
C      lc = 40
C#endif
C#ifdefC MMAG
C      lc = 0
C#endif
C#ifdefC LMFA
CC      lc = 0
C#endif
C#ifdefC LMF
C      lc = 10
C#endif
C#ifndef LMFGWD
      if (.not. lgors('ctrl lqp,1',w(osctrl))) lc = lc+2
C#endif
C#ifndef LMMC
      call mksym(lc,slabl,strn,w(osctrl),w(oslat),w(ossite),w(osarry))
C#endif

C --- Allocate permanent class arrays, maps and other initialization ---
      if (lc .ge. 20)
     .  call clsprm(1,w(osctrl),w(osarry),w(osham),w(ospot),w(osspec),
     .  w(oslat),w(osbz),w(osstr),sstrn)

C --- Read available class parameters from file ---
C#ifdefC LMCTL | LMCHK
C      call aiocls(.false.,0,w(osarry),w(osctrl),w(osham),w(ospot),
C     .  w(osspec),w(oslat),1,0)
C#elseif LM | TBE | LMGF | LMSHF | LMPLAN
C#ifndef TBE
      call aiocls(.false.,12,w(osarry),w(osctrl),w(osham),w(ospot),
     .  w(osspec),w(oslat),1,0)
C#endif
      call rdctrl2(recrd,recln0,nrecs,
     .  prgnam,wksize,vrsion,vn,vn2,T,
     .  slabl,osbz,osctrl,osham,ospot,oslat,osmix,osspec,ossite,osstr,
     .  osarry,osmove,sstrn)
C#ifndef TBE
      call aiocls(.false.,17,w(osarry),w(osctrl),w(osham),w(ospot),
     .  w(osspec),w(oslat),1,0)
      if (lc .ge. 20)
     .  call clsprp(1,w(osctrl),w(osarry),w(osham),w(ospot),w(osspec),
     .  w(oslat),w(osbz),w(osstr),sstrn)
C#endif
C#endif

C --- Optionally resize spheres ---
C#ifdefC RESIZE
C      if (cmdopt('--sfill',7,0,strn)) then
C        print *, '*** --sfill DEPRECATED .. see lmto.html'
C        call sfill(strn(8:),slabl,w(osctrl),w(osarry),w(oslat),
C     .    w(osspec),w(ossite))
C      else
C        call sfill('auto',slabl,w(osctrl),w(osarry),w(oslat),w(osspec),
C     .    w(ossite))
C      endif
C#endif

C ... Set default values for species data
      call defspc(w(osspec))

C ... Patch for now ... maybe replace
C#ifdef LM | TBE | LMGF | LMF
C#ifndef LMFGWD
      ltet = lgors('ctrl lmet,1',w(osctrl)) .or.
     .       lgors('ctrl ldos,4+2+1',w(osbz))
      call mkqp(w(osctrl),w(osbz),w(oslat),ltet,F,1,-2)
C#endif
C#endif

C ... quit after SHOW
      if (igets('ctrl quit',w(osctrl)) .eq. 1) then
        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif

C#ifdefC LMFA
C      call freeat(w(osctrl),w(osspec),w(osham),w(ospot))
C#endif

C#ifdefC LMCHK & FP
C      if (cmdopt('--fp',4,0,strn)) then
C        call fpchk(w(osspec),w(ossite))
C        call cexit(0,1)
C      endif
C#endif
C#ifdefC LMCHK | LMXBS | LMGPOL | LMSHF
C      call lmaux(prgnam,w(osctrl),w(osham),w(ospot),w(oslat),w(osspec),
C     .  w(ossite),w(osstr),w(osarry),w(osbz),slabl,sstrn,auxmod,wksize)
C#endif

C#ifdef LM
      call lmasa(prgnam,w(osbz),w(osctrl),w(osham),w(ospot),w(oslat),
     .  w(osmix),w(osspec),w(ossite),w(osstr),w(osarry),sstrn)
C#endif

C#ifdefC LMF
C      call lmfp(prgnam,w(osbz),w(osctrl),w(osham),w(ospot),w(oslat),
C     .  w(osmix),w(osspec),w(ossite),w(osstr),w(osmove),w(osarry),sstrn)
C#endif

C#ifdefC LMMC
C      if (cmdopt('--fit',5,0,strn)) then
C        call lmcfit('MCFIT',w(osctrl),w(osham),w(ospot),w(oslat),
C     .    w(osmix),w(osspec),w(ossite),w(osstr),w(osmove),w(osarry),
C     .    sstrn)
C      elseif (cmdopt('--atom',6,0,strn)) then
C        call lmca('MCA',w(osctrl),w(osham),w(ospot),w(oslat),
C     .    w(osmix),w(osspec),w(ossite),w(osstr),w(osmove),w(osarry),
C     .    sstrn)
C      elseif (cmdopt('--xbs',5,0,strn)) then
C        call lmxbs('MCXBS',w(osctrl),w(osham),w(ospot),w(oslat),
C     .    w(osmix),w(osspec),w(ossite),w(osstr),w(osmove),w(osarry),
C     .    sstrn)
C      elseif (cmdopt('--rho',5,0,strn)) then
C        call lmrho('MCXBS',w(osctrl),w(osham),w(ospot),w(oslat),
C     .    w(osmix),w(osspec),w(ossite),w(osstr),w(osmove),w(osarry),
C     .    sstrn)
C      else
C        call lmce('MCE',w(osctrl),w(osham),w(ospot),w(oslat),
C     .    w(osmix),w(osspec),w(ossite),w(osstr),w(osmove),w(osarry),
C     .    sstrn)
C      endif
C#endif

C#ifdefC LMSTR
C      call asastr(prgnam,w(osctrl),w(oslat),w(osspec),w(osstr),
C     .  w(osarry),w(osham))
C#endif

C#ifdefC LMCTL
C      i = fopnn('LOG')
C      j = 0
C      if (cmdopt('-mad',4,0,strn)) j = 1
C      if (cmdopt('-spin1',6,0,strn)) j = j+10
C      if (cmdopt('-spinf',6,0,strn)) j = j+20
C      if (cmdopt('-enu',4,0,strn)) j = j+100
C      call shoctl(w(osctrl),w(osspec),w(ospot),w(osarry),j,i)
C      call fexit(0,0,' ',0)
C#endif

C#ifdefC LMGF
C      i = igets('ctrl lcgf',w(osctrl))
CC ... Special branch for exchange coupling
C      if (i .ge. 10 .and. i .le. 13) then
C      call exasa(w(osbz),w(osctrl),w(osham),w(ospot),w(oslat),
C     .  w(osspec),w(ossite),w(osarry),sstrn,slabl)
C      else
CC ... This lmasa must have gfasa uncommented
C      call lmasa(prgnam,w(osbz),w(osctrl),w(osham),w(ospot),w(oslat),
C     .  w(osmix),w(osspec),w(ossite),w(osstr),w(osarry),sstrn)
C      endif
C#endif

C#ifdefC LMDOS
C      call asados(prgnam,w(osbz),w(osctrl),w(osham),w(oslat),
C     .  w(osspec),w(ossite),w(osarry))
C#endif

C#ifdefC TBE
C      call tbzint(prgnam,w(osbz),w(osctrl),w(osham),w(ospot),w(oslat),
C     .  w(osmix),w(osspec),w(ossite),w(osstr),w(osarry),sstrn)
C#endif

C#ifdefC LMIMP
C      call lmaux(prgnam,w(osctrl),w(osham),w(ospot),w(oslat),w(osspec),
C     .  w(ossite),w(osstr),w(osarry),w(osbz),slabl,sstrn,2**6,wksize)
C#endif

C#ifdefC MMAG
CC     Uncomment these lines for Monte-Carlo atomistic simulations
CC      call mcasim(w(osctrl),w(osarry),w(osspec),w(oslat),w(osmove),
CC     .  sstrn)
C      call mmag(w(osctrl),w(osham),w(osarry),w(osspec),w(oslat),
C     .  w(osmove),sstrn)
C#endif

C#ifdefC LMPLAN
C      call lmaux(prgnam,w(osctrl),w(osham),w(ospot),w(oslat),w(osspec),
C     .  w(ossite),w(osstr),w(osarry),w(osbz),slabl,sstrn,2**1,wksize)
C#endif
C#endif

C -------------- End of program -------------
 1000 continue
C#ifdefC MPI | MPIK
C#ifdefC MPE
C      if (procid .eq. master) i = fextg(ext)
C      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
C      i = 0
C      call skp2bl(ext,20,i)
C      ierr = MPE_FINISH_LOG(ext(2:i))
C#endifC
C      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C      if ( procid .eq. master ) then
C        call rx0(prgnam//' on '//shortname(procid)(1:namelen(procid)))
C      else
C        call fexit(0,0,' ',0)
C      endif
C#else
      call rx0(prgnam)
C#endif

      end
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
C     implicit none
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

        if (prgnam .eq. 'lmaqu') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'lmctl') then
          call info0(0,0,0,
     .      '%N%1f -spin1 -spinf -mad -enu')
        endif

        if (prgnam .eq. 'lmimp') then
          call info0(0,0,0,
     .      '%N%1f -rs -4 -5 -3s -4s -47u -5s (1 is required)')
        endif

        if (prgnam .eq. 'lmcor') then
          call info0(0,0,0,
     .      '%N%1f -findr')
        endif

        if (prgnam .eq. 'lmfit') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'lmplan') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'lmavgm') then
          call info0(0,0,0,
     .      '%N%1f -spin1')
        endif

        if (prgnam .eq. 'lmmix') then
          call info0(0,0,0,
     .      '%N%1f-fn=mix-file-name -bin2a or -a2bin')
        endif

        if (prgnam .eq. 'tbpg') then
          call info0(0,0,0,
     .      '%N%1f-cont -dumph')
        endif

        if (prgnam .eq. 'tbe') then
          call info0(0,0,0,
     .      '%N%1f-cont -dumph --rpos=fnam --wpos=fnam')
        endif

        if (prgnam .eq. 'tbbnd') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'tbfit') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'lmpg') then
          call info0(0,0,0,
     .      '%N%1f -map -onesp')
        endif

        if (prgnam .eq. 'lmgpol') then
          call info0(0,0,0,
     .      '%N%1f---')
        endif

        if (prgnam .eq. 'lmmc') then
          call info0(0,0,0,
     .      '%N%1f--atom --fit --rs --st'//
     .       '%N --atom invokes the free atom program'//
     .       '%N --fit  invokes the two-center fit'
     .      )
        endif

        if (prgnam .eq. 'lmstr') then
          call info0(0,0,0,
     .      '%N%1f--chk%11fcompares file str.ext with str1.ext'//
     .       '%N --plot[:con|:line[,v1x..z,v2x..z]|onec]'//
     .       '%N%16fplots envelope')
        endif

        if (prgnam .eq. 'lmshf') then
          call info0(0,0,0,
     .      '%N%1f-enu=expr   linearize pot. pars around ''expr''')
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

        if (prgnam .eq. 'lmdos' .or. prgnam .eq. 'tbdos') then
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
     .  '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
        endif

        if (prgnam .eq. 'lmscell') then
          call info0(0,0,0,
     .  '%N%1f--wsite=name'
     .  //'%N --wrsj[:fn=name][:scl=#]'
     .  //'%N --shorten'
     .  //'%N --first'
     .  //'%N --pl:expr'
     .  //'%N --sites:site-list'
     .  //'%N --sort:expr [expr] [expr]'
     .  //'%N --disp:fname:site-list')
C     . disp#rcut#expr-center-list#expr-neighbor-list
        endif

        if (prgnam .eq. 'lmxbs') then
          call info0(0,0,0,
     .      '%N%1f-shift=x1,x1,x3 -spec '//
     .      '-dup=d1,d2,d3[,expr] -bs=val -ss=val')
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
      end
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
      end


c      subroutine set_optio( optio_in )
c      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
c     .        iend,ichoos,nchoos,optio,nx ,optio_in
c      logical noerr, ioc,a2bin
c      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
c     .               iend,ichoos,nchoos,noerr,optio
c      optio=optio_in
c      end

      subroutine findctrlstart(nfilin)
      character(len=9):: strn, ccc
      do
        read(nfilin,"(a)",err=1010,end=1010) strn
        if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
      end
