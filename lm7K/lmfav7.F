Cgetarg...info...           structure ['', 'ctrl', 'v_sctrl']
      program lmfa
      use m_struc_def
      use m_rdctrl2_func, only: rdctrl2
      use m_globalvariables
      implicit none

      integer procid, master, mpipid, nproc


C Heap allocation
c      integer wksize
c
c      parameter(wksize= 80 000 000)
c
c      integer w(wksize)
C     Next two lines guarantee w is aligned along a d.p. boundary
c      double precision ws
c      equivalence (ws,w(1))
c      common /w/ w

C ... Controls for IO
      integer lstrn
      parameter (lstrn=20000)

      character prgnam*8, vrsion(2)*6, ext*100
C ... Structure for strings
      character*(lstrn) sstrn

C ... For structures
      logical lgors
      integer mxspec
CKi      integer v_sbz,v_sctrl,v_slat,v_sham,v_spot,v_smix,v_sspec,v_ssite,v_sarry,
CKi     *  v_smove,v_sstr,v_stb
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

ckino Jan.04.2012:        type(s_spec),pointer:: v_sspec(:) =>NULL()
      type(s_spec),allocatable:: v_sspec(:) 
ckino Jan.04.2012:        type(s_site),pointer:: v_ssite(:) =>NULL()
      type(s_site),allocatable:: v_ssite(:) 


      parameter (mxspec=256)
      character*8 slabl(mxspec)

C ... miscellaneous local variables
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     *igets,a2vec,NULLI
      double precision vn(2),vn2(2),dglob
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)




C ... Program-dependent name and help

      data prgnam /'LMFA'/


C ... Program-dependent cagetories


      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=120)
      character*8 alabl
c      character*1 recrd(0:1)
c      equivalence( recrd(0),w(1))
CKi      character*1,allocatable:: recrd(:)
      character:: recrd_(0:mxrecs*recln0-1)
      character*(mxrecs*recln0):: recrd
      equivalence ( recrd, recrd_ )


      real(8),parameter::    NULLR =-99999
cMay1 patch for MPI from y.Yoshimoto
#if MPI|MPIK
      include 'mpif.h'
      integer ierr
#endif

C -------------- First executable statement ---------------
C     Link in routines that cause potential library conflicts
c      call nada

c toksw was here. MOved to rdccat

C --- Version ---
      call show_programinfo(6)

      vrsion(1) = 'LM'
      vrsion(2) = ' '
      vn(1) = 7.00d0
      vn2(1) = 0
      vn(2) = 0.00d0
      vn2(2) = 0

      vrsion(2) = 'FP'
      vn(2) = 7.00d0
      vn2(2) = 0

      vn(1) = vn(1) + vn2(1)/1d5
      vn(2) = vn(2) + vn2(2)/1d5
      stdo = lgunit(1)
      i = dglob('stdo',dble(stdo),1)
#if MPI|MPIK
      call mpi_init(ierr)
#endif
      master = 0
      procid = mpipid(1)
      nproc  = mpipid(0)

C --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs))
c     *call lmhelp(prgnam,vn(2),vn2(2),wksize)
     *call lmhelp(prgnam,vn(2),vn2(2),0)
      if (cmdopt('--version',9,0,outs)) then
        write(*,'(f4.2)') vn(2)
        call cexit(0,1)
      endif

C --- Dynamic memory allocation and other initialization ---

C     Add version to variables table?
C     j = 0
C     swtmp = a2bin('T',ltet,1,0,' ',j,-1)
C     call addsyv('version',vn,i)
c      call headl2(prgnam,wksize,stdo)
      call headl2(prgnam,0,stdo)
      call finits(2,0,0,i)
c      call pshpr(0)
c      call wkinit(wksize)
c      call wkfast(T)
c      call poppr
      call poseof(fopn('LOG'))
      i = fextg(ext)
      call word(ext,1,i,j)
      if (ext(i:i) .eq. '.') i=i+1
      if (ext(j:j) .eq. '.') j=j-1
      if (j .ge. i) call ptenv('EXT='//ext(i:j))




C ... Abort with error message if ctrl file is missing (swtmp = .true.)
      swtmp = .false.

      if (cmdopt('--input',6,0,strn)) then
        if (nproc .gt. 0) call rx('--input not allowed with MPI')
      else
        if (fxst('CTRL') .ne. 1) then
          call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     *    i1mach(2))
          swtmp = .true.
        endif
      endif

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



C --- Set the top-level verbosity if spec'd from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif

C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
CKi      allocate( recrd( 0:mxrecs*recln0-1 ) )
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
c        stop 'xxxxxxxxxxxxxxxxx'

          i = 60
          if (cmdopt('--show',6,0,strn)) then
            i = 1
            call info0(i,0,0,' ------------------------ End '//
     *      'of input file ----------------------')
          endif
          call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     *    ' input file',nrecs,0)
          if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
        endif
      endif


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
ckino Jan.04.2012:       *slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
ckino Jan.04.2012:       *slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     *slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     *v_sarry,v_smove,v_stb,sstrn)

c

      print *,'end of rdctrl2 in imfav7'
C ... Input to scale MT radii (define RESIZE)




C --- Lattice setup ---

ctakao will add |1)
      call setcg(v_slat,8,12)
ckino Jan.04.2012:        call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)


C --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      call upacks('strn symg',i,j)
      if (j .ge. i) strn = sstrn(i:j+1)
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      lc = 20

Cgetarg       if (.not. lgors('ctrl lqp,1',v_sctrl)) lc = lc+2
      if ( .not. iand(1,int(v_sctrl%lqp)) .ne.0 ) lc = lc+2 



ctakao will add |1)
      print *,'goto mksym'
ckino Jan.04.2012:        call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
      print *,'end of mksym x'

C --- Allocate permanent class arrays, maps and other initialization ---
c      if (lc .ge. 20)
c     *  call clsprm(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
c     *  v_slat,v_sbz,v_sstr,sstrn)
c      print *,'end of clsprm'

C --- Read available class parameters from file ---

C     Discard contents of input file
Cki      deallocate(recrd)

C --- Optionally resize spheres ---

C$$$      if (cmdopt('--sfill',7,0,strn)) then
C$$$        print *, '*** --sfill DEPRECATED .. see lmto.html'
C$$$        call sfill(strn(8:),slabl,v_sctrl,v_sarry,v_slat,
C$$$     *    v_sspec,v_ssite)
C$$$      else
C$$$        call sfill('auto',slabl,v_sctrl,v_sarry,v_slat,v_sspec,
C$$$     *    v_ssite)
C$$$      endif


C ... Set default values for species data
      print *,'goto defspc'
ckino Jan.04.2012:        call defspc(v_sspec)
      call defspc(v_sspec)
      print *,'end of defspc'

C ... Patch for now ... maybe replace


C ... quit after SHOW
Cgetarg       if (igets('ctrl quit',v_sctrl) .eq. 1) then
      if ( int(v_sctrl%quit) .eq. 1 ) then 

        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif
      print *,'goto freeat'
ckino Jan.04.2012:        call freeat(v_sctrl,v_sspec,v_sham,v_spot)
      call freeat(v_sctrl,v_sspec,v_sham,v_spot)
      call rx0(prgnam)
Cgetarg       end
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
     *'%a [--OPTION] [-var-assign] [ext]')

      print 343
      print 344
  343 format(/' --h'/' --help',t17,'Print this message, and quit'
     */' --input',t17,
     *'List categories, tokens, and data program expects, and quit'
     */' --show',t17,
     *'Print control file after parsing by preprocessor,'/t17,
     *'and echo input data as read from the control file'
     */' --showp',t17,
     *'Same as --show, but quit after input parsed'
     */' --iactiv',t17,'(--no-iactiv) ',
     *'Turn on (off) interactive mode'/
     *t17,'This switch overrides input file setting',
     */' --pr#1[,#2...]',t17,
     *'Set the verbosity (stack) to values #1,#2, ...'
     */' --time=#1[,#2]',t17,
     *'Print timing info to # levels (#1=summary; #2=on-the-fly)'/
     */' -vnam=expr',t17,
     *'Define numerical variable "nam"; set to result of ''expr'''
     */' -cnam=strn',t17,
     *'Define character variable "nam"; set to ''strn'''
     *)


  344 format(
     */' --rpos=filnam',t17,
     *'After reading input file, read site positions from "filnam"'/
     *' --fixpos[:tol=#]',
     *' Adjust positions slightly, rendering them as'/
     *t17,'exactly as possible consistent with the symmetry group')

      if (.true.) then
        outs = '%N '//prgnam//'%a-specific options:'
        call strip(outs,i1,i2)
        call info0(0,0,0,outs(1:i2))

        if (prgnam .eq. 'lmfa') then
          call info0(0,0,0,
     *    '%N%1f ')
        endif

        if (prgnam .eq. 'lmfgwd') then
          call info0(0,0,0,
     *    '%N%1f ')
        endif

        if (prgnam .eq. 'lmf') then
          call info0(0,0,0,
     *    '%N%1f--rs=#1,#2,#3,#4,#5'//
     *    '%N%6f#1=0 start from atm file; 1 from rst file;'//
     *    ' 2 from rsta file'//
     *    '%N%11fadd 10 to shift sm-rho 1st iter'//
     *    '%N%11fadd 100 to rotate local rho 1st iter'//
     *    '%N%6f#2=1 save rst file'//
     *    '%N%6f(#3,#4,#5)=0 read (pos,E_f,pnu) from rst file'//
     *    '%N%N --band --wden --rdbasp --optbas --etot'//
     *    ' (see documentation)')
        endif

        if (prgnam .eq. 'lmmc') then
          call info0(0,0,0,
     *    '%N%1f--atom --fit --rs --st'//
     *    '%N --atom invokes the free atom program'//
     *    '%N --fit  invokes the two-center fit'
     *    )
        endif

        if (prgnam .eq. 'lm') then
          call info0(0,0,0,
     *    '%N%1f--rs=#1,#2 -mix=#1[,#2] -onesp -map[0]')
        endif

        if (prgnam .eq. 'lmgf') then
          call info0(0,0,0,
     *    '%N%1f-ef=#%N Exchange-mode-specific options:'//
     *    '%N  --sites[:pair]:site-list'//
     *    '%N Further switches for exchange printout mode:'//
     *    '%N --wrsj[:j00][:amom][:sscl][:[g]scl=#][:tol=#] --wmfj -rcut=#'
     *    )
        endif

        if (prgnam .eq. 'lmdos') then
          call info0(0,0,0,
     *    '%N%1f--dos:options  modifies number and kinds of dos '//
     *    'generated;%N%16fsee documentation (doc/lmto.html)')
        endif

        if (prgnam .eq. 'lmchk') then
          call info0(0,0,0,
     *    '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
     *    '[:pairs:pair-list]...'//
     *    '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
     *    '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
     *    '%N --findes'//
     *    '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
        endif

        if (prgnam .eq. 'lmscell') then
          call info0(0,0,0,
     *    '%N%1f--wsite=name'
     *    //'%N%1f--wsitex=name'
     *    //'%N --wrsj[:fn=name][:scl=#]'
     *    //'%N --shorten'
     *    //'%N --first'
     *    //'%N --pl:expr'
     *    //'%N --sites:site-list'
     *    //'%N --sort:expr [expr] [expr]'
     *    //'%N --disp:fname:site-list')
C     . disp#rcut#expr-center-list#expr-neighbor-list
        endif

        if (prgnam .eq. 'lmstr') then
          call info0(0,0,0,
     *    '%N%1f--chk%11fcompares file str.ext with str1.ext'//
     *    '%N --plot[:con|:line[,v1x..z,v2x..z]|onec] '//
     *    'plots envelope'//
     *    '%N --pltg[:con|:line[,v1x..z,v2x..z]|onec] '//
     *    'plots envelope, val-lap')
        endif

        if (prgnam .eq. 'lmxbs') then
          call info0(0,0,0,
     *    '%N%1f-shift=x1,x1,x3 -spec '//
     *    '-dup=d1,d2,d3[,expr] -bs=val -ss=val')
        endif

        if (prgnam .eq. 'lmctl') then
          call info0(0,0,0,
     *    '%N%1f -spin1 -spinf -mad -enu')
        endif

        if (prgnam .eq. 'lmpg') then
          call info0(0,0,0,
     *    '%N%1f -map -onesp')
        endif

        if (prgnam .eq. 'lmplan') then
          call info0(0,0,0,
     *    '%N%1f---')
        endif

        if (prgnam .eq. 'lmshf') then
          call info0(0,0,0,
     *    '%N%1f-enu=expr   linearize pot. pars around ''expr''')
        endif

        if (prgnam .eq. 'lmimp') then
          call info0(0,0,0,
     *    '%N%1f -rs -4 -5 -3s -4s -47u -5s (1 is required)')
        endif

        if (prgnam .eq. 'lmcor') then
          call info0(0,0,0,
     *    '%N%1f -findr')
        endif

        if (prgnam .eq. 'lmavgm') then
          call info0(0,0,0,
     *    '%N%1f -spin1')
        endif

        if (prgnam .eq. 'lmmix') then
          call info0(0,0,0,
     *    '%N%1f-fn=mix-file-name -bin2a or -a2bin')
        endif

        if (prgnam .eq. 'tbe') then
          call info0(0,0,0,
     *    '%N%1f--band --wpos=fnam -cont -dumph'//
     *    ' --st --md=# --mv=# --xyz=#'
     *    //'%N%3fThe last 4 switches apply to'//
     *    ' molecular dynamics simulations')
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
     *vn-vn2/1d5,wksize/1000)
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

