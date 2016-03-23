      program lmfa
      use m_struc_def
      use m_rdctrl2_func, only: rdctrl2
      use m_globalvariables
      implicit none
      integer procid, master, mpipid, nproc
      integer,parameter:: lstrn=10000
      character(lstrn):: sstrnmix,sstrnsymg='' !,jobid !header,
      character prgnam*8, vrsion(2)*6, ext*100
      logical lgors
      integer mxspec
      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl
      type(s_lat):: v_slat
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_mix):: v_smix
      type(s_array):: v_sarry
      type(s_move):: v_smove
      type(s_str):: v_sstr
      type(s_spec),allocatable:: v_sspec(:)
      type(s_site),allocatable:: v_ssite(:)
      parameter (mxspec=256)
      character*8 slabl(mxspec)
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     *igets,a2vec,NULLI
      double precision dglob !
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)
      data prgnam /'LMFA'/
      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=511)
      character*8 alabl
      character:: recrd_(0:mxrecs*recln0-1)
      character*(mxrecs*recln0):: recrd
      equivalence ( recrd, recrd_ )
      real(8),parameter::    NULLR =-99999
#if MPI|MPIK
      include 'mpif.h'
      integer ierr
#endif
      call show_programinfo(6)
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
      if (swtmp .or. cmdopt('--h',3,0,outs)) call lmhelp(prgnam)!,vn(2),vn2(2),0)
      call headl2(prgnam,0,stdo)
      call finits(2,0,0,i)
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
C --- Set the top-level verbosity if spec'd from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif
C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
      if (procid .eq. master) then
        if (.not.cmdopt('--input',7,0,strn)) then
          nfilin = fopna('CTRL',-1,1)
          call findctrlstart(nfilin)
          alabl = '#{}% ct '
          if (cmdopt('--show',6,0,strn)) alabl = '#{}% ctp'
          call rdfile(nfilin,alabl,recrd,mxrecs,strn,recln0,nrecs)
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
C --- Read recrd ---
      call rdctrl2(recrd,recln0,nrecs,prgnam,F, !,vrsion,vn,vn2
     *slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     *v_sarry,v_smove,sstrnmix,sstrnsymg) !,v_stb
      print *,'end of rdctrl2 in imfav7'
      call setcg(v_slat,8,12)
      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
C --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      lc = 20
      if ( .not. iand(1,int(v_sctrl%lqp)) .ne.0 ) lc = lc+2
      print *,'goto mksym:',trim(strn)
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
      print *,'end of mksym x'
C ... Set default values for species data
      print *,'goto defspc'
      call defspc(v_sspec)
      print *,'end of defspc'
      if ( int(v_sctrl%quit) .eq. 1 ) then
        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif
      print *,'goto freeat'
      call freeat(v_sctrl,v_sspec,v_sham,v_spot)
      call rx0(prgnam)
      end

      subroutine lmhelp(prgnam)!,vn,vn2,wksize)
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
      character prgnam*8
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
        if (prgnam .eq. 'mmag') then
          call info0(0,0,0,'%N%1f--cont --wrsj[:fn=name]')
        endif

      endif
      call fexit(0,0,' ',0)
      end subroutine lmhelp


      subroutine findctrlstart(nfilin)
      character(len=9):: strn, ccc
      do
        read(nfilin,"(a)",err=1010,end=1010) strn
        if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
      end subroutine findctrlstart

