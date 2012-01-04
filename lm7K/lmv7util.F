C For lmdos and lmchk. Single core mode only. No MPI modes.
      program lmfutil
C- Main program for lmf and so on.
      use m_struc_def
      use m_rdctrl2_func, only: rdctrl2
      use m_globalvariables
      implicit none
c      integer procid, master, mpipid, nproc

cC Heap allocation
c      integer wksize
c      parameter(wksize= 10 000 000)
c      integer w(wksize)
cC     Next two lines guarantee w is aligned along a d.p. boundary
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

ckino Jan.04.2012:        type(s_spec),pointer:: v_sspec(:) =>NULL()
      type(s_spec),allocatable:: v_sspec(:) 
ckino Jan.04.2012:        type(s_site),pointer:: v_ssite(:) =>NULL()
      type(s_site),allocatable:: v_ssite(:) 

      parameter (mxspec=256)
      character*8 slabl(mxspec)

C ... miscellaneous local variables
      character strn*1000,outs*20
      integer i,j,k,fxst,fadd,fopn,lgunit,i1mach,auxmod,lc,stdo,fextg,
     .igets,a2vec,NULLI
      double precision vn(2),vn2(2),dglob
      logical T,F,swtmp,cmdopt,ltet
      parameter (T=.true., F=.false., NULLI=-99999)

C ... Program-dependent name and help
#if  LMDOS
      data prgnam /'LMDOS'/
#endif
#if LMCHK
      data prgnam /'LMCHK'/ auxmod /1/
#endif

      integer:: nfilin,mxrecs,recln0,nrecs,fopna
      parameter (mxrecs=10000,recln0=120)
      character*8 alabl

c      character:: recrd(0:mxrecs*recln0-1)
c      character*(1000):: recrd
c      equivalence ( recrd, recrd_ )
      character*(mxrecs*recln0):: recrd

      real(8),parameter::    NULLR =-99999
      real(8)::pwmode,pwemin,pwemax

      integer::ic
      real(8):: bz_w,dval!test

C --- Version ---
      call show_programinfo(6)

      vrsion(1) = 'LM'
      vrsion(2) = ' '
      vn(1) = 7.00d0
      vn2(1) = 0
      vn(2) = 0.00d0
      vn2(2) = 0
      vn(1) = vn(1) + vn2(1)/1d5
      vn(2) = vn(2) + vn2(2)/1d5
      stdo = lgunit(1)
Changenglob      i = dglob('stdo',dble(stdo),1)
      globalvariables%stdo = dble(stdo); globalvariables%l_stdo = globalvariables%l_stdo +1; i = dble(stdo)
c      master = 0
c      procid = mpipid(1)
c      nproc  = mpipid(0)

C --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs))
c     .call lmhelp(prgnam,vn(2),vn2(2),wksize)
     .call lmhelp(prgnam,vn(2),vn2(2),0)
      if (cmdopt('--version',9,0,outs)) then
        write(*,'(f4.2)') vn(2)
        call cexit(0,1)
      endif

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
c      if (cmdopt('--input',6,0,strn)) then
c        if (nproc .gt. 0) call rx('--input not allowed with MPI')
c      else
      if (fxst('CTRL') .ne. 1) then
        call awrit0(' '//prgnam//'%a:%9pmissing ctrl file',' ',80,
     .  i1mach(2))
        swtmp = .true.
      endif
c      endif
c      call mpibc1(swtmp,1,1,.false.,'lmf','error')
c      if (swtmp) call cexit(-1,1)

C ... Set special file directory for temporary files
C     User may which to customize the directory
C     Default is to use the standard directory
C     sttmpd is located at the bottom of this file.
c      call sttmpd

C ... File logical units
      i = fadd('TMP',-1,4)
      i = fadd('BAND',-1,4)
c#if TBE
c      i = fadd('STRT',-1,4)
c      i = fadd('QMOM',-1,4)
c#endif
c#if LM | LMDOS | LMGF
#if LMDOS
      i = fadd('MOMS',-1,4)
      i = fadd('CLS',-1,4)
#endif

C --- Set the top-level verbosity if specified from cmd line ---
      if (cmdopt('--pr',4,0,outs)) then
        i = 4
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,j,k)
        if (i .eq. 1) call setpr(k)
      endif

C --- Input from ctrl file ---
C     recrd, nrecs are obtained.
      nrecs  = 0
Cki      allocate( recrd( 0:mxrecs*recln0-1 ) )
c      if (procid .eq. master) then
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
     .    'of input file ----------------------')
        endif
        call info2(i,0,1,' '//prgnam//'%a : %i lines read from'//
     .  ' input file',nrecs,0)
        if (cmdopt('--showp',7,0,strn)) call cexit(0,1)
      endif
c      endif

C      write(stdo,*)' === TEST out recrd === '
C      j = 0
C      do i = 1, nrecs
C        write(stdo,"(i2,1x,256a1)") procid,recrd(j:j+recln0-1)
C        j = j + recln0
C      enddo
C      write(stdo,*)' --- end of rdfile --- '
C      call rx('done')

C --- Read recrd ---
c tkotani found nmap is not zero cleared.
cccccccccccccccccc
c      print *,' uuuu000 111 nmap=',int(v_sctrl%nmap) 
c      v_sctrl%nmap=9999
cccccccccccccccccc

      call rdctrl2(recrd,recln0,nrecs,prgnam,vrsion,vn,vn2,F,
ckino Jan.04.2012:       .slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
ckino Jan.04.2012:       .slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     .slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     .v_sarry,v_smove,v_stb,sstrn)


cccccccccccccccccc
c      print *,' uuuu000 222 nmap=',int(v_sctrl%nmap) 
cccccccccccccccccc

C --- Lattice setup ---
      call setcg(v_slat,8,12)
ckino Jan.04.2012:        call lattic(v_slat,v_sctrl,v_ssite,v_sarry)
      call lattic(v_slat,v_sctrl,v_ssite,v_sarry)

C --- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      call upacks('strn symg',i,j)
      if (j .ge. i) strn = sstrn(i:j+1)
      if (cmdopt('--nosym',7,0,outs)) strn = ' '
      lc = 20
c#ifndef LMFGWD
Cgetarg       if (.not. lgors('ctrl lqp,1',v_sctrl)) lc = lc+2
c this is equivalent to
      if ( v_sctrl%lqp==0 ) lc = lc+2 
c#endif
ckino Jan.04.2012:        call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
c nmap is not zero cleared ---> takao think it is zero for all cases here.

C --- Allocate permanent class arrays, maps and other initialization ---
ckino Jan.04.2012:        call clsprm(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
      call clsprm(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .v_slat,v_sbz,v_sstr,sstrn)

C --- Read available class parameters from file ---
#if LMCHK
      call aiocls(.false.,0,v_sarry,v_sctrl,v_sham,v_spot,
ckino Jan.04.2012:       .v_sspec,v_slat,1,0)
     .v_sspec,v_slat,1,0)
      call aiocls(.false.,17,v_sarry,v_sctrl,v_sham,v_spot,
ckino Jan.04.2012:       .v_sspec,v_slat,1,0)
     .v_sspec,v_slat,1,0)
      if (lc .ge. 20)
ckino Jan.04.2012:       .call clsprp(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .call clsprp(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .v_slat,v_sbz,v_sstr,sstrn)
#endif

      if (cmdopt('--sfill',7,0,strn)) then
        call rx('no --sfill option in lm7K')
      endif

C ... Set default values for species data
ckino Jan.04.2012:        call defspc(v_sspec)
      call defspc(v_sspec)

C ... quit after SHOW
Cgetarg       if (igets('ctrl quit',v_sctrl) .eq. 1) then
      if ( int(v_sctrl%quit) .eq. 1 ) then 

        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif

#if LMCHK & FP
      if (cmdopt('--fp',4,0,strn)) then
ckino Jan.04.2012:          call fpchk(v_sspec,v_ssite)
ckino Jan.04.2012:          call fpchk(v_sspec,v_ssite)
        call fpchk(v_sspec,v_ssite)
        call cexit(0,1)
      endif
#endif

#if LMCHK
c      if (cmdopt('--findes',8,0,strn)) auxmod = 128
ckino Jan.04.2012:        call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
      call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
ckino Jan.04.2012:       .v_ssite,v_sstr,v_sarry,v_sbz,slabl,sstrn,auxmod,0) !wksize)
     .v_ssite,v_sstr,v_sarry,v_sbz,slabl,sstrn,auxmod,0) !wksize)
#endif

#if LMDOS
      call asados(prgnam,v_sbz,v_sctrl,v_sham,v_slat,
ckino Jan.04.2012:       .v_sspec,v_ssite,v_sarry)
ckino Jan.04.2012:       .v_sspec,v_ssite,v_sarry)
     .v_sspec,v_ssite,v_sarry)
#endif

C -------------- End of program -------------
 1000 continue
      call rx0(prgnam)
      end program lmfutil


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

        if (prgnam .eq. 'lmdos') then
          call info0(0,0,0,
     .    '%N%1f--dos:options  modifies number and kinds of dos '//
     .    'generated;%N%16fsee documentation (doc/lmto.html)')
        endif

        if (prgnam .eq. 'lmchk') then
          call info0(0,0,0,
     .    '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
     .    '[:pairs:pair-list]...'//
     .    '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
     .    '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
     .    '%N --findes'//
     .    '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
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


C$$$      subroutine sttmpd
C$$$C- Creates special tmpdir for saving temporary files
C$$$C  User may wish to customize this routine.
C$$$      character tmpdir*100
C$$$      integer fopnT,ifi
C$$$C ... for henry, lm-MPIK
C$$$C      integer i1,i2,nw
C$$$C      character*40 strn

C$$$C     return
C$$$C ... Set customization of temporary directory here, if desired
C$$$C     This is usual default (current working directory)
C$$$      tmpdir = ' '
C$$$C     call gtenv('HOME',tmpdir)
C$$$C     call gtenv('TMPDIR',tmpdir)
C$$$C     call getenv('HOME',tmpdir)
C$$$C     call getenv('TMPDIR',tmpdir)
C$$$C     tmpdir = '/home/tmp/'

C$$$C ... for henry, lm-MPIK specific
C$$$C      call getenv('HOME',strn)
C$$$C      call strip(strn,i1,i2)
C$$$C      call wrdsg(strn(i1:i2),0,'/',nw)
C$$$C      call wordg(strn,0,'/',nw,i1,i2)
C$$$C      strn = strn(i1:)
C$$$C      if (strn .eq. 'markv') strn = 'svan2'
C$$$C      call word(strn,1,i1,i2)
C$$$C      tmpdir = '/home/' // strn(i1:i2) // '/tmp'

C$$$C ... Set the directory
C$$$      ifi = fopnT(tmpdir,0,0,11)
C$$$C     debugging check
C$$$C      ifi = fopnT('tmp' ,-1,0,0)
C$$$C      call fshow
C$$$C      print *, ifi
C$$$C      write(ifi,*) 'test'
C$$$C      call rx('done')
C$$$Cgetarg       end
C$$$       end subroutine sttmpd 



      subroutine findctrlstart(nfilin)
      use m_globalvariables
      character(len=9):: strn, ccc
      do
        read(nfilin,"(a)",err=1010,end=1010) strn
        if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
Cgetarg       end
      end subroutine findctrlstart 

Cgetarg...info...           structure ['clsprm', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['clsprm', 'ham', 'sham']
Cgetarg...info...           structure ['clsprm', 'array', 'sarray']
Cgetarg...info...           structure ['clsprm', 'pot', 'spot']
Cgetarg...info...           structure ['clsprm', 'lat', 'slat']
Cgetarg...info...           structure ['clsprm', 'spec', 'sspec']
Cgetarg...info...           structure ['clsprm', 'bz', 'sbz']
Cgetarg...info...           use_to_add ['clsprm', 'm_struc_def']
      subroutine clsprm(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .sstr,sstrn)
      
      use m_struc_def  !Cgetarg

C- Allocate memory for, and try to read, class-related data
C ----------------------------------------------------------------------
Ci    mode :1s digit
Ci         : 1 allocates ASA arrays (see spot,sham,sarray)
Ci         : 2 allocates TBE arrays:  oifrlx
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nclass nl nspin nspec loptc nmap
Ci     Stored:    *
Ci     Passed to: lgors ioqpp aiocls clsprp
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oics oipc oclabl opos
Ci     Stored:    ogroup ormax ohave oifrlx ozos
Ci     Passed to: aiocls clsprp
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: *
Ci     Stored:    osop oivso
Ci     Passed to: aiocls clsprp
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: ovrmax orhrmx opnu oqnu opp oves
Ci     Stored:    opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp orhos
Ci                osop ogrrme ovintr opmpol oqpp
Ci     Passed to: ioqpp aiocls clsprp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rmt idmod z lmxa
Ci     Stored:    *
Ci     Passed to: aiocls clsprp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat
Ci     Stored:    *
Ci     Passed to: aiocls clsprp
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: ndos
Ci     Stored:    *
Ci     Passed to: clsprp
Ci
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci     Elts read: n map
Ci     Stored:    *
Ci     Passed to: *
Ci   sstrn :struct for global strings
Ci     Elts read: map
Ci     Stored:    *
Ci     Passed to: clsprp
Cu Updates
Cu   08 Nov 07 (J. Xu) LDA+U implementation; qpp is complex
Cu   30 Sep 04 Reads/writes relativistic ppar's
Cu   02 Apr 04 SO array redimensioned; make and pack pot->nrhos
Cu   18 Mar 03 Allocate space for relativistic potential parameters
Cu   07 Feb 03 SO array redimensioned
Cu   19 Feb 02 entry clsprp copies class data to doubly padded layers
Cu   28 Apr 98 allocate space for radial matrix elements of grad
C ----------------------------------------------------------------------
C     implicit none
      integer mode
      character*(*) sstrn
Cgetarg       double precision sctrl(1),sarray(1),sham(1),spot(1),sspec(1),
Cgetarg      .  slat(1),sbz(1),sstr(1)
      type(s_ctrl)::sctrl
      type(s_array)::sarray
      type(s_ham)::sham
      type(s_pot)::spot
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_bz)::sbz
      type(s_str)::sstr

C Local variables
      logical lgors,cmdopt
      character*80 strn,clabl*8
                           integer:: ogroup , oifrlx , oivso , oipc , opos , ozos , n0 
ckino Dec.27.2011:                             real(8),pointer :: rv_p_ormax(:) =>NULL()

ckino Dec.27.2011:                            real(8),pointer :: rv_p_osop(:) =>NULL()

ckino Dec.27.2011:                           real(8),pointer :: rv_p_opprel(:) =>NULL()

ckino Dec.27.2011:                          real(8),pointer :: rv_p_orhrmx(:) =>NULL()

ckino Dec.27.2011:                         real(8),pointer :: rv_p_orhos(:) =>NULL()

ckino Dec.27.2011:                        real(8),pointer :: rv_p_oqt(:) =>NULL()

ckino Dec.27.2011:                       real(8),pointer :: rv_p_oqpp(:) =>NULL()

ckino Dec.27.2011:                      real(8),pointer :: rv_p_oqnu(:) =>NULL()

ckino Dec.27.2011:                     real(8),pointer :: rv_p_oqc(:) =>NULL()

ckino Dec.27.2011:                    real(8),pointer :: rv_p_opp(:) =>NULL()

ckino Dec.27.2011:                   real(8),pointer :: rv_p_opnu(:) =>NULL()

ckino Dec.27.2011:                  real(8),pointer :: rv_p_opmpol(:) =>NULL()

ckino Dec.27.2011:                 real(8),pointer :: rv_p_ogrrme(:) =>NULL()

ckino Dec.27.2011:                real(8),pointer :: rv_p_obxc(:) =>NULL()

ckino Dec.27.2011:               integer,pointer :: iv_p_oinitc(:) =>NULL()

ckino Dec.27.2011:              real(8),pointer :: rv_p_ovrmax(:) =>NULL()

ckino Dec.27.2011:             real(8),pointer :: rv_p_ovintr(:) =>NULL()

ckino Dec.27.2011:            real(8),pointer :: rv_p_oves(:) =>NULL()

ckino Dec.27.2011:           real(8),pointer :: rv_p_ovdif(:) =>NULL()

ckino Dec.27.2011:          integer,pointer :: iv_p_oics(:) =>NULL()

ckino Dec.27.2011:         real(8),pointer :: rv_p_oclabl(:) =>NULL()

      parameter (n0=10)
      integer nbas,nclass,nspec,nl,nsp,nlspc,igets,nclasp,i1,i2,nmap,ic,
     .is,ib,iclbsj,idmod(n0),iprint,lmx,nclspp,nqpp,lnsph,isw,neula,
     .nbf,nrhos,nlibu
C     integer nmto,oppn
      double precision plat(9),z,xx,dval
C ... Heap

      integer,allocatable:: ics(:)
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg       call upack('ctrl nbas nclass nl nspin nspec',sctrl,nbas,nclass,nl,
Cgetarg      .  nsp,nspec)

      nbas=sctrl%nbas
      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec

Cgetarg...info...           integer :: nlibu , integer(8) :: sham%nlibu 1
Cgetarg       call upack('ham nlibu',sham,nlibu,0,0,0,0)
      
      nlibu=sham%nlibu

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg       call upack1('array nclasp',sarray,nclasp)
      
      nclasp=sarray%nclasp

      nclspp = max(2*nclasp-nclass,nspec)
      nlspc = nl*nsp*max(nclspp,nspec)

C ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
ckino Dec.27.2011:           allocate(spot%rv_p_opnu(abs(-nlspc)))
ckino 12月.27.2011:           allocate(spot%rv_p_opnu(abs(-nlspc)))
         allocate(spot%rv_a_opnu(abs(-nlspc)))
ckino Dec.27.2011:           if (-nlspc<0) spot%rv_p_opnu(:)=0.0d0
ckino 12月.27.2011:           if (-nlspc<0) spot%rv_p_opnu(:)=0.0d0
         if (-nlspc<0) spot%rv_a_opnu(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_oqnu(abs(-3*nlspc)))
ckino 12月.27.2011:           allocate(spot%rv_p_oqnu(abs(-3*nlspc)))
         allocate(spot%rv_a_oqnu(abs(-3*nlspc)))
ckino Dec.27.2011:           if (-3*nlspc<0) spot%rv_p_oqnu(:)=0.0d0
ckino 12月.27.2011:           if (-3*nlspc<0) spot%rv_p_oqnu(:)=0.0d0
         if (-3*nlspc<0) spot%rv_a_oqnu(:)=0.0d0

ckino Dec.27.2011:           allocate(rv_p_opp(abs(-6*nlspc)))
         allocate(spot%rv_a_opp(abs(-6*nlspc)))
ckino Dec.27.2011:           if (-6*nlspc<0) rv_p_opp(:)=0.0d0
         if (-6*nlspc<0) spot%rv_a_opp(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_opprel(abs(-32*nl*nl*max(nclspp,nspec))))
ckino 12月.27.2011:           allocate(spot%rv_p_opprel(abs(-32*nl*nl*max(nclspp,nspec))))
         allocate(spot%rv_a_opprel(abs(-32*nl*nl*max(nclspp,nspec))))
ckino Dec.27.2011:           if (-32*nl*nl*max(nclspp,nspec)<0) spot%rv_p_opprel(:)=0.0d0
ckino 12月.27.2011:           if (-32*nl*nl*max(nclspp,nspec)<0) spot%rv_p_opprel(:)=0.0d0
         if (-32*nl*nl*max(nclspp,nspec)<0) spot%rv_a_opprel(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_oqc(abs(-nclasp)))
ckino Dec.28.2011:           allocate(spot%rv_p_oqc(abs(-nclasp)))
ckino Dec.27.2011:           if (-nclasp<0) spot%rv_p_oqc(:)=0.0d0
ckino Dec.28.2011:           if (-nclasp<0) spot%rv_p_oqc(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_oqt(abs(-nclasp)))
ckino 12月.27.2011:           allocate(spot%rv_p_oqt(abs(-nclasp)))
         allocate(spot%rv_a_oqt(abs(-nclasp)))
ckino Dec.27.2011:           if (-nclasp<0) spot%rv_p_oqt(:)=0.0d0
ckino 12月.27.2011:           if (-nclasp<0) spot%rv_p_oqt(:)=0.0d0
         if (-nclasp<0) spot%rv_a_oqt(:)=0.0d0

ckino Dec.27.2011:           allocate(iv_p_oinitc(abs(-nclasp)))
ckino Dec.27.2011:           allocate(sarray%iv_p_ohave(abs(-nclasp)))
         allocate(sarray%iv_a_oipq(abs(-nclasp)))
ckino Dec.27.2011:           if (-nclasp<0) iv_p_oinitc(:)=0
ckino Dec.27.2011:           if (-nclasp<0) sarray%iv_p_ohave(:)=0
         if (-nclasp<0) sarray%iv_a_oipq(:)=0

c        call defi (ogroup,-2*nclasp)
ckino Dec.27.2011:           allocate(spot%rv_p_oves(abs(-nclspp)))
ckino 12月.27.2011:           allocate(spot%rv_p_oves(abs(-nclspp)))
         allocate(spot%rv_a_oves(abs(-nclspp)))
ckino Dec.27.2011:           if (-nclspp<0) spot%rv_p_oves(:)=0.0d0
ckino 12月.27.2011:           if (-nclspp<0) spot%rv_p_oves(:)=0.0d0
         if (-nclspp<0) spot%rv_a_oves(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_ovdif(abs(-nclspp)))
ckino 12月.27.2011:           allocate(spot%rv_p_ovdif(abs(-nclspp)))
         allocate(spot%rv_a_ovdif(abs(-nclspp)))
ckino Dec.27.2011:           if (-nclspp<0) spot%rv_p_ovdif(:)=0.0d0
ckino 12月.27.2011:           if (-nclspp<0) spot%rv_p_ovdif(:)=0.0d0
         if (-nclspp<0) spot%rv_a_ovdif(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_orhrmx(abs(-nclspp)))
ckino 12月.27.2011:           allocate(spot%rv_p_orhrmx(abs(-nclspp)))
         allocate(spot%rv_a_orhrmx(abs(-nclspp)))
ckino Dec.27.2011:           if (-nclspp<0) spot%rv_p_orhrmx(:)=0.0d0
ckino 12月.27.2011:           if (-nclspp<0) spot%rv_p_orhrmx(:)=0.0d0
         if (-nclspp<0) spot%rv_a_orhrmx(:)=0.0d0

ckino Dec.27.2011:           allocate(spot%rv_p_ovrmax(abs(-2*nclspp)))
ckino Dec.28.2011:           allocate(spot%rv_p_ovrmax(abs(-2*nclspp)))
         allocate(spot%rv_a_ovrmax(abs(-2*nclspp)))
ckino Dec.27.2011:           if (-2*nclspp<0) spot%rv_p_ovrmax(:)=0.0d0
ckino Dec.28.2011:           if (-2*nclspp<0) spot%rv_p_ovrmax(:)=0.0d0
         if (-2*nclspp<0) spot%rv_a_ovrmax(:)=0.0d0

C       call defdr(oqnu,-3*nlspc)
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg         call upack1('array oics',sarray,oics)
        
ckino Dec.27.2011:           iv_p_oics => sarray%iv_p_oics 


c takao sp2cls is replaced.
c        call sp2cls('spec rmt',sspec,w(oics),1,1,nclspp,ormax)
        allocate( ics(nclspp) )
ckino Dec.27.2011:           call icopy ( nclspp , sarray%iv_p_oics , 1 , ics , 1 ) 
ckino Dec.27.2011:           call icopy ( nclspp , sarray%iv_p_oics , 1 , ics , 1 ) 
         call icopy ( nclspp , sarray%iv_a_oics , 1 , ics , 1 ) 

        print *,' qqqqq ics=',ics
ckino Dec.27.2011:           allocate(sarray%rv_p_ormax(abs(nclspp)))
ckino 12月.27.2011:           allocate(sarray%rv_p_ormax(abs(nclspp)))
         allocate(sarray%rv_a_ormax(abs(nclspp)))
ckino Dec.27.2011:           if (nclspp<0) sarray%rv_p_ormax(:)=0.0d0
ckino 12月.27.2011:           if (nclspp<0) sarray%rv_p_ormax(:)=0.0d0
         if (nclspp<0) sarray%rv_a_ormax(:)=0.0d0

ckino Dec.27.2011:           call dcopy ( nclspp , sspec ( ics ( 1:nclspp ) ) %rmt , 1 , sarray%rv_p_ormax 
ckino 12月.27.2011:           call dcopy ( nclspp , sspec ( ics ( 1:nclspp ) ) %rmt , 1 , sarray%rv_p_ormax 
         call dcopy ( nclspp , sspec ( ics ( 1:nclspp ) ) %rmt , 1 , sarray%rv_a_ormax 
     .   , 1 ) 


Cgetarg...info...           integer :: ogroup , integer(8) :: sarray%ogroup 1
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg...info...           integer :: oinitc , integer(8) :: sarray%ohave 1
Cgetarg         call pack5('array ogroup ormax ohave',sarray,ogroup,ormax,
Cgetarg      .    oinitc,0,0)
        
c        sarray%ogroup=ogroup 
ckino Dec.27.2011:           sarray%rv_p_ormax => rv_p_ormax 

ckino Dec.27.2011:           sarray%iv_p_ohave => iv_p_oinitc 


Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: oqc , integer(8) :: spot%oqc 1
Cgetarg...info...           integer :: oqt , integer(8) :: spot%oqt 1
Cgetarg         call pack5('pot opnu oqnu oqc oqt',spot,opnu,oqnu,oqc,oqt,0)
        
ckino Dec.27.2011:           spot%rv_p_opnu => rv_p_opnu 

ckino Dec.27.2011:           spot%rv_p_oqnu => rv_p_oqnu 

ckino Dec.27.2011:           spot%rv_p_oqc => rv_p_oqc 

ckino Dec.27.2011:           spot%rv_p_oqt => rv_p_oqt 


Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg         call pack5('pot oves ovdif ovrmax orhrmx opp',spot,oves,ovdif,
Cgetarg      .    ovrmax,orhrmx,opp)
        
ckino Dec.27.2011:           spot%rv_p_oves => rv_p_oves 

ckino Dec.27.2011:           spot%rv_p_ovdif => rv_p_ovdif 

ckino Dec.27.2011:           spot%rv_p_ovrmax => rv_p_ovrmax 

ckino Dec.27.2011:           spot%rv_p_orhrmx => rv_p_orhrmx 

ckino Dec.27.2011:           spot%rv_p_opp => rv_p_opp 


C       Even in collinear case, bxc might be written to
ckino Dec.27.2011:           allocate(spot%rv_p_obxc(abs(-3*nclasp)))
ckino Dec.28.2011:           allocate(spot%rv_p_obxc(abs(-3*nclasp)))
         allocate(spot%rv_a_obxc(abs(-3*nclasp)))
ckino Dec.27.2011:           if (-3*nclasp<0) spot%rv_p_obxc(:)=0.0d0
ckino Dec.28.2011:           if (-3*nclasp<0) spot%rv_p_obxc(:)=0.0d0
         if (-3*nclasp<0) spot%rv_a_obxc(:)=0.0d0

Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg         call pack5('pot opprel obxc nrhos',spot,opprel,obxc,0,0,0)
        
ckino Dec.27.2011:           spot%rv_p_opprel => rv_p_opprel 

ckino Dec.27.2011:           spot%rv_p_obxc => rv_p_obxc 

        spot%nrhos=0 

Cgetarg         if (lgors('ctrl lncol,-1',sctrl)) then
c
c takao Kino's bug--> his replied that iand(-1,i)=i 
c         if ( iand(-1,int(sctrl%lncol)) .ne.0 ) then 

        if ( sctrl%lncol /=0 ) then 
Cgetarg...info...           integer :: neula , integer(8) :: sham%neula 1
Cgetarg...info...           integer :: nbf , integer(8) :: sham%nbf 1
Cgetarg           call upack2('ham neula nbf',sham,neula,nbf)
          neula=sham%neula
          nbf=sham%nbf
          nrhos = max(nl,neula,nbf)
ckino Dec.27.2011:             if (associated(spot%rv_p_obxc)) deallocate(spot%rv_p_obxc)
ckino Dec.28.2011:             if (associated(spot%rv_p_obxc)) deallocate(spot%rv_p_obxc)
           if (allocated(spot%rv_a_obxc)) deallocate(spot%rv_a_obxc)

ckino Dec.27.2011:             allocate(spot%rv_p_orhos(abs(2*3*4*nrhos*nclasp)))
ckino Dec.28.2011:             allocate(spot%rv_p_orhos(abs(2*3*4*nrhos*nclasp)))
ckino Dec.27.2011:             if (2*3*4*nrhos*nclasp<0) spot%rv_p_orhos(:)=0.0d0
ckino Dec.28.2011:             if (2*3*4*nrhos*nclasp<0) spot%rv_p_orhos(:)=0.0d0

ckino Dec.27.2011:             allocate(spot%rv_p_osop(abs(-nlspc*nsp*9)))
ckino 12月.27.2011:             allocate(spot%rv_p_osop(abs(-nlspc*nsp*9)))
           allocate(spot%rv_a_osop(abs(-nlspc*nsp*9)))
ckino Dec.27.2011:             if (-nlspc*nsp*9<0) spot%rv_p_osop(:)=0.0d0
ckino 12月.27.2011:             if (-nlspc*nsp*9<0) spot%rv_p_osop(:)=0.0d0
           if (-nlspc*nsp*9<0) spot%rv_a_osop(:)=0.0d0

ckino Dec.27.2011:             allocate(spot%rv_p_obxc(abs(-3*nclasp)))
ckino Dec.28.2011:             allocate(spot%rv_p_obxc(abs(-3*nclasp)))
           allocate(spot%rv_a_obxc(abs(-3*nclasp)))
ckino Dec.27.2011:             if (-3*nclasp<0) spot%rv_p_obxc(:)=0.0d0
ckino Dec.28.2011:             if (-3*nclasp<0) spot%rv_p_obxc(:)=0.0d0
           if (-3*nclasp<0) spot%rv_a_obxc(:)=0.0d0

Cgetarg...info...           integer :: nrhos , integer(8) :: spot%nrhos 1
Cgetarg...info...           integer :: orhos , integer(8) :: spot%orhos 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg           call pack5('pot nrhos orhos obxc osop',spot,
Cgetarg      .      nrhos,orhos,obxc,osop,0)
          spot%nrhos=nrhos 
ckino Dec.27.2011:             spot%rv_p_orhos => rv_p_orhos 

ckino Dec.27.2011:             spot%rv_p_obxc => rv_p_obxc 

ckino Dec.27.2011:             spot%rv_p_osop => rv_p_osop 

        endif

cki        ogrrme = 0
ckino Dec.27.2011:            nullify(spot%rv_p_ogrrme)
ckino Dec.28.2011:            nullify(spot%rv_p_ogrrme)
        if (allocated(spot%rv_a_ogrrme)) deallocate(spot%rv_a_ogrrme)
Cgetarg         if (igets('ctrl loptc',sctrl) .gt. 0)
Cgetarg      .    call defrr(ogrrme,-16*nl*nsp*nclasp)
         if ( int ( sctrl%loptc ) .gt. 0 ) then 
ckino Dec.27.2011:           allocate(spot%rv_p_ogrrme(abs(-16*nl*nsp*nclasp)))
ckino Dec.28.2011:           allocate(spot%rv_p_ogrrme(abs(-16*nl*nsp*nclasp)))
         allocate(spot%rv_a_ogrrme(abs(-16*nl*nsp*nclasp)))
ckino Dec.27.2011:           if (-16*nl*nsp*nclasp<0) spot%rv_p_ogrrme(:)=0.0d0
ckino Dec.28.2011:           if (-16*nl*nsp*nclasp<0) spot%rv_p_ogrrme(:)=0.0d0
         if (-16*nl*nsp*nclasp<0) spot%rv_a_ogrrme(:)=0.0d0
         endif 


Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg         call pack1('pot ogrrme',spot,ogrrme)
        
ckino Dec.27.2011:           spot%rv_p_ogrrme => rv_p_ogrrme 


Cgetarg         if (lgors('ctrl lsx,3',sctrl) .or.
Cgetarg      .      lgors('ctrl lscr,-1',sctrl)) then
        if ( iand(3,int(sctrl%lsx)) .ne.0 .or. iand(-1,int(sctrl%lscr)) .ne.0 
     .  ) then 

ckino Dec.27.2011:             allocate(spot%rv_p_ovintr(abs(-nclasp*(nl*nsp)**2)))
ckino Dec.28.2011:             allocate(spot%rv_p_ovintr(abs(-nclasp*(nl*nsp)**2)))
           allocate(spot%rv_a_ovintr(abs(-nclasp*(nl*nsp)**2)))
ckino Dec.27.2011:             if (-nclasp*(nl*nsp)**2<0) spot%rv_p_ovintr(:)=0.0d0
ckino Dec.28.2011:             if (-nclasp*(nl*nsp)**2<0) spot%rv_p_ovintr(:)=0.0d0
           if (-nclasp*(nl*nsp)**2<0) spot%rv_a_ovintr(:)=0.0d0

Cgetarg...info...           integer :: ovintr , integer(8) :: spot%ovintr 1
Cgetarg           call pack1('pot ovintr',spot,ovintr)
          
ckino Dec.27.2011:             spot%rv_p_ovintr => rv_p_ovintr 


        endif
Cgetarg         lnsph = isw(lgors('ctrl lasa,32',sctrl))
c         lnsph = isw ( iand(32,int(sctrl%lasa)) .ne.0 ) 
ctakao remove asa
        lnsph = isw ( .false.)

C   ... Allocate memory for, attempt to read qpp from file
        if (lnsph .ne. 0 .or. nlibu .ne. 0) then
ckino Dec.27.2011:             allocate(spot%rv_p_opmpol(abs(-(2*nl-1)*nl**2*3*nsp*nclasp)))
ckino 12月.27.2011:             allocate(spot%rv_p_opmpol(abs(-(2*nl-1)*nl**2*3*nsp*nclasp)))
           allocate(spot%rv_a_opmpol(abs(-(2*nl-1)*nl**2*3*nsp*nclasp)))
ckino Dec.27.2011:             if (-(2*nl-1)*nl**2*3*nsp*nclasp<0) spot%rv_p_opmpol(:)=0.0d0
ckino 12月.27.2011:             if (-(2*nl-1)*nl**2*3*nsp*nclasp<0) spot%rv_p_opmpol(:)=0.0d0
           if (-(2*nl-1)*nl**2*3*nsp*nclasp<0) spot%rv_a_opmpol(:)=0.0d0

          i1 = nl**2
          nqpp = (i1*(i1+1))/2
ckino Dec.27.2011:             allocate(spot%rv_p_oqpp(abs(2*nqpp*4*nsp*nbas)))
           allocate(spot%rv_a_oqpp(abs(2*nqpp*4*nsp*nbas)))
ckino Dec.27.2011:             if (2*nqpp*4*nsp*nbas<0) spot%rv_p_oqpp(:)=0.0d0
           if (2*nqpp*4*nsp*nbas<0) spot%rv_a_oqpp(:)=0.0d0

ckino Dec.27.2011:             call dvset ( spot%rv_p_oqpp , 1 , 1 , - 1d0 ) 
           call dvset ( spot%rv_a_oqpp , 1 , 1 , - 1d0 ) 

Cgetarg...info...           integer :: opmpol , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oqpp , integer(8) :: spot%oqpp 1
Cgetarg           call pack2('pot opmpol oqpp',spot,opmpol,oqpp)
ckino Dec.27.2011:             spot%rv_p_opmpol => rv_p_opmpol 

ckino Dec.27.2011:             spot%rv_p_oqpp => rv_p_oqpp 

          call ioqpp(.false.,sctrl,spot)
        endif

Cgetarg...info...           real(8) :: plat(9) , real(8) :: slat%plat 9
Cgetarg         call upack1('lat plat',slat,plat)
        
C$$$         i_copy_size=size(slat%plat) 
C$$$         call dcopy(i_copy_size,slat%plat,1,plat,1) 


C$$$C   --- Create maps ---
C$$$Cgetarg         nmap = igets('ctrl nmap',sctrl)
C$$$         nmap = int(sctrl%nmap) 



C$$$cccccccccccccccccc
C$$$         print *,' xxxxxx 222 nmap=',nmap
C$$$cccccccccccccccccc
 


C$$$        if (nmap .gt. 0) then
C$$$          call upacks('strn map',i1,i2)
C$$$Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
C$$$Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
C$$$Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
C$$$Cgetarg           call upack('array oipc oclabl opos',sarray,oipc,oclabl,opos,0,
C$$$Cgetarg      .      0)
        
C$$$           oipc=sarray%oipc
C$$$           oclabl=sarray%oclabl
C$$$           opos=sarray%opos

C$$$Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
C$$$Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
C$$$Cgetarg           call upack('pot ovrmax orhrmx',spot,ovrmax,orhrmx,0,0,0)
        
C$$$           ovrmax=spot%ovrmax
C$$$           orhrmx=spot%orhrmx

C$$$Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
C$$$Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
C$$$Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
C$$$Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
C$$$Cgetarg           call upack('pot opnu oqnu opp oves',spot,opnu,oqnu,opp,oves,0)
        
C$$$           opnu=spot%opnu
C$$$           oqnu=spot%oqnu
C$$$           opp=spot%opp
C$$$           oves=spot%oves

C$$$          do  10  ic = 1, nclass
C$$$            xx = dval(w(oclabl),ic)
C$$$            call r8tos8(xx,clabl)
C$$$            is = w(oics+ic-1)
C$$$            ib = iclbsj(ic,w(oipc),-nbas,1)
C$$$Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
C$$$Cgetarg             call upack('spec idmod',sspec,is,idmod,0,0,0)
        
C$$$             i_copy_size=size(sspec(is)%idmod) 
C$$$             call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 

C$$$            if (ib .lt. 0) then
C$$$              if (iprint() .ge. 20) print *,
C$$$     .          'CLSPRM: no basis corresponding to class ' //clabl
C$$$              goto 10
C$$$            endif
C$$$Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
C$$$Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
C$$$Cgetarg             call upack('spec z lmxa',sspec,is,z,lmx,0,0)
        
C$$$             z=sspec(is)%z
C$$$             lmx=sspec(is)%lmxa

C$$$            call wrmap(sstrn(i1:i2),nmap,clabl,ic,ib,z,plat,w(opos))
C$$$Cgetarg             call mapq(clabl,ic-1,nl,nsp,lmx,idmod,w(opnu),w(oqnu),
C$$$Cgetarg      .        w(opp),w(opprel),w(orhrmx),w(ovrmax),
C$$$Cgetarg      .        lgors('ctrl lves,1',sctrl),w(oves),w(oinitc))
C$$$             call mapq ( clabl , ic-1 , nl , nsp , lmx , idmod , w ( opnu 
C$$$     .       ) , w ( oqnu ) , w ( opp ) , w ( opprel ) , w ( orhrmx ) , w 
C$$$     .       ( ovrmax ) , iand(1,int(sctrl%lves)) .ne.0 , w ( oves ) , w ( 
C$$$     .       oinitc ) ) 

C$$$   10     continue

C$$$          if (cmdopt('-spin1',6,0,strn))
C$$$     .      call spinav(0,nclass,nl,nsp,w(opnu),w(oqnu))
C$$$          if (cmdopt('-map0',5,0,strn)) then
C$$$            call aiocls(.true.,0,sarray,sctrl,sham,spot,sspec,slat,1,
C$$$     .        nclass)
C$$$            call rx0('mapped moments written')
C$$$          endif
C$$$        endif
      endif

C$$$C --- Allocate tbe arrays ---
C$$$      if (mod(mode,10) .eq. 2) then
C$$$        call defi(oifrlx,-4*nbas)
C$$$Cgetarg...info...           integer :: oifrlx , integer(8) :: sarray%oifrlx 1
C$$$Cgetarg         call pack5('array oifrlx',sarray,oifrlx,0,0,0,0)
      
C$$$         sarray%oifrlx=oifrlx 

C$$$Cgetarg         if (lgors('ctrl lncol,-1',sctrl)) then
C$$$         if ( iand(-1,int(sctrl%lncol)) .ne.0 ) then 

C$$$          call rx('clsprm: is this needed?')
C$$$          call defdr(osop,-nlspc*nsp*9)
C$$$          call defi (oivso,-nclasp*nl)
C$$$CKi          call pack5('ham osop oivso',sham,osop,oivso,0,0,0)
C$$$Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
C$$$Cgetarg...info...           integer :: oivso , integer(8) :: spot%oivso 1
C$$$Cgetarg           call pack5('pot osop oivso',spot,osop,oivso,0,0,0)
      
C$$$           spot%osop=osop 
C$$$           spot%oivso=oivso 

C$$$        endif
C$$$      endif

C$$$Cgetarg       if (lgors('ctrl lstonr,-1',sctrl)) then
C$$$      if ( iand(-1,int(sctrl%lstonr(1))) .ne.0 ) then 

C$$$C        call defdr(oammx,-nclasp)
C$$$C        call defdr(ostni,-nclasp)
C$$$Cgetarg         call defdr(ozos,igets('bz ndos',sbz)*nl*nsp*nclasp)
C$$$         call defdr ( ozos , int(sbz%ndos) * nl * nsp * nclasp ) 

C$$$Cgetarg...info...           integer :: ozos , integer(8) :: sarray%ozos 1
C$$$Cgetarg         call pack5('array ozos',sarray,ozos,0,0,0,0)
      
C$$$         sarray%ozos=ozos 

C$$$      endif

C --- Misc variables ---
C      oips = igets('array oips',sarray)
C      nla = 0
C      do  20  ib = 1, nbas
C        is = w(oips+ib-1)
C        lmx = igetss('spec lmxa',is,sspec)
C        nla = nla + (lmx+1)**2
C   20 continue
C      call pack1('ctrl nla',sctrl,nla)

C     call shstru('array',sarray,0,0)
C     stop

      return

      entry clsprp(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .sstr,sstrn)
C- Copy class data to doubly padded case

Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg       call upack('ctrl nclass nl nspin nspec',sctrl,nclass,nl,nsp,nspec,
Cgetarg      .  0)
      
      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg       call upack1('array nclasp',sarray,nclasp)
      
      nclasp=sarray%nclasp

      if (nclasp .eq. nclass) return
      nclspp = max(2*nclasp-nclass,nspec)
C ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then

        i1 = nclasp-nclass
Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg         call upack('pot opnu oqnu',spot,opnu,oqnu,0,0,0)
        
ckino Dec.27.2011:           rv_p_opnu => spot%rv_p_opnu 

ckino Dec.27.2011:           rv_p_oqnu => spot%rv_p_oqnu 


ckino Dec.27.2011:           call dpscop ( spot%rv_p_opnu , spot%rv_p_opnu , nl * nsp * i1 , nl * nsp 
ckino 12月.27.2011:           call dpscop ( spot%rv_p_opnu , spot%rv_p_opnu , nl * nsp * i1 , nl * nsp 
         call dpscop ( spot%rv_a_opnu , spot%rv_a_opnu , nl * nsp * i1 , nl * nsp 
     .   * nclass + 1 , nl * nsp * nclasp + 1 , 1d0 ) 

ckino Dec.27.2011:           call dpscop ( spot%rv_p_oqnu , spot%rv_p_oqnu , 3 * nl * nsp * i1 , 3 * 
ckino 12月.27.2011:           call dpscop ( spot%rv_p_oqnu , spot%rv_p_oqnu , 3 * nl * nsp * i1 , 3 * 
         call dpscop ( spot%rv_a_oqnu , spot%rv_a_oqnu , 3 * nl * nsp * i1 , 3 * 
     .   nl * nsp * nclass + 1 , 3 * nl * nsp * nclasp + 1 , 1d0 ) 


Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg         call upack('pot oves ovdif ovrmax orhrmx',spot,oves,ovdif,
Cgetarg      .    ovrmax,orhrmx,0)
        
ckino Dec.27.2011:           rv_p_oves => spot%rv_p_oves 

ckino Dec.27.2011:           rv_p_ovdif => spot%rv_p_ovdif 

ckino Dec.27.2011:           rv_p_ovrmax => spot%rv_p_ovrmax 

ckino Dec.27.2011:           rv_p_orhrmx => spot%rv_p_orhrmx 


ckino Dec.27.2011:           call dpscop ( spot%rv_p_oves , spot%rv_p_oves , i1 , nclass + 1 , nclasp 
ckino 12月.27.2011:           call dpscop ( spot%rv_p_oves , spot%rv_p_oves , i1 , nclass + 1 , nclasp 
         call dpscop ( spot%rv_a_oves , spot%rv_a_oves , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

ckino Dec.27.2011:           call dpscop ( spot%rv_p_ovdif , spot%rv_p_ovdif , i1 , nclass + 1 , nclasp 
ckino 12月.27.2011:           call dpscop ( spot%rv_p_ovdif , spot%rv_p_ovdif , i1 , nclass + 1 , nclasp 
         call dpscop ( spot%rv_a_ovdif , spot%rv_a_ovdif , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

ckino Dec.27.2011:           call dpscop ( spot%rv_p_ovrmax , spot%rv_p_ovrmax , i1 , nclass + 1 , nclasp 
ckino Dec.28.2011:           call dpscop ( spot%rv_p_ovrmax , spot%rv_p_ovrmax , i1 , nclass + 1 , nclasp 
         call dpscop ( spot%rv_a_ovrmax , spot%rv_a_ovrmax , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

ckino Dec.27.2011:           call dpscop ( spot%rv_p_orhrmx , spot%rv_p_orhrmx , i1 , nclass + 1 , nclasp 
ckino 12月.27.2011:           call dpscop ( spot%rv_p_orhrmx , spot%rv_p_orhrmx , i1 , nclass + 1 , nclasp 
         call dpscop ( spot%rv_a_orhrmx , spot%rv_a_orhrmx , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 


      else
        call rxi('clsprp: bad mode',mode)
      endif

Cgetarg       end
      end subroutine clsprm 

C================================================================
c take subs/ioqpp.F here
Cgetarg...info...           structure ['ioqpp', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['ioqpp', 'pot', 'spot']
Cgetarg...info...           use_to_add ['ioqpp', 'm_struc_def']
      subroutine ioqpp(lio,sctrl,spot)
      
      use m_struc_def  !Cgetarg

C- File i/o for phi-phi, phi-dot, dot-dot products
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lio: true for write, false for read
Ci          <0 write
Cu Updates
Cu   08 Nov 07 (J. Xu) qpp is complex
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lio
Cgetarg       double precision sctrl(1),spot(1)
      type(s_ctrl)::sctrl
      type(s_pot)::spot

C ... Heap
C ... Local parameters
       integer:: i1 , i2 , nbas , nl , nsp , igets , ifi , fopna , rdm 
     . , ipr 
ckino Dec.27.2011:         real(8),pointer :: rv_p_oqpp(:) =>NULL()


      call getpr(ipr)
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbas nl nspin',sctrl,nbas,nl,nsp,0,0)
      
      nbas=sctrl%nbas
      nl=sctrl%nl
      nsp=sctrl%nspin

      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
Cgetarg       oqpp = igets('pot oqpp',spot)
ckino Dec.27.2011:         rv_p_oqpp => spot%rv_p_oqpp 


      ifi = fopna('qpp',-1,4+8)

      if (lio) then
ckino Dec.27.2011:           call ywrm ( 1 , 'lmasa' , 1 , ifi , ' ' , spot%rv_p_oqpp , 1 , i1 
         call ywrm ( 1 , 'lmasa' , 1 , ifi , ' ' , spot%rv_a_oqpp , 1 , i1 
     .   , i1 , i2 ) 

        if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
        if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
        call pshpr(0)
ckino Dec.27.2011:           if ( rdm ( ifi , 2 , i1 * i2 , ' ' , spot%rv_p_oqpp , i1 , i2 ) .lt. 
         if ( rdm ( ifi , 2 , i1 * i2 , ' ' , spot%rv_a_oqpp , i1 , i2 ) .lt. 
     .   0 ) then 

          if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
ckino Dec.27.2011:             call dvset ( spot%rv_p_oqpp , 1 , 1 , - 1d0 ) 
           call dvset ( spot%rv_a_oqpp , 1 , 1 , - 1d0 ) 

        endif
        call poppr
      endif

      call fclose(ifi)
Cgetarg       end
      end subroutine ioqpp 


C================================================================
c take subs/spinav.F here
      subroutine spinav(mode,nclass,nl,nsp,pnu,qnu)
C- Averages up+down spin moments + pp's for all classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 average spins
Ci         :1 do not average, but exchange spins
Ci   nclass:number of inequivalent classes
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Co Outputs :moments are spin-averaged
Ci   pnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
Ci   qnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
Co   nsp   :set to 1 on output (mode=0)
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   10 Jan 06 Added mode
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nclass,nl,nsp,ic
      double precision pnu(nl,nsp,nclass),qnu(3,nl,nsp,nclass)

      if (nsp .eq. 1) return

      if (mode .gt. 1) then
        call rx('spinav: bad mode')
      elseif (mode .eq. 1) then
        do  ic = 1, nclass
          call dswap(nl,pnu(1,2,ic),1,pnu(1,1,ic),1)
          call dswap(3*nl,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
        enddo
        return
      endif

      do  ic = 1, nclass
        call daxpy(nl,1d0,pnu(1,2,ic),1,pnu(1,1,ic),1)
        call dscal(nl,.5d0,pnu(1,1,ic),1)
        call daxpy(3*nl,1d0,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
      enddo

      do  ic = 2, nclass
        call dcopy(nl,pnu(1,1,ic),1,pnu(1,ic,1),1)
        call dcopy(3*nl,qnu(1,1,1,ic),1,qnu(1,1,ic,1),1)
      enddo
      nsp = 1
      end

Cgetarg...info...           structure ['aiocls', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['aiocls', 'pot', 'spot']
Cgetarg...info...           structure ['aiocls', 'array', 'sarray']
Cgetarg...info...           structure ['aiocls', 'spec', 'sspec']
Cgetarg...info...           structure ['aiocls', 'lat', 'slat']
Cgetarg...info...           use_to_add ['aiocls', 'm_struc_def']
c#define ONE_ATOM_UNIT 1
      subroutine aiocls(lio,mode,sarray,sctrl,sham,spot,sspec,slat,ic1,
     .ic2)
      
      use m_struc_def  !Cgetarg
      use m_globalvariables
C- File I/O atomic data for classes ic1..ic2
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio    F for read, T for write
Ci   mode   ones digit
Ci          1 use data from first class corresponding to same species,
Ci            if data from own class is missing
Ci          2 like 1, but attempt to read data from disk anyway
Ci          3 make no attempt to read from disk but copy info from
Ci            first class corresponding to same species, if data missing
Ci          4 Add 4 if to use default P,Q when not otherwise supplied
Ci          10s digit
Ci          1 assemble background rho
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oclabl ohave oics onrc
Ci     Stored:
Ci     Passed to:
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nl nspin zbak nbas nclass
Ci     Stored:    zbak
Ci     Passed to: lgors
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci          (not used now)
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opp osop ogrrme opnu oqnu orhrmx ovrmax ovintr opmpol
Ci                oves
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idmod p q
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: vol
Ci     Stored:
Ci     Passed to:
Ci   ic1,ic2: range of classes to read data
Cr Remarks
Cr   Right now, read always takes data from file if available
Cu Updates
Cu   09 Nov 07 Corrected sign of default moment (paioc2)
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   26 Apr 03 Added MPI calls
Cu   07 Feb 03 adjusted for redimensioned sop
Cu   30 May 02 Assign better default P
Cu   28 Apr 98 I/O of radial matrix elements of grad
Cu   28 Sep 00 Added setting default P,Q
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lio
      integer mode,ic1,ic2
Cgetarg       double precision sarray(1),sctrl(1),sspec(1),slat(1),sham(1),
Cgetarg      .  spot(1)
      type(s_array)::sarray
      type(s_ctrl)::sctrl
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_pot)::spot

C heap
C Local variables
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,lgors,
     .aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,
     .readpq,readpp,readso,readop,readva,readmp
                    integer::  onrc , n0 , nglob !olpp , olppr , olsop , olgrme , olva , olmp , opot  , ocor ,
ckino Dec.27.2011:                      real(8),pointer :: rv_p_osop(:) =>NULL()

ckino Dec.27.2011:                     real(8),pointer :: rv_p_opprel(:) =>NULL()

ckino Dec.27.2011:                    real(8),pointer :: rv_p_orhrmx(:) =>NULL()

ckino Dec.27.2011:                   real(8),pointer :: rv_p_oqnu(:) =>NULL()

ckino Dec.27.2011:                  real(8),pointer :: rv_p_omp(:) =>NULL()

ckino Dec.27.2011:                 real(8),pointer :: rv_p_opp(:) =>NULL()

ckino Dec.27.2011:                real(8),pointer :: rv_p_opnu(:) =>NULL()

ckino Dec.27.2011:               real(8),pointer :: rv_p_ogrrme(:) =>NULL()

ckino Dec.27.2011:              integer,pointer :: iv_p_ohave(:) =>NULL()

ckino Dec.27.2011:             real(8),pointer :: rv_p_ovrmax(:) =>NULL()

ckino Dec.27.2011:            real(8),pointer :: rv_p_ova(:) =>NULL()

ckino Dec.27.2011:           real(8),pointer :: rv_p_oves(:) =>NULL()

ckino Dec.27.2011:          integer,pointer :: iv_p_oics(:) =>NULL()

ckino Dec.27.2011:         real(8),pointer :: rv_p_oclabl(:) =>NULL()

      parameter(n0=10)
      integer ic,is,lmx,k,nl,nsp,ifi,jfi,fopn,lmxx,nspx,nrx,nr,
     .idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .nclspp,iprint,lgunit,mode0,mode00,lrel
c      integer mpipid,procid
      double precision rhrmx,vrmax(2),ves,z,rmxx,ax,qc,dq,vrmxx(2),
     .sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,dval,
     .zbak(2),dgets,pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100)

      real(8):: olpp(100),olppr(3000),olsop(500),olgrme(500),olva(100),olmp(3000),opot(3000),ocor(3000)

Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg       call upack('ctrl nl nspin nclass',sctrl,nl,nsp,nclass,0,0)
      
      nl=sctrl%nl
      nsp=sctrl%nspin
      nclass=sctrl%nclass

Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg       call upack('pot opp opprel osop ogrrme',spot,opp,opprel,osop,
Cgetarg      .  ogrrme,0)
      
ckino Dec.27.2011:         rv_p_opp => spot%rv_p_opp 

ckino Dec.27.2011:         rv_p_opprel => spot%rv_p_opprel 

ckino Dec.27.2011:         rv_p_osop => spot%rv_p_osop 

ckino Dec.27.2011:         rv_p_ogrrme => spot%rv_p_ogrrme 


Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg       call upack('pot opnu oqnu orhrmx ovrmax',spot,opnu,oqnu,orhrmx,
Cgetarg      .  ovrmax,0)
      
ckino Dec.27.2011:         rv_p_opnu => spot%rv_p_opnu 

ckino Dec.27.2011:         rv_p_oqnu => spot%rv_p_oqnu 

ckino Dec.27.2011:         rv_p_orhrmx => spot%rv_p_orhrmx 

ckino Dec.27.2011:         rv_p_ovrmax => spot%rv_p_ovrmax 


Cgetarg...info...           integer :: ova , integer(8) :: spot%ovintr 1
Cgetarg...info...           integer :: omp , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg       call upack('pot ovintr opmpol oves',spot,ova,omp,oves,0,0)
      
ckino Dec.27.2011:         rv_p_ova => spot%rv_p_ovintr 

ckino Dec.27.2011:         rv_p_omp => spot%rv_p_opmpol 

ckino Dec.27.2011:         rv_p_oves => spot%rv_p_oves 


Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: ohave , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg       call upack('array nclasp oclabl ohave oics',sarray,nclasp,oclabl,
Cgetarg      .  ohave,oics,0)
      
      nclasp=sarray%nclasp
ckino Dec.27.2011:         rv_p_oclabl => sarray%rv_p_oclabl 

ckino Dec.27.2011:         iv_p_ohave => sarray%iv_p_ohave 

ckino Dec.27.2011:         iv_p_oics => sarray%iv_p_oics 


Changenglob      lrel = nglob('lrel')
      lrel = globalvariables%lrel

      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)

C ... MPI: only master does sphere program
c      procid = mpipid(1)
c      if (procid .eq. 0) then

c      olppr = 1
c      call defdr(olpp,100)
c      call defdr(olppr,3000)
c      call defdr(olsop,500)
c      call defdr(olgrme,500)
c      call defdr(olva,100)
c      call defdr(olmp,3000)
c      call defdr(opot,3000)
c      call defdr(ocor,3000)
      i2 = ic2
      if (i2 .eq. 0) i2 = nclasp
c     call awrit2('%n:1i',' ',100,6,nclasp,w(oics))
      do  10  ic = ic1, i2
cki        is = w(oics+ic-1)
ckino Dec.27.2011:          is=sarray%iv_p_oics(ic)
ckino Dec.27.2011:          is=sarray%iv_p_oics(ic)
        is=sarray%iv_a_oics(ic)
ckino Dec.27.2011:           icmap = iclbsj ( is , sarray%iv_p_oics , - nclasp , 1 ) 
ckino Dec.27.2011:           icmap = iclbsj ( is , sarray%iv_p_oics , - nclasp , 1 ) 
         icmap = iclbsj ( is , sarray%iv_a_oics , - nclasp , 1 ) 

        if (icmap .eq. ic .or. mode00 .eq. 0) icmap = 0
Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
Cgetarg...info...           real(8) :: pdf(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: qdf(n0,2) , real(8) :: sspec%q 20
Cgetarg         call upack('spec lmxa idmod p q',sspec,is,lmx,idmod,pdf,qdf)
        
        lmx=sspec(is)%lmxa
        i_copy_size=size(sspec(is)%idmod) 
        call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
        i_copy_size=size(sspec(is)%p) 
        call dcopy(i_copy_size,sspec(is)%p,1,pdf,1) 
        i_copy_size=size(sspec(is)%q) 
        call dcopy(i_copy_size,sspec(is)%q,1,qdf,1) 


        call dpzero(pnuloc,100)
        call dpzero(qnuloc,100)

C   --- Open the atom file ---
        outs1 = ' '
        outs2 = ' '
        outs3 = ' '
ckino Dec.27.2011:           call r8tos8 ( dval ( sarray%rv_p_oclabl , ic ) , clabl ) 
ckino Dec.27.2011:           call r8tos8 ( dval ( sarray%rv_p_oclabl , ic ) , clabl ) 
         call r8tos8 ( dval ( sarray%rv_a_oclabl , ic ) , clabl ) 

        if (mode00 .ne. 3) then
c#if ONE_ATOM_UNIT
          ifi = fopn(clabl)
c#else
c         ifi = fopna(clabl,30+ic,0)
c#endif
        endif

C  --- Copy what is passed through to holding arrays ---
        k = nl*nsp
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 1 , ic , icmap , havepq , jc ) 
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 1 , ic , icmap , havepq , jc ) 
         call pvaioc ( sarray%iv_a_oipq , 1 , ic , icmap , havepq , jc ) 

        readpq = .not. havepq .or. mode00 .eq. 2 .and. jc .ne. ic
        readpq = readpq .and. mode00 .lt. 3
        if (havepq) then
ckino Dec.27.2011:             call dpscop ( spot%rv_p_opnu , pnuloc , k , 1 + ( jc - 1 ) * k , 1 
ckino 12月.27.2011:             call dpscop ( spot%rv_p_opnu , pnuloc , k , 1 + ( jc - 1 ) * k , 1 
           call dpscop ( spot%rv_a_opnu , pnuloc , k , 1 + ( jc - 1 ) * k , 1 
     .     , 1d0 ) 

ckino Dec.27.2011:             call dpscop ( spot%rv_p_oqnu , qnuloc , 3 * k , 1 + ( jc - 1 ) * 3 
ckino 12月.27.2011:             call dpscop ( spot%rv_p_oqnu , qnuloc , 3 * k , 1 + ( jc - 1 ) * 3 
           call dpscop ( spot%rv_a_oqnu , qnuloc , 3 * k , 1 + ( jc - 1 ) * 3 
     .     * k , 1 , 1d0 ) 

ckino Dec.27.2011:             call dpscop ( spot%rv_p_orhrmx , rhrmx , 1 , jc , 1 , 1d0 ) 
ckino 12月.27.2011:             call dpscop ( spot%rv_p_orhrmx , rhrmx , 1 , jc , 1 , 1d0 ) 
           call dpscop ( spot%rv_a_orhrmx , rhrmx , 1 , jc , 1 , 1d0 ) 

ckino Dec.27.2011:             call dpscop ( spot%rv_p_ovrmax , vrmax , 2 , 2 * jc - 1 , 1 , 1d0 
ckino Dec.28.2011:             call dpscop ( spot%rv_p_ovrmax , vrmax , 2 , 2 * jc - 1 , 1 , 1d0 
           call dpscop ( spot%rv_a_ovrmax , vrmax , 2 , 2 * jc - 1 , 1 , 1d0 
     .     ) 

          if (ic.ne.jc) call awrit0('%a pq,',outs1,len(outs1),0)
        endif
        k = 6*nl*nsp
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 2 , ic , icmap , havepp , jc ) 
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 2 , ic , icmap , havepp , jc ) 
         call pvaioc ( sarray%iv_a_oipq , 2 , ic , icmap , havepp , jc ) 

        readpp = .not. havepp .or. mode00 .eq. 2 .and. jc .ne. ic
        readpp = readpp .and. mode00 .lt. 3
        if (havepp) then
ckino Dec.27.2011:             call dpscop ( rv_p_opp , olpp , k , 1 + ( jc - 1 ) * k 
           call dpscop ( spot%rv_a_opp , olpp , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (lrel .eq. 2) then
            k = 4*nl*2*nl*2*2
ckino Dec.27.2011:               call dpscop ( spot%rv_p_opprel , olppr , k , 1 + ( jc - 1 ) 
ckino 12月.27.2011:               call dpscop ( spot%rv_p_opprel , olppr , k , 1 + ( jc - 1 ) 
             call dpscop ( spot%rv_a_opprel , olppr , k , 1 + ( jc - 1 ) 
     .       * k , 1 , 1d0 ) 

          endif
ckino Dec.27.2011:             call dpscop ( spot%rv_p_oves , ves , 1 , jc , 1 , 1d0 ) 
ckino 12月.27.2011:             call dpscop ( spot%rv_p_oves , ves , 1 , jc , 1 , 1d0 ) 
           call dpscop ( spot%rv_a_oves , ves , 1 , jc , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a pp,',outs1,len(outs1),0)
        endif
        k = (nl*nsp)**2
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 8 , ic , icmap , haveva , jc ) 
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 8 , ic , icmap , haveva , jc ) 
         call pvaioc ( sarray%iv_a_oipq , 8 , ic , icmap , haveva , jc ) 

cki        haveva = haveva .and. ova .gt. 0
ckino Dec.27.2011:          haveva = haveva .and. associated(spot%rv_p_ovintr)
ckino Dec.28.2011:          haveva = haveva .and. associated(spot%rv_p_ovintr)
        haveva = haveva .and. allocated(spot%rv_a_ovintr)
        readva = .not. haveva .or. mode00 .eq. 2 .and. jc .ne. ic
cki        readva = readva .and. ova .gt. 0 .and. mode00 .lt. 3
ckino Dec.27.2011:          readva = readva .and. associated(spot%rv_p_ovintr) .and. mode00 .lt. 3
ckino Dec.28.2011:          readva = readva .and. associated(spot%rv_p_ovintr) .and. mode00 .lt. 3
        readva = readva .and. allocated(spot%rv_a_ovintr) .and. mode00 .lt. 3
        if (haveva) then
ckino Dec.27.2011:             call dpscop ( spot%rv_p_ovintr , olva , k , 1 + ( jc - 1 ) * k 
ckino Dec.28.2011:             call dpscop ( spot%rv_p_ovintr , olva , k , 1 + ( jc - 1 ) * k 
           call dpscop ( spot%rv_a_ovintr , olva , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a va,',outs1,len(outs1),0)
        endif
        k = nl**2*(2*nl-1)*3*nsp
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 16 , ic , icmap , havemp , jc ) 
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 16 , ic , icmap , havemp , jc ) 
         call pvaioc ( sarray%iv_a_oipq , 16 , ic , icmap , havemp , jc ) 

cki        havemp = havemp .and. omp .gt. 0
ckino Dec.27.2011:          havemp = havemp .and. associated(spot%rv_p_opmpol)
ckino 12月.27.2011:          havemp = havemp .and. associated(spot%rv_p_opmpol)
        havemp = havemp .and. allocated(spot%rv_a_opmpol)
        readmp = .not. havemp .or. mode00 .eq. 2 .and. jc .ne. ic
cki        readmp = readmp .and. omp .gt. 0 .and. mode00 .lt. 3
ckino Dec.27.2011:          readmp = readmp .and. associated(spot%rv_p_opmpol) .and. mode00 .lt. 3
ckino 12月.27.2011:          readmp = readmp .and. associated(spot%rv_p_opmpol) .and. mode00 .lt. 3
        readmp = readmp .and. allocated(spot%rv_a_opmpol) .and. mode00 .lt. 3
        if (havemp) then
ckino Dec.27.2011:             call dpscop ( spot%rv_p_opmpol , olmp , k , 1 + ( jc - 1 ) * k 
ckino 12月.27.2011:             call dpscop ( spot%rv_p_opmpol , olmp , k , 1 + ( jc - 1 ) * k 
           call dpscop ( spot%rv_a_opmpol , olmp , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a mp,',outs1,len(outs1),0)
        endif
        k = nl*nsp*nsp*9
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 4 , ic , icmap , haveso , jc ) 
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 4 , ic , icmap , haveso , jc ) 
         call pvaioc ( sarray%iv_a_oipq , 4 , ic , icmap , haveso , jc ) 

cki        haveso = haveso .and. osop.gt. 0
ckino Dec.27.2011:          haveso = haveso .and. associated(spot%rv_p_osop)
ckino 12月.27.2011:          haveso = haveso .and. associated(spot%rv_p_osop)
        haveso = haveso .and. allocated(spot%rv_a_osop)
        readso = .not. haveso .or. mode00 .eq. 2 .and. jc .ne. ic
cki        readso = readso .and. osop .gt. 0 .and. mode00 .lt. 3
ckino Dec.27.2011:          readso = readso .and. associated(spot%rv_p_osop) .and. mode00 .lt. 3
ckino 12月.27.2011:          readso = readso .and. associated(spot%rv_p_osop) .and. mode00 .lt. 3
        readso = readso .and. allocated(spot%rv_a_osop) .and. mode00 .lt. 3
        if (haveso) then
ckino Dec.27.2011:             call dpscop ( spot%rv_p_osop , olsop , k , 1 + ( jc - 1 ) * 
ckino 12月.27.2011:             call dpscop ( spot%rv_p_osop , olsop , k , 1 + ( jc - 1 ) * 
           call dpscop ( spot%rv_a_osop , olsop , k , 1 + ( jc - 1 ) * 
     .     k , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a sop,',outs1,len(outs1),0)
        endif
        k = 16*nl*nsp
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 32 , ic , icmap , haveop , jc ) 
ckino Dec.27.2011:           call pvaioc ( sarray%iv_p_ohave , 32 , ic , icmap , haveop , jc ) 
         call pvaioc ( sarray%iv_a_oipq , 32 , ic , icmap , haveop , jc ) 

cki        haveop = haveop .and. ogrrme .gt. 0
ckino Dec.27.2011:          haveop = haveop .and. associated(spot%rv_p_ogrrme)
ckino Dec.28.2011:          haveop = haveop .and. associated(spot%rv_p_ogrrme)
        haveop = haveop .and. allocated(spot%rv_a_ogrrme)
        readop = .not. haveop .or. mode00 .eq. 2 .and. jc .ne. ic
cki        readop = readop .and. ogrrme .gt. 0 .and. mode00 .lt. 3
ckino Dec.27.2011:          readop = readop .and. associated(spot%rv_p_ogrrme) .and. mode00 .lt. 3
ckino Dec.28.2011:          readop = readop .and. associated(spot%rv_p_ogrrme) .and. mode00 .lt. 3
        readop = readop .and. allocated(spot%rv_a_ogrrme) .and. mode00 .lt. 3
        if (haveop) then
ckino Dec.27.2011:             call dpscop ( spot%rv_p_ogrrme , olgrme , k , 1 + ( jc - 1 ) 
ckino Dec.28.2011:             call dpscop ( spot%rv_p_ogrrme , olgrme , k , 1 + ( jc - 1 ) 
           call dpscop ( spot%rv_a_ogrrme , olgrme , k , 1 + ( jc - 1 ) 
     .     * k , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a opp,',outs1,len(outs1),0)
        endif

C   --- File WRITE ---
        if (lio) then
          lgen = .false.
          lpot = .false.
          lcor = .false.
C     ... Pick up GEN and POT, if available, to save again
          if (scat(iabs(ifi),'GEN:',':',.true.)) then
            lgen = aiogen(alabel,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,dq,
     .      vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          endif
          if (scat(iabs(ifi),'POT:',':',.true.)) then
            read(ifi,102) nr,nsp,a,rmax
  102       format(2i5,2f12.5)
            lpot = aiopot(nr,nsp,a,rmax,-99d0,opot,ifi)
          endif
          lcor = aiocor(nr,nsp,a,rmxx,ocor,sumec,sumtc,ifi)

          rewind ifi
          jfi = -ifi
          if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .    dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
          if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,idmod,
     .    nl,lmx,nsp,rhrmx,vrmax,jfi)
          if (havepp) sw = aiopar(clabl,lrel,olpp,olppr,ves,nl,
     .    lmx,nsp,jfi)
          if (haveva) sw = aiova(clabl,olva,nl,lmx,nsp,jfi)
          if (havemp) sw = aiomp(clabl,olmp,nl,2*nl-2,nsp,jfi)
          if (haveso) sw = aiosop(clabl,olsop,nl,lmx,nsp,jfi)
          if (haveop) sw = aiorme(clabl,olgrme,nl,nsp,jfi)
          if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,opot,jfi)
          if (lcor) lcor = aiocor(nr,nsp,a,rmxx,ocor,sumec,sumtc,jfi)

C   --- File READ ---
        else

C     ... Copy whatever is available on disk to holding arrays
          if (readpq .or. mode0 .ge. 4) then
            if (readpq) rewind ifi
            if (readpq) readpq = aiomom(clabl,pnuloc,qnuloc,idmod,
     .      nl,lmx,nsp,rhrmx,vrmax,ifi)
C           Couldn't read from atom file ; take default values
            if (readpq) call awrit0('%a pq,',outs2,len(outs2),0)
            if (mode0 .ge. 4 .and. .not. (readpq .or. havepq)) then
C             call dmcpy(pdf,n0,1,pnuloc,nl,1,nl,nsp)
              call paioc2(nsp,nl,n0,pdf,qdf,pnuloc,qnuloc)
              call awrit0('%a pq,',outs3,len(outs2),0)
              call dvset(vrmax,1,2,-.7d0)
              rhrmx = .1d0
              readpq = .true.
            endif
          endif
          if (readpp) then
            readpp = aiopar(clabl,lrel,olpp,olppr,ves,nl,lmx,nsp,
     .      ifi)
            if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
          endif
cki          if (readso .and. osop .gt. 0) then
ckino Dec.27.2011:            if (readso .and. associated(spot%rv_p_osop) ) then
ckino 12月.27.2011:            if (readso .and. associated(spot%rv_p_osop) ) then
          if (readso .and. allocated(spot%rv_a_osop) ) then
            readso = aiosop(clabl,olsop,nl,lmx,nsp,ifi)
            if (readso) call awrit0('%a so,',outs2,len(outs2),0)
          endif
cki          if (readop .and. ogrrme .gt. 0) then
ckino Dec.27.2011:            if (readop .and. associated(spot%rv_p_ogrrme)) then
ckino Dec.28.2011:            if (readop .and. associated(spot%rv_p_ogrrme)) then
          if (readop .and. allocated(spot%rv_a_ogrrme)) then
            readop = aiorme(clabl,olgrme,nl,nsp,ifi)
            if (readop) call awrit0('%a op,',outs2,len(outs2),0)
          endif
cki          if (readva .and. ova .gt. 0) then
ckino Dec.27.2011:            if (readva .and. associated(spot%rv_p_ovintr) ) then
ckino Dec.28.2011:            if (readva .and. associated(spot%rv_p_ovintr) ) then
          if (readva .and. allocated(spot%rv_a_ovintr) ) then
            readva = aiova(clabl,olva,nl,lmx,nsp,ifi)
            if (readva) call awrit0('%a va,',outs2,len(outs2),0)
          endif
cki          if (readmp .and. omp .gt. 0) then
ckino Dec.27.2011:            if (readmp .and. associated(spot%rv_p_opmpol) ) then
ckino 12月.27.2011:            if (readmp .and. associated(spot%rv_p_opmpol) ) then
          if (readmp .and. allocated(spot%rv_a_opmpol) ) then
            readmp = aiomp(clabl,olmp,nl,2*nl-2,nsp,ifi)
            if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
          endif

C     ... Update what parameters are available
cki          w(ohave+ic-1) = isw(haveop.or.readop)*32+
ckino Dec.27.2011:            sarray%iv_p_ohave(ic) = isw(haveop.or.readop)*32+
ckino Dec.27.2011:            sarray%iv_p_ohave(ic) = isw(haveop.or.readop)*32+
          sarray%iv_a_oipq(ic) = isw(haveop.or.readop)*32+
     .    isw(havemp.or.readmp)*16+
     .    isw(haveva.or.readva)*8 +
     .    isw(haveso.or.readso)*4 +
     .    isw(havepp.or.readpp)*2 +
     .    isw(havepq.or.readpq)*1 +
cki     .    w(ohave+ic-1) - bitand(w(ohave+ic-1),63)
ckino Dec.27.2011:       .    sarray%iv_p_ohave(ic) - bitand(sarray%iv_p_ohave(ic),63)
ckino Dec.27.2011:       .    sarray%iv_p_ohave(ic) - bitand(sarray%iv_p_ohave(ic),63)
     .    sarray%iv_a_oipq(ic) - bitand(sarray%iv_a_oipq(ic),63)

          k = nl*nsp
          if (havepq .or. readpq) then
ckino Dec.27.2011:               call dpscop ( pnuloc , spot%rv_p_opnu , k , 1 , 1 + ( ic - 1 ) * k 
ckino 12月.27.2011:               call dpscop ( pnuloc , spot%rv_p_opnu , k , 1 , 1 + ( ic - 1 ) * k 
             call dpscop ( pnuloc , spot%rv_a_opnu , k , 1 , 1 + ( ic - 1 ) * k 
     .       , 1d0 ) 

ckino Dec.27.2011:               call dpscop ( qnuloc , spot%rv_p_oqnu , 3 * k , 1 , 1 + ( ic - 1 ) 
ckino 12月.27.2011:               call dpscop ( qnuloc , spot%rv_p_oqnu , 3 * k , 1 , 1 + ( ic - 1 ) 
             call dpscop ( qnuloc , spot%rv_a_oqnu , 3 * k , 1 , 1 + ( ic - 1 ) 
     .       * 3 * k , 1d0 ) 

ckino Dec.27.2011:               call dpscop ( rhrmx , spot%rv_p_orhrmx , 1 , 1 , ic , 1d0 ) 
ckino 12月.27.2011:               call dpscop ( rhrmx , spot%rv_p_orhrmx , 1 , 1 , ic , 1d0 ) 
             call dpscop ( rhrmx , spot%rv_a_orhrmx , 1 , 1 , ic , 1d0 ) 

ckino Dec.27.2011:               call dpscop ( vrmax , spot%rv_p_ovrmax , 2 , 1 , 2 * ic - 1 , 1d0 
ckino Dec.28.2011:               call dpscop ( vrmax , spot%rv_p_ovrmax , 2 , 1 , 2 * ic - 1 , 1d0 
             call dpscop ( vrmax , spot%rv_a_ovrmax , 2 , 1 , 2 * ic - 1 , 1d0 
     .       ) 

          endif
          k = 6*nl*nsp
          if (havepp .or. readpp) then
ckino Dec.27.2011:               call dpscop ( olpp , rv_p_opp , k , 1 , 1 + ( ic - 1 ) 
             call dpscop ( olpp , spot%rv_a_opp , k , 1 , 1 + ( ic - 1 ) 
     .       * k , 1d0 ) 

            if (lrel .eq. 2) then
              k = 4*nl*2*nl*2*2
ckino Dec.27.2011:                 call dpscop ( olppr , spot%rv_p_opprel , k , 1 , 1 + ( ic - 
ckino 12月.27.2011:                 call dpscop ( olppr , spot%rv_p_opprel , k , 1 , 1 + ( ic - 
               call dpscop ( olppr , spot%rv_a_opprel , k , 1 , 1 + ( ic - 
     .         1 ) * k , 1d0 ) 

            endif
ckino Dec.27.2011:               call dpscop ( ves , spot%rv_p_oves , 1 , 1 , ic , 1d0 ) 
ckino 12月.27.2011:               call dpscop ( ves , spot%rv_p_oves , 1 , 1 , ic , 1d0 ) 
             call dpscop ( ves , spot%rv_a_oves , 1 , 1 , ic , 1d0 ) 

          endif
          k = (nl*nsp)**2
cki          if ((haveva .or. readva) .and. ova .gt. 0)
ckino Dec.27.2011:             if ( ( haveva .or. readva ) .and. associated ( spot%rv_p_ovintr ) ) call 
ckino Dec.28.2011:             if ( ( haveva .or. readva ) .and. associated ( spot%rv_p_ovintr ) ) call 
           if ( ( haveva .or. readva ) .and. allocated(spot%rv_a_ovintr ) ) call 
ckino Dec.27.2011:       .     dpscop ( olva  , spot%rv_p_ovintr , k , 1 , 1 + ( ic - 1 ) * k , 
ckino Dec.28.2011:       .     dpscop ( olva  , spot%rv_p_ovintr , k , 1 , 1 + ( ic - 1 ) * k , 
     .     dpscop ( olva  , spot%rv_a_ovintr , k , 1 , 1 + ( ic - 1 ) * k , 
     .     1d0 ) 

          k = nl**2*(2*nl-1)*3*nsp
cki          if ((havemp .or. readmp) .and. omp .gt. 0)
ckino Dec.27.2011:             if ( ( havemp .or. readmp ) .and. associated ( spot%rv_p_opmpol ) ) call 
ckino 12月.27.2011:             if ( ( havemp .or. readmp ) .and. associated ( spot%rv_p_opmpol ) ) call 
           if ( ( havemp .or. readmp ) .and. allocated(spot%rv_a_opmpol ) ) call 
ckino Dec.27.2011:       .     dpscop ( olmp  , spot%rv_p_opmpol , k , 1 , 1 + ( ic - 1 ) * k , 
ckino 12月.27.2011:       .     dpscop ( olmp  , spot%rv_p_opmpol , k , 1 , 1 + ( ic - 1 ) * k , 
     .     dpscop ( olmp  , spot%rv_a_opmpol , k , 1 , 1 + ( ic - 1 ) * k , 
     .     1d0 ) 

          k = nl*nsp*nsp*9
cki          if ((haveso .or. readso) .and. osop .gt. 0)
ckino Dec.27.2011:             if ( ( haveso .or. readso ) .and. associated ( spot%rv_p_osop ) ) 
ckino 12月.27.2011:             if ( ( haveso .or. readso ) .and. associated ( spot%rv_p_osop ) ) 
           if ( ( haveso .or. readso ) .and. allocated(spot%rv_a_osop ) ) 
ckino Dec.27.2011:       .     call dpscop ( olsop  , spot%rv_p_osop , k , 1 , 1 + ( ic - 1 
ckino 12月.27.2011:       .     call dpscop ( olsop  , spot%rv_p_osop , k , 1 , 1 + ( ic - 1 
     .     call dpscop ( olsop  , spot%rv_a_osop , k , 1 , 1 + ( ic - 1 
     .     ) * k , 1d0 ) 

          k = 16*nl*nsp
ckino Dec.27.2011:             if ( ( haveop .or. readop ) .and. associated ( spot%rv_p_ogrrme ) 
ckino Dec.28.2011:             if ( ( haveop .or. readop ) .and. associated ( spot%rv_p_ogrrme ) 
           if ( ( haveop .or. readop ) .and. allocated(spot%rv_a_ogrrme ) 
ckino Dec.27.2011:       .     ) call dpscop ( olgrme  , spot%rv_p_ogrrme , k , 1 , 1 + ( ic 
ckino Dec.28.2011:       .     ) call dpscop ( olgrme  , spot%rv_p_ogrrme , k , 1 , 1 + ( ic 
     .     ) call dpscop ( olgrme  , spot%rv_a_ogrrme , k , 1 , 1 + ( ic 
     .     - 1 ) * k , 1d0 ) 

        endif

        if (mode00 .lt. 3) call fclr(clabl,ifi)
        if (iprint() .gt. 40) then
          outs = ' '
          if (outs1 .ne. ' ') then
            call awrit1('%x '//clabl//'%a: copied '//outs1//'%a%b from '//
     .      'class %i',outs,len(outs),0,jc)
          endif
          if (outs2 .ne. ' ') then
            call awrit0('%x '//clabl//'%a: read '//outs2//'%a%b from '//
     .      'disk',outs,len(outs),0)
          endif
          if (outs3 .ne. ' ') then
            if (outs2 .eq.  ' ') then
              call awrit0('%x '//clabl//'%a: use defaults for: '//outs3//
     .        '%a%b',outs,len(outs),0)
            else
              call awrit0('%a; use defaults for: '//outs3//'%a%b',
     .        outs,len(outs),0)
            endif
          endif
          if (outs .eq. ' ')
     .    call awrit0(' '//clabl//'%a: nothing read',outs,len(outs),0)
          call awrit0(' aiocls class'//outs,' ',-len(outs),lgunit(1))
        endif
   10 continue
C     call awrit2('%n:1i',' ',100,6,nclasp,w(ohave))
c      call rlse(olpp)
C     End of MPI master-only branch
c      endif
      end subroutine aiocls 

      subroutine pvaioc(initc,mask,ic0,icmap,lhave,ic)

C- Find whether data avaiable either in class or mapped class
C     implicit none
      logical lhave
      integer initc(1),mask,ic0,icmap
      integer ic
      ic = ic0
      lhave = mod(initc(ic)/mask,2) .eq. 1
      if (.not. lhave .and. icmap .ne. 0) then
        lhave = mod(initc(icmap)/mask,2) .eq. 1
        ic = icmap
      endif
Cgetarg       end
      end subroutine pvaioc 

      subroutine paioc2(nsp,nl,n0,pat,qat,pnu,qnu)

C- Widget to copy pat,qat to pnu,qnu
C     implicit none
      integer n0,nl,nsp
      double precision pat(n0,2),qat(n0,2),pnu(nl,nsp),qnu(3,nl,nsp)
      integer i,il

      do  i = 1, nsp
        do  il = 1, nl
C       pnu(il,i) = int(pat(il,i)) + .5d0
          pnu(il,i) = pat(il,i)
          qnu(1,il,i) = qat(il,1)/nsp
          if (nsp .eq. 2) then
            if (pat(il,i) .eq. 0) pnu(il,i) = pnu(il,1)
            qnu(1,il,i) = qat(il,1)/nsp + qat(il,2)/2*dble(3-2*i)
          endif
          qnu(2,il,i) = 0d0
          qnu(3,il,i) = 0d0
        enddo
      enddo
Cgetarg       end
      end subroutine paioc2 


