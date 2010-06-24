Cchp1 cmdvar= integer,iv_p_,oinitc
Cchp1 cmdvar= integer,iv_p_,ohave
Cchp1 cmdvar= real(8),rv_p_,ovrmax
Cchp1 cmdvar= real(8),rv_p_,ova
Cchp1 cmdvar= real(8),rv_p_,ovintr
Cchp1 cmdvar= real(8),rv_p_,oves
Cchp1 cmdvar= real(8),rv_p_,ovdif
Cchp1 cmdvar= integer,iv_p_,oics
Cchp1 cmdvar= real(8),rv_p_,oclabl
C For lmdos and lmchk. Single core mode only. No MPI modes.
      program lmfutil
C- Main program for lmf and so on.
      use m_struc_def
      use m_rdctrl2_func, only: rdctrl2
      implicit none
c      integer procid, master, mpipid, nproc

C Heap allocation
      integer wksize
      parameter(wksize= 10 000 000)
      integer w(wksize)
C     Next two lines guarantee w is aligned along a d.p. boundary
      double precision ws
      equivalence (ws,w(1))
      common /w/ w

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

      type(s_spec),pointer:: v_sspec(:)
      type(s_site),pointer:: v_ssite(:)

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
      character:: recrd_(0:mxrecs*recln0-1)
      character*(1000):: recrd
      equivalence ( recrd, recrd_ )

      real(8),parameter::    NULLR =-99999
      real(8)::pwmode,pwemin,pwemax

      integer::ic
      real(8):: bz_w,dval!test

C --- Version ---
      vrsion(1) = 'LM'
      vrsion(2) = ' '
      vn(1) = 7.00d0
      vn2(1) = 0
      vn(2) = 0.00d0
      vn2(2) = 0
      vn(1) = vn(1) + vn2(1)/1d5
      vn(2) = vn(2) + vn2(2)/1d5
      stdo = lgunit(1)
      i = dglob('stdo',dble(stdo),1)
c      master = 0
c      procid = mpipid(1)
c      nproc  = mpipid(0)

C --- Help ---
      swtmp = .false.
      if (swtmp .or. cmdopt('--h',3,0,outs))
     .call lmhelp(prgnam,vn(2),vn2(2),wksize)
      if (cmdopt('--version',9,0,outs)) then
        write(*,'(f4.2)') vn(2)
        call cexit(0,1)
      endif

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
     .slabl,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,
     .v_sarry,v_smove,v_stb,sstrn)

cccccccccccccccccc
c      print *,' uuuu000 222 nmap=',int(v_sctrl%nmap) 
cccccccccccccccccc

C --- Lattice setup ---
      call setcg(v_slat,8,12)
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
      call mksym(lc,slabl,strn,v_sctrl,v_slat,v_ssite,v_sarry)
c nmap is not zero cleared ---> takao think it is zero for all cases here.

C --- Allocate permanent class arrays, maps and other initialization ---
      call clsprm(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .v_slat,v_sbz,v_sstr,sstrn)

C --- Read available class parameters from file ---
#if LMCHK
      call aiocls(.false.,0,v_sarry,v_sctrl,v_sham,v_spot,
     .v_sspec,v_slat,1,0)
      call aiocls(.false.,17,v_sarry,v_sctrl,v_sham,v_spot,
     .v_sspec,v_slat,1,0)
      if (lc .ge. 20)
     .call clsprp(1,v_sctrl,v_sarry,v_sham,v_spot,v_sspec,
     .v_slat,v_sbz,v_sstr,sstrn)
#endif

      if (cmdopt('--sfill',7,0,strn)) then
        call rx('no --sfill option in lm7K')
      endif

C ... Set default values for species data
      call defspc(v_sspec)

C ... quit after SHOW
Cgetarg       if (igets('ctrl quit',v_sctrl) .eq. 1) then
      if ( int(v_sctrl%quit) .eq. 1 ) then 

        call info0(0,0,0,' '//prgnam//'%a:  Q=SHOW encountered')
        call rx0(prgnam)
      endif

#if LMCHK & FP
      if (cmdopt('--fp',4,0,strn)) then
        call fpchk(v_sspec,v_ssite)
        call cexit(0,1)
      endif
#endif

#if LMCHK
      if (cmdopt('--findes',8,0,strn)) auxmod = 128
      call lmaux(prgnam,v_sctrl,v_sham,v_spot,v_slat,v_sspec,
     .v_ssite,v_sstr,v_sarry,v_sbz,slabl,sstrn,auxmod,wksize)
#endif

#if LMDOS
      call asados(prgnam,v_sbz,v_sctrl,v_sham,v_slat,
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
Cchp1       integer obxc,oinitc,opnu,oqnu,oqc,oqt,oves,ovdif,orhrmx,ovrmax,
Cchp1      .orhos,osop,ogroup,ovintr,opp,opprel,oifrlx,oivso,oipc,oclabl,
Cchp1      .opos,oics,ozos,ormax,opmpol,oqpp,ogrrme,n0
Cchp1        integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , ovdif 
Cchp1      . , orhrmx , ovrmax , orhos , osop , ogroup , ovintr , opp , opprel 
Cchp1      . , oifrlx , oivso , oipc , opos , oics , ozos , ormax , opmpol 
Cchp1      . , oqpp , ogrrme , n0 
Cchp1         integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , ovdif 
Cchp1      .  , orhrmx , ovrmax , orhos , osop , ogroup , ovintr , opp , opprel 
Cchp1      .  , oifrlx , oivso , oipc , opos , ozos , ormax , opmpol , oqpp 
Cchp1      .  , ogrrme , n0 
Cchp1          integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , orhrmx 
Cchp1      .   , ovrmax , orhos , osop , ogroup , ovintr , opp , opprel , oifrlx 
Cchp1      .   , oivso , oipc , opos , ozos , ormax , opmpol , oqpp , ogrrme 
Cchp1      .   , n0 
Cchp1           integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
Cchp1      .    ovrmax , orhos , osop , ogroup , ovintr , opp , opprel , oifrlx 
Cchp1      .    , oivso , oipc , opos , ozos , ormax , opmpol , oqpp , ogrrme 
Cchp1      .    , n0 
Cchp1            integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
Cchp1      .     ovrmax , orhos , osop , ogroup , opp , opprel , oifrlx , oivso 
Cchp1      .     , oipc , opos , ozos , ormax , opmpol , oqpp , ogrrme , n0 
Cchp1             integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
Cchp1      .      orhos , osop , ogroup , opp , opprel , oifrlx , oivso , oipc 
Cchp1      .      , opos , ozos , ormax , opmpol , oqpp , ogrrme , n0 
             integer:: obxc , opnu , oqnu , oqc , oqt , orhrmx , orhos , osop 
     .       , ogroup , opp , opprel , oifrlx , oivso , oipc , opos , ozos 
     .       , ormax , opmpol , oqpp , ogrrme , n0 
             integer,pointer :: iv_p_oinitc(:)

            real(8),pointer :: rv_p_ovrmax(:)

           real(8),pointer :: rv_p_ovintr(:)

          real(8),pointer :: rv_p_oves(:)

         real(8),pointer :: rv_p_ovdif(:)

        integer,pointer :: iv_p_oics(:)

       real(8),pointer :: rv_p_oclabl(:)

      parameter (n0=10)
      integer nbas,nclass,nspec,nl,nsp,nlspc,igets,nclasp,i1,i2,nmap,ic,
     .is,ib,iclbsj,idmod(n0),iprint,lmx,nclspp,nqpp,lnsph,isw,neula,
     .nbf,nrhos,nlibu
C     integer nmto,oppn
      double precision plat(9),z,xx,dval
C ... Heap
      integer w(1)
      common /w/ w

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
        call defdr(opnu,-nlspc)
        call defdr(oqnu,-3*nlspc)
        call defdr(opp,-6*nlspc)
        call defdr(opprel,-32*nl*nl*max(nclspp,nspec))
        call defdr(oqc,-nclasp)
        call defdr(oqt,-nclasp)
Cchp1         call defi (oinitc,-nclasp)
         allocate(iv_p_oinitc(abs(-nclasp)))
         if (-nclasp<0) iv_p_oinitc(:)=0

        call defi (ogroup,-2*nclasp)
Cchp1         call defdr(oves,-nclspp)
         allocate(rv_p_oves(abs(-nclspp)))
         if (-nclspp<0) rv_p_oves(:)=0.0d0

Cchp1         call defdr(ovdif,-nclspp)
         allocate(rv_p_ovdif(abs(-nclspp)))
         if (-nclspp<0) rv_p_ovdif(:)=0.0d0

        call defdr(orhrmx,-nclspp)
Cchp1         call defdr(ovrmax,-2*nclspp)
         allocate(rv_p_ovrmax(abs(-2*nclspp)))
         if (-2*nclspp<0) rv_p_ovrmax(:)=0.0d0

C       call defdr(oqnu,-3*nlspc)
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg         call upack1('array oics',sarray,oics)
        
Cchp1         oics=sarray%oics
         iv_p_oics => sarray%iv_p_oics 


c takao sp2cls is replaced.
c        call sp2cls('spec rmt',sspec,w(oics),1,1,nclspp,ormax)
        allocate( ics(nclspp) )
Cchp1         call icopy(nclspp,w(oics),1,ics,1)
         call icopy ( nclspp , iv_p_oics , 1 , ics , 1 ) 

        print *,' qqqqq ics=',ics
        call defdr(ormax,nclspp)
        call dcopy(nclspp,sspec(ics(1:nclspp))%rmt,1,w(ormax),1)

Cgetarg...info...           integer :: ogroup , integer(8) :: sarray%ogroup 1
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg...info...           integer :: oinitc , integer(8) :: sarray%ohave 1
Cgetarg         call pack5('array ogroup ormax ohave',sarray,ogroup,ormax,
Cgetarg      .    oinitc,0,0)
        
        sarray%ogroup=ogroup 
        sarray%ormax=ormax 
Cchp1         sarray%ohave=oinitc 
         sarray%iv_p_ohave => iv_p_oinitc 


Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: oqc , integer(8) :: spot%oqc 1
Cgetarg...info...           integer :: oqt , integer(8) :: spot%oqt 1
Cgetarg         call pack5('pot opnu oqnu oqc oqt',spot,opnu,oqnu,oqc,oqt,0)
        
        spot%opnu=opnu 
        spot%oqnu=oqnu 
        spot%oqc=oqc 
        spot%oqt=oqt 

Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg         call pack5('pot oves ovdif ovrmax orhrmx opp',spot,oves,ovdif,
Cgetarg      .    ovrmax,orhrmx,opp)
        
Cchp1         spot%oves=oves 
         spot%rv_p_oves => rv_p_oves 

Cchp1         spot%ovdif=ovdif 
         spot%rv_p_ovdif => rv_p_ovdif 

Cchp1         spot%ovrmax=ovrmax 
         spot%rv_p_ovrmax => rv_p_ovrmax 

        spot%orhrmx=orhrmx 
        spot%opp=opp 

C       Even in collinear case, bxc might be written to
        call defdr(obxc,-3*nclasp)
Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg         call pack5('pot opprel obxc nrhos',spot,opprel,obxc,0,0,0)
        
        spot%opprel=opprel 
        spot%obxc=obxc 
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
          call rlse(obxc)
          call defdr(orhos,2*3*4*nrhos*nclasp)
          call defdr(osop,-nlspc*nsp*9)
          call defdr(obxc,-3*nclasp)
Cgetarg...info...           integer :: nrhos , integer(8) :: spot%nrhos 1
Cgetarg...info...           integer :: orhos , integer(8) :: spot%orhos 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg           call pack5('pot nrhos orhos obxc osop',spot,
Cgetarg      .      nrhos,orhos,obxc,osop,0)
          spot%nrhos=nrhos 
          spot%orhos=orhos 
          spot%obxc=obxc 
          spot%osop=osop 
        endif

        ogrrme = 0
Cgetarg         if (igets('ctrl loptc',sctrl) .gt. 0)
Cgetarg      .    call defrr(ogrrme,-16*nl*nsp*nclasp)
        if ( int(sctrl%loptc) .gt. 0 ) call defrr ( ogrrme , -16 * nl 
     .  * nsp * nclasp ) 

Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg         call pack1('pot ogrrme',spot,ogrrme)
        
        spot%ogrrme=ogrrme 

Cgetarg         if (lgors('ctrl lsx,3',sctrl) .or.
Cgetarg      .      lgors('ctrl lscr,-1',sctrl)) then
        if ( iand(3,int(sctrl%lsx)) .ne.0 .or. iand(-1,int(sctrl%lscr)) .ne.0 
     .  ) then 

Cchp1           call defdr(ovintr,-nclasp*(nl*nsp)**2)
           allocate(rv_p_ovintr(abs(-nclasp*(nl*nsp)**2)))
           if (-nclasp*(nl*nsp)**2<0) rv_p_ovintr(:)=0.0d0

Cgetarg...info...           integer :: ovintr , integer(8) :: spot%ovintr 1
Cgetarg           call pack1('pot ovintr',spot,ovintr)
          
Cchp1           spot%ovintr=ovintr 
           spot%rv_p_ovintr => rv_p_ovintr 


        endif
Cgetarg         lnsph = isw(lgors('ctrl lasa,32',sctrl))
c         lnsph = isw ( iand(32,int(sctrl%lasa)) .ne.0 ) 
ctakao remove asa
        lnsph = isw ( .false.)

C   ... Allocate memory for, attempt to read qpp from file
        if (lnsph .ne. 0 .or. nlibu .ne. 0) then
          call defrr(opmpol,-(2*nl-1)*nl**2*3*nsp*nclasp)
          i1 = nl**2
          nqpp = (i1*(i1+1))/2
          call defrr(oqpp,2*nqpp*4*nsp*nbas)
          call dvset(w(oqpp),1,1,-1d0)
Cgetarg...info...           integer :: opmpol , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oqpp , integer(8) :: spot%oqpp 1
Cgetarg           call pack2('pot opmpol oqpp',spot,opmpol,oqpp)
          spot%opmpol=opmpol 
          spot%oqpp=oqpp 
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
        
        opnu=spot%opnu
        oqnu=spot%oqnu

        call dpscop(w(opnu),w(opnu),nl*nsp*i1,
     .  nl*nsp*nclass+1,nl*nsp*nclasp+1,1d0)
        call dpscop(w(oqnu),w(oqnu),3*nl*nsp*i1,
     .  3*nl*nsp*nclass+1,3*nl*nsp*nclasp+1,1d0)

Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg         call upack('pot oves ovdif ovrmax orhrmx',spot,oves,ovdif,
Cgetarg      .    ovrmax,orhrmx,0)
        
Cchp1         oves=spot%oves
         rv_p_oves => spot%rv_p_oves 

Cchp1         ovdif=spot%ovdif
         rv_p_ovdif => spot%rv_p_ovdif 

Cchp1         ovrmax=spot%ovrmax
         rv_p_ovrmax => spot%rv_p_ovrmax 

        orhrmx=spot%orhrmx

Cchp1         call dpscop(w(oves),w(oves),i1,nclass+1,nclasp+1,1d0)
         call dpscop ( rv_p_oves , rv_p_oves , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

Cchp1         call dpscop(w(ovdif),w(ovdif),i1,nclass+1,nclasp+1,1d0)
         call dpscop ( rv_p_ovdif , rv_p_ovdif , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

Cchp1         call dpscop(w(ovrmax),w(ovrmax),i1,nclass+1,nclasp+1,1d0)
         call dpscop ( rv_p_ovrmax , rv_p_ovrmax , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

        call dpscop(w(orhrmx),w(orhrmx),i1,nclass+1,nclasp+1,1d0)

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
      integer w(1)
      common /w/ w
C ... Local parameters
      integer i1,i2,nbas,nl,nsp,igets,ifi,oqpp,fopna,rdm,ipr

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
      oqpp = int(spot%oqpp) 

      ifi = fopna('qpp',-1,4+8)

      if (lio) then
        call ywrm(1,'lmasa',1,ifi,' ',w(oqpp),1,i1,i1,i2)
        if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
        if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
        call pshpr(0)
        if (rdm(ifi,2,i1*i2,' ',w(oqpp),i1,i2) .lt. 0) then
          if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
          call dvset(w(oqpp),1,1,-1d0)
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
      integer w(1)
      common /w/ w
C Local variables
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,lgors,
     .aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,
     .readpq,readpp,readso,readop,readva,readmp
Cchp1       integer opnu,oqnu,orhrmx,ovrmax,opp,opprel,oves,ova,omp,osop,olpp,
Cchp1      .olppr,olsop,olgrme,olva,olmp,opot,ocor,oics,ohave,oclabl,onrc,
Cchp1      .ogrrme,n0,nglob
Cchp1        integer:: opnu , oqnu , orhrmx , ovrmax , opp , opprel , oves 
Cchp1      . , ova , omp , osop , olpp , olppr , olsop , olgrme , olva , olmp 
Cchp1      . , opot , ocor , oics , ohave , onrc , ogrrme , n0 , nglob 
Cchp1         integer:: opnu , oqnu , orhrmx , ovrmax , opp , opprel , oves 
Cchp1      .  , ova , omp , osop , olpp , olppr , olsop , olgrme , olva , olmp 
Cchp1      .  , opot , ocor , ohave , onrc , ogrrme , n0 , nglob 
Cchp1          integer:: opnu , oqnu , orhrmx , ovrmax , opp , opprel , ova 
Cchp1      .   , omp , osop , olpp , olppr , olsop , olgrme , olva , olmp , 
Cchp1      .   opot , ocor , ohave , onrc , ogrrme , n0 , nglob 
Cchp1           integer:: opnu , oqnu , orhrmx , ovrmax , opp , opprel , omp 
Cchp1      .    , osop , olpp , olppr , olsop , olgrme , olva , olmp , opot , 
Cchp1      .    ocor , ohave , onrc , ogrrme , n0 , nglob 
Cchp1            integer:: opnu , oqnu , orhrmx , opp , opprel , omp , osop , 
Cchp1      .     olpp , olppr , olsop , olgrme , olva , olmp , opot , ocor , ohave 
Cchp1      .     , onrc , ogrrme , n0 , nglob 
            integer:: opnu , oqnu , orhrmx , opp , opprel , omp , osop , 
     .      olpp , olppr , olsop , olgrme , olva , olmp , opot , ocor , onrc 
     .      , ogrrme , n0 , nglob 
            integer,pointer :: iv_p_ohave(:)

           real(8),pointer :: rv_p_ovrmax(:)

          real(8),pointer :: rv_p_ova(:)

         real(8),pointer :: rv_p_oves(:)

        integer,pointer :: iv_p_oics(:)

       real(8),pointer :: rv_p_oclabl(:)

      parameter(n0=10)
      integer ic,is,lmx,k,nl,nsp,ifi,jfi,fopn,lmxx,nspx,nrx,nr,
     .idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .nclspp,iprint,lgunit,mode0,mode00,lrel
c      integer mpipid,procid
      double precision rhrmx,vrmax(2),ves,z,rmxx,ax,qc,dq,vrmxx(2),
     .sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,dval,
     .zbak(2),dgets,pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100)

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
      
      opp=spot%opp
      opprel=spot%opprel
      osop=spot%osop
      ogrrme=spot%ogrrme

Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg       call upack('pot opnu oqnu orhrmx ovrmax',spot,opnu,oqnu,orhrmx,
Cgetarg      .  ovrmax,0)
      
      opnu=spot%opnu
      oqnu=spot%oqnu
      orhrmx=spot%orhrmx
Cchp1       ovrmax=spot%ovrmax
       rv_p_ovrmax => spot%rv_p_ovrmax 


Cgetarg...info...           integer :: ova , integer(8) :: spot%ovintr 1
Cgetarg...info...           integer :: omp , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg       call upack('pot ovintr opmpol oves',spot,ova,omp,oves,0,0)
      
Cchp1       ova=spot%ovintr
       rv_p_ova => spot%rv_p_ovintr 

      omp=spot%opmpol
Cchp1       oves=spot%oves
       rv_p_oves => spot%rv_p_oves 


Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: ohave , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg       call upack('array nclasp oclabl ohave oics',sarray,nclasp,oclabl,
Cgetarg      .  ohave,oics,0)
      
      nclasp=sarray%nclasp
Cchp1       oclabl=sarray%oclabl
       rv_p_oclabl => sarray%rv_p_oclabl 

Cchp1       ohave=sarray%ohave
       iv_p_ohave => sarray%iv_p_ohave 

Cchp1       oics=sarray%oics
       iv_p_oics => sarray%iv_p_oics 


      lrel = nglob('lrel')

      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)

C ... MPI: only master does sphere program
c      procid = mpipid(1)
c      if (procid .eq. 0) then

      olppr = 1
      call defdr(olpp,100)
      call defdr(olppr,3000)
      call defdr(olsop,500)
      call defdr(olgrme,500)
      call defdr(olva,100)
      call defdr(olmp,3000)
      call defdr(opot,3000)
      call defdr(ocor,3000)
      i2 = ic2
      if (i2 .eq. 0) i2 = nclasp
c     call awrit2('%n:1i',' ',100,6,nclasp,w(oics))
      do  10  ic = ic1, i2
cki        is = w(oics+ic-1)
        is=iv_p_oics(ic)
Cchp1         icmap = iclbsj(is,w(oics),-nclasp,1)
         icmap = iclbsj ( is , iv_p_oics , - nclasp , 1 ) 

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
Cchp1         call r8tos8(dval(w(oclabl),ic),clabl)
         call r8tos8 ( dval ( rv_p_oclabl , ic ) , clabl ) 

        if (mode00 .ne. 3) then
c#if ONE_ATOM_UNIT
          ifi = fopn(clabl)
c#else
c         ifi = fopna(clabl,30+ic,0)
c#endif
        endif

C  --- Copy what is passed through to holding arrays ---
        k = nl*nsp
Cchp1         call pvaioc(w(ohave),1,ic,icmap,havepq,jc)
         call pvaioc ( iv_p_ohave , 1 , ic , icmap , havepq , jc ) 

        readpq = .not. havepq .or. mode00 .eq. 2 .and. jc .ne. ic
        readpq = readpq .and. mode00 .lt. 3
        if (havepq) then
          call dpscop(w(opnu),pnuloc,k,1+(jc-1)*k,1,1d0)
          call dpscop(w(oqnu),qnuloc,3*k,1+(jc-1)*3*k,1,1d0)
          call dpscop(w(orhrmx),rhrmx,1,jc,1,1d0)
Cchp1           call dpscop(w(ovrmax),vrmax,2,2*jc-1,1,1d0)
           call dpscop ( rv_p_ovrmax , vrmax , 2 , 2 * jc - 1 , 1 , 1d0 
     .     ) 

          if (ic.ne.jc) call awrit0('%a pq,',outs1,len(outs1),0)
        endif
        k = 6*nl*nsp
Cchp1         call pvaioc(w(ohave),2,ic,icmap,havepp,jc)
         call pvaioc ( iv_p_ohave , 2 , ic , icmap , havepp , jc ) 

        readpp = .not. havepp .or. mode00 .eq. 2 .and. jc .ne. ic
        readpp = readpp .and. mode00 .lt. 3
        if (havepp) then
          call dpscop(w(opp),w(olpp),k,1+(jc-1)*k,1,1d0)
          if (lrel .eq. 2) then
            k = 4*nl*2*nl*2*2
            call dpscop(w(opprel),w(olppr),k,1+(jc-1)*k,1,1d0)
          endif
Cchp1           call dpscop(w(oves),ves,1,jc,1,1d0)
           call dpscop ( rv_p_oves , ves , 1 , jc , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a pp,',outs1,len(outs1),0)
        endif
        k = (nl*nsp)**2
Cchp1         call pvaioc(w(ohave),8,ic,icmap,haveva,jc)
         call pvaioc ( iv_p_ohave , 8 , ic , icmap , haveva , jc ) 

cki        haveva = haveva .and. ova .gt. 0
        haveva = haveva .and. associated(rv_p_ova)
        readva = .not. haveva .or. mode00 .eq. 2 .and. jc .ne. ic
cki        readva = readva .and. ova .gt. 0 .and. mode00 .lt. 3
        readva = readva .and. associated(rv_p_ova) .and. mode00 .lt. 3
        if (haveva) then
Cchp1           call dpscop(w(ova),w(olva),k,1+(jc-1)*k,1,1d0)
           call dpscop ( rv_p_ova , w ( olva ) , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a va,',outs1,len(outs1),0)
        endif
        k = nl**2*(2*nl-1)*3*nsp
Cchp1         call pvaioc(w(ohave),16,ic,icmap,havemp,jc)
         call pvaioc ( iv_p_ohave , 16 , ic , icmap , havemp , jc ) 

        havemp = havemp .and. omp .gt. 0
        readmp = .not. havemp .or. mode00 .eq. 2 .and. jc .ne. ic
        readmp = readmp .and. omp .gt. 0 .and. mode00 .lt. 3
        if (havemp) then
          call dpscop(w(omp),w(olmp),k,1+(jc-1)*k,1,1d0)
          if (ic.ne.jc) call awrit0('%a mp,',outs1,len(outs1),0)
        endif
        k = nl*nsp*nsp*9
Cchp1         call pvaioc(w(ohave),4,ic,icmap,haveso,jc)
         call pvaioc ( iv_p_ohave , 4 , ic , icmap , haveso , jc ) 

        haveso = haveso .and. osop.gt. 0
        readso = .not. haveso .or. mode00 .eq. 2 .and. jc .ne. ic
        readso = readso .and. osop .gt. 0 .and. mode00 .lt. 3
        if (haveso) then
          call dpscop(w(osop),w(olsop),k,1+(jc-1)*k,1,1d0)
          if (ic.ne.jc) call awrit0('%a sop,',outs1,len(outs1),0)
        endif
        k = 16*nl*nsp
Cchp1         call pvaioc(w(ohave),32,ic,icmap,haveop,jc)
         call pvaioc ( iv_p_ohave , 32 , ic , icmap , haveop , jc ) 

        haveop = haveop .and. ogrrme .gt. 0
        readop = .not. haveop .or. mode00 .eq. 2 .and. jc .ne. ic
        readop = readop .and. ogrrme .gt. 0 .and. mode00 .lt. 3
        if (haveop) then
          call dpscop(w(ogrrme),w(olgrme),k,1+(jc-1)*k,1,1d0)
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
            lpot = aiopot(nr,nsp,a,rmax,-99d0,w(opot),ifi)
          endif
          lcor = aiocor(nr,nsp,a,rmxx,w(ocor),sumec,sumtc,ifi)

          rewind ifi
          jfi = -ifi
          if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .    dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
          if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,idmod,
     .    nl,lmx,nsp,rhrmx,vrmax,jfi)
          if (havepp) sw = aiopar(clabl,lrel,w(olpp),w(olppr),ves,nl,
     .    lmx,nsp,jfi)
          if (haveva) sw = aiova(clabl,w(olva),nl,lmx,nsp,jfi)
          if (havemp) sw = aiomp(clabl,w(olmp),nl,2*nl-2,nsp,jfi)
          if (haveso) sw = aiosop(clabl,w(olsop),nl,lmx,nsp,jfi)
          if (haveop) sw = aiorme(clabl,w(olgrme),nl,nsp,jfi)
          if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,w(opot),jfi)
          if (lcor) lcor = aiocor(nr,nsp,a,rmxx,w(ocor),sumec,sumtc,jfi)

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
            readpp = aiopar(clabl,lrel,w(olpp),w(olppr),ves,nl,lmx,nsp,
     .      ifi)
            if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
          endif
          if (readso .and. osop .gt. 0) then
            readso = aiosop(clabl,w(olsop),nl,lmx,nsp,ifi)
            if (readso) call awrit0('%a so,',outs2,len(outs2),0)
          endif
          if (readop .and. ogrrme .gt. 0) then
            readop = aiorme(clabl,w(olgrme),nl,nsp,ifi)
            if (readop) call awrit0('%a op,',outs2,len(outs2),0)
          endif
cki          if (readva .and. ova .gt. 0) then
          if (readva .and. associated(rv_p_ova) ) then
            readva = aiova(clabl,w(olva),nl,lmx,nsp,ifi)
            if (readva) call awrit0('%a va,',outs2,len(outs2),0)
          endif
          if (readmp .and. omp .gt. 0) then
            readmp = aiomp(clabl,w(olmp),nl,2*nl-2,nsp,ifi)
            if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
          endif

C     ... Update what parameters are available
cki          w(ohave+ic-1) = isw(haveop.or.readop)*32+
          iv_p_ohave(ic) = isw(haveop.or.readop)*32+
     .    isw(havemp.or.readmp)*16+
     .    isw(haveva.or.readva)*8 +
     .    isw(haveso.or.readso)*4 +
     .    isw(havepp.or.readpp)*2 +
     .    isw(havepq.or.readpq)*1 +
cki     .    w(ohave+ic-1) - bitand(w(ohave+ic-1),63)
     .    iv_p_ohave(ic) - bitand(iv_p_ohave(ic),63)

          k = nl*nsp
          if (havepq .or. readpq) then
            call dpscop(pnuloc,w(opnu),k,1,1+(ic-1)*k,1d0)
            call dpscop(qnuloc,w(oqnu),3*k,1,1+(ic-1)*3*k,1d0)
            call dpscop(rhrmx,w(orhrmx),1,1,ic,1d0)
Cchp1             call dpscop(vrmax,w(ovrmax),2,1,2*ic-1,1d0)
             call dpscop ( vrmax , rv_p_ovrmax , 2 , 1 , 2 * ic - 1 , 1d0 
     .       ) 

          endif
          k = 6*nl*nsp
          if (havepp .or. readpp) then
            call dpscop(w(olpp),w(opp),k,1,1+(ic-1)*k,1d0)
            if (lrel .eq. 2) then
              k = 4*nl*2*nl*2*2
              call dpscop(w(olppr),w(opprel),k,1,1+(ic-1)*k,1d0)
            endif
Cchp1             call dpscop(ves,w(oves),1,1,ic,1d0)
             call dpscop ( ves , rv_p_oves , 1 , 1 , ic , 1d0 ) 

          endif
          k = (nl*nsp)**2
cki          if ((haveva .or. readva) .and. ova .gt. 0)
Cchp1           if ((haveva .or. readva) .and. associated(rv_p_ova) ) 
Cchp1      .    call dpscop(w(olva),w(ova),k,1,1+(ic-1)*k,1d0)
           if ( ( haveva .or. readva ) .and. associated ( rv_p_ova ) ) call 
     .     dpscop ( w ( olva ) , rv_p_ova , k , 1 , 1 + ( ic - 1 ) * k , 
     .     1d0 ) 

          k = nl**2*(2*nl-1)*3*nsp
          if ((havemp .or. readmp) .and. omp .gt. 0)
     .    call dpscop(w(olmp),w(omp),k,1,1+(ic-1)*k,1d0)
          k = nl*nsp*nsp*9
          if ((haveso .or. readso) .and. osop .gt. 0)
     .    call dpscop(w(olsop),w(osop),k,1,1+(ic-1)*k,1d0)
          k = 16*nl*nsp
          if ((haveop .or. readop) .and. ogrrme .gt. 0)
     .    call dpscop(w(olgrme),w(ogrrme),k,1,1+(ic-1)*k,1d0)
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
      call rlse(olpp)
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


