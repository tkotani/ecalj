


      subroutine lmhelp(prgnam) !,vn,vn2,wksize)
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
c      double precision vn,vn2
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

        if (prgnam .eq. 'lmfa') then
          call info0(0,0,0,
     *    '%N%1f ')
        endif

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

c         if (prgnam .eq. 'lmdos') then
c            call info0(0,0,0,
c     .           '%N%1f--dos:options  modifies number and kinds of dos '//
c     .           'generated;%N%16fsee documentation (doc/lmto.html)')
c         endif

         if (prgnam .eq. 'lmchk') then
            call info0(0,0,0,
     .           '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
     .           '[:pairs:pair-list]...'//
     .           '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
     .           '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
     .           '%N --findes'//
     .           '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
         endif
      endif

c      if (nint(vn2) .ne. 0) then
c        ch = char(nint(vn2)+ichar('a')-1)
c      else
c        ch = ' '
c      endif
c      call info2(0,0,0,'%N version %,2d'//ch//'   worksize %i K',
c     .vn-vn2/1d5,wksize/1000)
      call fexit(0,0,' ',0)
      end subroutine lmhelp


c$$$      subroutine sttmpd
c$$$C- Creates special tmpdir for saving temporary files
c$$$C  User may wish to customize this routine.
c$$$      character tmpdir*100
c$$$      integer fopnT,ifi
c$$$C ... for henry, lm-MPIK
c$$$C      integer i1,i2,nw
c$$$C      character*40 strn
c$$$
c$$$C     return
c$$$C ... Set customization of temporary directory here, if desired
c$$$C     This is usual default (current working directory)
c$$$      tmpdir = ' '
c$$$C     call gtenv('HOME',tmpdir)
c$$$C     call gtenv('TMPDIR',tmpdir)
c$$$C     call getenv('HOME',tmpdir)
c$$$C     call getenv('TMPDIR',tmpdir)
c$$$C     tmpdir = '/home/tmp/'
c$$$
c$$$C ... for henry, lm-MPIK specific
c$$$C      call getenv('HOME',strn)
c$$$C      call strip(strn,i1,i2)
c$$$C      call wrdsg(strn(i1:i2),0,'/',nw)
c$$$C      call wordg(strn,0,'/',nw,i1,i2)
c$$$C      strn = strn(i1:)
c$$$C      if (strn .eq. 'markv') strn = 'svan2'
c$$$C      call word(strn,1,i1,i2)
c$$$C      tmpdir = '/home/' // strn(i1:i2) // '/tmp'
c$$$
c$$$C ... Set the directory
c$$$      ifi = fopnT(tmpdir,0,0,11)
c$$$
c$$$C     debugging check
c$$$C      ifi = fopnT('tmp' ,-1,0,0)
c$$$C      call fshow
c$$$C      print *, ifi
c$$$C      write(ifi,*) 'test'
c$$$C      call rx('done')
c$$$Cgetarg       end
c$$$      end subroutine sttmpd



      subroutine findctrlstart(nfilin)
! if we find 'ctrlstart', locate reading at the next line of ctrlstart.
!     this is useful if you like to use script, GWinput, ctrl in a file.
      implicit none
      integer:: nfilin
      character(len=9):: strn, ccc
      do
        read(nfilin,"(a)",err=1010,end=1010) strn
        if(strn == 'ctrlstart') return
      enddo
 1010 continue
      rewind(nfilin)
      end subroutine findctrlstart


