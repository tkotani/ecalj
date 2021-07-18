      subroutine lmhelp(prgnam)
      implicit none
      character prgnam*8
      integer i1,i2,stdo,lgunit
      character ch*1,outs*1000
      stdo = lgunit(1)
      call locase(prgnam)
      call info0(0,0,0,' usage:  '//prgnam//'%a [--OPTION] [-var-assign] [ext]')
      print 343
      print 344
  343 format(/' --help',t17,
     .'List categories, tokens, and data program expects, and quit'
     ./' --show',t17,
     .'Print control file after parsing by preprocessor,'/t17,
     .'and echo input data as read from the control file'
     ./' --showp',t17,
     .'Same as --show, but quit after input parsed'
     ./' --pr#1[,#2...]',t17,
     .'Set the verbosity (stack) to values #1,#2, ...'
     ./' --time=#1[,#2]',t17,
     .'Print timing info to # levels (#1=summary; #2=on-the-fly)'/
     ./' -vnam=expr',t17,
     .'Define numerical variable "nam"; set to result of ''expr''')
  344 format(
     ./' --rpos=filnam',t17,
     .'After reading input file, read site positions from "filnam" ==>NEED CHECK!'/
     .' --fixpos[:tol=#]',
     .' Adjust positions slightly, rendering them as'/
     .t17,'exactly as possible consistent with the symmetry group ==>NEED CHECK!')
      if (.true.) then
        outs = '%N '//prgnam//'%a-specific options:'
        call strip(outs,i1,i2)
        call info0(0,0,0,outs(1:i2))
        if (prgnam .eq. 'lmfa') then
          call info0(0,0,0,  '%N%1f ')
        endif
        if (prgnam .eq. 'lmfgwd') then
          call info0(0,0,0, '%N%1f ')
        endif
        if (prgnam .eq. 'lmf') then
          call info0(0,0,0,
     .    '%N%1f--rs=#1,#2,#3,#4,#5'//
     .    '%N%6f#1=0 start from atm file; 1 from rst file;'//
     .    ' 2 from rsta file'//
     .    '%N%11fadd 10 to shift sm-rho 1st iter'//
     .    '%N%11fadd 100 to rotate local rho 1st iter'//
     .    '%N%6f#2=1 save rst file'//
     .          '%N%6f(#3,#4,#5)=0 read (pos,E_f,pnu) from rst file')
          write(stdo,"(/a)") ' Other command-line-options:  Search "call cmdopt" in *.F'
        endif
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
      call fexit(0,0,' ',0)
      end subroutine lmhelp
