      subroutine lmhelp(prgnam)
      implicit none
      character prgnam*8
      integer i1,i2,stdo,lgunit
      character ch*1,outs*1000
      stdo = lgunit(1)
      call locase(prgnam)
      call info0(0,0,0,' usage:  '//prgnam//'%a [--OPTION] [-var-assign] [ext]')
      print 343
c      print 344
  343 format(/' --help',t17,
     .'List categories, tokens, and data program expects, and quit'
     ./' --show',t17,
     .'Print control file after parsing by preprocessor,'/t17,
     .'and echo input data as read from the control file'
c     ./' --showp',t17,
c     .'Same as --show, but quit after input parsed'
     ./' --pr=#1',t17,
     .'Set the verbosity (stack) to values #1'
     ./' --time=#1[,#2]',t17,
     .'Print timing info to # levels (#1=summary; #2=on-the-fly)'/
     ./' -vnam=expr',t17,
     .'Define numerical variable "nam"; set to result of ''expr''')
c  344 format(
c     ./' --rpos=filnam',t17,
c     .'After reading input file, read site positions from "filnam" ==>NEED CHECK!'/
c     .' --fixpos[:tol=#]',
c     .' Adjust positions slightly, rendering them as'/
c     .t17,'exactly as possible consistent with the symmetry group ==>NEED CHECK!')
c      if (.true.) then
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
           write(6,*)'--rs=#1,#2,#3,#4,#5'
           write(6,*)'#1: =1 from rst, =2 from rsta'
           write(6,*)'#2: =1 save to rst, =2 save to rsta'
           write(6,*)'#3: =1 read pos from ctrl even when we have rst'
           write(6,*)'#4: obsolate'
           write(6,*)'#5: =1 read P from ctrl even when we have rst'
c     .    '%N%11fadd 10 to shift sm-rho 1st iter'//
c     .    '%N%11fadd 100 to rotate local rho 1st iter'//
c     .    '%N%6f#2=1 save rst file 2 save rsta'//
c     .          '%N%6f(#3,#4,#5)=0 read (pos,E_f,pnu) from rst file')
          write(stdo,"(/a)") ' Other command-line-options:  Search "call cmdopt" in *.F'
        endif
c$$$        if (prgnam .eq. 'lmchk') then
c$$$            call info0(0,0,0,
c$$$     .           '%N%1f--shell[:v][:e][:r=#][:sites:site-list]'//
c$$$     .           '[:pairs:pair-list]...'//
c$$$     .           '%N%8f...[:tab[=#]][:disp=fnam][:nn][:fn=fnam]'//
c$$$     .           '%N --mino[:dxmx=#][:xtol=#][:maxit=#][:style=#]:list'//
c$$$     .           '%N --findes'//
c$$$     .           '%N --wpos=fnam'//'%N --angles[=#]'//'%N --terse')
c$$$        endif
c      endif
      call fexit0(0,' ')
      end subroutine lmhelp
