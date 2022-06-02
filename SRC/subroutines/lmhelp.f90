subroutine lmhelp(prgnam)
  use m_lgunit,only:stdo
  implicit none
  character prgnam*8
  integer :: i1,i2
  character ch*1,outs*1000
  call locase(prgnam)
  call info0(0,0,0,' usage:  '//prgnam//'%a [--OPTION] [-var-assign] [ext]')
  print 343
343 format(/' --help',t17, &
       'List categories, tokens, and data program expects, and quit' &
       /' --show',t17, &
       'Print control file after parsing by preprocessor,'/t17, &
       'and echo input data as read from the control file' &
       /' --pr=#1',t17, &
       'Set the verbosity (stack) to values #1' &
       /' --time=#1[,#2]',t17, &
       'Print timing info to # levels (#1=summary; #2=on-the-fly)'/ &
       /' -vnam=expr',t17, &
       'Define numerical variable "nam"; set to result of ''expr''')
  outs = '%N '//prgnam//'%a-specific options:'
  call strip(outs,i1,i2)
  call info0(0,0,0,outs(1:i2))
  if (prgnam == 'lmfa') then
     call info0(0,0,0,  '%N%1f ')
  endif
  if (prgnam == 'lmfgwd') then
     call info0(0,0,0, '%N%1f ')
  endif
  if (prgnam == 'lmf') then
     write(6,*)'--rs=#1,#2,#3,#4,#5'
     write(6,*)'#1: =1 from rst, =2 from rsta'
     write(6,*)'#2: =1 save to rst, =2 save to rsta'
     write(6,*)'#3: =1 read pos from ctrl even when we have rst'
     write(6,*)'#4: obsolate'
     write(6,*)'#5: =1 read P from ctrl even when we have rst'
     write(6,*)'lmf-MPIK --jobgw works as GW driver'
     write(stdo,"(/a)") ' Other command-line-options:  Search "call cmdopt" in *.F'
  endif
  call fexit0(0,' ')
end subroutine lmhelp
