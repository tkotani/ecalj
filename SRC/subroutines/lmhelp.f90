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
  if (prgnam == 'lmf') then
     write(stdo,*)'=== lmf-MPIK options=== (No --rs option now: 2022-6-20)'
     write(stdo,*)' --jobgw works as GW driver'
     write(stdo,*)' -v{foobar}=xxx replace {foobar} with xxx'
     write(stdo,*)' Read rst.* prior to atm.* file'
     write(stdo,*)' Use --quit=band, --quit=mkpot or --quit=dmat to surpress writing rst'
     write(stdo,*) ' Other command-line-options:  Search "call cmdopt" in *.F'
  endif
  call fexit0(0,' ')
end subroutine lmhelp
