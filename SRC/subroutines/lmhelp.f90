subroutine lmhelp(prgnam)
  use m_lgunit,only:stdo
  implicit none
  character prgnam*8
  integer :: i1,i2
  character ch*1,outs*1000
  call locase(prgnam)
  write(stdo,*)' usage:  '//trim(prgnam)//' [--OPTION] [-var-assign] [extension]'
  print 343
343 format(/' --help',t17, &
         'List categories, tokens, and data program expects, and quit' &
         /' --show',t17, 'Print control file after parsing by preprocessor,'/t17, &
         'and echo input data as read from the control file' &
         /' --pr=#1',t17, 'Set the verbosity (stack) to values #1' &
         /' --time=#1[,#2]',t17, 'Print timing info to # levels (#1=summary; #2=on-the-fly)'/ &
         /' -vnam=expr',t17,'Define numerical variable "nam"; set to result of ''expr''')
    write(stdo,*)' --jobgw       lmf-MPIK works as the GW driver (previous lmfgw-MPIK)'
    write(stdo,*)' --quit=band, --quit=mkpot or --quit=dmat: Stop points. Surpress writing rst'
    write(stdo,*)
    write(stdo,*)' NOTE: Read rst.* prior to atm.* file (No --rs options: 2022-6-20)'
    write(stdo,*)' NOTE: Other command-line-options => Search "call cmdopt" in SRC/*/*.f90'
    call fexit0(0,' ')
end subroutine lmhelp
