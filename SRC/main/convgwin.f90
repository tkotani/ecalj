program convgwin
  !- convert GWIN to GWIN_V2
  !r takao kotani Apr 2002
  ! ----------------------------------------------------------------------
  implicit integer (i-n)
  character(len=100)::record
  !      integer:: ifile_hanlde
  !      ifigwin  = ifile_handle()
  !      ifigwino = ifile_handle()
  open (newunit=ifigwino,file='GWIN')
  open (newunit=ifigwin, file='GWIN_V2')

  do i=1,20
     read(ifigwino,"(a)") record
     write(ifigwin,"(a)") record
  enddo

  do
     read(ifigwino,*,err=110)i1,i2,i3,i4
     i3=2
     write(ifigwin,"(4i5)")   i1,i2,i3,i4
  enddo
110 continue
  !      write(6,*)'1111111111111'
  write(ifigwin,"(a)") &
       "  atom   l    n  occ;unocc;Valence(1=yes,0=no) "
  do
     read(ifigwino,*,err=120)i1,i2,i3,i4,i5
     write(ifigwin,"(5i5)")   i1,i2,i3,i4,i5
     i3=2;i4=0;i5=0
     if(i4==2) i4=1
     if(i5==2) i5=1
     write(ifigwin,"(5i5)")   i1,i2,i3,i4,i5
  enddo
120 continue
  !      write(6,*)'22222222222222'
  write(ifigwin,"(a)") &
       '  atom   l    n  occ unocc   '// &
       'ForX0 ForSxc :CoreState(1=yes, 0=no)'
  do
     read(ifigwino,"(a)",end=130,err=130)record
     write(ifigwin,"(a)")record
  enddo
130 continue
  write(6,*)' GWIN is converted to GWIN_V2!'
END PROGRAM convgwin
