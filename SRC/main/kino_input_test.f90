program kino_input_test
  !--- This is a sample program to explain how to use keyvalue.f by Dr.H.kino.
  ! This program works with IN file which can be
  ! ======from here ===============================================
  ! int2  200
  ! realxxx -100.0
  ! intv  2 3 4 5 6
  ! realvv -9.0 -8.0 -2
  ! char2 sldkjflsdjfljsd
  ! <block2>
  ! re -----------------------------------
  ! test test ok
  ! t s e t
  ! 1
  ! 2
  ! abcd fgw kij
  ! </block2>
  ! ======to here ================================================
  ! Anyway, you can see how to use it by chaning this IN file and codes here.
  use m_keyvalue,only: getkeyvalue

  implicit none

  integer:: ikey,iflag,ivkey(10),file,i,ret
  real(8) :: rkey,rvkey(10),rvkey2(10)
  logical :: lflag,lkey

  character(200):: ckey,buf
  !        call input_substituteKeys('%d%s/#vnew.data',2,(/'%d','%s'/),
  !     i       (/'/home/kino','psnmame   '/),ret)

  integer(4):: nmbas=300, imbas(300) , imbas0(300),istat
  istat=-9999
  call getkeyvalue("GWinput","MagAtom", &
       imbas,nmbas,status=istat)
  write(6,*) ' nmbas istat=',nmbas,istat
  ! top2rx 2013.08.09 kino      stop
  call rx( '')

  write(6,*)' =============test1======================='
  call getkeyvalue("IN","int",ikey,default=999,status=ret )
  write(*,*) ret,ikey

  write(6,*)' =============test2======================='
  call getkeyvalue("IN","real",rkey,default=1430d0,status=ret)
  write(*,*) ret,rkey

  write(6,*)' =============test3======================='
  call getkeyvalue("IN","logical",lkey,status=ret )
  write(*,*) ret,lkey

  write(6,*)' =============test4======================='
  call getkeyvalue("IN","intv",ivkey,5,status=ret )
  write(*,*) ret,ivkey

  write(6,*)' =============test5======================='
  !        rvkey2(1:3)=(/10.0d0,10.0d0,10.0d0/)
  call getkeyvalue("IN","realvv", rvkey,3,status=ret, &
       default=(/10.0d0,10.0d0,10.0d0,10d0,10d0,10d0/))
  write(*,*) ret,rvkey

  write(6,*)' =============test6======================='
  call getkeyvalue("IN","<block>",unit=file,status=ret)
  write(*,*) ret,file
  do i=1,ret
     read(file,'(a)') buf
     write(*,'(a)') buf
  enddo
  close(file)

  write(6,*)' =============test7======================='
  call getkeyvalue("IN","char",ckey,status=ret)
  write(*,*) ret,ckey(:len_trim(ckey))

  ! top2rx 2013.08.09 kino      stop
  call rx( '')
end program kino_input_test

