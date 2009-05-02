      program kino_input_test
c--- This is a sample program to explain how to use keyvalue.f by Dr.H.kino.
c This program works with IN file which can be
c ======from here ===============================================
c int2  200
c realxxx -100.0
c intv  2 3 4 5 6
c realvv -9.0 -8.0 -2
c char2 sldkjflsdjfljsd
c <block2>
c re -----------------------------------
c test test ok
c t s e t
c 1
c 2
c abcd fgw kij
c </block2>
c ======to here ================================================
c Anyway, you can see how to use it by chaning this IN file and codes here.
      use keyvalue
      implicit none

      integer:: ikey,iflag,ivkey(10),file,i,ret
      real(8) :: rkey,rvkey(10),rvkey2(10)
      logical :: lflag,lkey

      character(200):: ckey,buf
c        call input_substituteKeys('%d%s/#vnew.data',2,(/'%d','%s'/),
c     i       (/'/home/kino','psnmame   '/),ret)

      integer(4):: nmbas=300, imbas(300) , imbas0(300),istat
      istat=-9999
      call getkeyvalue("GWinput","MagAtom", 
     &        imbas,nmbas,status=istat)
      write(6,*) ' nmbas istat=',nmbas,istat
      stop

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
c        rvkey2(1:3)=(/10.0d0,10.0d0,10.0d0/)
      call getkeyvalue("IN","realvv", rvkey,3,status=ret,
     &                 default=(/10.0d0,10.0d0,10.0d0,10d0,10d0,10d0/))
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

      stop
      end program

