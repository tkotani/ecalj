      module m_ext
      implicit none
      character(512),public::sname
      public:: m_ext_init
      contains
      subroutine m_ext_init()
      logical master
      integer:: ifi,ipos,fopn,i,na
      character*256:: sss,s222,argv
!! Get the extension of ctrlfile. It is stored to trim(sname).
c      if(initex) then
      na= iargc()
      do i = 1, na
         call getarg( i, argv )
c         write( 6, '( i5, x, a )' ) i, trim(argv)
         if(argv(1:1)/='-') then
            sname=trim(argv)
            exit
         endif
      enddo
c      print *,'file ext=',trim(sname)
c      stop 'xxxxxxxxxxxxxxxx'
c        ifi=fopn('temptempxxx') 
c        sss=''
c        INQUIRE (ifi, NAME=sss) !second run, 'sss' must contain extention, eg, '.copt'
c        print *,'sss=',sss
c        call fclose(ifi)
c        ipos= scan(sss,'.',back=.true.) !backward search of '.' from the last.
c        sname=sss(ipos+1:256)  !sname contains extention such as 'copt'.
c        open(newunit=ifi, file='temptempxxx.'//trim(sname)) 
c        close(ifi,status='delete')
c        initex=.false.
c      endif
      end subroutine
      end module

      character(5) function charnum5(num)
      integer(4) ::num
      charnum5 = char(48+mod(num/10000,10))//
     &           char(48+mod(num/1000,10))//
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
      
      character(4) function charnum4(num)
      integer(4) ::num
      charnum4 = 
     &           char(48+mod(num/1000,10))//
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
c
      character*20 function xts(num1,num2)
      integer(4) :: num,num1,num2
c      character*(*) add
      num = num2
      xts = char(48+mod(num,10))
      if(num>9)     xts = char(48+mod(num/10,10))//xts
      if(num>99)    xts = char(48+mod(num/100,10))//xts
      if(num>999)   xts = char(48+mod(num/1000,10))//xts
      if(num>9999)  xts = char(48+mod(num/10000,10))//xts
      if(num>99999) xts = char(48+mod(num/100000,10))//xts
      if(num>999999) call rx( ' xts:can not produce')
      xts ='.L'//xts
      num = num1
      xts = char(48+mod(num,10))//xts
      if(num>9)     xts = char(48+mod(num/10,10))//xts
      if(num>99)    xts = char(48+mod(num/100,10))//xts
      if(num>999)   xts = char(48+mod(num/1000,10))//xts
      if(num>9999)  xts = char(48+mod(num/10000,10))//xts
      if(num>99999) xts = char(48+mod(num/100000,10))//xts
      if(num>999999) call rx( ' xts:can not produce')
      end
      
      character*8 function charext(num)
      integer(4) ::num
      charext = char(48+mod(num,10))
      if(num>9)   charext= char(48+mod(num/10,10))//charext
      if(num>99)  charext= char(48+mod(num/100,10))//charext
      if(num>999) charext= char(48+mod(num/1000,10))//charext
      if(num>9999)charext= char(48+mod(num/10000,10))//charext
      if(num >99999) call rx( ' charext:can not produce')
      end

      character*9 function ftoa9(arg)
      real(8):: arg
      write(ftoa9,"(1x,f8.3)") arg
      end
