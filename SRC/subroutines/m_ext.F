      module m_ext
      character(512),public,protected::sname='???'
      public:: m_ext_init
      contains
      subroutine m_ext_init()
      logical master
      integer:: ifi,ipos,i,na
      character*256:: sss,s222,argv
!! Get the extension of ctrlfile. It is stored to trim(sname).
      na= iargc()
      do i = 1, na
         call getarg( i, argv )
c         write( 6, '( i5, x, a )' ) i, trim(argv)
         if(argv(1:1)/='-') then
            sname=trim(argv)
            exit
         endif
      enddo
      end subroutine
      end module
      
! use i2char instead !-----------------      
c$$$      character(9) function itoa(numin)
c$$$      integer(4) ::numin,num
c$$$      num=abs(numin)
c$$$                     itoa= char(48+mod(num,10))
c$$$      if(num>9)      itoa= char(48+mod(num/10,10))//itoa
c$$$      if(num>99)     itoa= char(48+mod(num/100,10))//itoa
c$$$      if(num>999)    itoa= char(48+mod(num/1000,10))//itoa
c$$$      if(num>9999)   itoa= char(48+mod(num/10000,10))//itoa
c$$$      if(num>99999)  itoa= char(48+mod(num/100000,10))//itoa
c$$$      if(num>999999) itoa= char(48+mod(num/1000000,10))//itoa
c$$$      if(num>9999999)itoa= char(48+mod(num/10000000,10))//itoa
c$$$      if(num>99999999) call rx( ' charext:can not produce')
c$$$      if(numin<0) itoa='-'//itoa
c$$$      end
c
c$$$      character(5) function charnum5(num)
c$$$      integer(4) ::num
c$$$      charnum5 = char(48+mod(num/10000,10))//
c$$$     &           char(48+mod(num/1000,10))//
c$$$     &           char(48+mod(num/100,10))//
c$$$     &           char(48+mod(num/10,10))//
c$$$     &           char(48+mod(num,10))
c$$$      end
      
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
      

      character*9 function ftoa9(arg)
      real(8):: arg
      write(ftoa9,"(1x,f8.3)") arg
      end

c-------------
      character*8 function xt(num)
      integer :: num
      if(num==0) then
        xt=''
        return
      endif
      xt = char(48+mod(num,10))
      if(num>9)     xt = char(48+mod(num/10,10))//xt
      if(num>99)    xt = char(48+mod(num/100,10))//xt
      if(num>999)   xt = char(48+mod(num/1000,10))//xt
      if(num>9999)  xt = char(48+mod(num/10000,10))//xt
      if(num>99999) xt = char(48+mod(num/100000,10))//xt
      if(num>999999) call rx( ' xt:can not produce')
      xt='.'//xt
      end
c-------------
      character*20 function xxt(num1,num2)
      integer :: num,num1,num2
      num = num2
      xxt = char(48+mod(num,10))
      if(num>9)     xxt = char(48+mod(num/10,10))//xxt
      if(num>99)    xxt = char(48+mod(num/100,10))//xxt
      if(num>999)   xxt = char(48+mod(num/1000,10))//xxt
      if(num>9999)  xxt = char(48+mod(num/10000,10))//xxt
      if(num>99999) xxt = char(48+mod(num/100000,10))//xxt
      if(num>999999) call rx( ' xxt:can not produce')
      xxt ='to'//xxt
      num = num1
      xxt = char(48+mod(num,10))//xxt
      if(num>9)     xxt = char(48+mod(num/10,10))//xxt
      if(num>99)    xxt = char(48+mod(num/100,10))//xxt
      if(num>999)   xxt = char(48+mod(num/1000,10))//xxt
      if(num>9999)  xxt = char(48+mod(num/10000,10))//xxt
      if(num>99999) xxt = char(48+mod(num/100000,10))//xxt
      if(num>999999) call rx( ' xxt:can not produce')
      end


      
      
