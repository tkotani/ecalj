      module m_ext
      implicit none
      character(512)::sname
      logical,private:: initex=.true.
      contains
      subroutine getext(master)
      logical master
      integer:: ifi,ipos,fopn
      character*256:: sss,s222
!! A trick to get the extension of ctrlfile. It is stored to trim(sname).
      if(initex) then
        ifi=fopn('temptempxxx') 
        sss=''
        INQUIRE (ifi, NAME=sss) !second run, 'sss' must contain extention, eg, '.copt'
        print *,'sss=',sss
        call fclose(ifi)
        ipos= scan(sss,'.',back=.true.) !backward search of '.' from the last.
        sname=sss(ipos+1:256)  !sname contains extention such as 'copt'.
        open(newunit=ifi, file='temptempxxx.'//trim(sname)) 
        close(ifi,status='delete')
        initex=.false.
      endif
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
