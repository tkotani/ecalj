      integer function ifile_handle()
!! find unused file handle
      implicit none
      integer:: i
      logical:: nexist
      integer,save:: irem=2001
      character*256::nnn
      ifile_handle=-999999
      do i=irem,9999
         inquire(unit=i,opened=nexist,name=nnn)
         if(.not.nexist) then
            ifile_handle=i
            irem=i+1
            return
         endif
c         print *,'mpipid i=',mpipid(1),i,trim(nnn)
      enddo
      do i=5001,irem
         inquire(unit=i,opened=nexist)
         if(.not.nexist) then
            ifile_handle=i
            irem=i
            return
         endif
      enddo
      call rx('ifile_handle: we did not find open file handle')
      end
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine rangedq(qin, qout) !renewed jun2015
      implicit none! qout is in [-0.5d0,0.5d0)
      intent(in)::       qin
      intent(out)::           qout
      integer :: ix
      real(8):: qin(3),qout(3),tol=1d-12 
      qout= qin-nint(qin)
      do ix=1,3
        if(qout(ix)>0.5d0-tol) qout(ix)=-0.5d0 !this is needed to distinguish 0.5d0 and -0.5d0.
      enddo 
      end
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      character(10) function i2char(numin)
!! convert num to char. See charnum3 to understand this.
      implicit none
      integer(4) ::num,itens,iten,numin
      num=abs(numin)
      if(num>=10**9) call rx('i2char:num>10**9')
      i2char=''
      do itens=8,1,-1 !itens mean itens+1's digit
        iten=10**itens
        if(num>=iten) i2char=trim(i2char)//char(48+mod(num/iten,10))
      enddo
      i2char = trim(i2char)//char(48+mod(num,10)) !1st digit
      if(numin<0) i2char = '-'//trim(i2char)//char(48+mod(num,10)) !1st digit
      end
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      character*8 function charext(num)
      integer(4) ::num
      charext = char(48+mod(num,10))
      if(num>9)   charext= char(48+mod(num/10,10))//charext
      if(num>99)  charext= char(48+mod(num/100,10))//charext
      if(num>999) charext= char(48+mod(num/1000,10))//charext
      if(num>9999)charext= char(48+mod(num/10000,10))//charext
      if(num >99999) call rx( ' charext:can not produce')
      end
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      integer function ichangesign(a,n)
      implicit none
      integer:: i,n
      real(8):: a(n)
      ichangesign=-1
      do i=1,n-1
        if(a(i)*a(i+1) <0) then
          ichangesign=i
          exit
        endif
      enddo
      end
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      character(3) function charnum3(num)
      integer(4) ::num
      charnum3 = 
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
cssssssssssssssssssssssssssssssssssssssssssssssssssss
      character(4) function charnum4(num)
      integer(4) ::num
      charnum4 = 
     &           char(48+mod(num/1000,10))//
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
cssssssssssssssssssssssssssssssssssssssssssssssssssss
      character*9 function ftoa9(arg)
      real(8):: arg
      write(ftoa9,"(1x,f8.3)") arg
      end
cssssssssssssssssssssssssssssssssssssssssssssssssssss
      character(15) function f2a(arg)
      real(8):: arg
      write(f2a,"(d13.5)") arg
      end
cssssssssssssssssssssssssssssssssssssssssssssssssssss
      character*20 function xts(num1,num2)
      integer(4) :: num,num1,num2
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
cssssssssssssssssssssssssssssssssssssssssssssssssssss
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
cssssssssssssssssssssssssssssssssssssssssssssssssssss
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
