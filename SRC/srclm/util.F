!! ----------------------------------------------------
      integer function ifile_handle()
!! find unused file handle
      implicit none
      integer:: i
      logical:: nexist
      integer,save:: irem=2001
      character*256::nnn
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

      character(10) function i2char(num)
!! convert num>0 to char. See charnum3 to understand this.
      implicit none
      integer(4) ::num,itens,iten,i
      i2char=''
      do itens=9,1,-1 !itens mean itens+1's digit
        iten=10**itens
        if(num>=iten) i2char=trim(i2char)//char(48+mod(num/iten,10))
      enddo
      i2char = trim(i2char)//char(48+mod(num,10)) !1st digit
c      print *,'num ',num,' i2char= ',i2char,' xxx',48+mod(num,10)
c      print *,'i2char ',char(48+mod(num,10)) , 'xxx',iii//char(48+mod(num,10)) 
      end
c$$$      character(3) function charnum3(num) !3digit.
c$$$      integer(4) ::num
c$$$      charnum3 = 
c$$$     &           char(48+mod(num/100,10))//
c$$$     &           char(48+mod(num/10,10))//
c$$$     &           char(48+mod(num,10))
c$$$      end
c$$$      character(26) function f2char(adat)
c$$$      real(8):: adat
c$$$      character(38):: bbb
c$$$      write(bbb,"(f38.16)")adat
c$$$      f2char=adjustl(bbb)
c$$$      end

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
