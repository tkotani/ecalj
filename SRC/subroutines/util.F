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
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      character(10) function i2char(num)
!! convert num>0 to char. See charnum3 to understand this.
      implicit none
      integer(4) ::num,itens,iten
      i2char=''
      do itens=9,1,-1 !itens mean itens+1's digit
        iten=10**itens
        if(num>=iten) i2char=trim(i2char)//char(48+mod(num/iten,10))
      enddo
      i2char = trim(i2char)//char(48+mod(num,10)) !1st digit
      end
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      character(3) function charnum3(num)
      integer(4) ::num
      charnum3 = 
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
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
