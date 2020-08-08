      implicit none
      integer(4):: i,irow,icol,j,iq,idat,ndat,nq,idd
      character(10000) ::iline
      character(10) :: label(20),flag='nul'
      real(8),allocatable:: dmat(:,:)

      character(3):: charnum3
      character(100):: ddd

      open(116,file='table.x')
c--------data plot
c data block is
c   ndat \times nq
      ndat=5
      nq=3
      idd=0
      allocate(dmat(ndat,nq))
      do 
        read(5,"(a)",end=9999) iline
        if(iline(1:1) =='!'.or.iline(1:1) =='#'.or.iline(1:1)=='%')
     &    cycle
        if(iline(1:1) =='$') then
          if(iline(2:2)/=' ') then
            read(iline(2:),* ) flag
            write(6,*)'flag=',flag
          else
            flag='data'
          endif
c          write(6,"(a)") trim(iline)
          if(trim(flag)=='data') then
             idd=idd+1
            open(16,file='sdat.'//charnum3(idd))
            ddd= '  #'//trim(iline(2:))
            write(16,"('# data: ',a)") trim(iline(2:))
            do idat=1,ndat 
              read(5,"(a)") iline
              read(iline,*) dmat(idat,1),dmat(idat,2),dmat(idat,3)
            enddo
c            do iq=1,nq
c            write(6,"(20f9.4)") iq, dmat(1:ndat,iq)
c            enddo

             write(116,"(7f9.2,$)") 
     &       dmat(4:5,1), dmat(3:5,2), dmat(2:4,3)
             write(116,"(a)") trim(ddd)
             
             cycle

            write(16,1122) 1, dmat(2,1)
            write(16,1122) 2, dmat(2,2)
            write(16,1122) 3, dmat(2,3)
            write(16,"(' = =')")
            write(16,1122) 1, dmat(3,1)
            write(16,1122) 2, dmat(3,2)
            write(16,1122) 3, dmat(3,3)
            write(16,"(' = =')")
            write(16,1122) 1, dmat(4,1)
            write(16,1122) 2, dmat(3,2)
            write(16,"(' = =')")
            write(16,1122) 2, dmat(4,2)
            write(16,1122) 3, dmat(3,3)
            close(16)
          endif
        endif
 1122   format(i5,20f9.2,a)
ccccccccccccc
c       if(trim(flag)=='Expt') then
c          write(6,*)'xxxxxxx=',flag
      enddo
 9999 continue
      end

      character(3) function charnum3(num)
      integer(4) ::num
      charnum3 = 
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
