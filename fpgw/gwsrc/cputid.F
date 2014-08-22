c Time routines -------------------------------------------------------------
c This dummy time routine
#ifdef ETIME_DUMMY
      real*4 function etime(ar)
      real*4 ar(2)
      ar = 0.0
      end
#endif
c
c       subroutine fdate(s)
c       character(*) s
c       s = ' '
c       end
c----------------------------------------------------------------------------





c----------------------------------------------------------------------
      double precision function cpumin(ifile,cpu0)
      implicit real*8(a-h,o-z)

c cpu time in minutes for RISC: mclock()/600.d0
c cpusek is the total cputime - cpu0

c     cpumin   = mclock()/600.d0 - cpu0
      cpusec   = (secnds(0.0) - cpu0)
      cpumin   = cpusec/60.d0
c     cpusec=0d0
c     cpumin=0d0
c----
      if(ifile .le. 0)then
        write(*,6400)cpusec,cpumin
      else
        write(ifile,6400)cpusec,cpumin
      end if
 6400 format(/' cpu time = ',f15.4,' secs ',f15.4,' mins ')

      return
      end
c----------------------------------------------------------------------
c     double precision function cpusec(cpu0)
c     implicit real*8(a-h,o-z)
c
c cpu time in seconds for RISC
c cpusek is the total cputime - cpu0
c
c     cpusec   = mclock()/100.d0 - cpu0
c     return
c     end
c----------------------------------------------------------------------
#define REALTIME
      subroutine cputid(ifile)
      implicit none
      include 'mpif.h'
c cpu time in seconds and minutes for RISC IBM: mclock()/100 (secs)
c                     for HP: secnds(0.0)
c                        DEC: etime(etw)
c ifile = file number where the cpu time information is
c         to be printed
c ifile = 0 ==> screen
      real*4 cpuetime, etw(2),cpulast, etime
      real*4 cpu0 /-1.d0/
      save cpu0
      character*(*):: message
      real(4):: cpusec,cpumax,cpumin
      integer:: ierr,rank ,ifile
#ifdef REALTIME
      logical,save::firsttime=1
      integer,save::i1
      integer:: i2,irate,imax
      real(4)::diff
#endif
Cstop2rx 2013.08.09 kino      if(ifile .lt. 0)stop 'cputid: negative unit number'
      if(ifile .lt. 0)call rx( 'cputid: negative unit number')
#ifdef REALTIME
      if (firsttime) then 
         call system_clock(i1)
         firsttime=.false.
         cpusec=0.0d0
         cpumin=0.0d0
      else
         call system_clock(i2,irate,imax)
         diff=i2-i1
         if (diff<0) diff=imax-i1+i2
         diff=diff/dble(irate)
         cpusec=diff
         cpumin=cpusec/60.d0
      endif
#else
      cpuetime = etime(etw)
      if (cpu0 .lt. 0.d0) cpu0 = cpuetime
      cpusec   = cpuetime - cpu0
      cpumin   = cpusec/60.d0
#endif
      if(ifile .eq. 0)then
        write(*,6400)cpusec,cpumin
      else
        write(ifile,6400)cpusec,cpumin
      end if
 6400 format(' CPU= ',f12.4,' secs',f7.1,' mins')
      return
!!
      entry cputid2(message,ifile)
Cstop2rx 2013.08.09 kino      if(ifile .lt. 0)stop 'cputid: negative unit number'
      if(ifile .lt. 0)call rx( 'cputid: negative unit number')
#ifdef REALTIME
      if (firsttime) then 
         call system_clock(i1)
         firsttime=.false.
         cpusec=0.0d0
         cpumin=0.0d0
      else
         call system_clock(i2,irate,imax)
         diff=i2-i1
         if (diff<0) diff=imax-i1+i2
         diff=diff/dble(irate)
         cpusec=diff
         cpumin=cpusec/60.d0
      endif
#else
      cpuetime = etime(etw)
      if (cpu0 .lt. 0.d0) cpu0 = cpuetime
      cpusec   = cpuetime - cpu0
      cpumin   = cpusec/60.d0
#endif
      if(ifile .eq. 0)then
        write(*,6401)trim(message),cpusec,cpumin
      else
        write(ifile,6401)trim(message),cpusec,cpumin
      end if
 6401 format(a,' CPU= ',f12.4,' secs',f7.1,' mins')


      entry mpi_cputid2(message,ifile)
Cstop2rx 2013.08.09 kino      if(ifile .lt. 0)stop 'cputid: negative unit number'
      if(ifile .lt. 0)call rx( 'cputid: negative unit number')
      if (firsttime) then
         call system_clock(i1)
         firsttime=.false.
         cpusec=0.0d0
      else
         call system_clock(i2,irate,imax)
         diff=i2-i1
         if (diff<0) diff=imax-i1+i2
         diff=diff/dble(irate)
         cpusec=diff
      endif

      call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call mpi_reduce(cpusec,cpumax,1,MPI_REAL4,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(cpusec,cpumin,1,MPI_REAL4,MPI_MIN,0,MPI_COMM_WORLD,ierr)

      if (rank.eq.0) then 
      if(ifile .eq. 0)then
        write(*,6402)trim(message),cpumax,cpumax/60.0d0 , cpumin,cpumin/60.0d0
      else
        write(ifile,6402)trim(message),cpumax,cpumax/60.d0, cpumin,cpumin/60.d0
      end if
      endif
 6402 format(a,', realtime nodeMAX=',f12.4,' secs (',f7.1,' mins), MIN=',f12.4,' secs (',f7.1,' mins)')

      return
      end



      module m_realtimediff
      implicit none
      integer,parameter:: nmax=10000
      integer:: i1(nmax)=nmax*0 
      character(30):: msg(nmax)=''
      real(8):: elapse(nmax)=nmax*0.0d0
      real(8):: start(nmax)=nmax*0.0d0
      end module m_realtimediff


      subroutine realtimediff(id,msg0)
      use m_realtimediff 
      implicit none
      integer,intent(in):: id
      character(*),intent(in)::msg0

      integer::mode
      real(8):: diff
      integer:: i2,imax
      character(20):: errmsg
#define USEMPIWTIME 
#ifdef USEMPIWTIME
      real(8):: mpi_wtime

      if (id<=0 .or. id>nmax) then 
          write(errmsg,"(2i5)") id,nmax
          call rx('fatal error, id<=0 or id>nmax, id,nmax='//trim(errmsg))
      endif

      mode=0
      if (len_trim(msg0)>0) mode=1 

      if (mode.eq.0) then
        i1(id)=i1(id)+1
        start(id)=mpi_wtime()
      else
        diff=mpi_wtime()
        msg(id)=msg0
        elapse(id)=elapse(id)+diff-start(id)
      endif
#else
      if (mode.eq.0) then 
        call system_clock(i1(id))
      else
        msg(id)=msg0
        call system_clock(i2,irate,imax)
        diff = i2-i1(id)
        if (i1(id)>i2) diff=imax-i1(id)+i2
        diff= dble(i2-i1(id))/dble(irate)  ! real time
       ! write(*,*)'mpi=',str_rank,trim(msg),diff
        elapse(id)=elapse(id)+diff  
      endif 
#endif

      end subroutine realtimediff

      subroutine print_realtimediff()
!      use m_mpi
      use m_realtimediff
      implicit none
      integer:: i,j,isize
      character(30):: str
         do i=1,nmax
           !if (i1(i).ne.0 .and. len_trim(msg(i))>0 ) then 
           if (i1(i).ne.0 ) then 
           str=msg(i)
           do j=1,30
             if (str(j:j).eq.' ') str(j:j)='_'
           enddo
           write(6,'(a,i4,1x,a,1x,F10.2,1x,a,1x,i5,1x,a)')'timediff,id=',i,'('//trim(str)//')',elapse(i),'sec' ,i1(i),'times'
!            write(6,*)'timediff,id=',i,trim(str),elapse(i),'sec', i1(i),'times'
           endif 
         enddo
      end subroutine print_realtimediff

