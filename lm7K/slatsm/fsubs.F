      subroutine cexit(pv,ps)
      implicit none
#if MPI|MPIK
      include 'mpif.h'
#endif
      integer:: pv,ps
      integer:: status,ierr
      if (ps.ne.0) then
#if MPI|MPIK
        if (pv.eq.0) then
          call MPI_finalized(status,ierr)
          if (status.eq.0) then
            call MPI_finalize(ierr)
          endif
        endif
#endif
        call exit(pv)
      endif
      
      end subroutine cexit


      

      subroutine fsystm(ps,res)
#ifdef HASIFPORT
      use ifport
#endif
      character(*) ps
      integer:: res
      res = system(ps)
#ifdef HASIFPORT
      if (res.ne.0) then
         write(*,*) 'system: return code=',res,ierrno()
      endif
#endif
      end subroutine fsystm

      subroutine cwrite(ps,i1,i2,newln)
      character(*) ps
      integer:: i1,i2,newln
      integer,save:: iii
      if (newln.ne.0) then
        write(*,'(a,$)') ps(i1+1:i2+1)
      else
        write(*,'(a)')   ps(i1+1:i2+1)
      endif
      end subroutine cwrite

      subroutine flushs(iout)
      integer iout
      if (i>=0) then
        call flush(6)
      else
        write(*,*) 'flushs(',iout,
     .  ') is called, but not executed in flushs()'
      endif
c       another option is  comitqq()
      end subroutine flushs

      integer function  bitor(i1,i2)
      bitor= ior(i1,i2)
      end function  bitor

      integer function bitand(i1,i2)
      bitand= iand(i1,i2)
      end function bitand 

      subroutine nlchar(ich,ps)
      integer ich
#if 1
      character(*) ps
c   '\n'=10 
      if (ich==1) then
        ps(1:1) =  char(10) 
      endif
#else
c  may be
      character ps
c   '\n'=10
      if (ich==1) then
        ps =  char(10)
      endif
#endif
c        intel ifort accepts both

      end subroutine nlchar


      subroutine locase(ps)
      character(*) ps
      integer::i,n,shift
      n=len_trim(ps)
      shift=-ichar('A')+ichar('a') 
      do i=1,n
        if ( ichar(ps(i:i)) >= ichar('A')
     .  .and. ichar(ps(i:i)) <= ichar('Z') ) then
          ps(i:i) = char( ichar(ps(i:i))+ shift )
        endif
      enddo
      end subroutine locase


      subroutine sectim(tsec, tusec)
#ifdef HASIFPORT
      use ifport
#endif
      integer tsec,tusec
      integer:: ierr
#if HASGETTIMEOFDAY
c ------ intel fortran
      integer(4):: v(2)
c            ---- not integer:: v(2) 
      call gettimeofday(v,ierr)
      tsec= v(1)
      tusec=v(2)
#elif HASCPUTIME
c ------- gfortran
      real:: x
      call cpu_time( x)
      tsec = int(x)
      tusec = int((x-tsec)*1000000.0)
#else
c#error "needs subroutine to get cpu time"
      tsec=0.0
      tusec=0.0
#endif

      end subroutine sectim

      subroutine gtenv(pnam,pval)
#ifdef HASIFPORT
      use ifport
#endif
      character(*) pnam, pval
      integer:: ret
#ifdef HASGETENVQQ
c      ----- intel fortran
      ret = getenvqq(pnam,pval)
#elif HASGETENVIRONMENTVARIABLE
c     ----- fortran2003std
      call get_environment_variable(pnam,pval) 
#elif HASGETENV
      call getenv(pnam,pval) 
#else
#error "needs subroutine to get environment variable"
#endif
      end subroutine gtenv


      subroutine ptenv(pnam)
#ifdef HASIFPORT
      use ifport
#endif
      character(*) pnam
      integer ret
#ifdef HASSETENVQQ
c --- intel fortran
      ret=  setenvqq(pnam)
      write(*,*) 'ptenv() is called with ',pnam
      write(*,*) 'ptenv() is supported.'
#else
      write(*,*) 'ptenv() is called with ',pnam
      write(*,*) 'ptenv() not supported, but continue.'
#endif
      end subroutine ptenv 



      subroutine mkdcon(dmach,d1mach)
      implicit none
      real(8) :: dmach(0:2), d1mach(0:4)
      real(8) eps,ubt,bubt
      real(8):: t,b,l,u

      real(8):: dlamch
c machine constant from lapack routine
c base
      b = dlamch('b')
c eps*base
      eps = dlamch('p')
c emin
      l = dlamch('m')
c emax
      u = dlamch('l')
c        write(*,*) b,eps,l,u

      t = int(1.0-(log(eps)/log(b)))
c         write(*,*) 't=',t

      dmach(0) = b**(1.-t)
      dmach(1) = 100.*b**(l+t)
      dmach(2) = (b**(u-t))/100.0

      d1mach(0) = (b**(l-1.))
      ubt = u*(1.0-(b**(-t)))
      bubt = b*b**(ubt-1.0)
      d1mach(1) = bubt
      d1mach(2) = (b**(-t))
      d1mach(3) = (b**(1.0-t))
      d1mach(4) = log(b)/log(10.0)
      end subroutine mkdcon


#if 0
      integer function unlink(filename)
      implicit none
       character(*):: filename
       integer:: istart,iend,unusedfid,i
       logical:: L      
       unlink=1
       istart=100; iend=999
      do i=istart,iend
        inquire(i,opened=L)
        if (.not.L) then
          unusedfid=i
          return
        endif
      enddo
      unusedfid=0
      if (unusedfid==0) return
      open(unusedfid,file=filename)
      close(unusedfid,status='delete')
      unlink=0
      end function unlink
#endif
