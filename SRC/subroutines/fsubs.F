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
#ifdef __INTEL_COMPILER
C#if  ( defined __INTEL_COMPILER && __INTEL_COMPILER!=1110)
      use ifport
C#endif
#endif
      character(*) ps
      integer:: res
      res = system(ps)
      if (res.ne.0) then
#ifdef __INTEL_COMPILER
         write(*,*) 'system: return code=',res,ierrno()
#else
         write(*,*) 'system: return code=',res
#endif
C fortran2008std
C     call  execute_command_line (ps,exitstat=res) 
      endif
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
      integer tsec,tusec
      integer:: ierr
      logical,save:: firsttime=.true.
      real :: diff
      integer,save:: i1
      integer:: i2,irate,imax
      real:: x
      call cpu_time( x)
      tsec = int(x)
      tusec = int((x-tsec)*1000000.0)
      end subroutine sectim

      subroutine gtenv(pnam,pval)
      character(*) pnam, pval
      integer:: ret
      call get_environment_variable(pnam,pval)
      end subroutine gtenv

      subroutine ptenv(pnam)
      character(*) pnam
      integer ret
      write(*,*) 'ptenv() is called with ',pnam
      write(*,*) 'ptenv() not supported, but continue.'
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
c$$$
c$$$
c$$$#if 0
c$$$      integer function unlink(filename)
c$$$      implicit none
c$$$       character(*):: filename
c$$$       integer:: istart,iend,unusedfid,i
c$$$       logical:: L      
c$$$       unlink=1
c$$$       istart=100; iend=999
c$$$      do i=istart,iend
c$$$        inquire(i,opened=L)
c$$$        if (.not.L) then
c$$$          unusedfid=i
c$$$          return
c$$$        endif
c$$$      enddo
c$$$      unusedfid=0
c$$$      if (unusedfid==0) return
c$$$      open(unusedfid,file=filename)
c$$$      close(unusedfid,status='delete')
c$$$      unlink=0
c$$$      end function unlink
c$$$#endif
