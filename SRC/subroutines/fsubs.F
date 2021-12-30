      subroutine cexit(pv,ps)
      implicit none
#if MPI|MPIK
      include 'mpif.h'
#endif
      integer:: pv,ps,i
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
      integer:: i1,i2,newln,i
      integer,save:: iii
      if (newln.ne.0) then
        write(*,'(a,$)') ps(i1+1:i2+1)
      else
        write(*,'(a)')   ps(i1+1:i2+1)
      endif
      end subroutine cwrite

      subroutine flushs(iout)
      integer iout,i
      if (i>=0) then
        call flush(6)
      else
        write(*,*) 'flushs(',iout,
     .  ') is called, but not executed in flushs()'
      endif
c       another option is  comitqq()
      end subroutine flushs

      integer function  bitor(i1,i2)
      integer::i1,i2
      bitor= ior(i1,i2)
      end function  bitor

      integer function bitand(i1,i2)
      integer:: i1,i2
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

c      subroutine gtenv(pnam,pval)
c      character(*) pnam, pval
c      integer:: ret
c      call get_environment_variable(pnam,pval)
c      end subroutine gtenv

c$$$      subroutine ptenv(pnam)
c$$$      character(*) pnam
c$$$      integer ret
c$$$      write(*,*) 'ptenv() is called with ',pnam
c$$$      write(*,*) 'ptenv() not supported, but continue.'
c$$$      end subroutine ptenv 

c$$$      subroutine mkdcon(dmach,d1mach)
c$$$      implicit none
c$$$      real(8) :: dmach(0:2), d1mach(0:4)
c$$$      real(8) eps,ubt,bubt
c$$$      real(8):: t,b,l,u
c$$$      real(8):: dlamch
c$$$c machine constant from lapack routine
c$$$c base
c$$$      b = dlamch('b')
c$$$c eps*base
c$$$      eps = dlamch('p')
c$$$c emin
c$$$      l = dlamch('m')
c$$$c emax
c$$$      u = dlamch('l')
c$$$c        write(*,*) b,eps,l,u
c$$$
c$$$      t = int(1.0-(log(eps)/log(b)))
c$$$c         write(*,*) 't=',t
c$$$      dmach(0) = b**(1.-t)
c$$$      dmach(1) = 100.*b**(l+t)
c$$$      dmach(2) = (b**(u-t))/100.0
c$$$      d1mach(0) = (b**(l-1.))
c$$$      ubt = u*(1.0-(b**(-t)))
c$$$      bubt = b*b**(ubt-1.0)
c$$$      d1mach(1) = bubt
c$$$      d1mach(2) = (b**(-t))
c$$$      d1mach(3) = (b**(1.0-t))
c$$$      d1mach(4) = log(b)/log(10.0)
c$$$      end subroutine mkdcon
      
c$$$      integer function unlink(filename)
c$$$      integer:: id
c$$$      character(*):: filename
c$$$      open(newunit=id,file=filename)
c$$$      close(id,status='delete')
c$$$      end function unlink
