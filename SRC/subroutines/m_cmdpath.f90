!>Get cmdpath, which is the full path to the command (fortran program) itself.
module m_cmdpath
   character(1024), public :: cmdpath = ''
   public setcmdpath, setcmdpathc !command path
   private
contains
   subroutine setcmdpathc(cname, nsize) bind(C) !Set cmdpath, which is the full path to the command (fortran program) itself.
      implicit none
      integer:: nsize, i
      character(1):: cname(nsize)
      character(1024):: string = ''
      forall (i=1:nsize) string(i:i) = cname(i)
      cmdpath = adjustl(string)
      print *, 'cmdpath=', trim(cmdpath)
   end subroutine setcmdpathc
   subroutine setcmdpath() !Set cmdpath, which is the full path to the command (fortran program) itself.
      use ISO_C_BINDING
      implicit none
      interface
         function readlink(path, buf, bufsize) bind(C)
            import
            integer(C_SIZE_T) :: readlink
            character(KIND=C_CHAR), intent(IN) :: path(*)
            character(KIND=C_CHAR) :: buf(*)
            integer(C_SIZE_T), value :: bufsize
         end function readlink
         function getpid() bind(C)
            import
            integer(C_SIZE_T) :: getpid
         end function getpid
      end interface
      integer :: pid, i, idx
      integer(C_SIZE_T) :: szret
      character(512) :: path
      character(KIND=C_CHAR) :: cbuf(256)
      if (len_trim(cmdpath) /= 0) return !already set
      cbuf = ''
      pid = getpid()
      write (path, '(i0)') pid
      path = '/proc/'//TRIM(path)//'/exe'
      szret = readlink(TRIM(path)//C_NULL_CHAR, cbuf, SIZE(cbuf, KIND=C_SIZE_T))
      if (szret == -1) stop 'Error reading link'
      path = ''
      do i = 1, SIZE(cbuf)
         if (cbuf(i) == C_NULL_CHAR) exit
         path(i:i) = cbuf(i)
      end do
      idx = INDEX(path, '/', BACK=.TRUE.)
      cmdpath = TRIM(path(:idx - 1))//'/'
   end subroutine setcmdpath
end module m_cmdpath
