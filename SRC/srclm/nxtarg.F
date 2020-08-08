      logical function nxtarg(strn)
C ----------------------------------------------------------------------
Co Outputs
Co  strn   :argument returned in string a
Co  cmdstr :returns false if iarg is outside the range of arguments
Cr Remarks
Cr  nxtarg
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) strn
C ... Local parameters
      integer iarg
      logical cmdstr
      common /nxalst/ iarg

      nxtarg = cmdstr(iarg,strn)
      iarg = iarg+1

      end

      integer function nxargs(jarg)
C- Query or set current argument index
C ----------------------------------------------------------------------
Co  strn   :argument returned in string a
Ci Inputs
Ci   jarg  :if >=0, set current argument index to jarg
Ci         :if <0,  return in nxargs current value of argument index 
Co Outputs
Co  nxarg  :current value of argument index
C ----------------------------------------------------------------------
      integer jarg
      integer iarg
      common /nxalst/ iarg

      if (jarg .ge. 0) then
        iarg = jarg
      endif
      nxargs = iarg

      end

      block data dnxarg
C     implicit none
C ... Local parameters
      integer iarg
      common /nxalst/ iarg
      data iarg /-1/

      end

