      subroutine rdm2(ifi,s,nr,nc)
C     implicit none
C- Alternative rdm using a2bin
      integer ifi,nr,nc
      double precision s(nr,nc)
C Local variables
      integer strmax,mxline,mr
      parameter (strmax=200,mxline=5000)
      character*(strmax) instr(mxline),a
      logical a2bin,retnr
      integer cp,i,j,nrecs
      common /krdm2/ instr

C --- Read entire file into instr ---
      do 1  i = 1, mxline
    1 instr(i) = ' '
      call rdfile(ifi,'#{}%',instr,mxline,a,strmax,nrecs)
      if (nrecs .eq. 0) then
        if (nr .eq. 0) return
        call rx('rdm2: file has no records')
      endif

C --- Prepare for two passes when nr is to be determined ---
    5 continue
      retnr = .false.
      mr = nr
      if (nr .eq. 0) then 
        retnr = .true.
        nr = mxline
      endif

C --- Convert nr*nc strings to binary, or until a2bin false ---
      cp = 0
      do  10  i = 0, nr-1
        do  10  j = 0, nc-1
          call skipbl(instr,strmax,cp)
          if (cp .ge. strmax*nrecs-1) goto 91
          if (.not. a2bin(instr,s,4,i+mr*j,' ',cp,strmax*nrecs)) goto 91
   10 continue
      return

C --- If requested, determine nr and make second pass ---
   91 if (.not. retnr) stop 'rdm: missing matrix'
      nr = i
      goto 5
      end

