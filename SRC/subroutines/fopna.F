      integer function lgunit(i)
      use m_ext,only: sname
      implicit none
      character(10):: i2char
      character*100 ext
      integer i, fopn, i1mach, fhndl,ifile_handle,procid
      integer,save:: lgunit1=0,lgunit2=0,lgunit3=0
C- Returns stdout for i=1, log for i=2, mlog for i=3 (MPI logfile)
      include "mpif.h"
      integer mpipid
      procid = mpipid(1) 
      lgunit = 6!i1mach(2)
      if (i .eq. 1) return
      if (i .eq. 2) then
         if(lgunit2==0) then
            lgunit2=ifile_handle()
            open(lgunit2,file='log.'//trim(sname),position='append')
         endif
         lgunit = lgunit2
      elseif (i .eq. 3) then
         if(lgunit3==0) then
            lgunit3=ifile_handle()
            open(lgunit3,file='mlog.'//trim(sname)//'_'//trim(i2char(procid)))
         endif
         lgunit = lgunit3
      endif
      end
      
      subroutine namcat(dest,fnam1,fnam2,lcase)
C- Concatenates two names, optionally lowering the case
      implicit none
      logical lcase
      character*(*) fnam1,fnam2,dest
      integer iprint,i1,i2,j1,j2
      call word(fnam1,1,i1,i2)
      call word(fnam2,1,j1,j2)
      i2 = max(i1-1,i2)
      j2 = max(j1,j2)
      if (i2 .ge. i1) then
        dest = fnam1(i1:i2)//fnam2(j1:j2)
      else
        dest = fnam2(j1:j2)
      endif
      if (i2-i1+j2-j1+2 .gt. len(dest) .and. iprint() .gt. 0) then
        print *, 'namcat: string "',fnam1(i1:i2),fnam2(j1:j2),
     .  '" was truncated to "',dest,'"'
      endif
      if (lcase) call locase(dest)
      end
      
