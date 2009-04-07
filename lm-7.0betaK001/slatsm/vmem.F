      subroutine vmem(o1,o2)
C- Virtual memory routines
C ----------------------------------------------------------------
Ci Inputs
Ci   o1,o2
Co Outputs
Co   Nothing
Cr Remarks
Cr
C ----------------------------------------------------------------
C Passed parameters
      integer o1,o2
C Local parameters
      logical lsave
      integer len,offst,fopn,fopno,ifi,i,ii
      save len,lsave,offst
C heap:
      integer w(1)
      common /w/ w

      offst = o1
      len = -o2
      lsave = (o2 .lt. 0)
      if (o2 .lt. 0) o2 = o1
      return

      entry vmems(i)
      if (.not. lsave) return
      ifi = fopn('TMP')
      rewind ifi
      do  10  ii = 1, i
   10 call vmem2(-ifi,w(offst),len)
      return

      entry vmemg(i)
      if (.not. lsave) return
      ifi = fopno('TMP')
      rewind ifi
      do  20  ii = 1, i
   20 call vmem2(ifi,w(offst),len)
      call fclose(ifi)

      end
      subroutine vmem2(ifi,w,len)
      integer ifi,len,w(len)

      if (ifi .gt. 0) read(ifi) w
      if (ifi .lt. 0) write(-ifi) w
      end

