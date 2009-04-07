      subroutine ioqpp(lio,sctrl,spot)
C- File i/o for phi-phi, phi-dot, dot-dot products
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lio: true for write, false for read
Ci          <0 write
Cu Updates
Cu   08 Nov 07 (J. Xu) qpp is complex
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lio
      double precision sctrl(1),spot(1)
C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
      integer i1,i2,nbas,nl,nsp,igets,ifi,oqpp,fopna,rdm,ipr

      call getpr(ipr)
      call upack('ctrl nbas nl nspin',sctrl,nbas,nl,nsp,0,0)
      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
      oqpp = igets('pot oqpp',spot)
      ifi = fopna('qpp',-1,4+8)

      if (lio) then
        call ywrm(1,'lmasa',1,ifi,' ',w(oqpp),1,i1,i1,i2)
        if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
        if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
        call pshpr(0)
        if (rdm(ifi,2,i1*i2,' ',w(oqpp),i1,i2) .lt. 0) then
          if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
          call dvset(w(oqpp),1,1,-1d0)
        endif
        call poppr
      endif

      call fclose(ifi)
      end

