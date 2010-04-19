Cgetarg...info...           structure ['ioqpp', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['ioqpp', 'pot', 'spot']
Cgetarg...info...           use_to_add ['ioqpp', 'm_struc_def']
      subroutine ioqpp(lio,sctrl,spot)
       
       use m_struc_def  !Cgetarg

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
Cgetarg       double precision sctrl(1),spot(1)
       type(s_ctrl)::sctrl
       type(s_pot)::spot

C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
      integer i1,i2,nbas,nl,nsp,igets,ifi,oqpp,fopna,rdm,ipr

      call getpr(ipr)
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbas nl nspin',sctrl,nbas,nl,nsp,0,0)
       
       nbas=sctrl%nbas
       nl=sctrl%nl
       nsp=sctrl%nspin

      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
Cgetarg       oqpp = igets('pot oqpp',spot)
       oqpp = int(spot%oqpp) 

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
Cgetarg       end
       end subroutine ioqpp 


