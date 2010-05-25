      logical function aiorme(alabl,radme,nl,nsp,ifi)
C- File I/O for matrix elements of radial gradients of w.f.
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   alabl,nl,nsp
Ci   ifi:  logical unit: positive for read, negative for write
Cio  radme(4,2,nl,2,nsp):  <g2_l grad g1_l'> with l' = l +/- 1
Cio        i = 1 for g1=phi,    g2=phidot
Cio            2 for g1=phidot, g2=phi
Cio            3 for g1=phi,    g2=phidot
Cio            4 for g1=phidot, g2=phidot
Cio        j = 1 for <g2 | grad g1> - (l+1) < g2 | 1/r g1 >
Cio            2 for <g2 | grad g1> +    l  < g2 | 1/r g1 >
Cio       ll = 1 for l' = l-1  for l' = l+1
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,nsp,ifi
      double precision radme(4,2,nl,2,nsp)
C Local parameters
      integer is1,nl2,nsp2,ipr
      logical scat,pars1v,lfdmp
      character s*72

      aiorme = .false.
      call getpr(ipr)
      if (ifi .gt. 0) then
C   ... return unless file has RGRAD category
        if (.not. scat(ifi,'RGRAD:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nl .ne. nl2 .or. nsp .ne. nsp2)
     .  call rxs('AIORME: mismatch in nl or nsp, class ',alabl)
        call dpzero(radme,4*2*nl*2*nsp)
        do  10  is1 = 1, nsp
          aiorme = lfdmp(radme(1,1,1,1,is1),4*2*nl*2,ifi)
   10   continue
        return
   18   continue
        print *, 'aiorme: (input skipped) bad syntax, class '//alabl
      else
        write(-ifi,21) alabl, nl, nsp
        do  20  is1 = 1, nsp
          aiorme = lfdmp(radme(1,1,1,1,is1),4*2*nl*2,ifi)
   20   continue
      endif

   21 format('RGRAD: ',a4,'  nl=',i1,'  nsp=',i1)
      end

