c$$$      logical function aiomp(alabl,pmpol,nl,lmxv,nsp,ifi)
c$$$C- File I/O for ASA multipole moments integrals
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs/Outputs
c$$$Ci   alabl,nl,lmxv,nsp
c$$$Ci   ifi:    logical unit: positive for read, negative for write
c$$$Cio  pmpol: intra-atomic density-density response matrix
c$$$Cr Remarks
c$$$Cr Bugs
c$$$Cr   Input not checked when file nsp mismatches passed nsp
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$C Passed parameters
c$$$      character*8 alabl
c$$$      integer nl,lmxv,nsp,ifi
c$$$      double precision pmpol(nl,nl,0:lmxv,3,nsp)
c$$$C Local parameters
c$$$      double precision x2(10),damax
c$$$      integer i,n,ip,k,a2vec,ix2(10),ipr,nl2,nsp2,l1,l2,lf,idamax
c$$$      logical scat,rdstrn,pars1v
c$$$      character*100 s
c$$$
c$$$      aiomp = .false.
c$$$      call getpr(ipr)
c$$$C --- File READ ---
c$$$      if (ifi .gt. 0) then
c$$$C   ... return unless file has PMPOL: category
c$$$        if (.not. scat(ifi,'PMPOL:',':',.true.)) return
c$$$C   ... read nl and nsp ... abort if missing
c$$$        backspace ifi
c$$$        read(ifi,'(a72)') s
c$$$        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
c$$$        if (nl .ne. nl2 .and. ipr .ge. 10)
c$$$     .  print *, 'aiomp (warning) mismatch in nl, class '//alabl
c$$$        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
c$$$        if (nsp .ne. nsp2 .and. ipr .ge. 10)
c$$$     .  print *, 'aiomp (warning) mismatch in nsp, class '//alabl
c$$$C   ... For now, abort
c$$$        if (nsp .ne. nsp2 .or. nl .ne. nl2) goto 18
c$$$        call dpzero(pmpol,nl*nl*(lmxv+1)*3*nsp)
c$$$C   ... read pmpol
c$$$        n = min(nsp,nsp2)
c$$$        do  100  i = 1, n
c$$$        do  10  k = 1, 3
c$$$
c$$$          do  120  l1 = 0, nl-1
c$$$          do  12  l2 = 0, nl-1
c$$$
c$$$   16       if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
c$$$            if (s .eq. ' ') goto 16
c$$$            ip = 0
c$$$C    ...  Abort if fail to read lmxv+1 numbers
c$$$            ip = a2vec(s,len(s),ip,4,' ',1,-2,-(lmxv+3),ix2,x2)
c$$$            if (ip .ne. lmxv+3) call rxs('AIOMP: failed to parse ',s)
c$$$            call dcopy(lmxv+1,x2(3),1,pmpol(l1+1,l2+1,0,k,i),nl**2)
c$$$   12     continue
c$$$ 120      continue
c$$$   10   continue
c$$$ 100    continue
c$$$        aiomp = .true.
c$$$        return
c$$$   18   continue
c$$$        print *, 'aiomp: (input skipped) bad syntax, class '//alabl
c$$$        return
c$$$
c$$$C --- File WRITE ---
c$$$      else
c$$$        write(-ifi,21) alabl, nl, nsp
c$$$        do  121  i = 1, nsp
c$$$        do  20  k = 1, 3
c$$$
c$$$          do  122  l1 = 0, nl-1
c$$$          do  22  l2 = 0, nl-1
c$$$C           x2(1) = damax(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl)
c$$$            ip = idamax(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl)
c$$$            x2(1) = pmpol(l1+1,l2+1,ip-1,k,i)
c$$$            if (x2(1) .gt. 9999d0 .or. x2(1) .lt. -999d0) then
c$$$              call dcopy(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl,x2,1)
c$$$C             call awrit4('%,4i%,4i %n:1g',' ',100,6,l1,l2,lmxv+1,x2)
c$$$              call awrit4('%,4i%,4i %n:1;11F',' ',100,-ifi,l1,l2,lmxv+1,
c$$$     .            x2)
c$$$            else
c$$$              write(-ifi,333) l1,l2,(pmpol(l1+1,l2+1,lf,k,i),lf=0,lmxv)
c$$$            endif
c$$$   22     continue
c$$$  122     continue
c$$$  333     format(2i4,10f12.6)
c$$$          if (i .ne. nsp .or. k .ne. 3) write(-ifi,333)
c$$$
c$$$   20   continue
c$$$  121   continue
c$$$        aiomp = .true.
c$$$      endif
c$$$   21 format('PMPOL: ',a4,'  nl=',i1,'  nsp=',i1)
c$$$
c$$$      end
c$$$
