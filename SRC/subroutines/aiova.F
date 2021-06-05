      logical function aiova(alabl,vintra,nl,lmax,nsp,ifi)
C- File I/O for ASA intra-atomic density-density response matrix
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   alabl,nl,lmax,nsp
Ci   ifi:    logical unit: positive for read, negative for write
Cio  vintra: intra-atomic density-density response matrix
Cr Remarks
Cr Bugs
Cr   Input not checked when file nsp mismatches passed nsp
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmax,nsp,ifi
      double precision vintra(nl,nl,nsp,nsp)
C Local parameters
      double precision x2(0:10)
      integer i,i2,l,n,ip,k,a2vec,ix2(7),ipr,nl2,nsp2,lmx
      logical scat,rdstrn,pars1v
      character*120 s

      aiova = .false.
      call getpr(ipr)
C --- File READ ---
      if (ifi .gt. 0) then
C   ... return unless file has VINTRA: category
        if (.not. scat(ifi,'VINTRA:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a120)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (lmax+1 .ne. nl2 .and. ipr .ge. 10)
     .  print *, 'aiova (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp .ne. nsp2 .and. ipr .ge. 10)
     .  print *, 'aiova (warning) mismatch in nsp, class '//alabl
        call dpzero(vintra,nl*nl*nsp*nsp)
        lmx = min(nl,nl2)-1
C   ... read vintra
        n = min(nsp,nsp2)
        do   i = 1, n
        do   l = 0, lmx
          if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
          ip = 0
          k = a2vec(s,len(s),ip,4,' ',1,-2,-(lmx+1)*n-1,ix2,x2)
C    ...  Abort if failed to read lmx+2 numbers
          if (k .ne. (lmx+1)*n+1) call rx('AIOVA: failed to parse '//s)
          call dcopy(lmx+1,x2(1),1,vintra(l+1,1,i,1),nl)
          if (n.eq.2) call dcopy(lmx+1,x2(lmx+2),1,vintra(l+1,1,i,2),nl)
          if (nsp2 .lt. nsp) then
            call dcopy(lmx+1,x2(1),1,vintra(l+1,1,2,2),nl)
            call dcopy(lmx+1,x2(1),1,vintra(l+1,1,1,2),nl)
            call dcopy(lmx+1,x2(1),1,vintra(l+1,1,2,1),nl)
          endif
        enddo
        enddo
        aiova = .true.
        return
   18   continue
        print *, 'aiova: (input skipped) bad syntax, class '//alabl
        return

C --- File WRITE ---
      else
        write(-ifi,21) alabl, lmax+1, nsp
        do  i = 1, nsp
        do  l = 0, lmax
          write(-ifi,333) l,((vintra(l+1,k,i,i2), k=1,1+lmax), i2=1,nsp)
 333      format(i4,8f12.6)
        enddo
        enddo
        aiova = .true.
      endif
   21 format('VINTRA: ',a4,'  nl=',i1,'  nsp=',i1)

      end

