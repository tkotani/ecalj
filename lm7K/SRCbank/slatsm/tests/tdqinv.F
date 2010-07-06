C test dqinv
C     To check many cases, eg
C     foreach n ( `mlist 150:188` )
C     a.out << END | grep -E 'overwritten|factor'
C
C     $n
C
C     END
C     end
      subroutine fmain
C     implicit none
      INTEGER   lda, ldb, ldc, ldw, n, i, j,
     .  opw, opw2, wksize, dinv, nlev, nmin
      logical a2bin,cmdopt
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,wksize=6000000)
      character*1 cs, outs*80
      double precision sa( lda, lda ), pa( lda, lda ), pw( lda, lda)
      double precision s_time, temp, s_mflops, q_time, q_mflops
      double precision cpusec,xx,errmx
      integer ipiv(lda)
      integer w(wksize)
      common /w/ w
      common /static/ sa, pa, pw

C --- Setup ---
      call finits(2,0,0,i)
      call wkinit(wksize)
      call nada
      n = 300
      cs = 'n'
      print *, 'cs= (n or s)?'
      read(*,'(a1)') cs
      call initqu(.true.)
      print *, 'test dqinv: cs=',cs
      call query('n=?',2,n)
      if (n .gt. lda) call rx('n gt lda')
      mflops = 2.0d-6 * n * n * n

      j = 6
      if (cmdopt('-nmin=',j,0,outs)) then
        if (.not. a2bin(outs,nmin,2,0,' ',j,72))
     .    call rxs('cannot parse',outs)
        call dqinv0(nmin,1)
      endif

C#ifdefC gfortran
C      nmin = 32
C#else
      call dqinv0(nmin,0)
C#endif

C --- Matrix inversion, dgefa,dgedi ---
      CALL init ( sa,lda,n, cs )
*     call ywrm(.false.,'a',1,6,'(8f16.10)',sa,lda,n,n)
      temp = cpusec()
      if (cs .eq. 's') then
        call dsifa(sa,lda,n,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')

        call dsidi(sa,lda,n,pw,xx,i,pw(1,2),1)
        do  16  i = 1, n
        do  16  j = 1, i
   16   sa(i,j) = sa(j,i)
      else
C        call dgefa(sa,lda,n,pw,i)
C        if (i .ne. 0) call rx(' matrix is singular')
C        call dgedi(sa,lda,n,pw,xx,pw(1,2),1)
        call dgetrf(n,n,sa,lda,ipiv,i)
        if (i .ne. 0) call rx(' matrix is singular')
        call dgetri(n,sa,lda,ipiv,pw,lda*lda,i)
        if (i .ne. 0) call rx(' matrix is singular')
      endif
      s_time = cpusec() - temp
      if (s_time .ne. 0) then
        s_mflops = mflops / s_time
      else
        s_mflops = 0
      endif
*     call ywrm(.false.,'a^-1',1,6,'(4f16.10)',sa,lda,n,l)

C --- Inversion by dqinv ---
      CALL init ( pa,lda,n, cs )
C ... (debugging) check usage of w
      ldw = min(n+5,lda)
      write(*,'('' using lda,ldw,n,nmin='',4i5)') lda,ldw,n,nmin
      call defrr(opw,  ldw*(n))
      call defrr(opw2, n*n)
      call dvset(w(opw),1,ldw*(n),-99d0)
      call dvset(w(opw2),1,n*n,-99d0)

      temp = cpusec()
C      call dqinv(cs,pa,lda,2,n,w(opw),ldw,i)
      i = dinv(cs,n,lda,pa)

      q_time = cpusec() - temp
      if (q_time .ne. 0) then
        q_mflops = mflops / q_time
      else
        q_mflops = 0
      endif

C --- Check quality of inversion ---
      CALL diff('compare dgetrf,i to  dqinv', sa, pa, lda, n )
      CALL init ( sa,lda,n, cs )
      temp = cpusec()
      call dgemm('N','N',n,n,n,1d0,sa,lda,pa,lda,0d0,pw,lda)
      temp = cpusec() - temp + 1d-8
      do  40  i = 1, n
   40 pw(i,i) = pw(i,i) - 1
      errmx = 0
      do 42 i = 1, n
      do 42 j = 1, n
   42 errmx = max(errmx,dabs(pw(i,j)))
      print 100, 'inv(a) * a - 1     errmax ',errmx,nint(mflops/temp)
  100 format(1X,a,':  errmx=',g10.3,2i4)

      print 110, s_time,s_mflops,q_time,q_mflops,s_time/(q_time+1d-6),n

  110 format(/1X,'Serial time: ',F7.3,'  Serial MFlops: ',F6.1,
     .       /1X,' dqinv time: ',F7.3,'   dqinv MFlops: ',F6.1
     .       /1X,'     factor: ',F7.3,'  for n =',i4)

*     call ywrm(.false.,'pw',1,6,'(16f8.2)',pw,ldw,n+1,n)

C ... Make sure pw not overwritten beyond ldw*(n)
      call wkchk('check if pw2 still intact')
      call dcopy(n*n,w(opw2),1,pw,1)
      j = 1
      do  50  i = 1, n*n
        if (pw(i,j) .ne. -99d0) call fexit2(-1,111,
     .    ' STOP: pw overwritten, i,j=%i %i',i,j)
   50 continue
      if (ldw .gt. n) then
        call dmcpy(w(opw),ldw,1,pw,lda,1,ldw,n)
        do  52  i = 1, ldw
        do  52  j = 1, n
            if (i .le. n .and. j .le. n) goto 52
            if (pw(i,j) .ne. -99d0) call fexit2(-1,111,
     .        ' STOP: pw overwritten, i,j=%i %i',i,j)
   52   continue
      endif
      end
      integer function dinv(cs,n,lda,a)
C- Inversion of a real matrix
C     implicit none
      character*1 cs
      integer n,lda
      double precision a(lda,n)
      integer ldw,ow,i
      integer w(1)
      common /w/ w

      ldw = n
      call defrr(ow,  ldw*(n+1))
      call dqinv(cs,a,lda,2,n,w(ow),ldw,i)
      call rlse(ow)
      dinv = i
      end

      SUBROUTINE init ( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double precision a( lda, * )
      character *1 cs
      real ran1
      call ran1in(1)

      do 10 i = 1, n
      do 10 j = 1, n
   10 a(i,j) = ran1()

      if (n .eq. 4) then
        a(1,1) = 16
        a(2,1) = 1
        a(3,1) = 2
        a(4,1) = 3
        a(1,2) = 3
        a(2,2) = 14
        a(3,2) = 2
        a(4,2) = 1
        a(1,3) = 1
        a(2,3) = 1
        a(3,3) = 12
        a(4,3) = 3
        a(1,4) = 3
        a(2,4) = 1
        a(3,4) = 2
        a(4,4) = 10
      endif

      if (cs .eq. 's') then
        do 20 i = 1, n
        do 20 j = 1, i
   20   a(i,j) = a(j,i)
      endif

      end
      SUBROUTINE diff (strn, sc, pc, ldc, n)
C- Compare the two arrays for differences
C     implicit none
      character*(*) strn
      integer ldc, n, m, i, j
      double precision sc( ldc, * ), pc( ldc, * ), errmx
      errmx = 0d0
      do 10 i = 1, n
      do 10 j = 1, n
   10 errmx = max(errmx,dabs(sc(i,j)-pc(i,j)))
      print 100, strn,errmx
  100 format(1X,a,':  errmx=',g10.3)
      end
C      subroutine ywrm(lbin,filel,icast,ifi,fmt,s,ns,nr,nc)
CC- Writes complex matrix to file ifi
CCr lbin: writes in binary mode
C      logical lbin
C      character*(*) filel
C      integer icast,nr,nc
C      double precision s(ns,nc,2)
C      character*(*) fmt, outs*10
C      integer i,j
C
C      if (lbin) then
C        if (filel .eq. ' ') then
C          write(ifi) nr,nc,icast
C        else
C          write(ifi) nr,nc,icast,len(filel)
C          write(ifi) filel(1:len(filel))
C        endif
C        call dpdump(s,nr*nc*mod(icast,10),-ifi)
C        return
C      endif
C      outs = ' '
C      if (icast .eq. 1)  outs = ' real'
C      if (icast .eq. 11) outs = ' symm'
C      if (icast .eq. 2)  outs = ' complex'
C      if (icast .eq. 12) outs = ' herm'
C      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)
Cc      if (filel .ne. ' ') call awrit0(filel,' ',len(filel),ifi)
C      if (filel .ne. ' ') print 333, filel
C  333 format('#',a)
C      do  10  i = 1, nr
C   10 write(ifi,fmt) (s(i,j,1), j=1,nc)
C      if (mod(icast,10) .gt. 1) then
C       write(ifi,'(1x)')
C      do  20  i = 1, nr
C   20 write(ifi,fmt) (s(i,j,2), j=1,nc)
C      endif
C      end
C      subroutine prm(a,lda,n)
C      integer lda,n
C      double precision a(lda,lda)
C
C      call ywrm(.false.,'from prm',1,6,'(4f12.6)',a,lda,n,n)
C      end
