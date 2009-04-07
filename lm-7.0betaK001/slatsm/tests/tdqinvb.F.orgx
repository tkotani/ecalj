C Tests dqinvb
      subroutine fmain
      implicit none
      INTEGER   lda, ldb, ldc, ldw, n, nb, i, j, imax, jmax, nnb,
     .  opw, opw2, wksize, j1,j2, nlev
      logical a2bin,cmdopt
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,wksize=6000000)
C     PARAMETER (lda=401,ldb=401,ldc=401,wksize=6000000)
      character cs*4,outs*80
      double precision sa(lda,lda),sb(ldb,ldb),
     .  pb(ldb,ldb),pw(lda,lda),pw2(ldb*ldb+ldb),piv(lda)
      double precision s_time, temp, s_mflops, q_time, q_mflops
      double precision xflops,xtim
      double precision cpusec,xx
      character*10 strn,trans*1
      integer w(wksize)
      common /w/ w
      common /static/ sa, sb, pb, pw, pw2

C --- Setup ---
      call finits(2,0,0,i)
      call initqu(.true.)
      call wkinit(lda*lda*4)
      call nada
      n = 400

C ... command-line switches -nlev=
      nlev = 2
      j = 6
      if (cmdopt('-nlev=',j,0,outs)) then
        if (.not. a2bin(outs,nlev,2,0,' ',j,72))
     .    call rxs('cannot parse',outs)
      endif

      cs = 's'
      print *, 'cs= (t, s, ts)?'
      read(*,'(a3)') cs
      print *, 'test dqinvb: cs=',cs
      call query('n=?',2,n)
      nb = n-1
      call query('nb=?',2,nb)
      nnb = max(n,nb)
      if (n .gt. lda .or. nb .gt. lda) call rx('n gt lda')
      mflops = 2.0d-6 * n * n * n

C --- Matrix inversion, dgefa,dgedi ---
      call ran1in(1)
      CALL init ( sa,lda,nnb, cs )
      CALL init ( sb,ldb,nnb, cs )

C      call ywrm(0,'a',1,6,'(8f16.10)',sa,lda,n,n)
C      call ywrm(0,'b',1,6,'(8f16.10)',sb,ldb,n,nb)
      temp = cpusec()
      if (cs(1:1) .eq. 's' .or. cs(2:2) .eq. 's' .or.
     .    cs(3:3) .eq. 's' .or. cs(4:4) .eq. 's') then
C        print *, 'using dsifa,sl ..'
C        call dsifa(sa,lda,n,pw,i)
C        if (i .ne. 0) call rx(' matrix is singular')
C        do  11  j = 1, nb
C   11   call dsisl(sa,lda,n,pw,sb(1,j))
C        strn = 'dsifa,sl'
        call dsytrf('L',n,sa,lda,piv,pw,lda*lda,i)
        if (i .ne. 0) call rx('matrix is singular')
        if (cs(1:1) .eq. 't') then
          do  24  j = 1, nb
          do  24  i = 1, n
   24     pb(i,j) = sb(j,i)
          call dsytrs('L',n,nb,sa,lda,piv,pb,ldb,i)
          if (i .ne. 0) call rx('bad call to dsytrs')
          do  28  j = 1, nb
          do  28  i = 1, n
   28     sb(j,i) = pb(i,j)
        else
          call dsytrs('L',n,nb,sa,lda,piv,sb,ldb,i)
          if (i .ne. 0) call rx('bad call to dsytrs')
        endif
        strn = 'dsytrf,rs'
      else
C        call dgefa(sa,lda,n,pw,i)
C        if (i .ne. 0) call rx(' matrix is singular')
C        if (cs(1:1) .eq. 't') then
C          do  14  j = 1, nb
C            do  16  i = 1, n
C   16       pb(i,j) = sb(j,i)
C            call dgesl(sa,lda,n,pw,pb(1,j),1)
C            do  18  i = 1, n
C   18       sb(j,i) = pb(i,j)
C   14     continue
C*         call ywrm(0,'b a^-1',1,6,'(5f16.10)',sb,ldb,nb,n)
C        else
C          do  12  j = 1, nb
C   12     call dgesl(sa,lda,n,pw,sb(1,j),0)
C*         call ywrm(0,'lpack a^-1 b',1,6,'(4f16.10)',sb,ldb,n,nb)
C        endif
C        strn = 'dgefa,sl'
        call dgetrf(n,n,sa,lda,piv,i)
        if (i .ne. 0) call rx(' matrix is singular')
        if (cs(1:1) .eq. 't') then
          do  14  j = 1, nb
          do  14  i = 1, n
   14     pb(i,j) = sb(j,i)
          call dgetrs('T',n,nb,sa,lda,piv,pb,ldb,i)
          if (i .ne. 0) call rx('bad call to dgetrs')
          do  18  j = 1, nb
          do  18  i = 1, n
            sb(j,i) = pb(i,j)
   18     continue
        else
          call dgetrs('N',n,nb,sa,lda,piv,sb,ldb,i)
          if (i .ne. 0) call rx('bad call to dgetrs')
        endif
        strn = 'dgetrf,rs'
      endif
      s_time = cpusec() - temp
      if (s_time .eq. 0d0) s_time = 1
      if (s_time .ne. 0) then
        s_mflops = mflops / s_time
      else
        s_mflops = 0
      endif

C --- Inversion-and-backsubstitution by dqinvb ---
      call ran1in(1)
      CALL init ( sa,lda,nnb, cs )
      CALL init ( pb,ldb,nnb, cs )
C ... (debugging) check usage of w
      ldw = min(n+5,lda)
      print *, 'using lda,ldw,n=',lda,ldw,n
      call defrr(opw,  ldw*n)
      call defrr(opw2, n*n)
      call dvset(w(opw),1,ldw*n,-99d0)
      call dvset(w(opw2),1,n*n,-99d0)
      call dvset(pw2,1,ldb*(ldb+1),-99d0)

      temp = cpusec()
      call dqinvb(cs,sa,lda,nlev,n,nb,pw,ldw,pw2,pb,ldb,i)
      q_time = cpusec() - temp
      if (q_time .eq. 0d0) q_time = 1
      q_mflops = mflops / q_time
      if (i .ne. 0) call rx('dqinvb failed to find inverse')
*     call ywrm(0,'dqinvb a^-1 b',1,6,'(4f16.10)',pb,ldb,n,nb)

C --- Check quality of inversion ---
      if (cs(1:1) .ne. 't') then
        CALL diff('compare '//strn//'to dqinvb', sb, pb, ldb, n, nb)
      else
        CALL diff('compare '//strn//'to dqinvb', sb, pb, ldb, nb, n)
      endif
      call ran1in(1)
      CALL init ( sa,lda,nnb, cs )
      if (cs(1:1) .ne. 't') then
        call dgemm('N','N',n,nb,n,1d0,sa,lda,pb,ldb,0d0,sb,ldb)
        CALL init ( pb,ldb,nnb, cs )
        CALL diff('compare a (a^-1 b) to b', sb, pb, ldb, n, nb)
      else
        call dgemm('N','N',nb,n,n,1d0,pb,ldb,sa,lda,0d0,sb,ldb)
        CALL init ( pb,ldb, nnb, cs )
        CALL diff('compare (b a^-1) a to b', sb, pb, ldb, nb, n)
      endif

      print 110, s_time, s_mflops, q_time, q_mflops, s_time/q_time,n,nb

  110 format(/1X,'Serial time: ',F7.3,'  Serial MFlops: ',F6.1,
     .       /1X,'dqinvb time: ',F7.3,'  dqinvb MFlops: ',F6.1,
     .       /1X,'     factor: ',F7.3,'  for n,nb =',2i4)

C ... Make sure pw not overwritten beyond block(1..n,1..n+1)
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

C ... Make sure pw2 not overwritten beyond (n+1)*nb
      imax = ldb*(ldb+1)
      do  60  i = (n+1)*nb+1, imax
        if (pw2(i) .ne. -99d0) call fexit2(-1,111,
     .    ' STOP: pw2 overwritten, i,j=%i %i',i,j)
   60 continue

*     call ywrm(0,'pw',1,6,'(16f8.2)',pw,ldw,n+1,n)

C ... Check that 'b' option works
      call ran1in(1)
      CALL init ( sa,lda,nnb, cs )
      CALL init ( pb,ldb,nnb, cs )
      call dcopy(ldb*nnb,pb,1,sb,1)
      call dqinvb(cs,sa,lda,nlev,n,nb,pw,ldw,pw2,sb,ldb,i)

      temp = cpusec()
      call word(cs,1,j1,j2)
      cs(max(j2,0)+1:max(j2,0)+1) = 'b'
      call dqinvb(cs,sa,lda,nlev,n,nb,pw,ldw,pw2,pb,ldb,i)
      xtim = cpusec() - temp
      if (xtim .eq. 0d0) xtim = 1

      print *, ' '
      if (cs(1:1) .ne. 't') then
        CALL diff('compare dqinvb to same w/b', sb, pb, ldb, n, nb)
      else
        CALL diff('compare dqinvb to same w/b', sb, pb, ldb, nb, n)
      endif

      print 111, xtim
  111 format(/1X,'     B time: ',F7.3,'       B MFlops: ',F6.1)

      end

      SUBROUTINE init ( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double precision a( lda, * )
      character *(*) cs
      real ran1

      do 10 i = 1, n
      do 10 j = 1, n
   10 a(i,j) = ran1()

      if (cs .eq. 's' .or. cs .eq. 'ts' .or.
     .    cs .eq. 'rs' .or. cs .eq. 'trs') then
C       print *, '... symmetrizing a'
        do 20 i = 1, n
        do 20 j = 1, i
   20   a(i,j) = a(j,i)
      endif

      end
      SUBROUTINE diff (strn, sc, pc, ldc, n, nb)
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      integer ldc, n, m, i, j, nb, ii,jj
      double precision sc( ldc, * ), pc( ldc, * ), errmx
      errmx = 0d0
      ii = 0
      jj = 0
      do 10 i = 1, n
      do 10 j = 1, nb
        if (dabs(sc(i,j)-pc(i,j)) .gt. errmx) then
          ii = i
          jj = j
        endif
   10 errmx = max(errmx,dabs(sc(i,j)-pc(i,j)))
      print 100, strn,n,nb,errmx,ii,jj
  100 format(1X,a,t30,'(',i4,' rows,',i4,' columns):  errmx=',g10.3,2i4)
      end
      subroutine ywrm(lbin,filel,icast,ifi,fmt,s,ns,nr,nc)
C- Writes complex matrix to file ifi
Cr lbin: writes in binary mode
      logical lbin
      character*(*) filel
      integer icast,nr,nc
      double precision s(ns,nc,2)
      character*(*) fmt, outs*10
      integer i,j

      if (lbin) then
        if (filel .eq. ' ') then
          write(ifi) nr,nc,icast
        else
          write(ifi) nr,nc,icast,len(filel)
          write(ifi) filel(1:len(filel))
        endif
        call dpdump(s,nr*nc*mod(icast,10),-ifi)
        return
      endif
      outs = ' '
      if (icast .eq. 1)  outs = ' real'
      if (icast .eq. 11) outs = ' symm'
      if (icast .eq. 2)  outs = ' complex'
      if (icast .eq. 12) outs = ' herm'
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)
c      if (filel .ne. ' ') call awrit0(filel,' ',len(filel),ifi)
      if (filel .ne. ' ') print 333, filel
  333 format('#',a)
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(i,j,1), j=1,nc)
      if (mod(icast,10) .gt. 1) then
       write(ifi,'(1x)')
      do  20  i = 1, nr
   20 write(ifi,fmt) (s(i,j,2), j=1,nc)
      endif
      end
      subroutine prm(a,lda,n1,n2)
      integer lda,n1,n2
      double precision a(lda,lda)
      call ywrm(.false.,' ',1,6,'(12f12.6)',a,lda,n1,n2)
      end
