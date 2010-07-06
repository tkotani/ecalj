C test yyqinv
      subroutine fmain
      implicit none
      INTEGER   lda, ldb, ldc, ldw, n, nb, i, j, i0, j0
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,ldw=901)
      character*1 cs
      double precision sa(lda,lda,2),sb(ldb,ldb,2),
     .  pb(ldb,ldb,2),pw(lda,lda,2)
      double precision s_time, temp, s_mflops, q_time, q_mflops
      double precision cpusec,xx,errmx
      common /static/ sa, sb, pb, pw

C --- Setup ---
      call finits(2,0,0,i)
      n = 5
      cs = 'n'
      call initqu(.true.)
      print *, 'test yyqnvb: cs=',cs
      call query('n=?',2,n)
      nb = n-1
      call query('nb=?',2,nb)
      if (n .gt. lda .or. nb .gt. lda) call rx('n gt lda')
      mflops = 2.0d-6 * n * n * n * 4
C ... (debugging) check usage of w
      call dvset(pw,1,ldw*ldw*2,-99d0)

C --- Matrix inversion, dgefa,dgedi ---
      call ran1in(1)
      CALL init ( sa,lda,n, cs )
      CALL init ( sb,ldb,n, cs )

*     call yprm(.false.,'a',2,6,'(8f16.10)',sa,lda,n,lda,n)
*     call yprm(.false.,'b',2,6,'(8f16.10)',sb,ldb,n,ldb,n)
      temp = cpusec()
      if (cs .eq. 'h') then
        call yyhifa(sa,sa(1,1,2),lda,n,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')
        stop 'not ready for ax=b'
        call yyhidi(sa,sa(1,1,2),lda,n,pw,xx,i,pw(1,2,1),1)
        do  16  i = 1, n
        do  16  j = 1, i
        sa(i,j,1) =  sa(j,i,1)
   16   sa(i,j,2) = -sa(j,i,2)
      else
        call yygefa(sa,sa(1,1,2),lda,n,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')
        do  12  j = 1, nb
   12   call yygesl(sa,sa(1,1,2),lda,n,pw,sb(1,j,1),sb(1,j,2),0)
      endif
      s_time = cpusec() - temp
      s_mflops = mflops / s_time
*     call yprm(.false.,'a^-1 b',2,6,'(5f16.10)',sb,ldb,n,lda,nb)

C --- Check a (a^-1 b) = b ---
      call ran1in(1)
      CALL init ( sa,lda,n, cs )
      call yygemm('N','N',n,n,n,1d0,sa,sa(1,1,2),lda,sb,sb(1,1,2),ldb,
     .  0d0,pw,pw(1,1,2),ldw)
*     call yprm(.false.,'a (a^-1 b)',2,6,'(8f16.10)',pw,ldw,n,ldw,n)
      CALL init ( sb,ldb,n, cs )

      do  40  j = 1, nb
      do  40  i = 1, n
      pw(i,j,1) = pw(i,j,1) - sb(i,j,1)
   40 pw(i,j,2) = pw(i,j,2) - sb(i,j,2)
      errmx = 0
      do  42  j = 1, nb
      do  42  i = 1, n
   42 errmx = max(errmx,dabs(pw(i,j,1)),dabs(pw(i,j,2)))
      print 100, 'a (inv(a) b) - b',errmx,nint(s_mflops)
  100 format(1X,a,':  errmx=',g10.3,2i4)

      end

      SUBROUTINE init ( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double precision a( lda, lda, 2)
      character *1 cs

      do 10 i = 1, n
      do 10 j = 1, n
      a(i,j,1) = ran1()
   10 a(i,j,2) = ran1()

      if (cs .eq. 'h') then
        do 20 i = 1, n
        do 20 j = 1, i
        a(i,j,1) =  a(j,i,1)
   20   a(i,j,2) = -a(j,i,2)
        do  22 i = 1, n
   22   a(i,i,2) = 0
      endif

      end
      SUBROUTINE diff (strn, sc, pc, ldc, n)
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      integer ldc, n, m, i, j
      double precision sc(ldc,ldc,2), pc(ldc,ldc,2), errmx
      errmx = 0d0
      do 10 i = 1, n
      do 10 j = 1, n
   10 errmx = max(errmx,dabs(sc(i,j,1)-pc(i,j,1)),
     .                  dabs(sc(i,j,2)-pc(i,j,2)))
      print 100, strn,errmx
  100 format(1X,a,':  errmx=',g10.3)
      end
      subroutine yprm(lbin,filel,icast,ifi,fmt,s,ns,nr,nsc,nc)
C- Writes complex matrix to file ifi
Cr lbin: writes in binary mode
      logical lbin
      character*(*) filel
      integer icast,nr,nc
      double precision s(ns,nsc,2)
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
