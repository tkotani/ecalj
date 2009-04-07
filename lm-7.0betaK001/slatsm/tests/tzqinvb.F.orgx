C Tests zqinvb
C To test each branch, use:
C
C  rm out
C  foreach cs (1 t th th4 th4l h h4 h4l 4 4l )
C    echo $cs'\n\n' | a.out >> out
C  end
      subroutine fmain
      implicit none
      INTEGER   lda, ldb, ldc, ldw, n, nb, i, j, imax, jmax,
     .  opw, opw2, wksize, j1,j2, ifi,fopng,rdm
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,wksize=6000000)
      character cs*5,outs*80
      double precision sa(lda,lda,2),sb(ldb,ldb,2),
     .  pb(ldb,ldb,2),pw(lda,lda,2),pw2(ldb*ldb+ldb)
      double precision stime, temp, smflops, qtime, qmflops
      double precision cpusec,xx,xtim,xflops
      integer w(wksize)
      common /w/ w
      common /static/ sa, sb, pb, pw, pw2

C --- Setup ---
      call finits(2,0,0,i)
      call initqu(.true.)
      call wkinit(lda*lda*4)
C#ifdefC RDA
C      cs = ' '
C      n = -1
C      nb = n
C      cs = 't'
C#else
      n = 300
      cs = 'th'
      print *, 'cs= (combination of t, h, 4, l)?'
      read(*,'(a4)') cs
      call query('n=?',2,n)
      nb = max(n-1,3)
      call query('nb=?',2,nb)
C#endif

       print *, 'test zqinvb: cs=',cs
      if (n .gt. lda .or. nb .gt. lda) call rx('n gt lda')

C --- Matrix inversion, yygefa,yygesl ---
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C      nb = 5
C#endif
      CALL init (sa,lda,n,cs)
      CALL init2(sb,ldb,max(n,nb),cs)
      mflops = 2.0d-6 * n * n * n * 4

C     call yprm(.false.,'a',2,6,'(8f16.10)',sa,lda,n,lda,n)
C     call yprm(.false.,'b',2,6,'(8f16.10)',sb,ldb,n,ldb,nb)
C     call yprm(.false.,'bt',2,6,'(8f16.10)',sb,ldb,nb,ldb,n)
      temp = cpusec()
      if (cs .eq. 'h' .and. .false.) then
        call yyhifa(sa,sa(1,1,2),lda,n,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')
        stop 'not ready for ax=b, symmetric case'
      else
        if (cs .eq. 'h') print *, '(NB): not using symmetric yyhifa'
        call yygefa(sa,sa(1,1,2),lda,n,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')
        if (cs(1:1) .eq. 't') then
          do  14  j = 1, nb
            do  16  i = 1, n
            pb(i,j,1) = sb(j,i,1)
   16       pb(i,j,2) = sb(j,i,2)
            call yygesl(sa,sa(1,1,2),lda,n,pw,pb(1,j,1),pb(1,j,2),1)
            do  18  i = 1, n
            sb(j,i,1) = pb(i,j,1)
            sb(j,i,2) = pb(i,j,2)
   18       continue
   14     continue
          stime = cpusec() - temp
*         call yprm(.false.,'b a^-1',2,6,'(5f16.10)',sb,ldb,nb,lda,n)
        else
          do  12  j = 1, nb
   12     call yygesl(sa,sa(1,1,2),lda,n,pw,sb(1,j,1),sb(1,j,2),0)
          stime = cpusec() - temp
*         call yprm(.false.,'a^-1 b',2,6,'(5f16.10)',sb,ldb,n,lda,nb)
        endif
      endif
      if (stime .ne. 0) then
        smflops = mflops / stime
      else
        smflops = 1
      endif

C --- Inversion-and-backsubstitution by zqinvb ---
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init (sa,lda,n,cs)
      CALL init2(pb,ldb,max(n,nb),cs)
C     call yprm(.false.,'a',2,6,'(10f15.10)',sa,lda,n,lda,n)
C     call yprm(.false.,'b',2,6,'(10f15.10)',pb,ldb,n,ldb,nb)
C ... (debugging) check usage of w
      ldw = min(n+5,lda)
C     lapack uses data as a 1d array; must use ldw=n to check w
      ldw = n
      write(*,'('' using lda,ldw,n,nb='',4i4)') lda,ldw,n,nb
      call defrr(opw,  ldw*(n+1))
      call defrr(opw2, n*n)
      call dvset(w(opw),1,ldw*(n+1),-99d0)
      call dvset(w(opw2),1,n*n,-99d0)
      call dvset(pw2,1,ldb*(ldb+1),-99d0)

      call ztoyy(sa,lda,lda,n,n,0,1)
      if (cs(1:1) .ne. 't') call ztoyy(pb,ldb,ldb,n,nb,0,1)
      if (cs(1:1) .eq. 't') call ztoyy(pb,ldb,ldb,nb,n,0,1)
C     call zprm('a',2,sa,lda,n,n)
C     if (cs(1:1) .ne. 't') call zprm('b',2,pb,ldb,n,nb)
C     if (cs(1:1) .eq. 't') call zprm('bt',2,pb,ldb,nb,n)

      temp = cpusec()
      call zqinvb(cs,sa,lda,n,nb,w(opw),ldw,pw2,pb,ldb,i)
C      if (cs(1:1) .ne. 't') call zprm('a^-1 b',2,pb,ldb,n,nb)
C      if (cs(1:1) .eq. 't') call zprm('b a^-1',2,pb,ldb,nb,n)

      qtime = cpusec() - temp
      if (cs(1:1) .ne. 't') call ztoyy(pb,ldb,ldb,n,nb,1,0)
      if (cs(1:1) .eq. 't') call ztoyy(pb,ldb,ldb,nb,n,1,0)

      if (qtime .ne. 0) then
        qmflops = mflops / qtime
      else
        qmflops = 1
      endif
      if (i .ne. 0) call rx('zqinvb failed to find inverse')
*     call yprm(.false.,'a^-1 b',2,6,'(5f16.10)',pb,ldb,n,lda,nb)

C --- Check quality of inversion ---
      if (cs(1:1) .ne. 't') then
        CALL diff('compare yygefa,sl to zqinvb', sb, pb, ldb, n, nb)
      else
        CALL diff('compare yygefa,sl to zqinvb', sb, pb, ldb, nb, n)
      endif
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init( sa,lda,n,cs)
      if (cs(1:1) .ne. 't') then
        call yygemm('N','N',n,nb,n,1d0,sa,sa(1,1,2),lda,pb,pb(1,1,2),
     .  ldb,0d0,sb,sb(1,1,2),ldb)
        CALL init2( pb,ldb,max(n,nb),cs)
        CALL diff('compare a (a^-1 b) to b', sb, pb, ldb, n, nb)
      else
        call yygemm('N','N',nb,n,n,1d0,pb,pb(1,1,2),ldb,
     .    sa,sa(1,1,2),lda,0d0,sb,sb(1,1,2),ldb)
        CALL init2(pb,ldb, max(n,nb),cs)
        CALL diff('compare (b a^-1) a to b', sb, pb, ldb, nb, n)
      endif

      print 110, stime, smflops, qtime, qmflops, qmflops/smflops,n,nb

  110 format(/1X,'Serial time: ',F7.3,'  Serial MFlops: ',F6.1,
     .       /1X,'zqinvb time: ',F7.3,'  zqinvb MFlops: ',F6.1,
     .       /1X,'     factor: ',F7.3,'  for n,nb =',2i4)

C ... Make sure pw not overwritten beyond ldw*(n+1)
      call wkchk('check if pw2 still intact')
      call dcopy(n*n,w(opw2),1,pw,1)
      j = 1
      do  50  i = 1, n*n
        if (pw(i,j,1) .ne. -99d0) call fexit2(-1,111,
     .    ' STOP: pw2 overwritten, i,j=%i %i',i,j)
   50 continue
      if (ldw .gt. n .or. .true.) then
        call wkchk('check if pw still intact')
        call dmcpy(w(opw),ldw,1,pw,lda,1,ldw,n+1)
        do  52  i = 1, ldw
        do  52  j = 1, n+1
            if (i .le. n .and. j .le. n+1) goto 52
            if (pw(i,j,1) .ne. -99d0) call fexit2(-1,111,
     .        ' STOP: pw overwritten, i,j=%i %i',i,j)
   52   continue
      endif

C ... Make sure pw2 not overwritten beyond (n+1)*nb
      imax = ldb*(ldb+1)
      do  60  i = (n+1)*nb+1, imax
        if (pw2(i) .ne. -99d0) call fexit2(-1,111,
     .    ' STOP: pw2 overwritten, i=%i',i,j)
   60 continue

C ... Check that 'b' option works
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init (sa,lda,n,cs)
      CALL init2(pb,ldb,max(n,nb),cs)
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init (sa,lda,n,cs)
      CALL init2(sb,ldb,max(n,nb),cs)
      call dcopy(ldb*max(n,nb),pb,1,sb,1)
      call ztoyy(sa,lda,lda,n,n,0,1)
      if (cs(1:1) .ne. 't') call ztoyy(sb,ldb,ldb,n,nb,0,1)
      if (cs(1:1) .eq. 't') call ztoyy(sb,ldb,ldb,nb,n,0,1)

C     call zprm('a',2,sa,lda,n,n)
C     if (cs(1:1) .ne. 't') call zprm('b',2,sb,ldb,n,nb)
C     if (cs(1:1) .eq. 't') call zprm('bt',2,sb,ldb,nb,n)
      call zqinvb(cs,sa,lda,n,nb,w(opw),ldw,pw2,sb,ldb,i)
C     if (cs(1:1) .ne. 't') call zprm('a^-1 b',2,sb,ldb,n,nb)
C     if (cs(1:1) .eq. 't') call zprm('b a^-1',2,sb,ldb,nb,n)
      if (cs(1:1) .ne. 't') call ztoyy(sb,ldb,ldb,n,nb,1,0)
      if (cs(1:1) .eq. 't') call ztoyy(sb,ldb,ldb,nb,n,1,0)


      if (cs(1:1) .ne. 't') call ztoyy(pb,ldb,ldb,n,nb,0,1)
      if (cs(1:1) .eq. 't') call ztoyy(pb,ldb,ldb,nb,n,0,1)
      temp = cpusec()
      call word(cs,1,j1,j2)
      cs(max(j2,0)+1:max(j2,0)+1) = 'b'
      call zqinvb(cs,sa,lda,n,nb,w(opw),ldw,pw2,pb,ldb,i)
      xtim = cpusec() - temp
      xflops = mflops / max(xtim,1d-6)
      if (cs(1:1) .ne. 't') call ztoyy(pb,ldb,ldb,n,nb,1,0)
      if (cs(1:1) .eq. 't') call ztoyy(pb,ldb,ldb,nb,n,1,0)

      print *, ' '
      if (cs(1:1) .ne. 't') then
        CALL diff('compare zqinvb to same w/b', sb, pb, ldb, n, nb)
      else
        CALL diff('compare zqinvb to same w/b', sb, pb, ldb, nb, n)
      endif

      print 111, xtim, xflops
  111 format(/1X,'     B time: ',F7.3,'       B MFlops: ',F6.1)

      end

      SUBROUTINE init(a,lda,n,cs)
C- Initialize arrays
      implicit none
      integer lda,n
      double precision a( lda, lda, 2)
      character*(*) cs
      integer i,j
      real ran1
C#ifdefC RDA
C      double precision pb(128*128*4)
C      integer rdm,ifi,fopng
C
C      if (n .eq. -1) then
C      n = 128
C      ifi = fopng('a',-1,0)
C      rewind ifi
C      i = rdm(ifi,1011,128*128*4,' ',pb,n,n)
C      if (i .le. 0) call rx('failed to read a from file')
C      call ymcpy(pb,n,1,n*n,a,lda,1,lda*lda,n,n)
C      return
C      endif
C#endif

      do 10 i = 1, n
      do 10 j = 1, n
      a(i,j,1) = ran1()
   10 a(i,j,2) = ran1()

      i = 0
      call chrpos(cs,'h',len(cs),i)
      if (i .lt. len(cs)) then
        do 20 i = 1, n
        do 20 j = 1, i
        a(i,j,1) =  a(j,i,1)
   20   a(i,j,2) = -a(j,i,2)
        do  22 i = 1, n
   22   a(i,i,2) = 0
      endif

      end
      SUBROUTINE init2( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double precision a( lda, lda, 2)
      character*(*) cs

      do 10 i = 1, n
      do 10 j = 1, n
      a(i,j,1) = ran1()
   10 a(i,j,2) = ran1()

      i = 0
      call chrpos(cs,'h',len(cs),i)
      if (i .lt. len(cs)) then
        do 20 i = 1, n
        do 20 j = 1, i
        a(i,j,1) =  a(j,i,1)
   20   a(i,j,2) = -a(j,i,2)
        do  22 i = 1, n
   22   a(i,i,2) = 0
      endif

      end
      SUBROUTINE diff (strn, sc, pc, ldc, n, nb)
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      integer ldc, n, m, i, j, nb
      double precision sc(ldc,ldc,2), pc(ldc,ldc,2), errmx
      errmx = 0d0
      do 10 j = 1, nb
      do 10 i = 1, n
   10 errmx = max(errmx,dabs(sc(i,j,1)-pc(i,j,1)),
     .                  dabs(sc(i,j,2)-pc(i,j,2)))
      print 100, strn,n,nb,errmx
  100 format(1X,a,t30,'(',i4,' rows,',i4,' columns):  errmx=',g10.3)
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
      subroutine yyprm(strn,icast,s,ofi,ns,nr,nc)
C ofi used only for kcplx=0
C ns,nr,nc are formal dimensions, not real ones
      implicit none
      integer icast,ofi,ns,nr,nc,ifi
*      double precision s(ns,nsc,2)
      double precision s(ns,nc)
      character*(20) fmt, fmt0, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      save fmt
      data fmt /'(9f15.10)'/
*     fmt = '(9f15.10)'
*     fmt = '(%9;10,10d)'
      ifi = fopna('out',29,0)

      if (icast .ne. 0) then
        call ywrm(0,' ',icast,ifi,fmt,s,ofi,ns,nr,nc)
      else
        call awrit2('%% rows %i cols %i',' ',80,ifi,nr,nc)
        do  10  i = 1, nr
   10   write(ifi,'(22i7)') (s(i,j), j=1,nc)
      endif

      call fclose(ifi)
      outs = ' prm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in prmx')
      return
      entry yprm0(fmt0)
      fmt = fmt0
      end
      subroutine zprm(strn,icast,s,ns,nr,nc)
      implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(10) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      fmt = '(9f15.10)'
      fmt = '(5f20.15)'
      outs = ' '
      if (icast .eq. 1)  outs = ' real'
      if (icast .eq. 11) outs = ' symm'
      if (icast .eq. 2)  outs = ' complex'
      if (icast .eq. 12) outs = ' herm'
      ifi = fopna('out',29,0)
      rewind ifi
C#ifdefC NOAWRITE
C      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
C#else
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)
C#endif
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(1,i,j), j=1,nc)
      if (mod(icast,10) .gt. 1) then
      write(ifi,*)
      do  20  i = 1, nr
   20 write(ifi,fmt) (s(2,i,j), j=1,nc)
      endif
      call fclose(ifi)
C#ifdefC NOAWRITE
C      outs = ' prm: wrote '//strn//' continue?'
C      print *, outs
C#else
      outs = ' zprm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
C#endif
      read(*,'(a80)') outs

C#ifdefC NOAWRITE
C      if (outs .eq. 'q') stop 'quit in zprm'
C#else
      if (outs .eq. 'q') call rx0('quit in zprm')
C#endif
      end
