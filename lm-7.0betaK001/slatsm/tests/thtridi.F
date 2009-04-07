C test timing of htridi
      subroutine fmain
      implicit none
      INTEGER   lda, ldb, ldc, ldw, n, i, j, i0, j0,
     .  opw, opw2, wksize
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,wksize=6000000)
      character*1 cs, outs*80
      double precision h(lda,lda,2), wk(lda,10), e(lda,2)
      double precision s_time, temp, s_mflops, q_time, q_mflops
      double precision z_time, z_mflops
      double precision cpusec,xx,errmx,abstol,dlamch
      integer ierr,nb,ilaenv,lwork,owork,orwork,oiwork,m,ifail,info
      integer w(wksize)
      common /w/ w
      common /static/ h

C --- Setup ---
      call finits(2,0,0,i)
      call wkinit(lda*lda*4)
      call nada
      n = 300
      call initqu(.true.)
      call query('n=?',2,n)
      if (n .gt. lda) call rx('n gt lda')
      mflops = 1.0d-6 * n * n * n * 16d0/3

C --- Eigenvalues and eigenvectors, htridi ---
      CALL init (h,lda,n,'h')
c     call yprm(.false.,'a',2,6,'(8f16.10)',h,lda,n,lda,n)
      temp = cpusec()
      call htridx(lda,n,h,h(1,1,2),e(1,1),wk,wk(1,4),wk(1,2))
      do  10  j = 1, n
   10 wk(j,1) = wk(j,1)**2
      call tqlrat(n,e(1,1),wk,ierr)
      call rxx(ierr.ne.0,'tqlrat cannot find all evals') 
      q_time = max(cpusec() - temp,1d-6)
      q_mflops = mflops / q_time

      CALL init (h,lda,n,'h')
      temp = cpusec()
      call htridi(lda,n,h,h(1,1,2),e(1,2),wk,wk,wk(1,2))
      do  20  j = 1, n
   20 wk(j,1) = wk(j,1)**2
      call tqlrat(n,e(1,2),wk,ierr)
      s_time = max(cpusec() - temp,1d-6)
      s_mflops = mflops / s_time

      print *, 'lowest, highest evals'
      do  30  j = 1, n, n-1
   30 print 345, e(j,1), e(j,2), e(j,1)-e(j,2)
  345 format(2f20.10,1pe9.1)
      call diff('compare evals',e,e(1,2),lda,n)

      print 110, s_time, s_mflops, q_time, q_mflops, s_time/q_time,n
  110 format(/1X,'htridi time: ',F7.3,'  htridi MFlops: ',F6.1,
     .       /1X,'htridx time: ',F7.3,'  htridx MFlops: ',F6.1
     .       /1X,'     factor: ',F7.3,'  for n =',i4)

C --- Lapack  ---
C#ifdef LAPACK
      print *
      print *, '... compare to lapack'
      CALL initz (h,lda,n,'h')
C     call zprm('h',h,lda,n,n)
      temp = cpusec()
      abstol=2*DLAMCH('S')
      nb = ilaenv(1, 'ZHEEVX', 'VIU', n, -1, -1, -1)
      lwork = (nb+1)*n
      call defcc(owork, lwork)
      call defrr(orwork, 7*n)
      call defi(oiwork, 5*n)
      call dpzero(e(1,2),n)
      call ZHEEVX( 'N', 'A','U', n, h, lda,
     .              0d0, 0d0, ! not used
     .              1, n, ! which eigenvalues to calculate
     .              abstol,
     .              m, e(1,2),
     .              w,   n, ! eigenvectors T(n,n)
     .              w(owork), lwork, w(orwork),
     .              w(oiwork), ifail, info )
      z_time = cpusec() - temp
      z_mflops = mflops / z_time
      print *, 'lowest, highest evals'
      do  40  j = 1, n, n-1
   40 print 345, e(j,1), e(j,2), e(j,1)-e(j,2)
      call diff('compare evals',e,e(1,2),lda,n)

      print 120, z_time, z_mflops, q_time, q_mflops,
     .  z_time/q_time,n
  120 format(/1X,'zheevx time: ',F7.3,'  zheevx MFlops: ',F6.1,
     .       /1X,'htridx time: ',F7.3,'  htridx MFlops: ',F6.1,
     .       /1X,'     factor: ',F7.3,'  for n =',i4)
      return
C#endif

      end

      SUBROUTINE init ( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double precision a( lda, lda, 2)
      character *1 cs
      real ran1
      call ran1in(1)

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
      SUBROUTINE initz ( a, lda, n, cs)
C- Initialize arrays
      integer lda, n
      double complex a( lda, lda)
      double precision xx1,xx2
      character *1 cs
      real ran1
      call ran1in(1)

      do 10 i = 1, n
      do 10 j = 1, n
      xx1 = ran1()
      xx2 = ran1()
   10 a(i,j) = dcmplx(xx1,xx2)

      if (cs .eq. 'h') then
        do 20 i = 1, n
        do 20 j = 1, i
   20   a(i,j) =  dconjg(a(j,i))
        do  22 i = 1, n
   22   a(i,i) = dble(a(i,i))
      endif
      end
      SUBROUTINE diff (strn, sc, pc, ldc, n)
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      integer ldc, n, i, j
      double precision sc(ldc,ldc,2), pc(ldc,ldc,2), errmx
      errmx = 0d0
      do 10 j = 1, 1
      do 10 i = 1, n
   10 errmx = max(errmx,dabs(sc(i,j,1)-pc(i,j,1)))
      print 100, strn,errmx
  100 format(1X,a,':  errmx=',1pe9.2)
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
      subroutine zprm(strn,s,ns,nr,nc)
      logical lbin
      integer icast,nr,nc
      double precision s(2,ns,nc)
      character*(10) fmt, strn*(*), outs*80
      integer i,j

      fmt = '(5f18.13)'
      ifi = 19
      open(ifi,file='out')
      rewind ifi
      write(ifi,234) nr, nc
  234 format('% rows',i5,' cols',i5,' complex')
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(1,i,j), j=1,nc)
      do  20  i = 1, nr
   20 write(ifi,fmt) (s(2,i,j), j=1,nc)
      close(ifi)
      outs = 'zprm: done writing data '//strn
      print *, outs
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx('quit in zprm')
      end

C# speedup factor on the SGI:
C# n    cs=n  cs=h
C  150  3.758 2.643
C  151  4.085 2.995
C  152  3.778 2.694
C  153  4.109 3.028
C  154  3.787 2.734
C  155  4.140 2.952
C  156  3.827 2.807
C  157  4.184 2.990
C  158  3.750 2.754
C  159  4.159 3.005
C  160  3.995 2.796
C  161  4.207 3.087
C  162  3.885 2.704
C  163  4.227 3.013
C  164  3.879 2.743
C  165  4.176 3.062
C  166  3.907 2.723
C  167  4.356 3.085
C  168  4.009 2.810
C  169  4.323 3.073
C  170  3.975 2.717
C  171  4.369 3.076
C  172  3.962 2.788
C  173  4.306 3.110
C  174  3.983 2.786
C  175  4.390 3.110
C  176  4.045 2.813
C  177  4.358 3.154
C  178  4.037 2.807
C  179  4.421 3.152
C  180  3.990 2.846
C  181  4.395 3.149
C  182  4.061 2.786
C  183  4.499 3.166
C  184  4.108 2.850
C  185  4.534 3.184
C  186  4.036 2.809
C  187  4.446 3.158
C  188  4.094 2.826
