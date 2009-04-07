C test yygemm
      subroutine fmain
      implicit none
      INTEGER   lda, ldb, ldc, ldw, n, nb, nl, i, j
      double precision mflops
      PARAMETER (lda=901,ldb=901,ldc=901,ldw=901)
      character*1 csa,csb,outs*80
      double precision sa(lda,lda,2),sb(ldb,ldb,2),
     .  pb(ldb,ldb,2),pw(lda,lda,2)
      double precision s_time, temp, s_mflops, q_time, q_mflops
      double precision cpusec,xx,errmx
      common /static/ sa, sb, pb, pw

C --- Setup ---
      call finits(2,0,0,i)
      n = 300
      csa = 'N'
      csb = 'N'
      outs(1:2) = '  '
      print *, 'csa,csb= (two char, eg TC or CN)'
      read(*,'(a2)') outs(1:2)
      if (outs(1:2) .eq. '  ') outs(1:2) = 'NN'
      csa = outs(1:1)
      csb = outs(2:2)
      call initqu(.true.)
      print *, 'test yygemm: csa,csb=',csa,csb
      call query('na=?',2,n)
      nb = n
      call query('nb=?',2,nb)
      nl = n-2
      call query('nl=?',2,nl)
      if (n .gt. lda .or. nb .gt. lda) call rx('n gt lda')
      mflops = 2.0d-6 * n * nb * n * 4

C --- Reference : yygemm0 ---
      call ran1in(1)
      CALL init ( sa,lda,max(n,nl), 'n' )
      CALL init ( sb,ldb,max(nb,nl), 'n' )
C      call yprm(.false.,'a',2,6,'(8f16.10)',sa,lda,n,lda,nl)
C      call yprm(.false.,'b',2,6,'(8f16.10)',sb,ldb,nl,ldb,nb)

      temp = cpusec()
      call yygemm0(csa,csb,n,nb,nl,1d0,sa,sa(1,1,2),lda,sb,sb(1,1,2),
     .  ldb,0d0,pw,pw(1,1,2),ldw)
      s_time = cpusec() - temp
      s_mflops = mflops / s_time
C      call yprm(.false.,'a*b',2,6,'(8f16.10)',pw,ldw,n,ldw,nb)

      temp = cpusec()
      call yygemm(csa,csb,n,nb,nl,1d0,sa,sa(1,1,2),lda,sb,sb(1,1,2),ldb,
     .  0d0,pb,pb(1,1,2),ldb)
C      call yprm(.false.,'a*b',2,6,'(8f16.10)',pb,ldb,n,ldb,nb)
      q_time = cpusec() - temp
      q_mflops = mflops / q_time

      CALL diff('compare yygemm0 to  yygemm', pw, pb, ldb, n, nb)

      print 110, s_time, s_mflops, q_time, q_mflops

  110 format(/1X,'yygemm0 time: ',F7.3,'  yygemm0 MFlops: ',F6.1,
     .       /1X,'yygemm  time: ',F7.3,'  yygemm  MFlops: ',F6.1)


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
      SUBROUTINE diff (strn, sc, pc, ldc, n, nb)
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      integer ldc, n, m, i, j, nb
      double precision sc(ldc,ldc,2), pc(ldc,ldc,2), errmx
      errmx = 0d0
      do  10  j = 1, nb
      do  10  i = 1, n
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

      subroutine yygemm0(transa,transb,m,n,k,alpha,ar,ai,lda,br,bi,ldb,
     .  beta,cr,ci,ldc)
C- Analog of zgemm, using real arithmetic
C ----------------------------------------------------------------
Ci Inputs:
Ci   ar,ai; br,bi; cr,ci: analog of a,b,c in dgemm
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      double precision   alpha, beta
      double precision   ar(lda,1), br(ldb,1), cr(ldc,1),
     .                   ai(lda,1), bi(ldb,1), ci(ldc,1)
C Local variables
      logical ls
      integer im,in
      double precision xx,s

      s = alpha

C --- Real part ---
      ls = transa .eq. 'C' .eqv. transb .eq. 'C'
      if (ls) then
        do  10  in = 1, n
        do  10  im = 1, m
   10   cr(im,in) = -cr(im,in)
      endif
      call dgemm(transa,transb,m,n,k,s,ai,lda,bi,ldb,beta,cr,ldc)

      if (ls) then
        do  12  in = 1, n
        do  12  im = 1, m
   12   cr(im,in) = -cr(im,in)
      endif
      call dgemm(transa,transb,m,n,k,s,ar,lda,br,ldb,1d0,cr,ldc)

C --  Imaginary part ---
      ls = transb .eq. 'C'
      if (ls) then
        do  14  in = 1, n
        do  14  im = 1, m
   14   ci(im,in) = -ci(im,in)
      endif

      call dgemm(transa,transb,m,n,k,s,ar,lda,bi,ldb,beta,ci,ldc)

      if (ls .neqv. transa .eq. 'C') then
        do  16  in = 1, n
        do  16  im = 1, m
   16   ci(im,in) = -ci(im,in)
      endif
      call dgemm(transa,transb,m,n,k,s,ai,lda,br,ldb,1d0,ci,ldc)

      if (transa .eq. 'C') then
        do  18  in = 1, n
        do  18  im = 1, m
   18   ci(im,in) = -ci(im,in)
      endif

      end
