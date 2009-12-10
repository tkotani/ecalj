C Tests yqinvb
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
      character cs*4,outs*80
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
      print *, 'cs= (combination of t, h, 4)?'
      read(*,'(a4)') cs
      call query('n=?',2,n)
      nb = n-1
      call query('nb=?',2,nb)
C#endif

       print *, 'test yqinvb: cs=',cs
      if (n .gt. lda .or. nb .gt. lda) call rx('n gt lda')

C --- Matrix inversion, yygefa,yygesl ---
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C      nb = 5
C#endif
      CALL init( sa,lda,n, cs )
      CALL init2( sb,ldb,max(n,nb), cs )
      mflops = 2.0d-6 * n * n * n * 4

C      call yprm(.false.,'a',2,6,'(8f16.10)',sa,lda,n,lda,n)
C      call yprm(.false.,'b',2,6,'(8f16.10)',sb,ldb,n,ldb,nb)
C      call yprm(.false.,'b',2,6,'(8f16.10)',sb,ldb,nb,ldb,n)
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

C --- Inversion-and-backsubstitution by yqinvb ---
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init( sa,lda,n, cs )
      CALL init2(pb,ldb,max(n,nb),cs)
C     call yprm(.false.,'a',2,6,'(9f8.2)',sa,lda,n,lda,n)
C     call yprm(.false.,'b',2,6,'(9f15.10)',pb,ldb,n,ldb*ldb,nb)
C ... (debugging) check usage of w
      ldw = min(n+5,lda)
      write(*,'('' using lda,ldw,n='',4i4)') lda,ldw,n
      call defrr(opw,  ldw*(n+1))
      call defrr(opw2, n*n)
      call dvset(w(opw),1,ldw*(n+1),-99d0)
      call dvset(w(opw2),1,n*n,-99d0)
      call dvset(pw2,1,ldb*(ldb+1),-99d0)

      temp = cpusec()
      call yqinvb(cs,sa,lda*lda,lda,n,nb,w(opw),ldw,pw2,
     .  pb,ldb*ldb,ldb,i)
      qtime = cpusec() - temp
      if (qtime .ne. 0) then
        qmflops = mflops / qtime
      else
        qmflops = 1
      endif
      if (i .ne. 0) call rx('yqinvb failed to find inverse')
*     call yprm(.false.,'a^-1 b',2,6,'(5f16.10)',pb,ldb,n,lda,nb)

C --- Check quality of inversion ---
      if (cs(1:1) .ne. 't') then
        CALL diff('compare yygefa,sl to yqinvb', sb, pb, ldb, n, nb)
      else
        CALL diff('compare yygefa,sl to yqinvb', sb, pb, ldb, nb, n)
      endif
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init( sa,lda,n, cs )
      if (cs(1:1) .ne. 't') then
        call yygemm('N','N',n,nb,n,1d0,sa,sa(1,1,2),lda,pb,pb(1,1,2),
     .  ldb,0d0,sb,sb(1,1,2),ldb)
        CALL init2( pb,ldb,max(n,nb), cs )
        CALL diff('compare a (a^-1 b) to b', sb, pb, ldb, n, nb)
      else
        call yygemm('N','N',nb,n,n,1d0,pb,pb(1,1,2),ldb,
     .    sa,sa(1,1,2),lda,0d0,sb,sb(1,1,2),ldb)
        CALL init2( pb,ldb, max(n,nb), cs )
        CALL diff('compare (b a^-1) a to b', sb, pb, ldb, nb, n)
      endif

      print 110, stime, smflops, qtime, qmflops, qmflops/smflops,n,nb

  110 format(/1X,'Serial time: ',F7.3,'  Serial MFlops: ',F6.1,
     .       /1X,'yqinvb time: ',F7.3,'  yqinvb MFlops: ',F6.1,
     .       /1X,'     factor: ',F7.3,'  for n,nb =',2i4)

C ... Make sure pw not overwritten beyond ldw*(n+1)
      call wkchk('check if pw2 still intact')
      call dcopy(n*n,w(opw2),1,pw,1)
      j = 1
      do  50  i = 1, n*n
        if (pw(i,j,1) .ne. -99d0) call fexit2(-1,111,
     .    ' STOP: pw overwritten, i,j=%i %i',i,j)
   50 continue
      if (ldw .gt. n) then
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
      CALL init ( sa,lda,n, cs )
      CALL init2( pb,ldb,max(n,nb), cs )
      call ran1in(1)
C#ifdefC RDA
C      n = -1
C#endif
      CALL init ( sa,lda,n, cs )
      CALL init2( sb,ldb,max(n,nb), cs )
      call dcopy(ldb*max(n,nb),pb,1,sb,1)
      call yqinvb(cs,sa,lda*lda,lda,n,nb,w(opw),ldw,pw2,
     .  sb,ldb*ldb,ldb,i)

      temp = cpusec()
      call word(cs,1,j1,j2)
      cs(max(j2,0)+1:max(j2,0)+1) = 'b'
      call yqinvb(cs,sa,lda*lda,lda,n,nb,w(opw),ldw,pw2,
     .  pb,ldb*ldb,ldb,i)

      xtim = cpusec() - temp
      xflops = mflops / max(xtim,1d-6)

      print *, ' '
      if (cs(1:1) .ne. 't') then
        CALL diff('compare yqinvb to same w/b', sb, pb, ldb, n, nb)
      else
        CALL diff('compare yqinvb to same w/b', sb, pb, ldb, nb, n)
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
C# speedup factor on the SGI:
C  n  nb  cs=n 
C 150 149 3.896
C 151 150 3.775
C 152 151 3.767
C 153 152 3.804
C 154 153 3.948
C 155 154 3.768
C 156 155 3.820
C 157 156 3.833
C 158 157 4.014
C 159 158 3.860
C 160 159 3.920
C 161 160 3.910
C 162 161 4.035
C 163 162 3.872
C 164 163 3.821
C 165 164 3.924
C 166 165 4.101
C 167 166 3.906
C 168 167 3.963
C 169 168 3.958
C 170 169 4.082
C 171 170 3.912
C 172 171 3.920
C 173 172 3.865
C 174 173 4.186
C 175 174 3.975
C 176 175 4.006
C 177 176 3.982
C 178 177 4.178
C 179 178 3.991
C 180 179 4.000
C 181 180 4.030
C 182 181 4.193
C 183 182 4.027
C 184 183 4.038
C 185 184 4.065
C 186 185 4.232
C 187 186 4.003
C 188 187 4.038

C 150 15  2.316
C 151 15  2.234
C 152 15  2.344
C 153 15  2.268
C 154 15  2.376
C 155 15  2.278
C 156 15  2.375
C 157 15  2.297
C 158 15  2.420
C 159 15  2.318
C 160 16  2.467
C 161 16  2.365
C 162 16  2.460
C 163 16  2.347
C 164 16  2.445
C 165 16  2.361
C 166 16  2.495
C 167 16  2.345
C 168 16  2.520
C 169 16  2.408
C 170 17  2.553
C 171 17  2.426
C 172 17  2.528
C 173 17  2.437
C 174 17  2.579
C 175 17  2.441
C 176 17  2.600
C 177 17  2.505
C 178 17  2.593
C 179 17  2.466
C 180 18  2.607
C 181 18  2.482
C 182 18  2.627
C 183 18  2.494
C 184 18  2.638
C 185 18  2.548
C 186 18  2.670
C 187 18  2.535
C 188 18  2.665
