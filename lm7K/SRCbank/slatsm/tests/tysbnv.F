      subroutine fmain
      implicit none
      INTEGER   lda, ldw, nn1, n, i, j, i0, j0,
     .  opw, opw2, wksize, nlev, ldap, na, nmax, nd, nn
      double precision mflops
      PARAMETER (lda=901,ldap=12,na=8,nn=24,wksize=6000000,nmax=nn*nn+1)
      integer offs(lda),ija(2,nmax)
      character*1 cs, outs*80
      double precision apr(ldap,ldap,na),api(ldap,ldap,na),
     .  sa( lda, lda, 2 ), pa( lda, lda, 2), pw( lda, lda, 2)
      double precision stime, temp, smflops, qtime, qmflops
      double precision cpusec,xx,errmx
      integer w(wksize)
      common /w/ w
      common /static/ sa, pa, pw

C --- Setup ---
      call finits(2,0,0,i)
      call wkinit(lda*lda*4)
      call nada
      call tc('on')
      nlev = 2
      cs = ' '

C --- Matrix inversion, yygefa,yygedi ---
      CALL init(sa,lda,apr,api,ija,offs,nn)
      nd = offs(nn+1)
      mflops = 2.0d-6 * nd * nd * nd * 4
C     call yprm(.false.,'a',2,6,'(9f8.2)',sa,lda,nd,lda,nd)
      temp = cpusec()
      if (cs .eq. 'h') then
        call yyhifa(sa,sa(1,1,2),lda,nd,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')
        call yyhidi(sa,sa(1,1,2),lda,nd,pw,xx,i,pw(1,2,1),pw(1,3,1),1)
        do  16  i = 1, nd
        do  16  j = 1, i
        sa(i,j,1) =  sa(j,i,1)
   16   sa(i,j,2) = -sa(j,i,2)
      else
        call yygefa(sa,sa(1,1,2),lda,nd,pw,i)
        if (i .ne. 0) call rx(' matrix is singular')
        call yygedi(sa,sa(1,1,2),lda,nd,pw,xx,pw(1,2,1),pw(1,3,1),1)
      endif
      stime = cpusec() - temp
      if (stime .eq. 0d0) stime = 1
      smflops = mflops / stime
*     call yprm(.false.,'a^-1',2,6,'(5f16.10)',sa,lda,nd,lda,nd)

C --- Inversion by yyqinv ---
      call dvset(pa,1,lda*lda*2,-99d0)
      CALL init(pa,lda,apr,api,ija,offs,nn)
C ... (debugging) check usage of w
      ldw = min(nd+5,lda)
      write(*,'('' using lda,ldw,nd,nlev='',4i4)') lda,ldw,nd,nlev
    1 call defcc(opw,  ldw*nd)
      call defrr(opw2, nd*nd)
      call dvset(w(opw),1,ldw*nd*2,-99d0)
      call dvset(w(opw2),1,nd*nd,-99d0)

      temp = cpusec()
C     call yyqinv(cs,pa,pa(1,1,2),lda,nlev,nd,w(opw),ldw,i)
      nn1 = 1
      call yysbnv(11,apr,api,ldap,ija,offs,nlev,nn1,nn,
     .  w(opw),ldw,pa,pa(1,1,2),lda,i)

      qtime = cpusec() - temp
      if (qtime .eq. 0d0) qtime = 1
      qmflops = mflops / qtime

C  ...Check explicitly hermetian
      if (cs .eq. 'h') then
        errmx = 0d0
        i0 = 0
        j0 = 0
        do 10 i = 1, nd
        do 10 j = 1, nd
          if (max(dabs(pa(i,j,1)-pa(j,i,1)),
     .            dabs(pa(i,j,2)+pa(j,i,2))) .gt. errmx) then
            i0 = i
            j0 = j
          endif
   10   errmx = max(errmx,dabs(pa(i,j,1)-pa(j,i,1)),
     .                    dabs(pa(i,j,2)+pa(j,i,2)))
        print 100, 'checking hermetian errmax ',errmx,i0,j0
  100   format(1X,a,':  errmx=',g10.3,' ',2i4)
      endif

C --- Check quality of inversion ---
      CALL diff('compare yygefa,di to yqsnv', sa, pa, lda, nd )
      CALL init (sa,lda,apr,api,ija,offs,nn)
      temp = cpusec()
      call yygemm('N','N',nd,nd,nd,1d0,sa,sa(1,1,2),lda,pa,pa(1,1,2),
     .  lda,0d0,pw,pw(1,1,2),lda)
      temp = cpusec() - temp
      if (temp .eq. 0d0) temp = 1
      do  40  i = 1, nd
   40 pw(i,i,1) = pw(i,i,1) - 1
*     call yprm(.false.,'a^1 * a - 1',2,6,'(16f8.2)',pw,ldw,nd,ldw,nd)
      errmx = 0
      do 42 i = 1, nd
      do 42 j = 1, nd
   42 errmx = max(errmx,dabs(pw(i,j,1)),dabs(pw(i,j,2)))
      print 100, 'inv(a) * a - 1   errmax ',errmx,nint(mflops/temp)

      print 110, stime, smflops, qtime, qmflops, stime/qtime,nd

  110 format(/1X,'Serial time: ',F7.3,'  Serial MFlops: ',F6.1,
     .       /1X,' yqsnv time: ',F7.3,'   yqsnv MFlops: ',F6.1
     .       /1X,'     factor: ',F7.3,'  for nd =',i4)

C ... Make sure pw not overwritten beyond ldw*(nd+1)
      call wkchk('check if pw2 still intact')
      call dcopy(nd*nd,w(opw2),1,pw,1)
      j = 1
      do  50  i = 1, nd*nd
        if (pw(i,j,1) .ne. -99d0) call fexit2(-1,111,
     .    ' STOP: pw overwritten, i,j=%i %i',i,j)
   50 continue
      if (ldw .gt. nd) then
        call dmcpy(w(opw),ldw,1,pw,lda,1,ldw,nd+1)
        do  52  i = 1, ldw
        do  52  j = 1, nd+1
            if (i .le. nd .and. j .le. nd+1) goto 52
            if (pw(i,j,1) .ne. -99d0) call fexit2(-1,111,
     .        ' STOP: pw overwritten, i,j=%i %i',i,j)
   52   continue
      endif
C     call ywrm(0,'pw',1,6,'(16f8.2)',pw,ldw*ldw,ldw,nd+2,nd+2)

      end

      SUBROUTINE init( a, lda, apr, api, ija, offs, nn)
C- Initialize arrays
      implicit none
      integer lda, ldap, na, na1, nn, nnn, offs(1), ija(2,1)
      parameter(ldap=12,na=8,na1=3,nnn=24)
      double precision a(lda,lda,2),apr(ldap,ldap,na),api(ldap,ldap,na)
C Local
      integer ipiax(nnn,nnn),i,j,ip,k
      double precision apl(ldap,ldap,na),pfun(2,ldap,na)
      real ran1


C ... set up apl,pfun
      call ran1in(5)
      do  100  i = 1, ldap*ldap*na
        apr(i,1,1) = int(2000*ran1()-1000)/100d0
  100 continue
      do  110  i = 1, 2*ldap*na
        pfun(i,1,1) = int(2000*ran1()-1000)/100d0
  110 continue
C ... set up ipiax.  Those < na1 go on diagonal
      do  120  i = 1, nn*nn
        ipiax(i,1) = nint(na*ran1()) + na1+1
        if (ipiax(i,1) .gt. na) ipiax(i,1) = 0
  120 continue
      do  122  i = 1, nn
        ipiax(i,i) = int(na1*ran1())+1
  122 continue

C#ifdefC DEBUG
      do  130  i = 1, nn
  130 write(6,'(30i2)') (ipiax(i,j), j=1,nn)
C#endif

C ... set up api
      call dpzero(api,ldap*ldap*na)
      do  8  i = 1, na
        if (i .le. na1) then
          call daxpy(ldap,1d0,pfun(1,1,i),2,apr(1,1,i),ldap+1)
          call dcopy(ldap,pfun(2,1,i),2,api(1,1,i),ldap+1)
        endif
    8 continue

C ... Make offs
      offs(1) = 0
      do  12  i = 1, nn
        offs(i+1) = offs(i) + nint(ldap*ran1())
   12 continue

      call awrit2(' offs %n:1i',' ',120,6,nn+1,offs)

C ... Assemble ija
      ija(1,1) = nn+2
      k = ija(1,1)-1
      do  30  i = 1, nn
      do  32  j = 1, nn
        ip = ipiax(i,j)
        if (ip .ne. 0 .and. i .ne. j) then
          k = k+1
          ija(2,k) = ip
          ija(1,k) = j
        endif
   32 continue
      ija(2,i) = ipiax(i,i)
      ija(1,i+1) = k+1
   30 continue

C ... Unpack a
      call yysp2a(1,1,nn,1,nn,apr,api,ldap,ija,offs,a,a(1,1,2),lda)
      end
      SUBROUTINE diff (strn, sc, pc, ldc, n)
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      integer ldc, n, i, j
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

