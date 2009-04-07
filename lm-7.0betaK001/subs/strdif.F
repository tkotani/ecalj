      logical function strdif(fn1,fn2)
C- Compares two structure constants files
C ----------------------------------------------------------------------
Ci Inputs
Ci   fn1   :file name for first strux
Ci   fn2   :file name for second strux
Co Outputs
Co   strdif:T, if files compare
Cl Local variables
Cl  n1,n2  :n?(1) nttab
Cl         :n?(2) nl
Cl         :n?(3) nbas
Cl         :n?(4) nkap
Cl         :n?(5) itral
Cl         :n?(6) lmaxw
Cl         :n?(7) nitab
Cl         :n?(8) lio
Cr Remarks
Cr
Cu Updates
Cu   07 Aug 06 Updated; works for 2-kappa case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character fn1*(*),fn2*(*)
C ... Local parameters
      logical skips
      integer n1(8),n2(8),ifi1,ifi2,i,ns,nttab,nitab,nl,
     .  nbas,nkap,niax,oalp1,oiax1,onpr1,os1,oalp2,oiax2,onpr2,os2,ik,
     .  ikn,it,ii,lio,itral,nds,nmto,lmaxw,getdig,nalp,ncplx,nkaps,nkapn
      parameter (niax=10)
      logical ldif,dcmp,iostr1,iostr
      integer w(1)
      common /w/ w
      double precision ckbas1,ckbas2,e1(2),e2(2),dval,errmx1,errmx2
      character out1*80, out2*80, cargs(8)*5
      data cargs /'nsite','nl','nbas','nkap','itral','lmaxw','nitab',
     .            'lio'/

      strdif = .false.

C --- Compare headers ---
C     read(ifi) lio,itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
      n1(8)=0
      ldif = iostr1(n1(8),fn1,n1(5),n1(1),n1(7),n1(2),n1(3),n1(4),n1(6),
     .              ckbas1,ifi1)
      if (.not. ldif) then
        out1 = ' strdif: could not read header, file '//fn1
        call info0(10,0,0,out1)
        return
      endif
      n2(8)=0
      ldif = iostr1(n2(8),fn2,n2(5),n2(1),n2(7),n2(2),n2(3),n2(4),n2(6),
     .              ckbas2,ifi2)
      if (.not. ldif) then
        out1 = ' strdif: could not read header, file '//fn2
        call info0(10,0,0,out1)
        return
      endif

C     skips is set to .true. if strux are too different to compare
      out1 = ' '
      out2 = ' '
      skips = .false.
      do  10  i = 1, 8
        if (n1(i) .ne. n2(i)) then
          call awrit0('%a  '//cargs(i),out1,80,0)
          call awrit0('%a  '//cargs(i),out2,80,0)
          call awrit1('%a1=%i',out1,80,0,n1(i))
          call awrit1('%a2=%i',out2,80,0,n2(i))
          if (i .ne. 8) skips = .true.
        endif
   10 continue

      if (.not. dcmp(ckbas1,ckbas2,1,1d-10,i,errmx1)) then
        call awrit1('%a  ckbas1=%1;5d',out1,80,0,ckbas1)
        call awrit1('%a  ckbas2=%1;5d',out2,80,0,ckbas2)
        skips = .true.
      endif

      errmx1 = 0
      errmx2 = 0

      if (out1 .ne. ' ') then
        call info0(10,0,0,
     .    ' STRDIF: structure constant header files differ')
        call strip(out1,ii,it)
        print *, out1(ii:it)
        call strip(out2,ii,it)
        print *, out2(ii:it)
      else
        call info5(10,0,0,'%9fheaders match ...%N'//
     .    '%9flio=%i  nbas=%i  nkap=%i  nttab=%i',
     .    n1(8),n1(3),n1(4),n1(1),0)
      endif

      if (skips) return

C --- Compare contents of strux ---
      out1 = ' '
      out2 = ' '
      lio = n1(8)
      nttab = n1(1)
      nl = n1(2)
      nbas = n1(3)
      nkap = n1(4)
      nitab = n1(7)
      lmaxw = n1(6)
      itral = n1(5)

      ldif = iostr(lio+8,fn1,nl,nbas,nkap,e1,itral,ckbas1,
     .  lmaxw,nitab,oalp1,oiax1,onpr1,os1)
      if (.not. ldif) then
        out1 = 'strdif: could not read strux, file'//fn1
        call info0(10,0,0,out1)
        return
      endif
      ldif = iostr(lio+8,fn2,nl,nbas,nkap,e2,itral,ckbas2,
     .  lmaxw,nitab,oalp2,oiax2,onpr2,os2)
      if (.not. ldif) then
        out1 = 'strdif: could not read strux, file'//fn2
        call info0(10,0,0,out1)
        return
      endif

C     Disentangle nkaps,nkapn
      nmto = getdig(mod(lio/10000,100),0,2)
      nkaps = nkap
      nkapn = 1
      if (nmto .eq. 1) then
        nkaps = 1
        nkapn = nkap
      endif

C     Use nds=nl**2 for now, until nds is passed
      nds = nl*nl
      ncplx = 1+getdig(mod(lio/100,100),0,2)
      nalp = nds*nbas*nkaps*nkaps*nkapn
      ns = nds**2*nitab*nkaps*nkaps*nkapn*ncplx

      if (.not. dcmp(w(oalp1),w(oalp2),nalp,1d-12*0,i,errmx1)) then
        call awrit2('%a  alp1(%i)=%1,6;6g',out1,80,0,i,dval(w(oalp1),i))
        call awrit2('%a  alp2(%i)=%1,6;6g',out2,80,0,i,dval(w(oalp2),i))
      endif
      if (.not. dcmp(w(os1),w(os2),ns,1d-12*0,i,errmx2)) then
C       s=s(nds*nds,nkaps*nkaps,nitab,nkapn)
C       Given i, find s(ii,ik,it) by: it=(i-1)
        ii = i-1
        ikn = ii/(ns/nkapn)
        ii = ii - ikn*(ns/nkapn)
        it = ii/(ns/nitab/nkapn)
        ii = ii - it*(ns/nitab/nkapn)
        ik = ii/(nds**2)
        ii = ii - ik*(nds**2)
        ii = ii+1
        ik = ik+1
        it = it+1
        ikn = ikn+1
        call awrit5('%a  s1(%i,kap=%i,it=%i,kn=%i)=%1;6g',out1,80,0,ii,
     .    it,ik,ikn,dval(w(os1),i))
        call awrit5('%a  s2(%i,kap=%i,it=%i,nk=%i)=%1;6g',out2,80,0,ii,
     .    it,ik,ikn,dval(w(os2),i))
      endif

C --- Printout ---
      if (max(errmx1,errmx2) .ne. 0d0) then
        call info2(10,0,0,' STRDIF: structure constants differ: '//
     .    'errmx(alp)=%1;4g; errmx(s)=%1;4g',errmx1,errmx2)
        call strip(out1,ii,it)
        print *, out1(ii:it)
        call strip(out2,ii,it)
        print *, out2(ii:it)
      else
        call info0(10,0,0,
     .    ' STRDIF: structure constant files are equivalent')
        strdif = .true.
      endif

      end
      logical function dcmp(a1,a2,n,tol,m,errmx)
C- compares two arrays for equality within a tolerance
C     implicit none
      integer n,m,i
      double precision a1(n),a2(n),tol,errmx

      errmx = 0d0
      m = 0
      do  10  i = 1, n
        if (dabs(a1(i)-a2(i)) .le. tol) goto 10
        if (dabs(a1(i)-a2(i)) .gt. errmx) then
          errmx = dabs(a1(i)-a2(i))
          m = i
        endif
   10 continue

      dcmp = m .eq. 0
      end

