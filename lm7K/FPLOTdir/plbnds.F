      program plbnds
C- Band plotter using bnds file in standard format for lmto programs.
C ----------------------------------------------------------------------
Cl Local variables
Cl   ifi bands file
Cl   ifm file logical unit for merge file
Cl   ifw file logical unit for new merge file to be created
Cu Updates
Cu   08 Jul 08 Read bnds file with q-dependent number of bands
Cu   07 Jun 08 printout 'mc' merge command for 3 color weights
Cu   01 May 07 plbnds used variables ymin,ymax when making fplot cmd
Cu   12 Jan 07 New switch -nocol to suppress using color weight info
Cu   06 Jun 06 can read files with color weights
C ----------------------------------------------------------------------
      implicit none
      integer nqmx,maxlin,louts,nouts
      parameter(nqmx=4500, maxlin=50, louts=140, nouts=80 )
      double precision xline(0:maxlin),xdum(2),ydum(2)
      double precision xq(0:nqmx),eq(nqmx)
      real(8),allocatable:: eb(:,:),ebw(:,:,:),eb2(:),ebw2(:,:)
      double precision dq,emin,ebot,ebotb,etopb,emax,htcm,
     .etik,xmin,xmax,ql1,ql2,ql3,e,h,x,q1,q2,q3,eferm,efermm,widcm,
     .etop,xq0,efnew,scl,syma(4),clip(4),fgs(0:3),sgs(0:4),qm1,qm2,
     .qm3,csw(12)
      integer ib,jb,kb,lb,nb,nbf,nbfil,nbm,iq,nq,nlab,iq1,nq1,i,ifi,ifm,
     .ifw,iline,nline,i1mach,isw,nql(maxlin),ib1,ib2,j,ltfrme(3),psym,
     .nargf,fxst,iouts,iiouts,ipr,icolwt,ncolwt
      integer NULLI
      parameter (NULLI=-99999)
C Heap allocation
      integer wksize
      parameter(wksize=2000000)
      integer w(wksize)
      common /w/ w

      logical cmdopt,a2bin,ls1,ls2,lpad,rdstrn,lcol3
      character title*120,xlabel*120,ylabel*120,font*12,ltype*120,
     .ll*120,s*80,outs(nouts)*(louts),prmed*40,ss*80,smerge*120,
     .xnfmt*20,ynfmt*20,datnam*10,scol3*512,fncol3*80
      character*10 lab(maxlin)
      integer fontsz
      integer fopna,fopnx,lgcmd,fopng,lenf,awrite
      integer iarg

C --- Setup and defaults ---
      xnfmt = '%;4d'
      ynfmt = '%;4d'
      if (cmdopt('-lt=',4,0,ltype)) then
        ltype(1:4) = '-lt '
      else
        ltype = ' '
      endif
      if (cmdopt('-dat=',4,0,s)) then
        datnam = s(6:)
      else
        datnam = 'dat'
      endif
      call dpzero(csw,8)

      ls1 = cmdopt('-spin1',6,0,s)
      ls2 = cmdopt('-spin2',6,0,s)
C     Prints out a hypothetical merge instruction if an independent
C     set of bands were made with a third set of color weights.      
      lcol3 = cmdopt('-col3',5,0,s) .or. cmdopt('--col3',6,0,s)
      scol3 = ' '
      if (lcol3) scol3 = 'mc ''-ff8.5,9f9.4:/(8x,9f9.4)'''

      if (ls1 .and. ls2) call rx('plbnds only spin1 or spin2 allowed')

      sgs(0) = 2
      sgs(1) = 0
      sgs(2) = 0
      sgs(3) = 0
      sgs(4) = -1
      fgs(0) = 3
      fgs(1) = 1
      fgs(2) = 1
      fgs(3) = 1
      ncolwt = 0
      iouts = 1
      call pshpr(0)
      call wkinit(wksize)
      call poppr
      call finits(2,0,0,i)
      do  19  i = 1, nouts
   19 outs(i) = ' '
C ... Check that all arguments before last are switches
      if (nargf() .le. 1) goto 999
      ifm = 0
      ifw = 0
      lpad = .false.
      smerge = ' '
      do  iarg = 1, nargf()-2
        call getarf(iarg,s)
        if (s(1:1) .ne. '-') goto 999
        if (s(1:7) .eq. '-merge=') smerge = s(8:)
        if (s(1:8) .eq. '-mergep=') then
          smerge = s(9:)
          lpad = .true.
        endif
      enddo
C  18 continue
C ... Last argument is file name : open file
C     First check whether name exists
      call getarf(nargf()-1,s)
      j = fopnx(s,172,-1,-1)
      if (j .ne. 0) then
        ifi = fopng(s,-1,1)
C     Else check whether extension exists
      else
        j = fxst('bnds')
        if (j .eq. 0) then
          call strip(s,i,j)
          print 345, s(i:j), s(i:j)
  345     format(' plbnds: failed to open file "',
     .    a,'" or file "bnds.',a,'".')
          goto 999
        endif
        ifi = fopna('bnds',-1,0)
      endif

C ... Open merge file, if specified, file unit ifm
C     First check whether name exists
      if (smerge .ne. ' ') then
        call wordg(smerge,0,', ',1,i,j)
        if (smerge(i:j) .ne. smerge) then
          ss = smerge(i:j)
          call wordg(smerge,0,', ',2,i,j)          
          ifw = fopng(smerge(i:j),-1,2)
          smerge = ss
        endif
        if (smerge .eq. s)
     .  call rx('merge and original file of same name')
        s = smerge
        j = fopnx(s,172,-1,-1)
        if (j .ne. 0) then
          ifm = fopng(s,-1,1)
C       Else check whether extension exists
        else
          s = 'bnds.'//smerge        
          j = fopnx(s,172,-1,-1)
          if (j .eq. 0) then
            call strip(smerge,i,j)
            print 345, smerge(i:j), smerge(i:j)
            goto 999
          endif
          ifm = fopng(s,-1,1)
        endif
      endif

C      j = fxst('bnds')
C      if (j .eq. 0) goto 999
C      ifi = fopna('bnds',-1,0)
      htcm = 16d0
C     gtoxcm = 6d0
      title = ' '
      xlabel = ' '
      prmed = ' ps'
      psym = 6
      scl = 1
      syma(1) = .2d0
      syma(2) = 8
      if (cmdopt('-ascii',6,0,prmed)) then
        print *, 'writing ascii'
        psym = 1
      endif
      call getpr(ipr)
      ylabel = 'Energy (Ry)'
      font = '/Times-Roman'
      fontsz = 24
      do  1  i = 1, maxlin
    1 lab(i) = ' '
      ltfrme(1) = 3
      ltfrme(2) = 0
      ltfrme(3) = 0

C --- Band file header, setup ---
C     For backwards compatibility, read 3rd argument as real; if so, icolwt=0
C     smerge used as a temporary character variable
      rewind ifi
      if (.not. rdstrn(ifi,smerge,120,.false.)) goto 999
      call words(smerge,i)
      rewind ifi
      if (i .lt. 3) then
        read(ifi,*,end=999,err=999) nbf,eferm
        x = 0
      else
        read(ifi,*,end=999,err=999) nbf,eferm, x
C       For color weights, icolwt must be a positive integer
        icolwt = 0
        if (abs(nint(x)-x) .eq. 0) icolwt = nint(x)
        if (icolwt .lt. 0) icolwt = 0
        call word(smerge,3,i,j)
        i = i-1
        call chrpos(smerge,'.',j,i)
        if (i .lt. j) icolwt = 0
      endif

      if (icolwt .eq. 2) then
        call info0(0,0,0,
     .  ' plbnds : bands file contains two sets of color weights')
      elseif (icolwt .ne. 0) then
        call info0(0,0,0,' plbnds : bands file contains color weights')
      endif
      if (icolwt .gt. 2) then
        call info0(0,0,0,' plbnds (warning) reset # color weights to 1')
        icolwt = 1
      endif
      if (icolwt .lt. 2) lcol3 = .false.
      if (lcol3) then
        call info0(0,0,0,' plbnds : merge third color weight')        
      endif
      ncolwt = icolwt
      if (cmdopt('-nocol',6,0,s) .or. cmdopt('--nocol',7,0,s)) then
        ncolwt = 0
        call info0(0,0,0,' plbnds : color weights suppressed')
      endif

C#ifdef F90
      allocate(eb(nbf,nqmx),eb2(nbf))
      if (icolwt .ne. 0) then
        allocate(ebw(nbf,nqmx,icolwt),ebw2(nbf,icolwt))
      endif
C#endif
      nb = nbf
      if (ifm .ne. 0) then
        read(ifm,*,end=999,err=999) nbm,efermm
        if (nbm .ne. nbf)
     .  call info2(0,0,0,' plbnds (warning) : original file has %i '//
     .  'bands while merged file has %i bands',nbf,nbm)
        nb = min(nbm,nbf)
        if (lpad) nb = max(nbm,nbf)
        if (eferm .ne. efermm)
     .  call info2(0,0,0,' plbnds (warning) : original file has ef='//
     .  '%d while merged file has ef=%d',eferm,efermm)
        if (ls1 .or. ls2)
     .  call rx('-merge not compatible with -spin1 or -spin2')
      endif
      if (ifw .ne. 0) then
        write(ifw,'(i5,f10.5,i5)') nb,eferm,0
      endif

      nq = 0
      nline = 0
      xq(0) = 0d0
      ebot = 1d10
      etop = -1d10
      xq0 = 0
      xline(0) = 0
      efnew = eferm
      if (cmdopt('-ef=',4,0,s)) then
        j = 4
        call rxx(.not.a2bin(s,efnew,4,0,' ',j,len(s)),
     .  'failed to parse '//s)
      endif
      if (cmdopt('-scl=',5,0,s)) then
        j = 5
        call rxx(.not.a2bin(s,scl,4,0,' ',j,len(s)),
     .  'failed to parse '//s)
      endif
      call awrit2(' plbnds: %i bands  fermi=%;6d',
     .s,80,0,nb,eferm)
      if (scl .ne. 1)
     .call awrit1('%a  scaled by %;6d',s,80,0,scl)
      if (efnew .ne. eferm)
     .call awrit1('%a  shifted to %;6d',s,80,0,efnew)
      call awrit0('%a',s,80,-i1mach(2))
C#ifndefC F90
C      call rxx(nb.gt.nbmx,'nb gt nbmx')
C#endif

C --- For each panel, do ---
   91 continue
      ebotb = 1d10
      etopb = -1d10
      read(ifi,*) nq1
C     Check for line match with merge file, if one is to be read
      if (ifm .ne. 0) then
        read(ifm,*) j
        if (nq1 .ne. j) 
     .  call info5(0,0,0,' plbnds (fatal) : mismatch panel %i: file'//
     .  ' has %i qp while merged file has %i qp',nline+1,nq1,j,0,0)
        if (nq1 .ne. j) goto 995        
        nq1 = nq1*2
      endif
C     If new merge file is to be created, write number of points on line
      if (ifw .ne. 0) then
        write(ifw,'(i6)') nq1
      endif
      if (nq1 .le. 0) goto 90
      if ((ls1 .or. ls2) .and. mod(nq1,2) .ne. 0) then
        call setpr(30)
        call rx('plbnds: spin pol case needs even number of bands')
      endif
      if (ls1 .or. ls2 .or. ifm .ne. 0) nq1 = nq1/2
C     Read energy bands along this line from bands file
      do  3  iq1 = 1, nq1
        nq = nq+1
        if (nq .gt. nqmx) stop 'increase nqmx in plbnds'
        nbfil = 0
        call rdbnd(ifi,q1,q2,q3,nbf,nbfil,-NULLI,eb(1,nq))
C       read(ifi,500) q1,q2,q3,(eb(i,nq),i=1,nbf)
        if (icolwt .ne. 0) then
C        Read associated file weights
          do  j = 1, icolwt
            call rdbnd(ifi,qm1,qm2,qm3,nbf,nbfil,0,ebw(1,nq,j))
C           read(ifi,500) qm1,qm2,qm3,(ebw(i,nq,j),i=1,nbf)
            if (qm1 .ne. q1) goto 995
            if (qm2 .ne. q2) goto 995
            if (qm3 .ne. q3) goto 995
          enddo
        endif
C       Read corresponding energy bands from merge file
        if (ifm .ne. 0) then
          read(ifm,*) qm1,qm2,qm3,(eb2(i),i=1,nbm)
          if (qm1 .ne. q1) goto 995
          if (qm2 .ne. q2) goto 995
          if (qm3 .ne. q3) goto 995
C       If only one spin is to be extracted, read the second spin
        else
          if (ls1 .or. ls2) then
            nbfil = 0
            call rdbnd(ifi,q1,q2,q3,nb,nbfil,-NULLI,eb2)
C           read(ifi,500) q1,q2,q3,(eb2(i),i=1,nb)
            if (icolwt .ne. 0) then
              do  j = 1, icolwt
                call rdbnd(ifi,q1,q2,q3,nb,nbfil,0,ebw2(1,j))
C               read(ifi,500) q1,q2,q3,(ebw2(i,j),i=1,nb)
              enddo
            endif
          endif
        endif
        if (lpad) then
          do  i= nbf+1,nb
            eb(i,nq) = -NULLI
          enddo
          do  i= nbm+1,nb
            eb2(i) = -NULLI
          enddo
        endif
C       If new merge file is to be created, write the bands
        if (ifw .ne. 0) then
          write(ifw,500) q1,q2,q3,(eb(i,nq),i=1,nb)
          write(ifw,500) q1,q2,q3,(eb2(i),i=1,nb)
        endif
C       If 2nd spin is to be extracted, poke bands into spin1
        if (ls2) then
          call dcopy(nb,eb2,1,eb(1,nq),1)
          if (icolwt .ne. 0) then
            do  j = 1, icolwt
              call dcopy(nb,ebw2(1,j),1,ebw(1,nq,j),1)
            enddo
          endif
        endif
  500   format(3f10.5/(10f8.4))
  501   format((10f8.4))
C       Shift and scale bands
        do  9  i = 1, nb
    9   eb(i,nq) = eb(i,nq)*scl + efnew - eferm*scl
C       Find upper and lower bounds (global and for this panel)
        do  2  i = 1, nb
          ebot = dmin1(ebot,eb(i,nq))
          etop = dmax1(etop,eb(i,nq))
          ebotb = dmin1(ebotb,eb(i,nq))
          etopb = dmax1(etopb,eb(i,nq))
    2   continue
        dq = 0
        if (iq1 .gt. 1) dq = dsqrt((q1-ql1)**2+(q2-ql2)**2+(q3-ql3)**2)
        xq(nq) = xq(nq-1)+dq
        if (ifm .ne. 0) then
          xq(nq+1) = xq(nq)
          nq = nq+1
        endif
        if (ipr .gt. 40)
     .  call awrit7(' panel %i iq=%i  dq=%;6d (%;6d %;6d %;6d) xq=%;6d',
     .  ' ',80,i1mach(2),nline,nq,dq,(q1-ql1),(q2-ql2),(q3-ql3),xq(nq))
        ql1 = q1
        ql2 = q2
        ql3 = q3
    3 continue
      if (ifm .ne. 0) nq1 = nq1*2
      nline = nline+1
      nql(nline) = nq1
      xline(nline) = xq(nq)
      call awrit5(' panel %i  nq=%i  ebot=%1;6d  etop=%1;6d'//
     .'  delta q=%1;6d',' ',80,i1mach(2),nline,nq1,ebotb,etopb,
     .xq(nq)-xq0)
      xq0 = xq(nq)
      goto 91

C --- Get dimensions for plot ---
   90 nline = nline-1
      if (ifw .ne. 0) call fclose(ifw)
      xmin = xq(1)
      xmax = xq(nq)
      call awrit6(' nq=%i  npan=%i  emin=%1;6d  '//
     .'ef=%1;6d  emax=%1;6d  sum dq=%1;6d',
     .' ',80,i1mach(2),nq,nline,ebot,efnew,etop,xmax-xmin)
      emin = ebot
      emax = etop
C     gtoxcm = 5
      htcm = 10
   93 print '('' emin, emax, width(cm), height(cm) ?'')'
      read(5,*)  emin,emax,widcm,htcm
      if (widcm+2d0 .gt. 25d0) then
        call awrit1(' Total width (%1;1dcm) cannot exceed 23cm',
     .  ' ',80,i1mach(2),widcm)
        goto 93
      endif
      call dcopy(4,0d0,0,clip,1)
      clip(2) = widcm/14d0
      clip(4) = htcm/14d0

      nlab = 0
      if (prmed .ne. '-ascii' .and.
     .(cmdopt('-lbl',4,0,s) .or. cmdopt('-tl',3,0,s))) then
        if (cmdopt('-tl=',4,0,s)) then
          title = s(5:)
        elseif (.not. cmdopt('-lbl=',5,0,s)) then
          call cwrite(' Title ? ',0,28,0)
          read(*,'(a80)') title
        endif
        s = ' '
        if (cmdopt('-lbl=',5,0,s)) then
          call wrdsg(s(6:len(s)),0,',',nlab)
          do  46  j = 1, nlab
            call wordg(s(6:len(s)),0,',',j,ib1,ib2)
   46     lab(j) = s(5+ib1:5+ib2)
        elseif (cmdopt('-lbl',4,0,s)) then
          call cwrite(' Symmetry labels ? ',0,28,0)
          read(*,'(a80)') s
          call words(s,nlab)
          do  45  j = 1, nlab
            call word(s,j,ib1,ib2)
   45     lab(j) = s(ib1:ib2)
        endif
      endif

C --- Make data files in format suitable for fplot ---
      if (cmdopt('-fplot',6,0,s)) then
        lgcmd = fopng('plot.plbnds',-1,0)
C   ... Assign line type to character variable
        if (ltype .ne. ' ') then
          write(lgcmd,'(''%char0 ltb="'',a,''"'')') trim(ltype(5:))
          ltype(5:) = '{ltb}'
        endif
        outs(1) = 'fplot '
        nq1 = 0
        do  40  iline = 1, nline+1
          call awrit5('%a'//
     .    '%?#n==1#%N%% var ymin=%;3d ymax=%;3d%N#%2j#'//
     .    ' -frme %;3d,%;3d,0,1 -p0 -x 0,1 -y {ymin},{ymax}'//
     .    ' -noxn -1p',outs(iouts),louts,0,iline,emin,emax,
     .    xline(iline-1)/(xmax-xmin),xline(iline)/(xmax-xmin))
          if (title .ne. ' ') then
            call strip(outs(1),ib1,ib2)
            outs(1)(ib2+1:) = ' -tl "' // title
            call strip(outs(1),ib1,ib2)
            outs(1)(ib2+1:ib2+1) = '"'
            title = ' '
          endif
C         make ib1,ib2 = indices to lowest, highest band for this panel
          ib1 = nb+1
          ib2 = -1
          do  30  ib = 1, nb
            ebotb = 1d10
            etopb = -1d10
            do  20  iq = 1+nq1, nql(iline)+nq1
              ebotb = dmin1(ebotb,eb(ib,iq))
              etopb = dmax1(etopb,eb(ib,iq))
   20       continue
            if (ebotb .lt. emax .and. etopb .gt. emin) then
              ib1 = min(ib1,ib)
              ib2 = max(ib2,ib)
            endif
   30     continue

C         Open band file for this panel
          call awrit3(' write file bnd%i.'//datnam//'%a, bands %i - %i',
     .    ' ',80,i1mach(2),iline,ib1,ib2)
          lenf = awrite('bnd%i.'//datnam//'%a',s,len(s),0,iline,0,0,0,0,
     .    0,0,0)
          open(10, file=s(1:lenf), status='UNKNOWN', err=99)

          if (lcol3) then
            if (cmdopt('-col3:',6,0,ss)) then
              fncol3 = trim(ss(7:)) // '/' // s
            elseif (cmdopt('--col3:',7,0,ss)) then
              fncol3 = trim(ss(8:)) // '/' // s
            else
              fncol3 = 'tmp/' // s
            endif
          endif

C         Two distinct spin channels: write them both out
          if (cmdopt('-fplot:s',8,0,ss)) then
            if (mod(nql(iline),2) .ne. 0)
     .      call rx('spin pol bands should have even nq')
            j = 2*(ib2-ib1+1)
            if (ncolwt .ne. 0) j = j*(ncolwt+1)
            call awrit2('%% rows %i cols %i',' ',80,10,nql(iline)/2,1+j)
            do  36  iq = 1+nq1, nql(iline)+nq1,2
              write(10,333)
     .        (xq(iq)-xq(1+nq1))/(xq(nql(iline)+nq1)-xq(1+nq1)),
     .        (eb(ib,iq),eb(ib,iq+1),ib=ib1,ib2)
              if (ncolwt .ne. 0) then
                do  j = 1, ncolwt
                  write(10,334) (ebw(ib,iq,j),ebw(ib,iq+1,j),ib=ib1,ib2)
                enddo
              endif
   36       continue

C         Write out bands and optional weights
          else
            j = ib2-ib1+1
            if (ncolwt .ne. 0) j = j*(ncolwt+1)
            call awrit2('%% rows %i cols %i',' ',80,10,nql(iline),1+j)
            do  35  iq = 1+nq1, nql(iline)+nq1
              write(10,333)
     .        (xq(iq)-xq(1+nq1))/(xq(nql(iline)+nq1)-xq(1+nq1)),
     .        (eb(ib,iq),ib=ib1,ib2)
              if (ncolwt .ne. 0) then
                do  j = 1, ncolwt
                  write(10,334) (ebw(ib,iq,j),ib=ib1,ib2)
                enddo
              endif
   35       continue
          endif
C         Build up mc command for color merge
          if (lcol3 .and. ncolwt .ne. 0) then
            j = ib2-ib1+1
C           if (ncolwt .ne. 0) j = j*(ncolwt+1)
c            print *, '-qr', s(1:lenf), ' -qr ' 'tmp/'//s(1:lenf)
            call awrit2('%a -qr '//s(1:lenf)//' -qr '//trim(fncol3)//
     .      ' -sub 1,nr,%i,%i -ccat -w '//s(1:lenf),scol3,len(scol3),
     .      0,2+j,1+2*j)
          endif

  333     format(f8.5,9f9.4:/(8x,9f9.4))
  334     format((8x,9f9.4))
          close(10)
          nq1 = nq1 + nql(iline)
          if (iline .gt. 1)call awrit1('%a -noyn',outs(iouts),louts,0,1)
          iouts = iiouts(outs,iouts,louts+64)
          if (iouts .gt. nouts) call rx('plbnd output string too long')
          ib = 1+(ib2-ib1+1)
          if (cmdopt('-fplot:s',lenf,0,ss)) ib = 1+2*(ib2-ib1+1)
          jb = 1+2*(ib2-ib1+1)
          if (cmdopt('-fplot:s',lenf,0,ss)) jb = 1+4*(ib2-ib1+1)
          kb = 1+3*(ib2-ib1+1)
          lb = 1+4*(ib2-ib1+1)
          if (cmdopt('-fplot:s',lenf,0,ss)) kb = 1+6*(ib2-ib1+1)

          call awrit2('%a -lt 2,bold=3,col=0,0,0,1,.5,.3,.5 '//
     .    '-tp 2~0,%1;4d,1,%1;4d -lt 1',
     .    outs(iouts),louts,0,efnew,efnew)
          iouts = iiouts(outs,iouts,louts+90)
          if (iouts .gt. nouts) call rx('plbnd output string too long')
          if (ncolwt .eq. 2 .and. lcol3) then
            call awrit8('%a -colsy 2:%i'//
     .      ' -colsw %i:%i -colsw2 %i:%i -colsw3 %i:%i '//ltype//
     .      '%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,
     .      0,ib,ib+1,jb,jb+1,kb,kb+1,lb,iline)
          elseif (ncolwt .eq. 2) then
            call awrit6('%a -colsy 2:%i'//
     .      ' -colsw %i:%i -colsw2 %i:%i '//ltype//
     .      '%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,
     .      0,ib,ib+1,jb,jb+1,kb,iline)
          elseif (ncolwt .ne. 0) then
            call awrit4('%a -colsy 2:%i'//
     .      ' -colsw %i:%i '//ltype//'%a -qr bnd%i.'//datnam//'%a',
     .      outs(iouts),louts,0,ib,ib+1,jb,iline)
          else
            call awrit2('%a -colsy 2:%i '//
     .      ltype//'%a -qr bnd%i.'//datnam//'%a',outs(iouts),louts,0,
     .      ib,iline)
          endif

C     ... Symmetry label, if there is one
          if (iline .le. nlab) then
            iouts = iiouts(outs,iouts,louts+30)
            call awrit3('%a -font %?#n#s#t#24'//
     .      ' -lblu %;2d,{ymin-(ymax-ymin)/60} cd '//
     .      lab(iline),outs(iouts),louts,0,isw(lab(iline).eq.'G'),
     .      0d0,emin-(emax-emin)/60)
          endif
          iouts = iiouts(outs,iouts,louts+10)
          if (iouts .gt. nouts) call rx('plbnd output string too long')
   40   continue
C   ... Last symmetry label, if there is one
        if (nline+2 .le. nlab) then
          iouts = iiouts(outs,iouts,louts+30)
          call awrit3('%a -font %?#n#s#t#24'//
     .    ' -lblu %;2d,{ymin-(ymax-ymin)/60} cd '//
     .    lab(nline+2),outs(iouts),louts,0,isw(lab(nline+2).eq.'G'),
     .    1d0,emin-(emax-emin)/60)
        endif
        do  25  ib = 1, iouts
          call skpblb(outs(ib),louts,i)
C         call cwrite(outs(ib),0,i,1)
          write(lgcmd,'(a)') outs(ib)(1:i+1)
   25   continue
        if (lcol3 .and. ncolwt .ne. 0) then
          print *, ' ... to merge another color from file tmp:'
          print '(a)', trim(scol3)
        endif
        print *, ' ... to plot, invoke:'
        print *, ' fplot -disp -f plot.plbnds'
        call fclose(lgcmd)
        call cexit(0,1)
   99   call rx(' file open error, file='//s(1:lenf))
      endif

C --- Plot frame ---
      call pltini(prmed(2:10),0,.false.)
      call pltstp(clip(1),clip(2),clip(3),clip(4))
      call pltstu(xmin,xmax,emin,emax,0d0)
      etik = 0.1d0
      if (dabs(emax-emin) .gt. 5d0) etik = 0.5d0
      if (dabs(emax-emin) .gt. 10d0) etik = 1d0
      if (dabs(emax-emin) .gt. 20d0) etik = 2d0
      call setfnt(font,fontsz)
      call frme(-1d30,xmax-xmin,1,0d0,   .6d0,0,0,xnfmt,xlabel,
     .-1d30,etik,     2,.025d0,.6d0,0,1,ynfmt,ylabel,
     .' ',xdum,xdum,xdum,title,ltfrme, fgs)
      do  4  iline = 1, nline
        call mve(xline(iline),emax)
        call drw(xline(iline),emin)
    4 continue

C --- Put in vertical lines and symmetry labels ---
      do  5  i = 0, nline+1
        ll = ' '
        ll = lab(i+1)
C       print *, ll
        call setfnt(font,fontsz)
        if (ll .eq. 'G') call setfnt('/Symbol     ',fontsz)
        if (i .eq. 0) then
          x = xmin
        elseif (i .eq. nline+1) then
          x = xmax
        else
          x = xline(i)
        endif
        h = fontsz / (72 * 2.54 * 2) * scl
        call pstr(x,emin-h,0d0,ll,0d0,0d0,0d0,'cc','u',1)
C       call mve(x,emin)
C       call pslabl('a',0,0d0,ll,'h','f',1,0d0)
    5 continue

C --- Plot points ---
      do  7  ib = 1, nb
        if (mod(ib,10) .eq. 1)
     .  call awrit1(' Plotting band %i',' ',80,i1mach(2),ib)
        do  8  iq = 1, nq
    8   eq(iq) = eb(ib,iq)
        if (cmdopt('-fast',5,0,s)) then
          call plsym(xq(1),eq,0d0,nq,psym,syma,sgs)
        else
          do  6  iq = 1, nq
            e = eb(ib,iq)
            if (e .ge. emin .and. e .le. emax) then
              call plsym(xq(iq),e,0d0,1,psym,syma,sgs)
            endif
    6     continue
        endif
    7 continue

C --- Draw in Fermi energy ---
      if (efnew .ge. emin .and. efnew .le. emax) then
        xdum(1) = xmin
        xdum(2) = xmax
        ydum(1) = efnew
        ydum(2) = efnew
        call plntyp(2,3,2d0,1d0,0d0,0d0)
        call plcrv(xdum,ydum,xdum,xdum,xdum,2,xdum(1),xdum(2),
     .  ydum(1)-0.1d0,ydum(2)+0.1d0,0,-2d0,csw,11)
      endif
      call pltdmp(0)
      call cexit(0,1)

C --- Error exit ---
  999 continue
      print '(a,a)',
     .' usage: plbnds [-switches] filename'
      print 998, 'Last arg can be true file name or "ext" with'//
     .' filename = bnds.ext'
      print 998, 'Switches:'
      print 998, '-ef=#     shifts fermi energy to #'
      print 998, '-scl=#    scales bands by #'
      print 998, '-spin1    plots only spin-1 bands'
      print 998, '-spin2    plots only spin-1 bands'
      print 998, '-fast     '
      print 998, '-lbl      plbnds prompts for symmetry labels'//
     .' and title'
      print 998, '-lbl=CCC.. supply labels through this switch'
      print 998, '-tl=string supply title through this switch'
      print 998, '-fplot[:s] generates an fplot command and files'//
     .' to use with fplot'
      print 998, '          use optional :s if file contains two spins'
      print 998, '-lt=...   (with fplot) uses specified fplot line type'
      print 998, '-col3[:dir] (with fplot) flag to merge separate '//
     .' bnds file with extra color'
      print 998, '-dat=name (with fplot) writes data files bnd*.name'
      print 998, '-merge=file2[,file3]'
      print 998, '-mergep=file2[,file3]'
      print 998, '          reads two input files, one for each spin'//
     .' (spin-pol case)'
      print 998, '          Optional file3 causes plbnds to write'//
     .' merged file to file3'
      print 998, '          -mergep pads channel with fewer bands'
  998 format(15x,a:a)
      call cexit(-1,1)

  995 continue
      call pshpr(40)
      call rx('mismatch between original and merged files')

      end
      integer function iiouts(outs,iln,itol)
      implicit none
      integer iln,itol,louts,i
      character*(*) outs(1)
      louts = len(outs(1))
      call skpblb(outs(iln),louts,i)
      if (i .gt. louts-itol) then
C       outs(iln)(i+2:i+2) = '\\'
        iln = iln+1
      endif
      iiouts = iln
      end
      subroutine rdbnd(ifi,q1,q2,q3,nbf,nbfil,ipad,eb)
C- Read bands for one qp from standard bands file
C  and pads bands from nbfil .. nbf
Ci nbfil number of bands to read
Ci       input nbfil 0 => attempt to read from file.
Ci       If missing, attempt to read nbf bands.
Co       On output, nbfil = number actually read
      implicit none
      integer ifi,nbf,ipad,nbfil
      double precision q1,q2,q3,eb(nbf)
      character strn*120
      integer nw,j,j1,j2,ix(4),a2vec,i
      double precision fac(4)
      logical rdstrn

C     Check whether format contains number of bands to read
      if (.not. rdstrn(ifi,strn,len(strn),.false.)) goto 99
      call words(strn,nw)
      if (nw .lt. 3) goto 99
      if (nw .eq. 3) then
        if (nbfil .eq. 0) nbfil = nbf
        j = 0
        j = a2vec(strn,len(strn),j,4,', ',2,-3,3,ix,fac)
        if (j .ne. 3) goto 99
      elseif (nw .eq. 4) then
        call word(strn,4,j1,j2)
        j = 0
        j = a2vec(strn,len(strn),j,4,', ',2,-3,4,ix,fac)
        if (j .ne. 4) goto 99
        if (nbfil .eq. 0) then
          nbfil = nint(fac(4))
          if (nbfil .ne. nint(fac(4))) then
            call rx('plbnds:  expected integer for 4th '//
     .      'argument, line: '//trim(strn))
          endif
        endif
      endif
      q1 = fac(1)
      q2 = fac(2)
      q3 = fac(3)

      read(ifi,501) (eb(i),i=1,nbfil)
  501 format((10f8.4))
      do  i = nbfil+1, nbf
        eb(i) = ipad
      enddo
      return
   99 continue
      call rx('plbnds:  failed to read qp from file, line:  '//
     .trim(strn))
      
      end
