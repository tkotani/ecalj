      program main
      implicit none
#if MPI|MPIK
      include 'mpif.h'
      integer ierr, nsize, id 
      character(14):: testbuf
      call mpi_init(ierr)
cstart--- check mpi. You can delete this. 
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      write(*,*) 'mpi: size and id=',nsize,id
cend-----
#endif      
      call fmain()
#if MPI|MPIK
      call mpi_finalize(ierr)
#endif
      end program main

      subroutine fmain
C - Band plotter using DOS file in standard format for lmto programs.
Cu Updates
Cu   07 Jun 06 Added -ref options
Cu   04 Nov 02 -lst2 w/out arguments automatically takes members
Cu                   of -lst and increments by 1 (for spin pol)
Cu   13 Jun 02 Added -ef= and -escl= command-line arguments
      implicit none
      integer nemx,ndinmx,ndmx,maxnum
      parameter(ndmx=25, maxnum=300)
      integer fopng,fopnx,i,j,k,j1,j2,id,idtop,ie,ifmt,is,ist,istart,
     .is1,is2,i1mach,awrite,ifi,ifr,md,nd,ne,ne1,ne2,nnum,num1,num2,
     .nrnde,nspin,npdif,lerr,iouts,iiouts,iarg,nargf,lgcmd,ner,mdi,
     .nchan,iscale
      integer num(maxnum),fontsz,ix(maxnum),a2vec,ltfrme(3)
      double precision cmhgt,dcm,dssum,dbot,dtop,dtopmx,
     .dxtik,dy,dytik,dytk,e1,e2,eferm,emax,
     .emin,estart,top,xmax,ymax,ypos1,ypos2,emun,
     .emux,efurm,dul,fgs(0:3),de,eshift,escl,eminr,emaxr,csw(0:11)
      double precision dmx(ndmx+1),xd(2),yd(2),xxv(10),wgt(maxnum)
      character*12 font,dc*1
      integer xnum,ynum
      logical cmdopt,lrat,l2chan,a2bin
      character*100 outs(58)*80,xnfmt*20,ynfmt*20,datnam*40
      character xlab*120,ylab*120,title*120,fn*120,ltdos*80
      character*512 strn,doslst,dosls2
      character*1 strn2(0:1)
      equivalence (strn,strn2)
C#ifdef F90
      real(8),allocatable:: e(:),dtot(:),din(:,:),d(:,:,:),dref(:)
C#elseC
C      parameter(nemx=5001,ndinmx=2000)
C      double precision e(nemx),dtot(nemx),dref(nemx)
C      double precision din(nemx,ndinmx),d(nemx,2,ndmx)
C#endif

C Heap allocation
      integer wksize
      parameter(wksize=100000)
C#ifndefC F90
C      common /static/ d,din,e,dtot
C#endif

      call pshpr(0)
      call wkinit(wksize)
      call poppr

C --- Setup and defaults ---
      call dpzero(csw,8)
      l2chan = .false.
      xnfmt = '%;4d'
      ynfmt = '%;4d'
      fgs(0) = 3
      fgs(1) = 1
      fgs(2) = 1
      fgs(3) = 1
      iouts = 1
      do  5  i = 1, 58
    5 outs(i) = ' '
      call finits(2,0,0,i)
      dtopmx = 9999d0
C ... multiple of dos to round to
      nrnde = 5
C ... trial tick spacings
      dytik = 10d0
      dxtik = .5d0
C ... page size in cm
      xmax = 18d0
      ymax = 25d0
C ... character fonts
      font = '/Times-Roman'
      fontsz = 24
C ... Frame line type
      ltfrme(1) = 3
      ltfrme(2) = 0
      ltfrme(3) = 0
C ... Generates an fplot command and appropriate files
      if (cmdopt('-fplot',6,0,strn)) then
        datnam = 'dat'
        if (cmdopt('-ext',4,0,strn)) then
          dc = strn(5:5)
          if (dc .eq. ' ') then
            call info0(0,0,0,
     .      ' pldos needs file name to accompany '//strn//'%a')
            goto 199
          endif
          datnam = strn(6:)
        endif
        ifi = fopng('dosp.'//datnam,10,0)
      endif

C ---- Input from file ifi -----
C     if (fxst('dos') .lt. 1) goto 199
C     ifi = fopna('dos',-1,0)
      do  7  iarg = 1, nargf()-1
        call getarf(iarg,strn)
        if (strn(1:1) .ne. '-') then
          ifi = fopng(strn,-1,1)
          goto 8
        endif
    7 continue
      goto 199
    8 continue
      print *

C --- Case using reference dos.  Integer flags this option ---
      ifr = 0
      iscale = 0
      if (cmdopt('-ref',4,0,strn)) then
        iscale = 0
        nchan = 0
        dc = strn(5:5)
        if (dc .eq. ' ') then
          call info0(0,0,0,
     .    ' pldos needs file name to accompany '//strn//'%a')
          goto 199
        endif
C   ... Return here to resume parsing for arguments
        j2 = 5-1
   10   continue
        j2 = j2+1
        if (strn(j2:j2) .eq. dc) goto 10
        j1 = min(len(strn),j2)
        call nwordg(strn,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then
          elseif (strn(j1:j1+2) .eq. 'fn=')  then
            if (j1+3 .gt. j2) call rx('pldos: bad file name')
            fn = strn(j1+3:j2)
            k = fopnx(strn(j1+3:j2),172,-1,-1)
            if (k .eq. 0) then
              call rxs2('pldos: missing file "',strn(j1+3:j2),'"')
            endif
            ifr = fopng(strn(j1+3:j2),-1,1)
          elseif (strn(j1:j1+4) .eq. 'chan=')  then
            if (j1+5 .gt. j2) call rx('pldos: bad list, chan=..')
            call mkils0(strn(j1+5:j2),nchan,ix)
            if (k .gt. 100) call rx('increase maxnum')
            call mkilst(strn(j1+5:j2),k,ix)
            if (k .lt. 0) call rx('pldos: bad list, chan=..')
          elseif (strn(j1:j2) .eq. 'scale')  then
            iscale = 1
          else
            call rxs('pldos: fail to recognize argument ',
     .      strn(j1:j2))
          endif
          goto 10
        endif
        if (ifr .eq. 0) then
          call info0(0,0,0,
     .    ' pldos: -ref option requires a file be specified')
          goto 199
        endif
        if (nchan .eq. 0) then
          nchan = 1
          ix(1) = 1
        endif
C#ifdef F90
        rewind ifr
        read(ifr,*) eminr,emaxr,ner,nd,nspin
        rewind ifr
        nemx = ner
        ndinmx = nd*nspin
        allocate(dref(nemx),din(nemx,ndinmx))
C#endif
        call dosio(din,nemx,ndinmx,ner,nd,eminr,emaxr,nspin,
     .  efurm,dul,ifmt,ifr)
C       Assemble reference DOS
        call dpzero(dref,ner)
        do  i = 1, nchan
          k = ix(i)
          if (k .gt. nd*nspin) call rxi(
     .    'pldos: reference channel outside file''s range 1 ..',
     .    nd*nspin)
          call daxpy(ner,1d0,din(1,k),1,dref,1)
        enddo
C#ifdef F90
        deallocate(din)
C#endif
        call fclose(ifr)
        ifr = 1
      endif

C#ifdef F90
      rewind ifi
      read(ifi,*) emun,emux,ne,nd,nspin
      rewind ifi
      nemx = ne
      ndinmx = nd*nspin
      allocate(e(nemx),dtot(nemx),din(nemx,ndinmx),d(nemx,2,ndmx))
C#endif
      call dosio(din,nemx,ndinmx,ne,nd,emun,emux,nspin,
     .efurm,dul,ifmt,ifi)
      call fclose(ifi)
      nd = nd*nspin

      if (ifr .ne. 0) then
        call info0(20,0,1,' pldos: '//
     .  ' verifying reference and dos files have equivalent mesh ...')
        call isanrg(ner,ne,ne,' ','ne in reference file',.true.)
        call fsanrg(eminr,emun,emun,0d-9,' ','emin in reference file',
     .  .true.)
        call fsanrg(emaxr,emux,emux,0d-9,' ','emax in reference file',
     .  .true.)
      endif

      ltdos = '1,bold=3,col=0,0,0'
      if (iscale .ne. 0) ltdos = '1,bold=3,col=1,0,0,colw=0,0,1'


      eshift = 0
      if (cmdopt('-ef=',4,0,strn)) then
        ie = 4
        if (.not. a2bin(strn,eshift,4,0,' ',ie,-1))
     .  call rxs('fplot: failed to parse argument  ',strn)
        eshift = efurm-eshift
      elseif (cmdopt('-shiftef=',9,0,strn)) then
        ie = 9
        if (.not. a2bin(strn,eshift,4,0,' ',ie,-1))
     .  call rxs('fplot: failed to parse argument  ',strn)
      endif
      emin = emun-eshift
      emax = emux-eshift
      eferm = efurm-eshift

      if (cmdopt('-escl=',6,0,strn)) then
        ie = 6
        if (.not. a2bin(strn,escl,4,0,' ',ie,-1))
     .  call rxs('fplot: failed to parse argument  ',strn)
        emin = emin*13.6
        emax = emax*13.6
        eferm = eferm*13.6
      endif

      call awrit5(' pldos:  %i dos read  ef=%d  emin=%d  emax=%d  ne=%i'
     .,' ',80,i1mach(2),nd,eferm,emin,emax,ne)
      if (ne .gt. nemx) stop '**** DIMENSION NEMX TOO SMALL'
      if (nd .gt. ndinmx) stop '**** DIMENSION NDINMX TOO SMALL'
      do  1  ie = 1, ne
    1 e(ie) = emin + (ie-1d0)*(emax-emin)/(ne-1d0)

C --- Select the next linear combination of dos ---
C     ifr should be zero, or 1 if reference dos is 1st dos
      md = 0
      dssum = 0

C ... Reference case: dref becomes the first DOS
      if (ifr .eq. 1) then
        call dcopy(ne,dref,1,d,1)
        call idmxmn(ne,dref,1,i,j)
        dmx(1) = dref(j)
        dssum = dssum+dmx(md+ifr)
      endif

   90 continue
      print *, ' '
      call awrit1(' dos %i:  Enter list of dos to be added:',strn,80,
     .0,md+1)

      if (cmdopt('-lst=',5,0,strn)) then
        call wordg(strn,1,'; ',md+1,is1,is2)
        doslst = ' '
        is1 = max(is1,6)
        if (is2 .ge. is1) doslst = strn(is1:is2)
      else
        call skpblb(strn,len(strn),i)
        call cwrite(strn,0,i+1,0)
        read(*,'(a512)') doslst
      endif
      if (cmdopt('-lst2=',6,0,strn)) then
        call wordg(strn,1,'; ',md+1,is1,is2)
        dosls2 = ' '
        is1 = max(is1,7)
        if (is2 .ge. is1) dosls2 = strn(is1:is2)
      elseif (cmdopt('-lst2',5,0,strn)) then
        dosls2 = '*'
      else
        dosls2 = ' '
      endif

      is = 0
      nnum = 0
      num2 = 0
      call skipbl(doslst,512,is)
      if (is .ge. 512) goto 32
      call mkilst(doslst,nnum,num)
   32 continue
      if (dosls2 .eq. '*') then
        if (nnum .eq. 0) goto 33
        num2 = nnum
        do  i = 1, nnum
          num(nnum+i) = num(i)+1
        enddo
        nnum = 2*nnum
        l2chan = .true.
        goto 33
      endif
      call skipbl(dosls2,512,is)
      if (is .ge. 512) goto 33
      call mkilst(dosls2,num2,num(1+nnum))
      nnum = nnum + num2
      if (num2 .gt. 0) l2chan = .true.
   33 continue

      if (num2 .ne. 0 .and. ifr .ne. 0)
     .call rx('pldos: -ref not set up for 2nd channel')
      num1 = nnum - num2
      if (cmdopt('-fplot',6,0,strn) .and. num2 .gt. 0) then
        ifi = fopng('dosp2.'//datnam,11,0)
      endif
      if (nnum .eq. 0) goto 98
      do  40  i = 1, nnum
        if (num(i) .gt. nd) call fexit2(-1,1,' pldos (abort) dos %i '/
     .  /'exceeds number read (%i)',num(i),nd)
   40 wgt(i) = 1d0
      if (cmdopt('-lst=',5,0,strn)) then
        call awrit3(
     .  ' combining %i dos: %n:1i',' ',80,i1mach(2),nnum,nnum,num)
        strn = ' '
      else
        i = awrite(' Weights for dos:%n:1i?',strn,80,0,
     .  nnum,num,0,0,0,0,0,0)
        call cwrite(strn,0,i+1,0)
        read(*,'(a512)') strn
      endif
      if (strn .eq. ' ' .or. strn .eq. '/') goto 42
      is = 0
      ist = a2vec(strn,len(strn),is,4,', ',2,3,nnum,ix,wgt)
      if (ist .le. 0) then
        print *, 'bad weights, try again'
        goto 90
      endif
      if (ist .lt. nnum)
     .call awrit2('w=%n:1d',' ',80,i1mach(2),nnum,wgt)
   42 continue
C       call awrit1('%Nsum w = %d',' ',80,i1mach(2),dsum(nnum,wgt,1))
      md = md+1
      if (md .gt. ndmx) call fexit(-1,1,
     .' Exit -1:  no. dos exceeds ndmx=%i',ndmx)
C   ... Accumulate the linear combination of dos
      call dpzero(d(1,1,md+ifr),ne)
      call dpzero(d(1,2,md+ifr),ne)
      do  44  i = 1, num1
   44 call daxpy(ne,wgt(i),din(1,num(i)),1,d(1,1,md+ifr),1)
      do  45  i = num1+1, nnum
   45 call daxpy(ne,-wgt(i),din(1,num(i)),1,d(1,2,md+ifr),1)

      if (cmdopt('-dif',4,0,strn)) then
        if (num2 .ne. 0) call rx('pldos: -dif not set up for num2')
        npdif = 2
        lrat = .false.
        if (strn(5:5) .eq. ':') then
          is = 5
          ist = a2vec(strn,len(strn),is,4,', ',2,2,2,ix,xxv)
          if (ist .lt. 0) goto 199
          npdif = nint(xxv(1))
          if (ist .gt. 1 .and. xxv(2) .gt. 0) lrat = .true.
        endif
        call dcopy(ne,d(1,1,md+ifr),1,dtot,1)
        call poldvm(e,dtot,ne,npdif,lrat,1d-6,lerr,d(1,1,md+ifr))
        if (lerr.ne.0) call rx('pldos: error in rational poly interp')
C   ... Integration by the trapezoidal rule
      elseif (cmdopt('-int',4,0,strn)) then
        if (num2 .ne. 0) call rx('pldos: -int not set up for num2')
        if (ifr .ne. 0) call rx('pldos: -int not set up for ifr')
        call dcopy(ne,d(1,1,md+ifr),1,dtot,1)
        d(1,1,md+ifr) = 0d0
        de = e(2) - e(1)
        do  47  ie = 2, ne
          d(ie,1,md+ifr) = d(ie-1,1,md+ifr) +
     .    (dtot(ie)+dtot(ie-1))*de/2
   47   continue
      endif
      dmx(md+ifr) = 0d0
      do  46  ie = 1, ne
   46 dmx(md+ifr) =
     .dmax1(dmx(md+ifr),abs(d(ie,1,md+ifr)),abs(d(ie,2,md+ifr)))
      dssum = dssum+dmx(md+ifr)
      goto 90

C --- Find the first nonzero dos ---
   98 continue
      mdi = md+ifr
      if (mdi .eq. 0) call rx('pldos: no dos to plot')
      do  80  ie = 1, ne
        istart = ie
        top = 0d0
        do  81  id = 1, mdi
   81   top = dmax1(top,d(ie,1,id))
        if (top .gt. 1d-5) goto 82
   80 continue
   82 estart = emin + (istart-1)*(emax-emin)/(ne-1d0)

C --- Do the plotting ---
      print *, ' '
      call awrit5(' %i plots:  emin=%d  estart=%d  emax=%d  ef=%d',
     .strn,72,i1mach(2),mdi,emin,estart,emax,eferm)
      call awrit3(' dmx=%n:1,2;2d  sum=%,2;2d',' ',80,
     .i1mach(2),mdi,dmx,dssum)
      dtopmx = 0d0
      do  87  id = 1, mdi
   87 dtopmx = dmax1(dtopmx,dmx(id))
      cmhgt = min(ymax/dssum*dtopmx,5d0)
      e1 = emin
      e2 = emax
C     print 335, 'Enter dmax, height of dmax (cm), e1, e2: '
C 335 format(1x,a)
      call cwrite(' Enter dmax, height of dmax (cm), e1, e2: ',0,41,0)
      read (5,*)  dtopmx,cmhgt,e1,e2
      ne1 = max0(nint((e1-emin)*ne/(emax-emin)+1),1)
      ne2 = min0(nint((e2-emin)*ne/(emax-emin)),ne)
      do  58  id = 1, mdi
        dmx(id) = 0d0
        do  58  ie = ne1, ne2
   58 dmx(id) = dmax1(dmx(id),abs(d(ie,1,id)),abs(d(ie,2,id)))
      if (e1 .ne. emin .or. e2 .ne. emax) call awrit4(
     .'dmx[%;3d,%;3d] =%n:1,2;2d',' ',80,i1mach(2),e1,e2,mdi,dmx)
      print *, ' '

      ypos1 = 0.97d0

C ... Write files for fplot
      if (cmdopt('-fplot',6,0,strn)) then

C       Scale columns 2...md+ifr
        if (ifr .eq. 1 .and. iscale .ne. 0) then
          do  ie = ne1, ne2
            if (d(ie,1,1) .ne. 0) then
              do  id = 2, mdi
                d(ie,1,id) = d(ie,1,id)/d(ie,1,1)
              enddo
            else
              do  id = 2, mdi
                d(ie,1,id) = 1
              enddo
            endif
          enddo
        endif

        ifi = fopng('dosp.'//datnam,-1,0)
        ypos1 = 1.4d0
        call awrit2('%% rows %i cols %i',strn,80,ifi,ne2-ne1+1,mdi+1)
        do  60  ie = ne1, ne2
          write(ifi,333) e(ie), (d(ie,1,id), id=1, mdi)
  333     format(f9.5,6f11.4:/(9x,6f11.4))
   60   continue
        call fclose(ifi)
        if (l2chan) then
          ifi = fopng('dosp2.'//datnam,-1,0)
          call awrit2('%% rows %i cols %i',strn,80,ifi,ne2-ne1+1,mdi+1)
          do  61  ie = ne1, ne2
            write(ifi,333) e(ie), (d(ie,2,id), id=1, mdi)
   61     continue
          call fclose(ifi)
        endif
      else
        call pltini('ps',0,.false.)
        call setfnt(font,fontsz)
      endif

      outs(1) = ' '
      if (cmdopt('-fplot',6,0,strn)) then
        lgcmd = fopng('plot.dos',-1,0)
        if (ltdos .ne. ' ')
     .  write(lgcmd,'(''%char0 ltdos="'',a,''"''/)') trim(ltdos)

        call awrit0('fplot',' ',-80,lgcmd)
      endif
      if (iscale .ne. 0) mdi = 1
      do  20  id = 1, mdi
C     Round off dtop to a multiple of nrnde
        idtop = int(dmx(id)*1.05d0/nrnde + 1d0)
        dtop = idtop*nrnde
        dtop = dmin1(dmax1(dtop,5d0),dtopmx)
        dbot = 0
        if (l2chan) dbot = -dtop
        dcm = dtop*cmhgt/dtopmx
        dy = dcm/ymax
        ypos2 = ypos1
        ypos1 = ypos2-dy
        if (ypos1 .lt. 0)
     .  print *, 'pldos (warning): frame past bottom of page ...'
        call awrit6(' plot %i  dmx=%1,3;3d  dtop=%1,3;3d  dy=%1,3;3d'//
     .  '  y=%1,3;3d %1,3;3d',' ',80,i1mach(2),
     .  id,dmx(id),dtop,dy,ypos1,ypos2)
        if (cmdopt('-fplot',6,0,strn)) then
C       First frame: set up variables
          if (iouts .eq. 1) then
            iouts = iouts+1
            call awrit5('%% var ymax=%1;3d dy=%1;3d dw=.00 ymax+=dy '//
     .      'emin=%1;3d emax=%1;3d ef=%d',
     .      outs(iouts),80,0,ypos2,ypos2-ypos1,e1,e2,eferm)
            iouts = iouts+1
          endif
C        call skpblb(outs(iouts),80,i)
C        if (i .ge. 0) iouts = iouts+1
          iouts = iouts+1
          call awrit3('%% var ymax-=dy+dw dy=%1;3d dmin=%1;3d dmax=%1;3d',
     .    outs(iouts),80,0,ypos2-ypos1,dbot,dtop)
          iouts = iouts+1
          call awrit1('%a -frme 0,1,{ymax-dy},{ymax} -p0 -x {emin},{emax}'
     .    //' -y {dmin},{dmax} -1p%?#n# -noxn##',outs(iouts),80,0,id-md)
          iouts = iiouts(outs,iouts,72)
          if (iscale .ne. 0) then
            call awrit1('%a -colsy %i -colsw 3'//' -lt {ltdos} -qr dosp.'
     .      //trim(datnam),outs(iouts),250,0,id+1)
          else
            call awrit1('%a -colsy %i -lt {ltdos} -qr dosp.'//
     .      trim(datnam),outs(iouts),250,0,id+1)
          endif
          iouts = iiouts(outs,iouts,72)
          if (iouts .gt. 58-1) call rx('pldos output string too long')
          if (l2chan) then
            call awrit1('%a -colsy %i -lt {ltdos} -qr dosp2.'//
     .      trim(datnam),outs(iouts),250,0,id+1)
            iouts = iiouts(outs,iouts,72)
            if (iouts .gt. 58-1) call rx('pldos output string too long')
          endif
          if (eferm .lt. e2 .and. eferm .gt. e1) then
            call awrit0('%a -lt 2,bold=3,col=0,0,0,2,.5,.05,.5 '//
     .      '-tp 2~{ef},{dmin},{ef},{dmax}',outs(iouts),250,0)
            iouts = iiouts(outs,iouts,72)
          endif
        else
          call pltstp(5d-2,0.98d0,ypos1,ypos2)
          call pltstu(e1,e2,0d0,dtop,0d0)
          dytk = dytik
          if (dtop .lt. 10.01d0) dytk = 1d0
          if (dtop .gt. 200d0) dytk = 50d0
          xnum = 0
          ynum = 0
          if (id .eq. mdi) then
            xnum = 1
            ynum = 1
          endif
          xlab = ' '
          ylab = ' '
          title = ' '
          call frme(-1d30,dxtik,1,.025d0,.6d0,0,xnum,xnfmt,xlab,
     .    -1d30,dytk ,1,.025d0,.6d0,0,ynum,ynfmt,ylab,
     .    ' ',xxv,xxv,xxv,title,ltfrme,fgs)
          call plntyp(1,3,2d0,0.5d0,0d0,0d0)
          call plcrv(e(ne1),d(ne1,1,id),e,e,e,ne2-ne1+1,e1,e2,0d0,dtop,0,
     .    -2d0,csw,11)
          if (eferm .lt. e2 .and. eferm .gt. e1) then
            xd(1) = eferm
            xd(2) = eferm
            yd(1) = 0
            yd(2) = dtop
            call plntyp(2,3,2d0,1d0,0d0,0d0)
            call plcrv(xd,yd,xd,xd,xd,2,xd(1)-0.1d0,xd(2)+0.1d0,yd(1),
     .      yd(2),0,-2d0,csw,11)
          endif
          print *, ' '
        endif
        if (iscale .ne. 0) goto 21
   20 continue
   21 continue
      if (cmdopt('-fplot',6,0,strn)) then
C       call awrit1('%a',outs,80,-i1mach(2),0)
        lgcmd = fopng('plot.dos',-1,0)
        do  25  is = 1, iouts
          call skpblb(outs(is),80,i)
C         call cwrite(outs(is),0,i,1)
          write(lgcmd,'(a)') outs(is)(1:i+1)
   25   continue
        call fclose(lgcmd)
        print *, ' ... invoke:'
        print *, ' fplot -disp -pr10 -f plot.dos'
      else
        call pltdmp(0)
      endif
C#ifdef F90
      deallocate(e,dtot,din,d)
C#endif
      call cexit(0,1)

   99 continue
      call rx('pldos encountered error in attempting to open file')
  199 print *, 'usage:  pldos ',
     .'[-fplot] [-lst=list] [-lst2|-lst2=list] [-dif:#[,rat]]'
      print *,
     .'              [-ref:fn=filename[:chan=ref-channel][:scale]]'
      print *, '              [-ext=fplot-data-file-extension]'
      print *, '              [-ef=#|-shiftef=#] [-escl=#] [-int] fname'
      call cexit(-1,1)
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
      subroutine dosio(dos,nemx,ndmx,ndos,nld,emin,emax,nspin,
     .eferm,del,ifmt,ifile)
C- I/O for dos, MSM's format
      implicit none
C Passed parameters
      integer ifile,ifmt,ndmx,ndos,nemx,nld,nspin
      double precision del,eferm,emax,emin
      double precision dos(nemx,1)
C Local parameters
      integer ie,ild,iedum

C --- Write branch ---
      if (ifile .lt. 0) then
        write(-ifile,760) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ifmt .eq. 0) then
          do  10  ie = 1, ndos
   10     write(-ifile,761) ie,(dos(ie,ild),ild=1,nld*nspin)
        else
          do  11  ild = 1, nld*nspin
   11     write(-ifile,762) (dos(ie,ild),ie=1,ndos)
        endif
  761   format(i5,6f12.5/(5x,6f12.5))
  760   format(2f10.5,3i5,2f10.5,i5)
  762   format(5f14.6)
      endif
C --- Read branch ---
      if (ifile .gt. 0) then
        read(ifile,760) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ndos .gt. nemx) call rxi('dosio: ndos.gt.nemx =',nemx)
        if (nld*nspin .gt. ndmx) call rxi('dosio: nld > ndmx=',ndmx)
        if (ifmt .eq. 0) then
          do  20  ie = 1, ndos
   20     read(ifile,761) iedum,(dos(ie,ild),ild=1,nld*nspin)
        elseif (ifmt .eq. 1) then
          do  21  ild = 1, nld*nspin
   21     read(ifile,762) (dos(ie,ild),ie=1,ndos)
        else
          call rx('dosio: bad fmt')
        endif
      endif
      end
