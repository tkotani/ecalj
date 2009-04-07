      subroutine aioxtn(sctrl,sspec,ssite,sarray,slat,slabl,pnuc,qnuc)
C- I/O moments and pp's for all classes, other formats
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbasp nclass nl nspin lrs
Ci     Stored:    lrs
Ci     Passed to: iors
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa rmt z
Ci     Stored:    name
Ci     Passed to: uspec dval spacks iors
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: pnu
Ci     Stored:    class spec
Ci     Passed to: spackv asars1 usite dval iors aioxt1
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: oclabl oics onrcp
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat vol opos awald nkd nkq odlv oqlv
Ci     Stored:    *
Ci     Passed to: ulat dval iors
Ci   pnuc  :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci          This array is class-based
Ci   qnuc  :energy-weighted moments of the sphere charges
Ci          This array is class-based
Co Outputs
Co   *When importing from atom files, output are updated
Co    atom files.
Co   *When importing from rst files, output is pnuc,qnuc
Cl Local variables
Cl         :
Cb Bugs
Cb   When importing from atom files,
Cb   file importing for layer case (nbasp>nbas) fails to generate
Cb   proper Madelung potential.
Cr Remarks
Cr
Cu Updates
Cu   04 Nov 04 Upgrade of functionality 
Cu   05 Mar 02 added to the restart file editor
Cu   01 Mar 02 First cut at import mode from rst files (ASA mode)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sctrl(1),sspec(1),ssite(1),slat(1),sarray(1)
      double precision pnuc(*),qnuc(*)
      character*8 slabl(*)
C ... Local parameters
      integer ic,ifi,fopna,fopnx,nr,lmx,nsp,idmod(9),nl,nbas,nbasp,
     .  nbasn,nclass,j1,j2,nkd,nkq,i,j,k,iors,mode,lgunit,stdo,lio,
     .  iosits,nglob,nspec,nati
      double precision a,rmax,vrmax(2),rhrmx(2),ves,dval,
     .  z,qc,dq,sumec,sumev,thrpv,ekin,utot,rhoeps,etot(26),
     .  pnl(2*9),qnl(3*2*9),pp(3*2*9),evrl(9*2),dsum,wk(10)
      double precision awald,alat,vol,plat(9),tau(3),emad,
     .  trumad,vmtz(2),dglob
      integer odlv,oqlv,onrcp,oves,ospot,oclabl,oics,ov,is,oipc,oips,
     .  oqt,ormax,omad,opos,ovrmax,ossite,osspec,oslat,osbz,nbasi,a2vec,
     .  nspeci,nit,irs,igets,opnui,oqnui,n0,opnus,oqnus,oipa,nlsts,nw,
     .  nlstd,ip,lsts(10000),lstd(10000),oz,olmx,oqcs,oqts,oamoms,owk
      integer opos1,opos2,opos3
      parameter (n0=10)
      character outs*80, clabl*8, alabel*4, fnam*80
      logical sw,cmdopt,lrel,lgen,lmom,lpp,lpot,lqtot,lves,a2bin
      logical aiogen,aiomom,aiopar,aiopot
      logical io4g,io4mm,io4pp,io4pt,io3smm,lstutt,lrs
      logical io4sg,io4smm,io4spp,io4spt
      logical io5sg,io5smm,io5spp,io5spt
      logical io47um,io47up
      character*32 jobid
C ... Heap
      integer w(1)
      common /w/ w

      call upack('array oclabl oics onrcp',sarray,oclabl,oics,onrcp,0,0)
      call upack('ctrl nbas nbasp nclass nl nspin',sctrl,nbas,nbasp,
     .  nclass,nl,nsp)
      call upack('lat alat plat',slat,alat,plat,0,0,0)
      nspec = nglob('nspec')
      call defi(oipc,nbasp)
      call spackv(10,'site class',ssite,1,nbasp,w(oipc))
C     call upack('lat alat vol opos',slat,alat,vol,opos,0,0)
      stdo = lgunit(1)

      call dpzero(vrmax,2)
c     if charges at all sites are found, lqtot=.true.
      lqtot = .true.
      lves = .false.
      call defrr(oqt,-nclass)
      call defrr(ormax,-nclass)
      call defrr(ovrmax,-nclass*2)
      call defrr(oves,-nclass)

      lrs = cmdopt('-rs',3,0,outs)

c   ... import from rst files
        if (lrs) then

          irs = igets('ctrl lrs',sctrl)

          call defi(oips,nbasp)
          call spackv(10,'site spec',ssite,1,nbasp,w(oips))
          call sp2cls('spec z',sspec,w(oips),1,1,nbasp,oz)
          call sp2cls('spec lmxa',sspec,w(oips),1,1,nbasp,olmx)

          call defdr(ospot,1000)
          call upot(w(ospot),-1,-000,0,0,0)
          call redfrr(ospot,nint(dval(w(ospot),1)))

          call defdr(osbz,1000)
          call ubz(w(osbz),-1,-000,0,0,0)
          call redfrr(osbz,nint(dval(w(osbz),1)))

          call defdr(oslat,1000)
          call ulat(w(oslat),-1,0,0,0,0)
          call redfrr(oslat,nint(dval(w(oslat),1)))

          call defrr(opnus,-n0*nsp*nbasp)
          call defrr(oqnus,-3*n0*nsp*nbasp)
          call defi(oipa,nbasp)
          call asars1(30,ssite,nbasp,nsp,nl,ov,w(opnus),w(oqnus),
     .      w(oipc),w(oipa),pnuc,qnuc)
          call rlse(oipa)


   99     print 301
  301     format(//
     .      ' Welcome to the restart file editor.')
  100     continue
          print '('' Option (enter ? for menu) : '',$)'
          outs = ' '
          read(*,'(a80)') outs
          call locase(outs)
          if (.false.) then
          elseif (outs .eq. ' ') then
            print '('' Enter ''''q'''' to exit, ''''a'''' to abort'')'
            goto 100
          elseif (outs(1:1) .eq. 'a') then
            call rx0('exiting aioxtn without rewriting rst file')
          elseif (outs(1:1) .eq. 'q') then
            goto 199
          elseif (outs(1:1) .eq. 'i') then
            goto 101
          elseif (outs(1:3) .eq. 'dup') then
            call words(outs,nw)
            if (nw .eq. 1) then
              print 315
  315         format(/' Enter the lower and upper range of sites ',
     .          'to duplicate: ',$)
              read(*,*) i,j
            elseif (nw .eq. 2 .or. nw .eq. 3) then
              ip = 3
              k = a2vec(outs,len(outs),ip,2,', ',2,-3,-2,lsts,lstd)
              if (k .lt. 0) then
                print 306, outs(1:ip+1)
                goto 100
              elseif (k .eq. 1) then
                i = lstd(1)
                j = i
              elseif (k .eq. 2) then
                i = lstd(1)
                j = lstd(2)
              endif
            endif
            if (j .lt. i) then
              print *, 'upper range < lower range ... skipping'
              goto 100
            endif

            nbasn = nbasp + j-i+1
C           duplicate ssite
            k = nint(ssite(1))
            call defrr(ossite,k*nbasn)
            call dcopy(k*nbasp,ssite,1,w(ossite),1)
            call dmscop(ssite,k,w(ossite),k,1,k,i,nbasp,1,j+1,1d0)
            call aioxt3(k,w(ossite),nbasp,ssite,nbasn)
C           duplicate pnu
            k = n0*nsp
            owk = opnus
            call defrr(opnus,k*nbasn)
            call dcopy(k*nbasp,w(owk),1,w(opnus),1)
            call dmscop(w(opnus),k,w(owk),k,1,k,i,nbasp,1,j+1,1d0)
            call aioxt3(k,w(owk),nbasp,w(opnus),nbasn)
C           duplicate qnu
            k = 3*n0*nsp
            owk = oqnus
            call defrr(oqnus,k*nbasn)
            call dcopy(k*nbasp,w(owk),1,w(oqnus),1)
            call dmscop(w(oqnus),k,w(owk),k,1,k,i,nbasp,1,j+1,1d0)
            call aioxt3(k,w(owk),nbasp,w(oqnus),nbasn)
C           Update nbas,nbasp
            nbas  = nbas  + j-i+1
            nbasp = nbasn
            i = dglob('nbasp',dble(nbasp),1)
            i = dglob('nbas',dble(nbas),1)
            call pack5('ctrl nbas nbasp',sctrl,nbas,nbasp,0,0,0)
C           Update ipc,ips,z,lmxa
            call defi(oipc,nbasp)
            call spackv(10,'site class',ssite,1,nbasp,w(oipc))
            call defi(oips,nbasp)
            call spackv(10,'site spec',ssite,1,nbasp,w(oips))
            call sp2cls('spec z',sspec,w(oips),1,1,nbasp,oz)
            call sp2cls('spec lmxa',sspec,w(oips),1,1,nbasp,olmx)

            print '('' basis enlarged to'',i5,'' atoms'')', nbasp
            goto 100

          elseif (outs(1:5) .eq. 'shift') then
            call words(outs,nw)
            if (nw .eq. 1) then
              print 311
  311         format(/' Enter the lower and upper range of sites ',
     .          'to shift: ',$)
              read(*,*) i,j
            elseif (nw .eq. 2 .or. nw .eq. 3) then
              ip = 5
              k = a2vec(outs,len(outs),ip,2,', ',2,-3,-2,lsts,lstd)
              if (k .lt. 0) then
                print 306, outs(1:ip+1)
                goto 100
              elseif (k .eq. 1) then
                i = lstd(1)
                j = i
              elseif (k .eq. 2) then
                i = lstd(1)
                j = lstd(2)
              endif
            endif
            if (j .lt. i) then
              print *, 'upper range < lower range ... skipping'
              goto 100
            endif
            print 314
  314       format(/' Enter dx dy dz: ',$)
            read(*,*) tau
            call defrr(opos,3*nbasp)
            call spackv(10,'site pos',ssite,1,nbasp,w(opos))
            call aioxt2(i,j,tau,w(opos))
            call spackv(11,'site pos',ssite,1,nbasp,w(opos))
            call rlse(opos)
            goto 100

          elseif (outs(1:1) .eq. 's') then


            call defrr(oqcs,nbasp)
            call defrr(oqts,nbasp)
            call defrr(oamoms,nbasp)
            call getq(nsp,n0,w(olmx),nbasp,w(oz),w(opnus),w(oqnus),
     .        w(oics),sspec,w(oqcs),w(oqts),w(oamoms))

            call defrr(opos,3*nbasp)
            call spackv(10,'site pos',ssite,1,nbasp,w(opos))
            call defrr(opos1,nbasp)
            call defrr(opos2,nbasp)
            call defrr(opos3,nbasp)
            call dmscop(w(opos1),1,w(opos),3,1,1,1,nbasp,1,1,1d0)
            call dmscop(w(opos2),1,w(opos),3,2,2,1,nbasp,1,1,1d0)
            call dmscop(w(opos3),1,w(opos),3,3,3,1,nbasp,1,1,1d0)

            outs = ' Site   Z      x         y         z'//
     .        '         dq     amom'
            print *
C            call arrprt(outs,'%,4i%:-3,4;1d%:-3,4;4d%:-2,4;4d',
            call arrprt(outs,
     .        '%,4i%;6,1D%;10,5D%;10,5D%;10,5D%:-3,4;4d%:-2,4;4d',
     .        'Idddddd',nbasp,0,2,0,'  | ',w,
     .        w(oz),w(opos1),w(opos2),w(opos3),w(oqts),w(oamoms),w)
            wk(1) = dsum(nbasp,w(oqts),1)
            wk(2) = dsum(nbasp,w(oamoms),1)
            call awrit2(' sum-of-sphere charges=%:,6;6d'//
     .        '  sum-of-sphere moments=%:,6;6d',' ',80,
     .        stdo,wk,wk(2))
            call rlse(oqcs)
            print *, ' '
            goto 100

          elseif (outs .eq. 'wsite') then
            lio = 1000*(2+4+8+16+32) + 1
            if (iosits(lio,3d0,0,'site',ifi,slabl,alat,plat,nbasp,nspec,
     .        sspec,ssite).lt. 0) call rx('failed to write ssite')
            goto 100

          elseif (outs .eq. '?') then
            print 310
  310       format(/' Select one of these options:'/
     .        t4,'i',t10,'to import a portion of a restart file'/
     .        t4,'dup',t10,
     .                   'to duplicate some site data, enlarging basis'/
     .        t4,'s',t10,'to show some site data'/
     .        t4,'shift',t10,'to shift some site positions'/
     .        t4,'wsite',t10,'to write site file'/
     .        t4,'q',t10,'to quit the editor'/
     .        )
            goto 100

          else
            print '(1x,''unrecognized option: '',a10)', outs
            goto 100

          endif

C     ... Import some pnu,qnu from a rst file
  101     continue
          fnam = ' '
          print 302
  302     format(/' Enter the full name of a restart file ',
     .            '(enter nothing to exit): ',$)
          read(*,'(a80)') fnam
          if (fnam .eq. ' ') goto 100
          ifi = fopnx(fnam,72,-1,-1)
          if (ifi .ne. 1) then
            print *, 'sorry, file doesn''t seem to exist.  try again.'
            goto 101
          endif
C         call pshpr(110)
          ifi = fopnx(fnam,72,0,-1)
C         call poppr

          i = iors(-1,sctrl,w,w,w(oslat),w,w(osbz),
     .      jobid,nbasi,nati,nspeci,w,w,nit,.false.,ifi)

          call word(fnam,1,j1,j2)
          print 303, fnam(j1:j2), nbasi, nspeci
  303     format(' file ',a,' contains nbas =',i4,' and nspec =',i4)

C         setup for file read ... much of the data will be discarded
          call pack1('ctrl lrs',sctrl,0)

          call defrr(ossite,1000)
          call usite(w(ossite),-1,0,1,0,0,0)
          call redfrr(ossite,nbasi*nint(dval(w(ossite),1)))
          call usite(w(ossite),-1,0,nbasi,0,0,0)
          call spackv(1,'site spec',w(ossite),1,nbasi,-1)

          call defrr(osspec,1000)
          call uspec(w(osspec),-1,0,1,0,0,0)
          call redfrr(osspec,nspeci*nint(dval(w(osspec),1)))
          call uspec(w(osspec),-1,-000,nspeci,0,0,0)
          do  i = 1, nspeci
            call spacks(1,'spec name',w(osspec),'zz',i,i)
          enddo
C         call defrr(opnui,-1*n0*nsp*nbasi)
          call defrr(oqnui,-3*n0*nsp*nbasi)

          i = iors(0,sctrl,w(ossite),w(osspec),w(oslat),w(ospot),
     .      w(osbz),jobid,nbasi,nati,nspeci,w,w(oqnui),nit,.false.,ifi)

          call pack1('ctrl lrs',sctrl,irs)
          call fclose(ifi)


  110     print 304, fnam(j1:j2)
  304     format(/' Enter list of sites from input file ',a,
     .      ' (? for help) : ',$)
          read(*,'(a80)') outs

          if (outs .eq. ' ') then
            goto 101

          elseif (outs .eq. '?') then
          print 305
  305     format(/
     .      ' If you enter a single number, you can poke data from'/
     .      ' this site into a group of sites.'/
     .      ' If you enter a list, the list of destination sites'/
     .      ' you enter following this will have to be the same length'/
     .      ' as the list you enter now.'/
     .      ' See documentation for site-list syntax.')

          goto 110
          endif

          call mkils0(outs,nlsts,ip)
          if (nlsts .lt. 0) then
            print 306, outs(1:ip+1)
  306       format(' Failed to parse string ',a,' ... try again.')
            goto 101
          endif
          call mkilst(outs,nlsts,lsts)
          print 307, outs(1:ip+1), nlsts
  307     format(' List ',a,' contains',i4,' sites.')


  120     print 308
  308     format(/' Enter list of destination sites : ',$)
          read(*,'(a80)') outs

          if (.false.) then
          elseif (outs .eq. ' ' .or.  outs .eq. '?') then
            print 312
  312       format(/' Enter ''q'' to quit'/
     .        ' Enter ''?'' to see this message'/
     .        ' Otherwise, enter a site-list ',
     .        ' (see documentation for syntax)')
            if (nlsts .eq. 1)  print 313
  313       format(2x,' Because you imported only one site from the',
     .        ' input file,'/
     .        2x,' you can enter an implicit list with z=#')
            goto 120
          elseif (outs .eq. 'q') then
            goto 99
          elseif (outs(1:2) .eq. 'z=') then
            if (nlsts .ne. 1) then
              print *, 'Sorry, z= option allowed only when you import',
     .          ' from one site'
              goto 120
            endif
            i = 2
            if (.not. a2bin(outs,wk,4,0,' ',i,len(outs)-1)) then
              print *, 'Sorry, couldn''t parse string ... try again'
              goto 120
            endif
            nlstd = 0
            do  i = 1, nbasp
              if (dval(w(oz),i) .eq. wk(1)) then
                nlstd = nlstd + 1
                lstd(nlstd) = i
              endif
            enddo
            call word(outs,1,j1,j2)
            print 307, outs(j1:j2), nlstd
            if (nlstd .eq. 0) then
              print *, 'Sorry, no sites with this z ... try again'
              goto 120
            endif

          else

            call mkils0(outs,nlstd,ip)
            if (nlstd .lt. 0) then
              print 306, outs(1:ip+1)
              goto 101
            endif
            call mkilst(outs,nlstd,lstd)
            print 307, outs(1:ip+1), nlstd
          endif

          if (nlsts .gt. 1 .and. nlsts .ne. nlstd) then
            print 309
  309       format(' Size of lists do not match ... try again.')
            goto 101
          endif


C         Poke imported site->pnu and qnui to pnus and qnus
          mode = 1
          if (nlsts .eq. 1) mode = 0
          call aioxt1(mode,nlstd,lsts,lstd,nsp,nbasi,w(ossite),
     .      w(oqnui),nbasp,w(opnus),w(oqnus))

          goto 101


C     ... Cleanup and exit rst editor
  199     continue
          call defi(oipa,nbasp)
          call asars1(31,ssite,nbasp,nsp,nl,ov,w(opnus),w(oqnus),
     .      w(oipc),w(oipa),pnuc,qnuc)
          call rlse(oipc)

          return

        endif

c     call pshpr(110)
    5 continue
      print 357
  357 format(/' aioxtn: read gen   mom   ppar  pot')
      do  10  ic = 1, nclass
        lgen = .false.
        lmom = .false.
        lpp  = .false.
        lpot = .false.
        rhrmx(1) = 0
        rhrmx(2) = 0
        call r8tos8(dval(w(oclabl),ic),clabl)
        is = w(oics+ic-1)
        call upack('spec lmxa rmt z',sspec,is,lmx,rmax,z,0)
        lstutt = cmdopt('-3s ',4,0,outs) .or.
     .           cmdopt('-4s ',4,0,outs) .or.
     .           cmdopt('-47u',4,0,outs) .or.
     .           cmdopt('-5s ',4,0,outs)

c   ... stuttgart lmtoxx don't lower case, and no extension
        if (lstutt) then
          ifi = fopnx(clabl,2,0,-1)
        else
          ifi = fopna(clabl,-1,0)
        endif
        rewind ifi

        if (cmdopt('-4 ',4,0,outs) .or. cmdopt('-3s ',4,0,outs)) then
          if (cmdopt('-4 ',4,0,outs)) then
          lgen = io4g(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .      qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          lmom = io4mm(alabel,pnl,qnl,idmod,nl,lmx,nsp,rhrmx,vrmax,ifi)
          lpp  = io4pp(alabel,pp,ves,nl,lmx,nsp,ifi)
          else
          lgen = io4sg(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .      qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          lmom = io3smm(alabel,pnl,qnl,idmod,nl,lmx,nsp,rhrmx,vrmax,ifi)
          lpp  = io4spp(alabel,pp,ves,nl,lmx,nsp,ifi)
          endif
          call defdr(ov,   nr*nsp)
          lpot = io4pt(nr,nsp,a,rmax,w(ov),ifi)
          clabl = alabel
        elseif (cmdopt('-5 ',3,0,outs)) then
          call rx('aioxtn not ready for v5')
        elseif (cmdopt('-4s',3,0,outs)) then
          alabel = clabl
          lgen = io4sg(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .      qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          lmom = io4smm(alabel,pnl,qnl,evrl,idmod,nl,lmx,nsp,vrmax,ifi)
          lpp  = io4spp(alabel,pp,ves,nl,lmx,nsp,ifi)
          call defdr(ov,   nr*nsp)
          lpot = io4spt(nr,nsp,a,rmax,w(ov),ifi)
          call fclose(ifi)
          ifi = fopna(clabl,-1,0)
        elseif (cmdopt('-47u',4,0,outs)) then
          nr = 0
          alabel = clabl
c         aiogen same as 5s
          lgen = io5sg(a,alabel,etot,lmx,lrel,nr,nsp,qc,dq,rmax,z,ifi)
          ekin  = etot(2)
          utot  = etot(3)
          sumec = etot(5)
          sumev = etot(6)
          rhoeps = etot(10)
          thrpv = 0
          lmom = io47um(alabel,lmx,nl,nsp,pnl,idmod,qnl,vrmax,ifi)
c     ... for now this doesn't mean anything, so we leave it out.
          lpp = io47up(alabel,lmx,nl,nsp,pp,ves,ifi)
          call defdr(ov,   nr*nsp)
c         aiopot same as 5s
          lpot = io5spt(nr,nsp,a,rmax,w(ov),ifi)
          call fclr(clabl,ifi)
          ifi = fopna(clabl,-1,0)
        elseif (cmdopt('-5s',3,0,outs)) then
          nr = 0
c     ... lgen picks up for us a,nr
          alabel = clabl
          lgen = io5sg(a,alabel,etot,lmx,lrel,nr,nsp,qc,dq,rmax,z,ifi)
          ekin  = etot(2)
          utot  = etot(3)
          sumec = etot(5)
          sumev = etot(6)
          rhoeps = etot(10)
          thrpv = 0
          lmom = io5smm(alabel,idmod,lmx,nl,nsp,pnl,qnl,vrmax,ifi)
          lpp  = .false.
c     ... for now this doesn't mean anything, so we leave it out.
c          lpp  = io5spp(alabel,lmx,nl,nsp,pp,ves,ifi)
          call defdr(ov,   nr*nsp)
          lpot = io5spt(nr,nsp,a,rmax,w(ov),ifi)
          if (lpot .and. lves) then
            call dpscop(w(oves),ves,1,ic,1,1d0)
            call daxpy(nr*nsp,-1d0,ves,0,w(ov),1)
          endif
          call fclr(clabl,ifi)
c          call fclose(ifi)
c          ifi = fopna(clabl,-1,0)
        else
          call rx('aioxtn: need conversion -4, -5 or -3s, -4s or -5s')
        endif
        if (lgen) then
          call dvset(w(oqt),ic,ic,dq)
          call dvset(w(ormax),ic,ic,rmax)
        else
          lqtot = .false.
        endif
        print 358, clabl, lgen,lmom,lpp,lpot
  358   format(2x,a8,5(l6))

        if (lstutt) then
          call fclose(ifi)
          ifi = fopna(clabl,-1,0)
          rewind ifi
        endif
        rewind ifi
        if (lgen) sw = aiogen(clabl,z,rmax,lmx,nsp,lrel,nr,a,qc,
     .    dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,-ifi)
        if (lmom) sw = aiomom(clabl,pnl,qnl,idmod,nl,lmx,nsp,rhrmx,
     .    vrmax,-ifi)
        if (lpp) sw = aiopar(clabl,0,pp,w,ves,nl,lmx,nsp,-ifi)
        if (lpot) sw = aiopot(nr,nsp,a,rmax,-99d0,w(ov),-ifi)
        if (lpot) call rlse(ov)
        call fclr(clabl,ifi)
   10 continue

c     sphere charges are known ... get madelung potential to shift pot
      if (lstutt .and. .not. lves) then
        clabl = 'no'
        if (lqtot) clabl = 'yes'

        print '(/'' ... were all sphere charges found? '',a)', clabl
        if (lqtot) then

        print '('' ... undo madelung shift in site potentials'')'
        call upack('lat alat vol opos',slat,alat,vol,opos,0,0)
        call upack('lat awald nkd nkq odlv oqlv',slat,awald,nkd,nkq,
     .    odlv,oqlv)
        call defdr(omad,nbasp**2)
        call madmat(nbasp,w(opos),awald,alat,vol,w(odlv),nkd,
     .    w(oqlv),nkq,w(omad))
        vmtz(1) = 0
        call madpot(nbasp,1,nclass,w(onrcp),w(oipc),w(oclabl),w(oqt),
     .    0d0,rhrmx,w(ormax),w(omad),w,w,0d0,.false.,w(ovrmax),w(oves),
     .    emad,trumad,vmtz)

        print '(/'' ... re-read files, shifting potentials by ves'')'
        lves = .true.
        goto 5
        endif
      endif


      end

      subroutine aioxt1(mode,nlst,lsts,lstd,nsp,nbasi,ssitei,qnui,
     .  nbas,pnu,qnu)
C- Poke P,Q parameters from source list into destination list
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 import data is for one site only (lsts(1))
Ci         :1 import Pi,Qi(lsts(1..nlst) to Pd,Qd(lstd(1..nlst)
Ci   nlst  :number of sites for which to distribute imported P,Q
Ci   lsts  :list of sites from which to distribute imported P,Q
Ci   lstd  :list of sites to which to distribute imported P,Q
Ci   nbasi :number of sites for which imported data is available
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ssitei:struct for site-specific information; see routine usite
Ci         :It contains pnu for imported data
Ci   qnui  :energy-weighted moments of the sphere charges to import
Co Outputs
Co   pnu   :boundary conditions.  In this routine array is site-based
Co   qnu   :energy-weighted moments of the sphere charges
Co         :In this routine array is site-based
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   28 Feb 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nlst,nsp,lsts(nlst),lstd(nlst),n0,nbasi,nbas
      parameter (n0=10)
      double precision ssitei(1)
      double precision qnui(3,n0,nsp,*)
      double precision pnu(n0,nsp,nbas),qnu(3,n0,nsp,nbas)
C ... Local parameters
      integer ib,ibs,ibd
      double precision ploc(n0,2)

      print '(''      source  dest'')'
      do  ib = 1, nlst
        if (mode .eq. 0) then
          ibs = lsts(1)
        else
          ibs = lsts(ib)
        endif
        ibd = lstd(ib)
        print 333, ib, ibs,ibd
  333   format(i4,2i6,$)
        if (ibs .gt. nbasi) then
          print 334, 'ib(s)'
  334     format(' *** ', a,' exceeds max val ... copy suppressed')
        elseif (ibd .gt. nbas) then
          print 334, 'dest'
        else
          print *, ' '
          call upack2('site pnu',ssitei,ibs,ploc)
          call dcopy(n0*nsp,ploc,1,pnu(1,1,ibd),1)
          call dcopy(3*n0*nsp,qnui(1,1,1,ibs),1,qnu(1,1,1,ibd),1)
        endif
      enddo
      end

c ---- file i/o atomic data, v4 ---
      logical function io4g(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     . qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
c- file i/o for atomic general data.  ifi>0 for read, <0 for write
c ----------------------------------------------------------------
ci inputs
ci
co outputs
co
cr remarks
cr
c ----------------------------------------------------------------
C     implicit none
c passed parameters
      integer lmx,nsp,nr,ifi
      double precision z,rmax,a,
     . qc,dq,vrmax(2),sumec,sumev,thrpv,ekin,utot,rhoeps,etot
      character*4 alabel
      logical lrel
c local parameters
      logical scat

      io4g = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'gen:',':',.true.)) return
        read(ifi,200) alabel,lmx,nsp,rmax,lrel,nr,a
        read(ifi,201) z,qc,dq,vrmax,rhoeps
        read(ifi,202) sumec,utot,ekin,sumev,etot,thrpv
        io4g = .true.
      else
        write(-ifi,'(''gen:'')')
        write(-ifi,100) alabel,lmx,nsp,rmax,lrel,nr,a
        write(-ifi,101) int(z),int(qc),dq,vrmax,rhoeps
        write(-ifi,102) sumec,utot,ekin,sumev,etot,thrpv
      endif
  100 format(3x,a4,'  lmx=',i1,'  nspin=',i1,'  rmax=',f9.6,'  rel=',l1,
     .       '  nr=',i4,'  a=',f5.3)
  101 format(3x,'z=',i2,  '  qc=',i2,'  qtot=',f9.6,'  vrmax=',2f9.6,
     .          '  rhoeps=',f12.6)
  102 format('   sumec=',f15.7,'  utot=',f15.7,'  ekin= ',f15.7/
     .       '   sumev=',f15.7,'  etot=',f15.7,'  thrpv=',f15.7)

  200 format(3x,a4,6x,i1,8x,i1,7x,f9.6,6x,l1,5x,i4,4x,f4.3)
  201 format(3x,2x,f2.0,5x,f2.0,7x,f9.6,8x,2f9.6,9x,f12.6)
  202 format(2(9x,f15.7,7x,f15.7,8x,f15.7/))

      end

      logical function io3smm(albl,pl,ql,idmod,nl,lmax,nsp,rh,vrmax,ifi)
c- file i/o for moments.  ifi>0 for read, <0 for write
c ----------------------------------------------------------------
ci inputs
ci   albl,nl,lmax,nsp
ci   ifi:  logical unit: positive for read, negative for write
ci   pl:   related to log derivatives (see atom)
ci   ql:   moments q (see atom) (ifi > 0)
co   rh,vrmax density and vxc at rmax (ifi>0)
co outputs
co   pl,ql:   moments q (see atom) (ifi < 0)
co   aiomom:true unless read error or category not found
co   rh,vrmax density and vxc at rmax (ifi<0)
cr remarks
cr   input for second spin defaults to that for first spin if not
cr   available for second spin.
cr   bug in reading vrmax for nsp=2
c ----------------------------------------------------------------
C     implicit none
c passed parameters
      integer nl,lmax,nsp,ifi
      character*4 albl
      integer idmod(0:nl-1)
      double precision pl(0:nl-1,nsp),ql(3,0:nl-1,nsp),rh,vrmax(2)
c local parameters
      double precision pi,dl,xx(4)
      integer i,l,ii,j,k,iprint
      logical scat,rdstrn
#if AIX
      character*72 s
#endif
      io3smm = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        io3smm = .true.
        backspace ifi
        rh = 0
        vrmax(2) = 0
        read(ifi,'(21x,f8.6)') vrmax(1)
        read(ifi,*)
        do  10  i = 1, nsp
        do  16  l = 0, lmax
          read(ifi,100,err=12) ii, (xx(j), j=1,4), k
          if (ii .ne. l) call rx('IO3SMM: bad l quantum number')
          idmod(l) = k
          pl(l,i) = xx(1)
          pl(l,nsp) = xx(1)
          do  15  ii = 1, 3
            ql(ii,l,nsp) = xx(1+ii)/nsp
            ql(ii,l,i) = xx(1+ii)
   15     continue
   16   continue
C --- Patch for bug in AIX err= ---
#if AIX
        if (i .eq. 1 .and. nsp .eq. 2) then
          if (.not. rdstrn(ifi,s,72,.false.)) goto 18
          if (s(1:5) .ne. '   0 ') goto 18
          backspace ifi
        endif
#endif
   10   continue
      else
        pi = 4*datan(1.d0)
        write(-ifi,'(''MOMNTS:  '',a4,''  rho,vrmax='',3f10.6)')
     .    albl, rh, vrmax
        write(-ifi,891)
  891   format(
     .  '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          dl = dtan(pi*(.5d0 - pl(l,i)))
          if (dabs(dl) .gt. 9999) dl = 0
          write(-ifi,100) l,pl(l,i),(ql(ii,l,i),ii=1,3),idmod(l),dl
   11   continue
      endif
      return

  100 format(i2,4f11.7,i3,f13.7,f11.7)

C --- handle read exception ---
#if AIX
   18 l = 0
      i = 2
#endif
   12 continue
      if (l .eq. 0 .and. i .eq. 2) then
        if (iprint() .ge. 20)
     .  print *, 'IO3SMM, ATOM=',albl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  40  l = 0, lmax
          do  40  ii = 1, 3
          ql(ii,l,1) = ql(ii,l,2)
   40   continue
      else
        io3smm = .false.
      endif
      end

      logical function io4mm(albl,pl,ql,idmod,nl,lmax,nsp,rh,vrmax,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs
Ci   albl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   pl:   related to log derivatives (see atom)
Ci   ql:   moments q (see atom) (ifi > 0)
Co   rh,vrmax density and Vxc at rmax (ifi>0)
Co Outputs
Co   pl,ql:   moments q (see atom) (ifi < 0)
Co   aiomom:true unless read error or category not found
Co   rh,vrmax density and Vxc at rmax (ifi<0)
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   Bug in reading vrmax for nsp=2
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,lmax,nsp,ifi
      character*4 albl
      integer idmod(0:nl-1)
      double precision pl(0:nl-1,nsp),ql(3,0:nl-1,nsp),rh,vrmax(2)
C Local parameters
      double precision pi,dl,xx(4)
      integer i,l,ii,j,k,iprint
      logical scat,rdstrn
#if AIX
      character*72 s
#endif

      io4mm = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        io4mm = .true.
        backspace ifi
        read(ifi,'(25x,3f10.6)') rh, vrmax
        read(ifi,*)
        do  10  i = 1, nsp
        do  16  l = 0, lmax
          read(ifi,100,err=12) ii, (xx(j), j=1,4), k
          if (ii .ne. l) call rx('IO4MM: bad l quantum number')
          idmod(l) = k
          pl(l,i) = xx(1)
          pl(l,nsp) = xx(1)
          do  15  ii = 1, 3
            ql(ii,l,nsp) = xx(1+ii)/nsp
            ql(ii,l,i) = xx(1+ii)
   15     continue
   16   continue
C --- Patch for bug in AIX err= ---
#if AIX
        if (i .eq. 1 .and. nsp .eq. 2) then
          if (.not. rdstrn(ifi,s,72,.false.)) goto 18
          if (s(1:5) .ne. '   0 ') goto 18
          backspace ifi
        endif
#endif
   10   continue
      else
        pi = 4*datan(1.d0)
        write(-ifi,'(''MOMNTS:  '',a4,''  rho,vrmax='',3f10.6)')
     .    albl, rh, vrmax
        write(-ifi,891)
  891   format(
     .  '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          dl = dtan(pi*(.5d0 - pl(l,i)))
          if (dabs(dl) .gt. 9999) dl = 0
          write(-ifi,100) l,pl(l,i),(ql(ii,l,i),ii=1,3),idmod(l),dl
   11   continue
      endif
      return

  100 format(i4,4f13.7,i4,f13.7)

C --- handle read exception ---
#if AIX
   18 l = 0
      i = 2
#endif
   12 continue
      if (l .eq. 0 .and. i .eq. 2) then
        if (iprint() .ge. 20)
     .  print *, 'IO4MM, ATOM=',albl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  40  l = 0, lmax
          do  40  ii = 1, 3
          ql(ii,l,1) = ql(ii,l,2)
   40   continue
      else
        io4mm = .false.
      endif
      end

      logical function io4pp(albl,pp,ves,nl,lmax,nsp,ifi)
C- File I/O for potential parameters.
C ----------------------------------------------------------------
Ci Inputs
Ci   albl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   ves:  electrostatic potential at rmax
Ci   pp:   potential parameters (ifi > 0) read from file ifi
Co Outputs
Co   pp:   potential parameters (ifi < 0) written to file -ifi
Co   aiopar:true unless read error or category not found
Cr Remarks
Cr   11 Apr 94 Added convention: reading 1/p=0 => set p to 0
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*4 albl
      integer nl,lmax,nsp,ifi
      double precision pp(6,0:nl-1,nsp),ves
C Local parameters
      double precision xx,x2(6)
      equivalence (xx,x2)
      integer i,l,ll,k,iprint,i1mach
      logical scat,sw
#if AIX
      logical rdstrn
      character*72 s
#endif

      sw = .false.
      io4pp = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        io4pp = .true.
        backspace ifi
        read(ifi,'(17x,f12.8)') ves
        read(ifi,*)
        do  10  i = 1, nsp
        do  16  l = 0, lmax
          read(ifi,20,err=12) ll, (x2(k), k=1,6)
          if (ll .ne. l) call rx('IO4PP: bad l quantum number')
C Map delta into sqrt(delta), preserving sign and 1/sqrt(p) into p
          x2(3) = dsign(1.d0,x2(3))*dsqrt(dabs(x2(3)))
          if (x2(4) .eq. 0) then
            if (iprint() .gt. 40)
     .        print *, 'io4pp: encountered 1/p=0 ... set p to 0'
          else
            x2(4) = 1/x2(4)**2
          endif
          do  15  k = 1, 6
          pp(k,l,i) = x2(k)
   15     pp(k,l,nsp) = x2(k)
   16     continue
C --- Patch for bug in AIX err= ---
#if AIX
        if (i .eq. 1 .and. nsp .eq. 2) then
          if (.not. rdstrn(ifi,s,72,.false.)) goto 18
          if (s(1:3) .ne. ' 0 ') goto 18
          backspace ifi
        endif
#endif
   10   continue
      else
        write(-ifi,21) albl, ves
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          xx = pp(3,l,i)
          if (pp(4,l,i) .eq. 0) then
            write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .        0d0, (pp(k,l,i), k=5,6)
            sw = .true.
          else
            write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .        1/dsqrt(pp(4,l,i)), (pp(k,l,i), k=5,6)
        endif
   11   continue
        if (iprint() .ge. 30  .and. sw)
     .    write(i1mach(2),*) 'io4pp: encountered p=0 ... wrote 1/p=0'
        io4pp = .true.
      endif
      return

C --- handle read exception ---
#if AIX
   18 l = 0
      i = 2
#endif
   12 continue
      if (l .eq. 0 .and. i .eq. 2) then
        if (iprint() .ge. 20)
     .  print *, 'IO4PP, ATOM=',albl,
     .  ':  spin 2 input missing; taken from spin 1 data'
      else
        io4pp = .false.
      endif

   20 format(i2,3f12.8,f12.7,3f12.8)
   21 format('PPAR:  ',a4,'  ves=',f12.8/ ' l',5x,'e_nu',10x,'C',8x,
     .       '+/-del',5x,'1/sqrt(p)',6x,'gam',9x,'alp')

      end

      logical function io4pt(nr,nsp,a,rmax,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   v, if file write
Co Outputs
Co   v, if file read
Cr Remarks
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, aiopot returns true only if the category is found,
Cr    the file's value of a and nr match input and rmax is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,v(nr,nsp)
C Local parameters
      integer i,isp,nr2,nsp2
      double precision a2,rmax2
      logical scat

      io4pt = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'POT:',':',.true.)) return
        read(ifi,102,err=15) nr2,nsp2,a2,rmax2
        if (nr .eq. 0) nr=nr2
        if (nsp2 .eq. 0) nsp=nsp2
        if (a2 .eq. 0) a=a2
        if (rmax .eq. 0) rmax=rmax2
        if (a2 .ne. a .or. nr .ne. nr2 .or.
     .      dabs(rmax2-rmax) .gt. .03) goto 15
        do  10  isp = 1, min0(nsp2,nsp)
c          read(ifi,101,err=15) (v(i,isp),i = 1,nr)
          read(ifi,101) (v(i,isp),i = 1,nr)
          do  10  i = 1, nr
          v(i,nsp) = v(i,isp)
   10   continue
        io4pt = .true.
   15   continue
      else
        write(-ifi,'(''POT:'')')
        write(-ifi,102) nr,nsp,a,rmax
        do  11  isp = 1, nsp
          write(-ifi,101) (v(i,isp),i = 1,nr)
   11   continue
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5)
      return
      end

C ---- File I/O atomic data, Stuttgart v4 ---
      logical function io4sg(alabel,z,rmax,lmx,nsp,lrel,nr,a,
     .  qc,dq,vrmax,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
C     implicit none
C Passed parameters
      integer lmx,nsp,nr,ifi
      double precision z,rmax,a,
     . qc,dq,vrmax(2),sumec,sumev,thrpv,ekin,utot,rhoeps,etot
      character*8 alabel
      logical lrel
C Local parameters
      logical scat

      io4sg = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'GEN:',':',.true.)) return
        read(ifi,200) alabel,lmx,nsp,rmax,lrel,nr,a
        read(ifi,201) z,qc,dq,vrmax,rhoeps
        read(ifi,202) sumec,utot,ekin,sumev,etot,thrpv
        io4sg = .true.
      else
        write(-ifi,'(''GEN:'')')
        write(-ifi,100) alabel,lmx,nsp,rmax,lrel,nr,a
        write(-ifi,101) int(z),int(qc),dq,vrmax,rhoeps
        write(-ifi,102) sumec,utot,ekin,sumev,etot,thrpv
      endif
  100 format(3X,A4,'  LMX=',I1,'  NSPIN=',I1,'  RMAX=',F9.6,'  REL=',L1,
     .       '  NR=',I3,'  A=',F4.3)
  101 format(3X,'Z=',I2,  '  QC=',I2,'  QTOT=',F9.6,'  VRMAX=',2F9.6,
     .          '  RHOEPS=',F12.6)
  102 format(
     .  '   SUMEC=',F15.6,'  UTOT=',F15.6,'  EKIN= ',F15.6/
     .  '   SUMEV=',F15.6,'  ETOT=',F15.6,'  THRPV=',F15.6)

  200 format(3X,A4,6X,I1,8X,I1,7X,F9.6,6X,L1,5X,I3,4X,F4.3)
  201 format(3X,2X,F2.0,5X,F2.0,7X,F9.6,8X,2F9.6,9X,F12.6)
  202 format(2(9X,F15.6,7X,F15.6,8X,F15.6/))

      end
      logical function io4smm(alabl,pl,ql,evrl,idmod,nl,lmax,nsp,vrmax,
     .        ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   pl:   related to log derivatives (see atom)
Ci   ql:   moments q (see atom) (ifi > 0)
Ci   evrl: eshift (if idmod.ne.0, ifi > 0)
Co Outputs
Co   pl,ql:   moments q (see atom) (ifi < 0)
Ci   evrl: eshift (if idmod.ne.0, ifi < 0)
Co   io4smm:true unless read error or category not found
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,lmax,nsp,ifi
      character*8 alabl
      integer idmod(0:nl-1)
      double precision pl(0:nl-1,nsp),ql(3,0:nl-1,nsp),evrl(0:nl-1,nsp),
     .                 vrmax
C Local parameters
      double precision pi,dl,xx(6)
      integer i,l,ii,j,k,iprint
      logical scat
      external scat,iprint

      io4smm = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        io4smm = .true.
        backspace ifi
        read(ifi,'(21x,f10.6)') vrmax
        read(ifi,*)
        do  10  i = 1, nsp
          do  10  l = 0, lmax
          read(ifi,100,err=12) ii, (xx(j), j=1,4), k, (xx(j), j=5,6)
          if (ii .ne. l) stop 'IO4SMM: bad l quantum number'
          idmod(l) = k
          evrl(l,i) = xx(6)
          evrl(l,nsp) = xx(6)
          pl(l,i) = xx(1)
          pl(l,nsp) = xx(1)
          do  15  ii = 1, 3
            ql(ii,l,nsp) = xx(1+ii)/nsp
            ql(ii,l,i) = xx(1+ii)
   15     continue
   10   continue
      else
        pi = 4*datan(1.d0)
        write(-ifi,'(''MOMNTS:  '',A4,''  VRMAX='',F10.6)') alabl, vrmax
        write(-ifi,891)
  891   format(
     .  ' L',6X,'PL',9X,'Q0',9X,'Q1',9X,'Q2',4X,' ID ',5X,'DL',9X,'ESH')
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          dl = dtan(pi*(.5d0 - pl(l,i)))
          write(-ifi,100) l,pl(l,i),(ql(ii,l,i),ii=1,3),idmod(l),dl,
     .                    evrl(l,i)
   11   continue
      endif
      return

  100 format(i2,4f11.7,i3,f12.7,f11.7)

C --- handle read exception ---
   12 continue
      if (l .eq. 0 .and. i .eq. 2) then
        if (iprint() .ge. 20)
     .  print *, 'IO4SMM, ATOM=',alabl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  40  l = 0, lmax
          do  40  ii = 1, 3
          ql(ii,l,1) = ql(ii,l,2)
   40   continue
      else
        io4smm = .false.
      endif
      end
      logical function io4spp(alabl,pp,ves,nl,lmax,nsp,ifi)
C- File I/O for potential parameters.
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl,nl,LMAX,NSP
Ci   IFI:  logical unit: positive for read, negative for write
Ci   ves:  electrostatic potential at rmax
Ci   pp:   potential parameters (ifi > 0) read from file ifi
Co Outputs
Co   pp:   potential parameters (ifi < 0) written to file -ifi
Co   io4spp:true unless read error or category not found
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmax,nsp,ifi
      double precision pp(6,0:nl-1,nsp),ves
C Local parameters
      double precision xx,x2(6)
      equivalence (xx,x2)
      integer i,l,ll,k,iprint
      logical scat
      external scat,iprint

      io4spp = .false.

      if (ifi .gt. 0) then
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        io4spp = .true.
        backspace ifi
        read(ifi,'(17x,f10.6)') ves
        read(ifi,*)
        do  10  i = 1, nsp
          do  10  l = 0, lmax
          read(ifi,20,err=12) ll, (x2(k), k=1,6)
          if (ll .ne. l) stop 'IO4SPP: bad l quantum number'
C Map delta into sqrt(delta), preserving sign and 1/sqrt(p) into p
#ifndef NOMAP_PP
          x2(3) = dsign(1.d0,x2(3))*dsqrt(dabs(x2(3)))
          x2(4) = 1/x2(4)**2
#endif
          do  15  k = 1, 6
          pp(k,l,i) = x2(k)
   15     pp(k,l,nsp) = x2(k)
   10   continue
      else
        write(-ifi,21) alabl, ves
        do  11  i = 1, nsp
          do  11  l = 0, lmax
          xx = pp(3,l,i)
c aek avoid 1/sqrt(p) to be infinity
          if ( pp(4,l,i).lt.1.d-99 ) then
          write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .                   0.0d0, (pp(k,l,i), k=5,6)
          else
          write(-ifi,20) l, (pp(k,l,i), k=1,2), xx**2*dsign(1.d0,xx),
     .                   1/dsqrt(pp(4,l,i)), (pp(k,l,i), k=5,6)
          endif
   11   continue
      endif
      return

C --- handle read exception ---
   12 continue
      if (l .eq. 0 .and. i .eq. 2) then
        if (iprint() .ge. 20)
     .  print *, 'IO4SPP, ATOM=',alabl,
     .  ':  spin 2 input missing; taken from spin 1 data'
      else
        io4spp = .false.
      endif

   20 format(i2,6f12.7)
   21 format('PPAR:  ',A4,'  VES=',F10.6/ ' L',8X,'ENU',9X,'C',8X,
     .       '+/-DEL',4X,'1/SQRT(P)',6X,'GAM',9X,'ALP')

      end
      logical function io4spt(nr,nsp,a,rmax,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: file logical unit, but >0 for read, <0 for write
Ci   nr,nsp,a,rmax
Ci   v, if file write
Co Outputs
Co   v, if file read
Cr Remarks
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, rmax,
Cr    followed by the potential.
Cr    On reading, io4spt returns true only if the category is found,
Cr    the file's value of a and nr match input and rmax is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ifi,nr,nsp
      double precision a,rmax,v(nr,nsp)
C Local parameters
      integer i,isp,nr2,nsp2
      double precision a2,rmax2
      logical scat

      io4spt = .false.
      if (ifi .gt. 0) then
        if (.not. scat(ifi,'POT:',':',.true.)) return
        read(ifi,102,err=15) nr2,nsp2,a2,rmax2
        if ( dabs(a2-a) .gt. 1.d-4 .or. nr .ne. nr2 .or.
     .      dabs(rmax2-rmax) .gt. .03) goto 15
        do  10  isp = 1, min0(nsp2,nsp)
c          read(ifi,101,err=15) (v(i,isp),i = 1,nr)
          read(ifi,101) (v(i,isp),i = 1,nr)
          do  10  i = 1, nr
          v(i,nsp) = v(i,isp)
   10   continue
        io4spt = .true.
   15   continue
      else
        write(-ifi,'(''POT:'')')
        write(-ifi,102) nr,nsp,a,rmax
        do  11  isp = 1, nsp
          write(-ifi,101) (v(i,isp),i = 1,nr)
   11   continue
      endif
  101 format(1p,5d16.9)
  102 format(2i5,2f12.5)
      return
      end

C ---- File I/O atomic data, Stuttgart v47+u ---
C ... aiogen is same as 5g
      logical function io47um(clabl,lmx,nl,nsp,pnl,idmod,qnl,vrmax,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------------
Ci Inputs:
Ci   clabl :name of the different inequivalent atom
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Co Outputs ( Inputs if file write: ifi < 0)
Cio  pnl   :principal quantum number
Cio  qnl   :moments
Cio  vrmax :potential at rmax=wsr
Cr Remarks:
Cr   aiomom:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   zeroth moment in atomic file is taken by convention to be
Cr   electronic charge q0 = m0 + p^gamma m2,  rather than m0.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nl,lmx,nsp,ifi
C-ise
      integer idmod(0:nl-1)
C-ise
      double precision pnl(0:nl-1,*),qnl(0:2,0:nl-1,*),vrmax
      character*4 clabl
C Local variables:
      integer i,iost,isp,l,l2,j,k,idm
      double precision pi,dl,xx(6)
      logical scat
      character chr
      character*72 messg
      parameter(pi=3.14159265358979324d0)
C External calls:
      external dcopy,dscal,errmsg,scat
C Intrinsic functions:
      intrinsic dtan

      io47um = .false.
      if (ifi.gt.0) then
C ----  Read moments and p's from atomic file
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        backspace ifi
        read(ifi,200,err=15,end=15,iostat=iost)vrmax
        read(ifi,*  ,err=15,end=15,iostat=iost)
        do isp = 1, 2
          do l = 0, lmx
           read(ifi,201,err=15,end=15,iostat=iost)chr
           if (chr.ge.'0'.and.chr.le.'9') then
             backspace ifi
             read(ifi,301,err=15,end=15,iostat=iost)
     .            l2,(xx(j),j=1,4),idm
             if (isp.le.nsp) then
               idmod(l) = idm
               pnl(l,isp) = xx(1)
               call dcopy(3,xx(2),1,qnl(0,l,isp),1)
             elseif (l.eq.l2) then
               pnl(l,1) = (pnl(l,1)+xx(1))*0.5d0
               do k=0,2
                 qnl(k,l,1)= qnl(k,l,1)+xx(k+2)
               enddo
             endif
           else
             if (isp.le.nsp) then
               write(messg,400)clabl,isp,l
               call errmsg(messg,1)
               goto 15
             else
               io47um = .true.
               return
             endif
           endif
          enddo
        enddo
        io47um = .true.
      else
C ----  write atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,vrmax
        do isp = 1, nsp
          do l = 0, lmx
            dl = dtan(pi*(0.5d0 - pnl(l,isp)))
            write(-ifi,301,err=15,iostat=iost)
C-ise     .       l,pnl(l,isp),(qnl(k,l,isp),k=0,2),0,dl
     .       l,pnl(l,isp),(qnl(k,l,isp),k=0,2),idmod(l),dl
          enddo
        enddo
      endif
      return

   15 continue
      if (iost.ne.0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi.gt.0 .and. l.eq.0 .and. isp.eq.2) then
        write(messg,401) clabl
        call errmsg(messg,1)
        call dcopy(  nl,pnl (0,1)  ,1,pnl (0,2),1)
        call dcopy(3*nl,qnl (0,0,1),1,qnl (0,0,2),1)
        call dscal(6*nl,0.5d0,qnl,1)
        io47um = .true.
      endif

200   format(21x,f10.6)
201   format(1x,a1)
300   format('MOMNTS:  ',a4,'  VRMAX=',f10.6,
     .  /' L',6x,'PL',9x,'Q0',9x,'Q1',9x,'Q2   IDMOD',5x,'DL')
301   format(i2,4f11.7,i3,f12.7,f11.7)
400   format(' IO47UM: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing moments.$')
401   format(' IO47UM: ATOM ',a4,':spin 2 input missing;',
     .       ' moments split.$')
402   format(' IO47UM: I/O error, iostat=',i3,'$')
      end

      logical function io47up(clabl,lmx,nl,nsp,pp,ves,ifi)
C- File I/O for moments.  ifi > 0 for read, < 0 for write
C ----------------------------------------------------------------------
Ci Inputs: (pp: if file write : ifi < 0)
Ci   ifi   :logical unit: positive for read, negative for write
Ci   clabl :name of the different inequivalent atom
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Cio Inputs/Outputs:
Cio  pp    :potential parameters
Cio  ves   :electrostatic potential at wsr
Cr Remarks:
Cr   aiopar:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer npp
      parameter(npp=6)
      integer nl,lmx,nsp,ifi
      double precision pp(npp,0:nl-1,nsp),ves
      character*4 clabl
C Local variables:
      integer i,isp,l,l2,k,iost
      double precision x2(6),tiny,vmax
      logical scat
      parameter(tiny=1.d-6,vmax=.999999d9)
      character  chr
      character*4 cdum
      character*80 messg
C External calls:
      external errmsg,scat
C Intrinsic functions:
      intrinsic dabs,dsign,dsqrt

      io47up = .false.

      if (ifi.gt.0) then
C ----- Read potential parameters from atomic file
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        backspace (ifi)
        read(ifi,500,err=15,end=15,iostat=iost)cdum,ves
        do isp = 1, nsp
          do l = 0, lmx
            read(ifi,200,err=15,end=15,iostat=iost)chr
            if ( chr.ge.'0'.and. chr.le.'9') then
              backspace ifi
              read(ifi,*  ,err=15,end=15,iostat=iost)l2,(x2(k),k=1,6)
            else
              write(messg,401)clabl,isp,l
              call errmsg(messg,1)
              goto 15
            endif
C --------- Map delta into sqrt(delta) preserving sign, 1/sqrt(p) into p
            x2(3) = dsqrt(dabs(x2(3)))*dsign(1.d0,x2(3))
C-ise            x2(4) = dmax1(x2(4),tiny)
C-ise            x2(4) = 1/x2(4)/x2(4)
            call dcopy(6,x2,1,pp(1,l,isp),1)
          enddo
        enddo
        io47up = .true.
      else
C ----- Write potential parameters to atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,ves
        do isp = 1, nsp
          do l = 0, lmx
            do k=1,6
              x2(k) = pp(k,l,isp)
            enddo
C-ise            if (dabs(x2(4)).le.1.d0/vmax/vmax)  x2(4) = 1.d0/vmax/vmax
            x2(3) = x2(3)*x2(3)*dsign(1.d0,x2(3))
C-ise            x2(4) = 1/dsqrt(x2(4))
            write(-ifi,301,err=15,iostat=iost)l
            do k=1,6
              if(x2(k).gt.-9.9999999d0.and.x2(k).lt.99.9999999d0) then
                write(-ifi,302,err=15,iostat=iost)x2(k)
              else
                write(-ifi,303,err=15,iostat=iost)x2(k)
              endif
            enddo
            write(-ifi,304,err=15,iostat=iost)
          enddo
        enddo
      endif
      return

15    continue
      if (iost.ne.0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi.gt.0 .and. l.eq.0 .and. isp.eq.2) then
        write(messg,400) clabl
        call errmsg(messg,1)
        call dcopy(6*(lmx+1),pp(1,0,1),1,pp(1,0,2),1)
        io47up = .true.
      endif

200   format(1x,a1)
300   format('PPAR:  ',a4,'  VES=',f10.6/ ' L',5x,'ENU',9x,'C',8x,
     .       '+/-DEL',7x,'P',8x,'GAMMA',6x,'ALPHA')
301   format(1x,i1,$)
302   format(f11.7,$)
303   format(e11.5e1,$)
304   format()
400   format(' IO47UP: ATOM ',a4,':spin 2 input missing; ',
     .          'taken from spin 1.$')
401   format(' IO47UP: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing potential parameters$')
402   format(' IO47UP: I/O error, iostat=',i3,'$')
500   format(7x,a4,6x,f10.6/,1x)
      end

C ---- File I/O atomic data, Stuttgart v56 ---
      logical function io5sg(a,clabl,etot,lmx,lrel,nr,nsp,qc,qtot,
     .                       wsr,z,ifi)
C- File I/O for atomic general data.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci Inputs/Outputs:
Cio  clabl :name of the different inequivalent atom
Cio  etot  :total energy of the atom
Cio  lmx   :maximum l-quantum number inside the sphere
Cio  lrel  :true, for relativistic calculation, false otherwise
Cio  nr    :number of mesh points
Cio  nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Cio  qc    :core electronic charge
Cio  qtot  :difference between spin up and spin down charge
Cio  wsr   :Wigner-Seitz radius
Cio  z     :nuclear charge
C ----------------------------------------------------------------
C     implicit none
C Passed variables:
      integer lmx,nsp,nr,ifi
      double precision a,etot(13,0:2),qc,qtot,wsr,z
      character*4 clabl
      logical lrel
C Local variables:
      integer i,iqc,iost,isp,iz,iz2,lmx2,lgunit,nr2,nsp2
      double precision a2,dsum,dum,wsr2,tiny
      logical scat
      character*500 messg
      parameter(tiny=1.d-5)
C External calls:
      external  dsum,errmsg,scat,lgunit
C Intrinsic functions:
      intrinsic idnint

      io5sg = .false.

      iz=idnint(z)

      if (ifi .gt. 0) then
        if (.not. scat(ifi,'GENERAL:',':',.true.)) return
        backspace ifi
        read(ifi,300,err=15,end=15,iostat=iost)lmx2,nsp2,wsr2,
     .                                         lrel,nr2,a2
        read(ifi,301,err=15,end=15,iostat=iost)iz2,iqc,qtot
        qc = iqc
        if (nsp.eq.2) then
          do  10  isp = 1, nsp2
            read(ifi,302,err=15,end=15,iostat=iost)i
            read(ifi,303,err=15,end=15,iostat=iost)(etot(i,isp),i=1,13)
   10     continue
        endif
        read(ifi,303,err=15,end=15,iostat=iost)(etot(i,0),i=1,13)
C        if (dabs(a-a2).gt.tiny .or.
C     .      lmx.ne.lmx2 .or. nr.ne.nr2 .or. nsp.ne.nsp2 .or.
C     .      dabs(wsr-wsr2).gt.tiny .or. iz.ne.iz2) then
C          write(messg,400)clabl,a,a2,lmx,lmx2,nr,nr2,nsp,
C     .                    nsp2,wsr,wsr2,iz,iz2
C          call errmsg(messg,-1)
C   ... Require match for these only wsr,z only
        if (lmx.ne.lmx2 .or. nsp.ne.nsp2 .or. dabs(wsr-wsr2).gt.tiny
     .      .or. iz.ne.iz2) then
          print *, 'io5sg: file mismatch'
          io5sg = .false.
          return
        else
          io5sg = .true.
        endif
        a = a2
        nr = nr2
      else
       if (-ifi.ne.lgunit(1)) then
         write(-ifi,300,err=15,iostat=iost)lmx,nsp,wsr,lrel,nr,a
         write(-ifi,301,err=15,iostat=iost)iz,idnint(qc),qtot
       endif
       if (nsp.eq.2) then
         do  20  isp = 1, nsp
           write(-ifi,302,iostat=iost)isp
           write(-ifi,303,iostat=iost)(etot(i,isp),i=1,13)
   20    continue
       endif
       write(-ifi,303,iostat=iost)(etot(i,0),i=1,13)
       io5sg = .true.
      endif
      return

C --- handle read exception ---
15    continue
      write(messg,401)iost
      call errmsg(messg,1)

300   format('GENERAL: LMX=',i1,'  NSPIN=',i1,'   WSR=',f9.6,'  REL=',
     .       l1,'  NR=',i3,'  A=',f4.3)
301   format(3x,'Z=',i2,  '  QC=',i2,'  QTOT=',f9.6)
302   format(/' SPIN',i2,':')
303   format(/' ETOT  =',f15.7,2x,'EKIN = ',f15.7,2x,'UTOT =',f15.7,
     .       /' EMAD  =',f15.7,2x,'SUMEC= ',f15.7,2x,'SUMEV=',f15.7,
     .       /' RHOV  =',f15.7,2x,'UNUCL= ',f15.7,2x,'UELEC=',f15.7,
     .       /' RHOEPS=',f15.7,2x,'REPSNL=',f15.7,2x,'RHOMU=',f15.7,
     .       /' RMUNL =',f15.7,2x)

400   format(' AIOGEN: difference between data in atomic file ',a4,
     .       '|        and actual values:',
     .       '|        A    ',': actual value:',f8.4 ,'   read:',f8.4,
     .       '|        LMX  ',': actual value:',i8   ,'   read:',i8,
     .       '|        NR   ',': actual value:',i8   ,'   read:',i8,
     .       '|        NSP  ',': actual value:',i8   ,'   read:',i8,
     .       '|        WSR  ',': actual value:',f8.4 ,'   read:',f8.4,
     .       '|        Z    ',': actual value:',i8   ,'   read:',i8,'$')
401   format(' AIOGEN: I/O error, iostat=',i3,'$')
      end
      logical function io5smm(clabl,idmenu,lmx,nl,nsp,pnl,qnl,vrmax,ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------------
Ci Inputs:
Ci   clabl :name of the different inequivalent atom
Ci   idmenu:if lmul=t and idmenu=1 => orbital collected to a single enu
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Co Outputs ( Inputs if file write: ifi < 0)
Cio  pnl   :principal quantum number
Cio  qnl   :moments
Cio  vrmax :potential at rmax=wsr
Cr Remarks:
Cr   aiomom:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   zeroth moment in atomic file is taken by convention to be
Cr   electronic charge q0 = m0 + p^gamma m2,  rather than m0.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer idmenu(0:*),nl,lmx,nsp,ifi
      double precision pnl(0:nl-1,*),qnl(0:2,0:nl-1,*),vrmax
      character*4 clabl
C Local variables:
      integer i,iost,isp,l,l2,j,k,ii
      double precision pi,dl,xx(6)
      logical scat
      character chr
      character*72 messg
      parameter(pi=3.14159265358979324d0)
C External calls:
      external dcopy,dscal,errmsg,scat
C Intrinsic functions:
      intrinsic dtan

      io5smm = .false.
      if (ifi.gt.0) then
C ----  Read moments and p's from atomic file
        if (.not. scat(ifi,'MOMNTS:',':',.true.)) return
        backspace ifi
        read(ifi,200,err=15,end=15,iostat=iost)vrmax
        read(ifi,*  ,err=15,end=15,iostat=iost)
        do  10  isp = 1, 2
          do  12  l = 0, lmx
           read(ifi,201,err=15,end=15,iostat=iost)chr
           if (chr.ge.'0'.and.chr.le.'9') then
             backspace ifi
             read(ifi,301,err=15,end=15,iostat=iost)
     .            l2,(xx(j),j=1,4),ii
             if (isp.le.nsp) then
               pnl(l,isp) = xx(1)
               call dcopy(3,xx(2),1,qnl(0,l,isp),1)
               idmenu(l) = ii
             elseif (l.eq.l2) then
               pnl(l,1) = (pnl(l,1)+xx(1))*0.5d0
               do  14  k=0,2
                 qnl(k,l,1)= qnl(k,l,1)+xx(k+2)
   14          continue
               idmenu(l) = ii
             endif
           else
             if (isp.le.nsp) then
               write(messg,400)clabl,isp,l
               call errmsg(messg,1)
               goto 15
             else
               io5smm = .true.
               return
             endif
           endif
   12    continue
   10   continue
        io5smm = .true.
      else
C ----  write atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,vrmax
        do  20  isp = 1, nsp
          do  22  l = 0, lmx
            dl = dtan(pi*(0.5d0 - pnl(l,isp)))
            write(-ifi,301,err=15,iostat=iost)
     .       l,pnl(l,isp),(qnl(k,l,isp),k=0,2),idmenu(l),dl
   22     continue
   20   continue
      endif
      return

   15 continue
      if (iost.ne.0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi.gt.0 .and. l.eq.0 .and. isp.eq.2) then
        write(messg,401) clabl
        call errmsg(messg,1)
        call dcopy(  nl,pnl (0,1)  ,1,pnl (0,2),1)
        call dcopy(3*nl,qnl (0,0,1),1,qnl (0,0,2),1)
        call dscal(6*nl,0.5d0,qnl,1)
        io5smm = .true.
      endif

200   format(21x,f10.6)
201   format(1x,a1)
300   format('MOMNTS:  ',a4,'  VRMAX=',f10.6,
     .  /' L',6x,'PL',9x,'Q0',9x,'Q1',9x,'Q2    IDMENU',5x,'DL')
301   format(i2,4f11.7,2x,i3,1x,f12.7,f11.7)
400   format(' AIOMOM: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing moments.$')
401   format(' AIOMOM: ATOM ',a4,':spin 2 input missing;',
     .       ' moments split.$')
402   format(' AIOMOM: I/O error, iostat=',i3,'$')
      end
      logical function io5spp(clabl,lmx,nl,nsp,pp,ves,ifi)
C- File I/O for moments.  ifi > 0 for read, < 0 for write
C ----------------------------------------------------------------------
Ci Inputs: (pp: if file write : ifi < 0)
Ci   ifi   :logical unit: positive for read, negative for write
Ci   clabl :name of the different inequivalent atom
Ci   lmx   :maximum l-quantum number inside the sphere
Ci   nl    :number of l's
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci Inputs/Outputs: (if file read : ifi > 0)
Cio  pp    :potential parameters
Cio  ves   :electrostatic potential at wsr
Cr Remarks:
Cr   aiopar:true unless read error or category not found
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer npp
      parameter(npp=11)
      integer nl,lmx,nsp,ifi
      double precision pp(npp,0:nl-1,nsp),ves
      character*4 clabl
C Local variables:
      integer i,isp,l,l2,k,iost
      double precision x2(7),tiny,vmax
      logical scat
      parameter(tiny=1.d-6,vmax=.999999d9)
      character  chr
      character*4 cdum
      character*80 messg
C External calls:
      external errmsg,scat
C Intrinsic functions:
      intrinsic dabs

      io5spp = .false.

      if (ifi.gt.0) then
C ----- Read potential parameters from atomic file
        if (.not. scat(ifi,'PPAR:',':',.true.)) return
        backspace (ifi)
        read(ifi,300,err=15,end=15,iostat=iost)cdum,ves
        do  10  isp = 1, nsp
          do  12  l = 0, lmx
            read(ifi,200,err=15,end=15,iostat=iost)chr
            if ( chr.ge.'0'.and. chr.le.'9') then
              backspace ifi
              read(ifi,*  ,err=15,end=15,iostat=iost)l2,(x2(k),k=1,7)
            else
              write(messg,401)clabl,isp,l
              call errmsg(messg,1)
              goto 15
            endif
            call dcopy(4,x2,1,pp(1,l,isp),1)
            pp( 9,l,isp)=x2(5)
            pp(10,l,isp)=x2(6)
            pp(11,l,isp)=x2(7)
   12     continue
   10   continue
        io5spp = .true.
      else
C ----- Write potential parameters to atomic file
        write(-ifi,*,err=15,iostat=iost)
        write(-ifi,300,err=15,iostat=iost)clabl,ves
        do  20  isp = 1, nsp
          do  22  l = 0, lmx
            do  24  k=1,4
              x2(k) = pp(k,l,isp)
   24       continue
            x2(5) = pp( 9,l,isp)
            x2(6) = pp(10,l,isp)
            x2(7) = pp(11,l,isp)
            write(-ifi,301,err=15,iostat=iost)l
            do  25  k=1,7
              if(x2(k).gt.-9.9999999d0.and.x2(k).lt.99.9999999d0) then
                write(-ifi,302,err=15,iostat=iost)x2(k)
              else
                write(-ifi,303,err=15,iostat=iost)x2(k)
              endif
   25       continue
            write(-ifi,304,err=15,iostat=iost)
   22     continue
   20   continue
      endif
      return

15    continue
      if (iost.ne.0) then
        write(messg,402)iost
        call errmsg(messg,1)
      endif
C --- handle read exception ---
      if (ifi.gt.0 .and. l.eq.0 .and. isp.eq.2) then
        write(messg,400) clabl
        call errmsg(messg,1)
        call dcopy(6*(lmx+1),pp(1,0,1),1,pp(1,0,2),1)
        io5spp = .true.
      endif

200   format(1x,a1)
300   format('PPAR:  ',a4,'  VES=',f10.6/ ' L',7x,'ENU',8x,'S',7x,
     .       '  O   ',3x,'   P     ',5x,' PHI ',5x,'PHIDOT',
     .       4x,' PGAMMA')
301   format(1x,i1,$)
302   format(f11.7,$)
303   format(e11.5e1,$)
304   format()
400   format(' AIOPAR: ATOM ',a4,':spin 2 input missing; ',
     .          'taken from spin 1.$')
401   format(' AIOPAR: ATOM ',a4,':spin ',i1,', l=',i1,
     .       '; missing potential parameters$')
402   format(' AIOPAR: I/O error, iostat=',i3,'$')
      end
      logical function io5spt(nr,nsp,a,wsr,v,ifi)
C- File I/O for cell potential.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   nr    :number of mesh points
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci Inputs:/Outputs (i if file write (ifi < 0), o if file read (ifi > 0))
Cio  v     :spherical potential
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio  wsr   :Wigner-Seitz radius
Co   aiopot:see remarks
Cr Remarks:
Cr    Format for potential in atomic file begins with
Cr    category 'POT:', followed by a line containing nr, nsp, a, wsr,
Cr    followed by the potential.
Cr    On reading, aiopot returns true only if the category is found,
Cr    the file's value of a and nr match input and wsr is
Cr    close to file's value and the potential is read without error.
Cr    spin-down potential copied from spin-up potential if s-down absent
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer ifi,nr,nsp
      double precision a,wsr,v(nr,nsp)
C Local variables:
      integer i,iost,isp,nr2,nsp2
      double precision a2,wsr2,tiny
      logical scat
      character*300 messg
      parameter(tiny=1.d-5)
C External call
      external errmsg,scat
C Intrinsic functions:
      intrinsic dabs,min0

      io5spt = .false.
      if (ifi .gt. 0) then
C ----- read file
        if (.not. scat(ifi,'POT:',':',.true.)) return
        read(ifi,301,err=15,end=15,iostat=iost)nr2,nsp2,a2,wsr2
        if (dabs(a2-a).gt.tiny .or. nr.ne.nr2 .or.
     .      dabs(wsr2-wsr).gt.tiny) then
          write(messg,400)nr,nr2,a,a2,wsr,wsr2
          call errmsg(messg,-1)
        else
          do  10  isp = 1, min0(nsp2,nsp)
            read(ifi,300,err=15,end=15,iostat=iost)(v(i,isp),i = 1,nr)
            do  12  i = 1, nr
   12       v(i,nsp) = v(i,isp)
   10     continue
          io5spt = .true.
        endif
      else
C ----- write file
        write(-ifi,'(''POT:'')',err=15,iostat=iost)
        write(-ifi,301,err=15,iostat=iost)nr,nsp,a,wsr
        do  20  isp = 1, nsp
          write(-ifi,300,err=15,iostat=iost)(v(i,isp),i=1,nr)
   20   continue
        io5spt = .true.
      endif
      return

C --- handle read exception ---
15    continue
      write(messg,401)iost
      call errmsg(messg,1)

300   format(1p,5d16.9)
301   format(2i5,2f12.5)
400   format(' AIOPOT: difference between data in atomic file ',
     .       '|        and actual values:',
     .       '|        NR ',': actual value:',i8  ,'   read:',i8,
     .       '|        A  ',': actual value:',f8.4,'   read:',f8.4,
     .       '|        WSR',': actual value:',f8.4,'   read:',f8.4,'$')
401   format(' AIOPOT: I/O error, iostat=',i3,'$')
      end
      subroutine aioxt2(i1,i2,tau,pos)
C- Shifts pos(1..3,i1..i2) by a const tau
C     implicit none
      integer i1,i2
      double precision tau(3),pos(3,*)
      integer i,j
      do  i = i1, i2
        do  j = 1, 3
          pos(j,i) = pos(j,i) + tau(j)
        enddo
      enddo
      end
      subroutine aioxt3(k,ssite1,n1,ssite2,n2)
      integer n1,n2,k
      double precision ssite1(k,n1),ssite2(k,n2)

      end

