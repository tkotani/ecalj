C$$$      subroutine supcel(mode,sctrl,slat,ssite,slabl,sspec,plx,nkd,dlat,
C$$$     .  nbx)
C$$$C- Generates basis vectors in a supercell
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   mode  :not used now.  Should be zero
C$$$Ci   sctrl :struct for program flow parameters; see routine uctrl
C$$$Ci     Elts read: nbas nspec
C$$$Ci     Stored:
C$$$Ci     Passed to:
C$$$Ci   slat  :struct for lattice information; see routine ulat
C$$$Ci     Elts read: alat plat
C$$$Ci     Stored:
C$$$Ci     Passed to:
C$$$Ci   ssite :struct for site-specific information; see routine usite
C$$$Ci     Elts read: spec pos
C$$$Ci     Stored:    pos spec pl
C$$$Ci     Passed to: spackv usite dval dpscop spack dvset dvprm
C$$$Ci   slabl :vector of species labels
C$$$Ci     Elts read:
C$$$Ci     Stored:
C$$$Ci     Passed to: baslst iosits
C$$$Ci   sspec :struct for species-specific information; see routine uspec
C$$$Ci     Elts read: z
C$$$Ci     Stored:    name
C$$$Ci     Passed to: spacks iosits
C$$$Ci   plx   :supercell primitive lattice vectors
C$$$Ci   nkd   :number of direct lattice vectors (lattic.f)
C$$$Ci   dlat  :direct lattice vectors
C$$$Co Outputs
C$$$Co   nbx   :number of sites in superlattice
C$$$Cl Local variables
C$$$Cl   qlx   :supercell reciprocal lattice
C$$$Cr Remarks
C$$$Cb Bugs
C$$$Cb   No check is made to see whether ther are enough dlat to encompass
C$$$Cb   all lattice vectors
C$$$Cu Updates
C$$$Cu   01 Aug 08 New --wsitex
C$$$Cu   19 Apr 03 Changed distortion to read from displacements file
C$$$Cu   12 Apr 03 Some bug fixes when partial list of sites saved
C$$$Cu   09 Apr 03 Added pairwise distortion; also writes original cell
C$$$Cu             site sdata first
C$$$Cu   24 Nov 02 R.S. pairwise hamiltonian may be created for supercell
C$$$Cu             from file rsj
C$$$Cu   26 Nov 02 supcel can restrict supercell sites to a specified list
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$C ... Passed parameters
C$$$      character*8 slabl(1)
C$$$      integer nkd,mode,nbx
C$$$      double precision sctrl(1),slat(1),ssite(1),sspec(1),dlat(3,nkd)
C$$$C ... Local parameters
C$$$      character*72 alabl*8,outs,sortex(5),out2,jfilnm,dc*1
C$$$      logical cmdopt,lsort,lfirst,latvec
C$$$      integer lgunit,i,j,k,m,ib,nbas,is,lok,ndsite,iosits,lio,ifi,nspec,
C$$$     .  ixv(36),imark(36),it(3),ip,parg,j1,j2,nexpr,nbxsav,a2vec,nmark,
C$$$     .  fopna,fxst,nfilj,nfilm,nw,stdo,fopn,rdm,nn,nc,iprint
C$$$      integer mxcsiz,nbasp,ifirsj,iorsj,nttab,nlst,ilst,nlst1
C$$$      integer ontab,oiax,ortab,mxcsz,ohpair,oaamom,oimap,oamom2,oips,oz,
C$$$     .  oilst,olst1,odisp,opos2,owk,oiprm,ossite,opos,odpos
C$$$      double precision sclrsj
C$$$      double precision plat(3,3),plx(3,3),qlx(3,3),qlat(3,3)
C$$$      double precision dd(3),tau(3),x(3),vol,volx,dval,tol,alat,rcut,xx
C$$$C ... Heap
C$$$      integer w(1)
C$$$      common /w/ w

C$$$      call upack1('lat alat',slat,alat)
C$$$      call upack1('lat plat',slat,plat)
C$$$      call upack1('ctrl nbas',sctrl,nbas)
C$$$      call upack1('ctrl nspec',sctrl,nspec)
C$$$      tol = 1d-6
C$$$      stdo = lgunit(1)

C$$$C --- Get list of lattice vectors to include ---
C$$$      nmark = 0
C$$$      if (cmdopt('-use:',5,0,outs)) then
C$$$        m = 5
C$$$        k = a2vec(outs,72,m,2,',; ',3,-3,36,ixv,imark)
C$$$        nmark = k/3
C$$$        if (3*nmark .ne. k) call rxi('SUPCEL: ''-use:'' takes'//
C$$$     .    ' multiples of 3 arguments, but read',k)

C$$$        do  10  ib = 1, nmark
C$$$C         Assemble vector as specified multiples of prim. latt. vec.
C$$$          call dpzero(x,3)
C$$$          call daxpy(3,dble(imark(3*ib-2)),plat(1,1),1,x,1)
C$$$          call daxpy(3,dble(imark(3*ib-1)),plat(1,2),1,x,1)
C$$$          call daxpy(3,dble(imark(3*ib-0)),plat(1,3),1,x,1)
C$$$C         Mark which of dlat it belongs to
C$$$          do  12  k = 1, nkd
C$$$            if (tol .gt. abs(x(1)-dlat(1,k))+
C$$$     .                   abs(x(2)-dlat(2,k))+
C$$$     .                   abs(x(3)-dlat(3,k))) then
C$$$              imark(ib) = k
C$$$              goto 14
C$$$            endif
C$$$   12   continue
C$$$        call rx('supcel: missing dlat')
C$$$   14   continue
C$$$   10   continue
C$$$      endif

C$$$C --- Optionally get a list (subset) of sites in original basis ---
C$$$      nlst = 0
C$$$      call defi(oilst,-nbas)
C$$$      if (cmdopt('--sites',7,0,outs)) then
C$$$        i = 8
C$$$        call sp2cls('site spec',ssite,0,1,1,nbas,oips)
C$$$        call sp2cls('spec z',sspec,0,1,1,nspec,oz)
C$$$C       call awrit2('%N Z= %n:1,3;3d',' ',100,6,nspec,w(oz))
C$$$        dc = outs(i:i)
C$$$        call baslst(0,11,outs(i:),j2,w(oips),nbas,slabl,w(oz),0,' ',x,
C$$$     .    nlst,w(oilst))
C$$$        if (nlst .eq. 0) w(oilst) = 0
C$$$      endif

C$$$C --- Supercell reciprocal lattice vectors ---
C$$$      write(stdo,'(/a)') ' Basis vectors of supercell:'
C$$$      call dinv33(plx,1,qlx,volx)
C$$$      call dinv33(plat,0,qlat,vol)
C$$$      if (dabs(volx) .lt. 1d-6)
C$$$     .   call rx('supercell vectors are linearly dependent')
C$$$      write(stdo,351)
C$$$  351 format(/15x,'plx',30x,'qlx')
C$$$      write(stdo,350) ((plx(m,k),m=1,3),(qlx(m,k),m=1,3),k=1,3)
C$$$  350 format(3f10.5,5x,3f10.5)

C$$$      call ivset(ixv,1,3,2)
C$$$      if (cmdopt('-shorten=',9,0,outs) .or.
C$$$     .    cmdopt('--shorten=',10,0,outs)) then
C$$$        ip = 9
C$$$        if (cmdopt('--shorten=',10,0,outs)) ip = 10
C$$$        call skipbl(outs,len(outs),ip)
C$$$        if (parg(' ',2,outs,ip,len(outs),', ',2,3,it,ixv) .lt. 0)
C$$$     .    call rxs2('supcel: failed to parse "',outs(1:ip+5),' ..."')
C$$$        print '(/''  supcel: shorten vectors with mode ='',3i2)',ixv
C$$$      endif

C$$$C --- Count supercell basis vectors (nbx) ---
C$$$C     and make table mapping supercell sites to original sites
C$$$      write (stdo,*) ' '
C$$$      nbx = 0
C$$$      ilst = 0
C$$$      do  20  ib = 1, nbas
C$$$        if (nlst .gt. 0) then
C$$$   11     if (ilst .ge. nlst) goto 20
C$$$          if (w(oilst+ilst) .lt. ib) then
C$$$            ilst = ilst+1
C$$$            goto 11
C$$$          endif
C$$$          if (w(oilst+ilst) .ne. ib) goto 20
C$$$        endif

C$$$        call upack('site spec pos',ssite,ib,is,tau,0,0)
C$$$        call spacks(0,'spec name',sspec,alabl,is,is)
C$$$        do  22  k = 1, nkd
C$$$          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C$$$          if (lok .eq. 1) write(stdo,120) ib,alabl,(dd(m),m=1,3)
C$$$  120     format(i6,2x,a,3f11.6)
C$$$   22   continue
C$$$   20 continue
C$$$C     Make table imap mapping supercell sites to original sites
C$$$      call defi(oimap,-nbx)
C$$$      nbx = 0
C$$$      ilst = 0
C$$$      do  24  ib = 1, nbas
C$$$        if (nlst .gt. 0) then
C$$$   21     if (ilst .ge. nlst) goto 24
C$$$          if (w(oilst+ilst) .lt. ib) then
C$$$            ilst = ilst+1
C$$$            goto 21
C$$$          endif
C$$$          if (w(oilst+ilst) .ne. ib) goto 24
C$$$        endif
C$$$        call upack('site spec pos',ssite,ib,is,tau,0,0)
C$$$        call spacks(0,'spec name',sspec,alabl,is,is)
C$$$        do  26  k = 1, nkd
C$$$          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C$$$          if (lok .eq. 1) w(oimap+nbx-1) = ib
C$$$   26   continue
C$$$   24 continue

C$$$C --- Printout ---
C$$$      write(stdo, 333) abs(vol), abs(volx), abs(volx/vol)
C$$$  333 format(/' cell vol=',f10.6,'     scell vol=',f12.6,
C$$$     .        '     ratio=',f7.2)
C$$$      write(stdo, 334) nbas, nbx, dble(nbx)/nbas
C$$$  334 format( ' cell nbas=',i9,'     scell nbas=',i11,
C$$$     .        '     ratio=',f7.2)
C$$$      if (nlst .eq. 0) then
C$$$        if (abs(dble(nbx)/nbas-abs(volx/vol)) .gt. 1d-5)
C$$$     .    write(stdo, 335)
C$$$      else
C$$$        if (abs(dble(nbx)/nlst-abs(volx/vol)) .gt. 1d-5)
C$$$     .    write(stdo, 335)
C$$$      endif
C$$$  335 format(' (warning) supercell has incommensurate lattice')

C$$$      lsort = cmdopt('-sort:',6,0,outs) .or. cmdopt('--sort:',7,0,outs)
C$$$      lfirst= cmdopt('-first',6,0,outs) .or. cmdopt('--first',7,0,outs)

C$$$C --- Set up expressions to sort ---
C$$$      i = 7
C$$$      if (cmdopt('--sort:',7,0,outs)) i = 8
C$$$      nexpr = 0
C$$$      do  48  k = 1, 3
C$$$        call word(outs(i:),k,j1,j2)
C$$$        if (j2 .ge. j1) then
C$$$          nexpr = nexpr + 1
C$$$          sortex(nexpr) = outs(6+j1:6+j2)
C$$$        else
C$$$          goto 49
C$$$        endif
C$$$   48 continue
C$$$   49 continue

C$$$      if (.not. (cmdopt('-wsite=',7,0,outs) .or.
C$$$     .           cmdopt('--wsite=',8,0,outs) .or.
C$$$     .           cmdopt('--wsitex=',9,0,outs))) return

C$$$C --- Get list of positions; sort them ---
C$$$      if (lsort .or. lfirst) then
C$$$      call defrr(opos,3*nbx)
C$$$      call defi(owk,-nbx)
C$$$      call defi(oiprm,nbx)
C$$$      nbxsav = nbx
C$$$      nbx = 0
C$$$      ilst = 0
C$$$      do  40  ib = 1, nbas
C$$$        if (nlst .gt. 0) then
C$$$   41     if (ilst .ge. nlst) goto 40
C$$$          if (w(oilst+ilst) .lt. ib) then
C$$$            ilst = ilst+1
C$$$            goto 41
C$$$          endif
C$$$          if (w(oilst+ilst) .ne. ib) goto 40
C$$$        endif

C$$$        call upack('site spec pos',ssite,ib,is,tau,0,0)
C$$$        call spacks(0,'spec name',sspec,alabl,is,is)
C$$$        do  42  k = 1, nkd
C$$$          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C$$$          if (lok .eq. 1) then
C$$$            call dmscop(w(opos),nbxsav,dd,1,1,1,1,3,nbx,1,1d0)
C$$$            w(owk+nbx-1) = nbx
C$$$C           case lfirst: if dd = tau, this was the original vector
C$$$            if (lfirst) then
C$$$              xx = (dd(1)-tau(1))**2+(dd(2)-tau(2))**2+(dd(3)-tau(3))**2
C$$$              x(1) = (dd(1)-tau(1))
C$$$              x(2) = (dd(2)-tau(2))
C$$$              x(3) = (dd(3)-tau(3))
C$$$              if (latvec(1,tol,qlx,x)) then
C$$$                w(owk+nbx-1) = ib-nbas
C$$$              endif
C$$$            endif
C$$$          endif
C$$$   42   continue
C$$$   40 continue

C$$$C     call prmx('unsorted pos',w(opos),nbx,nbx,3)

C$$$      if (lsort) then
C$$$        call defrr(opos2,3*nbx)
C$$$        call dcopy(3*nbx,w(opos),1,w(opos2),1)
C$$$        call mapdat(nexpr,sortex,' ',w(oiprm),nbx,3,w(opos2),w)
C$$$        call defrr(owk,nexpr*nbx)
C$$$        call dmcpy(w(opos2),nbx,1,w(owk),1,nexpr,nbx,nexpr)
C$$$C       call prmx('expressions',w(owk),nexpr,nexpr,nbx)
C$$$        call dcopy(3*nbx,w(opos),1,w(opos2),1)
C$$$        call dvheap(nexpr,nbx,w(owk),w(oiprm),0d0,1)
C$$$C       call msort(nbx,3,nexpr,w(oiprm),w(owk),w(opos2),w(opos))
C$$$C       call prmx('sorted pos',w(opos),nbx,nbx,3)
C$$$      else
C$$$        call ivheap(1,nbx,w(owk),w(oiprm),1)
C$$$C       call yprm('ordering',0,w(owk),0,nbx,nbx,1)
C$$$C       call yprm('iprm',0,w(oiprm),0,nbx,nbx,1)
C$$$      endif

C$$$      endif

C$$$      if (.not. (cmdopt('-wsite=',7,0,outs) .or.
C$$$     .           cmdopt('--wsite=',8,0,outs) .or.
C$$$     .           cmdopt('--wsitex=',9,0,outs))) return


C$$$C --- Read R.S. hamiltonian from file rsj ---
C$$$C     ifirsj is the file logical unit for real-space hamiltonian
C$$$C     ifirsj = 0 flags that no file will be written
C$$$C     Switch syntax:  --wrsj[:fn=name][:scl=#]
C$$$      ifirsj = 0
C$$$      if (cmdopt('--wrsj',6,0,jfilnm) .or.
C$$$     .    cmdopt('-wrsj',5,0,jfilnm)) then

C$$$        out2 = jfilnm(7:)
C$$$        if (jfilnm(1:5) .eq. '-wrsj') out2 = jfilnm(6:)
C$$$C       This is default name for supercell rsj file
C$$$        jfilnm = 'rsj2'

C$$$        dc = out2(1:1)
C$$$        sclrsj = 1
C$$$        if (dc .ne. ' ') then
C$$$C     ... Return here to resume parsing for arguments
C$$$          j2 = 0
C$$$   50     continue
C$$$          j2 = j2+1
C$$$          if (out2(j2:j2) .eq. dc) goto 50
C$$$          j1 = min(len(out2),j2)
C$$$          call nwordg(out2,0,dc//' ',1,j1,j2)
C$$$          if (j2 .ge. j1) then
C$$$            if (.false.) then
C$$$            elseif (out2(j1:j1+3) .eq. 'scl=')  then
C$$$              j = 0
C$$$              i = parg('scl=',4,out2(j1:),j,len(out2(j1:)),
C$$$     .          dc//' ',1,1,k,sclrsj)
C$$$              if (i .le. 0) goto 59

C$$$            elseif (out2(j1:j1+2) .eq. 'fn=')  then
C$$$              jfilnm = out2(j1+5:j2)
C$$$            else
C$$$              goto 59
C$$$            endif
C$$$            goto 50
C$$$   59       continue
C$$$            call rxs('supcel: failed to parse wrsj option', out2)
C$$$          endif
C$$$        endif

C$$$        call info0(20,1,0,
C$$$     .    ' ... Writing supercell r.s. hamiltonian to file '//jfilnm)

C$$$C   ... Setup call for file read
C$$$        mxcsiz = 0
C$$$        nbasp = nbas
C$$$        if (fxst('rsj') .ne. 1) then
C$$$          call rx('supcel: specified --wrsj but no file rsj')
C$$$        endif
C$$$        ifirsj = fopna('rsj',-1,0)
C$$$        call defrr(oaamom,nbas)
C$$$        nfilj = iorsj(0,ifirsj,nbas,alat,plat,w(oaamom),w,w,w,w,x)

C$$$C       This makes iax table and rtab suitable for iorsj
C$$$        call defrr(opos,3*nbas)
C$$$        call spackv(10,'site pos',ssite,1,nbas,w(opos))
C$$$        call pairs(nbas,nbasp,1d0,plat,x(1)/2*(1+1d-6),w(opos),
C$$$     .    -1,3,-1,w,nttab,ontab,oiax,mxcsiz)
C$$$        call defrr(ortab,3*nttab)
C$$$        call mkrtab(000,1d0,plat,w(opos),w(oiax),nttab,w(opos),w(ortab))

C$$$C       Read pair hamiltonian for entries matching iax table
C$$$        call defrr(ohpair,nttab)
C$$$        nfilm = iorsj(03,ifirsj,nbas,w,plat,w(oaamom),w(oiax),w(ortab),
C$$$     .    w(ontab),w(ohpair),x)

C$$$        if (nfilj.ne.nfilm) call info2(20,0,0,'%8f(warning) only'//
C$$$     .    ' matched %i out of %i file entries',nfilm,nfilj)
C$$$        if (nfilm.ne.nttab) then
C$$$          call info2(20,0,0,'%8f(warning) iax table has %i pairs'//
C$$$     .      ' but matched only %i file entries',nttab,nfilm)

C$$$C         Purge iax table of unused entries
C$$$          call psup3(nlst,w(oilst),w(oiax),nttab,nfilm,w(ortab))
C$$$          call symiax(1,plat,nbas,w,w,w,0,w(ontab),w(oiax),nttab,mxcsz)
C$$$          if (nttab .ne. nfilm) then
C$$$            call fexit2(-1,111,' Exit -1 : iorsj: bug in symiax:'//
C$$$     .        ' nttab=%i, nfilm=%i',nttab,nfilm)
C$$$          endif

C$$$C         Remake rtab
C$$$          call mkrtab(000,1d0,plat,w(opos),w(oiax),nttab,w(opos),
C$$$     .      w(ortab))

C$$$C         Re-read pair hamiltonian with new iax table
C$$$          nfilm = iorsj(03,ifirsj,nbas,w,plat,w(oaamom),w(oiax),
C$$$     .      w(ortab),w(ontab),w(ohpair),x)
C$$$          if (nttab .ne. nfilm) call rx('supcel: bug in iorsj')

C$$$        endif

C$$$C       Close input file rsj
C$$$        call fclose(ifirsj)

C$$$C       Open supercell file; write header
C$$$        call defrr(oamom2,nbx)
C$$$        call cpvprm(0,1,nbx,w(oimap),w(oaamom),w(oamom2))
C$$$        ifirsj = fopna(jfilnm,-1,0)
C$$$        i = iorsj(0,-ifirsj,nbx,alat,plx,w(oamom2),w,w,w,w,w)

C$$$      endif

C$$$C --- Create supercell ssite ---
C$$$      call defrr(ossite,1000)
C$$$      call usite(w(ossite),-1,0,1,0,0,0)
C$$$      ndsite = nint(dval(w(ossite),1))
C$$$      call redfrr(ossite,nbx*ndsite)
C$$$      call usite(w(ossite),-1,0,nbx,0,0,0)
C$$$      nspec = 0
C$$$      nbx = 0

C$$$      ilst = 0
C$$$      do  30  ib = 1, nbas
C$$$        if (nlst .gt. 0) then
C$$$   31     if (ilst .ge. nlst) goto 30
C$$$          if (w(oilst+ilst) .lt. ib) then
C$$$            ilst = ilst+1
C$$$            goto 31
C$$$          endif
C$$$          if (w(oilst+ilst) .ne. ib) goto 30
C$$$        endif

C$$$        call upack('site spec pos',ssite,ib,is,tau,0,0)
C$$$        call spacks(0,'spec name',sspec,alabl,is,is)
C$$$        nspec = max(nspec,is)
C$$$        do  32  k = 1, nkd
C$$$          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C$$$          if (lok .eq. 1) then

C$$$C       ... Copy parameters from site ib to site nbx
C$$$            call dpscop(ssite,w(ossite),ndsite,
C$$$     .        1+(ib-1)*ndsite,1+(nbx-1)*ndsite,1d0)
C$$$C           12 Apr 03 ssite(1) may not be struc size if 1st ib<>1
C$$$            call dvset(w(ossite),1,1,ssite)
C$$$            call pack5('site spec pos',w(ossite),nbx,is,dd,0,0)

C$$$C       ... Write rsj to supercell file
C$$$            if (ifirsj .ne. 0) then
C$$$            call psup2(ifirsj,ib,nbx,w(ontab),w(ortab),sclrsj,w(ohpair))
C$$$            endif

C$$$C       ... PL index
C$$$            if (cmdopt('--pl:',5,0,out2)) then
C$$$              call prjpos(2,ixv,w,plat,1,dd,x)
C$$$              call mapdat(1,out2(6:),' ',m,1,3,x,w)
C$$$              call spack(2,'site pl',ssite,1,m,j1,j2,k)
C$$$              call dvset(w(ossite),m+(nbx-1)*ndsite,m+(nbx-1)*ndsite,x)
C$$$            endif
C$$$          endif
C$$$   32   continue
C$$$   30 continue

C$$$C --- Permute site strux if sorted ---
C$$$      if (lsort .or. lfirst) then
C$$$        call defrr(owk,nbx*ndsite)
C$$$        call dvprm(ndsite,nbx,w(ossite),w(owk),w(oiprm),.true.)
C$$$        call dvset(w(ossite),1,1,dble(ndsite))
C$$$        call rlse(owk)
C$$$      endif

C$$$C --- Displacement around an atom ---
C$$$      call defrr(opos,3*nbx)
C$$$      call spackv(10,'site pos',w(ossite),1,nbx,w(opos))
C$$$      call defrr(odpos,-3*nbx)
C$$$      if (cmdopt('-disp',5,0,out2) .or.
C$$$     .    cmdopt('--disp',6,0,out2)) then
C$$$        dc = out2(6:6)
C$$$        if (out2(1:6) .eq. '--disp') dc = out2(7:7)
C$$$        call wrdsg(out2,1,dc//' ',nw)
C$$$        if (nw .lt. 3) call rx('supcel: --disp requires 2 arguments')

C$$$C   ... Read displacements file
C$$$        call wordg(out2,0,dc//' ',2,j1,j2)
C$$$        ifi = fopn(out2(j1:j2))
C$$$        rewind ifi
C$$$C       Determine size of file
C$$$        nn = 0
C$$$        nc = 0
C$$$        j = rdm(ifi,0,0,' ',xx,nn,nc)
C$$$        if (j .lt. 0) call rxs
C$$$     .    ('supcel: failed to read displacements file ',out2(j1:j2))
C$$$        if (nc .ne. 6) call rxs('supcel: displacements file '//
C$$$     .    'apparently does not have 6 columns: ',out2(j1:j2))
C$$$        call defrr(odisp,nn*nc)
C$$$        rewind ifi
C$$$        j = rdm(ifi,0,nn*nc,' ',w(odisp),nn,nc)
C$$$        call psup5(nn,w(odisp),rcut)
C$$$        call wordg(out2,0,dc//' ',3,i,j)
C$$$        call info(30,1,0,' ... Displacements from file '//
C$$$     .    out2(j1:j2)//' for site list='//out2(i:j)//'. rcut=%d',rcut,0)

C$$$        nlst1 = 0
C$$$        call defi(olst1,-nbx)
C$$$        call wordg(out2,0,dc//' ',3,j1,j2)
C$$$        if (j2 .lt. j1) call rx('supcel: --disp requires 2 arguments')
C$$$        call sp2cls('site spec',w(ossite),0,1,1,nbx,oips)
C$$$        call sp2cls('spec z',sspec,0,1,1,nspec,oz)
C$$$        call baslst(1,11,out2(j1-1:),j2,w(oips),nbx,slabl,w(oz),0,' ',
C$$$     .    x,nlst1,w(olst1))

C$$$C       This makes iax table and rtab
C$$$        mxcsiz = 0
C$$$C       call setpr(50)
C$$$        call pshpr(iprint()-10)
C$$$        call pairs(nbx,nbx,1d0,plx,rcut/2*(1+1d-6),w(opos),
C$$$     .    -1,3,-1,w,nttab,ontab,oiax,mxcsiz)
C$$$        call defrr(ortab,3*nttab)
C$$$        call mkrtab(000,1d0,plx,w(opos),w(oiax),nttab,w(opos),w(ortab))
C$$$        call poppr
C$$$        call ivset(ixv,1,3,2)
C$$$        call shorps(nttab,plx,ixv,w(ortab),w(ortab))

C$$$C       Get displacements
C$$$        call psup4(w(oiax),nttab,
C$$$     .    nlst1,w(olst1),w(ortab),w(opos),w(odisp),nn,w(odpos))

C$$$C       Add to ssite
C$$$        call daxpy(3*nbx,1d0,w(odpos),1,w(opos),1)
C$$$        call spackv(11,'site pos',w(ossite),1,nbx,w(opos))
C$$$      endif

C$$$C --- Write site file ---
C$$$      lio = 1000*(2+4+8+16+32) + 1
C$$$C     lio = lio+16000
C$$$      if (outs(1:7) .eq. '-wsite=') then
C$$$        out2 = outs(8:)
C$$$        outs = out2
C$$$      elseif (outs(1:8) .eq. '--wsite=') then
C$$$        out2 = outs(9:)
C$$$        outs = out2
C$$$      elseif (outs(1:9) .eq. '--wsitex=') then
C$$$        out2 = outs(10:)
C$$$        outs = out2
C$$$        lio = lio + 10
C$$$      endif
C$$$      if (iosits(lio,3d0,0,outs,ifi,slabl,alat,plx,nbx,nspec,sspec,
C$$$     .  w(ossite)).lt. 0) call rx('failed to write ssite')

C$$$C --- Write positions to file ---
C$$$      if (cmdopt('--wpos=',7,0,outs)) then
C$$$        call defrr(opos,3*nbx)
C$$$        call spackv(10,'site pos',w(ossite),1,nbx,w(opos))
C$$$        call iopos(.true.,0,outs(8:),nbx,w(opos))
C$$$      endif


C$$$      return
C$$$C 999 continue

C$$$      end

C$$$      subroutine psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C$$$C- Kernel called by supcel: to decide whether a vector should be added
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   nmark :?
C$$$Ci   imark :?
C$$$Ci   tau   :position vector for this site
C$$$Ci   k     :which of the lattice vectors in dlat to check
C$$$Ci   dlat  :a list of lattice vectors, of which dlat(k) is to be used here.
C$$$Ci   plx   :(passed to shorps; used if vectors are to be shortened)
C$$$Ci   qlx   :reciprocal lattice vectors to supercell
C$$$Ci   ixv   :(passed to shorps; used if vectors are to be shortened)
C$$$Ci         :0 suppresses shifts along plat(j)
C$$$Ci         :1 shifts to unit cell at origin (pos in 1st quadrant)
C$$$Ci         :2 shifts to minimize length of pos
C$$$Co Outputs
C$$$Co   lok   :0 do not add this vector to the list
C$$$Co         :1 add this vector to the list
C$$$Co   dd    :connecting vector
C$$$Co   nbx   :incremented by one if this vector should be included
C$$$Cl Local variables
C$$$Cl         :
C$$$Cr Remarks
C$$$Cr
C$$$Cu Updates
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$C ... Passed parameters
C$$$      integer nmark,imark(nmark),lok,nbx,ixv(3),k
C$$$      double precision dlat(3,k),tau(3),plx(3,3),qlx(3,3),dd(3)
C$$$C ... Local parameters
C$$$      double precision x(3)
C$$$      integer m
C$$$      logical cmdopt
C$$$      character*72 outs

C$$$      dd(1) = dlat(1,k)+tau(1)
C$$$      dd(2) = dlat(2,k)+tau(2)
C$$$      dd(3) = dlat(3,k)+tau(3)
C$$$      if (nmark .gt. 0) then
C$$$        lok = 0
C$$$        do  10  m = 1, nmark
C$$$   10   if (imark(m) .eq. k) lok = 1
C$$$      else
C$$$        lok = 1
C$$$        do  20  m = 1, 3
C$$$          x(m) = dd(1)*qlx(1,m)+dd(2)*qlx(2,m)+dd(3)*qlx(3,m)
C$$$          if (x(m).lt.-0.001d0 .or. x(m).gt.0.999d0) lok = 0
C$$$   20   continue
C$$$      endif
C$$$      if (lok .eq. 1) then
C$$$        nbx = nbx+1
C$$$        if (cmdopt('-shorten',8,0,outs)) call shorps(1,plx,ixv,dd,dd)
C$$$        if (cmdopt('--shorten',9,0,outs)) call shorps(1,plx,ixv,dd,dd)
C$$$C 120   format(2x,a,3f11.6)
C$$$      endif
C$$$      end

C$$$      subroutine psup2(ifi,ib,ib2,ntab,rtab,sclrsj,rsj)
C$$$C- Kernel called by supcel: write r.s. pairwise ham for this site
C$$$C  Note that jb2 of supercell is not known, so we just write zero
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   ifi   :file handle
C$$$Ci   ib    :site index of original cell
C$$$Ci   ib2   :site index of supercell
C$$$Ci   iax   :neighbor table containing pair information (pairc.f)
C$$$Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
C$$$Ci   rsj
C$$$Co Outputs
C$$$Cl Local variables
C$$$Cl         :
C$$$Cr Remarks
C$$$Cr
C$$$Cu Updates
C$$$Cu   24 Nov 02
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$C ... Passed parameters
C$$$      integer ifi,ib,ib2,ntab(ib+1)
C$$$      double precision sclrsj
C$$$C     integer niax
C$$$C     parameter (niax=10)
C$$$C     integer iax(niax,*)
C$$$      double precision rsJ(*),rtab(3,*)
C$$$C ... Local parameters
C$$$      integer it

C$$$      do  it = ntab(ib)+1, ntab(ib+1)
C$$$C        call awrit6('%,4i%,4i%;12,7D%;12,7D%;12,7D%;12,4D',' ',80,ifi,
C$$$C     .    ib2,0,rtab(1,it),rtab(2,it),rtab(3,it),sclrsj*rsj(it))
C$$$        call awrit6('%,6i%,6i%;12,7D%;12,7D%;12,7D %;6g',' ',80,ifi,
C$$$     .    ib2,0,rtab(1,it),rtab(2,it),rtab(3,it),sclrsj*rsj(it))
C$$$      enddo
C$$$      end

C$$$      subroutine psup3(nlst,ilst,iax,nttab,n,rtab)
C$$$C- Kernel called by supcel to purge pairs outside of site list
C$$$C  used when a subset of all sites is specified
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   nlst  :number of sites in site list
C$$$Ci   ilst  :site list
C$$$Ci   iax   :neighbor table containing pair information (pairc.f)
C$$$Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
C$$$Co Outputs
C$$$Co   n     :number of pairs after purging
C$$$Co         :NB: n is not set if nlst=0
C$$$Cr Remarks
C$$$Cr
C$$$Cu Updates
C$$$Cu   26 Nov 02 First created
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$C ... Passed parameters
C$$$      integer niax,nlst,ilst(nlst),nttab
C$$$      parameter (niax=10)
C$$$      integer iax(niax,nttab)
C$$$      double precision rtab(3,nttab)
C$$$C ... Local parameters
C$$$      integer i,ib,jb,low,n

C$$$      if (nlst .eq. 0) return
C$$$      n = 0
C$$$      do  i = 1, nttab
C$$$        ib = iax(1,i)
C$$$        low = 1
C$$$        call hunti(ilst,nlst,ib,0,low)
C$$$        if (low .ge. nlst) then
C$$$C         print *, 'i, iax',i,iax(1,i)
C$$$          iax(1,i) = 0
C$$$        elseif (ilst(low+1) .ne. ib) then
C$$$C         print *, 'i, iax',i,iax(1,i)
C$$$          iax(1,i) = 0
C$$$        else
C$$$          jb = iax(2,i)
C$$$          call hunti(ilst,nlst,jb,0,low)
C$$$          if (low .ge. nlst) then
C$$$C           print *, 'i, iax',i,iax(1,i),iax(2,i)
C$$$            iax(1,i) = 0
C$$$          elseif (ilst(low+1) .ne. jb) then
C$$$C           print *, 'i, iax',i,iax(1,i),iax(2,i)
C$$$            iax(1,i) = 0
C$$$          else
C$$$            n = n+1
C$$$          endif
C$$$        endif
C$$$      enddo

C$$$      call info2(20,0,0,' supcel: site-list further reduced pair'//
C$$$     .  ' table to %i pairs',n,0)

C$$$      end

C$$$      subroutine mapdat(nexpr,expr,sincl,iwk,nr,nc,dat,dat2)
C$$$C- Replace columns of dat with algebraic expressions of them
C$$$C     implicit none
C$$$C Passed parameters
C$$$      integer nexpr,nr,nc,iwk(nr)
C$$$      character*(*) expr(1), sincl
C$$$      double precision dat(nr,nc),dat2(*)
C$$$C Local variables
C$$$      integer ir,iv0,ival,i,j,ii,jr
C$$$      character*4 xn, outs*80
C$$$      logical a2bin,logi

C$$$      call numsyv(iv0)
C$$$      jr = 0
C$$$      do  20  ir = 1, nr
C$$$        call clrsyv(iv0)
C$$$        call addsyv('i',dble(ir),ival)

C$$$C   ... Load data table
C$$$        do  22  j = 1, nc
C$$$          ii = 1
C$$$          xn = 'x   '
C$$$          call bin2a('(i2)',0,0,j,2,0,4,xn,ii)
C$$$          call addsyv(xn,dat(ir,j),ival)
C$$$   22   continue

C$$$C   ... Exclude points if not satisfy sincl
C$$$        logi = .true.
C$$$        if (sincl .ne. ' ') then
C$$$          j = 0
C$$$          if (.not. a2bin(sincl,logi,0,0,' ',j,-1))
C$$$     .      call rx('mapdat:  error parsing sincl')
C$$$        endif

C$$$C   ... Put expressions of these vars into dat
C$$$        do  24  i = 1, nexpr
C$$$          j = 0
C$$$          if (.not. a2bin(expr(i),dat(ir,i),4,0,' ',j,-1)) then
C$$$            outs = expr(i)
C$$$            call skpblb(expr(i),len(expr(i)),ii)
C$$$            call fexit(-1,1,'MAPDAT: failed to parse expr: '
C$$$     .        //outs(1:j+1)//'<-->'//outs(j+2:ii+1),0)
C$$$          endif
C$$$   24   continue
C$$$        if (logi) then
C$$$          jr = jr+1
C$$$          iwk(jr) = ir
C$$$        endif
C$$$   20 continue
C$$$      if (jr .ne. nr) then
C$$$        call dcopy(nr*nexpr,dat,1,dat2,1)
C$$$        call rowpmt(nr,jr,nexpr,iwk,dat2,dat)
C$$$        nr = jr
C$$$      endif
C$$$      call clrsyv(iv0)

C$$$      end

C$$$      subroutine msort(nr,ncs,ncw,iprm,wk,s,ssort)
C$$$C- sort matrix according to increasing wk.
C$$$C     implicit none
C$$$      integer nr,ncs,ncw,iprm(nr)
C$$$      double precision wk(ncw,nr),s(nr,ncs),ssort(nr,ncs)
C$$$      integer ir,ic,ix

C$$$      call dvheap(ncw,nr,wk,iprm,0d0,1)
C$$$      do  20  ir = 1, nr
C$$$        ix = iprm(ir)
C$$$        do  30  ic = 1, ncs
C$$$   30   ssort(ir,ic) = s(ix,ic)
C$$$   20 continue
C$$$      end

C$$$      subroutine rowpmt(nf,nt,nc,ipmt,afrom,ato)
C$$$C     implicit none
C$$$      integer nf,nt,nc,ipmt(nt),i
C$$$      double precision afrom(nf,nc), ato(nt,nc)

C$$$      do  10  i = 1, nt
C$$$   10 call dcopy(nc,afrom(ipmt(i),1),nf,ato(i,1),nt)
C$$$      end

C$$$      subroutine psup4(iax,nttab,nlst1,lst1,rtab,pos,disp,nn,dpos)
C$$$C- Kernel called by supcel to displace pairs radially
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   iax   :neighbor table containing pair information (pairc.f)
C$$$Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
C$$$Ci   nlst1 :number of sites of "center" type
C$$$Ci   lst1  :list of sites of "center" type
C$$$Ci   nlst2 :number of sites of "neighbor" type
C$$$Ci   lst2  :list of sites of "neighbor" type
C$$$Ci   rtab  :site positions corresponding to entries in a neighbor table
C$$$Ci          relative to some center
C$$$Co Inputs/Outputs
C$$$Cio   pos  :On input, basis vectors
C$$$Cio        :On output, displaced basis vectors
C$$$Cu Updates
C$$$Cu   07 Apr 03 First created
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$C ... Passed parameters
C$$$      integer nlst1,lst1(nlst1),niax,nttab,nn
C$$$      double precision pos(3,*),rtab(3,*),dpos(3,*),disp(nn,6)
C$$$      parameter (niax=10)
C$$$      integer iax(niax,1)
C$$$C ... Local parameters
C$$$      integer iat,low1,i1,i2,ipr,count,stdo,lgunit,i,j
C$$$      double precision rhat(3),tol

C$$$      count = 0
C$$$      low1 = 0
C$$$      call getpr(ipr)
C$$$      ipr = 100
C$$$      stdo = lgunit(1)
C$$$      tol = 1d-6

C$$$C      if (ipr .ge. 50) write(stdo,330) nlst1
C$$$C  330 format(' supcel: add displacement around',i4,' site(s)'/
C$$$C     .  '  ia  ib',15x,'rb-ra',27x,'delta pos(ib)')
C$$$      if (ipr .ge. 50) write(stdo,330)
C$$$  330 format('  ia  ib',15x,'rb-ra',27x,'delta pos(ib)')

C$$$      do  iat = 1, nttab

C$$$C   ... If the center atom isn't in the list, skip this pair
C$$$        i1 = iax(1,iat)
C$$$        call hunti(lst1,nlst1,i1,0,low1)
C$$$        if (i1 .ne. lst1(low1+1)) goto 10

C$$$C   ... See if this pair is in any of in the disp list
C$$$        call dcopy(3,rtab(1,iat),1,rhat,1)
C$$$        do  i = 1, nn
C$$$          if (abs(rhat(1)-disp(i,1)) .lt. tol) then
C$$$            if (abs(rhat(2)-disp(i,2)) .lt. tol) then
C$$$C             We have a match!  Add to dpos
C$$$              if (abs(rhat(3)-disp(i,3)) .lt. tol) then
C$$$                i2 = iax(2,iat)
C$$$                dpos(1,i2) = dpos(1,i2) + disp(i,4)
C$$$                dpos(2,i2) = dpos(2,i2) + disp(i,5)
C$$$                dpos(3,i2) = dpos(3,i2) + disp(i,6)
C$$$                 if (ipr .ge. 50) write(stdo,333)
C$$$     .            i1,i2,(rtab(j,iat),j=1,3),(disp(i,j),j=4,6)
C$$$  333           format(2i4,3f11.6,2x,3f11.6)
C$$$                count = count+1
C$$$              endif
C$$$            endif
C$$$          endif
C$$$        enddo

C$$$   10   continue
C$$$      enddo
C$$$      call info(30,0,0,' ... Added displacements for %i site(s),'//
C$$$     .  ' %i pair(s)',nlst1,count)
C$$$      end
C$$$      subroutine psup5(nn,disp,rmax)
C$$$C- Determine largest length of position vectors
C$$$C     implicit none
C$$$      integer nn
C$$$      double precision disp(nn,3),rmax,ddot
C$$$      integer i

C$$$      rmax = 0
C$$$      do  i = 1, nn
C$$$        rmax = max(rmax,ddot(3,disp(i,1),nn,disp(i,1),nn))
C$$$      enddo
C$$$      rmax = dsqrt(rmax)
C$$$      end

