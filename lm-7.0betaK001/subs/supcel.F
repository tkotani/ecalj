      subroutine supcel(mode,sctrl,slat,ssite,slabl,sspec,plx,nkd,dlat,
     .  nbx)
C- Generates basis vectors in a supercell
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :not used now.  Should be zero
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nspec
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat
Ci     Stored:
Ci     Passed to:
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    pos spec pl
Ci     Passed to: spackv usite dval dpscop spack dvset dvprm
Ci   slabl :vector of species labels
Ci     Elts read:
Ci     Stored:
Ci     Passed to: baslst iosits
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z
Ci     Stored:    name
Ci     Passed to: spacks iosits
Ci   plx   :supercell primitive lattice vectors
Ci   nkd   :number of direct lattice vectors (lattic.f)
Ci   dlat  :direct lattice vectors
Co Outputs
Co   nbx   :number of sites in superlattice
Cl Local variables
Cl   qlx   :supercell reciprocal lattice
Cr Remarks
Cb Bugs
Cb   No check is made to see whether ther are enough dlat to encompass
Cb   all lattice vectors
Cu Updates
Cu   01 Aug 08 New --wsitex
Cu   19 Apr 03 Changed distortion to read from displacements file
Cu   12 Apr 03 Some bug fixes when partial list of sites saved
Cu   09 Apr 03 Added pairwise distortion; also writes original cell
Cu             site sdata first
Cu   24 Nov 02 R.S. pairwise hamiltonian may be created for supercell
Cu             from file rsj
Cu   26 Nov 02 supcel can restrict supercell sites to a specified list
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*8 slabl(1)
      integer nkd,mode,nbx
      double precision sctrl(1),slat(1),ssite(1),sspec(1),dlat(3,nkd)
C ... Local parameters
      character*72 alabl*8,outs,sortex(5),out2,jfilnm,dc*1
      logical cmdopt,lsort,lfirst,latvec
      integer lgunit,i,j,k,m,ib,nbas,is,lok,ndsite,iosits,lio,ifi,nspec,
     .  ixv(36),imark(36),it(3),ip,parg,j1,j2,nexpr,nbxsav,a2vec,nmark,
     .  fopna,fxst,nfilj,nfilm,nw,stdo,fopn,rdm,nn,nc,iprint
      integer mxcsiz,nbasp,ifirsj,iorsj,nttab,nlst,ilst,nlst1
      integer ontab,oiax,ortab,mxcsz,ohpair,oaamom,oimap,oamom2,oips,oz,
     .  oilst,olst1,odisp,opos2,owk,oiprm,ossite,opos,odpos
      double precision sclrsj
      double precision plat(3,3),plx(3,3),qlx(3,3),qlat(3,3)
      double precision dd(3),tau(3),x(3),vol,volx,dval,tol,alat,rcut,xx
C ... Heap
      integer w(1)
      common /w/ w

      call upack1('lat alat',slat,alat)
      call upack1('lat plat',slat,plat)
      call upack1('ctrl nbas',sctrl,nbas)
      call upack1('ctrl nspec',sctrl,nspec)
      tol = 1d-6
      stdo = lgunit(1)

C --- Get list of lattice vectors to include ---
      nmark = 0
      if (cmdopt('-use:',5,0,outs)) then
        m = 5
        k = a2vec(outs,72,m,2,',; ',3,-3,36,ixv,imark)
        nmark = k/3
        if (3*nmark .ne. k) call rxi('SUPCEL: ''-use:'' takes'//
     .    ' multiples of 3 arguments, but read',k)

        do  10  ib = 1, nmark
C         Assemble vector as specified multiples of prim. latt. vec.
          call dpzero(x,3)
          call daxpy(3,dble(imark(3*ib-2)),plat(1,1),1,x,1)
          call daxpy(3,dble(imark(3*ib-1)),plat(1,2),1,x,1)
          call daxpy(3,dble(imark(3*ib-0)),plat(1,3),1,x,1)
C         Mark which of dlat it belongs to
          do  12  k = 1, nkd
            if (tol .gt. abs(x(1)-dlat(1,k))+
     .                   abs(x(2)-dlat(2,k))+
     .                   abs(x(3)-dlat(3,k))) then
              imark(ib) = k
              goto 14
            endif
   12   continue
        call rx('supcel: missing dlat')
   14   continue
   10   continue
      endif

C --- Optionally get a list (subset) of sites in original basis ---
      nlst = 0
      call defi(oilst,-nbas)
      if (cmdopt('--sites',7,0,outs)) then
        i = 8
        call sp2cls('site spec',ssite,0,1,1,nbas,oips)
        call sp2cls('spec z',sspec,0,1,1,nspec,oz)
C       call awrit2('%N Z= %n:1,3;3d',' ',100,6,nspec,w(oz))
        dc = outs(i:i)
        call baslst(0,11,outs(i:),j2,w(oips),nbas,slabl,w(oz),0,' ',x,
     .    nlst,w(oilst))
        if (nlst .eq. 0) w(oilst) = 0
      endif

C --- Supercell reciprocal lattice vectors ---
      write(stdo,'(/a)') ' Basis vectors of supercell:'
      call dinv33(plx,1,qlx,volx)
      call dinv33(plat,0,qlat,vol)
      if (dabs(volx) .lt. 1d-6)
     .   call rx('supercell vectors are linearly dependent')
      write(stdo,351)
  351 format(/15x,'plx',30x,'qlx')
      write(stdo,350) ((plx(m,k),m=1,3),(qlx(m,k),m=1,3),k=1,3)
  350 format(3f10.5,5x,3f10.5)

      call ivset(ixv,1,3,2)
      if (cmdopt('-shorten=',9,0,outs) .or.
     .    cmdopt('--shorten=',10,0,outs)) then
        ip = 9
        if (cmdopt('--shorten=',10,0,outs)) ip = 10
        call skipbl(outs,len(outs),ip)
        if (parg(' ',2,outs,ip,len(outs),', ',2,3,it,ixv) .lt. 0)
     .    call rxs2('supcel: failed to parse "',outs(1:ip+5),' ..."')
        print '(/''  supcel: shorten vectors with mode ='',3i2)',ixv
      endif

C --- Count supercell basis vectors (nbx) ---
C     and make table mapping supercell sites to original sites
      write (stdo,*) ' '
      nbx = 0
      ilst = 0
      do  20  ib = 1, nbas
        if (nlst .gt. 0) then
   11     if (ilst .ge. nlst) goto 20
          if (w(oilst+ilst) .lt. ib) then
            ilst = ilst+1
            goto 11
          endif
          if (w(oilst+ilst) .ne. ib) goto 20
        endif

        call upack('site spec pos',ssite,ib,is,tau,0,0)
        call spacks(0,'spec name',sspec,alabl,is,is)
        do  22  k = 1, nkd
          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          if (lok .eq. 1) write(stdo,120) ib,alabl,(dd(m),m=1,3)
  120     format(i6,2x,a,3f11.6)
   22   continue
   20 continue
C     Make table imap mapping supercell sites to original sites
      call defi(oimap,-nbx)
      nbx = 0
      ilst = 0
      do  24  ib = 1, nbas
        if (nlst .gt. 0) then
   21     if (ilst .ge. nlst) goto 24
          if (w(oilst+ilst) .lt. ib) then
            ilst = ilst+1
            goto 21
          endif
          if (w(oilst+ilst) .ne. ib) goto 24
        endif
        call upack('site spec pos',ssite,ib,is,tau,0,0)
        call spacks(0,'spec name',sspec,alabl,is,is)
        do  26  k = 1, nkd
          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          if (lok .eq. 1) w(oimap+nbx-1) = ib
   26   continue
   24 continue

C --- Printout ---
      write(stdo, 333) abs(vol), abs(volx), abs(volx/vol)
  333 format(/' cell vol=',f10.6,'     scell vol=',f12.6,
     .        '     ratio=',f7.2)
      write(stdo, 334) nbas, nbx, dble(nbx)/nbas
  334 format( ' cell nbas=',i9,'     scell nbas=',i11,
     .        '     ratio=',f7.2)
      if (nlst .eq. 0) then
        if (abs(dble(nbx)/nbas-abs(volx/vol)) .gt. 1d-5)
     .    write(stdo, 335)
      else
        if (abs(dble(nbx)/nlst-abs(volx/vol)) .gt. 1d-5)
     .    write(stdo, 335)
      endif
  335 format(' (warning) supercell has incommensurate lattice')

      lsort = cmdopt('-sort:',6,0,outs) .or. cmdopt('--sort:',7,0,outs)
      lfirst= cmdopt('-first',6,0,outs) .or. cmdopt('--first',7,0,outs)

C --- Set up expressions to sort ---
      i = 7
      if (cmdopt('--sort:',7,0,outs)) i = 8
      nexpr = 0
      do  48  k = 1, 3
        call word(outs(i:),k,j1,j2)
        if (j2 .ge. j1) then
          nexpr = nexpr + 1
          sortex(nexpr) = outs(6+j1:6+j2)
        else
          goto 49
        endif
   48 continue
   49 continue

      if (.not. (cmdopt('-wsite=',7,0,outs) .or.
     .           cmdopt('--wsite=',8,0,outs) .or.
     .           cmdopt('--wsitex=',9,0,outs))) return

C --- Get list of positions; sort them ---
      if (lsort .or. lfirst) then
      call defrr(opos,3*nbx)
      call defi(owk,-nbx)
      call defi(oiprm,nbx)
      nbxsav = nbx
      nbx = 0
      ilst = 0
      do  40  ib = 1, nbas
        if (nlst .gt. 0) then
   41     if (ilst .ge. nlst) goto 40
          if (w(oilst+ilst) .lt. ib) then
            ilst = ilst+1
            goto 41
          endif
          if (w(oilst+ilst) .ne. ib) goto 40
        endif

        call upack('site spec pos',ssite,ib,is,tau,0,0)
        call spacks(0,'spec name',sspec,alabl,is,is)
        do  42  k = 1, nkd
          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          if (lok .eq. 1) then
            call dmscop(w(opos),nbxsav,dd,1,1,1,1,3,nbx,1,1d0)
            w(owk+nbx-1) = nbx
C           case lfirst: if dd = tau, this was the original vector
            if (lfirst) then
              xx = (dd(1)-tau(1))**2+(dd(2)-tau(2))**2+(dd(3)-tau(3))**2
              x(1) = (dd(1)-tau(1))
              x(2) = (dd(2)-tau(2))
              x(3) = (dd(3)-tau(3))
              if (latvec(1,tol,qlx,x)) then
                w(owk+nbx-1) = ib-nbas
              endif
            endif
          endif
   42   continue
   40 continue

C     call prmx('unsorted pos',w(opos),nbx,nbx,3)

      if (lsort) then
        call defrr(opos2,3*nbx)
        call dcopy(3*nbx,w(opos),1,w(opos2),1)
        call mapdat(nexpr,sortex,' ',w(oiprm),nbx,3,w(opos2),w)
        call defrr(owk,nexpr*nbx)
        call dmcpy(w(opos2),nbx,1,w(owk),1,nexpr,nbx,nexpr)
C       call prmx('expressions',w(owk),nexpr,nexpr,nbx)
        call dcopy(3*nbx,w(opos),1,w(opos2),1)
        call dvheap(nexpr,nbx,w(owk),w(oiprm),0d0,1)
C       call msort(nbx,3,nexpr,w(oiprm),w(owk),w(opos2),w(opos))
C       call prmx('sorted pos',w(opos),nbx,nbx,3)
      else
        call ivheap(1,nbx,w(owk),w(oiprm),1)
C       call yprm('ordering',0,w(owk),0,nbx,nbx,1)
C       call yprm('iprm',0,w(oiprm),0,nbx,nbx,1)
      endif

      endif

      if (.not. (cmdopt('-wsite=',7,0,outs) .or.
     .           cmdopt('--wsite=',8,0,outs) .or.
     .           cmdopt('--wsitex=',9,0,outs))) return


C --- Read R.S. hamiltonian from file rsj ---
C     ifirsj is the file logical unit for real-space hamiltonian
C     ifirsj = 0 flags that no file will be written
C     Switch syntax:  --wrsj[:fn=name][:scl=#]
      ifirsj = 0
      if (cmdopt('--wrsj',6,0,jfilnm) .or.
     .    cmdopt('-wrsj',5,0,jfilnm)) then

        out2 = jfilnm(7:)
        if (jfilnm(1:5) .eq. '-wrsj') out2 = jfilnm(6:)
C       This is default name for supercell rsj file
        jfilnm = 'rsj2'

        dc = out2(1:1)
        sclrsj = 1
        if (dc .ne. ' ') then
C     ... Return here to resume parsing for arguments
          j2 = 0
   50     continue
          j2 = j2+1
          if (out2(j2:j2) .eq. dc) goto 50
          j1 = min(len(out2),j2)
          call nwordg(out2,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (out2(j1:j1+3) .eq. 'scl=')  then
              j = 0
              i = parg('scl=',4,out2(j1:),j,len(out2(j1:)),
     .          dc//' ',1,1,k,sclrsj)
              if (i .le. 0) goto 59

            elseif (out2(j1:j1+2) .eq. 'fn=')  then
              jfilnm = out2(j1+5:j2)
            else
              goto 59
            endif
            goto 50
   59       continue
            call rxs('supcel: failed to parse wrsj option', out2)
          endif
        endif

        call info0(20,1,0,
     .    ' ... Writing supercell r.s. hamiltonian to file '//jfilnm)

C   ... Setup call for file read
        mxcsiz = 0
        nbasp = nbas
        if (fxst('rsj') .ne. 1) then
          call rx('supcel: specified --wrsj but no file rsj')
        endif
        ifirsj = fopna('rsj',-1,0)
        call defrr(oaamom,nbas)
        nfilj = iorsj(0,ifirsj,nbas,alat,plat,w(oaamom),w,w,w,w,x)

C       This makes iax table and rtab suitable for iorsj
        call defrr(opos,3*nbas)
        call spackv(10,'site pos',ssite,1,nbas,w(opos))
        call pairs(nbas,nbasp,1d0,plat,x(1)/2*(1+1d-6),w(opos),
     .    -1,3,-1,w,nttab,ontab,oiax,mxcsiz)
        call defrr(ortab,3*nttab)
        call mkrtab(000,1d0,plat,w(opos),w(oiax),nttab,w(opos),w(ortab))

C       Read pair hamiltonian for entries matching iax table
        call defrr(ohpair,nttab)
        nfilm = iorsj(03,ifirsj,nbas,w,plat,w(oaamom),w(oiax),w(ortab),
     .    w(ontab),w(ohpair),x)

        if (nfilj.ne.nfilm) call info2(20,0,0,'%8f(warning) only'//
     .    ' matched %i out of %i file entries',nfilm,nfilj)
        if (nfilm.ne.nttab) then
          call info2(20,0,0,'%8f(warning) iax table has %i pairs'//
     .      ' but matched only %i file entries',nttab,nfilm)

C         Purge iax table of unused entries
          call psup3(nlst,w(oilst),w(oiax),nttab,nfilm,w(ortab))
          call symiax(1,plat,nbas,w,w,w,0,w(ontab),w(oiax),nttab,mxcsz)
          if (nttab .ne. nfilm) then
            call fexit2(-1,111,' Exit -1 : iorsj: bug in symiax:'//
     .        ' nttab=%i, nfilm=%i',nttab,nfilm)
          endif

C         Remake rtab
          call mkrtab(000,1d0,plat,w(opos),w(oiax),nttab,w(opos),
     .      w(ortab))

C         Re-read pair hamiltonian with new iax table
          nfilm = iorsj(03,ifirsj,nbas,w,plat,w(oaamom),w(oiax),
     .      w(ortab),w(ontab),w(ohpair),x)
          if (nttab .ne. nfilm) call rx('supcel: bug in iorsj')

        endif

C       Close input file rsj
        call fclose(ifirsj)

C       Open supercell file; write header
        call defrr(oamom2,nbx)
        call cpvprm(0,1,nbx,w(oimap),w(oaamom),w(oamom2))
        ifirsj = fopna(jfilnm,-1,0)
        i = iorsj(0,-ifirsj,nbx,alat,plx,w(oamom2),w,w,w,w,w)

      endif

C --- Create supercell ssite ---
      call defrr(ossite,1000)
      call usite(w(ossite),-1,0,1,0,0,0)
      ndsite = nint(dval(w(ossite),1))
      call redfrr(ossite,nbx*ndsite)
      call usite(w(ossite),-1,0,nbx,0,0,0)
      nspec = 0
      nbx = 0

      ilst = 0
      do  30  ib = 1, nbas
        if (nlst .gt. 0) then
   31     if (ilst .ge. nlst) goto 30
          if (w(oilst+ilst) .lt. ib) then
            ilst = ilst+1
            goto 31
          endif
          if (w(oilst+ilst) .ne. ib) goto 30
        endif

        call upack('site spec pos',ssite,ib,is,tau,0,0)
        call spacks(0,'spec name',sspec,alabl,is,is)
        nspec = max(nspec,is)
        do  32  k = 1, nkd
          call psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
          if (lok .eq. 1) then

C       ... Copy parameters from site ib to site nbx
            call dpscop(ssite,w(ossite),ndsite,
     .        1+(ib-1)*ndsite,1+(nbx-1)*ndsite,1d0)
C           12 Apr 03 ssite(1) may not be struc size if 1st ib<>1
            call dvset(w(ossite),1,1,ssite)
            call pack5('site spec pos',w(ossite),nbx,is,dd,0,0)

C       ... Write rsj to supercell file
            if (ifirsj .ne. 0) then
            call psup2(ifirsj,ib,nbx,w(ontab),w(ortab),sclrsj,w(ohpair))
            endif

C       ... PL index
            if (cmdopt('--pl:',5,0,out2)) then
              call prjpos(2,ixv,w,plat,1,dd,x)
              call mapdat(1,out2(6:),' ',m,1,3,x,w)
              call spack(2,'site pl',ssite,1,m,j1,j2,k)
              call dvset(w(ossite),m+(nbx-1)*ndsite,m+(nbx-1)*ndsite,x)
            endif
          endif
   32   continue
   30 continue

C --- Permute site strux if sorted ---
      if (lsort .or. lfirst) then
        call defrr(owk,nbx*ndsite)
        call dvprm(ndsite,nbx,w(ossite),w(owk),w(oiprm),.true.)
        call dvset(w(ossite),1,1,dble(ndsite))
        call rlse(owk)
      endif

C --- Displacement around an atom ---
      call defrr(opos,3*nbx)
      call spackv(10,'site pos',w(ossite),1,nbx,w(opos))
      call defrr(odpos,-3*nbx)
      if (cmdopt('-disp',5,0,out2) .or.
     .    cmdopt('--disp',6,0,out2)) then
        dc = out2(6:6)
        if (out2(1:6) .eq. '--disp') dc = out2(7:7)
        call wrdsg(out2,1,dc//' ',nw)
        if (nw .lt. 3) call rx('supcel: --disp requires 2 arguments')

C   ... Read displacements file
        call wordg(out2,0,dc//' ',2,j1,j2)
        ifi = fopn(out2(j1:j2))
        rewind ifi
C       Determine size of file
        nn = 0
        nc = 0
        j = rdm(ifi,0,0,' ',xx,nn,nc)
        if (j .lt. 0) call rxs
     .    ('supcel: failed to read displacements file ',out2(j1:j2))
        if (nc .ne. 6) call rxs('supcel: displacements file '//
     .    'apparently does not have 6 columns: ',out2(j1:j2))
        call defrr(odisp,nn*nc)
        rewind ifi
        j = rdm(ifi,0,nn*nc,' ',w(odisp),nn,nc)
        call psup5(nn,w(odisp),rcut)
        call wordg(out2,0,dc//' ',3,i,j)
        call info(30,1,0,' ... Displacements from file '//
     .    out2(j1:j2)//' for site list='//out2(i:j)//'. rcut=%d',rcut,0)

        nlst1 = 0
        call defi(olst1,-nbx)
        call wordg(out2,0,dc//' ',3,j1,j2)
        if (j2 .lt. j1) call rx('supcel: --disp requires 2 arguments')
        call sp2cls('site spec',w(ossite),0,1,1,nbx,oips)
        call sp2cls('spec z',sspec,0,1,1,nspec,oz)
        call baslst(1,11,out2(j1-1:),j2,w(oips),nbx,slabl,w(oz),0,' ',
     .    x,nlst1,w(olst1))

C       This makes iax table and rtab
        mxcsiz = 0
C       call setpr(50)
        call pshpr(iprint()-10)
        call pairs(nbx,nbx,1d0,plx,rcut/2*(1+1d-6),w(opos),
     .    -1,3,-1,w,nttab,ontab,oiax,mxcsiz)
        call defrr(ortab,3*nttab)
        call mkrtab(000,1d0,plx,w(opos),w(oiax),nttab,w(opos),w(ortab))
        call poppr
        call ivset(ixv,1,3,2)
        call shorps(nttab,plx,ixv,w(ortab),w(ortab))

C       Get displacements
        call psup4(w(oiax),nttab,
     .    nlst1,w(olst1),w(ortab),w(opos),w(odisp),nn,w(odpos))

C       Add to ssite
        call daxpy(3*nbx,1d0,w(odpos),1,w(opos),1)
        call spackv(11,'site pos',w(ossite),1,nbx,w(opos))
      endif

C --- Write site file ---
      lio = 1000*(2+4+8+16+32) + 1
C     lio = lio+16000
      if (outs(1:7) .eq. '-wsite=') then
        out2 = outs(8:)
        outs = out2
      elseif (outs(1:8) .eq. '--wsite=') then
        out2 = outs(9:)
        outs = out2
      elseif (outs(1:9) .eq. '--wsitex=') then
        out2 = outs(10:)
        outs = out2
        lio = lio + 10
      endif
      if (iosits(lio,3d0,0,outs,ifi,slabl,alat,plx,nbx,nspec,sspec,
     .  w(ossite)).lt. 0) call rx('failed to write ssite')

C --- Write positions to file ---
      if (cmdopt('--wpos=',7,0,outs)) then
        call defrr(opos,3*nbx)
        call spackv(10,'site pos',w(ossite),1,nbx,w(opos))
        call iopos(.true.,0,outs(8:),nbx,w(opos))
      endif


      return
C 999 continue

      end

      subroutine psupcl(nmark,imark,tau,k,dlat,plx,qlx,ixv,lok,nbx,dd)
C- Kernel called by supcel: to decide whether a vector should be added
C ----------------------------------------------------------------------
Ci Inputs
Ci   nmark :?
Ci   imark :?
Ci   tau   :position vector for this site
Ci   k     :which of the lattice vectors in dlat to check
Ci   dlat  :a list of lattice vectors, of which dlat(k) is to be used here.
Ci   plx   :(passed to shorps; used if vectors are to be shortened)
Ci   qlx   :reciprocal lattice vectors to supercell
Ci   ixv   :(passed to shorps; used if vectors are to be shortened)
Ci         :0 suppresses shifts along plat(j)
Ci         :1 shifts to unit cell at origin (pos in 1st quadrant)
Ci         :2 shifts to minimize length of pos
Co Outputs
Co   lok   :0 do not add this vector to the list
Co         :1 add this vector to the list
Co   dd    :connecting vector
Co   nbx   :incremented by one if this vector should be included
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nmark,imark(nmark),lok,nbx,ixv(3),k
      double precision dlat(3,k),tau(3),plx(3,3),qlx(3,3),dd(3)
C ... Local parameters
      double precision x(3)
      integer m
      logical cmdopt
      character*72 outs

      dd(1) = dlat(1,k)+tau(1)
      dd(2) = dlat(2,k)+tau(2)
      dd(3) = dlat(3,k)+tau(3)
      if (nmark .gt. 0) then
        lok = 0
        do  10  m = 1, nmark
   10   if (imark(m) .eq. k) lok = 1
      else
        lok = 1
        do  20  m = 1, 3
          x(m) = dd(1)*qlx(1,m)+dd(2)*qlx(2,m)+dd(3)*qlx(3,m)
          if (x(m).lt.-0.001d0 .or. x(m).gt.0.999d0) lok = 0
   20   continue
      endif
      if (lok .eq. 1) then
        nbx = nbx+1
        if (cmdopt('-shorten',8,0,outs)) call shorps(1,plx,ixv,dd,dd)
        if (cmdopt('--shorten',9,0,outs)) call shorps(1,plx,ixv,dd,dd)
C 120   format(2x,a,3f11.6)
      endif
      end

      subroutine psup2(ifi,ib,ib2,ntab,rtab,sclrsj,rsj)
C- Kernel called by supcel: write r.s. pairwise ham for this site
C  Note that jb2 of supercell is not known, so we just write zero
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   ib    :site index of original cell
Ci   ib2   :site index of supercell
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   rsj
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Nov 02
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,ib,ib2,ntab(ib+1)
      double precision sclrsj
C     integer niax
C     parameter (niax=10)
C     integer iax(niax,*)
      double precision rsJ(*),rtab(3,*)
C ... Local parameters
      integer it

      do  it = ntab(ib)+1, ntab(ib+1)
C        call awrit6('%,4i%,4i%;12,7D%;12,7D%;12,7D%;12,4D',' ',80,ifi,
C     .    ib2,0,rtab(1,it),rtab(2,it),rtab(3,it),sclrsj*rsj(it))
        call awrit6('%,6i%,6i%;12,7D%;12,7D%;12,7D %;6g',' ',80,ifi,
     .    ib2,0,rtab(1,it),rtab(2,it),rtab(3,it),sclrsj*rsj(it))
      enddo
      end

      subroutine psup3(nlst,ilst,iax,nttab,n,rtab)
C- Kernel called by supcel to purge pairs outside of site list
C  used when a subset of all sites is specified
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlst  :number of sites in site list
Ci   ilst  :site list
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Co Outputs
Co   n     :number of pairs after purging
Co         :NB: n is not set if nlst=0
Cr Remarks
Cr
Cu Updates
Cu   26 Nov 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax,nlst,ilst(nlst),nttab
      parameter (niax=10)
      integer iax(niax,nttab)
      double precision rtab(3,nttab)
C ... Local parameters
      integer i,ib,jb,low,n

      if (nlst .eq. 0) return
      n = 0
      do  i = 1, nttab
        ib = iax(1,i)
        low = 1
        call hunti(ilst,nlst,ib,0,low)
        if (low .ge. nlst) then
C         print *, 'i, iax',i,iax(1,i)
          iax(1,i) = 0
        elseif (ilst(low+1) .ne. ib) then
C         print *, 'i, iax',i,iax(1,i)
          iax(1,i) = 0
        else
          jb = iax(2,i)
          call hunti(ilst,nlst,jb,0,low)
          if (low .ge. nlst) then
C           print *, 'i, iax',i,iax(1,i),iax(2,i)
            iax(1,i) = 0
          elseif (ilst(low+1) .ne. jb) then
C           print *, 'i, iax',i,iax(1,i),iax(2,i)
            iax(1,i) = 0
          else
            n = n+1
          endif
        endif
      enddo

      call info2(20,0,0,' supcel: site-list further reduced pair'//
     .  ' table to %i pairs',n,0)

      end

      subroutine mapdat(nexpr,expr,sincl,iwk,nr,nc,dat,dat2)
C- Replace columns of dat with algebraic expressions of them
C     implicit none
C Passed parameters
      integer nexpr,nr,nc,iwk(nr)
      character*(*) expr(1), sincl
      double precision dat(nr,nc),dat2(*)
C Local variables
      integer ir,iv0,ival,i,j,ii,jr
      character*4 xn, outs*80
      logical a2bin,logi

      call numsyv(iv0)
      jr = 0
      do  20  ir = 1, nr
        call clrsyv(iv0)
        call addsyv('i',dble(ir),ival)

C   ... Load data table
        do  22  j = 1, nc
          ii = 1
          xn = 'x   '
          call bin2a('(i2)',0,0,j,2,0,4,xn,ii)
          call addsyv(xn,dat(ir,j),ival)
   22   continue

C   ... Exclude points if not satisfy sincl
        logi = .true.
        if (sincl .ne. ' ') then
          j = 0
          if (.not. a2bin(sincl,logi,0,0,' ',j,-1))
     .      call rx('mapdat:  error parsing sincl')
        endif

C   ... Put expressions of these vars into dat
        do  24  i = 1, nexpr
          j = 0
          if (.not. a2bin(expr(i),dat(ir,i),4,0,' ',j,-1)) then
            outs = expr(i)
            call skpblb(expr(i),len(expr(i)),ii)
            call fexit(-1,1,'MAPDAT: failed to parse expr: '
     .        //outs(1:j+1)//'<-->'//outs(j+2:ii+1),0)
          endif
   24   continue
        if (logi) then
          jr = jr+1
          iwk(jr) = ir
        endif
   20 continue
      if (jr .ne. nr) then
        call dcopy(nr*nexpr,dat,1,dat2,1)
        call rowpmt(nr,jr,nexpr,iwk,dat2,dat)
        nr = jr
      endif
      call clrsyv(iv0)

      end

      subroutine msort(nr,ncs,ncw,iprm,wk,s,ssort)
C- sort matrix according to increasing wk.
C     implicit none
      integer nr,ncs,ncw,iprm(nr)
      double precision wk(ncw,nr),s(nr,ncs),ssort(nr,ncs)
      integer ir,ic,ix

      call dvheap(ncw,nr,wk,iprm,0d0,1)
      do  20  ir = 1, nr
        ix = iprm(ir)
        do  30  ic = 1, ncs
   30   ssort(ir,ic) = s(ix,ic)
   20 continue
      end

      subroutine rowpmt(nf,nt,nc,ipmt,afrom,ato)
C     implicit none
      integer nf,nt,nc,ipmt(nt),i
      double precision afrom(nf,nc), ato(nt,nc)

      do  10  i = 1, nt
   10 call dcopy(nc,afrom(ipmt(i),1),nf,ato(i,1),nt)
      end

      subroutine psup4(iax,nttab,nlst1,lst1,rtab,pos,disp,nn,dpos)
C- Kernel called by supcel to displace pairs radially
C ----------------------------------------------------------------------
Ci Inputs
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   nlst1 :number of sites of "center" type
Ci   lst1  :list of sites of "center" type
Ci   nlst2 :number of sites of "neighbor" type
Ci   lst2  :list of sites of "neighbor" type
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Co Inputs/Outputs
Cio   pos  :On input, basis vectors
Cio        :On output, displaced basis vectors
Cu Updates
Cu   07 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlst1,lst1(nlst1),niax,nttab,nn
      double precision pos(3,*),rtab(3,*),dpos(3,*),disp(nn,6)
      parameter (niax=10)
      integer iax(niax,1)
C ... Local parameters
      integer iat,low1,i1,i2,ipr,count,stdo,lgunit,i,j
      double precision rhat(3),tol

      count = 0
      low1 = 0
      call getpr(ipr)
      ipr = 100
      stdo = lgunit(1)
      tol = 1d-6

C      if (ipr .ge. 50) write(stdo,330) nlst1
C  330 format(' supcel: add displacement around',i4,' site(s)'/
C     .  '  ia  ib',15x,'rb-ra',27x,'delta pos(ib)')
      if (ipr .ge. 50) write(stdo,330)
  330 format('  ia  ib',15x,'rb-ra',27x,'delta pos(ib)')

      do  iat = 1, nttab

C   ... If the center atom isn't in the list, skip this pair
        i1 = iax(1,iat)
        call hunti(lst1,nlst1,i1,0,low1)
        if (i1 .ne. lst1(low1+1)) goto 10

C   ... See if this pair is in any of in the disp list
        call dcopy(3,rtab(1,iat),1,rhat,1)
        do  i = 1, nn
          if (abs(rhat(1)-disp(i,1)) .lt. tol) then
            if (abs(rhat(2)-disp(i,2)) .lt. tol) then
C             We have a match!  Add to dpos
              if (abs(rhat(3)-disp(i,3)) .lt. tol) then
                i2 = iax(2,iat)
                dpos(1,i2) = dpos(1,i2) + disp(i,4)
                dpos(2,i2) = dpos(2,i2) + disp(i,5)
                dpos(3,i2) = dpos(3,i2) + disp(i,6)
                 if (ipr .ge. 50) write(stdo,333)
     .            i1,i2,(rtab(j,iat),j=1,3),(disp(i,j),j=4,6)
  333           format(2i4,3f11.6,2x,3f11.6)
                count = count+1
              endif
            endif
          endif
        enddo

   10   continue
      enddo
      call info(30,0,0,' ... Added displacements for %i site(s),'//
     .  ' %i pair(s)',nlst1,count)
      end
      subroutine psup5(nn,disp,rmax)
C- Determine largest length of position vectors
C     implicit none
      integer nn
      double precision disp(nn,3),rmax,ddot
      integer i

      rmax = 0
      do  i = 1, nn
        rmax = max(rmax,ddot(3,disp(i,1),nn,disp(i,1),nn))
      enddo
      rmax = dsqrt(rmax)
      end

