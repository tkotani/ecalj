      subroutine lmaux(prgnam,sctrl,sham,spot,slat,sspec,ssite,sstr,
     .  sarray,sbz,slabl,sstrn,mode,wksize)
C- Auxilliary programs
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nclass nl nspec nspin lncol modep lpgf lrs,7 
Ci                omax1 omax2 wsrmax
Ci     Stored:    *
Ci     Passed to: asars supot asamad asvsph aioxtn
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oeula neula qss
Ci     Stored:    *
Ci     Passed to: asvsph
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opnu oqnu opp oqc oqt oaamom orhrmx oves ovrmax
Ci     Stored:    *
Ci     Passed to: asars supot asamad asvsph
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: avw alat plat plat2 nkd odlv nkq oqlv
Ci     Stored:    *
Ci     Passed to: asars supot asamad plana asvsph aioxtn
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: idmod lmxa z rham mxcst rmt a nr
Ci     Stored:    z rmt
Ci     Passed to: asars spackv asamad plana pp2enu asvsph shopol aioxtn
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    pos spec
Ci     Passed to: asars iosits spackv plana aioxtn
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci     Elts read: mxnbr rmax
Ci     Duplicate: mxnbr
Ci     Stored:    *
Ci     Passed to: *
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: ohave oics oclabl onrc npadl npadr nclasp opgfsl
Ci                ormax opos oipc oips
Ci     Stored:    ormax
Ci     Passed to: asamad asvsph aioxtn
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: asars
Ci   slabl :vector of species labels
Ci   sstrn :struct for global strings
Ci   mode  :a compound of bits, which are independent of each other
Ci         :  2**0 Show neighbors
Ci         :  2**1 Plane analysis
Ci         :  2**2 Generate input to xbs program
Ci         :  2**3 Shift moments, pp's to new linearization energy
Ci         :  2**4 Interpolate core to another mesh
Ci         :  2**5 Display poles of potential functions
Ci         :  2**6 Import data from other formats
Ci         :  2**7 Find empty spheres
Cu Updates
Cu   12 Aug 08 (L. Ke) empty sphere finder
Cu   04 Nov 04 Upgrade of rsta editor
Cu   26 Jan 03 Call to angtab changed
Cu   17 May 02 Modified MT radii scaling to lower priority for E.S.
Cu   23 Apr 02 Added option (--getwsr) to find MT radii
Cu   01 Mar 02 Updated Import data mode
Cu   05 Oct 01 Adapted mode 2**3 to work with lm v6.11
Cu   24 Nov 97 changed ovmin to run quickly
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,wksize
      character*(*) sstrn, prgnam*8
      character*8 slabl(*)
      double precision sarray(1),sctrl(1),sham(1),spot(1),slat(1),
     .  sspec(1),ssite(1),sstr(1),sbz(1)
C ... Local parameters
      character*120 outs,fnam*8
      integer NULLI
      logical cmdopt,T,F,swtmp
      parameter (T=.true., F=.false., NULLI=-99999)
      integer fopn,getdig,i,ip,j,k,m,ifi,iprint,lpgf,lpbc,lgunit,lncol,
     .  nbas,nclasp,nclass,nl,nlspc,nsp,modep(3),mxnbr,parg,nbasp,
     .  nbaspp,nkd,nkq,nspec,npadl,npadr,neul,nc,mxcsiz,nttab,igets,
     .  iosits,cmplat,ngrp,ival,stdo,irs(5),fopna,nclspp,bitand,igetss,
     .  ngmx,nsgrp
      integer oeold,oeula,oics,oidmod,oinitc,oipc,olmx,onrc,opnu,opold,
     .  opp,oqc,oqnu,oqt,orhrmx,ormax,oves,ovrmax,owk,owk2,oclabl,
     .  opos,opos2,oqlv,odlv,orham,oamsh,onrmsh,oalpha,oaamom,
     .  onpr,ormt,olock,olockc,os,oz,ozz,ormx,oip,oips,oips2,opgfsl,
     .  ontab,oiax,og,ossite,oag,ozc,mxclas,ormtc
      double precision xv(10),xx,alat,plat(3,3),plat2(9),facrmx,facrng,
     .  dval,avw,ekap(2),enu,rmaxs,qss(4),ckbas,cksumf,ehterm(4),
     .  qlat(9),emad,trumad,vmtz(2),omax1(3),omax2(3),wsrmax
      parameter (ngmx=48,mxclas=1000)
C ... Heap
      integer w(1)
      common /w/ w

      call upack('ctrl nbas nclass nl nspec nspin',sctrl,nbas,nclass,nl,
     .  nspec,nsp)
      call upack('ctrl lncol modep',sctrl,lncol,modep,0,0,0)
      lpgf = igets('ctrl lpgf',sctrl)
C     lpbc = 0 for pbc in 3 dimensions, 11 pgf padded geometry
      lpbc = 0
      if (lpgf .gt. 0 .and. lpgf .lt. 10) lpbc = 11
      call upack('array ohave oics oclabl onrc',sarray,oinitc,oics,
     .  oclabl,onrc,0)
      call upack('array npadl npadr nclasp opgfsl ormax',sarray,npadl,
     .  npadr,nclasp,opgfsl,ormax)
      call upack('array opos oipc oips',sarray,opos,oipc,oips,0,0)
      call upack('ham oeula neula qss',sham,oeula,neul,qss,0,0)
      call upack('pot opnu oqnu opp oqc oqt',spot,opnu,oqnu,opp,oqc,oqt)
      call upack('pot oaamom orhrmx oves ovrmax',spot,oaamom,orhrmx,
     .  oves,ovrmax,0)
      call upack('lat avw alat plat plat2',slat,avw,alat,plat,plat2,0)
      call upack('lat nkd odlv nkq oqlv',slat,nkd,odlv,nkq,oqlv,0)
      call upack2('str mxnbr rmax',sstr,mxnbr,rmaxs)
      call sp2cls('spec idmod',sspec,w(oics),nl,1,nclasp,oidmod)
      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmx)
      call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)
      call sp2cls('spec rham',sspec,w(oics),1,1,nclasp,orham)
      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas

      stdo = lgunit(1)

      j = 10
      if (cmdopt('--shorten',j-1,0,outs)) then
        call shorps(nbasp,plat,modep,w(opos),w(opos))
      endif

C ... Read from restart file
      if (cmdopt('--rs=',5,0,outs)) then
        irs(1) = igets('ctrl lrs,7',sctrl)
        if (irs(1) .gt. 0) then
          ifi = fopna('rsta',-1,0)
          call asars(0,sctrl,ssite,sspec,slat,spot,sbz,
     .      w(opnu),w(oqnu),.false.,ifi)
          call fclr('rsta',ifi)
C         call shoctl(sctrl,sspec,spot,sarray,F,stdo)
C         call rx('done')
        endif
      endif

C --- Neighbor tables and sphere overlaps ---
      if (getdig(mode,0,2) .ne. 0) then

      if (rmaxs .le. 0d0) then
        rmaxs = 2.7d0*avw
        call info5(30,0,0,'%1f'//
     .    'Use default rmaxs = %;3d a.u. = %;3d*avw = %;3d*alat',
     .    rmaxs,rmaxs/avw,rmaxs/alat,0,0)
      endif

C ... Get neighbor table iax for each atom in the cluster
      if (lpbc .eq. 0) then
        i = 3
        j = -1
      elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
        i = 2
        j = 1
      else
        call rx('ASASTR: not implemented for lpbc>1')
      endif
      mxcsiz = igets('str mxnbr',sstr)
      call pshpr(iprint()-20)
      call pairs(nbas,nbasp,alat,plat,rmaxs/2,w(opos),
     .  -1,i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
      call poppr

C --- Print out a few superlattice vectors ---
      j = 6
      if (cmdopt('--slat',j-1,0,outs)) then
      if (iprint() .ge. 10) then
        call info0(10,1,0,' LMCHK:  print multiples of plat%N'//
     .    '  i1  i2  i3%7fx%11fy%11fz%11flen')
        do  i = -2, 2
        do  j = -2, 2
        do  k = -2, 2
        xx = 0
        do  m = 1, 3
          xv(m) = i*plat(m,1) + j*plat(m,2) + k*plat(m,3)
          xx = xx + xv(m)**2
        enddo
        xx = dsqrt(xx)
        print 368, i,j,k, xv(1), xv(2), xv(3), xx
  368   format(3i4, 3f12.7, 1x, f12.5)
        enddo
        enddo
        enddo
      endif
      endif

C --- Find sphere overlaps ---
      j = 9
      if (cmdopt('--getwsr',j-1,0,outs)) then
        call info(10,1,0,' ... Make sphere radii',0,0)
C        xx = dglob('lrel',1d0,1)
C        xx = dglob('nsp',1d0,1)
C       Initial estimate for sphere radii: overlapping atom potentials
        call defrr(ozz,nspec)
        call spackv(10,'spec z',sspec,1,nspec,w(ozz))
        call defrr(ormt,nspec)
        call spackv(10,'spec rmt',sspec,1,nspec,w(ormt))
        call defi(olock,-nspec)
        do  i = 1, nspec
        call ivset(w(olock),i,i,bitand(igetss('spec mxcst',i,sspec),2))
        enddo
        if (lpbc .eq. 0) then
          i = 3
        elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
          i = 2
        else
          call rx('LMAUX: not implemented for lpbc>1')
        endif
        call makrm0(101,nspec,nbas,alat,plat,w(opos),slabl,w(oips),
     .    modep,w(olock),w(ozz),w(ormt))
C   ... Scale sphere radii satisfying constraints
        call upack('ctrl omax1 omax2 wsrmax',sctrl,omax1,omax2,wsrmax,0,
     .    0)
        call sclwsr(20,nbas,nbasp,nspec,alat,plat,w(opos),w(oips),modep,
     .    slabl,w(ozz),w(olock),1d0,wsrmax,omax1,omax2,w(ormt))
        call spackv(11,'spec rmt',sspec,1,nspec,w(ormt))
        nclspp = max(2*nclasp-nclass,nspec)
        call sp2cls('spec rmt',sspec,w(oics),1,1,nclspp,ormax)
        call pack1('array ormax',sarray,ormax)
      endif

C --- Show neighbors by shell ---
      outs = ' '
      j = 8
      if (cmdopt('--shell',j-1,0,outs)) then
        call shoshl(outs(j:),nbas,w(opos),plat,mxnbr,w(oz),slabl,
     .    w(oclabl),w(oips),w(oipc),w(oves),w(oeula),nclass)
      endif

C --- Show angles between neighbors ---
      j = 9
      if (cmdopt('--angles',j-1,0,outs)) then
        call shoang(outs(j:),nbas,w(opos),plat,mxnbr,slabl,w(oips))
      endif

C --- Check whether basis w(opos) differs from site file by translation ---
      j = 9
      if (cmdopt('--basis=',j-1,0,outs)) then
C       call upack('lat nsgrp osymgr',slat,ngrp,og,0,0,0)
        fnam = outs(j:)
        call info(20,1,0,' checking whether basis equivalent to file '
     .    //fnam//'...',0,0)
        j = iosits(8070,3d0,0,fnam,ifi,slabl,alat,plat,nbas,nspec,
     .    sspec,ssite)
        call defrr(opos2,3*nbas)
        call defi (oips2, nbas)
        call spackv(10,'site pos',ssite,1,nbas,w(opos2))
        call spackv(10,'site spec',ssite,1,nbas,w(oips2))

        call defrr(og,   9*ngmx)
        call symlat(plat,ngrp,w(og),j)
        j = cmplat(nbas,plat,plat,ngrp,w(og),w(oips),w(opos),w(oips2),
     .    w(opos2))
        call fexit(j,1,' Exit %i lmchk --basis= ...',j)
      endif

C ... Write positions in Cartesian coordinates and as multiples plat
      if (iprint() .ge. 50) then
      write(stdo,357)
  357 format(/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
C     qlat = (plat+)^-1
      call dinv33(plat,1,qlat,xx)
      do  i = 1, nbas
        call dpscop(w(opos),xv,3,3*i-2,1,1d0)
C       posp+ = (plat)^-1 pos+
        call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
        ip = ival(w(oips),i)
        print 345, i, slabl(ip), (xv(j),j=1,3), (xv(3+j),j=1,3)
  345   format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
      enddo
      endif

C --- Print overlaps, optionally minimize wrt spec'd sites ---
      outs = ' '
      i = 6
      swtmp = cmdopt('-mino',5,0,outs)
      swtmp = cmdopt('--mino',6,0,outs)
      if (swtmp) i = 7
      j = 1
      if (iprint() .lt. 30) j = 0
      call ovmin(outs(i:),nbas,nbasp,alat,plat,w(ormax),w(ormax),
     .  w(oclabl),w(oipc),modep,w(oz),ontab,oiax,w(opos),j)
C ... Write positions to file
      if (cmdopt('--wpos=',7,0,outs))
     .  call iopos(T,0,outs(8:),nbasp,w(opos))

C ... Inner products between Euler angles
      if (mod(lncol,2) .eq. 1) then
        call amagnc(nbas,nl,w(oipc),w(1),1,w(oqnu),w(oeula),neul,1,xv,
     .    w(oaamom))
        print '(1x)'
        call angtab(nbas,w(opos),alat,w(ormax),qss,w(oqlv),
     .    w(odlv),nkd,w(oipc),neul,w(oeula))
      endif

      endif

C --- Plane analysis branch ---
      if (getdig(mode,1,2) .ne. 0) then
        call upack('ctrl nbas nbasp nl',sctrl,nbas,nbasp,nl,0,0)
        call upack('pot opnu oqnu oqc oqt',spot,opnu,oqnu,oqc,oqt,0)
        call defrr(oves,nbaspp)
        call supot(1,sctrl,slat,spot)
        call asamad(sctrl,spot,sarray,slat,sspec,0,
     .    w(opnu),w(oqnu),0d0,w(oves),emad,trumad,vmtz)
        call plana(wksize,npadl,npadr,nbaspp,slabl,slat,sspec,ssite,
     .    w(oves),w(opnu),w(oqnu))
        return
      endif

C --- Generate input file to xbs program ---
      if (getdig(mode,2,2) .ne. 0) then
      ifi = fopn('XBS')
      facrmx = .5d0
      ip = 4
      if (cmdopt('-bs=',ip,0,outs)) then
        call skipbl(outs,len(outs),ip)
        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,facrmx)
      endif
      facrng = 1d0
      ip = 4
      if (cmdopt('-ss=',ip,0,outs)) then
        call skipbl(outs,len(outs),ip)
        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,facrng)
      endif
C ... Copy wsr*facrmx into rmax, wsr*facrng into rham (if nonzero)
      call defdr(ormx, nclass)
      call dpcopy(w(ormax),w(ormx),1,nclass,facrmx)
C ... Copy wsr*facrng into rham, if zero
      do  30  i = 1, nclass
      xx = dval(w(orham),i)
   30 if (xx .eq. 0d0 .or. xx .eq. NULLI)
     .  call dvset(w(orham),i,i,facrng*dval(w(ormax),i))
      if (iprint() .ge. 20) then
        call awrit2('%N ball size = %d * sphere size;  '//
     .    'def stick length = %d * sum sphere sizes',
     .    ' ',80,stdo,facrmx,facrng)
      endif
      if (cmdopt('-spec',5,0,outs) .or. cmdopt('--spec',6,0,outs)) then
        nc = nspec
        oip = oips
      else
        nc = nclass
        oip = oipc
      endif
      call ioxbs(ifi,nbas,nc,alat,plat,w(orham),w(ormx),w(oclabl),
     .  w(oip),w(oz),w(opos))
      call fclose(ifi)
      endif

C --- Shift pp's (and optionally) moments by enu ---
C     pp's are remade using the potential if available.
C     use -enu=val to shift all enu's to val.
C     Use -mom if to save shifted moments.  Potential NOT remade.
      if (getdig(mode,3,2) .ne. 0) then
      ip = 5
      if (cmdopt('-enu=',ip,0,outs)) then
        call skipbl(outs,len(outs),ip)
        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,enu)
        if (i .eq. -1) call rxs('LMSHF: failed to parse ',outs)
      else
        call rx('LMSHF: missing argument -enu=val')
      endif
      call supot(1,sctrl,slat,spot)
      call defdr(oeold,3*nl*nsp*nclass)
      call defdr(opold,  nl*nsp*nclass)
      call dcopy(3*nl*nsp*nclass,w(oqnu),1,w(oeold),1)
      call dcopy(  nl*nsp*nclass,w(opnu),1,w(opold),1)
      call pp2enu(sspec,nclass,w(oics),w(onrc),nsp,nl,w(ormax),avw,
     .  w(oqnu),w(opnu),w(oqnu),xv,w(oinitc),w(oclabl),enu,w(opp))
      if (.not. cmdopt('-mom',4,0,outs)) then
        oqnu = oeold
        opnu  = opold
      endif
      call asvsph(sctrl,slat,sarray,sspec,sham,spot,0d0,2,ehterm,i)
      call rx0('done')

      call sp2cls('spec a',sspec,w(oics),1,1,nclass,oamsh)
      call sp2cls('spec nr',sspec,w(oics),1,1,nclass,onrmsh)
      call prmx('a',w(oamsh),nclass,nclass,1)
      call mkppar(w(oclabl),nl,nsp,nclass,w(olmx),w(oz),w(ormax),avw,
     .  w(oamsh),w(onrmsh),w(opnu),w(oidmod),w(oves),w(oqnu),w(opp))
      call rx('implement call to write out new parms')
C      call aiomms(nclass,nl,w(olmx),nsp,w(oclabl),w(opnu),w(oqnu),
C     .  w(oidmod),w(opp),w(orhrmx),w(ovrmax),w(oves),.true.)
      endif

C --- Interpolate core to another mesh ---
      if (getdig(mode,4,2) .ne. 0) then
        call rx('patch clabl for call to coritp')
C       call coritp(nclass,nsp,w(oclabl),w(onrmsh),w(oamsh),w(ormax))
      endif

C --- Display poles of potential functions ---
      if (getdig(mode,5,2) .ne. 0) then
      ckbas = cksumf(w(opos),3*nbas)
      call iostr(8,'STR',nl,nbasp,1,ekap,0,ckbas,-1,i,oalpha,
     .  oiax,onpr,os)
      call defdr(owk,nlspc)
      call pptrns(0,nl,w(oipc),nclass,nsp,w(oalpha),1,w(opp),w(owk))
      call rlse(owk)
      call defdr(owk2,3*2*nl**2*nbasp*nsp)
      call shopol(nl,nclass,nsp,sspec,w(oics),w(opp))
      endif

C --- Import data in other formats ---
      if (getdig(mode,6,2) .ne. 0) then

        j = nint(ssite(1)) * nbasp
        call defrr(ossite,-100*j)
        call dcopy(j,ssite,1,w(ossite),1)
        call aioxtn(sctrl,sspec,w(ossite),sarray,slat,slabl,w(opnu),
     .    w(oqnu))

C       Output to restart file
        if (cmdopt('--rs=',5,0,outs)) then
          irs(2) = igets('ctrl lrs,8+16',sctrl)/8
          if (irs(2) .gt. 0) then
            ifi = fopna('rsta',-1,0)
            call asars(0,sctrl,w(ossite),sspec,slat,spot,sbz,
     .        w(opnu),w(oqnu),.false.,-ifi)
            call fclr('rsta',ifi)
          endif
        endif
      endif

C --- Empty sphere finder ---
      if (getdig(mode,7,2) .ne. 0) then
        call upack('lat nsgrp osymgr oag',slat,nsgrp,og,oag,0,0)
        call defrr(oz,nspec)
        call spackv(10,'spec z',sspec,1,nspec,w(oz))
        call defrr(ormt,nspec)
        call spackv(10,'spec rmt',sspec,1,nspec,w(ormt))

        call defi(olock,-nspec)
        do  i = 1, nspec
        call ivset(w(olock),i,i,bitand(igetss('spec mxcst',i,sspec),2))
        enddo
        call defrr(ozc,mxclas)
        call defrr(ormtc,mxclas)
        call defrr(olockc,mxclas)
        owk = onrc
        call defi(onrc,mxclas)
        call icopy(nclass,w(owk),1,w(onrc),1)

        if (nclass .ge. mxclas) call rx('lmaux: increase mxclas')
        call spec2c(nspec,nclass,w(oics),w(ormt),w(ormtc),w(oz),w(ozc),
     .    w(olock),w(olockc))
        call findes(sctrl,alat,w(oclabl),nbas,nclass,nl,w(onrc),mxclas,
     .    nsgrp,opos,oipc,plat,w(og),w(oag),w(olockc),w(ormtc),w(ozc))

      endif

      end

       subroutine spec2c(nspec,nclass,ics,rmts,rmtc,z,zc,lock,lockc)
C- Copy species data to class data
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspec
Ci   nclass:number of inequivalent classes
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   rmtc  :rmt by species
Ci   z     :z by species
Co Outputs
Co   rmtc  :rmt by class
Co   zc    :Z by class
co   lockc :lock by class
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   11 Aug 08 
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nspec,nclass,ics(nclass),lock(nspec),lockc(nclass)
      double precision rmts(nspec),rmtc(nclass),z(nspec),zc(nclass)
C ... Local parameters
      integer j,k

      do  k = 1, nclass    
        j = ics(k)
        rmtc(k) = rmts(j)
        zc(k) = z(j)
        lockc(k) = lock(j)
C       if (iprint() .gt. 60) write(*,310) k,rmtc(k)
      enddo

C 310 format(1x,'class ',I3,T15,'rmt = ',f10.7)

      end
        

