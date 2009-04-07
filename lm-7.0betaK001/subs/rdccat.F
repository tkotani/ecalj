#define LM 1
      subroutine rdccat(mxrecs,recln0,unit,offset,iosw,ido,pass2,vstrn,
     .  vn,prgnam,TF,toksw,a,nspec,nsite,nbas,nl,nsp,dalat,slabl,
     .  nclasp,pnu,qnu,pp,ics,dclabl,initc,ves,
     .  sbz,sctrl,sham,spot,slat,smix,soptic,sordn,sspec,ssite,
     .  sstr,smove,sgw,stb,sstrn)
C- Reads one category from ctrl file
C ----------------------------------------------------------------
Ci Inputs
Ci   mxrecs:maximum number of records (lines) input file may contain
Ci   recln0:record length of input file
Ci   unit  :file logical unit
Ci   offset:Input file stored in recrd(offset:)
Ci   iosw  : 0 do not attempt to read input file, but show what
Ci         :   tokens rdccat would attempt to read
Ci         : 1 attempt to read tokens from input file
Ci         : 2 print out value of tokens corresponding to
Ci         :   those that would be read by iosw=1
Ci   ido:  :-1 read in file only
Ci         : 0 read for all blocks
Ci         : 1 and up: read special block
Ci   pass2 :second pass (START category)
Ci   vstrn :ascii string naming calling program; used to version check
Ci   vn    :version number
Ci   prgnam:name of main program
Ci   TF    :work array holding true-false table for
Ci   toksw :a string of integers, eg '13022' that determine, for each
Ci          category and token, whether it is optional, required, to be
Ci          ignored, or whether an internal default is to be taken.  The
Ci          first digit corresponds to the category; the second to the
Ci          first token, the third to the second, etc.  Each digit must
Ci          be between 0 and 3:
Ci            0=>optional  1=>req'd  2=>ignored  3=>use internal default
Ci          The caller must order the digits just as they are ordered in
Ci          each category in rdccat.  The orderings and internal defaults
Ci          are located at the start of each category.
Ci    a    :a string of length recln0
Cio  nspec :number of species (passed except when reading 'STRUC' categ)
Cio  nsite :number of atoms (passed except when reading 'STRUC' categ)
Cio  nbas  :size of basis (passed except when reading 'STRUC' categ)
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   dalat :added AFTER input read in, so R/W scales with alat
Ci   slabl :vector of species labels
Ci   nclasp:number of classes
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         : pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   pp    :potential parameters (atomsr.f)
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   dclabl:class name, packed as a real number
Ci   initc :one's  digit >0 if moments are read
Ci          two's  digit >0 if pp's read
Ci          four's digit >0 if spin-orbit parms read
Ci   ves   :Electrostatic potential (read in START category)
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: n w
Ci     Stored:    n w nkabc lshft
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nl nspin lasa,128+256 ltb,1+48 mdprm lordn loptc
Ci     Stored:    quit nl lgen3 lasa nmap
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: qss:4 nmto
Ci     Stored:    *
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: *
Ci     Stored:    *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ldist dist alat plat avw nkdmx plat2
Ci     Stored:    plat alat ldist dist vol avw plat2
Ci   smix  :struct for charge mixing parameters; see routine umix
Ci     Elts read: b w wc mode mmix bv lxpot
Ci     Stored:    fn
Ci   soptic:struct for parameters for optical ME; see routine uoptic
Ci     Elts read: *
Ci     Stored:    *
Ci   sordn :struct for parameters for Order-N; see routine uordn
Ci     Elts read: *
Ci     Stored:    rmaxg
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z rmt a nr lmxb lmxa idxdn mxcst p rg rsma rfoca
Ci                rsmfa lmxl lfoca kmxt idmod
Ci     Stored:    z rmt a nr lmxa lmxb lmxf alpha mxcst pb1 pb2 p q
Ci                lmxl lfoca kmxt kmxv rsmv rg rsma rfoca rsmfa name
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: pos eula
Ci     Stored:    spec pos eula ndelta
Ci     Passed to: iosits ioposs
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci     Elts read: nkaps n mix
Ci     Stored:    *
Ci   smove :struct for dynamics information; see routine umove
Ci     Elts read: nmodt prmint:6
Ci     Stored:    *
Ci   sgw
Ci     Elts read: *
Ci     Stored:    *
Ci   sstrn :struct for global strings
Ci     Elts read: mix
Ci     Stored:    *
Cr Remarks
Cr   Embedded cluster scheme:
Cr   nspec corresponds to number of classes in cluster+medium
Cr   nbas corresponds to number of atoms in cluster+medium
Cr   plat is replaced by nlate . plate (see u_ord0)
Cl Local variables
Cl   optio  0 do not attempt to read input, but show what would be read
Cl          1 attempt to read input
Cl          2 shows contents corresponding to input
Cu Updates
Cu   01 Aug 08 New ctrl->pfloat
Cu   05 Jul 08 new ham->(pwemax,pwemin,pwmode,npwpad,oveps)
Cu   20 Jun 08 New spec->ehvl, str->rmaxg, str->ivl
Cu   04 Jun 08 (ATP) new MOL= switch; rework MD: category
Cu   05 Nov 07 some cleanup in preparation for new v7 input
Cu   28 Nov 07 GAMMA representation as as option GAMMA=2
Cu   06 Aug 06 Input for strux redsigned: 2-kappa input;
Cu             New RMAXS=; NMTO SIGMA->HCR, KAPS->EL, NKAPS->NEL
Cu   27 Jun 06 Added input for sham->pmax
Cu   21 Mar 06 Input PDEF for ctrl->defm
Cu   01 Feb 06 Input for spec->rcfa
Cu   21 Dec 05 cleaner --show : requires updated slatsm/rdfiln
Cu   22 Aug 05 Bug fix when reading both from SITE and pos files
Cu    1 Aug 05 Parameters RMINES, RMAXES, NESABC for ES finder
Cu   09 Jun 05 (ATP) New Parameters for tight-binding
Cu   27 Apr 05 (WRL) parameters for LDA+U
Cu   10 Jan 05 (ATP) Parameters for molecular dynamics
Cu   15 Dec 04 Optional MT radius R/A= in place of R= or R/W=
Cu   30 Sep 04 Optionally read ASA enu from ctrl file
Cu   16 Aug 04 Changes for extended local orbitals
Cu   13 May 04 Changes to accommodate input for molecules prog.
Cu             NB: OPTIONS GRCOR= changed to HAM GGA=
Cu   29 Mar 04 Read lat->FTMESH if lat->GMAX=0 (HAM)
Cu    2 Mar 04 Read lat->rpad (EWALD)
Cu    4 Jan 04 New token ham->rsstol
Cu   23 Sep 03 SX switches redesigned for compatiblity with GW
Cu    9 Jul 03 SO switch turned on automatically if lrel=2
Cu   21 May 03 Added input for sham->sigp
Cu    8 Apr 03 read BZ GETQP before NKABC
Cu   18 Mar 03 Read lrel as an integer, for fully rel. case
Cu   08 Feb 03 Switch for B-field
Cu   12 Jan 03 Split gw->lsig into gw->rdsig and gw->mksig
Cu   19 Jul 02 Added tokens to read r.s. sigma in HAM category
Cu    1 Jun 02 change XIPMX token to integer
Cu   20 May 02 Added HAM token PMIN
Cu   09 May 02 Added SPEC token FRZWF
Cu   23 Apr 02 Added SPEC tokens SCLWSR,OMAX1,OMAX2,WSRMAX
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   21 Aug 01 Extended to read input for local orbitals.
Cu   14 Aug 01 Extended iosits to read Kotani style input
Cu   26 Apr 01 Bug fix in bz->N.W input
Cu   17 Apr 01 Extended QASA token
Cu    3 Feb 01 bz->SAVDOS, DOSWT merged into SAVDOS; BZMP removed
Cu             New tokens ham->QASA  optic->ESCISS
Cu             New input for partial DOS
Cu   24 Jan 01 correct for lmax misalignment in second spin channel
Cu             when reading pnu,qnu
Cu   22 Dec 00 (ATP) MSTAT input redesigned
Cu             First attempt at generating FTMESH from gmax
Cu             New GW category
Cu             Some revised NMTO input
Cu             Rearranged SPEC category
Cm MPI
Cm   The master process reads the ctrl file using rdfiln and
Cm   broadcasts the file as a character string to all other
Cm   processes. Then each process reads the categories.
C ----------------------------------------------------------------
C     implicit none
C ... Control parameters
#if MPI | MPIK
      include "mpif.h"
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*120 strn
      logical mlog
      integer irec,fopng
#endif
      logical pass2
      integer recln,TF(0:42)
      integer unit,offset,ido,iosw,mxrecs,recln0
      character*(*) prgnam*8, vstrn*6, toksw, a*(*)
      double precision vn,vsn

C ... Structures and quantities extracted from them
      integer nspec,nsite,nbas,nsp,nl,nclasp
      character*(*) sstrn
      character slabl(*)*8
      double precision dalat
      double precision sbz(1),sctrl(1),slat(1),sspec(1),ssite(1),
     .  sham(1),spot(1),sstr(1),soptic(1),sordn(1),smove(30),smix(1),
     .  sgw(1),stb(1)
C ... Class-related arrays
      integer initc(1),ics(1)
      double precision pnu(nl,nsp,1),qnu(3,nl,nsp,1),pp(6,nl,nsp,1),
     .  ves(1),dclabl(1)

C ... Local variables
      integer n0
      logical cmdopt,T,F,sw,sw2,lgors,a2bin,parmxp
      character*8 alabl
      character*256 outs
      parameter (T=.true., F=.false., n0=10)
      integer a2vec,casti(5),dvec1(3),dvec2(3),i,i0,i1,i1mach,idmod(n0),
     .  ifi,igets,igetss,ii,iosits,iosw2,iprint,ipsw,isw,it(10),nglob,
     .  ivec(1000),j,jc,jm1,js,kmxt,lfoca,lgunit,lmxa,lmxax,lmxb,lmxl,
     .  lp1,n,nelti(5),nmto,nr,nrmx,offi(5),offs,partks,partkf,partok,
     .  prtkss,nkaph,nlmax,lpzi,bitand,nbasp,ncor
#if FTB
      parameter ( nrmx=1001 )
#else
      parameter ( nrmx=1501 )
#endif
      double precision plat(3,3)
      double precision dum,xx,dgets,dasum,z,rmt,alat,fpi,
     .  vol,avw,xv(2*n0),tmp(2*n0),tbalph(5,5),ddot,dglob,
     .  wk(100),nmass,pmass,emass
C     logical bittst
C     integer bit
C     bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) .eq. bit)

C ... Input is stored in record(offset:)
      character*1 recrd(0:1)
      common /w/ recrd

C --- Common blocks ---
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr

      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      data tbalph/.214286D0,.000001D0,.000001D0,.000001D0,.001D-3,
     .            .287234D0,.025822D0,.000001D0,.000001D0,.001D-3,
     .            .348485D0,.053030D0,.010714D0,.000001D0,.001D-3,
     .            .385057D0,.073209D0,.022481D0,.006069D0,.001D-3,
     .            .404761D0,.087927D0,.032886D0,.012257D0,4.2735D-3/
      data dvec1 /1,1,1/ dvec2 /0,0,0/

#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#endif

      recoff = offset

C --- Do according to switch iosw ---
      do  1000  iosw2 = iosw, iosw + iosw/3
      optio = iosw2 - 2*(iosw/3)

C     print *, 'ido=',ido,optio

      goto (1,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,
     .  170,180,190,200,210,220,230,240,250,260,270), ido+1
      i = recln(recln0)
      return

    1 continue
#if MPI | MPIK
      if (procid .eq. master) then
#endif
      rewind unit
      sw = F
      alabl = '#{}% ct '
      if (cmdopt('--show',6,0,a)) alabl = '#{}% ctp'
      if (optio .eq. 1) then
        rewind unit
        call rdfile(unit,alabl,recrd(recoff),mxrecs,a,reclen,nrecs)
      endif
      if (cmdopt('--showp',7,0,a)) call cexit(0,1)
      if (iprint() .ge. 60 .and. optio .eq. 1)
     .  call awrit1(' RDCCAT: read %i lines',' ',80,lgunit(1),nrecs)
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//
     .    ' Master process broadcasting %i '//
     .    'records, length %i',' ',256,lgunit(3),nrecs,reclen)
      endif
      endif
      call MPI_BCAST(nrecs,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(reclen,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(recrd(recoff),mxrecs*recln0+100,MPI_CHARACTER,
     .               master,MPI_COMM_WORLD,ierr)
      if (procid .ne. master) then
        call gettime(datim)
        if (mlog) then
          call awrit3(' rdccat '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' received %i records',
     .      ' ',256,lgunit(3),procid,numprocs,nrecs)
        endif
        if (iprint() .gt. 50 .and. mlog) then
          call awrit0(' rdccat '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' received this ctrl file:',' ',256,lgunit(3))
          do  irec = 0, nrecs-1
            write (lgunit(3),123) (recrd(recoff+i),
     .        i=reclen*irec,reclen*(irec+1)-1)
  123       format (120a1)
            call flush(lgunit(3))
          enddo
        endif
      endif
      return
#else
      return
#endif

C --- Header ---
C Token ordering: HEADER
   10 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading HEADER',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0',toksw,3,TF)
      call getcat(recrd,'HEADER ',' ',TF)
      if (.not. noerr) goto 1000
      if (optio .eq. 1) then
        ifi = lgunit(2)
        if (ifi .gt. 0 .and. iprint() .gt. 0) call headl2(prgnam,0,ifi)
        if (noerr) then
          call lstra('strn jobid',i,offs,j)
          call ustrn(tmp,-offs,1,offi,i1,reclen)
          i = 6
          call skipbl(recrd(catbeg),reclen,i)
          call skpblb(recrd(catbeg),reclen,j)
          call strncp(sstrn(offi(1):),recrd,1,catbeg+i+1,j)
          call ustrn(tmp,-offs,1,offi,i1,j-i+1)
C          outs = sstrn(offi(1):offi(1)+j-i)
C          outs = ' HEADER '// sstrn(offi(1):offi(1)+j-i)
#if MPI | MPIK
          if (procid .eq. master) then
#endif
          outs = ' HEADER '// sstrn(offi(1):offi(1)+j-i)
          call strip(outs,i,j)
          write(lgunit(1),'(1x,a)') outs(i:j)
C         call awrit0(' HEADER '//outs,' ',-len(outs),lgunit(1))
          if (ifi .gt. 0) write(ifi,'(1x,a)') outs(i:j)
#if MPI | MPIK
          endif
#endif
        endif
      endif
      goto 1000

C --- Version ---
C Token ordering: VERS; vstrn
   20 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading VERS',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('11',toksw,3,TF)
      call getcat(recrd,'VERS ',' ',TF)
      if (.not. noerr) goto 1000
      vsn = dble(int(vn*100))/100
      if (optio .eq. 0) then
        call awrit2(' ---- Input for program version '//vstrn//
     .    '%a %d.%i ---',' ',80,lgunit(1),vsn,nint((vn-vsn)*1d5))
      else
        i = partok(recrd(catbeg),vstrn,'-',xx,' ',-1,4,0,TF(1))
        if (noerr .and. xx .lt. int(vsn)) call fexit(-1,1,
     .    '  Exit -1 RDCCAT: control file earlier than v%d',vsn)
      endif
      goto 1000

C --- IO ---
C Token ordering: IO; SHOW HELP VERBOS WKP IACTIV TIM
   30 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading IO',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0000000',toksw,3,TF)
      call getcat(recrd,'IO ',' ',TF)
      if (.not. noerr) goto 1000
      sw = F
      i0 = partok(recrd(catbeg),'SHOW=','=',sw,' ',-1,0,0,TF(1))
      if (cmdopt('--show',6,0,a)) sw = T
      if (sw) iosw = 3
C      if (noerr .and. sw) iosw = 3
      i0 = partok(recrd(catbeg),'HELP=','=',sw,' ',-1,0,0,TF(2))
      if (cmdopt('--input',7,0,a)) sw = T
      if (noerr .and. sw) iosw = 0
      i0 = partok(recrd(catbeg),'VERBOS=','= ,',ivec,' ',5,2,0,TF(3))
C ... Override ivec w/ -pr commmand-line arg
      if ((cmdopt('--pr',4,0,a) .or. cmdopt('-pr',3,0,a)) .and.
     .    optio .eq. 1) then
        i = 4
        if (cmdopt('-pr',3,0,a)) i = 3
        i = a2vec(a,len(a),i,2,', ',2,2,5,it,ivec)
        if (i .lt. 0) call rxs('error parsing switch ',a)
        i0 = max(i0,i)
      endif
C     Copy verbosities to print stack
      if (i0 .ge. 1) then
        do  36  i = 0, 5
   36   call sprt(i,ivec(min(i+1,i0)))
      endif
      i0 = partok(recrd(catbeg),'WKP=','=',sw,' ',-1,0,0,TF(4))
      if (noerr .and. sw) call wkprnt(1)
      i0 = partok(recrd(catbeg),'IACTIV=','=',sw,' ',-1,0,0,TF(5))
      if (optio .eq. 1 .and. cmdopt('--iactiv=no',11,0,a)) then
        call initqu(F)
      elseif (optio .eq. 1 .and. cmdopt('--iactiv',7,0,a)) then
        call initqu(T)
      elseif (optio .eq. 1 .and. cmdopt('--no-iactiv ',7,0,a)) then
        call initqu(F)
      elseif (noerr) then
        call initqu(sw)
      endif

C      i0 = partok(recrd(catbeg),'OUTPUT=','=',256,outs,1,1,0,TF(5))
C      if (i0 .eq. 1 .and. outs(1:1) .ne. '*') then
C        i = fopng(outs,-1,0)
C        call i1macs(2,i)
C      endif

      ivec(2) = 0
      i0 = partok(recrd(catbeg),'TIM=','=',ivec,' ',2,2,0,TF(3))
C     Override ivec with '--time=' commmand-line arg
      if (optio .eq. 1 .and. cmdopt('--time=',7,0,a)) then
        i = 7
        i = a2vec(a,len(a),i,2,', ',2,2,2,it,ivec)
        if (i .lt. 0) call rxs('error parsing switch ',a)
        i0 = max(i0,i)
      endif
      if (i0 .ge. 1) call tcinit(ivec(2),ivec(1))
      return

C --- CONST ---
C Token ordering: CONST
   40 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading CONST',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0',toksw,3,TF)
      call getcat(recrd,'CONST ',
     .  'constants may declared for use in expressions',TF)
      if (.not. noerr) goto 1000
      i = 5
      if (optio .eq. 1) call parsyv(recrd(catbeg),catsiz-1,1999,0,i)
      if (optio .eq. 2) call shosyv(0,0,0,i1mach(2))
      goto 1000

C --- Struc ---
C Token ordering
C (0-8) STRUC NSPEC NBAS ALAT PLAT NL SHEAR NLFIT (ROT,DEFGRD,STRAIN)
C (9-13) LOADV DALAT SLAT FILE NBASP
   50 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading STRUC',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('11111002000202',toksw,3,TF)
      call getcat(recrd,'STRUC ','Lattice information',TF)
      if (.not. noerr) goto 1000
      nchoos = 2
      i0 = partok(recrd(catbeg),'NSPEC=','=',nspec,
     .  'Number of atom species to read',-1,2,0,TF(1))
      if (.not. noerr) then
        i0 = partok(recrd(catbeg),'NCLASS=','=',nspec,' ',-1,2,0,TF(1))
      endif
      nchoos = 0
      if (TF(12) .ne. 2) then
      if (optio .eq. 0) print *,'  ...    Read NBAS, ',
     .  'PLAT, and optionally ALAT from site file with this token:'
      i0 = 0
C     Initialize alat=-99 to flag no alat has been input
      alat = -99
      if (optio .lt. 2)
     .  i0 = partok(recrd(catbeg),'FILE=','=',8,outs,1,1,0,TF(12))
      if (i0 .eq. 1) then
        j = iosits(7000,3d0,0,outs(1:8),js,slabl,alat,plat,nbas,nspec,
     .    sspec,ssite)
        nsite = nbas
        call pack1('lat plat',slat,plat)
        if (alat .ne. -99) call pack1('lat alat',slat,alat)
        call fclose(js)
      endif
      if (optio .eq. 0) print *, '  ...    ',
     .  'Otherwise, read them from the ctrl file:'
      endif
      if (i0 .lt. 1) then
        i0 = partok(recrd(catbeg),'NBAS=','=',nbas,'Size of basis',
     .    -1,2,0,TF(2))
        nsite = nbas
        i0 = partks(recrd(catbeg),'PLAT=','=','lat plat',slat,0,
     .    'primitive lattice vectors, in units of alat',-9,4,0,TF(4))
      endif
      if (optio .eq. 1) call pack1('ctrl nbasp',sctrl,nbas)
      i0 = partks(recrd(catbeg),'NBASP=','=','ctrl nbasp',sctrl,0,
     .  'nbas + no. of point multipoles',-1,2,0,TF(13))
      call upack1('ctrl nbasp',sctrl,nbasp)
      nsite = nbasp

      if (alat .eq. -99)
     .  i0 = partks(recrd(catbeg),'ALAT=','=','lat alat',
     .  slat,0,'scaling of lattice vectors, in a.u.',-1,4,0,TF(3))
      i0 = partks(recrd(catbeg),'SLAT=','=','lat slat',slat,0,
     .  'supercell lattice vectors',-9,4,0,TF(11))
      i0 = partok(recrd(catbeg),'DALAT=','=',dalat,
     .  'added to alat after input is read',-1,4,0,TF(10))
      i0 = partks(recrd(catbeg),'NL=','=','ctrl nl',sctrl,0,
     .  'lmax+1 for basis and augmentation',-1,2,0,TF(5))
      nl = igets('ctrl nl',sctrl)
C ... lattice distortion or rotation
      nchoos = 4
      i0 = partks(recrd(catbeg),'SHEAR=','=','lat gam',slat,0,
     .  'voluming-conserving shear of PLAT (1=ideal)',-4,4,0,TF(6))
      call upack2('lat ldist dist',slat,j,xv)
      if (optio .ne. 2) then
        i0 = partok(recrd(catbeg),'ROT=','=',len(outs),outs,1,1,0,TF(8))
        if (i0 .eq. 1) then
          call a2rotm(outs,F,iprint()-10,xv)
          j=2
        endif
      endif
      if (optio .ne. 2 .and. j .eq. 0 .or. j .eq. 2) then
        i0 = partok(recrd(catbeg),'DEFGRD=','=',xv,' ',-9,4,0,TF(8))
        if (noerr) j=2
      endif
      if (optio .ne. 2 .and. j .eq. 0 .or. j .eq. 3) then
        i0 = partok(recrd(catbeg),'STRAIN=','=',xv,
     .       '6 Voigt strain matrix elements ',-6,4,0,TF(8))
        if (i0 .eq. 6 .or. optio .eq. 0 .or. optio .eq. 2) then
          i0 = partok(recrd(catbeg),'ALPHA=','=',xv(7),
     .      'amplitude of (Voigt) strain ',-1,4,0,TF(8))
        endif
        if (noerr) j=3
      endif
      nchoos = 0
      call pack2('lat ldist dist',slat,j,xv)

      i0 = partok(recrd(catbeg),'LOADV ',' ',sw,'When token is present,'
     .  //' vol and avw are loaded into variables table',0,0,0,TF(9))
      if (optio .eq. 1) then
        call upack('lat alat plat',slat,alat,plat,0,0,0)
        fpi = 16*datan(1d0)
        vol = alat**3*dabs(
     .        plat(1,1)*(plat(2,2)*plat(3,3) - plat(2,3)*plat(3,2)) +
     .        plat(1,2)*(plat(2,3)*plat(3,1) - plat(2,1)*plat(3,3)) +
     .        plat(1,3)*(plat(2,1)*plat(3,2) - plat(2,2)*plat(3,1)))
        avw = (3/(fpi*nbas)*vol)**(1d0/3d0)
        call pack2('lat vol avw',slat,vol,avw)
        if (sw) then
          call lodsyv('vol',1,vol,i)
          call lodsyv('avw',1,avw,i)
        endif
      endif

      goto 1000

C --- Options ---
C Token ordering: OPTIONS
C (1-7)   REL CCOR TPAN ADNF NSPIN FRZ INVIT
C (8-13)  NONCOL SDMOD SDPRM SS SO STONER
C (14-21) TWOC ELIN NEWREP NOHYB MTCOR ZBAK NRMIX Q=
C (22-29) SHARM EMCLUS SX SCR HF NSPH XCQS SAVVEC
C (30-32) (RMINES,RMAXES) NESABC PFLOAT
   60 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading OPTIONS',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('000200000000000000000000000222002',toksw,3,TF)
      call getcat(recrd,'OPTIONS ',' ',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'HF=','=','ctrl lcd,2',sctrl,0,
     .  'T for non-selfconsistent Harris',1,0,0,TF(26))
C     ASA-specific
      i0 = partks(recrd(catbeg),'CCOR=','=','ctrl lasa,4',sctrl,0,
     .  'Turn on combined correction',1,0,0,TF(2))
      i0 = partks(recrd(catbeg),'ADNF=','=','ctrl lham,4',sctrl,0,
     .  'Turn on automatic downfolding',1,0,0,TF(4))
C      i0 = partks(recrd(catbeg),'GRCOR=','=','ctrl lxcf,32',sctrl,0,
C     .  'Turn on gradient corrections',1,0,0,TF(3))
      i0 = partks(recrd(catbeg),'FRZ=','=','ctrl lcd,1',sctrl,0,
     .  'Freeze core',1,0,0,TF(6))
      i0 = partks(recrd(catbeg),'NSPH=','=','ctrl lasa,32',sctrl,0,
     .  'generate multipole moments in output density',-1,0,0,TF(27))
      i0 = partks(recrd(catbeg),'XCQS=','=','ctrl lcd,64',sctrl,0,
     .  'XC potential in reciprocal space',-1,0,0,TF(28))
      i0 = partks(recrd(catbeg),'SHARM=','=','ctrl lham,256',sctrl,0,
     .  'Use true spherical harmonics',1,0,0,TF(22))
      i0 = partks(recrd(catbeg),'NSPIN=','=','ctrl nspin',sctrl,0,
     .  'set to 2 for spin polarized calculations',1,2,0,TF(5))
      i0 = partks(recrd(catbeg),'PFLOAT=','=','ctrl pfloat',sctrl,0,
     .  'controls how band CG is determined in floating Pnu.'//
     .  '%N%10f0 Band CG found by traditional method'//
     .  '%N%10f1 Band CG found from true energy moment of density',
     .  1,2,0,TF(32))
      nsp = igets('ctrl nspin',sctrl)
      call isanrg(nsp,1,2,'rdccat','nsp',.true.)
#if LM
      if (optio .eq. 0) then
        print *,'  ... the following tokens for nspin=2:'
      elseif (nsp .eq. 1) then
        call ivset(TF(1),8,12,2)
      endif
      i0 = partks(recrd(catbeg),'NONCOL=','=','ctrl lncol,1',sctrl,0,
     .  'noncollinear magnetism',1,0,0,TF(8))
      i0 = partks(recrd(catbeg),'BFIELD=','=','ctrl lncol,8',sctrl,0,
     .  'applied magnetic field (requires NONCOL=t',1,0,0,TF(8))

      i0 = partks(recrd(catbeg),'SS=','=','ham qss',sham,0,'dir'//
     .'ection vector and angle of magnetic spin spiral',4,4,0,TF(11))
      if (i0 .ge. 3 .and. ((dgets('ham qss:4',sham) .ne. 0) .or. T))
     .  call lsets('ctrl lncol',sctrl,.true.,2)

      i0 = partks(recrd(catbeg),'REL=','=','ctrl lrel',sctrl,0,
     .  '0 for nonrelativistic Schrodinger equation'//
     .  '%N%10f1 for scalar relativistic Schrodinger equation'//
     .  '%N%10f2 for Dirac equation',
     .  1,2,0,TF(1))
C     Unpack and pack SO coupling
      j =   isw(lgors('ctrl lncol,4',sctrl))
     .  + 2*isw(lgors('ctrl lncol,32',sctrl))
     .  + 3*isw(lgors('ctrl lncol,64',sctrl))
      i0 = partok(recrd(catbeg),'SO=','=',j,
     .  'spin-orbit coupling (for REL=1)'//
     .  '%N%10f0 : no SO coupling (default)'//
     .  '%N%10f1 : Add L.S to hamiltonian'//
     .  '%N%10f2 : Add Lz.Sz only to hamiltonian'//
     .  '%N%10f3 : Like 2, but also compute <L.S-LzSz> by perturbation',
     .  -1,2,0,TF(12))
      call isanrg(j,0,3,' rdccat:','SO',.true.)
      if (optio .eq. 1) then
      if (j .eq. 1) then
        call lsets('ctrl lncol',sctrl,.true.,4)
      elseif (j .eq. 2) then
        call lsets('ctrl lncol',sctrl,.true.,32)
      elseif (j .eq. 3) then
        call lsets('ctrl lncol',sctrl,.true.,64)
      endif
      endif
C     Fully relativistic => spin-orbit coupling
      if (optio .eq. 1 .and. igets('ctrl lrel',sctrl) .eq. 2) then
        call lsets('ctrl lncol',sctrl,.true.,4)
        call lsets('ham lncol',sham,.true.,4)
      endif
C     If spin-orbit or SS, also turn on noncollinear
      call lsets('ctrl lncol',sctrl,isw(lgors('ctrl lncol,31',sctrl)),1)
      if (optio .eq. 0 .or. lgors('ctrl lncol,1',sctrl))
     .  i0 = partks(recrd(catbeg),'SDYN:',':, ',
     .  'ctrl lncol;force,16 sdmod?$ sdprm',
     .  sctrl,0,'Parameters for spin dynamics:'//
     .  '%N%10fforce=T calculate magnetic forces; move atoms according'
     .  //' to mode:'//
     .  '%N%10fsdmod= 0: output Euler angles from dens. matrix'//
     .  '%N%17f1: relax along force'//
     .  '%N%17f2: dynamics'//
     .  '%N%17fAdd 10   to mix angles independently of P,Q'//
     .  '%N%17fAdd 1000 to suppress self-const in P,Q'//
     .  '%N%10fsdprm= fscal,tau,Nose-etot0,max-theta,etol',
     .  -1,5,0,TF(9))

      i0 = partks(recrd(catbeg),'SCR=','=','ctrl lscr',sctrl,0,
     . 'Use scr to accelerate convergence:',-1,2,0,TF(25))
      if (optio .eq. 0 .and. TF(25) .ne. 2) then
        print 356
  356   format(
     .     10x,'0 do nothing'
     .    /10x,'1 Make ASA static response function',
     .          ' (see documentation)'
     .    /10x,'2 Use response to screen output q and ves'
     .    /10x,'4 Use model response to screen output q'
     .    /10x,'6 Use response to screen output ves only'
     .    /8x,'Add 1 to combine mode 1 with another mode',
     .    /8x,'Add 10*k to compute intra-site contribution to vbare',
     .          ' each kth iteration'
     .    /8x,'Add 100*k to compute response function',
     .          ' on every kth iteration')

      endif
      i0 = partks(recrd(catbeg),'SX=','=','ctrl lsx',sctrl,0,
     .  'Screened exchange:',-1,2,0,TF(24))
      if (optio .eq. 0 .and. TF(24) .ne. 2) then
        print 357
  357   format(
     .     10x,' 0 do nothing'
     .    /10x,' 1 Calculate SX Sigma'
     .    /10x,'11 Calculate SX Sigma w/ local W')
      endif
      i0 = prtkss(recrd(catbeg),'SXOPTS=','=','strn sxopt',sstrn,reclen,
     .  'Options for screened exchange',0,TF(24))
      i0 = partks(recrd(catbeg),'TWOC=','=','ctrl lham,3',sctrl,0,
     .  'two-center ASA hamiltonian',1,2,0,TF(14))

C      i0 = partks(recrd(catbeg),'GRCOR=','=','ctrl lxcf,32',sctrl,0,
C     .  'Turn on gradient corrections',1,0,0,TF(3))

      i0 = partks(recrd(catbeg),'TPAN=','=','ctrl lham,1',sctrl,0,
     .  'two-panel (molecule)',1,0,0,TF(3))

C     Unpack and pack GAMMA
      j = isw(lgors('ctrl lham,128',sctrl))
     .  + isw(lgors('ctrl lasa,512',sctrl))
      i0 = partok(recrd(catbeg),'GAMMA=','=',j,
     .  'gamma representation'//
     .  '%N%10f0 : alpha (shortest range, default)'//
     .  '%N%10f1 : gamma'//
     . '%N%09f>1 : (gamma(1)+gamma(nsp))/2'//
     .  '%N%10f  : Some noncollinear algorithms automatically'//
     .  '%N%10f  : choose this latter representation',
     .  -1,2,0,TF(14))
      call isanrg(j,0,2,' rdccat:','GAMMA',.true.)
      if (optio .eq. 1) then
      if (j .gt. 0) then
        call lsets('ctrl lham',sctrl,.true.,128)
      endif
      if (j .gt. 1) then
        call lsets('ctrl lasa',sctrl,.true.,512)
      endif
      endif
C      i0 = partks(recrd(catbeg),'GAMMA=','=','ctrl lham,128',sctrl,0,
C     .  'gamma rep''n',1,0,0,TF(14))
      i0 = partks(recrd(catbeg),'ELIN=','=','ctrl elin',sctrl,0,
     .  'energy to linearize for ccor in 2C hamiltonian ',1,4,0,TF(15))
      i0 = partks(recrd(catbeg),'NEWREP=','=','ctrl lham,8',sctrl,0,
     .  'Interactively transform tb-representation',1,0,0,TF(16))
      i0 = partks(recrd(catbeg),'NOHYB=','=','ctrl lham,16',sctrl,0,
     .  'turns off hybridisation ',1,0,0,TF(17))
      i0 = partks(recrd(catbeg),'SAVVEC=','=','ctrl lham,64',sctrl,0,
     .  'Save eigenvectors on disk',-1,0,0,TF(29))
      i0 = partks(recrd(catbeg),'ZBAK:',':, ','ctrl lasa;mtcor,64 '//
     .  'zbak?$',sctrl,0,'mtcor turns on Ewald MT correction%N%10f'//
     .  'zbak, first entry adds constant background charge to density'//
     .  '%N%10fOptional second number specifies background q for'//
     .  ' Ewald MT correction',
     .  1,5,0,TF(19))
      i0 = partks(recrd(catbeg),'NRMIX=','=','ctrl smalit',sctrl,0,
     .  'Sphere program, (1) max iter; (2) no. prior iter for Anderson'
     .  //' mixing ',2,2,0,TF(20))

      i0 = partks(recrd(catbeg),'RMINES=','=','ctrl rmines',sctrl,0,
     .  'Minimum MT radius when finding new ES',1,4,0,TF(30))
      i0 = partks(recrd(catbeg),'RMAXES=','=','ctrl rmaxes',sctrl,0,
     .  'Minimum MT radius when finding new ES',1,4,0,TF(30))
      call pvcat1('NESABC=','No. divisions '//
     .  'when searching for empty spheres','ctrl nesabc',sctrl,TF(31))
      i0 = partks(recrd(catbeg),'STONER=','=','ctrl lstonr',sctrl,0,
     .  'generalised Stoner rigid band calculation '//
     .  '%N%9f second argument is number of points; '//
     .  'third for graphical output',
     .  3,2,0,TF(13))
      nchoos = 3
      i0 = partks(recrd(catbeg),'Q=SHOW ',' ','ctrl quit,1',sctrl,0,
     .  'Stop after band pass',0,0,0,TF(21))
      i0 = partks(recrd(catbeg),'Q=ATOM ',' ','ctrl quit,2',sctrl,0,
     .  'Stop after band pass',0,0,0,TF(21))
      i0 = partks(recrd(catbeg),'Q=BAND ',' ','ctrl quit,4',sctrl,0,
     .  'Stop after band pass',0,0,0,TF(21))
      i0 = partks(recrd(catbeg),'Q=HAM ',' ','ctrl quit,4',sctrl,0,
     .  'Stop after band pass',0,0,0,TF(21))
      if (cmdopt('--quit=show',11,0,a)) call pack1('ctrl quit',sctrl,1)
      if (cmdopt('--quit=atom',11,0,a)) call pack1('ctrl quit',sctrl,2)
      if (cmdopt('--quit=band',11,0,a)) call pack1('ctrl quit',sctrl,4)
      if (cmdopt('--quit=ham',10,0,a)) call pack1('ctrl quit',sctrl,8)
C      nchoos = 2
C      i0 = partks(recrd(catbeg),'VDB=','=','ctrl lordn,4',sctrl,0,
C     .  'Vanderbilt order N scheme',-1,0,0,TF(22))
C      i0 = partks(recrd(catbeg),'EMCLUS=','=','ctrl lordn,8',sctrl,0,
C     .  'Embedded cluster order N scheme',-1,0,0,TF(23))
      nchoos = 0
#endif
      goto 1000

C --- Species (old CLASS) ---
C Inputs parameters from control file and optionally moments and pp's
C from disk.  Input must take the following order:
C (1) read atom label, and lmx from control file;
C (2) attempt to read princ. quant. #  CONF from control file;
C (3) attempt to read moments and pp's from atom file;
C (4) attempt to read idmod from control file;
C Step (1) must be first as this information is necessary for input
C to atom file; step (3) follows (2) to rank the priority of the
C atomic file higher wrt to CONF; step (4) follows (3) to rank the
C priority of the control file higher wrt IDMOD.
C Token ordering: CLASS;
C (1-10):  ATOM Z (R,R/W,R/A) EREF (NR,A) (LMX,LMXA) (P,Q) RHAM RFREE IDXDN
C (11-20)  ALPHA HCR DV IDMOD GROUP GRP2 AMASS MXCST C-HOLE PZ
C (21-30)  PB1,PB2 LMXPB LMXL (BAS,RSMH,EH) RSMG RFOCA RSMFA KMXA,RSMA QPOL I
C (31-36)  (SCLWSR,OMAX1,OMAX2) WSRMAX FRZWF (LXI,EXI) RINT RCUT
C (37-42)  (COLOUR,RADIUS) (IDU,UH,JH) RCFA VSO, (FRZQ1,FRZSO) EHVL
   70 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading SPEC',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('1111000000000000000000000000000000222202220',toksw,3,TF)

      avw = dgets('lat avw',slat)
      if (avw .eq. 0) avw = 1
      lpzi = 1
      j = TF(0)
      if (optio .eq. 1 .and. j .eq. 1) j = 0
      call getcat(recrd,'SPEC ','For species information'//
     .  ' (may also use CLASS)',j)
      if (.not. noerr) then
        call getcat(recrd,'CLASS ',' ',TF)
        if (.not. noerr) goto 1000
      endif

      if (optio .eq. 0) then
        print *,' '
        print *,
     .    '  ...    The next four tokens apply to the sphere resizer'
      endif

C     Global species input
      i0 = partks(recrd(catbeg),'SCLWSR=','= ,','ctrl sclwsr',sctrl,0,
     .  'SCLWSR>0 turns on sphere resizer and scales spheres trying'//
     .  '%N%10fto reach volume = SCLWSR * cell volume.'//
     .  '%N%10fAdd  10  to initially scale non-ES first;'//
     .  '%N%10f or  20  to scale ES independently.',
     .  1,4,0,TF(31))
      i0 = partks(recrd(catbeg),'OMAX1=','= ,','ctrl omax1',sctrl,0,
     .  'Sphere overlap constraints when adjusting MT radii',3,4,0,
     .  TF(31))
      i0 = partks(recrd(catbeg),'OMAX2=','= ,','ctrl omax2',sctrl,0,
     .  'Sphere overlap constraints of second type',3,4,0,TF(31))
      i0 = partks(recrd(catbeg),'WSRMAX=','= ,','ctrl wsrmax',sctrl,0,
     .  'Constrain maximum scaled sphere radius to WSRMAX',1,4,
     .  0,TF(32))

      if (optio .eq. 0) then
        print 382
  382   format(/'    ...  ',
     .    'The following tokens are input for each species.'/9x,
     .    'Data sandwiched between successive occurences of'/9x,
     .    'token ATOM= apply to one species.')
      endif

      jc = 0
      lmxax = 0
      nkaph = 1
      if (optio .eq. 2) nkaph = nglob('nkaph')
      do  71  j = 1, nspec
        jc = jc+1
        subsiz = catsiz
        i = partok(recrd(catbeg),'ATOM=','=',8,slabl(j),1,1,iend,TF(1))
C   ... i1 points to start of subcategory
        i1 = iend
        call subcat(recrd(catbeg),'ATOM=','=',i1)
        call upack('spec z rmt a nr',sspec,j,z,rmt,xx,nr)
        i0 = partok(recrd(catbeg),'Z=','=',z,'atomic number',-1,4,i1,
     .    TF(2))
C   ... Pick up R=, or R/W=, or R/A=
        nchoos = 3
        i0 = partok(recrd(catbeg),'R=','=',rmt,
     .    'rmax for augmentation',-1,4,i1,TF(3))
        if (.not. noerr .and. optio .ne. 2) then
          i0 = partok(recrd(catbeg),'R/W=','=',rmt,
     .      'rmax relative to average WS radius',-1,4,i1,TF(3))
          if (i0 .eq. 1) rmt = rmt*dgets('lat avw',slat)
        endif
        if (.not. noerr .and. optio .ne. 2) then
          i0 = partok(recrd(catbeg),'R/A=','=',rmt,
     .      'rmax relative to lattice constant',-1,4,i1,TF(3))
          if (i0 .eq. 1) rmt = rmt*dgets('lat alat',slat)
        endif
        nchoos = 0
C   ... Radial mesh parameters
        i0 = partok(recrd(catbeg),'NR=','=',nr,
     .    'number of radial mesh points',1,2,i1,TF(5))
        i0 = partok(recrd(catbeg),'A=','=',xx,
     .    'radial mesh point spacing',1,4,i1,TF(5))
        call pshpr(0)
        call rmesh(z,rmt,igets('ctrl lrel',sctrl),F,nrmx,xx,nr)
        call poppr
        call pack5('spec z rmt a nr',sspec,j,z,rmt,xx,nr)

        if (optio .eq. 1 .or. optio .eq. 0)
     .  call pack5('spec lmxa lmxb lmxf',sspec,j,nl-1,nl-1,2*nl-2,0)
        if (optio .eq. 1 .and. rmt .eq. 0 .and. TF(3) .ne. 2)
     .  call pack5('spec lmxa lmxb lmxf',sspec,j,-1,-1,-1,0)

        i0 = partks(recrd(catbeg),'LMX=','=','spec lmxb',sspec,j,
     .    'l-cutoff for basis',-1,2,i1,TF(6))
C       Flag that lmxb actually read for later checks
        sw = i0 .gt. 0
        lmxb = igetss('spec lmxb',j,sspec)
#ifndef LMMC
        call isanrg(lmxb,0,nl-1,' rdccat:','lmxb',.true.)
#endif
C       Floating orbitals have no lmxa
        if (rmt .eq. 0) then
        lmxa = -1
        else
        if (optio .eq. 1) call pack2('spec lmxa',sspec,j,lmxb)
        i0 = partks(recrd(catbeg),'LMXA=','=','spec lmxa',sspec,j,
     .    'l-cutoff for augmentation',-1,2,i1,TF(6))
        lmxa = igetss('spec lmxa',j,sspec)
        lmxax = max(lmxa,lmxax)
#ifndef LMMC
C       lmxb may not exceed lmxa
        if (sw) then
          if (lmxb .gt. lmxa) call fexit2(-1,1,' Exit -1 '//
     .      'species '//slabl(j)//'%a : LMX=%i exceeds LMXA=%i.'//
     .      '  Revise input so that LMX<=LMXA.',lmxb,lmxa)
        else
          lmxb = min(lmxb,lmxa)
          if (optio .eq. 1) call pack2('spec lmxb',sspec,j,lmxb)
        endif
#endif
        endif
        call pack5('spec lmxa lmxb',sspec,j,lmxa,lmxb,0,0)

C   ... Basis set for lmf
        i0 = partks(recrd(catbeg),'RSMH=','= ,','spec orbp;rsmh:1.10',
     .    sspec,j,'basis smoothing radii',lmxb+1,4,i1,TF(24))
        i0 = partks(recrd(catbeg),'EH=','= ,','spec orbp;eh:11.20',
     .    sspec,j,'basis kinetic energies',lmxb+1,4,i1,TF(24))
        i0 = partks(recrd(catbeg),'EHVL=','= ,','spec ehvl',
     .    sspec,j,'val-lap fit energies',lmxb+1,4,i1,TF(42))
        i0 = partks(recrd(catbeg),'RSMH2=','= ,',
     .    'spec orbp;rsmh2:21.30',sspec,j,
     .    'basis smoothing radii, second group',lmxb+1,4,i1,TF(24))
        i0 = partks(recrd(catbeg),'EH2=','= ,','spec orbp;eh:31.40',
     .    sspec,j,'basis kinetic energies, second group',lmxb+1,4,i1,
     .    TF(24))
        if (i0 .gt. 0 .and. optio .eq. 1) nkaph = 2
        i0 = partks(recrd(catbeg),'RS3=','= ,','spec rs3',sspec,j,
     .    'minimum smoothing radius for local orbital',1,4,i1,TF(24))
        i0 = partks(recrd(catbeg),'RHAM=','=','spec rham',sspec,j,
     .    'radius for range of hamiltonian',1,4,i1,TF(8))
C        i0 = partks(recrd(catbeg),'RFREE=','=','spec rsmfa',sspec,j,
C     .    'radius for ''FREE=T'' atomic calculation',1,4,i1,TF(9))
        i0 = partks(recrd(catbeg),'IDXDN=','= ,','spec idxdn',sspec,j,
     .    'downfolding index: 0, auto; 1, no dnf; 2, fold down;'//
     .    ' 3, neglect',1+lmxb,2,i1,TF(10))
        nchoos = 2
        i0 = partkF(recrd(catbeg),'HCR=','= ,','spec hcr',sspec,j,
     .    'Hard sphere radii for structure constants',
     .    lmxb+1,104,i1,TF(12),1d0)
        if (.not. noerr) then
        if (rmt .eq. 0) rmt = 1
        i0 = partkF(recrd(catbeg),'HCR/R=','= ,','spec hcr',sspec,j,
     .    'Hard sphere radii for structure constants, units of Rmax',
     .    lmxb+1,104,i1,TF(12),rmt)
        endif
C       Nothing read in ... set hcr to default
        if (.not. noerr .and. optio .eq. 1) then
          call dvset(xv,1,n0,.7d0*rmt)
          call pack5('spec hcr',sspec,j,xv,0,0,0)
        endif
        nchoos = 0
        i0 = partks(recrd(catbeg),'ALPHA=','= ,','spec alpha',sspec,j,
     .    'screening parameters for structure constants',-lmxb-1,4,i1,
     .    TF(11))
C       Find l-cutoff for TB alphas, to choose appropriate alpha set
C       (doesn't belong here ...)
        lp1 = lmxb+1
        if (optio .eq. 1 .and. TF(10) .ne. 2) then
          call upack2('spec idxdn',sspec,j,ivec)
   75     if (lp1 .gt. 1 .and. ivec(lp1) .ge. 4) then
            lp1 = lp1-1
            goto 75
          endif
        endif
        if (i0 .eq. 0) then
          call dpzero(xv,n0)
          call dcopy(5,tbalph(1,lp1),1,xv,1)
          call pack2('spec alpha',sspec,j,xv)
        endif
        i0 = partks(recrd(catbeg),'DV=','=','spec dv',sspec,j,
     .    'extra constant potential shift for this species',1,4,i1,
     .    TF(13))
        i0 = partks(recrd(catbeg),'IDMOD=','= ,','spec idmod',sspec,j,
     .    'idmod=0 floats P to band CG, 1 freezes P, 2 freezes enu',
     .    lmxa+1,2,i1,TF(14))
        i0 = partks(recrd(catbeg),'GROUP=','=','spec group',sspec,j,' ',
     .    1,2,i1,TF(15))
        i0 = partks(recrd(catbeg),'GRP2=', '=','spec grp2',sspec,j,' ',
     .    1,2,i1,TF(16))

        i0 = partks(recrd(catbeg),'MIX=','=','spec mxcst,1',sspec,j,
     .    '(ASA) set nonzero to suppress self-consistency of classes in'
     .    //' this spec',1,2,i1,TF(18))

        i = bitand(igetss('spec mxcst',j,sspec),2)
        sw = i .ne. 0
        i0 = partok(recrd(catbeg),'CSTRMX=','=',sw,'Exclude this'
     .    //' species when automatically resizing sphere radii'//
     .    ' (SCLWSR>0)',1,0,i1,TF(32))
        call pack2('spec mxcst',sspec,j,
     .    igetss('spec mxcst',j,sspec)-i+2*isw(sw))

        i = bitand(igetss('spec mxcst',j,sspec),4)
        sw = i .ne. 0
        i0 = partok(recrd(catbeg),'FRZWF=','=',sw,'Freeze augmentation'
     .    //' wave functions for this species',1,0,i1,TF(33))
        call pack2('spec mxcst',sspec,j,
     .    igetss('spec mxcst',j,sspec)-i+4*isw(sw))

        if (optio .eq. 0 .and. TF(38) .lt. 2) then
          print *,' '
          print *, '  ...    The next three tokens are for LDA+U'
        endif

        i0 = partks(recrd(catbeg),'IDU=','= ,','spec idu',sspec,j,
     .    'idu=0,1,2,3 LDA+U, 0 nothing, 1 AMF, 2 FLL, 3 mixed ',
     .    4,2,i1,TF(38))

        i0 = partks(recrd(catbeg),'UH=','= ,','spec uh',sspec,j,
     .    'Hubbard U for LDA+U',
     .    4,4,i1,TF(38))

        i0 = partks(recrd(catbeg),'JH=','= ,','spec jh',sspec,j,
     .    'Exchange parameter J for LDA+U',
     .    4,4,i1,TF(38))

C   ... linear-response-specific
        if (optio .eq. 1 .and. z .eq. 0) then
          call spacks(1,'spec pb1',sspec,'        ',j,j)
          call spacks(1,'spec pb2',sspec,'        ',j,j)
        endif
        i0 = partks(recrd(catbeg),'PBAS=','=','spec pb1',sspec,j,
     .    'product basis for GW',1,1,i1,TF(21))
        i0 = partks(recrd(catbeg),'PBAS2=','=','spec pb2',sspec,j,
     .    'second product basis for GW',1,1,i1,TF(21))
        i0 = partks(recrd(catbeg),'LMXPB=','=','spec lmxpb',sspec,j,
     .    'l cutoff for product basis',1,2,i1,TF(22))

C   ... tbe specific
        i0 = partks(recrd(catbeg),'QPOL=','=','spec qpol',sspec,j,
     .    '(tbe) up to ten polarisability parameters',10,4,i1,TF(29))
        i0 = partks(recrd(catbeg),'I=','=','spec stni',sspec,j,
     .    'TB Stoner I (l=2)',1,4,i1,TF(30))
        i0 = partks(recrd(catbeg),'FRZQ1=','=','spec iq1',sspec,j,
     .    '(tbfit) switches to freeze site energies',nl*nsp,2,i1,TF(41))
        if (lgors('ctrl lncol,4',sctrl) .and. TF(40) .ne. 2
     .      .and. lmxa .ge. 1 .or. optio .eq. 0) then
          i0 = partks(recrd(catbeg),'VSO=','=','spec vso:2',sspec,j,
     .      'TB spin-orbit parameters',-lmxb,4,i1,TF(40))
          i0 = partks(recrd(catbeg),'FRZSO=','=','spec ivso:2',sspec,j,
     .      '(tbfit) switches to freeze TB spin-orbit parameters',
     .      lmxa,2,i1,TF(41))
        endif

#if FTB
C    ... Basis
        i0 = partks(recrd(catbeg),'BAS:',':, ','spec ntorb orbp',
     .    sspec,j,'Defines basis:'//
     .    '%N%12fntorb=number of orbital types'//
     .    '%N%12forbp=orbital parameters lmin,lmax,rsm,e for each type',
     .    -1,5,i1,TF(24))
#endif

C   ... File read initial sphere boundary conditions P,Q
        i0 = partks(recrd(catbeg),'P=','= ,','spec p',sspec,j,
     .  'Starting log der. parameters for each l',1+lmxa,4,i1,TF(7))
        i0 = partks(recrd(catbeg),'PZ=','= ,','spec pz',sspec,j,
     .  'Starting semicore log der. parameters'//
     .    '%N%10fAdd 10 to attach Hankel tail'//
     .    '%N%10fAdd 20 to include perturbatively',
     .    1+lmxa,4,i1,TF(20))
        if (i0 .gt. 0) then
C         print *, j,lmxa,slabl(j)
          call upack('spec pz',sspec,j,wk,0,0,0)
          if (dasum(lmxb+1,wk,1) .ne. 0) lpzi = max(lpzi,2)
        endif
        i0 = partks(recrd(catbeg),'Q=','= ,','spec q',sspec,j,
     .  'Starting sphere charges for each l channel',1+lmxa,4,i1,TF(7))

C   ... Set up default P,Q in absence of explicit specification
        nsp = igets('ctrl nspin',sctrl)
        if (optio .eq. 1) then
C         Default P
          call dpzero(tmp,2*n0)
          call dpzero(xv, 2*n0)
          call upack('spec p',sspec,j,tmp,0,0,0)
          call pshpr(iprint()-30)
          call defpq(z,lmxa,nsp,tmp,xv)
          call poppr
          call pack2('spec p',sspec,j,tmp)
C         Default q
          if (i0 .eq. 0) then
            call dpzero(xv,2*n0)
            call pshpr(iprint()-30)
            call defpq(z,lmxa,1,tmp,xv)
            call poppr
            call pack2('spec q',sspec,j,xv)
          endif
        endif
        if (nsp .eq. 2 .or. optio .eq. 0)
     .    i0 = partks(recrd(catbeg),'MMOM=','= ,','spec q:11',sspec,j,
     .    'Starting mag. moms for each l channel',1+lmxa,4,i1,TF(7))

C   ... Defaults for:
C       lmxl = l-cutoff for numerical rep'sn of density in sphere
C       lfoca = mode for treating core
C       kmxt = kmax for expansion of envelope wf tails
C       kmxv = cutoff to expand smoothed potential
        if (optio .le. 1) then
          i0 = 1
          if (z .le. 8) i0 = 0
          call pack5('spec lmxl lfoca kmxt kmxv',sspec,j,lmxa,i0,3,15)
C         rdccat no longer sets defaults: need set after rescaling of rmt
C         call pack2('spec rsmv',sspec,j,rmt*.5d0)
        endif

        if (lmxa .ne. -1) then
        i0 = partks(recrd(catbeg),'LMXL=','=','spec lmxl',sspec,j,
     .    'lmax for which to accumulate rho,V in sphere',1,2,i1,TF(23))
        endif
        i0 = partks(recrd(catbeg),'RSMG=','=','spec rg',sspec,j,
     .    'smoothing for projection of charge in sphere.'//
     .    '%N%10finput<0 => choose default * -input',-1,4,i1,TF(25))

        i0 = partks(recrd(catbeg),'KMXA=','=','spec kmxt',sspec,j,
     .    'k-cutoff for projection of wave functions in sphere.',1,2,i1,
     .    TF(28))
        i0 = partks(recrd(catbeg),'RSMA=','=','spec rsma',sspec,j,
     .    'smoothing for projection of wave functions in sphere.'//
     .    '%N%10finput<0 => choose default * -input',-1,4,i1,TF(28))

        i0 = partks(recrd(catbeg),'LFOCA=','=','spec lfoca',sspec,j,
     .    'switch 0, 1 or 2 for FOCA',-1,2,i1,TF(26))
        i0 = partks(recrd(catbeg),'RFOCA=','=','spec rfoca',sspec,j,
     .    'smoothing for core tail.  input<0 => choose default * -input'
     .    ,-1,4,i1,TF(26))

        i0 = partks(recrd(catbeg),'RSMFA=','=','spec rsmfa',sspec,j,
     .    'smoothing for free atom.  input<0 => choose default * -input'
     .    ,-1,4,i1,TF(27))

        i0 = partks(recrd(catbeg),'RCFA=','= ,','spec rcfa',sspec,j,
     .    'Cutoff radius for renormalization of free atom density.'//
     .    '%N%10fOptional 2nd argument = width'//
     .    '%N%10fRCFA<0 => renormalize potential instead of density',
     .    2,4,i1,TF(39))

C   ... molecules-specific.  lmxb will change, so use TF(34) for now
        i0 = partks(recrd(catbeg),'LMXB=','=','spec lmxpb',sspec,j,
     .    'compound l-cutoffs for basis (molecules)',-1,2,i1,TF(34))
        i0 = partks(recrd(catbeg),'LXI=','=','spec lxi',sspec,j,
     .    'l-cutoffs for interstitial product basis',-1,2,i1,TF(34))
        i = igetss('spec lxi',j,sspec)
        call lx2vec(i,0,ii,ivec)
        i0 = partks(recrd(catbeg),'EXI=','=','spec exi',sspec,j,
     .    'Hankel energies for interstitial product basis',-ii,4,i1,
     .    TF(34))
        i0 = partks(recrd(catbeg),'RINT=','=','spec rint',sspec,j,
     .    'range of interstitial product basis',-1,4,i1,TF(35))
        i0 = partks(recrd(catbeg),'RCUT=','=','spec rcut',sspec,j,
     .    'radius over which true,smooth Hankels differ',-1,4,i1,TF(36))

C       rdccat no longer sets defaults: need set after rescaling of rmt
C       Default values of rg,rsma,rfoca,rsmfa
C        call upack('spec rg rsma rfoca rsmfa',sspec,j,rg,rsma,rfoca,
C     .    rsmfa)
C        if (rg    .eq. 0) rg    = -1
C        if (rg    .lt. 0) rg    = -rg*0.25d0*rmt
C        if (rsma  .eq. 0) rsma  = -1
C        if (rsma  .lt. 0) rsma  = -rsma*0.4d0*rmt
C        if (rfoca .eq. 0) rfoca = -1
C        if (rfoca .lt. 0) rfoca = -rfoca*0.4d0*rmt
C        if (rsmfa .eq. 0) rsmfa = -1
C        if (rsmfa .lt. 0) rsmfa = -rsmfa*0.5d0*rmt
C        call pack5('spec rg rsma rfoca rsmfa',sspec,j,rg,rsma,rfoca,
C     .    rsmfa)

C       Some sanity checks
#ifndef LMMC
        call upack('spec lmxl lfoca kmxt',sspec,j,lmxl,lfoca,kmxt,0)
        call isanrg(lmxl,min(0,lmxa),max(0,lmxa),'rdccat','lmxl',T)
        call isanrg(lfoca,0,2,'rdccat','lfoca',T)
        call isanrg(kmxt,2,25,'rdccat','kmxt',T)
#endif

        if (optio .eq. 0 .and. TF(19) .lt. 2) then
          print *,' '
          print *, '  ...    The next two tokens are for core hole'
        endif
        i0 = partks(recrd(catbeg),'C-HOLE=','=','spec coreh',sspec,j,
     .    'channel for core hole',1,1,i1,TF(19))
        i0 = partks(recrd(catbeg),'C-HQ=','= ,','spec coreq',sspec,j,
     .    'Charge in core hole.  '//
     .    'Optional 2nd entry is moment of core hole:'//
     .    '%N%12fQ(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2'//
     .    '%N%12fQ(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2',2,4,i1,TF(19))

        i0 = partks(recrd(catbeg),'EREF=','=','spec eref',sspec,j,
     .    'Reference energy subtracted from total energy',1,4,i1,TF(4))

        nchoos = 2
        call upack('spec mass',sspec,j,xx,0,0,0)
        noerr = .false.
        if (xx .eq. 0) then
          i0 = partok(recrd(catbeg),'ISOTOPE=','=',xx,
     .      'relative atomic mass (for dynamics)',-1,4,i1,
     .      TF(17))
          if (i0 .eq. 1) then
            ncor = 0
            i0 = partok(recrd(catbeg),'NCOR=','=',ncor,' ',-1,2,i1,F)
C       ... Work out atomic mass from neutron, proton and electron masses
            nmass = 1836.153d0/2d0
            pmass = 1838.684d0/2d0
            emass = 0.5d0
            xx = z*pmass + (xx-z)*nmass + ncor*emass
            call pack5('spec mass',sspec,j,xx,0,0,0)
            noerr = .true.
          endif
        endif
        if (.not. noerr) then
          i0 = partks(recrd(catbeg),'AMASS=','=','spec mass',sspec,j,
     .      'nuclear mass in a.u. (for dynamics)',1,4,i1,TF(17))
        endif
        nchoos = 0

        i0 = partks(recrd(catbeg),'COLOUR=','=','spec colxbs',sspec,j,
     .    'colour for xbs',-3,4,i1,TF(37))
        i0 = partks(recrd(catbeg),'RADIUS=','=','spec radxbs',sspec,j,
     .    'radius for xbs',-1,4,i1,TF(37))

        iend = i1
   71 continue
C ... Get the maximum L-cutoff set global variable nkaph and mxorb
      nlmax = (lmxax+1)**2
* fix nkaph by redesign of uspecb
      if (optio .eq. 1) then
        nkaph = nkaph+lpzi-1
        xx = dglob('nkaph',dble(nkaph),1)
        xx = dglob('nlmax',dble(nlmax),1)
        xx = dglob('mxorb',dble(nkaph)*nlmax,1)
      endif
C ... Copy species labels to array
      call spacks(1,'spec name',sspec,slabl,1,nspec)
      if (nl .lt. lmxax+1 .and. iprint() .gt. 5)
     .  call awrit2(' rdccat (warning):  default NL(=%i)'//
     .  ' less than max LMXA+1(=%i)',' ',80,lgunit(1),nl,lmxax+1)
C ... Update nl
      if (optio .eq. 1) then
        nl = lmxax+1
        call pack1('ctrl nl',sctrl,nl)
      endif

C     jc = jc+1
C      if (jc .lt. nspec) then
C        if (iprint().gt.30) call awrit1('rdccat: retaining %i species',
C     .    ' ',recl0,i1mach(2),nspec)
C        nspec = jc
C        stop
C      endif
C      if (optio .eq. 2) call shstru('spec',sspec,1,nspec)
C      pause
      goto 1000

C --- Site ---
C Token ordering: SITE;
C (1-12) POS ROT RELAX DELTA PL PLV VSHFT DPOS MODE SID FILE BFIELD
C (13-14)   V0 PM
   80 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading SITE',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('11000000000022',toksw,3,TF)
      call getcat(recrd,'SITE ','Site-specific data',TF)
      if (.not. noerr) goto 1000
      if (TF(11) .ne. 2) then
      if (optio .eq. 0) print *,'  ...    ',
     .  'You may read SITE category from site file with following token'
      i0 = 0
      if (optio .lt. 2)
     .  i0 = partok(recrd(catbeg),'FILE=','=',8,outs,1,1,0,TF(11))
      if (i0 .eq. 1) then
        ii = (8+32)*1000 + 140
        if (tf(4) .lt. 2) ii = ii + (64)*1000
        call upack1('lat plat',slat,plat)
        j = iosits(ii,3d0,0,outs(1:8),j,slabl,alat,plat,nbas,nspec,
     .    sspec,ssite)
        goto 89
      endif
      if (optio .eq. 0) print *, '  ...    ',
     .  'Otherwise, read it from the ctrl file:'
      endif

      sw = F
      i0 = partok(recrd(catbeg),'MODE=','=',sw,'Set to T if input '//
     .  'positions are fractions of plat',1,0,0,TF(9))
C ... Site pos associated quantities for all basis atoms
      call spack(2,'site spec',ssite,1,offi,casti,nelti,i)
      do  81  j = 1, nbas
        jm1 = j-1
        subsiz = catsiz
        if (optio .eq. 2) then
          i = ssite(offi(1)+nint(ssite(1))*jm1)
          alabl = slabl(i)
        endif
        i = partok(recrd(catbeg),'ATOM=','=',8,alabl,1,1,iend,TF)
        i1 = iend
        call subcat(recrd(catbeg),'ATOM=','=',i1)
        if (optio .eq .1) then
          call tokmat(alabl,slabl,nspec,8,' ',i,i0,F)
          if (i .lt. 0) call rxs('rdctrl: Category SITE referred to'//
     .      ' nonexistent species: ',alabl)
          i = i+1
          ssite(offi(1)+nint(ssite(1))*jm1) = i
        endif
        i0 = partks(recrd(catbeg),'POS=','=','site pos',ssite,j,
     .    'Atom coordinates, in units of alat',-3,4,i1,TF(1))
C       If positions are input as fractions of plat, scale to Cartesian
        if (sw .and. optio .eq. 1) then
          call spack(2,'site pos',ssite,1,ivec,tmp(1),tmp(2),tmp(3))
          i = ivec(1)+nint(ssite(1))*jm1
          call dcopy(3,ssite(i),1,tmp,1)
          call upack1('lat plat',slat,plat)
          do  84  ii = 1, 3
   84     ssite(i+ii-1)= ddot(3,tmp,1,plat(ii,1),3)
        endif
        if (optio .ne. 2) then
          i0 = partok(recrd(catbeg),'DPOS=','=',tmp,
     .      'shift in atom coordinates added to pos',-3,4,i1,TF(8))
          if (i0 .eq. 3) then
            call upack2('site pos',ssite,j,tmp(4))
            do  82  i = 1, 3
   82       tmp(i) = tmp(i)+tmp(3+i)
            call pack2('site pos',ssite,j,tmp)
          endif
        endif
        i0 = partks(recrd(catbeg),'V0=','=','site vel',ssite,j,
     .    'initial velocity for molecular dynamics',-3,4,i1,TF(13))

        if (optio .ne. 2) i0 =
     .    partok(recrd(catbeg),'ROT=','=',len(outs),outs,1,1,i1,TF(2))
        if (i0 .eq. 1 .and. optio .ne. 2) then
          call numsyv(i0)
          call upack2('site pos',ssite,j,xv)
          call lodsyv('ib',0,dble(j+1),ii)
          call lodsyv('x', 0,xv,ii)
          call lodsyv('y', 0,xv(2),ii)
          call lodsyv('z', 0,xv(3),ii)
          call a2rotm(outs,F,iprint()-10,tmp)
          call clrsyv(i0)
          call rm2eua(tmp,xv,xv(2),xv(3))
          call pack2('site eula',ssite,j,xv)
        endif
        if (optio .eq. 0 .and. TF(2) .ne. 2) then
          print 346
  346     format(10x,'(non-collinear fields): rotation of local ',
     .      'coordinates.  Syntax:'/10x,'ROT=rot1[,rot2...] with rotj=',
     .      '(x,y,z)angle.  (x,y,z) may be x: y: or z:')
        elseif (optio .eq. 2 .and. TF(2) .ne. 2) then
          call upack2('site eula',ssite,j,xv)
          call awrit3('%4pROT= : Euler alpha = %1;6d  beta = %1;6d'//
     .      '  gamma = %1;6d',' ',80,i1mach(2),xv,xv(2),xv(3))
        endif
        i0 = partks(recrd(catbeg),'RELAX=','=','site relax',ssite,j,
     .    'relax site positions (lattice dynamics) or '//
     .    'Euler angles (spin dynamics)',3,2,i1,TF(3))
        i0 = partks(recrd(catbeg),'DELTA=','=','site delta',ssite,j,
     .    '(tbe) vector of on-site energy shifts',nl*nsp,4,i1,TF(4))
        call pack2('site ndelta',ssite,j,nl*nsp)
        i0 = partks(recrd(catbeg),'PL=','=','site pl',ssite,j,'(pgf)'//
     .    ' Assign principal layer number to this site',1,2,i1,TF(5))
        i0 = partks(recrd(catbeg),'PLV=','=','site plv',ssite,j,
     .    '(pgf) Assign PL potential index to this site',1,2,i1,TF(6))
        i0 = partks(recrd(catbeg),'VSHFT=','=','site vshft',ssite,j,
     .    'Constant potential shift for this site',1,2,i1,TF(7))
        i0 = partks(recrd(catbeg),'SID=','=','site sid',ssite,j,'(rgf)'
     .    //' Site ID, used in the merging of coincident clusters',1,2,
     .    i1,TF(10))
        iend = i1

   81 continue

      call upack1('ctrl nbasp',sctrl,nbasp)
      if (nbasp .gt. nbas .and. TF(14) .lt. 2) then
        iend = 5
        do  86  j = nbas+1, nbasp
          subsiz = catsiz
          i0 = partok(recrd(catbeg),'PM ',' ',0,' ',0,1,iend,TF(14))
          i1 = iend
          call subcat(recrd(catbeg),'PM ',' ',i1)
          i0 = partks(recrd(catbeg),'POS=','= ,','site pos',ssite,j,
     .      'Point charge coordinates, in units of alat',-3,4,i1,TF(14))
          i0 = partks(recrd(catbeg),'Q=','=','site mpole',ssite,j,
     .      'Point monopole moment',-1,4,i1,TF(14))
          i0 = partks(recrd(catbeg),'P=','= ,','site dpole',ssite,j,
     .      'Point dipole moment',-3,4,i1,TF(14))
          iend = i1
   86   continue
      endif
   89 continue
      if (optio .eq. 1 .and. cmdopt('--rpos=',7,0,outs)) then
        call upack1('ctrl nbasp',sctrl,nbasp)
        call ioposs(.false.,-1,outs(8:),nbasp,ssite)
      endif
      goto 1000

C --- Symmetry group ---
C Token ordering: SYMGRP
   90 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading SYMGRP',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0',toksw,3,TF)
      call getcat(recrd,'SYMGRP ',
     .  'for explicit specification of symmetry group',TF)
      if (.not. noerr) goto 1000
      if (optio .eq. 2) then
        call upacks('strn symg',i,j)
        if (j .lt. i) goto 1000
      endif
      i0 = prtkss(recrd(catbeg),'SYMGRP ',' ','strn symg',sstrn,catsiz,
     .  'Generators for symmetry group',0,TF(0)+200)
C      if (i0 .eq. 0 .and. optio .eq. 1) then
C        call lstra('strn symg',tmp,offs,i)
C        call ustrn(tmp,-offs,1,offi,tmp,1)
C        sstrn(offi(1):offi(1)) = ' '
C      endif
      goto 1000

C --- Structure constants ---
C Token ordering: STR;
C (1-9) RMAX MODE (NEL,EL,EW2) MXNBR SHOW IINV EQUIV LMAXW,DRWATS
C (10-15) NBISI NALF NCUPL NDUST ADEC WZTCF (These 6 for TCF)
C (16)  RMAXG,IVL
  100 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading STR',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('00000000002222222',toksw,3,TF)
C               0123456789012345
      call getcat(recrd,'STR ','Real-space structure constants',TF)
      if (.not. noerr) goto 1000
      avw = dgets('lat avw',slat)
      if (avw .eq. 0) avw = 1
C ... Pick up RMAXS=, or RMAX=, or RMAX/A=
      nchoos = 2
      i0 = partkf(recrd(catbeg),'RMAXS=','=','str rmax',sstr,0,
     .  'Radial cutoff for strux, in a.u.',-1,4,0,TF(1),1d0)
      if (.not. noerr) then
      i0 = partkf(recrd(catbeg),'RMAX=','=','str rmax',sstr,0,
     .  'Radial cutoff for strux, in units of avw',-1,4,
     .  0,TF(1),avw)
      endif

      i0 = partks(recrd(catbeg),'RVL/R=','=','str rmaxg',sstr,0,
     .  'Radial cutoff for val-lap, in units of RMAX',-1,4,0,TF(16))
      i0 = partks(recrd(catbeg),'VLFUN=','=','str ivl',sstr,0,
     .  'Functions for val-lap basis'//
     .  '%N%10f0 G0 + G1'//
     .  '%N%10f1 G0 + Hsm'//
     .  '%N%10f2 G0 + Hsm-dot',1,2,0,TF(16))
      nchoos = 0
      i0 = partks(recrd(catbeg),'MODE=','=','str amode',sstr,0,
     .  '0 for 2nd gen (def); 1 for SSSW; 2 for NMTO; '//
     .  '3 SSSW and val-lap basis',1,2,0,TF(2))
C     For now, turn off OKA conventions for SSSW
c     if (igets('str amode',sstr) .eq. 1) then
      if (igets('str amode',sstr) .eq. 1 .or. 
     .  igets('str amode',sstr) .eq. 3) then
        call pack2('str loka',sstr,0,0)
      endif
      i0 = partks(recrd(catbeg),'NEL=','=','str nkaps',sstr,0,
     .  'Number of SSSW or NMTO energies',-1,2,0,TF(3))
      n = igets('str nkaps',sstr)
C ... Pick up LMTO energies for structure constants
      nchoos = 2
      i0 = partkf(recrd(catbeg),'EL=','=','str kaps',sstr,0,
     .  'SSSW or NMTO energies, units avw^-2',-n,4,0,TF(3),1d0)
      if (.not. noerr) then
      i0 = partkf(recrd(catbeg),'EW2=','=','str kaps',sstr,0,
     .  'SSSW or NMTO energies, units avw^-2',-n,4,0,TF(3),1/avw**2)
      endif
      nchoos = 0
      i0 = partks(recrd(catbeg),'MXNBR=','=','str mxnbr',sstr,0,
     .  'Max # of nbrs: for dimensioning arrays ',1,2,0,TF(4))
      i0 = partks(recrd(catbeg),'SHOW=','=','str lshow,1',sstr,0,' ',1,
     .  0,0,TF(5))
      i0 = partks(recrd(catbeg),'EQUIV=','=','str lequiv,1',sstr,0,' ',
     .  1,0,0,TF(7))

      i0 = partks(recrd(catbeg),'IINV:',':, ','str iinv;ncut:1.1'//
     .  ' iinv;nit:2.2 iinv;tol:3.3',sstr,0,
     .  'Iterative inversion of screened strux:'//
     .  '%N%12fncut=number of sites for inner block'//
     .  '%N%12fnit=number of iterations'//
     .  '%N%12ftol=tolerance in errors',
     .  -1,5,0,TF(6))

      i0 = partks(recrd(catbeg),'LMAXW=','=','str lmaxw',sstr,0,
     .  'l-cutoff for Watson sphere, used to help localize strux',
     .  -1,2,0,TF(8))
      i0 = partks(recrd(catbeg),'DELRW=','=','str drwats',sstr,0,
     .  'padding beyond cluster for Watson sphere',-1,4,
     .  0,TF(8))
      if (optio .eq. 0) call info(0,1,0,
     .  '   ... The following tokens are used for the 2-center fit'//
     .  '(TCF)',0,0)
      i0 = partks(recrd(catbeg),'NBISI=','=','str nbisi',sstr,0,
     .  'mesh parameters for TCF ',-3,2,0,TF(10))
      i0 = partks(recrd(catbeg),'NALF=','=','str nalf',sstr,0,
     .  'polynomial order for d-dependence of TCF fit coffs ',
     .  1,2,0,TF(11))
      i0 = partks(recrd(catbeg),'NCUPL=','=','str ncupl',sstr,0,
     .  'reduced polynomial order in TCF ',1,2,0,TF(12))
      i0 = partks(recrd(catbeg),'NDUST=','=','str ndust',sstr,0,
     .  'no. iterations to improve on TCF ',1,2,0,TF(13))
      i0 = partks(recrd(catbeg),'ADEC=','=','str adec',sstr,0,
     .  'log spacing in distance for TCF ',1,4,0,TF(14))
      i0 = partks(recrd(catbeg),'WZTCF=','=','str wztcf',sstr,0,
     .  'weighting of z-points in TCF ',1,4,0,TF(14))

      goto 1000

C --- GW ---
C Token ordering: GW;
C (1-10)  GCUTB GCUTX QOFFP NKABC NBAND NIME DELRE DELTA DELTAW PBTOL
C (11-18) LMXPB MKSIG FAST GSMEAR ECUTS
  110 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading GW',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('00000000000000000000000000000000000',toksw,3,TF)
      call getcat(recrd,'GW ',' ',TF)
      if (.not. noerr) goto 1000
      if (optio .eq. 0) call info(0,0,0,
     .  '   ... The following tokens are used for GW job=-1',0,0)
      i0 = partks(recrd(catbeg),'MKSIG=','=','gw mksig',sgw,0,
     .  'Make full sigma(k) for self-consistent calculations:'//
     .  '%N%10f0 do not create self-energy matrix sigma'//
     .  '%N%10f1 Make sigma_nn''(E_f)'//
     .  '%N%10f2 Make sigma_nn''(E_n/2+E_n''/2)'//
     .  '%N%10f3 Make sigma_nn''(E_n)/2 + sigma_nn''(E_n'')/2',
     .  1,2,0,TF(12))
      call pvcat1('NKABC=','No. qp divisions '//
     .  'along each of 3 lattice vectors','gw nkabc',sgw,TF(4))
      i0 = partks(recrd(catbeg),'GCUTB=','=','gw gcutb',sgw,0,
     .  'G-vector cutoff for basis functions',1,4,0,TF(1))
      i0 = partks(recrd(catbeg),'GCUTX=','=','gw gcutx',sgw,0,
     .  'G-vector cutoff for response function',1,4,0,TF(2))
      i0 = partks(recrd(catbeg),'ECUTS=','=','gw ecuts',sgw,0,
     .  '(self-consistency) max energy for which to calc. sigma',1,4,0,
     .  TF(18))
C      i0 = partks(recrd(catbeg),'FAST=','=','gw lgw,2',sgw,0,
C     .  'Make response function and sigma using Faleev''s '//
C     .  '%N%10ffast energy integration',
C     .  1,2,0,TF(13))
      i0 = partks(recrd(catbeg),'NBAND=','=','gw nband',sgw,0,
     .  'No. bands to include in calc. response function',1,2,0,TF(5))
      i0 = partks(recrd(catbeg),'NIME=','=','gw nime',sgw,0,
     .  'Number of imaginary energy points',1,2,0,TF(6))
      i0 = partks(recrd(catbeg),'DELRE=','=','gw delre',sgw,0,
     .  'Energy mesh spacing on real axis.'//
     .  '%N%10fSecond term is coff for nonlinear part',2,4,0,TF(7))
      i0 = partks(recrd(catbeg),'DELTA=','=','gw deltax',sgw,0,
     .  'Imaginary energy added to denominator in integration of x0.'//
     .  '%N%10fNot used in Faleev''s fast integration.'//
     .  '%N%10fNB: Delta should be negative.',1,4,0,TF(8))
      i0 = partks(recrd(catbeg),'DELTAW=','=','gw deltaw',sgw,0,
     . 'Width for finite difference in diff. of sigma for Z factor',
     .  1,4,0,TF(9))
      i0 = partks(recrd(catbeg),'GSMEAR=','=','gw gsmear',sgw,0,
     . 'Broadening in pole of GF for computation of sigma',
     .  1,4,0,TF(14))
      i0 = partks(recrd(catbeg),'PBTOL=','=','gw pbtol',sgw,0,
     .  'retain prod-basis with overlap eigenvalues>PBTOL',1,4,0,TF(10))
      i0 = partks(recrd(catbeg),'QOFFP=','=','gw qoffp',sgw,0,'k-point'
     .  //' offset parameter for special BZ integration for xi',1,4,0,
     .  TF(3))
C      i0 = partks(recrd(catbeg),'LMXPB=','=','gw lmxpb',sgw,0,
C     .  'default l-cutoff in product basis',1,2,0,TF(11))
      goto 1000

C --- Brillouin Zone ---
C Token ordering: BZ;
C (1-10)  NKABC GETQP PUTQP TETRA BZJOB METAL DOSWT N W RANGE
C (11-17) NPTS SAVDOS DOS EFMAX NEVMX ZVAL NOINV
C (18-25) (COND,V) INVIT IO FSMOM MULL S_EMSH DMAT (EF0,DEF)
  120 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading BZ',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('11000000000000000000000200',toksw,3,TF)
      call getcat(recrd,'BZ ',
     .  'for Brillouin zone and energy integrations',TF)
      if (.not. noerr) goto 1000
      nchoos = 2
      i0 = partks(recrd(catbeg),'GETQP=','=',
     .  'bz lio,1',sbz,0,'Read qp from disk',1,0,0,TF(2))
      if (.not. noerr .or. .not. lgors('bz lio,1',sbz))
     .  call pvcat1('NKABC=','No. qp along each of 3 lattice vectors',
     .  'bz nkabc',sbz,TF(1))
      nchoos = 0
      call pvcat1('BZJOB=','0 centers BZ mesh at origin, 1 centers '
     .  //'off origin','bz lshft',sbz,TF(5))
      i0 = partks(recrd(catbeg),'METAL=','=','bz lmet',sbz,0,
     . '0 assume insulator; 1 save evecs on disk; '//
     . '2 use wgt from prior iter%N%10f3 always make two band passes;'//
     . ' 4 BZ integration with 3-point scheme',1,2,0,TF(6))
C ... set logical metal in ctrl structure
      call lsets('ctrl lmet',sctrl,isw(lgors('bz lmet,-1',sbz)),1)
      i0 = partks(recrd(catbeg),'TETRA=','=','ctrl lmet,2',sctrl,0,
     .  'Tetrahedron integration',1,0,0,TF(4))
      i0 = partks(recrd(catbeg),'N=','=','bz n',sbz,0,
     .  'N>0: Polynomial order for Methfessel-Paxton sampling%N%10f'//
     .  'N=0: Conventional Gaussian sampling%N%10f'//
     .  'N<0: Broadening by Fermi-Dirac distribution%N%10f'//
     .  'To be used in conjunction with W= ; see next',
     .  1,2,0,TF(8))
      nchoos = 2
      i0 = partks(recrd(catbeg),'W=','=','bz w',sbz,0,
     .  'N>=0: Line broadening for sampling integration%N%10f'//
     .  'N<0 : Temperature for Fermi distribution (Ry)',1,4,0,TF(9))
      if (.not. noerr .and. optio .le. 1) then
        i0 = partks(recrd(catbeg),'N.W=','=','bz w',sbz,0,
     .  'Specify N and W simultaneously by this token',1,4,0,TF(9))
        call upack2('bz n w',sbz,i,dum)
        if (dum .lt. -1) then
          i = -int(-dum)
          dum = -dum+i
        else
          i = int(dum)
          dum = dum-i
        endif
        call pack2('bz n w',sbz,i,dum)
      endif
      i0 = partks(recrd(catbeg),'EF0=','=','bz ef',sbz,0,
     .  'Initial guess at Fermi energy',1,4,0,TF(25))
      i0 = partks(recrd(catbeg),'DELEF=','=','bz def',sbz,0,
     .  'Initial uncertainty in Fermi energy',1,4,0,TF(25))
      i0 = partks(recrd(catbeg),'ZBAK=','=','ctrl zbak',sctrl,0,
     .  'Background charge',1,4,0,TF(25))
      nchoos = 0
      i0 = partks(recrd(catbeg),'SAVDOS=','=','ctrl ldos,7',sctrl,0,
     .        'Choose some combination of the following:'
     .//'%N%10f1 Write DOS to directly disk (NPTS and DOS also needed)'
     .//'%N%10f2 Write weights for partial DOS'
     .//'%N%10f4 Same as (2), but weights m-resolved',
     .  1,2,0,TF(12))
      i0 = partks(recrd(catbeg),'DOS=','=','bz dosw',sbz,0,
     .  'Energy window over which DOS accumulated',-2,4,0,TF(13))
      i0 = partks(recrd(catbeg),'NPTS=','=','bz ndos',sbz,0,
     .  'No. DOS points (sampling integration, and lmdos)',1,2,0,TF(11))
      i0 = partks(recrd(catbeg),'PUTQP=','=','bz lio,2',sbz,0,
     .  'Write qp to disk',1,0,0,TF(3))
#if LM
      i0 = partks(recrd(catbeg),'EFMAX=','=','bz efmax',sbz,0,
     .  'Find evecs up to efmax',1,4,0,TF(14))
      i0 = partks(recrd(catbeg),'NEVMX=','=','bz nevmx',sbz,0,
     .  'Find at most nevmx eigenvectors'//
     .  '%N%10fIf NEVMX=0, program uses internal default'//
     .  '%N%10fIf NEVMX<0, no eigenvectors are generated',
     .  1,2,0,TF(15))
      i0 = partks(recrd(catbeg),'ZVAL=','=','bz zval',sbz,0,'Number of'
     .  //' electrons to accumulate in BZ integration',-1,4,0,TF(16))
      i0 = partks(recrd(catbeg),'NOINV=','=','ctrl lqp,1',sctrl,0,
     .  'Suppress automatic inclusion of inversion symmetry for BZ',
     .  1,0,0,TF(17))
      i0 = partks(recrd(catbeg),'FSMOM=','=','bz fsmom',sbz,0,'Fixed'//
     .  '-spin moment (fixed-spin moment method)',-1,4,0,TF(21))
      i0 = partks(recrd(catbeg),'DMAT=','=','bz lio,8',sbz,0,
     .  'Calculate density matrix',1,0,0,TF(24))
C     This needs to be cleaned up to handle more general cases.
C     Wait for rewrite of partok. See mkcond.f for different modes.
C     cond=1: ballistic conductivity.  vec=direction vector
C     cond=2: v_i . v_i.  1 or 2 elts of vec should be 1,marking i,j
      i0 = partks(recrd(catbeg),'COND:',':, ',
     .  'bz lcond;cond:1.1 lcond;vec:2?$',
     .  sbz,0,'Calculate ballistic conductivity',-1,5,0,TF(18))
      i0 = partks(recrd(catbeg),'INVIT=','=','ctrl lqp,2',sctrl,0,
     .  'Use inverse iteration for diagonalization',1,0,0,TF(19))
C      i0 = partks(recrd(catbeg),'BZMP=','=','ctrl ldos,8',sctrl,0,
C     .  'Brillouin zone map: NOS vs k ',-1,0,0,TF(21))
      i0 = partks(recrd(catbeg),'MULL=','=','bz lmull',sbz,0,
     .  'Mulliken population analysis ',-1,2,0,TF(22))
      i0 = partks(recrd(catbeg),'EMESH=','=','bz semsh',sbz,0,
     .  'Mesh for energy contour integration'
     .  //'%N%10fentry 1   2    3    4    5   6    7    8    9    10'
     .  //'%N%16fnz mode ebot etop  ----  depend on mode ---- '
     .  //'(see documentation)',
     .  10,4,0,TF(23))
      if (optio .eq. 0 .and. TF(23) .ne. 2) then
        write(*,753)
  753   format(10x,'modes 0-2',14x,'Im(z)'
     .        /10x,'modes 10,11',13x,'ecc eps'
     .        /10x,'Add  +100,200',11x,'        nzne vne delne delse'
     .        /10x,'Add  +300,400',11x,'        delz'
     .    )
      endif
#endif
      goto 1000

C --- Ewald sums ---
C Token ordering: (0-3) EWALD NKDMX AS TOL
  130 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading EWALD',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0000',toksw,3,TF)
      call getcat(recrd,'EWALD ','for Ewald summations',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'AS=','=','lat as',
     .  slat,0,'Ewald smoothing parameter',-1,4,0,TF(2))
      i0 = partks(recrd(catbeg),'TOL=','=','lat tol',
     .  slat,0,'Ewald tolerance',-1,4,0,TF(3))
      i0 = partks(recrd(catbeg),'NKDMX=','=','lat nkdmx',
     .  slat,0,'Ewald tolerance',-1,2,0,TF(1))
      i0 = partks(recrd(catbeg),'RPAD=','=','lat rpad',
     .  slat,0,'Scale rcutoff by rpad when lattice vectors padded '//
     .  'in oblong geometries',-1,4,0,TF(3))
      call pack1('lat nkqmx',slat,igets('lat nkdmx',slat))
      goto 1000

C --- Mixing ---
C Token ordering: MIX; (1-8)  MODE AMODE XIPMX BETV TOL BETA (UMIX,TOLU,NITU)
  140 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading MIX',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('000000002',toksw,3,TF)
      call getcat(recrd,'MIX ',
     .  'for density mixing in self-consistency cycle',TF)
      if (.not. noerr) goto 1000

      i0 = prtkss(recrd(catbeg),'MODE=','=; ','strn mix',sstrn,
     .  catsiz,'Mixing rules for charge mixing.  Syntax:',0,TF(1)+100)
      if (optio .eq. 0 .and. TF(1) .ne. 2) print 345
  345 format(10x,
     .'A[nmix][,b=beta][,bv=betv][,n=nit][,w=w1,w2][,nam=fn][,k=nkill]',
     .  '[;...]  or'/
     .10x,'B[nmix][,b=beta][,bv=betv][,wc=wc][,n=#][,w=w1,w2][,nam=fn]',
     .  '[,k=nkill]')
C ... Sanity check
      if (optio .gt. 0) then
        call upack('mix b w wc mode mmix',smix,xv(4),xv,xv(5),i,j)
        call upack1('mix bv',smix,xv(6))
        xv(7) = 0
        call spacks(1,'mix fn',smix,alabl,0,0)
        if (i0 .ge. 0) call pshpr(0)
        call upacks('strn mix',i1,ii)
        if (.not. parmxp(-1,sstrn(i1:ii),ii-i1+1,i.ne.0,j,xv,xv(4),
     .    xv(7),alabl,xv(5),-1,xv(6),0d0)) call rx(
     .    'RDCCAT: parse in parmxp failed')
        if (i0 .ge. 0) call poppr
      endif
#ifndef LMMC
      i0 = partks(recrd(catbeg),'CONV=','=','ctrl tol:3.3',sctrl,0,
     .  'tolerance in energy change from prior iteration for'//
     .  ' self-consistency',-1,4,0,TF(5))
      i0 = partks(recrd(catbeg),'CONVC=','=','ctrl tol:1.1',sctrl,0,
     .  'tolerance in output-input charge for self-consistency',
     .  -1,4,0,TF(5))
#else
C      i0 = partks(recrd(catbeg),'BETASP=','=','mix b',smix,0,
C     .  'mixing parameter for spheres',-1,4,0,TF(5))
C      i0 = partks(recrd(catbeg),'BETAI=','=','mix bv',smix,0,
C     .  'mixing parameter for interstitial',-1,4,0,TF(5))
      i0 = partks(recrd(catbeg),'CONV=','=','ctrl tol:1.1',sctrl,0,
     .  'tolerance in energy change from prior iteration for'//
     .  ' self-consistency',-1,4,0,TF(5))
      i0 = partks(recrd(catbeg),'QTOLSP=','=','ctrl tol:2.2',sctrl,0,
     .  'tolerance in sphere charge from prior iteration for'//
     .  ' self-consistency',-1,4,0,TF(5))
      i0 = partks(recrd(catbeg),'QTOLI=','=','ctrl tol:3.3',sctrl,0,
     .  'tolerance in interstitial charge for self-consistency',
     .  -1,4,0,TF(5))
#endif
#if LM
      i0 = prtkss(recrd(catbeg),'AMODE=','=; ','strn amix',sstrn,
     .  catsiz,'Mixing rules for Euler angles.  Syntax as in MODE',0,
     .  TF(1)+100)
      i0 = partks(recrd(catbeg),'XIPMX=','=','mix lxpot,3',smix,0,
     .  'mix potential independently of charge:'//
     .  '%N%10fXIPMX=1: mix vin and v(qmix)'//
     .  '%N%10fXIPMX=2: mix vin and v(qout)',1,2,0,TF(3))
      if (lgors('mix lxpot,3',smix)) call lsets('ctrl lves',sctrl,T,1)
      if (lgors('mix lxpot,1',smix) .and. lgors('mix lxpot,2',smix))
     .  call rxi('rdccat: illegal value, token XIPMX=',
     .  igets('mix lxpot',smix))

      i0 = partks(recrd(catbeg),'UMIX=','=','mix umix',smix,0,
     .  'Mixing parameter for densmat in LDA+U',1,4,0,TF(8))

      i0 = partks(recrd(catbeg),'TOLU=','=','mix tolu',smix,0,
     .  'Tolerance for densmat in LDA+U',1,4,0,TF(8))

      i0 = partks(recrd(catbeg),'NITU=','=','mix nitu',smix,0,
     .  'Max number of LDA+U iterations of densmat',1,2,0,TF(8))


#endif
      goto 1000

C --- Start ---
C Token ordering: START;
C (1-10)  CNTROL BEGMOM NIT CNVG * RDVES V DV P Q
  150 if ( pass2 .and. optio .eq. 2) goto 1000
      continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading START',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('00000000000',toksw,3,TF)
      call getcat(recrd,'START ',' ',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'NIT=','=','ctrl maxit',sctrl,0,
     .  'maximum number of iterations in self-consistency cycle',-1,2,0,
     .  TF(3))
      i0 = partks(recrd(catbeg),'BEGMOM=','=','ctrl lasa,3',sctrl,0,
     .        '0 start from file pot pars'//
     .  '%N%10f1 Create potential and pot pars from P,Q'//
     .  '%N%10f2 Create pot pars from file potential'//
     .  '%N%10f3 Create potential from restart file',
     .  -1,2,0,TF(2))
      i0 = partks(recrd(catbeg),'CNVG=','=','ctrl tol',sctrl,0,' ',-1,4,
     .  0,TF(4))
C      i0 = partks(recrd(catbeg),'FREE=','=','ctrl lasa,8',sctrl,0,' ',
C     .  -1,0,0,TF(5))
      if (.not. pass2) goto 1000
C     The following checks if to read input from control file
      i0 = partok(recrd(catbeg),'CNTROL=','=',sw,
     .  'if F, the remainder of this category is ignored',-1,0,0,TF(1))
      sw = noerr .and. sw
      if (.not. sw) goto 1000
C     And whether to read VES from the control file
      sw2 = .false.
      i0 = partok(recrd(catbeg),'RDVES=','=',sw2,
     .  'Read Ves(RMT) from this category',-1,0,0,TF(6))
      if (sw2) call shftpp(nclasp,nl*nsp,pp,ves,ves,T,F)
      i1 = iend
      outs = ' '
  151 continue
        call pshpr(0)
        call getcls(recrd,slabl,nspec,F,j)
        call poppr
        j = j+1
C   ... If not one of species, possibly a class member
        if (j .eq. 0) then
          i0 = partok(recrd(catbeg),'ATOM=','=',8,alabl,1,1,i1,F)
          if (i0 .le. 0) goto 156
          call s8tor8(alabl,xx)
          do  158  ii = 1, nclasp
            if (xx .eq. dclabl(ii)) then
              j = ii
              call awrit0('%a,'//alabl,outs,len(outs),0)
              goto 159
            endif
  158     continue
          if (iprint() .ge. 30) call awrit0(' rdccat: (warning) class "'
     .      //alabl//'%a" not in class table',' ',80,lgunit(1))
        else
          call awrit0('%a'//slabl(j)//'%a,',outs,len(outs),0)
        endif
  159   if (j .eq. 0) goto 156
        js = ics(j)
        i1 = iend
C iend now first char after ATOM=nam, j=class index, js=spec index
        if (sw2) then
          dum = 0
          i0 = partok(recrd(catbeg),'V=','=',ves(j),' ',1,4,i1,TF(7))
          i0 = partok(recrd(catbeg),'DV=','=',dum,' ',1,4,i1,TF(8))
          ves(j) = ves(j) + dum
        endif
        call upack('spec lmxa z idmod',sspec,js,lmxa,z,idmod,0)
        i0 = partok(recrd(catbeg),'P=','= ,',pnu(1,1,j),' ',
     .    -(lmxa+1)*nsp,4,i1,TF(9))
        if (noerr) then
C         24 Jan 01 correct for lmax misalignment in second spin channel
          if (lmxa+1 .lt. nl .and. nsp .eq. 2) then
            call upack('spec p',sspec,j,tmp,0,0,0)
            call dcopy(lmxa+1,pnu(lmxa+2,1,j),1,tmp(1+n0),1)
            call dcopy(lmxa+1,pnu(1,1,j),     1,tmp,1)
            call dmcpy(tmp,n0,1,pnu(1,1,j),nl,1,nl,2)
          endif
          call config(pnu(1,1,j),-1,z,ivec,n)
          do  154  i = 1, nsp
          do  154  lp1 = 1, lmxa+1
            ii = ivec(lp1)
            if (pnu(lp1,i,j) .lt. 0) ii = ii-1
            if (dabs(pnu(lp1,i,j)) .gt. 1) ii=0
            if (.not. cmdopt('-keepnu',7,0,a))
     .        pnu(lp1,i,j) = dabs(pnu(lp1,i,j)) + ii
C            if (ii .ne. 0 .and. iprint() .gt. 30) then
C              call awrit2(' rdccat: set pnu(%i,%i,'//slabl(j),
C     .          a,len(a),0,lp1,i)
C              call awrit2('%a) = %1;6d  idmod=%i',
C     .          a,-len(a),-i1mach(2),pnu(lp1,i,j),idmod(lp1))
C            endif
  154     continue

C         30 Sep 04 read enu as well as pnu.  Used when IDMOD=2
          call dcopy(nl,pp(1,1,1,j),6,tmp,1)
          call dcopy(nl,pp(1,1,nsp,j),6,tmp(1+nl*(nsp-1)),1)
          i0 = partok(recrd(catbeg),'ENU=','= ,',tmp,' ',
     .      -(lmxa+1)*nsp,4,i1,TF(9))
          if (i0 .gt. 0 .and. sw2) then
            do  i = 1, (lmxa+1)*nsp
              tmp(i) = tmp(i)-ves(j)
            enddo
          endif
          call dcopy(nl,tmp,1,pp(1,1,1,j),6)
          call dcopy(nl,tmp(1+nl*(nsp-1)),1,pp(1,1,nsp,j),6)

          i0 = partok(recrd(catbeg),'Q=','= ,',qnu(1,1,1,j),' ',
     .      -3*(lmxa+1)*nsp,4,i1,TF(10))
          if (i0 .ne. 3*(lmxa+1)*nsp) then
            iend = i1+1
          else
C           24 Jan 01 correct for lmax misalignment in second spin chan
            if (lmxa+1 .lt. nl .and. nsp .eq. 2) then
              call dpzero(wk,6*n0)
              call dcopy(i0/2,qnu(1,lmxa+2,1,j),1,wk(1+3*n0),1)
              call dcopy(i0/2,qnu(1,1,1,j),     1,wk,1)
              call dmcpy(wk,3*n0,1,qnu(1,1,1,j),3*nl,1,3*nl,2)
            endif
            initc(j) = 2*(initc(j)/2) + 1
          endif
          goto 151
        endif

  156 iend = subsiz
      if (subsiz .eq. catsiz) goto 152
      goto 151
  152 continue
      if (sw2) call shftpp(nclasp,nl*nsp,pp,ves,ves,F,T)
C ... Printout
      if (iosw .eq. 3) then
        call awrit0('%23o%0p%N rdccat: read P,Q for %a',outs,80,
     .    -lgunit(1))
        print *, '   PNU, PPARS'
        do  157  j = 1, nclasp
          call r8tos8(dclabl(j),alabl)
          js = ics(j)
          call upack2('spec lmxa',sspec,js,lmxa)
          if (mod(initc(j),2) .eq. 1) then
            print 335, 'class ', alabl, (pnu(lp1,1,j),lp1=1,lmxa+1)
            if (nsp .eq. 2)
     .      print 335,' spin ', '2       ',(pnu(lp1,2,j),lp1=1,lmxa+1)
            print 335,'      ', 'q       ',(qnu(1,lp1,1,j),lp1=1,lmxa+1)
            if (nsp .eq. 2)
     .      print 335,' spin ', '2       ',(qnu(1,lp1,2,j),lp1=1,lmxa+1)
  335       format(4x,a,a,7f8.4)
          else
            print 335, 'class ', alabl // '  missing pnu'
          endif
          if (mod(initc(j)/2,2) .eq. 1) then
            print 336,'  enu       ',      (pp(1,lp1,1,j), lp1=1,lmxa+1)
            if (nsp .eq. 2)
     .      print 335, ' spin ', '2       ',(pp(1,lp1,2,j),lp1=1,lmxa+1)
            print 336,'    c       ',      (pp(2,lp1,1,j), lp1=1,lmxa+1)
            if (nsp .eq. 2)
     .      print 335, ' spin ', '2       ',(pp(2,lp1,2,j),lp1=1,lmxa+1)
            print 336,'srdel       ',      (pp(3,lp1,1,j), lp1=1,lmxa+1)
            if (nsp .eq. 2)
     .      print 335, ' spin ', '2       ',(pp(3,lp1,2,j),lp1=1,lmxa+1)
  336       format(6x,a,7f8.4)
          else
            print 335, 'class ', alabl // '  missing ppars'
          endif
          print *, ' '
  157   continue
      endif
      goto 1000

C --- Master  ---
C Token ordering: MASTER; JOB
  160 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading MASTER',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('01',toksw,3,TF)
      dum = 0d0
      call getsyv('JOB',dum,i0)
      j = nint(dum)
       if (cmdopt('-job=',5,0,a)) then
        i = 5
        if (.not. a2bin(a,j,2,0,' ',i,-1))
     .    call rxs2('RDCCAT: could not parse "',a(1:50),
     .    '%a" in MASTER category')
        call lodsyv('job',1,dble(j),i)
      endif
      if (i0 .ne. 0) TF(0) = 1
      call getcat(recrd,'MASTER ','defines job-dependent variables',TF)
      if (.not. noerr) goto 1000
C ... Get symbolic variables preceding JOB
      if (catsiz .eq. 0) goto 1000
      subsiz = catsiz
      i1 = 7
      call subcat(recrd(catbeg),'JOB ',' ',i1)
      if (optio .eq. 1) call parsyv(recrd(catbeg),subsiz-1,1999,0,i1)
C ... Get symbolic variables specific to JOB
  165 continue
        subsiz = catsiz
        i0 = partok(recrd(catbeg),'JOB ',' ',i,' ',-1,2,iend,TF(1))
        i1 = iend
        call subcat(recrd(catbeg),'JOB ',' ',i1)
        if (optio .eq. 1 .and. j .ne. i .and.
     .      (j .ne. 0 .or. i .ne. 1)) goto 165
        if (optio .eq. 1) call parsyv(recrd(catbeg),subsiz-1,1999,0,i1)
        iend = i1
      goto 1000

C --- Hamiltonian parameters ---
C Token ordering:
C 1..10 KMESH ELIND FORCES XCFUN FRZ FTMESH TOL NMTO,KMTO,EKAP EWALD VMTZ
C 11..20  QASA (PMIN,PMAX) RDSIG RSRNGE,SIGTOL SIGP ALFSI DABC DQVAL EBAS
C 21..26  UDIAG PWMODE PWEMAX PWEMIN PWPAD OVTOL
  170 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading HAM',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('000000000002022222222222222',toksw,3,TF)
      call getcat(recrd,'HAM ','Parameters defining hamiltonian',TF)
      if (.not. noerr) goto 1000

      nchoos = 2
C     Define mesh implicitly through gmax
      i0 = partks(recrd(catbeg),'GMAX=','=','lat gmax',
     .  slat,0,'Energy cutoff for FT mesh',-1,4,0,TF(6))
      xx = dgets('lat gmax',slat)
C     Define mesh explicitly
      if (.not. noerr .or. xx .eq. 0) then
        call pvcat1('FTMESH=','No. divisions for charge density mesh',
     .    'lat nabc',slat,TF(6))
      endif

      i0 = partks(recrd(catbeg),'TOL=','=','lat tolft',
     .  slat,0,'FT mesh tolerance',-1,4,0,TF(7))

      i0 = partks(recrd(catbeg),'FRZWF=','=','ctrl lbas,16',sctrl,0,
     . 'Freeze augmentation wave functions for all species',1,0,0,TF(5))

      i0 = partks(recrd(catbeg),'FORCES=','=','ctrl lfrce',sctrl,0,
     .  'controls the ansatz for density shift in force calculation.'//
     .  '%N%10f-1 no force%10f0 no shift'//
     .  '%N%10f 1 free-atom shift  12 screened core+nucleus',
     .  -1,2,0,TF(3))

      i0 = partks(recrd(catbeg),'ELIND=','=','ham elind',sham,0,
     .  'Lindhard energy for model screening',1,4,0,TF(2))

      i = mod(igets('ctrl lxcf',sctrl),10)
      j = igets('ctrl lxcf',sctrl) - i
      i0 = partok(recrd(catbeg),'XCFUN=','= ,',i,
     .  'Chooses local exchange correlation functional:'//
     .  '%N%10f 1 for Ceperly-Alder'//
     .  '%N%10f 2 for Barth-Hedin (ASW fit)'//
     .  '%N%10f 3 for PW91 %N%10f 4 for PBE',
     .  -1,2,0,TF(4))
C     Pack LDA switch into 1's digit of lxcf
      if (i0 .eq. 1) then
        call pack1('ctrl lxcf',sctrl,i+j)
      endif
      i = igets('ctrl lxcf',sctrl)/10
      i0 = partok(recrd(catbeg),'GGA=','= ,',i,
     .  'Chooses GGA functional:'//
     .  '%N%10f 0 for LSDA, 1 for LMH, 2 for PW91,'//
     .  ' 3 for PBE, 4 for PBE with Becke exchange',-1,2,0,TF(4))
C     Pack GGA switch into 10's digit of lxcf
      if (i0 .eq. 1) then
        call pack1('ctrl lxcf',sctrl,
     .              10*i+mod(igets('ctrl lxcf',sctrl),10))
      endif
      i0 = partks(recrd(catbeg),'EBAS=','=','ham kmto',
     .  sham,0,'lmto basis kinetic energies (molecules)',3,4,0,TF(20))
      if (optio .ne. 0 .and. i0 .gt. 0) then
        call pack1('ham nmto',sham,i0)
      endif
      i0 = partks(recrd(catbeg),'ALFSI=','=','ham alfsi',sham,0,
     .  'coefficient of artificial term (0 turns off)',1,4,0,TF(17))
      i0 = partks(recrd(catbeg),'DABC=','=','ham dabc',sham,0,
     .  'Spacing for real-space mesh (molecules)',3,4,0,TF(18))
      if (i0 .lt. 3) then
        call upack1('ham dabc',sham,xv)
        if (i0 .eq. 1) xv(2) = xv(1)
        if (i0 .eq. 1 .or. i0 .eq. 2) xv(3) = xv(2)
        call pack1('ham dabc',sham,xv)
      endif
      i0 = partks(recrd(catbeg),'DQVAL=','=','ctrl zbak',sctrl,0,
     .  'Total charge (molecule)',1,4,0,TF(19))

      i0 = partks(recrd(catbeg),'RDSIG=','=','ham lsig',sham,0,
     .  'Add self-energy to local exchange correlation functional:'//
     .  '%N%10f 1s digit:'//
     .  '%N%10f   0 do not read Sigma'//
     .  '%N%10f   1 to add sigma to ham.'//
     .  '%N%10f   Add 2 to symmetrize Sigma(T)'//
     .  '%N%10f   Add 4 to include Re(Sigma(T)) only',
     .  1,2,0,TF(13))
      if (optio .eq. 0 .and. tf(13) .ne. 2) call info0(0,0,0,
     .    '%10f 10s digit:'//
     .  '%N%10f   0 simple interpolation'//
     .  '%N%10f   1 approx high and low sigma by diagonal (see sigp)'//
     .  '%N%10f   3 interpolate sigma by LDA evecs'//
     .  '%N%10f     In this case use 100s digit for:'//
     .  ' no. interpolation points'//
     .  '%N%10f Add 10000 to flag sigma stored without symops'//
     .  '%N%10f Add 20000 to use minimum neighbor table'//
     .  '%N%10f Add 40000 to allow file qp mismatch')

      i0 = partks(recrd(catbeg),'RSRNGE=','=','ham rsrnge',sham,0,
     .  'Maximum range in connecting vectors for r.s. sigma '//
     .  '(units of alat)',1,4,0,TF(14))

      i0 = partks(recrd(catbeg),'RSSTOL=','=','ham rsstol',sham,0,
     .  'Max tolerance in Bloch sum error for r.s. sigma ',1,4,0,TF(14))

      i0 = partks(recrd(catbeg),'SIGP:',':, ',
     .  'ham sigp;mode:1.1 sigp;nmin:2.2 sigp;emin:3.3 sigp;nmax:4.4 '//
     .  'sigp;emax:5.5 sigp;a:6.6 sigp;b:7.7 sigp;efit:8.8',
     .  sham,0,'Parameters for approximating sigma:'//
     .  '%N%10fmode: specifies how sigii is to be constrained:',
     .  -1,5,0,TF(15))
      if (optio .eq. 0 .and. tf(15) .ne. 2) then
        call info0(0,0,0,
     .    '%15f 0 sigii > a + b*e'//
     .  '%N%15f 1 sigii = a + b*e'//
     .  '%N%15f 2 sigii > a  and  sigii < b'//
     .  '%N%15f 3 sigii = a + b*e'//
     .  '%N%10fnmin: approx levels 1..nmin'//
     .  '%N%10femin: approx levels e<emin (not used if nmin>=0)'//
     .  '%N%10fnmax: approx levels nmax+1..'//
     .  '%N%10femax: approx levels e>emax (not used if nmax>0)'//
     .  '%N%10fa and b: '//
     .    'constraint parameters whose meaning depends on mode'//
     .  '%N%10fefit (mode 3): fit sigii for E(lda)>efit')
      endif

      i0 = partks(recrd(catbeg),'NMTO=','=','ham nmto',sham,0,
     .  'order of polynomial approximation for NMTO'//
     .  '%N%10f(For 3rd generation, require NMTO>=2)',-1,2,0,TF(8))
      if (optio .eq. 0) then
        nmto = 99
      else
        nmto = igets('ham nmto',sham)
        if (i0.gt.0 .and. nmto.gt.0) call pack1('ctrl lgen3',sctrl,1)
        if (nmto .ne. 0) call isanrg(nmto,2,5,'rdccat','nmto',.true.)
      endif
      if (nmto .ne. 0) i0 = partks(recrd(catbeg),'KMTO=','=','ham kmto',
     .  sham,0,'corresponding NMTO kinetic energies',-nmto,4,0,TF(8))

      if (nmto .ne. 0) i0 = partks(recrd(catbeg),'EWALD=','=',
     .  'ctrl lgen3,2',sctrl,0,'make strux by Ewald summation',-1,0,0,
     .  TF(9))


      i0 = partks(recrd(catbeg),'UDIAG=','=','ham udiag',sham,0,
     .  'nonzero => diagonal-only LDA+U',1,2,0,TF(21))

      i0 = partks(recrd(catbeg),'VMTZ=','=','pot vmtz0',spot,0,
     .  'Muffin-tin zero defining wave functions',1,4,0,TF(10))

      i = igets('ctrl lasa,128+256',sctrl)/128
      j = igets('ctrl lasa',sctrl) - 128*i
      i0 = partok(recrd(catbeg),'QASA=','= ,',i,
     .  '0 => Methfessel conventions for 2nd gen ASA moments '//
     .             'Q0..2'//
     .  '%N%10f1 => Q2 = coff to phidot**2 - p phi**2 in sphere'//
     .  '%N%10f2 => Q1,Q2 accumulated as coffs to <phi*phidot> '//
     .                 'and <phidot**2>'//
     .  '%N%10f3 => 1+2 (Stuttgart conventions)',-1,2,0,TF(11))
      call pack1('ctrl lasa',sctrl,j+128*mod(i,4))


      i0 = partks(recrd(catbeg),'PMIN=','= ,','ham pmin',sham,0,
     .  'global minimum in fractional part of potential '//
     .  'functions.%N%10fOne # for each l does the following:'//
     .  '%N%10f0: no minimum constraint'//
     .  '%N%10f#: with #<1, floor of fractional P is #'//
     .  '%N%10f1: use free-electron value as minimum',
     .  10,2,0,TF(12))
      i0 = partks(recrd(catbeg),'PMAX=','= ,','ham pmax',sham,0,
     .  'global maximum in fractional part of potential '//
     .  'functions.%N%10fOne # for each l does the following:'//
     .  '%N%10f0: no maximum constraint'//
     .  '%N%10f#: with #<1, ceiling of fractional P is #',
     .  10,2,0,TF(12))

      i0 = partks(recrd(catbeg),'OVEPS=','= ,','ham oveps',sham,0,
     .  'Diagonalize hamiltonian in reduced hilbert space,'//
     .  '%N%10fdiscarding part with evals of overlap < OVEPS',
     .  1,4,0,TF(26))

      i0 = partks(recrd(catbeg),'PWMODE=','= ,','ham pwmode',sham,0,
     .  'Controls APW addition to LMTO basis'//
     .  '%N%10f1s digit:'//
     .  '%N%10f0: LMTO basis only'//
     .  '%N%10f1: Mixed LMTO+PW'//
     .  '%N%10f2: PW basis only'//
     .  '%N%10f10s digit:'//
     .  '%N%10f0: PW basis fixed'//
     .  '%N%10f1: PW basis q-dependent',
     .  1,2,0,TF(22))
      if (igets('ham pwmode',sham) .ne. 0) TF(23) = 1
      i0 = partks(recrd(catbeg),'PWEMAX=','= ,','ham pwemax',sham,0,
     .  'Energy cutoff in PW part of basis (Ry)',1,4,0,TF(23))
      i0 = partks(recrd(catbeg),'PWEMIN=','= ,','ham pwemin',sham,0,
     .  'Minimum-Energy bound in PW part of basis (Ry)',1,4,0,TF(24))
      i0 = partks(recrd(catbeg),'NPWPAD=','= ,','ham npwpad',sham,0,
     .  'Overrides default padding of variable basis dimension',1,2,0,
     .  TF(25))
      if (igets('ham pwmode',sham) .ne. 0) TF(23) = 1

      goto 1000

C --- Plane  ---
C Token ordering: PLANE; NORMAL X Y
  180 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading PLANE',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0000',toksw,3,TF)
      call upack1('lat plat2',slat,xv)
      call getcat(recrd,'PLANE ','Parameters defining a plane',TF)
      if (.not. noerr) goto 1000
      i0 = partok(recrd(catbeg),'NORMAL=','=',xv(7),
     .  ' ',-3,4,0,TF(1))
      i0 = partok(recrd(catbeg),'X=','=',xv,' ',-3,4,0,TF(2))
      i0 = partok(recrd(catbeg),'Y=','=',xv(4),' ',-3,4,0,TF(3))
      call pack1('lat plat2',slat,xv)
      goto 1000

C --- TB ---
C Note at 6.10b switch 2^13 (NOMAD) was freed up
C Token ordering: TB;
C (1-11)  OVLP CRYSF OVCF ADDES RMAXH RMAXR FORCES FIJ MD 3PV EVDISC
C (12-23) PAIR TRH RHO U NOMAD UAVG CGMIN CGTOL * * * *
C (24-31) START (TSTEP,TS/FS) NITEQU (TEMP,T/K) SYRUP NIT DONLY
C (31-36) FMODE NBFIT EBFIT RFIT LAM SCL
  190 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading TB',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
C    (a ruler)  01234567890123456789012345678901234567890
      i = ipsw('10000000000000000000000000000000000002222',toksw,3,TF)
      call getcat(recrd,'TB ',' ',TF)
      if (.not. noerr) goto 1000
#if LM
      i0 = partks(recrd(catbeg),'OVLP=','=','ctrl ltb,1',sctrl,0,
     .  'Non orthogonal tight-binding ',-1,0,0,TF(1))
      i0 = partks(recrd(catbeg),'CRYSF=','=','ctrl ltb,2',sctrl,0,
     .  'Crystal field terms in hamiltonian ',-1,0,0,TF(2))
      i0 = partks(recrd(catbeg),'OVCF=','=','ctrl ltb,4',sctrl,0,
     .  'Crystal field terms in overlap ',-1,0,0,TF(3))
      i0 = partks(recrd(catbeg),'ADDES=','=','ctrl ltb,8',sctrl,0,
     .  'add ebarLL'' * sLL'' to hLL''',-1,0,0,TF(4))
      i0 = partks(recrd(catbeg),'RMAXH=','=','str rmax',sstr,0,
     .  'Hamiltonian cut-off length in units of a',-1,4,0,TF(5))
C      if (optio .eq. 1) rmaxr = rmaxh
C      i0 = partok(recrd(catbeg),'RMAXR=','=',rmaxr,
C     .  'Density matrix cut-off',-1,4,0,TF(6))
      i0 = partks(recrd(catbeg),'FORCES=','=','ctrl ltb,16',sctrl,0,
     .  'calculate forces ',-1,0,0,TF(7))
      i0 = partks(recrd(catbeg),'FIJ=','=','ctrl ltb,32',sctrl,0,
     .  'to get forces by atom pair, e.g., for stresses ',-1,0,0,TF(8))
      i0 = partks(recrd(catbeg),'3PV=','=','ctrl ltb,128',sctrl,0,
     .  'calculate pressure',-1,0,0,TF(10))
      i0 = partks(recrd(catbeg),'EVDISC=','=','ctrl ltb,256',sctrl,0,
     .  'can be F '//
     .  'for insulators or to save space for metals ',-1,0,0,TF(11))
      i0 = partks(recrd(catbeg),'PAIR=','=','ctrl ltb,512',sctrl,0,
     .  'pair potential only',-1,0,0,TF(12))
      i0 = partks(recrd(catbeg),'TRH=','=','ctrl ltb,1024',sctrl,0,
     .  'calculate local projection of band energies',-1,0,0,TF(13))
      i0 = partks(recrd(catbeg),'RHO=','=','ctrl ltb,2048',sctrl,0,
     .  'calculate local projection of charges',-1,0,0,TF(14))
      i0 = partks(recrd(catbeg),'U1=','=','ctrl ltb,4096',sctrl,0,
     .  'for electrostatics with L>=0 (first model)',-1,0,0,TF(15))
      i0 = partks(recrd(catbeg),'UL=','=','ctrl ltb,2^15',sctrl,0,
     .  'for electrostatics with L>=0',-1,0,0,TF(15))
      i0 = partks(recrd(catbeg),'TBU=','=','ctrl ltb,2^13',sctrl,0,
     .  'tight-binding+U',-1,0,0,TF(15))
      if (optio .eq. 0)
     . print*, '  ... if U1=T, UL=T or TBU=T,'//
     .         ' the next two tokens are used:'
      i0 = partks(recrd(catbeg),'NOUAVG=','=','ctrl ltb,16384',sctrl,0,
     .  'use individual U_l from Q=; don''t average them',-1,0,0,TF(17))
      i0 = partks(recrd(catbeg),'IODEL=','=','ctrl ltb,2^16',sctrl,0,
     .  'Attempt to read increments from disk',-1,0,0,TF(15))
      if (lgors('ctrl ldos,8',sctrl)) call lsets('ctrl ltb',sctrl,T,64)
      i0 = partks(recrd(catbeg),'GAMMA=','=','ctrl ltb,2^17',sctrl,0,
     .  'do gamma-point only',-1,0,0,TF(15))
      i0 = partks(recrd(catbeg),'MOL=','=','ctrl ltb,2^18',sctrl,0,
     .  'molecule: no PBCs',-1,0,0,TF(15))
      if (lgors('ctrl ltb,2^17',sctrl) .or.
     .    lgors('ctrl ltb,2^18',sctrl)) then
        call pack1('bz nkabc',sbz,dvec1)
        call pack1('bz lshft',sbz,dvec2)
      endif

      if (optio .eq. 0) print *,
     .  '  ... the next six tokens are used for tight-binding fitting'
      i0 = partks(recrd(catbeg),'FMODE=','=','tb fmode',stb,0,
     .  'Fit mode: 0 fit params, 1 fit in range of RFIT',1,2,0,TF(31))
      i0 = partks(recrd(catbeg),'NBFIT=','=','tb nbfit',stb,0,
     .  'Min and max band number to fit',-2,2,0,TF(32))
      i0 = partks(recrd(catbeg),'EBFIT=','=','tb ebfit',stb,0,
     .  'Energy range to fit eigenvalues',-2,4,0,TF(33))
      i = igets('tb fmode',stb)
      if (i .eq. 1 .or. optio .eq. 0) TF(34) = 1
      i0 = partks(recrd(catbeg),'RFIT=','=','tb rmfit',stb,0,
     .  '(fmode=1) range for fitting Slater-Koster parameters',
     .  -2,4,0,TF(34))
      i0 = partks(recrd(catbeg),'LAM=','=','tb alam',stb,0,
     .  'Initial value of Levenberg-Marquardt lambda variable',
     .  -1,4,0,TF(35))
      i0 = partks(recrd(catbeg),'SCL=','=','tb alsc',stb,0,
     .  'Scale factor for Levenberg-Marquardt lambda variable',
     .  -1,4,0,TF(36))
 

C      i0 = partok(recrd(catbeg),'CGTOL=','=',cgtol,
C     .  'Tolerance for O(N) conjugate gradients',-1,4,0,TF(19))

C ... Sanity checks
      call rxx(igets('ctrl ltb,1+4',sctrl).eq.4,
     .  'RDCCAT: OVCF=T but OVLP=F')
      call rxx(igets('ctrl ltb,1+8',sctrl).eq.8,
     .  'RDCCAT: ADDES=T but OVLP=F')
C ... Turn on forces for whatever switches require them
      if (dgets('ctrl mdprm',sctrl) .gt. 0 .or.
     .    lgors('ctrl ltb,32',sctrl)) call lsets('ctrl ltb',sctrl,T,16)
#endif

C      i0 = partok(recrd(catbeg),'EAM=','=',switch(63),
C     .  'Include EAM-like contribution to total energy',-1,0,0,TF(31))
      goto 1000

C --- Map ---
C Token ordering: MAP; MAP
  200 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading MAP',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('01',toksw,3,TF)
      call getcat(recrd,'MAP ',' ',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'MAP ',' ','ctrl nmap,1',sctrl,0,' ',1,
     .  0,0,TF)
      if (i0 .eq. 1) then
        if (lgors('ctrl nmap,1',sctrl).or.cmdopt('-map',4,0,outs)) then
          call lstra('strn map',i,offs,j)
          call ustrn(tmp,-offs,1,i,casti,catsiz)
          call strncp(sstrn(i:),recrd(catbeg),1,1,catsiz)
          sstrn(i:i+iend-1) = ' '
          call pack1('ctrl nmap',sctrl,catsiz/reclen)
        endif
      endif
      goto 1000

C --- Planar Green's function ---
C Token ordering: (PGF,MODE) (PLATL,PLATR) GFOPTS
  210 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading PGF',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('010',toksw,3,TF)
      call getcat(recrd,'PGF ','   planar Green''s function',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'MODE=','=','ctrl lpgf',sctrl,0,
     .      '0: do nothing'//
     .'%N%10f1: diagonal layer GF'//
     .'%N%10f2: left- and right-bulk bulk GF'//
     .'%N%10f3: find k(E) for left bulk'//
     .'%N%10f4: find k(E) for right bulk'//
     .'%N%10f5: Calculate current',
     .  1,2,0,TF(0))
      if (optio.ne.0 .and. .not.lgors('ctrl lpgf,-1',sctrl)) goto 1000
      i0 = partks(recrd(catbeg),'SPARSE=','=','ctrl lpgf:2',sctrl,0,
     .  '1: Calculate GF with LU decomposition',
     .  1,2,0,TF(0))
      i0 = partks(recrd(catbeg),'PLATL=','=','lat platl:7',slat,0,
     .  'third lattice vector of left bulk',3,4,0,TF(1))
      i0 = partks(recrd(catbeg),'PLATR=','=','lat platr:7',slat,0,
     .  'third lattice vector of right bulk',3,4,0,TF(1))
      i0 = prtkss(recrd(catbeg),'GFOPTS=','=; ','strn gfopt',sstrn,
     .  catsiz,'Switches governing execution of PGF',0,TF(2)+100)
      goto 1000

C --- TBFIT ---
C Token ordering: TBFIT; FMODE NBFIT EBFIT RFIT
  220 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading TBFIT',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('00000',toksw,3,TF)
      call getcat(recrd,'TBFIT ','For fitting of TB parameters to'//
     .  ' energy bands',TF)
      if (.not. noerr) goto 1000
      call rx('tbfit input not implemented')
C      i0 = partok(recrd(catbeg),'FMODE=','=',fmode,
C     .  'Fit mode: 0 fit params, 1 fit in range of RFIT',1,2,0,TF(1))
C      i0 = partok(recrd(catbeg),'NBFIT=','=',nbfit,
C     .  'Min and max band number to fit',-2,2,0,TF(2))
C      i0 = partok(recrd(catbeg),'EBFIT=','=',ebfit,
C     .  'Energy range to fit eigenvalues',-2,4,0,TF(3))
CC     if (fmode .eq. 1 .or. optio .eq. 0) TF(4) = 1
C      i0 = partok(recrd(catbeg),'RFIT=','=',rmfit,'(fmode=1) range '//
C     .  'for fitting Slater-Koster parameters',-2,4,0,TF(4))
      goto 1000

C --- Dynamics ---
C Token ordering: DYN;
C (1-6)  HAM KT TS NTS NITEQU START
C (7-14) TOL TS0 MX MI NSEQ MODET CT NTHERM ---- for micromagnetics
C (15-)  MSTAT MD NIT
  230 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading DYN',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('000000000000000000',toksw,3,TF)
      call getcat(recrd,'DYN  ','parameters for dynamics',TF)
      if (.not. noerr) goto 1000

      i0 = partks(recrd(catbeg),'MSTAT:',':, ',
     .  'ctrl mdprm;relax:1.1 mdprm;hess:2.2 mdprm;xtol:3.3 '//
     .  'mdprm;gtol:4.4 mdprm;step:5.5 mdprm;nkill:6.6',
     .  sctrl,0,'Parameters for molecular statics:'//
     .  '%N%10frelax=0 (no relaxation) 4 (conjugate gradients) 5 '//
     .  '(Fletcher-Powell)%N%16f6 (Broyden) '//
     .  '%N%10fhess=T  read hessian matrix'//
     .  '%N%10fxtol=#   convergence criterion in displacements'//
     .  '%N%19f>0: use length;  <0: use max val; =0: do not use'//
     .  '%N%10fgtol=#   convergence criterion in gradients'//
     .  '%N%19f>0: use length;  <0: use max val; =0: do not use'//
     .  '%N%10fstep=#  step length'//
     .  '%N%10fnkill=# Remove hessian after # iter',
     .  -1,5,0,TF(15))

      if (dgets('ctrl mdprm',sctrl) .gt. 100 .or. optio .eq. 0) then
        if (optio .eq. 0)
     .    print*, '  ... if MSTAT with lattice relaxation,'//
     .    ' the next token is used:'
        i0 = partks(recrd(catbeg),'PDEF=','= ,','ctrl defm',sctrl,0,
     .    'lattice deformation modes',-6,4,0,TF(15))
      endif

      if (dgets('ctrl mdprm',sctrl) .lt. 4 .or. optio .eq. 0) then
        i0 = partks(recrd(catbeg),'MD:',':, ',
     .  'ctrl mdprm;dyn:1.1 mdprm;tstep:2.2 mdprm;temp:3.3 '//
     .  'mdprm;taup:4.4 mdprm;time:5.5 mdprm;taub:6.6',
     .  sctrl,0,'Parameters for molecular dynamics:'//
     .  '%N%10fdyn=0 (no MD) 1 (NVE) 2 (NVT) 3 (NPT) '//
     .  '%N%10ftstep=#  time step (a.u.)'//
     .  '%N%10ftemp=#   temperature (a.u.)'//
     .  '%N%10ftaup=# thermostat relaxation time (a.u.)'//
     .  '%N%10ftime=# total MD time (a.u.)'//
     .  '%N%10ftaub=# barostat relaxation time (a.u.)',-1,5,0,TF(15))
      if (optio .eq. 0 .and. tf(15) .ne. 2) print
     .  '(10x,''NB: 1 deg.K = 6.3333e-6 a.u.; 1 fs = 20.67098 a.u.'')'

        if (i0 .eq. 1 .and. dgets('ctrl mdprm',sctrl) .gt. 3)
     .    call rx('RDCCAT: illegal "dyn" parameter in token "MD"')
      endif
      if (optio .eq. 0 .and. tf(16) .ne. 2) print
     .  '(10x,''NB: 1 deg.K = 6.333e-6 a.u.; 1 sec = 20.671e15 a.u.'')'
      i0 = partks(recrd(catbeg),'NIT=','=','ctrl nitmv',sctrl,0,
     . 'maximum number of relaxation steps (MSTAT or MD)',-1,2,0,TF(17))

C ... Make sure 'HAM' is there; give up if not
      if (optio .ne. 2)
     .  i0 = partok(recrd(catbeg),'HAM ',' ',1,outs,1,1,0,TF(1))
C ... Find the beginning and the end, and poke into sstrn
      if (optio .ne. 2 .and. i0 .gt. 0) then
C       Get offs(1) = index to mmham
        call lstra('strn mmham',tmp,offs,i)
        i1 = 0
        call subcat(recrd(catbeg),'HAM ',' ',i1)
        i1 = subsiz-1
        i = i1 + 5
        call subcat(recrd(catbeg),'ENDH ',' ',i1)
C       Allocate string of length 1 to get offset
        call ustrn(tmp,-offs,1,offi,i1,1)
C       copy to sstrn, compressing spaces; return length in j
        call cpstr(recrd(catbeg),subsiz,101,'}',i,j,sstrn(offi(1):))
C       Reallocate string of proper length
        call ustrn(tmp,-offs,1,offi,i1,j)
      endif
      subsiz = catsiz
C      this part isn't implemented yet
C      call u_sdyn(s_sdyn,tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))
C      if (optio .eq. 1) it(2) = 0
C      call u_bsi(s_sdyn,i1,ntherm,ivec,tmp2,tmp2(4),tmp2(5),j,it,it(2))
      i0 = partks(recrd(catbeg),'KT=','=','move kt',smove,0,
     .  'Temperature, in a.u.',1,4,0,TF(2))
      i0 = partks(recrd(catbeg),'NTHERM=','=','move nmodt',smove,0,
     .  'Number of thermostats',-1,2,0,TF(14))
C     Number of thermostats is number of modes, coffs to read
      j = -igets('move nmodt',smove)
      if (optio .eq. 0) j = 3
      i0 = partks(recrd(catbeg),'MODET=','=','move modt',smove,0,
     .  'thermostat mode(s)',j,2,i1,TF(12))
      i0 = partks(recrd(catbeg),'CT=','=','move ct',smove,0,
     .  'corresponding thermostat coefficient(s)',j,4,i1,TF(13))
      i0 = partks(recrd(catbeg),'TS=','=','move ts',smove,0,
     .  'time step in a.u.',-1,4,0,TF(3))
      i0 = partks(recrd(catbeg),'TTOT=','=','move tstot',smove,0,
     .  'duration of total simulation, in a.u.',-1,4,0,TF(4))
      i0 = partks(recrd(catbeg),'TEQU=','=','move tsequ',smove,0,
     .  'equilibration time, a.u.',-1,4,0,TF(5))
      if (optio .eq. 0 .and. TF(7) .ne. 2) print *,
     .  '  ... the following control the Bulirsch-Stoer integration:'
      i0 = partks(recrd(catbeg),'NEW=','=','move prmint:2.2,1',smove,0,
     .  'start new SD run',-1,0,0,TF(6))
      i0 = partks(recrd(catbeg),'BSTOL=','=','move prmint:4',smove,0,
     .  'tolerance in numerical integration',-1,4,0,TF(7))
      i0 = partks(recrd(catbeg),'TS0=','=','move prmint:3',smove,0,
     .  'minimum time step in units of TS',-1,4,0,TF(8))
      i0 = partks(recrd(catbeg),'MX=','=','move prmint:5',smove,0,
     .  'maximum order of rational function extrapolation',-1,2,0,TF(9))
      i0 = partks(recrd(catbeg),'MI=','=','move prmint:6',smove,0,
     .  'number of midpoint rules to invoke',-1,2,0,TF(10))
      j = -11
      if (optio .ne. 0) j = igets('move prmint:6',smove)
      if (j .gt. 11) call rxi('rdccat (dyn): mi > 11, mi=',j)
      i0 = partks(recrd(catbeg),'NSEQ=','=','move prmint:7.17',smove,
     .  0,'sequence of no. midpoint divisions',-j,2,0,TF(11))
      goto 1000

C --- Order-N Embedded cluster ---
C Token ordering: ORDERN; (1-5) MODE PLATE RG EMB EFREE
  240 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading ORDERN',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('011100',toksw,3,TF)
      call getcat(recrd,'ORDERN ',' ',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'MODE=','=','ctrl lordn',sctrl,0,
     .  '0: do nothing;  1: embedded cluster',
     .  1,2,0,TF(0))
      if (igets('ctrl lordn',sctrl).eq.0 .and. optio .ne. 0) goto 1000

      i0 = partks(recrd(catbeg),'PLATE=','=','lat plate',slat,0,
     .  'subcluster primitive lattice vectors, in units of alat',
     .  -9,4,0,TF(1))

      i0 = prtkss(recrd(catbeg),'EMB=','=; ','strn gemb',sstrn,
     .  catsiz,'File names for embedding gf',0,TF(4)+100)

      i0 = partks(recrd(catbeg),'EFREE=','=','ordn efre',sordn,0,
     .  'Lowest state for which the free electron propagate'//
     .  '%N(approximately the MT zero)',
     .  -1,4,0,TF(5))

      nchoos = 2
      i0 = partks(recrd(catbeg),'RG=','=','ordn rmaxg',sordn,0,
     .  'Ranges for Green''s function:'//
C     .  '%N%10f(1): range for direct inversion'//
C     .  '%N%10f(2): range for iterative inversion'//
C     .  '%N%10f(3): range for assumed form for G'
     .  ' ',
     .  1,4,0,TF(3))
      if (.not. noerr .and. optio .ne. 2) then
        i0 = partks(recrd(catbeg),'RG/W=','=','ordn rmaxg',sordn,0,
     .    'same as RG, but relative to a.v.w.',1,4,0,TF(3))
        if (i0 .eq. 1) then
          call spack(2,'ordn rmaxg',sordn,offi,casti,nelti,0,0)
          call dscal(1,dgets('lat avw',slat),sordn(offi(1)),1)
        endif
      endif
C      call upack2('lat plat plate',slat,plat,plate)
C      call dqinvb('t',plate,3,0,3,3,wk,3,wk,plat,3,i)
C      print 350, ((plat(i,j),i=1,3),(plate(i,j),i=1,3),j=1,3)
C  350 format(3f10.5,5x,3f10.5)

      call rxx(i.eq.0,'rdccat: plat is singular')
      goto 1000

C --- CMD ---
C Token ordering: CMD
  250 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading CMD',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('0',toksw,3,TF)
      call getcat(recrd,'CMD ',
     .  'Contents are appended to command-line arguments',TF)
      if (.not. noerr) goto 1000
      if (optio .ge. 1) call acmdop(recrd(catbeg+6),catsiz-6,optio-1)
      goto 1000

C --- Green's function ---
C Token ordering: GF MODE GFOPT
  260 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading GF',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('010',toksw,3,TF)
      call getcat(recrd,'GF ','   Green''s function',TF)
      if (.not. noerr) goto 1000
      i0 = partks(recrd(catbeg),'MODE=','=','ctrl lcgf',sctrl,0,
     .      ' 0: do nothing;'//
     .'%N%10f10: generate J(q)'//
     .'%N%10f11: generate derivative properties of J(q)'//
     .'%N%10f12: generate g(R,R'')',
     .  1,2,0,TF(1))
      i0 = prtkss(recrd(catbeg),'GFOPTS=','=; ','strn gfopt',sstrn,
     .  catsiz,'Switches governing execution of GF',0,TF(2)+100)
      goto 1000

C --- Optics ---
C Token ordering:
C OPTIC MODE NE,WINDOW OCRNG UNRNG PARTIAL AXES CLS ESCISS
  270 continue
#if MPI | MPIK
      if (mlog) then
        call gettime(datim)
        call awrit2(' rdccat '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' reading OPTICS',' ',256,lgunit(3),procid,numprocs)
      endif
#endif
      i = ipsw('011000000',toksw,3,TF)
      call getcat(recrd,'OPTICS ','for optics calculations',TF)
      if (.not. noerr) goto 1000
      if (tf(1) .ne. 2) then
        i0 = partks(recrd(catbeg),'MODE=','=','ctrl loptc',sctrl,0,
     .        ' 0: do nothing;'//
     .   '%N%10f-1: generate joint DOS'//
     .   '%N%11f1: generate linear eps_2'//
     .   '%N%11f2: generate linear eps_2 with on-the-fly sampling',
     .    1,2,0,TF(1))
C     .   '%N%11f2: generate second harmonic eps'//
C     .   '%N%11f   1 and 2 allowed in combination',
C     .    1,2,0,TF(1))
        if (igets('ctrl loptc',sctrl).eq.0 .and. optio .ne. 0) goto 1000
        i0 = partks(recrd(catbeg),'NPTS=','=','optic ne',soptic,0,
     .    'No. energy points',1,2,0,TF(2))
        i0 = partks(recrd(catbeg),'WINDOW=','=','optic window',soptic,0,
     .    'Energy window over which to calc. eps',2,4,0,TF(2))
        i0 = partks(recrd(catbeg),'FILBND=','=','optic ocrng',soptic,0,
     .    'Occ. energy bands from which to calc. eps',2,2,0,TF(3))
        i0 = partks(recrd(catbeg),'EMPBND=','=','optic unrng',soptic,0,
     .    'Unocc energy bands from which to calc. eps',2,2,0,TF(4))
        i0 = partks(recrd(catbeg),'PART=','=','optic mode,1',soptic,0,
     .    'Band-to-band decomposition of epsilon',-1,0,0,TF(5))
C       NLO axes.  Should i1 below really be zero?
        i0 = partks(recrd(catbeg),'CHI2:',':, ','optic nchi2 axes',
     .    soptic,0,'For second harmonic generation:'//
     .    '%N%12fnchi2=number of (abc) sets for which to calc. chi2'//
     .    '%N%12faxes=a,b,c for each of the nchi2 sets',
     .    -1,5,0,TF(6))
        i0 = partks(recrd(catbeg),'ESCISS=','=','optic esciss',soptic,0,
     .    'Scissors operator (energy added to unoccupied bands)',
     .      1,4,0,TF(8))

        endif

      goto 1000

 1000 continue

      end

      subroutine pvcat1(token,descr,sname,struc,TF)
C- Reads a token of length 3, padding 2nd and third with first
C     implicit none
      character*(*) token,descr,sname,struc
      integer ivec(3),TF,partok,i0

C ... Input is stored in recrd(offset:)
      character*1 recrd(0:1)
      common /w/ recrd
C ... Common blocks
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr

      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      call upack1(sname,struc,ivec)
      i0 = partok(recrd(catbeg),token,'= ,',ivec,descr,3,2,0,TF)
      if (i0 .eq. 1) ivec(2) = ivec(1)
      if (i0 .eq. 1) ivec(3) = ivec(1)
C     if (i0 .eq. 1 .or. i0 .eq. 2) ivec(3) = ivec(1)
      call pack1(sname,struc,ivec)
      end

      integer function partkf(categ,token,term,sname,struc,
     .  is,chr,count,cast,istart,iopt,fac)
C- See partks
C  Should supersede partks; need add fac
Ci fac  scaling factor when nchoos>0
Ci cast = 100 + real cast: expand to count
C     implicit none
C ... Passed parameters
      integer cast,is,count,istart,iopt
      character*(*) sname,categ(0:*)*1,term,chr,token
      double precision struc(1)
      double precision fac
C ... Local parameters
      integer partks,i0,lcast,i
      double precision x1(100),x2(100),x3(100),x4(100),x5(100)
C iolib common block
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr
      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      lcast = mod(cast,100)

      if (optio .eq. 2 .and. fac .ne. 1 .and. fac .ne. 0) then
        if (optio .eq. 2) i0 = iabs(count)
        if (is .ne. 0) then
          call upack(sname,struc,is,x1,x2,x3,x4)
        else
          call upack(sname,struc,x1,x2,x3,x4,x5)
        endif
        call dscal(i0,1/fac,x1,1)
        if (is .ne. 0) then
          call pack5(sname,struc,is,x1,x2,x3,x4)
        else
          call pack5(sname,struc,x1,x2,x3,x4,x5)
        endif
      endif

      i0 = partks(categ,token,term,sname,struc,
     .  is,chr,count,lcast,istart,iopt)
      partkf = i0

      if (cast .ge. 100 .and. i0 .lt. count .and. noerr) then
        if (lcast .ne. 4) call rx('partkf not ready')
        if (is .ne. 0) then
          call upack(sname,struc,is,x1,x2,x3,x4)
        else
          call upack(sname,struc,x1,x2,x3,x4,x5)
        endif
        do  i = i0+1, count
          x1(i) = x1(i0)
        enddo
        i0 = count
        if (is .ne. 0) then
          call pack5(sname,struc,is,x1,x2,x3,x4)
        else
          call pack5(sname,struc,x1,x2,x3,x4,x5)
        endif
      endif

C ... Scale
      if ((noerr .or. optio .eq. 2) .and. fac.ne.1 .and. fac.ne.0) then
        if (optio .eq. 2) i0 = iabs(count)
        if (is .ne. 0) then
          call upack(sname,struc,is,x1,x2,x3,x4)
        else
          call upack(sname,struc,x1,x2,x3,x4,x5)
        endif
        call dscal(i0,fac,x1,1)
        if (is .ne. 0) then
          call pack5(sname,struc,is,x1,x2,x3,x4)
        else
          call pack5(sname,struc,x1,x2,x3,x4,x5)
        endif
      endif

      end

