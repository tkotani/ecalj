Cgetarg...info...           structure ['lmaux', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['lmaux', 'array', 'sarray']
Cgetarg...info...           structure ['lmaux', 'ham', 'sham']
Cgetarg...info...           structure ['lmaux', 'pot', 'spot']
Cgetarg...info...           structure ['lmaux', 'lat', 'slat']
Cgetarg...info...           structure ['lmaux', 'str', 'sstr']
Cgetarg...info...           structure ['lmaux', 'spec', 'sspec']
Cgetarg...info...           structure ['lmaux', 'site', 'ssite']
Cgetarg...info...           use_to_add ['lmaux', 'm_struc_def']
      subroutine lmaux(prgnam,sctrl,sham,spot,slat,sspec,ssite,sstr,
     .  sarray,sbz,slabl,sstrn,mode,wksize)
       
       use m_struc_def  !Cgetarg

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
      implicit none
C ... Passed parameters
      integer mode,wksize
      character*(*) sstrn, prgnam*8
      character*8 slabl(*)
Cgetarg       double precision sarray(1),sctrl(1),sham(1),spot(1),slat(1),
Cgetarg      .  sspec(1),ssite(1),sstr(1),sbz(1)
       type(s_array)::sarray
       type(s_ctrl)::sctrl
       type(s_ham)::sham
       type(s_pot)::spot
       type(s_lat)::slat
       type(s_spec)::sspec(*)
       type(s_site)::ssite(*)
       type(s_str)::sstr
       type(s_bz)::sbz

       type(s_site),allocatable::ssite2(:)

C ... Local parameters
      character*120 outs,fnam*8
      integer NULLI
      logical cmdopt,T,F,swtmp
      parameter (T=.true., F=.false., NULLI=-99999)
      integer fopn,getdig,i,ip,j,k,m,ifi,iprint,lpbc,lgunit,lncol,
     .  nbas,nclasp,nclass,nl,nlspc,nsp,modep(3),mxnbr,parg,nbasp,
     .  nbaspp,nkd,nkq,nspec,npadl,npadr,neul,nc,mxcsiz,nttab,igets,
     .  iosits,cmplat,ngrp,ival,stdo,irs(5),fopna,nclspp,bitand,igetss,
     .  ngmx,nsgrp
      integer oeold,oeula,oinitc,oipc,olmx,onrc,opnu,opold,
     .  opp,oqc,oqnu,oqt,orhrmx,ormax,oves,ovrmax,owk,owk2,oclabl,
     .  opos,opos2,oqlv,odlv,orham,oamsh,onrmsh,oalpha,oaamom,
     .  onpr,ormt,olock,olockc,os,oz,ozz,ormx,oip,oips,oips2,opgfsl,
     .  ontab,oiax,og,oag,ozc,mxclas,ormtc,oics
      double precision xv(10),xx,alat,plat(3,3),plat2(9),facrmx,facrng,
     .  dval,avw,ekap(2),enu,rmaxs,qss(4),ckbas,cksumf,ehterm(4),
     .  qlat(9),emad,trumad,vmtz(2),omax1(3),omax2(3),wsrmax
      parameter (ngmx=48,mxclas=1000)

      integer:: i_copy_size, i_spackv, i_spacks
C ... Heap
      integer w(1)
      common /w/ w

ctakao
      integer:: ifx
      integer,allocatable:: ics(:)
      integer,allocatable:: lmxa(:)
c      real(8),allocatable::idmod(:,:), z(:),amsh(:),rmax(:), nrmsh(:)
      real(8),allocatable:: z(:),rmax(:)

      print *,' lmaux:' 

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbas nclass nl nspec nspin',sctrl,nbas,nclass,nl,
Cgetarg      .  nspec,nsp)
       
       nbas=sctrl%nbas
       nclass=sctrl%nclass
       nl=sctrl%nl
       nspec=sctrl%nspec
       nsp=sctrl%nspin

Cgetarg...info...           integer :: lncol , integer(8) :: sctrl%lncol 1
Cgetarg...info...           integer :: modep(3) , integer(8) :: sctrl%modep 3
Cgetarg       call upack('ctrl lncol modep',sctrl,lncol,modep,0,0,0)
       
       lncol=sctrl%lncol
       i_copy_size=size(sctrl%modep) 
       call i8icopy(i_copy_size,sctrl%modep,1,modep,1) 

Cgetarg       lpgf = igets('ctrl lpgf',sctrl)
c       lpgf = int(sctrl%lpgf(1)) 

C     lpbc = 0 for pbc in 3 dimensions, 11 pgf padded geometry
      lpbc = 0
c      if (lpgf .gt. 0 .and. lpgf .lt. 10) lpbc = 11
Cgetarg...info...           integer :: oinitc , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: onrc , integer(8) :: sarray%onrc 1
Cgetarg       call upack('array ohave oics oclabl onrc',sarray,oinitc,oics,
Cgetarg      .  oclabl,onrc,0)
       
       oinitc=sarray%ohave
       oics=sarray%oics
       oclabl=sarray%oclabl
       onrc=sarray%onrc
ccc

Cgetarg...info...           integer :: npadl , integer(8) :: sarray%npadl 1
Cgetarg...info...           integer :: npadr , integer(8) :: sarray%npadr 1
Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: opgfsl , integer(8) :: sarray%opgfsl 1
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg       call upack('array npadl npadr nclasp opgfsl ormax',sarray,npadl,
Cgetarg      .  npadr,nclasp,opgfsl,ormax)
       
       npadl=sarray%npadl
       npadr=sarray%npadr
       nclasp=sarray%nclasp
       opgfsl=sarray%opgfsl
       ormax=sarray%ormax

Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
Cgetarg...info...           integer :: oips , integer(8) :: sarray%oips 1
Cgetarg       call upack('array opos oipc oips',sarray,opos,oipc,oips,0,0)
       
       opos=sarray%opos
       oipc=sarray%oipc
       oips=sarray%oips

Cgetarg...info...           integer :: oeula , integer(8) :: sham%oeula 1
Cgetarg...info...           integer :: neul , integer(8) :: sham%neula 1
Cgetarg...info...           real(8) :: qss(4) , real(8) :: sham%qss 4
Cgetarg       call upack('ham oeula neula qss',sham,oeula,neul,qss,0,0)
       
       oeula=sham%oeula
       neul=sham%neula
       i_copy_size=size(sham%qss) 
       call dcopy(i_copy_size,sham%qss,1,qss,1) 

Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg...info...           integer :: oqc , integer(8) :: spot%oqc 1
Cgetarg...info...           integer :: oqt , integer(8) :: spot%oqt 1
Cgetarg       call upack('pot opnu oqnu opp oqc oqt',spot,opnu,oqnu,opp,oqc,oqt)
       
       opnu=spot%opnu
       oqnu=spot%oqnu
       opp=spot%opp
       oqc=spot%oqc
       oqt=spot%oqt

Cgetarg...info...           integer :: oaamom , integer(8) :: spot%oaamom 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg       call upack('pot oaamom orhrmx oves ovrmax',spot,oaamom,orhrmx,
Cgetarg      .  oves,ovrmax,0)
       
       oaamom=spot%oaamom
       orhrmx=spot%orhrmx
       oves=spot%oves
       ovrmax=spot%ovrmax

Cgetarg...info...           real(8) :: avw , real(8) :: slat%avw 1
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: plat2(9) , real(8) :: slat%plat2 9
Cgetarg       call upack('lat avw alat plat plat2',slat,avw,alat,plat,plat2,0)
       
       avw=slat%avw
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       i_copy_size=size(slat%plat2) 
       call dcopy(i_copy_size,slat%plat2,1,plat2,1) 

Cgetarg...info...           integer :: nkd , integer(8) :: slat%nkd 1
Cgetarg...info...           integer :: odlv , integer(8) :: slat%odlv 1
Cgetarg...info...           integer :: nkq , integer(8) :: slat%nkq 1
Cgetarg...info...           integer :: oqlv , integer(8) :: slat%oqlv 1
Cgetarg       call upack('lat nkd odlv nkq oqlv',slat,nkd,odlv,nkq,oqlv,0)
       
       nkd=slat%nkd
       odlv=slat%odlv
       nkq=slat%nkq
       oqlv=slat%oqlv

Cgetarg...info...           integer :: mxnbr , integer(8) :: sstr%mxnbr 1
Cgetarg...info...           real(8) :: rmaxs , real(8) :: sstr%rmax 1
Cgetarg       call upack2('str mxnbr rmax',sstr,mxnbr,rmaxs)
       
       mxnbr=sstr%mxnbr
       rmaxs=sstr%rmax

ctakao
c      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmx) --->lmxa
c      call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)
      allocate(ics(1:nclasp),lmxa(nclasp),z(nclasp))
      call icopy(nclasp,w(oics),1,ics,1)
      lmxa(1:nclasp) = sspec(ics(1:nclasp))%lmxa
      z   (1:nclasp) = sspec(ics(1:nclasp))%z
      print *,' lmxa=',lmxa
      print *,' z   =',z

      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas

      stdo = lgunit(1)

      j = 10
      if (cmdopt('--shorten',j-1,0,outs)) then
        call shorps(nbasp,plat,modep,w(opos),w(opos))
      endif

C ... Read from restart file
      if (cmdopt('--rs=',5,0,outs)) then
Cgetarg         irs(1) = igets('ctrl lrs,7',sctrl)
         irs ( 1 ) = iand(7,int(sctrl%lrs)) 

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
Cgetarg       mxcsiz = igets('str mxnbr',sstr)
       mxcsiz = int(sstr%mxnbr) 

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
      ifx=0
      if (cmdopt('--getwsr',j-1,0,outs)) then
        call info(10,1,0,' ... Make sphere radii',0,0)
C        xx = dglob('lrel',1d0,1)
C        xx = dglob('nsp',1d0,1)
C       Initial estimate for sphere radii: overlapping atom potentials
        call defrr(ozz,nspec)
Cgetarg         call spackv(10,'spec z',sspec,1,nspec,w(ozz))
         i_copy_size=1; 
         do i_spackv=1,nspec 
         call spackv_array_copy_r8_r8('u',sspec(i_spackv)%z,i_copy_size,i_spackv+1-1,w(ozz))
         enddo

        call defrr(ormt,nspec)
Cgetarg         call spackv(10,'spec rmt',sspec,1,nspec,w(ormt))
         i_copy_size=1; 
         do i_spackv=1,nspec 
         call spackv_array_copy_r8_r8('u',sspec(i_spackv)%rmt,i_copy_size,i_spackv+1-1,w(ormt))
         enddo

        call defi(olock,-nspec)
        do  i = 1, nspec
Cgetarg         call ivset(w(olock),i,i,bitand(igetss('spec mxcst',i,sspec),2))
         call ivset ( w ( olock ) , i , i , bitand ( int(sspec(i)%mxcst) 
     .   , 2 ) ) 

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
Cgetarg...info...           real(8) :: omax1(3) , real(8) :: sctrl%omax1 3
Cgetarg...info...           real(8) :: omax2(3) , real(8) :: sctrl%omax2 3
Cgetarg...info...           real(8) :: wsrmax , real(8) :: sctrl%wsrmax 1
Cgetarg         call upack('ctrl omax1 omax2 wsrmax',sctrl,omax1,omax2,wsrmax,0,
Cgetarg      .    0)
         
         i_copy_size=size(sctrl%omax1) 
         call dcopy(i_copy_size,sctrl%omax1,1,omax1,1) 
         i_copy_size=size(sctrl%omax2) 
         call dcopy(i_copy_size,sctrl%omax2,1,omax2,1) 
         wsrmax=sctrl%wsrmax

        call sclwsr(20,nbas,nbasp,nspec,alat,plat,w(opos),w(oips),modep,
     .    slabl,w(ozz),w(olock),1d0,wsrmax,omax1,omax2,w(ormt))
Cgetarg         call spackv(11,'spec rmt',sspec,1,nspec,w(ormt))
         i_copy_size=1; 
         do i_spackv=1,nspec 
         call spackv_array_copy_r8_r8('p',sspec(i_spackv)%rmt,i_copy_size,i_spackv+1-1,w(ormt))
         enddo

        nclspp = max(2*nclasp-nclass,nspec)
c takao
c        call sp2cls('spec rmt',sspec,ics,1,1,nclspp,ormax)
        allocate(rmax(nclspp))
        print *,' zzzz nclspp=',nclspp
        do i=1,nclspp
          rmax(i) =sspec(ics(i))%rmt
        enddo
        call defdr(ormax,nclspp)
        call dcopy(nclspp,rmax,1,w(ormax),1)

Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg         call pack1('array ormax',sarray,ormax)
         sarray%ormax=ormax 
         ifx=1
      endif

c-------
      if(ifx==0) then
        allocate(rmax(nclasp))
        call dcopy(nclasp,w(ormax),1,rmax,1)
      endif
      
c


C --- Show neighbors by shell ---
      outs = ' '
      j = 8
      if (cmdopt('--shell',j-1,0,outs)) then
        call shoshl(outs(j:),nbas,w(opos),plat,mxnbr,z,slabl,
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
Cgetarg         call spackv(10,'site pos',ssite,1,nbas,w(opos2))
         i_copy_size=size(ssite(1)%pos)
         do i_spackv=1,nbas 
         call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos2))
         enddo

Cgetarg         call spackv(10,'site spec',ssite,1,nbas,w(oips2))
         i_copy_size=1; 
         do i_spackv=1,nbas 
         call spackv_array_copy_i8_i('u',ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,w(oips2))
         enddo


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
      call ovmin(outs(i:),nbas,nbasp,alat,plat,rmax,rmax,
     .  w(oclabl),w(oipc),modep,z,ontab,oiax,w(opos),j)
C ... Write positions to file
      if (cmdopt('--wpos=',7,0,outs))
     .  call iopos(T,0,outs(8:),nbasp,w(opos))

C ... Inner products between Euler angles
      if (mod(lncol,2) .eq. 1) then
        call amagnc(nbas,nl,w(oipc),w(1),1,w(oqnu),w(oeula),neul,1,xv,
     .    w(oaamom))
        print '(1x)'
        call angtab(nbas,w(opos),alat,rmax,qss,w(oqlv),
     .    w(odlv),nkd,w(oipc),neul,w(oeula))
      endif

      endif

C$$$C --- Plane analysis branch ---
C$$$      if (getdig(mode,1,2) .ne. 0) then
C$$$Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
C$$$Cgetarg...info...           integer :: nbasp , integer(8) :: sctrl%nbasp 1
C$$$Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
C$$$Cgetarg         call upack('ctrl nbas nbasp nl',sctrl,nbas,nbasp,nl,0,0)
         
C$$$         nbas=sctrl%nbas
C$$$         nbasp=sctrl%nbasp
C$$$         nl=sctrl%nl

C$$$Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
C$$$Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
C$$$Cgetarg...info...           integer :: oqc , integer(8) :: spot%oqc 1
C$$$Cgetarg...info...           integer :: oqt , integer(8) :: spot%oqt 1
C$$$Cgetarg         call upack('pot opnu oqnu oqc oqt',spot,opnu,oqnu,oqc,oqt,0)
         
C$$$         opnu=spot%opnu
C$$$         oqnu=spot%oqnu
C$$$         oqc=spot%oqc
C$$$         oqt=spot%oqt

C$$$        call defrr(oves,nbaspp)
C$$$        call supot(1,sctrl,slat,spot)
C$$$        call asamad(sctrl,spot,sarray,slat,sspec,0,
C$$$     .    w(opnu),w(oqnu),0d0,w(oves),emad,trumad,vmtz)
C$$$        call plana(wksize,npadl,npadr,nbaspp,slabl,slat,sspec,ssite,
C$$$     .    w(oves),w(opnu),w(oqnu))
C$$$        return
C$$$      endif

C --- Generate input file to xbs program ---
C$$$      if (getdig(mode,2,2) .ne. 0) then
C$$$      ifi = fopn('XBS')
C$$$      facrmx = .5d0
C$$$      ip = 4
C$$$      if (cmdopt('-bs=',ip,0,outs)) then
C$$$        call skipbl(outs,len(outs),ip)
C$$$        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,facrmx)
C$$$      endif
C$$$      facrng = 1d0
C$$$      ip = 4
C$$$      if (cmdopt('-ss=',ip,0,outs)) then
C$$$        call skipbl(outs,len(outs),ip)
C$$$        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,facrng)
C$$$      endif
C$$$C ... Copy wsr*facrmx into rmax, wsr*facrng into rham (if nonzero)
C$$$      call defdr(ormx, nclass)
C$$$      call dpcopy(w(ormax),w(ormx),1,nclass,facrmx)
C$$$C ... Copy wsr*facrng into rham, if zero
C$$$      do  30  i = 1, nclass
C$$$      xx = dval(rham,i)
C$$$   30 if (xx .eq. 0d0 .or. xx .eq. NULLI)
C$$$     .  call dvset(rham,i,i,facrng*dval(w(ormax),i))
C$$$      if (iprint() .ge. 20) then
C$$$        call awrit2('%N ball size = %d * sphere size;  '//
C$$$     .    'def stick length = %d * sum sphere sizes',
C$$$     .    ' ',80,stdo,facrmx,facrng)
C$$$      endif
C$$$      if (cmdopt('-spec',5,0,outs) .or. cmdopt('--spec',6,0,outs)) then
C$$$        nc = nspec
C$$$        oip = oips
C$$$      else
C$$$        nc = nclass
C$$$        oip = oipc
C$$$      endif
C$$$      call ioxbs(ifi,nbas,nc,alat,plat,rham,w(ormx),w(oclabl),
C$$$     .  w(oip),z,w(opos))
C$$$      call fclose(ifi)
C$$$      endif

C$$$C --- Shift pp's (and optionally) moments by enu ---
C$$$C     pp's are remade using the potential if available.
C$$$C     use -enu=val to shift all enu's to val.
C$$$C     Use -mom if to save shifted moments.  Potential NOT remade.
C$$$      if (getdig(mode,3,2) .ne. 0) then
C$$$      ip = 5
C$$$      if (cmdopt('-enu=',ip,0,outs)) then
C$$$        call skipbl(outs,len(outs),ip)
C$$$        i = parg(' ',4,outs,ip,len(outs),' ',1,1,j,enu)
C$$$        if (i .eq. -1) call rxs('LMSHF: failed to parse ',outs)
C$$$      else
C$$$        call rx('LMSHF: missing argument -enu=val')
C$$$      endif
C$$$      call supot(1,sctrl,slat,spot)
C$$$      call defdr(oeold,3*nl*nsp*nclass)
C$$$      call defdr(opold,  nl*nsp*nclass)
C$$$      call dcopy(3*nl*nsp*nclass,w(oqnu),1,w(oeold),1)
C$$$      call dcopy(  nl*nsp*nclass,w(opnu),1,w(opold),1)
C$$$      call pp2enu(sspec,nclass,ics,w(onrc),nsp,nl,w(ormax),avw,
C$$$     .  w(oqnu),w(opnu),w(oqnu),xv,w(oinitc),w(oclabl),enu,w(opp))
C$$$      if (.not. cmdopt('-mom',4,0,outs)) then
C$$$        oqnu = oeold
C$$$        opnu  = opold
C$$$      endif
C$$$      call asvsph(sctrl,slat,sarray,sspec,sham,spot,0d0,2,ehterm,i)
C$$$      call rx0('done')
C$$$ctakao
C$$$c      call sp2cls('spec a',sspec,ics,1,1,nclass,oamsh)
C$$$c      call sp2cls('spec nr',sspec,ics,1,1,nclass,onrmsh)
C$$$      allocate(amsh(nclass),nrmsh(nclass))
C$$$      do i=1,nclass
C$$$        amsh(i) =sspec(ics(i))%a 
C$$$        nrmsh(i)=sspec(ics(i))%nr
C$$$        print *,' amsh nrmsh=',amsh(i),nrmsh(i)
C$$$      enddo

C$$$      call prmx('a',amsh,nclass,nclass,1)
C$$$      call mkppar(w(oclabl),nl,nsp,nclass,lmxa,z,w(ormax),avw,
C$$$     .  amsh,nrmsh,w(opnu),idmod,w(oves),w(oqnu),w(opp))
C$$$      call rx('implement call to write out new parms')
C$$$C      call aiomms(nclass,nl,lmxa,nsp,w(oclabl),w(opnu),w(oqnu),
C$$$C     .  idmoc,w(opp),w(orhrmx),w(ovrmax),w(oves),.true.)
C$$$      endif

C --- Interpolate core to another mesh ---
      if (getdig(mode,4,2) .ne. 0) then
        call rx('patch clabl for call to coritp')
C       call coritp(nclass,nsp,w(oclabl),nrmsh,amsh,w(ormax))
      endif

C$$$C --- Display poles of potential functions ---
C$$$      if (getdig(mode,5,2) .ne. 0) then
C$$$      ckbas = cksumf(w(opos),3*nbas)
C$$$      call iostr(8,'STR',nl,nbasp,1,ekap,0,ckbas,-1,i,oalpha,
C$$$     .  oiax,onpr,os)
C$$$      call defdr(owk,nlspc)
C$$$      call pptrns(0,nl,w(oipc),nclass,nsp,w(oalpha),1,w(opp),w(owk))
C$$$      call rlse(owk)
C$$$      call defdr(owk2,3*2*nl**2*nbasp*nsp)
C$$$      call shopol(nl,nclass,nsp,sspec,ics,w(opp))
C$$$      endif

C$$$C --- Import data in other formats ---
C$$$      if (getdig(mode,6,2) .ne. 0) then

C$$$ctakao
C$$$        allocate( ssite2(nbasp))
C$$$        do i=1,nbasp
C$$$        call usite_init(ssite2(i))
C$$$        enddo
C$$$        ssite2(1:nbasp) = ssite(1:nbasp)
C$$$c        j = nint(ssite(1)%size) * nbasp
C$$$c        call defrr(ossite,-100*j)
C$$$c        call dcopy(j,ssite,1,ssite2,1)
C$$$        call aioxtn(sctrl,sspec,ssite2,sarray,slat,slabl,w(opnu),
C$$$     .    w(oqnu))

C$$$C       Output to restart file
C$$$        if (cmdopt('--rs=',5,0,outs)) then
C$$$Cgetarg           irs(2) = igets('ctrl lrs,8+16',sctrl)/8
C$$$           irs ( 2 ) = iand(8+16,int(sctrl%lrs)) /8 

C$$$          if (irs(2) .gt. 0) then
C$$$            ifi = fopna('rsta',-1,0)
C$$$            call asars(0,sctrl,ssite2,sspec,slat,spot,sbz,
C$$$     .        w(opnu),w(oqnu),.false.,-ifi)
C$$$            call fclr('rsta',ifi)
C$$$          endif
C$$$        endif
C$$$        deallocate(ssite2)
C$$$      endif

C --- Empty sphere finder ---
      if (getdig(mode,7,2) .ne. 0) then
Cgetarg...info...           integer :: nsgrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: og , integer(8) :: slat%osymgr 1
Cgetarg...info...           integer :: oag , integer(8) :: slat%oag 1
Cgetarg         call upack('lat nsgrp osymgr oag',slat,nsgrp,og,oag,0,0)
         
         nsgrp=slat%nsgrp
         og=slat%osymgr
         oag=slat%oag

        call defrr(oz,nspec)
Cgetarg         call spackv(10,'spec z',sspec,1,nspec,z)
         i_copy_size=1; 
         do i_spackv=1,nspec 
         call spackv_array_copy_r8_r8('u',sspec(i_spackv)%z,i_copy_size,i_spackv+1-1,z)
         enddo

        call defrr(ormt,nspec)
Cgetarg         call spackv(10,'spec rmt',sspec,1,nspec,w(ormt))
         i_copy_size=1; 
         do i_spackv=1,nspec 
         call spackv_array_copy_r8_r8('u',sspec(i_spackv)%rmt,i_copy_size,i_spackv+1-1,w(ormt))
         enddo

        call defi(olock,-nspec)
        do  i = 1, nspec
Cgetarg         call ivset(w(olock),i,i,bitand(igetss('spec mxcst',i,sspec),2))
         call ivset ( w ( olock ) , i , i , bitand ( int(sspec(i)%mxcst) 
     .   , 2 ) ) 

        enddo
        call defrr(ozc,mxclas)
        call defrr(ormtc,mxclas)
        call defrr(olockc,mxclas)
        owk = onrc
        call defi(onrc,mxclas)
        call icopy(nclass,w(owk),1,w(onrc),1)

        if (nclass .ge. mxclas) call rx('lmaux: increase mxclas')
        call spec2c(nspec,nclass,ics,w(ormt),w(ormtc),z,w(ozc),
     .    w(olock),w(olockc))
        call findes(sctrl,alat,w(oclabl),nbas,nclass,nl,w(onrc),mxclas,
     .    nsgrp,opos,oipc,plat,w(og),w(oag),w(olockc),w(ormtc),w(ozc))

      endif
ctakao
      deallocate(ics,lmxa,z)
Cgetarg       end
       end subroutine lmaux 


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

Cgetarg       end
       end subroutine spec2c 

        

