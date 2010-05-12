Cdelw1 var ['opnu', '1', 'opnu', '-nlspc', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oqnu', '1', 'oqnu', '-3*nlspc', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['opp', '1', 'opp', '-6*nlspc', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['opprel', '1', 'opprel', '-32*nl*nl*max(nclspp,nspec)', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oqc', '1', 'oqc', '-nclasp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oqt', '1', 'oqt', '-nclasp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oinitc', '1', 'oinitc', '-nclasp', 'integer', 'iv_w_', 'clsprm']
Cdelw1 var ['ogroup', '1', 'ogroup', '-2*nclasp', 'integer', 'iv_w_', 'clsprm']
Cdelw1 var ['oves', '1', 'oves', '-nclspp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['ovdif', '1', 'ovdif', '-nclspp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['orhrmx', '1', 'orhrmx', '-nclspp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['ovrmax', '1', 'ovrmax', '-2*nclspp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oics', '-1', 'undef', '0', 'undef', 'x', 'clsprm']
Cdelw1 var ['ormax', '1', 'ormax', 'nclspp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['ovintr', '2', 'ovintr', '-nclasp*(nl*nsp)**2', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['opmpol', '2', 'opmpol', '-(2*nl-1)*nl**2*3*nsp*nclasp', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oqpp', '2', 'oqpp', '2*nqpp*4*nsp*nbas', 'real(8)', 'rv_w_', 'clsprm']
Cdelw1 var ['oqpp', '-1', 'undef', '0', 'undef', 'x', 'ioqpp']
Cdelw1 undel [['undef', 'wref'], 'oics', 'clsprm']
Cdelw1 undel [['undef', 'wref'], 'oqpp', 'ioqpp']
Cdelw1 undel [['wref'], 'ogroup', 'clsprm']
Cdelw1 undel [['wref'], 'ormax', 'clsprm']
Cdelw1 undel [['wref'], 'oinitc', 'clsprm']
Cdelw1 undel [['wref'], 'opnu', 'clsprm']
Cdelw1 undel [['wref'], 'oqnu', 'clsprm']
Cdelw1 undel [['wref'], 'oqc', 'clsprm']
Cdelw1 undel [['wref'], 'oqt', 'clsprm']
Cdelw1 undel [['wref'], 'oves', 'clsprm']
Cdelw1 undel [['wref'], 'ovdif', 'clsprm']
Cdelw1 undel [['wref'], 'ovrmax', 'clsprm']
Cdelw1 undel [['wref'], 'orhrmx', 'clsprm']
Cdelw1 undel [['wref'], 'opp', 'clsprm']
Cdelw1 undel [['wref'], 'ovintr', 'clsprm']
Cdelw1 undel [['wref'], 'opmpol', 'clsprm']
Cdelw1 undel [['wref'], 'oqpp', 'clsprm']
Cgetarg...info...           structure ['clsprm', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['clsprm', 'ham', 'sham']
Cgetarg...info...           structure ['clsprm', 'array', 'sarray']
Cgetarg...info...           structure ['clsprm', 'pot', 'spot']
Cgetarg...info...           structure ['clsprm', 'lat', 'slat']
Cgetarg...info...           structure ['clsprm', 'spec', 'sspec']
Cgetarg...info...           structure ['clsprm', 'bz', 'sbz']
Cgetarg...info...           use_to_add ['clsprm', 'm_struc_def']
      subroutine clsprm(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .  sstr,sstrn)

      use m_struc_def  !Cgetarg

C- Allocate memory for, and try to read, class-related data
C ----------------------------------------------------------------------
Ci    mode :1s digit
Ci         : 1 allocates ASA arrays (see spot,sham,sarray)
Ci         : 2 allocates TBE arrays:  oifrlx
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nclass nl nspin nspec loptc nmap
Ci     Stored:    *
Ci     Passed to: lgors ioqpp aiocls clsprp
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oics oipc oclabl opos
Ci     Stored:    ogroup ormax ohave oifrlx ozos
Ci     Passed to: aiocls clsprp
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: *
Ci     Stored:    osop oivso
Ci     Passed to: aiocls clsprp
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: ovrmax orhrmx opnu oqnu opp oves
Ci     Stored:    opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp orhos
Ci                osop ogrrme ovintr opmpol oqpp
Ci     Passed to: ioqpp aiocls clsprp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rmt idmod z lmxa
Ci     Stored:    *
Ci     Passed to: aiocls clsprp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat
Ci     Stored:    *
Ci     Passed to: aiocls clsprp
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: ndos
Ci     Stored:    *
Ci     Passed to: clsprp
Ci
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci     Elts read: n map
Ci     Stored:    *
Ci     Passed to: *
Ci   sstrn :struct for global strings
Ci     Elts read: map
Ci     Stored:    *
Ci     Passed to: clsprp
Cu Updates
Cu   08 Nov 07 (J. Xu) LDA+U implementation; qpp is complex
Cu   30 Sep 04 Reads/writes relativistic ppar's
Cu   02 Apr 04 SO array redimensioned; make and pack pot->nrhos
Cu   18 Mar 03 Allocate space for relativistic potential parameters
Cu   07 Feb 03 SO array redimensioned
Cu   19 Feb 02 entry clsprp copies class data to doubly padded layers
Cu   28 Apr 98 allocate space for radial matrix elements of grad
C ----------------------------------------------------------------------
C     implicit none
      integer mode
      character*(*) sstrn
Cgetarg       double precision sctrl(1),sarray(1),sham(1),spot(1),sspec(1),
Cgetarg      .  slat(1),sbz(1),sstr(1)
      type(s_ctrl)::sctrl
      type(s_array)::sarray
      type(s_ham)::sham
      type(s_pot)::spot
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_bz)::sbz
      type(s_str)::sstr

C Local variables
      logical lgors,cmdopt
      character*80 strn,clabl*8
Cdelw1 do not change  [['wref'], 'oinitc', 'clsprm']
Cdelw1 do not change  [['wref'], 'opnu', 'clsprm']
Cdelw1 do not change  [['wref'], 'oqnu', 'clsprm']
Cdelw1 do not change  [['wref'], 'oqc', 'clsprm']
Cdelw1 do not change  [['wref'], 'oqt', 'clsprm']
Cdelw1 do not change  [['wref'], 'oves', 'clsprm']
Cdelw1 do not change  [['wref'], 'ovdif', 'clsprm']
Cdelw1 do not change  [['wref'], 'orhrmx', 'clsprm']
Cdelw1 do not change  [['wref'], 'ovrmax', 'clsprm']
Cdelw1 do not change  [['wref'], 'ogroup', 'clsprm']
Cdelw1 do not change  [['wref'], 'ovintr', 'clsprm']
Cdelw1 do not change  [['wref'], 'opp', 'clsprm']
Cdelw1 do not change  [['undef', 'wref'], 'oics', 'clsprm']
Cdelw1 do not change  [['wref'], 'ormax', 'clsprm']
Cdelw1 do not change  [['wref'], 'opmpol', 'clsprm']
Cdelw1 do not change  [['wref'], 'oqpp', 'clsprm']
Cdelw1       integer obxc,oinitc,opnu,oqnu,oqc,oqt,oves,ovdif,orhrmx,ovrmax,
Cdelw1      .  orhos,osop,ogroup,ovintr,opp,opprel,oifrlx,oivso,oipc,oclabl,
Cdelw1      .  opos,oics,ozos,ormax,opmpol,oqpp,ogrrme,n0
      integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , ovdif 
     . , orhrmx , ovrmax , orhos , osop , ogroup , ovintr , opp , oifrlx 
     . , oivso , oipc , oclabl , opos , oics , ozos , ormax , opmpol 
     . , oqpp , ogrrme , n0
c       real(8) ,allocatable :: rv_w_opprel(:)

      parameter (n0=10)
      integer nbas,nclass,nspec,nl,nsp,nlspc,igets,nclasp,i1,i2,nmap,ic,
     .  is,ib,iclbsj,idmod(n0),iprint,lmx,nclspp,nqpp,lnsph,isw,neula,
     .  nbf,nrhos,nlibu
C     integer nmto,oppn
      double precision plat(9),z,xx,dval
C ... Heap
      integer w(1)
      common /w/ w

      integer,allocatable:: ics(:)
      integer:: i_dummy=-999
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg       call upack('ctrl nbas nclass nl nspin nspec',sctrl,nbas,nclass,nl,
Cgetarg      .  nsp,nspec)

      nbas=sctrl%nbas
      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec

Cgetarg...info...           integer :: nlibu , integer(8) :: sham%nlibu 1
Cgetarg       call upack('ham nlibu',sham,nlibu,0,0,0,0)

      nlibu=sham%nlibu

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg       call upack1('array nclasp',sarray,nclasp)

      nclasp=sarray%nclasp

      nclspp = max(2*nclasp-nclass,nspec)
      nlspc = nl*nsp*max(nclspp,nspec)

C ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
Cdelw1 do not change opnu because of ['wref']
        call defdr(opnu,-nlspc)
Cdelw1 do not change oqnu because of ['wref']
        call defdr(oqnu,-3*nlspc)
Cdelw1 do not change opp because of ['wref']
        call defdr(opp,-6*nlspc)
cCdelw1         call defdr(opprel,-32*nl*nl*max(nclspp,nspec))
c         allocate(rv_w_opprel(32*nl*nl*max(nclspp,nspec)))
c         rv_w_opprel(:)=0.0d0

Cdelw1 do not change oqc because of ['wref']
        call defdr(oqc,-nclasp)
Cdelw1 do not change oqt because of ['wref']
        call defdr(oqt,-nclasp)
Cdelw1 do not change oinitc because of ['wref']
        call defi (oinitc,-nclasp)
Cdelw1 do not change ogroup because of ['wref']
        call defi (ogroup,-2*nclasp)
Cdelw1 do not change oves because of ['wref']
        call defdr(oves,-nclspp)
Cdelw1 do not change ovdif because of ['wref']
        call defdr(ovdif,-nclspp)
Cdelw1 do not change orhrmx because of ['wref']
        call defdr(orhrmx,-nclspp)
Cdelw1 do not change ovrmax because of ['wref']
        call defdr(ovrmax,-2*nclspp)
C       call defdr(oqnu,-3*nlspc)
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg         call upack1('array oics',sarray,oics)

Cdelw1 warning(1) , probably  oics  is not defined yet at linenumber= 143
Cdelw1 do not change oics because of ['undef', 'wref']
        oics=sarray%oics

c takao sp2cls is replaced.
c        call sp2cls('spec rmt',sspec,w(oics),1,1,nclspp,ormax)
        allocate( ics(nclspp) )
Cdelw1 warning(1) , probably  oics  is not defined yet at linenumber= 148
Cdelw1 do not change oics because of ['undef', 'wref']
        call icopy(nclspp,w(oics),1,ics,1)
        print *,' qqqqq ics=',ics
Cdelw1 do not change ormax because of ['wref']
        call defdr(ormax,nclspp)
Cdelw1 do not change ormax because of ['wref']
        call dcopy(nclspp,sspec(ics(1:nclspp))%rmt,1,w(ormax),1)

Cgetarg...info...           integer :: ogroup , integer(8) :: sarray%ogroup 1
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg...info...           integer :: oinitc , integer(8) :: sarray%ohave 1
Cgetarg         call pack5('array ogroup ormax ohave',sarray,ogroup,ormax,
Cgetarg      .    oinitc,0,0)

Cdelw1 do not change ogroup because of ['wref']
        sarray%ogroup=ogroup
Cdelw1 do not change ormax because of ['wref']
        sarray%ormax=ormax
Cdelw1 do not change oinitc because of ['wref']
        sarray%ohave=oinitc

Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: oqc , integer(8) :: spot%oqc 1
Cgetarg...info...           integer :: oqt , integer(8) :: spot%oqt 1
Cgetarg         call pack5('pot opnu oqnu oqc oqt',spot,opnu,oqnu,oqc,oqt,0)

Cdelw1 do not change opnu because of ['wref']
        spot%opnu=opnu
Cdelw1 do not change oqnu because of ['wref']
        spot%oqnu=oqnu
Cdelw1 do not change oqc because of ['wref']
        spot%oqc=oqc
Cdelw1 do not change oqt because of ['wref']
        spot%oqt=oqt

Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg         call pack5('pot oves ovdif ovrmax orhrmx opp',spot,oves,ovdif,
Cgetarg      .    ovrmax,orhrmx,opp)

Cdelw1 do not change oves because of ['wref']
        spot%oves=oves
Cdelw1 do not change ovdif because of ['wref']
        spot%ovdif=ovdif
Cdelw1 do not change ovrmax because of ['wref']
        spot%ovrmax=ovrmax
Cdelw1 do not change orhrmx because of ['wref']
        spot%orhrmx=orhrmx
Cdelw1 do not change opp because of ['wref']
        spot%opp=opp

C       Even in collinear case, bxc might be written to
c        call defdr(obxc,-3*nclasp)
Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg         call pack5('pot opprel obxc nrhos',spot,opprel,obxc,0,0,0)

        spot%opprel=i_dummy !opprel
c         spot%obxc=obxc
c         spot%nrhos=0

Cgetarg         if (lgors('ctrl lncol,-1',sctrl)) then
c
c takao Kino's bug--> his replied that iand(-1,i)=i
c         if ( iand(-1,int(sctrl%lncol)) .ne.0 ) then

c         if ( sctrl%lncol /=0 ) then
Cgetarg...info...           integer :: neula , integer(8) :: sham%neula 1
Cgetarg...info...           integer :: nbf , integer(8) :: sham%nbf 1
Cgetarg           call upack2('ham neula nbf',sham,neula,nbf)
c          neula=sham%neula
c          nbf=sham%nbf
c          nrhos = max(nl,neula,nbf)
c          call rlse(obxc)
c          call defdr(orhos,2*3*4*nrhos*nclasp)
c          call defdr(osop,-nlspc*nsp*9)
c          call defdr(obxc,-3*nclasp)
Cgetarg...info...           integer :: nrhos , integer(8) :: spot%nrhos 1
Cgetarg...info...           integer :: orhos , integer(8) :: spot%orhos 1
Cgetarg...info...           integer :: obxc , integer(8) :: spot%obxc 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg           call pack5('pot nrhos orhos obxc osop',spot,
Cgetarg      .      nrhos,orhos,obxc,osop,0)
c          spot%nrhos=nrhos
c          spot%orhos=orhos
c          spot%obxc=obxc
c          spot%osop= osop
c        endif
        spot%osop=i_dummy !osop


        ogrrme = 0
Cgetarg         if (igets('ctrl loptc',sctrl) .gt. 0)
Cgetarg      .    call defrr(ogrrme,-16*nl*nsp*nclasp)
c         if ( int(sctrl%loptc) .gt. 0 ) call defrr ( ogrrme , -16 * nl
c     .   * nsp * nclasp )

Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg         call pack1('pot ogrrme',spot,ogrrme)

        spot%ogrrme=ogrrme

Cgetarg         if (lgors('ctrl lsx,3',sctrl) .or.
Cgetarg      .      lgors('ctrl lscr,-1',sctrl)) then
c         if ( iand(3,int(sctrl%lsx)) .ne.0 .or. iand(-1,int(sctrl%lscr)) .ne.0
c     .   ) then
c         if ( iand(-1,int(sctrl%lscr)) .ne.0  ) then
        if ( int(sctrl%lscr) .ne.0  ) then
Cdelw1 do not change ovintr because of ['wref']
          call defdr(ovintr,-nclasp*(nl*nsp)**2)
Cgetarg...info...           integer :: ovintr , integer(8) :: spot%ovintr 1
Cgetarg           call pack1('pot ovintr',spot,ovintr)
Cdelw1 do not change ovintr because of ['wref']
          spot%ovintr=ovintr

        endif
Cgetarg         lnsph = isw(lgors('ctrl lasa,32',sctrl))
c         lnsph = isw ( iand(32,int(sctrl%lasa)) .ne.0 )
ctakao remove asa
        lnsph = isw ( .false.)

C   ... Allocate memory for, attempt to read qpp from file
        if (lnsph .ne. 0 .or. nlibu .ne. 0) then
Cdelw1 do not change opmpol because of ['wref']
          call defrr(opmpol,-(2*nl-1)*nl**2*3*nsp*nclasp)
          i1 = nl**2
          nqpp = (i1*(i1+1))/2
Cdelw1 do not change oqpp because of ['wref']
          call defrr(oqpp,2*nqpp*4*nsp*nbas)
Cdelw1 do not change oqpp because of ['wref']
          call dvset(w(oqpp),1,1,-1d0)
Cgetarg...info...           integer :: opmpol , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oqpp , integer(8) :: spot%oqpp 1
Cgetarg           call pack2('pot opmpol oqpp',spot,opmpol,oqpp)
Cdelw1 do not change opmpol because of ['wref']
          spot%opmpol=opmpol
Cdelw1 do not change oqpp because of ['wref']
          spot%oqpp=oqpp
          call ioqpp(.false.,sctrl,spot)
        endif

Cgetarg...info...           real(8) :: plat(9) , real(8) :: slat%plat 9
Cgetarg         call upack1('lat plat',slat,plat)

C$$$         i_copy_size=size(slat%plat)
C$$$         call dcopy(i_copy_size,slat%plat,1,plat,1)


C$$$C   --- Create maps ---
C$$$Cgetarg         nmap = igets('ctrl nmap',sctrl)
C$$$         nmap = int(sctrl%nmap)



C$$$cccccccccccccccccc
C$$$         print *,' xxxxxx 222 nmap=',nmap
C$$$cccccccccccccccccc



C$$$        if (nmap .gt. 0) then
C$$$          call upacks('strn map',i1,i2)
C$$$Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
C$$$Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
C$$$Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
C$$$Cgetarg           call upack('array oipc oclabl opos',sarray,oipc,oclabl,opos,0,
C$$$Cgetarg      .      0)

C$$$           oipc=sarray%oipc
C$$$           oclabl=sarray%oclabl
C$$$           opos=sarray%opos

C$$$Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
C$$$Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
C$$$Cgetarg           call upack('pot ovrmax orhrmx',spot,ovrmax,orhrmx,0,0,0)

C$$$           ovrmax=spot%ovrmax
C$$$           orhrmx=spot%orhrmx

C$$$Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
C$$$Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
C$$$Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
C$$$Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
C$$$Cgetarg           call upack('pot opnu oqnu opp oves',spot,opnu,oqnu,opp,oves,0)

C$$$           opnu=spot%opnu
C$$$           oqnu=spot%oqnu
C$$$           opp=spot%opp
C$$$           oves=spot%oves

C$$$          do  10  ic = 1, nclass
C$$$            xx = dval(w(oclabl),ic)
C$$$            call r8tos8(xx,clabl)
C$$$            is = w(oics+ic-1)
C$$$            ib = iclbsj(ic,w(oipc),-nbas,1)
C$$$Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
C$$$Cgetarg             call upack('spec idmod',sspec,is,idmod,0,0,0)

C$$$             i_copy_size=size(sspec(is)%idmod)
C$$$             call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1)

C$$$            if (ib .lt. 0) then
C$$$              if (iprint() .ge. 20) print *,
C$$$     .          'CLSPRM: no basis corresponding to class ' //clabl
C$$$              goto 10
C$$$            endif
C$$$Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
C$$$Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
C$$$Cgetarg             call upack('spec z lmxa',sspec,is,z,lmx,0,0)

C$$$             z=sspec(is)%z
C$$$             lmx=sspec(is)%lmxa

C$$$            call wrmap(sstrn(i1:i2),nmap,clabl,ic,ib,z,plat,w(opos))
C$$$Cgetarg             call mapq(clabl,ic-1,nl,nsp,lmx,idmod,w(opnu),w(oqnu),
C$$$Cgetarg      .        w(opp),w(opprel),w(orhrmx),w(ovrmax),
C$$$Cgetarg      .        lgors('ctrl lves,1',sctrl),w(oves),w(oinitc))
C$$$             call mapq ( clabl , ic-1 , nl , nsp , lmx , idmod , w ( opnu
C$$$     .       ) , w ( oqnu ) , w ( opp ) , w ( opprel ) , w ( orhrmx ) , w
C$$$     .       ( ovrmax ) , iand(1,int(sctrl%lves)) .ne.0 , w ( oves ) , w (
C$$$     .       oinitc ) )

C$$$   10     continue

C$$$          if (cmdopt('-spin1',6,0,strn))
C$$$     .      call spinav(0,nclass,nl,nsp,w(opnu),w(oqnu))
C$$$          if (cmdopt('-map0',5,0,strn)) then
C$$$            call aiocls(.true.,0,sarray,sctrl,sham,spot,sspec,slat,1,
C$$$     .        nclass)
C$$$            call rx0('mapped moments written')
C$$$          endif
C$$$        endif
      endif

C$$$C --- Allocate tbe arrays ---
C$$$      if (mod(mode,10) .eq. 2) then
C$$$        call defi(oifrlx,-4*nbas)
C$$$Cgetarg...info...           integer :: oifrlx , integer(8) :: sarray%oifrlx 1
C$$$Cgetarg         call pack5('array oifrlx',sarray,oifrlx,0,0,0,0)

C$$$         sarray%oifrlx=oifrlx

C$$$Cgetarg         if (lgors('ctrl lncol,-1',sctrl)) then
C$$$         if ( iand(-1,int(sctrl%lncol)) .ne.0 ) then

C$$$          call rx('clsprm: is this needed?')
C$$$          call defdr(osop,-nlspc*nsp*9)
C$$$          call defi (oivso,-nclasp*nl)
C$$$CKi          call pack5('ham osop oivso',sham,osop,oivso,0,0,0)
C$$$Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
C$$$Cgetarg...info...           integer :: oivso , integer(8) :: spot%oivso 1
C$$$Cgetarg           call pack5('pot osop oivso',spot,osop,oivso,0,0,0)

C$$$           spot%osop=osop
C$$$           spot%oivso=oivso

C$$$        endif
C$$$      endif

C$$$Cgetarg       if (lgors('ctrl lstonr,-1',sctrl)) then
C$$$      if ( iand(-1,int(sctrl%lstonr(1))) .ne.0 ) then

C$$$C        call defdr(oammx,-nclasp)
C$$$C        call defdr(ostni,-nclasp)
C$$$Cgetarg         call defdr(ozos,igets('bz ndos',sbz)*nl*nsp*nclasp)
C$$$         call defdr ( ozos , int(sbz%ndos) * nl * nsp * nclasp )

C$$$Cgetarg...info...           integer :: ozos , integer(8) :: sarray%ozos 1
C$$$Cgetarg         call pack5('array ozos',sarray,ozos,0,0,0,0)

C$$$         sarray%ozos=ozos

C$$$      endif

C --- Misc variables ---
C      oips = igets('array oips',sarray)
C      nla = 0
C      do  20  ib = 1, nbas
C        is = w(oips+ib-1)
C        lmx = igetss('spec lmxa',is,sspec)
C        nla = nla + (lmx+1)**2
C   20 continue
C      call pack1('ctrl nla',sctrl,nla)

C     call shstru('array',sarray,0,0)
C     stop

ccc #error, have return with len(w_varlist)>0 at line 417
Cdelw1 w_varlist remains: opnu oqnu opp opprel oqc oqt oinitc ogroup oves ovdif orhrmx ovrmax ormax ovintr opmpol oqpp
Cdelw1 w_varlistundel: [oics ogroup ormax oinitc opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp ovintr opmpol oqpp]
Cdelw1 w_varlist (undel), remains: [opnu oqnu opp oqc oqt oinitc ogroup oves ovdif orhrmx ovrmax ormax ovintr opmpol oqpp]
Cdelw1 w_varlist (del), remains: [opprel]
Cdelw1 not deallocate oqpp because of [wref]
Cdelw1 not deallocate opmpol because of [wref]
Cdelw1 not deallocate ovintr because of [wref]
Cdelw1 not deallocate ormax because of [wref]
Cdelw1 not deallocate ovrmax because of [wref]
Cdelw1 not deallocate orhrmx because of [wref]
Cdelw1 not deallocate ovdif because of [wref]
Cdelw1 not deallocate oves because of [wref]
Cdelw1 not deallocate ogroup because of [wref]
Cdelw1 not deallocate oinitc because of [wref]
Cdelw1 not deallocate oqt because of [wref]
Cdelw1 not deallocate oqc because of [wref]
Cdelw1 not deallocate opp because of [wref]
Cdelw1 not deallocate oqnu because of [wref]
Cdelw1 not deallocate opnu because of [wref]
Cdelw1       return
c       if (allocated(rv_w_opprel)) deallocate(rv_w_opprel)
      return


      entry clsprp(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .  sstr,sstrn)
C- Copy class data to doubly padded case

Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg       call upack('ctrl nclass nl nspin nspec',sctrl,nclass,nl,nsp,nspec,
Cgetarg      .  0)

      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg       call upack1('array nclasp',sarray,nclasp)

      nclasp=sarray%nclasp

c#error, have return with len(w_varlist)>0 at line 440
Cdelw1 w_varlist remains: oqpp opmpol ovintr ormax ovrmax orhrmx ovdif oves ogroup oinitc oqt oqc opprel opp oqnu opnu
Cdelw1 w_varlistundel: [oics ogroup ormax oinitc opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp ovintr opmpol oqpp]
Cdelw1 w_varlist (undel), remains: [oqpp opmpol ovintr ormax ovrmax orhrmx ovdif oves ogroup oinitc oqt oqc opp oqnu opnu]
Cdelw1 w_varlist (del), remains: [opprel]
Cdelw1 not deallocate opnu because of [wref]
Cdelw1 not deallocate oqnu because of [wref]
Cdelw1 not deallocate opp because of [wref]
Cdelw1 not deallocate oqc because of [wref]
Cdelw1 not deallocate oqt because of [wref]
Cdelw1 not deallocate oinitc because of [wref]
Cdelw1 not deallocate ogroup because of [wref]
Cdelw1 not deallocate oves because of [wref]
Cdelw1 not deallocate ovdif because of [wref]
Cdelw1 not deallocate orhrmx because of [wref]
Cdelw1 not deallocate ovrmax because of [wref]
Cdelw1 not deallocate ormax because of [wref]
Cdelw1 not deallocate ovintr because of [wref]
Cdelw1 not deallocate opmpol because of [wref]
Cdelw1 not deallocate oqpp because of [wref]
Cdelw1       if (nclasp .eq. nclass) return
      if ( nclasp .eq. nclass ) then
        if (allocated(rv_w_opprel)) deallocate(rv_w_opprel)
        return
      endif

      nclspp = max(2*nclasp-nclass,nspec)
C ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then

        i1 = nclasp-nclass
Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg         call upack('pot opnu oqnu',spot,opnu,oqnu,0,0,0)

Cdelw1 do not change opnu because of ['wref']
        opnu=spot%opnu
Cdelw1 do not change oqnu because of ['wref']
        oqnu=spot%oqnu

Cdelw1 do not change opnu because of ['wref']
        call dpscop(w(opnu),w(opnu),nl*nsp*i1,
     .    nl*nsp*nclass+1,nl*nsp*nclasp+1,1d0)
Cdelw1 do not change oqnu because of ['wref']
        call dpscop(w(oqnu),w(oqnu),3*nl*nsp*i1,
     .    3*nl*nsp*nclass+1,3*nl*nsp*nclasp+1,1d0)

Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg...info...           integer :: ovdif , integer(8) :: spot%ovdif 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg         call upack('pot oves ovdif ovrmax orhrmx',spot,oves,ovdif,
Cgetarg      .    ovrmax,orhrmx,0)

Cdelw1 do not change oves because of ['wref']
        oves=spot%oves
Cdelw1 do not change ovdif because of ['wref']
        ovdif=spot%ovdif
Cdelw1 do not change ovrmax because of ['wref']
        ovrmax=spot%ovrmax
Cdelw1 do not change orhrmx because of ['wref']
        orhrmx=spot%orhrmx

Cdelw1 do not change oves because of ['wref']
        call dpscop(w(oves),w(oves),i1,nclass+1,nclasp+1,1d0)
Cdelw1 do not change ovdif because of ['wref']
        call dpscop(w(ovdif),w(ovdif),i1,nclass+1,nclasp+1,1d0)
Cdelw1 do not change ovrmax because of ['wref']
        call dpscop(w(ovrmax),w(ovrmax),i1,nclass+1,nclasp+1,1d0)
Cdelw1 do not change orhrmx because of ['wref']
        call dpscop(w(orhrmx),w(orhrmx),i1,nclass+1,nclasp+1,1d0)

      else
        call rxi('clsprp: bad mode',mode)
      endif

Cgetarg       end
Cdelw1 w_varlist remains: opnu oqnu opp opprel oqc oqt oinitc ogroup oves ovdif orhrmx ovrmax ormax ovintr opmpol oqpp
Cdelw1 w_varlistundel: [oics ogroup ormax oinitc opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp ovintr opmpol oqpp]
Cdelw1 w_varlist (undel), remains: [opnu oqnu opp oqc oqt oinitc ogroup oves ovdif orhrmx ovrmax ormax ovintr opmpol oqpp]
Cdelw1 w_varlist (del), remains: [opprel]
Cdelw1 not deallocate oqpp because of [wref]
Cdelw1 not deallocate opmpol because of [wref]
Cdelw1 not deallocate ovintr because of [wref]
Cdelw1 not deallocate ormax because of [wref]
Cdelw1 not deallocate ovrmax because of [wref]
Cdelw1 not deallocate orhrmx because of [wref]
Cdelw1 not deallocate ovdif because of [wref]
Cdelw1 not deallocate oves because of [wref]
Cdelw1 not deallocate ogroup because of [wref]
Cdelw1 not deallocate oinitc because of [wref]
Cdelw1 not deallocate oqt because of [wref]
Cdelw1 not deallocate oqc because of [wref]
Cdelw1 not deallocate opp because of [wref]
Cdelw1 not deallocate oqnu because of [wref]
Cdelw1 not deallocate opnu because of [wref]
      if (allocated(rv_w_opprel)) deallocate(rv_w_opprel)

      end subroutine clsprm


C================================================================
c take subs/ioqpp.F here
Cgetarg...info...           structure ['ioqpp', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['ioqpp', 'pot', 'spot']
Cgetarg...info...           use_to_add ['ioqpp', 'm_struc_def']
      subroutine ioqpp(lio,sctrl,spot)

      use m_struc_def  !Cgetarg

C- File i/o for phi-phi, phi-dot, dot-dot products
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lio: true for write, false for read
Ci          <0 write
Cu Updates
Cu   08 Nov 07 (J. Xu) qpp is complex
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lio
Cgetarg       double precision sctrl(1),spot(1)
      type(s_ctrl)::sctrl
      type(s_pot)::spot

C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
Cdelw1 do not change  [['undef', 'wref'], 'oqpp', 'ioqpp']
      integer i1,i2,nbas,nl,nsp,igets,ifi,oqpp,fopna,rdm,ipr

      call getpr(ipr)
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbas nl nspin',sctrl,nbas,nl,nsp,0,0)

      nbas=sctrl%nbas
      nl=sctrl%nl
      nsp=sctrl%nspin

      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
Cgetarg       oqpp = igets('pot oqpp',spot)
Cdelw1 warning(1) , probably  oqpp  is not defined yet at linenumber= 526
Cdelw1 do not change oqpp because of ['undef', 'wref']
      oqpp = int(spot%oqpp)

      ifi = fopna('qpp',-1,4+8)

      if (lio) then
Cdelw1 warning(1) , probably  oqpp  is not defined yet at linenumber= 531
Cdelw1 do not change oqpp because of ['undef', 'wref']
        call ywrm(1,'lmasa',1,ifi,' ',w(oqpp),1,i1,i1,i2)
        if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
        if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
        call pshpr(0)
Cdelw1 warning(1) , probably  oqpp  is not defined yet at linenumber= 536
Cdelw1 do not change oqpp because of ['undef', 'wref']
        if (rdm(ifi,2,i1*i2,' ',w(oqpp),i1,i2) .lt. 0) then
          if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
Cdelw1 warning(1) , probably  oqpp  is not defined yet at linenumber= 538
Cdelw1 do not change oqpp because of ['undef', 'wref']
          call dvset(w(oqpp),1,1,-1d0)
        endif
        call poppr
      endif

      call fclose(ifi)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oqpp]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine ioqpp


C================================================================
c take subs/spinav.F here
      subroutine spinav(mode,nclass,nl,nsp,pnu,qnu)
C- Averages up+down spin moments + pp's for all classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 average spins
Ci         :1 do not average, but exchange spins
Ci   nclass:number of inequivalent classes
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Co Outputs :moments are spin-averaged
Ci   pnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
Ci   qnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
Co   nsp   :set to 1 on output (mode=0)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Jan 06 Added mode
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nclass,nl,nsp,ic
      double precision pnu(nl,nsp,nclass),qnu(3,nl,nsp,nclass)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (nsp .eq. 1) return

      if (mode .gt. 1) then
        call rx('spinav: bad mode')
      elseif (mode .eq. 1) then
        do  ic = 1, nclass
          call dswap(nl,pnu(1,2,ic),1,pnu(1,1,ic),1)
          call dswap(3*nl,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
        enddo
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return
      endif

      do  ic = 1, nclass
        call daxpy(nl,1d0,pnu(1,2,ic),1,pnu(1,1,ic),1)
        call dscal(nl,.5d0,pnu(1,1,ic),1)
        call daxpy(3*nl,1d0,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
      enddo

      do  ic = 2, nclass
        call dcopy(nl,pnu(1,1,ic),1,pnu(1,ic,1),1)
        call dcopy(3*nl,qnu(1,1,1,ic),1,qnu(1,1,ic,1),1)
      enddo
      nsp = 1
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

