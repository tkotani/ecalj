Cdelw1 var ['olpp', '1', 'olpp', '100', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['olppr', '1', 'olppr', '3000', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['olsop', '1', 'olsop', '500', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['olgrme', '1', 'olgrme', '500', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['olva', '1', 'olva', '100', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['olmp', '1', 'olmp', '3000', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['opot', '1', 'opot', '3000', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['ocor', '1', 'ocor', '3000', 'real(8)', 'rv_w_', 'aiocls']
Cdelw1 var ['oics', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['oclabl', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['ohave', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['opnu', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['oqnu', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['orhrmx', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['ovrmax', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['opp', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['opprel', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['oves', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['ova', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['omp', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['osop', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 var ['ogrrme', '-1', 'undef', '0', 'undef', 'x', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'oics', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'oclabl', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'ohave', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'opnu', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'oqnu', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'orhrmx', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'ovrmax', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'opp', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'opprel', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'oves', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'ova', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'omp', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'osop', 'aiocls']
Cdelw1 undel [['undef', 'wref'], 'ogrrme', 'aiocls']
Cdelw1 undel [['wref'], 'onrc', 'aiocls']
Cdelw1 undel [['wref'], 'olppr', 'aiocls']
Cgetarg...info...           structure ['aiocls', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['aiocls', 'pot', 'spot']
Cgetarg...info...           structure ['aiocls', 'array', 'sarray']
Cgetarg...info...           structure ['aiocls', 'spec', 'sspec']
Cgetarg...info...           structure ['aiocls', 'lat', 'slat']
Cgetarg...info...           use_to_add ['aiocls', 'm_struc_def']
#define ONE_ATOM_UNIT 1
      subroutine aiocls(lio,mode,sarray,sctrl,sham,spot,sspec,slat,ic1,
     .  ic2)
       
       use m_struc_def  !Cgetarg

C- File I/O atomic data for classes ic1..ic2
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio    F for read, T for write
Ci   mode   ones digit
Ci          1 use data from first class corresponding to same species,
Ci            if data from own class is missing
Ci          2 like 1, but attempt to read data from disk anyway
Ci          3 make no attempt to read from disk but copy info from
Ci            first class corresponding to same species, if data missing
Ci          4 Add 4 if to use default P,Q when not otherwise supplied
Ci          10s digit
Ci          1 assemble background rho
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oclabl ohave oics onrc
Ci     Stored:
Ci     Passed to:
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nl nspin zbak nbas nclass
Ci     Stored:    zbak
Ci     Passed to: lgors
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci          (not used now)
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opp osop ogrrme opnu oqnu orhrmx ovrmax ovintr opmpol
Ci                oves
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idmod p q
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: vol
Ci     Stored:
Ci     Passed to:
Ci   ic1,ic2: range of classes to read data
Cr Remarks
Cr   Right now, read always takes data from file if available
Cu Updates
Cu   09 Nov 07 Corrected sign of default moment (paioc2)
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   26 Apr 03 Added MPI calls
Cu   07 Feb 03 adjusted for redimensioned sop
Cu   30 May 02 Assign better default P
Cu   28 Apr 98 I/O of radial matrix elements of grad
Cu   28 Sep 00 Added setting default P,Q
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lio
      integer mode,ic1,ic2
Cgetarg       double precision sarray(1),sctrl(1),sspec(1),slat(1),sham(1),
Cgetarg      .  spot(1)
       type(s_array)::sarray
       type(s_ctrl)::sctrl
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_ham)::sham
       type(s_pot)::spot

C heap
      integer w(1)
      common /w/ w
C Local variables
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,lgors,
     .  aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,
     .        readpq,readpp,readso,readop,readva,readmp
Cdelw1 do not change  [['undef', 'wref'], 'opnu', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'oqnu', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'orhrmx', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'ovrmax', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'opp', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'opprel', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'oves', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'ova', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'omp', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'osop', 'aiocls']
Cdelw1 do not change  [['wref'], 'olppr', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'oics', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'ohave', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'oclabl', 'aiocls']
Cdelw1 do not change  [['wref'], 'onrc', 'aiocls']
Cdelw1 do not change  [['undef', 'wref'], 'ogrrme', 'aiocls']
Cdelw1       integer opnu,oqnu,orhrmx,ovrmax,opp,opprel,oves,ova,omp,osop,olpp,
Cdelw1      .  olppr,olsop,olgrme,olva,olmp,opot,ocor,oics,ohave,oclabl,onrc,
Cdelw1      .  ogrrme,n0,nglob
       integer:: opnu , oqnu , orhrmx , ovrmax , opp , opprel , oves 
     . , ova , omp , osop  , oics , ohave , oclabl , onrc , ogrrme  !, olppr
     . , n0 , nglob 
       real(8) ,allocatable :: rv_w_olpp(:)
       real(8) ,allocatable :: rv_w_olsop(:)
       real(8) ,allocatable :: rv_w_olgrme(:)
       real(8) ,allocatable :: rv_w_olva(:)
       real(8) ,allocatable :: rv_w_olmp(:)
       real(8) ,allocatable :: rv_w_opot(:)
       real(8) ,allocatable :: rv_w_ocor(:)

      parameter(n0=10)
      integer ic,is,lmx,k,nl,nsp,ifi,jfi,fopn,lmxx,nspx,nrx,nr,
     .  idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .  nclspp,iprint,lgunit,mode0,mode00,lrel
      integer mpipid,procid
      double precision rhrmx,vrmax(2),ves,z,rmxx,ax,qc,dq,vrmxx(2),
     .  sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,dval,
     .  zbak(2),dgets,pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100)
ctakao
      real(8),allocatable:: rv_w_olppr

Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg       call upack('ctrl nl nspin nclass',sctrl,nl,nsp,nclass,0,0)
       
       nl=sctrl%nl
       nsp=sctrl%nspin
       nclass=sctrl%nclass

Cgetarg...info...           integer :: opp , integer(8) :: spot%opp 1
Cgetarg...info...           integer :: opprel , integer(8) :: spot%opprel 1
Cgetarg...info...           integer :: osop , integer(8) :: spot%osop 1
Cgetarg...info...           integer :: ogrrme , integer(8) :: spot%ogrrme 1
Cgetarg       call upack('pot opp opprel osop ogrrme',spot,opp,opprel,osop,
Cgetarg      .  ogrrme,0)
       
Cdelw1 warning(1) , probably  opp  is not defined yet at linenumber= 111
Cdelw1 do not change opp because of ['undef', 'wref']
       opp=spot%opp
Cdelw1 warning(1) , probably  opprel  is not defined yet at linenumber= 112
Cdelw1 do not change opprel because of ['undef', 'wref']
       opprel=spot%opprel
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 113
Cdelw1 do not change osop because of ['undef', 'wref']
       osop=spot%osop
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 114
Cdelw1 do not change ogrrme because of ['undef', 'wref']
       ogrrme=spot%ogrrme

Cgetarg...info...           integer :: opnu , integer(8) :: spot%opnu 1
Cgetarg...info...           integer :: oqnu , integer(8) :: spot%oqnu 1
Cgetarg...info...           integer :: orhrmx , integer(8) :: spot%orhrmx 1
Cgetarg...info...           integer :: ovrmax , integer(8) :: spot%ovrmax 1
Cgetarg       call upack('pot opnu oqnu orhrmx ovrmax',spot,opnu,oqnu,orhrmx,
Cgetarg      .  ovrmax,0)
       
Cdelw1 warning(1) , probably  opnu  is not defined yet at linenumber= 123
Cdelw1 do not change opnu because of ['undef', 'wref']
       opnu=spot%opnu
Cdelw1 warning(1) , probably  oqnu  is not defined yet at linenumber= 124
Cdelw1 do not change oqnu because of ['undef', 'wref']
       oqnu=spot%oqnu
Cdelw1 warning(1) , probably  orhrmx  is not defined yet at linenumber= 125
Cdelw1 do not change orhrmx because of ['undef', 'wref']
       orhrmx=spot%orhrmx
Cdelw1 warning(1) , probably  ovrmax  is not defined yet at linenumber= 126
Cdelw1 do not change ovrmax because of ['undef', 'wref']
       ovrmax=spot%ovrmax

Cgetarg...info...           integer :: ova , integer(8) :: spot%ovintr 1
Cgetarg...info...           integer :: omp , integer(8) :: spot%opmpol 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg       call upack('pot ovintr opmpol oves',spot,ova,omp,oves,0,0)
       
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 133
Cdelw1 do not change ova because of ['undef', 'wref']
       ova=spot%ovintr
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 134
Cdelw1 do not change omp because of ['undef', 'wref']
       omp=spot%opmpol
Cdelw1 warning(1) , probably  oves  is not defined yet at linenumber= 135
Cdelw1 do not change oves because of ['undef', 'wref']
       oves=spot%oves

Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: ohave , integer(8) :: sarray%ohave 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg       call upack('array nclasp oclabl ohave oics',sarray,nclasp,oclabl,
Cgetarg      .  ohave,oics,0)
       
       nclasp=sarray%nclasp
Cdelw1 warning(1) , probably  oclabl  is not defined yet at linenumber= 145
Cdelw1 do not change oclabl because of ['undef', 'wref']
       oclabl=sarray%oclabl
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 146
Cdelw1 do not change ohave because of ['undef', 'wref']
       ohave=sarray%ohave
Cdelw1 warning(1) , probably  oics  is not defined yet at linenumber= 147
Cdelw1 do not change oics because of ['undef', 'wref']
       oics=sarray%oics

      lrel = nglob('lrel')

      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)

C ... MPI: only master does sphere program
      procid = mpipid(1)
      if (procid .eq. 0) then

Cdelw1 warning(1) , probably  olppr  is not defined yet at linenumber= 158
Cdelw1 do not change olppr because of ['wref']
c      olppr = 1
Cdelw1       call defdr(olpp,100)
       allocate(rv_w_olpp(100))
Ctakao_ZeroClear_NotRequiered        if (100<0) rv_w_olpp(:)=0.0d0

Cdelw1 do not change olppr because of ['wref']
c      call defdr(olppr,3000)
       allocate(rv_w_olppr(3000))
c

Cdelw1       call defdr(olsop,500)
       allocate(rv_w_olsop(500))
Ctakao_ZeroClear_NotRequiered        if (500<0) rv_w_olsop(:)=0.0d0

Cdelw1       call defdr(olgrme,500)
       allocate(rv_w_olgrme(500))
Ctakao_ZeroClear_NotRequiered        if (500<0) rv_w_olgrme(:)=0.0d0

Cdelw1       call defdr(olva,100)
       allocate(rv_w_olva(100))
Ctakao_ZeroClear_NotRequiered        if (100<0) rv_w_olva(:)=0.0d0

Cdelw1       call defdr(olmp,3000)
       allocate(rv_w_olmp(3000))
Ctakao_ZeroClear_NotRequiered        if (3000<0) rv_w_olmp(:)=0.0d0

Cdelw1       call defdr(opot,3000)
       allocate(rv_w_opot(3000))
Ctakao_ZeroClear_NotRequiered        if (3000<0) rv_w_opot(:)=0.0d0

Cdelw1       call defdr(ocor,3000)
       allocate(rv_w_ocor(3000))
Ctakao_ZeroClear_NotRequiered        if (3000<0) rv_w_ocor(:)=0.0d0

      i2 = ic2
      if (i2 .eq. 0) i2 = nclasp
c     call awrit2('%n:1i',' ',100,6,nclasp,w(oics))
      do  10  ic = ic1, i2
Cdelw1 1 is digit, dropped
Cdelw1 warning(1) , probably  oics  is not defined yet at linenumber= 171
Cdelw1 do not change oics because of ['undef', 'wref']
        is = w(oics+ic-1)
Cdelw1 warning(1) , probably  oics  is not defined yet at linenumber= 172
Cdelw1 do not change oics because of ['undef', 'wref']
        icmap = iclbsj(is,w(oics),-nclasp,1)
        if (icmap .eq. ic .or. mode00 .eq. 0) icmap = 0
Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
Cgetarg...info...           real(8) :: pdf(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: qdf(n0,2) , real(8) :: sspec%q 20
Cgetarg         call upack('spec lmxa idmod p q',sspec,is,lmx,idmod,pdf,qdf)
         
         lmx=sspec(is)%lmxa
         i_copy_size=size(sspec(is)%idmod) 
         call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
         i_copy_size=size(sspec(is)%p) 
         call dcopy(i_copy_size,sspec(is)%p,1,pdf,1) 
         i_copy_size=size(sspec(is)%q) 
         call dcopy(i_copy_size,sspec(is)%q,1,qdf,1) 


        call dpzero(pnuloc,100)
        call dpzero(qnuloc,100)

C   --- Open the atom file ---
        outs1 = ' '
        outs2 = ' '
        outs3 = ' '
Cdelw1 warning(1) , probably  oclabl  is not defined yet at linenumber= 196
Cdelw1 do not change oclabl because of ['undef', 'wref']
        call r8tos8(dval(w(oclabl),ic),clabl)
        if (mode00 .ne. 3) then
#if ONE_ATOM_UNIT
          ifi = fopn(clabl)
#else
         ifi = fopna(clabl,30+ic,0)
#endif
        endif

C  --- Copy what is passed through to holding arrays ---
        k = nl*nsp
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 207
Cdelw1 do not change ohave because of ['undef', 'wref']
        call pvaioc(w(ohave),1,ic,icmap,havepq,jc)
        readpq = .not. havepq .or. mode00 .eq. 2 .and. jc .ne. ic
        readpq = readpq .and. mode00 .lt. 3
        if (havepq) then
Cdelw1 warning(1) , probably  opnu  is not defined yet at linenumber= 211
Cdelw1 do not change opnu because of ['undef', 'wref']
          call dpscop(w(opnu),pnuloc,k,1+(jc-1)*k,1,1d0)
Cdelw1 warning(1) , probably  oqnu  is not defined yet at linenumber= 212
Cdelw1 do not change oqnu because of ['undef', 'wref']
          call dpscop(w(oqnu),qnuloc,3*k,1+(jc-1)*3*k,1,1d0)
Cdelw1 warning(1) , probably  orhrmx  is not defined yet at linenumber= 213
Cdelw1 do not change orhrmx because of ['undef', 'wref']
          call dpscop(w(orhrmx),rhrmx,1,jc,1,1d0)
Cdelw1 warning(1) , probably  ovrmax  is not defined yet at linenumber= 214
Cdelw1 do not change ovrmax because of ['undef', 'wref']
          call dpscop(w(ovrmax),vrmax,2,2*jc-1,1,1d0)
          if (ic.ne.jc) call awrit0('%a pq,',outs1,len(outs1),0)
        endif
        k = 6*nl*nsp
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 218
Cdelw1 do not change ohave because of ['undef', 'wref']
        call pvaioc(w(ohave),2,ic,icmap,havepp,jc)
        readpp = .not. havepp .or. mode00 .eq. 2 .and. jc .ne. ic
        readpp = readpp .and. mode00 .lt. 3
        if (havepp) then
Cdelw1 warning(1) , probably  opp  is not defined yet at linenumber= 222
Cdelw1 do not change opp because of ['undef', 'wref']
Cdelw1 warning(1) , probably  opp  is not defined yet at linenumber= 222
Cdelw1 do not change opp because of ['undef', 'wref']
Cdelw1           call dpscop(w(opp),w(olpp),k,1+(jc-1)*k,1,1d0)
           call dpscop ( w ( opp ) , rv_w_olpp , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (lrel .eq. 2) then
            k = 4*nl*2*nl*2*2
Cdelw1 do not change olppr because of ['wref']
Cdelw1 warning(1) , probably  opprel  is not defined yet at linenumber= 225
Cdelw1 do not change opprel because of ['undef', 'wref']
            call dpscop(w(opprel),rv_w_olppr,k,1+(jc-1)*k,1,1d0)
          endif
Cdelw1 warning(1) , probably  oves  is not defined yet at linenumber= 227
Cdelw1 do not change oves because of ['undef', 'wref']
          call dpscop(w(oves),ves,1,jc,1,1d0)
          if (ic.ne.jc) call awrit0('%a pp,',outs1,len(outs1),0)
        endif
        k = (nl*nsp)**2
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 231
Cdelw1 do not change ohave because of ['undef', 'wref']
        call pvaioc(w(ohave),8,ic,icmap,haveva,jc)
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 232
Cdelw1 do not change ova because of ['undef', 'wref']
        haveva = haveva .and. ova .gt. 0
        readva = .not. haveva .or. mode00 .eq. 2 .and. jc .ne. ic
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 234
Cdelw1 do not change ova because of ['undef', 'wref']
        readva = readva .and. ova .gt. 0 .and. mode00 .lt. 3
        if (haveva) then
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 236
Cdelw1 do not change ova because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 236
Cdelw1 do not change ova because of ['undef', 'wref']
Cdelw1           call dpscop(w(ova),w(olva),k,1+(jc-1)*k,1,1d0)
           call dpscop ( w ( ova ) , rv_w_olva , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a va,',outs1,len(outs1),0)
        endif
        k = nl**2*(2*nl-1)*3*nsp
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 240
Cdelw1 do not change ohave because of ['undef', 'wref']
        call pvaioc(w(ohave),16,ic,icmap,havemp,jc)
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 241
Cdelw1 do not change omp because of ['undef', 'wref']
        havemp = havemp .and. omp .gt. 0
        readmp = .not. havemp .or. mode00 .eq. 2 .and. jc .ne. ic
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 243
Cdelw1 do not change omp because of ['undef', 'wref']
        readmp = readmp .and. omp .gt. 0 .and. mode00 .lt. 3
        if (havemp) then
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 245
Cdelw1 do not change omp because of ['undef', 'wref']
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 245
Cdelw1 do not change omp because of ['undef', 'wref']
Cdelw1           call dpscop(w(omp),w(olmp),k,1+(jc-1)*k,1,1d0)
           call dpscop ( w ( omp ) , rv_w_olmp , k , 1 + ( jc - 1 ) * k 
     .     , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a mp,',outs1,len(outs1),0)
        endif
        k = nl*nsp*nsp*9
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 249
Cdelw1 do not change ohave because of ['undef', 'wref']
        call pvaioc(w(ohave),4,ic,icmap,haveso,jc)
        haveso = haveso .and. osop.gt. 0
        readso = .not. haveso .or. mode00 .eq. 2 .and. jc .ne. ic
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 252
Cdelw1 do not change osop because of ['undef', 'wref']
        readso = readso .and. osop .gt. 0 .and. mode00 .lt. 3
        if (haveso) then
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 254
Cdelw1 do not change osop because of ['undef', 'wref']
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 254
Cdelw1 do not change osop because of ['undef', 'wref']
Cdelw1           call dpscop(w(osop),w(olsop),k,1+(jc-1)*k,1,1d0)
           call dpscop ( w ( osop ) , rv_w_olsop , k , 1 + ( jc - 1 ) * 
     .     k , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a sop,',outs1,len(outs1),0)
        endif
        k = 16*nl*nsp
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 258
Cdelw1 do not change ohave because of ['undef', 'wref']
        call pvaioc(w(ohave),32,ic,icmap,haveop,jc)
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 259
Cdelw1 do not change ogrrme because of ['undef', 'wref']
        haveop = haveop .and. ogrrme .gt. 0
        readop = .not. haveop .or. mode00 .eq. 2 .and. jc .ne. ic
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 261
Cdelw1 do not change ogrrme because of ['undef', 'wref']
        readop = readop .and. ogrrme .gt. 0 .and. mode00 .lt. 3
        if (haveop) then
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 263
Cdelw1 do not change ogrrme because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 263
Cdelw1 do not change ogrrme because of ['undef', 'wref']
Cdelw1           call dpscop(w(ogrrme),w(olgrme),k,1+(jc-1)*k,1,1d0)
           call dpscop ( w ( ogrrme ) , rv_w_olgrme , k , 1 + ( jc - 1 ) 
     .     * k , 1 , 1d0 ) 

          if (ic.ne.jc) call awrit0('%a opp,',outs1,len(outs1),0)
        endif

C   --- File WRITE ---
        if (lio) then
          lgen = .false.
          lpot = .false.
          lcor = .false.
C     ... Pick up GEN and POT, if available, to save again
          if (scat(iabs(ifi),'GEN:',':',.true.)) then
            lgen = aiogen(alabel,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,dq,
     .        vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
          endif
          if (scat(iabs(ifi),'POT:',':',.true.)) then
            read(ifi,102) nr,nsp,a,rmax
  102       format(2i5,2f12.5)
Cdelw1             lpot = aiopot(nr,nsp,a,rmax,-99d0,w(opot),ifi)
             lpot = aiopot ( nr , nsp , a , rmax , - 99d0 , rv_w_opot , ifi 
     .       ) 

          endif
Cdelw1           lcor = aiocor(nr,nsp,a,rmxx,w(ocor),sumec,sumtc,ifi)
           lcor = aiocor ( nr , nsp , a , rmxx , rv_w_ocor , sumec , sumtc 
     .     , ifi ) 


          rewind ifi
          jfi = -ifi
          if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .      dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
          if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,idmod,
     .      nl,lmx,nsp,rhrmx,vrmax,jfi)
Cdelw1 do not change olppr because of ['wref']
Cdelw1 do not change olppr because of ['wref']
Cdelw1           if (havepp) sw = aiopar(clabl,lrel,w(olpp),w(olppr),ves,nl,
Cdelw1      .      lmx,nsp,jfi)
           if ( havepp ) sw = aiopar ( clabl , lrel , rv_w_olpp , 
     .       rv_w_olppr , ves , nl , lmx , nsp , jfi ) 

Cdelw1           if (haveva) sw = aiova(clabl,w(olva),nl,lmx,nsp,jfi)
           if ( haveva ) sw = aiova ( clabl , rv_w_olva , nl , lmx , nsp 
     .     , jfi ) 

Cdelw1           if (havemp) sw = aiomp(clabl,w(olmp),nl,2*nl-2,nsp,jfi)
           if ( havemp ) sw = aiomp ( clabl , rv_w_olmp , nl , 2 * nl - 
     .     2 , nsp , jfi ) 

Cdelw1           if (haveso) sw = aiosop(clabl,w(olsop),nl,lmx,nsp,jfi)
           if ( haveso ) sw = aiosop ( clabl , rv_w_olsop , nl , lmx , nsp 
     .     , jfi ) 

Cdelw1           if (haveop) sw = aiorme(clabl,w(olgrme),nl,nsp,jfi)
           if ( haveop ) sw = aiorme ( clabl , rv_w_olgrme , nl , nsp , 
     .     jfi ) 

Cdelw1           if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,w(opot),jfi)
           if ( lpot ) sw = aiopot ( nr , nsp , a , rmax , - 99d0 , rv_w_opot 
     .     , jfi ) 

Cdelw1           if (lcor) lcor = aiocor(nr,nsp,a,rmxx,w(ocor),sumec,sumtc,jfi)
           if ( lcor ) lcor = aiocor ( nr , nsp , a , rmxx , rv_w_ocor , 
     .     sumec , sumtc , jfi ) 


C   --- File READ ---
        else

C     ... Copy whatever is available on disk to holding arrays
          if (readpq .or. mode0 .ge. 4) then
            if (readpq) rewind ifi
            if (readpq) readpq = aiomom(clabl,pnuloc,qnuloc,idmod,
     .        nl,lmx,nsp,rhrmx,vrmax,ifi)
C           Couldn't read from atom file ; take default values
            if (readpq) call awrit0('%a pq,',outs2,len(outs2),0)
            if (mode0 .ge. 4 .and. .not. (readpq .or. havepq)) then
C             call dmcpy(pdf,n0,1,pnuloc,nl,1,nl,nsp)
              call paioc2(nsp,nl,n0,pdf,qdf,pnuloc,qnuloc)
              call awrit0('%a pq,',outs3,len(outs2),0)
              call dvset(vrmax,1,2,-.7d0)
              rhrmx = .1d0
              readpq = .true.
            endif
          endif
          if (readpp) then
Cdelw1 do not change olppr because of ['wref']
Cdelw1 do not change olppr because of ['wref']
Cdelw1             readpp = aiopar(clabl,lrel,w(olpp),w(olppr),ves,nl,lmx,nsp,
Cdelw1      .        ifi)
             readpp = aiopar ( clabl , lrel , rv_w_olpp , rv_w_olppr , ves 
     .       , nl , lmx , nsp , ifi ) 

            if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
          endif
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 323
Cdelw1 do not change osop because of ['undef', 'wref']
          if (readso .and. osop .gt. 0) then
Cdelw1             readso = aiosop(clabl,w(olsop),nl,lmx,nsp,ifi)
             readso = aiosop ( clabl , rv_w_olsop , nl , lmx , nsp , ifi ) 

            if (readso) call awrit0('%a so,',outs2,len(outs2),0)
          endif
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 327
Cdelw1 do not change ogrrme because of ['undef', 'wref']
          if (readop .and. ogrrme .gt. 0) then
Cdelw1             readop = aiorme(clabl,w(olgrme),nl,nsp,ifi)
             readop = aiorme ( clabl , rv_w_olgrme , nl , nsp , ifi ) 

            if (readop) call awrit0('%a op,',outs2,len(outs2),0)
          endif
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 331
Cdelw1 do not change ova because of ['undef', 'wref']
          if (readva .and. ova .gt. 0) then
Cdelw1             readva = aiova(clabl,w(olva),nl,lmx,nsp,ifi)
             readva = aiova ( clabl , rv_w_olva , nl , lmx , nsp , ifi ) 

            if (readva) call awrit0('%a va,',outs2,len(outs2),0)
          endif
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 335
Cdelw1 do not change omp because of ['undef', 'wref']
          if (readmp .and. omp .gt. 0) then
Cdelw1             readmp = aiomp(clabl,w(olmp),nl,2*nl-2,nsp,ifi)
             readmp = aiomp ( clabl , rv_w_olmp , nl , 2 * nl - 2 , nsp , 
     .       ifi ) 

            if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
          endif

C     ... Update what parameters are available
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 341
Cdelw1 do not change ohave because of ['undef', 'wref']
          w(ohave+ic-1) = isw(haveop.or.readop)*32+
     .                    isw(havemp.or.readmp)*16+
     .                    isw(haveva.or.readva)*8 +
     .                    isw(haveso.or.readso)*4 +
     .                    isw(havepp.or.readpp)*2 +
     .                    isw(havepq.or.readpq)*1 +
     .                    w(ohave+ic-1) - bitand(w(ohave+ic-1),63)

          k = nl*nsp
          if (havepq .or. readpq) then
Cdelw1 warning(1) , probably  opnu  is not defined yet at linenumber= 351
Cdelw1 do not change opnu because of ['undef', 'wref']
            call dpscop(pnuloc,w(opnu),k,1,1+(ic-1)*k,1d0)
Cdelw1 warning(1) , probably  oqnu  is not defined yet at linenumber= 352
Cdelw1 do not change oqnu because of ['undef', 'wref']
            call dpscop(qnuloc,w(oqnu),3*k,1,1+(ic-1)*3*k,1d0)
Cdelw1 warning(1) , probably  orhrmx  is not defined yet at linenumber= 353
Cdelw1 do not change orhrmx because of ['undef', 'wref']
            call dpscop(rhrmx,w(orhrmx),1,1,ic,1d0)
Cdelw1 warning(1) , probably  ovrmax  is not defined yet at linenumber= 354
Cdelw1 do not change ovrmax because of ['undef', 'wref']
            call dpscop(vrmax,w(ovrmax),2,1,2*ic-1,1d0)
          endif
          k = 6*nl*nsp
          if (havepp .or. readpp) then
Cdelw1 warning(1) , probably  opp  is not defined yet at linenumber= 358
Cdelw1 do not change opp because of ['undef', 'wref']
Cdelw1 warning(1) , probably  opp  is not defined yet at linenumber= 358
Cdelw1 do not change opp because of ['undef', 'wref']
Cdelw1             call dpscop(w(olpp),w(opp),k,1,1+(ic-1)*k,1d0)
             call dpscop ( rv_w_olpp , w ( opp ) , k , 1 , 1 + ( ic - 1 ) 
     .       * k , 1d0 ) 

            if (lrel .eq. 2) then
              k = 4*nl*2*nl*2*2
Cdelw1 do not change olppr because of ['wref']
Cdelw1 warning(1) , probably  opprel  is not defined yet at linenumber= 361
Cdelw1 do not change opprel because of ['undef', 'wref']
              call dpscop(rv_w_olppr,w(opprel),k,1,1+(ic-1)*k,1d0)
            endif
Cdelw1 warning(1) , probably  oves  is not defined yet at linenumber= 363
Cdelw1 do not change oves because of ['undef', 'wref']
            call dpscop(ves,w(oves),1,1,ic,1d0)
          endif
          k = (nl*nsp)**2
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 366
Cdelw1 do not change ova because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 366
Cdelw1 do not change ova because of ['undef', 'wref']
Cdelw1           if ((haveva .or. readva) .and. ova .gt. 0)
Cdelw1      .      call dpscop(w(olva),w(ova),k,1,1+(ic-1)*k,1d0)
           if ( ( haveva .or. readva ) .and. ova .gt. 0 ) call dpscop ( 
     .     rv_w_olva , w ( ova ) , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 

          k = nl**2*(2*nl-1)*3*nsp
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 369
Cdelw1 do not change omp because of ['undef', 'wref']
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 369
Cdelw1 do not change omp because of ['undef', 'wref']
Cdelw1           if ((havemp .or. readmp) .and. omp .gt. 0)
Cdelw1      .      call dpscop(w(olmp),w(omp),k,1,1+(ic-1)*k,1d0)
           if ( ( havemp .or. readmp ) .and. omp .gt. 0 ) call dpscop ( 
     .     rv_w_olmp , w ( omp ) , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 

          k = nl*nsp*nsp*9
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 372
Cdelw1 do not change osop because of ['undef', 'wref']
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 372
Cdelw1 do not change osop because of ['undef', 'wref']
Cdelw1           if ((haveso .or. readso) .and. osop .gt. 0)
Cdelw1      .      call dpscop(w(olsop),w(osop),k,1,1+(ic-1)*k,1d0)
           if ( ( haveso .or. readso ) .and. osop .gt. 0 ) call dpscop ( 
     .     rv_w_olsop , w ( osop ) , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 

          k = 16*nl*nsp
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 375
Cdelw1 do not change ogrrme because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 375
Cdelw1 do not change ogrrme because of ['undef', 'wref']
Cdelw1           if ((haveop .or. readop) .and. ogrrme .gt. 0)
Cdelw1      .      call dpscop(w(olgrme),w(ogrrme),k,1,1+(ic-1)*k,1d0)
           if ( ( haveop .or. readop ) .and. ogrrme .gt. 0 ) call dpscop 
     .     ( rv_w_olgrme , w ( ogrrme ) , k , 1 , 1 + ( ic - 1 ) * k , 1d0 
     .     ) 

        endif

        if (mode00 .lt. 3) call fclr(clabl,ifi)
        if (iprint() .gt. 40) then
        outs = ' '
        if (outs1 .ne. ' ') then
          call awrit1('%x '//clabl//'%a: copied '//outs1//'%a%b from '//
     .      'class %i',outs,len(outs),0,jc)
        endif
        if (outs2 .ne. ' ') then
          call awrit0('%x '//clabl//'%a: read '//outs2//'%a%b from '//
     .      'disk',outs,len(outs),0)
        endif
        if (outs3 .ne. ' ') then
          if (outs2 .eq.  ' ') then
            call awrit0('%x '//clabl//'%a: use defaults for: '//outs3//
     .        '%a%b',outs,len(outs),0)
          else
            call awrit0('%a; use defaults for: '//outs3//'%a%b',
     .      outs,len(outs),0)
          endif
        endif
        if (outs .eq. ' ')
     .    call awrit0(' '//clabl//'%a: nothing read',outs,len(outs),0)
        call awrit0(' aiocls class'//outs,' ',-len(outs),lgunit(1))
       endif
   10 continue
C     call awrit2('%n:1i',' ',100,6,nclasp,w(ohave))
Cdelw1 rlse name= olpp old_list= olpp olppr olsop olgrme olva olmp opot ocor 
Cdelw1 rlse name= olpp new_list= (None)
Cdelw1 not deallocate olppr because of [wref]
Cdelw1       call rlse(olpp)
       if (allocated(rv_w_ocor)) deallocate(rv_w_ocor)
       if (allocated(rv_w_opot)) deallocate(rv_w_opot)
       if (allocated(rv_w_olmp)) deallocate(rv_w_olmp)
       if (allocated(rv_w_olva)) deallocate(rv_w_olva)
       if (allocated(rv_w_olgrme)) deallocate(rv_w_olgrme)
       if (allocated(rv_w_olsop)) deallocate(rv_w_olsop)
       if (allocated(rv_w_olpp)) deallocate(rv_w_olpp)
c       call rlse(olppr)

C     End of MPI master-only branch
      endif

C --- MPI broadcast everything passed out of aiocls ---
C     sarray and contents of sarray pointers
Cdelw1 warning(1) , probably  ohave  is not defined yet at linenumber= 411
Cdelw1 do not change ohave because of ['undef', 'wref']
      call mpibc1(w(ohave),nclasp,2,.false.,'aiocls','initc')
C     sctrl and contents
      call mpibc1(zbak,2,4,.false.,'aiocls','zbak')
C     spot and and contents of spot pointers
Cdelw1 warning(1) , probably  opp  is not defined yet at linenumber= 415
Cdelw1 do not change opp because of ['undef', 'wref']
      call mpibc1(w(opp),6*nl*nsp*nclasp,4,.false.,'aiocls','pp')
      if (lrel .eq. 2) then
        k = 4*nl*2*nl*2*2
Cdelw1 warning(1) , probably  opprel  is not defined yet at linenumber= 418
Cdelw1 do not change opprel because of ['undef', 'wref']
        call mpibc1(w(opprel),k*nclasp,4,.false.,'aiocls','pp')
      endif
C     if (igets('ctrl lncol',sctrl) .gt. 0) then
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 421
Cdelw1 do not change osop because of ['undef', 'wref']
      if (osop .gt. 0) then
Cdelw1 warning(1) , probably  osop  is not defined yet at linenumber= 422
Cdelw1 do not change osop because of ['undef', 'wref']
       call mpibc1(w(osop),9*nl*nsp*nsp*nclasp,4,.false.,'aiocls','sop')
      endif
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 424
Cdelw1 do not change ogrrme because of ['undef', 'wref']
      if (ogrrme .gt. 0) then
Cdelw1 warning(1) , probably  ogrrme  is not defined yet at linenumber= 425
Cdelw1 do not change ogrrme because of ['undef', 'wref']
       call mpibc1(w(ogrrme),8*nl*nsp*nclasp,4,.false.,'aiocls','grrme')
      endif
Cdelw1 warning(1) , probably  opnu  is not defined yet at linenumber= 427
Cdelw1 do not change opnu because of ['undef', 'wref']
      call mpibc1(w(opnu),nl*nsp*nclasp,4,.false.,'aiocls','pnu')
Cdelw1 warning(1) , probably  oqnu  is not defined yet at linenumber= 428
Cdelw1 do not change oqnu because of ['undef', 'wref']
      call mpibc1(w(oqnu),3*nl*nsp*nclasp,4,.false.,'aiocls','qnu')
Cdelw1 warning(1) , probably  orhrmx  is not defined yet at linenumber= 429
Cdelw1 do not change orhrmx because of ['undef', 'wref']
      call mpibc1(w(orhrmx),nclasp,4,.false.,'aiocls','rhrmx')
      nclspp = 2*nclasp-nclass
Cdelw1 warning(1) , probably  ovrmax  is not defined yet at linenumber= 431
Cdelw1 do not change ovrmax because of ['undef', 'wref']
      call mpibc1(w(ovrmax),2*nclspp,4,.false.,'aiocls','vrmax')
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 432
Cdelw1 do not change ova because of ['undef', 'wref']
      if (ova .gt. 0) then
        k = nclasp*(nl*nsp)**2
Cdelw1 warning(1) , probably  ova  is not defined yet at linenumber= 434
Cdelw1 do not change ova because of ['undef', 'wref']
        call mpibc1(w(ova),k,4,.false.,'aiocls','vintr')
      endif
C     if (lgors('ctrl lasa,32',sctrl)) then
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 437
Cdelw1 do not change omp because of ['undef', 'wref']
      if (omp .gt. 0) then
        k = (2*nl-1)*nl**2*3*nsp*nclasp
Cdelw1 warning(1) , probably  omp  is not defined yet at linenumber= 439
Cdelw1 do not change omp because of ['undef', 'wref']
        call mpibc1(w(omp),k,4,.false.,'aiocls','pmpol')
      endif
Cdelw1 warning(1) , probably  oves  is not defined yet at linenumber= 441
Cdelw1 do not change oves because of ['undef', 'wref']
      call mpibc1(w(oves),nclspp,4,.false.,'aiocls','ves')

C --- Assemble background rho (rmt) ---
Cgetarg       if (mod(mode/10,10) .ne. 0 .and. lgors('ctrl lasa,64',sctrl)) then
ctakao assumes lasa64=F (no asa mode) Thus this is always false.
c       if ( mod ( mode/10 , 10 ) .ne. 0 .and. iand(64,int(sctrl%lasa)) .ne.0 
c     . ) then 
       if (.false.) then


Cgetarg...info...           real(8) :: zbak(2) , real(8) :: sctrl%zbak 2
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg         call upack('ctrl zbak nbas nclass',sctrl,zbak,nbas,nclass,0,0)
         
         i_copy_size=size(sctrl%zbak) 
         call dcopy(i_copy_size,sctrl%zbak,1,zbak,1) 
         nbas=sctrl%nbas
         nclass=sctrl%nclass

Cgetarg...info...           integer :: onrc , integer(8) :: sarray%onrc 1
Cgetarg         call upack1('array onrc',sarray,onrc)
         
         onrc=sarray%onrc

C        call awrit2('%n:1d',' ',100,6,nclass,w(orhrmx))
C        call awrit2('%n:1i',' ',100,6,nclasp,w(onrc))
        if (zbak(2) .eq. 0) then
          do  15  ic = 1, nclass
Cdelw1 warning(1) , probably  oclabl  is not defined yet at linenumber= 470
Cdelw1 do not change oclabl because of ['undef', 'wref']
            z = dval(w(oclabl),ic)
            call r8tos8(z,clabl)
Cdelw1 warning(1) , probably  orhrmx  is not defined yet at linenumber= 472
Cdelw1 do not change orhrmx because of ['undef', 'wref']
            z = dval(w(orhrmx),ic)
            if (iprint() .ge. 30 .and. z .eq. 0) print *,
     .        ' (AIOCLS, mtcor): rho(rmax) for class ',clabl,' is 0'
Cdelw1 1 is digit, dropped
            zbak(2) = zbak(2) + (w(onrc+ic-1)*z)/nbas
   15     continue
Cgetarg           zbak(2) = zbak(2) * dgets('lat vol',slat)
           zbak ( 2 ) = zbak ( 2 ) * slat%vol 

        endif
Cgetarg...info...           real(8) :: zbak(2) , real(8) :: sctrl%zbak 2
Cgetarg         call pack1('ctrl zbak',sctrl,zbak)
         
         i_copy_size=size(sctrl%zbak) 
         call dcopy(i_copy_size,zbak,1,sctrl%zbak,1) 

      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oics oclabl ohave opnu oqnu orhrmx ovrmax opp opprel oves ova omp osop ogrrme onrc olppr]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine aiocls 

      subroutine pvaioc(initc,mask,ic0,icmap,lhave,ic)

C- Find whether data avaiable either in class or mapped class
C     implicit none
      logical lhave
      integer initc(1),mask,ic0,icmap
      integer ic
      ic = ic0
      lhave = mod(initc(ic)/mask,2) .eq. 1
      if (.not. lhave .and. icmap .ne. 0) then
        lhave = mod(initc(icmap)/mask,2) .eq. 1
        ic = icmap
      endif
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvaioc 

      subroutine paioc2(nsp,nl,n0,pat,qat,pnu,qnu)

C- Widget to copy pat,qat to pnu,qnu
C     implicit none
      integer n0,nl,nsp
      double precision pat(n0,2),qat(n0,2),pnu(nl,nsp),qnu(3,nl,nsp)
      integer i,il

      do  i = 1, nsp
      do  il = 1, nl
C       pnu(il,i) = int(pat(il,i)) + .5d0
        pnu(il,i) = pat(il,i)
        qnu(1,il,i) = qat(il,1)/nsp
        if (nsp .eq. 2) then
          if (pat(il,i) .eq. 0) pnu(il,i) = pnu(il,1)
          qnu(1,il,i) = qat(il,1)/nsp + qat(il,2)/2*dble(3-2*i)
        endif
        qnu(2,il,i) = 0d0
        qnu(3,il,i) = 0d0
      enddo
      enddo
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine paioc2 


