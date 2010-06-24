Cchp1 cmdvar= real(8),rv_p_,ogrrme
Cchp1 cmdvar= integer,iv_p_,ohave
Cchp1 cmdvar= real(8),rv_p_,ovrmax
Cchp1 cmdvar= real(8),rv_p_,ovintr
Cchp1 cmdvar= real(8),rv_p_,oves
Cchp1 cmdvar= real(8),rv_p_,ovdif
Cchp1 cmdvar= integer,iv_p_,oics
Cchp1 cmdvar= real(8),rv_p_,oqp
      subroutine clsprm(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .sstr,sstrn)

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
Cchp1       integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , ovdif 
Cchp1      ., orhrmx , ovrmax , orhos , osop , ogroup , ovintr , opp , oifrlx 
Cchp1      ., oivso , oipc , oclabl , opos , oics , ozos , ormax , opmpol 
Cchp1      ., oqpp , ogrrme , n0
Cchp1        integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , ovdif 
Cchp1      . , orhrmx , ovrmax , orhos , osop , ogroup , ovintr , opp , oifrlx 
Cchp1      . , oivso , oipc , oclabl , opos , ozos , ormax , opmpol , oqpp 
Cchp1      . , ogrrme , n0 
Cchp1         integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , oves , orhrmx 
Cchp1      .  , ovrmax , orhos , osop , ogroup , ovintr , opp , oifrlx , oivso 
Cchp1      .  , oipc , oclabl , opos , ozos , ormax , opmpol , oqpp , ogrrme 
Cchp1      .  , n0 
Cchp1          integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
Cchp1      .   ovrmax , orhos , osop , ogroup , ovintr , opp , oifrlx , oivso 
Cchp1      .   , oipc , oclabl , opos , ozos , ormax , opmpol , oqpp , ogrrme 
Cchp1      .   , n0 
Cchp1           integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
Cchp1      .    ovrmax , orhos , osop , ogroup , opp , oifrlx , oivso , oipc 
Cchp1      .    , oclabl , opos , ozos , ormax , opmpol , oqpp , ogrrme , n0 
Cchp1            integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
Cchp1      .     orhos , osop , ogroup , opp , oifrlx , oivso , oipc , oclabl 
Cchp1      .     , opos , ozos , ormax , opmpol , oqpp , ogrrme , n0 
            integer:: obxc , oinitc , opnu , oqnu , oqc , oqt , orhrmx , 
     .      orhos , osop , ogroup , opp , oifrlx , oivso , oipc , oclabl 
     .      , opos , ozos , ormax , opmpol , oqpp , n0 
            real(8),pointer :: rv_p_ogrrme(:)

           real(8),pointer :: rv_p_ovrmax(:)

          real(8),pointer :: rv_p_ovintr(:)

         real(8),pointer :: rv_p_oves(:)

        real(8),pointer :: rv_p_ovdif(:)

       integer,pointer :: iv_p_oics(:)

c       real(8) ,allocatable :: pprel_rv(:)

      parameter (n0=10)
      integer nbas,nclass,nspec,nl,nsp,nlspc,igets,nclasp,i1,i2,nmap,ic,
     .is,ib,iclbsj,idmod(n0),iprint,lmx,nclspp,nqpp,lnsph,isw,neula,
     .nbf,nrhos,nlibu
C     integer nmto,oppn
      double precision plat(9),z,xx,dval
C ... Heap
      integer w(1)
      common /w/ w

      integer,allocatable:: ics(:)
      integer:: i_dummy=-999

      nbas=sctrl%nbas
      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec


      nlibu=sham%nlibu


      nclasp=sarray%nclasp

      nclspp = max(2*nclasp-nclass,nspec)
      nlspc = nl*nsp*max(nclspp,nspec)

C ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
        call defdr(opnu,-nlspc)
        call defdr(oqnu,-3*nlspc)
        call defdr(opp,-6*nlspc)
cCdelw1         call defdr(opprel,-32*nl*nl*max(nclspp,nspec))
c         allocate(pprel_rv(32*nl*nl*max(nclspp,nspec)))
c         pprel_rv(:)=0.0d0

        call defdr(oqc,-nclasp)
        call defdr(oqt,-nclasp)
        call defi (oinitc,-nclasp)
        call defi (ogroup,-2*nclasp)
Cchp1         call defdr(oves,-nclspp)
         allocate(rv_p_oves(abs(-nclspp)))
         if (-nclspp<0) rv_p_oves(:)=0.0d0

Cchp1         call defdr(ovdif,-nclspp)
         allocate(rv_p_ovdif(abs(-nclspp)))
         if (-nclspp<0) rv_p_ovdif(:)=0.0d0

        call defdr(orhrmx,-nclspp)
Cchp1         call defdr(ovrmax,-2*nclspp)
         allocate(rv_p_ovrmax(abs(-2*nclspp)))
         if (-2*nclspp<0) rv_p_ovrmax(:)=0.0d0

C       call defdr(oqnu,-3*nlspc)

Cchp1         oics=sarray%oics
         iv_p_oics => sarray%iv_p_oics 


c takao sp2cls is replaced.
c        call sp2cls('spec rmt',sspec,w(oics),1,1,nclspp,ormax)
        allocate( ics(nclspp) )
Cchp1         call icopy(nclspp,w(oics),1,ics,1)
         call icopy ( nclspp , iv_p_oics , 1 , ics , 1 ) 

        print *,' qqqqq ics=',ics
        call defdr(ormax,nclspp)
        call dcopy(nclspp,sspec(ics(1:nclspp))%rmt,1,w(ormax),1)


        sarray%ogroup=ogroup
        sarray%ormax=ormax
Cchp1         sarray%ohave=oinitc
         sarray%iv_p_ohave => iv_p_oinitc 



        spot%opnu=opnu
        spot%oqnu=oqnu
        spot%oqc=oqc
        spot%oqt=oqt


Cchp1         spot%oves=oves
         spot%rv_p_oves => rv_p_oves 

Cchp1         spot%ovdif=ovdif
         spot%rv_p_ovdif => rv_p_ovdif 

Cchp1         spot%ovrmax=ovrmax
         spot%rv_p_ovrmax => rv_p_ovrmax 

        spot%orhrmx=orhrmx
        spot%opp=opp

C       Even in collinear case, bxc might be written to
c        call defdr(obxc,-3*nclasp)

        spot%opprel=i_dummy !opprel
c         spot%obxc=obxc
c         spot%nrhos=0

c
c takao Kino's bug--> his replied that iand(-1,i)=i
c         if ( iand(-1,int(sctrl%lncol)) .ne.0 ) then

c         if ( sctrl%lncol /=0 ) then
c          neula=sham%neula
c          nbf=sham%nbf
c          nrhos = max(nl,neula,nbf)
c          call rlse(obxc)
c          call defdr(orhos,2*3*4*nrhos*nclasp)
c          call defdr(osop,-nlspc*nsp*9)
c          call defdr(obxc,-3*nclasp)
c          spot%nrhos=nrhos
c          spot%orhos=orhos
c          spot%obxc=obxc
c          spot%osop= osop
c        endif
        spot%osop=i_dummy !osop


cki        ogrrme = 0
          nullify(rv_p_ogrrme)
c         if ( int(sctrl%loptc) .gt. 0 ) call defrr ( ogrrme , -16 * nl
c     .   * nsp * nclasp )


Cchp1         spot%ogrrme=ogrrme
         spot%rv_p_ogrrme => rv_p_ogrrme 


c         if ( iand(3,int(sctrl%lsx)) .ne.0 .or. iand(-1,int(sctrl%lscr)) .ne.0
c     .   ) then
c         if ( iand(-1,int(sctrl%lscr)) .ne.0  ) then
        if ( int(sctrl%lscr) .ne.0  ) then
Cchp1           call defdr(ovintr,-nclasp*(nl*nsp)**2)
           allocate(rv_p_ovintr(abs(-nclasp*(nl*nsp)**2)))
           if (-nclasp*(nl*nsp)**2<0) rv_p_ovintr(:)=0.0d0

Cchp1           spot%ovintr=ovintr
           spot%rv_p_ovintr => rv_p_ovintr 


        endif
c         lnsph = isw ( iand(32,int(sctrl%lasa)) .ne.0 )
ctakao remove asa
        lnsph = isw ( .false.)

C   ... Allocate memory for, attempt to read qpp from file
        if (lnsph .ne. 0 .or. nlibu .ne. 0) then
          call defrr(opmpol,-(2*nl-1)*nl**2*3*nsp*nclasp)
          i1 = nl**2
          nqpp = (i1*(i1+1))/2
          call defrr(oqpp,2*nqpp*4*nsp*nbas)
          call dvset(w(oqpp),1,1,-1d0)
          spot%opmpol=opmpol
          spot%oqpp=oqpp
          call ioqpp(.false.,sctrl,spot)
        endif


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
c       if (allocated(pprel_rv)) deallocate(pprel_rv)
      return


      entry clsprp(mode,sctrl,sarray,sham,spot,sspec,slat,sbz,
     .sstr,sstrn)
C- Copy class data to doubly padded case


      nclass=sctrl%nclass
      nl=sctrl%nl
      nsp=sctrl%nspin
      nspec=sctrl%nspec


      nclasp=sarray%nclasp

c#error, have return with len(w_varlist)>0 at line 440
      if ( nclasp .eq. nclass ) then
        if (allocated(pprel_rv)) deallocate(pprel_rv)
        return
      endif

      nclspp = max(2*nclasp-nclass,nspec)
C ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then

        i1 = nclasp-nclass

        opnu=spot%opnu
        oqnu=spot%oqnu

        call dpscop(w(opnu),w(opnu),nl*nsp*i1,
     .  nl*nsp*nclass+1,nl*nsp*nclasp+1,1d0)
        call dpscop(w(oqnu),w(oqnu),3*nl*nsp*i1,
     .  3*nl*nsp*nclass+1,3*nl*nsp*nclasp+1,1d0)


Cchp1         oves=spot%oves
         rv_p_oves => spot%rv_p_oves 

Cchp1         ovdif=spot%ovdif
         rv_p_ovdif => spot%rv_p_ovdif 

Cchp1         ovrmax=spot%ovrmax
         rv_p_ovrmax => spot%rv_p_ovrmax 

        orhrmx=spot%orhrmx

Cchp1         call dpscop(w(oves),w(oves),i1,nclass+1,nclasp+1,1d0)
         call dpscop ( rv_p_oves , rv_p_oves , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

Cchp1         call dpscop(w(ovdif),w(ovdif),i1,nclass+1,nclasp+1,1d0)
         call dpscop ( rv_p_ovdif , rv_p_ovdif , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

Cchp1         call dpscop(w(ovrmax),w(ovrmax),i1,nclass+1,nclasp+1,1d0)
         call dpscop ( rv_p_ovrmax , rv_p_ovrmax , i1 , nclass + 1 , nclasp 
     .   + 1 , 1d0 ) 

        call dpscop(w(orhrmx),w(orhrmx),i1,nclass+1,nclasp+1,1d0)

      else
        call rxi('clsprp: bad mode',mode)
      endif

      if (allocated(pprel_rv)) deallocate(pprel_rv)

      end subroutine clsprm


C================================================================
c take subs/ioqpp.F here
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
      type(s_ctrl)::sctrl
      type(s_pot)::spot

C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
      integer i1,i2,nbas,nl,nsp,igets,ifi,oqpp,fopna,rdm,ipr

      call getpr(ipr)

      nbas=sctrl%nbas
      nl=sctrl%nl
      nsp=sctrl%nspin

      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
      oqpp = int(spot%oqpp)

      ifi = fopna('qpp',-1,4+8)

      if (lio) then
        call ywrm(1,'lmasa',1,ifi,' ',w(oqpp),1,i1,i1,i2)
        if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
        if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
        call pshpr(0)
        if (rdm(ifi,2,i1*i2,' ',w(oqpp),i1,i2) .lt. 0) then
          if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
          call dvset(w(oqpp),1,1,-1d0)
        endif
        call poppr
      endif

      call fclose(ifi)
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

      if (nsp .eq. 1) return

      if (mode .gt. 1) then
        call rx('spinav: bad mode')
      elseif (mode .eq. 1) then
        do  ic = 1, nclass
          call dswap(nl,pnu(1,2,ic),1,pnu(1,1,ic),1)
          call dswap(3*nl,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
        enddo
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
      end

