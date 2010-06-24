Cchp1 cmdvar= real(8),rv_p_,og
Cchp1 cmdvar= real(8),rv_p_,osymgr
Cchp1 cmdvar= real(8),rv_p_,odlv
Cchp1 cmdvar= real(8),rv_p_,oag
Cchp1 cmdvar= real(8),rv_p_,ormax
Cchp1 cmdvar= integer,iv_p_,owk
Cchp1 cmdvar= integer,iv_p_,onrc
Cchp1 cmdvar= integer,iv_p_,oips
Cchp1 cmdvar= real(8),rv_p_,orhrmx
Cchp1 cmdvar= real(8),rv_p_,oqt
Cchp1 cmdvar= real(8),rv_p_,oqnu
Cchp1 cmdvar= real(8),rv_p_,oqc
Cchp1 cmdvar= real(8),rv_p_,opp
Cchp1 cmdvar= real(8),rv_p_,opnu
Cchp1 cmdvar= integer,iv_p_,oinitc
Cchp1 cmdvar= integer,iv_p_,ohave
Cchp1 cmdvar= real(8),rv_p_,ovrmax
Cchp1 cmdvar= real(8),rv_p_,oves
Cchp1 cmdvar= real(8),rv_p_,oqlv
Cchp1 cmdvar= integer,iv_p_,oics
Cchp1 cmdvar= real(8),rv_p_,oclabl
      subroutine lmaux(prgnam,sctrl,sham,spot,slat,sspec,ssite,sstr,
     .sarray,sbz,slabl,sstrn,mode,wksize)

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
     .nbas,nclasp,nclass,nl,nlspc,nsp,modep(3),mxnbr,parg,nbasp,
     .nbaspp,nkd,nkq,nspec,neul,nc,mxcsiz,nttab,igets, !npadl,npadr,
     .iosits,cmplat,ngrp,ival,stdo,irs(5),fopna,nclspp,bitand,igetss,
     .ngmx,nsgrp
Cchp1       integer:: oeold , oeula , oinitc , oipc , olmx , onrc , opnu 
Cchp1      ., opold , opp , oqc , oqnu , oqt , orhrmx , ormax , oves , ovrmax 
Cchp1      ., owk , owk2 , oclabl , opos , oqlv , odlv , orham , oamsh , 
Cchp1 cki     .onrmsh , oalpha , oaamom , onpr , os , ormx , oip , oips , opgfsl 
Cchp1      .onrmsh , oalpha  , onpr , os , ormx , oip , oips , opgfsl 
Cchp1      ., ontab , oiax , og , oag , mxclas , oics
Cchp1        integer:: oeold , oeula , oinitc , oipc , olmx , onrc , opnu 
Cchp1      . , opold , opp , oqc , oqnu , oqt , orhrmx , ormax , oves , ovrmax 
Cchp1      . , owk , owk2 , opos , oqlv , odlv , orham , oamsh , onrmsh , 
Cchp1      . oalpha , onpr , os , ormx , oip , oips , opgfsl , ontab , oiax 
Cchp1      . , og , oag , mxclas , oics 
Cchp1         integer:: oeold , oeula , oinitc , oipc , olmx , onrc , opnu 
Cchp1      .  , opold , opp , oqc , oqnu , oqt , orhrmx , ormax , oves , ovrmax 
Cchp1      .  , owk , owk2 , opos , oqlv , odlv , orham , oamsh , onrmsh , 
Cchp1      .  oalpha , onpr , os , ormx , oip , oips , opgfsl , ontab , oiax 
Cchp1      .  , og , oag , mxclas 
Cchp1          integer:: oeold , oeula , oinitc , oipc , olmx , onrc , opnu 
Cchp1      .   , opold , opp , oqc , oqnu , oqt , orhrmx , ormax , oves , ovrmax 
Cchp1      .   , owk , owk2 , opos , odlv , orham , oamsh , onrmsh , oalpha 
Cchp1      .   , onpr , os , ormx , oip , oips , opgfsl , ontab , oiax , og 
Cchp1      .   , oag , mxclas 
Cchp1           integer:: oeold , oeula , oinitc , oipc , olmx , onrc , opnu 
Cchp1      .    , opold , opp , oqc , oqnu , oqt , orhrmx , ormax , ovrmax , 
Cchp1      .    owk , owk2 , opos , odlv , orham , oamsh , onrmsh , oalpha , 
Cchp1      .    onpr , os , ormx , oip , oips , opgfsl , ontab , oiax , og , 
Cchp1      .    oag , mxclas 
Cchp1            integer:: oeold , oeula , oinitc , oipc , olmx , onrc , opnu 
Cchp1      .     , opold , opp , oqc , oqnu , oqt , orhrmx , ormax , owk , owk2 
Cchp1      .     , opos , odlv , orham , oamsh , onrmsh , oalpha , onpr , os , 
Cchp1      .     ormx , oip , oips , opgfsl , ontab , oiax , og , oag , mxclas 
Cchp1             integer:: oeold , oeula , oipc , olmx , onrc , opnu , opold , 
Cchp1      .      opp , oqc , oqnu , oqt , orhrmx , ormax , owk , owk2 , opos , 
Cchp1      .      odlv , orham , oamsh , onrmsh , oalpha , onpr , os , ormx , oip 
Cchp1      .      , oips , opgfsl , ontab , oiax , og , oag , mxclas 
Cchp1              integer:: oeold , oeula , oipc , olmx , onrc , opold , opp , 
Cchp1      .       oqc , oqnu , oqt , orhrmx , ormax , owk , owk2 , opos , odlv 
Cchp1      .       , orham , oamsh , onrmsh , oalpha , onpr , os , ormx , oip , 
Cchp1      .       oips , opgfsl , ontab , oiax , og , oag , mxclas 
Cchp1               integer:: oeold , oeula , oipc , olmx , onrc , opold , oqc , 
Cchp1      .        oqnu , oqt , orhrmx , ormax , owk , owk2 , opos , odlv , orham 
Cchp1      .        , oamsh , onrmsh , oalpha , onpr , os , ormx , oip , oips , opgfsl 
Cchp1      .        , ontab , oiax , og , oag , mxclas 
Cchp1                integer:: oeold , oeula , oipc , olmx , onrc , opold , oqnu , 
Cchp1      .         oqt , orhrmx , ormax , owk , owk2 , opos , odlv , orham , oamsh 
Cchp1      .         , onrmsh , oalpha , onpr , os , ormx , oip , oips , opgfsl , 
Cchp1      .         ontab , oiax , og , oag , mxclas 
Cchp1                 integer:: oeold , oeula , oipc , olmx , onrc , opold , oqt , 
Cchp1      .          orhrmx , ormax , owk , owk2 , opos , odlv , orham , oamsh , onrmsh 
Cchp1      .          , oalpha , onpr , os , ormx , oip , oips , opgfsl , ontab , oiax 
Cchp1      .          , og , oag , mxclas 
Cchp1                  integer:: oeold , oeula , oipc , olmx , onrc , opold , orhrmx 
Cchp1      .           , ormax , owk , owk2 , opos , odlv , orham , oamsh , onrmsh , 
Cchp1      .           oalpha , onpr , os , ormx , oip , oips , opgfsl , ontab , oiax 
Cchp1      .           , og , oag , mxclas 
Cchp1                   integer:: oeold , oeula , oipc , olmx , onrc , opold , ormax 
Cchp1      .            , owk , owk2 , opos , odlv , orham , oamsh , onrmsh , oalpha 
Cchp1      .            , onpr , os , ormx , oip , oips , opgfsl , ontab , oiax , og 
Cchp1      .            , oag , mxclas 
Cchp1                    integer:: oeold , oeula , oipc , olmx , onrc , opold , ormax 
Cchp1      .             , owk , owk2 , opos , odlv , orham , oamsh , onrmsh , oalpha 
Cchp1      .             , onpr , os , ormx , oip , opgfsl , ontab , oiax , og , oag , 
Cchp1      .             mxclas 
Cchp1                     integer:: oeold , oeula , oipc , olmx , opold , ormax , owk , 
Cchp1      .              owk2 , opos , odlv , orham , oamsh , onrmsh , oalpha , onpr , 
Cchp1      .              os , ormx , oip , opgfsl , ontab , oiax , og , oag , mxclas 
Cchp1                      integer:: oeold , oeula , oipc , olmx , opold , ormax , owk2 
Cchp1      .               , opos , odlv , orham , oamsh , onrmsh , oalpha , onpr , os , 
Cchp1      .               ormx , oip , opgfsl , ontab , oiax , og , oag , mxclas 
Cchp1                       integer:: oeold , oeula , oipc , olmx , opold , owk2 , opos , 
Cchp1      .                odlv , orham , oamsh , onrmsh , oalpha , onpr , os , ormx , oip 
Cchp1      .                , opgfsl , ontab , oiax , og , oag , mxclas 
Cchp1                        integer:: oeold , oeula , oipc , olmx , opold , owk2 , opos , 
Cchp1      .                 odlv , orham , oamsh , onrmsh , oalpha , onpr , os , ormx , oip 
Cchp1      .                 , opgfsl , ontab , oiax , og , mxclas 
Cchp1                         integer:: oeold , oeula , oipc , olmx , opold , owk2 , opos , 
Cchp1      .                  orham , oamsh , onrmsh , oalpha , onpr , os , ormx , oip , opgfsl 
Cchp1      .                  , ontab , oiax , og , mxclas 
                         integer:: oeold , oeula , oipc , olmx , opold , owk2 , opos , 
     .                   orham , oamsh , onrmsh , oalpha , onpr , os , ormx , oip , opgfsl 
     .                   , ontab , oiax , mxclas 
                         real(8),pointer :: rv_p_og(:)

                        real(8),pointer :: rv_p_odlv(:)

                       real(8),pointer :: rv_p_oag(:)

                      real(8),pointer :: rv_p_ormax(:)

                     integer,pointer :: iv_p_owk(:)

                    integer,pointer :: iv_p_onrc(:)

                   integer,pointer :: iv_p_oips(:)

                  real(8),pointer :: rv_p_orhrmx(:)

                 real(8),pointer :: rv_p_oqt(:)

                real(8),pointer :: rv_p_oqnu(:)

               real(8),pointer :: rv_p_oqc(:)

              real(8),pointer :: rv_p_opp(:)

             real(8),pointer :: rv_p_opnu(:)

            integer,pointer :: iv_p_oinitc(:)

           real(8),pointer :: rv_p_ovrmax(:)

          real(8),pointer :: rv_p_oves(:)

         real(8),pointer :: rv_p_oqlv(:)

        integer,pointer :: iv_p_oics(:)

       real(8),pointer :: rv_p_oclabl(:)

      real(8) ,allocatable :: pos2_rv(:)
      real(8) ,allocatable :: rmt_rv(:)
      integer ,allocatable :: lock_iv(:)
      real(8) ,allocatable :: lockc_rv(:)
      real(8) ,allocatable :: z_rv(:)
      real(8) ,allocatable :: zz_rv(:)
      integer ,allocatable :: ips2_iv(:)
      real(8) ,allocatable :: zc_rv(:)
      real(8) ,allocatable :: rmtc_rv(:)

      double precision xv(10),xx,alat,plat(3,3),facrmx,facrng, !,plat2(9)
     .dval,avw,ekap(2),enu,rmaxs,qss(4),ckbas,cksumf,ehterm(4),
     .qlat(9),emad,trumad,vmtz(2),omax1(3),omax2(3),wsrmax
      parameter (ngmx=48,mxclas=1000)

      integer:: i_copy_size, i_spackv, i_spacks
C ... Heap
      integer w(1)
      common /w/ w

ctakao
      integer:: ifx,w_dummy=1
      integer,allocatable:: ics(:)
      integer,allocatable:: lmxa(:)
c      real(8),allocatable::idmod(:,:), z(:),amsh(:),rmax(:), nrmsh(:)
      real(8),allocatable:: z(:),rmax(:)

      print *,' lmaux:'


      nbas=sctrl%nbas
      nclass=sctrl%nclass
      nl=sctrl%nl
      nspec=sctrl%nspec
      nsp=sctrl%nspin


c       lncol=sctrl%lncol
      i_copy_size=size(sctrl%modep)
      call i8icopy(i_copy_size,sctrl%modep,1,modep,1)

c       lpgf = int(sctrl%lpgf(1))

C     lpbc = 0 for pbc in 3 dimensions, 11 pgf padded geometry
      lpbc = 0
c      if (lpgf .gt. 0 .and. lpgf .lt. 10) lpbc = 11

Cchp1       oinitc=sarray%ohave
       iv_p_oinitc => sarray%iv_p_ohave 

Cchp1       oics=sarray%oics
       iv_p_oics => sarray%iv_p_oics 

Cchp1       oclabl=sarray%oclabl
       rv_p_oclabl => sarray%rv_p_oclabl 

Cchp1       onrc=sarray%onrc
       iv_p_onrc => sarray%iv_p_onrc 

ccc


c       npadl=sarray%npadl
c       npadr=sarray%npadr
      nclasp=sarray%nclasp
c       opgfsl=sarray%opgfsl
Cchp1       ormax=sarray%ormax
       rv_p_ormax => sarray%rv_p_ormax 



      opos=sarray%opos
      oipc=sarray%oipc
Cchp1       oips=sarray%oips
       iv_p_oips => sarray%iv_p_oips 



c       oeula=sham%oeula
c       neul=sham%neula
c       i_copy_size=size(sham%qss)
c       call dcopy(i_copy_size,sham%qss,1,qss,1)


Cchp1       opnu=spot%opnu
       rv_p_opnu => spot%rv_p_opnu 

Cchp1       oqnu=spot%oqnu
       rv_p_oqnu => spot%rv_p_oqnu 

Cchp1       opp=spot%opp
       rv_p_opp => spot%rv_p_opp 

Cchp1       oqc=spot%oqc
       rv_p_oqc => spot%rv_p_oqc 

Cchp1       oqt=spot%oqt
       rv_p_oqt => spot%rv_p_oqt 



cki      oaamom=spot%oaamom
Cchp1       orhrmx=spot%orhrmx
       rv_p_orhrmx => spot%rv_p_orhrmx 

Cchp1       oves=spot%oves
       rv_p_oves => spot%rv_p_oves 

Cchp1       ovrmax=spot%ovrmax
       rv_p_ovrmax => spot%rv_p_ovrmax 



      avw=slat%avw
      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
c       i_copy_size=size(slat%plat2)
c       call dcopy(i_copy_size,slat%plat2,1,plat2,1)


      nkd=slat%nkd
Cchp1       odlv=slat%odlv
       rv_p_odlv => slat%rv_p_odlv 

      nkq=slat%nkq
Cchp1       oqlv=slat%oqlv
       rv_p_oqlv => slat%rv_p_oqlv 



      mxnbr=sstr%mxnbr
      rmaxs=sstr%rmax

ctakao
c      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmx) --->lmxa
c      call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)
      allocate(ics(1:nclasp),lmxa(nclasp),z(nclasp))
Cchp1       call icopy(nclasp,w(oics),1,ics,1)
       call icopy ( nclasp , iv_p_oics , 1 , ics , 1 ) 

      lmxa(1:nclasp) = sspec(ics(1:nclasp))%lmxa
      z   (1:nclasp) = sspec(ics(1:nclasp))%z
      print *,' lmxa=',lmxa
      print *,' z   =',z

      nbasp = nbas !+ npadl + npadr
      nbaspp = 2*nbasp - nbas

      stdo = lgunit(1)

      j = 10
      if (cmdopt('--shorten',j-1,0,outs)) then
        call shorps(nbasp,plat,modep,w(opos),w(opos))
      endif

C ... Read from restart file
      if (cmdopt('--rs=',5,0,outs)) then
        irs ( 1 ) = iand(7,int(sctrl%lrs))

        if (irs(1) .gt. 0) then
          ifi = fopna('rsta',-1,0)
Cchp1           call asars(0,sctrl,ssite,sspec,slat,spot,sbz,
Cchp1      .    w(opnu),w(oqnu),.false.,ifi)
Cchp1            call asars ( 0 , sctrl , ssite , sspec , slat , spot , sbz , 
Cchp1      .     rv_p_opnu , w ( oqnu ) , .false. , ifi ) 
            call asars ( 0 , sctrl , ssite , sspec , slat , spot , sbz , 
     .      rv_p_opnu , rv_p_oqnu , .false. , ifi ) 


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
        mxcsiz = int(sstr%mxnbr)

        call pshpr(iprint()-20)
c      call pairs(nbas,nbasp,alat,plat,rmaxs/2,w(opos),
c     .  -1,i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
        call pairs(nbas,nbasp,alat,plat,rmaxs/2,w(opos),
     .  -1,i,j,w_dummy,nttab,ontab,oiax,mxcsiz)
        call poppr

C --- Print out a few superlattice vectors ---
        j = 6
        if (cmdopt('--slat',j-1,0,outs)) then
          if (iprint() .ge. 10) then
            call info0(10,1,0,' LMCHK:  print multiples of plat%N'//
     .      '  i1  i2  i3%7fx%11fy%11fz%11flen')
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
  368             format(3i4, 3f12.7, 1x, f12.5)
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
          allocate(zz_rv(nspec))
Ctakao_ZeroClear_NotRequiered          if (nspec<0) zz_rv(:)=0.0d0

          i_copy_size=1;
          do i_spackv=1,nspec
            call spackv_array_copy_r8_r8 ( 'u' , sspec ( i_spackv ) %z , 
     .      i_copy_size , i_spackv + 1 - 1 , zz_rv )

          enddo

          allocate(rmt_rv(nspec))
Ctakao_ZeroClear_NotRequiered          if (nspec<0) rmt_rv(:)=0.0d0

          i_copy_size=1;
          do i_spackv=1,nspec
            call spackv_array_copy_r8_r8 ( 'u' , sspec ( i_spackv ) %rmt 
     .      , i_copy_size , i_spackv + 1 - 1 , rmt_rv )

          enddo

          allocate(lock_iv(abs(-nspec)))
          if (-nspec<0) lock_iv(:)=0

          do  i = 1, nspec
            call ivset ( lock_iv , i , i , bitand ( int ( sspec ( i ) 
     .      %mxcst ) , 2 ) )


          enddo
          if (lpbc .eq. 0) then
            i = 3
          elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
            i = 2
          else
            call rx('LMAUX: not implemented for lpbc>1')
          endif
Cchp1           call makrm0 ( 101 , nspec , nbas , alat , plat , w ( opos ) , 
Cchp1      .    slabl , w ( oips ) , modep , lock_iv , zz_rv , rmt_rv 
Cchp1      .    )
           call makrm0 ( 101 , nspec , nbas , alat , plat , w ( opos ) , 
     .     slabl , iv_p_oips , modep , lock_iv , zz_rv , rmt_rv ) 


C   ... Scale sphere radii satisfying constraints

          i_copy_size=size(sctrl%omax1)
          call dcopy(i_copy_size,sctrl%omax1,1,omax1,1)
          i_copy_size=size(sctrl%omax2)
          call dcopy(i_copy_size,sctrl%omax2,1,omax2,1)
          wsrmax=sctrl%wsrmax

Cchp1           call sclwsr ( 20 , nbas , nbasp , nspec , alat , plat , w ( opos 
Cchp1      .    ) , w ( oips ) , modep , slabl , zz_rv , lock_iv , 1d0 
Cchp1      .    , wsrmax , omax1 , omax2 , rmt_rv )
           call sclwsr ( 20 , nbas , nbasp , nspec , alat , plat , w ( opos 
     .     ) , iv_p_oips , modep , slabl , zz_rv , lock_iv , 1d0 , wsrmax 
     .     , omax1 , omax2 , rmt_rv ) 


          i_copy_size=1;
          do i_spackv=1,nspec
            call spackv_array_copy_r8_r8 ( 'p' , sspec ( i_spackv ) %rmt 
     .      , i_copy_size , i_spackv + 1 - 1 , rmt_rv )

          enddo

          nclspp = max(2*nclasp-nclass,nspec)
c takao
c        call sp2cls('spec rmt',sspec,ics,1,1,nclspp,ormax)
          allocate(rmax(nclspp))
          print *,' zzzz nclspp=',nclspp
          do i=1,nclspp
            rmax(i) =sspec(ics(i))%rmt
          enddo
Cchp1           call defdr(ormax,nclspp)
           allocate(rv_p_ormax(abs(nclspp)))
           if (nclspp<0) rv_p_ormax(:)=0.0d0

Cchp1           call dcopy(nclspp,rmax,1,w(ormax),1)
           call dcopy ( nclspp , rmax , 1 , rv_p_ormax , 1 ) 


Cchp1           sarray%ormax=ormax
           sarray%rv_p_ormax => rv_p_ormax 

          ifx=1
        endif

c-------
        if(ifx==0) then
          allocate(rmax(nclasp))
Cchp1           call dcopy(nclasp,w(ormax),1,rmax,1)
           call dcopy ( nclasp , rv_p_ormax , 1 , rmax , 1 ) 

        endif

c


C --- Show neighbors by shell ---
        outs = ' '
        j = 8
        if (cmdopt('--shell',j-1,0,outs)) then
Cchp1           call shoshl(outs(j:),nbas,w(opos),plat,mxnbr,z,slabl,
Cchp1      .    w(oclabl),w(oips),w(oipc),w(oves),w(oeula),nclass)
Cchp1            call shoshl ( outs ( j: ) , nbas , w ( opos ) , plat , mxnbr 
Cchp1      .     , z , slabl , rv_p_oclabl , w ( oips ) , w ( oipc ) , w ( oves 
Cchp1      .     ) , w ( oeula ) , nclass ) 
Cchp1             call shoshl ( outs ( j: ) , nbas , w ( opos ) , plat , mxnbr 
Cchp1      .      , z , slabl , rv_p_oclabl , w ( oips ) , w ( oipc ) , rv_p_oves 
Cchp1      .      , w ( oeula ) , nclass ) 
             call shoshl ( outs ( j: ) , nbas , w ( opos ) , plat , mxnbr 
     .       , z , slabl , rv_p_oclabl , iv_p_oips , w ( oipc ) , rv_p_oves 
     .       , w ( oeula ) , nclass ) 



        endif

C --- Show angles between neighbors ---
        j = 9
        if (cmdopt('--angles',j-1,0,outs)) then
Cchp1           call shoang(outs(j:),nbas,w(opos),plat,mxnbr,slabl,w(oips))
           call shoang ( outs ( j: ) , nbas , w ( opos ) , plat , mxnbr 
     .     , slabl , iv_p_oips ) 

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
          allocate(pos2_rv(3*nbas))
Ctakao_ZeroClear_NotRequiered          if (3*nbas<0) pos2_rv(:)=0.0d0

          allocate(ips2_iv(nbas))
Ctakao_ZeroClear_NotRequiered          if (nbas<0) ips2_iv(:)=0

          i_copy_size=size(ssite(1)%pos)
          do i_spackv=1,nbas
            call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .      , i_copy_size , i_spackv + 1 - 1 , pos2_rv )

          enddo

          i_copy_size=1;
          do i_spackv=1,nbas
            call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %spec 
     .      , i_copy_size , i_spackv + 1 - 1 , ips2_iv )

          enddo


Cchp1           call defrr(og,   9*ngmx)
           allocate(rv_p_og(abs(9*ngmx)))
           if (9*ngmx<0) rv_p_og(:)=0.0d0

Cchp1           call symlat(plat,ngrp,w(og),j)
           call symlat ( plat , ngrp , rv_p_og , j ) 

Cchp1           j = cmplat ( nbas , plat , plat , ngrp , w ( og ) , w ( oips 
Cchp1      .    ) , w ( opos ) , ips2_iv , pos2_rv )
Cchp1            j = cmplat ( nbas , plat , plat , ngrp , w ( og ) , iv_p_oips 
Cchp1      .     , w ( opos ) , ips2_iv , pos2_rv ) 
            j = cmplat ( nbas , plat , plat , ngrp , rv_p_og , iv_p_oips 
     .      , w ( opos ) , ips2_iv , pos2_rv ) 



          call fexit(j,1,' Exit %i lmchk --basis= ...',j)
        endif

C ... Write positions in Cartesian coordinates and as multiples plat
        if (iprint() .ge. 50) then
          write(stdo,357)
  357     format(/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .    'pos (multiples of plat)')
C     qlat = (plat+)^-1
          call dinv33(plat,1,qlat,xx)
          do  i = 1, nbas
            call dpscop(w(opos),xv,3,3*i-2,1,1d0)
C       posp+ = (plat)^-1 pos+
            call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
Cchp1             ip = ival(w(oips),i)
             ip = ival ( iv_p_oips , i ) 

            print 345, i, slabl(ip), (xv(j),j=1,3), (xv(3+j),j=1,3)
  345       format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
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
Cchp1         call ovmin(outs(i:),nbas,nbasp,alat,plat,rmax,rmax,
Cchp1      .  w(oclabl),w(oipc),modep,z,ontab,oiax,w(opos),j)
         call ovmin ( outs ( i: ) , nbas , nbasp , alat , plat , rmax 
     .   , rmax , rv_p_oclabl , w ( oipc ) , modep , z , ontab , oiax 
     .   , w ( opos ) , j ) 

C ... Write positions to file
        if (cmdopt('--wpos=',7,0,outs))
     .  call iopos(T,0,outs(8:),nbasp,w(opos))

C ... Inner products between Euler angles
c      if (mod(lncol,2) .eq. 1) then
c        call amagnc(nbas,nl,w(oipc),w(1),1,w(oqnu),w(oeula),neul,1,xv,
c     .    w(oaamom))
c        print '(1x)'
c        call angtab(nbas,w(opos),alat,rmax,qss,w(oqlv),
c     .    w(odlv),nkd,w(oipc),neul,w(oeula))
c      endif

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

        nsgrp=slat%nsgrp
Cchp1         og=slat%osymgr
         rv_p_og => slat%rv_p_osymgr 

Cchp1         oag=slat%oag
         rv_p_oag => slat%rv_p_oag 


        allocate(z_rv(nspec))
Ctakao_ZeroClear_NotRequiered          if (nspec<0) z_rv(:)=0.0d0

        i_copy_size=1;
        do i_spackv=1,nspec
          call spackv_array_copy_r8_r8('u',sspec(i_spackv)%z,i_copy_size,i_spackv+1-1,z)
        enddo

        allocate(rmt_rv(nspec))
Ctakao_ZeroClear_NotRequiered          if (nspec<0) rmt_rv(:)=0.0d0

        i_copy_size=1;
        do i_spackv=1,nspec
          call spackv_array_copy_r8_r8 ( 'u' , sspec ( i_spackv ) %rmt 
     .    , i_copy_size , i_spackv + 1 - 1 , rmt_rv )

        enddo

        allocate(lock_iv(abs(-nspec)))
        if (-nspec<0) lock_iv(:)=0

        do  i = 1, nspec
          call ivset ( lock_iv , i , i , bitand ( int ( sspec ( i ) 
     .    %mxcst ) , 2 ) )


        enddo
        allocate(zc_rv(mxclas))
Ctakao_ZeroClear_NotRequiered          if (mxclas<0) zc_rv(:)=0.0d0

        allocate(rmtc_rv(mxclas))
Ctakao_ZeroClear_NotRequiered          if (mxclas<0) rmtc_rv(:)=0.0d0

        allocate(lockc_rv(mxclas))
Ctakao_ZeroClear_NotRequiered          if (mxclas<0) lockc_rv(:)=0.0d0

cki        owk = onrc
        iv_p_owk => iv_p_onrc
Cchp1         call defi(onrc,mxclas)
         allocate(iv_p_onrc(abs(mxclas)))
         if (mxclas<0) iv_p_onrc(:)=0

Cchp1         call icopy(nclass,w(owk),1,w(onrc),1)
Cchp1          call icopy ( nclass , w ( owk ) , 1 , iv_p_onrc , 1 ) 
          call icopy ( nclass , iv_p_owk , 1 , iv_p_onrc , 1 ) 



        if (nclass .ge. mxclas) call rx('lmaux: increase mxclas')
        call spec2c ( nspec , nclass , ics , rmt_rv , rmtc_rv , 
     .  z , zc_rv , lock_iv , lockc_rv )

Cchp1         call findes ( sctrl , alat , w ( oclabl ) , nbas , nclass , nl 
Cchp1      .  , w ( onrc ) , mxclas , nsgrp , opos , oipc , plat , w ( og ) 
Cchp1      .  , w ( oag ) , lockc_rv , rmtc_rv , zc_rv )
Cchp1          call findes ( sctrl , alat , rv_p_oclabl , nbas , nclass , nl 
Cchp1      .   , w ( onrc ) , mxclas , nsgrp , opos , oipc , plat , w ( og ) 
Cchp1      .   , w ( oag ) , lockc_rv , rmtc_rv , zc_rv ) 
Cchp1           call findes ( sctrl , alat , rv_p_oclabl , nbas , nclass , nl 
Cchp1      .    , iv_p_onrc , mxclas , nsgrp , opos , oipc , plat , w ( og ) 
Cchp1      .    , w ( oag ) , lockc_rv , rmtc_rv , zc_rv ) 
Cchp1            call findes ( sctrl , alat , rv_p_oclabl , nbas , nclass , nl 
Cchp1      .     , iv_p_onrc , mxclas , nsgrp , opos , oipc , plat , w ( og ) 
Cchp1      .     , rv_p_oag , lockc_rv , rmtc_rv , zc_rv ) 
            call findes ( sctrl , alat , rv_p_oclabl , nbas , nclass , nl 
     .      , iv_p_onrc , mxclas , nsgrp , opos , oipc , plat , rv_p_og , 
     .      rv_p_oag , lockc_rv , rmtc_rv , zc_rv ) 






      endif
ctakao
      deallocate(ics,lmxa,z)
      if (allocated(lockc_rv)) deallocate(lockc_rv)
      if (allocated(rmtc_rv)) deallocate(rmtc_rv)
      if (allocated(zc_rv)) deallocate(zc_rv)
      if (allocated(z_rv)) deallocate(z_rv)
      if (allocated(ips2_iv)) deallocate(ips2_iv)
      if (allocated(pos2_rv)) deallocate(pos2_rv)
      if (allocated(lock_iv)) deallocate(lock_iv)
      if (allocated(rmt_rv)) deallocate(rmt_rv)
      if (allocated(zz_rv)) deallocate(zz_rv)

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

      end subroutine spec2c



