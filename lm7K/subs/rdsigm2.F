#define USEOPTIONAL 1
      subroutine seneinterp(nbas,nsp,ndimh,slat,sham,sbz,ifis,rsrnge,
     &nk1,nk2,nk3,sfz) !,qin,isp,sene)

      use m_struc_def  !Cgetarg
      use m_hft2rs 
      use m_globalvariables
      use m_sigflg
C- takao. obtain sene at given q (qin).
C
Ci sfz
Co sene
C ------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,ifis,ndimh,lwsig
C     integer osig(nbas),otau(nbas),oppi(nbas)
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_bz)::sbz

C     double precision ssite(1),sspec(1)
      double precision rsrnge
C ... Local parameters
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
     .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .ndhrs,j1,k1,k2,k3,iq1,nl,nspc,mxorb,nqsig,modsgp,nmin,nmax,
     .fopna,kcplx,lrot,iprint,stdo,lssym,ledit,niax,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      parameter (niax=10)
c      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
                 integer:: ohrss 
ckino Dec.14.2011:                   integer,pointer :: iv_p_ontabs(:) =>NULL()
ckino 12月.27.2011:                   integer,allocatable :: iv_a_ontabs(:) 

ckino Dec.14.2011:                  integer,pointer :: iv_p_oiaxs(:) =>NULL()
ckino Dec.27.2011:                  integer,allocatable :: iv_a_oiaxs(:) 

ckino Dec.22.2011:                 real(8),pointer :: rv_p_opos(:) =>NULL()

ckino Dec.22.2011:                real(8),pointer :: rv_p_og(:) =>NULL()

ckino Dec.22.2011:               real(8),pointer :: rv_p_oag(:) =>NULL()

ckino Dec.22.2011:              real(8),pointer :: rv_p_ohrs(:) =>NULL()

ckino Dec.22.2011:           integer,pointer :: iv_p_oistab(:) =>NULL()

ckino Dec.22.2011:            integer,pointer :: iv_p_oiprmb(:) =>NULL()

ckino Dec.22.2011:           real(8),pointer :: rv_p_oqsig(:) =>NULL()

        real(8),pointer :: rv_p_oqp(:) =>NULL()

ckino Dec.22.2011:         integer,pointer :: iv_p_ooffh(:) =>NULL()

      integer ,allocatable :: gstar_iv(:)
      integer ,allocatable :: ipq_iv(:)
      real(8) ,allocatable :: wgt_rv(:)
      complex(8) ,allocatable :: wk_zv(:)
      complex(8) ,allocatable :: sigm_zv(:)
      complex(8) ,allocatable :: sigm2_zv(:)
      real(8) ,allocatable :: delt_rv(:)
      integer ,allocatable :: istb2_iv(:)

      double precision tmp(3)
      real(8),allocatable:: evls(:),evlz(:),sigii(:)
      complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
C     For offset q mesh
      integer is(3),lshft(3),ifac(3),lqoffo
      double precision rb(3,3),qb(3,3),qoffi(3),qoffo(3)
C MPI
      logical mlog
      integer procid,mpipid,master
C ... for mixing sigma files: Use alf1*sigm + alf2*sigm1
      double precision alf(2)
C     integer os
      double precision plat(3,3),qp(3),tolq,rsstol,sigp(10)  ,qlat(3,3)
      double precision emin,emax,asig,bsig,efit,qoff(3),rotm(3,3),ddot
      equivalence (emin,sigp(3)),(emax,sigp(5))
      equivalence (asig,sigp(6)),(bsig,sigp(7)),(efit,sigp(8))
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      double precision qk
      integer jj1,jj2,jj3,k
C.....titus
C      integer ib,ib1,ib2,nlmaa,lidim,ldh,osfz1,nsp1
C      integer nkap0,n0H
C      parameter  (nkap0=3,n0H=5)
C      integer offH(n0H,nkap0,nbas)
C.....titus
C ... Heap
       integer ::iwdummy 


      real(8):: eseavr
      integer:: napw_in
      complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp)


      real(8):: qin(3)
ckino delete integer(4) def.      integer(4):: ndhamx,debugmode
      integer:: ndhamx,debugmode
      complex(8):: sene(ndimh,ndimh)

      complex(8):: zv_dummy(1)
      real(8):: rv_dummy(1)

      integer:: i_copy_size,nhrss
      real(8),allocatable::hrss(:)

      logical:: l_dummy_isanrg, isanrg
ccccccccccccccccccccccccccccccc
c      integer:: oss
ccccccccccccccccccccccccccccccc
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .(jj2*ifac(2)-1)*qb(k,2) +
     .(jj3*ifac(3)-1)*qb(k,3)


c
      print *
      print *,'seneinterp:'
      mode=12 !lrsig=12
      lwsig=0



C --- Setup and printout ---
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

C     tolq is allowed fuzz in qp
      tolq = 1d-6
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      ledit = 0
      lrsig = mod(mode,10)
      lwsig = 0
      moditp = mod(mode/10,10)
      alf(1) = 1
      alf(2) = 0
      nscnts = 0
      lqoffo = 0
      lfbzout = .false.

C     Switch flagging whether input file has symops or not
      hreal = 0
      if (lrsig .ge. 4) hreal = 1
      if (lrsig .eq. 0) return

      ltrans = mod(mode/1000,10) .eq. 0
      kcplx = 1
      if (mod(mode/1000,10) .eq. 1) kcplx = 0
      lphase = mod(mode/1000,10) .eq. 0
      lssym = mod(mode/10000,10)
Ckino isanrg is logical function,       call isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)
      l_dummy_isanrg=isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)

C ... Command-line switches for input sigma and defaults
#ifdef USEOPTIONAL
       call sigflg (dc= ' ' , sio='rsig' , strn=' ' 
     . , lfbz=lfbzin  , lonesp=lonesp , qpoff=qoffi  ) 
#else
       call sigflg ( ' ' , 'rsig' , ' ' , iwdummy , iwdummy , iwdummy 
     . , lfbzin , iwdummy , lonesp , qoffi , iwdummy ) 
#endif

      if (cmdopt('--rsig',6,0,outs)) then
        out2 = outs(7:)
        dc = out2(1:1)
#ifdef USEOPTIONAL
         call sigflg ( dc=dc , sio='rsig' , strn=out2 
     .   , lfbz=lfbzin  ,lonesp= lonesp , qpoff=qoffi  ) 
#else
         call sigflg ( dc , 'rsig' , out2 , iwdummy , iwdummy , iwdummy 
     .   , lfbzin , iwdummy , lonesp , qoffi , iwdummy ) 
#endif

      endif
#if 1
      if ( lfbzin .or. mod(lssym,2).eq.1 ) then
        lfbzin = .true.
      else
        lfbzin = .false.
      endif
#else
      lfbzin = lfbzin .or. mod(lssym,2)
#endif

      if (cmdopt('--mixsig=',9,0,outs)) then
        j = 0
        i = parg('--mixsig=',4,outs,j,len(outs),
     .  ', ',2,2,ifac,alf)
        if (i .ne. 1 .and. i .ne. 2) call rx
     .  ('rdsigm: failed to parse arguments to --mixsig=')
      endif

      call info(10,1,0,' RDSIGM: read file sigm and create '//
     .'%?#n#REAL#COMPLEX# sigma(R) by FT ...',hreal,0)

      i_copy_size=size(sham%ldham)
      call i8icopy(i_copy_size,sham%ldham,1,ldham,1)
ckino Dec.22.2011:         iv_p_ooffh => sham%iv_p_ooffh 

ckino Dec.22.2011:         iv_p_oiprmb => sham%iv_p_oindxo 

      i_copy_size=size(sham%sigp)
      call dcopy(i_copy_size,sham%sigp,1,sigp,1)
      rsstol=sham%rsstol

      modsgp = nint(sigp(1))
      nmin   = nint(sigp(2))
      nmax   = nint(sigp(4))
      call info2(30,0,0,
     .'%9fSigm will be approximated by:  '//
     .'%?#(n==0)#Simple bloch sum##%-1j'//
     .'%?#(n==1)#diagonal Sigma for high and low states##%-1j'//
     .'%?#(n==2)#Perturbation for higher kappa blocks##%-1j'//
     .'%?#(n==3)#Interpolation from known points##%-1j'//
     .' ',moditp,0)
      if (moditp .eq. 1 .or. moditp .eq. 2) then
        call info5(30,0,0,
     .  '%9fApproximate sigma '//
     .  '%?#(n<0)#for energies E(lda)<%d; and %-2j#%-1j#'//
     .  '%?#(n>0)#for states up to %-1jn=%i; and ##%j'//
     .  '%?#(n<=0)#for energies E(lda)>%d%-2j#%-1j#'//
     .  '%?#(n>0)#for states above %-1jn=%i##%j',
     .  nmin,emin,nmax,emax,0)
        call info5(30,0,0,'%9fFor high states '//
     .  '%?#(n==0)#Sigii > %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
     .  '%?#(n==1)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
     .  '%?#(n==3)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
     .  '%?#(n==2)#%,3;4d < Sigii < %,3;4d%-2j##%-1j'//
     .  ' ',modsgp,asig,bsig,0,0)
        if (nmax .gt. 0) call info5(30,0,0,
     .  '%9fLinear fit to diagonal sigma for E(lda)>%d'//
     .  '%?#(n==3)# and n<%i'//
     .  ' ',efit,modsgp,nmax,0,0)
        if (nmax .eq. 0) call info5(30,0,0,
     .  '%9fLinear fit to diagonal sigma for E(lda)>%d'//
     .  '%?#(n==3)# and E(lda)<%d'//
     .  ' ',efit,modsgp,emax,0,0)
Ckino isanrg is logical function,         call isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
        l_dummy_isanrg=isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
      endif
      if (lfbzin .or. ddot(3,qoffi,1,qoffi,1).ne.0) then
        call info5(30,0,0,
     .  '%?#(n==1)#%9fFile sigm in saved in FBZ: '//
     .  'symmetrization suppressed##'//
     .  '%?#(n==1)#%N%9fFile k-mesh is offset by:%3:2g#%j#'//
     .  ' ',isw(lfbzin),isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi,0,0)
      endif

      if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
        call info5(30,0,0,
     .  '%9fUse for sigma: %d*(file sigm)%?#n# + '//
     .  '%d*(file sigm1)',alf(1),isw(alf(2).ne.0),alf(2),0,0)
      endif

      call tcn('seneinterp')
C --- Read sigma(orbital basis) from file ---

C ... Read k-mesh parameters
      if (procid .eq. master) then
C       Require file contains sigm(orbital basis)
        call iosigh(0,nscnts,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
     .  lshft(3),ifis)
        if (nscnts .ne. 0 .and. nscnts .ne. 5) then
          call iosigh(2,0,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
     .    lshft(3),ifis)
        endif
      endif
      call mpibc1(i,1,2,.false.,'rdsigm','nsp')
      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      if (i .eq. nsp) then
        lsplts = .false.
      elseif (i .gt. nsp) then
        call rx(
     .  'rdsigm: sigm file spin polarized but calculation is not')
      else
        lsplts = .true.
      endif

C     Check for consistency in remaining ndimh, nqp
C     and set file sigma file pointer past header
      if (procid .eq. master) then
        rewind ifis
        if (lfbzin) nqp = nk1*nk2*nk3
        call iosigh(3,nscnts,i,ndimh,nk1,nk2,nk3,nqp,
     .  lshft(1),lshft(2),lshft(3),ifis)
      endif
      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      call info5(30,0,0,'%9fsigm file has %i irreducible QP: nk ='//
     .' ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)
C     Sanity checks
      if (lsplts) call info0(30,0,0,
     .'%9f(warning) sigm file not spin pol .. splitting spins')

      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
ckino Dec.22.2011:         rv_p_opos => slat%rv_p_opos 



      nsgrp=slat%nsgrp
ckino Dec.22.2011:         iv_p_oistab => slat%iv_p_oistab 

ckino Dec.22.2011:         rv_p_og => slat%rv_p_osymgr 

ckino Dec.22.2011:         rv_p_oag => slat%rv_p_oag 


C     Get leading dimension of hrs
ckino Dec.22.2011:         call offshp ( iv_p_oiprmb , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs 
       call offshp ( sham%iv_p_oindxo , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs 
     . ) 

C     Setup for FFT
      call fftz30(nk1,nk2,nk3,k1,k2,k3) !takao comment out because k1 k2 k3 are supplied at the begining of routine.
      if (nk1.ne.k1 .or. nk2.ne.k2 .or. nk3.ne.k3) call rx(
     .'rdsigm: not ready for FFT w/ dimensions ne no div.')
C ... Number of group operations for input file sigma
      if (mod(lssym,2) .eq. 0) nsgrps = nsgrp
      if (mod(lssym,2) .eq. 1) nsgrps = 1
      if (lfbzin) nsgrps = 1
C ... Make is,ifac,qb,qlat,qoff
      do  8  i = 1, 3
        llshft(i) = lshft(i) .ne. 0
    8 continue
      call pshpr(0)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
      call poppr
      do  i = 1, 3
        qoff(i) = qk(i,1,1,1) + qoffi(i)
      enddo

C ... Setup for r.s. hamiltonian: allocate and create iax,ntab; allocate hrs
      i = 1000 + 10*hreal + 0
C      print *, '!!'; i = 1000*0 + 10*hreal + 0
      if (mod(lssym,4) .ge. 2) i = 10*hreal + 0
      if (ltrans) i = i + 10000
      print *,' go to hft2rs i=',i
          call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , 1 , 
ckino Dec.22.2011:       .    nsp , nbas , rv_p_og , rv_p_oag , nsgrp , rsrnge , iv_p_ooffh 
ckino Dec.22.2011:       .    nsp , nbas , slat%rv_p_osymgr , rv_p_oag , nsgrp , rsrnge , iv_p_ooffh 
ckino Dec.22.2011:       .    nsp , nbas , slat%rv_p_osymgr , slat%rv_p_oag , nsgrp , rsrnge , iv_p_ooffh 
     .    nsp , nbas , slat%rv_p_osymgr , slat%rv_p_oag , nsgrp , rsrnge , sham%iv_p_ooffh 
cki     .    , 1 , nbas , ndimh , ndimh , w , plat , rv_p_opos , iv_p_ontabs , 
ckino Dec.14.2011:       .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , rv_p_opos , iv_p_ontabs , 
ckino Dec.22.2011:       .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , rv_p_opos , iv_a_ontabs , 
ckino 12月.27.2011:       .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , slat%rv_p_opos , sham%iv_a_ontabs , 
     .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , slat%rv_p_opos , sham%iv_a_ontabs , 
cki     .    iv_p_oiaxs , ndhrs , w ) 
ckino Dec.14.2011:       .    iv_p_oiaxs , ndhrs , rv_dummy ) 
ckino Dec.27.2011:       .    sham%iv_a_oiaxs , ndhrs , rv_dummy ) 
     .    sham%iv_a_oiaxs , ndhrs , rv_dummy ) 

ccccccccccccccccccccccccccccccccc
c          print *,' size iv_p_ontabls=',size(iv_p_ontabs)
cccccccccccccccccccccccccccccccccccccc


ckino Dec.14.2011:         nttabs = ival ( iv_p_ontabs , nbas + 1 ) 
ckino 12月.27.2011:         nttabs = ival ( sham%iv_a_ontabs , nbas + 1 ) 
       nttabs = ival ( sham%iv_a_ontabs , nbas + 1 ) 

C     Allocate memory for hrs
      if (hreal .eq. 1) then
ckino Dec.22.2011:           allocate(rv_p_ohrs(ndhrs**2*nttabs*nsp))
         allocate(sham%rv_p_ohrs(ndhrs**2*nttabs*nsp))
ckino Dec.22.2011:           rv_p_ohrs(:)=0d0
         sham%rv_p_ohrs(:)=0d0
      else
ckino Dec.22.2011:           allocate(rv_p_ohrs(2*ndhrs**2*nttabs*nsp))
         allocate(sham%rv_p_ohrs(2*ndhrs**2*nttabs*nsp))
ckino Dec.22.2011:           rv_p_ohrs(:)=0d0
         sham%rv_p_ohrs(:)=0d0
      endif

C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
ckino Dec.22.2011:         allocate(rv_p_oqsig(abs(3*mxkp)))
       allocate(sham%rv_p_oqsig(abs(3*mxkp)))
ckino Dec.22.2011:         if (3*mxkp<0) rv_p_oqsig(:)=0.0d0
       if (3*mxkp<0) sham%rv_p_oqsig(:)=0.0d0

       allocate(rv_p_oqp(abs(3*mxkp)))
       if (3*mxkp<0) rv_p_oqp(:)=0.0d0

      allocate(gstar_iv(abs(-mxkp-1)))
      if (-mxkp-1<0) gstar_iv(:)=0

      gstar_iv = - 2

      allocate(ipq_iv(mxkp))

      allocate(wgt_rv(abs(-mxkp)))
      if (-mxkp<0) wgt_rv(:)=0.0d0

      call info(20,1,0,
     .' q-points in full BZ where sigma calculable ...',0,0)
        call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
ckino Dec.22.2011:       .  , 0 , ipq_iv , rv_p_oqsig , wgt_rv , nqsig , mxkp , 0 , 0 ) 
     .  , 0 , ipq_iv , sham%rv_p_oqsig , wgt_rv , nqsig , mxkp , 0 , 0 ) 




      sham%nqsig=nqsig
ckino Dec.22.2011:         sham%rv_p_oqsig => rv_p_oqsig 


      call dpzero ( wgt_rv , mxkp )

      call info2(20,0,0,
     .' Irr. qp for which sigma is calculated ...',0,0)
ckino Dec.22.2011:          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og 
        call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_p_osymgr 
     .  , nsgrps , ipq_iv , rv_p_oqp , wgt_rv , nqps , mxkp , gstar_iv 
     .  , 0 ) 



Ckino isanrg is logical function,       call isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
      l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
C     call prmx('q for which sigma is read',w(oqp),3,3,nqp)
      call info2(30,0,0,
     .'%?#(n==1)#%1fFile sigm has k-offset=%3:2g '//
     .'relative to the above#%j#'//
     .' ',isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi)

C ... titus
CC     List of sites
C      call getoffH(w(ooffH),offH,n0H,nkap0,nbas)
C      ib1 = 1
C      ib2 = nbas
C      lidim = offH(4,1,nbas+1)
C
C      print*,'lidim',lidim
C
C      do 100 ib = ib1, ib2
C      rewind ifis
C      call iosigh(0,nscnts,nsp1,ndimh,nk1,nk2,nk3,nqp,
C     .  lshft(1),lshft(2),lshft(3),ifis)
C      nlmaa = offH(4,1,ib+1) - offH(4,1,ib)
C      print*,'nlmaa',nlmaa
C      print*,'offH(4,1,ib)',offH(4,1,ib)
C      print*,'offH(4,1,ib+1)',offH(4,1,ib+1)
C ... titus

C --- Generate hrs = sigma(T) from file sigma(k) ---
C ... Loop over spins, keeping sig in full BZ for only one spin
      if (procid .eq. master) then
        do  isp = 1, nsp

C     If sigma file not spin polarized, use sigma from spin 1
          if (isp .eq. 2 .and. lsplts) then
            call iosigh(3,0,1,ndimh,nk1,nk2,nk3,nqp,
     .      lshft(1),lshft(2),lshft(3),ifis)
          endif

C     If to mix in some other sigma file, set up ifis2
          if (alf(2) .ne. 0) then
            if (isp .eq. 1 .or. isp .eq. 2 .and. lsplts) then
              ifis2 = fopna('sigm1',-1,4)
              call iosigh(3,0,nsp,ndimh,nk1,nk2,nk3,nqp,
     .        lshft(1),lshft(2),lshft(3),ifis2)
            endif
          endif




cccccccccccccccccccccccccccccccccccccccccccc
          goto 8888
cccccccccccccccccccccccccccccccccccccccccccc





C ... File sigma in irr BZ -> sigma in full BZ by rotations
c      call defcc(osfz,k1*k2*k3*ndimh**2)

C ... titus
C      call defcc(osfz,k1*k2*k3*ndimh*nlmaa)
C      call defcc(osfz1,k1*k2*k3*ndimh*nlmaa)
C ... titus
          allocate(wk_zv(ndimh**2))

Changenglob          nl = nglob('nl')
          nl = globalvariables%nl
          nspc = 1
          allocate(sigm_zv(ndimh**2))

ccccccccccccccccccccccccc
ctakaox
c      call defcc(oss,ndimh**2)
ccccccccccccccccccccccccc
cccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          do  iq1 = 1, nqps
            read(ifis) qp,  eseavr !june2009  add eseavr
            if(debugmode()>0) print *,' uuu2 iq1 eseavr=',iq1,eseavr
            call dpscop(qp,tmp,3,1,1,1d0)
             call dpsadd ( tmp , rv_p_oqp , 3 , 1 , 3 * iq1 - 2 , - 1d0 ) 

            call dpsadd(tmp,qoffi,3,1,1,-1d0)

            !tmp= qp - w((oqp(3*(iq1-1)+1:3*(iq1-1)+3))+qoffi)
            if (.not. latvec(1,tolq,plat,tmp)) then
               call dpscop ( rv_p_oqp , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 ) 

              if (lssym .ge. 4) then
                call info5(20,0,0,' rdsigm (warning) qp %i : '//
     .          'expected%3;8,4D read%3;8,4D',iq1,tmp,qp,0,0)
                call dpscop(tmp,qp,3,1,1,1d0)
              else
                print 456, iq1,tmp,qp
  456           format(' rdsigm: error on reading qp no',i4/
     .          ' Expected qp=',3f12.6/
     .          '     File qp=',3f12.6)
                call rx(' incompatible q-mesh')
              endif
            endif
C       File read sigma
            print *,' nscnts =',nscnts
            if (nscnts == 5) call dpdump(tmp,1,ifis)
            call dpdump ( sigm_zv , ndimh * * 2 * 2 , ifis )


C        if (isp .eq. 2) then
C          print *, iq1,isp
C         call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
C        endif

cccccccccccccccccccccccccccccccccccccc
ctakaox !see also bndfp.F and hambls.F
c        if(iq1==1) then
c          print *,'reading ovlmat'
c          ifio = fopna('ovlmat',-1,4)
c        endif
c        call dpdump(w(oss),ndimh**2*2,ifio)
c        print *,'read ovlmat for iq1=',iq1
c        call divctakao(w(osigm),w(oss),ndimh)
ccccccccccccccccccccccccccccccccccccc
c

c ... scalced sigm. alf(1:2)
            if (alf(1) .ne. 1d0) then
              call dscal ( ndimh * * 2 * 2 , alf ( 1 ) , sigm_zv , 1 )

            endif
            if (alf(2) .ne. 0d0) then
              read(ifis2) qp
              call dpscop(qp,tmp,3,1,1,1d0)
               call dpsadd ( tmp , rv_p_oqp , 3 , 1 , 3 * iq1 - 2 , - 1d0 ) 

              call dpsadd(tmp,qoffi,3,1,1,-1d0)
              if (abs(tmp(1))+abs(tmp(2))+abs(tmp(3)) .gt. tolq) then
                 call dpscop ( rv_p_oqp , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 ) 

                if (lssym .ge. 4) then
                  call info5(20,0,0,' rdsigm (warning) qp %i : '//
     .            'expected%3;8,4D read%3;8,4D (file sigm1)',iq1,tmp,qp,
     .            0,0)
                else
                  call rx(' file sigm1 has incompatible irr mesh')
                endif
              endif
C         call dpsadd(qp,w(oqp),3,1,3*iq1-2,1d0)
              allocate(sigm2_zv(ndimh**2))

              call dpdump ( sigm2_zv , ndimh * * 2 * 2 , ifis2 )

              call daxpy ( ndimh * * 2 * 2 , alf ( 2 ) , sigm2_zv , 1 , 
     .        sigm_zv , 1 )

              if (allocated(sigm2_zv)) deallocate(sigm2_zv)

            endif

c ... rotate sigm. w(osigm) at qp ---> w(ofbz) at stars of qp.
C       call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
             write(6,"(a,13f13.5)")'goto hamfbk3 xx input qp=', qp
c        napw_in=0 ! for test

            i_copy_size=size(slat%qlat)
            call dcopy(i_copy_size,slat%qlat,1,qlat,1)


            call hamfb3k ( qp , iq1 , nk1 , nk2 , nk3 , k1 , k2 , k3 , ipq_iv 
     .      , napw_in , ndimh , ndimh , ndimh , qb , plat , qlat , ifac , 
     .      gstar_iv , sigm_zv , sfz ( 1 , 1 , 1 , 1 , 1 , isp ) )

            if(debugmode()>0) write(6,"(a,3f13.5)")'end of hamfbk3'
c
c        i = 100*kcplx + 00
c        if (lphase) i = i+40
c        print *,' hamfb3: i =',i
c        call hamfb3(nbas,nl,w(ooffH),w(oiprmb),i,w(opos),iq1,nk1,
c     .    nk2,nk3,k1,k2,k3,w(oipq),w(oistab),w(og),w(oag),
c     .    w(ogstar),ifac,ndimh,ndimh,ndimh,nspc,qb,w(osigm),w(owk),
c     .    w(owk),w(osfz))

          enddo





ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 8888     continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc








ccccccccccccccccccccccccccccccccccccccccccc
ctakaox !see also bndfp.F and hambls.F
c      print *,' end of iq1 loop'
c      call fclose(ifio)
c      call rlse(oss)
c      print *,' end of iq1 loop xxxx'
cccccccccccccccccccccccccccccccccccccccccccc


C      print *, 'one element of sigma(k)'
C      call pvtrod(w(osfz),k1,k2,k3,ndimh)

C ... FT sfz and copy sfz(T) to hrs
          i = 10*hreal + 1
          if (ltrans) i = i + 10000
C ... titus
C      call hft2rs(i,nk1,nk2,nk3,k1,k2,k3,qoff,isp,nsp,nbas,w(og),w(oag)
C     .  ,nsgrps,rsrnge,w(ooffH),ib,ib,nlmaa,ndimh,w(osfz),plat,
C     .  w(opos),ontabs,
C     .  oiaxs,ndhrs,w(ohrs))
C ... titus
               call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , isp 
ckino Dec.22.2011:       .         , nsp , nbas , rv_p_og , rv_p_oag , nsgrps , rsrnge , iv_p_ooffh 
ckino Dec.22.2011:       .         , nsp , nbas , slat%rv_p_osymgr , rv_p_oag , nsgrps , rsrnge , iv_p_ooffh 
ckino Dec.22.2011:       .         , nsp , nbas , slat%rv_p_osymgr , slat%rv_p_oag , nsgrps , rsrnge , iv_p_ooffh 
     .         , nsp , nbas , slat%rv_p_osymgr , slat%rv_p_oag , nsgrps , rsrnge , sham%iv_p_ooffh 
     .         , 1 , nbas , ndimh , ndimh , sfz ( 1 , 1 , 1 , 1 , 1 , isp ) 
ckino Dec.14.2011:       .         , plat , rv_p_opos ,iv_p_ontabs , iv_p_oiaxs , ndhrs , rv_p_ohrs ) 
ckino Dec.14.2011:       .         , plat , rv_p_opos ,iv_p_ontabs , iv_a_oiaxs , ndhrs , rv_p_ohrs ) 
ckino Dec.22.2011:       .         , plat , rv_p_opos ,iv_a_ontabs , iv_a_oiaxs , ndhrs , rv_p_ohrs ) 
ckino Dec.22.2011:       .         , plat , slat%rv_p_opos ,iv_a_ontabs , iv_a_oiaxs , ndhrs , rv_p_ohrs ) 
ckino Dec.27.2011:       .         , plat , slat%rv_p_opos ,iv_a_ontabs , sham%iv_a_oiaxs , ndhrs , sham%rv_p_ohrs ) 
ckino 12月.27.2011:       .         , plat , slat%rv_p_opos ,sham%iv_a_ontabs , sham%iv_a_oiaxs , ndhrs , sham%rv_p_ohrs ) 
     .         , plat , slat%rv_p_opos ,sham%iv_a_ontabs , sham%iv_a_oiaxs , ndhrs , sham%rv_p_ohrs ) 






C      print *, 'print out sigma(T)'
C      call pvtrof(w(ohrs),ndhrs,nsp,nk1,nk2,nk3,w(oiaxs),1,
C     .  nttabs,w(osfz),ndimh)
C      stop

C ... End loop over spins
c      call rlse(osfz)
        enddo
C ... titus
C 100  continue
C     exit of MPI master loop
      endif
      if (allocated(sigm_zv)) deallocate(sigm_zv)
      if (allocated(wk_zv)) deallocate(wk_zv)
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq_iv)) deallocate(ipq_iv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)



      sham%ndhrs=ndhrs
ckino Dec.14.2011:         sham%iv_p_ontabs => iv_p_ontabs 
ckino Dec.14.2011:         sham%iv_a_ontabs => iv_a_ontabs 
ckino 12月.27.2011:        call move_alloc(from= iv_a_ontabs ,to= sham%iv_a_ontabs )

ckino Dec.14.2011:         sham%iv_p_oiaxs => iv_p_oiaxs 
ckino Dec.14.2011:         sham%iv_a_oiaxs => iv_a_oiaxs 
ckino Dec.27.2011:        call move_alloc(from= iv_a_oiaxs ,to= sham%iv_a_oiaxs )

ckino Dec.22.2011:         sham%rv_p_ohrs => rv_p_ohrs 



C ... Check that FT replicates sigma(k) at file qp to tolerance rsstol
C     rsstol = 1d-6
C     print *, '!!'
C     if (alf(1) .ne. 1 .or. alf(2) .ne. 0 .or. .true.) then
      if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
        call info0(20,0,0,' skipping check of Bloch summed sigma'//
     .  ' (sigma scaled) ')
      else
        j1 = 0
        if (lssym .ge. 4) j1 = 2
        if (procid .eq. master) then
          call dpzero(tmp,3)
               call chksgr ( j1 , ltrans , kcplx , plat , nsp , ndimh , ifis 
ckino Dec.14.2011:       .         , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.22.2011:       .         , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.22.2011:       .         , hreal , sham%iv_p_oindxo , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.27.2011:       .         , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .         , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .         , rsstol , i , rv_p_oqp , nbas , 0 , rotm , iwdummy ) 





        endif
      endif

      if ( .not. ( cmdopt ( '--wsig' , 6 , 0 , outs ) .or.cmdopt ( 
     .'-wsig' , 5 , 0 , outs ) ) ) then
         if (associated(rv_p_oqp)) deallocate(rv_p_oqp)

      endif


C ... Broadcast R.S. sigma before symmetrization
      if (hreal .eq. 1) then
ckino Dec.22.2011:           call mpibc1 ( rv_p_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog 
         call mpibc1 ( sham%rv_p_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog 
     .   , 'rdsigm' , 'hrs' ) 

      else
ckino Dec.22.2011:           call mpibc1 ( rv_p_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4 
         call mpibc1 ( sham%rv_p_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4 
     .   , mlog , 'rdsigm' , 'hrs' ) 

      endif

C --- Symmetrize hrs ---
C     Best to use nsgrp here instead of nsgrps.
C     However, neighbor table was generated for nsgrps, and may be
C     increase if nsgrp>nsgrps.  This branch will fail in that case.
      if (mod(lrsig,4) .ge. 2 .and. nsgrp .gt. 1) then

cC       Allocate memory for hrs
c        if (hreal .eq. 1) then
c          call defdr(ohrss,ndhrs**2*nttabs*nsp)
cCchp1           call dcopy(ndhrs**2*nttabs*nsp*1,w(ohrs),1,w(ohrss),1)
c           call dcopy ( ndhrs * * 2 * nttabs * nsp * 1 , rv_p_ohrs , 1 , 
c     .     w ( ohrss ) , 1 ) 
c        else
c          call defdc(ohrss,ndhrs**2*nttabs*nsp)
cCchp1           call dcopy(ndhrs**2*nttabs*nsp*2,w(ohrs),1,w(ohrss),1)
c           call dcopy ( ndhrs * * 2 * nttabs * nsp * 2 , rv_p_ohrs , 1 , 
c     .     w ( ohrss ) , 1 ) 
c        endif
ckino Dec.22.2011:          nhrss= size(rv_p_ohrs)
        nhrss= size(sham%rv_p_ohrs)
        allocate(hrss(nhrss))
ckino Dec.22.2011:          hrss = rv_p_ohrs
        hrss = sham%rv_p_ohrs

Changenglob        mxorb = nglob('mxorb')
        mxorb = globalvariables%mxorb
C       10s digit distributes ri-rj pairs to avg (ri-rj),(rj-ri) pairs
C       Use 10*2 for hermitian matrices, 10*1 for symmetric ones
        i = 1-hreal + 10*2 + 100*(1-hreal) + 100000
        allocate(istb2_iv(nsgrp*nbas))

ckino Dec.22.2011:           call istbpm ( iv_p_oistab , nbas , nsgrp , istb2_iv ) 
         call istbpm ( slat%iv_p_oistab , nbas , nsgrp , istb2_iv ) 


Changenglob        nl = nglob('nl')
        nl = globalvariables%nl
ckino Dec.22.2011:                call rsmsym ( i , plat , mxorb , iv_p_oiprmb , ndimh , nbas , 
              call rsmsym ( i , plat , mxorb , sham%iv_p_oindxo , ndimh , nbas , 
ckino Dec.14.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_p_ontabs , iv_p_oiaxs 
ckino Dec.14.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_p_ontabs , iv_a_oiaxs 
ckino Dec.22.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_a_ontabs , iv_a_oiaxs 
ckino Dec.27.2011:       .        slat%rv_p_opos , nl , nsp , 1 , nttabs , iv_a_ontabs , sham%iv_a_oiaxs 
ckino 12月.27.2011:       .        slat%rv_p_opos , nl , nsp , 1 , nttabs , sham%iv_a_ontabs , sham%iv_a_oiaxs 
     .        slat%rv_p_opos , nl , nsp , 1 , nttabs , sham%iv_a_ontabs , sham%iv_a_oiaxs 
ckino Dec.22.2011:       .        , rv_p_og , istb2_iv , nsgrp , ndhrs , hrss , rv_p_ohrs ) 
ckino Dec.22.2011:       .        , slat%rv_p_osymgr , istb2_iv , nsgrp , ndhrs , hrss , rv_p_ohrs ) 
     .        , slat%rv_p_osymgr , istb2_iv , nsgrp , ndhrs , hrss , sham%rv_p_ohrs ) 



        if (allocated(istb2_iv)) deallocate(istb2_iv)
c        call rlse(ohrss)
        deallocate(hrss)

C        Debugging
C        ifis2 = fopna('out',-1,0)
C        print *, 'hello',ndhrs,nttabs
C        call ywrm(0,'sig',2,ifis2,'(9f20.10)',w(ohrs),1,ndhrs**2,
C     .    ndhrs**2,nttabs*nsp)
C        call rx0('done')

C   ... Check how well symmetrized FT replicates sigma(k) at file qp
        if (procid .eq. master) then
          rsstol = 0
          if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
            call info0(20,0,0,' skipping check of symmetrized sigma'//
     .      ' (sigma scaled) ')
          else
                call chksgr ( 0 , ltrans , kcplx , plat , nsp , ndimh , ifis 
ckino Dec.14.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.22.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.22.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.27.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .          , rsstol , i , 0 , nbas , 0 , rotm , iwdummy ) 




          endif
        endif

C   ... Broadcast R.S. sigma after symmetrization
        if (hreal .eq. 1) then
ckino Dec.22.2011:             call mpibc1 ( rv_p_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog 
           call mpibc1 ( sham%rv_p_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog 
     .     , 'rdsigm' , 'hrs' ) 

        else
ckino Dec.22.2011:             call mpibc1 ( rv_p_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4 
           call mpibc1 ( sham%rv_p_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4 
     .     , mlog , 'rdsigm' , 'hrs' ) 

        endif

      endif

C ... lwsig=0 modes that internally transform sigma, save and exit
      if (lwsig .eq. 0) then
        if (cmdopt('--wsig',6,0,outs) .or. cmdopt('-wsig',5,0,outs)) then
          if (procid .eq. master) then
            ifis2 = fopna('sigm2',-1,4)
            rsstol = 0
            out2 = outs(7:)
            if (outs(1:5) .eq. '-wsig') out2 = outs(6:)
            dc = out2(1:1)
            call sigflg(dc,'wsig',out2,ledit,lwsig,lrot,lfbzout,lnwmsh,
     .      lonesp,qoffo,rots)

            if (lonesp .eq. 1) then
              if (nsp .eq. 1) call rx(
     .        'onesp option nonsensical unless sigma file is spin pol')
ckino Dec.22.2011:                     call siged ( 1 , nbas , nsp , ndhrs , plat , rv_p_opos , ndimh 
                   call siged ( 1 , nbas , nsp , ndhrs , plat , slat%rv_p_opos , ndimh 
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_p_oiaxs , rv_p_ohrs 
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_a_oiaxs , rv_p_ohrs 
ckino Dec.22.2011:       .             , iv_p_oiprmb , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs 
ckino Dec.22.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs 
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
ckino 12月.27.2011:       .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
     .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
     .             ) 





            endif

            if (ledit .eq. 1) then
ckino Dec.22.2011:                     call siged ( 0 , nbas , nsp , ndhrs , plat , rv_p_opos , ndimh 
                   call siged ( 0 , nbas , nsp , ndhrs , plat , slat%rv_p_opos , ndimh 
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_p_oiaxs , rv_p_ohrs 
ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_a_oiaxs , rv_p_ohrs 
ckino Dec.22.2011:       .             , iv_p_oiprmb , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs 
ckino Dec.22.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs 
ckino Dec.27.2011:       .             , sham%iv_p_oindxo , hreal , iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
ckino 12月.27.2011:       .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
     .             , sham%iv_p_oindxo , hreal , sham%iv_a_ontabs , sham%iv_a_oiaxs , sham%rv_p_ohrs 
     .             ) 





            endif

C       Handle case q-mesh for file output differs from file input
C       ltmp = true in cases where file symops artificially reduced
C        ltmp = mod(lssym,2) .eq. 1 .and. lwsig .ne. 3
C        if (lnwmsh .or. ltmp) then
            if (lnwmsh) then

              i_copy_size=size(sbz%nkabc)
              call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
              i_copy_size=size(sbz%lshft)
              call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)
              nqp=sbz%nkp
               rv_p_oqp => sbz%rv_p_oqp 



              mxkp = nk1*nk2*nk3
               allocate(rv_p_oqp(abs(3*mxkp)))
               if (3*mxkp<0) rv_p_oqp(:)=0.0d0

              allocate(gstar_iv(abs(-mxkp-1)))
              if (-mxkp-1<0) gstar_iv(:)=0

              gstar_iv = - 2

              allocate(ipq_iv(mxkp))

              allocate(wgt_rv(abs(-mxkp)))
              if (-mxkp<0) wgt_rv(:)=0.0d0

              do  i = 1, 3
                llshft(i) = lshft(i) .ne. 0
              enddo
              call dpzero ( wgt_rv , mxkp )

              if (lfbzout) then
                call info0(20,1,0,
     .          ' Creating sigma for new mesh, full BZ (file sigm2) ...')
                  call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
     .            , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 ) 



              else
                call info0(20,1,0,
     .          ' Creating sigma for new mesh (file sigm2) ...')
ckino Dec.22.2011:                    call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og 
                  call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_p_osymgr 
     .            , nsgrp , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , gstar_iv 
     .            , 0 ) 




              endif

            elseif (lfbzout) then
cki              oqp = oqsig
ckino Dec.22.2011:                rv_p_oqp => rv_p_oqsig
              rv_p_oqp => sham%rv_p_oqsig
              nqp = nqsig
              call info2(20,1,0,
     .        ' Creating sigma for full BZ, %i qp (file sigm2) ...',nqp,0)
            else
              call info2(20,0,0,
     .        ' Writing sigma for irr BZ, %i qp (file sigm2) ...',nqp,0)
            endif

C   ... Setup rotation matrix
            if (rots .ne. ' ') then
              call a2rotm(rots,.false.,0,rotm)
              if (iprint() .ge. 20) then
                call info0(20,0,0,' Rotate sigma by rotation matrix:')
                write (stdo,350) ((rotm(i,j),j=1,3),i=1,3)
  350           format(3f11.6)
              endif
              lrot = lrot + 1
            endif

C   ... Setup phase matrix
            if (lrot .ge. 2) then
              allocate(delt_rv(abs(-3*nbas)))
              if (-3*nbas<0) delt_rv(:)=0.0d0

              call info0(0,0,0,' Phase shift of sigma.  '//
     .        'Translation vectors from shorps file:')
              call iopos ( .false. , - 1 , 'shorps' , nbas , delt_rv )

            endif
C       call prmx('q for which sigm is written',w(oqp),3,3,nqp)

            call iosigh(0,nscnts,nsp,ndimh,nk1,nk2,nk3,nqp,
     .      lshft(1),lshft(2),lshft(3),-ifis2)
            if (lqoffo .ne. 0) then
C         call prmx('qp',w(oqp),3,3,nqp)
              do  i = 1, nqp
                 call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 * 
     .           i - 2 , 1 , 1d0 ) 

              enddo
C         call prmx('qp',w(oqp),3,3,nqp)
            endif
                call chksgr ( 1 , ltrans , kcplx , plat , nsp , ndimh , ifis2 
ckino Dec.14.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.22.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.22.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs 
ckino Dec.27.2011:       .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .          , hreal , sham%iv_p_oindxo , nttabs , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs 
     .          , rsstol , nqp , rv_p_oqp , nbas , lrot , rotm , delt_rv ) 





            call fclose(ifis2)
            call rx0('done writing sigma, file sigm2')
C      elseif (lwsig .eq. 1 .or. lwsig .eq. 2) then
C        if (lfbzout) call info0(20,0,0,' ')
C        if (lfbzout) goto 99
          endif
        endif
      endif

C     Broadcast entire sham structure (shouldn't be needed)
C     call mpibc1(sham,i,4,mlog,'rdsigm','sham')
C     Broadcast neighbor table (shouldn't be needed)
C     call mpibc1(w(ontabs),nbas+1,2,mlog,'rdsigm','ntabs')
C     call mpibc1(w(oiaxs),niax*nttabs,2,mlog,'rdsigm','iaxs')

C ... Transformations done by caller: write qp list to disk
   99 continue
      if (lwsig .ne. 0) then

        call info0(30,0,0,' ')

        if (lnwmsh) then

          i_copy_size=size(sbz%nkabc)
          call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
          i_copy_size=size(sbz%lshft)
          call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)
          nqp=sbz%nkp
           rv_p_oqp => sbz%rv_p_oqp 


          call info2(20,0,0,
     .    ' Transform sigma on new k-mesh:  '//
     .    'nk=%3:1i  shft=%3:1l',nkxyz,lshft)
        endif

        if (lfbzout) then
          mxkp = nk1*nk2*nk3
          call info(20,0,0,' rdsigm: use full BZ (%i q-points)...',mxkp,0)
           allocate(rv_p_oqp(abs(3*mxkp)))
           if (3*mxkp<0) rv_p_oqp(:)=0.0d0

          allocate(ipq_iv(mxkp))

          allocate(wgt_rv(abs(-mxkp)))
          if (-mxkp<0) wgt_rv(:)=0.0d0


          i_copy_size=size(slat%plat)
          call dcopy(i_copy_size,slat%plat,1,plat,1)

C   ... Make is,ifac,qb,qlat
          do   i = 1, 3
            llshft(i) = lshft(i) .ne. 0
          enddo
          call pshpr(0)
          call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
          call poppr
            call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
     .      , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 ) 



C       call pack2('bz nkp oqp',sbz,nqp,oqp)
          if (allocated(delt_rv)) deallocate(delt_rv)
          if (allocated(wgt_rv)) deallocate(wgt_rv)

        endif

C ... Write points to disk, shifting by qoffo
        if (lqoffo .ne. 0) then
C       call prmx('qp',w(oqp),3,3,nqp)
          do  i = 1, nqp
             call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 * 
     .       i - 2 , 1 , 1d0 ) 

          enddo
C       call prmx('qp',w(oqp),3,3,nqp)
        endif

        call info0(30,0,0,' RDSIGM: writing data to file QPTS ...')
        ifiz = fopna('qpts',-1,0)
          call getqp ( 1 , - ifiz , nqp , nkxyz , lshft , 0 , rv_p_oqp 
     .    , iwdummy , iwdummy ) 



        if (lfbzout) then
          if (allocated(ipq_iv)) deallocate(ipq_iv)

        endif

      endif

      call tcx('seneinterp')

      if (allocated(gstar_iv)) deallocate(gstar_iv)

      end subroutine seneinterp





      subroutine rdsigm2(nbas,nsp,ndimh,slat,sham,sbz,ifis, 
     &nk1,nk2,nk3,ldim,qsmesh,sfz)

      use m_struc_def  !Cgetarg
      use m_globalvariables
      use m_hamindex
C- takao. Expand self-energy (read by ifis) to all the q point on mesh. In developing.
Ci ifis:  file hundle for self-energy file sigm. only at irreducible q point.
Co complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp):  self-energy (\Sigma-Vxc) for  all the q points on mesh.
Co real(8):: qsmesh(3,nk1,nk2,nk3)
c  We have to clean up this routine. The purpose of this routine is "read sigm file and expand it in full BZ".
c  Not do more than that. (in future, we do scaling of simga in bndfp.F.
c  Many un-used local variables are contained.
c  Especially qsmesh (regular q mesh for self-energy.) is very problematic. It should be given at a place, and then
c  it should be used somewhere else.
c
c  original text below are not so meanigful. but kept for hints.
c
C- Read sigm(k) from file and generate sigm(R)
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 do not read self-energy; just exit
Ci         :1 read and FT sigma(k) to make sigma(T)
Ci         :2 symmetrize sigma(T)
Ci         :4 take the real part of sigma(T) only
Ci         :10s digit (used here only for printout)
Ci         :0 Simple bloch sum of sigma
Ci         :1 approx high- and low-energy sigma with diagonal
Ci         :2 perturbation approach, diagonalizing exactly only
Ci         :  the first kappa block.
Ci         :3 Linear interpolation of sigma from neighboring p
Ci         :  sig(q) = sum_i wt_i A+_i sig(qpi) A_i
Ci         :  where A_i = z^LDA_qi (z_q^LDA)^-1
Ci         :100s digit specifies number of interpolation points
Ci         :    for interpolation mode 3.
Ci         :  0 -> use default (4 points)
Ci         :1000s digit specifies ASA conventions:
Ci         :      poke sigm(k)_RL,R'L' -> sigm(T)_RL,R'L'
Ci         :      (10000s digit of hft2rs)
Ci         :10000s digit
Ci         :  1 specifies file sigm stored with no group operations
Ci         :    sigm is stored at k-points in the full BZ.
Ci         :    Equivalent to --rsig:fbz
Ci         :  2 do not force R.S. neighbor table to conform to
Ci         :    symmetry
Ci         :  4 do not force file qp to match those generated by
Ci              bzmesh.  Use qp from bzmesh.
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ndimh :hamiltonian dimension
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat opos nsgrp oistab osymgr oag
Ci     Stored:
Ci     Passed to:
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: nkabc lshft nkp oqp
Ci         : (used only when rdsigm writes sigma on a new mesh)
Ci   ifis  :file logical unit for sigma
Ci   rsrnge:maximum length for connecting vectors in sigma(R)
Cio Inputs/Outputs:
Cio  sham  :struct for parameters defining hamiltonian; see routine uham
Cio    Elts read: ldham ooffH oindxo sigp rsstol
Cio    Stored:    nqsig oqsig ndhrs ontabs oiaxs ohrs
Cio Outputs:
Co   lwsig :0  No special transformations of sigma
Co         :1  Mode transforms sigma from orbital into LDA basis.
Co         :   rdsigm reads sigma from file sigm in orbital basis;
Co         :   Calling program generates and stores transformed sigma
Co         :2  Similar to lwsig=1, except
Co             low- and high- energy blocks replaced by diagonal parts
Co         :-1 Mode transforms sigm from LDA to orbital basis
Co         :   (Inverse operation of lwsig=1 or lwsig=2).
Co         :   This mode requires both sigma and LDA eigenvectors go be
Co         :   stored on disk in files 'sigm' and 'evec'
Co         :   rdsigm generates and stores transformation in file 'sigm2'
Co         :3  Returns lwsig=3 to flag calling program.  It should generate and
Co         :   store LDA eigenvalues and eigenvectors. No sigm file is read.
Co         :4  Returns lwsig=4 to flag calling program.  It should generate and
Co         :   store eigenvalues and eigenvectors.
Co         :5  Writes sigm as in lwsig=2 with low- and high- energy blocks
Co         :   replaced by diagonal parts, but writes sigm in orbital bas.
Cl Local variables
Cl   sigp  :parameters for approximating self-energy sigma.  sigma
Cl         :is approximated by its diagonal part sigii for energies
Cl         :below a low-energy cutoff (specified nmin or emin) and
Cl         :above a low-energy cutoff (specified nmax or emax).
Cl         : arg 1: specifies how to set diagonal part sigii
Cl         :        for states above the high-energy cutoff nmax or emax
Cl         :        0 constrain sigii to be > asig+bsig*e
Cl         :        1 constrain sigii to be = asig+bsig*e
Cl         :        2 constrain sigii to be > asig and < bsig
Cl         :        3 constraint same as case 1.
Cl         :          arg1=3 differs in that the least-squares fit to
Cl         :          sigii (for informational purposes only, to help
Cl         :          estimate asig and bsig) is done for states between
Cl         :          efit and nmax or emax
Cl         : arg 2: nmin : sigma for states 1..nmin are approximated by sigii
Cl         : arg 3: emin : (used only if nmin<0)
Cl         :             : sigma for levels e<emin are approximated by sigii
Cl         : arg 4: nmax : sigma for levels i>nmax are approximated by
Cl         :               sigii AND constrained according to mode
Cl         : arg 5: emax : (used only if nmax<=0)
Cl         :             : sigma for levels e<emax are approximated by
Cl         :               sigii AND constrained according to mode
Cl         : arg 6: asig : constraint used to approximate
Cl         :               sigii = asig + E * bsig  or
Cl         :               asig < sigii < bsig
Cl         : arg 7: bsig : constraint used to approximate
Cl         :               sigii = asig + E * bsig  or
Cl         :               asig < sigii < bsig
Cl         : arg 8: efit : (mode 3) energy minimium
Cl                         for fitting asig and bsig (not used here)
Cl   Note  :sigp takes a similar, but slightly different meaning when
Cl         :used in conjuction with sigm(LDA)->sigm(orbital), i.e. lwsig=-1
Cl         :In this case, sigp applies to the new basis, whose
Cl         :eigenvalues and eigenvectors are stored in 'evec.'
Cl         :Arguments are used in rotevs, with the following relations:
Cl         : arg1 (mode) has same meaning as before
Cl         : arg2 (nmin) takes meaning of lcore.
Cl         :       nmin>0 => sigm does not contain this block,  but the
Cl         :                 new basis does.  Use in emin (arg3) in
Cl         :                 place of sigm for this block .
Cl         :       nmin<0 => new basis does not contain this block
Cl         : arg3 (emin) If nmin>0, use emin for diag. sigma, this block
Cl         : arg4, arg5 (nmax,emax) could play the role of lhigh in
Cl         :       rotevs.  However, they are not used.  lhigh is
Cl         :       determined from constraint ndims+lccore+lhigh = ndimz
Cl
Cl   hreal :0 allow sigm(R) to be complex
Cl         :1 assume sigm(R) is real
Cl  ltrans :specifies indexing of sigm(T)_RL,R'L'L
Cl         : ltrans = F taken from sigm(k)_RL,R'L'
Cl         : ltrans = T taken from sigm(k)_RL,R'L'
Cl  kcplx  :0 sigm(k) has real, imaginary separated
Cl         :1 sigm(k) is in complex*16 format:
Cl         :2 sigm(k) has real, imaginary separated by columns
Cl  lphase :T if phase convention phi = q * [(g R_j + a) - R_i]
Cl         :  for rotations should be scaled by -1
Cl  lssym  :10000s digit mode
Cl  lfbzin :flags whether input self-energy file has suppressed symops
Cl  nscnts :compound of switches containing file contents of sigma
Cl  qoffi  :k-mesh offset for input sigma file
Cl  qoffo  :k-mesh offset for generated sigma file
Cr Remarks
Cb Bugs
Cb   Routine should be able to accomodate offset qp case (ifac<>0)
Cb   Routine should be able to accomodate case when a different
Cb   list of irreducible qp is used than the one generated by
Cb   bzmesh.
Cu Updates
Cu   27 Jun 08 Redesigned transformation of sigma to new basis
Cu   15 May 07 Parallelized symmetrizer
Cu   24 Jan 07 New option --wsig:onesp
Cu   20 Jan 07 Bug fix, --mixsig, spin polarized case
Cu   16 Jan 07 New option lwsig=-1; updated iosigh
Cu   24 Jul 06 MPI broadcast of sigma
Cu   24 Feb 05 Rotation of sigma matrix now in hrs
Cu             Switch lwsig to flag transformed sigm to be generated
Cu   27 Jan 05 New 40000s digit mode (allow qp mismatch)
Cu   20 Dec 04 New --wsig:rot and --wsig:phase switches
Cu             New 10000s digit mode
Cu   22 Nov 04 New --wsig:newkp switch
Cu   25 Sep 04 rdsigm can read sigm(q) file for offset q-mesh
Cu             rdsigm can write sigm(q) for mesh sbz->nkabc
Cu    5 Apr 04 bug fix for --wsig switch.  New --mixsig switch
Cu   15 Feb 04 Enable rdsigm to read a linear combination of
Cu             data from file ifis and data from file 'sigm1'
Cu    4 Jan 04 rdsigm checks bloch sum against ham->rsstol
Cu   10 Oct 03 rdsigm can spin-split a non-spin-polarized sigma
Cu   23 Sep 03 rdsigm can now read sigm fit to ASA hamiltonian
Cu   11 Jul 03 enable write of sigma to full BZ
Cu   24 May 03 Redesigned switches in accord with new interpolation
Cu   09 May 03 Added symmetrization of sigma
Cu   10 Jan 03 Some adaptations for Sergey's interpolation
Cu   14 Aug 02 Added option to orthogonalize sigm.
Cu   27 Jul 02 first created
C ------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,ifis,ndimh,lwsig
C     integer osig(nbas),otau(nbas),oppi(nbas)
      type(s_lat)::slat
      type(s_ham)::sham
      type(s_bz)::sbz

C     double precision ssite(1),sspec(1)
      double precision rsrnge
C ... Local parameters
      logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
     .latvec,lfbzin,lfbzout
      integer parg,isw,lonesp
      character outs*80,out2*80,dc*1,rots*120
      integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
     .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
     .ndhrs,j1,k1,k2,k3,iq1,nl,nspc,mxorb,nqsig,modsgp,nmin,nmax,
     .fopna,kcplx,lrot,iprint,stdo,lssym,ledit,niax,nscnts,
     .ndims,ndimz,iq,n123(4),lcore,lhigh
      parameter (niax=10)
c      equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
                integer:: ohrss , osigm2 , odelt , oistb2 
ckino Dec.14.2011:                  integer,pointer :: iv_p_ontabs(:)
ckino 12月.27.2011:                  integer,allocatable :: iv_a_ontabs(:)

ckino Dec.14.2011:                 integer,pointer :: iv_p_oiaxs(:) 
ckino Dec.27.2011:                 integer,allocatable :: iv_a_oiaxs(:) 

ckino Dec.22.2011:                real(8),pointer :: rv_p_opos(:) 

ckino Dec.22.2011:               real(8),pointer :: rv_p_og(:) 

ckino Dec.22.2011:              real(8),pointer :: rv_p_oag(:)

ckino Dec.22.2011:             real(8),pointer :: rv_p_ohrs(:) 

ckino Dec.22.2011:           integer,pointer :: iv_p_oistab(:) 

ckino Dec.22.2011:           integer,pointer :: iv_p_oiprmb(:)

ckino Dec.22.2011:          real(8),pointer :: rv_p_oqsig(:) 

ckino Dec.22.2011:         integer,pointer :: iv_p_ooffh(:) 

      integer ,allocatable :: gstar_iv(:)
      integer ,allocatable :: ipq_iv(:)
      real(8) ,allocatable :: qp_rv(:)
      real(8) ,allocatable :: wgt_rv(:)
      complex(8) ,allocatable :: wk_zv(:)
      complex(8) ,allocatable :: sigm_zv(:)

      double precision tmp(3)
      real(8),allocatable:: evls(:),evlz(:),sigii(:)
      complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
C     For offset q mesh
      integer is(3),lshft(3),ifac(3) !,lqoffo
      double precision rb(3,3),qb(3,3) !,qoffi(3),qoffo(3)
C MPI
      logical mlog
      integer procid,mpipid,master
C ... for mixing sigma files: Use alf1*sigm + alf2*sigm1
      double precision alf(2)
C     integer os
      double precision plat(3,3),qp(3),tolq,rsstol,sigp(10)  ,qlat(3,3)
      double precision emin,emax,asig,bsig,efit,rotm(3,3),ddot !,qoff(3)
      equivalence (emin,sigp(3)),(emax,sigp(5))
      equivalence (asig,sigp(6)),(bsig,sigp(7)),(efit,sigp(8))
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
c      double precision qk
      integer jj1,jj2,jj3,k
C.....titus
C      integer ib,ib1,ib2,nlmaa,lidim,ldh,osfz1,nsp1
C      integer nkap0,n0H
C      parameter  (nkap0=3,n0H=5)
C      integer offH(n0H,nkap0,nbas)
C.....titus
C ... Heap
       integer ::iwdummy 


      real(8):: eseavr,qsmesh(3,nk1,nk2,nk3)
      integer:: i1,i2,i3,ikt,ldim,napw_in,debugmode
      complex(8)::sfz(nk1,nk2,nk3,ndimh,ndimh,nsp)
      integer:: i_copy_size

      logical:: isanrg, l_dummy_isanrg
ccccccccccccccccccccccccccccccc
c      integer:: oss
ccccccccccccccccccccccccccccccc
c      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
c     .                    (jj2*ifac(2)-1)*qb(k,2) +
c     .                    (jj3*ifac(3)-1)*qb(k,3)

c
      print *
      print *,'rdsigm2:'
      mode=12
      lwsig=0

C --- Setup and printout ---
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

C     tolq is allowed fuzz in qp
      tolq = 1d-6
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      ledit = 0
      lrsig = 2  !mod(mode,10)
      lwsig = 0
      moditp = 1 !mod(mode/10,10)
      alf(1) = 1
      alf(2) = 0
      nscnts = 0
c      lqoffo = 0
      lfbzout = .false.

C     Switch flagging whether input file has symops or not
      hreal = 0
c      if (lrsig .ge. 4) hreal = 1
c      if (lrsig .eq. 0) return
c
c      ltrans = mod(mode/1000,10) .eq. 0
      ltrans=.true.

c      kcplx = 1
c      if (mod(mode/1000,10) .eq. 1) kcplx = 0
      kcplx = 0

c      lphase = mod(mode/1000,10) .eq. 0
c      lssym = mod(mode/10000,10)

      lphase = .true.
      lssym  = 0
Ckino isanrg is logical function, c      call isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)
c      l_dummy_isanrg=isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)

C ... Command-line switches for input sigma and defaults
c      call sigflg(' ','rsig',' ',w,w,w,lfbzin,w,
c     .  lonesp,qoffi,w)

c      if (cmdopt('--rsig',6,0,outs)) then
c        out2 = outs(7:)
c        dc = out2(1:1)
c        call sigflg(dc,'rsig',out2,w,w,w,lfbzin,w,lonesp,qoffi,w)
c      endif
      lonesp = 0
c      qoffi  = 0d0
      lfbzin = .false.

c#if 1
c      if ( lfbzin .or. mod(lssym,2).eq.1 ) then
c      lfbzin = .true.
c      else
c      lfbzin = .false.
c      endif
c#else
c      lfbzin = lfbzin .or. mod(lssym,2)
c#endif

c$$$      if (cmdopt('--mixsig=',9,0,outs)) then
c$$$        j = 0
c$$$        i = parg('--mixsig=',4,outs,j,len(outs),
c$$$     .        ', ',2,2,ifac,alf)
c$$$        if (i .ne. 1 .and. i .ne. 2) call rx
c$$$     .    ('rdsigm: failed to parse arguments to --mixsig=')
c$$$      endif
c      call info(10,1,0,' RDSIGM: read file sigm and create '//
c     .  '%?#n#REAL#COMPLEX# sigma(R) by FT ...',hreal,0)

c      call upack('ham ldham ooffH oindxo sigp rsstol',sham,ldham,ooffH,
c     .  oiprmb,sigp,rsstol)
c      modsgp = nint(sigp(1))
c      nmin   = nint(sigp(2))
c      nmax   = nint(sigp(4))


c$$$c info section
c$$$      call info2(30,0,0,
c$$$     .  '%9fSigm will be approximated by:  '//
c$$$     .  '%?#(n==0)#Simple bloch sum##%-1j'//
c$$$     .  '%?#(n==1)#diagonal Sigma for high and low states##%-1j'//
c$$$     .  '%?#(n==2)#Perturbation for higher kappa blocks##%-1j'//
c$$$     .  '%?#(n==3)#Interpolation from known points##%-1j'//
c$$$     .  ' ',moditp,0)
c$$$      if (moditp .eq. 1 .or. moditp .eq. 2) then
c$$$        call info5(30,0,0,
c$$$     .    '%9fApproximate sigma '//
c$$$     .    '%?#(n<0)#for energies E(lda)<%d; and %-2j#%-1j#'//
c$$$     .    '%?#(n>0)#for states up to %-1jn=%i; and ##%j'//
c$$$     .    '%?#(n<=0)#for energies E(lda)>%d%-2j#%-1j#'//
c$$$     .    '%?#(n>0)#for states above %-1jn=%i##%j',
c$$$     .    nmin,emin,nmax,emax,0)
c$$$        call info5(30,0,0,'%9fFor high states '//
c$$$     .    '%?#(n==0)#Sigii > %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
c$$$     .    '%?#(n==1)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
c$$$     .    '%?#(n==3)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
c$$$     .    '%?#(n==2)#%,3;4d < Sigii < %,3;4d%-2j##%-1j'//
c$$$     .    ' ',modsgp,asig,bsig,0,0)
c$$$        if (nmax .gt. 0) call info5(30,0,0,
c$$$     .    '%9fLinear fit to diagonal sigma for E(lda)>%d'//
c$$$     .    '%?#(n==3)# and n<%i'//
c$$$     .    ' ',efit,modsgp,nmax,0,0)
c$$$        if (nmax .eq. 0) call info5(30,0,0,
c$$$     .    '%9fLinear fit to diagonal sigma for E(lda)>%d'//
c$$$     .    '%?#(n==3)# and E(lda)<%d'//
c$$$     .    ' ',efit,modsgp,emax,0,0)
Ckino isanrg is logical function, c$$$        call isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
c$$$        l_dummy_isanrg=isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
c$$$      endif
c$$$      if (lfbzin .or. ddot(3,qoffi,1,qoffi,1).ne.0) then
c$$$        call info5(30,0,0,
c$$$     .  '%?#(n==1)#%9fFile sigm in saved in FBZ: '//
c$$$     .  'symmetrization suppressed##'//
c$$$     .  '%?#(n==1)#%N%9fFile k-mesh is offset by:%3:2g#%j#'//
c$$$     .  ' ',isw(lfbzin),isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi,0,0)
c$$$        endif
c$$$      if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
c$$$        call info5(30,0,0,
c$$$     .    '%9fUse for sigma: %d*(file sigm)%?#n# + '//
c$$$     .    '%d*(file sigm1)',alf(1),isw(alf(2).ne.0),alf(2),0,0)
c$$$      endif

C --- Read sigma(orbital basis) from file ---
      call tcn('rdsigm2')

      if (procid .eq. master) then
C C       Require file contains sigm(orbital basis)
        call iosigh(0,nscnts,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
     .  lshft(3),ifis)
C         if (nscnts .ne. 0 .and. nscnts .ne. 5) then
C         call iosigh(2,0,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
C      .    lshft(3),ifis)
C         endif
      endif
      print *,'nscnts=',nscnts
c      stop 'xxxxxxxxxx'
      lshft=0
      call mpibc1(i,1,2,.false.,'rdsigm','nsp')
      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
c
      if (i == nsp) then
        lsplts = .false.
      elseif (i > nsp) then
        call rx(
     .  'rdsigm: sigm file spin polarized but calculation is not')
      else
        lsplts = .true.
      endif
      if (lsplts) call info0(30,0,0,
     .'%9f(warning) sigm file not spin pol .. splitting spins')

C     Check for consistency in remaining ndimh, nqp
C     and set file sigma file pointer past header
      if (procid .eq. master) then
        rewind ifis
c        if (lfbzin) nqp = nk1*nk2*nk3
        call iosigh(3,nscnts,i,ndimh,nk1,nk2,nk3,nqp,
     .  lshft(1),lshft(2),lshft(3),ifis)
      endif
      call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
      call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
      call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
      call info5(30,0,0,'%9fsigm file has %i irreducible QP: nk ='//
     .' ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)

C     Sanity checks

      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
ckino Dec.22.2011:         rv_p_opos => slat%rv_p_opos 



      nsgrp=slat%nsgrp
ckino Dec.22.2011:         iv_p_oistab => slat%iv_p_oistab 

ckino Dec.22.2011:         rv_p_og => slat%rv_p_osymgr 

ckino Dec.22.2011:         rv_p_oag => slat%rv_p_oag 


C     Get leading dimension of hrs
c      call offsHp(w(oiprmb),1,nbas,2,0,ndimh,0,ndhrs)
C     Setup for FFT
      call fftz30(nk1,nk2,nk3,k1,k2,k3) !takao comment out because k1 k2 k3 are supplied at the begining of routine.
      if (nk1.ne.k1 .or. nk2.ne.k2 .or. nk3.ne.k3) call rx(
     .'rdsigm: not ready for FFT w/ dimensions ne no div.')
C ... Number of group operations for input file sigma
c      if (mod(lssym,2) .eq. 0) nsgrps = nsgrp
c      if (mod(lssym,2) .eq. 1) nsgrps = 1
c      if (lfbzin) nsgrps = 1
      nsgrps = nsgrp

C ... Make is,ifac,qb,qlat,qoff
      do i = 1, 3
        llshft(i) = lshft(i) .ne. 0
      enddo
      call pshpr(0)
      call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
      call poppr

c.takao qmesh
      do i1=1,nk1
        do i2=1,nk2
          do i3=1,nk3
            qsmesh(:,i1,i2,i3) = (i1*ifac(1)-1)*qb(:,1) +
     &      (i2*ifac(2)-1)*qb(:,2) +
     &      (i3*ifac(3)-1)*qb(:,3)
          enddo
        enddo
      enddo

c      do  i = 1, 3
c        qoff(i) = qk(i,1,1,1) + qoffi(i)
c      enddo


C ... Setup for r.s. hamiltonian: allocate and create iax,ntab; allocate hrs
c      i = 1000 + 10*hreal + 0
C      print *, '!!'; i = 1000*0 + 10*hreal + 0
c      if (mod(lssym,4) .ge. 2) i = 10*hreal + 0
c      if (ltrans) i = i + 10000
c      i=11000
c      print *,' go to hft2rs i=',i
c      call hft2rs(i,nk1,nk2,nk3,k1,k2,k3,qoff,1,nsp,nbas,w(og),w(oag),
c     .  nsgrp,rsrnge,w(ooffH),1,nbas,ndimh,ndimh,w,plat,w(opos),ontabs,
c     .  oiaxs,ndhrs,w)
c      nttabs = ival(w(ontabs),nbas+1)
C     Allocate memory for hrs
c      if (hreal .eq. 1) then
c        call defdr(ohrs,-ndhrs**2*nttabs*nsp)
c      else
c      call defdc(ohrs,-ndhrs**2*nttabs*nsp)
c      endif

C --- Get list of irreducible k-points, and ipq and gstar arrays ---
      mxkp = nk1*nk2*nk3
ckino Dec.22.2011:         allocate(rv_p_oqsig(abs(3*mxkp)))
       allocate(sham%rv_p_oqsig(abs(3*mxkp)))
ckino Dec.22.2011:         if (3*mxkp<0) rv_p_oqsig(:)=0.0d0
       if (3*mxkp<0) sham%rv_p_oqsig(:)=0.0d0

      allocate(qp_rv(3*mxkp))

      allocate(gstar_iv(abs(-mxkp-1)))
      if (-mxkp-1<0) gstar_iv(:)=0

      gstar_iv = - 2

      allocate(ipq_iv(mxkp))

      allocate(wgt_rv(abs(-mxkp)))
      if (-mxkp<0) wgt_rv(:)=0.0d0

      call info(20,1,0,
     .' q-points in full BZ where sigma calculable ...',0,0)
        call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy 
ckino Dec.22.2011:       .  , 0 , ipq_iv , rv_p_oqsig , wgt_rv , nqsig , mxkp , 0 , 0 ) 
     .  , 0 , ipq_iv , sham%rv_p_oqsig , wgt_rv , nqsig , mxkp , 0 , 0 ) 




      sham%nqsig=nqsig
ckino Dec.22.2011:         sham%rv_p_oqsig => rv_p_oqsig 


      call dpzero ( wgt_rv , mxkp )

      call info2(20,0,0,
     .' Irr. qp for which sigma is calculated ...',0,0)
ckino Dec.22.2011:         call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og 
       call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , slat%rv_p_osymgr 
     . , nsgrps , ipq_iv , qp_rv , wgt_rv , nqps , mxkp , gstar_iv , 
     . 0 ) 


Ckino isanrg is logical function,       call isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
      l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
C     call prmx('q for which sigma is read',w(oqp),3,3,nqp)
c      call info2(30,0,0,
c     .  '%?#(n==1)#%1fFile sigm has k-offset=%3:2g '//
c     .  'relative to the above#%j#'//
c     .  ' ',isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi)

C ... titus
CC     List of sites
C      call getoffH(w(ooffH),offH,n0H,nkap0,nbas)
C      ib1 = 1
C      ib2 = nbas
C      lidim = offH(4,1,nbas+1)
C
C      print*,'lidim',lidim
C
C      do 100 ib = ib1, ib2
C      rewind ifis
C      call iosigh(0,nscnts,nsp1,ndimh,nk1,nk2,nk3,nqp,
C     .  lshft(1),lshft(2),lshft(3),ifis)
C      nlmaa = offH(4,1,ib+1) - offH(4,1,ib)
C      print*,'nlmaa',nlmaa
C      print*,'offH(4,1,ib)',offH(4,1,ib)
C      print*,'offH(4,1,ib+1)',offH(4,1,ib+1)
C ... titus


C --- Generate hrs = sigma(T) from file sigma(k) ---
C ... Loop over spins, keeping sig in full BZ for only one spin
      if (procid .eq. master) then

        do  isp = 1, nsp
C     If sigma file not spin polarized, use sigma from spin 1
          if (isp .eq. 2 .and. lsplts) then
            call iosigh(3,0,1,ndimh,nk1,nk2,nk3,nqp,
     .      lshft(1),lshft(2),lshft(3),ifis)
          endif
C     If to mix in some other sigma file, set up ifis2
          if (alf(2) .ne. 0) then
            if (isp .eq. 1 .or. isp .eq. 2 .and. lsplts) then
              ifis2 = fopna('sigm1',-1,4)
              call iosigh(3,0,nsp,ndimh,nk1,nk2,nk3,nqp,
     .        lshft(1),lshft(2),lshft(3),ifis2)
            endif
          endif

C ... File sigma in irr BZ -> sigma in full BZ by rotations
c      call defcc(osfz,k1*k2*k3*ndimh**2)

C ... titus
C      call defcc(osfz,k1*k2*k3*ndimh*nlmaa)
C      call defcc(osfz1,k1*k2*k3*ndimh*nlmaa)
C ... titus

          allocate(wk_zv(ndimh**2))
Changenglob          nl = nglob('nl')
          nl = globalvariables%nl
          nspc = 1
          allocate(sigm_zv(ndimh**2))

          do  iq1 = 1, nqps
            read(ifis) qp,  eseavr !june2009  add eseavr
            if(debugmode()>0) print *,' uuu2 iq1 eseavr=',iq1,eseavr
            call dpscop(qp,tmp,3,1,1,1d0)
            call dpsadd ( tmp , qp_rv , 3 , 1 , 3 * iq1 - 2 , - 1d0 )

c        call dpsadd(tmp,qoffi,3,1,1,-1d0)

            !tmp= qp - w((oqp(3*(iq1-1)+1:3*(iq1-1)+3))+qoffi)
            if (.not. latvec(1,tolq,plat,tmp)) then
              call dpscop ( qp_rv , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 )

              if (lssym .ge. 4) then
                call info5(20,0,0,' rdsigm (warning) qp %i : '//
     .          'expected%3;8,4D read%3;8,4D',iq1,tmp,qp,0,0)
                call dpscop(tmp,qp,3,1,1,1d0)
              else
                print 456, iq1,tmp,qp
  456           format(' rdsigm: error on reading qp no',i4/
     .          ' Expected qp=',3f12.6/
     .          '     File qp=',3f12.6)
                call rx(' incompatible q-mesh')
              endif
            endif
C       File read sigma
c        print *,' nscnts =',nscnts
c        if (nscnts == 5) call dpdump(tmp,1,ifis)
            call dpdump ( sigm_zv , ndimh * * 2 * 2 , ifis )


C        if (isp .eq. 2) then
C          print *, iq1,isp
C         call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
C        endif
C
cccccccccccccccccccccccccccccccccccccc
ctakaox !see also bndfp.F and hambls.F
c        if(iq1==1) then
c          print *,'reading ovlmat'
c          ifio = fopna('ovlmat',-1,4)
c        endif
c        call dpdump(w(oss),ndimh**2*2,ifio)
c        print *,'read ovlmat for iq1=',iq1
c        call divctakao(w(osigm),w(oss),ndimh)
ccccccccccccccccccccccccccccccccccccc

c ... scalced sigm. alf(1:2)
            if (alf(1) .ne. 1d0) then
              call dscal ( ndimh * * 2 * 2 , alf ( 1 ) , sigm_zv , 1 )

            endif
c$$$        if (alf(2) .ne. 0d0) then
c$$$          read(ifis2) qp
c$$$          call dpscop(qp,tmp,3,1,1,1d0)
c$$$          call dpsadd(tmp,w(oqp),3,1,3*iq1-2,-1d0)
c$$$c          call dpsadd(tmp,qoffi,3,1,1,-1d0)
c$$$          if (abs(tmp(1))+abs(tmp(2))+abs(tmp(3)) .gt. tolq) then
c$$$            call dpscop(w(oqp),tmp,3,3*iq1-2,1,1d0)
c$$$            if (lssym .ge. 4) then
c$$$              call info5(20,0,0,' rdsigm (warning) qp %i : '//
c$$$     .          'expected%3;8,4D read%3;8,4D (file sigm1)',iq1,tmp,qp,
c$$$     .          0,0)
c$$$            else
c$$$              call rx(' file sigm1 has incompatible irr mesh')
c$$$            endif
c$$$          endif
c$$$C         call dpsadd(qp,w(oqp),3,1,3*iq1-2,1d0)
c$$$          call defcc(osigm2,ndimh**2)
c$$$          call dpdump(w(osigm2),ndimh**2*2,ifis2)
c$$$          call daxpy(ndimh**2*2,alf(2),w(osigm2),1,w(osigm),1)
c$$$          call rlse(osigm2)
c$$$        endif
c ... rotate sigm. w(osigm) at qp ---> w(ofbz) at stars of qp.
C       call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
            write(6,"(a,13f13.5)")' Goto hamfb3k  xxx input qp=', qp
            ikt = getikt(qp)
            napw_in= napwk(ikt)

            i_copy_size=size(slat%qlat)
            call dcopy(i_copy_size,slat%qlat,1,qlat,1)

            call hamfb3k ( qp , iq1 , nk1 , nk2 , nk3 , k1 , k2 , k3 , ipq_iv 
     .      , napw_in , ndimh , ndimh , ndimh , qb , plat , qlat , ldim , 
     .      ifac , gstar_iv , sigm_zv , sfz ( 1 , 1 , 1 , 1 , 1 , isp 
     .      ) )

            if(debugmode()>0) write(6,"(a,3f13.5)")'end of hamfbk3'
          enddo
          deallocate(sigm_zv)
          deallocate(wk_zv)
C ... End loop over spins
        enddo
      endif ! exit of MPI master loop
      call tcx('rdsigm2')
      if (allocated(wgt_rv)) deallocate(wgt_rv)
      if (allocated(ipq_iv)) deallocate(ipq_iv)
      if (allocated(qp_rv)) deallocate(qp_rv)
      if (allocated(gstar_iv)) deallocate(gstar_iv)

      end subroutine rdsigm2




