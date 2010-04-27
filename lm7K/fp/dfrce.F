Cgetarg...info...           structure ['dfrce', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['dfrce', 'lat', 'slat']
Cgetarg...info...           structure ['dfrce', 'spec', 'sspec']
Cgetarg...info...           structure ['dfrce', 'site', 'ssite']
Cgetarg...info...           structure ['pvdf1', 'lat', 'slat']
Cgetarg...info...           structure ['pvdf1', 'site', 'ssite']
Cgetarg...info...           structure ['pvdf1', 'spec', 'sspec']
Cgetarg...info...           structure ['pvdf4', 'lat', 'slat']
Cgetarg...info...           structure ['pvdf4', 'site', 'ssite']
Cgetarg...info...           structure ['pvdf4', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['pvdf1', 'm_struc_def']
Cgetarg...info...           use_to_add ['dfrce', 'm_struc_def']
Cgetarg...info...           use_to_add ['pvdf4', 'm_struc_def']
      subroutine dfrce(ssite,sspec,slat,sctrl,k1,k2,k3,nvl,
     .  orhoat,orhat1,elind,qmom,smrho,smrout,dfh)
       use m_struc_def  !Cgetarg
#if MPE 
      use m_events !Event numbers: include "events.ins"
#endif

C-Correction to force theorem, Harris functional
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: pvdf4 pvdf2 rhomom pvdf1 smvxcm
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl z p pz lmxa a nr rmt nxi exi chfa rsmfa rg
Ci     Stored:    *
Ci     Passed to: pvdf4 pvdf2 rhomom pvdf1 gtpcor corprm smvxcm
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc ng ogv okv vol alat plat qlat
Ci     Stored:    *
Ci     Passed to: pvdf4 pvdf2 pvdf1 smvxcm
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: lfrce
Ci     Stored:    *
Ci     Passed to: *
Ci   k1..3 :dimensions smrho
Ci   nvl   :sum of local (lmxl+1)**2, lmxl = density l-cutoff
Ci   orhoat:vector of offsets containing site density
Ci   orhat1:pointer to local densities
Ci   elind :Lindhard parameter, used for Lindhard screening
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   smrho :smooth density on uniform mesh
Ci   smrho :smooth (input) density that generated the hamiltonian
Cio  smrout:smooth (output) density that the hamiltonian generated
Co Outputs
Co   dfh   :correction to the HF force
Cl Local variables
Cl    job  :describes which ansatz for charge shift is used for correction
Cl         :<=0  do not calculate correction to force
Cl         :  1  shift in free-atom density
Cl         :  2  shift in core+nuclear density
Cl         :+10  to screen the rigid shift by the Lindhard function
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   18 Dec 03 adapted to modified smvxc
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   17 Sep 01 Adapted for local orbitals
Cu   21 Jun 00 spin polarized
Cu   18 Jun 98 adapted from nfp dfrce.f
Cu   16 Jun 98 MvS parallelized for SGI
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI | MPIK
      include "mpif.h"
#if MPE
      include "mpef.h"
#endif 
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog,cmdopt
#endif
      integer nvl,k1,k2,k3,orhoat(3,1),orhat1(3,1)
Cgetarg       double precision ssite(1),sspec(1),slat(1),sctrl(1),dfh(3,1),
Cgetarg      .  qmom(1)
       real(8):: dfh(3,1) , qmom(1) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_ctrl)::sctrl

      double complex smrho(k1,k2,k3,*),smrout(k1,k2,k3,*)
C Local variables
      integer nbas,job,n1,n2,n3,ng,iprint,nsp,nglob,ib,is,lmxl,iv0,nlm,
     .  ip,m,i,ngabc(3),ltop,nspec,nlmtop,igets,igetss,nn,lgunit,stdo
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer ogv,okv,ocdv,oceps,ocdn,ocnomi,osmro,odvxc,ovxcp,ovxcm,
     .  ocdvx,oqmout,ocvin,ocdn0,oyl,og2,og,oiv,ocs,osn,owk1,owk2,owk3
      double precision vol,plat(3,3),qlat(3,3),alat,vsum,pi,tpiba,elind,
     .  fes1(3),fes2(3),fxc(3),c,avgdf(3)
      integer npmx,mp
      parameter (npmx=32)
      integer oicdn(npmx),oicdn0(npmx),oicdv(npmx),oics(npmx),oisn(npmx)
#if MPI | MPIK
      integer, dimension(:), allocatable :: bproc, iiv0
      integer pid,jb
#endif
      integer::ibini,ibend
C ... for SGI parallel
C$    integer nbmx,mp_my_threadnum,mp_numthreads
C$     parameter (nbmx=512)
C$    integer iiv0(nbmx)
C     double precision dqsmo,dval
      character*40 strn
C ... Heap
      integer w(1)
      common /w/ w

Cgetarg       job = igets('ctrl lfrce',sctrl)
       job = int(sctrl%lfrce) 

      if (job .le. 0) return
      call tcn('dfrce')

#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      if (mlog) then
        do  pid = 0, numprocs-1
          call MPI_BCAST(shortname(pid),10,MPI_CHARACTER,pid,
     .                   MPI_COMM_WORLD,ierr) 
          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,
     .                   MPI_COMM_WORLD,ierr) 
        enddo
      endif
#endif

C --- Setup ---
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat nabc ng ogv okv vol',slat,ngabc,ng,ogv,okv,vol)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       ng=slat%ng
       ogv=slat%ogv
       okv=slat%okv
       vol=slat%vol

Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: qlat(3,3) , real(8) :: slat%qlat 9
Cgetarg       call upack('lat alat plat qlat',slat,alat,plat,qlat,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       i_copy_size=size(slat%qlat) 
       call dcopy(i_copy_size,slat%qlat,1,qlat,1) 

      c = 1000
      stdo = lgunit(1)
      nsp  = nglob('nsp')
      nbas = nglob('nbas')
      nspec= nglob('nspec')
      nn   = k1*k2*k3

C ... Arrays needed for pvdf1
      call defcc(oceps, ng)
      call defcc(ocnomi,ng)
      call defcc(ocvin, ng)
      call defcc(ocdvx, ng*nsp)

C ... Set up for vectorized Y_lm and gaussians
      ltop = 0
      do   is = 1, nspec
Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        ltop = max0(ltop,lmxl)
      enddo
      nlmtop = (ltop+1)**2
      call defrr(oyl,  ng*nlmtop)
      call defrr(og2,  ng)
      call defrr(og,   ng*3)
      call suylg(ltop,alat,ng,w(ogv),w(og),w(og2),w(oyl))
      call rlse(og)
      call defi(oiv, ng*3)
      call suphs0(plat,ng,w(ogv),w(oiv))

C --- Make ves(rhoin,q) ---
      call defcc(osmro, nn)
      call defrr(ocs, ng)
      call defrr(osn, ng)
      call dpcopy(smrho,w(osmro),1,2*nn,1d0)
      if (nsp.eq.2) call dpadd(w(osmro),smrho(1,1,1,2),1,2*nn,1d0)
      call fftz3(w(osmro),n1,n2,n3,k1,k2,k3,1,0,-1)
      call gvgetf(ng,1,w(okv),k1,k2,k3,w(osmro),w(ocvin))
      call pvdf4(ssite,sspec,slat,qmom,ng,
     .  w(og2),w(oyl),w(ocs),w(osn),w(oiv),qlat,w(ocvin))
C     call gvputf(ng,1,w(okv),k1,k2,k3,w(ocvin),w(osmro))
C     call zprm3('ves(input rho)',0,w(osmro),k1,k2,k3)
      call rlse(osmro)

C --- Make dVxc(in)/dn ---
      call defcc(odvxc, nn*nsp)
      call defcc(osmro, nn*nsp)
      call defcc(ovxcp, nn*nsp)
      call defcc(ovxcm, nn*nsp)
      call defcc(owk1,  nn*nsp)
      call defcc(owk2,  nn*nsp)
      call defcc(owk3,  nn*nsp)
      call dpcopy(smrho,w(osmro),1,2*nn*nsp,1d0)
      call pvdf2(nbas,nsp,ssite,sspec,slat,n1,n2,n3,k1,k2,k3,w(osmro),
     .  w(ovxcp),w(ovxcm),w(owk1),w(owk2),w(owk3),w(odvxc))
      call rlse(ovxcp)

C --- cdvx = FFT ((n0_out-n0_in) dVxc/dn) ---
C     Use total n0_out-n0_in but keep vxc spin polarized
      call dpzero(w(osmro),2*nn)
      do  i = 1, nsp
        call dpadd(w(osmro),smrout(1,1,1,i),1,2*nn,1d0)
        call dpadd(w(osmro),smrho(1,1,1,i),1,2*nn,-1d0)
      enddo
C     call zprm3('drho',0,w(osmro),k1,k2,k3)
      call pvdf3(n1,n2,n3,k1,k2,k3,nsp,w(osmro),w(odvxc))
      call fftz3(w(odvxc),n1,n2,n3,k1,k2,k3,nsp,0,-1)
      call gvgetf(ng,nsp,w(okv),k1,k2,k3,w(odvxc),w(ocdvx))

C --- Cnomi = (n0_out(q) - n0_in(q)) ---
      call fftz3(w(osmro),n1,n2,n3,k1,k2,k3,1,0,-1)
c     call zprm3('FFT smrout-smrin',w(osmro),k1,k2,k3*nsp)
      call gvgetf(ng,1,w(okv),k1,k2,k3,w(osmro),w(ocnomi))
C     dqsmo = vol*dval(w(ocnomi),1)

      call rlse(odvxc)

C ... Debugging slot smrho(out) for out-in
C      print *, '*** debugging ... subs smrout for out-in'
C      call dpcopy(smrout,w(osmro),1,2*nn,1d0)
C      call fftz3(w(osmro),n1,n2,n3,k1,k2,k3,1,0,-1)
C      call gvgetf(ng,1,w(okv),k1,k2,k3,w(osmro),w(ocnomi))
C      call zprm3('rho-out(q)',w(osmro),k1,k2,k3)

C --- Multipole moments of the output density ---
      call defrr (oqmout,   nvl)
      call pshpr(0)
      call rhomom(nbas,ssite,sspec,orhat1,w(oqmout),vsum)
      call poppr
      call dpadd(w(oqmout),qmom,1,nvl,-1d0)

C --- Lindhard dielectric function ---
      if (job .gt. 10) then
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
C        call lindxx(123,n1,n2,n3,k1,k2,k3,ng,w(okv),w(oceps),w(ogv),
C     .    tpiba,elind,w,w,w,w)
        call lindsc(3,ng,w(ogv),tpiba,elind,w(oceps))
      endif

C --- For each site, get correction to force ---
#if ! (MPI | MPIK)
      if (iprint() .ge. 30) then
        strn = 'shift in free-atom density'
        if (job .eq. 11) strn = 'screened shift in free-atom density'
        if (job .eq. 12) strn = 'screened shift in core+nuclear density'
        write(stdo,201) strn
      endif
#endif
  201 format(/' Harris correction to forces: ',a/
     .  '  ib',9x,'delta-n dVes',13x,'delta-n dVxc',15x,'total')
C  201   format(/' Harris correction to forces:'/
C     .    '  ib',11x,'dn0 dVes',15x,'dnloc dVes',15x,
C     .    'dn0 dVxc',15x,'total')

C ... Setup array iiv0 = (vector of iv0 for parallel); allocate work arrays
#if MPI | MPIK
      iv0 = 0
      allocate(iiv0(1:nbas), stat=ierr)
      do  12  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        nlm = (lmxl+1)**2
        iiv0(ib) = iv0
        iv0 = iv0+nlm
   12 continue
#endif 
C$    iv0 = 0
C$    if (nbas .gt. nbmx) call rxi('dfrce: need nbmx at least',nbas)
C$    do  12  ib = 1, nbas
C$      is = igetss('site spec',ib,ssite)
C$      lmxl = igetss('spec lmxl',is,sspec)
C$      nlm = (lmxl+1)**2
C$      iiv0(ib) = iv0
C$      iv0 = iv0+nlm
C$ 12 continue
      mp = 1
C$    mp=mp_numthreads()
C$    if(npmx.lt.mp) call rxi('rsif_q: increase npmx, needed',mp)
      do  14  ip = 1, min(nbas,mp)
        call defcc(oicdn(ip),  ng)
        call defcc(oicdn0(ip), ng*nsp)
        call defcc(oicdv(ip),  ng)
        call defrr(oics(ip), ng)
        call defrr(oisn(ip), ng)
   14 continue

C ... Estimate shift in density for each site, and corresponding force
      ip = 1
      iv0 = 0
C$DOACROSS LOCAL (ib,is,lmxl,nlm,i,m,iv0,ip,ocdn,ocdn0,ocdv,
C$&               ocs,osn,fes1,fes2,fxc)
C$&        SHARED (job,ng,iiv0,ocnomi,ocvin,ocdvx)
C$&        MP_SCHEDTYPE=RUNTIME
#if MPI | MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_START_DFRCE,procid,"dfrce")
#endif 
      allocate (bproc(0:numprocs), stat=ierr)
      call dstrbp(nbas,numprocs,1,bproc(0))
      ibini = bproc(procid)
      ibend = bproc(procid+1)-1
#else
      ibini=1
      ibend=nbas
#endif 

      do  ib = ibini, ibend
#if MPI | MPIK
        if (mlog .and. ib .eq. bproc(procid)) then
          call gettime(datim)
          call awrit4(' dfrce '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting atoms %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        iv0 = iiv0(ib)
#endif 

Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        if (lmxl .eq. -1) goto 20

        nlm = (lmxl+1)**2
C$      iv0 = iiv0(ib)
C$      ip  = mp_my_threadnum()+1
        ocs = oics(ip)
        osn = oisn(ip)
        ocdn = oicdn(ip)
        ocdn0 = oicdn0(ip)
        ocdv = oicdv(ip)

        call pvdf1(job,ssite,sspec,slat,nsp,ib,iv0,qmom,w(oqmout),
     .    ng,w(ogv),w(og2),w(oyl),w(ocs),w(osn),w(oiv),qlat,
     .    0,w(ocnomi),w(oceps),w(ocdn0),w(ocdn),w(ocdv),
     .    w(ocdvx),w(ocvin),orhoat(1,ib),fes1,fes2,fxc)

C        call gvputf(ng,1,w(okv),k1,k2,k3,w(ocdn),smpot)
C        call fftz3(smpot,n1,n2,n3,k1,k2,k3,1,0,1)
C        call zprm3('dn',smpot,k1,k2,k3)
C
        do  i = 1, 3
          dfh(i,ib) = -(fes1(i) + fes2(i) + fxc(i))
        enddo

#if ! (MPI | MPIK)
        if (iprint() .ge. 30)
     .    write(stdo,200) ib,(c*(fes1(m)+fes2(m)),m=1,3),
     .    (c*fxc(m),m=1,3),(c*dfh(m,ib),m=1,3)
  200   format(i4,3f8.2,1x,3f8.2,1x,3f8.2:1x,3f8.2)
#endif

        iv0 = iv0+nlm
   20 continue
      enddo
#if MPI | MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_DFRCE,procid,"dfrce")
      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
#endif 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
      ierr = MPE_LOG_EVENT(EVENT_START_BCAST,procid,"broadcast")
#endif 
      do  pid = 0, numprocs-1
        ib = bproc(pid)
        jb = bproc(pid+1) - ib
        call MPI_BCAST(dfh(1,ib),3*jb,MPI_DOUBLE_PRECISION,pid,
     .                 MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit6(' dfrce '//datim//' Process %i of %i on '
     .          //shortname(procid)(1:namelen(procid))//
     .          ' bcast dfh(1-3,%i-%i) %i d.p. numbers'//
     .          ' from process %i on '
     .          //shortname(pid)(1:namelen(pid)),' ',
     .          256,lgunit(3),procid,numprocs,bproc(pid),bproc(pid+1)-1,
     .          3*jb,pid)
        endif
      enddo
      deallocate (bproc, stat=ierr)
      deallocate (iiv0, stat=ierr)
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_BCAST,procid,"broadcast")
#endif 
#endif 
      call dpzero(avgdf,3)
      do  23  ib = 1, nbas
      do  23  i = 1, 3
   23 avgdf(i) = avgdf(i) + dfh(i,ib)/nbas

C ... MPI printout
#if MPI | MPIK
      if (iprint() .ge. 30) then
        strn = 'shift in free-atom density'
        if (job .eq. 11) strn = 'screened shift in free-atom density'
        if (job .eq. 12) strn = 'screened shift in core+nuclear density'
        write(stdo,201) strn
        do  26  ib = 1, nbas
   26   write(stdo,202) ib, (c*dfh(i,ib),i=1,3)
  202   format(i4,50x,3f8.2)
      endif
#endif

C ... Shift all forces to make avg correction zero
      do  24  ib = 1, nbas
      do  24  i = 1, 3
   24 dfh(i,ib) = dfh(i,ib) - avgdf(i)
      if (iprint() .ge. 30) write(stdo,331) (c*avgdf(m),m=1,3)
  331 format(' shift forces to make zero average correction:',8x,3f8.2)
      call rlse(oceps)
      call tcx('dfrce')
Cgetarg       end
       end subroutine dfrce 


      subroutine pvdf1(job,ssite,sspec,slat,nsp,ib,iv0,qmom,qmout,ng,gv,
     .  g2,yl,cs,sn,iv,qlat,kmax,cnomin,ceps,cdn0,cdn,cdv,cdvxc,cvin,
     .  orhoat,fes1,fes2,fxc)
       
       use m_struc_def  !Cgetarg

C- Estimate shift in local density for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite,sspec,slat
Ci   ng,gv,kmax
Ci   orhoat
Ci   job: 1  shift in free-atom density
Ci        2  shift in core+nuclear density
Ci      +10  to screen the rigid shift by the response function
Ci   ib      which site is being shifted
Ci   iv0     offset to qmom
Ci   qmom,qmout moments of input and output densities
Ci   cnomin  difference betw. smoothed output and input density n0
Ci   cvin    electrostatic potential of input density Ves[n0~_in]
Ci   ceps    response function
Ci   cdvxc   dVxc/dn (nout-nin)
Co Outputs
Co   cdn0:   Job 1:  shift in valence part of the free atom density
Co           Job 12: shift in atom density (1/eps - 1)
Co   cdn:    Job 1:  dn^(u) where dn is the unscreened shift in
Co           in the free-atom density.
Co           Job 12: dn^(u) 1/eps where dn is unscreened shift in
Co           the charge density.  Local density approximated
Co   NB:     In all cases, the local part of density is approximated
Co           by a gaussian of the equivalent multipole moment.
Co   cdv:    shift in the electrostatic potential
Co   fes1,fes2,fxc
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,nsp,iv0,kmax,ib,job,orhoat(3),iv(ng,3)
Cgetarg       double precision ssite(1),sspec(1),slat(1),qmom(1),qmout(1),
Cgetarg      .  gv(ng,3),tau(3),fes1(3),fes2(3),fxc(3),g2(ng),yl(ng,1),
Cgetarg      .  cs(ng),sn(ng),qlat(3,3)
       real(8):: qmom(1) , qmout(1) , gv(ng,3) , tau(3) , fes1(3) , 
     . fes2(3) , fxc(3) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex cdn0(ng,nsp),cdn(ng),cdv(ng),ceps(ng),
     .  cnomin(ng),cdvxc(ng,nsp),cvin(ng)
C ... Local parameters
      integer ig,ilm,l,lmxl,m,nlm,nlmx,k,is,jv0,jb,nbas,nglob,js,ll,n0,
     .  nrmx
      parameter (nlmx=64, nrmx=1501, n0=10)
      integer lmxa,nr,nxi,ie,ixi,job0,kcor,lcor,lfoc,i,
     .  ngabc(3),n1,n2,n3,nlml
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision pi,alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg,
     .  vol,y0,z,v(3),df(0:20),feso(3),qcor(2),gpot0(nlmx,3),fesdn(3),
     .  fesgg(3),pnu(n0,2),pnz(n0,2),a,rmt,qloc,rsmfa,exi(n0),hfc(n0,2),
     .  qfat,gam,qall,qc,qval,qg,e,aa,q0(3),sum
      double precision rwgt(nrmx),cc,gamf,cfoc,cvol
C     parameter (k0=3)
C     double complex gkl(0:k0,nlmx)
      double complex tpia,cxx,phase,gc0,xc0,cof(nlmx)
C ... Heap
      integer w(1)
      common /w/ w
      data q0 /0d0,0d0,0d0/

      call tcn('pvdf1')
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat nabc alat vol',slat,ngabc,alat,vol,0,0)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       alat=slat%alat
       vol=slat%vol

      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      nbas = nglob('nbas')
C     tpiba = 2*pi/alat
      tpia = 2*pi*dcmplx(0d0,-1d0)/alat
      job0 = mod(job,10)

      call dpzero(cdn,2*ng)
      call dpzero(cdn0,2*ng*nsp)
      cdv(1) = 0d0
      call dpzero(fes1,3)
      call dpzero(fxc,3)
      call dpzero(fesdn,3)
      call dpzero(gpot0,nlmx*3)
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg       call upack('site spec pos',ssite,ib,is,tau,0,0)
       
       is=ssite(ib)%spec
       i_copy_size=size(ssite(ib)%pos) 
       call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

      call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)

C --- Unscreened rigid charge density shift, job 1, in cdn0 ---
      if (job0 .eq. 1) then
C       call upack('site spec pos',ssite,ib,is,tau,0,0)
c        call upack('spec z p pz lmxa lmxl',sspec,is,z,pnu,pnz,lmxa)
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: pnu(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: pnz(n0,2) , real(8) :: sspec%pz 20
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg         call upack('spec z p pz lmxa',sspec,is,z,pnu,pnz,lmxa)
         
         z=sspec(is)%z
         i_copy_size=size(sspec(is)%p) 
         call dcopy(i_copy_size,sspec(is)%p,1,pnu,1) 
         i_copy_size=size(sspec(is)%pz) 
         call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1) 
         lmxa=sspec(is)%lmxa

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg         call upack('spec a nr rmt lmxl',sspec,is,a,nr,rmt,lmxl)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt
         lmxl=sspec(is)%lmxl

Cgetarg...info...           integer :: nxi , integer(8) :: sspec%nxi 1
Cgetarg...info...           real(8) :: exi(n0) , real(8) :: sspec%exi 10
Cgetarg...info...           real(8) :: hfc(n0,2) , real(8) :: sspec%chfa 20
Cgetarg...info...           real(8) :: rsmfa , real(8) :: sspec%rsmfa 1
Cgetarg         call upack('spec nxi exi chfa rsmfa',sspec,is,nxi,exi,hfc,rsmfa)
         
         nxi=sspec(is)%nxi
         i_copy_size=size(sspec(is)%exi) 
         call dcopy(i_copy_size,sspec(is)%exi,1,exi,1) 
         i_copy_size=size(sspec(is)%chfa) 
         call dcopy(i_copy_size,sspec(is)%chfa,1,hfc,1) 
         rsmfa=sspec(is)%rsmfa

        call gtpcor(sspec,is,kcor,lcor,qcor)
        if (nr .gt. nrmx) call rx('dfrce: nr gt nrmx')
        call radwgt(rmt,a,nr,rwgt)
        nlml = (lmxl+1)**2
        call radsum(nr,nr,nlml,nsp,rwgt,w(orhoat(1)),qloc)
        call radsum(nr,nr,nlml,nsp,rwgt,w(orhoat(2)),sum)
        qloc = (qloc-sum)/y0
        qfat = 0d0
        do  12  i  = 1, nsp
        do  12  ie = 1, nxi
          gam  = 0.25d0*rsmfa**2
          qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
          qfat = qfat + hfc(ie,i)*qall
   12   continue
        call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
        qg = qval+qsc-qfat-qloc
C   ... Shift in free atom density
        do  14    i = 1, nsp
        do  14  ixi = 1, nxi
        e = exi(ixi)
C       do  14  ig = 1, ng
C         v(1) = gv(ig,1)*tpiba
C         v(2) = gv(ig,2)*tpiba
C         v(3) = gv(ig,3)*tpiba
C         v2 = v(1)**2+v(2)**2+v(3)**2
C         aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
C         scalp = -alat*(tau(1)*v(1)+tau(2)*v(2)+tau(3)*v(3))
C         phase = dcmplx(dcos(scalp),dsin(scalp))
C         cdn0(ig) = cdn0(ig) + hfc(ixi,i)*aa*phase*y0/vol
C     ... Vectorized version
          cc = -4d0*pi*hfc(ixi,i)*y0/vol
          do  15  ig = 1, ng
          aa = cc*dexp(gam*(e-g2(ig)))/(e-g2(ig))
   15     cdn0(ig,i) = cdn0(ig,i) + aa*dcmplx(cs(ig),sn(ig))
   14   continue

C   ... Add gaussian to conserve local charge
C     ... Add gaussian to conserve local charge.  If density corresponds
C         to the free-atom density, qfat+qloc = qval+qsc; then qg=0

C       do  16  ig = 1, ng
C         v(1) = gv(ig,1)*tpiba
C         v(2) = gv(ig,2)*tpiba
C         v(3) = gv(ig,3)*tpiba
C         v2 = v(1)**2+v(2)**2+v(3)**2
C         scalp = -alat*(tau(1)*v(1)+tau(2)*v(2)+tau(3)*v(3))
C         phase = dcmplx(dcos(scalp),dsin(scalp))
C         cdn0(ig) = cdn0(ig) + qg*phase*dexp(-gam*v2)/vol
C  16   continue
C   ... Vectorized version
        cc = qg/vol/nsp
        do  16  i = 1, nsp
        do  16  ig = 1, ng
   16   cdn0(ig,i)=cdn0(ig,i)+cc*dcmplx(cs(ig),sn(ig))*dexp(-gam*g2(ig))
      endif

C --- Coefficients defining local valence + core density ---
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg       call upack('site spec pos',ssite,ib,is,tau,0,0)
       
       is=ssite(ib)%spec
       i_copy_size=size(ssite(ib)%pos) 
       call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg       call upack('spec lmxl rg',sspec,is,lmxl,rg,0,0)
       
       lmxl=sspec(is)%lmxl
       rg=sspec(is)%rg

      call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
      nlm = (lmxl+1)**2
      if (nlm .gt. nlmx) call rxi('pvdf1: increase nlmx to',nlm)
      ilm = 0
      cxx = dcmplx(0d0,1d0)
      do  20  l = 0, lmxl
      cxx = cxx*dcmplx(0d0,-1d0)
      do  22  m = -l,l
      ilm = ilm+1
   22 cof(ilm) = cxx*qmom(ilm+iv0)*4*pi/df(2*l+1)
   20 continue
C     cof(1) = cof(1) + 4*pi*y0*(qcorg+qsc-z)
      cof(1) = cof(1) + 4*pi*y0*(qcorg-z)

C --- Shift in n0, ves~ for list of G vectors ---
      gam = 0.25d0*rg*rg
      gamf = 0.25d0*rfoc*rfoc
      cfoc = -4d0*pi*y0*cofh/vol
      cvol = 1d0/vol
      do  30  ig = 2, ng

        v(1) = gv(ig,1)
        v(2) = gv(ig,2)
        v(3) = gv(ig,3)

C   ... Accumulate unscreened smoothed core+nuclear density

CC       Old, serial version
C        call gkl_ft(v,rg,0d0,tau,alat,kmax,nlm,k0,cy,gkl)
C        do  32  ilm = 1, nlm
C   32   cdn(ig) = cdn(ig) + cof(ilm)*gkl(0,ilm)/vol
C
CC       This part for (grad g) ves(in)
C        do  33  k = 1, 3
C        cxx = tpia*v(k)*cvin(ig)
C        do  33  ilm = 1, nlm
C   33   gpot0(ilm,k) = gpot0(ilm,k) + dconjg(cxx)*gkl(0,ilm)

C   ... Vectorized version (absorb -i**l)
        phase = dcmplx(cs(ig),sn(ig))
        gc0 = phase*dexp(-gam*g2(ig))*cvol
        xc0 = dcmplx(0d0,1d0)*dconjg(tpia*cvin(ig))*gc0*vol
        ilm = 0
        do  32  l = 0, lmxl
          xc0 = xc0*dcmplx(0d0,-1d0)
          do  33  m = -l, l
            ilm = ilm+1
            cdn(ig) = cdn(ig) + yl(ig,ilm)*cof(ilm)*gc0
            gpot0(ilm,1) = gpot0(ilm,1) + yl(ig,ilm)*gv(ig,1)*xc0
            gpot0(ilm,2) = gpot0(ilm,2) + yl(ig,ilm)*gv(ig,2)*xc0
            gpot0(ilm,3) = gpot0(ilm,3) + yl(ig,ilm)*gv(ig,3)*xc0
   33     continue
   32   continue

C   ... Accumulate unscreened foca density
        aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
        cdn(ig) = cdn(ig) + cfoc*aa*phase
C       A slow, unvectorized version
C       call hkl_ft(v,rfoc,ceh,tau,alat,kmax,1,k0,cy, gkl)
C       cdn(ig) = cdn(ig) + cofh*gkl(0,1)/vol

C   ... Make the screened shift in input density n0~
C       Job 1: cdn0 = (valence part of) cdn^u ; cdn = cdn^u
        if (job0 .eq. 1) then
          cdn(ig) = cdn(ig) + (cdn0(ig,1) + cdn0(ig,nsp))/(3-nsp)
          if (job .gt. 10) cdn(ig) = cdn(ig) / ceps(ig)
C       Job 12: cdn0 = cdn^u (1/eps - 1); cdn = cdn^s = cdn^u / eps
        elseif (job .eq. 12) then
          do  i = 1, nsp
            cdn0(ig,i) = cdn(ig) * (1/ceps(ig)-1) / nsp
          enddo
          cdn(ig) = cdn(ig) / ceps(ig)
        else
          call rxi('dfrce: nonsensical job',job)
        endif

C   ... Electrostatic potential shift = 1/eps dv [n0~]
C       g2 = tpiba*tpiba*(gv(ig,1)**2+gv(ig,2)**2+gv(ig,3)**2)
        cdv(ig)  = cdn(ig) * (8*pi/g2(ig))

C       fes1 = (n0_out - n0_in) d ves[n0~]
C       fxc  = dVxc/dn (nout-nin) d n0~
        do  36  k = 1, 3
          fes1(k) = fes1(k) + dconjg(cnomin(ig)) * tpia*v(k)*cdv(ig)
          do  i = 1, nsp
          fxc(k)  = fxc(k)  + dconjg(cdvxc(ig,i)) * tpia*v(k)*cdn0(ig,i)
          enddo
C         fesdn(k)= fesdn(k) + dconjg(cvin(ig))  * tpia*v(k)*cdn0(ig,i)
   36   continue

   30 continue

      do  37  k = 1, 3
C     fesdn(k)  = fesdn(k)*vol
      fxc(k)  = fxc(k)*vol
   37 fes1(k) = fes1(k)*vol

C --- Integral of grad g (output-input local charge) ves~ ---
      call dpzero(fesgg,3)
      do  39  k = 1, 3
      do  39  ilm = 1, nlm
        l = ll(ilm)
        gpot0(ilm,k) = gpot0(ilm,k)*4*pi/df(2*l+1)

C       fesgg(k) = fesgg(k) + qmom(iv0+ilm)*gpot0(ilm,k)
        fesgg(k) = fesgg(k) + qmout(iv0+ilm)*gpot0(ilm,k)
   39 continue

C      print 339, 'n0(out-in) * g dves ',fes1
C      print 339, 'd(g) qmom(out-in) ves[n0~]',fesgg
C      print 339, 'n0~(out-in) * dvxc   ',fxc
C  339 format(a,6p,3f8.2)

C --- Integral of dves~ (output-input local charge) for all sites ---
      call dpzero(fes2,3)
      call dpzero(feso,3)
      jv0 = 0
      do  40  jb = 1, nbas
Cgetarg...info...           integer :: js , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg         call upack('site spec pos',ssite,jb,js,tau,0,0)
         
         js=ssite(jb)%spec
         i_copy_size=size(ssite(jb)%pos) 
         call dcopy(i_copy_size,ssite(jb)%pos,1,tau,1) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack('spec lmxl rg',sspec,js,lmxl,rg,0,0)
         
         lmxl=sspec(js)%lmxl
         rg=sspec(js)%rg

        nlm = (lmxl+1)**2

C ... For this jb, mesh density for all G vectors
        if (nlm .gt. nlmx) call rxi('pvdf1: increase nlmx to',nlm)
        call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
        call dpzero(gpot0,nlmx*3)
        gam = 0.25d0*rg*rg
        do  50  ig = 2, ng

C          v(1) = gv(ig,1)
C          v(2) = gv(ig,2)
C          v(3) = gv(ig,3)
C          call gkl_ft(v,rg,0d0,tau,alat,0,nlm,k0,cy,gkl)
C          do  55  k = 1, 3
C          cxx = tpia*v(k)*cdv(ig)
C          do  55  ilm = 1, nlm
C   55     gpot0(ilm,k) = gpot0(ilm,k) + dble(dconjg(cxx)*gkl(0,ilm))

C ... This is the vectorized version
          aa = dexp(-gam*g2(ig))
          gc0 = dcmplx(0d0,1d0)*aa*
     .          dconjg(tpia*cdv(ig))*dcmplx(cs(ig),sn(ig))
          ilm = 0
          do  55  l = 0, lmxl
            gc0 = gc0*dcmplx(0d0,-1d0)
            do  56  m = -l,l
              ilm = ilm+1
              gpot0(ilm,1) = gpot0(ilm,1)+dble(gc0)*yl(ig,ilm)*gv(ig,1)
              gpot0(ilm,2) = gpot0(ilm,2)+dble(gc0)*yl(ig,ilm)*gv(ig,2)
              gpot0(ilm,3) = gpot0(ilm,3)+dble(gc0)*yl(ig,ilm)*gv(ig,3)
   56       continue
   55     continue

C          print 357, ig, cxx, dconjg(cxx)*gkl(0,1)
C  357     format(i4,1p,4e18.8)

   50   continue

C   ... Multiply factors into gpot0, accumulate force
        ilm = 0
        do  60  l = 0, lmxl
          do  62  m = -l, l
          ilm = ilm+1
          do  64  k = 1, 3
          gpot0(ilm,k) = gpot0(ilm,k)*4*pi/df(2*l+1)

          feso(k) = feso(k) + qmom(jv0+ilm)*gpot0(ilm,k)
          fes2(k) = fes2(k) + qmout(jv0+ilm)*gpot0(ilm,k)
   64   continue
   62 continue
   60 continue

        jv0 = jv0+nlm
   40 continue

C      print 339, 'qmom(in) dv         ',feso
C      print 339, 'qmom(out-in) dv     ',fes2

      call dpadd(fes2,fesgg,1,3,-1d0)
      call tcx('pvdf1')

Cgetarg       end
       end subroutine pvdf1 


      subroutine pvdf2(nbas,nsp,ssite,sspec,slat,n1,n2,n3,k1,k2,k3,
     .  smrho,vxcp,vxcm,wk1,wk2,wk3,dvxc)
      use m_struc_def
C- Makes derivative of smoothed xc potential wrt density.
C     implicit none
C ... Passed parameters
      integer nbas,nsp,n1,n2,n3,k1,k2,k3
Cgetarg       double precision ssite(1),sspec(1),slat(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex vxcp(k1,k2,k3,nsp),vxcm(k1,k2,k3,nsp),
     .               dvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp),
     .               wk1(k1,k2,k3,nsp),wk2(k1,k2,k3,nsp),
     .               wk3(k1,k2,k3,nsp)
C ... Local parameters
      integer i1,i2,i3,i,nn
      double precision fac,dmach,f1,f2,f,alfa,dfdr,rrho,dvdr,
     .  rmusm(2),rvmusm(2),rvepsm(2),repsm(2),repsmx(2),repsmc(2),
     .  fcexc0(2),fcex0(2),fcec0(2),fcvxc0(2)

      fac = dmach(1)**(1d0/3d0)
      alfa = 2d0/3d0
      nn = k1*k2*k3
      call pshpr(0)

C ... Add fac (rho+ + rho-)/2 into rho+, rho- for spin pol case,
C     Add fac * rho into rho if not spin polarized
      if (nsp .eq. 1) then
        call dpcopy(smrho,smrho,1,nn*2,1d0+fac)
      else
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rrho = smrho(i1,i2,i3,1) + smrho(i1,i2,i3,2)
          smrho(i1,i2,i3,1) = smrho(i1,i2,i3,1) + rrho*fac/2
          smrho(i1,i2,i3,2) = smrho(i1,i2,i3,2) + rrho*fac/2
        enddo
        enddo
        enddo
      endif

C ... vxcp = vxc (smrho+drho)
      call dpzero(vxcp, nn*2*nsp)
      call dpzero(wk1, nn*2*nsp)
      call dpzero(wk2, nn*2*nsp)
      call dpzero(wk3, nn*2)
      call smvxcm(ssite,sspec,slat,nbas,0,k1,k2,k3,smrho,
     .  vxcp,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... Replace fac*rho with -fac*rho
      if (nsp .eq. 1) then
        call dpcopy(smrho,smrho,1,nn*2,(1d0-fac)/(1d0+fac))
      else
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rrho = (smrho(i1,i2,i3,1) + smrho(i1,i2,i3,2))/(1d0+fac)
          smrho(i1,i2,i3,1) = smrho(i1,i2,i3,1) - rrho*fac
          smrho(i1,i2,i3,2) = smrho(i1,i2,i3,2) - rrho*fac
        enddo
        enddo
        enddo
      endif

C ... vxcm = vxc (smrho-drho)
      call dpzero(vxcm, nn*2*nsp)
      call smvxcm(ssite,sspec,slat,nbas,0,k1,k2,k3,smrho,
     .  vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... Restore rho+, rho-
      if (nsp .eq. 1) then
        call dpcopy(smrho,smrho,1,nn*2,1/(1d0-fac))
      else
        do  i3 = 1, n3
        do  i2 = 1, n2
        do  i1 = 1, n1
          rrho = (smrho(i1,i2,i3,1) + smrho(i1,i2,i3,2))/(1d0-fac)
          smrho(i1,i2,i3,1) = smrho(i1,i2,i3,1) + rrho*fac/2
          smrho(i1,i2,i3,2) = smrho(i1,i2,i3,2) + rrho*fac/2
        enddo
        enddo
        enddo
      endif

C ... Overwrite vxcp with df/drho
      do  i = 1, nsp
        do  i1 = 1, nn
          rrho = (smrho(i1,1,1,1)+smrho(i1,1,1,nsp))/(3-nsp)
          if (rrho .gt. 0) then
            f1 = vxcm(i1,1,1,i)*(rrho*(1-fac))**alfa
            f2 = vxcp(i1,1,1,i)*(rrho*(1+fac))**alfa
            dfdr = (f2-f1)/(2d0*fac*rrho)
            vxcp(i1,1,1,i) = dfdr
          else
            vxcp(i1,1,1,i) = 0
          endif
        enddo
      enddo

C ... vxcm = vxc (smrho)
      call dpzero(vxcm, nn*2*nsp)
      call smvxcm(ssite,sspec,slat,nbas,0,k1,k2,k3,smrho,
     .  vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .  rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... dvxc/drho into dvxc
      do  i = 1, nsp
        do  i1 = 1, nn
          rrho = (smrho(i1,1,1,1)+smrho(i1,1,1,nsp))/(3-nsp)
          if (rrho .gt. 0) then
            f = vxcm(i1,1,1,i) * rrho**alfa
            dvdr = (vxcp(i1,1,1,i) - alfa*f/rrho) / rrho**alfa
            dvxc(i1,1,1,i) = dvdr
          else
            dvxc(i1,1,1,i) = 0
          endif
        enddo
      enddo

C     call zprm3('d vxc / dn',0,dvxc,k1,k2,k3*nsp)
      call poppr

Cgetarg       end
       end subroutine pvdf2 


      subroutine pvdf3(n1,n2,n3,k1,k2,k3,nsp,deln0,dvxc)

C- Overwrites dvxc with (nout-nin)*dvxc
C     implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,nsp
      double complex deln0(k1,k2,k3),dvxc(k1,k2,k3,nsp)
C ... Local parameters
      integer i1,i2,i3,i

      do  10  i  = 1, nsp
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1
   10 dvxc(i1,i2,i3,i) = dvxc(i1,i2,i3,i)*deln0(i1,i2,i3)

C     call zprm3('dvxc/dn * (nout-nin)',0,dvxc,k1,k2,k3*nsp)

Cgetarg       end
       end subroutine pvdf3 


      subroutine pvdf4(ssite,sspec,slat,qmom,ng,g2,yl,cs,sn,iv,qlat,cv)
       
       use m_struc_def  !Cgetarg

C- Makes smoothed ves from smoothed density and qmom, incl nuc. charge
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl rg
Ci     Stored:    *
Ci     Passed to: corprm
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc vol
Ci     Stored:    *
Ci     Passed to: *
Ci   qmom  :multipole moments of on-site densities (rhomom.f)
Ci   ng    :number of G-vectors
Ci   g2    :square of G-vectors
Ci   yl    :spherical harmonics 
Ci   cs    :vector of cosines for the ng vectors
Ci   sn    :vector of sines for the ng vectors
Ci   iv
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Co Outputs
Co   cv    :local gaussian density added to cv
Co         :estatatic potential make from density
Cr Remarks
Cr   Local charge consists of a sum of gaussians that compensate for
Cr   the difference in multipole moments of true and smooth local charge
Cr   and a contribution from the smooth core charge.
Cr     g(qmpol) + g(qcore-z) + h(ncore)
Cr
Cr   Adapted from vesgcm to make strictly FT ves(nloc)
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ng,iv(ng,3)
Cgetarg       double precision ssite(1),sspec(1),qmom(1),slat(1),g2(ng),
Cgetarg      .  yl(ng,1),cs(ng),sn(ng),qlat(3,3)
       real(8):: qmom(1) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex cv(ng)
C ... Local parameters
      integer ig,ib,ilm,is,iv0,l,lmxl,m,nbas,nlm,nlmx,nglob,n1,n2,n3,
     .  ngabc(3),lfoc
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      parameter (nlmx=64)
      double precision tau(3),df(0:20),pi,y0,vol,rg,qcorg,qcorh,qsc,
     .  cofg,cofh,ceh,rfoc,z,q0(3),gam,gamf,cfoc,cvol,aa
      double complex cof(nlmx),cfac,phase
C      parameter (k0=3)
C      double precision gv(ng,3),v(3)
C      double complex gkl(0:k0,nlmx)
      data q0 /0d0,0d0,0d0/

      call tcn('pvdf4')
      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      nbas = nglob('nbas')
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat nabc vol',slat,ngabc,vol,0,0,0)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       vol=slat%vol


C --- FT of gaussian density, all sites, for list of G vectors ---
      iv0 = 0
      do  10  ib = 1, nbas
Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: tau(3) , real(8) :: ssite%pos 3
Cgetarg         call upack('site spec pos',ssite,ib,is,tau,0,0)
         
         is=ssite(ib)%spec
         i_copy_size=size(ssite(ib)%pos) 
         call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack('spec lmxl rg',sspec,is,lmxl,rg,0,0)
         
         lmxl=sspec(is)%lmxl
         rg=sspec(is)%rg

        if (lmxl .eq. -1) goto 10

        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
        nlm = (lmxl+1)**2
        if (nlm .gt. nlmx) call rxi('pvdf4: increase nlmx to',nlm)
        ilm = 0
        cfac = dcmplx(0d0,1d0)
        do  20  l = 0, lmxl
        cfac = cfac*dcmplx(0d0,-1d0)
        do  21  m = -l, l
          ilm = ilm+1
          cof(ilm) = cfac*qmom(ilm+iv0)*4*pi/df(2*l+1)
   21     continue
   20   continue
        cof(1) = cof(1) + 4*pi*y0*(qcorg-z)

        gam = 0.25d0*rg*rg
        gamf = 0.25d0*rfoc*rfoc
        cfoc = -4d0*pi*y0*cofh/vol
        cvol = 1d0/vol
        do  30  ig = 1, ng
          phase = dcmplx(cs(ig),sn(ig))
          aa = dexp(-gam*g2(ig))*cvol
          do  32  ilm = 1, nlm
   32     cv(ig) = cv(ig) + aa*yl(ig,ilm)*cof(ilm)*phase
C     ... Add foca hankel part
          aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
          cv(ig) = cv(ig) + cfoc*aa*phase
   30   continue

        iv0 = iv0+nlm
   10 continue

C --- Potential is 8pi/G**2 * density; overwrite cv with potential ---
      cv(1) = (0d0,0d0)
      do  40  ig = 2, ng
   40 cv(ig) = (8*pi)*cv(ig)/g2(ig)

      call tcx('pvdf4')
Cgetarg       end
       end subroutine pvdf4 


