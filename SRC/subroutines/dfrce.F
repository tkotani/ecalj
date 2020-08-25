      subroutine dfrce (job, ssite , sspec , slat , k1 , k2 , k3 !sctrl , 
     .  , nvl , sv_p_orhoat , sv_p_orhat1 , elind , qmom , smrho , smrout 
     .     , dfh )
      use m_work,only: rv_a_ogv,iv_a_okv
      use m_struc_def  !Cgetarg
      use m_globalvariables
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
Ci   smrout:smooth (output) density that the hamiltonian generated
Co Outputs
Co   dfh   :correction to the HF force
Cl Local variables
Cl    job  :describes which ansatz for charge shift is used for correction
Cl         :<=0  do not calculate correction to force
Cl         :  1  shift in free-atom density
Cl         :  2  shift in core+nuclear density
Cl         :+10  to screen the rigid shift by the Lindhard function
Cr Remarks
!! Density
!!  sv_p_orhoat: input atomic density that generaed Hamiltonian
!!  sv_p_orhat1: new atomic density that the Hamiltonian generated
!!  smrho:  input  density that generaed Hamiltonian
!!  smrout: output density that the Hamiltonian generated
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   18 Dec 03 adapted to modified smvxc
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   17 Sep 01 Adapted for local orbitals
Cu   21 Jun 00 spin polarized
Cu   18 Jun 98 adapted from nfp dfrce.f
Cu   16 Jun 98 MvS parallelized for SGI
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
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
      double precision starttime, endtime
      logical mlog,cmdopt
#endif
      integer:: nvl , k1 , k2 , k3
      type(s_rv1) :: sv_p_orhat1(3,*)
      type(s_rv1) :: sv_p_orhoat(3,*)
      real(8):: dfh(3,*) , qmom(*)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
c      type(s_ctrl)::sctrl

      double complex smrho(k1,k2,k3,*),smrout(k1,k2,k3,*)
C Local variables
      integer nbas,job,n1,n2,n3,ng,iprint,nsp,nglob,ib,is,lmxl,iv0,nlm,
     .ip,m,i,ngabc(3),ltop,nspec,nlmtop,igets,igetss,nn,lgunit,stdo
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
c        integer:: ocdv , ocdn , ocdn0 , ocs , osn
ckino Dec.9.2011:          integer,pointer :: iv_p_okv(:) =>NULL()

ckino Dec.9.2011:         real(8),pointer :: rv_p_ogv(:) =>NULL()

      complex(8) ,allocatable :: ceps_zv(:)
      complex(8) ,allocatable :: cnomi_zv(:)
      complex(8) ,allocatable :: smro_zv(:)
      complex(8) ,allocatable :: dvxc_zv(:)
      complex(8) ,allocatable :: vxcp_zv(:)
      complex(8) ,allocatable :: vxcm_zv(:)
      complex(8) ,allocatable :: cdvx_zv(:)
      real(8) ,allocatable :: qmout_rv(:)
      complex(8) ,allocatable :: cvin_zv(:)
      real(8) ,allocatable :: yl_rv(:)
      real(8) ,allocatable :: g2_rv(:)
      real(8) ,allocatable :: g_rv(:)
      integer ,allocatable :: iv_iv(:)
      complex(8) ,allocatable :: wk1_zv(:)
      complex(8) ,allocatable :: wk2_zv(:)
      complex(8) ,allocatable :: wk3_zv(:)

      double precision vol,plat(3,3),qlat(3,3),alat,vsum,pi,tpiba,elind,
     .fes1(3),fes2(3),fxc(3),c,avgdf(3)
      integer mp
c      integer npmx,mp
c      parameter (npmx=32)
c      integer oicdn(npmx),oicdn0(npmx),oicdv(npmx),oics(npmx),oisn(npmx)
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
      real(8),allocatable:: cs_(:),sn_(:)

c      job = int(sctrl%lfrce)
c      if (job .le. 0) return
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
     .    MPI_COMM_WORLD,ierr)
          call MPI_BCAST(namelen(pid),1,MPI_INTEGER,pid,
     .    MPI_COMM_WORLD,ierr)
        enddo
      endif
#endif

C --- Setup ---
      ngabc=slat%nabc
      ng=slat%ng
      vol=slat%vol
      alat=slat%alat
      plat=slat%plat
      qlat=slat%qlat
      c = 1000
      stdo = lgunit(1)
      nsp  = globalvariables%nsp
      nbas = globalvariables%nbas
      nspec= globalvariables%nspec
      nn   = k1*k2*k3

C ... Arrays needed for pvdf1
      allocate(ceps_zv(ng))
      allocate(cnomi_zv(ng))
      allocate(cvin_zv(ng))
      allocate(cdvx_zv(ng*nsp))

C ... Set up for vectorized Y_lm and gaussians
      ltop = 0
      do   is = 1, nspec
        lmxl = int(sspec(is)%lmxl)
        ltop = max0(ltop,lmxl)
      enddo
      nlmtop = (ltop+1)**2
      allocate(yl_rv(ng*nlmtop))
      allocate(g2_rv(ng))
      allocate(g_rv(ng*3))
      call suylg ( ltop , alat , ng , rv_a_ogv , g_rv , g2_rv , yl_rv)
      if (allocated(g_rv)) deallocate(g_rv)
      allocate(iv_iv(ng*3))
      call suphs0 ( plat , ng , rv_a_ogv , iv_iv )

C --- Make ves(rhoin,q) ---
      allocate(smro_zv(nn))
      allocate(cs_(ng), sn_(ng))
c      call dpcopy ( smrho , smro_zv , 1 , 2 * nn , 1d0 )
c     if ( nsp.eq.2 ) call dpadd ( smro_zv , smrho ( 1 , 1 , 1 , 2 ) , 1 , 2 * nn , 1d0 )
      call dcopy(2*nn, smrho,1,smro_zv,1)
      if(nsp==2) call daxpy(2*nn, 1d0, smrho (1,1,1,2), 1, smro_zv, 1)
      call fftz3 ( smro_zv , n1 , n2 , n3 , k1 , k2 , k3 , 1 , 0 , - 1 )
      call gvgetf ( ng , 1 , iv_a_okv , k1 , k2 , k3 , smro_zv , cvin_zv )
      call pvdf4 ( ssite , sspec , slat , qmom , ng , g2_rv , yl_rv , cs_ , sn_ , iv_iv , qlat , cvin_zv )
      if (allocated(smro_zv)) deallocate(smro_zv)
      deallocate(cs_,sn_)

C --- Make dVxc(in)/dn ---
      allocate(dvxc_zv(nn*nsp))
      allocate(smro_zv(nn*nsp))
      allocate(vxcp_zv(nn*nsp))
      allocate(vxcm_zv(nn*nsp))
      allocate(wk1_zv(nn*nsp))
      allocate(wk2_zv(nn*nsp))
      allocate(wk3_zv(nn*nsp))
      call dpcopy ( smrho , smro_zv , 1 , 2 * nn * nsp , 1d0 )
      call pvdf2 ( nbas , nsp , ssite , sspec , slat , n1 , n2 , n3 
     ., k1 , k2 , k3 , smro_zv , vxcp_zv , vxcm_zv , wk1_zv 
     ., wk2_zv , wk3_zv , dvxc_zv )
      deallocate(wk3_zv)
      deallocate(wk2_zv)
      deallocate(wk1_zv)
      deallocate(vxcm_zv)
      deallocate(vxcp_zv)

C --- cdvx = FFT ((n0_out-n0_in) dVxc/dn) ---
C     Use total n0_out-n0_in but keep vxc spin polarized
      call dpzero ( smro_zv , 2 * nn )
      do  i = 1, nsp
c        call dpadd ( smro_zv , smrout ( 1 , 1 , 1 , i ) , 1 , 2*nn , 1d0 )
c        call dpadd ( smro_zv , smrho ( 1 , 1 , 1 , i ) , 1 , 2*nn , - 1d0 )
        call daxpy (2*nn,  1d0, smrout( 1 , 1 , 1 , i ), 1, smro_zv,1)
        call daxpy (2*nn, -1d0, smrho ( 1 , 1 , 1 , i ), 1, smro_zv,1)
      enddo
C     call zprm3('drho',0,w(osmro),k1,k2,k3)
      call pvdf3 ( n1 , n2 , n3 , k1 , k2 , k3 , nsp , smro_zv , dvxc_zv )
      call fftz3 ( dvxc_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp , 0 , - 1 )
      call gvgetf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , dvxc_zv , cdvx_zv )

C --- Cnomi = (n0_out(q) - n0_in(q)) ---
      call fftz3 ( smro_zv , n1 , n2 , n3 , k1 , k2 , k3 , 1 , 0 , - 1 )
c     call zprm3('FFT smrout-smrin',w(osmro),k1,k2,k3*nsp)
      call gvgetf ( ng , 1 , iv_a_okv , k1 , k2 , k3 , smro_zv , cnomi_zv  )
      if (allocated(smro_zv)) deallocate(smro_zv)
      if (allocated(dvxc_zv)) deallocate(dvxc_zv)
C ... Debugging slot smrho(out) for out-in
C      print *, '*** debugging ... subs smrout for out-in'
C      call dpcopy(smrout,w(osmro),1,2*nn,1d0)
C      call fftz3(w(osmro),n1,n2,n3,k1,k2,k3,1,0,-1)
C      call gvgetf(ng,1,w(okv),k1,k2,k3,w(osmro),w(ocnomi))
C      call zprm3('rho-out(q)',w(osmro),k1,k2,k3)

C --- Multipole moments of the output density ---
      allocate(qmout_rv(nvl))
      call pshpr(0)
      call rhomom ( nbas , ssite , sspec , sv_p_orhat1 , qmout_rv , vsum )
      call poppr
c      call dpadd ( qmout_rv , qmom , 1 , nvl , - 1d0 )
      qmout_rv = qmout_rv- qmom(1:nvl)
C --- Lindhard dielectric function ---
      if (job .gt. 10) then
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
        call lindsc ( 3 , ng , rv_a_ogv , tpiba , elind , ceps_zv )
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
     .'  ib',9x,'delta-n dVes',13x,'delta-n dVxc',15x,'total')
C  201   format(/' Harris correction to forces:'/
C     .    '  ib',11x,'dn0 dVes',15x,'dnloc dVes',15x,
C     .    'dn0 dVxc',15x,'total')

C ... Setup array iiv0 = (vector of iv0 for parallel); allocate work arrays
#if MPI | MPIK
      iv0 = 0
      allocate(iiv0(1:nbas), stat=ierr)
      do  12  ib = 1, nbas
        is = int(ssite(ib)%spec)
        lmxl = int(sspec(is)%lmxl)
        nlm = (lmxl+1)**2
        iiv0(ib) = iv0
        iv0 = iv0+nlm
   12 continue
#endif
      mp = 1

C ... Estimate shift in density for each site, and corresponding force
      ip = 1
      iv0 = 0
C$DOACROSS LOCAL (ib,is,lmxl,nlm,i,m,iv0,ip,ocdn,ocdn0,ocdv,
C$&               ocs,osn,fes1,fes2,fxc)
C$&        SHARED (job,ng,iiv0,ocnomi,ocvin,ocdvx)
C$&        MP_SCHEDTYPE=RUNTIME
#if MPI | MPIK
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
     .    //shortname(procid)(1:namelen(procid))//
     .    ' starting atoms %i to %i',' ',256,lgunit(3),
     .    procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        iv0 = iiv0(ib)
#endif
        is = int(ssite(ib)%spec)
        lmxl = int(sspec(is)%lmxl)
        if (lmxl .eq. -1) goto 20
        nlm = (lmxl+1)**2
        call pvdf1 ( job , ssite , sspec , slat , nsp , ib , iv0 , qmom
     .    , qmout_rv , ng , rv_a_ogv , g2_rv , yl_rv , iv_iv , qlat , 0 
     .    , cnomi_zv , ceps_zv , cdvx_zv , cvin_zv , sv_p_orhoat ( 1 , 
     .    ib ) , fes1 , fes2 , fxc )
        do  i = 1, 3
          dfh(i,ib) = -(fes1(i) + fes2(i) + fxc(i))
        enddo
#if ! (MPI | MPIK)
        if (iprint() .ge. 30)
     .  write(stdo,200) ib,(c*(fes1(m)+fes2(m)),m=1,3),
     .  (c*fxc(m),m=1,3),(c*dfh(m,ib),m=1,3)
  200   format(i4,3f8.2,1x,3f8.2,1x,3f8.2:1x,3f8.2)
#endif
        iv0 = iv0+nlm
   20   continue
      enddo
#if MPI | MPIK
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do  pid = 0, numprocs-1
        ib = bproc(pid)
        jb = bproc(pid+1) - ib
        call MPI_BCAST(dfh(1,ib),3*jb,MPI_DOUBLE_PRECISION,pid,
     .  MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit6(' dfrce '//datim//' Process %i of %i on '
     .    //shortname(procid)(1:namelen(procid))//
     .    ' bcast dfh(1-3,%i-%i) %i d.p. numbers'//
     .    ' from process %i on '
     .    //shortname(pid)(1:namelen(pid)),' ',
     .    256,lgunit(3),procid,numprocs,bproc(pid),bproc(pid+1)-1,
     .    3*jb,pid)
        endif
      enddo
      deallocate (bproc, stat=ierr)
      deallocate (iiv0, stat=ierr)
#endif
      call dpzero(avgdf,3)
      do  23  ib = 1, nbas
      do  23  i = 1, 3
        avgdf(i) = avgdf(i) + dfh(i,ib)/nbas
   23 continue

C ... MPI printout
#if MPI | MPIK
      if (iprint() .ge. 30) then
        strn = 'shift in free-atom density'
        if (job .eq. 11) strn = 'screened shift in free-atom density'
        if (job .eq. 12) strn = 'screened shift in core+nuclear density'
        write(stdo,201) strn
        do  26  ib = 1, nbas
          write(stdo,202) ib, (c*dfh(i,ib),i=1,3)
   26   continue
  202   format(i4,50x,3f8.2)
      endif
#endif

C ... Shift all forces to make avg correction zero
      do  24  ib = 1, nbas
      do  24  i = 1, 3
        dfh(i,ib) = dfh(i,ib) - avgdf(i)
   24 continue
      if (iprint() .ge. 30) write(stdo,331) (c*avgdf(m),m=1,3)
  331 format(' shift forces to make zero average correction:',8x,3f8.2)
      if (allocated(qmout_rv)) deallocate(qmout_rv)
      if (allocated(iv_iv)) deallocate(iv_iv)
      if (allocated(g2_rv)) deallocate(g2_rv)
      if (allocated(yl_rv)) deallocate(yl_rv)
      if (allocated(cdvx_zv)) deallocate(cdvx_zv)
      if (allocated(cvin_zv)) deallocate(cvin_zv)
      if (allocated(cnomi_zv)) deallocate(cnomi_zv)
      if (allocated(ceps_zv)) deallocate(ceps_zv)
      call tcx('dfrce')
      end subroutine dfrce




c      subroutine pvdf1(job,ssite,sspec,slat,nsp,ib,iv0,qmom,qmout,ng,gv,
c     .g2,yl,cs,sn,iv,qlat,kmax,cnomin,ceps,cdn0,cdn,cdv,cdvxc,cvin,
c     .orhoat,fes1,fes2,fxc)
      subroutine pvdf1 ( job , ssite , sspec , slat , nsp , ib , iv0 
     . , qmom , qmout , ng , gv , g2 , yl , iv , qlat , kmax , cnomin 
     . , ceps , cdvxc , cvin , sv_p_orhoat , fes1 , fes2 , fxc )
      use m_struc_def  !Cgetarg
      use m_globalvariables
c need to modify texts.
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
      implicit none
C ... Passed parameters
      integer:: ng , nsp , iv0 , kmax , ib , job , iv(ng,3),i_copy_size
      type(s_rv1) :: sv_p_orhoat(3)

      real(8):: qmom(*) , qmout(*) , gv(ng,3) , tau(3) , fes1(3) , 
     .fes2(3) , fxc(3) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex cdn0(ng,nsp),cdn(ng),cdv(ng),ceps(ng),
     .cnomin(ng),cdvxc(ng,nsp),cvin(ng)
C ... Local parameters
      integer ig,ilm,l,lmxl,m,nlm,nlmx,k,is,jv0,jb,nbas,nglob,js,ll,n0,
     .nrmx
      parameter (nlmx=64, nrmx=1501, n0=10)
      integer lmxa,nr,nxi,ie,ixi,job0,kcor,lcor,lfoc,i,
     .ngabc(3),n1,n2,n3,nlml
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision pi,alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg,
     .vol,y0,z,v(3),df(0:20),feso(3),qcor(2),gpot0(nlmx,3),fesdn(3),
     .fesgg(3),pnu(n0,2),pnz(n0,2),a,rmt,qloc,rsmfa,exi(n0),hfc(n0,2),
     .qfat,gam,qall,qc,qval,qg,e,aa,q0(3),sum
      double precision rwgt(nrmx),cc,gamf,cfoc,cvol
C     parameter (k0=3)
C     double complex gkl(0:k0,nlmx)
      double complex tpia,cxx,phase,gc0,xc0,cof(nlmx)
C ... Heap
      data q0 /0d0,0d0,0d0/

      call tcn('pvdf1')

c      i_copy_size=size(slat%nabc)
c      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ngabc=slat%nabc
      alat=slat%alat
      vol=slat%vol

      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
Changenglob      nbas = nglob('nbas')
      nbas = globalvariables%nbas
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

      is=ssite(ib)%spec
c      i_copy_size=size(ssite(ib)%pos)
c      call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)
      tau=ssite(ib)%pos(1:3)
      call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)

C --- Unscreened rigid charge density shift, job 1, in cdn0 ---
      if (job0 .eq. 1) then
C       call upack('site spec pos',ssite,ib,is,tau,0,0)
c        call upack('spec z p pz lmxa lmxl',sspec,is,z,pnu,pnz,lmxa)

        z=sspec(is)%z
        i_copy_size=size(sspec(is)%p)
        call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
        i_copy_size=size(sspec(is)%pz)
        call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)
        lmxa=sspec(is)%lmxa


        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        lmxl=sspec(is)%lmxl


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
        call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 1 )%v 
     .   , qloc )

        call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 2 )%v 
     .   , sum )

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
            cdn0(ig,i) = cdn0(ig,i) + aa*dcmplx(cs(ig),sn(ig))
   15     continue
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
          cdn0(ig,i)=cdn0(ig,i)+cc*dcmplx(cs(ig),sn(ig))*dexp(-gam*g2(ig))
   16   continue
      endif

C --- Coefficients defining local valence + core density ---

      is=ssite(ib)%spec
      i_copy_size=size(ssite(ib)%pos)
      call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)


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
          cof(ilm) = cxx*qmom(ilm+iv0)*4*pi/df(2*l+1)
   22   continue
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
        fes1(k) = fes1(k)*vol
   37 continue

C --- Integral of grad g (output-input local charge) ves~ ---
      fesgg=0d0 !call dpzero(fesgg,3)
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

        js=ssite(jb)%spec
        i_copy_size=size(ssite(jb)%pos)
        call dcopy(i_copy_size,ssite(jb)%pos,1,tau,1)


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
     .    dconjg(tpia*cdv(ig))*dcmplx(cs(ig),sn(ig))
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
   64       continue
   62     continue
   60   continue
        jv0 = jv0+nlm
   40 continue
C      print 339, 'qmom(in) dv         ',feso
C      print 339, 'qmom(out-in) dv     ',fes2
c      call dpadd(fes2,fesgg,1,3,-1d0)
      fes2=fes2-fesgg
      call tcx('pvdf1')
      end subroutine pvdf1


      subroutine pvdf2(nbas,nsp,ssite,sspec,slat,n1,n2,n3,k1,k2,k3,
     .smrho,vxcp,vxcm,wk1,wk2,wk3,dvxc)
      use m_struc_def
C- Makes derivative of smoothed xc potential wrt density.
      implicit none
C ... Passed parameters
      integer nbas,nsp,n1,n2,n3,k1,k2,k3
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      complex(8):: vxcp(k1,k2,k3,nsp),vxcm(k1,k2,k3,nsp),
     .dvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp),
     .wk1(k1,k2,k3,nsp),wk2(k1,k2,k3,nsp),
     .wk3(k1,k2,k3,nsp)
C ... Local parameters
      integer i1,i2,i3,i,nn
      double precision fac,dmach,f1,f2,f,alfa,dfdr,rrho,dvdr,
     .rmusm(2),rvmusm(2),rvepsm(2),repsm(2),repsmx(2),repsmc(2),
     .fcexc0(2),fcex0(2),fcec0(2),fcvxc0(2)

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
     .vxcp,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

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
     .vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

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
      do i = 1, nsp
        do i3=1,n3
          do i2=1,n2
            do i1=1,n1
              rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
              if (rrho .gt. 0) then
                f1 = vxcm(i1,i2,i3,i)*(rrho*(1-fac))**alfa
                f2 = vxcp(i1,i2,i3,i)*(rrho*(1+fac))**alfa
                dfdr = (f2-f1)/(2d0*fac*rrho)
                vxcp(i1,i2,i3,i) = dfdr
              else
                vxcp(i1,i2,i3,i) = 0d0
              endif
            enddo
          enddo
        enddo
      enddo
C$$$      do  i = 1, nsp
C$$$        do i1 = 1,nn
C$$$          rrho = (smrho(i1,1,1,1)+smrho(i1,1,1,nsp))/(3-nsp)
C$$$          if (rrho .gt. 0) then
C$$$            f1 = vxcm(i1,1,1,i)*(rrho*(1-fac))**alfa
C$$$            f2 = vxcp(i1,1,1,i)*(rrho*(1+fac))**alfa
C$$$            dfdr = (f2-f1)/(2d0*fac*rrho)
C$$$            print *,' f f f=',f1,f2,dfdr,vxcp(1,1,1,1)
C$$$            vxcp(i1,1,1,i) = dfdr
C$$$          else
C$$$            vxcp(i1,1,1,i) = 0d0
C$$$          endif
C$$$        enddo
C$$$      enddo
C$$$      endif

C ... vxcm = vxc (smrho)
      call dpzero(vxcm, nn*2*nsp)
      call smvxcm(ssite,sspec,slat,nbas,0,k1,k2,k3,smrho,
     .vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm,
     .rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,f)

C ... dvxc/drho into dvxc
      do  i = 1, nsp
        do i3=1,n3
          do i2=1,n2
            do i1=1,n1
              rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
              if (rrho .gt. 0) then
                f = vxcm(i1,i2,i3,i) * rrho**alfa
                dvdr = (vxcp(i1,i2,i3,i) - alfa*f/rrho) / rrho**alfa
                dvxc(i1,i2,i3,i) = dvdr
              else
                dvxc(i1,i2,i3,i) = 0d0
              endif
            enddo
          enddo
        enddo
      enddo
C     call zprm3('d vxc / dn',0,dvxc,k1,k2,k3*nsp)
      call poppr
      end subroutine pvdf2



      subroutine pvdf3(n1,n2,n3,k1,k2,k3,nsp,deln0,dvxc)
C- Overwrites dvxc with (nout-nin)*dvxc
      implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,nsp
      double complex deln0(k1,k2,k3),dvxc(k1,k2,k3,nsp)
C ... Local parameters
      integer i1,i2,i3,i

      do  10  i  = 1, nsp
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1
        dvxc(i1,i2,i3,i) = dvxc(i1,i2,i3,i)*deln0(i1,i2,i3)
   10 continue

C     call zprm3('dvxc/dn * (nout-nin)',0,dvxc,k1,k2,k3*nsp)

      end subroutine pvdf3


      subroutine pvdf4(ssite,sspec,slat,qmom,ng,g2,yl,cs,sn,iv,qlat,cv)
      use m_struc_def  !Cgetarg
      use m_globalvariables
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
      implicit none
C ... Passed parameters
      integer ng,iv(ng,3),i_copy_size
      real(8):: qmom(1) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex cv(ng)
C ... Local parameters
      integer ig,ib,ilm,is,iv0,l,lmxl,m,nbas,nlm,nlmx,nglob,n1,n2,n3,
     .ngabc(3),lfoc
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      parameter (nlmx=64)
      double precision tau(3),df(0:20),pi,y0,vol,rg,qcorg,qcorh,qsc,
     .cofg,cofh,ceh,rfoc,z,q0(3),gam,gamf,cfoc,cvol,aa
      double complex cof(nlmx),cfac,phase
C      parameter (k0=3)
C      double precision gv(ng,3),v(3)
C      double complex gkl(0:k0,nlmx)
      data q0 /0d0,0d0,0d0/

      call tcn('pvdf4')
      call stdfac(20,df)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
Changenglob      nbas = nglob('nbas')
      nbas = globalvariables%nbas

c      i_copy_size=size(slat%nabc)
c      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      ngabc=slat%nabc
      vol=slat%vol


C --- FT of gaussian density, all sites, for list of G vectors ---
      iv0 = 0
      do  10  ib = 1, nbas

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pos)
        call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)


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
            cv(ig) = cv(ig) + aa*yl(ig,ilm)*cof(ilm)*phase
   32     continue
C     ... Add foca hankel part
          aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
          cv(ig) = cv(ig) + cfoc*aa*phase
   30   continue

        iv0 = iv0+nlm
   10 continue

C --- Potential is 8pi/G**2 * density; overwrite cv with potential ---
      cv(1) = (0d0,0d0)
      do  40  ig = 2, ng
        cv(ig) = (8*pi)*cv(ig)/g2(ig)
   40 continue

      call tcx('pvdf4')
      end subroutine pvdf4


