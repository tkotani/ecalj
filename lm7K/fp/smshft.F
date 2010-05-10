Cdelw1 var ['ocgr', '0', 'ocgr', 'ng*nsp', 'complex(8)', 'zv_w_', 'smshft']
Cdelw1 var ['okv', '-1', 'undef', '0', 'undef', 'x', 'smshft']
Cdelw1 var ['ocgs', '0', 'ocgs', 'ng*nsp', 'complex(8)', 'zv_w_', 'smshft']
Cdelw1 var ['ocwk', '0', 'ocwk', 'ng', 'complex(8)', 'zv_w_', 'smshft']
Cdelw1 var ['ogv', '-1', 'undef', '0', 'undef', 'x', 'smshft']
Cdelw1 undel [['undef', 'wref'], 'okv', 'smshft']
Cdelw1 undel [['undef', 'wref'], 'ogv', 'smshft']
Cdelw1 undel [['call'], 'orhoat', 'pvsms1']
Cgetarg...info...           structure ['smshft', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['smshft', 'lat', 'slat']
Cgetarg...info...           structure ['smshft', 'site', 'ssite']
Cgetarg...info...           structure ['smshft', 'spec', 'sspec']
Cgetarg...info...           structure ['smshft', 'ham', 'sham']
Cgetarg...info...           structure ['pvsms1', 'lat', 'slat']
Cgetarg...info...           structure ['pvsms1', 'site', 'ssite']
Cgetarg...info...           structure ['pvsms1', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['smshft', 'm_struc_def']
Cgetarg...info...           use_to_add ['pvsms1', 'm_struc_def']
      subroutine smshft(ssite,sspec,slat,sctrl,sham,orhoat,smrho)
       
       use m_struc_def  !Cgetarg

C- Estimate the smooth density for a shift in atomic positions.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos pos0
Ci     Duplicate: spec spec
Ci     Stored:    pos
Ci     Passed to: pvsms1 rhgcmp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z p pz lmxa lmxl a nr rmt nxi exi chfa rsmfa
Ci     Duplicate: z p pz lmxa
Ci     Stored:    *
Ci     Passed to: pvsms1 gtpcor rhgcmp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat qlat nabc ng ogv okv vol
Ci     Duplicate: alat vol
Ci     Stored:    *
Ci     Passed to: pvsms1 symsmr rhgcmp
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: lfrce
Ci     Stored:    *
Ci     Passed to: *
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: elind
Ci     Stored:    *
Ci     Passed to: *
Ci   orhoat:vector of offsets containing site density
Co Outputs
Co   smrho :a perturbation is added to smrho, depending on job
Cr Remarks
Cr   job describes which ansatz for charge shift is used for correction
Cr     <=0  do not calculate correction to force
Cr       1  shift in free-atom density
Cr       2  shift in core+nuclear density
Cr     +10  to screen the rigid shift by the Lindhard function
Cr   (job taken from ctrl->lfrce)
Cu Updates
Cu   17 Sep 01 Adapted for local orbitals.  Altered argument list
Cu    3 Jul 00 Adapted from nfp smshft.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer orhoat(3,1)
Cgetarg       double precision ssite(1),sspec(1),slat(1),sctrl(1),sham(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_ctrl)::sctrl
       type(s_ham)::sham

      double complex smrho(*)
C ... Local parameters
      integer i,ib,igets,igetss,iprint,is,k1,k2,k3,kcor,lcor,lmxa,n0,n1,
     .  n2,n3,nbas,ng,ngabc(3),nglob,nsp
Cdelw1 do not change  [['undef', 'wref'], 'ogv', 'smshft']
Cdelw1 do not change  [['undef', 'wref'], 'okv', 'smshft']
Cdelw1       integer ocgr,ocgs,ocwk,ogv,okv
       integer:: ogv , okv 
       complex(8) ,allocatable :: zv_w_ocgr(:)
       complex(8) ,allocatable :: zv_w_ocgs(:)
       complex(8) ,allocatable :: zv_w_ocwk(:)

      integer kmax,job
      parameter (n0=10)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,dgets,elind,pi,qc,qsc,qcor(2),plat(3,3),
     .  qlat(3,3),qv,qval,tpiba,vol,z,pnu(n0,2),pnz(n0,2)
C ... Heap
      integer w(1)
      common /w/ w

C     stdo = lgunit(1)
Cgetarg       job = igets('ctrl lfrce',sctrl)
       job = int(sctrl%lfrce) 

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [okv ogv]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (job .le. 0) return

      call tcn('smshft')

C --- Setup and printout ---
      nsp  = nglob('nsp')
      nbas = nglob('nbas')
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: qlat(3,3) , real(8) :: slat%qlat 9
Cgetarg       call upack('lat alat plat qlat',slat,alat,plat,qlat,0,0)
       
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       i_copy_size=size(slat%qlat) 
       call dcopy(i_copy_size,slat%qlat,1,qlat,1) 

Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat nabc ng ogv okv vol',slat,ngabc,ng,ogv,okv,vol)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       ng=slat%ng
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 111
Cdelw1 do not change ogv because of ['undef', 'wref']
       ogv=slat%ogv
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 112
Cdelw1 do not change okv because of ['undef', 'wref']
       okv=slat%okv
       vol=slat%vol

      call fftz30(n1,n2,n3,k1,k2,k3)
C     call zprm3('input smrho',0,smrho,k1,k2,k3)

C ... Hold on to original smrho (cgr)
Cdelw1       call defcc (ocgr,  ng*nsp)
       allocate(zv_w_ocgr(ng*nsp))
       if (ng*nsp<0) zv_w_ocgr(:)=0.0d0

      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 121
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 121
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1       call gvgetf(ng,nsp,w(okv),k1,k2,k3,smrho,w(ocgr))
       call gvgetf ( ng , nsp , w ( okv ) , k1 , k2 , k3 , smrho , zv_w_ocgr 
     . ) 


C --- Shift in unscreened density at the two positions ---
Cdelw1       call defcc (ocgs, ng*nsp)
       allocate(zv_w_ocgs(ng*nsp))
       if (ng*nsp<0) zv_w_ocgs(:)=0.0d0

Cdelw1       call defcc (ocwk, ng)
       allocate(zv_w_ocwk(ng))
       if (ng<0) zv_w_ocwk(:)=0.0d0

      kmax = 0
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 127
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 127
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call pvsms1(ssite,sspec,slat,nbas,nsp,kmax,
Cdelw1      .  ng,w(ogv),orhoat,w(ocwk),w(ocgs),job)
       call pvsms1 ( ssite , sspec , slat , nbas , nsp , kmax , ng , 
     . w ( ogv ) , orhoat , zv_w_ocwk , zv_w_ocgs , job ) 

Cdelw1 rlse name= ocwk old_list= ocgr ocgs ocwk 
Cdelw1 rlse name= ocwk new_list= ocgr ocgs 
Cdelw1       call rlse(ocwk)
       if (allocated(zv_w_ocwk)) deallocate(zv_w_ocwk)


C ... Debugging: print unscreened shift in pseudo core density
C      call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocgs),smrho)
C      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C      call zprm3('unscreened local density',0,smrho,k1,k2,k3)

C --- Screened shift ---
      if (job .gt. 10) then
C       Compute elind if not given
        qval = 0d0
        do  12  ib = 1, nbas
Cgetarg           is = igetss('site spec',ib,ssite)
           is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: pnu(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: pnz(n0,2) , real(8) :: sspec%pz 20
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg           call upack('spec z p pz lmxa',sspec,is,z,pnu,pnz,lmxa)
           
           z=sspec(is)%z
           i_copy_size=size(sspec(is)%p) 
           call dcopy(i_copy_size,sspec(is)%p,1,pnu,1) 
           i_copy_size=size(sspec(is)%pz) 
           call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1) 
           lmxa=sspec(is)%lmxa

          if (lmxa .eq. -1) goto 12
          call gtpcor(sspec,is,kcor,lcor,qcor)
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
          qval = qval+qv
   12   continue
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
Cgetarg         elind = dgets('ham elind',sham)
         elind = sham%elind 

        if (elind .lt. 0d0) elind = -(3*pi**2*qval/vol)**.66666d0*elind
Cdelw1         if (nsp.eq.2) call dsumdf(ng*2,1d0,w(ocgs),0,1,w(ocgs),ng*2,1)
         if ( nsp.eq.2 ) call dsumdf ( ng * 2 , 1d0 , zv_w_ocgs , 0 , 
     .   1 , zv_w_ocgs , ng * 2 , 1 ) 

C        call lindxx(122,n1,n2,n3,k1,k2,k3,ng,w(okv),w(ocgs),w(ogv),
C     .    tpiba,elind,w,w,w,w)
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 171
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 171
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1         call lindsc(2,ng,w(ogv),tpiba,elind,w(ocgs))
         call lindsc ( 2 , ng , w ( ogv ) , tpiba , elind , zv_w_ocgs 
     .   ) 

C        call dscal(2*ng,.001d0,w(ocgs),1)

Cdelw1         if (nsp.eq.2) call dsumdf(ng*2,.5d0,w(ocgs),0,1,w(ocgs),ng*2,1)
         if ( nsp.eq.2 ) call dsumdf ( ng * 2 , .5d0 , zv_w_ocgs , 0 , 
     .   1 , zv_w_ocgs , ng * 2 , 1 ) 

      endif

C ... Debugging: show delta smrho
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 178
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 178
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1       call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocgs),smrho)
       call gvputf ( ng , nsp , w ( okv ) , k1 , k2 , k3 , zv_w_ocgs 
     . , smrho ) 

      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C      call zprm3('screened delta smrho',0,smrho,k1,k2,k3)
C      print *, 'returning with delta smrho'
C      return

C --- Add shift to smrho, ensuring no shift in <rho> ---
      do  i = 1, nsp
Cdelw1         call dvset(w(ocgs),1+ng*(i-1),1+ng*(i-1),0d0)
         call dvset ( zv_w_ocgs , 1 + ng * ( i - 1 ) , 1 + ng * ( i - 
     .   1 ) , 0d0 ) 

      enddo
Cdelw1       call dpadd(w(ocgr),w(ocgs),1,ng*2*nsp,1d0)
       call dpadd ( zv_w_ocgr , zv_w_ocgs , 1 , ng * 2 * nsp , 1d0 ) 

Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 189
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 189
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1       call gvputf(ng,nsp,w(okv),k1,k2,k3,w(ocgr),smrho)
       call gvputf ( ng , nsp , w ( okv ) , k1 , k2 , k3 , zv_w_ocgr 
     . , smrho ) 

      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)

C --- Symmetrize the shifted density ---
      call symsmr(slat,nsp,k1,k2,k3,smrho)

      if (iprint() .gt. 100)
     .  call zprm3('shifted smrho',0,smrho,k1,k2,k3*nsp)

Cdelw1 rlse name= ocgr old_list= ocgr ocgs 
Cdelw1 rlse name= ocgr new_list= (None)
Cdelw1       call rlse(ocgr)
       if (allocated(zv_w_ocgs)) deallocate(zv_w_ocgs)
       if (allocated(zv_w_ocgr)) deallocate(zv_w_ocgr)

      call tcx('smshft')

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [okv ogv]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine smshft 


      subroutine pvsms1(ssite,sspec,slat,nbas,nsp,kmax,ng,
     .  gv,orhoat,cwk,cg,job)
       
       use m_struc_def  !Cgetarg

C- Shift in smoothed density according to job.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos pos0
Ci     Stored:    pos
Ci     Passed to: rhgcmp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxl z p pz lmxa a nr rmt nxi exi chfa rsmfa
Ci     Stored:    name
Ci     Passed to: spacks gtpcor rhgcmp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat vol
Ci     Stored:    *
Ci     Passed to: rhgcmp
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   cy    :Normalization constants for spherical harmonics
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   orhoat:vector of offsets containing site density
Ci   job :describes which ansatz for charge shift is used for correction
Ci         <=0  do not calculate correction to force
Ci           1  shift in free-atom density
Ci           2  shift in core+nuclear density
Ci         +10  to screen the rigid shift by the Lindhard function
Co Outputs
Co  cg   coefficients to FT of shifted density
Cl Local variables
Cl  qloc   :difference in true and smoothed local charge.
Cl         :If the charge is assembled from overlapping atom-centered
Cl         :densities, qloc is the difference between the smoothed
Cl         :and head densities.
Cr Remarks
Cr   Shift of the "free atom densities."
Cr   The table below shows the densities and corresponding charges, and
Cr   parameters that hold their representations (true and smoothed
Cr   approximate forms):
Cr      density     charge    reps'n     smooth -reps'n
Cr      rho(smH)     qfat    cofh,ceh    already smooth
Cr      rho1-rho2    qloc     rhoat      qkl
Cr      rhoc         qc
Cr   This routine constructs the following difference
Cr   rhat(final) - rhat(initial)  positions, in the smooth reps'n, where
Cr      rhat = rho(smH) + qg * g(r)
Cr   where
Cr       qg = qval+qsc-qfat-qloc
Cr
Cr   In the special case rho is assembled from a superposition of
Cr   free-atom densities, and rho1-rho2 = rhoval(free-atm)-rho(smH)
Cr   (see ovlcor.f).  Thus in this case:
Cr      rho(free-atom) = rho(smH) + rho1-rho2 + rhoc
Cr   with the corresponding integrated charges
Cr       qval=z-qc     = qfat     + qloc      - qc
Cr   Thus in this special case qg=0: the only shift comes from rho(smH).
Cr   Because the local density (which contains the remaining part of
Cr   the free-atom density) will automatically be shifted, it follows
Cr   that the shifted smooth density will correspond to the
Cr   smooth sum-of-FA densities constructed at the shifted positions.
Cr
Cr   In the general case, qg is not zero.  By shifting the a gaussian
Cr   along with the sm-Hankels, the integrated total density of charge
Cr   shifted (local density + mesh density) is neutral.
Cr
Cr   Improvements: if the tail density were also shifted inside each
Cr   augmentation sphere, the total density would correspond exactly
Cr   to the sum-of-FA densities at the shifted positions, when the
Cr   starting density is also a sum-of-FA densities.
Cr
Cr   Shift of the "core + valence densities."
Cr
Cu Updates
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
Cdelw1 do not change  [['call'], 'orhoat', 'pvsms1']
      integer nbas,nsp,ng,job,orhoat(3,1)
Cgetarg       double precision ssite(1),sspec(1),slat(1),gv(ng,3)
       real(8):: gv(ng,3) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex cg(ng,nsp),cwk(ng)
C Local variables
      integer ib,i,is,iv0,kmax,lmxl,igetss,lmxa,nr,nrmx,n0,
     .  nlml,nxi,ie,ixi,ig,ipr,iprint,kcor,lcor,stdo,lgunit
      parameter (nrmx=1501, n0=10)
      double precision a,aa,alat,df(0:20),e,
     .  exi(n0),gam,hfc(n0,2),pi,pnew(3),pnu(n0,2),pnz(n0,2),pold(3),
     .  pp,qall,qc,qcor(2),qsc,qfat,qg,qloc,qval,rmt,
     .  rsmfa,rwgt(nrmx),scalp,sum,tpiba,v(3),v2,vol,volsp,y0,z
      character*35 strn,spid*8
      double complex phase
C ... Heap
      integer w(1)
      common /w/ w

C ... Setup
      call tcn('pvsms1')
      call stdfac(20,df)
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat alat vol',slat,alat,vol,0,0,0)
       
       alat=slat%alat
       vol=slat%vol

      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      tpiba = 2*pi/alat
      ipr = iprint()
      stdo = lgunit(1)
      volsp = vol*nsp

C --- For each site, accumulate shift in density ---
      call dpzero(cg,2*ng*nsp)
      iv0 = 0
      strn = 'free atom densities'
      if (job .eq. 11) strn = 'screened free atom densities'
      if (job .eq. 12) strn = 'screened core+multipole densities'
      if (ipr .ge. 30) write(stdo,339) strn
  339 format(/' smshft:  add shifted ',a/'   site',16x,'old pos',22x,
     .  'new pos',14x,'shift')

      do  10  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg         call spacks(0,'spec name',sspec,spid,is,is)
         do i_spacks=is,is 
         call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
         enddo

Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        nlml = (lmxl+1)**2
        if (lmxl .eq. -1) goto 10

Cgetarg...info...           integer :: is , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: pnew(3) , real(8) :: ssite%pos 3
Cgetarg...info...           real(8) :: pold(3) , real(8) :: ssite%pos0 3
Cgetarg         call upack('site spec pos pos0',ssite,ib,is,pnew,pold,0)
         
         is=ssite(ib)%spec
         i_copy_size=size(ssite(ib)%pos) 
         call dcopy(i_copy_size,ssite(ib)%pos,1,pnew,1) 
         i_copy_size=size(ssite(ib)%pos0) 
         call dcopy(i_copy_size,ssite(ib)%pos0,1,pold,1) 

        pp = alat*dsqrt((pnew(1)-pold(1))**2 + (pnew(2)-pold(2))**2
     .                + (pnew(3)-pold(3))**2)
        if (ipr .ge. 30) write(stdo,340) ib,spid,pold,pnew,pp/alat
  340   format(i4,':',a,f8.5,2f9.5,2x,3f9.5,2x,f9.6)
C       Skip this site if shift is negligible
        if (pp .le. 1d-6) goto 18

C   --- Shift in mesh density, job 1 ---
        if (mod(job,10) .eq. 1) then
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: pnu(n0,2) , real(8) :: sspec%p 20
Cgetarg...info...           real(8) :: pnz(n0,2) , real(8) :: sspec%pz 20
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg           call upack('spec z p pz lmxa',sspec,is,z,pnu,pnz,lmxa)
           
           z=sspec(is)%z
           i_copy_size=size(sspec(is)%p) 
           call dcopy(i_copy_size,sspec(is)%p,1,pnu,1) 
           i_copy_size=size(sspec(is)%pz) 
           call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1) 
           lmxa=sspec(is)%lmxa

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg           call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
           
           a=sspec(is)%a
           nr=sspec(is)%nr
           rmt=sspec(is)%rmt

Cgetarg...info...           integer :: nxi , integer(8) :: sspec%nxi 1
Cgetarg...info...           real(8) :: exi(n0) , real(8) :: sspec%exi 10
Cgetarg...info...           real(8) :: hfc(n0,2) , real(8) :: sspec%chfa 20
Cgetarg...info...           real(8) :: rsmfa , real(8) :: sspec%rsmfa 1
Cgetarg           call upack('spec nxi exi chfa rsmfa',sspec,is,nxi,exi,hfc,
Cgetarg      .      rsmfa)
           
           nxi=sspec(is)%nxi
           i_copy_size=size(sspec(is)%exi) 
           call dcopy(i_copy_size,sspec(is)%exi,1,exi,1) 
           i_copy_size=size(sspec(is)%chfa) 
           call dcopy(i_copy_size,sspec(is)%chfa,1,hfc,1) 
           rsmfa=sspec(is)%rsmfa

          call gtpcor(sspec,is,kcor,lcor,qcor)
          if (nr .gt. nrmx) call rx('dfrce: nr gt nrmx')
          call radwgt(rmt,a,nr,rwgt)
Cdelw1 1 is digit, dropped
          call radsum(nr,nr,nlml,nsp,rwgt,w(orhoat(1,ib)),qloc)
Cdelw1 2 is digit, dropped
          call radsum(nr,nr,nlml,nsp,rwgt,w(orhoat(2,ib)),sum)
          qloc = (qloc-sum)/y0
          qfat = 0d0
          do  12  i  = 1, nsp
          do  12  ie = 1, nxi
            gam  = 0.25d0*rsmfa**2
            qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
            qfat = qfat + hfc(ie,i)*qall
   12     continue
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
C         Excess sphere charge.  See Remarks above.
          qg = qval+qsc-qfat-qloc

C     ... Shift in smoothed free atom density
          do  14    i = 1, nsp
          do  14  ixi = 1, nxi
          e = exi(ixi)
          do  15 ig = 1, ng
            v(1) = gv(ig,1)*tpiba
            v(2) = gv(ig,2)*tpiba
            v(3) = gv(ig,3)*tpiba
            v2 = v(1)**2+v(2)**2+v(3)**2
            aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
            scalp = -alat*(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
            phase = dcmplx(dcos(scalp),dsin(scalp))
            scalp = -alat*(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
            phase = phase - dcmplx(dcos(scalp),dsin(scalp))
            cg(ig,i) = cg(ig,i) + hfc(ixi,i)*aa*phase*y0/vol
   15     continue
   14     continue

C     ... Add gaussian to conserve local charge; see Remarks
          do  16   i = 1, nsp
          do  16  ig = 1, ng
            v(1) = gv(ig,1)*tpiba
            v(2) = gv(ig,2)*tpiba
            v(3) = gv(ig,3)*tpiba
            v2 = v(1)**2+v(2)**2+v(3)**2
            scalp = -alat*(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
            phase = dcmplx(dcos(scalp),dsin(scalp))
            scalp = -alat*(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
            phase = phase - dcmplx(dcos(scalp),dsin(scalp))
            cg(ig,i) = cg(ig,i) + qg*phase*dexp(-gam*v2)/volsp
   16     continue

C   --- Shift in mesh density, job 12 ---
        elseif (job .eq. 12) then

C     ... Core + valence at old position
          call dpzero(cwk,ng*2)
Cgetarg...info...           real(8) :: pold(3) , real(8) :: ssite%pos 3
Cgetarg           call pack2('site pos',ssite,ib,pold)
           
           i_copy_size=size(ssite(ib)%pos) 
           call dcopy(i_copy_size,pold,1,ssite(ib)%pos,1) 

          call rhgcmp(131,ib,ib,ssite,sspec,slat,orhoat,kmax,ng,cwk)
          call dscal(ng*2,-1d0,cwk,1)
C     ... Core + valence at new position
Cgetarg...info...           real(8) :: pnew(3) , real(8) :: ssite%pos 3
Cgetarg           call pack2('site pos',ssite,ib,pnew)
           
           i_copy_size=size(ssite(ib)%pos) 
           call dcopy(i_copy_size,pnew,1,ssite(ib)%pos,1) 

          call rhgcmp(131,ib,ib,ssite,sspec,slat,orhoat,kmax,ng,cwk)

C     ... Add to cg
          do  i = 1, nsp
            call daxpy(ng*2,1d0/nsp,cwk,1,cg(1,i),1)
          enddo

        else
          call rxi('smshft: bad job:',job)
        endif

   18   iv0 = iv0+nlml
   10 continue

      call tcx('pvsms1')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhoat]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvsms1 


