subroutine smshft (job,ppnew,ppold)
  use m_density,only:  sv_p_orhoat=>orhoat,smrho=>osmrho !input and output
  use m_supot,only:  iv_a_okv,rv_a_ogv
  use m_lmfinit,only: ctrl_lfrce,ham_elind, &
       lat_alat,nsp,nbas, sspec=>v_sspec, ssite=>v_ssite,n0
  use m_lattic,only: lat_qlat, lat_vol,lat_plat
  use m_supot,only: lat_nabc, lat_ng
  use m_lgunit,only:stdo
  !- Estimate the smooth density for a shift in atomic positions.
  ! nput  orhoat:vector of offsets containing site density
  ! nputoutput
  !   smrho :a perturbation is added to smrho, depending on job
  !r Remarks
  !r   job describes which ansatz for charge shift is used for correction
  !r     <=0  do not calculate correction to force
  !r       1  shift in free-atom density
  !r       2  shift in core+nuclear density
  !r     +10  to screen the rigid shift by the Lindhard function
  !u Updates
  !u   17 Sep 01 Adapted for local orbitals.  Altered argument list
  !u    3 Jul 00 Adapted from nfp smshft.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: i,ib,igets,igetss,iprint,is,k1,k2,k3,kcor,lcor,lmxa,n1,n2,n3,ng,ngabc(3),i_copy_size
  complex(8) ,allocatable :: cgr_zv(:,:)
  complex(8) ,allocatable :: cgs_zv(:,:),cgs_zvv(:,:)
  complex(8) ,allocatable :: cwk_zv(:)
  integer :: kmax,job
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  double precision :: alat,dgets,elind,pi,qc,qsc,qcor(2),plat(3,3), &
       qlat(3,3),qv,qval,tpiba,vol,z,pnu(n0,2),pnz(n0,2)
  real(8)::ppnew(3,nbas),ppold(3,nbas)
  if (job <= 0) return
  call tcn('smshft')
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
  ngabc=lat_nabc
  ng=lat_ng
  vol=lat_vol
  call fftz30(n1,n2,n3,k1,k2,k3)
  ! ... Hold on to original smrho (cgr)
  allocate(cgr_zv(ng,nsp))
  call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
  call gvgetf(ng , nsp , iv_a_okv , k1 , k2 , k3 , smrho , cgr_zv  )
  ! --- Shift in unscreened density at the two positions ---
  allocate(cgs_zv(ng,nsp),cgs_zvv(ng,nsp),cwk_zv(ng))
  kmax = 0
  call pvsms1 ( ssite , sspec ,  nbas , nsp , kmax , ng ,&
       rv_a_ogv , sv_p_orhoat , cwk_zv , cgs_zv , job ,ppnew,ppold)
  if (allocated(cwk_zv)) deallocate(cwk_zv)
  ! --- Screened shift ---
  if (job > 10) then !Compute elind if not given
     qval = 0d0
     do  12  ib = 1, nbas
        is  =ssite(ib)%spec
        z   =sspec(is)%z
        pnu =ssite(ib)%pnu
        pnz =ssite(ib)%pz
        lmxa=sspec(is)%lmxa
        if (lmxa == -1) goto 12
        call gtpcor(sspec,is,kcor,lcor,qcor)
        call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
        qval = qval+qv
12   enddo
     pi = 4d0*datan(1d0)
     tpiba = 2*pi/alat
     elind = ham_elind
     if (elind < 0d0) elind = -(3*pi**2*qval/vol)**.66666d0*elind
     cgs_zvv(:,1)= cgs_zv(:,1) + cgs_zv(:,2)
     cgs_zvv(:,2)= cgs_zv(:,1) - cgs_zv(:,2)
     call lindsc( 2 , ng , rv_a_ogv , tpiba , elind , cgs_zvv )
  endif
  ! ... Debugging: show delta smrho
  call gvputf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , cgs_zv , smrho )
  call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
  !     --- Add shift to smrho, ensuring no shift in <rho> ---
  cgs_zv(1,1:nsp)=0d0 !real part zero bug? 2022-5-27
  !    do  i = 1, nsp
  !      call dvset(cgs_zv,1 + ng*( i - 1 ), 1 + ng*(i-1 ) , 0d0 )<=== real part zero bug? 2022-5-27
  !    enddo
  cgr_zv = cgr_zv + cgs_zv
  call gvputf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , cgr_zv , smrho )
  call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
  ! --- Symmetrize the shifted density ---
  call symsmrho(smrho)
  if (allocated(cgs_zv)) deallocate(cgs_zv,cgs_zvv)
  if (allocated(cgr_zv)) deallocate(cgr_zv)
  call tcx('smshft')
end subroutine smshft
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine pvsms1 ( ssite , sspec ,  nbas , nsp , kmax &
     , ng , gv , sv_p_orhoat , cwk , cg , job, ppnew,ppold )
  use m_struc_def
  use m_lmfinit,only:lat_alat,n0,slabl
  use m_lattic,only: lat_vol!,rv_a_opos
  use m_lgunit,only:stdo
  !- Shift in smoothed density according to job.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos pos0
  !i     Stored:    pos
  !i     Passed to: rhgcmp
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl z p pz lmxa a nr rmt nxi exi chfa rsmfa
  !i     Stored:    name
  !i     Passed to: spacks gtpcor rhgcmp
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat vol
  !i     Stored:    *
  !i     Passed to: rhgcmp
  !i   nbas  :size of basis
  !i   nsp   :number of spin channels
  !i   cy    :Normalization constants for spherical harmonics
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   orhoat:vector of offsets containing site density
  !i   job :describes which ansatz for charge shift is used for correction
  !i         <=0  do not calculate correction to force
  !i           1  shift in free-atom density
  !i           2  shift in core+nuclear density
  !i         +10  to screen the rigid shift by the Lindhard function
  !o Outputs
  !o  cg   coefficients to FT of shifted density
  !l Local variables
  !l  qloc   :difference in true and smoothed local charge.
  !l         :If the charge is assembled from overlapping atom-centered
  !l         :densities, qloc is the difference between the smoothed
  !l         :and head densities.
  !r Remarks
  !r   Shift of the "free atom densities."
  !r   The table below shows the densities and corresponding charges, and
  !r   parameters that hold their representations (true and smoothed
  !r   approximate forms):
  !r      density     charge    reps'n     smooth -reps'n
  !r      rho(smH)     qfat    cofh,ceh    already smooth
  !r      rho1-rho2    qloc     rhoat      qkl
  !r      rhoc         qc
  !r   This routine constructs the following difference
  !r   rhat(final) - rhat(initial)  positions, in the smooth reps'n, where
  !r      rhat = rho(smH) + qg * g(r)
  !r   where
  !r       qg = qval+qsc-qfat-qloc
  !r
  !r   In the special case rho is assembled from a superposition of
  !r   free-atom densities, and rho1-rho2 = rhoval(free-atm)-rho(smH)
  !r   (see ovlcor.f).  Thus in this case:
  !r      rho(free-atom) = rho(smH) + rho1-rho2 + rhoc
  !r   with the corresponding integrated charges
  !r       qval=z-qc     = qfat     + qloc      - qc
  !r   Thus in this special case qg=0: the only shift comes from rho(smH).
  !r   Because the local density (which contains the remaining part of
  !r   the free-atom density) will automatically be shifted, it follows
  !r   that the shifted smooth density will correspond to the
  !r   smooth sum-of-FA densities constructed at the shifted positions.
  !r
  !r   In the general case, qg is not zero.  By shifting the a gaussian
  !r   along with the sm-Hankels, the integrated total density of charge
  !r   shifted (local density + mesh density) is neutral.
  !r
  !r   Improvements: if the tail density were also shifted inside each
  !r   augmentation sphere, the total density would correspond exactly
  !r   to the sum-of-FA densities at the shifted positions, when the
  !r   starting density is also a sum-of-FA densities.
  !r
  !r   Shift of the "core + valence densities."
  !r
  !u Updates
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  ! ----------------------------------------------------------------------
  implicit none
  integer:: nbas , nsp , ng , job,i_copy_size
  type(s_rv1) :: sv_p_orhoat(3,1)
  real(8):: gv(ng,3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  complex(8):: cg(ng,nsp),cwk(ng)
  integer :: ib,i,is,iv0,kmax,lmxl,igetss,lmxa,nr,nlml,nxi,ie,ixi,ig,ipr,iprint,kcor,lcor
  integer,parameter:: nrmx=1501
  real(8) :: a,aa,alat,df(0:20),e,exi(n0),gam,hfc(n0,2),pnew(3),pnu(n0,2),pnz(n0,2),pold(3), &
       pp,qall,qc,qcor(2),qsc,qfat,qg,qloc,qval,rmt,rsmfa,rwgt(nrmx),scalp,summ,tpiba,v(3),&
       v2,vol,volsp,z,ppnew(3,nbas),ppold(3,nbas)
  character(35) :: strn
  character:: spid*8
  complex(8):: phase,img=(0d0,1d0)
  real(8),parameter:: pi = 4d0*datan(1d0), y0=1d0/dsqrt(4d0*pi)
  call tcn('pvsms1')
  call stdfac(20,df)
  alat=lat_alat
  vol=lat_vol
  tpiba = 2*pi/alat
  ipr = iprint()
  volsp = vol*nsp
  ! --- For each site, accumulate shift in density ---
  cg = 0d0
  iv0 = 0
  strn = 'free atom densities'
  if (job == 11) strn = 'screened free atom densities'
  if (job == 12) strn = 'screened core+multipole densities'
  if (ipr >= 30) write(stdo,339) strn
339 format(/' smshft:  add shifted ',a/'   site',16x,'old pos',22x,'new pos',14x,'shift')
  do  10  ib = 1, nbas
     is = int(ssite(ib)%spec)
     spid = slabl(is) !sspec(is)%name
     lmxl = sspec(is)%lmxl
     nlml = (lmxl+1)**2
     if (lmxl == -1) cycle
     is  = ssite(ib)%spec
     pnew= ppnew(:,ib) !ssite(ib)%pos !rv_a_opos(:,ib) !
     pold= ppold(:,ib) !ssite(ib)%pos0
     pp = alat*dsqrt(sum((pnew-pold)**2))
     if(ipr>=30) write(stdo,340) ib,spid,pold,pnew,pp/alat
340  format(i4,':',a,f8.5,2f9.5,2x,3f9.5,2x,f9.6)
     !       Skip this site if shift is negligible
     if (pp <= 1d-6) goto 18
     !   --- Shift in mesh density, job 1 ---
     if (mod(job,10)==1) then
        z=sspec(is)%z
        pnu=ssite(ib)%pnu
        pnz=ssite(ib)%pz
        lmxa=sspec(is)%lmxa
        a  =sspec(is)%a
        nr  =sspec(is)%nr
        rmt =sspec(is)%rmt
        nxi =sspec(is)%nxi
        exi =sspec(is)%exi
        hfc =sspec(is)%chfa
        rsmfa=sspec(is)%rsmfa
        gam  = 0.25d0*rsmfa**2
        call gtpcor(sspec,is,kcor,lcor,qcor)
        if (nr > nrmx) call rx('dfrce: nr gt nrmx')
        call radwgt(rmt,a,nr,rwgt)
        call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 1 , ib )%v, qloc )
        call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 2 , ib )%v, summ )
        qloc = (qloc-summ)/y0
        qfat = 0d0
        do   i  = 1, nsp
           do   ie = 1, nxi
              qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
              qfat = qfat + hfc(ie,i)*qall
           enddo
        enddo
        call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
        qg = qval+qsc-qfat-qloc !Excess sphere charge.  See Remarks above.
        !     ... Shift in smoothed free atom density
        do    i = 1, nsp
           do  ixi = 1, nxi
              e = exi(ixi)
              do ig = 1, ng
                 v(:) = gv(ig,:)*tpiba
                 v2 = sum(v**2)
                 aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
                 phase = exp(-img*alat*sum(pnew*v)) - exp(-img*alat*sum(pold*v))
                 cg(ig,i) = cg(ig,i) + hfc(ixi,i)*aa*phase*y0/vol
              enddo
           enddo
        enddo
14      continue
        !     ... Add gaussian to conserve local charge; see Remarks
        do    i = 1, nsp
           do    ig = 1, ng
              v(:) = gv(ig,:)*tpiba
              v2 = sum(v**2)
              phase = exp(-img*alat*sum(pnew*v)) - exp(-img*alat*sum(pold*v))
              cg(ig,i) = cg(ig,i) + qg*phase*dexp(-gam*v2)/volsp
           enddo
        enddo
        !   --- Shift in mesh density, job 12 ---
     elseif (job == 12) then  !     ... Core + valence at old position
        call rx('pvsms1: not supporting HAM_FORCES=12 now for simplicity')
        ! cwk=0d0
        ! ssite(ib)%pos=pold
        ! call rhgcmp(131, ib, ib, ssite, sspec ,  sv_p_orhoat, kmax , ng , cwk )
        ! cwk=-cwk !call dscal(ng*2,-1d0,cwk,1)
        ! !     ... Core + valence at new position
        ! ssite(ib)%pos=pnew
        ! call rhgcmp(131, ib, ib, ssite , sspec,  sv_p_orhoat, kmax , ng , cwk )
        ! do i=1,nsp
        !    cg(:,i)= cwk/nsp + cg(:,i)
        ! enddo
     else
        call rxi('smshft: bad job:',job)
     endif
18   continue
     iv0 = iv0+nlml
10 enddo
  call tcx('pvsms1')
end subroutine pvsms1


