      subroutine smshft (job)
      use m_density,only:  sv_p_orhoat=>orhoat,smrho=>osmrho !input and output
      use m_supot,only:  iv_a_okv,rv_a_ogv
      use m_lmfinit,only: ctrl_lfrce,ham_elind,
     &     lat_alat,nsp,nbas, sspec=>v_sspec, ssite=>v_ssite
      use m_lattic,only: lat_qlat, lat_vol,lat_plat
      use m_supot,only: lat_nabc, lat_ng
      use m_lgunit,only:stdo
C- Estimate the smooth density for a shift in atomic positions.
!input  orhoat:vector of offsets containing site density
!inputoutput 
!   smrho :a perturbation is added to smrho, depending on job
Cr Remarks
Cr   job describes which ansatz for charge shift is used for correction
Cr     <=0  do not calculate correction to force
Cr       1  shift in free-atom density
Cr       2  shift in core+nuclear density
Cr     +10  to screen the rigid shift by the Lindhard function
Cu Updates
Cu   17 Sep 01 Adapted for local orbitals.  Altered argument list
Cu    3 Jul 00 Adapted from nfp smshft.f
C ----------------------------------------------------------------------
      implicit none
      integer i,ib,igets,igetss,iprint,is,k1,k2,k3,kcor,lcor,lmxa,n0,n1,n2,n3,ng,ngabc(3),i_copy_size
      complex(8) ,allocatable :: cgr_zv(:,:)
      complex(8) ,allocatable :: cgs_zv(:,:),cgs_zvv(:,:)
      complex(8) ,allocatable :: cwk_zv(:)
      integer kmax,job
      parameter (n0=10)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision alat,dgets,elind,pi,qc,qsc,qcor(2),plat(3,3),
     .     qlat(3,3),qv,qval,tpiba,vol,z,pnu(n0,2),pnz(n0,2)
      if (job .le. 0) return
      call tcn('smshft')
      alat=lat_alat
      plat=lat_plat
      qlat=lat_qlat
      ngabc=lat_nabc
      ng=lat_ng
      vol=lat_vol
      call fftz30(n1,n2,n3,k1,k2,k3)
C ... Hold on to original smrho (cgr)
      allocate(cgr_zv(ng,nsp))
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,-1)
      call gvgetf(ng , nsp , iv_a_okv , k1 , k2 , k3 , smrho , cgr_zv  )
C --- Shift in unscreened density at the two positions ---
      allocate(cgs_zv(ng,nsp),cgs_zvv(ng,nsp),cwk_zv(ng))
      kmax = 0
      call pvsms1 ( ssite , sspec ,  nbas , nsp , kmax , ng ,!slat ,
     .  rv_a_ogv , sv_p_orhoat , cwk_zv , cgs_zv , job )
      if (allocated(cwk_zv)) deallocate(cwk_zv)
C --- Screened shift ---
      if (job .gt. 10) then !Compute elind if not given
        qval = 0d0
        do  12  ib = 1, nbas
          is  =ssite(ib)%spec
          z   =sspec(is)%z
          pnu =sspec(is)%p
          pnz =sspec(is)%pz
          lmxa=sspec(is)%lmxa
          if (lmxa .eq. -1) goto 12
          call gtpcor(sspec,is,kcor,lcor,qcor)
          call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qv,qsc)
          qval = qval+qv
   12   continue
        pi = 4d0*datan(1d0)
        tpiba = 2*pi/alat
        elind = ham_elind
        if (elind .lt. 0d0) elind = -(3*pi**2*qval/vol)**.66666d0*elind
        cgs_zvv(:,1)= cgs_zv(:,1) + cgs_zv(:,2)
        cgs_zvv(:,2)= cgs_zv(:,1) - cgs_zv(:,2)
        call lindsc( 2 , ng , rv_a_ogv , tpiba , elind , cgs_zvv )
      endif
C ... Debugging: show delta smrho
      call gvputf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , cgs_zv , smrho )
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C     --- Add shift to smrho, ensuring no shift in <rho> ---
      cgs_zv(1,1:nsp)=0d0 !real part zero bug? 2022-5-27
c    do  i = 1, nsp     
c      call dvset(cgs_zv,1 + ng*( i - 1 ), 1 + ng*(i-1 ) , 0d0 )<=== real part zero bug? 2022-5-27
c    enddo
      cgr_zv = cgr_zv + cgs_zv 
      call gvputf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , cgr_zv , smrho )
      call fftz3(smrho,n1,n2,n3,k1,k2,k3,nsp,0,1)
C --- Symmetrize the shifted density ---
      call symsmrho(smrho) 
      if (allocated(cgs_zv)) deallocate(cgs_zv,cgs_zvv)
      if (allocated(cgr_zv)) deallocate(cgr_zv)
      call tcx('smshft')
      end subroutine smshft
!!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine pvsms1 ( ssite , sspec ,  nbas , nsp , kmax 
     . , ng , gv , sv_p_orhoat , cwk , cg , job )
      use m_struc_def          
      use m_lmfinit,only:lat_alat
      use m_lattic,only: lat_vol
      use m_lgunit,only:stdo
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
      implicit none
      integer:: nbas , nsp , ng , job,i_copy_size
      type(s_rv1) :: sv_p_orhoat(3,1)
      real(8):: gv(ng,3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      double complex cg(ng,nsp),cwk(ng)
      integer ib,i,is,iv0,kmax,lmxl,igetss,lmxa,nr,nrmx,n0,
     .nlml,nxi,ie,ixi,ig,ipr,iprint,kcor,lcor
      parameter (nrmx=1501, n0=10)
      double precision a,aa,alat,df(0:20),e,
     .exi(n0),gam,hfc(n0,2),pi,pnew(3),pnu(n0,2),pnz(n0,2),pold(3),
     .pp,qall,qc,qcor(2),qsc,qfat,qg,qloc,qval,rmt,
     .rsmfa,rwgt(nrmx),scalp,summ,tpiba,v(3),v2,vol,volsp,y0,z
      character*35 strn,spid*8
      double complex phase
      complex(8):: img=(0d0,1d0)
      call tcn('pvsms1')
      call stdfac(20,df)
      alat=lat_alat
      vol=lat_vol
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      tpiba = 2*pi/alat
      ipr = iprint()
      volsp = vol*nsp
C --- For each site, accumulate shift in density ---
      cg = 0d0
      iv0 = 0
      strn = 'free atom densities'
      if (job .eq. 11) strn = 'screened free atom densities'
      if (job .eq. 12) strn = 'screened core+multipole densities'
      if (ipr .ge. 30) write(stdo,339) strn
  339 format(/' smshft:  add shifted ',a/'   site',16x,'old pos',22x,'new pos',14x,'shift')
      do  10  ib = 1, nbas
        is = int(ssite(ib)%spec)
        spid = sspec(is)%name
        lmxl = sspec(is)%lmxl
        nlml = (lmxl+1)**2
        if (lmxl .eq. -1) goto 10
        is  = ssite(ib)%spec
        pnew= ssite(ib)%pos
        pold= ssite(ib)%pos0
        pp = alat*dsqrt((pnew(1)-pold(1))**2 + (pnew(2)-pold(2))**2+(pnew(3)-pold(3))**2)
        if(ipr>=30) write(stdo,340) ib,spid,pold,pnew,pp/alat
  340   format(i4,':',a,f8.5,2f9.5,2x,3f9.5,2x,f9.6)
C       Skip this site if shift is negligible
        if (pp .le. 1d-6) goto 18
C   --- Shift in mesh density, job 1 ---
        if (mod(job,10)==1) then
          z=sspec(is)%z
          i_copy_size=size(sspec(is)%p)
          call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
          i_copy_size=size(sspec(is)%pz)
          call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)
          lmxa=sspec(is)%lmxa
          a  =sspec(is)%a
          nr  =sspec(is)%nr
          rmt =sspec(is)%rmt
          nxi =sspec(is)%nxi
          exi =sspec(is)%exi
          hfc =sspec(is)%chfa
          rsmfa=sspec(is)%rsmfa
          call gtpcor(sspec,is,kcor,lcor,qcor)
          if (nr .gt. nrmx) call rx('dfrce: nr gt nrmx')
          call radwgt(rmt,a,nr,rwgt)
          call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 1 , ib )%v, qloc )
          call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 2 , ib )%v, summ )
          qloc = (qloc-summ)/y0
          qfat = 0d0
          do   i  = 1, nsp
          do   ie = 1, nxi
            gam  = 0.25d0*rsmfa**2
            qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
            qfat = qfat + hfc(ie,i)*qall
         enddo
         enddo
         call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
          qg = qval+qsc-qfat-qloc !Excess sphere charge.  See Remarks above.
C     ... Shift in smoothed free atom density
          do    i = 1, nsp
          do  ixi = 1, nxi
            e = exi(ixi)
            do ig = 1, ng
              v(:) = gv(ig,:)*tpiba
              v2 = sum(v**2) !v(1)**2+v(2)**2+v(3)**2
              aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
c              scalp = -alat*sum(pnew*v) !(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
c              phase = dcmplx(dcos(scalp),dsin(scalp))
c              scalp = -alat*sum(pold*v) !(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
c              phase = phase - dcmplx(dcos(scalp),dsin(scalp))
              phase = exp(-img*alat*sum(pnew*v)) - exp(-img*alat*sum(pold*v))
              cg(ig,i) = cg(ig,i) + hfc(ixi,i)*aa*phase*y0/vol
          enddo
          enddo
          enddo
   14     continue
C     ... Add gaussian to conserve local charge; see Remarks
          do    i = 1, nsp
          do    ig = 1, ng
            v(:) = gv(ig,:)*tpiba
            v2 = sum(v**2) !v(1)**2+v(2)**2+v(3)**2
c            scalp = -alat*sum(pnew*v) !(pnew(1)*v(1)+pnew(2)*v(2)+pnew(3)*v(3))
c            phase = dcmplx(dcos(scalp),dsin(scalp))
c            scalp = -alat*sum(pold*v) !(pold(1)*v(1)+pold(2)*v(2)+pold(3)*v(3))
!     phase - dcmplx(dcos(scalp),dsin(scalp))
            phase = exp(-img*alat*sum(pnew*v)) - exp(-img*alat*sum(pold*v))
            cg(ig,i) = cg(ig,i) + qg*phase*dexp(-gam*v2)/volsp
          enddo
          enddo
C   --- Shift in mesh density, job 12 ---
        elseif (job .eq. 12) then
C     ... Core + valence at old position
          cwk=0d0
          ssite(ib)%pos=pold
          call rhgcmp(131, ib, ib, ssite, sspec ,  sv_p_orhoat, kmax , ng , cwk )
          cwk=-cwk !call dscal(ng*2,-1d0,cwk,1)
C     ... Core + valence at new position
          ssite(ib)%pos=pnew
          call rhgcmp(131, ib, ib, ssite , sspec,  sv_p_orhoat, kmax , ng , cwk )
C     ... Add to cg
          cg(:,i)= cwk/nsp + cg(:,i)
c          do  i = 1, nsp
c            call daxpy(ng*2,1d0/nsp,cwk,1,cg(1,i),1)
c          enddo
        else
          call rxi('smshft: bad job:',job)
        endif
   18   iv0 = iv0+nlml
   10 continue
      call tcx('pvsms1')
      end subroutine pvsms1


