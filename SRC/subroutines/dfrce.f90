module m_dfrce
  public dfrce
  contains
subroutine dfrce (job, sv_p_orhoat , sv_p_orhat1 ,  qmom , smrho , smrout , dfh )
  use m_lmfinit,only: nvl=>pot_nlml
  use m_supot,only: rv_a_ogv,iv_a_okv
  use m_struc_def
  use m_lmfinit,only:lat_alat,nsp,nbas,nspec,ssite=>v_ssite,sspec=>v_sspec
  use m_lattic,only: lat_qlat, lat_vol,lat_plat
  use m_supot,only: lat_nabc,k1,k2,k3,lat_ng
  use m_lgunit,only:stdo
  ! Correction to force theorem, Harris functional
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: pvdf4 pvdf2 rhomom pvdf1 smvxcm
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl z p pz lmxa a nr rmt nxi exi chfa rsmfa rg
  !i     Stored:    *
  !i     Passed to: pvdf4 pvdf2 rhomom pvdf1 gtpcor corprm smvxcm
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: nabc ng ogv okv vol alat plat qlat
  !i     Stored:    *
  !i     Passed to: pvdf4 pvdf2 pvdf1 smvxcm
  !i   sctrl :struct for program flow parameters; see routine uctrl
  !i     Elts read: lfrce
  !i     Stored:    *
  !i     Passed to: *
  !i   k1..3 :dimensions smrho
  !i   nvl   :sum of local (lmxl+1)**2, lmxl = density l-cutoff
  !i   orhoat:vector of offsets containing site density
  !i   orhat1:pointer to local densities
  !i   elind :Lindhard parameter, used for Lindhard screening
  !i   qmom  :multipole moments of on-site densities (rhomom.f)
  !i   smrho :smooth density on uniform mesh
  !i   smrho :smooth (input) density that generated the hamiltonian
  !i   smrout:smooth (output) density that the hamiltonian generated
  !o Outputs
  !o   dfh   :correction to the HF force
  !l Local variables
  !l    job  :describes which ansatz for charge shift is used for correction
  !l         :<=0  do not calculate correction to force
  !l         :  1  shift in free-atom density
  !l         :  2  shift in core+nuclear density
  !l         :+10  to screen the rigid shift by the Lindhard function
  !r Remarks
  !! Density
  !!  sv_p_orhoat: input atomic density that generaed Hamiltonian
  !!  sv_p_orhat1: new atomic density that the Hamiltonian generated
  !!  smrho:  input  density that generaed Hamiltonian
  !!  smrout: output density that the Hamiltonian generated
  !u Updates
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   18 Dec 03 adapted to modified smvxc
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   17 Sep 01 Adapted for local orbitals
  !u   21 Jun 00 spin polarized
  !u   18 Jun 98 adapted from nfp dfrce.f
  !u   16 Jun 98 MvS parallelized for SGI
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  type(s_rv1) :: sv_p_orhat1(3,*)
  type(s_rv1) :: sv_p_orhoat(3,*)
  real(8):: dfh(3,nbas) , qmom(*)
  double complex smrho(k1,k2,k3,*),smrout(k1,k2,k3,*)
  integer :: job,n1,n2,n3,ng,iprint,ib,is,lmxl,iv0,nlm, &
       ip,m,i,ngabc(3),ltop,nlmtop,igets,igetss,nn
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
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
  double precision :: vol,plat(3,3),qlat(3,3),alat,vsum,pi,tpiba,elind=0d0, &
       fes1(3),fes2(3),fxc(3),c,avgdf(3)
  integer::ibini,ibend,iiv0(nbas)
  character(40) :: strn
  real(8),allocatable:: cs_(:),sn_(:)
  call tcn('dfrce')
  ! --- Setup ---
  dfh=0d0
  ngabc=lat_nabc
  ng=lat_ng
  vol=lat_vol
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
  c = 1000
  nn   = k1*k2*k3
  ! ... Arrays needed for pvdf1
  allocate(ceps_zv(ng))
  allocate(cnomi_zv(ng))
  allocate(cvin_zv(ng))
  allocate(cdvx_zv(ng*nsp))
  ! ... Set up for vectorized Y_lm and gaussians
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
  ! --- Make ves(rhoin,q) ---
  allocate(smro_zv(nn))
  allocate(cs_(ng), sn_(ng))
  call dcopy(2*nn, smrho,1,smro_zv,1)
  if(nsp==2) call daxpy(2*nn, 1d0, smrho (1,1,1,2), 1, smro_zv, 1)
  call fftz3 ( smro_zv , n1 , n2 , n3 , k1 , k2 , k3 , 1 , 0 , - 1 )
  call gvgetf ( ng , 1 , iv_a_okv , k1 , k2 , k3 , smro_zv , cvin_zv )
  call pvdf4 ( ssite , sspec ,  qmom , ng , g2_rv , yl_rv , cs_ , sn_ , iv_iv , qlat , cvin_zv )
  if (allocated(smro_zv)) deallocate(smro_zv)
  deallocate(cs_,sn_)
  ! --- Make dVxc(in)/dn ---
  allocate(dvxc_zv(nn*nsp))
  allocate(smro_zv(nn*nsp))
  allocate(vxcp_zv(nn*nsp))
  allocate(vxcm_zv(nn*nsp))
  allocate(wk1_zv(nn*nsp))
  allocate(wk2_zv(nn*nsp))
  allocate(wk3_zv(nn*nsp))
  call dpcopy ( smrho , smro_zv , 1 , 2 * nn * nsp , 1d0 )
  call pvdf2 ( nbas , nsp , ssite , sspec ,  n1 , n2 , n3 & ! & slat ,
  , k1 , k2 , k3 , smro_zv , vxcp_zv , vxcm_zv , wk1_zv &
       , wk2_zv , wk3_zv , dvxc_zv )
  deallocate(wk3_zv)
  deallocate(wk2_zv)
  deallocate(wk1_zv)
  deallocate(vxcm_zv)
  deallocate(vxcp_zv)
  ! --- cdvx = FFT ((n0_out-n0_in) dVxc/dn) ---
  !     Use total n0_out-n0_in but keep vxc spin polarized
  call dpzero ( smro_zv , 2 * nn )
  do  i = 1, nsp
     call daxpy (2*nn,  1d0, smrout( 1 , 1 , 1 , i ), 1, smro_zv,1)
     call daxpy (2*nn, -1d0, smrho ( 1 , 1 , 1 , i ), 1, smro_zv,1)
  enddo
  call pvdf3 ( n1 , n2 , n3 , k1 , k2 , k3 , nsp , smro_zv , dvxc_zv )
  call fftz3 ( dvxc_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp , 0 , - 1 )
  call gvgetf ( ng , nsp , iv_a_okv , k1 , k2 , k3 , dvxc_zv , cdvx_zv )
  ! --- Cnomi = (n0_out(q) - n0_in(q)) ---
  call fftz3 ( smro_zv , n1 , n2 , n3 , k1 , k2 , k3 , 1 , 0 , - 1 )
  call gvgetf ( ng , 1 , iv_a_okv , k1 , k2 , k3 , smro_zv , cnomi_zv  )
  if (allocated(smro_zv)) deallocate(smro_zv)
  if (allocated(dvxc_zv)) deallocate(dvxc_zv)
  ! ... Debugging slot smrho(out) for out-in
  !      print *, '*** debugging ... subs smrout for out-in'
  !      call dpcopy(smrout,w(osmro),1,2*nn,1d0)
  !      call fftz3(w(osmro),n1,n2,n3,k1,k2,k3,1,0,-1)
  !      call gvgetf(ng,1,w(okv),k1,k2,k3,w(osmro),w(ocnomi))
  !      call zprm3('rho-out(q)',w(osmro),k1,k2,k3)

  ! --- Multipole moments of the output density ---
  allocate(qmout_rv(nvl))
  call pshpr(0)
  call rhomom (sv_p_orhat1 , qmout_rv , vsum )
  call poppr
  qmout_rv = qmout_rv- qmom(1:nvl)
  ! --- Lindhard dielectric function ---
  if (job > 10) then
     pi = 4d0*datan(1d0)
     tpiba = 2*pi/alat
     call lindsc ( 3 , ng , rv_a_ogv , tpiba , elind , ceps_zv )
  endif
  ! --- For each site, get correction to force ---
  if (iprint() >= 30) then
     strn = 'shift in free-atom density'
     if (job == 11) strn = 'screened shift in free-atom density'
     if (job == 12) strn = 'screened shift in core+nuclear density'
     write(stdo,201) strn
  endif
201 format(/' Harris correction to forces: ',a/ &
       '  ib',9x,'delta-n dVes',13x,'delta-n dVxc',15x,'total')
  ! ... Setup array iv0 offset of qmom to ibas
  iv0 = 0
  do ib = 1, nbas
     is = ssite(ib)%spec
     lmxl = sspec(is)%lmxl
     nlm = (lmxl+1)**2
     iiv0(ib) = iv0
     iv0 = iv0+nlm
  enddo
  ibini=1
  ibend=nbas
  do ib = ibini, ibend
     is   = ssite(ib)%spec
     lmxl = sspec(is)%lmxl
     if (lmxl == -1) goto 20
     nlm = (lmxl+1)**2
     call pvdf1 ( job , ssite , sspec , nsp , ib , iiv0(ib), qmom &
          , qmout_rv , ng , rv_a_ogv , g2_rv , yl_rv , iv_iv , qlat , 0 &
          , cnomi_zv , ceps_zv , cdvx_zv , cvin_zv , sv_p_orhoat ( 1 , &
          ib ) , fes1 , fes2 , fxc )
     do  i = 1, 3
        dfh(i,ib) = -(fes1(i) + fes2(i) + fxc(i))
     enddo
     ! if ! (MPI | MPIK)
     if (iprint() >= 30) &
          write(stdo,200) ib,(c*(fes1(m)+fes2(m)),m=1,3), &
          (c*fxc(m),m=1,3),(c*dfh(m,ib),m=1,3)
200  format(i4,3f8.2,1x,3f8.2,1x,3f8.2:1x,3f8.2)
     ! endif
20   continue
  enddo
  avgdf=0d0
  do  ib = 1, nbas
     do   i = 1, 3
        avgdf(i) = avgdf(i) + dfh(i,ib)/nbas
     enddo
  enddo
  ! ... Shift all forces to make avg correction zero
  do  ib = 1, nbas
     do  i = 1, 3
        dfh(i,ib) = dfh(i,ib) - avgdf(i)
     enddo
  enddo
  if (iprint() >= 30) write(stdo,331) (c*avgdf(m),m=1,3)
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

subroutine pvdf1 ( job , ssite , sspec ,  nsp , ib , iv0 & ! & slat ,
  , qmom , qmout , ng , gv , g2 , yl , iv , qlat , kmax , cnomin &
       , ceps , cdvxc , cvin , sv_p_orhoat , fes1 , fes2 , fxc )
  use m_struc_def  !Cgetarg
  use m_lmfinit,only: nbas
  use m_lmfinit,only:lat_alat,pnux=>pnu,pzx=>pz
  use m_lattic,only: lat_vol,rv_a_opos
  use m_supot,only: lat_nabc
  ! need to modify texts.
  !- Estimate shift in local density for one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite,sspec,slat
  !i   ng,gv,kmax
  !i   orhoat
  !i   job: 1  shift in free-atom density
  !i        2  shift in core+nuclear density
  !i      +10  to screen the rigid shift by the response function
  !i   ib      which site is being shifted
  !i   iv0     offset to qmom
  !i   qmom,qmout moments of input and output densities
  !i   cnomin  difference betw. smoothed output and input density n0
  !i   cvin    electrostatic potential of input density Ves[n0~_in]
  !i   ceps    response function
  !i   cdvxc   dVxc/dn (nout-nin)
  !o Outputs
  !o   cdn0:   Job 1:  shift in valence part of the free atom density
  !o           Job 12: shift in atom density (1/eps - 1)
  !o   cdn:    Job 1:  dn^(u) where dn is the unscreened shift in
  !o           in the free-atom density.
  !o           Job 12: dn^(u) 1/eps where dn is unscreened shift in
  !o           the charge density.  Local density approximated
  !o   NB:     In all cases, the local part of density is approximated
  !o           by a gaussian of the equivalent multipole moment.
  !o   cdv:    shift in the electrostatic potential
  !o   fes1,fes2,fxc
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer:: ng , nsp , iv0 , kmax , ib , job , iv(ng,3),i_copy_size
  type(s_rv1) :: sv_p_orhoat(3)

  real(8):: qmom(*) , qmout(*) , gv(ng,3) , tau(3) , fes1(3) , &
       fes2(3) , fxc(3) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  !      type(s_lat)::slat

  double complex cdn0(ng,nsp),cdn(ng),cdv(ng),ceps(ng), &
       cnomin(ng),cdvxc(ng,nsp),cvin(ng)
  ! ... Local parameters
  integer :: ig,ilm,l,lmxl,m,nlm,nlmx,k,is,jv0,jb,js,ll,n0, &
       nrmx
  parameter (nlmx=64, nrmx=1501, n0=10)
  integer :: lmxa,nr,nxi,ie,ixi,job0,kcor,lcor,lfoc,i, &
       ngabc(3),n1,n2,n3,nlml
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  double precision :: pi,alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg, &
       vol,y0,z,v(3),df(0:20),feso(3),qcor(2),gpot0(nlmx,3),fesdn(3), &
       fesgg(3),pnu(n0),pnz(n0),a,rmt,qloc,rsmfa,exi(n0),hfc(n0,2), &
       qfat,gam,qall,qc,qval,qg,e,aa,q0(3),sum
  double precision :: rwgt(nrmx),cc,gamf,cfoc,cvol
  !     parameter (k0=3)
  !     double complex gkl(0:k0,nlmx)
  double complex tpia,cxx,phase,gc0,xc0,cof(nlmx)
  ! ... Heap
  data q0 /0d0,0d0,0d0/
  call tcn('pvdf1')
  ngabc=lat_nabc
  alat=lat_alat
  vol=lat_vol
  call stdfac(20,df)
  pi = 4d0*datan(1d0)
  y0 = 1d0/dsqrt(4d0*pi)
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
  tau=rv_a_opos(:,ib) !ssite(ib)%pos(1:3)
  call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
  ! --- Unscreened rigid charge density shift, job 1, in cdn0 ---
  if (job0 == 1) then
     z=sspec(is)%z
     pnu=pnux(1:n0,1,is) !sspec(is)%p
     pnz=pzx(1:n0,1,is)  !sspec(is)%pz
     lmxa=sspec(is)%lmxa
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     lmxl=sspec(is)%lmxl
     nxi=sspec(is)%nxi
     exi=sspec(is)%exi
     hfc=sspec(is)%chfa
     rsmfa=sspec(is)%rsmfa
     gam  = 0.25d0*rsmfa**2
     call gtpcor(sspec,is,kcor,lcor,qcor)
     if (nr > nrmx) call rx('dfrce: nr gt nrmx')
     call radwgt(rmt,a,nr,rwgt)
     nlml = (lmxl+1)**2
     call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 1 )%v , qloc )
     call radsum ( nr , nr , nlml , nsp , rwgt , sv_p_orhoat( 2 )%v , sum )
     qloc = (qloc-sum)/y0
     qfat = 0d0
     do  i  = 1, nsp
        do  ie = 1, nxi
           qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
           qfat = qfat + hfc(ie,i)*qall
        enddo
     enddo
     call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
     qg = qval+qsc-qfat-qloc
     !   ... Shift in free atom density
     do  141    i = 1, nsp
        do  142  ixi = 1, nxi
           e = exi(ixi)
           cc = -4d0*pi*hfc(ixi,i)*y0/vol
           do  15  ig = 1, ng
              aa = cc*dexp(gam*(e-g2(ig)))/(e-g2(ig))
              cdn0(ig,i) = cdn0(ig,i) + aa*dcmplx(cs(ig),sn(ig))
15         enddo
142     enddo
141  enddo
     !   ... Add gaussian to conserve local charge
     !     ... Add gaussian to conserve local charge.  If density corresponds
     !         to the free-atom density, qfat+qloc = qval+qsc; then qg=0
     cc = qg/vol/nsp
     do   i = 1, nsp
        do   ig = 1, ng
           cdn0(ig,i)=cdn0(ig,i)+cc*dcmplx(cs(ig),sn(ig))*dexp(-gam*g2(ig))
        enddo
     enddo
  endif
  ! --- Coefficients defining local valence + core density ---
  is=ssite(ib)%spec
  tau=rv_a_opos(:,ib) !ssite(ib)%pos
  lmxl=sspec(is)%lmxl
  rg=sspec(is)%rg
  call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
  nlm = (lmxl+1)**2
  if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
  ilm = 0
  cxx = dcmplx(0d0,1d0)
  do  20  l = 0, lmxl
     cxx = cxx*dcmplx(0d0,-1d0)
     do  22  m = -l,l
        ilm = ilm+1
        cof(ilm) = cxx*qmom(ilm+iv0)*4*pi/df(2*l+1)
22   enddo
20 enddo
  !     cof(1) = cof(1) + 4*pi*y0*(qcorg+qsc-z)
  cof(1) = cof(1) + 4*pi*y0*(qcorg-z)
  ! --- Shift in n0, ves~ for list of G vectors ---
  gam = 0.25d0*rg*rg
  gamf = 0.25d0*rfoc*rfoc
  cfoc = -4d0*pi*y0*cofh/vol
  cvol = 1d0/vol
  do  30  ig = 2, ng
     v = gv(ig,:)
     !   ... Accumulate unscreened smoothed core+nuclear density
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
33      enddo
32   enddo
     !   ... Accumulate unscreened foca density
     aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
     cdn(ig) = cdn(ig) + cfoc*aa*phase
     !   ... Make the screened shift in input density n0~
     !       Job 1: cdn0 = (valence part of) cdn^u ; cdn = cdn^u
     if (job0 == 1) then
        cdn(ig) = cdn(ig) + (cdn0(ig,1) + cdn0(ig,nsp))/(3-nsp)
        if (job > 10) cdn(ig) = cdn(ig) / ceps(ig)
        !       Job 12: cdn0 = cdn^u (1/eps - 1); cdn = cdn^s = cdn^u / eps
     elseif (job == 12) then
        do  i = 1, nsp
           cdn0(ig,i) = cdn(ig) * (1/ceps(ig)-1) / nsp
        enddo
        cdn(ig) = cdn(ig) / ceps(ig)
     else
        call rxi('dfrce: nonsensical job',job)
     endif
     !   ... Electrostatic potential shift = 1/eps dv [n0~]
     !       g2 = tpiba*tpiba*(gv(ig,1)**2+gv(ig,2)**2+gv(ig,3)**2)
     cdv(ig)  = cdn(ig) * (8*pi/g2(ig))
     do  36  k = 1, 3
        fes1(k) = fes1(k) + dconjg(cnomin(ig)) * tpia*v(k)*cdv(ig)
        do  i = 1, nsp
           fxc(k)  = fxc(k)  + dconjg(cdvxc(ig,i)) * tpia*v(k)*cdn0(ig,i)
        enddo
36   enddo
30 enddo
  do  37  k = 1, 3
     fxc(k)  = fxc(k)*vol
     fes1(k) = fes1(k)*vol
37 enddo
  ! --- Integral of grad g (output-input local charge) ves~ ---
  fesgg=0d0 !call dpzero(fesgg,3)
  do  k = 1, 3
     do  ilm = 1, nlm
        l = ll(ilm)
        gpot0(ilm,k) = gpot0(ilm,k)*4*pi/df(2*l+1)
        fesgg(k) = fesgg(k) + qmout(iv0+ilm)*gpot0(ilm,k)
     enddo
  enddo
  !      print 339, 'n0(out-in) * g dves ',fes1
  !      print 339, 'd(g) qmom(out-in) ves[n0~]',fesgg
  !      print 339, 'n0~(out-in) * dvxc   ',fxc
  !  339 format(a,6p,3f8.2)

  ! --- Integral of dves~ (output-input local charge) for all sites ---
  call dpzero(fes2,3)
  call dpzero(feso,3)
  jv0 = 0
  do  40  jb = 1, nbas
     js=ssite(jb)%spec
     tau=rv_a_opos(:,jb) !ssite(jb)%pos
     lmxl=sspec(js)%lmxl
     rg=sspec(js)%rg
     nlm = (lmxl+1)**2
     ! ... For this jb, mesh density for all G vectors
     if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
     call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
     call dpzero(gpot0,nlmx*3)
     gam = 0.25d0*rg*rg
     do  50  ig = 2, ng
        aa = dexp(-gam*g2(ig))
        gc0 = dcmplx(0d0,1d0)*aa* &
             dconjg(tpia*cdv(ig))*dcmplx(cs(ig),sn(ig))
        ilm = 0
        do  55  l = 0, lmxl
           gc0 = gc0*dcmplx(0d0,-1d0)
           do  56  m = -l,l
              ilm = ilm+1
              gpot0(ilm,1) = gpot0(ilm,1)+dble(gc0)*yl(ig,ilm)*gv(ig,1)
              gpot0(ilm,2) = gpot0(ilm,2)+dble(gc0)*yl(ig,ilm)*gv(ig,2)
              gpot0(ilm,3) = gpot0(ilm,3)+dble(gc0)*yl(ig,ilm)*gv(ig,3)
56         enddo
55      enddo
50   enddo
     !   ... Multiply factors into gpot0, accumulate force
     ilm = 0
     do  60  l = 0, lmxl
        do  62  m = -l, l
           ilm = ilm+1
           do  64  k = 1, 3
              gpot0(ilm,k) = gpot0(ilm,k)*4*pi/df(2*l+1)
              feso(k) = feso(k) + qmom(jv0+ilm)*gpot0(ilm,k)
              fes2(k) = fes2(k) + qmout(jv0+ilm)*gpot0(ilm,k)
64         enddo
62      enddo
60   enddo
     jv0 = jv0+nlm
40 enddo
  fes2=fes2-fesgg
  call tcx('pvdf1')
end subroutine pvdf1


subroutine pvdf2(nbas,nsp,ssite,sspec,n1,n2,n3,k1,k2,k3, smrho,vxcp,vxcm,wk1,wk2,wk3,dvxc)
  use m_struc_def
  use m_smvxcm,only: smvxcm
  !- Makes derivative of smoothed xc potential wrt density.
  implicit none
  ! ... Passed parameters
  integer :: nbas,nsp,n1,n2,n3,k1,k2,k3
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  complex(8):: vxcp(k1,k2,k3,nsp),vxcm(k1,k2,k3,nsp), &
       dvxc(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp), &
       wk1(k1,k2,k3,nsp),wk2(k1,k2,k3,nsp), &
       wk3(k1,k2,k3,nsp)
  ! ... Local parameters
  integer :: i1,i2,i3,i,nn
  double precision :: fac,dmach,f1,f2,f,alfa,dfdr,rrho,dvdr, &
       rmusm(2),rvmusm(2),rvepsm(2),repsm(2),repsmx(2),repsmc(2), &
       fcexc0(2),fcex0(2),fcec0(2),fcvxc0(2),ff(1,1)
  fac = dmach(1)**(1d0/3d0)
  alfa = 2d0/3d0
  nn = k1*k2*k3
  call pshpr(0)
  ! ... Add fac (rho+ + rho-)/2 into rho+, rho- for spin pol case,
  !     Add fac * rho into rho if not spin polarized
  if (nsp == 1) then
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
  ! ... vxcp = vxc (smrho+drho)
  call dpzero(vxcp, nn*2*nsp)
  call dpzero(wk1, nn*2*nsp)
  call dpzero(wk2, nn*2*nsp)
  call dpzero(wk3, nn*2)
  call smvxcm(ssite,sspec,nbas,0,k1,k2,k3,smrho, &
  vxcp,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm, &
       rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,ff)
  ! ... Replace fac*rho with -fac*rho
  if (nsp == 1) then
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
  ! ... vxcm = vxc (smrho-drho)
  call dpzero(vxcm, nn*2*nsp)
  call smvxcm(ssite,sspec,nbas,0,k1,k2,k3,smrho,&
  vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm, &
       rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,ff)
  ! ... Restore rho+, rho-
  if (nsp == 1) then
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
  ! ... Overwrite vxcp with df/drho
  do i = 1, nsp
     do i3=1,n3
        do i2=1,n2
           do i1=1,n1
              rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
              if (rrho > 0) then
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
  ! ... vxcm = vxc (smrho)
  call dpzero(vxcm, nn*2*nsp)
  call smvxcm(ssite,sspec,nbas,0,k1,k2,k3,smrho,&
       vxcm,dvxc,wk1,wk2,wk3,repsm,repsmx,repsmc,rmusm, &
       rvmusm,rvepsm,fcexc0,fcex0,fcec0,fcvxc0,ff)
  ! ... dvxc/drho into dvxc
  do  i = 1, nsp
     do i3=1,n3
        do i2=1,n2
           do i1=1,n1
              rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
              if (rrho > 0) then
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
  call poppr
end subroutine pvdf2

subroutine pvdf3(n1,n2,n3,k1,k2,k3,nsp,deln0,dvxc)
  !- Overwrites dvxc with (nout-nin)*dvxc
  implicit none
  integer :: n1,n2,n3,k1,k2,k3,nsp
  double complex deln0(k1,k2,k3),dvxc(k1,k2,k3,nsp)
  integer :: i1,i2,i3,i
  do  i  = 1, nsp
     do  i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              dvxc(i1,i2,i3,i) = dvxc(i1,i2,i3,i)*deln0(i1,i2,i3)
           enddo
        enddo
     enddo
  enddo
end subroutine pvdf3

subroutine pvdf4(ssite,sspec,qmom,ng,g2,yl,cs,sn,iv,qlat,cv) !slat,
  use m_struc_def
  use m_lmfinit,only: nbas
  use m_lattic,only: lat_vol,rv_a_opos
  use m_supot,only: lat_nabc
  !- Makes smoothed ves from smoothed density and qmom, incl nuc. charge
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl rg
  !i     Stored:    *
  !i     Passed to: corprm
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: nabc vol
  !i     Stored:    *
  !i     Passed to: *
  !i   qmom  :multipole moments of on-site densities (rhomom.f)
  !i   ng    :number of G-vectors
  !i   g2    :square of G-vectors
  !i   yl    :spherical harmonics
  !i   cs    :vector of cosines for the ng vectors
  !i   sn    :vector of sines for the ng vectors
  !i   iv
  !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
  !o Outputs
  !o   cv    :local gaussian density added to cv
  !o         :estatatic potential make from density
  !r Remarks
  !r   Local charge consists of a sum of gaussians that compensate for
  !r   the difference in multipole moments of true and smooth local charge
  !r   and a contribution from the smooth core charge.
  !r     g(qmpol) + g(qcore-z) + h(ncore)
  !r
  !r   Adapted from vesgcm to make strictly FT ves(nloc)
  !u Updates
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ng,iv(ng,3),i_copy_size
  real(8):: qmom(1) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex cv(ng)
  integer :: ig,ib,ilm,is,iv0,l,lmxl,m,nlm,nlmx,nglob,n1,n2,n3, ngabc(3),lfoc
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  parameter (nlmx=64)
  double precision :: tau(3),df(0:20),pi,y0,vol,rg,qcorg,qcorh,qsc, &
       cofg,cofh,ceh,rfoc,z,q0(3),gam,gamf,cfoc,cvol,aa
  double complex cof(nlmx),cfac,phase
  data q0 /0d0,0d0,0d0/
  call tcn('pvdf4')
  call stdfac(20,df)
  pi = 4d0*datan(1d0)
  y0 = 1d0/dsqrt(4d0*pi)
  ngabc=lat_nabc
  vol=lat_vol
  ! --- FT of gaussian density, all sites, for list of G vectors ---
  iv0 = 0
  do  10  ib = 1, nbas
     is=ssite(ib)%spec
     tau=rv_a_opos(:,ib) !ssite(ib)%pos
     lmxl=sspec(is)%lmxl
     rg=sspec(is)%rg
     if (lmxl == -1) goto 10
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
     nlm = (lmxl+1)**2
     if (nlm > nlmx) call rxi('pvdf4: increase nlmx to',nlm)
     ilm = 0
     cfac = dcmplx(0d0,1d0)
     do  20  l = 0, lmxl
        cfac = cfac*dcmplx(0d0,-1d0)
        do  21  m = -l, l
           ilm = ilm+1
           cof(ilm) = cfac*qmom(ilm+iv0)*4*pi/df(2*l+1)
21      enddo
20   enddo
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
32      enddo
        !     ... Add foca hankel part
        aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
        cv(ig) = cv(ig) + cfoc*aa*phase
30   enddo
     iv0 = iv0+nlm
10 enddo
  ! --- Potential is 8pi/G**2 * density; overwrite cv with potential ---
  cv(1) = (0d0,0d0)
  do  40  ig = 2, ng
     cv(ig) = (8*pi)*cv(ig)/g2(ig)
40 enddo
  call tcx('pvdf4')
end subroutine pvdf4

end module m_dfrce
