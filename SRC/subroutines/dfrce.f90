!>Correction to force theorem, Harris functional
module m_dfrce 
  use m_lmfinit,only: nsp,nbas,nspec,ispec,spec_a,rmt_i=>rmt,&
       nr_i=>nr,lmxa_i=>lmxa,lmxl_i=>lmxl,spec_z=>z,rg_i=>rg
  use m_ll,only:ll
  use m_fatom,only: sspec
  public dfrce
contains
  subroutine dfrce(job,orhoat,orhoat_out,qmom,smrho,smrout, dfh)
    use m_lmfinit,only: nvl=>pot_nlml,lat_alat
    use m_supot,only: rv_a_ogv,iv_a_okv
    use m_struc_def
    use m_lattic,only: lat_qlat, lat_vol,lat_plat
    use m_supot,only: lat_ng,n1,n2,n3
    use m_lgunit,only:stdo
    use m_rhomom,only: rhomom
    !i Inputs
    !i   n1..3 :dimensions smrho
    !i   nvl   :sum of local (lmxl+1)**2, lmxl = density l-cutoff
    !i   orhoat:vector of offsets containing site density
    !i   orhoat_out:pointer to local densities
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
    !!  orhoat: input atomic density that generaed Hamiltonian
    !!  orhoat_out: new atomic density that the Hamiltonian generated
    !!  smrho:  input  density that generaed Hamiltonian
    !!  smrout: output density that the Hamiltonian generated
    !u Updates
    !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
    implicit none
    type(s_rv1) :: orhoat_out(3,*)
    type(s_rv1) :: orhoat(3,*)
    real(8):: dfh(3,nbas) , qmom(*)
    double complex smrho(n1,n2,n3,*),smrout(n1,n2,n3,*)
    integer :: job,ng,iprint,ib,is,lmxl,iv0,nlm, ip,m,i,ltop,nlmtop,igets,igetss,nn
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
    integer ,allocatable :: iv_iv(:,:)
    complex(8) ,allocatable :: wn1_zv(:)
    complex(8) ,allocatable :: wn2_zv(:)
    complex(8) ,allocatable :: wn3_zv(:)
    double precision :: vol,plat(3,3),qlat(3,3),alat,vsum,pi,tpiba,elind=0d0, fes1(3),fes2(3),fxc(3),c,avgdf(3)
    integer::ibini,ibend,iiv0(nbas)
    character(40) :: strn
    real(8),allocatable:: cs_(:),sn_(:)
    call tcn('dfrce')
    dfh=0d0
    ng=lat_ng
    vol=lat_vol
    alat=lat_alat
    plat=lat_plat
    qlat=lat_qlat
    c = 1000
    nn   = n1*n2*n3
    ! ... Arrays needed for pvdf1
    allocate(ceps_zv(ng))
    allocate(cnomi_zv(ng))
    allocate(cvin_zv(ng))
    allocate(cdvx_zv(ng*nsp))
    ! ... Set up for vectorized Y_lm and gaussians
    ltop = 0
    do   is = 1, nspec
       lmxl = lmxl_i(is)
       ltop = max0(ltop,lmxl)
    enddo
    nlmtop = (ltop+1)**2
    allocate(yl_rv(ng*nlmtop))
    allocate(g2_rv(ng))
    allocate(g_rv(ng*3))
    call suylg ( ltop , alat , ng , rv_a_ogv , g_rv , g2_rv , yl_rv)
    if (allocated(g_rv)) deallocate(g_rv)
    allocate(iv_iv(ng,3))
    iv_iv = nint(matmul(rv_a_ogv,plat))
    ! --- Make ves(rhoin,q) ---
    allocate(smro_zv(nn))
    allocate(cs_(ng), sn_(ng))
    call dcopy(2*nn, smrho,1,smro_zv,1)
    if(nsp==2) call daxpy(2*nn, 1d0, smrho (1,1,1,2), 1, smro_zv, 1)
    call fftz3 ( smro_zv , n1 , n2 , n3 , n1 , n2 , n3 , 1 , 0 , - 1 )
    call gvgetf ( ng , 1 , iv_a_okv , n1 , n2 , n3 , smro_zv , cvin_zv )
    call pvdf4 (  qmom , ng , g2_rv , yl_rv , cs_ , sn_ , iv_iv , qlat , cvin_zv )
    if (allocated(smro_zv)) deallocate(smro_zv)
    deallocate(cs_,sn_)
    ! --- Make dVxc(in)/dn ---
    allocate(dvxc_zv(nn*nsp))
    allocate(smro_zv(nn*nsp))
    allocate(vxcp_zv(nn*nsp))
    allocate(vxcm_zv(nn*nsp))
    allocate(wn1_zv(nn*nsp))
    allocate(wn2_zv(nn*nsp))
    allocate(wn3_zv(nn*nsp))
    call dpcopy ( smrho , smro_zv , 1 , 2 * nn * nsp , 1d0 )
    call pvdf2 ( nbas , nsp, n1 , n2 , n3 , smro_zv , vxcp_zv , vxcm_zv , wn1_zv, wn2_zv , wn3_zv , dvxc_zv )
    deallocate(wn3_zv)
    deallocate(wn2_zv)
    deallocate(wn1_zv)
    deallocate(vxcm_zv)
    deallocate(vxcp_zv)
    ! --- cdvx = FFT ((n0_out-n0_in) dVxc/dn) ---
    !     Use total n0_out-n0_in but keep vxc spin polarized
    call dpzero ( smro_zv , 2 * nn )
    do  i = 1, nsp
       call daxpy (2*nn,  1d0, smrout( 1 , 1 , 1 , i ), 1, smro_zv,1)
       call daxpy (2*nn, -1d0, smrho ( 1 , 1 , 1 , i ), 1, smro_zv,1)
    enddo
    call pvdf3 ( n1 , n2 , n3 , nsp , smro_zv , dvxc_zv )
    call fftz3 ( dvxc_zv , n1 , n2 , n3 , n1 , n2 , n3 , nsp , 0 , - 1 )
    call gvgetf ( ng , nsp , iv_a_okv , n1 , n2 , n3 , dvxc_zv , cdvx_zv )
    ! --- Cnomi = (n0_out(q) - n0_in(q)) ---
    call fftz3 ( smro_zv , n1 , n2 , n3 , n1 , n2 , n3 , 1 , 0 , - 1 )
    call gvgetf ( ng , 1 , iv_a_okv , n1 , n2 , n3 , smro_zv , cnomi_zv  )
    if (allocated(smro_zv)) deallocate(smro_zv)
    if (allocated(dvxc_zv)) deallocate(dvxc_zv)
    ! ... Debugging slot smrho(out) for out-in
    !      print *, '*** debugging ... subs smrout for out-in'
    !      call dpcopy(smrout,w(osmro),1,2*nn,1d0)
    !      call fftz3(w(osmro),n1,n2,n3,n1,n2,n3,1,0,-1)
    !      call gvgetf(ng,1,w(okv),n1,n2,n3,w(osmro),w(ocnomi))
    !      call zprm3('rho-out(q)',w(osmro),n1,n2,n3)

    ! --- Multipole moments of the output density ---
    allocate(qmout_rv(nvl))
    call pshpr(0)
    call rhomom (orhoat_out , qmout_rv , vsum )
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
       is = ispec(ib)
       lmxl = lmxl_i(is)
       nlm = (lmxl+1)**2
       iiv0(ib) = iv0
       iv0 = iv0+nlm
    enddo
    ibini=1
    ibend=nbas
    do ib = ibini, ibend
       is   = ispec(ib)
       lmxl = lmxl_i(is)
       if (lmxl == -1) goto 20
       nlm = (lmxl+1)**2
       call pvdf1 ( job , nsp , ib , iiv0(ib), qmom &
            , qmout_rv , ng , rv_a_ogv , g2_rv , yl_rv , iv_iv , qlat , 0 &
            , cnomi_zv , ceps_zv , cdvx_zv , cvin_zv , orhoat ( 1 , &
            ib ) , fes1 , fes2 , fxc )
       do  i = 1, 3
          dfh(i,ib) = -(fes1(i) + fes2(i) + fxc(i))
       enddo
       if (iprint() >= 30) &
            write(stdo,200) ib,(c*(fes1(m)+fes2(m)),m=1,3), &
            (c*fxc(m),m=1,3),(c*dfh(m,ib),m=1,3)
200    format(i4,3f8.2,1x,3f8.2,1x,3f8.2:1x,3f8.2)
20     continue
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
  subroutine pvdf1(job,nsp,ib,iv0,qmom, qmout,ng,gv,g2,yl,iv,qlat,kmax,cnomin,ceps,cdvxc,cvin , orhoat, fes1,fes2,fxc)
    use m_struc_def 
    use m_lmfinit,only:lat_alat,pnuall,pnzall
    use m_lattic,only: lat_vol,rv_a_opos
    use m_supot,only: n1,n2,n3
    use m_hansr,only:corprm
    ! need to modify texts.
    !- Estimate shift in local density for one site
    ! ----------------------------------------------------------------------
    !i Inputs
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
    integer:: ng , nsp , iv0 , kmax , ib , job , iv(ng,3),i_copy_size
    type(s_rv1) :: orhoat(3)
    real(8):: qmom(*) , qmout(*) , gv(ng,3) , tau(3) , fes1(3) , &
         fes2(3) , fxc(3) , g2(ng) , yl(ng,1) , cs(ng) , sn(ng) , qlat(3,3)
    double complex cdn0(ng,nsp),cdn(ng),cdv(ng),ceps(ng), &
         cnomin(ng),cdvxc(ng,nsp),cvin(ng)
    integer :: ig,ilm,l,lmxl,m,nlm,nlmx,k,is,jv0,jb,js,n0, nrmx
    parameter (nlmx=64, nrmx=1501, n0=10)
    integer :: lmxa,nr,nxi,ie,ixi,job0,kcor,lcor,lfoc,i, nlml
    double precision :: alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg, &
         vol,z,v(3),df(0:20),feso(3),qcor(2),gpot0(nlmx,3),fesdn(3), &
         fesgg(3),pnu(n0),pnz(n0),a,rmt,qloc,exi(n0),hfc(n0,2), &
         qfat,gam,qall,qc,qval,qg,e,aa,q0(3),ssum
    double precision :: rwgt(nrmx),cc,gamf,cfoc,cvol,rsmfa
    complex(8):: tpia,cxx,gc0,xc0,cof(nlmx),phase(ng),img=(0d0,1d0)
    real(8),parameter:: pi = 4d0*datan(1d0),tpi=2d0*pi,y0 = 1d0/dsqrt(4d0*pi)
    data q0 /0d0,0d0,0d0/
    call tcn('pvdf1')
    alat=lat_alat
    vol=lat_vol
    call stdfac(20,df)
    tpia = 2*pi*dcmplx(0d0,-1d0)/alat
    job0 = mod(job,10)
    call dpzero(cdn,2*ng)
    call dpzero(cdn0,2*ng*nsp)
    cdv(1) = 0d0
    call dpzero(fes1,3)
    call dpzero(fxc,3)
    call dpzero(fesdn,3)
    gpot0=0d0
    is=ispec(ib) 
    tau=rv_a_opos(:,ib) 
    phase = exp(-img*tpi*sum(q0*tau)) * exp(-img*tpi*matmul(tau, matmul(qlat, transpose(iv))))
    ! --- Unscreened rigid charge density shift, job 1, in cdn0 ---
    if (job0 == 1) then
       z=spec_z(is)
       pnu=pnuall(1:n0,1,ib)
       pnz=pnzall(1:n0,1,ib)
       lmxa=lmxa_i(is)
       a=  spec_a(is)
       nr= nr_i(is)
       rmt=rmt_i(is)
       lmxl=lmxl_i(is)
       nxi=sspec(is)%nxi
       exi=sspec(is)%exi
       hfc=sspec(is)%chfa
       rsmfa=sspec(is)%rsmfa
       gam  = 0.25d0*rsmfa**2
       call gtpcor(is,kcor,lcor,qcor)
       if (nr > nrmx) call rx('dfrce: nr gt nrmx')
       call radwgt(rmt,a,nr,rwgt)
       nlml = (lmxl+1)**2
       call radsum ( nr , nr , nlml , nsp , rwgt , orhoat( 1 )%v , qloc )
       call radsum ( nr , nr , nlml , nsp , rwgt , orhoat( 2 )%v , ssum )
       qloc = (qloc-ssum)/y0
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
                cdn0(ig,i) = cdn0(ig,i) + aa*phase(ig) !dcmplx(cs(ig),sn(ig))
15           enddo
142       enddo
141    enddo
       !   ... Add gaussian to conserve local charge
       !     ... Add gaussian to conserve local charge.  If density corresponds
       !         to the free-atom density, qfat+qloc = qval+qsc; then qg=0
       cc = qg/vol/nsp
       do   i = 1, nsp
          do   ig = 1, ng
             cdn0(ig,i)=cdn0(ig,i)+cc*phase(ig)*dexp(-gam*g2(ig))
          enddo
       enddo
    endif
    ! --- Coefficients defining local valence + core density ---
    is=ispec(ib)
    tau=rv_a_opos(:,ib) 
    lmxl=lmxl_i(is)
    rg=rg_i(is)
    call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
    nlm = (lmxl+1)**2
    if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
    ilm = 0
    cxx = dcmplx(0d0,1d0)
    do  20  l = 0, lmxl
       cxx = cxx*dcmplx(0d0,-1d0)
       do  22  m = -l,l
          ilm = ilm+1
          cof(ilm) = cxx*qmom(ilm+iv0)*4d0*pi/df(2*l+1)
22     enddo
20  enddo
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
       gc0 = phase(ig)*dexp(-gam*g2(ig))*cvol
       xc0 = dcmplx(0d0,1d0)*dconjg(tpia*cvin(ig))*gc0*vol
       ilm = 0
       do  32  l = 0, lmxl
          xc0 = xc0*dcmplx(0d0,-1d0)
          do  33  m = -l, l
             ilm = ilm+1
             cdn(ig) = cdn(ig) + yl(ig,ilm)*cof(ilm)*gc0
             gpot0(ilm,:) = gpot0(ilm,:) + yl(ig,ilm)*gv(ig,:)*xc0
33        enddo
32     enddo
       !   ... Accumulate unscreened foca density
       aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
       cdn(ig) = cdn(ig) + cfoc*aa*phase(ig)
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
       fes1(:) = fes1(:) + dconjg(cnomin(ig)) * tpia*v(:)*cdv(ig)
       do  i = 1, nsp
          fxc(:)  = fxc(:)  + dconjg(cdvxc(ig,i)) * tpia*v(:)*cdn0(ig,i)
       enddo
30  enddo
    fxc  = fxc*vol
    fes1 = fes1*vol
    ! --- Integral of grad g (output-input local charge) ves~ ---
    fesgg=0d0 !call dpzero(fesgg,3)
    do  k = 1, 3
       do  ilm = 1, nlm
          l = ll(ilm)
          gpot0(ilm,k) = gpot0(ilm,k)*4d0*pi/df(2*l+1)
          fesgg(k) = fesgg(k) + qmout(iv0+ilm)*gpot0(ilm,k)
       enddo
    enddo
    !      print 339, 'n0(out-in) * g dves ',fes1
    !      print 339, 'd(g) qmom(out-in) ves[n0~]',fesgg
    !      print 339, 'n0~(out-in) * dvxc   ',fxc
    !  339 format(a,6p,3f8.2)

    ! --- Integral of dves~ (output-input local charge) for all sites ---
    fes2=0d0
    feso=0d0 
    jv0 = 0
    do  40  jb = 1, nbas
       js=ispec(jb) 
       tau=rv_a_opos(:,jb)
       lmxl=lmxl_i(js)
       rg=rg_i(js)
       nlm = (lmxl+1)**2
       ! ... For this jb, mesh density for all G vectors
       if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
       !     call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
       phase = exp(-img*tpi*sum(q0*tau)) * exp(-img*tpi*matmul(tau, matmul(qlat, transpose(iv))))
       gpot0=0d0
       gam = 0.25d0*rg*rg
       do  50  ig = 2, ng
          aa = dexp(-gam*g2(ig))
          gc0 = dcmplx(0d0,1d0)*aa* &
               dconjg(tpia*cdv(ig))*phase(ig) !dcmplx(cs(ig),sn(ig))
          ilm = 0
          do  55  l = 0, lmxl
             gc0 = gc0*dcmplx(0d0,-1d0)
             do  56  m = -l,l
                ilm = ilm+1
                gpot0(ilm,:) = gpot0(ilm,:)+dble(gc0)*yl(ig,ilm)*gv(ig,:)
56           enddo
55        enddo
50     enddo
       !   ... Multiply factors into gpot0, accumulate force
       ilm = 0
       do  60  l = 0, lmxl
          do  62  m = -l, l
             ilm = ilm+1
             gpot0(ilm,:) = gpot0(ilm,:)*4d0*pi/df(2*l+1)
             feso(:) = feso(:) + qmom(jv0+ilm)*gpot0(ilm,:)
             fes2(:) = fes2(:) + qmout(jv0+ilm)*gpot0(ilm,:)
62        enddo
60     enddo
       jv0 = jv0+nlm
40  enddo
    fes2=fes2-fesgg
    call tcx('pvdf1')
  end subroutine pvdf1
  subroutine pvdf2(nbas,nsp,n1,n2,n3, smrho,vxcp,vxcm,wn1,wn2,wn3,dvxc)
    use m_struc_def
    use m_smvxcm,only: smvxcm
    !- Makes derivative of smoothed xc potential wrt density.
    implicit none
    integer :: nbas,nsp,n1,n2,n3
    complex(8):: vxcp(n1,n2,n3,nsp),vxcm(n1,n2,n3,nsp), &
         dvxc(n1,n2,n3,nsp),smrho(n1,n2,n3,nsp), &
         wn1(n1,n2,n3,nsp),wn2(n1,n2,n3,nsp), &
         wn3(n1,n2,n3,nsp)
    ! ... Local parameters
    integer :: i1,i2,i3,i,nn
    double precision :: fac,dmach,f1,f2,f,alfa,dfdr,rrho,dvdr, &
         rmusm(nsp),rvmusm(nsp),rvepsm(nsp),repsm(nsp),repsmx(nsp),repsmc(nsp), ff(1,1)
    !fcexc0(nsp),fcex0(nsp),fcec0(nsp),fcvxc0(nsp),
    fac = dmach(1)**(1d0/3d0)
    alfa = 2d0/3d0
    nn = n1*n2*n3
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
    call dpzero(wn1, nn*2*nsp)
    call dpzero(wn2, nn*2*nsp)
    call dpzero(wn3, nn*2)
    call smvxcm(0,smrho, vxcp,dvxc,wn1,wn2,wn3,repsm,repsmx,repsmc,rmusm, &
         rvmusm,rvepsm,ff) !,fcexc0,fcex0,fcec0,fcvxc0
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
    call smvxcm(0,smrho,  vxcm,dvxc,wn1,wn2,wn3,repsm,repsmx,repsmc,rmusm, &
         rvmusm,rvepsm,ff) !,fcexc0,fcex0,fcec0,fcvxc0
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
    call smvxcm(0,smrho, vxcm,dvxc,wn1,wn2,wn3,repsm,repsmx,repsmc,rmusm, &
         rvmusm,rvepsm,ff) !,fcexc0,fcex0,fcec0,fcvxc0
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
  subroutine pvdf3(n1,n2,n3,nsp,deln0,dvxc)
    !- Overwrites dvxc with (nout-nin)*dvxc
    implicit none
    integer :: n1,n2,n3,nsp
    double complex deln0(n1,n2,n3),dvxc(n1,n2,n3,nsp)
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
  subroutine pvdf4(qmom,ng,g2,yl,cs,sn,iv,qlat,cv)
    use m_struc_def
    use m_lattic,only: lat_vol,rv_a_opos
    use m_supot,only: n1,n2,n3
    use m_hansr,only:corprm
    !- Makes smoothed ves from smoothed density and qmom, incl nuc. charge
    ! ----------------------------------------------------------------------
    !i Inputs
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
    double complex cv(ng)
    integer :: ig,ib,ilm,is,iv0,l,lmxl,m,nlm,nlmx,lfoc
    parameter (nlmx=64)
    double precision :: tau(3),df(0:20),vol,rg,qcorg,qcorh,qsc, &
         cofg,cofh,ceh,rfoc,z,q0(3),gam,gamf,cfoc,cvol,aa
    complex(8):: cof(nlmx),cfac,phase(ng),img=(0d0,1d0)
    real(8),parameter:: pi = 4d0*datan(1d0),tpi = 2d0*pi,y0 = 1d0/dsqrt(4d0*pi)
    data q0 /0d0,0d0,0d0/
    call tcn('pvdf4')
    call stdfac(20,df)
    vol=lat_vol
    ! --- FT of gaussian density, all sites, for list of G vectors ---
    iv0 = 0
    do  10  ib = 1, nbas
       is=ispec(ib)
       tau=rv_a_opos(:,ib) 
       lmxl=lmxl_i(is)
       rg=rg_i(is)
       if (lmxl == -1) goto 10
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       !     call suphas(q0,tau,ng,iv,n1,n2,n3,qlat,cs,sn)
       phase = exp(-img*tpi*sum(q0*tau)) * exp(-img*tpi*matmul(tau, matmul(qlat, transpose(iv))))
       nlm = (lmxl+1)**2
       if (nlm > nlmx) call rxi('pvdf4: increase nlmx to',nlm)
       ilm = 0
       cfac = dcmplx(0d0,1d0)
       do  20  l = 0, lmxl
          cfac = cfac*dcmplx(0d0,-1d0)
          do  21  m = -l, l
             ilm = ilm+1
             cof(ilm) = cfac*qmom(ilm+iv0)*4d0*pi/df(2*l+1)
21        enddo
20     enddo
       cof(1) = cof(1) + 4*pi*y0*(qcorg-z)
       gam = 0.25d0*rg*rg
       gamf = 0.25d0*rfoc*rfoc
       cfoc = -4d0*pi*y0*cofh/vol
       cvol = 1d0/vol
       do  30  ig = 1, ng
          aa = dexp(-gam*g2(ig))*cvol
          do  32  ilm = 1, nlm
             cv(ig) = cv(ig) + aa*yl(ig,ilm)*cof(ilm)*phase(ig)
32        enddo
          !     ... Add foca hankel part
          aa = dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))
          cv(ig) = cv(ig) + cfoc*aa*phase(ig)
30     enddo
       iv0 = iv0+nlm
10  enddo
    ! --- Potential is 8pi/G**2 * density; overwrite cv with potential ---
    cv(1) = (0d0,0d0)
    do  40  ig = 2, ng
       cv(ig) = (8*pi)*cv(ig)/g2(ig)
40  enddo
    call tcx('pvdf4')
  end subroutine pvdf4
  subroutine suylg(ltop,alat,ng,gv,g,g2,yl)
    use m_ropyln,only: ropyln
    !- Set up vectors g, g2, yl from list of vectors gv
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ltop  :l-cutoff for YL
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !o Outputs
    !o   g     :gv scaled by (2 pi / alat)
    !o   g2    :square of g
    !o   yl    :YL(g)
    !u Updates
    !u   30 May 00 adapted from nfp su_ylg
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ltop,ng
    double precision :: alat,gv(ng,3),g(ng,3),yl(ng,1),g2(ng)
    integer :: i
    double precision :: pi,tpiba
    ! ... Make (2*pi/alat)*gv in g
    pi = 4d0*datan(1d0)
    tpiba = 2d0*pi/alat
    do  i = 1, ng
       g(i,1) = tpiba*gv(i,1)
       g(i,2) = tpiba*gv(i,2)
       g(i,3) = tpiba*gv(i,3)
    enddo
    ! ... Make the yl's and g2
    call ropyln(ng,g(1,1),g(1,2),g(1,3),ltop,ng,yl,g2)
  end subroutine suylg
  subroutine lindsc(job,ng,gv,tpiba,elind,cn)! Make screened density, using a Lindhard response function
    !i Inputs
    !i   job  1s digit specifies what cn contains on output;
    !i          0 cn is untouched
    !i          1 cn is overwritten by eps^-1 cn
    !i          2 cn is overwritten by (eps^-1-1) cn
    !i          3 cn is overwritten by eps^-1.  Input cn is irrelevant
    !i   ng    :number of G-vectors
    !i   gv    :list of dimensionless reciprocal lattice vectors (gvlist.f)
    !i   tpiba :2*pi/lattice constant
    !i   elind :Lindhard screening parameter
    !o Outputs
    !o   cn:     overwritten by
    !o           cn untouched                (job = 0)
    !o           eps^-1 (input cn)           (job = 1)
    !o           (eps^-1 - 1)(input cn)      (job = 2)
    !o           eps                         (job = 3)
    !o           eps^-1                      (job = 4)
    !r Remarks
    !r   The Thomas-Fermi dielectric response is:
    !r     eps_TF(q) = 1 + 4 (k_F a_0) / pi (q a_0)^2
    !r               = 1 + 4 k_F / pi q^2 (Rydberg units)
    !r   The Lindhard function is
    !r     eps(q) = 1 + 4 k_F / pi q^2 * [...] , where
    !r              [...] = 1/2 + (1-x*x)/4x ln((1+x)/(1-x)) and
    !r              x = q / 2 k_F.
    !u Updates
    !u   28 Oct 01 routine revamped.  Original was ridiculously complicated.
    ! ----------------------------------------------------------------------
    !     implicit none
    integer :: job,ng
    double precision :: gv(ng,3),tpiba,elind
    double complex cn(ng)
    ! Local variables
    integer :: i
    double precision :: g2,xx,yy,eps,pi
    logical:: l_dummy_isanrg,isanrg
    pi = 4d0*datan(1d0)
    ! ino isanrg is logical function,       call isanrg(job,0,4,'lindsc:','job', .true.)
    l_dummy_isanrg=isanrg(job,0,4,'lindsc:','job', .true.)
    if (job == 0) return
    cn(1) = 0
    ! ... Early exit if elind is zero => eps = 1
    if (elind == 0) then
       if (job == 1) then
       elseif (job == 2) then
          call dpzero(cn,ng*2)
       elseif (job == 3) then
          do  i = 2, ng
             cn(i) = 1
          enddo
       endif
       return
    endif
    do  22  i = 2, ng
       g2 = tpiba*tpiba*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
       xx = sqrt(g2/elind)/2
       yy = 0.5d0 + (1-xx**2)/(4*xx)*dlog(dabs((1+xx)/(1-xx)))
       eps = 1 + 4*dsqrt(elind)/(pi*g2)*yy
       if (job == 1) then
          cn(i) = cn(i) / eps
       elseif (job == 2) then
          cn(i) = cn(i) * (1/eps-1)
       elseif (job == 3) then
          cn(i) = eps
       elseif (job == 4) then
          cn(i) = 1 / eps
       endif
22  enddo
  end subroutine lindsc
end module m_dfrce
