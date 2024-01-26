!>Correction to force theorem, Harris functional
module m_dfrce 
  use m_lmfinit,only: nsp,nbas,nspec,ispec,spec_a,rmt_i=>rmt,nr_i=>nr,lmxa_i=>lmxa,lmxl_i=>lmxl,spec_z=>z,rg_i=>rg,nlmxlx
  use m_ll,only:ll
  use m_fatom,only: sspec
  public dfrce
contains
  subroutine dfrce(job,orhoat,orhoat_out,qmom,smrho,smrout, dfh)
    use m_lmfinit,only: nvl=>pot_nlml,alat=>lat_alat
    use m_supot,only:  gv=>rv_a_ogv,iv_a_okv
    use m_struc_def
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,plat=>lat_plat
    use m_supot,only:  ng=>lat_ng,n1,n2,n3
    use m_lgunit,only: stdo
    use m_rhomom,only: rhomom
    use m_ropyln,only: ropyln     !- Set up vectors g, g2, yl from list of vectors gv
    !i Inputs
    !i   n1..3 :dimensions smrho
    !i   nvl   :sum of local (lmxl+1)**2, lmxl = density l-cutoff
    !i   orhoat:vector of offsets containing site density
    !i   orhoat_out:pointer to local densities
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
    !i   smrho :smooth density on uniform mesh
    !i   smrho :smooth (input) density that generated the hamiltonian
    !i   smrout:smooth (output) density that the hamiltonian generated
    !o Outputs
    !o   dfh   :correction to the HF force
    !l Local variables
    !l    job  :describes which ansatz for charge shift is used for correction
    !l         :  0  do not calculate correction to force
    !l         :  1   shift in free-atom density
    !l         :  >1 no shift
    !r Remarks
    !! Density
    !!  orhoat:    input atomic density that generaed Hamiltonian
    !!  orhoat_out: new atomic density that the Hamiltonian generated
    !!  smrho:  input  density that generaed Hamiltonian
    !!  smrout: output density that the Hamiltonian generated
    implicit none
    type(s_rv1) :: orhoat_out(3,*)
    type(s_rv1) :: orhoat(3,*)
    real(8):: dfh(3,nbas), qmom(nlmxlx,nbas)
    complex(8):: smrho(n1,n2,n3,nsp),smrout(n1,n2,n3,nsp)
    integer :: job,iprint,is,lmxl, ip,i,lmxlx,nlmtop,nn,ib,m,nlm
    complex(8) ,allocatable :: ceps_zv(:)
    complex(8) ,allocatable :: cnomi_zv(:)
    complex(8) ,allocatable :: smro_zv(:,:,:,:),smav(:,:,:)
    complex(8) ,allocatable :: dvxc_zv(:,:,:,:)
    complex(8) ,allocatable :: vxcp_zv(:)
    complex(8) ,allocatable :: vxcm_zv(:)
    complex(8) ,allocatable :: cdvx_zv(:)
    real(8) ,allocatable :: qmout_rv(:,:)
    complex(8) ,allocatable :: cv(:)
    real(8) ,allocatable :: yl(:,:)
    real(8) ,allocatable :: g2(:)
    real(8) ,allocatable :: g_rv(:,:)
    integer ,allocatable :: iv(:,:)
    complex(8) ,allocatable :: wn1_zv(:)
    complex(8) ,allocatable :: wn2_zv(:)
    complex(8) ,allocatable :: wn3_zv(:)
    real(8):: vsum,tpiba, fes1(3),fes2(3),fxc(3),c=1000d0,avgdf(3),qloc,a,rmt
    real(8),parameter:: pi = 4d0*datan(1d0),tpi = 2d0*pi,y0 = 1d0/dsqrt(4d0*pi)
    integer::ibini,ibend ,i1,i2,i3,nr
    character(40) :: strn
    call tcn('dfrce')
    write(stdo,*)'dfrce job=',job
    dfh=0d0
    nn   = n1*n2*n3
    allocate(cnomi_zv(ng),cv(ng),cdvx_zv(ng*nsp))
    lmxlx= maxval(lmxl_i)
    nlmtop = (lmxlx+1)**2
    allocate(yl(ng,nlmtop),g2(ng),g_rv(ng,3),iv(ng,3))
    tpiba = 2d0*pi/alat
    g_rv(:,:) = tpiba*gv(:,:)
    call ropyln(ng,g_rv(1,1),g_rv(1,2),g_rv(1,3),lmxlx,ng,yl,g2) !Make the yl's and g2
    if (allocated(g_rv)) deallocate(g_rv)
    iv = nint(matmul(gv,plat))
    allocate(smav,source=sum(smrho(:,:,:,:),dim=4))
    call fftz3 ( smav , n1 , n2 , n3 , n1 , n2 , n3 , 1 , 0 , - 1 )
    call gvgetf ( ng , 1 , iv_a_okv , n1 , n2 , n3 , smav , cv ) 
    deallocate(smav)
    pvdf4:block
      use m_lattic,only: rv_a_opos
!      use m_hansr,only:corprm
      integer :: ig,ilm,l,lfoc
      real(8):: tau(3),df(0:20),rg,qcorg,qcorh,qsc,cofg,cofh,ceh,rfoc,z,gam,gamf,cfoc,cvol
      complex(8):: cof(nlmxlx),cfac,phase(ng),img=(0d0,1d0)
      call stdfac(20,df)
      ibloop: do ib = 1, nbas !FT of gaussian density, all sites, for list of G vectors ---
         is  = ispec(ib)
         tau = rv_a_opos(:,ib) 
         lmxl= lmxl_i(is)
         rg=rg_i(is)
         if (lmxl == -1) cycle
         call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
         phase = exp(-img*tpi*matmul(tau, matmul(qlat, transpose(iv))))
         do ilm=1,(lmxl+1)**2 !l = 0, lmxl
            l=ll(ilm)
            cof(ilm) = (-img)**l*qmom(ilm,ib)*4d0*pi/df(2*l+1)
         enddo
         nlm = (lmxl+1)**2
         gam = 0.25d0*rg*rg
         gamf = 0.25d0*rfoc*rfoc
         cfoc = -4d0*pi*y0*cofh/vol
         cvol = 1d0/vol
         do ig = 1, ng
            cv(ig) = cv(ig) + cvol*dexp(-gam*g2(ig))*sum(yl(ig,1:nlm)*cof(1:nlm)*phase(ig)) &
                 +            cfoc*dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))*phase(ig)  !Add foca hankel part
         enddo
      enddo ibloop
      cv(1) = 0d0
      cv(2:ng) = (8*pi)*cv(2:ng)/g2(2:ng) ! --- Potential is 8pi/G**2 * density; overwrite cv with potential ---
    endblock pvdf4
    allocate(dvxc_zv(n1,n2,n3,nsp))
    call pvdf2(n1 , n2 , n3 , smrho, dvxc_zv ) !Make dVxc(in)/dn 
    allocate(smro_zv,mold=smrho) 
    smro_zv(:,:,:,1) = smrout(:,:,:,1)-smrho(:,:,:,1)
    if(nsp==2) smro_zv(:,:,:,1)=smro_zv(:,:,:,1)+ smrout(:,:,:,2)-smrho(:,:,:,2)
    forall(i=1:nsp) dvxc_zv(:,:,:,i) = dvxc_zv(:,:,:,i)*smro_zv(:,:,:,1) !- Overwrites dvxc with (nout-nin)*dvxc
    call fftz3 ( dvxc_zv , n1 , n2 , n3 , n1 , n2 , n3 , nsp , 0 , - 1 )
    call gvgetf ( ng , nsp , iv_a_okv , n1 , n2 , n3 , dvxc_zv , cdvx_zv ) !cdvx = FFT ((n0_out-n0_in) dVxc/dn) Use total n0_out-n0_in but keep vxc spin polarized
    call fftz3 ( smro_zv , n1 , n2 , n3 , n1 , n2 , n3 , 1 , 0 , - 1 )
    call gvgetf ( ng , 1 , iv_a_okv , n1 , n2 , n3 , smro_zv , cnomi_zv  ) !Cnomi = (n0_out(q) - n0_in(q)) ---
    deallocate(smro_zv,dvxc_zv)
    allocate(qmout_rv(nlmxlx,nbas))
    call pshpr(0)
    call rhomom (orhoat_out , qmout_rv , vsum ) !Multipole moments of the output density ---
    call poppr
    qmout_rv = qmout_rv- qmom
    if(iprint()>=30)write(stdo,201) merge('shift in free-atom density','no shift of atomic density',job==1)
201 format(/' Harris correction to forces: ',a/'  ib',9x,'delta-n dVes',13x,'delta-n dVxc',15x,'total')
    do ib = 1,nbas
       is   = ispec(ib)
       lmxl = lmxl_i(is)
       if (lmxl == -1) cycle
       nlm = (lmxl+1)**2
       nr= nr_i(is)
       block
         real(8):: rwgt(nr)
         call radwgt(rmt_i(is),spec_a(is),nr,rwgt)
         call radsum(nr,nr,nlm,nsp,rwgt,orhoat(1,ib)%v-orhoat(2,ib)%v, qloc)
         call pvdf1(job,ib,qmom,qmout_rv,ng,gv,g2,yl,iv,qlat,0,cnomi_zv,cdvx_zv,cv,qloc, fes1,fes2,fxc)
         dfh(:,ib) = -(fes1(:) + fes2(:) + fxc(:))
         if(iprint()>=30)&
              write(stdo,"(i4,3f8.2,1x,3f8.2,1x,3f8.2:1x,3f8.2)")ib,(c*(fes1(m)+fes2(m)),m=1,3),(c*fxc(m),m=1,3),(c*dfh(m,ib),m=1,3)
       endblock
    enddo
    avgdf = sum(dfh,dim=2)/nbas
    forall(ib = 1:nbas) dfh(:,ib) = dfh(:,ib) - avgdf !Shift all forces to make avg correction zero
    if(iprint() >= 30) write(stdo,"(' shift forces to make zero average correction:',8x,3f8.2)") (c*avgdf(m),m=1,3)
    deallocate(qmout_rv,iv,g2,yl,cdvx_zv,cv,cnomi_zv)
    call tcx('dfrce')
  end subroutine dfrce
  subroutine pvdf1(job,ib,qmom, qmout,ng,gv,g2,yl,iv,qlat,kmax,cnomin,cdvxc,cvin,qloc, fes1,fes2,fxc)
    use m_struc_def 
    use m_lmfinit,only:lat_alat
    use m_density,only: pnuall,pnzall
    use m_lattic,only: lat_vol,rv_a_opos
    use m_supot,only: n1,n2,n3
!    use m_hansr,only:corprm
    implicit none
    intent(in)::   job,ib,qmom, qmout,ng,gv,g2,yl,iv,qlat,kmax,cnomin,cdvxc,cvin,qloc
    intent(out)::                                                                               fes1,fes2,fxc
    !- Estimate shift in local density for one site
    !i   ng,gv,kmax,qloc
    !i   job: 1  shift in free-atom density
    !i        /=1  no shift 
    !i   ib      which site is being shifted
    !i   qmom,qmout moments of input and output densities
    !i   cnomin  difference betw. smoothed output and input density n0
    !i   cvin    electrostatic potential of input density Ves[n0~_in]
    !i   ceps    response function
    !i   cdvxc   dVxc/dn (nout-nin)
    !l internal
    !o   cdv:    shift in the electrostatic potential
    !o   fes1,fes2,fxc
    !l   cdn0:  Job 1:  shift in valence part of the free atom density
    !l           Job/=1: =0
    !l   cdn:   Job 1:  dn^(u) where dn is the unscreened shift in the free-atom density.
    !l           Job/=1: =no shift
    !o   NB:     In all cases, the local part of density is approximated a gaussian of the equivalent multipole moment.
    integer:: ng , kmax , ib , job , iv(ng,3),i_copy_size
    type(s_rv1) :: orhoat(3)
    real(8):: qmom(nlmxlx,nbas) , qmout(nlmxlx,nbas) , gv(ng,3) , tau(3) , fes1(3) , &
         fes2(3) , fxc(3) , g2(ng) , yl(ng,1)  , qlat(3,3)
    complex(8):: cdn0(ng,nsp),cdn(ng),cdv(ng), cnomin(ng),cdvxc(ng,nsp),cvin(ng)
    integer :: ig,ilm,l,lmxl,m,nlm,nlmx,k,is,jb,js,n0
    parameter (nlmx=64, n0=10)
    integer :: lmxa,nxi,ie,ixi,kcor,lcor,lfoc,i, nlml
    real(8):: alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg, &
         vol,z,v(3),df(0:20),qcor(2),gpot0(nlmx,3),fesdn(3), & !,feso(3)
         fesgg(3),pnu(n0),pnz(n0),rmt,qloc,exi(n0),hfc(n0,2), &
         qfat,gam,qall,qc,qval,qg,e,q0(3),ssum, cc,gamf,cfoc,cvol
    complex(8):: tpia,cxx,gc0,xc0,cof(nlmx),phase(ng),img=(0d0,1d0)
    real(8),parameter:: pi = 4d0*datan(1d0),tpi=2d0*pi,y0 = 1d0/dsqrt(4d0*pi)
    data q0 /0d0,0d0,0d0/
    call tcn('pvdf1')
    alat=lat_alat
    vol=lat_vol
    call stdfac(20,df)
    tpia = 2*pi*dcmplx(0d0,-1d0)/alat
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
    if (job == 1) then
       z=spec_z(is)
       pnu=pnuall(1:n0,1,ib)
       pnz=pnzall(1:n0,1,ib)
       lmxa=lmxa_i(is)
       nxi=sspec(is)%nxi
       exi=sspec(is)%exi
       hfc=sspec(is)%chfa
       gam = 0.25d0*sspec(is)%rsmfa**2
       call gtpcor(is,kcor,lcor,qcor)
       qfat = 0d0
       do  i  = 1, nsp
          do  ie = 1, nxi
             qall = -4d0*pi*y0*dexp(gam*exi(ie))/exi(ie)
             qfat = qfat + hfc(ie,i)*qall
          enddo
       enddo
       call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor,qc,qval,qsc)
       qg = qval+qsc-qfat-qloc/y0
       !   ... Shift in free atom density
       do  141    i = 1, nsp
          do  142  ixi = 1, nxi
             e = exi(ixi)
             cc = -4d0*pi*hfc(ixi,i)*y0/vol
             do  15  ig = 1, ng
                cdn0(ig,i) = cdn0(ig,i) + cc*dexp(gam*(e-g2(ig)))/(e-g2(ig))*phase(ig)
15           enddo
142       enddo
141    enddo
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
       do  m = -l,l
          ilm = ilm+1
          cof(ilm) = cxx*qmom(ilm,ib)*4d0*pi/df(2*l+1)
       enddo
20  enddo
    ! --- Shift in n0, ves~ for list of G vectors ---
    gam = 0.25d0*rg*rg
    gamf = 0.25d0*rfoc*rfoc
    cfoc = -4d0*pi*y0*cofh/vol
    cvol = 1d0/vol
    igloop:do  30  ig = 2, ng
       v = gv(ig,:)
       !   ... Accumulate unscreened smoothed core+nuclear density
       gc0 = phase(ig)*dexp(-gam*g2(ig))*cvol
       xc0 = dcmplx(0d0,1d0)*dconjg(tpia*cvin(ig))*gc0*vol
       ilm = 0
       do  l = 0, lmxl
          xc0 = xc0*dcmplx(0d0,-1d0)
          do   m = -l, l
             ilm = ilm+1
             cdn(ig) = cdn(ig) + yl(ig,ilm)*cof(ilm)*gc0
             gpot0(ilm,:) = gpot0(ilm,:) + yl(ig,ilm)*gv(ig,:)*xc0
          enddo
       enddo
       cdn(ig) = cdn(ig) + cfoc*dexp(gamf*(ceh-g2(ig)))/(ceh-g2(ig))*phase(ig)
       if(job==1) then !Make the shift in input density n0~
          cdn(ig)=cdn(ig)+sum(cdn0(ig,1:nsp)) !cdn0 = (valence part of) cdn^u ; cdn = cdn^u
          fxc(:) = fxc(:) +  tpia*v(:)*sum(dconjg(cdvxc(ig,:))*cdn0(ig,:))
       endif   
       cdv(ig)  = cdn(ig) * (8*pi/g2(ig))
       fes1(:) = fes1(:) + dconjg(cnomin(ig)) * tpia*v(:)*cdv(ig)
30  enddo igloop
    fxc  = fxc*vol
    fes1 = fes1*vol
    ! --- Integral of grad g (output-input local charge) ves~ ---
    fesgg=0d0 !call dpzero(fesgg,3)
    do  k = 1, 3
       do  ilm = 1, nlm
          l = ll(ilm)
          gpot0(ilm,k) = gpot0(ilm,k)*4d0*pi/df(2*l+1)
          fesgg(k) = fesgg(k) + qmout(ilm,ib)*gpot0(ilm,k)
       enddo
    enddo
    ! --- Integral of dves~ (output-input local charge) for all sites ---
    fes2=0d0
    do  40  jb = 1, nbas
       js=ispec(jb) 
       tau=rv_a_opos(:,jb)
       lmxl=lmxl_i(js)
       rg=rg_i(js)
       nlm = (lmxl+1)**2
       ! ... For this jb, mesh density for all G vectors
       if (nlm > nlmx) call rxi('pvdf1: increase nlmx to',nlm)
       phase = exp(-img*tpi*sum(q0*tau)) * exp(-img*tpi*matmul(tau, matmul(qlat, transpose(iv))))
       gpot0=0d0
       gam = 0.25d0*rg*rg
       do  50  ig = 2, ng
          gc0 = dcmplx(0d0,1d0)*dexp(-gam*g2(ig))* dconjg(tpia*cdv(ig))*phase(ig)
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
             fes2(:) = fes2(:) + qmout(ilm,jb)*gpot0(ilm,:)
62        enddo
60     enddo
40  enddo
    fes2=fes2-fesgg
    call tcx('pvdf1')
  end subroutine pvdf1
  subroutine pvdf2(n1,n2,n3, smrho,dvxc)!- Makes derivative of smoothed xc potential wrt density.
    use m_struc_def
    use m_smvxcm,only: smvxcm
    implicit none
    integer :: n1,n2,n3,i,i1,i2,i3
    complex(8):: vxcp(n1,n2,n3,nsp),vxcm(n1,n2,n3,nsp),vxc0(n1,n2,n3,nsp), &
         dvxc(n1,n2,n3,nsp),smrho(n1,n2,n3,nsp),dsmrho(n1,n2,n3,nsp), wn1(n1,n2,n3,nsp),wn2(n1,n2,n3,nsp),wn3(n1,n2,n3,nsp)
    real(8):: fac,dmach,f1,f2,f,alfa,dfdr,rrho,dvdr, &
         rmusm(nsp),rvmusm(nsp),rvepsm(nsp),repsm(nsp),repsmx(nsp),repsmc(nsp), ff(1,1)
    fac = dmach(1)**(1d0/3d0)
    alfa = 2d0/3d0
    call pshpr(0)
    dsmrho(:,:,:,1) = sum(smrho(:,:,:,:),dim=4)/nsp*fac ! small separation for numerial deivative
    if(nsp==2) dsmrho(:,:,:,2) = dsmrho(:,:,:,1)
    vxcp=0d0; vxcm=0d0; vxc0=0d0; wn1=0d0; wn2=0d0; wn3=0d0
    call smvxcm(0,smrho+dsmrho,  vxcp,dvxc,wn1,wn2,wn3,repsm,repsmx,repsmc,rmusm, rvmusm,rvepsm,ff) ! ... vxcp = vxc (smrho+drho)
    call smvxcm(0,smrho-dsmrho,  vxcm,dvxc,wn1,wn2,wn3,repsm,repsmx,repsmc,rmusm, rvmusm,rvepsm,ff)! ... vxcm = vxc (smrho-drho)
    call smvxcm(0,smrho,         vxc0,dvxc,wn1,wn2,wn3,repsm,repsmx,repsmc,rmusm, rvmusm,rvepsm,ff) 
    wn1=0d0
    do i = 1, nsp
       do concurrent(i1=1:n1,i2=1:n2,i3=1:n3)
          rrho = (smrho(i1,i2,i3,1)+smrho(i1,i2,i3,nsp))/(3-nsp)
          if(rrho>0) then
             dvxc(i1,i2,i3,i)=((vxcp(i1,i2,i3,i)*(rrho*(1+fac))**alfa - vxcm(i1,i2,i3,i)*(rrho*(1-fac))**alfa)/(2d0*fac*rrho) &
                  - alfa*vxc0(i1,i2,i3,i)*rrho**alfa/rrho)/rrho**alfa
          else
             dvxc(i1,i2,i3,i)=wn1(i1,i2,i3,1)
          endif
       enddo
    enddo
    call poppr
  end subroutine pvdf2
end module m_dfrce
