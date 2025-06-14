!> Make the potential at atomic sites and augmentation matrices.
module m_locpot 
  use m_ftox
  use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxb_i=>lmxb,lmxl_i=>lmxl,spec_a,kmxt_i=>kmxt,rg_i=>rg,nlmxlx
  use m_ll,only:ll
  use m_lmfinit,only: lmxax,nsp,nbas,phispinsym
  use m_MPItk,only: master_mpi
  use m_lgunit,only:stdo
  use m_vxcatom,only: vxcnsp
  public locpot
  
  real(8),allocatable,public :: rotp(:,:,:,:,:) !rotation matrix
  real(8),public:: sumt0,sumec,sumtc
  
  private
  real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi  !integer:: idipole
contains
    subroutine locpot(job,novxc,orhoat,qmom,vval,gpot0, osig,otau,oppi,ohsozz,ohsopm,phzdphz,hab,vab,sab, & 
       vvesat,rhoexc,rhoex,rhoec,rhovxc, valvef,xcore, sqloc,sqlocc,saloc, qval,qsc )
    use m_density,only: v0pot,v1pot,pnzall,pnuall !output
    use m_lmfinit,only: lxcf,lhh,nkapii,nkaphh,lmaxu,lldau,n0,nppn,nrmx,nkap0,nlmx,nbas,nsp,lso,ispec,frzwfa,lmxax
    use m_lmfinit,only: slabl,idu,coreh,ham_frzwf,rsma,alat,v0fix,jnlml,vol,qbg=>zbak,lpzex
    use m_uspecb,only:uspecb
    use m_ldau,only: vorb !input. 'U-V_LDA(couter term)' of LDA+U
    use m_fatom,only:sspec
    use m_struc_def
    implicit none
    intent(in)::    job,novxc,orhoat,qmom,vval,gpot0
    !i Inputs
    !i   novxc
    !i   orhoat:vector of offsets containing site density
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
    !i   vval  :electrostatic potential at MT boundary; needed to computed matrix elements of local orbitals.
    !i   gpot0 :integrals of local gaussians * phi0~
    !i         :phi0~ is the estatic potential of the interstitial
    !i   job   := 1 make core and augmentation matrices, :=0 not
    !i rhobg  :compensating background density
    !i lso    :if nonzero, calculate LzSz matrix elements
    !i lmaxu  : max l for U
    !i vorb   : orbital dependent potential matrices
    !i lldau  :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat beginning at dmats(*,lldau(ib))
    !o Outputs
    !o   ... the following are summed over all spheres
    !o   vvesat:integral of valence density times electrostatic potential
    !o   rhoexc:integral of density times xc energy density
    !o   rhoex :integral of density times exch. energy density
    !o   rhoec :integral of density times corr. energy density
    !o   rhovxc:integral of density times xc potential
    !o   valvef:integral (rho1*(v1-2Z/r) - rho2*vsm) ??
    !o   xcore :integral rhoc*(v1-2Z/r)
    !o   sqloc :total valence charge rho1-rho2 in sphere
    !o   saloc :total valence magnetic moment in sphere +core moment for sites with core hole
    !o   sqlocc:total core charge qcor1-qcor2 in sphere
    !o   qval  :nominal total valence charge qv-z
    !o   qsc   :semicore charge from local orbitals
    !o   osig  :augmentation overlap integrals; see augmat
    !o   otau  :augmentation kinetic energy integrals; see augmat
    !o   oppi  :augmentation kinetic + potential integrals; see augmat
    !o   phzdphz  :phz dphz
    !o   hab   :integrals of the ham. with true w.f.  See Remarks in augmat
    !o   vab   :integrals of the pot. with true w.f.  See Remarks in augmat
    !o   sab   :integrals of    unity with true w.f.  See Remarks in augmat
    !l Local variables
    !l   lfltwf:T  update potential used to define basis
    !i   idu  :idu(l+1,ibas)>01 => this l has a nonlocal U matrix
    ! ----------------------------------------------------------------------
    integer:: ib,job,ibx,ir,isp,l,lm, kcor,lcor, i,nglob,ipr,iprint,j1,is,lmxl,lmxa,nr,lmxb,kmax,lfoc,nrml,nlml,&
         ifi, lsox,lmxh, lh(nkap0),nkapi,nkaph,k ,ifivesint
    type(s_rv1) :: orhoat(3,nbas)
    type(s_cv5) :: oppi(3,nbas)
    type(s_sblock):: ohsozz(3,nbas),ohsopm(3,nbas)
    type(s_rv4) :: otau(3,nbas), osig(3,nbas)
    real(8):: qmom(nlmxlx,nbas) , vval(nlmxlx,nbas),rhoexc(nsp),rhoex(nsp),rhoec(nsp),rhovxc(nsp), &
         qval,sqloc,sqlocc,saloc, valvef,vvesat,xcore,gpot0(nlmxlx,nbas),phzdphz(nppn,n0,nsp,nbas),rhobg,&
         eh(n0,nkap0),rsmh(n0,nkap0), ehl(n0),rsml(n0),z,a,rmt,qc,ceh,rfoc, &
         qcorg,qcorh,qsc,cofg,cofh,qsca,rg,qv, qloc,qlocc,xcor, aloc,alocc,rs3,vmtz,qcor(2),qc0,qsc0, ov0mean,pmean
    real(8),target::hab(3,3,n0,nsp,nbas),vab(3,3,n0,nsp,nbas),sab(3,3,n0,nsp,nbas)
    character spid*8
    logical :: lfltwf,cmdopt0,readov0,v0write,novxc
    real(8),pointer:: pnu(:,:),pnz(:,:)
    real(8),allocatable:: wk(:),efg(:,:),zz(:)
    logical,save:: secondcall=.false.
    character*20:: strib
    character strn*120
    call tcn('locpot')
    rhobg=qbg/vol
    ipr = iprint()
    if (ipr >= 30) write(stdo,"('locpot:')")
    k = nrmx*nlmx*nsp
    allocate(efg(5,nbas),zz(nbas))
    xcore   = 0d0
!    if(master_mpi) open(newunit=ifivesint,file='vesintloc',form='formatted',status='unknown')
    if(allocated(rotp)) deallocate(rotp)
    allocate(rotp(0:lmxax,nsp,2,2,nbas))
    ibblock: block
      real(8):: valvs(nbas),valvt(nbas)
      real(8):: qloc(nbas),aloc(nbas),qlocc(nbas),alocc(nbas) 
      real(8):: rhexc(nsp,nbas),rhex(nsp,nbas),rhec(nsp,nbas),rhvxc(nsp,nbas)
      real(8):: xcor(nbas),qv(nbas),qsca(nbas)
      valvs=0d0;valvt=0d0;  qloc=0d0;aloc=0d0;qlocc=0d0;alocc=0d0
      rhexc=0d0;rhex=0d0;rhec=0d0;rhvxc=0d0 ;xcor=0d0;qv=0d0;qsca=0d0
      sumtc=0d0; sumec=0d0; sumt0=0d0
      ibloop: do  ib = 1, nbas
         is=ispec(ib)
         pnu=>pnuall(:,:,ib)
         pnz=>pnzall(:,:,ib)
         z=   z_i(is)
         qc=  sspec(is)%qc
         rg=  rg_i(is)
         a=   spec_a(is)
         nr=  nr_i(is)
         rmt= rmt_i(is)
         lmxa=lmxa_i(is)
         if(lmxa == -1) cycle ! floating orbital
         lmxl=lmxl_i(is)
         lmxb=lmxb_i(is)
         kmax=kmxt_i(is)
         spid=slabl(is) 
         zz(ib)=z
         nlml = (lmxl+1)**2
         lfltwf = (.not.frzwfa(is)).and.(.not.ham_frzwf).and.job==1 ! modify b.c. of Rad.wave func.
         j1 = jnlml(ib) 
         call corprm(is, qcorg,qcorh,qsca(ib),cofg,cofh,ceh,lfoc,rfoc,z)
         call gtpcor(is, kcor,lcor,qcor) !qcor(1:2) is meaningful only when kcor/=0 
         call atqval(lmxa,pnu,pnz,z,kcor,lcor,qcor, qc0,qv(ib),qsc0)
         if(qsc0 /= qsca(ib) .OR. qc /= qc0-qsc0) then
            if(ipr>0)write(stdo,ftox)' is=',is,'qsc0=',ftof(qsc0),'qsca',ftof(qsca(ib)),'qc',ftof(qc),'qc0',ftof(qc0)
            call rxs('problem in locpot -- possibly low LMXA or orbital mismatch, species ',spid)
         endif
         if(ipr>= 20) write(stdo,"(' site',i3,'  z=',f7.3,'  rmt=',f8.5,'  nr=',i3,'   a=',f5.3, &
              '  nlml=',i2,'  rg=',f5.3,'  Vfloat=',l1)") ib,z,rmt,nr,a,nlml,rg,lfltwf
         if(ipr>=30.and. kcor/=0 .and. sum(abs(qcor))/=0 ) write(stdo,ftox)&
              ' core hole: kcor=',kcor,'lcor=',lcor,'qcor amom=',ftof(qcor)
         locpt2augmat: block
           real(8)::rhol1(nr,nlml,nsp),rhol2(nr,nlml,nsp),v1(nr,nlml,nsp),v2(nr,nlml,nsp),v1es(nr,nlml,nsp),v2es(nr,nlml,nsp),&
                gpotb(nlml),rofi(nr),rwgt(nr),v1out(nr,nlml,nsp)
           !real(8),pointer ::rho1(:,:,:),rho2(:,:,:),rhoc(:,:)
           !real(8),allocatable ::rho1(:,:,:),rho2(:,:,:),rhoc(:,:) !nr,nlml,nsp),rhoc(nr,nsp)
           !allocate(rho1,source=reshape(orhoat(1,ib)%v,[nr,nlml,nsp]))
           !allocate(rho2,source=reshape(orhoat(2,ib)%v,[nr,nlml,nsp]))
           !allocate(rhoc,source=reshape(orhoat(3,ib)%v,[nr,nsp]))
           call radmsh(rmt,a,nr,rofi)
           call radwgt(rmt,a,nr,rwgt)
           if(cmdopt0('--wrhomt'))call wrhomt('rhoMT.','density',ib,orhoat(1,ib)%v,rofi,nr,nlml,nsp) ! Write true density to file rhoMT.ib
           call locpt2(ib,z,rmt,rg,a,nr,nsp,cofg,cofh, & ! Make potential and energy terms at this site ---
                ceh,rfoc,lfoc,nlml,qmom(:,ib),vval(:,ib),rofi,rwgt, orhoat(1,ib)%v,orhoat(2,ib)%v,orhoat(3,ib)%v,   gpot0(:,ib), &
                rhol1,rhol2,v1,v2,v1es,v2es,&
                valvs(ib),rhexc(:,ib),rhex(:,ib),rhec(:,ib),rhvxc(:,ib),&
                valvt(ib),xcor(ib) ,qloc(ib),qlocc(ib),aloc(ib),alocc(ib),gpotb, rhobg,efg(1,ib),ifivesint,lxcf,    v1out) 
           ! write density 1st(true) component and counter components.
           if(cmdopt0('--density') .AND. master_mpi .AND. secondcall) then
              write(stdo,"(' TotalValenceChange diff in MT; ib,\int(rho2-rho1)=',i5,f13.5)") ib,qloc(ib)
              write(strib,'(i10)') ib
              open(newunit=ibx,file='rho1MT.'//trim(adjustl(strib)))
              do isp=1,nsp
                 do ir=1,nr
                    write(ibx,"(i2,d13.6,256f13.6)")isp,rofi(ir),(rhol1(ir,lm,isp),lm=1,nlml)
                 enddo
              enddo
              close(ibx)
              open(newunit=ibx,file='rho2MT.'//trim(adjustl(strib)))
              do isp=1,nsp
                 do ir=1,nr
                    write(ibx,"(i2,d13.6,256f13.6)")isp,rofi(ir),(rhol2(ir,lm,isp),lm=1,nlml)
                 enddo
              enddo
              close(ibx)
           endif
           if(cmdopt0('--wpotmt'))call wrhomt('vtrue.','potential',ib,v1,rofi,nr,nlml,nsp)! Write true potential to file vtrue.ib
           if(lfltwf) v0pot(ib)%v(1:nr,1:nsp) = y0*v1out(1:nr,1,1:nsp) ! Update the potential used to define radial basis set
           phispinsymB: block ! spin averaged oV0 to generate phi and phidot. takaoAug2019
             !phispinsym= cmdopt0('--phispinsym')
             if(phispinsym) then
                if(master_mpi.AND.nsp==2)write(6,*) 'locpot: --phispinsym mode: use spin-averaged potential for phi and phidot'
                do ir =1,nr
                   v0pot(ib)%v(ir,:)= sum([(v0pot(ib)%v(ir,isp),isp=1,nsp)])/nsp
                enddo   
             endif
           endblock phispinsymB
           v0fixblock: block ! experimental case --v0fix
             character charext*8
             real(8):: ov0(nr)
             if(v0fix) then
                inquire(file='v0pot.'//trim(charext(ib)),exist=readov0)
                write(6,*)'v0fixmode=',readov0,ib,nr
                if(.not.readov0) call rx('no v0pot files')
                open(newunit=ifi,file='v0pot.'//trim(charext(ib)),form='unformatted')
                read(ifi) ov0(1:nr)
                close(ifi)
                forall(ir=1:nr) v0pot(ib)%v(ir,:)= ov0(ir)
             endif
           endblock v0fixblock
           if(master_mpi .AND. nsp==2)then
              do l=0,lmxa
                 write(6,"(' ibas l=',2i3,' pnu(1:nsp) pnz(1:nsp)=',4f10.5)") ib,l,pnu(l+1,1:nsp),pnz(l+1,1:nsp)
              enddo
           endif
           v1pot(ib)%v(1:nr,1:nsp) = y0*v1out(1:nr,1,1:nsp) ! Store the potential used in mkrout to calculate the core
           if(lfoc==0) xcore = xcore + xcor(ib)
           if(kcor/=0.and.(dabs(qcor(2)-alocc(ib))>0.01d0).and.ipr>=10) & !  Check for core moment mismatch ; add to total moment
                write(stdo,ftox) ' (warning) core moment mismatch spec=',is,'input file=',qcor(2),'core density=',alocc
           ! Make augmentation matrices sig, tau, ppi ---
           if (job==1) then !     ... Smooth Hankel tails for local orbitals
              rsmh= 0d0
              eh  = 0d0
              call uspecb(is,rsmh,eh)
              nkapi=nkapii(is)
              nkaph=nkaphh(is)
              block
                if(lpzex(is)==1) rsml=rsmh(:,nkaph)
                if(lpzex(is)==1)  ehl=  eh(:,nkaph)
              endblock
              if (novxc) then !  ... Use effective potentials with modified xc
                 v1=v1es 
                 v2=v2es 
              endif
              !          if(idipole/=0) call adddipole(v1,rofi,nr,nlml,nsp,idipole,ssite(ib)%pos(idipole)*alat)
              !             call adddipole(v2,rofi,nr,nlml,nsp,idipole,ssite(ib)%pos(idipole)*alat)
              lsox = merge(1, lso, .NOT. novxc .AND. cmdopt0('--socmatrix') )
              lmxh = lmxb !MTO l of basis minimum
              if (ipr >= 20) write(stdo,"('     potential shift to crystal energy zero:',f12.6)") y0*(gpot0(1,ib)-gpotb(1))
              augmatblock: block
                use m_gaugm,only:  gaugm
                use m_augmat,only: vlm2us,momusl
                use m_potpus,only: potpus
                integer :: k,nlma,nlmh,i, lxa(0:kmax),kmax1
                real(8):: v0(nr,nsp),rsmaa,&
                     vdif(nr,nsp),sodb(3,3,n0,nsp,2), vum((lmxa+1)**2,nlml,3,3,nsp),qum((lmxa+1)**2,(lmxl+1),3,3,nsp),  &
                     fh(nr*(lmxh+1)*nkap0),   xh(nr*(lmxh+1)*nkap0),   vh((lmxh+1)*nkap0),   dh((lmxh+1)*nkap0), &
                     fp(nr*(lmxa+1)*(kmax+1)),xp(nr*(lmxa+1)*(kmax+1)),vp((lmxa+1)*(kmax+1)),dp((lmxa+1)*(kmax+1))
                complex(8):: vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,3,3,2,0:lmaxu)
                real(8),pointer:: hab_(:,:,:,:),sab_(:,:,:,:),vab_(:,:,:,:)
                rsmaa=rsma(is)
                nlma = (lmxa+1)**2
                lxa=lmxa
                nlmh = (lmxh+1)**2
                v0 = v0pot(ib)%v
                hab_=>hab(:,:,:,:,ib)
                sab_=>sab(:,:,:,:,ib)
                vab_=>vab(:,:,:,:,ib)
                vdif(1:nr,1:nsp) = y0*v1(1:nr,1,1:nsp) - v0(1:nr,1:nsp) !vdif= extra part of spherical potential for deterimning radial function
                call potpus(z,rmt,lmxa,v0,vdif,a,nr,nsp,lsox,pnu,pnz,ehl,rsml,rs3,vmtz, & !hab,vab,sab and phzdphz, and rotp
                     phzdphz(:,:,:,ib),hab_,vab_,sab_,sodb,rotp(0:lmxa,:,:,:,ib)) 
                call momusl(z,rmt,lmxa,pnu,pnz,rsml,ehl,lmxl,nlml,a,nr,nsp,rofi,rwgt,v0,v1, qum,vum)!Moments and potential integrals of ul*ul, ul*sl, sl*s
                call fradhd(nkaph,eh,rsmh,lhh(:,is),lmxh,nr,rofi, fh,xh,vh,dh) !head
                call fradpk(kmax,rsmaa,lmxa,             nr,rofi, fp,xp,vp,dp) !tail
                if(lldau(ib)>0)call vlm2us(lmaxu,rmt,idu(:,is),lmxa, count( idu(:,ispec(1:ib-1))>0 ), & !offset to the Ublock for ib
                     vorb,phzdphz(:,:,:,ib),rotp(0:lmxa,:,:,:,ib), vumm)!LDA+U: vumm of (u,s,gz) from vorb for phi.
                kmax1=kmax+1
                call gaugm(nr,nsp,lsox,rofi,rwgt,lmxa,lmxl,nlml,v2,gpot0(1:nlml,ib)-gpotb(1:nlml),hab_,vab_,sab_,sodb,qum,vum,& !...Pkl*Pkl !tail x tail
                     lmaxu,vumm,lldau(ib),idu(:,is),  lmxa,nlma,nlma,& !lmxa=lcutoff for augmentation
                     kmax1,kmax1,lmxa,lxa, fp,xp,vp,dp,&
                     kmax1,kmax1,lmxa,lxa, fp,xp,vp,dp,&
                     osig(1,ib)%v, otau(1,ib)%v, oppi(1,ib)%cv, ohsozz(1,ib)%sdiag, ohsopm(1,ib)%soffd)
                call gaugm(nr,nsp,lsox,rofi,rwgt,lmxa,lmxl,nlml,v2,gpot0(1:nlml,ib)-gpotb(1:nlml),hab_,vab_,sab_,sodb,qum,vum,& !...Hsm*Pkl! head x tail
                     lmaxu,vumm,lldau(ib),idu(:,is),  lmxh, nlmh,nlma, &!lmxh=lcutoff for basis (lmxh<=lmxa is assumed)
                     nkaph,nkapi,lmxh,lhh(:,is),fh,xh,vh,dh,&
                     kmax1,kmax1,lmxa,lxa,      fp,xp,vp,dp,&
                     osig(2,ib)%v, otau(2,ib)%v, oppi(2,ib)%cv, ohsozz(2,ib)%sdiag, ohsopm(2,ib)%soffd)
                call gaugm(nr,nsp,lsox,rofi,rwgt,lmxa,lmxl,nlml,v2,gpot0(1:nlml,ib)-gpotb(1:nlml),hab_,vab_,sab_,sodb,qum,vum,& !...Hsm*Hsm! head x head
                     lmaxu,vumm,lldau(ib),idu(:,is),  lmxh,nlmh,nlmh,&
                     nkaph,nkapi,lmxh,lhh(:,is),fh,xh,vh,dh,&
                     nkaph,nkapi,lmxh,lhh(:,is),fh,xh,vh,dh,&
                     osig(3,ib)%v, otau(3,ib)%v, oppi(3,ib)%cv, ohsozz(3,ib)%sdiag, ohsopm(3,ib)%soffd)
              endblock augmatblock
           endif
         endblock locpt2augmat
      enddo ibloop
      vvesat = sum(valvs)
      valvef=  sum(valvt)
      sqloc  =  sum(qloc)
      saloc  =  sum(aloc)
      if(kcor/=0) saloc = sum(aloc) + qcor(2)
      sqlocc =  sum(qlocc)
      rhoexc =  sum(rhexc,dim=2)
      rhoex  =  sum(rhex,dim=2)
      rhoec  =  sum(rhec,dim=2)
      rhovxc =  sum(rhvxc,dim=2)
      qval = sum(qv)+sum(qsca)
      qsc  = sum(qsca)
    endblock ibblock
    if(cmdopt0('--density') .AND. master_mpi) secondcall= .TRUE. 
!    if(master_mpi) close(ifivesint)
    if (ipr > 40) call elfigr(nbas,stdo,zz,efg)!  Electric field gradient
    deallocate(efg,zz)!,rhol1,rhol2,v1,v2,v1es,v2es)
    call tcx('locpot')
  end subroutine locpot
  subroutine locpt2(ib,z,rmt,rg,a,nr,nsp,cofg,cofh,ceh,rfoc,lfoc, & !- Makes the potential at one site, and associated energy terms.
       nlml,qmom,vval,rofi,rwgt,rho1,rho2,rhoc,gpot0,&
       rhol1,rhol2,v1,v2,v1es,v2es,&
       vvesat,rhoexc,rhoex,rhoec,rhovxc, valvef,xcore,qloc, & 
       qlocc,aloc,alocc,gpotb,rhobg,efg,ifivesint,lxcfun  , v1out) 
    use m_hansmr,only: hansmr,hansmronly
    use m_hansr,only:  hansr
    !i Inputs
    !i   z     :nuclear charge
    !i   rmt   :augmentation radius
    !i   rg    :smoothing radius for compensating gaussians used to correct the multipole moments of local smooth density
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
    !i         :cofg = Y0 * qcorg
    !i   cofh  :coefficient to Hankel part of pseudocore density (corprm)
    !i   ceh   :energy of hankel function to fit core tail
    !i   rfoc  :smoothing radius for hankel head fitted to core tail
    !i   lfoc  :switch specifying treatment of core density.
    !i         :0 => core w.f. have val,slo = 0 at rmt
    !i         :1 => core included explicitly with valence
    !ixxx         :2 => core included perturbatively
    !i   nlml  :L-cutoff for rho1,rho2
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
    !i   vval  :boundary condition for estat potential at MT boundary.
    !i         :See remarks.
    !i   rofi  :radial mesh points for density and potential
    !i   rwgt  :radial mesh weights to integrate function on radial mesh
    !i         :Integral f(r) dr = sum_i f_i * wt_i
    !i   rho1  :local valence density = sum_ilm rho1(ilm) * r**2  Y_L(ilm),
    !i         :on radial mesh.
    !i   rho2  :local smoothed valence density, defined as rho1
    !i         :Local atomic valence density is rho1-rho2
    !i   rhoc  :core density times 4*pi*r*r
    !i   rhobg: compensating background density
    !o Outputs
    !o   vvesat:integral of valence density times electrostatic potential
    !o   rhoexc:integral of valence density times xc energy density
    !o   rhoex :integral of valence density times exch. energy density
    !o   rhoec :integral of valence density times corr. energy density
    !o   rhovxc:integral of valence density times xc potential
    !o   valvef:integral valence density times full potential:
    !o         :used to compute double-counting for kinetic energy.
    !o         : = vefv1 - vefv2; defined in Local variables below
    !o   xcore :integral rhoc*(v1-2Z/r)
    !o   qloc  :total valence charge in sphere rho1 - rho2
    !o   aloc  :total valence magnetic moment in sphere
    !o   qlocc :total core charge in sphere qcor1 - qcor2
    !o   alocc :total core magnetic moment in sphere
    !o   gpotb :integrals of gaussians(radius rg) times smooth ves.
    !o         :Here smooth ves is the electrostatic potential of the
    !o         :compensated sm density rhol2, w/ b.c. ves(rmax)=0
    !o         :This is a local analog of gpot0 generated for the smooth
    !o         :density.  Here gpotb is generated in a Y_lm expansion
    !o         :of the gaussian orbitals
    !o   rhol1 :full true electron density, rho1 + rhoc
    !o   rhol2 :full smooth density, i.e. uncompensated rho2 plus
    !o         :compensating gaussians + pseudocore charge
    !o   v1    :true potential at this site, excluding nuclear contribution
    !o         :but including background contribution; see Remarks
    !o         :See Remarks concerning boundary condition at rmt
    !o   v2    :Total smooth potential including background potential
    !o         : = ves[rhol2] + ...
    !o         :   ... lfoc=0 : vxc(rho2)
    !o         :       lfoc=1 : vxc(rho2+sm-rhoc)
    !o         :Apart from differences in l-truncation in the XC potential
    !o         :(i.e. true and local sm. densities are different),
    !o         :v1(nr,1,1)-2*z/y0/rmt = v2(nr,1,1)
    !o         :See Remarks concerning boundary condition at rmt
    !o  efg    :l=2 potential at nucleus
    !l Local variables
    !l   rhol1 :full electron density, rho1 + rhoc
    !l   rhol2 :full smooth compensated density rho2 + gval + gcor + gnuc
    !l   rhochs:Part of smooth core density contains Hankel tails
    !l         :Used to make vxc[n2] and exc[n2]
    !l         :rhochs = srfpi*cofh*xi(0)*r*r
    !l   vefv1 :int (valence true density)  * (total true potential):
    !l         : = int n1 * v1~
    !l   vefv2 :int (valence smooth density)* (total smooth potential):
    !l         : = int (n2*v2~) +  sum_L qmom_L gpotb_L
    !l   qcor1 :true core charge inside rmt
    !l   qcor2 :smoothed core charge inside rmt
    !l   rhonsm:nuclear density smoothed into gaussian of width rg
    !l   rhocsm:core density smoothed into gaussian of width rg
    !r Remarks
    !r   On boundary conditions for the estat potential vval at the MT
    !r   boundary.  In the original formulation (see Springer book) the
    !r   the total hamiltonian and energy did not depend on the choice of
    !r   vval, as the contribution from true and smooth potentials exactly
    !r   cancel.  However, this is no longer true when local orbitals are
    !r   present that have no corresponding smooth part.
    !r
    !r   The potential due to uniform background density rhobg is added
    !r                 vbg= -(4pi/3)*(rmt^2-r^2)*rhobg to v1 and v2
    !r   in spherical components .
    implicit none
    integer:: nr,nsp,lfoc,nlml,ib, ipr,iprint,i,isp,ilm,l,lxcfun,nglob,nrml, ifivesint
    real(8):: rofi(nr),rwgt(nr), qmom(nlml),vval(nlmxlx),gpotb(nlml),gpot0(nlml), &
         rho1(nr,nlml,nsp),rhol1(nr,nlml,nsp), rho2(nr,nlml,nsp),rho2s(nr,nlml,nsp),rhol2(nr,nlml,nsp), &
         v1(nr,nlml,nsp),v1es(nr,nlml,nsp), v2(nr,nlml,nsp),v2es(nr,nlml,nsp), &
         wk(nr,nlml,nsp),rhoc(nr,nsp), efg(5), z,rmt,rg,a,cofg,cofh,ceh,rfoc,xcore,qloc,qlocc, &
         aloc,alocc,rhoexc(nsp),rhoex(nsp),rhoec(nsp),rhovxc(nsp), valvef,vvesat, rhobg,&
         rhochs(nr),rhonsm(nr),df(0:20),cof(nlml), rhocsm(nr),tmp(nsp),xill(nr),&
         ag2,cof0,fac,qv1,qv2,qcor1,qcor2, r,rep1(nsp),rep2(nsp),rep1x(nsp),rep2x(nsp),rep1c(nsp),rep2c(nsp), &
         rhves1,rhves2,rmu1(nsp),rmu2(nsp),rvs1,rvs2, rvsm(nsp),rvtr(nsp),samh,sfac,sgpotb,sum1,sum2,sumg,sumh,top, &
         ves1,vales1,vales2,vcpn1,vcpn2,vefv1,vefv2,vesc1,vesc2, &
         vesn1,vesn2,vnucl,vsum,vtr,a1,a2,smrhoc, qs(nsp),ves1int,ves2int, w2(nsp),fl(1,1,1),gnu,gg(nr), &
         v1out(nr,nlml,nsp),dEdQ(nlml)
    real(8),parameter:: pi = 4d0*datan(1d0),srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi,pi4=4d0*pi
    logical:: debug=.false.,topl
    call tcn('locpt2')
    call stdfac(20,df)
    ipr = iprint()
    ag2  = 1d0/rg**2
    nrml = nr*nlml
    gg   =  rofi**2 * dexp(-ag2*rofi**2) ! ... Renormalize gaussian
    if(abs(pi4*(ag2/pi)**1.5d0*sum(rwgt(2:nr)*gg(2:nr))-1d0)>1d-4)&
         write(stdo,ftox)' locpot (warning): large gaussian, integral=',ftod(sumg)
    fac  = 1d0/(sum(rwgt(2:nr)*gg(2:nr)))
    sfac = fac/(pi4*(ag2/pi)**1.5d0)
    do i=1,nr
      call hansmr(rofi(i),ceh,1d0/rfoc,xill(i),0) !Smooth nuc. and core rho, sm Hankel portion, true & smooth core q
    enddo  
    rhochs(:) = srfpi*cofh*xill(:)*rofi(:)**2   ! n^c_sH,a(r) smoothedcore density smHankel
    if(ipr>=20) then
       sumh  = sum(rwgt*rhochs)                           !
       samh = -y0*cofh*pi4*dexp(ceh*rfoc*rfoc*0.25d0)/ceh !total sm core charge (smHamkel)
       if(dabs(samh)>1d-6)write(stdo,ftox)'    sm core charge in MT=',ftof(sumh),'= total-spillout=',ftof(samh),'-',ftof(samh-sumh)
    endif
    rhol1 = rho1
    rhol1(:,1,:)= rhol1(:,1,:)+ y0*rhoc(:,:) ! True electron density = rhol1 -2Z/r = rho1 + y0*(rhoc -2Z/r) !1st component of Eq.(24)
    !     gval : qmom * compensating gaussians gg = \sum_L QaL^v GaL 
    do isp = 1, nsp
       do ilm = 1, nlml
          l = ll(ilm)
          rhol2(:,ilm,isp)=rho2(:,ilm,isp)+ qmom(ilm) *pi4/df(2*l+1)*sfac*(ag2/pi)**1.5d0 *(2d0*ag2)**l *rofi(:)**l*gg(:)/nsp
       enddo
       rhol2(:,1,isp)= rhol2(:,1,isp) + y0/nsp*rhochs ! rhol2 = n^c_sH,a + Gaussians(qmom) + rho2.  Eq.(30)
       ! caution: Eq(30) has typo. wrong: n^Zc_2,a+... ==> right: n^c_sH,a+...
    enddo
    if(nsp==2) then
       aloc = srfpi*sum(rwgt*((rho1(:,1,1)-rho1(:,1,2))-(rho2(:,1,1)-rho2(:,1,2))))
       alocc= sum(rwgt*(rhoc(:,1)-rhoc(:,2)))
    else
       aloc = 0d0
       alocc= 0d0
    endif
    rhototal: block
      real(8):: rhol1t(nr,nlml),rhol2t(nr,nlml),rmax
      rhol1t=sum(rhol1(:,:,1:nsp),dim=3)
      rhol2t=sum(rhol2(:,:,1:nsp),dim=3)
      rhol1t(:,1)=rhol1t(:,1)+srfpi*rhobg*rofi(:)**2 ! ... Add background density to spherical parts of rhol1 and rhol2
      rhol2t(:,1)=rhol2t(:,1)+srfpi*rhobg*rofi(:)**2
      qv1   = srfpi*sum(rwgt*sum(rho1(:,1,1:nsp),dim=2)) ! valence charges 1st
      qv2   = srfpi*sum(rwgt*sum(rho2(:,1,1:nsp),dim=2))
      qloc  = qv1-qv2
!      qcor1 = sum(rwgt*sum(rhoc(:,1:nsp),dim=2))          ! core charge 1st
!      qcor2 = sum(rwgt*(rhochs(:)+srfpi*cofg*fac*gg(:))) ! smcore+couter gaussian charge 2nd.
!      if(abs(qcor1- qcor2)>1d-6) call rx1('locpt2: qlocc/=0 strange core charge =',qcor1-qcor2)
      sum1 = srfpi*sum(rwgt*rhol1t(:,1)) - z !MT charge of \bar{n}^ZcV of 1st component  Eq.(29) See JPSJ
      sum2 = srfpi*sum(rwgt*rhol2t(:,1))     !MT charge of \bar{n}^ZcV of 2nd component  Eq.(30)
      if(dabs(sum1-sum2)>1d-6) call rx1('locpt2: sphere not neutral: charge =',sum1-sum2)
      Getv1out: block     
        call poinsp(z,vval,nlml,v1out,rofi,rhol1t,wk,nr,rvs1,rhves1,  vnucl,vsum)
        ! v1out is with the b.c. vval detemined by \bar{n0}^Zcv. This is needed for pnunew (set energy at the center of gravity).
        if (nsp == 2) v1out(:,:,2)=v1out(:,:,1)
        call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rhol1,lxcfun,rep1,rep1x,rep1c,rmu1,v1out,fl,qs) !v1out is for radial basis and pnu.
        xcore = sum([(sum(rwgt(2:nr)*rhoc(2:nr,isp)*(y0*v1out(2:nr,1,isp)-2d0*z/rofi(2:nr))),isp=1,nsp)]) ! Vin*rhoc
      endblock Getv1out
      v1esv2esgpotb: block
        call poinsp(0d0,[(0d0,i=1,nlml)],nlml,v2es,rofi,rhol2t,wk,nr,rvs2,rhves2, vnucl,vsum)! Ves[rhol2], Ves=0 at MTboundaries (vval=0).
        rmax=rofi(nr)
        do ilm = 1,nlml
           l = ll(ilm)
           gpotb(ilm) = pi4/df(2*l+1)*sfac*(ag2/pi)**1.5d0*(2d0*ag2)**l*sum(rwgt*v2es(:,ilm,1)*rofi(:)**l* gg(:) )
           dEdQ(ilm)  = (gpot0(ilm) - gpotb(ilm))*rmax**l  ! {\cal Q}^V_aL*rmax**l in Eq.(35) 
        enddo !NOTE     dEdQ gives the boundariy condition at MT
        call poinsp(z,  dEdQ,nlml,v1es,rofi,rhol1t,wk,nr,rvs1,rhves1, vnucl,vsum)! v1es-2*z(1/rofi(2:nr)-1/rmt) is the es part of 1st comp. of Eq.(34).
        call poinsp(0d0,dEdQ,nlml,v2es,rofi,rhol2t,wk,nr,rvs2,rhves2, vnucl,vsum)!                         v2es is the es part of 2nd comp. of Eq.(34).
        do ilm = 1, nlml ! gpotb = integrals of compensating gaussians times the estatic 2nd component v2es
           l = ll(ilm) 
           gpotb(ilm) = pi4/df(2*l+1)*sfac*(ag2/pi)**1.5d0*(2d0*ag2)**l*sum(rwgt*v2es(:,ilm,1)*rofi(:)**l* gg(:) )
        enddo
!        if(master_mpi) then
!           ves1int = pi4*(sum(rwgt*y0*v1es(:,1,1)*rofi(:)**2) - z*rofi(nr)**2)
!           ves2int = pi4*sum(rwgt*y0*v2es(:,1,1)*rofi(:)**2)
!           write(ifivesint,ftox)ib,ftof(ves1int-ves2int),ftof(ves1int),ftof(ves2int),' ! ib estatic average vesint1-vesint2 ves1int ves2int'
!        endif
      endblock v1esv2esgpotb
      efg(1:5)= merge(v1es(5,5:9,1)/rofi(5)**2,0d0,nlml >= 9 .AND. z > 0.01) !electric field at nucleus
      vnucl = 2d0*srfpi*sum(rwgt(2:nr)*rhol1t(2:nr,1)*(1d0/rofi(2:nr)-1d0/rmt)) + 2d0*z/rmt + y0*dEdQ(1) != v1es+vcore at ir=0 without 2z/r. Note b.c.
      ! Estatic term of 1st comp. of Eq.34. is given as  v1es +  2d0*(-z/r+z/rmt) =  Ves(rhol1t,zero at rmt)+y0*dEdQ(1) + 2d0*(-z/r+z/rmt)  
      vvesat = rvs1-z*vnucl - rvs2  ! density \times electrostatic potential (z-z self-interaction removed).
    endblock rhototal
    if(nsp==2) v1es(:,:,2)=v1es(:,:,1)
    if(nsp==2) v2es(:,:,2)=v2es(:,:,1)
    forall(isp=1:nsp) v1(:,:,isp) =v1es(:,:,1) 
    forall(isp=1:nsp) v2(:,:,isp) =v2es(:,:,1)
    AddVxc:block
      call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rhol1,lxcfun,rep1,rep1x,rep1c,rmu1,v1,fl,qs) !add xc to v1
      rho2s = rho2
      forall(isp=1:nsp) rho2s(1:nr,1,isp)=rho2s(1:nr,1,isp) + merge(y0/nsp*rhochs(1:nr),0d0,lfoc==1)    !add xc to v2
      call vxcnsp(0,a,rofi,nr,rwgt,nlml,nsp,rho2s,lxcfun,rep2,rep2x,rep2c,rmu2,v2,fl,qs)
    endblock AddVxc
    vefv: block 
      real(8):: rvsm(nlml,nsp),rvtr(nlml,nsp)
      do concurrent (ilm = 1:nlml) ! --- Integrals involving the full nonspherical potential ---
         rvtr(ilm,:) = [(sum(rwgt(2:nr)*rho1(2:nr,ilm,isp)*v1(2:nr,ilm,isp)),isp=1,nsp)]
         if(ilm==1) rvtr(ilm,:)=rvtr(ilm,:)- srfpi*2d0*z*[(sum(rwgt(2:nr)*rho1(2:nr,ilm,isp)/rofi(2:nr)),isp=1,nsp)]
         rvsm(ilm,:)= [(sum(rwgt(2:nr)*rho2(2:nr,ilm,isp)*v2(2:nr,ilm,isp)),isp=1,nsp)]
      enddo
      vefv1=sum(rvtr)
      vefv2=sum(rvsm)
      rhoexc = rep1  - rep2
      rhoex  = rep1x - rep2x
      rhoec  = rep1c - rep2c
      rhovxc = rmu1  - rmu2
      valvef = vefv1 - vefv2  ! v1*rho1_val  - v2*rho2_val
      if (ipr>=40) then
        if(nsp == 1) write(stdo,"(/' ilm',09x,'rho*vtrue',07x,'rho*vsm')")
        if(nsp == 2) write(stdo,"(/' ilm',19x,'rho*vtrue',30x,'rho*vsm'/13x, &
             'spin1',7x,'spin2',7x,'tot',11x,'spin1',7x,'spin2',7x,'tot')")
        do ilm= 1,nlml
          do isp=1,nsp
            if(dmax1(dabs(rvsm(ilm,isp)),dabs(rvtr(ilm,isp)))>1d-6) then
              if(nsp==1)write(stdo,"(i4,3x,2f15.6)") ilm,rvtr(ilm,isp),rvsm(ilm,isp)
              if(nsp==2)write(stdo,"(i4,3x,3f12.6,2x,3f12.6,2x)")&
                   ilm,rvtr(ilm,1:nsp),sum(rvtr(ilm,:)),rvsm(ilm,1:nsp),sum(rvsm(ilm,:))
            endif
          enddo
        enddo
      endif
    endblock vefv
    if(ipr>=20) then
       write(stdo,"(' local terms:     true',11x,'smooth',9x,'local')")
       write(stdo,"(' rhoeps:  ',3f15.6/' rhomu:   ',3f15.6)") sum(rep1),sum(rep2),sum(rhoexc),rmu1(1),rmu2(1),rhovxc(1)
       if(nsp==2) write(stdo,"(' spin2:   ',3f15.6/' total:   ',3f15.6)")rmu1(2),rmu2(2),rhovxc(2),sum(rmu1),sum(rmu2),sum(rhovxc)
       write(stdo,"(' val*vef  ',3f15.6/' val chg: ',3f15.6)") vefv1,vefv2,valvef,qv1,qv2,qloc
       if(nsp==2) write(stdo,"(' val mmom: ',f15.6,'  core mmom:',f11.6)") aloc,alocc !       write(stdo,"(' core chg:',3f15.6)") qcor1 !,qcor2,qlocc
    endif
    call tcx('locpt2')
  end subroutine locpt2
  subroutine poinsp(z,vval,nlm,v,rofi,rho,rho0,nr,rhoves,rhves1, vnucl,vsum) !- Solves non-spherical poisson Equation inside sphere
    use m_ll,only:ll
    !i   z     :nuclear charge
    !i   vval  :boundary conditions: potential for channel ilm = vval(ilm)
    !i   nlm   :L-cutoff for density
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   rofi  :radial mesh points
    !i   rho   :density, represented as sum_ilm rho(r,ilm) Y_L(ilm)
    !i   rho0  :work array, holding spherical charge density times 4*pi*r*r
    !i   nr    :number of radial mesh points
    !o   v     :electrostatic potential, satisfying boundary c. vval above
    !o         :It does not include nuclear contribution -2z/r
    !o   rhoves:electrostatic energy
    !o   rhves1:contribution to electrostatic energy from nonspherical rho
    !o   vnucl :potential at the nucleus
    !o   vsum  :integral over that potential which is zero at rmax.
    !r Remarks
    !r   Poisson's equation is d2u/dr2 = u*l*(l+1)/r**2 - 8pi*rho/r
    implicit none
    integer :: nlm,nr,nsp,ir,ilm,l,lp1
    real(8) :: z,a,b,rhoves,rhves1,vnucl,vsum,rho(nr,nlm),vval(nlm),v(nr,nlm),rofi(nr),        vhrho(2),rho0(nr)
    real(8) :: ea,a2b4,rpb,rmax,fllp1,df,r, drdi,srdrdi,g,x,f,y2,y3,y4,vnow,vhom,alfa,sum,wgt
    real(8),parameter::fpi = 16d0*datan(1d0),srfpi = dsqrt(fpi),y0 = 1d0/srfpi,atepi = 2d0*fpi
    nsp = 1
    a= log(rofi(3)/rofi(2)-1d0) !rofi(i)=b(e^(a*(i-1))-1). We get a and b
    b= rofi(3)/(exp(2d0*a)-1d0)
    a2b4 = a*a/4d0
    rmax = rofi(nr)
    rho0(1:nr) = rho(1:nr,1)*srfpi
    call poiss0(z,rofi,rho0,nr,vval(1)*y0,v,vhrho,vsum,nsp) !Call poiss0 for l=0 to get good pot for small r 
    rhoves = vhrho(1)
    vnucl = v(1,1)
    v(1:nr,1) = v(1:nr,1)*srfpi
    rhves1 = 0d0
    ilmloop: do  80  ilm = 2, nlm
       l = ll(ilm)
       lp1 = l+1
       fllp1 = l*(l+1d0)
       ! --- Numerov for inhomogeneous solution --- ---
       v(1,ilm) = 0d0
       df = 0d0
       do ir = 2, 3
          r = rofi(ir)
          drdi = a*(r+b)
          srdrdi = dsqrt(drdi)
          g = (r**lp1)/srdrdi
          v(ir,ilm) = r**l
          x = fllp1*drdi*drdi/(r*r) + a2b4
          f = g*(1d0-x/12d0)
          if (ir == 2) y2 = -atepi*rho(2,ilm)*drdi*srdrdi/r
          if (ir == 3) y3 = -atepi*rho(3,ilm)*drdi*srdrdi/r
          df = f-df
       enddo
       do ir=4,nr
          r = rofi(ir)
          drdi = a*(r+b)
          srdrdi = dsqrt(drdi)
          y4 = -atepi*drdi*srdrdi*rho(ir,ilm)/r
          df = df+g*x+(y4+10d0*y3+y2)/12d0
          f = f+df
          x = fllp1*drdi*drdi/(r*r) + a2b4
          g = f/(1d0-x/12d0)
          v(ir,ilm) = g*srdrdi/r
          y2 = y3
          y3 = y4
       enddo
       ! --- Add homogeneous solution --- ---
       vnow = v(nr,ilm)
       vhom = rmax**l
       alfa = (vval(ilm)-vnow)/vhom
       sum = 0d0
       do ir = 2, nr
          r = rofi(ir)
          wgt = 2*(mod(ir+1,2)+1)
          if (ir == nr) wgt = 1d0
          v(ir,ilm) = v(ir,ilm) + alfa*(r**l)
          sum = sum+wgt*(r+b)*rho(ir,ilm)*v(ir,ilm)
       enddo
       rhves1 = rhves1 + a*sum/3d0
       v(1,ilm) = 0d0
80  enddo ilmloop
    rhoves = rhoves + rhves1
  end subroutine poinsp
  subroutine elfigr(nc,stdo,z,efg1)    !- Computation of electric field gradient at Nucleus
    !i Inputs
    !i   nc    :number of classes or sites
    !i   stdo  :standard output
    !i   z     :nuclear charge
    ! o Inputs/Outputs
    ! o  efg1  :input:  l=2 part of electric field at nucleus
    ! o        :output: Electric field gradient
    implicit none
    integer :: nc,stdo
    real(8) :: z(nc),efg1(5,*)
    integer :: i,ifi,ic,ifesn,j,n,ier,imax,ixx,iyy
    real(8) :: v(3,3),vi(3,3),d(3),e(3),e2(3),tau(2,3)
    real(8) :: tr(3,3),ti(3,3),da(3)
    real(8) :: conv1,emmsfe,emmssn,Qfe,Qsn,pi,f0,s3,conv2,conv3, ax,bx,cx,dx,ex,dmax,eta,split
    data conv1 /162.1d0/
    data emmsfe,emmssn /4.8d-8,7.97d-8/
    data Qfe,Qsn /0.21d-28,-0.08d-28/
    pi = 4d0*datan(1d0)
    f0 = -dsqrt(15d0/16d0/pi)
    s3 = dsqrt(3d0)
    conv2 = conv1*3d0
    do  i = 1, 1
       ifi = stdo
       if(i == 2) ifi=71
       write(ifi,'(/," Electric Field Gradient:")')
       write(ifi,'(" Site ",5x,"Tensor axes",5x,"esu/cm^2",4x,"V/m^2 ",4x,"eta",4x,"line splitting")')
       write(ifi,'("      ",5x,"           ",5x," x10^13 ",4x,"x10^19",4x,"   ",4x,"    (mm/s)")')
    enddo
    do  ic = 1, nc
       if (z(ic) > 0.01d0) then
          ifesn = 1
          if (dabs(z(ic)-26d0) < 0.01d0) then
             conv3 = conv2*1.0d19*0.5d0*Qfe/emmsfe
          else if (dabs(z(ic)-50d0) < 0.01d0) then
             conv3 = conv2*1.0d19*0.5d0*Qsn/emmssn
          else
             ifesn = 0
             conv3 = 0d0
          endif
          ax = efg1(1,ic)
          bx = efg1(2,ic)
          cx = efg1(3,ic)
          dx = efg1(4,ic)
          ex = efg1(5,ic)
          v(1,1) =  2d0*ex-2d0/s3*cx
          v(2,2) = -2d0*ex-2d0/s3*cx
          v(3,3) = 4d0/s3*cx
          v(1,2) = 2d0*ax
          v(1,3) = 2d0*dx
          v(2,3) = 2d0*bx
          v(2,1) = v(1,2)
          v(3,1) = v(1,3)
          v(3,2) = v(2,3)
          do  i = 1, 3
             do  j = 1, 3
                v(i,j) = f0*v(i,j)
                vi(i,j) = 0d0
                tr(i,j) = 0d0
                ti(i,j) = 0d0
             enddo
             tr(i,i) = 1d0
          enddo
          n = 3
          diag: block
            complex(8):: tt(3,3),oo(3,3),vv(3,3),img=(0d0,1d0)
            integer:: nmx=3,nev
            oo(:,1) = [complex(8):: 1d0,0d0,0d0]
            oo(:,2) = [complex(8):: 0d0,1d0,0d0]
            oo(:,3) = [complex(8):: 0d0,0d0,1d0]
            vv=v+img*vi
            call diagcv(oo,vv,tt,n,d,nmx,1d99,nev)
            tr=dreal(tt)
          endblock diag
          do  i = 1, 3
             da(i) = dabs(d(i))
          enddo
          dmax = 0d0
          imax = 0
          do  i = 1, 3
             if (da(i) > dmax) then
                dmax = da(i)
                imax = i
             endif
          enddo
          !  d(imax) is field gradient (Vzz)
          ixx = mod(imax,3)+1
          iyy = mod(imax+1,3)+1
          eta = 0d0
          if(dabs(d(imax)) > 1.d-2) eta = dabs((d(ixx)-d(iyy))/d(imax))
          split = conv3*da(imax)*dsqrt(1d0+eta**2/3d0)
          do  i = 1, 3
             if(i == 1) then
                write(ifi,'(i4,3x,3f6.3,2x,f8.2,2x,f8.2,5x,f6.3,5x,f8.5)') ic,(tr(j,i),j=1,3),conv1*d(i),conv2*d(i),eta,split
             else
                write(ifi,'(4x,3x,3f6.3,2x,f8.2,2x,f8.2,5x)') (tr(j,i),j=1,3),conv1*d(i),conv2*d(i)
             endif
          enddo
          write(ifi,'(1x)')
       endif
    enddo
  end subroutine elfigr
  !$$$!!---------------------
  !$$$      subroutine adddipole(v,rofi,nr,nlml,nsp, idipole, basr)
  !$$$!! v=v+x (or y or z). to calculate dipole matrix of Hamiltonian
  !$$$      integer:: nr,nlml,nsp,idipole,isp,nrmx
  !$$$      real(8):: v(nr,nlml,nsp),rofi(nr),basr
  !$$$      real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),srfpi3 = dsqrt(4d0*pi/3d0)
  !$$$      do isp=1,nsp
  !$$$        v(1:nr,1,isp) = v(1:nr,1,isp) + basr * srfpi
  !$$$        if(idipole==1) v(1:nr,4,isp) = v(1:nr,4,isp)+rofi(1:nr)* srfpi3 !x
  !$$$        if(idipole==2) v(1:nr,2,isp) = v(1:nr,2,isp)+rofi(1:nr)* srfpi3 !y
  !$$$        if(idipole==3) v(1:nr,3,isp) = v(1:nr,3,isp)+rofi(1:nr)* srfpi3 !z
  !$$$      enddo
  !$$$      end
  subroutine wrhomt(filnam,descr,ib,rhol,rofi,nr,nlml,nsp)! Write augmented charge density or potential to file for 1 sphere
    !i   rhol  :sphere density, tabulated on a radial mesh
    !i         :rl = full charge density * r**2, written as rl = sum_ilm rhol(ilm) Y_L(ilm)
    !i   rofi  :radial mesh points
    !i   nr    :number of radial mesh points
    !i   nlml  :L-cutoff for charge density on radial mesh
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !Radial mesh rofi and sphere density rhol are written to filnam.ib
    integer :: ib,nr,nlml,nsp,jfi,iprint
    real(8) :: rofi(nr), rhol(nr,nlml,nsp)
    character :: descr*(*),filnam*(*)
    character(8) :: xtxx
    character(256):: fnam
    fnam=trim(filnam)//xtxx(ib)
    if(iprint()>30) write(6,"(a)")' Writing Sphere '//trim(descr)//' to '//trim(fnam)
    open(newunit=jfi,file=trim(fnam),form='unformatted')
    write (jfi) nr,1,1,0
    write (jfi) rofi
    write (jfi) nr,nlml*nsp,1,0,nsp
    write (jfi) rhol
    close(jfi)
  end subroutine wrhomt
end module m_locpot
