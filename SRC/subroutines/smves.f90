!>smooth part of electrostatic potential.
module m_smves 
  use m_lmfinit,only: rmt_i=>rmt,lmxl_i=>lmxl,spec_a,rg_i=>rg,nlmxlx,nbas
  use m_ll,only:ll
  private
  public smves
contains
  subroutine smves(qmom,gpot0,vval,hpot0,smrho,smpot,vconst,smq,qsmc,f,rhvsm0,rhvsm,zsum,vrmt,qbg) ! Electrostatic potential of the 0th component (represented by PlaneWave + Gaussians + smHankels)
    use m_supot,only: iv_a_okv, rv_a_ogv
    use m_lgunit,only:stdo
    use m_lmfinit,only: rv_a_ocy,nsp,ispec
    use m_lattic,only: vol=>lat_vol
    use m_supot,only: ng=>lat_ng,n1,n2,n3
    use m_MPItk,only: master_mpi
    use m_ext,only: sname
    use m_esmsmves,only: esmsmves
!    use m_hansr,only:corprm
    use m_vesgcm,only: vesgcm
    use m_ftox
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
    !i   qmom  :multipole moments (rhomom.f)!   NOTE:2023-10-5 Full moment Q_aL^Zc in Eq.(25).
    !i   smrho :smooth density on real-space mesh
    !i   qbg   : back ground charge
    !io  vconst:constant potential to be added to total
    !io        vconst is the average estat at the MT boundary.
    !o Outputs (see also Remarks)
    !o   gpot0 :integrals of compensating gaussians  G_RL * ves0 !ves0 is the electro static part of V of Eq.(34)
    !           Note that ves made from two contribution ves0_n0 and ves_(n^c_sH +Gaissians). See Eq.(28) (typo n^Zc_0==>n^c_sm,a).
    !o         :For accuracy, integral is split into 
    !o         :  G_RL*ves0_n0 (vesgcm) +  G_RL*ves0_(n^c_sH,a+Gaussians) (ugcomp)
    !o   vval  :coffs to YL expansion of es potential at MT boundary. This is needed to know eigenvalues for pnunew.
    !o   hpot0 : similar with gpot0 but for smHankels
    !o   smpot : ves0 (incluing mesh expansion of contributions from gaussians and smHankels).
    !o   smq   : integral of smooth density n0
    !o   qsmc  : \int n^c_sH
    !o   f     : electrostatic contribution to force.
    !o   rhvsm : electrostatic energy of 0th comp + vconst*smq
    !o   vrmt  : electrostatic potential at rmt, with G=0 term in smpot=0
    implicit none
    integer:: ib , igetss , ilm , ipr , iprint , is ,  lfoc, &
         lgunit , lmxl , m , nlm , j1 , j2 , j3,iwdummy, ifivsmconst
    real(8):: qmom(nlmxlx,nbas) , f(3,nbas) , gpot0(nlmxlx,nbas) , vval(nlmxlx,nbas) , hpot0(nbas) &
         , vrmt(nbas),qsmc,smq,rhvsm,vconst,zsum,zvnsm,qbg
    real(8):: ceh,cofg,cofh,dgetss,hsum,qcorg,qcorh,qsc, &
         rfoc,rmt,s1,s2,sbar,sumx,sum1,sum2,u00,u0g,ugg,usm,vbar,z,R,eint, rhvsm0
    complex(8):: smrho(n1,n2,n3,nsp),smpot(n1,n2,n3,nsp)
    complex(8) ,allocatable :: cg1_zv(:), cgsum_zv(:), cv_zv(:)
    real(8),parameter::pi = 4d0*datan(1d0), srfpi = dsqrt(4d0*pi), y0= 1d0/srfpi
    call tcn('smves')
    ipr = iprint()
    if(nsp==2) smrho(:,:,:,1)=smrho(:,:,:,1)+smrho(:,:,:,2)!Electrostatics is only on total density
    allocate(cv_zv(ng),cg1_zv(ng),cgsum_zv(ng))
    f=0d0
    call fftz3(smrho,n1,n2,n3,n1,n2,n3,1,0,-1) ! FT of smooth density to reciprocal space
    call vesft ( ng,rv_a_ogv,iv_a_okv,cv_zv,smrho,smpot,u00 )! Estatic potential of smooth density without gaussians
    call vesgcm (qmom,ng,rv_a_ogv,iv_a_okv,cv_zv,cg1_zv,cgsum_zv ,&
         smpot,f,gpot0,hpot0,qsmc,zsum,vrmt ) !Add estatic potential of  gaussians and smHankels to smpot
    if(ipr>=40)write(stdo,"(/' after vesgcomp: forces are:'/(i4,3f12.6))")(ib,(f(m,ib),m=1,3),ib=1,nbas)
    ESMsuppliedfromMOBATA: block
      call esmsmves(qmom, ng,rv_a_ogv,iv_a_okv,cv_zv,cg1_zv,cgsum_zv,smrho, qbg, smpot,f,gpot0,hpot0,qsmc,zsum,vrmt )
    endblock ESMsuppliedfromMOBATA
    call mshvmt(ng,rv_a_ogv,iv_a_okv, cv_zv, smpot,vval )
    call symvvl(vval,vrmt) !Compute e.s. potential at MT boundary vrmt
    deallocate(cgsum_zv,cg1_zv,cv_zv)
    vbar = 0d0
    sbar = 0d0
    do  ib = 1, nbas
       rmt = rmt_i(ispec(ib))
       vbar = vbar + rmt**2 * vrmt(ib) 
       sbar = sbar + rmt**2
    enddo
    vbar = vbar/sbar ! vbar =avg v(RMT) and optionally added to vconst
    vconst = -vbar 
    if(ipr>=20) write (stdo,"(' smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=',f9.6)") vconst
    if(master_mpi) then
       open(newunit=ifivsmconst,file='vessm.'//trim(sname))
       write(ifivsmconst,"(d23.15,a)") vconst, '!-(averaged electro static potential at MTs)'
       close(ifivsmconst)
    endif
    ! ... Adjust vbar, vval, gpot0 by vconst
    do  ib = 1, nbas
       lmxl = lmxl_i(ispec(ib))
       if (lmxl > -1) then
          nlm = (lmxl+1)**2
          vrmt(ib) = vrmt(ib) + vconst
          vval(1,ib) = vval(1,ib) + vconst/y0
          gpot0(1,ib) = gpot0(1,ib) + vconst/y0
       endif
    enddo
    if (ipr >= 40) then
       write (stdo,"(' average electrostatic potential at MT boundaries after shift')")
       write(stdo, "(a)") ' Site    ves'
       do ib=1,nbas
          write(stdo,"(i4, f12.6)")ib,vrmt(ib)
       enddo
    endif
    ! ... Back transform of density and potential to real-space mesh
    call fftz3(smrho,n1,n2,n3,n1,n2,n3,1,0,1)
    call fftz3(smpot,n1,n2,n3,n1,n2,n3,1,0,1)
    smrho(:,:,:,1)=smrho(:,:,:,1)+qbg/vol !Add background to smrho
    if (qbg /= 0) then
       R = (3d0/pi/4d0*vol)**(1d0/3d0)
       eint = qbg*2*9d0/10d0/R
       if(ipr>=30) write(stdo,ftox)' cell interaction energy from homogeneous background (q=',ftof(qbg),') is ',ftof(eint)
    endif
    smq = dreal(sum(smrho(:,:,:,1)))               *vol/(n1*n2*n3)  !Integral n0
    s1  = dreal(sum(smrho(:,:,:,1)*smpot(:,:,:,1)))*vol/(n1*n2*n3)  !Integral n0*smpot
    call ugcomp(qmom,gpot0,hpot0,ugg,f) !Gaussian and smHankel parts added
    if(ipr>=50)write (stdo,"(/' after ugcomp: forces are'/(i4,3f12.6))")(ib,(f(m,ib),m=1,3),ib=1,nbas)
    ! --- Collect energy terms; make zvnuc for smooth problem ---
    zvnsm = 0d0
    rhvsm0= s1 + vconst*smq !\int smpot*n1 !smpot is the 0th component of the electrostatic potential part in Eq.(34).
    rhvsm = s1 + vconst*smq
    sumx = 0d0
    do  ib = 1, nbas
       is = ispec(ib)
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       lmxl = lmxl_i(is)
       if (lmxl > -1) then
          nlm = (lmxl+1)**2
          hsum = -srfpi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh ! hsum = integral of charge in sm. Hankel
          hpot0(ib) = hpot0(ib) + vconst*hsum
          zvnsm = zvnsm  + cofh*hpot0(ib) !+ (qcorg-z)*y0*gpot0(iv0+1)
          do  ilm = 1, nlm
             rhvsm = rhvsm + qmom(ilm,ib)*gpot0(ilm,ib)
             sumx = sumx   + qmom(ilm,ib)*gpot0(ilm,ib)
          enddo
       endif
    enddo
    rhvsm = rhvsm + zvnsm
    usm = 0.5d0*rhvsm
    if (ipr >= 30) write (stdo,"('   smooth rhoves',f14.6,'   charge',f13.6)") usm,smq
    smrho(:,:,:,1)=smrho(:,:,:,1)-qbg/vol! ... subtract background
    smq=smq-qbg
    if (nsp == 2) then!     Restore spin 1 density, copy potential to second spin channel
       smrho(:,:,:,1)= smrho(:,:,:,1)-smrho(:,:,:,2)
       smpot(:,:,:,2)= smpot(:,:,:,1)
    endif
    call tcx('smves')
  end subroutine smves
  subroutine mshvmt(ng,gv, kv,cv,smpot,vval)!- Makes potential at MT surfaces given potential on a uniform mesh
    use m_struc_def
    use m_lattic,only:lat_plat,rv_a_opos
    use m_lmfinit,only:lat_alat,ispec
    use m_supot,only: n1,n2,n3
    use m_ropyln,only: ropyln
    use m_ropbes,only: ropbes
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !i   kv    :indices for gather/scatter operations (gvlist.f)
    !i   cv    :work array
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
    !i   smpot :estat potential
    !o Outputs
    !o   vval  :coffs to YL expansion of es potential at MT boundary for each site, computed from the mesh density.
    !r Remarks
    !r   A PW exp(i.q.r) has a one-center expansion at radius r
    !r      sum_L C_L Y_L(r) where C_L = 4 pi i^l j_l(|rq|) Y_L(q)
    !r   Routine symvvl symmetrizes the vval generated here.
    !b Bugs
    !b   Possible to make ves for sites with lmxl=-1, which tells
    !b   value of ves at point.  However, vval doesn't have the
    !b   space allocated.  So skip for now
    !u Updates
    !u   01 Jul 05 handle sites with lmxl=-1
    !u   22 Aug 01 Newly created.
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ng,kv(ng,3)
    real(8):: gv(ng,3),vval(nlmxlx,nbas)
    double complex smpot(n1,n2,n3),cv(ng)
    integer :: i,ib,is,lmxx,nlmx,lmxl,nlm,m,ilm,l,ipr
    real(8) :: alat,pi,tpiba,tau(3),rmt,fac,plat(3,3)
    double complex vvali,fprli
    parameter (lmxx=6, nlmx=(lmxx+1)**2)
    real(8),allocatable:: phil(:,:),yl(:,:)
    real(8),allocatable:: gv2(:,:),agv(:),cgp(:),sgp(:)
    call tcn('mshvmt')
    allocate(phil(ng,0:lmxx),yl(ng,nlmx))
    allocate(gv2(ng,3),agv(ng),cgp(ng),sgp(ng))
    call getpr(ipr)
    pi = 4d0*datan(1d0)
    alat=lat_alat
    plat =lat_plat
    tpiba = 2*pi/alat
    call gvgetf(ng,1,kv,n1,n2,n3,smpot,cv)
    call dpcopy(gv,gv2,1,3*ng,tpiba)
    call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxx,ng,yl,agv) !YL(G)*G**l, agv=|G| for each g ---
    do  i = 1, ng
       agv(i) = sqrt(agv(i))
    enddo
    do  ib = 1, nbas
       is=ispec(ib)
       tau=rv_a_opos(:,ib) 
       rmt=rmt_i(is)
       lmxl=lmxl_i(is)
       if (lmxl == -1) goto 10
       nlm = (lmxl+1)**2
       if (nlm > nlmx) call rxi('mshvmt: increase nlmx to',nlm)
       !       Add a negligibly small amount to rmt to handle case rmt=0
       rmt = rmt+1d-32
       !   --- j_l(|rmt*q|)/rmt**l for each G and l=0..lmax ---
       !       Does not evolve correctly in the correct large r limit
       call ropbes(agv,rmt**2,lmxl,cgp,sgp,phil,ng)
       tau=alat*tau !call dscal(3,alat,tau,1)
       do  i = 1, ng
          fac = -sum(tau*gv2(i,:)) !+tau(2)*gv2(i,2)+tau(3)*gv2(i,3))
          cgp(i) = dcos(fac)
          sgp(i) = dsin(fac)
       enddo
       !   --- Sum_G 4*pi*(i*rmt)**l j_l(|rmt*G|)/(rmt*G)**l YL(G) G**l ---
       ilm = 0
       fprli = 4*pi
       do  l  = 0, lmxl
          do  m = -l, l
             ilm = ilm+1
             vval(ilm,ib) = fprli*sum((phil(2:ng,l)*yl(2:ng,ilm))*(cv(2:ng)*dcmplx(cgp(2:ng),-sgp(2:ng))))
          enddo
          fprli = fprli*(0d0,1d0)*rmt
       enddo
10     continue
    enddo
    deallocate(phil,yl)
    deallocate(gv2,agv,cgp,sgp)
    call tcx('mshvmt')
  end subroutine mshvmt
  subroutine symvvl(vval,vrmt) !Symmetrizes the potential at the MT boundary.
    use m_mksym,only: symops,ag,ngrp,ipc=>iclasst
    use m_struc_def
    use m_lattic,only:lat_plat,rv_a_opos
    use m_lgunit,only:stdo
    use m_lmfinit,only:ispec
    !o Outputs
    ! o  vval  :On input,  unsymmetrized potential
    ! o        :On output, elements of potential for sites in the same
    ! o        :class are symmetrized.
    !o Outputs
    !o   vrmt  :spherical average of potential (i.e. Y0*vval(l=0)) returned
    !o         :for each site.
    !r Remarks
    !r   This routine symmetrizes any vector of the same structure as vval.
    !b Bugs
    !b   Possible to make ves for sites with lmxl=-1, which tells
    !b   value of ves at point.  However, vval doesn't have the
    !b   space allocated.  So skip for now
    !u Updates
    !u   01 Jul 05 handle sites with lmxl=-1
    !u   23 Aug 01 Newly created.
    ! ----------------------------------------------------------------------
    implicit none
    real(8):: vval(nlmxlx,nbas),vrmt(nbas)
    integer :: ic,ib,ilm,mxint,nclass,ipa(nbas),nrclas
    integer :: ips(nbas),lmxl(nbas) !ipc(nbas),
    real(8) :: pos(3,nbas),posc(3,nbas),plat(3,3),pi,y0
    integer:: igetss,nlml ,ipr,nn,ibas
    real(8) ,allocatable :: qwk_rv(:)
    real(8) ,allocatable :: sym_rv(:)
    call tcn('symvvl')
    call getpr(ipr)
    plat=lat_plat
    do ibas=1,nbas
       ips(ibas)  = ispec(ibas) 
       pos(:,ibas)= rv_a_opos(:,ibas)
    enddo
    nclass = maxval(ipc)
    do  ib = 1, nbas
       lmxl(ib)=lmxl_i(ips(ib))
    enddo
    do  ic = 1, nclass
       !   ... Make nrclas,ipa,posc !       call psymr0(lmxl,ic,nbas,ipc,pos,posc,ipa,nrclas)
       nrclas = 0
       do  ib = 1, nbas
          if (ipc(ib) == ic) then
             nrclas = nrclas+1
             ipa(nrclas) = ib
             posc(:,nrclas) = pos(:,ib)
          endif
       enddo
       if (nrclas > 0) then
          ib = findloc([(ipc(ibas)==ic,ibas=1,nbas)],dim=1,value=.true.) 
          if (lmxl(ib) > -1) then
             nlml = (lmxl(ib)+1)**2
             if (ipr >= 50) write(stdo,800) ic,nrclas,nlml
800          format(' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)
             call symqmp ( nrclas,nlml,nlml,plat,posc,ngrp,symops,ag,ipa,vval,nn )!,ipa,sym_rv,vval,nn )
          endif
       endif
    enddo
    ! ... Extract vrmt = l=0 term for each site, and printout
    pi = 4d0*datan(1d0)
    y0 = 1d0/dsqrt(4d0*pi)
    do  ib = 1, nbas
       if (lmxl(ib) == -1) cycle
       vrmt(ib) = vval(1,ib)*y0
    enddo
    write(stdo,"(/' site class  ilm      vval',6x,'ves(rmax)')")
    do  ib = 1, nbas
       if (lmxl(ib) == -1) cycle
       nlml = (lmxl(ib)+1)**2
       ic = ipc(ib)
       if(ib==findloc([(ipc(ibas)==ic,ibas=1,nbas)],dim=1,value=.true.).or.ipr>60) then ! (ib == iclbas(ic,ipc)) then !findloc(ipc,value=ic) ) then !
          do ilm = 1, nlml
             if (ilm == 1) then;                    write(stdo,"(i4,2i6,2f12.6)") ib,ic,ilm,vval(ilm,ib),vrmt(ib)
             elseif(dabs(vval(ilm,ib))>1d-6.AND.ipr>=50) then; write(stdo,"(10x,i6,f12.6)") ilm,vval(ilm,ib); endif
          enddo
       endif
    enddo
    call tcx('symvvl')
  end subroutine symvvl
  subroutine ugcomp(qmom,gpot0,hpot0,ugg,f) !Part of the smooth estatic energy from compensating G's alone.
    use m_lgunit,only:stml
    use m_smhankel,only:hhugbl,hgugbl,ggugbl
    use m_lattic,only: rv_a_opos
    use m_lmfinit,only: ispec
!    use m_hansr,only:corprm
    !i   qmom  :multipole moments of on-site densities (rhomom.f). defined at Eq.(25)
    ! o Inputs/Outputs
    ! o  Let n0  = smooth potential without compensating gaussians
    ! o      n0~ = smooth potential with compensating gaussians
    ! o    phi0  = ves[n0]
    ! o    phi0~ = ves[n0~]
    ! o    g_RL  = gaussian in RL channel
    ! o    h_R   = l=0 sm hankel in RL channel, for core density
    ! o  Then:
    ! o  gpot0 :On input, integrals g_RL * phi0
    ! o        :On output, integrals g_RL * phi0~
    ! o  hpot0 :On input, integrals h_R * phi0
    ! o        :On output, integrals h_R * phi0~
    !o Outputs
    !o   ugg   :electrostatic energy integral [n0~-n0]*[phi0~-phi0]
    !i   f     :contribution to forces is added
    !r Remarks
    !u Updates
    !u   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
    !u   15 Feb 02 (ATP) Added MPI parallelization
    !u   22 Apr 00 Adapted from nfp ugcomp
    ! ----------------------------------------------------------------------
    implicit none
    real(8):: qmom(nlmxlx,nbas),gpot0(nlmxlx,nbas),f(3,nbas),hpot0(nbas),ugg
    integer :: ndim,ndim0,i,ib,ilm1,ilm2,is,jb,js,l1,l2, &
         lfoc1,lfoc2,lmax1,lmax2,m,nlm1,nlm2
    parameter (ndim=49, ndim0=2)
    real(8) :: ceh1,ceh2,cof1,cof2,cofg1,cofg2,cofh1,cofh2,fpi, &
         pi,qcorg1,qcorg2,qcorh1,qcorh2,qsc1,qsc2,qm1,qm2,rg1,rg2,rh1, &
         rh2,srfpi,y0,z1,z2
    real(8) :: df(0:20),ff(3),tau1(3),tau2(3)
    double complex s(ndim,ndim),ds(ndim,ndim,3),s0(ndim0,ndim0),ds0(ndim0,ndim0,3)
    integer :: nlmx,npmx,nbmx
    parameter (nlmx=64, nbmx=256)
    real(8) :: xf(3,nbas),xhpot0(nbas), xgpot0(nlmxlx,nbas),xugg
    integer:: ibini,ibend
    call tcn('ugcomp')
    call stdfac(20,df)
    pi = 4d0*datan(1d0)
    fpi = 4d0*pi
    srfpi = dsqrt(fpi)
    y0 = 1d0/srfpi
    ! --- Loop over sites where charge lump making pot is centered ---
    ugg = 0d0
    xugg=0d0   
    xgpot0=0d0 
    xf=0d0     
    xhpot0=0d0
    do ib=1,nbas
       is=ispec(ib) 
       tau1=rv_a_opos(:,ib) 
       lmax1=lmxl_i(is)
       if (lmax1 ==-1) cycle
       rg1=rg_i(is)
       call corprm(is,qcorg1,qcorh1,qsc1,cofg1,cofh1,ceh1,lfoc1, rh1,z1)
       nlm1 = (lmax1+1)**2
       do  jb = 1, nbas!  Loop over sites where charge lump sees the potential
          js=ispec(jb)
          tau2=rv_a_opos(:,jb) 
          lmax2=lmxl_i(js)
          if(lmax2==-1) cycle
          rg2=rg_i(js)
          call corprm(js,qcorg2,qcorh2,qsc2,cofg2,cofh2,ceh2, lfoc2,rh2,z2)
          nlm2 = (lmax2+1)**2
          if (nlm1 > ndim) call rxi('ugcomp: ndim < nlm1=',nlm1)
          if (nlm2 > ndim) call rxi('ugcomp: ndim < nlm2=',nlm2)
          call ggugbl(tau1,tau2,rg1,rg2,nlm1,nlm2,ndim,ndim,s,ds)  !gaussian-gaussian
          ff = 0d0
          do  ilm1 = 1, nlm1
             l1 = ll(ilm1)
             qm1 = qmom(ilm1,ib)
             cof1 = qm1*fpi/df(2*l1+1)
             do  ilm2 = 1, nlm2
                l2 = ll(ilm2)
                qm2 = qmom(ilm2,jb)
                cof2 = qm2*fpi/df(2*l2+1)
                xugg = xugg + cof1*cof2*s(ilm1,ilm2)
                xgpot0(ilm2,jb) = xgpot0(ilm2,jb) + s(ilm1,ilm2)*cof1*fpi/df(2*l2+1)
                ff = ff + 0.5d0*cof1*cof2*ds(ilm1,ilm2,1:3) !Forces
             enddo
          enddo
          if (lfoc1 > 0 .OR. lfoc2 > 0) then ! Additional h*h, h*g, g*h terms for foca ---
             call hhugbl(0,tau1,tau2,[rh1],[rh2],[ceh1],[ceh2],1,1,ndim0,ndim0, s0,ds0)!smHankel-smHankel
             xugg = xugg + cofh1*s0(1,1)*cofh2
             xhpot0(jb) = xhpot0(jb) + cofh1*s0(1,1) !smooth Hankel part
             ff = ff + 0.5d0*cofh1*cofh2*ds0(1,1,1:3)
             call hgugbl(tau1,tau2,rh1,rg2,ceh1,1,nlm2,ndim,ndim, s,ds)!gaussian-smHankel
             do  ilm2 = 1, nlm2
                l2 = ll(ilm2)
                qm2 = qmom(ilm2,jb)
                cof2 = qm2*fpi/df(2*l2+1)
                xugg = xugg + cofh1*s(1,ilm2)*cof2
                ff = ff + 0.5d0*cofh1*cof2*ds(1,ilm2,1:3)
                xgpot0(ilm2,jb) = xgpot0(ilm2,jb) + s(1,ilm2)*cofh1*fpi/df(2*l2+1)
             enddo
             call hgugbl(tau2,tau1,rh2,rg1,ceh2,1,nlm1,ndim,ndim, s,ds) !gaussian-smHamel
             do  ilm1 = 1, nlm1
                l1 = ll(ilm1)
                qm1 = qmom(ilm1,ib)
                cof1 = qm1*fpi/df(2*l1+1)
                xugg = xugg + cof1*s(1,ilm1)*cofh2
                ff = ff - 0.5d0*cof1*cofh2*ds(1,ilm1,1:3)
                xhpot0(jb) = xhpot0(jb) + cof1*s(1,ilm1) !Gaussian part
             enddo
          endif
          if (jb /= ib) then
             xf(1:3,ib) = xf(1:3,ib) - ff
             xf(1:3,jb) = xf(1:3,jb) + ff
          endif
       enddo
    enddo
    f(1:3,:) = f(1:3,:) + xf(1:3,:)
    hpot0(:) = hpot0(:) + xhpot0(:)
    gpot0(:,:) = gpot0(:,:) + xgpot0(:,:)
    ugg = ugg + xugg
    call tcx('ugcomp')
  end subroutine ugcomp
  subroutine vesft(ng,gv,kv,cv,smrho,smpot,ssum)!Make electrostatic potential of density given in recip space
    use m_lmfinit,only:alat=>lat_alat
    use m_lattic,only: vol=>lat_vol
    use m_supot,only: n1,n2,n3
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !i   kv    :indices for gather/scatter operations (gvlist.f)
    !i   cv    :work array holding smrho and smpot in glist form
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
    !i   smrho :FT of smooth density on uniform mesh
    !o Outputs
    !o   smpot :FT of smooth electrostatic potential
    !o   sum   :integral pot*density
    implicit none
    integer :: ng,i, kv(ng,3)
    real(8):: gv(ng,3),ssum, tpiba,g2
    complex(8):: smrho(n1,n2,n3),smpot(n1,n2,n3),cv(ng),ccc(ng)
    real(8),parameter:: pi = 4d0*datan(1d0), pi8  = 8d0*pi
    call tcn('vesft')
    tpiba=2d0*pi/alat
    call gvgetf(ng,1,kv,n1,n2,n3,smrho,cv) ! ... Gather density coefficients
    ccc= [complex(8):: 0d0, ((pi8/(tpiba**2*sum(gv(i,:)**2)))*cv(i),i=2,ng)]
    ssum = vol*sum(ccc(2:ng)*cv(2:ng))
    call gvputf(ng,1,kv,n1,n2,n3,ccc,smpot)! smpot(G) =8pi/G**2 smrho(G)
    call tcx('vesft')
  end subroutine vesft
  subroutine symqmp(nrclas,nlml,nlmx,plat,posc,ngrp,g,ag,ipa,qmp,nn)  !- Symmetrize multipole moments for a single class
    !i Inputs
    !i   nrclas:number of atoms in the ith class
    !i   nlml  :L-cutoff for charge density on radial mesh
    !i   nlmx  :dimensions sym: sym is generated for ilm=1..nlmx
    !i   plat  :primitive lattice vectors, in units of alat
    !i   posc  :work array holding basis vectors for this class
    !i   ngrp  :number of group operations
    !i   g     :point group operations
    !i   ag    :translation part of space group
    !i   qwk   :work array of dimension nlml
    !i   ipa   :ipa(1..nrclas) = table of offsets to qmp corresponding
    !i         :to each member of the class
    ! o Inputs/Outputs
    ! o  qmp   :On input,  unsymmetrized multipole moments
    ! o        :On output, multipole moments are symmetrized
    !o Outputs
    !o   sym   :symmetry projectors for each member of the class
    !o   nn    :number of elements symmetrized
    implicit none
    integer :: nrclas,nlmx,ipa(nrclas),nlml,ngrp,nn,ia,ilm,ixx,iyy(1)
    real(8) :: plat(3,3),posc(3,nrclas),g(3,3,ngrp),ag(3,ngrp),sym(nlmx,nlmx,nrclas),qwk(nlml),qmp(nlmxlx,nbas),wgt,xx,qlat(3,3)
    call tcn('symqmp')
    if (nlml > nlmx) call rxi('symqmp: increase nlmx to',nlml)
    call dinv33(plat,1,qlat,xx)
    call symprj(nrclas,nlmx,ngrp,ixx,iyy,g,ag,plat,qlat,posc,sym)! ... Make the symmetry projectors
    qwk=0d0
    do  ia = 1, nrclas
       call pxsmr1(1d0,1,nlml,1,sym(1,1,ia),qmp(:,ipa(ia)),qwk,nn) ! ... Accumulate symmetrized qmpol on first site
    enddo
    wgt = nrclas
    do  ia = 1, nrclas ! ... Rotate and copy to all sites in class
       qmp(:,ipa(ia))=0d0
       call pysmr1(wgt,1,nlml,1,sym(1,1,ia),qwk,qmp(:,ipa(ia)),nn)
    enddo
    nn = 0
    do  ilm = 1, nlml
       if(dabs(qmp(ilm,ipa(1))) > 1d-6) nn = nn+1
    enddo
    call tcx('symqmp')
  end subroutine symqmp
end module m_smves
