!>smooth part of electrostatic potential.
module m_smves 
  use m_lmfinit,only: rmt_i=>rmt,lmxl_i=>lmxl,spec_a,rg_i=>rg
  use m_ll,only:ll
  private
  public smves
contains
  subroutine smves(qmom,gpot0,vval,hpot0,sgp0,smrho,smpot,vconst,smq,qsmc,f,rhvsm,zvnsm,zsum,vrmt,qbg) ! Electrostatic potential of the 0th component (represented by PlaneWave + Gaussians + smHankels)
    use m_supot,only: iv_a_okv, rv_a_ogv
    use m_lmfinit,only: rv_a_ocy,nsp,stdo,nbas,ispec
    use m_lattic,only: vol=>lat_vol
    use m_supot,only: ng=>lat_ng,n1,n2,n3
    use m_MPItk,only: master_mpi
    use m_ext,only: sname
    use m_esmsmves,only: esmsmves
    use m_hansr,only:corprm
    use m_vesgcm,only: vesgcm
    use m_ftox
    !i   nbas  :size of basis
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
    !i   smrho :smooth density on real-space mesh
    !i   qbg   : back ground charge
    !io  vconst:constant potential to be added to total
    !io        vconst is the average estat at the MT boundary.
    !o Outputs (see also Remarks)
    !o   gpot0 :integrals of compensating gaussians g_RL * phi0~
    !o         :For accuracy, integral is split into
    !o         :g_RL phi0 (vesgcm) + g_RL phi [n0~-n0] (ugcomp)
    !o         :vesgcm projects g_RL to the mesh to do the integral
    !o         :ugcomp does its integrals analytically (structure constants)
    !o         :NB: There is a local analog of gpot0 generated in locpt2.
    !o   vval  :coffs to YL expansion of es potential at MT boundary
    !o   hpot0 :integrals of semicore smooth Hankels * phi0~
    !o   sgp0  :sgp0 = sum_RL integral qmom_RL g_RL phi0~
    !o         :     = integral [n0~-n0 ] phi0~
    !o   smpot :smooth potential phi0~ (includes compensating gaussians)
    !o   smq   :integral of smooth density n0
    !o   qsmc  :pseudocore charge
    !o   f     :electrostatic contribution to force.
    !o   rhvsm :integral n0~ [phi0~ + vconst]
    !o         :(electrostatic energy of sm. density n0~) + vconst*smq
    !o   zvnsm :integral (qcorg-z + rhoc) phi0~
    !o   vrmt  :electrostatic potential at rmt, with G=0 term in smpot=0
    !l Local variables
    !l   u00   :integral n0 phi[n0] = n0 phi0
    !l   u0g   :integral n0 [phi0~-phi0]
    !l   ugg   :integral [n0~-n0] [phi0~-phi0]
    !l         :Note: ugg is not used.
    !r Remarks
    !r  The total density is a sum of three terms,
    !r
    !r    n0(mesh) + sum_RL (n_RL(r) - n0_RL(r))
    !r
    !r  The first term is the smooth density on a mesh of points; the
    !r  second is the true density and is defined on a radial mesh for each
    !r  sphere; the last is the 1-center expansion of the smooth density on
    !r  the radial mesh.  (Note: because of l-truncation, n0_R(r) is not
    !r  identical to the one-center expansion of n0(mesh).  The sum of the
    !r  three terms converges rapidly with l because errors in n_R(r) are
    !r  mostly canceled by errors in n0_R(r).)
    !r
    !r  We add and subtract a set of compensating gaussian orbitals
    !r
    !r    n0 + sum_RL Q_RL g_RL + sum_RL (n_RL(r) - n0_RL(r) - Q_RL g_RL)
    !r
    !r  which render the integral of the local part (the last 3 terms)
    !r  zero in each RL channel.  The g_RL must be localized enough that
    !r  their spillout beyond the MT radius is negligible.
    !r
    !r  We define
    !r
    !r    n0~ = n0 + compensating gaussians sum_RL Q_RL g_RL
    !r
    !r  In the interstitial, the electrostatic potential of n0~ is the true
    !r  estat potential.  The potential of n0 is called phi0 and the
    !r  potential of n0~ is called phi0~.  The total electrostatic energy
    !r  is computed as
    !r
    !r    the electrostatic energy of  n0~ + integral n0*vconst +
    !r    the electrostatic energy of (neutral) local parts
    !r
    !r  vconst may either be passed as an input (mode=0) or it is
    !r  generated here as the average ves(RMT).
    !r  This routine computes the estat potential and energy from the
    !r  first two terms.  Some variables used in smves and its subroutines:
    !r    Let n0  = smooth density without the compensating sum_RL Q_RL g_RL
    !r        n0~ = n0 + sum_RL Q_RL g_RL
    !r      phi0  = ves[n0]
    !r      phi0~ = ves[n0~]
    !r      g_RL  = gaussian in RL channel
    !r      h_R   = l=0 sm hankel in RL channel, (represents core densities)
    !r    qmom_RL = multipole moment in RL channel of (n_R(r) - n0_R(r))
    !r              so that int n_RL(r)-n0_RL(r) = qmom_RL * g_RL(r)
    !r      gpot0 = vector of integrals g_RL * phi0~
    !r            =  integral g_RL * (phi0 = phi[n0])
    !r              +integral g_RL * (phi0~-phi0 = phi[n0~-n0])
    !r               The integral is partitioned to minimize mesh errors.
    !r               The first part is done by projecting g_RL to a mesh
    !r               and integrating the product g_RL*phi0 on the mesh
    !r               The second is done analytically by structure constants
    !r      hpot0 = integrals h_R * phi0~ (contributions from core)
    !r            = integrals h_R * (phi0 = phi[n0])
    !r             +integrals h_R * (phi0~-phi0 = phi[n0~-n0])
    !r       u00   :integral n0 phi[n0] = integral n0 phi0
    !r       u0g   :integral n0 [phi0~-phi0]
    !r       sgp0  :integral [n0~-n0] phi0~
    !r   Therefore :u00 + u0g + sgp0 = integral n0~ phi0~
    !r       smq   :integral n0
    !r       vconst:constant potential to be added to total.
    !r             :It is computed from average (v(RMT))
    !r       rhvsm :u00 + u0g + sgp0 + vconst*smq
    !r             := integral n0~ phi0~ + vconst*smq
    !r       zvnsm :integral core density * phi0~
    !r
    !r  Subroutines called by smves:
    !r    vesft    computes the electrostatic potential of n0 = phi0
    !r             (i.e. without the compensating gaussians).  This
    !r             is pretty trivial, since nabla^2 -> G^2 in G-space
    !r
    !r    vesgcm   1. makes the first term in gpot0
    !r                = integral g_RL * (phi0 = phi[n0])
    !r             2. makes the first term in hpot0
    !r             3. adds ves[n0~-n0] to the mesh estat potential
    !r
    !r    ugcomp   1. makes the second term in gpot0
    !r             2. makes the second term in hpot0
    !u Updates
    !b Bugs
    !b   Possible to make vval(l=0) for sites with lmxl=-1, which tells
    !b   value of ves at point.  However, vval doesn't have the
    !b   space allocated.  So skip for now.
    !u Updates
    !u   01 Jul 05 handle sites with lmxl=-1
    !u   19 Sep 02 (WRL) Added background term
    !u   24 Aug 01 Extended to calc vval.  Altered argument list.
    !u   20 Apr 01 Generates vrmt
    !u   21 Jun 00 spin polarized
    !u   22 Apr 00 Adapted from nfp ves_smooth.f
    ! ----------------------------------------------------------------------
    implicit none
    integer:: ib , igetss , ilm , ipr , iprint , is , iv0 , lfoc, &
         lgunit , lmxl , m , nlm , j1 , j2 , j3,iwdummy, ifivsmconst
    real(8):: qmom(1) , f(3,nbas) , gpot0(1) , vval(1) , hpot0(nbas) &
         , vrmt(nbas),qsmc,smq,rhvsm,sgp0,vconst,zsum,zvnsm,qbg
    real(8):: ceh,cofg,cofh,dgetss,hsum,qcorg,qcorh,qsc, &
         rfoc,rmt,s1,s2,sbar,sumx,sum1,sum2,u00,u0g,ugg,usm,vbar,z,R,eint
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
    call esmsmves(qmom, ng,rv_a_ogv,iv_a_okv,cv_zv,cg1_zv,cgsum_zv, & !ESM method supplied from M.OBATA
         smrho, qbg, smpot,f,gpot0,hpot0,qsmc,zsum,vrmt )
    call mshvmt(ng,rv_a_ogv,iv_a_okv, cv_zv, smpot,vval )
    call symvvl(vval,vrmt) !Compute e.s. potential at MT boundary vrmt
    deallocate(cgsum_zv,cg1_zv,cv_zv)
    ! --- Make 
    vbar = 0d0
    sbar = 0d0
    do  ib = 1, nbas
       rmt = rmt_i(ispec(ib))
       vbar = vbar + rmt**2 * vrmt(ib) 
       sbar = sbar + rmt**2
    enddo
    vbar = vbar/sbar ! vbar =avg v(RMT) and optionally added to vconst
    vconst = -vbar 
    if(ipr>=20) write (stdo,232) vconst
232 format(' smves: Add vconst to Ele.Static Pot. so that avaraged Ves at Rmt is zero: vconst=',f9.6)
    if(master_mpi) then
       open(newunit=ifivsmconst,file='vessm.'//trim(sname))
       write(ifivsmconst,"(d23.15,a)") vconst, '!-(averaged electro static potential at MTs)'
       close(ifivsmconst)
    endif
    ! ... Adjust vbar, vval, gpot0 by vconst
    iv0 = 0
    do  ib = 1, nbas
       lmxl = lmxl_i(ispec(ib))
       if (lmxl > -1) then
          nlm = (lmxl+1)**2
          vrmt(ib) = vrmt(ib) + vconst
          vval(1+iv0) = vval(1+iv0) + vconst/y0
          gpot0(1+iv0) = gpot0(1+iv0) + vconst/y0
          iv0 = iv0 + nlm
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
       if(ipr>=30) write(stdo,ftox)' cell interaction energy from homogeneous'// &
            ' background (q=',ftof(qbg),') is ',ftof(eint)
    endif
    smq = dreal(sum(smrho(:,:,:,1)))               *vol/(n1*n2*n3)  !Integral n0
    s1  = dreal(sum(smrho(:,:,:,1)*smpot(:,:,:,1)))*vol/(n1*n2*n3)  !Integral n0*phi0~
    u0g = s1 - u00
    call ugcomp(qmom,gpot0,hpot0,ugg,f) 
    if(ipr>=50)write (stdo,"(/' after ugcomp: forces are'/(i4,3f12.6))")(ib,(f(m,ib),m=1,3),ib=1,nbas)
    if(ipr>=50)write(stdo,"(' u00,u0g,ugg=',3f14.6)") u00,u0g,ugg
    ! --- Collect energy terms; make zvnuc for smooth problem ---
    zvnsm = 0d0
    rhvsm = u00 + u0g + vconst*smq
    sumx = 0d0
    iv0 = 0
    do  ib = 1, nbas
       is = ispec(ib)
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       lmxl = lmxl_i(is)
       if (lmxl > -1) then
          nlm = (lmxl+1)**2
          hsum = -srfpi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh ! hsum = integral of charge in sm. Hankel
          hpot0(ib) = hpot0(ib) + vconst*hsum
          zvnsm = zvnsm + (qcorg-z)*y0*gpot0(iv0+1) + cofh*hpot0(ib)
          do  ilm = 1, nlm
             rhvsm = rhvsm + qmom(iv0+ilm)*gpot0(iv0+ilm)
             sumx = sumx + qmom(iv0+ilm)*gpot0(iv0+ilm)
          enddo
          iv0 = iv0+nlm
       endif
    enddo
    sgp0 = sumx
    usm = 0.5d0*(rhvsm+zvnsm)
    if (ipr >= 30) write (stdo,"('   smooth rhoves',f14.6,'   charge',f13.6)") usm,smq
    smrho(:,:,:,1)=smrho(:,:,:,1)-qbg/vol! ... subtract background
    smq=smq-qbg
    if (nsp == 2) then!     Restore spin 1 density, copy potential to second spin channel
       smrho(:,:,:,1)= smrho(:,:,:,1)-smrho(:,:,:,2)
       smpot(:,:,:,2)=smpot(:,:,:,1)
    endif
    call tcx('smves')
  end subroutine smves
  subroutine mshvmt(ng,gv, kv,cv,smpot,vval)
    use m_struc_def
    use m_lattic,only:lat_plat,rv_a_opos
    use m_lmfinit,only:lat_alat,nbas,ispec
    use m_supot,only: n1,n2,n3
    use m_ropyln,only: ropyln
    use m_ropbes,only: ropbes
    !- Makes potential at MT surfaces given potential on a uniform mesh
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   ng    :number of G-vectors
    !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
    !i   kv    :indices for gather/scatter operations (gvlist.f)
    !i   cv    :work array
    !i   n1,n2,n3 dimensions of smrho,smpot for smooth mesh density
    !i   smpot :estat potential
    !o Outputs
    !o   vval  :coffs to YL expansion of es potential at MT boundary
    !o         :for each site, computed from the mesh density.
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
    real(8):: gv(ng,3),vval(1)
    double complex smpot(n1,n2,n3),cv(ng)
    integer :: i,ib,is,lmxx,nlmx,iv0,lmxl,nlm,m,ilm,l,ipr
    double precision :: alat,pi,tpiba,tau(3),rmt,fac,plat(3,3)
    double complex vvali,fprli
!    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    parameter (lmxx=6, nlmx=(lmxx+1)**2)
    double precision,allocatable:: phil(:,:),yl(:,:)
    double precision,allocatable:: gv2(:,:),agv(:),cgp(:),sgp(:)
    call tcn('mshvmt')
    allocate(phil(ng,0:lmxx),yl(ng,nlmx))
    allocate(gv2(ng,3),agv(ng),cgp(ng),sgp(ng))
    call getpr(ipr)
    pi = 4d0*datan(1d0)
    alat=lat_alat
    plat =lat_plat
!    ngabc=lat_nabc
    tpiba = 2*pi/alat
    call gvgetf(ng,1,kv,n1,n2,n3,smpot,cv)
    ! --- YL(G)*G**l, agv=|G| for each g ---
    call dpcopy(gv,gv2,1,3*ng,tpiba)
    call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxx,ng,yl,agv)
    do  i = 1, ng
       agv(i) = sqrt(agv(i))
    enddo
    iv0 = 0
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
       !       call dpzero(vval(iv0+1),nlm)
       ilm = 0
       fprli = 4*pi
       do  l  = 0, lmxl
          do  m = -l, l
             ilm = ilm+1
             vvali = 0
             do  i = 2, ng
                vvali = vvali + (phil(i,l)*yl(i,ilm))*(cv(i)*dcmplx(cgp(i),-sgp(i)))
             enddo
             vval(ilm+iv0) = fprli*vvali
          enddo
          fprli = fprli*(0d0,1d0)*rmt
       enddo
       iv0 = iv0 + nlm
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
    use m_lmfinit,only:ispec,nbas
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
    real(8):: vval(1),vrmt(nbas)
    integer :: ic,ib,ilm,mxint,nclass,ipa(nbas),nrclas,iv0
    integer :: ips(nbas),lmxl(nbas) !ipc(nbas),
    double precision :: pos(3,nbas),posc(3,nbas),plat(3,3),pi,y0
    integer:: igetss,nlml ,ipr,jpr,nn,ibas
    real(8) ,allocatable :: qwk_rv(:)
    real(8) ,allocatable :: sym_rv(:)
    call tcn('symvvl')
    call getpr(ipr)
    plat=lat_plat
!    ngrp=lat_nsgrp
    do ibas=1,nbas
       ips(ibas)  = ispec(ibas) 
       pos(:,ibas)= rv_a_opos(:,ibas)
    enddo
    nclass = maxval(ipc)
    do  ib = 1, nbas
       lmxl(ib)=lmxl_i(ips(ib))
    enddo
    do  ic = 1, nclass
       !   ... Make nrclas,ipa,posc
       call psymr0(lmxl,ic,nbas,ipc,pos,posc,ipa,nrclas)
       if (nrclas > 0) then
          ib = findloc([(ipc(ibas)==ic,ibas=1,nbas)],dim=1,value=.true.) !iclbas(ic,ipc) !findloc(ipc,value=ic) !
          if (lmxl(ib) > -1) then
             nlml = (lmxl(ib)+1)**2
             if (ipr >= 50) write(stdo,800) ic,nrclas,nlml
800          format(' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)
             allocate(qwk_rv(nlml))
             allocate(sym_rv(nlml*nlml*nrclas))
             call symqmp ( nrclas,nlml,nlml,plat,posc,ngrp,symops,ag,qwk_rv,ipa,sym_rv,vval,nn )
             if (allocated(sym_rv)) deallocate(sym_rv)
             if (allocated(qwk_rv)) deallocate(qwk_rv)
          endif
       endif
    enddo
    ! ... Extract vrmt = l=0 term for each site, and printout
    pi = 4d0*datan(1d0)
    y0 = 1d0/dsqrt(4d0*pi)
    if (ipr >= 45) write(stdo,221)
221 format(/' site class  ilm      vval',6x,'ves(rmax)')
    iv0 = 0
    do  ib = 1, nbas
       if (lmxl(ib) == -1) goto 10
       nlml = (lmxl(ib)+1)**2
       vrmt(ib) = vval(1+iv0)*y0
       ic = ipc(ib)
       jpr = 0
       if (ipr > 60) jpr = 2
       if( ib == findloc([(ipc(ibas)==ic,ibas=1,nbas)],dim=1,value=.true.)) then ! (ib == iclbas(ic,ipc)) then !findloc(ipc,value=ic) ) then !
          if (ipr >= 45) jpr = 1
          if (ipr >= 50) jpr = 2
       endif
       if (jpr > 0) then
          do  ilm = 1, nlml
             if (ilm == 1) then
                write(stdo,650) ib,ic,ilm,vval(ilm+iv0),vrmt(ib)
             elseif (dabs(vval(ilm+iv0)) > 1d-6  .AND. jpr > 1) then
                write(stdo,651)    ilm,vval(ilm+iv0)
             endif
650          format(i4,2i6,2f12.6)
651          format(10x,i6,f12.6)
          enddo
       endif
       iv0 = iv0 + nlml
10     continue
    enddo
    call tcx('symvvl')
  end subroutine symvvl
  subroutine ugcomp(qmom,gpot0,hpot0,ugg,f) !Part of the smooth estatic energy from compensating G's alone.
    use m_lgunit,only:stml
    use m_smhankel,only:hhugbl,hgugbl,ggugbl
    use m_lattic,only: rv_a_opos
    use m_lmfinit,only: ispec,nbas
    use m_hansr,only:corprm
    !i   qmom  :multipole moments of on-site densities (rhomom.f)
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
    real(8):: qmom(*),gpot0(*),f(3,nbas),hpot0(nbas),ugg
    integer :: ndim,ndim0,i,ib,ilm1,ilm2,is,iv0,jb,js,jv0,nvl,l1,l2, &
         lfoc1,lfoc2,lmax1,lmax2,m,nlm1,nlm2
    parameter (ndim=49, ndim0=2)
    double precision :: ceh1,ceh2,cof1,cof2,cofg1,cofg2,cofh1,cofh2,fpi, &
         pi,qcorg1,qcorg2,qcorh1,qcorh2,qsc1,qsc2,qm1,qm2,rg1,rg2,rh1, &
         rh2,srfpi,y0,z1,z2
    double precision :: df(0:20),ff(3),tau1(3),tau2(3)
    double complex s(ndim,ndim),ds(ndim,ndim,3),s0(ndim0,ndim0), &
         ds0(ndim0,ndim0,3)!,wk(ndim0,ndim0),dwk(ndim0,ndim0,3)
    integer :: nlmx,npmx,ip,mp,nbmx
    parameter (nlmx=64, npmx=1, nbmx=256)
    double precision :: xf(3,nbas,npmx),xhpot0(nbas,npmx), xgpot0(nlmx*nbas,npmx),xugg(npmx)
    integer:: ibini,ibend
    call tcn('ugcomp')
    call stdfac(20,df)
    pi = 4d0*datan(1d0)
    fpi = 4d0*pi
    srfpi = dsqrt(fpi)
    y0 = 1d0/srfpi
    mp = 1
    if (npmx < mp) call rxi('ugcomp: increase npmx, needed',mp)
    ! --- Loop over sites where charge lump making pot is centered ---
    ugg = 0d0
    iv0 = 0
    ip = 1
    xugg=0d0   !call dpzero(xugg, mp)
    xgpot0=0d0 !call dpzero(xgpot0, nlmx*nbas*mp)
    xf=0d0     !, 3*nbas*mp)
    xhpot0=0d0
    ibini=1
    ibend=nbas
    do ib=ibini,ibend
       is=ispec(ib) 
       tau1=rv_a_opos(:,ib) 
       lmax1=lmxl_i(is)
       rg1=rg_i(is)
       call corprm(is,qcorg1,qcorh1,qsc1,cofg1,cofh1,ceh1,lfoc1, rh1,z1)
       nlm1 = (lmax1+1)**2
       if (lmax1 > -1) then
          jv0 = 0
          do  jb = 1, nbas!  Loop over sites where charge lump sees the potential
             js=ispec(jb)
             tau2=rv_a_opos(:,jb) 
             lmax2=lmxl_i(js)
             rg2=rg_i(js)
             if (lmax2 > -1) then
                call corprm(js,qcorg2,qcorh2,qsc2,cofg2,cofh2,ceh2, lfoc2,rh2,z2)
                nlm2 = (lmax2+1)**2
                if (nlm1 > ndim) call rxi('ugcomp: ndim < nlm1=',nlm1)
                if (nlm2 > ndim) call rxi('ugcomp: ndim < nlm2=',nlm2)
                call ggugbl(tau1,tau2,rg1,rg2,nlm1,nlm2,ndim,ndim,s,ds)  !gaussian-gaussian
                ff = 0d0
                do  ilm1 = 1, nlm1
                   l1 = ll(ilm1)
                   qm1 = qmom(iv0+ilm1)
                   if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
                   cof1 = qm1*fpi/df(2*l1+1)
                   do  ilm2 = 1, nlm2
                      l2 = ll(ilm2)
                      qm2 = qmom(jv0+ilm2)
                      if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
                      cof2 = qm2*fpi/df(2*l2+1)
                      xugg(ip) = xugg(ip) + cof1*cof2*s(ilm1,ilm2)
                      xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip) + s(ilm1,ilm2)*cof1*fpi/df(2*l2+1)
                      ff = ff + 0.5d0*cof1*cof2*ds(ilm1,ilm2,1:3) !Forces
                   enddo
                enddo
                if (lfoc1 > 0 .OR. lfoc2 > 0) then ! Additional h*h, h*g, g*h terms for foca ---
                   call hhugbl(0,tau1,tau2,[rh1],[rh2],[ceh1],[ceh2],1,1,ndim0,ndim0, s0,ds0)!smHankel-smHankel
                   xugg(ip) = xugg(ip) + cofh1*s0(1,1)*cofh2
                   xhpot0(jb,ip) = xhpot0(jb,ip) + cofh1*s0(1,1)
                   ff = ff + 0.5d0*cofh1*cofh2*ds0(1,1,1:3)
                   call hgugbl(tau1,tau2,rh1,rg2,ceh1,1,nlm2,ndim,ndim, s,ds)!gaussian-smHankel
                   do  ilm2 = 1, nlm2
                      l2 = ll(ilm2)
                      qm2 = qmom(jv0+ilm2)
                      if (ilm2 == 1) qm2 = qm2 + y0*(qcorg2-z2)
                      cof2 = qm2*fpi/df(2*l2+1)
                      xugg(ip) = xugg(ip) + cofh1*s(1,ilm2)*cof2
                      ff = ff + 0.5d0*cofh1*cof2*ds(1,ilm2,1:3)
                      xgpot0(jv0+ilm2,ip) = xgpot0(jv0+ilm2,ip) + s(1,ilm2)*cofh1*fpi/df(2*l2+1)
                   enddo
                   call hgugbl(tau2,tau1,rh2,rg1,ceh2,1,nlm1,ndim,ndim, s,ds) !gaussian-smHamel
                   do  ilm1 = 1, nlm1
                      l1 = ll(ilm1)
                      qm1 = qmom(iv0+ilm1)
                      if (ilm1 == 1) qm1 = qm1 + y0*(qcorg1-z1)
                      cof1 = qm1*fpi/df(2*l1+1)
                      xugg(ip) = xugg(ip) + cof1*s(1,ilm1)*cofh2
                      ff = ff - 0.5d0*cof1*cofh2*ds(1,ilm1,1:3)
                      xhpot0(jb,ip) = xhpot0(jb,ip) + cof1*s(1,ilm1)
                   enddo
                endif
                if (jb /= ib) then
                   xf(1:3,ib,ip) = xf(1:3,ib,ip) - ff
                   xf(1:3,jb,ip) = xf(1:3,jb,ip) + ff
                endif
                jv0 = jv0+nlm2
             endif
          enddo
          iv0 = iv0+nlm1
       endif
    enddo
    nvl = iv0
    do 80 ip = 1, mp
       do ib = 1, nbas
          f(1:3,ib) = f(1:3,ib) + xf(1:3,ib,ip)
          hpot0(ib) = hpot0(ib) + xhpot0(ib,ip)
       enddo
       gpot0(1:nvl) = gpot0(1:nvl) + xgpot0(1:nvl,ip)
       ugg = ugg + xugg(ip)
80  enddo
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
subroutine symqmp(nrclas,nlml,nlmx,plat,posc,ngrp,g,ag,qwk,ipa,sym,qmp,nn)  !- Symmetrize multipole moments for a single class
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
  !u Updates
  !u   23 Aug 01 adapted from psymql
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nrclas,nlmx,ipa(nrclas),nlml,ngrp,nn
  double precision :: plat(3,3),posc(3,nrclas),g(3,3,ngrp),ag(3,ngrp)
  double precision :: sym(nlmx,nlmx,nrclas),qwk(nlml),qmp(*)
  integer :: ia,ilm,ixx,iyy(1)
  double precision :: wgt,xx,qlat(3,3)
  call tcn('symqmp')
  if (nlml > nlmx) call rxi('symqmp: increase nlmx to',nlml)
  call dinv33(plat,1,qlat,xx)
  ! ... Make the symmetry projectors
  call symprj(nrclas,nlmx,ngrp,ixx,iyy,g,ag,plat,qlat,posc,sym)
  ! ... Accumulate symmetrized qmpol on first site
  call dpzero(qwk, nlml)
  do  ia = 1, nrclas
     call pxsmr1(1d0,1,nlml,1,sym(1,1,ia),qmp(1+ipa(ia)),qwk,nn)
  enddo
  ! ... Rotate and copy to all sites in class
  wgt = nrclas
  do  ia = 1, nrclas
     call dpzero(qmp(1+ipa(ia)), nlml)
     call pysmr1(wgt,1,nlml,1,sym(1,1,ia),qwk,qmp(1+ipa(ia)),nn)
  enddo
  nn = 0
  do  ilm = 1, nlml
     if (dabs(qmp(ilm+ipa(1))) > 1d-6) nn = nn+1
  enddo
  call tcx('symqmp')
end subroutine symqmp
end module m_smves
