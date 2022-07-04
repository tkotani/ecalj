subroutine mshvmt(nbas,ssite,sspec,ng,gv, kv,cv,k1,k2,k3,smpot,vval)
  use m_struc_def
  use m_lattic,only:lat_plat,rv_a_opos
  use m_lmfinit,only:lat_alat
  use m_supot,only: lat_nabc
  use m_ropyln,only: ropyln
  use m_ropbes,only: ropbes
!- Makes potential at MT surfaces given potential on a uniform mesh
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    class spec pos
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: rmt lmxl
  !i     Stored:    *
  !i     Passed to: *
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat plat nabc nsgrp osymgr oag
  !i     Stored:    *
  !i     Passed to: *
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   cv    :work array
  !i   k1,k2,k3 dimensions of smrho,smpot for smooth mesh density
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
  integer :: k1,k2,k3,nbas,ng,kv(ng,3)
  real(8):: gv(ng,3) , vval(1)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex smpot(k1,k2,k3),cv(ng)
  integer :: i,ib,is,lmxx,nlmx,iv0,lmxl,nlm,ngabc(3), n1,n2,n3,m,ilm,l,ipr
  double precision :: alat,pi,tpiba,tau(3),rmt,fac,plat(3,3)
  double complex vvali,fprli
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
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
  ngabc=lat_nabc
  tpiba = 2*pi/alat
  call gvgetf(ng,1,kv,k1,k2,k3,smpot,cv)
  ! --- YL(G)*G**l, agv=|G| for each g ---
  call dpcopy(gv,gv2,1,3*ng,tpiba)
  call ropyln(ng,gv2(1,1),gv2(1,2),gv2(1,3),lmxx,ng,yl,agv)
  do  i = 1, ng
     agv(i) = sqrt(agv(i))
  enddo
  iv0 = 0
  do  ib = 1, nbas
     is=ssite(ib)%spec
     tau=rv_a_opos(:,ib) !ssite(ib)%pos
     rmt=sspec(is)%rmt
     lmxl=sspec(is)%lmxl
     if (lmxl == -1) goto 10
     nlm = (lmxl+1)**2
     if (nlm > nlmx) call rxi('mshvmt: increase nlmx to',nlm)
     !       Add a negligibly small amount to rmt to handle case rmt=0
     rmt = rmt+1d-32
     !   --- j_l(|rmt*q|)/rmt**l for each G and l=0..lmax ---
     !       Does not evolve correctly in the correct large r limit
     call ropbes(agv,rmt**2,lmxl,cgp,sgp,phil,ng)
     ! Patch for now
     !       do  i = 1, ng
     !         call besslr(agv(i)**2*rmt**2,0,0,phil(i,0),yl)
     !       enddo

     !   ... Phases exp(-i.G.tau), fast version
     !       call suphs0(plat,ng,gv,gv2)
     !       call dinv33(plat,1,qlat,fac)
     !       call dpzero(q,3)
     !       call suphas(q,tau,ng,gv2,n1,n2,n3,qlat,cgp,sgp)
     !   ... Phases calculated straightforwardly.  Fast enough not to matter.
     call dscal(3,alat,tau,1)
     do  i = 1, ng
        fac = -(tau(1)*gv2(i,1)+tau(2)*gv2(i,2)+tau(3)*gv2(i,3))
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
              vvali = vvali + (phil(i,l)*yl(i,ilm))* &
                   (cv(i)*dcmplx(cgp(i),-sgp(i)))
           enddo
           vval(ilm+iv0) = fprli*vvali
        enddo
        fprli = fprli*(0d0,1d0)*rmt
     enddo
     !   ... Printout
     !        if (ipr .gt. 0) then
     !          do  ilm = 1, nlm
     !            if (ilm .eq. 1) then
     !              write(stdo,650) ib,ilm,vval(ilm+iv0)
     !            elseif (dabs(vval(ilm+iv0)) .gt. 1d-6) then
     !              write(stdo,651)    ilm,vval(ilm+iv0)
     !            endif
     !  650              format(i4,i6,2f12.6)
     !  651                     format(4x,i6,f12.6)
     !          enddo
     !        endif
     iv0 = iv0 + nlm
10   continue
  enddo
  deallocate(phil,yl)
  deallocate(gv2,agv,cgp,sgp)
  call tcx('mshvmt')
end subroutine mshvmt

subroutine symvvl(nbas,ssite,sspec,vval,vrmt)
  use m_mksym,only: rv_a_osymgr,rv_a_oag,lat_nsgrp,ipc=>iclasst
  use m_struc_def
  use m_lattic,only:lat_plat,rv_a_opos
  use m_lgunit,only:stdo
  !- Symmetrizes the potential at the MT boundary.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl rmt
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: plat nsgrp osymgr oag alat nabc
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
  integer :: nbas
  real(8):: vval(1) , vrmt(nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  integer :: ic,ib,ilm,mxint,nclass,ipa(nbas),nrclas,iv0
  integer :: ips(nbas),lmxl(nbas) !ipc(nbas),
  double precision :: pos(3,nbas),posc(3,nbas),plat(3,3),pi,y0
  integer:: igetss , nlml ,ipr , jpr , ngrp , nn , iclbas,ibas
  real(8) ,allocatable :: qwk_rv(:)
  real(8) ,allocatable :: sym_rv(:)
  call tcn('symvvl')
  call getpr(ipr)
  plat=lat_plat
  ngrp=lat_nsgrp
  do ibas=1,nbas
!     ipc(ibas)  = ssite(ibas)%class
     ips(ibas)  = ssite(ibas)%spec
     pos(:,ibas)= rv_a_opos(:,ibas) !ssite(i_spackv)%pos
  enddo
  nclass = mxint(nbas,ipc)
  do  ib = 1, nbas
     lmxl ( ib ) = int(sspec(ips(ib))%lmxl)
  enddo
  do  ic = 1, nclass
     !   ... Make nrclas,ipa,posc
     call psymr0(lmxl,ic,nbas,ipc,pos,posc,ipa,nrclas)
     if (nrclas > 0) then
        ib = iclbas(ic,ipc)
        if (lmxl(ib) > -1) then
           nlml = (lmxl(ib)+1)**2
           if (ipr >= 50) write(stdo,800) ic,nrclas,nlml
800        format(' Symmetry class',i3,'   nrclas=',i3,'   nlml=',i3)
           allocate(qwk_rv(nlml))
           allocate(sym_rv(nlml*nlml*nrclas))
           call symqmp ( nrclas , nlml , nlml , plat , posc , ngrp , rv_a_osymgr &
                , rv_a_oag , qwk_rv , ipa , sym_rv , vval , nn )
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
     if (ib == iclbas(ic,ipc)) then
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
650        format(i4,2i6,2f12.6)
651        format(10x,i6,f12.6)
        enddo
     endif
     iv0 = iv0 + nlml
10   continue
  enddo
  call tcx('symvvl')
end subroutine symvvl


