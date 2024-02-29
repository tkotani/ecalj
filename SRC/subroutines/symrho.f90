module m_symrhoat
  use m_lmfinit,only: nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxl_i=>lmxl,spec_a
  public symrhoat
  private
  contains
subroutine symrhoat(sv_p_orhoat, qbyl, hbyl, f)!- Symmetrize charge density and related quantities
  use m_supot,only: rv_a_ogv , iv_a_oips0 , zv_a_obgv ,iv_a_okv
  use m_struc_def
  use m_lmfinit,only: nbas,nsp,lfrce
  use m_lgunit,only:stdo
  !i   lfrce    :>0 symmetrize forces
  !i Inputs/Outputs
  ! o  orhoat:vector of offsets containing site density
  ! o        :Symmetrized on output
  ! o  qbyl  :site- and l-decomposed charges
  ! o        :Symmetrized on output
  ! o  hbyl  :site- and l-decomposed one-electron energies
  ! o        :Symmetrized on output
  ! o  f     :forces
  ! o        :Symmetrized on output
  !r Remarks
  !u Updates
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   19 Jun 00 Packaged from nfp symrat.f and symsmr.f
  implicit none
  type(s_rv1) :: sv_p_orhoat(3,*)
  real(8):: f(*) , qbyl(*) , hbyl(*)
  integer :: iprint
  call tcn('symrhoat')
  if(iprint()>=10) write(stdo,*)' Symmetrize density..'
  call symrat (sv_p_orhoat , qbyl , hbyl , f )
  if ( iprint ( ) > 60 ) call prrhat (sv_p_orhoat )
  call tcx('symrhoat')
end subroutine symrhoat
subroutine symrat(sv_p_orhoat , qbyl , hbyl , f ) !Symmetrize the atomic charge densities and the forces.
  use m_struc_def
  use m_mksym,only: oistab , symops, ag,ngrp,ipc_iv=>iclasst
  use m_lattic,only: lat_qlat,lat_plat,rv_a_opos
  use m_lgunit,only:stdo
  use m_lmfinit,only: ispec,lfrce,nbas,nsp,n0
  !i   nbas  :size of basis
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   lfrce    :>0 symmetrize forces
  !i Inputs/Outputs
  ! o  orhoat:vector of offsets containing site density
  ! o        :Symmetrized on output
  ! o  qbyl  :site- and l-decomposed charges
  ! o        :Symmetrized on output
  ! o  hbyl  :site- and l-decomposed one-electron energies
  ! o        :Symmetrized on output
  ! o  f     :forces
  ! o        :Symmetrized on output
  implicit none
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  real(8):: f(3,nbas),qbyl(n0,nsp,nbas),hbyl(n0,nsp,nbas)
  integer:: ib0,ic,ipr,iprint,is,lmxa,lmxl,nclass,nlml,nlmx,nr,nrclas,igetss,ival,ibas
  integer ,allocatable :: ipa_iv(:)
  integer ,allocatable :: ips_iv(:)
  real(8) ,allocatable :: pos_rv(:)
  real(8) ,allocatable :: pos0_rv(:,:)
  real(8) ,allocatable :: rho_rv(:)
  real(8) ,allocatable :: sym_rv(:)
  double precision :: plat(3,3),qlat(3,3)
  call tcn('symrat')
  ipr = iprint()
  plat=lat_plat
  qlat=lat_qlat
!  ngrp=lat_nsgrp
  allocate(ips_iv(nbas))
  allocate(pos0_rv(3,nbas))
  do ibas=1,nbas
     pos0_rv(:,ibas)= rv_a_opos(:,ibas) 
  enddo
  nclass = maxval(ipc_iv )
  allocate(ipa_iv(nbas))
  allocate(pos_rv(3*nbas))
  if (iprint() >= 40) write(stdo,"(a)")'                qbyl        hbyl        ebar   '
  do  ic = 1, nclass
     call psymr0 ( - 2,ic,nbas,ipc_iv,pos0_rv,pos_rv, ipa_iv,nrclas )
     if (nrclas <=0) cycle
     ib0 = ipa_iv(1) !ival ( ipa_iv,1 )
     is = ispec(ib0) 
     lmxl=lmxl_i(is)
     lmxa=lmxa_i(is)
     nr=nr_i(is)
     nlml = (lmxl+1)**2
     !   ... Make the projectors; make to at least to l=1 for forces
     nlmx = max0(nlml,4)
     allocate(sym_rv(nlmx*nlmx*nrclas))
     call symprj(nrclas,nlmx,ngrp,nbas,oistab,symops,ag,plat,qlat,pos_rv,sym_rv )
     !   ... Apply the projectors to rhoat
     if (lmxl > -1) then
        allocate(rho_rv(nr*nlml*nsp))
        call psymr1 ( nrclas,ipa_iv,nr,nlml,nsp,nlmx,sym_rv, rho_rv,sv_p_orhoat,1 )
        call psymr1 ( nrclas,ipa_iv,nr,nlml,nsp,nlmx,sym_rv, rho_rv,sv_p_orhoat,2 )
        !   ... Symmetrize site charges and eigval sum
        call psymrq ( nrclas,nsp,ipa_iv,lmxa,qbyl,hbyl ) !write qbyl
        deallocate(rho_rv)
     endif
     !   ... Symmetrize the forces
     if ( lfrce /= 0 ) call psymrf ( nrclas,ipa_iv,nlmx,sym_rv, f )
     deallocate(sym_rv)
  enddo
  call tcx('symrat')
end subroutine symrat
subroutine psymrf(nrclas,ipa,nlmx,s,f)
  implicit none
  integer :: nrclas,nlmx,ipa(nrclas)
  double precision :: s(nlmx,nlmx,nrclas),f(3,*)
  integer :: ia,ib
  double precision :: x(3)
  x = 0d0
  do  ia = 1, nrclas
     ib = ipa(ia)
     x(1)= x(1)+s(4,4,ia)*f(1,ib)+s(4,2,ia)*f(2,ib)+s(4,3,ia)*f(3,ib)
     x(2)= x(2)+s(2,4,ia)*f(1,ib)+s(2,2,ia)*f(2,ib)+s(2,3,ia)*f(3,ib)
     x(3)= x(3)+s(3,4,ia)*f(1,ib)+s(3,2,ia)*f(2,ib)+s(3,3,ia)*f(3,ib)
  enddo
  do  ia = 1, nrclas
     ib = ipa(ia)
     f(1,ib) = (s(4,4,ia)*x(1)+s(2,4,ia)*x(2)+s(3,4,ia)*x(3))*nrclas
     f(2,ib) = (s(4,2,ia)*x(1)+s(2,2,ia)*x(2)+s(3,2,ia)*x(3))*nrclas
     f(3,ib) = (s(4,3,ia)*x(1)+s(2,3,ia)*x(2)+s(3,3,ia)*x(3))*nrclas
  enddo
end subroutine psymrf
subroutine psymrq(nrclas,nsp,ipa,lmxa,qbyl,hbyl)!- Symmetrize l-decomposed site charges and eval sums
  use m_lgunit,only:stdo
  implicit none
  integer :: nrclas,nsp,lmxa,ipa(nrclas),n0
  parameter (n0=10)
  double precision :: qbyl(n0,nsp,*),hbyl(n0,nsp,*)
  integer :: ia,ib,iprint,l,isp,icopy_size
  double precision :: qsum(n0,2),hsum(n0,2),fac
  call dpzero(qsum,2*n0)
  call dpzero(hsum,2*n0)
  fac = 1d0/nrclas
  do  ia = 1, nrclas
     ib = ipa(ia)
     do  isp = 1, nsp
        do  l = 0, lmxa
           qsum(l+1,isp) = qsum(l+1,isp) + qbyl(l+1,isp,ib)*fac
           hsum(l+1,isp) = hsum(l+1,isp) + hbyl(l+1,isp,ib)*fac
        enddo
     enddo
  enddo
  if (iprint() >= 40) then
     !         write(stdo,"(a)")'                qbyl        hbyl        ebar   '
     !        write(stdo,"(a)")'   ebar=hbyl/qbyl= center of gravity of occpied states'
     !        write(stdo,770) (ipa(ia),ia = 1,nrclas)
     !  770   format(' symmetrized qbyl,hbyl for class containing atom site ib=',20i3)
     !  770   format(' atom site ib=',255i3)
     if (nsp == 1) write(stdo,780)(ipa(1),l,qsum(l+1,1),hsum(l+1,1), hsum(l+1,1)/qsum(l+1,1), l=0,lmxa)
     if (nsp == 2) write(stdo,781)(ipa(1),l,(qsum(l+1,isp),hsum(l+1,isp),hsum(l+1,isp)/qsum(l+1,isp),isp=1,nsp), l=0,lmxa)
780  format('ib=',i3,' l=',i2,3f12.6)
781  format('ib=',i3,' l=',i2,3f12.6,'   spin 2',3f12.6)
     write(stdo,*)'-----------------------------------------------'
  endif
  do  ia = 1, nrclas
     ib = ipa(ia)
     do  isp = 1, nsp
        do  l = 0, lmxa
           qbyl(l+1,isp,ib) = qsum(l+1,isp)
           hbyl(l+1,isp,ib) = hsum(l+1,isp)
        enddo
     enddo
  enddo
end subroutine psymrq
subroutine psymr1 ( nrclas,ipa,nr,nlml,nsp,nlmx,sym,rho,sv_p_orhoat,icmp )
  use m_struc_def, only: s_rv1
  use m_lgunit,only:stdo
  !- Symmetrize density for one class of atoms
  implicit none
  integer :: nrclas,nsp
  integer:: ipa(nrclas),nlmx,nr,nlml,icmp
  type(s_rv1) :: sv_p_orhoat(3,*)
  double precision :: sym(nlmx,nlmx,nrclas),rho(nr,nlml,nsp)
  integer :: ia,ib,iprint,nn
  double precision :: wgt
  call dpzero(rho, nr*nlml*nsp)
  do  ia = 1, nrclas
     ib = ipa(ia)
     call pxsmr1 ( 1d0,nr,nlml,nsp,sym ( 1,1,ia ),sv_p_orhoat( icmp,ib )%v, rho,nn )
  enddo
  wgt = nrclas
  do  ia = 1, nrclas
     ib = ipa(ia)
     call dpzero ( sv_p_orhoat( icmp,ib )%v,nr * nlml * nsp )
     call pysmr1 ( wgt,nr,nlml,nsp,sym ( 1,1,ia ),rho,sv_p_orhoat( icmp,ib )%v,nn )
  enddo
end subroutine psymr1
subroutine prrhat(sv_p_orhoat )
  use m_struc_def
  use m_lmfinit,only: nsp,nbas, ispec
  use m_lgunit,only:stdo
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  integer:: ib , lmxl , nr , nlml , is , igetss
  real(8) ,allocatable :: rofi_rv(:)
  real(8) ,allocatable :: rwgt_rv(:)
  double precision :: a,rmt
  do  ib = 1, nbas
     is = ispec(ib) 
     lmxl=lmxl_i(is)
     a=spec_a(is)
     nr=nr_i(is)
     rmt=rmt_i(is)
     if (lmxl == -1) goto 10
     allocate(rofi_rv(nr))
     allocate(rwgt_rv(nr))
     call radmsh ( rmt , a , nr , rofi_rv )
     call radwgt ( rmt , a , nr , rwgt_rv )
     nlml = (lmxl+1)**2
     write(stdo,200) ib,rmt,nr,nlml
200  format(/' Density at site',i3,'   rmt=',f8.4,'   nr=',i5,'   nlml=',i3)
     call prlrho('true density',  nr,nlml,nsp,rofi_rv,rwgt_rv,sv_p_orhoat(1,ib )%v )
     call prlrho('smooth density',nr,nlml,nsp,rofi_rv,rwgt_rv,sv_p_orhoat(2,ib )%v )
     if (allocated(rwgt_rv)) deallocate(rwgt_rv)
     if (allocated(rofi_rv)) deallocate(rofi_rv)
10   continue
  enddo
end subroutine prrhat
subroutine prlrho(str,nr,nlm,nsp,rofi,rwgt,rho)!- Print info about site density
  use m_ll,only:ll
  use m_lgunit,only:stdo
  implicit none
  integer :: nr,nlm,nsp
  double precision :: rho(nr,nlm,nsp),rofi(nr),rwgt(nr)
  character*(*) str
  integer :: ilm,l,itop,ibot,i
  double precision :: xx,pi,srfpi,top,bot,sum
  pi = 4d0*datan(1d0)
  srfpi = dsqrt(4d0*pi)
  write(stdo,601) str
  do  ilm = 1, nlm
     l = ll(ilm)
     top = -1d10
     bot = 1d10
     sum = 0d0
     itop = 0
     ibot = 0
     do  i = 1, nr
        xx = (rho(i,ilm,1)+rho(i,ilm,nsp))/(3-nsp)
        if (xx > top) then
           itop = i
           top = xx
        endif
        if (xx < bot) then
           ibot = i
           bot = xx
        endif
        sum = sum + rwgt(i)*xx*(rofi(i)+1d-32)**l
     enddo
     xx = dmax1(dabs(top),dabs(bot))
     if (xx > 1d-6) then
        if (ilm == 1) then
           write(stdo,600) ilm,bot,ibot,top,itop,sum,srfpi*sum
        else
           write(stdo,600) ilm,bot,ibot,top,itop,sum
        endif
     endif
600  format(i4,f15.8,i5,f15.8,i5,f15.8,f12.6)
601  format(/' prlrho: ',a/'  ilm',6x,'rhomin    pnt',7x,'rhomax    pnt       moment')
  enddo
end subroutine prlrho
end module m_symrhoat

subroutine symsmrho(smrho) !- Symmetrize the smooth charge density
  use m_supot,only:  iv_a_oips0, iv_a_okv,rv_a_ogv,zv_a_obgv
  use m_supot,only: lat_ng,n1,n2,n3
  use m_mksym,only: ngrp
  use m_lmfinit,only: nsp
  use m_lgunit,only:stdo
  implicit none
  double complex smrho(n1,n2,n3,nsp)
  integer:: ng,isp
  complex(8) ,allocatable :: csym_zv(:), cv_zv(:)
  call tcn('symsmrho')
  if (ngrp > 1) then
     ng=lat_ng
     allocate(cv_zv(ng),csym_zv(ng))
     call fftz3(smrho,n1,n2,n3,n1,n2,n3,nsp,0,-1)
     do  isp = 1, nsp
        call gvgetf( ng,1,iv_a_okv,n1,n2,n3,smrho(1,1,1,isp), cv_zv )
        call gvsym ( ng,rv_a_ogv,iv_a_oips0,zv_a_obgv,cv_zv,   csym_zv )
        csym_zv=csym_zv - cv_zv
        call gvaddf( ng,iv_a_okv,n1,n2,n3,csym_zv, smrho(1,1,1,isp) )
     enddo
     call fftz3(smrho,n1,n2,n3,n1,n2,n3,nsp,0,1)
     deallocate(csym_zv,cv_zv)
  else
     write(stdo,*)' Smooth density not symmetrized (ngrp=1)'
  endif
  call tcx('symsmrho')
end subroutine symsmrho
subroutine gvsym(ng,gv,ips0,bgv,c,csym)  !- Symmetrize a function c, given in the form of a list
  !i   gv,ng   Lattice vectors, and number
  !i   ips0    pointer to first vector in star of this vector; see sgvsym
  !i   bgv     phase factor sum; see sgvsym
  !i   c       unsymmetrized function
  !o Outputs:
  !o   csym    symmetrized function
  implicit none
  integer :: ng,ips0(ng)
  double precision :: gv(ng,3)
  double complex bgv(ng),c(ng),csym(ng)
  integer :: i,j,i0,kstar,ipr,iprint,nstar   ! ... Sum up coefficients for first vector in each star
  csym=0d0
  do i = 1, ng
     csym(ips0(i)) = csym(ips0(i)) + bgv(i)*c(i)
  enddo
  do i0 = 1, ng
     if (ips0(i0) == i0) then
        kstar = 0
        do i = i0, ng
           if (ips0(i) == i0) kstar = kstar+1
        enddo
        csym(i0) = csym(i0)/kstar
     endif
  enddo
  do i = 1, ng
     csym(i) = csym(ips0(i))* dconjg(bgv(i))
  enddo
end subroutine gvsym
