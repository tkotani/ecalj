module m_supot
  use m_struc_def,only: s_rv1
  complex(8) ,allocatable,protected ::  zv_a_obgv (:)
  integer,   allocatable,protected ::  iv_a_oips0 (:), iv_a_okv (:,:)
  real(8),   allocatable,protected ::  rv_a_ogv (:,:)
  real(8),protected:: lat_gmax
  integer,protected:: lat_ng           
  integer,protected,pointer:: n1,n2,n3
  integer,protected,target,private::  ngabc(3)
contains
  subroutine m_supot_init()! Initialization for G vectors bgv,ips0,gv,kv !See gvlst2 and sgvsym
    use m_lmfinit,only : lcd4,nsp,alat=>lat_alat,ftmesh,gmax=>lat_gmaxin
    use m_lgunit,only:stdo
    use m_mksym,only:  ngrp,symops,ag
    use m_lattic,only: plat=>lat_plat,rv_a_opos,qlat=>lat_qlat,vol=>lat_vol, awald=>lat_awald,nkd=>lat_nkd, nkq=>lat_nkq
    use m_shortn3,only: mshsiz
    use m_ftox
    implicit none
    integer:: ngmx, ng, iprint,ig,j1,j2,j3,ierr
    real(8):: xx,qpg(3),gg,qlat1(3,3),gmax2,gs(3),tpiba,  tol=1d-8
    call tcn('m_supot_init')
    ngabc=ftmesh
    n1=>ngabc(1)
    n2=>ngabc(2)
    n3=>ngabc(3)
    if(lcd4) then ! --- Setup for FFT charge density, potential representation ---
       call mshsiz(alat,plat,gmax, ngabc,ngmx) !return n1 n2 n3 (=ngabc) satisfying gmax    !write(stdo,ftox)' 000 gmax ngmx=',ftof(gmax),ngmx
       gvblock: block !Make list of lattice vectors within cutoff
         use m_shortn3,only: gvlst2
         real(8):: ogv(ngmx,3)
         integer:: okv(ngmx,3),ixx(1)
         call gvlst2(alat, plat, [0d0,0d0,0d0], n1,n2,n3, 0d0,gmax,[0],8+1000, ngmx, ng, okv, ogv, ixx)  !+1000 for symmetry cheker for sgvsym
         allocate(rv_a_ogv,source=ogv(1:ng,1:3))
         allocate(iv_a_okv,source=okv(1:ng,1:3))     !print *,'ogv(1:ng,1:3)',ogv(1,1:3),okv(1,1:3)
       endblock gvblock
       lat_ng = ng
       lat_gmax = gmax
       allocate(iv_a_oips0(ng),source=0)
       allocate(zv_a_obgv(ng),source=(0d0,0d0))
       call sgvsym(ngrp, symops , ag , ng , rv_a_ogv , iv_a_oips0 , zv_a_obgv, ierr)
       if(ierr/=0) call rxi('m_supot_init: sgvsym error ierr=',ierr)
    endif
    call tcx('m_supot_init')
  end subroutine m_supot_init
  pure subroutine sgvsym(ngrp,g,ag,ng,gv,ips0,bgv,ierr) !- Setup for symmetrization of a function in Fourier representation.
    use m_lgunit,only:stdo
    use m_ftox
    implicit none
    intent(in) ::        ngrp,g,ag,ng,gv
    intent(out)::                        ips0,bgv,ierr
    !i   ngrp,g,ag   space group
    !i   ng,gv       list of reciprocal lattice vectors
    !o Outputs:
    !o   ips0        pointer to first vector in star of this vector
    !o   bgv         phase factor sum; see Remarks
    !r Remarks:
    !r   The reciprocal lattice vectors are assumed to be sorted by length
    !r   Adapted from nfp su_gvsym.f
    !r
    !r   Symmetrized f means  f(G) = f(g(G)).
    !r   Let G be all r.l.v and G' be just first members of each star.  Then
    !r
    !r   f(r) = sum_G f(G) exp(i G.r)
    !r        = sum_G' sum_(G in star of G')  f(G') exp(i G.r)
    !r        = sum_G' f(G') exp(i G.r) sum_(G in star of G') exp(i(G-G').r)
    !r        = sum_G' f(G') exp(i G.r) bgv(star)
    integer :: ngrp,ng,ips0(ng),i,i00,irep,i0,nstar,k,j,j0,iprint,ksum,kstar,jx,jg,jjg,ierr
    real(8)::df,scalp,gg0,gg,fac,vv,v(3),diffmin,g(3,3,1),ag(3,1),gv(ng,3)
    real(8),parameter:: tpi = 8d0*datan(1d0),tol=1d-3,tol3=1d-3
    complex(8):: bgv(ng)
    complex(8),parameter:: img=(0d0,1d0)
    ips0 = 0
    bgv = 0d0
    i00 = 1
    ierr=-1
    do  irep = 1, ng+1
       i0 = findloc([(ips0(i)==0,i=i00,ng)],value=.true.,dim=1)+i00-1
       if(i0==i00-1) goto 81 
       nstar = irep !   --- Apply all point ops, find in list, add to phase sum ---
       do k = 1, ngrp  ! ... Find G' = g(k) G; j0 is index to G'
          v = matmul(g(:,:,k),gv(i0,:))
          j0 = findloc([(sum((v-gv(j,:))**2) < tol,j=i0,ng)],value=.true.,dim=1) + i0-1
          if(j0==i0-1) return
          !       if(j0==i0-1) call rxi('SGVSYM: cannot find mapped vector in list:',i0)
          ips0(j0) = i0
          bgv(j0)  = bgv(j0) + exp(img*tpi*sum(gv(j0,:)*ag(:,k)))
       enddo
       i00 = i0
    enddo
    ierr=-2
    return
    !  call rxi('SGVSYM: this cannot happen. irep=',irep)
81  continue
    !  if(iprint() >=20) write(stdo,"(' SGVSYM: ',i0,' symmetry stars found for ',i0,' reciprocal lattice vectors')") nstar,ng
    ksum = 0
    do 40  i0 = 1, ng ! --- Multiply phase sums by (star order)/(group order) ---
       if (ips0(i0) == i0) then
          kstar = 0
          gg0 = sum(gv(i0,:)**2)
          do i = i0, ng
             if(ips0(i) == i0) kstar = kstar+1
             if(dabs(sum(gv(i,:)**2)-gg0) > tol) exit
          enddo
          ksum = ksum+kstar
          fac = dble(kstar)/dble(ngrp)
          do i = i0, ng
             if(ips0(i) == i0) bgv(i) = bgv(i)*fac
             if(dabs(sum(gv(i,:)**2)-gg0) > tol) exit
          enddo
       endif
40  enddo
    ierr=0
    !  if (ksum /= ng) call rxi('SGVSYM error, ksum=',ksum)
  endsubroutine sgvsym
endmodule m_supot

