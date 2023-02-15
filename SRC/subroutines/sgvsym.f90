subroutine sgvsym(ngrp,g,ag,ng,gv,ips0,bgv)
  use m_lgunit,only:stdo
  use m_ftox
  !- Setup for symmetrization of a function in Fourier representation.
  ! ----------------------------------------------------------------------
  !i Inputs:
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
  !b Bugs
  !b   This routine copies su_gvsym.f conventions, which use inv(g)
  !u Updates
  !u    7 Sep 98 Adapted from nfp su_gvsym.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ngrp,ng,ips0(ng)
  double precision :: g(3,3,1),ag(3,1),gv(ng,3)
  complex(8):: bgv(ng)
  integer :: i,i00,irep,i0,nstar,k,j,j0,iprint,ksum,kstar
  real(8)::df,scalp,gg0,gg,fac,vv,v(3),diffmin
  integer:: jx,jg,jjg
  real(8),parameter:: tpi = 8d0*datan(1d0),tol=1d-3,tol3=1d-3
  ips0 = 0
  bgv = 0d0
  ! --- Main loop: look for next unclassified vector ---
  i00 = 1
  do  20  irep = 1, ng+1
     i0 = 0
     do  22  i = i00, ng
        i0 = i
        if (ips0(i) == 0) goto 80
22   enddo
     goto 81
80   continue
     !   --- Apply all point ops, find in list, add to phase sum ---
     nstar = irep
     do  30  k = 1, ngrp  ! ... Find G' = g(k) G; j0 is index to G'
        v = matmul(g(:,:,k),gv(i0,:))
        do  32  j = i0, ng
           j0 = j
           if (sum((v-gv(j,:))**2) < tol) goto 70
32      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        write(6,ftox) i0,'absv=',ftof(sum(v**2)),' v=',ftof(v)
        do  jx = 1, ng
           write(6,ftox) jx,' absg=',ftof(sum(gv(jx,:)**2)),'gv',ftof(gv(jx,:))
        enddo
        write(stdo,"('--- igvec',i6,' igrp',i4,' v=',3f9.5,' gv**2',f12.8)")i0,k,v,sum(v**2)
         diffmin=9999
         do jg=1,ng
            if(diffmin>sum((v-gv(jg,:))**2)) then
               jjg=jg
               diffmin = sum((v-gv(jg,:))**2)
            endif   
         enddo
         write(stdo,"(' nearest :',i5,13x,3f9.5,' gv**2 ',f12.8,' diff=',d8.2)") &
              jjg, gv(jjg,:), sum(gv(jjg,:)**2),diffmin
         call rxi('SGVSYM: cannot find mapped vector in list:',i0)
!!!!!!!!!!!!!!!!!!!!!!!!!         
70      continue
        ips0(j0) = i0
        scalp = gv(j0,1)*ag(1,k)+gv(j0,2)*ag(2,k)+gv(j0,3)*ag(3,k)
        bgv(j0) = bgv(j0) + cdexp(dcmplx(0d0,tpi*scalp))
30   enddo
     i00 = i0
20 enddo
  call rxi('SGVSYM: this cannot happen. irep=',irep)
81 continue

  !      if (iprint() .ge. 20) call awrit2(' SGVSYM: %i symmetry stars'//
  !     .' found for %i reciprocal lattice vectors',
  !     .' ',80,stdo,nstar,ng)
  if (iprint() >= 20) write(stdo,"(' SGVSYM: ',i0,' symmetry stars found for ', &
       i0,' reciprocal lattice vectors')") nstar,ng

  ! --- Multiply phase sums by (star order)/(group order) ---
  ksum = 0
  do  40  i0 = 1, ng
     if (ips0(i0) == i0) then
        kstar = 0
        gg0 = gv(i0,1)**2+gv(i0,2)**2+gv(i0,3)**2
        do  42  i = i0, ng
           if (ips0(i) == i0) kstar = kstar+1
           gg = gv(i,1)**2+gv(i,2)**2+gv(i,3)**2
           if (dabs(gg-gg0) > tol) goto 78
42      enddo
78      continue
        ksum = ksum+kstar
        fac = dble(kstar)/dble(ngrp)
        do  44  i = i0, ng
           if (ips0(i) == i0) bgv(i) = bgv(i)*fac
           gg = gv(i,1)**2+gv(i,2)**2+gv(i,3)**2
           if (dabs(gg-gg0) > tol) goto 79
44      enddo
79      continue
     endif
40 enddo
  if (ksum /= ng) call rxi('SGVSYM error, ksum=',ksum)
  ! --- Printout ---
  if (iprint() < 60) return
  j = min(ng,300)
  if (iprint() >= 102) j = ng
  nstar = 0
  write(stdo,252)
252 format(//'   elt     no           vector                phase')
  do  50  i0 = 1, j
     if (ips0(i0) == i0) then
        nstar = nstar+1
        vv = dsqrt(gv(i0,1)**2+gv(i0,2)**2+gv(i0,3)**2)
        write(stdo,250) nstar,vv
250     format(/'   Star',i6,'  length =',f8.4)
        kstar = 0
        do  56  i = i0, ng
           if (ips0(i) == i0) then
              kstar = kstar+1
              write(stdo,251) kstar,i,gv(i,1),gv(i,2),gv(i,3),bgv(i)
251           format(i5,i8,2x,3f7.2,3x,2f7.2)
           endif
56      enddo
     endif
50 enddo
end subroutine sgvsym

