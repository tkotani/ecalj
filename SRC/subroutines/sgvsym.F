      subroutine sgvsym(ngrp,g,ag,ng,gv,ips0,bgv)
C- Setup for symmetrization of a function in Fourier representation.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ngrp,g,ag   space group
Ci   ng,gv       list of reciprocal lattice vectors
Co Outputs:
Co   ips0        pointer to first vector in star of this vector
Co   bgv         phase factor sum; see Remarks
Cr Remarks:
Cr   The reciprocal lattice vectors are assumed to be sorted by length
Cr   Adapted from nfp su_gvsym.f
Cr
Cr   Symmetrized f means  f(G) = f(g(G)).
Cr   Let G be all r.l.v and G' be just first members of each star.  Then
Cr
Cr   f(r) = sum_G f(G) exp(i G.r)
Cr        = sum_G' sum_(G in star of G')  f(G') exp(i G.r)
Cr        = sum_G' f(G') exp(i G.r) sum_(G in star of G') exp(i(G-G').r)
Cr        = sum_G' f(G') exp(i G.r) bgv(star)
Cb Bugs
Cb   This routine copies su_gvsym.f conventions, which use inv(g)
Cu Updates
Cu    7 Sep 98 Adapted from nfp su_gvsym.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ngrp,ng,ips0(ng)
      double precision g(3,3,1),ag(3,1),gv(ng,3)
      double complex bgv(ng)
C ... Local parameters
      integer i,i00,irep,i0,nstar,k,j,j0,iprint,ksum,kstar,stdo,lgunit
      double precision tpi,df,scalp,gg0,gg,fac,vv,v(3)

      stdo = lgunit(1)
      tpi = 8d0*datan(1d0)

      do  10  i = 1, ng
        ips0(i) = 0
        bgv(i) = (0d0,0d0)
   10 continue

C --- Main loop: look for next unclassified vector ---
      i00 = 1
      do  20  irep = 1, ng+1
        i0 = 0
        do  22  i = i00, ng
          i0 = i
          if (ips0(i) .eq. 0) goto 80
   22   continue
        goto 81
  80    continue

C   --- Apply all point ops, find in list, add to phase sum ---
        nstar = irep
        do  30  k = 1, ngrp
C     ... Find G' = g(k) G; j0 is index to G'
          v(1) = g(1,1,k)*gv(i0,1)+g(1,2,k)*gv(i0,2)+g(1,3,k)*gv(i0,3)
          v(2) = g(2,1,k)*gv(i0,1)+g(2,2,k)*gv(i0,2)+g(2,3,k)*gv(i0,3)
          v(3) = g(3,1,k)*gv(i0,1)+g(3,2,k)*gv(i0,2)+g(3,3,k)*gv(i0,3)
          do  32  j = i0, ng
            df = (v(1)-gv(j,1))**2+(v(2)-gv(j,2))**2+(v(3)-gv(j,3))**2
            j0 = j
            if (df .lt. 1d-8) goto 70
   32     continue
          write(stdo,601) i0,k,gv(i0,1),gv(i0,2),gv(i0,3),v
  601     format(' ---- vec',i6,'   op',i3,2x,3f8.4,'  ->',3f8.4)
          call rxi('SGVSYM: cannot find mapped vector in list:',i0)
  70      continue
          ips0(j0) = i0
          scalp = gv(j0,1)*ag(1,k)+gv(j0,2)*ag(2,k)+gv(j0,3)*ag(3,k)
          bgv(j0) = bgv(j0) + cdexp(dcmplx(0d0,tpi*scalp))
   30   continue
        i00 = i0
   20 continue
      call rxi('SGVSYM: this cannot happen. irep=',irep)
  81  continue

      if (iprint() .ge. 20) call awrit2(' SGVSYM: %i symmetry stars'//
     .' found for %i reciprocal lattice vectors',
     .' ',80,stdo,nstar,ng)

C --- Multiply phase sums by (star order)/(group order) ---
      ksum = 0
      do  40  i0 = 1, ng
        if (ips0(i0) .eq. i0) then
          kstar = 0
          gg0 = gv(i0,1)**2+gv(i0,2)**2+gv(i0,3)**2
          do  42  i = i0, ng
            if (ips0(i) .eq. i0) kstar = kstar+1
            gg = gv(i,1)**2+gv(i,2)**2+gv(i,3)**2
            if (dabs(gg-gg0) .gt. 1d-3) goto 78
   42     continue
   78     continue
          ksum = ksum+kstar
          fac = dble(kstar)/dble(ngrp)
          do  44  i = i0, ng
            if (ips0(i) .eq. i0) bgv(i) = bgv(i)*fac
            gg = gv(i,1)**2+gv(i,2)**2+gv(i,3)**2
            if (dabs(gg-gg0) .gt. 1d-3) goto 79
   44     continue
   79     continue
        endif
   40 continue
      if (ksum .ne. ng) call rxi('SGVSYM error, ksum=',ksum)

C --- Printout ---
      if (iprint() .lt. 60) return
      j = min(ng,300)
      if (iprint() .ge. 102) j = ng
      nstar = 0
      write(stdo,252)
  252 format(//'   elt     no           vector                phase')
      do  50  i0 = 1, j
        if (ips0(i0) .eq. i0) then
          nstar = nstar+1
          vv = dsqrt(gv(i0,1)**2+gv(i0,2)**2+gv(i0,3)**2)
          write(stdo,250) nstar,vv
  250     format(/'   Star',i6,'  length =',f8.4)
          kstar = 0

          do  56  i = i0, ng
            if (ips0(i) .eq. i0) then
              kstar = kstar+1
              write(stdo,251) kstar,i,gv(i,1),gv(i,2),gv(i,3),bgv(i)
  251         format(i5,i8,2x,3f7.2,3x,2f7.2)
            endif
   56     continue
        endif
   50 continue

      end

