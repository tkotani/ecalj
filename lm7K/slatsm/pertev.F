c takao modified this based on delw_2.2.
      subroutine pertev(iopt,ld,dh,ds,wk,de,e0,z)
C ----------------------------------------------------------------------
Ci Inputs
Ci    iopt   0  1st order perturbation for eigenvalues dh only
Ci           1  include -e0 ds in perturbation
Ci          10s digit
Ci           1  perturbation also for eigenvetors
Ci           2  Renormalize eigenvectors
Ci           3  Orthornnormalize eigenvectors
Ci         100s digit
Ci           1  include 2nd order pert. correction for evals
Ci              (warning) not properly checked; may be wrong.
Ci              ... so not implemented
Ci   ld         dimension of hamiltonian
Ci   e0         eigenvalues of unperturbed H
Ci   dh,ds      Perturbing hamiltonian, and -e0 ds
Ci   z          eigenvectors of unperturbed H
Ci   wk         work array of dimension dh
Ci   de         work array of dimension e0
Co Outputs
Co   e0,z       depending on iopt
Cr Remarks
Cr   dz_k = sum_m z_m A_mk / (E_m - E_k)
Cr   A_mk = <z_m dH z_k>
C ----------------------------------------------------------------------
C     implicit none
      integer ld,iopt
      double precision e0(ld),dh(ld,ld,2),z(ld,ld,2),wk(ld,ld,2),de(ld)
      double precision ds(ld,ld,2)
      integer:: l2 , k , m , i , j 
      complex(8) ,allocatable :: zv_w_ow(:)

      double precision tol,xx,ddot,dr,di

      tol = 1d-12
      l2 = ld**2

C --- First order correction for eigenvalues ---
      call pvpert(iopt,ld,e0,dh,ds,z,z,wk,de)
      if (mod(iopt/10,10) .eq. 0) then
        call daxpy(ld,1d0,de,1,e0,1)
        return
      endif

C --- Eigenvectors, 1st order perturbation theory; overwrite dh ---
C ... hang onto dh for 2nd order perturbation correction
      if (mod(iopt/100,10) .eq. 1) then
         allocate(zv_w_ow(l2))
         call dcopy ( l2 * 2 , dh , 1 , zv_w_ow , 1 ) 
      endif
C ... Make A_mk = <z_m dH z_k>
      call yygemm('C','N',ld,ld,ld,1d0,z,z(1,1,2),ld,
     .  wk,wk(1,1,2),ld,0d0,dh,dh(1,1,2),ld)
C ... Make A_mk/(em-ek) and  sum_k |A_mk|/(em-ek)
c      call dpzero(wk,ld)
      wk=0d0
      do  40  k = 1, ld
      do  40  m = 1, ld
        if (abs(e0(m)-e0(k)) .lt. tol) then
          dh(m,k,1) = 0
          dh(m,k,2) = 0
        else
C          dh(m,k,1) = dh(m,k,1) / (e0(m)+de(m)-e0(k)-de(k))
C          dh(m,k,2) = dh(m,k,2) / (e0(m)+de(m)-e0(k)-de(k))
C          wk(m,1,1) = wk(m,1,1) +
C     .      (dh(m,k,1)**2 + dh(m,k,2)**2) / (e0(m)-e0(k))
          dh(m,k,1) = dh(m,k,1) / (e0(m)-e0(k))
          dh(m,k,2) = dh(m,k,2) / (e0(m)-e0(k))
        endif
   40 continue
C     call yprm('amk/(em-ek)',2,dh,ld*ld,ld,ld,ld)
C     call yprm('amk akm /(em-ek)',1,wk,ld*1,ld,ld,1)
C ... dz_k = z_m A_mk/(em-ek)
      call yygemm('N','N',ld,ld,ld,1d0,z,z(1,1,2),ld,
     .  dh,dh(1,1,2),ld,0d0,wk,wk(1,1,2),ld)
C     call yprm('dz',2,wk,ld*ld,ld,ld,ld)

C ... Restore w to dh; save z0 in w; add pert. correction to z
      if (mod(iopt/100,10) .eq. 1) then
        call dcopy ( l2 * 2 , zv_w_ow , 1 , dh , 1 ) 
        call dcopy ( l2 * 2 , z , 1 , zv_w_ow , 1 ) 
      endif
      call daxpy(l2*2,1d0,wk,1,z,1)
C     call yprm('new z',2,z,ld*ld,ld,ld,ld)

C ... Orthogonalize and/or normalize z
      if (mod(iopt/10,10) .gt. 1) then
        do  50  j = 1, ld
          if (mod(iopt/10,10) .gt. 2) then
          do  55  i = 1, j-1
           call yydotc(ld,z(1,i,1),z(1,i,2),1,z(1,j,1),z(1,j,2),1,dr,di)
           call yyaxpy(ld,-dr,-di,z(1,i,1),z(1,i,2),1,z(1,j,1),z(1,j,2),
     .                 1,.true.)
   55     continue
          endif
          xx = dsqrt(ddot(ld,z(1,j,1),1,z(1,j,1),1) +
     .               ddot(ld,z(1,j,2),1,z(1,j,2),1))
          call dscal(ld,1/xx,z(1,j,1),1)
          call dscal(ld,1/xx,z(1,j,2),1)
   50   continue
C     call yprm('after orthogonalization z',2,z,ld*ld,ld,ld,ld)
      endif

C ... Add perturbation to e0 
      call daxpy(ld,1d0,de,1,e0,1)
C     call prmx('1st order e',e0,ld,ld,1)

C --- Second order correction for eigenvalues ---
C      if (mod(iopt/100,10) .eq. 1) then
C        print *, 'second order'
Cc        call pvpert(iopt,ld,e0,dh,ds,w(ow),z,wk,de)
C        call rlse(ow)
C      endif
       if (allocated(zv_w_ow)) deallocate(zv_w_ow)

      end

      subroutine pvpert(iopt,ld,e0,dh,ds,z0,z,wk,de)
C- First order perturbation theory, correction to eigenvalues
C  iopt=1: treat dh -e0 ds as the perturbation
C  1st order pert theory, z and z0 are the same
C  return in wk:  dh z0 - ds z0 e0
C     implicit none
      integer ld
      double precision e0(ld),wk(ld,ld,2),de(ld)
      double precision dh(ld,ld,2),ds(ld,ld,2)
      double precision z(ld,ld,2),z0(ld,ld,2)
      integer ndgmx
      parameter (ndgmx=3)
      integer i,j,low,high,iwk(ndgmx),ndg,k,iopt
      double precision tol,xxv(ndgmx,2),xx
      tol = 1d-12

C --- dh - eps ds * z ---
      xx = 0
      if (mod(iopt,10) .eq. 1) then
        call yygemm('N','N',ld,ld,ld,1d0,ds,ds(1,1,2),ld,
     .    z,z(1,1,2),ld,0d0,wk,wk(1,1,2),ld)
        do  12  j = 1, ld
        do  12  i = 1, ld
        wk(i,j,1) = -wk(i,j,1)*e0(j)
   12   wk(i,j,2) = -wk(i,j,2)*e0(j)
        xx = 1
      endif
      call yygemm('N','N',ld,ld,ld,1d0,dh,dh(1,1,2),ld,
     .  z,z(1,1,2),ld,xx,wk,wk(1,1,2),ld)
C     call yprm('dh - e0 ds * z',2,wk,ld*ld,ld,ld,ld)

C --- First order eigenvalue correction ---
c      call dpzero(de,ld)
      de=0d0
      do  10  i = 1, ld
      do  10  j = 1, ld
   10 de(i) = de(i) + z0(j,i,1)*wk(j,i,1) + z0(j,i,2)*wk(j,i,2)
C     call prmx('pvpert de, 1st order',de,ld,ld,1)

C --- Permute z,hz,de,e0 for any degeneracies that were lifted ---
      low = 0
      high = 0
      i = 0
   20 i = i+1
        call huntx(e0,ld,e0(i)+tol,0,high)
        low = i
        high = min(high,ld)
C   ... Degeneracies between low,high
        if (high .ne. low) then
          ndg = high-low+1
          if (ndg .gt. ndgmx)
     .      call rxi('pertev: increase ndgmx; need',ndg)
          do  21  k = low, high
   21     xxv(k-low+1,1) = e0(k) + de(k)
          call dvheap(1,ndg,xxv,iwk,tol,1)
          do  22  j = 1, ld
            do  23  k = low, high
            xxv(iwk(k-low+1),1) = z(j,k,1)
   23       xxv(iwk(k-low+1),2) = z(j,k,2)
            do  24  k = low, high
            z(j,k,1) = xxv(k-low+1,1)
   24       z(j,k,2) = xxv(k-low+1,2)
            do  25  k = low, high
            xxv(iwk(k-low+1),1) = wk(j,k,1)
   25       xxv(iwk(k-low+1),2) = wk(j,k,2)
            do  26  k = low, high
            wk(j,k,1) = xxv(k-low+1,1)
   26       wk(j,k,2) = xxv(k-low+1,2)
   22     continue
C     ... Not necessary to sort e0, but do anyway
          do  28  k = low, high
          xxv(iwk(k-low+1),1) = e0(k)
   28     xxv(iwk(k-low+1),2) = de(k)
          do  29  k = low, high
          e0(k) = xxv(k-low+1,1)
   29     de(k) = xxv(k-low+1,2)
        endif
      if (i .lt. ld) goto 20
C     call prmx('de, sorted 1st order',de,ld,ld,1)
C     call yprm('z, permuted',2,z,ld*ld,ld,ld,ld)
C     call yprm('dh z, permuted',2,wk,ld*ld,ld,ld,ld)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

