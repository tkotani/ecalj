module  m_hansrz
  use m_hansr,only: hansr,hanr
contains
subroutine hansrz(rsml,lmin,lmax,e,rsq,nrx,nr,job,xi,fi)
  !- Vector of smoothed Hankel functions, any e
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   rsml    : vector of l-dependent smoothing radii of Hankels
  !i           : EITHER must be specified for lmin..lmax
  !i           : OR     rsmh(0) = const < 0. Implies rsmh(l) = -const for all l
  !i   nrx,lmin: dimensions of xi and fi
  !i   e,lmax  : energy and lmax to generate fns
  !i   rsq,nr  : vector of points r**2, and number of points.
  !i   job     : 1s digit nonzero, scale xi and fi by r**l
  !i           : 10s digit nonzero, input rsq is ordered as rsq(i-1) =< rsq(i)
  !o Outputs
  !o   xi   : e>0 smoothed Neumann fns for: xi(1..nr, lmin..lmax)
  !o        : e<0 smoothed Hankel fns for: xi(1..nr, lmin..lmax)
  !o   fi   : e>0 proportional to unsmoothed Bessel fns:
  !o        :     fi(1..nr, lmin..lmax) = e^{l+1/2} * j_l
  !o        : e<0 not referenced
  !o        : xi and fi are the radial part/r**l, so the solid Hankel is
  !o        : hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
  !l Local variables
  !l   n1,n2: Assuming points are sorted: points are evaluated as follows
  !l        : 1..n1-1  are evaluated by power series expansion
  !l        : n1..n2-1 are evaluated by explicit generation of l=-1,0
  !l                    and upward recursion for higher l.
  !l        : n2..nr   asymtotic form: sm-H has become regular H.
  !l        : If the points are not sorted, they are grouped into
  !l        : three bins, with n1-1, n2-n1, and nr-n2+1 points in them
  !l        : idx is a permutation index that keeps track of the grouping
  !l Local auto-arrays:
  !l   idx  : permutation index; not used if 10s digit of job set.
  !l   wrsq : reordered rsq array: rsq(i) = wrsq(idx(i));
  !l        : not used if 10s digit of job set.
  !l   gs   : gaussians for the intermediate points
  !l   wk   : double precision work array of length nrx (for hansr2)
  !r Remarks
  !r   Unlike hansr, the program scales xi(*,-1) and fi(*,-1) in the
  !r   same way as it is done for other l's if mod(job,10)>0.
  !r   This is solely done for the sake of consistency.
  !r
  !r   Points are partitioned into three length scales:
  !r     r<rc1 are evaluated by a polynomial expansion.
  !r     rc1<r<rc2 are evaluated from error functions and the higher l's
  !r       by upward recursion.
  !r     rc2<r are approximated with unsmoothed Hankels.
  !r   The relative error should be less than parameter 'tol', except
  !r   in a narrow region for r~rc1 and l>6, where the precision degrades
  !r   somewhat, worsening with higher l.  For all cases tested, the
  !r   relative error continued to be ~<10^-13 for l<=9.
  !u Updates
  !u   11 May 07 (S. Lozovoi) adapted from hansr.f
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nrx,nr,lmin,lmax,job
  double precision :: rsq(nrx),rsml(lmin:lmax),e
  double precision :: xi(nrx,lmin:lmax),fi(nrx,lmin:lmax)
  ! Local variables
  integer :: n0,ir,n1,n2,n2p,jerf,il,ill,llmax
  double precision :: a,rsm,rsm0,y0,a2,tol
  parameter (n0=10, tol=1d-15)
  double precision :: rc1,rc2,akap,akl,rl,rl0,rm2
  double precision :: fbess(-1:n0),fneum(-1:n0),rsmx(-1:n0)
  logical :: lsort,lscal,lpos
  ! Work arrays as automatic arrays
  integer :: idx(nrx)
  double precision :: whan(nrx,lmin:lmax),wbes(nrx,lmin:lmax)
  double precision :: wk(nrx),wrsq(nrx),gs(nrx)
  !     double precision, allocatable :: whan(:,:),wbes(:,:)

  lscal = mod(job,10) .ne. 0
  lsort = mod(job/10,10) .ne. 0
  lpos = e .gt. 0

  !     if (.not. lsort)
  !    . allocate (whan(nr,lmin:lmax),wbes(nr,lmin:lmax))

  ! --- Check lmax ---
  if (lmax < 0) call rx('hansrz: lmax is negative')
  if (lmax > n0) call rx('hansrz: increase parameter n0')

  ! ... Handle negative smoothing radii
  if (rsml(0) < 0d0) then
     rsm = - rsml(0)
     call dvset(rsmx(lmin),1,lmax-lmin+1,rsm)
  else
     call dcopy(lmax-lmin+1,rsml(lmin),1,rsmx(lmin),1)
  endif

  y0 = 1/dsqrt(16*datan(1d0))
  akap = dsqrt(dabs(e))
  ! Start big loop over smoothing radii
  rsm0 = -1d0
  n2p = nr + 1
  do  ill = lmax, lmin, -1
     rsm = rsmx(ill)
     if (dabs(rsm-rsm0) < tol) goto 60
     rsm0 = rsm
     llmax = max0(ill,0)

     ! --- Handle case rsm=0 ---
     if (rsm < 1d-12) then
        n2 = 1
        if (lpos) then
           do  ir = 1, nr
              rl = rsq(ir)
              call besslr(e*rl,0,lmin,llmax,fbess(lmin),fneum(lmin))
              rm2 = 1.d0/rl
              rl = dsqrt(rl)
              if (lmin == -1) then
                 fi(ir,-1) = fbess(lmin)
                 whan(ir,-1) = fneum(lmin) * rl
              endif
              do il = lmin, llmax
                 rl = rl * rm2
                 fi(ir,il) = fbess(il)
                 whan(ir,il) = fneum(il) * rl
              enddo
           enddo
        else
           call hanr(rsq(n2),lmin,llmax,nrx,nr+1-n2,e,whan(n2,lmin))
        endif
        goto 60
     endif

     ! --- Find cutoffs rc2 (negligible smoothing) and rc1 (power series) ---
     a = 1/rsm
     a2 = a*a
     !   ... For r>rc2 approximate smooth Hankels with normal ones
     rc2 = akap/(2*a)
     rc2 = ((rc2 + dsqrt(rc2**2 - dlog(tol)))/a)**2
     !   ... This rc1 generates a relative precision of ~10^-15 for r~rc1
     !       and machine precision for r>>rc1 or r<<rc1.
     !       For l>6 and r close to rc1, the precision degrades somewhat.
     rc1 = (rsm*(1.4d0+dble(llmax)/20))**2

     ! --- Partition the mesh: find n1, n2 and make idx, wrsq
     !       n1 is offset to block rc1<r<rc2,  n2 offset to block r>rc2
     !       idx is a map of original list, separating into the three groups
     !       wrsq is a table of r**2 for permuted list of points

     call rsort(nr,rsq,rc1,rc2,lsort,n1,n2,idx,wrsq)

     !   ... For debugging
     !       call awrit8(' hansrz: sort=%l rc1=%,2;2d rc2=%,2;2d'//
     !    .    ' lmax=%i nr=%i (%i pwr, %i smooth, %i asym)',
     !    .    ' ',120,6,lsort,dsqrt(rc1),dsqrt(rc2),llmax,nr,
     !    .    n1-1,n2-n1,nr-n2+1)


     ! --- Setup for the energy-independent gs, points n1..n2 ---
     if (lsort) then
        do  ir = max0(n1,1), n2-1
           gs(ir) = y0*dexp(-rsq(ir)*a2)
        enddo
     else
        do  ir = max0(n1,1), n2-1
           gs(ir) = y0*dexp(-wrsq(ir)*a2)
        enddo
     endif
     !       endif

     jerf = 0
     !   ... Case calculate points in original order (already sorted)
     if (lsort) then
        !     ... Power series for points within rc1
        if (n1 > 0) &
             call hansz1(rsq(1),lmin,llmax,nrx,n1-1,e,rsm, &
             dsqrt(rc1),whan(1,lmin))
        !     ... Normal evaluation of smoothed Hankels
        call hansz2(jerf,rsq(n1),lmin,llmax,nrx,n2-n1,e,rsm,gs(n1), &
             wk(n1),whan(n1,lmin))
        !     ... Asymtotic case, r>>rsm
        if (n2p > n2) then
           if (lpos) then
              do  ir = n2, n2p - 1
                 rl = rsq(ir)
                 call besslr(e*rl,0,lmin,llmax,fbess(lmin),fneum(lmin))
                 !               call bessll(e*rl,lmin,llmax,fbess(lmin),fneum(lmin))
                 rm2 = 1.d0/rl
                 rl = dsqrt(rl)
                 if (lmin == -1) then
                    fi(ir,lmin) = fbess(lmin)
                    whan(ir,lmin) = fneum(lmin) * rl
                 endif
                 do il = 0, llmax
                    rl = rl * rm2
                    fi(ir,il) = fbess(il)
                    whan(ir,il) = fneum(il) * rl
                 enddo
              enddo
           else
              call hanr(rsq(n2),lmin,llmax,nrx,n2p-n2,e,whan(n2,lmin))
           endif
        endif
        !     ... if e > 0, make also unsmoothed bessel functions  ---
        ! we have to do this only once, in addition bessels from n2 to nr already exist
        if (lpos .AND. llmax == lmax) then
           if (n2 > 1) then
              do  ir = 1, n2-1
                 !               call bessll(e*rsq(ir),lmin,llmax,
                 call besslr(e*rsq(ir),0,lmin,llmax, &
                      fbess(lmin),fneum(lmin))
                 do il = lmin, llmax
                    fi(ir,il) = fbess(il)
                 enddo
              enddo
           endif
        endif

        !   ... Case calculated points are not sorted
     else
        !     ... Power series for points within rc1
        if (n1 > 0) &
             call hansz1(wrsq(1),lmin,llmax,nrx,n1-1,e,rsm,dsqrt(rc1), &
             whan(1,lmin))
        !     ... Normal evaluation of smoothed Hankels
        call hansz2(jerf,wrsq(n1),lmin,llmax,nrx,n2-n1,e,rsm,gs(n1), &
             wk(n1),whan(n1,lmin))
        !     ... Asymtotic case, r>>rsm
        if (nr >= n2) then
           if (lpos) then
              do  ir = n2, nr
                 rl = wrsq(ir)
                 call besslr(e*rl,0,lmin,llmax,fbess(lmin),fneum(lmin))
                 !               call bessll(e*rl,lmin,llmax,fbess(lmin),fneum(lmin))
                 rm2 = 1.d0/rl
                 rl = dsqrt(rl)
                 if (lmin == -1) then
                    wbes(ir,lmin) = fbess(lmin)
                    whan(ir,lmin) = fneum(lmin) * rl
                 endif
                 do il = 0, llmax
                    rl = rl * rm2
                    wbes(ir,il) = fbess(il)
                    whan(ir,il) = fneum(il) * rl
                 enddo
              enddo
           else
              call hanr(wrsq(n2),lmin,llmax,nrx,nr+1-n2,e,whan(n2,lmin))
           endif
        endif
        !     ... if e > 0, make also unsmoothed bessel functions  ---
        if (lpos .AND. llmax == lmax) then
           if (n2 > 1) then
              do  ir = 1, n2-1
                 !               call bessll(e*wrsq(ir),lmin,llmax,
                 call besslr(e*wrsq(ir),0,lmin,llmax, &
                      fbess(lmin),fneum(lmin))
                 do il = lmin, llmax
                    wbes(ir,il) = fbess(il)
                 enddo
              enddo
           endif
           !     ... Poke bessels into fi, with the original ordering of points
           do  il = lmin, llmax
              do  ir = 1, nr
                 fi(ir,il) = wbes(idx(ir),il)
              enddo
           enddo
        endif
     endif
60   continue
     n2p = n2
     !     ... for given l, poke hankels into xi
     if (lsort .OR. rsm < 1d-12) then
        do  ir = 1, nr
           xi(ir,ill) = whan(ir,ill)
        enddo
     else
        do  ir = 1, nr
           xi(ir,ill) = whan(idx(ir),ill)
        enddo
     endif
     !       deallocate (whan,wbes)

     ! --- End big loop over smoothing radii
  enddo

  ! --- Scale bessels by k^(2l+1)
  if (lpos) then
     akl = akap
     if (lmin == -1) call dscal(nr,1d0/akl,fi(1,-1),1)
     do  il = 0, lmax
        call dscal(nr,akl,fi(1,il),1)
        akl = akl * e
     enddo
  endif

  ! --- Scale by r**l if job nonzero ---
  !     if (.not. lscal) return
  if (lscal) then
     do  ir = 1, nr
        rl0 = dsqrt(rsq(ir))
        if (lmin == -1) then
           xi(ir,-1) = xi(ir,-1)/rl0
           if (lpos) fi(ir,-1) = fi(ir,-1)/rl0
        endif
        if (lmax >= 1) then
           rl = rl0
           do  il = 1, lmax
              xi(ir,il) = xi(ir,il)*rl
              if (lpos) fi(ir,il) = fi(ir,il)*rl
              rl = rl*rl0
           enddo
        endif
     enddo
  endif

end subroutine hansrz
subroutine hansz1(rsq,lmin,lmax,nrx,nr,e,rsm,rmax,xi)
  !- Vector of smoothed hankel (e < 0 ) or neumann functions (e > 0)
  !- for l=lmin...lmax  by power series expansion.
  ! ---------------------------------------------------------------
  !i Inputs
  !i   rsq,nr : vector of points r**2, and number of points.
  !i   nrx    : dimensions xi; nrx must be gt nr.
  !i   e,rsm  : smoothing radius and energy
  !i   lmin   : starting l for which to evaluate xi (must be 0 or 1).
  !i   lmax   : highest l for which to evaluate xi (must be 0 or 1).
  !i   rmax:  points rsq are less than rmax**2.
  !o Outputs:
  !o   xi(1..nr,lmin:lmax) : radial part of smoothed hankel or neumann
  !o                         functions divided by r**l
  !r Remarks
  !r   This routine is intended for evaluation of smoothed hankels
  !r   for small r (r<rsm or so).
  !r   hansz1 tries to evaluate the polynomial in-line for a 14th order
  !r   polynomial, and a 20th order.  Failing that, it evaluates the
  !r   polynomial to whatever order is needed (up to nmax) to bring
  !r   the convergence to a relative precision of 'tol'.
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: nrx,nr,lmin,lmax,nmax,nm1,nm2
  double precision :: rsq(nrx),e,xi(nrx,lmin:lmax),rsm,rmax
  ! Local vairables
  parameter (nmax=40,nm1=14,nm2=20)
  double complex zerfc
  double precision :: derfc,cof0(-1:20),cofl(0:nmax),tol
  double precision :: a,a2,add,akap,cc,fac,rhs,ta,ta2l,y0,r2max,x
  integer :: i,il,l,ir,m,nmaxl
  logical :: lpos
  parameter (tol=1d-20)

  ! --- Setup ---
  if (lmax < lmin .OR. nr == 0) return
  if (lmin < -1 .OR. lmin > 0) call rx('hansz1: bad lmin')
  if (lmax < 0) call rx('hansz1: bad lmax')

  lpos = (e .gt. 0d0)
  y0 = 1/dsqrt(16*datan(1d0))
  a = 1/rsm
  ta = a+a
  a2 = a*a
  akap = dsqrt(dabs(e))
  cc = 4d0*y0*a*dexp(e/(ta*ta))
  r2max = rmax**(2*nm1)

  ! --- 0 order coefficients ---
  if (lpos) then
     cof0(-1) = dimag(zerfc(dcmplx(0d0,akap/ta))) / akap
  else
     cof0(-1) = derfc(akap/ta) / akap
  endif
  fac = cc
  do  il = 0, lmax
     cof0(il) = (e*cof0(il-1) + fac) / (2*il + 1)
     fac = fac * (2*a2)
  enddo

  ! --- For each l, generate xi(*,l) by power series ---
  ta2l = (2*a2)
  if (lmin == -1) cc = cc/ta2l
  do  20  l = lmin, lmax
     rhs = cc*ta2l
     add = cof0(l)
     cofl(0) = add
     !   --- Coffs to polynomial of order nm1 ---
     do  21  i = 1, nm1
        add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
        cofl(i) = add
        rhs = -rhs*a2/i
21   enddo
     !   ... Use it if it is accurate enough
     if (dabs(add*r2max) < dabs(cof0(l))*tol) then
        do  51  ir = 1, nr
           x = rsq(ir)
           xi(ir,l) = (((((((((((((cofl(14)*x+ &
                cofl(13))*x+cofl(12))*x+cofl(11))*x+cofl(10))*x+cofl(9))*x+ &
                cofl(8))*x+cofl(7))*x+cofl(6))*x+cofl(5))*x+cofl(4))*x+ &
                cofl(3))*x+cofl(2))*x+cofl(1))*x+cofl(0)
51      enddo
        ! ---   Coffs to polynomial of order nm2 ---
        !   ... Use it if it is accurate enough
     else
        do  22  i = nm1+1, nm2
           add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
           cofl(i) = add
           rhs = -rhs*a2/i
22      enddo
        if (dabs(add*r2max) < dabs(cof0(l))*tol) then
           do  52  ir = 1, nr
              x = rsq(ir)
              xi(ir,l) = (((((((((((((((((((cofl(20)*x+cofl(19))*x+ &
                   cofl(18))*x+cofl(17))*x+cofl(16))*x+cofl(15))*x+cofl(14))*x+ &
                   cofl(13))*x+cofl(12))*x+cofl(11))*x+cofl(10))*x+cofl(9))*x+ &
                   cofl(8))*x+cofl(7))*x+cofl(6))*x+cofl(5))*x+cofl(4))*x+ &
                   cofl(3))*x+cofl(2))*x+cofl(1))*x+cofl(0)
52         enddo
        else
           ! ---   Polynomial to nmaxl ---
           do  23  i = nm2+1, nmax
              add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
              cofl(i) = add
              nmaxl = i
              if (dabs(add*r2max) < dabs(cof0(l))*tol) goto 24
              rhs = -rhs*a2/i
23         enddo
           !       print 333, tol
           ! 333   format(' hansz1:  not converged to tol=',1pe8.1)
           print '(/1x,a,e8.1,a,i3)','hansz1: only ', &
                dabs(add*r2max/cof0(l)), &
                ' accuracy is achieved for l = ',l
           call rx1('hansz1: not converged to tol = %;4g',tol)
24         continue
           do  53  ir = 1, nr
              xi(ir,l) = cofl(nmaxl)
53         enddo
           do    m = nmaxl, 1, -1
              do    ir = 1, nr
                 xi(ir,l) = xi(ir,l)*rsq(ir) + cofl(m-1)
              enddo
           enddo
        endif
     endif
     ta2l = ta2l*(ta*a)
20 enddo

end subroutine hansz1
subroutine hansz2(jerf,rsq,lmin,lmax,nrx,nr,e,rsm,wk,wk2,xi)
  !- Vector of smoothed hankel (e < 0) or neumann (e > 0) functions
  !- for l=lmin...lmax.
  ! ---------------------------------------------------------------
  !i Inputs
  !i   jerf  : 1 - approximate erfc as a ratio of polinomials (e < 0)
  !i           0 - calculate erfc directly
  !i   rsq,nr: vector of points r**2, and number of points.
  !i   lmin  : starting l for which to evaluate xi (must be 0 or -1).
  !i   lmax  : highest l for which to evaluate xi (must be 0 or > 0)
  !i   e,rsm : smoothing radius and energy
  !i   nrx   : leading dimension of xi
  !i   wk    : array containing y0*dexp(-(r/rsm)**2)
  !i   wk2   : a work array of length nr.
  !o Outputs:
  !o   xi    : generated for points ir=1..nr and lmin..lmax
  !o   wk2   : (2/rsm**2)**(lmax)*4/rsm*dexp(-(akap*rsm/2)**2)*wk(ir)
  !o           (can be used to generate xi to higher l)
  !r Remarks
  !r   xi is the radial part divided by r**l.
  !r   xi is evaluated by upward recursion for l>lmin+2.
  !r   um and up are u_{+} and u_{-}, respectively, from the JMP paper
  !r             (ie the order is opposite)
  ! ---------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: jerf,nrx,nr,lmin,lmax
  double precision :: rsq(nrx),e,xi(nrx,lmin:lmax),rsm,wk(nr),wk2(nr)
  ! Local variables
  double precision :: sre,r2,a2,xx,ra,h0,arsm,earsm
  double precision :: akap,a,r,um,up,x,facgl,facdu,dudr
  integer :: il,ir
  logical :: lpos
  ! ... erfc(x) is evaluated as a ratio of polynomials,
  !     to a relative precision of <10^-15 for x<5.
  !     Different polynomials are used for x<1.3 and x>1.3.
  !     Numerators and denominators are t,b respectively.
  double precision :: w,f1,f2, &
       t10,t11,t12,t13,t14,t15,t16,t17,b11,b12,b13,b14,b15,b16,b17,b18, &
       t20,t21,t22,t23,t24,t25,t26,t27,b21,b22,b23,b24,b25,b26,b27,b28
  parameter ( &
       t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0, &
       t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0, &
       t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0, &
       t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0, &
       t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0, &
       t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0, &
       t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0, &
       t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0, &
       b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0, &
       b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0, &
       b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0, &
       b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0, &
       b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0, &
       b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0, &
       b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0, &
       b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0)
  double precision :: derfc,aexp
  double complex zerfc,zexp,vm,vp
  !     double precision upm
  double precision :: rzero,facgl0,am1
  integer :: n1
  parameter (rzero = 1d-40)
  ! ... f1(w=x-1/2) is erfc(x) for 0<x<1.3, if xx is y0*dexp(-x*x)
  f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)* &
       w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)* &
       w+b13)*w+b12)*w+b11)*w+1)
  ! ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
  f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)* &
       w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)* &
       w+b23)*w+b22)*w+b21)*w+1)
  ! ... upm(w) = u_{+}(r); upm(-w) = u_{-}(r) => h(r) = (upm(r) - upm(-r))/(2*r)
  !     upm(w) = dexp(-w*akap) * derfc(arsm-a*w)

  ! --- Setup ---
  if (lmax < lmin .OR. nr == 0) return
  if (lmin < -1 .OR. lmin > 0) call rx('hansz2: bad lmin')
  if (lmax < 0) call rx('hansz2: bad lmax')

  ! ... energy
  lpos = (e .gt. 0)
  if (lpos) then
     akap = dsqrt(e)
     !       zi = dcmplx(0d0,1d0)
  else
     akap = dsqrt(-e)
  endif
  if (lpos .AND. jerf == 1) &
       call rx('hansz2:  can''t approximate erfc with polynomials at'// &
       ' e > 0')

  a = 1/rsm
  a2 = 2*a*a
  arsm = akap*rsm*.5d0
  !      if (lpos) then
  !        earsm = dexp(arsm**2)*.5d0
  !      else
  !        earsm = dexp(-arsm**2)*.5d0
  !      endif
  earsm = dexp(e*rsm**2*0.25d0)*.5d0
  facdu = 8*a*earsm
  facgl = facdu
  if (lmin == 0) facgl = facgl * a2

  ! ... If the mesh starts in the orygin, set xi(1) to lim xi(r -> 0)
  n1 = 1
  if (rsq(1) < rzero) then
     n1 = 2
     !       facgl0 = facdu / dsqrt(16*datan(1d0))
     facgl0 = facdu * wk(1)
     if (lpos) then
        am1 = dimag(zerfc(dcmplx(0d0,arsm))) / akap
     else
        am1 = derfc(arsm) / akap
     endif
     if (lmin == -1) xi(1,-1) = am1
     xi(1,0) = facgl0 + e*am1
     if (lmax >= 1) then
        do  il = 1, lmax
           xx = 2*il + 1
           facgl0 = facgl0 * a2
           xi(1,il) = (e*xi(1,il-1) + facgl0) / xx
        enddo
     endif
  endif

  ! --- xi(*,lmin), xi(*,lmin+1) ---
  !  Approximate erfc as a ratio of polynomials
  if (jerf == 1) then
     do  ir = n1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        ra = r*a
        sre = akap*r
        h0 = dexp(-sre)/r
        xx = earsm*wk(ir)/r

        ! ---   Evaluate um,up ---
        x = ra - arsm
        if (x > 1.3d0) then
           um = h0 - f2(x-2d0)
        elseif (x > 0) then
           um = h0 - f1(x-.5d0)
        elseif (x > -1.3d0) then
           um = f1(-x-.5d0)
        else
           um = f2(-x-2d0)
        endif
        ! ...   Evaluation of up assumes x gt 0
        x = ra + arsm
        if (x > 1.3d0) then
           up = f2(x-2d0)
        else
           up = f1(x-.5d0)
        endif

        !   ... xi(0) = um - up
        xi(ir,0) = um - up
        !   ... xi(-1) = (um + up)*r/akap
        if (lmin == -1) then
           xi(ir,-1) = (um + up)*r/akap
           !   ... xi(1)
        elseif (lmax >= 1) then
           dudr = facdu*wk(ir) - sre*(um+up)
           xi(ir,1) = (xi(ir,0) - dudr)/r2
        endif
        wk2(ir) = facgl*wk(ir)
     enddo

     ! ... up and um directly from derfc;
  else
     do  ir = n1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        ra = r*a
        if (lpos) then
           zexp = cdexp(dcmplx(0d0,akap*r))
           vm = (1d0 - zerfc(dcmplx( ra,arsm))) * zexp
           vp = (1d0 - zerfc(dcmplx(-ra,arsm))) / zexp
           !   ... am1 = xi(-1) = (vm + vp) * i/(2*akap)
           am1 = -dimag(vm + vp) * 0.5d0/akap
           !   ... xi(0) = (vm - vp)/(2*r)
           xi(ir,0) = dreal(vm - vp) * (0.5d0/r)
        else
           aexp = dexp(akap*r)
           um = derfc(arsm-ra) / aexp
           up = derfc(arsm+ra) * aexp
           !   ... am1 = xi(-1) = (um + up)/(2*akap)
           am1 = (um + up) * 0.5d0/akap
           !   ... xi(0) = (um - up)/(2*r)
           xi(ir,0) = (um - up) * (0.5d0/r)
        endif

        !   ... xi(-1)
        if (lmin == -1) then
           xi(ir,-1) = am1
           !   ... xi(1)
        elseif (lmax >= 1) then
           dudr = facdu*wk(ir) + e*am1
           xi(ir,1) = (xi(ir,0) - dudr)/r2
        endif
        wk2(ir) = facgl*wk(ir)
     enddo
  endif

  ! --- xi(ir,l) for l>1 by upward recursion ---
  if (lmax > lmin+1) then
     facgl = a2
     do  il = lmin+2, lmax
        xx = 2*il-1
        do  ir = n1, nr
           xi(ir,il) = (xx*xi(ir,il-1) - e*xi(ir,il-2) - &
                wk2(ir)) / rsq(ir)
           wk2(ir) = facgl*wk2(ir)
        enddo
     enddo
  endif

end subroutine hansz2


end module m_hansrz
