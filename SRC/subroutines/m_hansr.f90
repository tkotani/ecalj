!>Smooth Hankel functions in real space. JMP39: https://doi.org/doi:10.1063/1.532437.
module m_hansr 
  use m_lmfinit,only: z_i=>z,rmt_i=>rmt,lmxb_i=>lmxb,lfoca_i=>lfoca,rfoca_i=>rfoca,rg_i=>rg
  ! JMP39:
  ! Bott, E., M. Methfessel, W. Krabs, and P. C. Schmidt.
  ! “Nonsingular Hankel Functions as a New Basis for Electronic Structure Calculations.”
  ! Journal of Mathematical Physics 39, no. 6 (June 1, 1998): 3393–3425.
  ! https://doi.org/doi:10.1063/1.532437.
  public hansr,hanr,hansmd, hansmr,corprm ! hansmr is equilvaent to hansr except numerical accuracy problem. See note of hansmr
  private
contains
  subroutine hansr(rsm,lmn,lmx,nxi,lxi,exi,rsq,nrx,nr,idx,job,xi) !- Vector of smoothed Hankel functions, set of negative e's
    !i Inputs
    !i   rsm     smoothing radius of smoothed Hankel
    !i   nrx,lmx dimensions xi
    !i   nxi,exi,lxi:number of energies, energies and lmax to generate fns
    !i   rsq,nr  vector of points r**2, and number of points.
    !i   job      1s digit nonzero, scale xi by r**l
    !i           10s digit nonzero, input rsq is already sorted
    !i   idx     integer work array of length 2*nrx.
    !i           Not needed if 10s digit of job set.
    !o Outputs
    !o   xi      smoothed Hankels for: xi(1..nr, 0..lxi(ie), 1..nxi)
    !o           xi is the radial part/r**l, so the solid Hankel is
    !o           hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
    !l Local variables
    !l   n1,n2 : Assuming points are sorted: points are evaluated as follows
    !l         : 1..n1-1  are evaluated by power series expansion
    !l         : n1..n2-1 are evaluated by explicit generation of l=-1,0
    !l                    and upward recursion for higher l.
    !l         : n2..nr   asymtotic form: sm-H has become regular H.
    !l         : If the points are not sorted, they are grouped into
    !l         : three bins, with n1-1, n2-n1, and nr-n2+1 points in them
    !l         : idx is a permutation index that keeps track of the grouping
    !r Remarks
    !NOTE: r is divided into domains [0,rc1),[rc1,rc2] [rc2,\infty]
    !NOTE: hansmr is equivalent to hansr but dividing r into two length scales.
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
    !u   11 May 07 (S. Lozovoi) small bug fixes; similarity with hansmz
    ! ---------------------------------------------------------------
    implicit none
    integer :: nrx,nr,lmn,lmx,idx(nrx,2),nxi,lxi(nxi),job
    real(8) :: rsq(nr),e,exi(nxi)
    real(8) :: xi(nrx,lmn:lmx,nxi),wk(nrx,4+lmx-lmn)
    integer :: ir,j,l,k,n0,n1,n2,ie,lmax
    real(8) :: a,rsm,y0,a2,emin,tol
    real(8) :: rc1,rc2,akap,rl,rl0
    parameter (tol=1d-15)
    logical :: ltmp,lsort,lscal
    lscal = mod(job,10) .ne. 0
    lsort = mod(job/10,10) .ne. 0
    ! --- Check lmx; handle case rsm=0 ---
    ltmp = rsm .lt. 1d-12
    lmax = -1
    do  5  ie = 1, nxi
       if (lxi(ie) > lmx) call rx('hansr: lxi gt lmx')
       if (exi(ie) > 0)   call rx('hansr: exi gt 0')
       if (ltmp) then
          call hanr(rsq,lmn,lxi(ie),nrx,nr,exi(ie),xi(1,lmn,ie))
       endif
       lmax = max(lmax,lxi(ie))
5   enddo
    if (ltmp) goto 60
    ! --- Find cutoffs rc2 (negligible smoothing) and rc1 (power series) ---
    emin = 0
    do  10  ie = 1, nxi
       emin = min(emin,exi(ie))
10  enddo
    akap = dsqrt(-emin)
    a = 1/rsm
    a2 = a*a
    y0 = 1/dsqrt(16*datan(1d0))
    ! ... For r>rc2 approximate smooth Hankels with normal ones
    rc2 = akap/(2*a)
    rc2 = ((rc2 + dsqrt(rc2**2 - dlog(tol)))/a)**2
    ! ... This rc1 generates a relative precision of ~10^-15 for r~rc1
    !     and machine precision for r>>rc1 or r<<rc1.
    !     For l>6 and r close to rc1, the precision degrades somewhat.
    rc1 = (rsm*(1.4d0+dble(lmax)/20))**2
    ! --- Separate the small from the large ---
    n0 = 0
    n1 = 0
    n2 = nr+1
    ! ... Case points already sorted.  Find n1,n2.
    if (lsort) then
       n1 = 1
       if (nr == 1) then
          if (rsq(1) < rc1) n1 = 2
          if (rsq(1) > rc2) n2 = 1
       else
          if (rsq(1) >= rc1) then
             n1 = 1
          else
             call huntx(rsq,nr,rc1,n1)
             n1 = n1+1
          endif
          if (rsq(nr) <= rc2) then
             n2 = nr+1
          else
             n2 = nr
             call huntx(rsq,nr,rc2,n2)
             n2 = n2+1
          endif
       endif

       ! ... Case points not sorted (iwk, wk(3) required now.)
    else
       !     On output, lsort is true if points already sorted.
       lsort = .true.
       do  12  ir = 1, nr
          !     n1 is offset to block rc1<r<rc2,  n2 offset to block r>rc2
          !     idx is a map of original list, separating into the three groups
          !     wk(*,3) is a table of r**2 for permuted list of points
          if (rsq(ir) < rc2) then
             if (rsq(ir) < rc1) then
                n0 = n0+1
                wk(n0,3) = rsq(ir)
                idx(ir,1) = n0
             else
                n1 = n1+1
                idx(ir,1) = n1
                idx(n1,2) = ir
             endif
          else
             n2 = n2-1
             wk(n2,3) = rsq(ir)
             idx(ir,1) = n2
          endif
          if (ir == 1) goto 12
          if (rsq(ir) < rsq(ir-1)) lsort = .FALSE. 
12     enddo
       ! ... Now we can poke wk(*,3) for the n1 intermediate points
       if ( .NOT. lsort .OR. .TRUE. ) then
          do  14  j = 1, n1
             k = idx(j,2)
             idx(k,1) = idx(k,1)+n0
             wk(n0+j,3) = rsq(k)
14        enddo
       endif
       n1 = n0+1
    endif
    ! --- Setup for the energy-independent wk, points n1..n2 ---
    if (lsort) then
       do  20  ir = n1, n2-1
          wk(ir,1) = y0*dexp(-rsq(ir)*a2)
20     enddo
    else
       do  22  ir = n1, n2-1
          wk(ir,1) = y0*dexp(-wk(ir,3)*a2)
22     enddo
    endif
    ! --- Start loop over energies ---
    do  40  ie = 1, nxi
       e = exi(ie)
       akap = dsqrt(-e)
       lmax = lxi(ie)
       !   ... Case calculate points in original order (already sorted)
       if (lsort) then
          !     ... Power series for points within rc1
          call hansr1(rsq(1),lmn,lmax,nrx,n1-1,e,rsm,dsqrt(rc1), xi(1,lmn,ie))
          !     ... Normal evaluation of smoothed Hankels
          if (n1 <= nr) call hansr2(rsq(n1),lmn,lmax,nrx,n2-n1,e,rsm,wk(n1,1),wk(n1,2),xi(n1,lmn,ie))
          !     ... Asymtotic case, r>>rsm
          if(n2 <= nr) call hanr(rsq(n2),lmn,lmax,nrx,nr+1-n2,e,xi(n2,lmn,ie))
          !   ... Case calculated points in sorted
       else
          !     ... Power series for points within rc1
          call hansr1(wk(1,3),lmn,lmax,nrx,n1-1,e,rsm,dsqrt(rc1), wk(1,4))
          !     ... Normal evaluation of smoothed Hankels
          call hansr2(wk(n1,3),lmn,lmax,nrx,n2-n1,e,rsm,wk(n1,1), wk(n1,2),wk(n1,4))
          !     ... Asymtotic case, r>>rsm
          if(n2 <= nr) call hanr(wk(n2,3),lmn,lmax,nrx,nr+1-n2,e,wk(n2,4))
          !     ... Poke into xi(lmn:lmax), with the original ordering of points
          do   l = lmn, lmax
             do    ir = 1, nr
                j = idx(ir,1)
                xi(ir,l,ie) = wk(idx(ir,1),4+l-lmn)
             enddo
          enddo
       endif
40  enddo
    ! --- Scale by r**l if job nonzero ---
60  continue
    if ( .NOT. lscal) return
    do   ir = 1, nr
       rl0 = dsqrt(rsq(ir))
       do   ie = 1, nxi
          rl = rl0
          do   l = 1, lxi(ie)
             xi(ir,l,ie) = xi(ir,l,ie)*rl
             rl = rl*rl0
          enddo
       enddo
    enddo
  end subroutine hansr
  subroutine hansr1(rsq,lmin,lmax,nrx,nr,e,rsm,rmax,xi)!Power series for points within rc1. Vector of smoothed hankel functions for l=0...lmax, negative e
    !  by power series expansion.
    !i Inputs
    !i   rsq,nr vector of points r**2, and number of points.
    !i   nrx    dimensions xi; nrx must be gt nr.
    !i   e,rsm  smoothing radius and energy
    !i   lmin   starting l for which to evaluate xi (must be 0 or 1).
    !i   lmax   highest l for which to evaluate xi (must be 0 or 1).
    !i   rmax:  points rsq are less than rmax**2.
    !o Outputs:
    !o   xi(1..nr,lmin:lmax)
    !r Remarks
    !r   xi is the radial part divided by r**l.
    !r   This routine is intended for evaluation of smoothed hankels
    !r   for small r (r<rsm or so).
    !r   hansr1 tries to evaluate the polynomial in-line for a 14th order
    !r   polynomial, and a 20th order.  Failing that, it evaluates the
    !r   polynomial to whatever order is needed to bring the convergence
    !r   to a relative precision of 'tol'.
    ! ---------------------------------------------------------------
    implicit none
    integer :: nrx,nr,lmin,lmax,nmax,nm1,nm2
    real(8) :: rsq(nrx),e,xi(nrx,lmin:lmax),rsm,rmax
    parameter (nmax=80,nm1=14,nm2=20)
    real(8) :: cof0(-1:20),cofl(0:nmax),tol
    real(8) :: a,a2,add,akap,al,cc,fac,rhs,ta,ta2l,y0,r2max,x,derfc
    integer :: i,l,ir,m,nmaxl
    parameter (tol=1d-20)
    if (lmax < lmin .OR. nr == 0) return
    if (lmin < -1 .OR. lmin > 0) call rx('hansr1: bad lmin')
    y0 = 1/dsqrt(16*datan(1d0))
    a = 1/rsm
    ta = a+a
    a2 = a*a
    akap = dsqrt(-e)
    cc = 4d0*y0*a*dexp(e/(ta*ta))
    r2max = rmax**(2*nm1)
    ! --- 0 order coefficient ---
    fac = derfc(akap/ta)
    cof0(-1) = fac/akap
    cof0(0)  = cc - akap*fac
    al = cof0(0)
    rhs = cc*(2*a2)
    fac = 1d0
    do  10  l = 1, lmax
       al = -(e*al + rhs) / (2*l*(2*l+1))
       rhs = -rhs*a2/l
       fac = -2d0*fac*l
       cof0(l) = fac*al
10  enddo
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
21     enddo
       !   ... Use it if it is accurate enough
       if (dabs(add*r2max) < cof0(l)*tol) then
          do  51  ir = 1, nr
             x = rsq(ir)
             xi(ir,l) = (((((((((((((cofl(14)*x+ &
                  cofl(13))*x+cofl(12))*x+cofl(11))*x+cofl(10))*x+cofl(9))*x+ &
                  cofl(8))*x+cofl(7))*x+cofl(6))*x+cofl(5))*x+cofl(4))*x+ &
                  cofl(3))*x+cofl(2))*x+cofl(1))*x+cofl(0)
51        enddo
          ! ---   Coffs to polynomial of order nm2 ---
          !   ... Use it if it is accurate enough
       else
          do  22  i = nm1+1, nm2
             add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
             cofl(i) = add
             rhs = -rhs*a2/i
22        enddo
          if (dabs(add*r2max) < cof0(l)*tol) then
             do  52  ir = 1, nr
                x = rsq(ir)
                xi(ir,l) = (((((((((((((((((((cofl(20)*x+cofl(19))*x+ &
                     cofl(18))*x+cofl(17))*x+cofl(16))*x+cofl(15))*x+cofl(14))*x+ &
                     cofl(13))*x+cofl(12))*x+cofl(11))*x+cofl(10))*x+cofl(9))*x+ &
                     cofl(8))*x+cofl(7))*x+cofl(6))*x+cofl(5))*x+cofl(4))*x+ &
                     cofl(3))*x+cofl(2))*x+cofl(1))*x+cofl(0)
52           enddo
          else
             ! ---   Polynomial to nmaxl ---
             do  23  i = nm2+1, nmax
                add = -(e*add + rhs) / ( 2*i*(2*i+(l+l+1)) )
                cofl(i) = add
                nmaxl = i
                if (dabs(add*r2max) < cof0(l)*tol) goto 24
                rhs = -rhs*a2/i
23           enddo
             print 333, tol
333          format(' hansr1 (warning):  not converged to tol=',1pe8.1)
             call rx('hansr1 failed to converge')
24           continue
             do  53  ir = 1, nr
                xi(ir,l) = cofl(nmaxl)
53           enddo
             do    m = nmaxl, 1, -1
                do    ir = 1, nr
                   xi(ir,l) = xi(ir,l)*rsq(ir) + cofl(m-1)
                enddo
             enddo
          endif
       endif
       ta2l = ta2l*(2*a2)
20  enddo
  end subroutine hansr1
  subroutine hansr2(rsq,lmin,lmax,nrx,nr,e,rsm,wk,wk2,xi) !Normal evaluation of smoothed Hankels !Vector of smoothed hankel functions for l=0...lmax, negative e.
    !i Inputs
    !i   rsq,nr:vector of points r**2, and number of points.
    !i   lmin:  starting l for which to evaluate xi (must be 0 or -1).
    !i   lmax:  highest l for which to evaluate xi.
    !i   e,rsm: smoothing radius and energy
    !i   nrx:   leading dimension of xi
    !i   wk:    array containing y0*dexp(-(r/rsm)**2)
    !i   wk2:   a work array of length nr.
    !o Outputs:
    !o   xi:    generated for points ir=1..nr and lmin..lmax
    !o   wk2:   (2/rsm**2)**(lmax)*4/rsm*dexp(-(akap*rsm/2)**2)*wk(ir)
    !o          (can be used to generate xi to higher l)
    !r Remarks
    !r   xi is the radial part divided by r**l.
    !r   xi is evaluated by upward recursion for l>lmin+2.
    ! ---------------------------------------------------------------
    implicit none
    integer :: nrx,nr,lmin,lmax
    real(8) :: rsq(nrx),e,xi(nrx,lmin:lmax),rsm,wk(nr),wk2(nr)
    real(8) :: sre,r2,xx,ra,h0,arsm,earsm
    real(8) :: akap,a,r,um,up,x,facgl,facdu,dudr
    integer :: l,ir
    ! ... erfc(x) is evaluated as a ratio of polynomials,
    !     to a relative precision of <10^-15 for x<5.
    !     Different polynomials are used for x<1.3 and x>1.3.
    !     Numerators and denominators are t,b respectively.
    real(8) :: w,f1,f2, &
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
    ! ... f1(w=x-1/2) is erfc(x) for 0<x<1.3, if xx is y0*dexp(-x*x)
    f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)* &
         w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)* &
         w+b13)*w+b12)*w+b11)*w+1)
    ! ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
    f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)* &
         w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)* &
         w+b23)*w+b22)*w+b21)*w+1)
    ! --- Setup ---
    if (lmax < lmin .OR. nr == 0) return
    if (lmin < -1 .OR. lmin > 0) call rx('hansr2: bad lmin')
    a = 1/rsm
    akap = dsqrt(-e)
    arsm = akap*rsm/2
    earsm = dexp(-arsm**2)/2
    facgl = (2*a**2)*8*a*earsm
    ! ... uncomment the following for upward recursion from l=1...
    if (lmin == -1) facgl = facgl/(2*a**2)
    facdu = 8*a*earsm
    ! --- xi(*,lmin), xi(*,lmin+1) ---
    do  20  ir = 1, nr
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
20  enddo
    ! --- xi(ir,l) for l>1 by upward recursion ---
    facgl = 2*a**2
    do  30  l = lmin+2, lmax
       xx = 2*l-1
       do  7  ir = 1, nr
          xi(ir,l) = (xx*xi(ir,l-1) - e*xi(ir,l-2) - wk2(ir))/rsq(ir)
          wk2(ir) = facgl*wk2(ir)
7      enddo
30  enddo
  end subroutine hansr2
  subroutine hanr(rsq,lmin,lmax,nrx,nr,e,xi) !!     ... Asymtotic case, r>>rsm
    !- Vector of unsmoothed hankel functions for l=0...lmax, negative e.
    ! ---------------------------------------------------------------
    !i Inputs
    !i   rsq,nr  vector of points r**2, and number of points.
    !i   e       energy of Hankel.
    !i   lmin:lmax generate xi from lmin:lmax.  lmin must be -1 or 0.
    !i   nrx     leading dimension of xi.
    !o Outputs
    !o   xi      radial Hankel functions/r**l for: xi(1..nr, lmin:lmax)
    !o           Solid hankel function is hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
    !o           Energy derivative is    hlp(ilm) = x(l-1)/2*cy(ilm)*yl(ilm)
    ! ---------------------------------------------------------------
    implicit none
    integer :: nrx,nr,lmin,lmax,l,ir
    real(8) :: rsq(nr),e,xi(nrx,lmin:lmax),sre,r2,r,h0,xx,akap
    if (lmin /= 0 .AND. lmin /= -1) call rx('hanr: input lmin must be -1 or 0')
    if (lmax < lmin .OR. nr <= 0) return
    akap = dsqrt(-e)
    ! --- Make xi(lmin), xi(lmin+1) ---   ! ... xi(-1) for lmax=-1 only
    if (lmin == -1 .AND. lmax == -1) then
       do   ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,-1) = -h0*sre/e
       enddo
    elseif (lmin == -1) then
       do    ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,0) = h0
          xi(ir,-1) = -h0*sre/e
       enddo
       ! ... xi(0) for lmax=0 only
    elseif (lmax == 0) then
       do   ir = 1, nr
          r2 = rsq(ir)
          r = dsqrt(r2)
          xi(ir,0) = dexp(-akap*r)/r
       enddo
    else
       do   ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,0) = h0
          xi(ir,1) = h0*(1d0+sre)/rsq(ir)
       enddo
    endif
    ! --- xi(*,lmin+2:lmax) by upward recursion ---
    do  30  l = lmin+2, lmax
       xx = 2*l-1
       do  ir = 1, nr
          xi(ir,l) = (xx*xi(ir,l-1) - e*xi(ir,l-2))/rsq(ir)
       enddo
30  enddo
  end subroutine hanr
  subroutine hansmd(mode,r,e,rsm,lmax,hs,dhs,ddhs,hsp,dhsp,ddhsp) !Value and some derivatives of smoothed radial Hankel functions
    !i Inputs
    !i   mode  :tells hansmd what derivatives to make.
    !i         :1s digit concerns 2nd radial derivative
    !i         :0 make neither 1st or 2nd radial derivative.
    !i         :>0 make 1st and second radial derivative:
    !i         :1 ddhs = radial part of Laplacian, 1/r d^2 (r*h_l) / dr^2
    !i         :2 ddhs = 1/r d^2 (r*h_l) / dr^2  - l(l+1)/r^2 h_l
    !i         :  NB: ddhs = laplacian of 3-dimensional hs_l YL
    !i         :3 ddhs = d^2 (h_l) / dr^2
    !i         :1s digit concerns energy derivative
    !i         :0 make none of hsp,dhsp,ddhsp
    !i         :1 make all  of hsp,dhsp,ddhsp
    !i   r     :radius
    !i   e     :hankel energy
    !i   rsm   :hankel smoothing radius
    !i   lmax  :make function values for l between (0:lmax)
    !o  Outputs:
    !o   hs    : function values of the radial sm. hankel h(e,r,0:lmax)
    !o         : A solid hankel is H=h*Y_L, where Y_L are the spherical
    !o         : harmonics for unit radius (no scaling by r**l)
    !o   dhs   : radial derivative of hs
    !o   ddhs  : radial part of Laplacian of hs, i.e. 1/r d^2 (r h) /dr^2
    !o         : OR some other second derivative (see mode)
    !o   hsp   : energy derivative of hs
    !o   dhsp  : mixed energy + radial derivative of hs
    !o   ddhsp : 3-d laplacian of hsp YL
    !r Remarks
    !r  See J. Math. Phys. 39, 3393 (1998).
    !r    For radial derivative, see JMP 39, 3393, Eq. 4.7
    !r      h'  = l/r h_l - h_l+1
    !r    Second radial derivative:
    !r      h'' = l(l-1)/r^2 xi_l - (2l+1)/r xi_l+1 + xi_l+2
    !r      1/r d^2/dr^2 (r*h_l) = l(l+1)/r^2 h_l - (2l+3)/r h_l+1 + h_l+2
    !r    Energy derivative:  see JMP 39, 3393, Eq. 7.5
    !r      hp_l = r/2 h_l-1  Special case l=0: hp_0 = 1/2 h_-1 ?
    !r    Mixed energy + radial derivative:
    !r      hp'_l = h(l-1)*l/2 - h(l)*r/2
    !r    Mixed energy + kinetic energy
    !r      hp'' = -(2l+3)/2 h_l + h_l+1*r/2
    !r
    !r  Note connection with hansmr, which makes xi(l) = h(l) / r^l
    !u Updates
    !u   28 Aug 04 Also generate ddhsp
    !u   16 Jun 04 First created
    ! ---------------------------------------------------------------
    implicit none
    integer :: mode,lmax,idx,l,mode0,mode1
    real(8) :: r,e,rsm, hs(0:lmax),dhs(0:lmax),ddhs(0:lmax),hsp(0:lmax),dhsp(0:lmax),ddhsp(0:lmax), xi(-1:lmax+2)
    if (lmax < 0) return
    mode0 = mod(mode,10)
    mode1 = mod(mode/10,10)
    !call hansr(rsm,-1,lmax+2,1,lmax+2,e,r**2,1,1,idx,wk,11,xi)
    call hansr(rsm,-1,lmax+2,1,[lmax+2],[e],[r**2],1,1,[idx],11,xi)
    do  54  l = 0, lmax
       hs(l)   = xi(l)
       if (mode0 /= 0) then
          dhs(l)  = xi(l)*l/r - xi(l+1)
          if (mode0 == 1) ddhs(l) = xi(l)*l*(l+1)/r**2 - (2*l+3)/r*xi(l+1) + xi(l+2)
          if (mode0 == 2) ddhs(l) =                    - (2*l+3)/r*xi(l+1) + xi(l+2)
          if (mode0 == 3) ddhs(l) = xi(l)*l*(l-1)/r**2 - (2*l+1)/r*xi(l+1) + xi(l+2)
       endif
       if (mode1 /= 0) then
          hsp(l)   = xi(l-1)*r/2
          dhsp(l)  = (xi(l-1)*l - xi(l)*r)/2
          ddhsp(l) = - (2*l+3)*xi(l)/2 + xi(l+1)*r/2
       endif
54  enddo
    if(mode1 /= 0) hsp(0) = xi(-1)/2
  end subroutine hansmd
  subroutine hansmr(r,e,a,xi,lmax) !Smoothed hankel functions for l=0...lmax, negative e.
    !NOTE: Except numerical minor differences, hansr is equivalent to hansr in m_hansr.
    !  hansmr may have numerical problem for rms<1d-9 (tailsm.f90).
    !      hansmr is usef for core fitting parts, while hansr is for valence part.
    !      Probably,
    !      hansr: r is divided into three section. Less smoothness but numericall good overall
    !      hansmr: better smoothness but less numericall problematic probably when rsm<1d-9
    !      (Thus we need special treatements as in tailsm.f90). 2022-6-29
    ! ---------------------------------------------------------------
    !o  Outputs: xi(0:lmax)
    !o
    !r  xi is the radial part divided by r**l.
    !r  A solid smoothed hankel is from xi as:
    !r  hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
    !r
    !r  See J. Math. Phys. 39, 3393 (1998).
    !r  xi(l)= 2/sqrt(pi) * 2^l int_a^inf u^2l dexp(-r^2*u^2+kap2/4u^2) du
    implicit none
    integer :: lmax,l,n,nmax
    real(8) :: r,e,a,xi(0:lmax),a0(0:40),a2,add,akap,al,cc,dudr,ema2r2,fac, &
         gl,r2n,radd,rfac,rhs,rlim,srpi,sum,ta,ta2l,tol,u,uminus,uplus,w,derfc
    parameter (nmax=1000, tol=1d-20, srpi=1.77245385090551602729817d0)
    if (e > 0d0) call rx('hansmr: e gt 0')
    if (lmax > 40) call rx('hansmr: lmax gt 40')
    rlim = 1.5d0/a
    akap = dsqrt(dabs(e))
    ta = a+a
    a2 = a*a
    cc = 4d0*a2*a*dexp(e/(ta*ta))/srpi
    if (a*r > 10d0) then !call bessl if exp(-a*a*r*r) is very small ---
       call bessl(e*r*r,lmax,a0,xi)
       rfac = r
       do  l = 0, lmax
          rfac = rfac*(1d0/(r*r))
          xi(l) = rfac*xi(l)
       enddo
       return
    endif
    ! --- Power series for small r ---
    if (r > rlim) goto 90
    a0(0) = cc/(ta*a) - akap*derfc(akap/ta)
    rhs = cc
    fac = 1d0
    al = a0(0)
    do  l = 1, lmax
       al = -(e*al+rhs)/(2*l*(2*l+1))
       rhs = -rhs*a2/l
       fac = -2d0*fac*l
       a0(l) = fac*al
    enddo
    ta2l = 1d0
    do  200  l = 0, lmax
       rhs = cc*ta2l
       sum = a0(l)
       add = sum
       r2n = 1d0
       do n = 1, nmax
          add = -(e*add+rhs)/( 2*n*(2*n+(l+l+1)) )
          r2n = r2n*(r*r)
          radd = add*r2n
          sum = sum+radd
          if (dabs(radd) < tol) goto 22
          rhs = -rhs*a2/n
       enddo
       print *, 'hansmr (warning): power series did not converge'
22     continue
       xi(l) = sum
       ta2l = ta2l*(2d0*a2)
200 enddo
    return
    ! --- Big r: make xi0,xi1 explicitly; the higher l by recursion ---
90  continue
    ema2r2 = dexp(-a*a*r*r)
    uminus = derfc(akap/ta-r*a)*dexp(-akap*r)
    uplus = derfc(akap/ta+r*a)*dexp(+akap*r)
    u = .5d0*(uminus-uplus)
    w = -.5d0*(uminus+uplus)
    dudr = akap*w + ta*dexp(e/(ta*ta))*ema2r2/srpi
    xi(0) = u/r
    if (lmax >= 1) then
       xi(1) = (u/r - dudr)/(r*r)
       gl = cc*ema2r2
       if (lmax >= 2) then
          do  l = 2, lmax
             xi(l) = ((2*l-1)*xi(l-1) -e*xi(l-2) - gl)/(r*r)
             gl = 2d0*a2*gl
          enddo
       endif
    endif
  end subroutine hansmr
  subroutine corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc, rfoc,z) !Returns parameters for Zc part of Eq.(28) TK.JPSJ034702
    use m_lmfinit,only: pnux=>pnusp,pzx=>pzsp,sspec=>v_sspec,n0
    !i  is: species index
    !i       pnusp, pzsp
    !o Outputs
    !o   cofg  :coefficient to Gaussian part of pseudocore density assigned so that pseudocore charge = true core charge
    !o   cofh  :coefficient to smHankel part of pseudocore density. Hankel contribution is determined by inputs (qcorh,ceh,rfoc)
    !           and should accurately represent the true core density for r>rmt
    !o   qcorg :charge in the gaussian part
    !o   qcorh :charge in the Hankel part
    !o   qsc   :number of electrons in semicore treated by local orbitals
    !o   lfoc  :switch specifying treatment of core density.
    !o          0 => val,slo = 0 at sphere boundary
    !o          1 => core tails included explicitly with valence
    !o   rfoc :smoothing radius for hankel head fitted to core tail
    !o   z     :nuclear charge
    !r Remarks
    !r   qcorg and qcorh are the charges in the Gaussian and smHankels. The hankel part is used when the core is allowed to spill out of
    !r   the augmentation sphere.
    !r
    !r   cofg and cofh are the coefficients in front of the standard gaussian and smoothed hankel functions for l=0.
    !r   That is: the pseudocore density is 
    !r      cofg*g0(rg;r)*Y0 + cofh*h0(rfoca;r)*Y0        (1) (See 0th part Eq.(28).
    !r   ceh and rfoc are the energy and sm.-radius for the hankel part.
    !r   cofg is set so that qc = integral of eq. 1 above.
    !r
    !r   For lfoc=0 there is no Hankel part; qc carried entirely by Gausian
    !r   For lfoc>0 there is no Hankel part; Gaussian carries difference between qc and charge in Hankel part.
    !r
    !r   To add to the radial density 4*pi*r**2*rho_true, multiply cofg,cofh by srfpi.
    !l Local variables
    !l    ccof :coefficient for core tail, for a smoothed Hankel.
    !l          ccof is differs from spec->ctail because ctail is constructed for an unsmoothed Hankel.
    implicit none
    integer :: is,i_copy_size
    real(8):: qcorg , qcorh , qsc , cofg , cofh , ceh , rfoc , z
    integer:: lfoc , lmxb , l,isp
    real(8):: pnu(n0),pz(n0),ccof,q0,q1,qc,rmt,rsm,x0(0:n0), xi(0:n0),dgetss
    real(8),parameter:: fpi = 16d0*datan(1d0), srfpi = dsqrt(fpi), y0 = 1d0/srfpi
    lfoc=lfoca_i(is)
    rfoc=rfoca_i(is)
    lmxb=lmxb_i(is)
    z=   z_i(is)
    rmt = rmt_i(is)
    qc=  sspec(is)%qc
    ccof=sspec(is)%ctail
    ceh= sspec(is)%etail
    pnu= pnux(1:n0,1,is) 
    pz = pzx(1:n0,1,is)  
    if ( rfoc <= 1d-5 ) rfoc = rg_i(is)
    qsc = 0
    isp=1 !we assme int pz(:,1)=pz(:,2) int pnu as well
    do  l = 0, lmxb
       if (int(pz(l+1)) /= 0) then
          if (int(mod(pz(l+1),10d0)) < int(pnu(l+1))) qsc = qsc + 4*l+2
       endif
    enddo
    if (ccof /= 0) then ! ... Scale smoothed hankel coeff for exact spillout charge
       call hansmr(rmt,0d0,1/rfoc,x0,1)
       call hansmr(rmt,ceh,1/rfoc,xi,1)
       q1 = srfpi/ceh*(-dexp(rfoc**2/4*ceh) - rmt**3*(xi(1)-dexp(rfoc**2/4*ceh)*x0(1))) !q1 = spillout charge in sm. Hankel
       rsm = 0.05d0
       call hansmr(rmt,0d0,1/rsm,x0,1)
       call hansmr(rmt,ceh,1/rsm,xi,1)
       q0 = srfpi/ceh*(-dexp(rsm**2/4*ceh) - rmt**3*(xi(1)-dexp(rsm**2/4*ceh)*x0(1))) !q0 = spillout charge in ordinary Hankel
       q0 = q0*y0
       q1 = q1*y0
       ccof = ccof*q0/q1
    endif
    qcorg = qc
    qcorh = 0d0
    if (lfoc > 0) then
       qcorh = -ccof*dexp(ceh*rfoc*rfoc/4d0)/ceh ! Set gaussian and smhankel charges
       qcorg = qc-qcorh
    endif
    cofh = -y0*qcorh*ceh*dexp(-ceh*rfoc*rfoc/4d0)! Coeffients to the the gaussian and smhankel
    cofg = y0*qcorg
  end subroutine corprm
end module m_hansr
