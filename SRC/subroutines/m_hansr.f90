!>Smooth Hankel functions in real space. JMP39: https://doi.org/doi:10.1063/1.532437.
module m_hansmr 
  ! JMP39:
  ! Bott, E., M. Methfessel, W. Krabs, and P. C. Schmidt.
  ! “Nonsingular Hankel Functions as a New Basis for Electronic Structure Calculations.”
  ! Journal of Mathematical Physics 39, no. 6 (June 1, 1998): 3393–3425.
  ! https://doi.org/doi:10.1063/1.532437.
  public hansmr,hansmronly !hansmr is equilvaent to hansr except numerical accuracy problem. See note.
  logical:: hansmronly=.true. !new test 2023-10-31. If =.false., we use hansr instead of hansmr in places (recover Mark's original setting at 2009).
  private
contains
  subroutine hansmr(r,e,a,xi,lmax) !Smoothed hankel functions for l=0...lmax, negative e. a=1/rsm
    !o  Outputs: xi(0:lmax)
    !r   xi is the radial part divided by r**l.
    !r   A solid smoothed hankel is from xi as   hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
    !r   See J. Math. Phys. 39, 3393 (1998).
    !r    xi(l)= 2/sqrt(pi) * 2^l int_a^inf u^2l dexp(-r^2*u^2+kap2/4u^2) du
    !  hansmr is usef for core fitting parts, while hansr is for valence part.
    !      Probably,
    !      hansr: r is divided into three section. Less smoothness but numericall good overall
    !      hansmr: better smoothness but less numericall problematic probably when rsm<1d-9
    !      (Thus we need special treatements as in tailsm.f90). 2022-6-29
    !Except numerical minor differences, hansr is equivalent to hansr in m_hansr. See note in hansr

    !TK thinks that the current version of hansmr is equivalent with hansr? (2023-11-1)
    implicit none
    integer :: lmax,l,n,nmax
    real(8) :: r,e,a,xi(0:lmax),a0(0:40),a2,add,akap,al,cc,dudr,ema2r2,fac, &
         gl,r2n,radd,rfac,rhs,rlim,srpi,sum,ta,ta2l,tol,u,uminus,uplus,w,derfc
    parameter (nmax=1000, tol=1d-20, srpi=1.77245385090551602729817d0)
    if (e > 0d0  ) call rx('hansmr: e gt 0')
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
    elseif(r<=rlim) then    ! --- Power series for small r ---
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
       do  l = 0, lmax
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
22        continue
          xi(l) = sum
          ta2l = ta2l*(2d0*a2)
       enddo
    else     ! --- Big r: make xi0,xi1 explicitly; the higher l by recursion ---
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
    endif
  endsubroutine hansmr
endmodule m_hansmr
module m_hansr !hansmr is equilvaent to hansr except numerical accuracy problem. See note.
  ! JMP39:
  ! Bott, E., M. Methfessel, W. Krabs, and P. C. Schmidt.
  ! “Nonsingular Hankel Functions as a New Basis for Electronic Structure Calculations.”
  ! Journal of Mathematical Physics 39, no. 6 (June 1, 1998): 3393–3425.
  ! https://doi.org/doi:10.1063/1.532437.
  public hansr 
  logical:: hansmronly=.true. !new test 2023-10-31
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
    !r     rc1<r<rc2 are evaluated from error functions and the higher l's by upward recursion.
    !r     rc2<r are approximated with unsmoothed Hankels.
    !r   The relative error should be less than parameter 'tol', except a narrow region for r~rc1 and l>6,
    !r    where the precision degrades somewhat, worsening with higher l.  For all cases tested, the relative error continued to be ~<10^-13 for l<=9.
    implicit none
    integer :: nrx,nr,lmn,lmx,idx(nrx,2),nxi,lxi(nxi),job
    real(8) :: rsq(nr),e,exi(nxi)
    real(8) :: xi(nrx,lmn:lmx,nxi),wk(nrx,4+lmx-lmn)
    integer :: ir,j,l,k,n0,n1,n2,ie,lmax
    real(8) :: a,rsm,y0,a2,emin,tol
    real(8) :: rc1,rc2,akap,rl,rl0,rc20
    parameter (tol=1d-15)
    logical :: ltmp,lsort,lscal
    call rx('unusedxxxxxxxxxxxxx')
    lscal = mod(job,10) .ne. 0
    lsort = mod(job/10,10) .ne. 0
    ! --- Check lmx; handle case rsm=0 ---
    ltmp = rsm .lt. 1d-12
    lmax = -1
    do ie = 1, nxi
       if(lxi(ie) > lmx) call rx('hansr: lxi gt lmx')
       if(exi(ie) > 0)   call rx('hansr: exi gt 0')
       if(ltmp) call hanr(rsq,lmn,lxi(ie),nrx,nr,exi(ie),xi(1,lmn,ie))
       lmax = max(lmax,lxi(ie))
    enddo
    if (ltmp) goto 60
    emin = minval(exi(1:nxi)) ! --- Find cutoffs rc2 (negligible smoothing) and rc1 (power series) ---    
    akap = dsqrt(-emin)
    a = 1/rsm
    a2 = a*a
    y0 = 1/dsqrt(16*datan(1d0))
    rc20 = akap/(2*a) ! ... For r>rc2 approximate smooth Hankels with normal ones
    rc2 = ((rc20 + dsqrt(rc20**2 - dlog(tol)))/a)**2
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
    do 40  ie = 1, nxi
       e = exi(ie)
       akap = dsqrt(-e)
       lmax = lxi(ie)
       if (lsort) then !   ... Case calculate points in original order (already sorted)
          call hansr1(rsq(1),lmn,lmax,nrx,n1-1,e,rsm,dsqrt(rc1), xi(1,lmn,ie)) !     ... Power series for points within rc1
          if(n1 <= nr) call hansr2(rsq(n1),lmn,lmax,nrx,n2-n1,e,rsm,wk(n1,1),wk(n1,2),xi(n1,lmn,ie)) !     ... Normal evaluation of smoothed Hankels
          if(n2 <= nr) call hanr(rsq(n2),lmn,lmax,nrx,nr+1-n2,e,xi(n2,lmn,ie)) !     ... Asymtotic case, r>>rsm
       else !   ... Case calculated points in sorted
          call hansr1(wk(1,3),lmn,lmax,nrx,n1-1,e,rsm,dsqrt(rc1), wk(1,4)) !     ... Power series for points within rc1
          call hansr2(wk(n1,3),lmn,lmax,nrx,n2-n1,e,rsm,wk(n1,1), wk(n1,2),wk(n1,4)) !     ... Normal evaluation of smoothed Hankels
          if(n2 <= nr) call hanr(wk(n2,3),lmn,lmax,nrx,nr+1-n2,e,wk(n2,4)) !     ... Asymtotic case, r>>rsm
          do   l = lmn, lmax
             do    ir = 1, nr
                j = idx(ir,1)
                xi(ir,l,ie) = wk(idx(ir,1),4+l-lmn) !     ... Poke into xi(lmn:lmax), with the original ordering of points
             enddo
          enddo
       endif
40  enddo
60  continue
    if(.NOT.lscal) return
    do    ir = 1, nr ! --- Scale by r**l if job nonzero ---
       do ie = 1, nxi
          xi(ir,1:lxi(ie),ie) = [(xi(ir,l,ie)*(rsq(ir)**.5)**l,l=1,lxi(ie))]
       enddo
    enddo
  endsubroutine hansr
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
    integer:: nrx,nr,lmin,lmax,l,ir
    real(8):: rsq(nrx),e,xi(nrx,lmin:lmax),rsm,wk(nr),wk2(nr),sre,r2,xx,ra,arsm,earsm,akap,a,r,um,up,x,facgl,facdu,dudr,erfcee
    if(lmax < lmin .OR. nr == 0) return
    if(lmin < -1 .OR. lmin > 0) call rx('hansr2: bad lmin')
    a = 1/rsm
    akap = dsqrt(-e)
    arsm = akap*rsm/2
    earsm = dexp(-arsm**2)/2
    facgl = (2*a**2)*8*a*earsm
    if (lmin == -1) facgl = facgl/(2*a**2)
    facdu = 8*a*earsm
    do  ir = 1, nr !xi(*,lmin), xi(*,lmin+1) ---
       r2 = rsq(ir)
       r = dsqrt(r2)
       ra = r*a
       sre = akap*r
       xx = earsm*wk(ir)/r
       x = ra-arsm
       if(x>0 ) um=dexp(-sre)/r-xx*erfcee(x) ! ---   Evaluate um,up ---
       if(x<=0) um=xx*erfcee(x)
       up= xx*erfcee(ra + arsm) !assumes x gt 0
       xi(ir,0) = um - up !   ... xi(0) = um - up
       if (lmin == -1) then !   ... xi(-1) = (um + up)*r/akap
          xi(ir,-1) = (um + up)*r/akap
       elseif (lmax >= 1) then !   ... xi(1)
          dudr = facdu*wk(ir) - sre*(um+up)
          xi(ir,1) = (xi(ir,0) - dudr)/r2
       endif
       wk2(ir) = facgl*wk(ir)
    enddo
    facgl = 2*a**2
    do l = lmin+2, lmax !xi(ir,l) for l>1 by upward recursion ---
       xx = 2*l-1
       xi(1:nr,l) = (xx*xi(1:nr,l-1) - e*xi(1:nr,l-2) - wk2(1:nr))/rsq(1:nr)
       wk2(1:nr) = facgl*wk2(1:nr)
    enddo
  end subroutine hansr2
  subroutine hanr(rsq,lmin,lmax,nrx,nr,e,xi) !!     ... Asymtotic case, r>>rsm - Vector of unsmoothed hankel functions for l=0...lmax, negative e.
    !i Inputs
    !i   rsq,nr  vector of points r**2, and number of points.
    !i   e       energy of Hankel.
    !i   lmin:lmax generate xi from lmin:lmax.  lmin must be -1 or 0.
    !i   nrx     leading dimension of xi.
    !o Outputs
    !o   xi      radial Hankel functions/r**l for: xi(1..nr, lmin:lmax)
    !o           Solid hankel function is hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
    !o           Energy derivative is    hlp(ilm) = x(l-1)/2*cy(ilm)*yl(ilm)
    implicit none
    integer :: nrx,nr,lmin,lmax,l,ir
    real(8) :: rsq(nr),e,xi(nrx,lmin:lmax),sre,r2,r,h0,xx,akap,srer(nr),rr(nr)
    if(lmin/=0 .AND.lmin/= -1) call rx('hanr: input lmin must be -1 or 0')
    if(lmax<lmin.OR.nr<= 0) return
    rr(1:nr)= dsqrt(rsq(1:nr))
    srer(1:nr)=dsqrt(-e)*rr(1:nr)
    if(lmin<=-1.and.-1<=lmax) xi(1:nr,-1)= -dexp(-srer)/rr*srer/e
    if(lmin<=0 .and. 0<=lmax) xi(1:nr,0) =  dexp(-srer)/rr
    if(lmin<=1 .and. 1<=lmax) xi(1:nr,1) =  dexp(-srer)/rr*(1d0+srer)/rsq
    do l = lmin+2, lmax !xi(*,lmin+2:lmax) by upward recursion ---
       xi(1:nr,l) = ((2*l-1)*xi(1:nr,l-1) - e*xi(1:nr,l-2))/rsq(1:nr)
    enddo
  endsubroutine hanr
endmodule m_hansr

subroutine corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc, rfoc,z) !Returns parameters for Zc part of Eq.(28) TK.JPSJ034702
  use m_lmfinit,only: pnux=>pnusp,pzx=>pzsp,n0
  use m_lmfinit,only: z_i=>z,rmt_i=>rmt,lmxb_i=>lmxb,lfoca_i=>lfoca,rfoca_i=>rfoca,rg_i=>rg
  use m_fatom,only:sspec
  use m_hansmr,only:hansmr
  !i  is: species index
  !o Outputs
  !o   cofg  :coefficient to Gaussian part of pseudocore density assigned so that pseudocore charge = true core charge
  !o   cofh  :coefficient to smHankel part for n^c_sH,a (See Eq.. Hankel contribution is determined by inputs (qcorh,ceh,rfoc)
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
  !r   cofh is the coefficients as 
  !      n_sH,a = cofh*h0(rfoca;r)*Y0 !Eq(23)
  !    cofg = y0* \int_0^rmt (core(r) - n_sH,a(r)) dr = y0*(qc-qcorh) because core(r) and n_sH,a(r) are the same for r>rmt
  !
  !r   ceh and rfoc are the energy and sm.-radius for the hankel part.
  !r   cofg is set so that qc = integral of eq. 1 above.
  !r
  !r   For lfoc=0 there is no Hankel part; qc carried entirely by Gausian
  !r   For lfoc>0 there is no Hankel part; Gaussian carries difference between qc and charge in Hankel part.
  !r
  !r   To add to the radial density 4*pi*r**2*rho_true, multiply cofg,cofh by srfpi.
  !
  !l Local variables
  implicit none
  integer :: is,i_copy_size
  real(8):: qcorg , qcorh , qsc , cofg , cofh , ceh , rfoc , z
  integer:: lfoc , lmxb , l,isp
  real(8):: pnu(n0),pz(n0),ccof,q0,q1,qc,rmt,rsm,x0(0:n0), xi(0:n0),dgetss
  real(8),parameter:: fpi = 16d0*datan(1d0), srfpi = dsqrt(fpi), y0 = 1d0/srfpi
  lfoc=lfoca_i(is)
  rfoc=rfoca_i(is)
  lmxb=lmxb_i(is)
  z =   z_i(is)
  rmt = rmt_i(is)
  qc  = sspec(is)%qc
  ccof= sspec(is)%ctail
  ceh=  sspec(is)%etail
  pnu= pnux(1:n0,1,is) 
  pz = pzx(1:n0,1,is)  
  if ( rfoc <= 1d-5 ) rfoc = rg_i(is)  !we assme int pz(:,1)=pz(:,2) int pnu as well
  qsc = sum([(4*l+2,l=0,lmxb)], mask=[(pz(l+1)>0d0.and.floor(mod(pz(l+1),10d0))<floor(pnu(l+1)),l=0,lmxb)])
  if(ccof /= 0) then ! ... Correct smHankel coefficient ccof to reproduce exact spillout charge
     !       ccof differs from spec->ctail because ctail is constructed for an unsmoothed Hankel.
     call hansmr(rmt,0d0,1/rfoc,x0,1)
     call hansmr(rmt,ceh,1/rfoc,xi,1)
     q1 = srfpi/ceh*(-dexp(rfoc**2/4*ceh) - rmt**3*(xi(1)-dexp(rfoc**2/4*ceh)*x0(1))) !q1 = spillout charge in smHankel
     rsm = 0.05d0
     call hansmr(rmt,0d0,1/rsm,x0,1)
     call hansmr(rmt,ceh,1/rsm,xi,1)
     q0 = srfpi/ceh*(-dexp(rsm**2/4*ceh) - rmt**3*(xi(1)-dexp(rsm**2/4*ceh)*x0(1))) !q0 = spillout charge in ordinary Hankel
     q0 = q0*y0
     q1 = q1*y0
     ccof = ccof*q0/q1
  endif
  qcorh = merge(-ccof*dexp(ceh*rfoc*rfoc/4d0)/ceh, 0d0,lfoc>0) ! hankel charge
  qcorg = merge(qc-qcorh, qc,lfoc>0)                           ! counter charge  
  cofh = -y0*qcorh*ceh*dexp(-ceh*rfoc*rfoc/4d0)! Coeffients of the smHankel to reproduce cores.
  cofg =  y0*qcorg                             ! charge for Y0 
endsubroutine corprm
