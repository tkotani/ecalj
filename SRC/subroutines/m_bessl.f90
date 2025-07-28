module m_bessl
  public :: bessl, wronkj, radkj
  private
  !C A table of (2l-1)!!
  real(8), parameter :: fac2l(0:24) = [ 0.1000000000000000D+01, & !l=0
    0.1000000000000000D+01, 0.3000000000000000D+01, 0.1500000000000000D+02, 0.1050000000000000D+03, & !l=1-4
    0.9450000000000000D+03, 0.1039500000000000D+05, 0.1351350000000000D+06, 0.2027025000000000D+07, & !l=5-8
    0.3445942500000000D+08, 0.6547290750000000D+09, 0.1374931057500000D+11, 0.3162341432250000D+12, & !l=9-12
    0.7905853580625000D+13, 0.2134580466768750D+15, 0.6190283353629375D+16, 0.1918987839625106D+18, & !l=13-16
    0.6332659870762850D+19, 0.2216430954766998D+21, 0.8200794532637891D+22, 0.3198309867728778D+24, & !l=17-20
    0.1311307045768799D+26, 0.5638620296805834D+27, 0.2537379133562626D+29, 0.1192568192774434D+31]   !l=21-24
  !$acc declare copyin(fac2l)
contains
subroutine bessl(y,lmax,fi,gi)! Spherical Bessel  and  Neumann Hankel functions
  !$acc routine seq
  !i Inputs:   
  !i   y     :y = r**2 or (i*r)**2= -r**2  !i=cmplex(0,1d0)
  !i   lmax  :maximum l
  ! TKotani think output for give y is (We use Methfessel's conventions)
  !     fi_l =   1/r**l     j_l(r) 
  !     gi_l =   1/r**(l+1) n_l(r)     for y>0
  !     gi_l =  -1/r**(l+1) h^1_l(i*r) for y<0 !(h^1 means 1st kind)
  ! ,where j_l, n_l, h^1_l are standard sperical bessel functions.
  !r  *Special cases:
  !r     fi(y,l) -> 1/(2l+1)!!     for y->0
  !r     gi(y,l) -> (2l-1)!!       for y->0
  !r     gi(y,l=0) = exp(i*akap*r) for y<=0    akap=sqrt(e) Im(akap)>=0
  !r     gi(y,l=0) = cos(akap*r)   for y>0     akap=sqrt(e) Im(akap)=0
  !b   For x > 40 this algorithm is numerically unstable !!!
  ! NOTE:   For r > 40 this algorithm is numerically unstable !!!
  ! xxx Check this and through away following old memos below. !2023-jan
  !
  implicit none
  integer :: lmax, i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,nf,tlp1,ll1,ll2,nlmax
  real(8) :: y,fi(0:lmax),gi(0:lmax), dt,dt2,exppr,my,srmy,g1,t
  ! --- A table of fac2l(l)=(2l-1)!!  data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
  real(8) :: dum(lmax*2+2)
  real(8),parameter:: tol=1d-15
  if (y == 0) then ! --- Case akap=0 ---
     do l=0,lmax
        fi(l) = 1/fac2l(l+1)
        gi(l) = fac2l(l)
     enddo   
     return
  endif
  lmx = max(lmax,2)
  my = -y
  ! --- Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
  tlp1 = lmx+lmx+1
  dt = 1d0
  t = 1d0
  i = 0
  do  20  k = 1, 1000
     if (dabs(dt) < tol) goto 21
     i = i+2
     dt2 = i+tlp1
     dt = dt*my/(i*dt2)
     t = t+dt
20 enddo
  ! call rx('BESSL: series not convergent 111')
21 continue
  dum(1) = t/fac2l(lmx+1)
  
  ! --- Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1) ---
  tlp1 =  tlp1-2
  dt = 1d0
  t = 1d0
  i = 0
  do  30  k = 1, 1000
     if (dabs(dt) < tol) goto 31
     i = i+2
     dt2 = i+tlp1
     dt = dt*my/(i*dt2)
     t = t+dt
30 enddo
  ! call rx('BESSL: series not convergent 222')
31 continue
  dum(2) = t/fac2l(lmx)
  
  ! --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
  nf = 2* lmx + 1
  do  k = 3, 2*lmx + 2
     nf = nf-2
     dum(k) = nf*dum(k-1) - y*dum(k-2)
  enddo
  do  k = 0, lmax ! --- Get fi and gi from dum ---
    fi(k) = dum(lmx+1-k)
  enddo
 
  !--- gi parts --------------------------   
  if( y>=0d0 )then
    do k = 0, lmax ! --- Get fi and gi from dum ---
      gi(k) = dum(lmx+2+k)*(-1)**k !  ... n_l(x) = j_{-l-1}*(-1)^{l+1} !Neumann function
    enddo
  else !For E<0, use Hankel functions rather than Neumann functions ---
     srmy = dsqrt(-y)
     gi(0) = 1d0
     g1 = 1d0+srmy
     if (lmax >= 1) gi(1) = g1
     if (lmax >= 2) then
        tlp1 = 1
        do   l = 2, lmax
           tlp1 = tlp1+2
           gi(l) = tlp1*gi(l-1) - y*gi(l-2)
        enddo
     endif
     exppr = 1d0/dexp(srmy)
     do  l = 0, lmax
        gi(l) = gi(l)*exppr
     enddo
  endif
end subroutine bessl

subroutine wronkj(e1,e2,r,lmax,fkk,fkj,fjk,fjj)! Wronskians for hankels and bessels on one sphere.
  !$acc routine seq
  !  fxy is continuous as e1-->e1.
  !  fkk,fjj are symmetric in e1,e2.  fkj(e1,e2)=fjk(e2,e1).
  implicit real*8 (a-h,p-z), integer(o)
  integer:: lmax,lp1,l
  dimension fkk(*),fkj(*),fjk(*),fjj(*),ak1(200),aj1(200), &
  ak2(200),aj2(200),dk2(200),dj2(200),dk1(200),dj1(200)
  ! ------ first: special case e1=e2=0 -------------
  if(dabs(e1) <= 1.d-6 .AND. dabs(e2) <= 1.d-6) then
     r3=r*r*r
     rk=-1.d0
     rj=1.d0/r
     do 20 l=0,lmax
        lp1=l+1
        rk=rk*(2*l-1)/r
        rj=rj*r/(2*l+1)
        fkk(lp1)=rk*rk*r3/(2*l-1)
        fjj(lp1)=-rj*rj*r3/(2*l+3)
        fkj(lp1)=-0.5d0*rj*rk*r3
        fjk(lp1)=fkj(lp1)
20   enddo
     return
  endif
  ! -------- case e1.ne.e2 --------------------
  if(dabs(e1-e2) > 1.d-6) then
     efac=1.d0/(e2-e1)
     call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,0)
     call radkj(e2,r,lmax,ak2,aj2,dk2,dj2,0)
     do 10 lp1=1,lmax+1
        fkk(lp1)=efac*r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
        fjj(lp1)=efac*r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
        fkj(lp1)=efac*r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))-efac
        fjk(lp1)=efac*r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))+efac
10   enddo
  else
     ! --------- case e1.eq.e2 but not zero ---------
     call radkj(e1,r,lmax,ak1,aj1,dk1,dj1,0)
     call radkj(e1,r,lmax,ak2,aj2,dk2,dj2,1)
     do 11 lp1=1,lmax+1
        fkk(lp1)=r*r*(ak1(lp1)*dk2(lp1)-dk1(lp1)*ak2(lp1))
        fjj(lp1)=r*r*(aj1(lp1)*dj2(lp1)-dj1(lp1)*aj2(lp1))
        fkj(lp1)=r*r*(ak1(lp1)*dj2(lp1)-dk1(lp1)*aj2(lp1))
        fjk(lp1)=r*r*(aj1(lp1)*dk2(lp1)-dj1(lp1)*ak2(lp1))
11   enddo
  endif
  return
end subroutine wronkj
! --------------------------------------------------------
subroutine radkj(e,r,lmax,ak,aj,dk,dj,job)
  !$acc routine seq
  !  radial parts of spherical hankels and bessels.
  !  job=0: makes values and slopes. job=1: makes energy derivatives.
  implicit real*8 (a-h,p-z), integer(o)
  integer:: lmax,l,lp1,job
  dimension ak(*),aj(*),dk(*),dj(*),phi(200),psi(200),php(200),psp(200) ! MIZUHO-IR
  er2=e*r*r
  if(job == 0) then
     call bessl(er2,lmax+1,phi,psi)
     rl=1.d0/r
     do 10 l=0,lmax
        lp1=l+1
        rl=rl*r
        ak(lp1)=psi(lp1)/(rl*r)
        aj(lp1)=phi(lp1)*rl
        dk(lp1)=(l*psi(lp1)-psi(l+2))/(rl*r*r)
        dj(lp1)=(l*phi(lp1)-er2*phi(l+2))*rl/r
10   enddo
  else
     call bessl(er2,lmax+2,phi,psi)
     do 11 lp1=1,lmax+2
        php(lp1)=-0.5d0*r*r*phi(lp1+1)
        psp(lp1)=0.5d0*((2*lp1-1)*psi(lp1)-psi(lp1+1))/e
11   enddo
     rl=1.d0/r
     do 12 l=0,lmax
        lp1=l+1
        rl=rl*r
        ak(lp1)=psp(lp1)/(rl*r)
        aj(lp1)=php(lp1)*rl
        dk(lp1)=(l*psp(lp1)-psp(l+2))/(rl*r*r)
        dj(lp1)=(l*php(lp1)-er2*php(l+2)-r*r*phi(l+2))*rl/r
12   enddo
  endif
  return
end subroutine radkj
! real(8) function fac2l(i)
!   !$acc routine seq
!   !C A table of (2l-1)!!
!   !     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
!   integer:: i
!   fac2l = fac2lx(i)
! end function fac2l
end module m_bessl
  ! Followings are old memo: Removes this old memo when you are definite for the definition above (or corrected).
  ! ----------------
  !i Inputs: We use Methfessel's conventions
  !i   y     :y = e*r**2 = z**2 = (kappa*r)**2 = -(akap*r)**2, Im kappa>=0
  !o Outputs:
  !o   fi    :proportional to (Bessel function) / r**l.  Constant
  !o         :of proportionality depends on conventions; see Remarks.
  !r         :fi is evaluated in a power series expansion.
  !o   gi    :proportional to :
  !o         :Neumann function * r**(l+1) for y>0
  !o         :Hankel function  * r**(l+1) for y<0.
  !
  !o         :See Remarks for constant of proportionality.
  !o         :For y>0 gi is evaluated in a power series expansion.
  !o         :gi(l=0) = cos(sqrt(y))
  !o         :For y<0:
  !o         :gi(l=0) = exp(-sqrt(-y)) = exp(- kappa r)
  !o         :gi(l=1) = (1+sqrt(-y))*exp(-sqrt(-y)) = (1-i kappa r)gi(l=0)
  !o         :gi(l)   = (2l+1)*gi(l-1) - y*gi(l-2)
  !r Remarks:
  !r  *Bessel and Hankel functions from fi and gi.
  !r   Let j, n, h be spherical Bessel, Neumann and Hankel functions.
  !r   We use Im kappa > 0.
  !r   Conventions for Hankel functions vary somewhat and are defined below.
  !r
  !r  *Jackson, and Morse and Feshback conventions: z = kappa * r
  !r     h = j + i*n,  h = Hankel function of the first kind h^1, and
  !r     h^1,2(z) = j(z) +/- i n(z)
  !r     h^1_0(z) = exp(iz)/iz
  !r     h^1_1(z) = exp(iz)/iz * (-1 - i/z)
  !r   i h^1_l(z) = gi_l(z**2) / z^(l+1) (computed if y<0)
  !r              ->(2l-1)!! / z^(l+1) for z -> 0
  !r       j_0(z) = sin(z)/z
  !r       j_l(z) = fi_l(z**2) * z^l
  !r       n_0(z) = -cos(z)/z
  !r     Limiting cases for z = kappa*r -> 0:
  !r       n_l(z) ->  -(2l-1)!!/z^(l+1) ( 1 - z^2/(2(2l-1) + ...)
  !r       j_l(z) ->   z^l / (2l+1)!!   ( 1 - z^2/(2(2l+3) + ...)
  !r       h^1_l(z) -> -i (2l-1)!!/z^(l+1)
  !r     Limiting cases for z = kappa*r -> infty:
  !r       h^1_l(z) ->  exp(i(z-l*pi/2))/iz
  !r       n_l(z)   -> -cos(z-l*pi/2)/z
  !r       j_l(z)   ->  sin(z-l*pi/2)/z
  !r
  !r  *Methfessel's conventions (label as,an,ak,aj) have the properties:
  !
  !r     (1)  functions ak and aj, and an are real for e<=0
  !r     (2)  cases e .ne. 0 and e .eq. 0 same equations
  !r     Define akap = -i kappa = sqrt(-e):  akap is real and >0 for e<0.
  !r       ak_0 = exp(-akap r)/r;  ak_1 = (1 + akap*r) ak_0 / r
  !r     Relation to standard conventions:
  !r       ak_l = -i kappa^(l+1) h_l = -i (-i akap)^(l+1) h_l
  !r            = gi_l / r^(l+1)
  !r       aj_l = 1/kappa^l j_l      = 1/(-i akap)^l j_l
  !r            = fi_l r^l
  !r     Define Neumann function as:
  !r       an = ah - i kap e^l aj = ah + akap e^l aj
  !r     Solid Hankels, Bessels are defined as (CHECK)
  !r       H_L = Y_L(-grad) ak(r);   J_L = E^(-l) Y_L (-grad) aj(r)
  !r     Limiting cases for energy->0:
  !r       ak -> (2l-1)!!/r^(l+1)      aj -> r^l/(2l+1)!!
  !r
  !r  *Generation of fi and gi:  fi (and gi for y>0) are calculated for
  !r   lmax and lmax-1 by an expansion in powers of x^2=y:
  !r
  !r                     (-x^2/2)^k
  !r       fi =  Sum_k  --------------  = dum(lmx+1-l)
  !r                    k! (2l+1+2k)!!
  !r
  !r       gi = dum(lmx+2+l), y>0
  !r
  !r     For y<0:
  !r       gi(l=0) = exp(-sqrt(-y))
  !r       gi(l=1) = (1d0+sqrt(-y))*exp(-sqrt(-y))
  !r       gi(l)   = (2l+1)*gi(l-1) - y*gi(l-2)
  !r
  !r     The remaining functions are obtained by the recursion formula:
  !r
  !r       j_{l-2}(x)=(2l-1)j_{l-1}(x)/x-j{l}(x)   l=lmx-2,lmx-3,..,-lmx-1
  !r
  !r     => dum(k) = (2*lmx+5-2*k)*dum(k-1)-y*dum(k-2)   k=3,4,...,2*lmx+2
  !r
  !r     and the Neumann function are given through:
  !r
  !r           n_l(x) = j_{-l-1}*(-1)^{l+1}
  !r
  !u Updates
  !u   23 Jul 08 bug fix: besslr doesn't make fi/gi(lmax+1) when lmax=0
  !u   19 May 04 Changed loka from logical to integer
  ! ----------------------------------------------------------------------
