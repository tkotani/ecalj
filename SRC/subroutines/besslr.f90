subroutine bessl(y,lmax,fi,gi)! Spherical Bessel Neumann Hankel functions
  !i Inputs:   
  !i   y     :y = r**2 or (i*r)**2= -r**2  !i=cmplex(0,1d0)
  !i   lmin  :minimum l
  !i   lmax  :maximum l
  ! TKotani think output for give y is (We use Methfessel's conventions)
  !     fi_l =   1/r**l     j_l(r) 
  !     gi_l =   1/r**(l+1) n_l(r) for y>0
  !     gi_l =  -1/r**(l+1) h^1_l(i*r) for y<0 !(h^1 means 1st kind)
  ! ,where j_l,n_l,h^1_l are standard sperical bessel functions.
  ! NOTE:   For r > 40 this algorithm is numerically unstable !!!
  ! xxx Check this and through away following old memos. !2023-jan
  !
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
  !r  *Special cases:
  !r     fi(y,l) -> 1/(2l+1)!!     for y->0
  !r     gi(y,l) -> (2l-1)!!       for y->0
  !r     gi(y,l=0) = exp(i*akap*r) for y<=0    akap=sqrt(e) Im(akap)>=0
  !r     gi(y,l=0) = cos(akap*r)   for y>0     akap=sqrt(e) Im(akap)=0
  !b   For x > 40 this algorithm is numerically unstable !!!
  !u Updates
  !u   23 Jul 08 bug fix: besslr doesn't make fi/gi(lmax+1) when lmax=0
  !u   19 May 04 Changed loka from logical to integer
  ! ----------------------------------------------------------------------
  implicit none
  integer,parameter :: lmin=0
  integer :: lmax
  double precision :: y,fi(lmin:lmax),gi(lmin:lmax)
  integer :: i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,nf,tlp1,ll1,ll2,nlmax
  parameter (nlmax=20)
  double precision :: dt,dt2,exppr,my,srmy,g1,t,tol, &
       dum(nlmax*4+2),fac2l(-nlmax:nlmax*2+3)
  parameter(tol=1.d-15)
  lmx = max0(lmax,2)
  if (lmx > nlmax+nlmax)call rxi(' BESSL : lmax gt nlmax*2, lmax=',lmx)
  ! --- A table of fac2l(l)=(2l-1)!!
  !     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
  fac2l(0) = 1d0
  do  10  l = 1, lmx+1
     fac2l(l) = fac2l(l-1) * (l+l-1)
10 enddo
  do  11  l = -1, lmin, -1
     fac2l(l) = fac2l(l+1) / (l+l+1)
11 enddo
  ! --- Case akap=0 ---
  if (y == 0) then
     do  12  l = lmin, lmax
        fi(l) = 1/fac2l(l+1)
        gi(l) = fac2l(l)
12   enddo
     goto 100
  endif
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
  call rx('BESSL: series not convergent')
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
  call rx('BESSL: series not convergent')
31 continue
  dum(2) = t/fac2l(lmx)
  ! --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
  ll1 = lmx + lmx + 1
  ll2 = ll1 + 1
  nf = ll1
  do  k = 3, ll2
     nf = nf-2
     dum(k) = nf*dum(k-1) - y*dum(k-2)
  enddo
  lmxp1 = lmx+1
  lmxp2 = lmx+2
  isn = (-1)**lmin
  do   k = lmin, lmax ! --- Get fi and gi from dum ---
     j1 = lmxp1-k
     j2 = lmxp2+k
     fi(k) = dum(j1)
     gi(k) = dum(j2)*isn !  ... n_l(x) = j_{-l-1}*(-1)^{l+1} !Neumann function
     isn = -isn
  enddo
! --- For E<0, use Hankel functions rather than Neumann functions ---
  if ( y < 0d0) then
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
     do  l = lmin, lmax
        gi(l) = gi(l)*exppr
     enddo
  endif
100 continue
end subroutine bessl
