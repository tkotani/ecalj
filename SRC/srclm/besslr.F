      subroutine besslr(y,loka,lmin,lmax,fi,gi)
C- Radial part of Bessel functions, standard definitions
C ----------------------------------------------------------------------
Ci Inputs:
Ci   y     :y = e*r**2 = z**2 = (kappa*r)**2 = -(akap*r)**2, Im kappa>=0
Ci   loka  :0 Methfessel's conventions
Ci         :1 Andersen conventions from 2nd generation LMTO
Ci   lmin  :minimum l
Ci   lmax  :maximum l
Co Outputs:
Co   fi    :proportional to (Bessel function) / r**l.  Constant
Co         :of proportionality depends on conventions; see Remarks.
Cr         :fi is evaluated in a power series expansion.
Co   gi    :proportional to :
Co         :Neumann function * r**(l+1) for y>0
Co         :Hankel function  * r**(l+1) for y<0.
Co         :See Remarks for constant of proportionality.
Co         :For y>0 gi is evaluated in a power series expansion.
Co         :gi(l=0) = cos(sqrt(y))
Co         :For y<0:
Co         :gi(l=0) = exp(-sqrt(-y)) = exp(i kappa r)
Co         :gi(l=1) = (1+sqrt(-y))*exp(-sqrt(-y)) = (1-i kappa r)gi(l=0)
Co         :gi(l)   = (2l+1)*gi(l-1) - y*gi(l-2)
Cr Remarks:
Cr  *Bessel and Hankel functions from fi and gi.
Cr   Let j, n, h be spherical Bessel, Neumann and Hankel functions.
Cr   We use Im kappa > 0.
Cr   Conventions for Hankel functions vary somewhat and are defined below.
Cr   See radhjz for generation of j,n,h following different conventions.
Cr
Cr  *Jackson, and Morse and Feshback conventions: z = kappa * r
Cr     h = j + i*n,  h = Hankel function of the first kind h^1, and
Cr     h^1,2(z) = j(z) +/- i n(z)
Cr     h^1_0(z) = exp(iz)/iz
Cr     h^1_1(z) = exp(iz)/iz * (-1 - i/z)
Cr   i h^1_l(z) = gi_l(z**2) / z^(l+1) (computed if y<0)
Cr              ->(2l-1)!! / z^(l+1) for z -> 0
Cr       j_0(z) = sin(z)/z
Cr       j_l(z) = fi_l(z**2) * z^l
Cr       n_0(z) = -cos(z)/z
Cr     Limiting cases for z = kappa*r -> 0:
Cr       n_l(z) ->  -(2l-1)!!/z^(l+1) ( 1 - z^2/(2(2l-1) + ...)
Cr       j_l(z) ->   z^l / (2l+1)!!   ( 1 - z^2/(2(2l+3) + ...)
Cr       h^1_l(z) -> -i (2l-1)!!/z^(l+1)
Cr     Limiting cases for z = kappa*r -> infty:
Cr       h^1_l(z) ->  exp(i(z-l*pi/2))/iz
Cr       n_l(z)   -> -cos(z-l*pi/2)/z
Cr       j_l(z)   ->  sin(z-l*pi/2)/z
Cr
Cr  *Gunnarsson's conventions (PRB 27, 7144, 1983).
Cr   Somewhat confusing.  Gunnarsson chooses  k = -kappa, Im k <= 0.
Cr   Apparently his definitions correspond to -1 * standard h^1.
Cr     h_l = j_l(k r) - i n_l(k r)     see eqn prior to eq. A8
Cr   Taking his A15 as a definition deriving from Andersen, we have
Cr   i h_l = -gi(z**2) / (kr)^(l+1)
Cr         -> - (2l-1)!! / (kr)^(l+1) for k->0  A10 (sign error fixed)
Cr     h_0 = -exp(-ikr)/ikr
Cr     h_1 = h_0 (1 - ikr) / kr
Cr     j_0 = sin(kr)/kr
Cr     Functions of odd order are imaginary for e<0
Cr     h_l(kr) = h^1_l(-kr)
Cr             = i gi_l((kr)**2) / (kr)^(l+1)
Cr     j_l(kr) = sin(kr)/kr
Cr     Limiting cases for energy->0:
Cr     h_l -> i (2l-1)!!/(kr)^(l+1)  j_l -> (kr)^l / (2l+1)!!
Cr     Wronskian: {h,j} = -1 / i kap
Cr
Cr  *Methfessel's conventions (label as an,ak,aj) have the properties:
Cr     (1)  functions ak and aj, and an are real for e<=0
Cr     (2)  cases e .ne. 0 and e .eq. 0 same equations
Cr     Define akap = -i kappa = sqrt(-e):  akap is real and >0 for e<0.
Cr       ak_0 = exp(-akap r)/r;  ak_1 = (1 + akap*r) ak_0 / r
Cr     Relation to standard conventions:
Cr       ak_l = -i kappa^(l+1) h_l = -i (-i akap)^(l+1) h_l
Cr            = gi_l / r^(l+1)
Cr       aj_l = 1/kappa^l j_l      = 1/(-i akap)^l j_l
Cr            = fi_l r^l
Cr     Define Neumann function as:
Cr       an = ah - i kap e^l aj = ah + akap e^l aj
Cr     Solid Hankels, Bessels are defined as (CHECK)
Cr       H_L = Y_L(-grad) ak(r);   J_L = E^(-l) Y_L (-grad) aj(r)
Cr     Limiting cases for energy->0:
Cr       ak -> (2l-1)!!/r^(l+1)      aj -> r^l/(2l+1)!!
Cr
Cr  *Andersen conventions (label as K,J)
Cr       K_l = -1/(2l-1)!! (kappa avw)^(l+1) i h_l(kr)
Cr           =  1/(2l-1)!! (avw)^(l+1) ak_l
Cr           =  gi(OKA) / (r/w)^(l+1)
Cr       J_l = 1/2*(2l-1)!! (kappa avw)^(-l) j_l(kr)
Cr           = 1/2*(2l-1)!! (avw)^(-l) aj
Cr           = fi(OKA) * (r/w)^l
Cr     Define Neumann function as:
Cr       N_l = K_l - i kap e^l 2/((2l-1)!!)^2 J_l
Cr     avw is some arbitrary length scale, e.g. average WSR.
Cr     By setting loka=.true. the fi and gi are rescaled as follows:
Cr       fi -> fi(OKA) = fi * (2l-1)!!/2
Cr       gi -> gi(OKA) = gi / (2l-1)!!
Cr     Limiting cases for energy->0
Cr       H_l -> (r/w)^(-l-1)         J_l -> (r/w)^l/(2(2l+1))
Cr     NB: in exact MTO paper, Andersen makes definitions j~ and n~:
Cr       n_l~(kappa r) =  (kappa r)^l+1/(2l-1)!! n_l(kappa r)
Cr                     -> 1 - (kappa r)^2/2(2l+3) + ... as kappa r -> 0
Cr       j_l~(kappa r) =  (2l+1)!!/(kappa r)^l j_l(kappa r)
Cr                     -> 1 - (kappa r)^2/2(2l+3) + ... as kappa r -> 0
Cr
Cr
Cr  *Generation of fi and gi:  fi (and gi for y>0) are calculated for
Cr   lmax and lmax-1 by an expansion in powers of x^2=y:
Cr
Cr                     (-x^2/2)^k
Cr       fi =  Sum_k  --------------  = dum(lmx+1-l)
Cr                    k! (2l+1+2k)!!
Cr
Cr       gi = dum(lmx+2+l), y>0
Cr
Cr     For y<0:
Cr       gi(l=0) = exp(-sqrt(-y))
Cr       gi(l=1) = (1d0+sqrt(-y))*exp(-sqrt(-y))
Cr       gi(l)   = (2l+1)*gi(l-1) - y*gi(l-2)
Cr
Cr     The remaining functions are obtained by the recursion formula:
Cr
Cr       j_{l-2}(x)=(2l-1)j_{l-1}(x)/x-j{l}(x)   l=lmx-2,lmx-3,..,-lmx-1
Cr
Cr     => dum(k) = (2*lmx+5-2*k)*dum(k-1)-y*dum(k-2)   k=3,4,...,2*lmx+2
Cr
Cr     and the Neumann function are given through:
Cr
Cr           n_l(x) = j_{-l-1}*(-1)^{l+1}
Cr
Cr  *Special cases:
Cr     fi(y,l) -> 1/(2l+1)!!     for y->0
Cr     gi(y,l) -> (2l-1)!!       for y->0
Cr     gi(y,l=0) = exp(i*akap*r) for y<=0    akap=sqrt(e) Im(akap)>=0
Cr     gi(y,l=0) = cos(akap*r)   for y>0     akap=sqrt(e) Im(akap)=0
Cr   As mentioned above
Cr     fi(OKA) = fi * (2l-1)!!/2  and
Cr     gi(OKA) = gi / (2l-1)!!    so
Cr     J(E,r) = fi(OKA) * (r/w)^l     -> (r/w)^l/(2(2l+1)) for e->0
Cr     K(E,r) = gi(OKA) / (r/w)^(l+1) ->
Cb Bugs
Cb   Program never checked for lmax < -1
Cb   For x > 40 this algorithm is numerically unstable !!!
Cu Updates
Cu   23 Jul 08 bug fix: besslr doesn't make fi/gi(lmax+1) when lmax=0
Cu   19 May 04 Changed loka from logical to integer
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer loka
      integer lmin,lmax
      double precision y,fi(lmin:lmax),gi(lmin:lmax)
C Local variables:
      integer i,isn,j1,j2,k,l,lmx,lmxp1,lmxp2,nf,tlp1,ll1,ll2,nlmax
      parameter (nlmax=20)
      double precision dt,dt2,exppr,my,srmy,g1,t,tol,
     .dum(nlmax*4+2),fac2l(-nlmax:nlmax*2+3)
      logical lhank
      parameter(tol=1.d-15)
      parameter(lhank = .true.)
C Intrinsic functions:
      intrinsic dabs,dexp,dsqrt,max0

      lmx = max0(lmax,2)
      if (lmin .gt. 0) call rx(' BESSL : lmin gt 0')
      if (lmx .gt. nlmax+nlmax) then
        call rxi(' BESSL : lmax gt nlmax*2, lmax=',lmx)
      endif

C --- A table of fac2l(l)=(2l-1)!!
c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      fac2l(0) = 1d0
      do  10  l = 1, lmx+1
        fac2l(l) = fac2l(l-1) * (l+l-1)
   10 continue
      do  11  l = -1, lmin, -1
        fac2l(l) = fac2l(l+1) / (l+l+1)
   11 continue

C --- Case akap=0 ---
      if (y .eq. 0) then
        do  12  l = lmin, lmax
          fi(l) = 1/fac2l(l+1)
          gi(l) = fac2l(l)
   12   continue
        goto 100
      endif
      my = -y

C --- Get dum(1) = j_{lmx}(x)/x^{lmx} = fi(lmx)
      tlp1 = lmx+lmx+1
      dt = 1d0
      t = 1d0
      i = 0
      do  20  k = 1, 1000
        if (dabs(dt) .lt. tol) goto 21
        i = i+2
        dt2 = i+tlp1
        dt = dt*my/(i*dt2)
        t = t+dt
   20 continue
      call rx('BESSLR: series not convergent')
   21 continue
      dum(1) = t/fac2l(lmx+1)

C --- Get dum(2) = j_{lmx-1}(x)/x^{lmx-1} = fi(lmx-1) ---
      tlp1 =  tlp1-2
      dt = 1d0
      t = 1d0
      i = 0
      do  30  k = 1, 1000
        if (dabs(dt) .lt. tol) goto 31
        i = i+2
        dt2 = i+tlp1
        dt = dt*my/(i*dt2)
        t = t+dt
   30 continue
      call rx('BESSLR: series not convergent')
   31 continue
      dum(2) = t/fac2l(lmx)

C --- Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
      ll1 = lmx + lmx + 1
      ll2 = ll1 + 1
      nf = ll1
      do  40  k = 3, ll2
        nf = nf-2
        dum(k) = nf*dum(k-1) - y*dum(k-2)
   40 continue

C --- Get fi and gi from dum ---
      lmxp1 = lmx+1
      lmxp2 = lmx+2
      isn = (-1)**lmin
      do  50  k = lmin, lmax
        j1 = lmxp1-k
        j2 = lmxp2+k
        fi(k) = dum(j1)
c   ... n_l(x) = j_{-l-1}*(-1)^{l+1}
        gi(k) = dum(j2)*isn
        isn = -isn
   50 continue

C --- For E<0, use Hankel functions rather than Neumann functions ---
      if (lhank .and. y .lt. 0d0) then
        srmy = dsqrt(-y)
        gi(0) = 1d0
        g1 = 1d0+srmy
        if (lmax .ge. 1) gi(1) = g1
        if (lmax .ge. 2) then
          tlp1 = 1
          do  62  l = 2, lmax
            tlp1 = tlp1+2
            gi(l) = tlp1*gi(l-1) - y*gi(l-2)
   62     continue
        endif
        if (lmin .le. -1) then
          gi(-1) = (gi(0) - g1)/y
          tlp1 = 1
          if (lmin .le. -2) then
            do  64  l = -2, lmin,-1
              tlp1  = tlp1-2
              gi(l) = ((l+l+3)*gi(l+1) - gi(l+2))/y
   64       continue
          endif
        endif
        exppr = 1d0/dexp(srmy)
        do  66  l = lmin, lmax
          gi(l) = gi(l)*exppr
   66   continue
      endif

C --- Scaling to Andersen's 2nd generation LMTO conventions ---
  100 continue
      if (loka .eq. 1) then
        do  68  l = lmin, lmax
          fi(l) = fi(l)*fac2l(l)*0.5d0
          gi(l) = gi(l)/fac2l(l)
   68   continue
      endif

      end
      subroutine bessl2(y,lmin,lmax,fi,gi)
C- Radial part of Bessel functions, Andersen's definitions
C  See besslr for definitions.
C  For y=0, fi(l) = 0.5/(2l+1), gi(l) = 1
C     implicit none
C Passed variables:
      integer lmin,lmax
      double precision y,fi(lmin:lmax),gi(lmin:lmax)

      call besslr(y,1,lmin,lmax,fi,gi)
      end
      subroutine besslm(y,lmax,fi,gi)
C- Radial part of Bessel functions, standard definitions
C  See besslr for definitions.
C  For y=0, fi(l) = (2l-1)!!, gi(l) = 1/(2l+1)!!
C     implicit none
C Passed variables:
      integer lmax
      double precision y,fi(0:lmax),gi(0:lmax)

      call besslr(y,0,0,lmax,fi,gi)
      end
      subroutine bessl(y,lmax,fi,gi)
C- Radial part of Bessel functions, standard definitions
C  See besslr for definitions.
C  For y=0, fi(l) = (2l-1)!!, gi(l) = 1/(2l+1)!!
C     implicit none
C Passed variables:
      integer lmax
      double precision y,fi(0:lmax),gi(0:lmax)

      call besslr(y,0,0,lmax,fi,gi)
      end
C      subroutine fmain
C      implicit none
C      integer lmin,lmax,il
C      double precision fi(-3:10),gi(-3:10),y
C
C      fi = -99d0
C      gi = -99d0
C      lmin = -1
C      lmax =  2
C      y = 1.1d0
C
C      call besslr(y,0,lmin,lmax,fi(lmin),gi(lmin))
C      print *, sngl(fi(lmin:lmax+1))
C      print *, sngl(gi(lmin:lmax+1))
C
C      end

