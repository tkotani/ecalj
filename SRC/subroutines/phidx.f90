module m_phidx
  private
  public phidx
  contains
subroutine phidx(job,z,l,v,hcr,vmtz,rofi,nr,nptdif,tol, &
     e,val,slo,nn,g,gp,phi,dphi,phip,dphip,p,phia,phiap,dla,dlap)
  use m_lgunit,only:stdo
  !- Generate potential parameters for a prescribed energy or b.c.
  ! ----------------------------------------------------------------
  !i Inputs:
  !i   job:   1s digit specifies boundary conditions
  !i          0, boundary conditions specified val,slo,nn (see Remarks)
  !i          1, same as 0, but also g and e assumed to be generated
  !i          2, boundary conditions specified by energy e (see Remarks)
  !i          10s digit
  !i          1, set dphip to satisfy Wronskian condition
  !i   z:     nuclear charge
  !i   l:     l quantum number for this g
  !i   v:     spherical potential on shifted logarithmic mesh
  !i   vmtz   flat potential used to generate dla,dlap
  !i          Not used if hcr=0
  !i   hcr:   hard sphere radius.  If nonzero, dla,dlap are generated
  !i   rofi:  list of points
  !i   nr:    number of mesh points
  !i   nptdif:2 or 4 for 3- or 5- point differentiation
  !i         :You may also set nptdif=0.  Then quantities related to
  !i         :energy differences are not calculated (dlap,phip,dphip,p)
  !i   tol:   precision to which wave function is integrated
  ! o Inputs/Outputs:
  ! o         Which subset of these quantities are needed for input and
  ! o         which quantities phidx alters depends on job; see Remarks.
  ! o
  ! o  e:     On input (job=1,2), energy eigenvalue
  ! o         On output (job=0), energy eigenvalue
  ! o  val:   On input (job=0,1), val(1)=value of g(r)=r*phi(r) at rmax
  ! o         On output (job=0,2), val(1)=value of normalized g(r) at rmax
  ! o         Also on output, val(2..1+nptdif) = energy derivatives of val
  ! o  slo:   On input (job=0,1), slo(1)=radial derivative of g(r) at rmax
  ! o         On output (job=0,2), slo(1)=der. of normalized g(r) at rmax
  ! o         Also on output, slo(2..1+nptdif) = energy derivatives of slo
  ! o  nn:    On input (job=0,1), number of nodes
  ! o         On output (job=2), number of nodes
  ! o  g:     On input (job=1) Wave function times r
  ! o         (assumed normalized so that int (g*g) dr = 1)
  ! o         On output (job=0,2) normalized wave function times r
  !o  Outputs:
  !o   gp:    first nptdif energy derivatives to G
  !o   phi:   wave function at rmax, i.e. g/rmax
  !o   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
  !o   phip:  energy derivative of wave function at rmax
  !o   dphip: energy derivative of slope of wave function at rmax
  !o   p:     <gp**2> (potential parameter)
  !o   phia:  (hcr>0) value of phi at hcr boundary, i.e. g(hcr)/hcr
  !o   phiap: (hcr>0) energy derivatives of phia
  !o   dla:   (hcr>0) hcr * logarithmic derivative of phi0 at hcr boundary
  !o                  where phi0 is back-extrapolated wave function
  !o          (hcr=0) not calculated
  !o   dlap:  (hcr>0) energy derivatives of dla
  !o          (hcr=0) not calculated
  !r Remarks:
  !r   This version makes parameters related to wave function g(r)/r
  !r   defined by potential v(r).
  !r
  !r   Boundary conditions are specified in one of two ways:
  !r     job=0   val,slo,nn are specified.  On output,
  !r             val,slo are renormalized so that val=g(nr), with <gg>=1
  !r             Here energy eigenvalue e is calculated
  !r             are assumed to correspond with g.)
  !r     job=2   the energy eigenvalue e is specified.
  !r             val,slo, and nn are calculated.
  !r     job=1   Assumes that all quantities val,slo,nn,e,g have been
  !r             generated, and calculates none of them.
  !b Bugs
  !b   Not really a bug, but phidx returns redundant information in
  !b   the following variables:
  !b      phi   = val(1)/rmax
  !b      dphi  = (slo(1) - phi)/rmax
  !b      phip  = vali(1)/rmax
  !b      dphip = (sloi(1) - phip)/rmax
  !u Updates
  !u   21 Jul 04 possibly set dphip to satisfy Wronskian condition
  !u   19 Dec 01 Return val,slo,phiap,dlap for nptdif energy derivatives
  !u    7 May 01 Allow nptdif=0
  !u   29 Nov 00 printout when phidx fails to converge, softer tole
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: job,l,nr,nn,nptdif
  double precision :: z,e,vmtz,hcr,val(*),slo(*),phi,dphi,phip,dphip, &
       dla,dlap(*),p,phia,phiap(*),tol,v(nr),rofi(nr), &
       g(2*nr),gp(2*nr,4),dlax(1)
  ! Local variables
  integer :: nre,i,iprint,job0,job1
  double precision :: rmax,eb1,eb2,dele,ddde,sum,a,b,aold,dmach,tola, &
       sloi(5),vali(5),phiai(4),dlai(4),ei(4),de1,de2,del1,del2,tole
  parameter (tole=1d-10)
  ! External calls
!  external dfphi,gintsr,iprint,makdla,rseq,rsq1,rx,lgunit

  ! ... Setup: extract a and b
  job0 = mod(job,10)
  job1 = mod(job/10,10)
  !      stdo = lgunit(1)
  rmax = rofi(nr)
  dele = .002d0
  a = log(rofi(nr)/rofi(nr-1))
  tola = 8*dmach(1)
  do   i = 1, 100
     aold = a
     a = log(rofi(nr)/rofi(nr-1)/(1-exp(-(nr-1)*a))*(1-exp(-(nr-2)*a)))
     if (i > 95) write(stdo,'(i4,1p,2e15.5)') i,a,a-aold
     if (abs(a-aold) <= tola) goto 1
  enddo
  call rx('phidx failed to determine ''a'' parameter')
1 continue
  b = rmax/(dexp(a*(nr-1)) - 1d0)

  ! --- Find energy e corresponding to val,slo ---
  if (job0 == 0) then
     eb1 = -20d0
     eb2 =  20d0
     e = (eb1+eb2)/2
     !       This generates g, normalized to <g g> = 1
     call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre)
     !       Scale val, slo to correspond to normalization <g g> = 1
     val(1) = val(1)/dsqrt(sum)
     slo(1) = slo(1)/dsqrt(sum)

     ! --- Find val,slo corresponding to energy e ---
  elseif (job0 == 2) then
     !       Initial estimate for val,slo
     call rsq1(0,e,l,z,v,nr,g,val,slo,nn,a,b,rofi,nr)
     !       Adjust slope iteratively until ei(slope)-e is within tole
     ei(1) = e
     eb1 = e-.1d0
     eb2 = e+.1d0
     do  i = 1, 5+5
        call rseq(eb1,eb2,ei,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre)
        if (iprint() > 0 .AND. i > 5) then
           !            call awrit7(' PHIDX  Z=%d  l=%i  nod=%i  bc=%;4g %;4g'//
           !     .      '  e(bc)=%;4g  e(bc)-e=%;4g',' ',80,stdo,z,l,nn,val(1),
           !     .      slo(1),ei(1),ei(1)-e)
           write(stdo,"(' PHIDX  Z=',f9.4,'  l nod=',i0,x,i0, &
                ' val slo=',2f9.4,' e(bc)=',f9.4,' e(bc)-e=',f9.4)") z,l,nn,val(1), &
                slo(1),ei(1),ei(1)-e
        endif

        if (abs(ei(1)-e) < tole) goto 2
        slo(1) = slo(1) + (ei(1)-e) * val(1) / g(nr)**2
     enddo
     call rx('phidx failed to converge')
2    continue
     val(1) = val(1)/dsqrt(sum)
     slo(1) = slo(1)/dsqrt(sum)
  elseif (job0 /= 1) then
     call rx('phidx: bad job')
  endif
  if (hcr /= 0) call makdla(e-vmtz,l,hcr,slo(1),val(1),rmax,phia,dla)

  if (nptdif /= 0) then
     ddde = -rmax/g(nr)**2
     ei(1) = 1
     ei(2) = -1
     ei(3) = 1.5d0
     ei(4) = -1.5d0
     eb1 = e-.1d0
     eb2 = e+.1d0

     ! cccccccccccccccccccccccccccccccccc
     ! choice 2 takao. For deep semicore, to avoid error in rseq (Mg dimer in 10\A cubic cell).
     ei=ei/2d0 !In future, we have to use better solution.
     !      you may need to use ei/3.0 or better algolism
     !  I had a problem that it results in warning rseq,
     !  because two exponential solution makes huge changes due to slight difference of energy
     !  and node number can not be the same.
     !    ATOM= Mg Z= 12 R= 3.000
     !     RSMH=   1.500 1.500 1.500 1.500 EH=  -1.0 -1.0 -1.0 -1.0
     !     RSMH2=  1.500 1.500 1.500 1.500 EH2= -2.0 -2.0 -2.0 -2.0
     !     PZ=0,12.9 P=0,3.3     KMXA={kmxa}  LMX=3 LMXA=4
     ! cccccccccccccccccccccccccccccccccc

     do  10  i = 1, nptdif
        sloi(i) = slo(1) + dele*ei(i)*ddde*val(1)/rmax
        ! cccccccccccccccc
        !          print *,' nptdif sloi val=',i,sloi(i),val(1)
        ! cccccccccccccccc
        ei(i) = e + dele*ei(i)
        call rseq(eb1,eb2,ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i), &
             sum,a,b,rofi,nr,nre)
        vali(i) = val(1)/dsqrt(sum)
        sloi(i) = sloi(i)/dsqrt(sum)
        if (hcr /= 0) call makdla(ei(i)-vmtz,l,hcr,sloi(i),vali(i), &
             rmax,phiai(i),dlai(i))

10   enddo

     de1  = (ei(1) - ei(2))/2
     del1 = (ei(1) + ei(2))/2 - e
     de2  = (ei(3) - ei(4))/2
     del2 = (ei(3) + ei(4))/2 - e
     !     Energy derivatives of value and slope
     call dfphi(de1,del1,de2,del2,1,val,vali,nptdif.eq.4)
     call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif.eq.4)
     !     Energy derivatives of dla
     if (hcr /= 0) then
        call dfphi(de1,del1,de2,del2,1,dlax,dlai,nptdif.eq.4)
        call dfphi(de1,del1,de2,del2,1,dlax,phiai,nptdif.eq.4)
        do  12  i = 1, nptdif
           dlap(i)  = dlai(i)
           phiap(i) = phiai(i)
12      enddo
     endif
     !     Energy derivatives of g
     call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif.eq.4)
     !     p = integral <gp gp>
     call gintsr(gp,gp,a,b,nr,z,e,l,v,rofi,p)
  endif

  !     phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
  phi = val(1)/rmax
  dphi = (slo(1) - phi)/rmax
  if (nptdif /= 0) then
     phip = vali(1)/rmax
     dphip = (sloi(1) - phip)/rmax
  endif
  !     Copy energy derivatives sloi to slo(2..)
  if (nptdif /= 0) then
     call dcopy(nptdif,sloi,1,slo(2),1)
     call dcopy(nptdif,vali,1,val(2),1)
  endif
  !     Set dphip from Wronskian condition:
  !     phi*dphip - dphi*phip = -1/rmax**2 =>
  !     dphip = (dphi*phip - 1/rmax**2)/phi
  if (nptdif /= 0 .AND. job1 /= 0) then
     !        print *, dphip,(dphi*phip - 1/rmax**2)/phi,
     !     . dphip - (dphi*phip - 1/rmax**2)/phi
     dphip = (dphi*phip - 1/rmax**2)/phi
  endif
  if (iprint() >= 111) write(stdo,749) phi,dphi,phip,dphip
749 format(' PHIDOT:  phi,phip,phip,dphip=',4f12.6)
end subroutine phidx


subroutine makdla(kbak,l,hcr,slo,val,rmax,phia,dla)
  !- Gets logarithmic derivative of wave function at hard core radius
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   kbak  :kinetic energy for back extrapolation
  !i   l     :angular momentum
  !i   hcr   :hard-core screening radius a, in atomic units
  !i   slo   :slope of g at rmax, where g = rmax * radial w.f.
  !i   val   :value of g at rmax, where g = rmax * radial w.f.
  !i   rmax  :radius at which val,slo are evaluated, in atomic units
  !o Outputs:
  !o   phia  :value of radial w.f. at a
  !o         :If a=rmax, phia = val/rmax
  !o   dla   :logarithmic derivative at a, i.e. a/phi*(dphi/dr)_a
  !o         :If a=rmax, dla = log deriv of phi = rmax*slo/val - 1
  !r Remarks:
  !r This was adapted from the Stuttgart third-generation LMTO package.
  ! ----------------------------------------------------------------------
  implicit none
  integer :: l,lmax,lx
  double precision :: kbak,hcr,slo,val,rmax,phia,dla
  integer :: nlmax,lmin
  parameter(nlmax=20)
  double precision :: er2,fi(0:nlmax+1),gi(0:nlmax+1),wn,wj,dlr,dj,dn, &
       sigma,phi,rdphia,rdfi,rdgi, fac2l(0:l+1)
  lmin=0
  lmax=l+1
  fac2l(0) = 1d0
  do  lx = 1, lmax
     fac2l(lx) = fac2l(lx-1) * (lx+lx-1)
  enddo
  er2 = kbak*rmax*rmax
  call bessl(er2,l+1,fi(0),gi(0)) !bessl2->bessl
  fi(lmin:lmax) = fi(lmin:lmax)*fac2l(lmin:lmax)*0.5d0 !Andersen factor
  gi(lmin:lmax) = gi(lmin:lmax)/fac2l(lmin:lmax)       !Andersen factor
  !     phi,dlr are value and logarithmic derivative r/phi dphi/dr at rmax
  !     free dlr=dj
  phi =  val/rmax
  dlr =  rmax*slo/val - 1
  !     dj,dn are the logarithmic derivatives of Bessel and Hankels
  dj  = l-fi(l+1)/fi(l)/(l+l+1)*er2
  dn  = l-gi(l+1)/gi(l)*(l+l+1)
  !     wj,wn are the amounts of Bessel and Hankel making up phi0
  wj  = (dlr-dj)*fi(l)*phi
  wn  = (dlr-dn)*gi(l)*phi
  sigma = hcr/rmax
  er2 = kbak*hcr**2
  call bessl(er2,l+1,fi(0),gi(0)) !call bessl2(er2,0,l+1,fi(0),gi(0))
  fi(lmin:lmax) = fi(lmin:lmax)*fac2l(lmin:lmax)*0.5d0 !Andersen factor
  gi(lmin:lmax) = gi(lmin:lmax)/fac2l(lmin:lmax)       !Andersen factor
  rdgi = l*gi(l) - gi(l+1)*(l+l+1)
  rdfi = l*fi(l) - fi(l+1)/(l+l+1)*er2
  phia   = 2d0*(wn*fi(l)*sigma**l - wj*gi(l)*sigma**(-l-1))
  rdphia = 2d0*(wn*rdfi*sigma**l  - wj*rdgi*sigma**(-l-1))
  dla    = rdphia/phia
end subroutine makdla

subroutine gintsr(g1,g2,a,b,nr,z,e,l,v,rofi,sum)
  !- Integrate inner product of two wave equations
  ! ----------------------------------------------------------------
  !i   g1,g2 :First and second radial wave functions
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   nr    :number of radial mesh points
  !i   z     :nuclear charge
  !i   e     :energy
  !i   l     :l quantum number of g1,g2
  !i   v     :spherical potential
  !i   rofi  :radial mesh points
  !o Outputs:
  !o   sum   :inner product
  !r Remarks:
  !r   Uses Simpson's rule
  ! ----------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nr,l
  double precision :: a,b,z,e,sum,g1(nr,2),g2(nr,2),v(nr),rofi(nr)
  ! ... Local parameters
  integer :: i,ir
  double precision :: fllp1,c,r
  double precision :: tmc,fi
  !     Speed of light, or infinity in nonrelativistic case
  common /cc/ c
  tmc(i,r) = c - (v(i) - 2d0*z/r - e)/c
  fi(i,r) = (r+b)*(g1(i,1)*g2(i,1)*(1 + fllp1/(tmc(i,r)*r)**2) &
       + g1(i,2)*g2(i,2))
  fllp1 = l*(l+1)
  sum = 0d0
  do  10  ir = 2, nr-1, 2
     sum = sum + fi(ir,rofi(ir))
10 enddo
  sum = 2*sum
  do  11  ir = 3, nr-2, 2
     sum = sum + fi(ir,rofi(ir))
11 enddo
  sum = (2*sum + fi(nr,rofi(nr)))*a/3
end subroutine gintsr
subroutine dfphi(de1,del1,de2,del2,nr,g,gp,fivep)
  !- Numerically differentiates g using five point finite difference
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   de1   :see Remarks
  !i   del1  :see Remarks
  !i   de2   :see Remarks
  !i   del2  :see Remarks
  !i   nr    :number of mesh points
  !i   g     :normalized wave function times r
  !i   fivep :if true, five-point formula; else three-point formula
  ! o Inputs/Outputs
  ! o  gp    :On input,
  ! o        :contains g at de1+del1, -de1+del1, de2+del2, -de2+del2
  ! o        :On output, contains energy derivatives of g:
  ! o        :two derivatives or four derivatives, depending on fivep
  !r Remarks
  !r   A simple five-point estimate for the numerical differentiation for
  !r   the first four derivatives of phi is based on integration of phi
  !r   at four energies de1+del1, -de1+del1, de2+del2, -de2+del2.  The
  !r   offset del1 and del2 are present because (for consistency's sake)
  !r   RSEQ is called which takes evenly spaced increments in slopes
  !r   about the average.  But the deviations del1 and del2 are only
  !r   nonzero to second order in the increment de, and this can be
  !r   exploited to obtain accurate five point estimates without having
  !r   to solve for the four simultaneous equations.
  !r
  !r   A three point for energy differences about Enu of 0, ep = de+del,
  !r   em = -de+del gives
  !r
  !r   (1)  gp =  2 del/(de^2-del^2) g(0) +
  !r              (de-del)/(de+del)/(2 de) g(ep) -
  !r              (de+del)/(de-del)/(2 de) g(em) +
  !r
  !r              1/6 gppp (de^2-del^2) + 1/12 gpppp del (de^2-del^2) + ...
  !r
  !r
  !r   (2) gpp = -2/(de^2-del^2) g(0) +
  !r              1/(de+del)/de g(ep) +
  !r              1/(de-del)/de g(em) -
  !r
  !r              2/3 gppp del + 1/12 gpppp del (de^2 + 3 del^2) + ...
  !r
  !r
  !r   The gppp term in (1) can be knocked out by taking two three point
  !r   formulas in linear combination (de1^2-del1^2) and (de2^2-del2^2)
  !r   leaving only fourth and higher order terms.  Because del is of
  !r   order de^2, the fourth order term is of the same order as the
  !r   fifth order term and there is no advantage in eliminating it.
  !r   Also from the difference between this more accurate (5-point)
  !r   estimate for gp, an estimate for gppp can be made as
  !r
  !r   (3) gppp = -6 ( gp(five point) - gp(three point) ) /(de1^2-del1^2)
  !r
  !r             + 1/2 gpppp del1
  !r
  !r   which is again accurate to order de^2.  Once gppp is known to this
  !r   order the term proportional to gppp in three point estimate for
  !r   gpp can be subtracted out directly and the gpppp term can be
  !r   eliminated by taking three point formulas (with the gppp term
  !r   removed) in linear combinations (de1^2+3*del1^2) and
  !r   (de2^2+3*del2^2), and finally the fourth order derivative can be
  !r   estimated from the difference in the five-point estimate for gpp
  !r   and the three point estimate.
  !  ----------------------------------------------------------------
  !     implicit none
  !     Passed parameters
  integer :: nr
  logical :: fivep
  double precision :: de1,del1,de2,del2,g(nr),gp(nr,4)
  !     Local parameters
  integer :: i
  double precision :: gp5p,gpp5p,gppp,gpp3p,gpp32
  double precision :: xx1,xx2,xx3,xx4,gp3p,w01,w11,w21,w02,w12,w22, &
       w01d,w11d,w21d,w02d,w12d,w22d,wp1d,wp2d,gpppp

  ! --- Constants common to 3-point and 5-point formulas ---
  w01 = 2*del1/(de1**2-del1**2)
  w11 = (de1-del1)/(de1+del1)/(2*de1)
  w21 = (de1+del1)/(de1-del1)/(2*de1)
  w01d = -2/(de1**2-del1**2)
  w11d = 1/(de1+del1)/de1
  w21d = 1/(de1-del1)/de1

  if ( .NOT. fivep) goto 20
  if (dabs(del1/de1) > .1 .OR.  dabs(del1/de1) > .1) then
     !       if (iprint() .ge. 20) print *, 'dfphi:  large del; use 3 point'
     goto 20
  endif

  ! --- Extra constants for 5-point formula ---
  xx1 = de1**2 - del1**2
  xx2 = de2**2 - del2**2
  xx3 = de1**2 + 3*del1**2
  xx4 = de2**2 + 3*del2**2
  w02 = 2*del2/(de2**2-del2**2)
  w12 = (de2-del2)/(de2+del2)/(2*de2)
  w22 = (de2+del2)/(de2-del2)/(2*de2)
  wp1d = 2d0*del1/3
  w02d = -2/(de2**2-del2**2)
  w12d = 1/(de2+del2)/de2
  w22d = 1/(de2-del2)/de2
  wp2d = 2d0*del2/3

  do  10  i = 1, nr

     ! Three point formula for gp; store in temporary gp3p
     gp3p = w01*g(i) + w11*gp(i,1) - w21*gp(i,2)

     ! Five point estimate for gp
     gp5p = (xx2*gp3p - xx1*(w02*g(i) + w12*gp(i,3) - w22*gp(i,4))) &
          /(xx2-xx1)

     ! Difference between five point and three point gives estimate for gppp
     gppp = -6/xx1*(gp5p - gp3p)

     ! Three point estimates for gpp with correction for gppp
     gpp3p = w01d * g(i) + w11d * gp(i,1) + w21d*gp(i,2) &
          -wp1d * gppp
     gpp32 = w02d * g(i) + w12d * gp(i,3) + w22d*gp(i,4) &
          -wp2d * gppp

     ! Five point estimate for gpp with correction for gppp
     gpp5p = (gpp3p*xx4 - gpp32*xx3) / (xx4 - xx3)

     ! Difference between five point and three point gives est for gpppp
     gpppp = -12/xx3*(gpp5p - gpp3p)

     gp(i,1) = gp5p
     gp(i,2) = gpp5p
     gp(i,3) = gppp
     gp(i,4) = gpppp
10 enddo
  return

  ! Three point formulae:  only gp, gpp calculated
20 continue
  do  30  i = 1, nr
     gp3p    = w01*g(i)  + w11*gp(i,1)  - w21*gp(i,2)
     gp(i,2) = w01d*g(i) + w11d*gp(i,1) + w21d*gp(i,2)
     gp(i,1) = gp3p
30 enddo

end subroutine dfphi

end module m_phidx
