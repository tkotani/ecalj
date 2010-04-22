      subroutine hansmz(r,e,rsml,xi,fi,lmax)
C- Smoothed hankel (e < 0) or neumann functions (e > 0) for l=0...lmax  
C ----------------------------------------------------------------------
Ci Inputs
Ci   r    : point at which the function is evaluated 
Ci   e    : energy   
Ci   rsml : vector of l-dependent smoothing radii of Hankels
Ci        : EITHER must be specified for 0..lmax
Ci        : OR     rsmh(0) = const < 0. Implies rsmh(l) = -const for all l
Ci   lmax : highest l for which to evaluate xi (must be positive or 0) 
Co Outputs:
Co   xi(0:lmax): radial part of smoothed hankel or neumann function 
Co             : divided by r**l
Co   fi(0:lmax): radial part of unsmoothed bessel function * k^{2l+1)/r**l
Co             : returned if e > 0, otherwise not referenced (see Remarks)
Cr Remarks
Cr   1) xi is evaluated by power series expansion for small r (r < rc1) and 
Cr   by upward recursion for large r (rc1 < r < rc2). Power series is 
Cr   terminated  at whatever order is needed to bring the convergence to 
Cr   a relative precision of 'tol'. If r > rc2, xi is replaced with the
Cr   unsmoothed function.  
Cr   2) For e > 0, hs_l = ns_l + i*e^{l+1/2}*j_l
Cu Updates
Cu   04 Jul 08 (S. Lozovoi) Bug fix
Cu   11 May 07 (S. Lozovoi) First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmax
      double precision r,e,xi(0:lmax),fi(0:lmax),rsml(0:lmax)
C ... Local parameters
      integer n0
      parameter (n0=10)
      integer il,ill,n,nmax
      double precision a0(0:n0),rsmx(0:n0),wxi(0:n0),wfi(0:n0)
      double precision rsm,rsm0,a,akrs
      double precision tol,srpi  
      double precision rm2,a2,ta,ta2,ta2l,add,akap,al,cc,dudr,ema2r2,
     .  rhs,rc1,rc2,sum,r2n,u,w,uminus,uplus,radd,gl,rfac,aexp,y0                                                            
      double precision derfc
      double complex zerfc,zexp,vm,vp
      logical lpos,lbes
      parameter (nmax=1000, tol=1d-20, srpi=1.77245385090551602729817d0)

      if (lmax .gt. n0) call rx('hansmz: lmax gt n0')
      if (lmax .lt. 0) call rx('hansmz: bad lmax')

      lbes = .false.
      lpos = (e .gt. 0d0)
      akap  = dsqrt(dabs(e))

c ---- return zero if exponential part very small -------
c      if (akap*r.gt.80d0) then
c        do  27  l = 0, ill
c  27    xi(l)=0d0
c        return
c        endif

C ... Handle negative smoothing radii
      if (rsml(0) .lt. 0d0) then
        call dvset(rsmx,1,lmax+1,-rsml(0))
      else
        call dcopy(lmax+1,rsml,1,rsmx,1)
      endif

C Start big loop over smoothing radii
      rsm0 = -1d0
      do  ill = lmax, 0, -1
        rsm = rsmx(ill)
        if (dabs(rsm-rsm0) .gt. tol) then 
          rsm0 = rsm
C ... For r > rc2 approximate smooth Hankels with normal ones
          akrs = akap*rsm*0.5d0
          rc2 = (akrs + dsqrt(akrs**2 - dlog(tol)))*rsm
C ... This rc1 generates a relative precision of ~10^-15 for r~rc1
C     and machine precision for r>>rc1 or r<<rc1.
C     For l>6 and r close to rc1, the precision degrades somewhat.
          rc1 = (1.4d0+dble(ill)/20)*rsm

C ... Make smooth hankels unless negligible smoothing or r > rc2 
          if (rsm .ge. 1d-12 .and. r .le. rc2) then  

C --- Setup ---
            a = 1.d0/rsm
            ta = a+a
            a2 = a*a
            ta2 = ta*a
            y0 = dexp(e/(ta*ta))/srpi
            cc = 4d0*a2*a*y0

C --- Power series for small r ---
            if (r .lt. rc1) then 

C ... 0 order coefficients 
              if (lpos) then
                al = dimag(zerfc(dcmplx(0d0,akrs)))/akap
              else
                al = derfc(akrs)/akap
              endif
              rhs = cc/ta2
              do  il = 0, ill
                al = (e*al + rhs)/(2*il+1)
                a0(il) = al  
                rhs = rhs*ta2
              enddo
C ... recursion in n 
              ta2l = 1d0
              do  il = 0, ill
                rhs = cc*ta2l
                sum = a0(il)
                add = sum
                r2n = 1d0
                do  n = 1, nmax
                  add = -(e*add + rhs)/(2*n*(2*n+(il+il+1)))
                  r2n = r2n*(r*r)
                  radd = add*r2n
                  sum = sum+radd
                  if (dabs(radd) .lt. tol) goto 22
                  rhs = -rhs*a2/n
                enddo   
                call rx(' hansmz: power series did not converge')
   22           continue
                wxi(il) = sum
                ta2l = ta2l*ta2
              enddo    

C --- rc1 < r < rc2 : smoothed hankels by upward recursion ---
            else
              rm2 = 1d0/(r*r)
              ema2r2 = dexp(-a2*r*r)
C ... make wxi(0) and wxi(1) explicitly
              if (lpos) then
                zexp = cdexp(dcmplx(0d0,akap*r))
                vm = (1d0 - zerfc(dcmplx( r*a,akrs)))*zexp
                vp = (1d0 - zerfc(dcmplx(-r*a,akrs)))/zexp
                wxi(0) = dreal(vm - vp) * (0.5d0/r)
                if (ill .ge. 1) then 
                  w = -dimag(vm + vp)*0.5d0
                  dudr = akap*w + ta*y0*ema2r2
                  wxi(1) = (wxi(0) - dudr)*rm2
                endif 
              else
                aexp = dexp(akap*r)
                uminus = derfc(akrs - r*a)/aexp 
                uplus = derfc(akrs + r*a)*aexp
                u = .5d0 * (uminus - uplus)/r
                wxi(0) = u
                if (ill .ge. 1) then 
                  w = -.5d0*(uminus + uplus)
                  dudr = akap*w + ta*y0*ema2r2 
                  wxi(1) = (u - dudr)*rm2
                endif 
              endif
C ... higher l by recursion 
              if (ill .ge. 2) then
                gl = cc * ema2r2
                do  il = 2, ill
                  wxi(il) = ((2*il-1)*wxi(il-1) - e*wxi(il-2) - gl)*rm2 
                  gl = ta2*gl
                enddo  
              endif 
            endif
C --- r > rc2: make unsmoothed hankels  ---
          else
            rm2 = 1d0/(r*r)
            if (.not. lbes) then
              lbes = .true.
              call besslr(e*r*r,0,0,lmax,wfi,fi)
              rfac = r
              do  il = 0, ill
                rfac = rfac*rm2 
                fi(il) = rfac*fi(il)
              enddo 
            endif
            call dcopy(ill+1,fi,1,wxi,1)
          endif 
        endif
        xi(ill) = wxi(ill)
 
C --- End big loop over smoothing radii
      enddo

C --- if e > 0, make also unsmoothed bessel functions  ---
      if (lpos) then
c ... make j unless already calculated (fi is a dummy array here)
        if (.not. lbes) call besslr(e*r*r,0,0,lmax,wfi,fi)
C ... scale j by k^(2l+1)
        rfac = akap
        do  il = 0, lmax
          fi(il) = wfi(il)*rfac
          rfac = rfac*e 
        enddo
      endif

      end

