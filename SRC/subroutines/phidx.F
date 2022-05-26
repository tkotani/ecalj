      subroutine phidx(job,z,l,v,hcr,vmtz,rofi,nr,nptdif,tol,
     .e,val,slo,nn,g,gp,phi,dphi,phip,dphip,p,phia,phiap,dla,dlap)
      use m_lgunit,only:stdo
C- Generate potential parameters for a prescribed energy or b.c.
C ----------------------------------------------------------------
Ci Inputs:
Ci   job:   1s digit specifies boundary conditions
Ci          0, boundary conditions specified val,slo,nn (see Remarks)
Ci          1, same as 0, but also g and e assumed to be generated
Ci          2, boundary conditions specified by energy e (see Remarks)
Ci          10s digit
Ci          1, set dphip to satisfy Wronskian condition
Ci   z:     nuclear charge
Ci   l:     l quantum number for this g
Ci   v:     spherical potential on shifted logarithmic mesh
Ci   vmtz   flat potential used to generate dla,dlap
Ci          Not used if hcr=0
Ci   hcr:   hard sphere radius.  If nonzero, dla,dlap are generated
Ci   rofi:  list of points
Ci   nr:    number of mesh points
Ci   nptdif:2 or 4 for 3- or 5- point differentiation
Ci         :You may also set nptdif=0.  Then quantities related to
Ci         :energy differences are not calculated (dlap,phip,dphip,p)
Ci   tol:   precision to which wave function is integrated
Cio Inputs/Outputs:
Cio         Which subset of these quantities are needed for input and
Cio         which quantities phidx alters depends on job; see Remarks.
Cio
Cio  e:     On input (job=1,2), energy eigenvalue
Cio         On output (job=0), energy eigenvalue
Cio  val:   On input (job=0,1), val(1)=value of g(r)=r*phi(r) at rmax
Cio         On output (job=0,2), val(1)=value of normalized g(r) at rmax
Cio         Also on output, val(2..1+nptdif) = energy derivatives of val
Cio  slo:   On input (job=0,1), slo(1)=radial derivative of g(r) at rmax
Cio         On output (job=0,2), slo(1)=der. of normalized g(r) at rmax
Cio         Also on output, slo(2..1+nptdif) = energy derivatives of slo
Cio  nn:    On input (job=0,1), number of nodes
Cio         On output (job=2), number of nodes
Cio  g:     On input (job=1) Wave function times r
Cio         (assumed normalized so that int (g*g) dr = 1)
Cio         On output (job=0,2) normalized wave function times r
Co  Outputs:
Co   gp:    first nptdif energy derivatives to G
Co   phi:   wave function at rmax, i.e. g/rmax
Co   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
Co   phip:  energy derivative of wave function at rmax
Co   dphip: energy derivative of slope of wave function at rmax
Co   p:     <gp**2> (potential parameter)
Co   phia:  (hcr>0) value of phi at hcr boundary, i.e. g(hcr)/hcr
Co   phiap: (hcr>0) energy derivatives of phia
Co   dla:   (hcr>0) hcr * logarithmic derivative of phi0 at hcr boundary
Co                  where phi0 is back-extrapolated wave function
Co          (hcr=0) not calculated
Co   dlap:  (hcr>0) energy derivatives of dla
Co          (hcr=0) not calculated
Cr Remarks:
Cr   This version makes parameters related to wave function g(r)/r
Cr   defined by potential v(r).
Cr
Cr   Boundary conditions are specified in one of two ways:
Cr     job=0   val,slo,nn are specified.  On output,
Cr             val,slo are renormalized so that val=g(nr), with <gg>=1
Cr             Here energy eigenvalue e is calculated
Cr             are assumed to correspond with g.)
Cr     job=2   the energy eigenvalue e is specified.
Cr             val,slo, and nn are calculated.
Cr     job=1   Assumes that all quantities val,slo,nn,e,g have been
Cr             generated, and calculates none of them.
Cb Bugs
Cb   Not really a bug, but phidx returns redundant information in
Cb   the following variables:
Cb      phi   = val(1)/rmax
Cb      dphi  = (slo(1) - phi)/rmax
Cb      phip  = vali(1)/rmax
Cb      dphip = (sloi(1) - phip)/rmax
Cu Updates
Cu   21 Jul 04 possibly set dphip to satisfy Wronskian condition
Cu   19 Dec 01 Return val,slo,phiap,dlap for nptdif energy derivatives
Cu    7 May 01 Allow nptdif=0
Cu   29 Nov 00 printout when phidx fails to converge, softer tole
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer job,l,nr,nn,nptdif
      double precision z,e,vmtz,hcr,val(*),slo(*),phi,dphi,phip,dphip,
     .dla,dlap(*),p,phia,phiap(*),tol,v(nr),rofi(nr),
     .g(2*nr),gp(2*nr,4)
C Local variables
      integer nre,i,iprint,job0,job1
      double precision rmax,eb1,eb2,dele,ddde,sum,a,b,aold,dmach,tola,
     .sloi(5),vali(5),phiai(4),dlai(4),ei(4),de1,de2,del1,del2,tole
      parameter (tole=1d-10)
C External calls
      external dfphi,gintsr,iprint,makdla,rseq,rsq1,rx,lgunit

C ... Setup: extract a and b
      job0 = mod(job,10)
      job1 = mod(job/10,10)
c      stdo = lgunit(1)
      rmax = rofi(nr)
      dele = .002d0
      a = log(rofi(nr)/rofi(nr-1))
      tola = 8*dmach(1)
      do   i = 1, 100
        aold = a
        a = log(rofi(nr)/rofi(nr-1)/(1-exp(-(nr-1)*a))*(1-exp(-(nr-2)*a)))
        if (i .gt. 95) write(stdo,'(i4,1p,2e15.5)') i,a,a-aold
        if (abs(a-aold) .le. tola) goto 1
      enddo
      call rx('phidx failed to determine ''a'' parameter')
    1 continue
      b = rmax/(dexp(a*(nr-1)) - 1d0)

C --- Find energy e corresponding to val,slo ---
      if (job0 .eq. 0) then
        eb1 = -20d0
        eb2 =  20d0
        e = (eb1+eb2)/2
C       This generates g, normalized to <g g> = 1
        call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre)
C       Scale val, slo to correspond to normalization <g g> = 1
        val(1) = val(1)/dsqrt(sum)
        slo(1) = slo(1)/dsqrt(sum)

C --- Find val,slo corresponding to energy e ---
      elseif (job0 .eq. 2) then
C       Initial estimate for val,slo
        call rsq1(0,e,l,z,v,nr,g,val,slo,nn,a,b,rofi,nr)
C       Adjust slope iteratively until ei(slope)-e is within tole
        ei(1) = e
        eb1 = e-.1d0
        eb2 = e+.1d0
        do  i = 1, 5+5
          call rseq(eb1,eb2,ei,tol,z,l,nn,val,slo,v,g,sum,a,b,rofi,nr,nre)
          if (iprint() .gt. 0 .and. i .gt. 5) then
c            call awrit7(' PHIDX  Z=%d  l=%i  nod=%i  bc=%;4g %;4g'//
c     .      '  e(bc)=%;4g  e(bc)-e=%;4g',' ',80,stdo,z,l,nn,val(1),
c     .      slo(1),ei(1),ei(1)-e)
            write(stdo,"(' PHIDX  Z=',f9.4,'  l nod=',i0,x,i0,
     &      ' val slo=',2f9.4,' e(bc)=',f9.4,' e(bc)-e=',f9.4)") z,l,nn,val(1),
     .      slo(1),ei(1),ei(1)-e
          endif

          if (abs(ei(1)-e) .lt. tole) goto 2
          slo(1) = slo(1) + (ei(1)-e) * val(1) / g(nr)**2
        enddo
        call rx('phidx failed to converge')
    2   continue
        val(1) = val(1)/dsqrt(sum)
        slo(1) = slo(1)/dsqrt(sum)
      elseif (job0 .ne. 1) then
        call rx('phidx: bad job')
      endif
      if (hcr .ne. 0) call makdla(e-vmtz,l,hcr,slo,val,rmax,phia,dla)

      if (nptdif .ne. 0) then
        ddde = -rmax/g(nr)**2
        ei(1) = 1
        ei(2) = -1
        ei(3) = 1.5d0
        ei(4) = -1.5d0
        eb1 = e-.1d0
        eb2 = e+.1d0

cccccccccccccccccccccccccccccccccccc
c choice 2 takao. For deep semicore, to avoid error in rseq (Mg dimer in 10\A cubic cell).
        ei=ei/2d0 !In future, we have to use better solution.
c      you may need to use ei/3.0 or better algolism
c  I had a problem that it results in warning rseq,
c  because two exponential solution makes huge changes due to slight difference of energy
c  and node number can not be the same.
c    ATOM= Mg Z= 12 R= 3.000
c     RSMH=   1.500 1.500 1.500 1.500 EH=  -1.0 -1.0 -1.0 -1.0
c     RSMH2=  1.500 1.500 1.500 1.500 EH2= -2.0 -2.0 -2.0 -2.0
c     PZ=0,12.9 P=0,3.3     KMXA={kmxa}  LMX=3 LMXA=4
cccccccccccccccccccccccccccccccccccc

        do  10  i = 1, nptdif
          sloi(i) = slo(1) + dele*ei(i)*ddde*val(1)/rmax
cccccccccccccccccc
c          print *,' nptdif sloi val=',i,sloi(i),val(1)
cccccccccccccccccc
          ei(i) = e + dele*ei(i)
          call rseq(eb1,eb2,ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i),
     .    sum,a,b,rofi,nr,nre)
          vali(i) = val(1)/dsqrt(sum)
          sloi(i) = sloi(i)/dsqrt(sum)
          if (hcr .ne. 0) call makdla(ei(i)-vmtz,l,hcr,sloi(i),vali(i),
     .    rmax,phiai(i),dlai(i))

   10   continue

        de1  = (ei(1) - ei(2))/2
        del1 = (ei(1) + ei(2))/2 - e
        de2  = (ei(3) - ei(4))/2
        del2 = (ei(3) + ei(4))/2 - e
C     Energy derivatives of value and slope
        call dfphi(de1,del1,de2,del2,1,val,vali,nptdif.eq.4)
        call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif.eq.4)
C     Energy derivatives of dla
        if (hcr .ne. 0) then
          call dfphi(de1,del1,de2,del2,1,dla,dlai,nptdif.eq.4)
          call dfphi(de1,del1,de2,del2,1,phia,phiai,nptdif.eq.4)
          do  12  i = 1, nptdif
            dlap(i)  = dlai(i)
            phiap(i) = phiai(i)
   12     continue
        endif
C     Energy derivatives of g
        call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif.eq.4)
C     p = integral <gp gp>
        call gintsr(gp,gp,a,b,nr,z,e,l,v,rofi,p)
      endif

C     phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
      phi = val(1)/rmax
      dphi = (slo(1) - phi)/rmax
      if (nptdif .ne. 0) then
        phip = vali(1)/rmax
        dphip = (sloi(1) - phip)/rmax
      endif

C     Copy energy derivatives sloi to slo(2..)
      if (nptdif .ne. 0) then
        call dcopy(nptdif,sloi,1,slo(2),1)
        call dcopy(nptdif,vali,1,val(2),1)
      endif

C     Set dphip from Wronskian condition:
C     phi*dphip - dphi*phip = -1/rmax**2 =>
C     dphip = (dphi*phip - 1/rmax**2)/phi
      if (nptdif .ne. 0 .and. job1 .ne. 0) then
C        print *, dphip,(dphi*phip - 1/rmax**2)/phi,
C     . dphip - (dphi*phip - 1/rmax**2)/phi
        dphip = (dphi*phip - 1/rmax**2)/phi
      endif

      if (iprint() .ge. 111) write(stdo,749) phi,dphi,phip,dphip
  749 format(' PHIDOT:  phi,phip,phip,dphip=',4f12.6)

      end

