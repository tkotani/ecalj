      subroutine strxq(mode,e,q,p,nlma,nlmh,ndim,alat,vol,awald,nkd,nkq,
     .  dlv,qlv,cg,indxcg,jcg,s,sd)
C- One-center expansion coefficents to j of Bloch summed h (strux)
C ----------------------------------------------------------------
Cr  Onsite contribution is not contained in the bloch sum in the case of p=0.
cr  See job=10 for hsmq.
Ci Inputs:
Ci   mode  :1's digit (not implemented)
Ci         :1: calculate s only
Ci         :2: calculate sd only
Ci         :any other number: calculate both s and sdot
Ci   e     :energy of Hankel function.  e must be <=0
Ci   q     :Bloch wave number
Ci   p     :position of Hankel function center;
Ci         :structure constants are for expansion about the origin
Ci   nlma  :Generate coefficients S_R'L',RL for L' < nlma
Ci   nlmh  :Generate coefficients S_R'L',RL for L  < nlmh
Ci   ndim  :leading dimension of s,sdot
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   vol   :cell volume
Ci   awald :Ewald smoothing parameter
Ci   nkq   :number of direct-space lattice vectors
Ci   nkq   :number of reciprocal-space lattice vectors
Ci   dlv   :direct-space lattice vectors, units of alat
Ci   qlv   :reciprocal lattice vectors, units of 2pi/alat
Ci   cg    :Clebsch Gordon coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L quantum number for the C.G. coefficients (scg.f)
Co Outputs
Co   s      :structure constant matrix S_R'L',RL
Co   sd     :Energy derivative of s
C ----------------------------------------------------------------
c1.  Bloch phase.  For translation vectors T, it's sum_T exp(+i q T)
c
c2.  Methfessel's definitions of Hankels and Bessel functions:
c
c    h_0 = Re e^(ikr)/r and j = sin(kr)/kr, k = sqrt(e), Im k >=0.
c    H_L = Y_L(-grad) h(r);   J_L = E^(-l) Y_L (-grad) j(r)
c
c    They are related to the usual n_l and j_l by factors (I think)
c       H_L  =  (i k)^(l+1) n_l (kr) Y_L   (E < 0)
c       J_L  =  (i k)^(-l)  j_l (kr) Y_L   (E < 0)
c
c   which explains how the energy-dependence is extracted out.
c   Also cases e .ne. 0 and e .eq. 0 have the same equations.
c
Cr Expansion Theorem: H_{RL}(r) = H_L(r-R)
Cr   H_{RL}(E,r) = J_{R'L'}(E,r) * S_{R'L',RL}
Cr   S_R'L',RL = 4 pi Sum_l" C_{LL'L"} (-1)^l (-E)^(l+l'-l")/2 H_L"(E,R-R')
c ---
      implicit none
C ... Passed parameters
      integer mode,ndim,nlma,nlmh
      integer indxcg(1),jcg(1),nkd,nkq
      double precision p(3),q(3),alat,awald,vol,e,cg(1),
     .  dlv(1),qlv(1)
      double complex s(ndim,nlmh),sd(ndim,nlmh)
C ... Local parameters
      integer lmxx,nrxmx,nlm0
c      parameter (lmxx=12,nlm0=(lmxx+1)**2,nrxmx=2000)
      double precision fpi,p1(3),sp
c     .  wk((lmxx*2+10)*nrxmx),yl(nrxmx*(lmxx+1)**2),efac(0:lmxx)
      real(8),allocatable :: wk(:),yl(:),efac(:)
      complex(8),allocatable :: dl(:),dlp(:)
      integer(4),allocatable :: sig(:)
      double complex phase,sumx,sud !dl(nlm0),dlp(nlm0)
      integer icg,icg1,icg2,ii,indx,ipow,l,lmax,nrx,nlm,
     .  ilm,ilma,la,ilmb,lh !sig(0:lmxx),
      logical ldot
      integer(4) :: job
#ifdef COMMONLL
      integer(4):: ll(51**2)
      common/llblock/ll
#else
      integer(4)::ll
#endif

C ... Setup
      ldot = .false.
      lmax = ll(nlma)+ll(nlmh)
      nlm = (lmax+1)**2
      nrx  = max(nkd,nkq)
      fpi  = 16d0*datan(1d0)
c
      if (nlma > ndim) call rxi('strxq: increase ndim: need',nlma)
c      if (lmax .gt. lmxx) call rxi('strxq: increase lmxx: need',lmax)
c      if (nrx .gt. nrxmx) call rxi('strxq: increase nrxmx: need',nrx)
      lmxx = lmax
      nlm0 =(lmxx+1)**2
      nrxmx= nrx
      allocate( wk((lmxx*2+10)*nrxmx),yl(nrxmx*(lmxx+1)**2),
     &  efac(0:lmxx),sig(0:lmxx),dl(nlm0),dlp(nlm0))

c ???
c      allocate( wk((lmax*2+10)*nrxmx),yl(nrxmx*(lmax+1)**2),
c     &  efac(0:lmax),sig(0:lmax),dl(nlm0),dlp(nlm0))

C --- Reduced structure constants ---
      call shortn(p,p1,dlv,nkd)
      sp = fpi/2*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
      phase = dcmplx(dcos(sp),dsin(sp))
      job = 0
      if( sum(abs(p))<1d-10 ) job = 10
c      print *," goto hsm e=",e
      if (e .lt. 0) then
        call hsmq(1,0,lmax,e,0d0,job,q,p1,nrx,nlm0,wk,yl,
     .    awald,alat,qlv,nkq,dlv,nkd,vol,dl,dlp)
      else
        call hsmqe0(lmax,0d0,job,q,p1,nrx,nlm0,wk,yl,
     .    awald,alat,qlv,nkq,dlv,nkd,vol,dl)
        ldot = .false.
      endif
c      print *," end of hsm"

C ... Put in phase to undo shortening
      if (sp .ne. 0d0) then
        do  20  ilm = 1, nlm
          dl(ilm) = phase*dl(ilm)
          if (ldot) dlp(ilm) = phase*dl(ilm)
   20   continue
      endif

c      print *," end of do 20"
C --- Combine with Clebsch-Gordan coefficients ---
C ... efac(l)=(-e)**l; sig(l)=(-)**l
      efac(0) = 1
      sig(0) = 1
      do  1  l = 1, lmax
        efac(l) = -e*efac(l-1)
        sig(l) = -sig(l-1)
    1 continue
c      print *," goto do11"
      do  11  ilma = 1, nlma
        la = ll(ilma)
        do  14  ilmb = 1, nlmh
          lh = ll(ilmb)
          ii = max0(ilma,ilmb)
          indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          sumx = 0d0
          sud = 0d0
c      print *," xxx 2"
c      print *," xxx =",ilma,ilmb,ll(ilma),ll(ilmb)
          if (ldot) then
            do  16  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sumx = sumx + cg(icg)*efac(ipow)*dl(ilm)
              sud = sud + cg(icg)*efac(ipow)*(dlp(ilm)+ipow*dl(ilm)/e)
   16       continue
          else
            do  15  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sumx  = sumx + cg(icg)*efac(ipow)*dl(ilm)
   15       continue
          endif
c      print *," xxx 3"
          s(ilma,ilmb) = fpi*sig(lh)*dconjg(sumx)
          if (ldot) sd(ilma,ilmb) = fpi*dconjg(sud)*sig(lh)
   14   continue
   11 continue
      if (allocated(wk))deallocate(wk)
      if (allocated(yl))deallocate(yl)
      if (allocated(efac))deallocate(efac)
      if (allocated(sig))deallocate(sig)
      if (allocated(dl))deallocate(dl)
      if (allocated(dlp))deallocate(dlp)
      end

      subroutine shortn(p,p1,dlat,nkd)
C- Get p1 = shortest vector such that p1-p is a lattice vector.
C ----------------------------------------------------------------
Ci Inputs
Ci   p     :vector to shorten
Ci   dlat  :lattice vectors, sorted by increasing length
Ci   nkd   :number of dlat
Co Outputs
Co   p1    :shortened vector
Cr Remarks
Cr   A slightly skewed norm is used to make result unique.
Cr   The first vector in the list must be the zero vector.
c ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkd
      double precision p(3),p1(3),dlat(3,nkd)
C ... Local parameters
      integer irep,k0,k
      double precision anrm2,x,y,z,p2,critk0,dd,crit
      anrm2(x,y,z) = x*x*1.00001d0 + y*y*1.00002d0 + z*z*1.00003d0
     .                -x*0.000004d0 - y*0.000003d0 - z*0.000002d0

      p1(1) = p(1)
      p1(2) = p(2)
      p1(3) = p(3)
      do  88  irep = 1, 20
        p2 = anrm2(p1(1),p1(2),p1(3))
        k0 = 1
        critk0 = p2
        do  52  k = 1, nkd
          dd = dlat(1,k)**2 + dlat(2,k)**2 + dlat(3,k)**2
C         Only valid if dlat is sorted:
          if (dd .gt. p2*4d0) goto 53
          crit = anrm2(p1(1)+dlat(1,k),p1(2)+dlat(2,k),p1(3)+dlat(3,k))
          if (crit .lt. critk0) then
            k0 = k
            critk0 = crit
          endif
  52    continue
  53    if (k0 .eq. 1) return
        p1(1) = p1(1) + dlat(1,k0)
        p1(2) = p1(2) + dlat(2,k0)
        p1(3) = p1(3) + dlat(3,k0)
   88 continue
      call rx('shortn: shortest vector not found')
      end
