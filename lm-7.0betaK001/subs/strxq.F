      subroutine strxq(job,e,q,p,nlma,nlmh,ndim,alat,vol,awald,nkd,nkq,
     .  dlv,qlv,cg,indxcg,jcg,s,sd)
C- One-center expansion coefficents to j of Bloch summed h (strux)
C ----------------------------------------------------------------
Ci Inputs:
Ci   job   :1's digit
Ci         :0,1 calculate s only
Ci         :2   calculate sd only (not implemented)
Ci         :any other number: calculate both s and sdot
Ci         :10's digit
Ci         :0 use standard phase convention, i.e. sum_T exp(+i q T)
Ci         :1 scale standard phase convention by exp(-i q . p)
Ci   e     :energy of Hankel function.
Ci   q     :Bloch wave number
Ci   p     :position of Hankel function center;
Ci         :structure constants are for expansion about the origin
Ci   nlma  :Generate coefficients S_R'L',RL for L' < nlma
Ci   nlmh  :Generate coefficients S_R'L',RL for L  < nlmh
Ci   ndim  :leading dimension of s,sdot
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   vol   :cell volume
Ci   awald :Ewald smoothing parameter
Ci   nkd   :number of direct-space lattice vectors
Ci   nkq   :number of reciprocal-space lattice vectors
Ci   dlv   :direct-space lattice vectors, units of alat
Ci   qlv   :reciprocal lattice vectors, units of 2pi/alat
Ci   cg    :Clebsch Gordon coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L quantum number for the C.G. coefficients (scg.f)
Co Outputs
Co   s      :structure constant matrix S_R'L',RL
Co   sd     :Energy derivative of s
Cr Remarks
Cr   This code uses Methfessel's definitions for Hankels and Bessels;
Cr   see besslr.f
Cr
Cr   Expansion Theorem: H_{RL}(r) = H_L(r-R)
Cr   H_{RL}(E,r) = J_{R'L'}(E,r) * S_{R'L',RL}
Cr   S_R'L',RL = 4 pi Sum_l" C_{LL'L"} (-1)^l (-E)^(l+l'-l")/2 H_L"(E,R-R')
Cu Updates
Cu  15 Aug 00 extended to e>0; added 10s digit job
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,ndim,nlma,nlmh
      integer indxcg(1),jcg(1),nkd,nkq
      double precision p(3),q(3),alat,awald,vol,e,cg(1),dlv(1),qlv(1)
      double complex s(ndim,nlmh),sd(ndim,nlmh)
C ... Local parameters
      integer lmxx,nlm0,nrxmx,jobh
      parameter (lmxx=12,nlm0=(lmxx+1)**2,nrxmx=2000)
      double precision efac(0:lmxx),fpi,wk((lmxx*2+10)*nrxmx),
     .  yl(nrxmx*(lmxx+1)**2),p1(3),sp
      double complex dl(nlm0),dlp(nlm0),phase,sum,sud
      integer icg,icg1,icg2,ii,indx,ipow,l,lmax,ll,sig(0:lmxx),nrx,nlm,
     .  ilm,ilma,la,ilmb,lh
      logical ldot

C ... Setup
      call tcn('strxq')
      ldot = mod(job,10) .ge. 2
      jobh = 10
      if (mod(job/10,10) .ne. 0) jobh = 1010
      lmax = ll(nlma)+ll(nlmh)
      nlm = (lmax+1)**2
      nrx  = max(nkd,nkq)
      fpi  = 16d0*datan(1d0)
      if (nlma .gt. ndim) call rxi('strxq: increase ndim: need',nlma)
      if (lmax .gt. lmxx) call rxi('strxq: increase lmxx: need',lmax)
      if (nrx .gt. nrxmx) call rxi('strxq: increase nrxmx: need',nrx)

C --- Reduced structure constants ---
      call shortn(p,p1,dlv,nkd)
      sp = fpi/2*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))

      if (e .ne. 0) then
        call hsmq(1,0,lmax,e,0d0,jobh,q,p1,nrx,nlm0,wk,yl,
     .    awald,alat,qlv,nkq,dlv,nkd,vol,dl,dlp)
      else
        call hsmqe0(lmax,0d0,jobh,q,p1,nrx,nlm0,wk,yl,
     .    awald,alat,qlv,nkq,dlv,nkd,vol,dl)
        ldot = .false.
      endif

C ... Put in phase to undo shortening
      if (sp .ne. 0d0) then
        phase = dcmplx(dcos(sp),dsin(sp))
        call zscal(nlm,phase,dl,1)
        if (ldot) call zscal(nlm,phase,dlp,1)
      endif

C --- Combine with Clebsch-Gordan coefficients ---
C ... efac(l)=(-e)**l; sig(l)=(-)**l
      efac(0) = 1
      sig(0) = 1
      do  1  l = 1, lmax
      efac(l) = -e*efac(l-1)
    1 sig(l) = -sig(l-1)
      do  11  ilma = 1, nlma
        la = ll(ilma)
        do  14  ilmb = 1, nlmh
          lh = ll(ilmb)
          ii = max0(ilma,ilmb)
          indx = (ii*(ii-1))/2 + min0(ilma,ilmb)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          sum = 0d0
          sud = 0d0
          if (ldot) then
            do  16  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sum = sum + cg(icg)*efac(ipow)*dl(ilm)
              sud = sud + cg(icg)*efac(ipow)*(dlp(ilm)+ipow*dl(ilm)/e)
   16       continue
          else
            do  15  icg = icg1, icg2
              ilm  = jcg(icg)
              ipow = (la+lh-ll(ilm))/2
              sum  = sum + cg(icg)*efac(ipow)*dl(ilm)
   15       continue
          endif
          s(ilma,ilmb) = fpi*sig(lh)*dconjg(sum)
          if (ldot) sd(ilma,ilmb) = fpi*dconjg(sud)*sig(lh)
   14   continue
   11 continue

      call tcx('strxq')

      end

