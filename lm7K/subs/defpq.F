      subroutine defpq(z,lmax,nsp,p,q)
C- Generate default values for P and Q
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   lmax  :l-cutoff
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Inputs/Outputs (see Remarks)
Cio  P     :Methfessel's potential function defining the logarithmic
Cio        :derivative: if D is the logarithmic derivative,
Cio        :P  = .5 - atan(D)/pi + (princ.quant.number)
Cio        :On output, if the initial p is zero, it is set to a default
Cio  Q     :sphere charges by l
Cr Remarks
Cr   defpq chooses default values for any P and Q not specified.
Cr   * For each l , if input P_l is zero, a default value is selected.
Cr   * If the sum of Q is nonzero, Q is not altered.  Otherwise
Cr     a set of default Q are selected depending on the input P.
Cu Updates
Cu   24 Oct 07 Handle NULL input values
Cu   02 Mar 01 Handles fractional Z case
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmax,nsp,n0
      parameter (n0=10)
      double precision z,p(n0,2),q(n0,2)
C ... Local parameters
      logical lchange
      integer l,k0(n0),q0(n0),q02(n0),iz,iz2,kl,lmx,lgunit,iprint,lmin !,k02(n0)
      double precision pdef(0:3,0:100),dasum,w1,w2
      integer NULLI
      parameter (NULLI=-99999)

      data ((pdef(l,iz),l=0,3),iz=0,9) /
     .1.50D+00,  2.50D+00,  3.50D+00,  4.50D+00,
     .1.76D+00,  2.30D+00,  3.16D+00,  4.11D+00,
     .1.91D+00,  2.16D+00,  3.12D+00,  4.10D+00,
     .2.63D+00,  2.46D+00,  3.18D+00,  4.11D+00,
     .2.80D+00,  2.52D+00,  3.21D+00,  4.13D+00,
     .2.77D+00,  2.64D+00,  3.21D+00,  4.12D+00,
     .2.89D+00,  2.71D+00,  3.18D+00,  4.12D+00,
     .2.92D+00,  2.79D+00,  3.18D+00,  4.11D+00,
     .2.92D+00,  2.86D+00,  3.18D+00,  4.11D+00,
     .2.92D+00,  2.89D+00,  3.18D+00,  4.10D+00/
      data ((pdef(l,iz),l=0,3),iz=10,19) /
!     .3.17D+00,  2.92D+00,  3.12D+00,  4.10D+00, !Ne
     .2.92D+00,  2.92D+00,  3.12D+00,  4.10D+00, !Ne
     .3.70D+00,  3.38D+00,  3.19D+00,  4.12D+00,
     .3.76D+00,  3.50D+00,  3.22D+00,  4.12D+00,
     .3.81D+00,  3.59D+00,  3.24D+00,  4.13D+00,
     .3.88D+00,  3.68D+00,  3.23D+00,  4.13D+00,
     .3.92D+00,  3.78D+00,  3.20D+00,  4.12D+00,
     .3.94D+00,  3.84D+00,  3.18D+00,  4.11D+00,
     .3.94D+00,  3.89D+00,  3.18D+00,  4.11D+00,
     .3.94D+00,  3.92D+00,  3.18D+00,  4.10D+00, !Ar !18
     .4.70D+00,  4.36D+00,  3.24D+00,  4.12D+00/
      data ((pdef(l,iz),l=0,3),iz=20,29) /
     .4.72D+00,  4.40D+00,  3.38D+00,  4.13D+00,
     .4.69D+00,  4.40D+00,  3.53D+00,  4.13D+00,
     .4.66D+00,  4.40D+00,  3.62D+00,  4.13D+00,
     .4.65D+00,  4.40D+00,  3.68D+00,  4.13D+00,
     .4.65D+00,  4.40D+00,  3.73D+00,  4.13D+00,
     .4.65D+00,  4.41D+00,  3.77D+00,  4.13D+00,
     .4.66D+00,  4.41D+00,  3.80D+00,  4.13D+00,
     .4.66D+00,  4.41D+00,  3.83D+00,  4.13D+00,
     .4.67D+00,  4.41D+00,  3.86D+00,  4.12D+00,
     .4.69D+00,  4.42D+00,  3.88D+00,  4.12D+00/
      data ((pdef(l,iz),l=0,3),iz=30,39) /
     .4.70D+00,  4.50D+00,  3.92D+00,  4.11D+00,
     .4.75D+00,  4.59D+00,  4.20D+00,  4.12D+00,
     .4.85D+00,  4.70D+00,  4.20D+00,  4.12D+00,
     .4.90D+00,  4.78D+00,  4.20D+00,  4.12D+00,
     .4.90D+00,  4.84D+00,  4.20D+00,  4.11D+00,
!     .5.14D+00,  4.89D+00,  4.20D+00,  4.11D+00, ! Br
!     .5.16D+00,  4.92D+00,  4.20D+00,  4.10D+00, !Kr
     .4.90D+00,  4.89D+00,  4.20D+00,  4.11D+00, !Br takao feb2011, Br2 dimer requires 4s as valence.  !It can be treated LO. Then PZ=4.9 P=5.2 is required
     .4.90D+00,  4.92D+00,  4.20D+00,  4.10D+00, !Kr 36
     .5.71D+00,  5.33D+00,  4.26D+00,  4.12D+00,
     .5.73D+00,  5.36D+00,  4.40D+00,  4.13D+00,
     .5.70D+00,  5.35D+00,  4.55D+00,  4.14D+00/
      data ((pdef(l,iz),l=0,3),iz=40,49) /
     .5.67D+00,  5.35D+00,  4.64D+00,  4.14D+00,
     .5.66D+00,  5.35D+00,  4.70D+00,  4.14D+00,
     .5.65D+00,  5.35D+00,  4.75D+00,  4.14D+00,
     .5.65D+00,  5.35D+00,  4.79D+00,  4.14D+00,
     .5.65D+00,  5.36D+00,  4.82D+00,  4.14D+00,
     .5.65D+00,  5.36D+00,  4.85D+00,  4.13D+00,
     .5.65D+00,  5.35D+00,  4.88D+00,  4.13D+00,
     .5.69D+00,  5.37D+00,  4.90D+00,  4.12D+00,
     .5.80D+00,  5.47D+00,  4.93D+00,  4.11D+00,
     .5.85D+00,  5.58D+00,  5.18D+00,  4.12D+00/
      data ((pdef(l,iz),l=0,3),iz=50,59) /
     .5.88D+00,  5.70D+00,  5.20D+00,  4.12D+00,
     .5.90D+00,  5.78D+00,  5.20D+00,  4.12D+00,
     .5.90D+00,  5.84D+00,  5.20D+00,  4.12D+00,
!     .6.15D+00,  5.89D+00,  5.20D+00,  4.11D+00, !iodine
!     .6.15D+00,  5.92D+00,  5.20D+00,  4.10D+00, !Xe
     .5.90D+00,  5.89D+00,  5.20D+00,  4.11D+00, !I 53 iodine feb2011 takao
     .5.90D+00,  5.92D+00,  5.20D+00,  4.10D+00, !Xe 54
     .6.71D+00,  5.92D+00,  5.31D+00,  4.13D+00,
     .6.69D+00,  6.27D+00,  5.52D+00,  4.15D+00,
     .6.61D+00,  6.25D+00,  5.65D+00,  4.19D+00,
     .6.56D+00,  6.25D+00,  5.66D+00,  4.28D+00,
     .6.57D+00,  6.25D+00,  5.66D+00,  4.38D+00/
      data ((pdef(l,iz),l=0,3),iz=60,69) /
     .6.57D+00,  6.25D+00,  5.65D+00,  4.48D+00,
     .6.57D+00,  6.25D+00,  5.64D+00,  4.55D+00,
     .6.58D+00,  6.25D+00,  5.63D+00,  4.61D+00,
     .6.61D+00,  6.25D+00,  5.62D+00,  4.65D+00,
     .6.63D+00,  6.25D+00,  5.60D+00,  4.69D+00,
     .6.65D+00,  6.25D+00,  5.59D+00,  4.73D+00,
     .6.67D+00,  6.25D+00,  5.56D+00,  4.76D+00,
     .6.69D+00,  6.26D+00,  5.54D+00,  4.78D+00,
     .6.70D+00,  6.28D+00,  5.52D+00,  4.80D+00,
     .6.73D+00,  6.32D+00,  5.47D+00,  4.82D+00/
      data ((pdef(l,iz),l=0,3),iz=70,79) /
     .6.74D+00,  6.34D+00,  5.44D+00,  4.85D+00,
     .6.73D+00,  6.35D+00,  5.53D+00,  4.93D+00,
     .6.72D+00,  6.36D+00,  5.61D+00,  5.14D+00,
     .6.71D+00,  6.37D+00,  5.68D+00,  5.14D+00,
     .6.71D+00,  6.36D+00,  5.73D+00,  5.14D+00,
     .6.70D+00,  6.37D+00,  5.78D+00,  5.14D+00,
     .6.70D+00,  6.37D+00,  5.81D+00,  5.14D+00,
     .6.70D+00,  6.37D+00,  5.85D+00,  5.14D+00,
     .6.71D+00,  6.36D+00,  5.87D+00,  5.13D+00,
     .6.73D+00,  6.36D+00,  5.90D+00,  5.12D+00/
      data ((pdef(l,iz),l=0,3),iz=80,89) /
     .6.86D+00,  6.38D+00,  5.92D+00,  5.11D+00,
     .6.90D+00,  6.59D+00,  6.14D+00,  5.12D+00,
     .6.93D+00,  6.70D+00,  6.16D+00,  5.12D+00,
     .6.94D+00,  6.79D+00,  6.17D+00,  5.12D+00,
     .6.96D+00,  6.84D+00,  6.17D+00,  5.12D+00,
!     .7.15D+00,  6.89D+00,  6.18D+00,  5.12D+00, !At
!     .7.17D+00,  6.92D+00,  6.13D+00,  5.10D+00, !Rn
     .6.96D+00,  6.89D+00,  6.18D+00,  5.12D+00, !At
     .6.96D+00,  6.92D+00,  6.13D+00,  5.10D+00, !Rn
     .7.81D+00,  6.98D+00,  6.19D+00,  5.10D+00,
     .7.72D+00,  6.96D+00,  6.47D+00,  5.13D+00,
     .7.67D+00,  6.96D+00,  6.58D+00,  5.16D+00/
      data ((pdef(l,iz),l=0,3),iz=90,100) /
     .7.62D+00,  6.95D+00,  6.61D+00,  5.22D+00,
     .7.60D+00,  6.95D+00,  6.60D+00,  5.31D+00,
     .7.56D+00,  6.94D+00,  6.58D+00,  5.41D+00,
     .7.53D+00,  6.94D+00,  6.56D+00,  5.51D+00,
     .7.51D+00,  6.94D+00,  6.54D+00,  5.58D+00,
     .7.52D+00,  6.94D+00,  6.53D+00,  5.65D+00,
     .7.54D+00,  6.94D+00,  6.52D+00,  5.70D+00,
     .7.57D+00,  6.94D+00,  6.50D+00,  5.74D+00,
     .7.60D+00,  6.94D+00,  6.48D+00,  5.77D+00,
     .7.64D+00,  6.95D+00,  6.46D+00,  5.80D+00,
     .7.68D+00,  6.95D+00,  6.44D+00,  5.82D+00/

      call fsanrg(z,0d0,100d0,0d0,' ','atomic number z',.true.)
      iz = nint(max(z,0d0))
      if (iz .gt. z) then !takao added -1d-12 for safe.
        iz2 = max(iz-1,0)
        if (iz2 .eq. iz) iz2 = iz+1
      else
        iz2 = iz+1
      endif
      w2 = (z-iz)/(iz2-iz)
      w1 = 1-w2

      if(iprint()>59) print *,'conf:xxx iz=',iz,w1,w2,iz,iz2


      lmx = min(lmax,3)
      if (lmax+1 .gt. n0) call rx('defpq: increase n0')
      lchange = .false.

C --- Set default values for P ---
      do  l = 0, lmx
        if (p(l+1,1) .eq. 0 .or. p(l+1,1) .eq. NULLI) then
          p(l+1,1) = pdef(l,iz)
          p(l+1,nsp) = pdef(l,iz)
          lchange = .true.
        endif
C       Poke p(1) into p(nsp), if integer parts don't match
        if (int(p(l+1,1)) .ne. int(p(l+1,nsp))) then
          p(l+1,nsp) = p(l+1,1)
          lchange = .true.
        endif
      enddo
C ... Fill in the P's for l>3
      do  l = 4, lmax
        if (p(l+1,1) .eq. NULLI .or. p(l+1,1) .eq. 0) then
          p(l+1,1) = l+1 + .1d0
          p(l+1,nsp) = l+1 + .1d0
          lchange = .true.
        endif
        if (int(p(l+1,1)) .ne. int(p(l+1,nsp))) then
          p(l+1,nsp) = p(l+1,1)
          lchange = .true.
        endif
      enddo

C --- Set default values for Q ---
      lmin = lmx+1              !Initially no default values to set
C     If all values are zero, Q has not been set.  Set all q
      if (dasum(lmx+1,q(1,1),1)+dasum(lmx+1,q(1,nsp),1) .eq. 0) then
        lmin = 0
C     If encounter NULLI for any l<lmx, set all q for that l and higher
      else
        do  l = 0, lmx
          if (q(l+1,1) .ne. NULLI .and. q(l+1,nsp) .ne. NULLI) cycle
          lmin = l; exit
        enddo
      endif

C     Atomic moments and configuration
c      call atmoms(k0,q0,dble(iz))
c      call atmoms(k0,q02,dble(iz2))
      k0=0
c$$$      k02=0
c$$$      call atmoms(k02,q02,dble(iz2))
      call atmoms(k0,q0,dble(iz)) !takao think this is correct order to obtain k0. Nov29 2010

c$$$ Sakakibara found this is wrong for Na0.5Co2O4
c$$$!! Sanity check for virtual crystal approx. takao Nov29 2010.
c$$$      if(abs(z-int(z))>1d-10 .and. sum(abs(k0-k02))>0) then
c$$$        call rx('For such a case for Virtual crystal approx.'
c$$$     &   //' Core configulation should be the same. Need to set P and Q')
c$$$      endif

C ... Make default values
      do  l = lmin, lmx
        kl = int(p(l+1,1))
C       If P has P. Q. N. higher it is empty
        if (kl .gt. k0(l+1)) then
          q0(l+1)  = 0
          q02(l+1) = 0
C       If P has P. Q. N. lower it is filled
        elseif (kl .lt. k0(l+1)) then
          q0(l+1)  = 2*(2*l+1)
          q02(l+1) = 2*(2*l+1)
        endif
C       Copy default to q
        q(l+1,1)   = (w1*q0(l+1) + w2*q02(l+1))/nsp
        q(l+1,nsp) = (w1*q0(l+1) + w2*q02(l+1))/nsp
        lchange = .true.
      enddo
      do  l = lmx+1, lmax
        q(l+1,1)   = 0
        q(l+1,nsp) = 0
      enddo
  100 continue
      if (iprint() .ge. 100 .or. (iprint() .ge. 50 .and. lchange)) then
        call awrit6(' DEFPQ, z=%d: P=%n;6,2D%?#n==2#;%n;6,2D##',
     .  ' ',80,lgunit(1),z,lmax+1,p,nsp,lmax+1,p(1,2))
        call awrit6('        z=%d: Q=%n;6,2D%?#n==2#;%n;6,2D##',
     .  ' ',80,lgunit(1),z,lmax+1,q,nsp,lmax+1,q(1,2))
      endif
      end
      subroutine atmoms(k0,q0,z)
C- Generate a default atomic configuration for specified nuclear charge
C ----------------------------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Co Outputs:
Co   k0    :modified core configuration
Co   q0    :modified moments configuration
Cr Remarks:
Cr   Adapted from Stuttgart code lmto56
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0(0:3),q0(0:3)
      double precision z
C ... Local parameters
      integer iz,iq0at(0:3,0:100),l,li,na,k,rb,cs,fr
      parameter (li=3,na=11,k=19,rb=37,cs=55,fr=87)
C Intrinsic functions:
      intrinsic dble,idnint

      data ((iq0at(l,iz),l=0,3),iz=0,49)
     ./0,0, 0, 0,   1,0, 0, 0,   2,0, 0, 0,   1,0, 0, 0,   2,0, 0, 0,
     .2,1, 0, 0,   2,2, 0, 0,   2,3, 0, 0,   2,4, 0, 0,   2,5, 0, 0,
     .2,6, 0, 0,   1,0, 0, 0,   2,0, 0, 0,   2,1, 0, 0,   2,2, 0, 0,
     .2,3, 0, 0,   2,4, 0, 0,   2,5, 0, 0,   2,6, 0, 0,   1,0, 0, 0,
     .2,0, 0, 0,   2,0, 1, 0,   2,0, 2, 0,   2,0, 3, 0,   1,0, 5, 0,
     .2,0, 5, 0,   2,0, 6, 0,   2,0, 7, 0,   2,0, 8, 0,   1,0,10, 0,
     .2,0,10, 0,   2,1,10, 0,   2,2,10, 0,   2,3,10, 0,   2,4,10, 0,
     .2,5,10, 0,   2,6,10, 0,   1,0, 0, 0,   2,0, 0, 0,   2,0, 1, 0,
     .2,0, 2, 0,   1,0, 4, 0,   1,0, 5, 0,   2,0, 5, 0,   1,0, 7, 0,
     .1,0, 8, 0,   0,0,10, 0,   1,0,10, 0,   2,0,10, 0,   2,1,10, 0/
      data ((iq0at(l,iz),l=0,3),iz=50,99)
     ./2,2,10, 0,   2,3,10, 0,   2,4,10, 0,   2,5,10, 0,   2,6,10, 0,
     .1,0, 0, 0,   2,0, 0, 0,   2,0, 1, 0,   2,0, 1, 1,   2,0, 1, 2,
     .2,0, 1, 3,   2,0, 1, 4,   2,0, 1, 5,   2,0, 0, 7,   2,0, 1, 7,
     .2,0, 1, 8,   2,0, 1, 9,   2,0, 1,10,   2,0, 1,11,   2,0, 1,12,
     .2,0, 0,14,   2,0, 1,14,   2,0, 2,14,   2,0, 3,14,   2,0, 4,14,
     .2,0, 5,14,   2,0, 6,14,   2,0, 7,14,   1,0, 9,14,   1,0,10,14,
     .2,0,10,14,   2,1,10,14,   2,2,10,14,   2,3,10,14,   2,4,10,14,
     .2,5,10,14,   2,6,10,14,   1,0, 0, 0,   2,0, 0, 0,   2,0, 1, 0,
     .2,0, 2, 0,   2,0, 1, 2,   2,0, 1, 3,   2,0, 1, 4,   2,0, 0, 6,
     .2,0, 0, 7,   2,0, 1, 7,   2,0, 0, 9,   2,0, 0,10,   2,0, 0,11/

      iz=idnint(z)

      call icopy(4,iq0at(0,iz),1,q0,1)

      if     (iz .ge. fr) then
        k0(0) = 7
        k0(1) = 7
        k0(2) = 6
        k0(3) = 5
      elseif (iz .ge. cs) then
        k0(0) = 6
        k0(1) = 6
        k0(2) = 5
        k0(3) = 4
      elseif (iz .ge. rb) then
        k0(0) = 5
        k0(1) = 5
        k0(2) = 4
        k0(3) = 4
      elseif (iz .ge. k ) then
        k0(0) = 4
        k0(1) = 4
        k0(2) = 3
        k0(3) = 4
      elseif (iz .ge. na) then
        k0(0) = 3
        k0(1) = 3
        k0(2) = 3
        k0(3) = 4
      elseif (iz .ge. li) then
        k0(0) = 2
        k0(1) = 2
        k0(2) = 3
        k0(3) = 4
      else
        k0(0) = 1
        k0(1) = 2
        k0(2) = 3
        k0(3) = 4
      endif

      end

