      subroutine atfold(isw,adnf,nl,nsp,vmtz,pp,sspec)
C- Automatic down-folding, 2nd generation LMTO
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw:  :1 intermediate waves mapped to lower
Ci         :2 higher waves mapped to 'neglected'
Ci         :3 combination of 1 and 2
Ci         :5 all waves mapped to lower
Ci   adnf  :T automatic downfolding turned ON; see Remarks.
Ci         :F automatic downfolding turned OFF.  Initially assign
Ci          any zero values of idxdn to lower set (see Remarks)
Ci   nl    :(global maximum l) + 1
Ci   nsp   :number of spins, needed for dimensioning ppar
Ci   vmtz  :muffin-tin zero (asamad.f)
Ci   pp    :potential parameters (atomsr.f)
Ci          NB: although pp's are class-based, atfold assumes the
Ci          first nspec classes correspond to the nspec species.
Cio Inputs/Outputs
Cio   sspec :struct containing species-specific information.
Cio          Those spec->idxdn input as zero are altered; see Remarks
Cr Remarks
Cr   There is an index in idxdn for each species and each l-channel.
Cr   It indicates in which set the orbital belongs:
Cr       1 - lower
Cr       2 - intermediate
Cr       3 - higher
Cr   If any of these are 0 on input, it will be automatically
Cr   assigned to a set using the following rules, provided adnf=T.
Cr   If |V^0 - vmtz|s^2 > 1  AND  |tau - 1| < tolh1
Cr      AND  (C - enu)sqrt p^gamma  >  tolh2  ==>  higher
Cr   If |gamma - gamma^safe|  >  |P^-1 - gamma^safe|
Cr      AND  (C - enu)sqrt p^gamma  >  toli   ==>  intermediate
Cr   Otherwise ==>  lower
Cr   You can make up your own downfolding rules.
Cu Updates
Cu   10 Apr 02 Redimensionsed idxdn to accomodate larger lmxb
Cu   14 Sep 99 atfold rewrote to species-based idxdn, with new arguments
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical adnf
      integer isw,nl,nsp
      double precision pp(6,0:nl-1,nsp,1),avw,vmtz,sspec(1)
C Local parameters
      integer n0,nkap0
      parameter(nkap0=3,n0=10)
      integer is,ic,l,isp,lgunit,iprint,idxdn(0:n0*nkap0),nspec,nglob,
     .  lmxb
      double precision tolh1,tolh2,toli,Gamma,work(5),rmax,dglob
      parameter ( tolh1 = 0.05d0, tolh2 = 0.3d0, toli = 0.25d0 )
      logical hcond(0:3),icond(0:2)

      avw = dglob('avw',0d0,0)
      nspec = nglob('nspec')
C     choose arbitrarily :  spins, classes not distinguished
      isp = 1
      do  1  is = 1, nspec
      call upack('spec idxdn rmt lmxb',sspec,is,idxdn,rmax,lmxb,0)
C     Assume representative ic = is
      ic = is
C ... Any idxdn for l>lmxb is automatically '4'
      do  3  l = lmxb+1, nl-1
    3 idxdn(l) = 4
C ... If no automatic downfolding, initial zero's to 1
      if (.not. adnf) then
        do  4  l = 0, lmxb
    4   if (idxdn(l) .eq. 0) idxdn(l) = 1
        if (isw .eq. 0) goto 1
      endif

      if (iprint() .ge. 40 .and. is .eq. 1) print 333
  333 format(' Atfold: checks for automatic downfolding'/
     .  ' spc l    |g - g_s|  |P^-1 - g_s| (C-enu)s(p) |V0-vmtz|s^2',
     .  ' |tau-1|  Conditions')

      do  2  l  = 0, lmxb

C   --- V^0 ---
        work(1) = pp(2,l,isp,ic) - pp(3,l,isp,ic)**2/pp(5,l,isp,ic)

C   --- tau ---
        Gamma = pp(3,l,isp,ic)**2/pp(5,l,isp,ic)**2
        work(2) = 2*(2*l+1)**2*(2*l+3)*(avw/rmax)**(2*l+1)
        work(2) = work(2)/Gamma/rmax**2

C   --- ( C - enu ) sqrt p^gamma ---
        work(3) = (pp(2,l,isp,ic) - pp(1,l,isp,ic))
     .                 *dsqrt(pp(4,l,isp,ic))

C   --- Gamma^safe ---
        work(4) = (rmax/avw)**(2*l+1)/(4*(l+1)*(2*l+1))

C   --- Inverse unscreened potential function at enu (2nd order) ---
        work(5) =  pp(3,l,isp,ic)**2 /
     .            (pp(1,l,isp,ic) - pp(2,l,isp,ic)) + pp(5,l,isp,ic)

        icond(1) = dabs(pp(5,l,isp,ic) - work(4)) .gt.
     .             dabs(work(5) - work(4))
        icond(2) =  work(3) .gt. toli
        hcond(1) =  abs( work(1) - vmtz )*rmax**2 .gt. 1d0
        hcond(2) =  work(3) .gt. tolh2
        hcond(3) =  dabs( work(2) -  1d0 ) .lt. tolh1
        if (iprint() .ge. 40) print 334, ic, l,
     .    dabs(pp(5,l,isp,ic)-work(4)),dabs(work(5)-work(4)),work(3),
     .    dabs(work(1)-vmtz)*rmax**2,dabs(work(2)-1d0),
     .    icond(1), icond(2), idxdn(l) .eq. 3 .or. idxdn(l) .eq. 0,
     .    hcond(1), hcond(2), hcond(3),
     .    idxdn(l) .eq. 3 .or. idxdn(l) .eq. 0
  334   format(2i3,5f12.5,2x,3l1,2x,4l1)

C   ... shift to intermediate if icond(0) is true
        icond(0) = icond(1) .and. icond(2) .and. idxdn(l) .eq. 0
        icond(0) = icond(0) .or.
     .    (hcond(1) .and. hcond(2) .and. hcond(3) .and. isw .ge. 2)
        icond(0) = icond(0) .or. idxdn(l) .gt. 2 .and. isw .ge. 2
        icond(0) = icond(0) .and. mod(isw,2) .eq. 0

C   ... shift to high if hcond(0) is true
        hcond(0) = hcond(1) .and. hcond(2) .and. hcond(3) .and.
     .             idxdn(l) .eq. 0 .and. isw .lt. 2

        if (hcond(0)) then
C --- higher set ---
          idxdn(l) = 3
          if (iprint() .ge. 10) write (lgunit(1),10) ic, l
          if (iprint() .ge. 10) write (lgunit(2),10) ic, l
          goto  2
        endif
        if (icond(0)) then
C --- intermediate set ---
          idxdn(l) = 2
          if (iprint() .ge. 10) write (lgunit(1),20) ic, l
          if (iprint() .ge. 10) write (lgunit(2),20) ic, l
          goto 2
        endif
C --- lower set ---
        if ( idxdn(l) .eq. 2 .and. mod(isw,2) .eq. 1 .or.
     .    idxdn(l) .gt. 2 .and. isw .eq. 5) then
          if (iprint() .ge. 10) write (lgunit(1),30) ic, l
          if (iprint() .ge. 10) write (lgunit(2),30) ic, l
        endif
        if ( idxdn(l) .eq. 0 ) idxdn(l) = 1
    2 continue
      call pack2('spec idxdn',sspec,is,idxdn)
    1 continue
   10 format(1x,'Atfold:  species',i4,',  l = ',i1,' to higher set')
   20 format(1x,'Atfold:  species',i4,',  l = ',i1,
     .                                           ' to intermediate set')
   30 format(1x,'Atfold:  species',i4,',  l = ',i1,' to lower set')
      end

