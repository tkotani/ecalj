      subroutine shoctl(sctrl,sspec,spot,sarray,iopt,ifi)
C- Writes P, Q to stdout and log in format suitable for ctrl file
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nl nspin
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa
Ci     Stored:
Ci     Passed to:
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opnu oqnu oves
Ci     Stored:
Ci     Passed to:
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: ohave oics oclabl nclasp
Ci     Stored:
Ci     Passed to:
Ci   iopt  :1s digit   0 -> don't print out ves(rmt)
Ci         :          >0 -> do    print out ves(rmt)
Ci         :10s digit  1 spin-average moments
Ci         :10s digit  2 spin-flip moments
Ci         :100s digit 0 -> don't print out enu
Ci         :          >0 -> do    print out enu
Ci   ifi   :file logical unit
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Sep 04 Optionally writes enu
Cu   23 Jan 03 Extended options to average spins
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,iopt
      double precision spot(1),sarray(1),sspec(1),sctrl(1)
C ... Local parameters
      logical swves,swenu
      integer nclasp,nl,nsp
      integer oics,oinitc,oclabl,olmx,opnu,oqnu,oves,opp
C ... Heap
      integer w(1)
      common /w/ w

      swves = mod(iopt,10) .ne. 0
      swenu = mod(iopt/100,10) .ne. 0

      call upack2('ctrl nl nspin',sctrl,nl,nsp)
      CALL upack('array ohave oics oclabl nclasp',sarray,oinitc,
     .  oics,oclabl,nclasp,0)
      call upack('pot opnu oqnu oves opp',spot,opnu,oqnu,oves,opp,0)

      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmx)

      if (mod(iopt/10,10) .ne. 0) then
        call spinav(mod(iopt/10,10)-1,nclasp,nl,nsp,w(opnu),w(oqnu))
      endif

      call pvshoc(nclasp,w(oclabl),nl,nsp,w(olmx),w(opnu),w(oqnu),
     .  w(opp),swves,swenu,w(oves),w(oinitc),ifi)
      call rlse(olmx)

      end
      subroutine pvshoc(nclass,dclabl,nl,nsp,lmx,pl,ql,pp,swves,swenu,
     .  ves,initc,ifi)
C ----------------------------------------------------------------
Ci Inputs
Ci   clabl,nl,nsp,lmx,pnu,emom
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical swves,swenu
      integer nclass,nl,nsp,lmx(0:nclass),ifi,initc(0:nclass)
      double precision pl(0:nl-1,nsp,0:nclass),pp(6,0:nl-1,nsp,0:1),
     .  ql(3,0:nl-1,nsp,0:nclass),ves(0:1),dclabl(0:1)
C Local parameters
      integer ic,i,l,isp
      character clabl*8,strn*40

      do  10  ic = 0, nclass-1
        call r8tos8(dclabl(ic),clabl)
        if (mod(initc(ic),2) .ne. 1) call info0(10,0,0,
     .    ' shoctl (warning) no moments for class '//clabl)
        if (lmx(ic) .ge. 3)
     .  write(ifi,23) clabl, ((pl(l,isp,ic), l=0,lmx(ic)),isp=1,nsp)
        if (lmx(ic) .eq. 2)
     .  write(ifi,22) clabl, ((pl(l,isp,ic), l=0,lmx(ic)),isp=1,nsp)
        if (lmx(ic) .eq. 1)
     .  write(ifi,21) clabl, ((pl(l,isp,ic), l=0,lmx(ic)),isp=1,nsp)
        if (lmx(ic) .eq. 0)
     .  write(ifi,21) clabl, ((pl(l,isp,ic), l=0,0),isp=1,nsp)
   23   format('        ATOM=',A8,' P=',20(4f11.7:,/24x))
   22   format('        ATOM=',A8,' P=',20(3f11.7:,/24x))
   21   format('        ATOM=',A8,' P=',20(2f11.7:,/24x))
        write(ifi,20) (((ql(i,l,isp,ic), i=1,3), l=0,lmx(ic)),isp=1,nsp)
   20   format(20x,'  Q=',20(3f11.7:,/24x))
        call awrit1('(20x,''ENU='',20(%if11.7:/24x))',
     .    strn,len(strn),0,lmx(ic)+1)
        if (swenu) write(ifi,strn)
     .    ((pp(1,l,isp,ic), l=0,lmx(ic)),isp=1,nsp)
        if (swves) write(ifi,43) ves(ic)
   10 continue
   43 format(20x,'  V=',f11.7)
      end

