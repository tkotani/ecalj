      subroutine vtoq(nbas,nl,nsp,nclass,ipc,nrclas,emad,frzves,
     .  lmx,clabl,vrmax,rhrmx,rmax,ves,mad,z,pnu,wk,dq,qnu)
C- Adjusts moments to conform to Madelung potential or energy
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas,nl,nsp,nclass,ipc,nrclas,rmax,dq,ves,mad,z,pnu
Ci   emad when frzves = 1
Ci   ves, when frzves = 2
Co Outputs
Co   wk:  inverse M of Madelung matrix
Co   dq:  total charges in each sphere
Co   qnu: adjusted moments
Co   ves  when frzves = 1
Cr Remarks
Cr   Case frzves=1:
Cr     Idea is to readjust charges, restoring to
Cr     given Madelung energy.  Minimizing (DQ)**2 subject to
Cr     constraints dE = 2 * V . DQ and sum DQ = 0 gives to first order
Cr     DQ = DE (V - Vbar) / (V - Vbar) . (V - Vbar)
Cr     DQ is calculated repeatedly until given emad is found
Cr   Case frzves=2:
Cr     is more stringent than frzves=1 because each total sphere
Cr     charge is determined by the input Madelung potential
Cr     (but shifted by a constant to make system neutral)
Cr   Case frzves=3:
Cr     calculates dq as in frzves=2, but qnu untouched
Cr   Case frzves=4:
Cr     calculates wk(j) such that sum wk(j) V(j) = total charge
Cr   Uses old ipc
C ----------------------------------------------------------------
C     implicit none
C Passed Parameters
      character*(*) clabl
      integer nbas,nl,nclass,nsp,ipc(1),nrclas(1),frzves,lmx(1)
      double precision ves(nclass),mad(nbas,nbas),wk(nclass,nclass),
     .  qnu(3,nl,nsp,nclass),pnu(nl,nsp,nclass),rmax(1),dq(1),z(1),
     .  emad,vrmax(1),rhrmx

C Heap:
      integer w(1)
      common /w/ w

C Local Parameters
Cdelw       integer ib, ic, iclbas, oiwk, odwk, n, ib0, ic0, il, j,
Cdelw      .  oqc, odelq, isp, iprint
       integer:: ib , ic , iclbas , n , ib0 , ic0 , il , j , isp , iprint 
       integer ,allocatable :: iv_w_oiwk(:)
       real(8) ,allocatable :: rv_w_odwk(:)
       real(8) ,allocatable :: rv_w_oqc(:)
       real(8) ,allocatable :: rv_w_odelq(:)

      double precision det,qc,delq,qtot,sumqsq,syschg,vconst,sum,
     .  emad0,vmtz(2),vdotv,vbar,trumad

      if (iprint() .ge. 20) print *
      goto (1,2,2,2), frzves
      return

C --- Case frzves = 1 ---
    1 continue
      if (emad .lt. 0) stop 'Vtoq: bad emad'
C --  Calculate Madelung potential and energy of this qnu  --
Cdelw       call defdr(odelq,nclass)
       allocate(rv_w_odelq(nclass)) 

Cdelw       call defdr(oqc,nclass)
       allocate(rv_w_oqc(nclass)) 

Cdelw       call getq(nsp,nl,lmx,nclass,z,pnu,qnu,0,w,w(oqc),dq,w(odelq))
       call getq ( nsp , nl , lmx , nclass , z , pnu , qnu , 0 , w , 
     . rv_w_oqc , dq , rv_w_odelq ) 

   66 continue
      call pshprt(iprint()-10)
      vmtz(1) = 0
      call madpot(nbas,1,nclass,nrclas,ipc,clabl,dq,0d0,rhrmx,rmax,
     .  mad,w,w,0d0,.false.,vrmax,ves,emad0,trumad,vmtz)
      call popprt

C --  Subtract off avg V, make (V - Vbar) . (V - Vbar) and dq --
      vbar = 0
      vdotv = 0
      do  60  ic = 1, nclass
   60 vbar  = vbar  + nrclas(ic)*ves(ic)
      vbar = vbar/nbas
      do  62  ic = 1, nclass
   62 ves(ic) = ves(ic) - vbar
      do  64  ic = 1, nclass
   64 vdotv =vdotv + nrclas(ic)*ves(ic)**2
      call daxpy(nclass,(emad-emad0)/vdotv,ves,1,dq,1)

C --  Shift charge and repeat until emad-emad = 0 --
      if (dabs(emad0-emad) .gt. 1d-8) goto 66
      goto 32


C --- Case frzves = 2,3,4 ---
    2 continue

C --  make the matrix M such that Q = M V --
      call dpzero(wk,nclass**2)
      do  10  ic0 = 1, nclass
        ib0 = iclbas(ic0-1,ipc)
        do  20  ib = 1, nbas
          ic = ipc(ib)
          wk(ic0,ic) = wk(ic0,ic) + 2*mad(ib0,ib)
   20   continue
        wk(ic0,ic0) = wk(ic0,ic0) + 2/rmax(ic0)
   10 continue

Cdelw       call defi(oiwk,nclass)
       allocate(iv_w_oiwk(nclass)) 

Cdelw       call defdr(odwk,nclass)
       allocate(rv_w_odwk(nclass)) 

Cdelw       call dgefa(wk,nclass,nclass,w(oiwk),n)
       call dgefa ( wk , nclass , nclass , iv_w_oiwk , n ) 

      if (n .ne. 0) stop 'Vtoq: madelung matrix singular'
Cdelw       call dgedi(wk,nclass,nclass,w(oiwk),det,w(odwk),1)
       call dgedi ( wk , nclass , nclass , iv_w_oiwk , det , rv_w_odwk 
     . , 1 ) 

Cdelw rlse name= oiwk old_list= odelq oqc oiwk odwk 
Cdelw rlse name= oiwk new_list= odelq oqc 
Cdelw       call rlse(oiwk)
       if (allocated(rv_w_odwk)) deallocate(rv_w_odwk)
       if (allocated(iv_w_oiwk)) deallocate(iv_w_oiwk)


C -- Find wk(ic) s.t. sum_ic wk(ic) v(ic) = tot Q --
      if (frzves .eq. 4) then
        do  50  ic = 1, nclass
          sum = 0
          do  52  j = 1, nclass
   52     sum = sum + wk(j,ic)*nrclas(j)
          wk(ic,1) = sum
   50   continue
Cdelw error in return 113 ['return']
Cdelw found return, but w_varlist remains: odelq oqc 
        return
      endif

C -- Find the shift in ves to make system neutral --
      call dmpy(wk,nclass,1,ves,nclass,1,dq,nclass,1,nclass,1,nclass)
      syschg = 0
      vconst = 0
      do  40  ic = 1, nclass
        sum = 0
        do  42  j = 1, nclass
   42   sum = sum + wk(ic,j)
        vconst = vconst + nrclas(ic)*sum
   40 syschg = syschg + nrclas(ic)*dq(ic)
      vconst = syschg / vconst
      do  45  ic = 1, nclass
   45 ves(ic) = ves(ic) - vconst

C -- Find the q corresponding to shifted ves --
      call dmpy(wk,nclass,1,ves,nclass,1,dq,nclass,1,nclass,1,nclass)

      if (frzves .eq. 3) goto 70

C --  Adjust moments to fit dq --
   32 continue
      if (nsp .eq. 2) print *, 'need check vtoq for nsp=2'

      if (iprint() .gt. 20) print *,
     .  'Class  New DQ     Old DQ   Difference      By l channel ...'

      do  30  ic = 1, nclass
        call getq(nsp,nl,nl-1,1,z(ic),pnu(1,1,ic),qnu(1,1,1,ic),0,w,
     .    qc,qtot,delq)

        sumqsq = 0
        do  35  isp = 1, nsp
        do  35  il = 1, nl
   35   sumqsq = sumqsq + qnu(1,il,isp,ic)**2

        if (iprint() .gt. 20)  print 334, ic, dq(ic), qtot, dq(ic)-qtot,
     .  (((dq(ic)-qtot)*qnu(1,il,isp,ic)**2/sumqsq, isp=1,nsp), il=1,nl)
  334   format(i4,7f11.6:/37x,4f11.6)

        do  37  isp = 1, nsp
          do  37  il = 1, nl
          qnu(1,il,isp,ic) = qnu(1,il,isp,ic) +
     .      (dq(ic)-qtot)*qnu(1,il,isp,ic)**2/sumqsq
   37   continue
   30 continue

   70 continue
      if (iprint() .gt. 20 .and. (frzves .ge. 1 .and. frzves .le. 3))
     .  print 335, vconst, syschg
  335 format(' Vtoq:  pot adjusted by', f9.6,
     .       ' to adjust for net charge', f9.6)
Cdelw w_varlist remains: odelq oqc 
      end

