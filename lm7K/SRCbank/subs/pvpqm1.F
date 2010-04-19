      subroutine pvpqm1(mode,sspec,nclass,nclasp,nsp,nbas,lpgf,lves,
     .  lqpp,leula,neul,ics,ipc,cnst,nqpp,qppo,qppn,vin,vout,eulao,
     .  eulan,vcnsto,vconst,n,xold,xnew)
C- Creates list of extra parameters to be included in self-consistency
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy to mix array
Ci         :1 copy from mix array
Ci   sspec :struct for species-specific information; see routine uspec
Ci   nclass:number of inequivalent classes
Ci   nclasp:number of classes in padded basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   lpgf  :modes and flags for layer GF
Ci   lves  :treatment of electrostatics
Ci   lqpp  :include nonspherical density coefficients in mix
Ci   leula :if include Euler angles in mixing
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   cnst  :work array holding constaints
Ci   nqpp  :number of elements for nonspherical charge density
Ci   ... if mode=0, the following are copied to xold
Ci   qppo  :input  nonspherical density coffs
Ci   vin   :input  electrostatic potential
Ci   eulao :input  Euler angles
Ci   vcnsto :input  potential shifts (not used)
Cio Inputs/Outputs
Cio  ... if mode=0, the following are input and copied to xnew
Cio  ... if mode=1, the following are read from xnew
Cio  qppn  :output nonspherical density coffs
Cio  vout  :output electrostatic potential
Cio  eulan :output Euler angles
Cio  vconst :output potential shifts (not used)
Co   xold  :if mode=0, elements from qppo,vin,eulao, etc are copied to
Co         : a single vector xold suitable for mixing.
Co         :if mode=1, xold is not used.
Co   xnew  :if mode=0, elements from qppn,vout,eulan, etc are copied to
Co         : a single vector xnew suitable for mixing.
Co         :if mode=1, xnew is copied into these arrays.
Co   n     :total number of elements to be included in mix
Cl Local variables
Cl         :
Cr Remarks
Cr  This returns 'extra' mixing parameters into a single vector.
Cb  Bugs
Cb  Bugs: qpp should be symmetrized, and based by class
Cu Updates
Cu   08 Nov 07 (J. Xu) qpp is complex
Cu   27 May 02 mixing constraint is just 1s bit of mxcst
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical leula,lves,lqpp
      integer mode,lpgf,nsp,nbas,nclass,nqpp,nclasp,neul,cnst(0:nclasp),
     .  ipc(nbas),ics(nclasp),n
      double precision eulao(nbas,neul,3),eulan(nbas,neul,3),
     .  sspec,vcnsto(-7:nbas),vconst(-7:nbas),
     .  vin(nclasp),vout(nclasp)
      double complex qppo(nqpp,4,nsp,nbas),qppn(nqpp,4,nsp,nbas)
      double precision xold(1),xnew(1)
C ... Local parameters
      logical lcnst
      integer i,ie,ic,ib,owk,k,lgunit,iprint,bitand
      character*80 outs
C ... Heap
      integer w(1)
      common /w/ w

C ... Count how many elements are in the extra mix
      outs = ' pvpqm1: elements to include in the mixing::'
      n = 0
      if (leula) then
        n = n + 3*neul*nbas
        call awrit0('%a%b Euler-angles,',outs,80,0)
      endif
      if (lqpp) then
        n = n + nqpp*4*nsp*nbas
        call awrit0('%a%b nonspherical moments,',outs,80,0)
      endif
      if (lves) then
        n = n + nclasp
        call awrit0('%a%b l=0 potentials,',outs,80,0)
      endif
      if (mode .eq. 0 .and. n .gt. 0 .and. iprint() .ge. 31)
     .  call awrit0('%a%b ',outs,-len(outs),-lgunit(1))
C     if (lpgf .eq. 1) n = n+1
      lcnst = cnst(0) .gt. 0

C --- Copy to mix array ---
      if (mode .eq. 0) then

C   ... Make cnst
        call sp2cls('spec mxcst',sspec,ics,1,1,nclasp,owk)
        do  i = 1, nclasp
          cnst(i) =  bitand(w(owk+i-1),1)
        enddo
        cnst(0) = 1
        if (lpgf .eq. 2) then
          call ivset(cnst,1,nclass+1,1)
        elseif (lpgf .ne. 0) then
          call ivset(cnst,nclass+2,nclasp+1,1)
        endif
        call rlse(owk)

C   --- Accumulate any parameters into mixing arrays ---
        i = 0
C        if (lpgf .eq. 1) then
C          i = i+1
C          xold(i) = vcnsto(-2)
C          xnew(i) = vconst(-2)
C        endif

C   ... Euler angles
        if (leula) then
          do  10  ie = 1, neul
            do  12  ib = 1, nbas
              ic = ipc(ib)
              if (lcnst .and. cnst(ic) .ne. 0) goto 12
              xold(i+1) = eulao(ib,ie,1)
              xold(i+2) = dcos(eulao(ib,ie,2))
              xold(i+3) = eulao(ib,ie,3)
              xnew(i+1) = eulan(ib,ie,1)
              xnew(i+2) = dcos(eulan(ib,ie,2))
              xnew(i+3) = eulan(ib,ie,3)
              i = i+3
   12       continue
   10     continue
        endif

C   ... qpp for each site
        if (lqpp) then
          k = nqpp*4*nsp
          do  20  ib = 1, nbas
            call dcopy(k,qppo(1,1,1,ib),2,xold(i+1),1)
            call dcopy(k,qppn(1,1,1,ib),2,xnew(i+1),1)
            if (xold(i+1) .eq. -1) xold(i+1) = 0
            i = i+k
   20     continue
        endif

C   ... l=0 electrostatic potential
        if (lves) then
          do  30  ic = 1, nclasp
            i = i+1
            xold(i) = vin(ic)
            xnew(i) = vout(ic)
   30     continue
        endif

        if (i .ne. n) call rx('bug in pvpqm1')
C        call prmx('xold',xold,n,n,1)
C        call prmx('xnew',xnew,n,n,1)

C --- Copy from mix array ---
      else

        i = 0

C        if (lpgf .eq. 1) then
C          i = i+1
C          vconst(-2) = xnew(i)
C        endif

        if (leula) then
          do  110  ie = 1, neul
            do  112  ib = 1, nbas
              ic = ipc(ib)
              if (lcnst .and. cnst(ic) .ne. 0) goto 112
              eulan(ib,ie,1) = xnew(i+1)
              eulan(ib,ie,2) = dacos(max(min(xnew(i+2),1d0),-1d0))
              eulan(ib,ie,3) = xnew(i+3)
              i = i+3
  112       continue
  110      continue
        endif

        if (lqpp) then
          k = nqpp*4*nsp
          do  120  ib = 1, nbas
            call dcopy(k,xnew(i+1),1,qppn(1,1,1,ib),2)
            i = i+k
  120     continue
        endif

        if (lves) then
          do  130  ic = 1, nclasp
            i = i+1
            vout(ic) = xnew(i)
  130     continue
        endif

      endif

      end

