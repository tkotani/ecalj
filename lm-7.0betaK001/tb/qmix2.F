      subroutine qmix2(lc,nbas,nsp,ltb,dclabl,ipc,it,itmax,cnvg,mmix,
     .                 nkill,nelts,beta,betav,tm,tj1,tj2,qmpol,qits,
     .                 mmom,mmom0,qdits,a1,a2,rms)
C- Mixing multpole moments for TB-L (l-independent U)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lc: switch, if T mix only l=0 monopole moments
Ci   nbas,nsp,ltb,it,itmax,cnvg,mmix,nkill,beta,betav,tm
Ci   beta and betav mix multipoles and magnetic moments respectively
Ci   qmpol: multipole moments from tbmpol
Ci   mmom: magnetic moments by site
Ci   qits, qdits, tj1, tj2, a1, a2: work arrays
Co Outputs:
Co   qmpol (and mmom) mixed, rms
Cr Remarks
Cr   If U is independent of l (but site dependent), i.e., Uav=T then
Cr   the most efficient is to mix the Q_RL and magnetic moments.
Cr
Cr   Mixing works like this:
Cr   At some iteration, i, we have a density (represented here simply
Cr   by multipoles, rather like the ASA) rho_i. This makes a hamiltonian
Cr   and the Schrodinger equation (bndtb, tbfrce) acts like a black box
Cr   to produce a function f(rho_i). The Anderson mixing consults all
Cr   previous rho_i and f(rho_i) and proposes a new density rho_i+1
Cr   which is fed back into the black box until rho_i+i is equal to
Cr   within an rms difference to rho_i. 
Cr   
Cr   qits is a work array which holds all previous rho_i, f(rho_i) and
Cr   rho_i+1, keeping just the previous mmix iterations. At each entry 
Cr   to qmix these are rolled back in the array and iteration mmix+1 is
Cr   dicarded so that after rolling the array is structured like this:
Cr
Cr                  i=0     i=1    2    3   ...
Cr   qits(..,i,1)  empty*   rho_i from previous iterations**
Cr   amix notation          x_0   x_1  x_2 ...
Cr
Cr    * mixed Q_RL will be put in here after amix
Cr   ** note i=1 is most recent, i=2 next most recent, etc.     
Cr      if it=1 then for i=1 these are zero or Q_RL from disc (restart)
Cr
Cr   qits(..,i,2)  i=0     i=1    2    3   ...
Cr                 empty*   f(rho_i) from previous iterations**
Cr   amix notation f(x_0) f(x_1) f(x_2) ...
Cr
Cr    * the most recent rho_i from tbmpol is copied into here before
Cr      mixing
Cr   ** note i=1 is most recent, i=2 next most recent, etc.
Cr
Cr   amix requires as input f(x_i) and d(x_i) = f(x_i) - x_i
Cr   hence the call to daxpy in constructing the work array, a.
Cr
Cr   If nsp=2, then qits(1,.) holds just the spin up charge, the spin
Cr   down charge is maintained in qdits and tacked onto the
Cr   end of the mixing work array, "a".
Cr
Cu Updates
Cu   qmix2 mixes the charges and the magnetic moments separately.
Cb Bugs
Cb   No doubt keeping both qits and a is redundant, maybe someone clever
Cb   can save memory by getting rid of qits and building "a" directly
Cb   from qmpol.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      logical lc
      integer nbas,nsp,ltb,it,itmax,mmix,nkill,ipc(1),nelts
      double precision qmpol(9,nbas),qits(9,nbas,0:mmix+1,2),dclabl(1),
     .                 mmom(nbas),mmom0(nbas),qdits(nbas,0:mmix+1,2),
     .                 a1(nelts,0:mmix+1,2),a2(nbas,0:mmix+1,2)
      double precision cnvg,beta,betav,tm,tj1(1),tj2(1),rms
C Local Variables
      integer neltsm,nmix1,nmix2,amix,npmix,i,ipr,iprint,
     .        ido,ib,ic,ierr,jump,i1mach
      integer onorm,okpvt
      double precision rms1,rms2,b,qmp(9),dqtot,d1mach,dabs,dsum
      character clabl*8, outs*20
      logical IO, kill, bittst, cmdopt
C Local iteration count
      integer LOCIT
      save LOCIT

C Heap
      integer w(1)
      common /w/ w

      IO = bittst(ltb,2**16)
      jump = 1
      if (lc) then
        jump = 9
      endif
      neltsm = 0
      if (nsp .eq. 2) then
        neltsm = nbas
      endif

C --- "zero'th iteration" ---
      if (it .eq. 0) then
        ierr = -1
        if (IO) then
          call ioqm(nsp,9*nbas,neltsm,qmpol,mmom,ierr)
        endif
        if (ierr .eq. -1) then
          call dcopy(9*nbas,0d0,0,qmpol,1)
          if (nsp .eq. 2) then
            call dcopy(nbas,mmom0,1,mmom,1)
          endif
        endif
        return
      endif

C --- kill "mix files" according to nkill or rms (see parmxp) ---
      kill = .true.
      if (cmdopt('--nomixkill',11,0,outs)) then
        kill = .false.
      endif

      if (it .eq. 1) then
        LOCIT = 0
      endif
      if ( ( nkill .lt. 0 .or. rms .lt. 0d0 .or.
     .     ( nkill .gt. 0 .and. mod(it,nkill) .eq. 0 ) )
     .    .and. kill ) then
        LOCIT = 1
      else
        LOCIT = LOCIT + 1
      endif
      npmix = min(LOCIT-1,mmix)

C --- check for charge neutrality ---
      dqtot = dsum(nbas,qmpol,9)
      if (dabs(dqtot) .gt. d1mach(3)) then
        if (iprint() .gt. 40 .or. (iprint() .ge. 30 .and.
     .      dabs(dqtot) .gt. 1d-4)) then
          call awrit1(' QMIX: input qtot=%;2e',' ',120,i1mach(2),dqtot)
        endif
      endif

C --- Roll back previous iterations ---
      do  i = mmix, 0, -1
        call dcopy(nelts,qits(1,1,i,1),jump,qits(1,1,i+1,1),jump)
        call dcopy(nelts,qits(1,1,i,2),jump,qits(1,1,i+1,2),jump)
        call dcopy(neltsm,qdits(1,i,1),1,qdits(1,i+1,1),1)
        call dcopy(neltsm,qdits(1,i,2),1,qdits(1,i+1,2),1)
      enddo

C --- Copy new Q_RL and mmom for this iteration from input qmpol, mmom
      if (it .eq. 1) then
        ierr = -1
        if (IO) then
          call pshprt(0)
          call ioqm(nsp,9*nbas,neltsm,qits(1,1,1,1),qdits(1,1,1),ierr)
          call popprt
        endif
        if (ierr .eq. -1) then
          call dcopy(9*nbas,0d0,0,qits(1,1,1,1),1)
          call dcopy(neltsm,mmom0,1,qdits(1,1,1),1)
        endif
      endif
      call dcopy(9*nbas,qmpol,1,qits(1,1,0,2),1)
      call dcopy(neltsm,mmom,1,qdits(1,0,2),1)

C --- Build work arrays for amix ---
      do  i = 0, npmix
        call dcopy(nelts,qits(1,1,i+1,1),jump,a1(1,i,2),1)
        call dcopy(nelts,qits(1,1,i,2),jump,a1(1,i,1),1)
        call dcopy(neltsm,qdits(1,i+1,1),1,a2(1,i,2),1)
        call dcopy(neltsm,qdits(1,i,2),1,a2(1,i,1),1)
        if (i .ne. 0) then
          call daxpy(nelts,-1d0,a1(1,i,2),1,a1(1,i,1),1)
          call daxpy(neltsm,-1d0,a2(1,i,2),1,a2(1,i,1),1)
        endif
      enddo

C --- Mix; don't chatter about it ---
      call pshprt(0)
      ipr = iprint()
      ido = 0
      call defdr(onorm,mmix**2)
      call defi(okpvt,mmix)
      b = beta
      nmix1 = amix(nelts,npmix,mmix,ido,b,ipr,tm,w(onorm),w(okpvt),
     .             a1,tj1,rms1)
      if (b .gt. 1d-4) then
        rms = rms1
      endif
      if (nsp .eq. 2) then
        b = betav
        nmix2 = amix(neltsm,npmix,mmix,ido,b,ipr,tm,w(onorm),w(okpvt),
     .               a2,tj2,rms2)
        if (betav .gt. 1d-4 .and. beta .lt. 1d-4) then
          rms = rms2
        endif
      endif
      call popprt
      call rlse(onorm)

C --- Get new Q_RL (and mmom) from work array ---
      call dcopy(nelts,a1(1,0,2),1,qits(1,1,0,1),jump)
      call dcopy(9*nbas,0d0,0,qmpol,1)
      call dcopy(nelts,a1(1,0,2),1,qmpol,jump)
      call dcopy(neltsm,a2(1,0,2),1,qdits(1,0,1),1)
      call dcopy(neltsm,a2(1,0,2),1,mmom,1)

C --- check for charge neutrality ---
      dqtot = dsum(nbas,qmpol,9)
      if (dabs(dqtot) .gt. d1mach(3)) then
        if (iprint() .gt. 40 .or. (iprint() .ge. 30 .and.
     .      dabs(dqtot) .gt. 1d-4)) then
          call awrit1(' QMIX: adding q=%;2e to conserve charge',' ',
     .      120,i1mach(2),-dqtot)
        endif
        dqtot = dqtot / nbas
        call daxpy(nbas,dqtot,-1d0,0,qmpol,9)
      endif

C --- write moments to disc ---
      if (it .eq. itmax .or. rms .lt. cnvg) then
        if (iprint() .gt. 40) then
          print *, ' '
          print *,' QMIX: writing moments to disc..'
        endif
        ierr = 1
        call ioqm(nsp,9*nbas,neltsm,qmpol,mmom,ierr)
      endif

C --- Printout ---
      if (iprint() .lt. 10) return
      print 100
      if (beta .gt. 1d-6) then
        call awrit6(
     .  ' Iteration %,3i: %,4i charges; mixed %i of %i, beta=%#4d '
     .  //'rms diff: %g',' ',90,i1mach(2),it,nelts,nmix1,npmix,beta,
     .  rms1)
        if (nmix1 .gt. 0) write (*,110) (tj1(i),i=1,nmix1)
      endif
      if (nsp .eq. 2 .and. betav .gt. 1d-6) then
        if (beta .gt. 1d-6) then
          call awrit5('                %,4i spins;'//
     .                '   mixed %i of %i, beta=%#4d '
     .    //'rms diff: %g',' ',90,i1mach(2),neltsm,nmix2,npmix,
     .    betav,rms2)
        else
          call awrit6(
     .  ' Iteration %,3i: %,4i spins;   mixed %i of %i, beta=%#4d '
     .  //'rms diff: %g',' ',90,i1mach(2),it,nelts,nmix1,npmix,betav,
     .  rms2)
        endif
        if (nmix2 .gt. 0) write (*,110) (tj2(i),i=1,nmix2)
      endif
      if (iprint() .lt. 40) return
      do  ib = 1, nbas
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
        call awrit1(' Atom %i '//clabl//'%cmultipole moments:',
     .        ' ',180,i1mach(2),ib)
        call dcopy(9,qits(1,ib,1,1),1,qmp,1)
        qmp(1) = qmp(1) + qdits(ib,1,1)
        call awrit3('        Q(in) %d, %3:1d, %5:1d',' ',180,
     .               i1mach(2),qmp,qmp(2),qmp(5))
        call dcopy(9,qits(1,ib,0,2),1,qmp,1)
        qmp(1) = qmp(1) + qdits(ib,0,2)
        call awrit3('       Q(out) %d, %3:1d, %5:1d',' ',180,
     .               i1mach(2),qmp,qmp(2),qmp(5))
        call dcopy(9,qits(1,ib,0,1),1,qmp,1)
        qmp(1) = qmp(1) + qdits(ib,0,1)
        call awrit3('     Q(mixed) %d, %3:1d, %5:1d',' ',180,
     .               i1mach(2),qmp,qmp(2),qmp(5))
        if (nsp .eq. 2) then
          call awrit3('    magnetic moment in %d, out %d, mixed %d',' ',
     .                120,i1mach(2),qits(1,ib,1,1)-qdits(ib,1,1),
     .                qits(1,ib,0,2)-qdits(ib,0,2),mmom(ib))
        endif
      enddo
  100 format(' QMIX mixing multipole moments:')
  110 format(28x,' t_j :',10f8.4)
      end
