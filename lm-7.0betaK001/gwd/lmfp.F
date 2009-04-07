      subroutine lmfp(prgnam,sbz,sctrl,sham,spot,slat,smix,
     .  sspec,ssite,sstr,smove,sarray,sstrn)
C- LM-FP self-consistency loop
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:nam of calling program
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: nevmx lmet
Ci     Stored:    nevmx
Ci     Passed to: iors bndfp
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nl nspec nspin lpgf lcgf lscr lsx zbak
Ci                maxit lrs,78+163 lbas,1 lfrce nitmv mdprm ltb tol:33
Ci     Stored:    mdprm ltb
Ci     Passed to: supot subasi suham rlxstp iors smshft bndfp nwit relax
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: ldham ehf ehk seref eterms
Ci     Passed to: subasi suham smshft bndfp
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: oorhat osmrho osmpot
Ci     Passed to: supot suham rdovfa iors bndfp
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat opos gam nsgrp oistab oag osymgr
Ci     Stored:    gam
Ci     Passed to: supot suham rdovfa iors smshft bndfp ioden
Ci   smix  :struct for charge mixing parameters; see routine umix
Ci     Elts read: b bv w wc nsave mmix
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Passed to: praugm subasi suham rdovfa iors smshft bndfp ioden
Ci                relax
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    pos pos0
Ci     Passed to: spackv suham rlxstp rdovfa iors smshft bndfp ioden
Ci                cppos relax
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci         :(not used for now)
Ci   smove :struct for dynamics information; see routine umove
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: npadl npadr
Ci     Stored:    *
Ci     Passed to: suham
Ci   sstrn :struct for global strings
Ci     Elts read: jobid
Ci     Stored:    *
Co Outputs
Co   See Remarks
Cl Local variables
Cl   lmaxu :max l for a U (used for dimensioning)
Ci w(olldau):lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat beginning at dmats(*,lldau(ib))
Cl   nlibu :total number of U blocks
Cl   irlxsh:counter for shear relaxations.  irlxsh=0 => new step
Cr Remarks
Cr  This is the entry point for the self-consistent FP band program
Cu Updates
Cu   05 Jul 08 Setup for new PW addition to basis
Cu   04 Jul 08 New restart file editor
Cu   20 Jun 06 Repackaged MPI
Cu   21 Mar 06 First cut at shear relaxations
Cu   08 Mar 06 Relaxation restores pos at minimum g when not convgd
Cu   08 Jan 06 can write to restart file rst.#, where # = iter-no.
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   07 Jul 05 rst file version 1.04
Cu   27 Apr 05 LDA+U added (Lambrecht)
Cu   26 Mar 05 Added switch --shorten=no to suppress pos shortening
Cu   23 Feb 05 Bug fix: forces correspondence betw/ pos and site->pos
Cu             after file read of positions.
Cu   11 Jan 05 energy convergence set to ehk when sigma included
Cu   21 Dec 04 Add option to rotate local density on file read
Cu             and to shorten basis vectors after file read
Cu   06 Sep 03 1st cut at automatic optimization of wave function
Cu    9 Jan 03 Undoes lattice shear when writing basis vectors
Cu   21 May 02 Writes restart file after smshft when moving atoms
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   08 Jun 01 Revised call to nwit
Cu   15 Feb 01 added density I/O; arguments to bndfp changed.
Cu   17 Jun 00 alpha version.  No relaxations yet.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters:
      character sstrn*(*), prgnam*8
      double precision sarray(1),sbz(1),sctrl(1),slat(1),smix(1),
     .  sspec(1),ssite(1),sham(1),spot(1),sstr(1),smove(1)
C ... Heap
      integer w(1)
      common /w/ w
C ... Local variables
      integer procid,master,mpipid,nproc
      logical mlog
      logical lgors,cmdopt,bittst
      integer fopna,fopng,i,ifi,igets,iors,ipr,irs(5),iscr,isw,ix(5),j,
     .  k,lcgf,leks,lfrce,lgunit,lpgf,lpnu,lrel,lrout,lsx,nbas,nat,
     .  nbaspp,ndham,nevmx,nglob,nit1,nl,npadl,npadr,nsp,nspec,numq,
     .  stdo,pdim,lsc
      integer oevl,ofrc,oftot,oorhat,opos,opos2,osmpot,osmrho,owk,op,o
      double precision plat(3,3),qlat(3,3),qbg,dgets,xv(10),fptol,umix
      character strn*120, fileid*68, alabl*8, flg*3
C For mixing.  Default parameters dmxp:
C 1(I) mixing scheme; 2 beta; 3 wc; 4,5 wgt; 6(I)nsave 7(I) mmix;
C 8(I) nkill; 9 betv; 10 rmscst.  11-20 are outputs:
C 11 rmsdel 12 rms2 13 nmix 14 actual broy 15 actual beta 16-24 tj
C 25 1 if wt ne 0, 10 if wa ne 0, 11 if all nonzero
C 27..29: hold parms for static parms block regular mixing
C 30..32: hold parms for static parms block Euler angle mixing
C 33 : Lindhard screening parameter
      double precision dmxp(33)
C ... for iterations
      logical lhf,lbin,a2bin
      integer maxit,iter
      double precision seref,etot(2),amom,qdiff,qtol,etol,alat
      equivalence (qdiff,dmxp(11))
C ... for relaxation
      logical xyzfrz(3),lshr,ltmp
      integer icom,natrlx,nvrelx,oindrx,ow,ltb,itrlx,nm,irlxsh,
     .  nitrlx,bitor,ng
      integer oistab,oag,og,oshr
      double precision mdprm(6),gam(4),gam1,bstim,rhosig,pletot(6,2),
     .  plat0(9),dist0(9),dist(9)
      parameter (nm=3)
C ... for LDA+U
      integer nlibu,lmaxu
      integer odmatu,olldau,ovorb,odmato
      double precision tolu
      data irlxsh /0/ dist0 /9*0d0/

C     parameter (T=.true., F=.false.)
#if MPI
      integer:: ierr
      include "mpif.h"
#endif

      call tcn('lmfp')
      etot(1) = 0
C     call wkfast(.false.)

C      if (cmdopt('--rdbasp',8,0,strn)) then
C        fileid = 'basp'
C        if (strn(9:12) .eq. ':fn=') then
C          fileid = strn(13:)
C        else
C        endif
C        call strip(fileid,i,j)
C        ifi = fopna(fileid(1:j),-1,0)
C        rewind ifi
C        call upack1('ctrl nspec',sctrl,nspec)
C        if (.not. ioorbp(111,2,1,nspec,sspec,k,ifi))
C     .    call rxs2('lmfp: failed to find BASIS: token in file "',
C     .    fileid(1:j),'"')
C        call fclr(' ',ifi)
C      endif

C ... MPI-specific
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
ckino. Mark says nproc <= nbas should be for MPI.
#if MPI
       call upack('ctrl nbas',sctrl,nbas,0,0,0,0)
c   mpipic() is defined in slatsm
       if ( nbas < nproc ) then
         if (procid == master ) then
          write(6,*) 
        write(6,"(a,/,a,/,a,/,a)")
     &   '# For lmf-MPI, process>NBAS is not allowed.',
     &   '# If you remove this sanity check, it will fail',
     &   '# at hsibl.F and smhsbl.F.',
     &   '# You may need a better MPI version for efficent use of CPUs.'
         endif
c        call mpi_abort(MPI_COMM_WORLD,-1,ierr)
         call mpi_finalize(ierr)
         call exit(-1)
c             ---------- '-1'  is the return code of the program
         stop
c             ---------- may be stop is unnecessary.
      endif
#endif

C -------------------------- Basis optimization -------------------
      call defdr(o,1)
    2 continue
      call rlse(o)
      call defdr(o,1)
      if (cmdopt('--optbas',8,0,strn)) then
C       call wkdbg2
C       No self-consistency
        call pack1('bz nevmx',sbz,-1)
C       No forces or dynamics
        call dpzero(mdprm,6)
        call pack2('ctrl lfrce mdprm',sctrl,0,mdprm)
        call lmfopb(strn(9:),slat,sspec,etot(1)-dgets('ham seref',sham))
C       call wkdbg2

C -------------------------- Total energy mode -------------------
      elseif (cmdopt('--etot',6,0,strn)) then
C       No forces or dynamics
        call dpzero(mdprm,6)
        call pack2('ctrl lfrce mdprm',sctrl,0,mdprm)
C       Suppress writing output density
        call lsets('ctrl lrs',sctrl,.false.,8)
        call lsets('ctrl lrs',sctrl,.false.,16)
C       Exactly one iteration
        call pack1('ctrl maxit',sctrl,1)
C       Suppress mixing of output density
        call lstra('strn mix',i,j,k)
        call ustrn(i,-j,1,k,i,4)
        sstrn(k:k+3) = 'none'
C       call upacks('strn mix',i,j)
C       print *, sstrn(i:j)
      endif

C -------------------------- Unpack & initialization -------------------
      call getpr(ipr)
      call upack('ctrl nbas nl nspec nspin',sctrl,nbas,nl,nspec,nsp,0)
      call upack('array npadl npadr',sarray,npadl,npadr,0,0,0)
      call upack('ctrl lpgf lcgf lscr lsx',sctrl,lpgf,lcgf,iscr,lsx,0)
C     call upack('ctrl lncol',sctrl,lncol,0,0,0,0)
      nat = nglob('nat')

      qbg    = dgets('ctrl zbak',sctrl)
      maxit  = igets('ctrl maxit',sctrl)
      lrel   = isw(lgors('ctrl lrel,-1',sctrl))
      lhf    = lgors('ctrl lcd,2',sctrl)
      if (lhf) maxit = 1
C     nbasp  = nbas +    npadl + npadr
      nbaspp = nbas + 2*(npadl + npadr)
      stdo   = lgunit(1)
C     stdl   = lgunit(2)
      call setcc(lrel)

      irs(1) =     igets('ctrl lrs,7',sctrl)
     .       +     8*isw(lgors('ctrl lrs,256',sctrl))
      irs(2) =     igets('ctrl lrs,8+16',sctrl)/8
      irs(3) = isw(lgors('ctrl lrs,32',sctrl))
      irs(4) = isw(lgors('ctrl lrs,64',sctrl))
      irs(5) = isw(lgors('ctrl lrs,128',sctrl))
C ... --rs=3 => always read from atom file
      if (igets('ctrl lrs,3',sctrl) .eq. 3) irs(1) = 0
C ... This shouldn't be ... for backwards compatibility
C     if (lhf) irs(1) = 0

C     Sanity checks: most ASA "extras" are not implemented here
      call isanrg(lcgf, 0,0,'lmfp:','lcgf', .true.)
      call isanrg(lpgf, 0,0,'lmfp:','lpgf', .true.)
      call isanrg(iscr, 0,0,'lmfp:','lscr', .true.)
      call isanrg(lsx,  0,0,'lmfp:','lsx',  .true.)
C     call isanrg(lncol,0,0,'lmfp:','lncol',.true.)

C ... Printout properties of species
      if (ipr .ge. 30) then
C       call pr_basis (sspec,0)
        call praugm(sspec,0)
      endif

C ... Setup for no screening transformation
      if (igets('ctrl lbas,1',sctrl) .ne. 0) then

C       Shorten site positions
        call upack2('lat plat opos',slat,plat,opos)
        if (.not. cmdopt('--shorten=no',12,0,strn)) then
C       unpack from site structure
        call defrr(opos2,3*nbas)
        call spackv(10,'site pos',ssite,1,nbas,w(opos2))
        ix(1) = 2
        ix(2) = 2
        ix(3) = 2
        call info0(50,1,0,' lmfp : shortening basis vectors ... ')
        call shorps(nbas,plat,ix,w(opos2),w(opos))
        call spackv(11,'site pos',ssite,1,nbaspp,w(opos))
C       Debugging printout
C       call prmx('starting basis vectors',w(opos2),3,3,nbas)
C       call prmx('shortened basis vectors',w(opos),3,3,nbas)
        call rlse(opos2)
        endif
      endif

C ... Setup for charge mixing
      call dpzero(dmxp,33)
      call upack('mix b bv w wc',smix,dmxp(2),dmxp(9),dmxp(4),dmxp(3),0)
      dmxp(6) = igets('mix nsave',smix)
      dmxp(7) = igets('mix mmix',smix)
      call parms0(0,0,0d0,0)

C ... Allocate memory for forces
      lfrce = igets('ctrl lfrce',sctrl)
      if (lfrce .ne. 0) then
        numq = 1
        if (igets('bz lmet',sbz) .eq. 4) numq = 3
        call defrr(ofrc,3*nbas*numq)
      endif

C ... Relaxation setup
      itrlx = 1
C     nstack = 0
C     Initial shear was already folded into plat
      call upack2('lat gam opos',slat,gam,opos)
      gam1 = gam(4)
      gam(4) = 1
      call pack2('lat gam ldist',slat,gam,0)

      call upack('ctrl nitmv mdprm ltb',sctrl,nitrlx,mdprm,ltb,0,0)
      lshr = nint(mdprm(1)) .gt. 100
      if (nint(mdprm(1)) .eq. 0) nitrlx = 0
      if (nint(mdprm(1)) .gt. 0 .and. nint(mdprm(1)) .lt. 4) then
        call rx('lmf not set up for MD yet')
      endif
      if (nitrlx .gt. 0) then
        call pack1('ctrl mdprm',sctrl,mdprm)
        call pack1('ctrl ltb',sctrl,bitor(ltb,16))
        call defi(oindrx,6*nbas)
C       Next lines in case lattice relaxation
        if (lshr) then
          if (abs(gam(4)-1) .gt. 1d-10) call rx('lmfp: '//
     .        'use of SHEAR= incompatible w/ lattice relaxation')
          call upack('lat plat0',slat,plat0,0,0,0,0)
        endif
        call rlxstp(sctrl,ssite,natrlx,nvrelx,w(oindrx),xyzfrz,pdim)
        call rlse(oindrx)
        icom = 0
        if (nvrelx .ne. 0) then
          call defi(oindrx,2*natrlx)
          call defdr(ow,nvrelx*nvrelx)
          call defdr(op,pdim)
        endif
        alat = dgets('lat alat',slat)
        if (procid .eq. master) then
        ifi = fopna('bsmv',-1,0)
        call defrr(opos2,3*nbas)
        j = 1
        call ivset(ix,1,3,j)
        call shorps(nbas,plat,ix,w(opos),w(opos2))
        call iobsm0(0,bstim,0d0,0d0,nbas,alat,w(opos2),ifi)
        call rlse(opos2)
        endif
      endif

C ... Re-entry for shear distortion
      call defrr(oshr,1)
    4 continue

C ... Potential setup
      call info0(50,0,0,' lmfp : potential setup ... ')
      call supot(0,sctrl,slat,spot)

C ... Setup of hamiltonian, augmentation
      if (ipr .ge. 50) then
        call info0(50,0,0,' lmfp : basis setup ... ')
      else
        call info0(30,0,0,' ')
      endif
      call subasi(sctrl,sspec,sham)
      call suham(sctrl,slat,sspec,ssite,sarray,sham,spot,sstrn)

C ... Quit if --quit=ham given
      if (igets('ctrl quit',sctrl) .eq. 8) then
        call rx0('quit = ham')
      endif

C --- Setup for iterations in a self-consistency cycle ---
C ... Unpack or allocate some permanent arrays
      call upack('pot oorhat osmrho osmpot',spot,oorhat,osmrho,osmpot,0,
     .  0)
      call defrr (oftot,   3*nbas)
      call upack1('ham ndham',sham,ndham)
      call defrr (oevl,    ndham*2)

C ... Set various switches
C     Whether forces, and how to calculate non Helman-Feynman corr.
      lfrce = igets('ctrl lfrce',sctrl)
C     Maximum number of eigenvalues
      nevmx = igets('bz nevmx',sbz)
C     Whether to evaluate output density and/or KS energy
      lrout = 1
      leks = 1
      j = 6
      if (cmdopt('-leks=',j,0,strn)) then
        if (.not. a2bin(strn,leks,2,0,' ',j,72))
     .    call rxs('failed to parse',strn)
      endif
      if (nevmx .eq. -1) then
        lrout = 0
        leks = 0
        call pack1('bz nevmx',sbz,nevmx)
      endif
C     Whether to float pnu's
      lpnu = 1
C     Sanity checks
      if (lrout .eq. 0 .and. lfrce .ne. 0) then
        write(stdo,333) 'when forces sought'
  333   format('lmfp (fatal): output density required ',a/
     .    '      To make output density turn off HF=t and/or NEVMX<0')
        call rx('incompatible input')
      endif
C     Sanity checks
      if (lrout .eq. 0 .and. cmdopt('--etot',6,0,strn)) then
        write(stdo,333) 'with --etot switch.'
        call rx('incompatible input')
      endif
      if (lrout.eq.0 .and. maxit.gt.1 .and. ipr.ge.20) then
        call awrit1('%N lmfp (warning): %i iterations sought but no'//
     .    ' output rho ... do 1 iteration',' ',80,stdo,maxit)
        maxit = 1
      endif

C... LDA+U initialization
      call defi(olldau,-nbas)
      call defrr(ovorb,1)
      call defrr(odmatu,1)
      call defrr(odmato,1)
C     Check for LDA+U ... return nlibu > 0 if any U blocks
      call suldau(nbas,sspec,ssite,nlibu,lmaxu,w(olldau))
      if (nlibu .gt. 0) then
        call rlse(ovorb)
        i = nsp*nlibu*(lmaxu*2+1)**2
        call defcc(ovorb,-i)
        call defcc(odmatu,-i)
        call defcc(odmato,-i)
C       need group info to symmetrize site density matrix
        call upack('lat nsgrp oistab osymgr',slat,ng,oistab,og,0,0)
C       defaults
        call upack('mix umix tolu',smix,umix,tolu,0,0,0)
        if (umix .eq. 0) umix = 1
C       if (tolu .eq. 0d0) tolu = 1d-4
C       initialize density matrix for LDA+U
        call sudmtu(nbas,nsp,nlibu,lmaxu,ssite,sspec,0,w(olldau),
     .    ng,w(og),w(oistab),w(odmatu),w(ovorb))
      endif
C     end LDA+U  initialization section

C ... Invoke the restart editor
      if (cmdopt('--rsedit',8,0,strn)) then
        call rsedit(strn(9:),1,sctrl,ssite,sspec,slat,spot,sbz,
     .    nbas,nat,nspec)
        call rx0('lmfp from rsedit')
      endif

C ... Invoke the response function editor
      if (cmdopt('--chimedit',10,0,strn)) then
        call chimedit(strn(11:),1,sctrl,ssite,sspec,slat,spot,sbz,
     .    nbas,nat,nspec)
        call rx0('lmfp from chimedit')
      endif

C ---------------- Re-entry point for a new iteration ---------------
      iter = 1
    5 continue

C --- Read restart file or overlap free atom densities ---
C     irs(1) tells what to read and whether to invoke smshft.
C     4s' bit of irs(1) -> invoke smshft after file read.
C     8s' bit of irs(1) -> rotate local density after file read
C     0+1's bits irs(1)     action
C           0              read from atom file
C           1              read from binary rst file
C           2              read from ascii rsta file
C           3              read nothing (data already input)
   10 continue
C     Harris-Foulkes -> always overlap free-atom densities
      if (irs(1) .eq. 0) then
        call rdovfa(nbas,nspec,ssite,sspec,slat,spot,qbg,w(oorhat))
        nit1 = 0
      elseif (mod(irs(1),4) .ge. 1 .and. mod(irs(1),4) .le. 2) then
        lbin = .not. bittst(irs(1),2)
        k = -1
        if (procid .eq. master) then
          if (lbin) ifi = fopna('rst',-1,4)
          if (.not. lbin) ifi = fopna('rsta',-1,0)
        endif
        call mpibc1(ifi,1,2,mlog,'lmfp','ifi')
        k = iors(1,sctrl,ssite,sspec,slat,spot,sbz,fileid,nbas,nat,
     .    nspec,w(oorhat),w,nit1,lbin,ifi)
        if (k .lt. 0 .and. nproc .gt. 0) then
C         irs(1) = irs(1) - mod(irs(1),4)
          call rx('MPI: rst read failed. Restart with --rs=0')
        endif
        if (procid .eq. master) then
          call fclose(ifi)
        endif
        call mpibc1(k,1,2,mlog,'lmfp','k')

C   ... Write positions array from site structure
        call upack('lat opos',slat,opos,0,0,0,0)
        call spackv(10,'site pos',ssite,1,nbaspp,w(opos))
        call mpibc1(slat,nint(slat(1)),4,.false.,'lmfp','slat')
        if (k .lt. 0) then
          irs(1) = 0
          goto 10
        endif
        if (mod(irs(1),8) .ge. 4) then
C         If no force switch set, use default
          k = igets('ctrl lfrce',sctrl)
          if (k .eq. 0) call pack5('ctrl lfrce',sctrl,1,0,0,0,0)
          call smshft(ssite,sspec,slat,sctrl,sham,w(oorhat),w(osmrho))
C         Restore force switch
          if (k .eq. 0) call pack5('ctrl lfrce',sctrl,k,0,0,0,0)
        endif
        if (mod(irs(1),16) .ge. 8) then
C          call upack('lat dist',slat,xv,0,0,0,0)
C          call pvsms2(ssite,sspec,xv,nbas,nsp,w(oorhat))
          irs(1) = irs(1)-8
          i = igets('ctrl lrs',sctrl)
          i = i-256
          call pack5('ctrl lrs',sctrl,i,0,0,0,0)
        endif
      endif

C ... Write positions after file read, and repack
      if (ipr .ge. 50) then
      write(stdo,357) 'Basis, after reading restart file'
  357 format(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
      call dinv33(plat,1,qlat,xv)
      do  i = 1, nbas
        call upack('site spec pos',ssite,i,j,xv,0,0)
        call spacks(0,'spec name',sspec,alabl,j,j)
        call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
        write(stdo,345) i, alabl, (xv(j),j=1,3), (xv(3+j),j=1,3)
  345   format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
      enddo
      endif

C --- Optionally re-shorten basis vectors ---
      if (cmdopt('--shorps',8,0,strn)) then
        call defrr(opos2,3*nbas)
        call spackv(10,'site pos',ssite,1,nbas,w(opos2))
        ix(1) = 2
        ix(2) = 2
        ix(3) = 2
        call shorps(-nbas,plat,ix,w(opos2),w(opos))
        call info0(20,1,-1,
     .    ' lmfp  : write shortening vectors to file shorps ...')
        call iopos(.true.,1,'shorps',nbas,w(opos))
        call shorps(nbas,plat,ix,w(opos2),w(opos))
        call spackv(11,'site pos',ssite,1,nbaspp,w(opos))
C       Debugging printout
C       call prmx('starting basis vectors',w(opos2),3,3,nbas)
C       call prmx('shortened basis vectors',w(opos),3,3,nbas)
        call rlse(opos2)
      endif

C     Hang on to previous site density matrix for this iteration
      if (nlibu .gt. 0) then
        i = nsp*nlibu*(lmaxu*2+1)**2
        call zcopy(i,w(odmatu),1,w(odmato),1)
        call dpzero(w(odmatu),2*i)
      endif

C --- Make and diagonalize hamiltonian, make new charge density ---
      if (maxit .eq. 0) call info0(20,1,0,
     .  ' lmfp  : zero iterations sought ... no band pass')

      call bndfp(nbas,nsp,nlibu,lmaxu,w(olldau),
     .  ssite,sspec,slat,sctrl,sham,spot,sbz,sstrn,
     .  ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit,w(oevl),w(ofrc),
     .  w(odmatu),w(ovorb))

C ... check convergence of dmatu and update it and vorb if necessary
      if (nlibu .gt. 0 .and. maxit .gt. 0 .and. lrout .gt. 0) then
        call chkdmu(nbas,nsp,nlibu,lmaxu,ssite,sspec,sham,0,w(odmatu),
     .    w(odmato),w(ovorb),tolu,umix,w(olldau),ng,w(og),w(oistab))
      endif

C ... basis optimization: just extract etot(1) and return to opt.
      if (cmdopt('--optbas',8,0,strn)) then
        call upack2('ham ehf ehk',sham,etot(1),etot(2))
        goto 2
      endif

C --- Write smoothed charge density for contour plotting ---
      if (procid .eq. master) then
      if (cmdopt('--wden',6,0,strn)) then
        call ioden(strn(7:),slat,ssite,sspec,w(oorhat),w(osmrho))
      endif

C --- Write restart file (skip if --quit=band) ---
      if (igets('ctrl quit',sctrl) .ne. 4) then
C     Suppress saving rst file in the middle of a shear (irlxsh > 0)
      if (irs(2).gt.0 .and. (lrout.gt.0 .or. maxit .eq. 0) .and.
     .    irlxsh .eq. 0) then
C       call p_dyn_ef(s_dyn, ef0,def,esmear)
        lbin = irs(2) .ne. 2
        if (lbin) fileid = 'rst'
        if (.not. lbin) fileid = 'rsta'
        if (irs(2) .eq. 3) then
          call word(fileid,1,i,j)
          j = j+1
          fileid(j:j) = '.'
          call bin2a(' ',0,0,iter,2,0,len(fileid),fileid,j)
          if (lbin) ifi = fopng(fileid,-1,8+4)
          if (.not. lbin) ifi = fopng(fileid,-1,8)
          call info0(10,1,-1,' lmfp:  writing to restart file '//fileid)
        else
          if (lbin) ifi = fopna(fileid,-1,4)
          if (.not. lbin) ifi = fopna(fileid,-1,0)
        endif
        call upacks('strn jobid',i,j)
        fileid = 'lmfp:  ' // sstrn(i:j)
        k = iors(1,sctrl,ssite,sspec,slat,spot,sbz,
     .    fileid,nbas,nat,nspec,w(oorhat),w,iter,lbin,-ifi)
        call fclose(ifi)
      endif
      endif
      endif
      if (cmdopt('--window=',9,0,strn))
     .  call rx0('lmf : early exit (--window option)')

C --- Add to save file; decide on next iteration ---
      if (maxit .gt. 0) then
      call upack2('ham ehf ehk',sham,etot(1),etot(2))
C ... Subtract reference energy
      seref = dgets('ham seref',sham)
      etot(1) = etot(1) - seref
      if (etot(2) .ne. 0) etot(2) = etot(2) - seref
      amom    = dgets('ham eterms:15',sham)
C     The desired tolerances in q,e
      qtol = dgets('ctrl tol',sctrl)
      etol = dgets('ctrl tol:3',sctrl)
      if (procid .eq. master) then
        rhosig = dgets('ham eterms:19',sham)
        i = 0
        if (rhosig .ne. -99 .and. rhosig .ne. 0) i = 10
        call nwit(igets('ctrl nvario',sctrl),iter,maxit,
     .    lhf.or.irs(1).eq.0.and.iter.eq.1,leks+i,etol,qtol,qdiff,
     .    'cxhi',amom,etot,lsc)
      endif
      call mpibc1(lsc,1,2,mlog,'lmfp','lsc')
      if (lsc .eq. 2 .and. .not. lhf .and. maxit .gt. 1) lsc = 3
      if (lsc .eq. 1 .and. lrout .gt. 0  .or. lsc .eq. 3) then
        call query('max it=',2,maxit)
        if (iter .ge. maxit) lsc = 1
        if (iter .lt. maxit) lsc = 3
      endif
      if (igets('ctrl quit',sctrl) .eq. 4)
     .  call rx0('lmf : exit (--quit=band)')
      iter = iter+1
C     Continue iterations toward self-consistency
      if (lsc .gt. 2) then
        irs(1) = 3
        goto 5
      endif

C ... Reset quantities for iterations towards self-consistency
      if (nvrelx .gt. 0 .and. nitrlx .gt. 0) then
        iter = 1
        dmxp(11) = 0
      endif

C     call upack('lat nsgrp oistab oag osymgr',slat,ng,oistab,oag,
C    .  og,0)
C     call shoist(w(oistab),nbas,w(oag),w(og),ng)

C --- Molecular statics ---
      if (nitrlx .gt. 0 .and. lsc .le. 2) then
        call cppos(1,nbas,ssite)
C       Buglet: this eats up a little memory.  Never released.
        call defrr(opos,3*nbas)
        call spackv(10,'site pos',ssite,1,nbas,w(opos))
        call upack1('ctrl mdprm',sctrl,mdprm)
        if (lshr) then
          call grdepl(nvrelx,w(oindrx),0.01d0,etot,irlxsh,pletot,dist)
          if (irlxsh .ne. 0) then
            call grdep2(1,nvrelx,w(oindrx),dist0,dist)
            goto 98
          else
            call relax(prgnam,sctrl,ssite,sspec,itrlx,w(oindrx),natrlx,
     .                 nvrelx,pletot,w(op),w(ow),0,0d0,dist0,icom)
            call dpzero(dist,6)
            call grdep2(1,nvrelx,w(oindrx),dist0,dist)
            dist(7) = 1
          endif

        else
        call relax(prgnam,sctrl,ssite,sspec,itrlx,w(oindrx),natrlx,
     .             nvrelx,w(ofrc),w(op),w(ow),0,0d0,w(opos),icom)
        endif

C       Restore lattice symmetry to machine precision
        if (cmdopt('--fixpos',8,0,strn)) then
          call upack('lat nsgrp oistab oag osymgr',slat,ng,oistab,oag,
     .      og,0)
C         call shoist(w(oistab),nbas,w(oag),w(og),ng)
          j = 8+1
          if (strn(9:13) .eq. ':tol=') then
            j = 13
          endif
          if (strn(9:9) .ne. ':' .or.
     .      .not. a2bin(strn,fptol,4,0,' ',j,len(strn))) fptol = 1d-5
          call fixpos(w(opos),nbas,fptol,ng,plat,w(og),w(oag),w(oistab))
        endif

C       Write updated positions to bsmv file
        if (procid .eq. master .and. .not. lshr) then
          ifi = fopna('bsmv',-1,0)
          call poseof(ifi)
          bstim = bstim+1
          call defrr(opos2,3*nbas)
          j = 1
          call ivset(ix,1,3,j)
          call shorps(nbas,plat,ix,w(opos),w(opos2))
          call iobsmv(0,bstim,0d0,0d0,nbas,alat,w(opos2),-ifi)
          call rlse(opos2)
          call fclose(ifi)
        endif
C       repack updated positions in site structure
        call spackv(11,'site pos',ssite,1,nbas,w(opos))

C        if (icom .eq. 1) then
C          print *, '!!'
C        endif

C   ... Exit when relaxation converged or maximum number of iterations
        if (icom .eq. 1) then
          if (procid .eq. master) then
            k = igets('ctrl nvario',sctrl)
            flg = 'C67'
            call nwitsv(1+2,k,flg,nsp,amom,etot)
          endif
          call tcx('lmfp')
          call fexit(0,111,
     .    ' LMFP: relaxation converged after %i iteration(s)',itrlx)
        else
          call query('proceed with next relaxation step',-1,0)
        endif

C   ... Restore minimum gradient positions if this is last step
C        if (itrlx .eq. nitrlx .and. icom .eq. -1) then
        if (itrlx .eq. nitrlx) then

          if (.not. lshr) then
            call info0(20,1,0,' lmfp: restore positions for minimum g')

C           call prmx('initial positions',w(opos),3,3,nbas)
            call prelx1(1,nm,lshr,natrlx,nvrelx,w(oindrx),w(op),w(opos))
C           call prmx('minimum-g positions',w(opos),3,3,nbas)
          else
            call prelx1(1,nm,lshr,natrlx,nvrelx,w(oindrx),w(op),dist0)
            call dpzero(dist,6)
            call grdep2(1,nvrelx,w(oindrx),dist0,dist)
            call info2(20,0,0,
     .        ' lmfp : strain of minimum gradient:'//
     .        '%N   PDEF=%6;8,4D'//
     .        '%N STRAIN=%6;8,4D',
     .        dist0,dist)
          endif

C         Repack updated positions in site structure
          call spackv(11,'site pos',ssite,1,nbas,w(opos))

        endif

C   ... New density after atom shifts
C       If explicitly told to read from atom files after atom movmment
        if (igets('ctrl lrs,3',sctrl) .eq. 3) then
          irs(1) = 0

C       Else, use self-consistent
        else if (.not. lshr) then
          irs(1) = 3
          call smshft(ssite,sspec,slat,sctrl,sham,w(oorhat),w(osmrho))
        endif

C   ... Write restart file (to include new positions)
        if (procid .eq. master .and. .not. lshr) then
          ifi = fopna('rst',-1,4)
          call upacks('strn jobid',i,j)
          fileid = 'lmfp:  ' // sstrn(i:j)
          k = iors(1,sctrl,ssite,sspec,slat,spot,sbz,
     .      fileid,nbas,nat,nspec,w(oorhat),w,iter,.true.,-ifi)
          call fclose(ifi)
        endif

C   ... Write positions to file
        if (cmdopt('--wpos=',7,0,strn) .or.
     .      cmdopt('--wpos:mode1:',13,0,strn)) then
          call upack1('lat gam',slat,gam)
          call defrr(owk, 3*nbas)
          call rdistn(w(opos),w(owk),nbas,gam(1),gam(2),gam(3),1/gam1)
          call iopos(.true.,0,strn(8:),nbas,w(owk))
C         call fclr(strn(8:),-1)
          call rlse(owk)
        endif

C   ... Remove mixing file
        if (procid .eq. master) then
          call info0(20,0,0,' Delete mixing and band weights files ...')
        ifi = fopna('mixm',-1,4)
        call dfclos(ifi)
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
        endif
C       reset mixing block
        call parms0(0,0,0d0,0)

C   ... Exit when maximum number of iterations encountered
        if (itrlx .eq. nitrlx) then
          if (procid .eq. master) then
            call tcx('lmfp')
            call fexit(1,111,
     .    ' LMFP: relaxation incomplete after %i iteration(s)',nitrlx)
          else
            call tcx('lmfp')
            call fexit(1,111,' ',0)
          endif
        endif
        itrlx = itrlx+1

        if (lshr) then
          goto 98
        else
          goto 5
        endif
      endif
      endif

C ... Write positions to file
      if (cmdopt('--wpos=',7,0,strn) .or.
     .    cmdopt('--wpos:mode1:',13,0,strn)) then
        call defrr(opos,3*nbas)
        call spackv(10,'site pos',ssite,1,nbas,w(opos))
        call upack1('lat gam',slat,gam)
        call defrr(owk, 3*nbas)
        call rdistn(w(opos),w(owk),nbas,gam(1),gam(2),gam(3),1/gam1)
        if (cmdopt('--wpos:mode1:',13,0,strn)) then
          call iopos(.true.,1,strn(14:),nbas,w(owk))
C         call fclr(strn(14:),-1)
        else
          call iopos(.true.,0,strn(8:),nbas,w(owk))
C         call fclr(strn(8:),-1)
        endif
        call rlse(opos)
      endif

      call tcx('lmfp')
C     if (maxit .eq. 0) call rx0(' zero iterations sought ... quitting')
      return

C --- Setup to start calculation at new shear ---
   98 continue

      if (procid .eq. master) then
        call info0(20,0,0,' Delete mixing and band weights files ...')
        ifi = fopna('mixm',-1,4)
        call dfclos(ifi)
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
      endif

C     Restore plat, pos to their undistorted state:
C     undo original transformation = P P_0^-1
      call defrr(opos2,3*nbas)
      call dinv33(plat,0,xv,xv(10))
      call dgemm('N','N',3,nbas,3,1d0,xv,3,w(opos),3,0d0,w(opos2),
     .  3)
C     Simultaneously pack in lat->pos and site->pos
      call upack('lat opos',slat,opos,0,0,0,0)
      call dgemm('N','N',3,nbas,3,1d0,plat0,3,w(opos2),3,0d0,
     .  w(opos),3)
C     call prmx('positions for plat0',w(opos),3,3,nbas)
      call spackv(11,'site pos',ssite,1,nbas,w(opos))
      call cppos(1,nbas,ssite)
C     New shear
      call pack5('lat plat ldist dist',slat,plat0,3,dist,0,0)
C     A little memory leakage rel to 1st pass, but not so serious
      call rlse(oshr)
      call lattic(slat,sctrl,ssite,sarray)
      call defrr(oshr,1)
      call upack('lat plat',slat,plat,0,0,0,0)
C     Remake qp
      ltmp = lgors('ctrl lmet,1',sctrl) .or.
     .       lgors('ctrl ldos,4+2+1',sbz)
      call mkqp(sctrl,sbz,slat,ltmp,.false.,1,-2)

C ... Write restart file (to include new positions)
      if (procid .eq. master .and. irlxsh .eq. 0) then
        ifi = fopna('rst',-1,4)
        call upacks('strn jobid',i,j)
        fileid = 'lmfp:  ' // sstrn(i:j)
        k = iors(1,sctrl,ssite,sspec,slat,spot,sbz,
     .    fileid,nbas,nat,nspec,w(oorhat),w,iter,.true.,-ifi)
        call fclose(ifi)
      endif

C     Decide on what density to use
      if (igets('ctrl lrs,3',sctrl) .eq. 3) then
        irs(1) = 0
C     Else, use file density
      else
        irs(1) =     igets('ctrl lrs,7',sctrl)
      endif
      goto 4

      end

      subroutine cppos(ib1,ib2,ssite)
C- Copy site positions to p0 for a range of sites
C     implicit none
      integer ib1,ib2
      double precision pos(3),ssite(1)
      integer ib

      do  ib = ib1, ib2
        call spackv(0,'site pos', ssite,ib,ib,pos)
        call spackv(1,'site pos0',ssite,ib,ib,pos)
      enddo
      end

      subroutine grdepl(nvrelx,indrlx,alpha,etot,irlxsh,pletot,dist)
C-
C     implicit none
      integer irlxsh,nvrelx,indrlx(nvrelx)
      double precision pletot(6,2),etot,dist(9),alpha
      double precision grad,vec1(6)
      integer iv,ipm,ipv

      call info5(30,1,0,' GRDEPL: point %i of %i for grad shear: '//
     .  'etot=%d',irlxsh,2*nvrelx,etot,0,0)

C   3 continue

C     Get index for current shear and store energy for that shear
      if (irlxsh .gt. 0) then
        iv = (irlxsh-1)/2 + 1
        ipv = iv
C       ipv = indrlx(iv)
        ipm = mod((irlxsh-1),2) + 1
        pletot(ipv,ipm) = etot
      endif

C     If this is last point, form gradient and exit
      if (irlxsh .eq. 2*nvrelx) then
        do  iv = 1, nvrelx
C         ipv = indrlx(iv)
          ipv = iv
          grad = (pletot(ipv,1) - pletot(ipv,2))/(2*alpha)
          pletot(ipv,1) = grad
        enddo
        irlxsh = 0
        return
      endif

C     Get shear index for next shear and whether + or -
      irlxsh = irlxsh+1
      iv = (irlxsh-1)/2 + 1
      ipv = indrlx(iv)
      ipm = mod((irlxsh-1),2) + 1
      if (ipv .lt. 1 .or. ipv .gt. 6)
     .  call rx('grdepl: something wrong with indrlx')
C     Make new shear
      call dvset(vec1,1,6,alpha)
      if (ipm .eq. 2) call dvset(vec1,1,6,-alpha)
      call dpzero(dist,6)
      call grdep2(iv,iv,indrlx,vec1,dist)
      dist(7) = 1

C     goto 3

      end


