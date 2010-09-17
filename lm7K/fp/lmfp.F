      subroutine lmfp(prgnam,sbz,sctrl,sham,spot,slat,smix,
     .sspec,ssite,sstr,smove,sarray,sstrn)

      use m_struc_def  !Cgetarg
      use m_struc_func, only: mpibc1_s_lat
      use m_globalvariables
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
      type(s_array)::sarray
      type(s_bz)::sbz
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_mix)::smix
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
      type(s_ham)::sham
      type(s_pot)::spot
      type(s_str)::sstr
      type(s_move)::smove

C ... Heap
       integer ::iwdummy 

C ... Local variables
      integer procid,master,mpipid,nproc
      logical mlog
      logical lgors,cmdopt,bittst
      integer fopna,fopng,i,ifi,igets,iors,ipr,irs(5),iscr,isw,ix(5),j,
     .k,lcgf,leks,lfrce,lgunit,lpnu,lrel,lrout,nbas,nat,
     .nbaspp,ndham,nevmx,nglob,nit1,nl,nsp,nspec,numq, !npadl,npadr,lsx,
     .stdo,pdim,lsc
          integer:: o 
          real(8),pointer :: rv_p_opos(:)

         type(s_rv1),pointer :: sv_p_oorhat(:)

        complex(8),pointer :: zv_p_osmrho(:)

       complex(8),pointer :: zv_p_osmpot(:)

c      real(8) ,allocatable :: evl_rv(:)
      real(8) ,allocatable :: frc_rv(:)
      real(8) ,allocatable :: ftot_rv(:)
      real(8) ,allocatable :: pos2_rv(:)
      real(8) ,allocatable :: wk_rv(:)
      real(8) ,allocatable :: p_rv(:)

      double precision plat(3,3),qlat(3,3),qbg,xv(10),fptol,umix
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
      integer:: icom , natrlx , nvrelx , ltb , itrlx , nm , irlxsh 
     ., nitrlx , bitor , ng
      integer ,allocatable :: indrx_iv(:)
      real(8) ,allocatable :: w_rv(:)

cki, for rlse and allocate oindrx
      integer,allocatable:: iv_tmp(:)

         real(8),pointer :: rv_p_og(:)

        real(8),pointer :: rv_p_oag(:)

       integer,pointer :: iv_p_oistab(:)

      real(8) ,allocatable :: shr_rv(:)

      double precision mdprm(6),gam(4),gam1,bstim,rhosig,pletot(6,2),
     .plat0(9),dist0(9),dist(9)
      parameter (nm=3)
C ... for LDA+U
      integer nlibu,lmaxu
c      integer:: odmatu , ovorbdmat , odmato
      complex(8),allocatable::vorbdmat(:),dmatu(:),dmato(:)
      integer ,allocatable :: lldau_iv(:)

      double precision tolu
      data irlxsh /0/ dist0 /9*0d0/

c takao test
      real(8):: pnu(20),bz_w
      logical :: llmfgw
#if MPI
      integer:: ierr
      include "mpif.h"
#endif

C     parameter (T=.true., F=.false.)

      call tcn('lmfp')

cccccccccccccc
c       call upack('spec p ',sspec,1,pnu,0,0,0)
c       print *,'lmfp ini pnu=',i,pnu(1:20)
cccccccccccccc

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

      nbas=sctrl%nbas

c   mpipic() is defined in slatsm
      if ( nbas < nproc ) then
        if (procid == master ) then
          write(6,*)
          write(6,"(a,/,a,/,a,/,a)")
     &    '# For lmf-MPI, process>NBAS is not allowed.',
     &    '# If you remove this sanity check, it will fail',
     &    '# at hsibl.F and smhsbl.F.',
     &    '# You may need a better MPI version for efficent use of CPUs.'
        endif
c        call mpi_abort(MPI_COMM_WORLD,-1,ierr)
        call mpi_finalize(ierr)
        call exit(-1)
c             ---------- '-1'  is the return code of the program
        stop
c             ---------- may be stop is unnecessary.
      endif
#endif

ctakao: optbas mode is removed
C -------------------------- Basis optimization -------------------
C$$$      call defdr(o,1)
C$$$    2 continue
C$$$      call rlse(o)
C$$$      call defdr(o,1)
C$$$      if (cmdopt('--optbas',8,0,strn)) then
C$$$C       call wkdbg2
C$$$C       No self-consistency
C$$$Cgetarg         call pack1('bz nevmx',sbz,-1)
C$$$         sbz%nevmx=-1
C$$$C       No forces or dynamics
C$$$        call dpzero(mdprm,6)
C$$$Cgetarg...info...           real(8) :: mdprm(6) , real(8) :: sctrl%mdprm 6
C$$$Cgetarg         call pack2('ctrl lfrce mdprm',sctrl,0,mdprm)
C$$$         sctrl%lfrce=0
C$$$         i_copy_size=size(sctrl%mdprm)
C$$$         call dcopy(i_copy_size,mdprm,1,sctrl%mdprm,1)
C$$$ctakao
C$$$c        call lmfopb(strn(9:),slat,sspec,etot(1)-dgets('ham seref',sham))
C$$$        call lmfopb(strn(9:),slat,sspec,etot(1)-sham%seref)
C$$$C       call wkdbg2

C -------------------------- Total energy mode -------------------
c      elseif (cmdopt('--etot',6,0,strn)) then
      if (cmdopt('--etot',6,0,strn)) then
C       No forces or dynamics
        call dpzero(mdprm,6)

        sctrl%lfrce=0
        i_copy_size=size(sctrl%mdprm)
        call dcopy(i_copy_size,mdprm,1,sctrl%mdprm,1)

C       Suppress writing output density

        call lsets_bitop_i8(sctrl%lrs,1,.false.,8)


        call lsets_bitop_i8(sctrl%lrs,1,.false.,16)

C       Exactly one iteration

        sctrl%maxit=1

C       Suppress mixing of output density
c        call lstra('strn mix',i,j,k)
c        call ustrn(i,-j,1,k,i,4)
c        sstrn(k:k+3) = 'none'
c       call packs('strn mix',sstrn,'none')
c       call upacks('strn mix',i,j)
c       print *, sstrn(i:j)
c       stop 'xxxxxxxxxxxxxxxxxx'
      endif

C -------------------------- Unpack & initialization -------------------
      call getpr(ipr)

      nbas=sctrl%nbas
      nl=sctrl%nl
      nspec=sctrl%nspec
      nsp=sctrl%nspin


c       npadl=sarray%npadl
c       npadr=sarray%npadr


c       i_copy_size=size(sctrl%lpgf)
c       call i8icopy(i_copy_size,sctrl%lpgf,1,lpgf,1)
c       lcgf=sctrl%lcgf
      iscr=sctrl%lscr
c       lsx=sctrl%lsx

C     call upack('ctrl lncol',sctrl,lncol,0,0,0,0)
Changenglob      nat = nglob('nat')
      nat = globalvariables%nat

      qbg = sctrl%zbak(1)

      maxit = int(sctrl%maxit)

c       lrel = isw ( iand(-1,int(sctrl%lrel)) .ne.0 )
c nim-hrkn <nim-hrkn@foo>  2010-05-07 15:45:00
      lrel = isw ( int(sctrl%lrel) .ne.0 )

      lhf = iand(2,int(sctrl%lcd)) .ne.0

      if (lhf) maxit = 1
C     nbasp  = nbas +    npadl + npadr
      nbaspp = nbas !+ 2*(npadl + npadr)
      stdo   = lgunit(1)
C     stdl   = lgunit(2)
      call setcc(lrel)

      irs ( 1 ) = iand(7,int(sctrl%lrs)) + 8 * isw ( iand(256,int(sctrl%lrs)) .ne.0 
     .)

      irs ( 2 ) = iand(8+16,int(sctrl%lrs)) /8

      irs ( 3 ) = isw ( iand(32,int(sctrl%lrs)) .ne.0 )

      irs ( 4 ) = isw ( iand(64,int(sctrl%lrs)) .ne.0 )

      irs ( 5 ) = isw ( iand(128,int(sctrl%lrs)) .ne.0 )

C ... --rs=3 => always read from atom file
      if ( iand(3,int(sctrl%lrs)) .eq. 3 ) irs ( 1 ) = 0

C ... This shouldn't be ... for backwards compatibility
C     if (lhf) irs(1) = 0

C     Sanity checks: most ASA "extras" are not implemented here
c      call isanrg(lcgf, 0,0,'lmfp:','lcgf', .true.)
c      call isanrg(lpgf, 0,0,'lmfp:','lpgf', .true.)
      call isanrg(iscr, 0,0,'lmfp:','lscr', .true.)
c      call isanrg(lsx,  0,0,'lmfp:','lsx',  .true.)
C     call isanrg(lncol,0,0,'lmfp:','lncol',.true.)

C ... Printout properties of species
      if (ipr .ge. 30) then
C       call pr_basis (sspec,0)
        call praugm(sspec,0)
      endif

C ... Setup for no screening transformation
      if ( iand(1,int(sctrl%lbas)) .ne. 0 ) then


C       Shorten site positions

        i_copy_size=size(slat%plat)
        call dcopy(i_copy_size,slat%plat,1,plat,1)
         rv_p_opos => slat%rv_p_opos 


        if (.not. cmdopt('--shorten=no',12,0,strn)) then
C       unpack from site structure
          allocate(pos2_rv(3*nbas))

          i_copy_size=size(ssite(1)%pos)
          do i_spackv=1,nbas
            call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .      , i_copy_size , i_spackv + 1 - 1 , pos2_rv )

          enddo

          ix(1) = 2
          ix(2) = 2
          ix(3) = 2
          call info0(50,1,0,' lmfp : shortening basis vectors ... ')
           call shorps ( nbas , plat , ix , pos2_rv , rv_p_opos ) 


          i_copy_size=size(ssite(1)%pos)
          do i_spackv=1,nbaspp
             call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .       , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

          enddo

C       Debugging printout
C       call prmx('starting basis vectors',w(opos2),3,3,nbas)
C       call prmx('shortened basis vectors',w(opos),3,3,nbas)
          if (allocated(pos2_rv)) deallocate(pos2_rv)

        endif
      endif

C ... Setup for charge mixing
      call dpzero(dmxp,33)

      dmxp(2)=smix%b
      dmxp(9)=smix%bv
c      i_copy_size=size(smix%w)
c      call dcopy(i_copy_size,smix%w,1,dmxp(4),1)
      dmxp(4:5)= smix%w(1:2)
      dmxp(3)= smix%wc

      dmxp ( 6 ) = int(smix%nsave)

      dmxp ( 7 ) = int(smix%mmix)

      call parms0(0,0,0d0,0)

C ... Allocate memory for forces
      lfrce = int(sctrl%lfrce)

      if (lfrce .ne. 0) then
        numq = 1
        if ( int(sbz%lmet) .eq. 4 ) numq = 3

        allocate(frc_rv(3*nbas*numq))

      endif

C ... Relaxation setup
      itrlx = 1
C     nstack = 0
C     Initial shear was already folded into plat

      i_copy_size=size(slat%gam)
      call dcopy(i_copy_size,slat%gam,1,gam,1)
       rv_p_opos => slat%rv_p_opos 


      gam1 = gam(4)
      gam(4) = 1

      i_copy_size=size(slat%gam)
      call dcopy(i_copy_size,gam,1,slat%gam,1)
      slat%ldist=0



      nitrlx=sctrl%nitmv
      i_copy_size=size(sctrl%mdprm)
      call dcopy(i_copy_size,sctrl%mdprm,1,mdprm,1)
c       ltb=sctrl%ltb

      lshr = nint(mdprm(1)) .gt. 100
      if (nint(mdprm(1)) .eq. 0) nitrlx = 0
      if (nint(mdprm(1)) .gt. 0 .and. nint(mdprm(1)) .lt. 4) then
        call rx('lmf not set up for MD yet')
      endif
      if (nitrlx .gt. 0) then

        i_copy_size=size(sctrl%mdprm)
        call dcopy(i_copy_size,mdprm,1,sctrl%mdprm,1)


c         sctrl%ltb=bitor(ltb,16)

        allocate(indrx_iv(6*nbas))

C       Next lines in case lattice relaxation
        if (lshr) then
          if (abs(gam(4)-1) .gt. 1d-10) call rx('lmfp: '//
     .    'use of SHEAR= incompatible w/ lattice relaxation')

          i_copy_size=size(slat%plat0)
          call dcopy(i_copy_size,slat%plat0,1,plat0,1)

        endif
        call rlxstp ( sctrl , ssite , natrlx , nvrelx , indrx_iv , 
     .  xyzfrz , pdim )

cki         if (allocated(indrx_iv)) deallocate(indrx_iv)

        icom = 0
        if (nvrelx .ne. 0) then
Cki The new indrx_iv uses the content of the old indrx_iv.
          isize_iv_tmp=size(indrx_iv)
          allocate( iv_tmp(isize_iv_tmp) )
          iv_tmp(:)= indrx_iv(:)
          deallocate(indrx_iv)
          allocate(indrx_iv(2*natrlx))
          isize_iv_tmp=min(isize_iv_tmp,2*natrlx)
          indrx_iv(:isize_iv_tmp)=iv_tmp(:isize_iv_tmp)
          deallocate(iv_tmp)

          allocate(w_rv(nvrelx*nvrelx))

          allocate(p_rv(pdim))

        endif
        alat = slat%alat

        if (procid .eq. master) then
          ifi = fopna('bsmv',-1,0)
          allocate(pos2_rv(3*nbas))

          j = 1
          call ivset(ix,1,3,j)
           call shorps ( nbas , plat , ix , rv_p_opos , pos2_rv ) 


          call iobsm0 ( 0 , bstim , 0d0 , 0d0 , nbas , alat , pos2_rv 
     .    , ifi )

          if (allocated(pos2_rv)) deallocate(pos2_rv)

        endif
      endif

C ... Re-entry for shear distortion
      allocate(shr_rv(1))

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

      call suham(sctrl,slat,sspec,ssite,sarray,sham,spot,sstrn,sbz)



C ... Quit if --quit=ham given
      if ( int(sctrl%quit) .eq. 8 ) then

        call rx0('quit = ham')
      endif

C --- Setup for iterations in a self-consistency cycle ---
C ... Unpack or allocate some permanent arrays
      sv_p_oorhat => spot%sv_p_oorhat 
      zv_p_osmrho => spot%zv_p_osmrho 
      zv_p_osmpot => spot%zv_p_osmpot 

      allocate(ftot_rv(3*nbas))
      ndham=sham%ndham
c      allocate(evl_rv(ndham*2))

C ... Set various switches
C     Whether forces, and how to calculate non Helman-Feynman corr.
      lfrce = int(sctrl%lfrce)

C     Maximum number of eigenvalues
      nevmx = int(sbz%nevmx)

C     Whether to evaluate output density and/or KS energy
      lrout = 1
      leks = 1
      j = 6
      if (cmdopt('-leks=',j,0,strn)) then
        if (.not. a2bin(strn,leks,2,0,' ',j,72))
     .  call rxs('failed to parse',strn)
      endif
      if (nevmx .eq. -1) then
        lrout = 0
        leks = 0

        sbz%nevmx=nevmx

      endif
C     Whether to float pnu's
      lpnu = 1
C     Sanity checks
      if (lrout .eq. 0 .and. lfrce .ne. 0) then
        write(stdo,333) 'when forces sought'
  333   format('lmfp (fatal): output density required ',a/
     .  '      To make output density turn off HF=t and/or NEVMX<0')
        call rx('incompatible input')
      endif
C     Sanity checks
      if (lrout .eq. 0 .and. cmdopt('--etot',6,0,strn)) then
        write(stdo,333) 'with --etot switch.'
        call rx('incompatible input')
      endif
      if (lrout.eq.0 .and. maxit.gt.1 .and. ipr.ge.20) then
        call awrit1('%N lmfp (warning): %i iterations sought but no'//
     .  ' output rho ... do 1 iteration',' ',80,stdo,maxit)
        maxit = 1
      endif

C... LDA+U initialization
      allocate(lldau_iv(nbas))
      lldau_iv(:)=0
C     Check for LDA+U ... return nlibu > 0 if any U blocks
      call suldau ( nbas , sspec , ssite , nlibu , lmaxu , lldau_iv )
      if (nlibu .gt. 0) then
        i = nsp*nlibu*(lmaxu*2+1)**2
      else
        i = 1 
      endif
      allocate(vorbdmat(i),dmatu(i),dmato(i))
      vorbdmat = 0d0
      dmatu = 0d0
      dmato = 0d0
      if (nlibu > 0) then
C       need group info to symmetrize site density matrix
        ng=slat%nsgrp
        iv_p_oistab => slat%iv_p_oistab 
        rv_p_og => slat%rv_p_osymgr 
        umix=smix%umix
        tolu=smix%tolu
        if (umix .eq. 0) umix = 1
C       if (tolu .eq. 0d0) tolu = 1d-4
C       initialize density matrix for LDA+U
          call sudmtu ( nbas , nsp , nlibu , lmaxu , ssite , sspec , 0 
     .    , lldau_iv , ng , rv_p_og , iv_p_oistab , dmatu , vorbdmat ) 
      endif
C     end LDA+U  initialization section


C ... Invoke the restart editor
      if (cmdopt('--rsedit',8,0,strn)) then
        call rsedit(strn(9:),1,sctrl,ssite,sspec,slat,spot,sbz,
     .  nbas,nat,nspec)
        call rx0('lmfp from rsedit')
      endif

C ... Invoke the response function editor
      if (cmdopt('--chimedit',10,0,strn)) then
        call chimedit(strn(11:),1,sctrl,ssite,sspec,slat,spot,sbz,
     .  nbas,nat,nspec)
        call rx0('lmfp from chimedit')
      endif

cccccccccccccccccc
c      do  i = 1, nbas
c        call upack2('site pnu',ssite,i,pnu)
c        print *,'check before reentry pnu=',i,pnu(1:20)
c      enddo

cccccccccccccccccc


C ---------------- Re-entry point for a new iteration ---------------
      iter = 1
    5 continue

cccccccccccccccccccccccccccccccccccccccccccccc
ctakao
c      do  i = 1, nbas
c        call upack2('site pnu',ssite,i,pnu)
c        print *,'check before 10 write pnu=',i,pnu(1:20)
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccc
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
         call rdovfa ( nbas , nspec , ssite , sspec , slat , spot , qbg 
     .   , sv_p_oorhat ) 

        nit1 = 0
      elseif (mod(irs(1),4) .ge. 1 .and. mod(irs(1),4) .le. 2) then
        lbin = .not. bittst(irs(1),2)
        k = -1
        if (procid .eq. master) then
          if (lbin) ifi = fopna('rst',-1,4)
          if (.not. lbin) ifi = fopna('rsta',-1,0)
        endif
        call mpibc1(ifi,1,2,mlog,'lmfp','ifi')
ctakao BZ_W is

        bz_w=sbz%w

c        print *,'yyyyyyy goto iors osmrho=',spot%osmrho
          k = iors ( 1 , sctrl , ssite , sspec , slat , spot , sbz , fileid 
     .    , nbas , nat , nspec , sv_p_oorhat , iwdummy , nit1 , lbin , 
     .    ifi ) 


c        print *,'yyyyyyyyyy end of iors'
ctakao

        sbz%w=bz_w

        if (k .lt. 0 .and. nproc .gt. 0) then
C         irs(1) = irs(1) - mod(irs(1),4)
          call rx('MPI: rst read failed. Restart with --rs=0')
        endif
        if (procid .eq. master) then
          call fclose(ifi)
        endif
        call mpibc1(k,1,2,mlog,'lmfp','k')

C   ... Write positions array from site structure

         rv_p_opos => slat%rv_p_opos 


        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbaspp
           call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .     , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

        enddo

cki        call mpibc1(slat,nint(slat%size),4,.false.,'lmfp','slat')
        call mpibc1_s_lat(slat,.false.,'lmfp','slat')
        if (k .lt. 0) then
          irs(1) = 0
          goto 10
        endif
        if (mod(irs(1),8) .ge. 4) then
C         If no force switch set, use default
          k = int(sctrl%lfrce)

          if ( k .eq. 0 ) then

            sctrl%lfrce=1
          endif

            call smshft ( ssite , sspec , slat , sctrl , sham , sv_p_oorhat 
     .      , zv_p_osmrho ) 


C         Restore force switch
          if ( k .eq. 0 ) then

            sctrl%lfrce=k
          endif

        endif
        if (mod(irs(1),16) .ge. 8) then
C          call upack('lat dist',slat,xv,0,0,0,0)
C          call pvsms2(ssite,sspec,xv,nbas,nsp,w(oorhat))
          irs(1) = irs(1)-8
          i = int(sctrl%lrs)

          i = i-256

          sctrl%lrs=i

        endif
      endif

cccccccccccccccccccccccccccccccccccccccccccccc
ctakao
c      do  i = 1, nbas
c        call upack2('site pnu',ssite,i,pnu)
c        print *,'check write pnu=',i,pnu(1:20)
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccc

C ... Write positions after file read, and repack
      if (ipr .ge. 50) then
        write(stdo,357) 'Basis, after reading restart file'
  357   format(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
        call dinv33(plat,1,qlat,xv)
        do  i = 1, nbas

          j=ssite(i)%spec
          i_copy_size=size(ssite(i)%pos)
          call dcopy(i_copy_size,ssite(i)%pos,1,xv,1)

          do i_spacks=j,j
            call spacks_copy('u',sspec(i_spacks)%name,j,j,alabl,i_spacks)
          enddo

          call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
          write(stdo,345) i, alabl, (xv(j),j=1,3), (xv(3+j),j=1,3)
  345     format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
        enddo
      endif

C --- Optionally re-shorten basis vectors ---
      if (cmdopt('--shorps',8,0,strn)) then
        allocate(pos2_rv(3*nbas))

        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbas
          call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .    , i_copy_size , i_spackv + 1 - 1 , pos2_rv )

        enddo

        ix(1) = 2
        ix(2) = 2
        ix(3) = 2
         call shorps ( - nbas , plat , ix , pos2_rv , rv_p_opos ) 


        call info0(20,1,-1,
     .  ' lmfp  : write shortening vectors to file shorps ...')
         call iopos ( .true. , 1 , 'shorps' , nbas , rv_p_opos ) 

         call shorps ( nbas , plat , ix , pos2_rv , rv_p_opos ) 


        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbaspp
           call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .     , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

        enddo

C       Debugging printout
C       call prmx('starting basis vectors',w(opos2),3,3,nbas)
C       call prmx('shortened basis vectors',w(opos),3,3,nbas)
        if (allocated(pos2_rv)) deallocate(pos2_rv)

      endif

C     Hang on to previous site density matrix for this iteration
      if (nlibu .gt. 0) then
c        i = nsp*nlibu*(lmaxu*2+1)**2
c        call zcopy(i,w(odmatu),1,w(odmato),1)
c        call dpzero(w(odmatu),2*i)
        dmato=dmatu
        dmatu=0d0
      endif

C --- Make and diagonalize hamiltonian, make new charge density ---
      if (maxit .eq. 0) call info0(20,1,0,
     .' lmfp  : zero iterations sought ... no band pass')

      llmfgw=.false. !takao now just lmfgw mode by lmfgw switch.
      if(trim(prgnam)=='LMFGWD') llmfgw=.true.

      call bndfp ( nbas , nsp , nlibu , lmaxu , lldau_iv , ssite 
     ., sspec , slat , sctrl , sham , spot , sbz , sstrn , ndham , 
     . leks , lrout , lfrce , lpnu , dmxp , iter , maxit  , !evl_rv, 
     . frc_rv , dmatu, vorbdmat, llmfgw )


C ... check convergence of dmatu and update it and vorbdmat if necessary
      if (nlibu .gt. 0 .and. maxit .gt. 0 .and. lrout .gt. 0) then
          call chkdmu ( nbas , nsp , nlibu , lmaxu , ssite , sspec , sham 
     .    , 0 , dmatu , dmato, vorbdmat, tolu , umix 
     .    , lldau_iv , ng , rv_p_og , iv_p_oistab ) 



      endif

ctakao optbas mode is removed
C ... basis optimization: just extract etot(1) and return to opt.
C$$$      if (cmdopt('--optbas',8,0,strn)) then
C$$$Cgetarg...info...           real(8) :: etot(2) , real(8) :: sham%ehf 1
C$$$Cgetarg...info...           real(8) :: etot(2) , real(8) :: sham%ehk 1
C$$$Cgetarg         call upack2('ham ehf ehk',sham,etot(1),etot(2))
C$$$         etot(1)=sham%ehf
C$$$         etot(2)=sham%ehk
C$$$        goto 2
C$$$      endif

C --- Write smoothed charge density for contour plotting ---
      if (procid .eq. master) then
        if (cmdopt('--wden',6,0,strn)) then
            call ioden ( strn ( 7: ) , slat , ssite , sspec , sv_p_oorhat 
     .      , zv_p_osmrho ) 


        endif

C --- Write restart file (skip if --quit=band) ---
        if ( int(sctrl%quit) .ne. 4 ) then

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
              k = iors ( 1 , sctrl , ssite , sspec , slat , spot , sbz , fileid 
     .        , nbas , nat , nspec , sv_p_oorhat , iwdummy , iter , lbin , 
     .        - ifi ) 


            call fclose(ifi)
c
c            takao save rst file option
c            if(cmdopt('--saveallrst',12,0,strn) then
c               call fsystm('cp "//rst)
            
          endif
        endif
      endif
      if (cmdopt('--window=',9,0,strn))
     .call rx0('lmf : early exit (--window option)')

C --- Add to save file; decide on next iteration ---
      if (maxit .gt. 0) then

        etot(1)=sham%ehf
        etot(2)=sham%ehk

C ... Subtract reference energy
        seref = sham%seref

        etot(1) = etot(1) - seref
        if (etot(2) .ne. 0) etot(2) = etot(2) - seref
        amom = sham%eterms(15)

C     The desired tolerances in q,e
        qtol = sctrl%tol(1)

        etol = sctrl%tol(3)

        if (procid .eq. master) then
          rhosig = sham%eterms(19)

          i = 0
          if (rhosig .ne. -99 .and. rhosig .ne. 0) i = 10
          call nwit ( int(sctrl%nvario) , iter , maxit , lhf.or.irs ( 1 
     .    ) .eq.0.and.iter.eq.1 , leks+i , etol , qtol , qdiff , 'cxhi' 
     .    , amom , etot , lsc )

        endif
        call mpibc1(lsc,1,2,mlog,'lmfp','lsc')
        if (lsc .eq. 2 .and. .not. lhf .and. maxit .gt. 1) lsc = 3
        if (lsc .eq. 1 .and. lrout .gt. 0  .or. lsc .eq. 3) then
          call query('max it=',2,maxit)
          if (iter .ge. maxit) lsc = 1
          if (iter .lt. maxit) lsc = 3
        endif
        if ( int(sctrl%quit) .eq. 4 ) call rx0 ( 'lmf : exit (--quit=band)' 
     .  )

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
           allocate(rv_p_opos(abs(3*nbas)))
           if (3*nbas<0) rv_p_opos(:)=0.0d0

          i_copy_size=size(ssite(1)%pos)
          do i_spackv=1,nbas
             call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .       , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

          enddo


          i_copy_size=size(sctrl%mdprm)
          call dcopy(i_copy_size,sctrl%mdprm,1,mdprm,1)

          if (lshr) then
            call grdepl ( nvrelx , indrx_iv , 0.01d0 , etot , irlxsh , 
     .      pletot , dist )

            if (irlxsh .ne. 0) then
              call grdep2 ( 1 , nvrelx , indrx_iv , dist0 , dist )

              goto 98
            else
              call relax ( prgnam , sctrl , ssite , sspec , itrlx , indrx_iv 
     .        , natrlx , nvrelx , pletot , p_rv , w_rv , 0 , 0d0 , dist0 
     .        , icom )

              call dpzero(dist,6)
              call grdep2 ( 1 , nvrelx , indrx_iv , dist0 , dist )

              dist(7) = 1
            endif

          else
             call relax ( prgnam , sctrl , ssite , sspec , itrlx , indrx_iv 
     .       , natrlx , nvrelx , frc_rv , p_rv , w_rv , 0 , 0d0 , rv_p_opos 
     .       , icom ) 


          endif

C       Restore lattice symmetry to machine precision
          if (cmdopt('--fixpos',8,0,strn)) then

            ng=slat%nsgrp
             iv_p_oistab => slat%iv_p_oistab 

             rv_p_oag => slat%rv_p_oag 

             rv_p_og => slat%rv_p_osymgr 


C         call shoist(w(oistab),nbas,w(oag),w(og),ng)
            j = 8+1
            if (strn(9:13) .eq. ':tol=') then
              j = 13
            endif
            if (strn(9:9) .ne. ':' .or.
     .      .not. a2bin(strn,fptol,4,0,' ',j,len(strn))) fptol = 1d-5
                call fixpos ( rv_p_opos , nbas , fptol , ng , plat , rv_p_og 
     .          , rv_p_oag , iv_p_oistab ) 




          endif

C       Write updated positions to bsmv file
          if (procid .eq. master .and. .not. lshr) then
            ifi = fopna('bsmv',-1,0)
            call poseof(ifi)
            bstim = bstim+1
            allocate(pos2_rv(3*nbas))

            j = 1
            call ivset(ix,1,3,j)
             call shorps ( nbas , plat , ix , rv_p_opos , pos2_rv ) 


            call iobsmv ( 0 , bstim , 0d0 , 0d0 , nbas , alat , pos2_rv 
     .      , - ifi )

            if (allocated(pos2_rv)) deallocate(pos2_rv)

            call fclose(ifi)
          endif
C       repack updated positions in site structure
          i_copy_size=size(ssite(1)%pos)
          do i_spackv=1,nbas
             call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .       , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

          enddo


C        if (icom .eq. 1) then
C          print *, '!!'
C        endif

C   ... Exit when relaxation converged or maximum number of iterations
          if (icom .eq. 1) then
            if (procid .eq. master) then
              k = int(sctrl%nvario)

              flg = 'C67'
              call nwitsv(1+2,k,flg,nsp,amom,etot)
            endif
            call tcx('lmfp')
            call fexit(0,111,
     .      ' LMFP: relaxation converged after %i iteration(s)',itrlx)
          else
            call query('proceed with next relaxation step',-1,0)
          endif

C   ... Restore minimum gradient positions if this is last step
C        if (itrlx .eq. nitrlx .and. icom .eq. -1) then
          if (itrlx .eq. nitrlx) then

            if (.not. lshr) then
              call info0(20,1,0,' lmfp: restore positions for minimum g')

C           call prmx('initial positions',w(opos),3,3,nbas)
               call prelx1 ( 1 , nm , lshr , natrlx , nvrelx , indrx_iv , p_rv 
     .         , rv_p_opos ) 


C           call prmx('minimum-g positions',w(opos),3,3,nbas)
            else
              call prelx1 ( 1 , nm , lshr , natrlx , nvrelx , indrx_iv , 
     .        p_rv , dist0 )

              call dpzero(dist,6)
              call grdep2 ( 1 , nvrelx , indrx_iv , dist0 , dist )

              call info2(20,0,0,
     .        ' lmfp : strain of minimum gradient:'//
     .        '%N   PDEF=%6;8,4D'//
     .        '%N STRAIN=%6;8,4D',
     .        dist0,dist)
            endif

C         Repack updated positions in site structure
            i_copy_size=size(ssite(1)%pos)
            do i_spackv=1,nbas
               call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .         , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

            enddo


          endif

C   ... New density after atom shifts
C       If explicitly told to read from atom files after atom movmment
          if ( iand(3,int(sctrl%lrs)) .eq. 3 ) then

            irs(1) = 0

C       Else, use self-consistent
          else if (.not. lshr) then
            irs(1) = 3
              call smshft ( ssite , sspec , slat , sctrl , sham , sv_p_oorhat 
     .        , zv_p_osmrho ) 


          endif

C   ... Write restart file (to include new positions)
          if (procid .eq. master .and. .not. lshr) then
            ifi = fopna('rst',-1,4)
            call upacks('strn jobid',i,j)
            fileid = 'lmfp:  ' // sstrn(i:j)
              k = iors ( 1 , sctrl , ssite , sspec , slat , spot , sbz , fileid 
     .        , nbas , nat , nspec , sv_p_oorhat , iwdummy , iter , .true. 
     .        , - ifi ) 


            call fclose(ifi)
          endif

C   ... Write positions to file
          if (cmdopt('--wpos=',7,0,strn) .or.
     .    cmdopt('--wpos:mode1:',13,0,strn)) then

            i_copy_size=size(slat%gam)
            call dcopy(i_copy_size,slat%gam,1,gam,1)

            allocate(wk_rv(3*nbas))

             call rdistn ( rv_p_opos , wk_rv , nbas , gam ( 1 ) , gam ( 2 
     .       ) , gam ( 3 ) , 1 / gam1 ) 


            call iopos ( .true. , 0 , strn ( 8: ) , nbas , wk_rv )

C         call fclr(strn(8:),-1)
            if (allocated(wk_rv)) deallocate(wk_rv)

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
     .        ' LMFP: relaxation incomplete after %i iteration(s)',nitrlx)
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
     .cmdopt('--wpos:mode1:',13,0,strn)) then
         allocate(rv_p_opos(abs(3*nbas)))
         if (3*nbas<0) rv_p_opos(:)=0.0d0

        i_copy_size=size(ssite(1)%pos)
        do i_spackv=1,nbas
           call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .     , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

        enddo


        i_copy_size=size(slat%gam)
        call dcopy(i_copy_size,slat%gam,1,gam,1)

        allocate(wk_rv(3*nbas))

         call rdistn ( rv_p_opos , wk_rv , nbas , gam ( 1 ) , gam ( 2 
     .   ) , gam ( 3 ) , 1 / gam1 ) 


        if (cmdopt('--wpos:mode1:',13,0,strn)) then
          call iopos ( .true. , 1 , strn ( 14: ) , nbas , wk_rv )

C         call fclr(strn(14:),-1)
        else
          call iopos ( .true. , 0 , strn ( 8: ) , nbas , wk_rv )

C         call fclr(strn(8:),-1)
        endif
        if (allocated(wk_rv)) deallocate(wk_rv)
         if (associated(rv_p_opos)) deallocate(rv_p_opos)


      endif

      call tcx('lmfp')
C     if (maxit .eq. 0) call rx0(' zero iterations sought ... quitting')
Cki#error, have return with len(w_varlist)>0 at line 1250
      if (allocated(lldau_iv)) deallocate(lldau_iv)
c      if (allocated(evl_rv)) deallocate(evl_rv)
      if (allocated(ftot_rv)) deallocate(ftot_rv)
      if (allocated(shr_rv)) deallocate(shr_rv)
      if (allocated(p_rv)) deallocate(p_rv)
      if (allocated(w_rv)) deallocate(w_rv)
      if (allocated(indrx_iv)) deallocate(indrx_iv)
      if (allocated(frc_rv)) deallocate(frc_rv)
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
      allocate(pos2_rv(3*nbas))

      call dinv33(plat,0,xv,xv(10))
       call dgemm ( 'N' , 'N' , 3 , nbas , 3 , 1d0 , xv , 3 , rv_p_opos 
     . , 3 , 0d0 , pos2_rv , 3 ) 


C     Simultaneously pack in lat->pos and site->pos

       rv_p_opos => slat%rv_p_opos 


       call dgemm ( 'N' , 'N' , 3 , nbas , 3 , 1d0 , plat0 , 3 , pos2_rv 
     . , 3 , 0d0 , rv_p_opos , 3 ) 


C     call prmx('positions for plat0',w(opos),3,3,nbas)
      i_copy_size=size(ssite(1)%pos)
      do i_spackv=1,nbas
         call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .   , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

      enddo

      call cppos(1,nbas,ssite)
C     New shear

      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,plat0,1,slat%plat,1)
      slat%ldist=3
      i_copy_size=size(slat%dist)
      call dcopy(i_copy_size,dist,1,slat%dist,1)

C     A little memory leakage rel to 1st pass, but not so serious
      if (allocated(pos2_rv)) deallocate(pos2_rv)
      if (allocated(lldau_iv)) deallocate(lldau_iv)
c      if (allocated(evl_rv)) deallocate(evl_rv)
      if (allocated(ftot_rv)) deallocate(ftot_rv)
      if (allocated(shr_rv)) deallocate(shr_rv)
c      call rlse(ovorbdmat)
      deallocate(vorbdmat,dmato,dmatu)

      call lattic(slat,sctrl,ssite,sarray)
      allocate(shr_rv(1))


      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)

C     Remake qp
      ltmp = iand(1,int(sctrl%lmet)) .ne.0 .or. iand(4+2+1,int(sctrl%ldos)) .ne.0

      call mkqp(sctrl,sbz,slat,ltmp,.false.,1,-2)

C ... Write restart file (to include new positions)
      if (procid .eq. master .and. irlxsh .eq. 0) then
        ifi = fopna('rst',-1,4)
        call upacks('strn jobid',i,j)
        fileid = 'lmfp:  ' // sstrn(i:j)
          k = iors ( 1 , sctrl , ssite , sspec , slat , spot , sbz , fileid 
     .    , nbas , nat , nspec , sv_p_oorhat , iwdummy , iter , .true. 
     .    , - ifi ) 


        call fclose(ifi)
      endif

C     Decide on what density to use
      if ( iand(3,int(sctrl%lrs)) .eq. 3 ) then

        irs(1) = 0
C     Else, use file density
      else
        irs ( 1 ) = iand(7,int(sctrl%lrs))

      endif
      goto 4

      if (allocated(shr_rv)) deallocate(shr_rv)
      if (allocated(p_rv)) deallocate(p_rv)
      if (allocated(w_rv)) deallocate(w_rv)
      if (allocated(indrx_iv)) deallocate(indrx_iv)
      if (allocated(frc_rv)) deallocate(frc_rv)

      end subroutine lmfp


      subroutine cppos(ib1,ib2,ssite)

      use m_struc_def  !Cgetarg

C- Copy site positions to p0 for a range of sites
C     implicit none
      integer ib1,ib2
      real(8):: pos(3)
      type(s_site)::ssite(*)

      integer ib

      do  ib = ib1, ib2
        i_copy_size=size(ssite(1)%pos)
        do i_spackv=ib,ib
          call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,1,pos)
        enddo

        i_copy_size=size(ssite(1)%pos0)
        do i_spackv=ib,ib
          call spackv_array_copy_r8_r8('p',ssite(i_spackv)%pos0,i_copy_size,1,pos)
        enddo

      enddo
      end subroutine cppos


      subroutine grdepl(nvrelx,indrlx,alpha,etot,irlxsh,pletot,dist)

C-
C     implicit none
      integer irlxsh,nvrelx,indrlx(nvrelx)
      double precision pletot(6,2),etot,dist(9),alpha
      double precision grad,vec1(6)
      integer iv,ipm,ipv

      call info5(30,1,0,' GRDEPL: point %i of %i for grad shear: '//
     .'etot=%d',irlxsh,2*nvrelx,etot,0,0)

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
     .call rx('grdepl: something wrong with indrlx')
C     Make new shear
      call dvset(vec1,1,6,alpha)
      if (ipm .eq. 2) call dvset(vec1,1,6,-alpha)
      call dpzero(dist,6)
      call grdep2(iv,iv,indrlx,vec1,dist)
      dist(7) = 1

C     goto 3

      end subroutine grdepl



