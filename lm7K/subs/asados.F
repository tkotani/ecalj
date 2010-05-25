      subroutine asados(prgnam,sbz,sctrl,sham,slat,sspec,ssite,sarray)

      use m_struc_def  !Cgetarg

C- Partial density-of-states generator
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:string labelling name of caller (for printout)
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: n w efmax lcond dosw nkabc nkp ntet oidtet owtkp oqp
Ci                ostar oipq
Ci     Stored:    *
Ci     Passed to: mkqp
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nclass nl nspec lncol lham ldos nsp
Ci     Stored:    *
Ci     Passed to: lgors subasi mkqp
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: subasi
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: avw alat plat
Ci     Stored:    npgrp
Ci     Passed to: mkqp
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z lmxa
Ci     Stored:    *
Ci     Passed to: subasi
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp oics oipc ormax oclabl
Ci     Stored:    *
Ci     Passed to: *
Ci   Optional command-line switches:
Ci     --dos:option1:option2...
Ci       Options are:
Ci       wtfn=name  read dos weights from file `name'
Ci       cls        is equivalent to wtfn=CLS
Ci       tbdos      uses moments file with tbe conventions; also moments
Ci                  file is named 'BAND'
Ci       mode=#     makes dos-related quantities:
Ci                  #=0 makes dos
Ci                  #=1 (ballistic) conductivity,
Ci                      1/2 | grad_k E(k) . vec
Ci                      In this mode, vec must also be specified.
Ci                  #=2 (diffusive) conductivity, or actually
Ci                      grad_1 E(k) . grad_2 E(k) where directions
Ci                      1 and 2 must be specifed by vec and vec2 below
Ci       idos       generate energy-integrated dos
Ci       fpdos      tells asados to run in 'fp' mode.  It only
Ci                  affects printout associating dos channels with
Ci                  site and class names
Ci       npts=#     number of energy points.  If not specified with
Ci                  command-line --dos, user is prompted for input.
Ci       window=#,# energy window over which data is to be computed
Ci                  #,# are minimum, maximum energies.  If not
Ci                  specified with command-line --dos, user is
Ci                  prompted for input.
Ci       vec=#,#,#  direction vector for conductivity; see mode=#
Ci       vec2=#,#,# second direction vector for `diffusive'
Ci                  conductivity; see mode= above
Ci       totdos     compute total dos by adding weights from all
Ci                  partial dos
Ci       bands=list compute contribution to dos from a prescribed
Ci                  list of bands
Ci       classes=list generate dos for a specified list of classes.
Ci                  list syntax follows standard class-list syntax
Ci                  as described in lmto.html.
Co Outputs
Co    DOS is generated and written to file 'dos'
Cl Local variables
Cu Updates
Cu   08 Jun 06  Bug fix for noncollinear case; new --mull argument
Cu    8 Jun 05  (ATP) handles empirical TB case
Cu   13 Sep 03  Revised printout of channels
Cu   18 Jan 02  Revised printout of identification of channels
Cu              with site orbitals
Cu   07 Feb 01 *several revisions to deal with the following cases:
Cu              different modes for DOS-related quantities
Cu              integration over a subset of bands
Cu              possibility to integrate quantities for total DOS
Cu              unification of nfp,etb,asa
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) prgnam*8
      type(s_array)::sarray
      type(s_bz)::sbz
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
      type(s_ham)::sham

C ... Local variables
      character*120 outs,sopts,dc*1,wtfn
      logical bittst,T,F,ltet,bzmap,lgors,lmet,cmdopt,tbdos,TBU,totdos,
     .lcond,lrange,lidos,lfp
      integer:: ogstar , oics , oidtet , oipc , oqp , ormax , owtkp 
     ., oipq , oclabl
      real(8) ,allocatable :: dos_rv(:)
      real(8) ,allocatable :: doswt_rv(:)
      real(8) ,allocatable :: eband_rv(:)
      real(8) ,allocatable :: ew_rv(:)
      integer ,allocatable :: nq_iv(:)
      integer ,allocatable :: nw_iv(:)
      real(8) ,allocatable :: vx_rv(:)
      real(8) ,allocatable :: vy_rv(:)
      real(8) ,allocatable :: wk_rv(:)
      real(8) ,allocatable :: x_rv(:)
      real(8) ,allocatable :: y_rv(:)
      real(8) ,allocatable :: z_rv(:)
      integer ,allocatable :: chan_iv(:)

      integer awrite,fopnn,fopno,i,iblst(1000),ichds,icond,il,imax,
     .iprint,igets,isw,iv(10),j,j1,j2,l,nband,nbandx,ldos, !lham,
     .lncol,lstyle,m,n1,n2,nbas,nblst,nchan,nchan2,nchds,nclasp,
     .nclass,ndum,nevmx,nfilem,nfilqp,nfstg,nkp,nkxyz(3),nl,nlf,nlo,
     .nlstc,norder,np0,npgrp,npln,npts,nqmx,nsp,nspc,nspd,nspec,nspx,
     .ntet,nwin,nwmx,parg,scrwid,stdo,moddos,nsite,nsitmx,lmxch,nll,
     .ib,iclbas,ic,ival,is,lmax,igetss,lopt,nchmx,nglob
      parameter (nsitmx = 1024)
C     This is a site or class list
      integer lsite(nsitmx)
      double precision vmtz(2),xx,efermi(2),efmax,emin,emax,drange(2),
     .width,wtbzm,cvec(0:6),del,ddot,avw,alat,plat(3,3)
      character*40 modstr(0:4),strn*120,strn2*120
      character*8 clabl
      parameter (T=.true., F=.false., scrwid=80)
C ... Heap
      integer w(1)
      common /w/ w
      data modstr /'dos','sigma(ballistic)','sigma(v*v)','*','*'/
ctakao
      integer,allocatable:: ics(:)
      real(8),allocatable:: z(:)

C --- Setup ---
      stdo = nglob('stdo')

      nbas=sctrl%nbas
      nclass=sctrl%nclass
      nl=sctrl%nl
      nspec=sctrl%nspec


      lncol=sctrl%lncol
c       lham=sctrl%lham
      ldos=sctrl%ldos


      nclasp=sarray%nclasp
      oics=sarray%oics
      oipc=sarray%oipc
      ormax=sarray%ormax


      oclabl=sarray%oclabl


      avw=slat%avw
      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)


      norder=sbz%n
      width=sbz%w
      efmax=sbz%efmax
      i_copy_size=size(sbz%lcond)
      call dcopy(i_copy_size,sbz%lcond,1,cvec,1)
      i_copy_size=size(sbz%dosw)
      call dcopy(i_copy_size,sbz%dosw,1,drange,1)

      norder = isign(1,norder) * mod(iabs(norder),100)
C     For now... until input in rdccat is patched
      call dcopy(3,cvec(1),1,cvec(4),1)
C ... Class-dependent arrays

c     call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)
      allocate(ics(nclasp),z(nclasp))
      call icopy(nclasp,w(oics),1,ics,1)
      do i=1,nclasp
        z(i)=sspec(ics(i))%z
      enddo

c      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmx)
C     call sp2cls('spec idxdn',sspec,w(oics),nl,1,nclass,oidxdn)
C ... Other initialization
      lmet = iand(1,int(sctrl%lmet)) .ne.0

      ltet = iand(2,int(sctrl%lmet)) .ne.0

      icond = nint(cvec(0))
      bzmap = bittst(ldos,8)
      call subasi(sctrl,sspec,sham)
      lfp = .false.

C --- Check for --pdos or --mull command-line argument ---
C     Tells lmdos that channels are oriented by site, not class.
C     It does not affect calculation, only the identification
C     of the channels printed at the end
C     ldos is not set properly when asados is called, so make assumptions
      nsite = 0
      lsite(1) = -1
      nll = nl
      if (ldos .eq. 0) ldos = 2
      if (cmdopt('--pdos',6,0,strn).or.cmdopt('--mull',6,0,strn2)) then
        lopt = 0
        if (cmdopt('--mull',6,0,strn)) then
          strn = strn2
          lopt = 1
        endif
        nchmx = nbas*nl**2
        allocate(chan_iv(nchmx))
Ctakao_ZeroClear_NotRequiered          if (nchmx<0) chan_iv(:)=0

        call sumlst ( lopt , nchmx , nbas , 1 , ssite , sspec , strn 
     .  ( 7: ) , moddos , nsite , lsite , lmxch , nchan , chan_iv 
     .  , nll )

        if (allocated(chan_iv)) deallocate(chan_iv)

        if (moddos .eq. 0 .or. moddos .eq. 3) ldos = 0
        if (moddos .eq. 1 .or. moddos .eq. 4) ldos = 2
        if (moddos .eq. 2 .or. moddos .eq. 5) ldos = 4
      endif
      nlo = nl
      if (ldos .eq. 0) nlo = 1
      if (ldos .eq. 2) nlo = nll
      if (ldos .eq. 4) nlo = nll*nll

C --- Parse DOS options ---
C ... Defaults
      nlstc = 0
      nchds = -1
      lstyle = 1
      wtfn = ' '
      TBU = .false.
      tbdos = .false.
      totdos = .false.
      iblst(1) = 0
      nblst = 0
      npts = 0
C     this will be set to .true. if set explicitly
      lrange = .false.
C     this will be set to .true. if integrated DOS sought.
      lidos = .false.
C     This may be set to .true. by switch
      lcond = .false.

C      print *, 'LCOND'
C      LCOND = .true.

      if (cmdopt('--dos',5,0,sopts)) then
        dc = sopts(6:6)
        if (dc .ne. ' ') then
C   ... Return here to resume parsing for arguments
          j2 = 5
   10     continue
          j2 = j2+1
          if (sopts(j2:j2) .eq. dc) goto 10
          j1 = min(len(sopts),j2)
          call nwordg(sopts,0,dc//' ',1,j1,j2)
          if (j2 .ge. j1) then
            if (.false.) then
            elseif (sopts(j1:j1+4) .eq. 'wtfn=')  then
              if (j1+5 .gt. j2) call rx('ASADOS: bad file name')
              wtfn = sopts(j1+5:j2)
C         cls switch just sets file name to 'cls'
            elseif (sopts(j1:j2) .eq. 'cls')  then
              wtfn = 'CLS'
              lfp = .true.
C         generate integrated DOS
            elseif (sopts(j1:j2) .eq. 'idos')  then
              lidos = .true.
C         tbdos has different dimensions in moments file
            elseif (sopts(j1:j2) .eq. 'tbdos')  then
              tbdos = .true.
C         tbdos has different dimensions in moments file
            elseif (sopts(j1:j2) .eq. 'tbu')  then
              tbdos = .true.
              TBU = .true.
C         fp mode has site-dependent dos channels (for printout)
            elseif (sopts(j1:j2) .eq. 'fpdos')  then
              lfp = .true.
C         Calculate conductivity or some quantity other than dos
            elseif (sopts(j1:j1+4) .eq. 'mode=')  then
              m = 0
              i = parg('mode=',2,sopts(j1:),m,len(sopts(j1:)),
     .        dc//' ',1,1,iv,icond)
              if (i .le. 0) goto 999
              call isanrg(icond,0,2,' ','dos mode',.true.)
C           lcond = .true.
C         Number of energy points
            elseif (sopts(j1:j1+4) .eq. 'npts=')  then
              m = 0
              i = parg('npts=',2,sopts(j1:),m,len(sopts(j1:)),
     .        dc//' ',1,1,iv,npts)
              if (i .le. 0) goto 999
C         DOS window
            elseif (sopts(j1:j1+6) .eq. 'window=')  then
              m = 0
              i = parg('window=',4,sopts(j1:),m,len(sopts(j1:)),
     .        ', '//dc,2,2,iv,drange)
              if (i .le. 1) goto 999
              lrange = .true.
C         First direction vector (icond=1,2)
            elseif (sopts(j1:j1+3) .eq. 'vec=')  then
              m = 0
              i = parg('vec=',4,sopts(j1:),m,len(sopts(j1:)),
     .        ', '//dc,2,3,iv,cvec(1))
              if (i .le. 0) goto 999
C         Second direction vector (icond=2)
            elseif (sopts(j1:j1+4) .eq. 'vec2=')  then
              m = 0
              i = parg('vec2=',4,sopts(j1:),m,len(sopts(j1:)),
     .        ', '//dc,2,3,iv,cvec(4))
              if (i .le. 0) goto 999
C         Generate total DOS (no partial decomposition into channels)
            elseif (sopts(j1:j2) .eq. 'totdos')  then
              totdos = .true.
C         DOS for a prescribed list of bands only
            elseif (sopts(j1:j1+5) .eq. 'bands=') then
              if (j1+6 .gt. j2) call rx('ASADOS: bad list')
              call mkils0(sopts(j1+6:j2),nblst,iblst)
              if (nblst .gt. 100) call rx('increase size of iblst')
              call mkilst(sopts(j1+6:j2),nblst,iblst)
C         DOS for a prescribed list of classes only
            elseif (sopts(j1:j1+6) .eq. 'classes') then
              if (nsite .ne. 0) call rx('ASADOS: specification of class '
     .        //'list not compatible with --pdos')
              j1 = j1+7
              m = j1
              i = parg('style=',2,sopts,m,len(sopts),sopts(j1:j1)//' ',1,
     .        1,iv,lstyle)
              if (i .lt. 0) then
                goto 999
              elseif (i .eq. 0) then
                j1 = j1+1
              else
                j1 = m+2
              endif
              call clist(lstyle,sopts(j1:j2),w(oclabl),z,nclass,nlstc,
     .        lsite)
            else
              goto 999
            endif
            goto 10
          endif
        endif

      endif

C --- Post-options initialization ---
C     lcond=T when contet branch is needed
      lcond = lcond .or. icond .ne. 0 .or. totdos .or. nblst .ne. 0

C     List of irreducible qp, tetrahedra, ipq
      call mkqp(sctrl,sbz,slat,ltet,lcond,1,0)

      i_copy_size=size(sbz%nkabc)
      call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
      nkp=sbz%nkp
      ntet=sbz%ntet
      oidtet=sbz%oidtet
      owtkp=sbz%owtkp


C ... bzmap setup
      if (bzmap) then
        nsp = int(sctrl%nspin)

C   ... Read in BZ plane(s) from file BZPL
        nfilqp = fopno('BZPL')
        rewind nfilqp
        nwmx = 0
        nkp = 0
        npln = 0
        nqmx = 0
C   ... Determine number of planes and energy windows for BZ map
   95   read (nfilqp,*,err=98,end=98) nwin,np0
        nwmx = max(nwin,nwmx)
        do  96  i = 1, nwin
          read (nfilqp,*,err=98,end=98) xx,xx
   96   continue
        do  97  i = 1, np0
          read (nfilqp,*,err=98,end=98) xx,xx,xx,xx,xx,n1,
     .    xx,xx,xx,xx,xx,n2,xx
          nkp = nkp + n1*n2
          npln = npln + 1
          nqmx = max(n1*n2,nqmx)
   97   continue
        goto 95
   98   call rxx(nkp .eq. 0,prgnam//': empty or badly formed BZPL file')
C   ... Get BZ plane(s) and BZ sampling parameters
        allocate(nw_iv(npln))
Ctakao_ZeroClear_NotRequiered          if (npln<0) nw_iv(:)=0

        allocate(nq_iv(2*npln))
Ctakao_ZeroClear_NotRequiered          if (2*npln<0) nq_iv(:)=0

        allocate(vx_rv(3*npln))
Ctakao_ZeroClear_NotRequiered          if (3*npln<0) vx_rv(:)=0.0d0

        allocate(vy_rv(3*npln))
Ctakao_ZeroClear_NotRequiered          if (3*npln<0) vy_rv(:)=0.0d0

        allocate(x_rv(2*npln))
Ctakao_ZeroClear_NotRequiered          if (2*npln<0) x_rv(:)=0.0d0

        allocate(y_rv(2*npln))
Ctakao_ZeroClear_NotRequiered          if (2*npln<0) y_rv(:)=0.0d0

        allocate(z_rv(npln))
Ctakao_ZeroClear_NotRequiered          if (npln<0) z_rv(:)=0.0d0

        allocate(ew_rv(abs(-2*npln*nwmx)))
        if (-2*npln*nwmx<0) ew_rv(:)=0.0d0

        call defdr(oqp, 3*nkp)
        call getbzp ( nfilqp , npln , nwmx , nkp , nw_iv , ew_rv 
     .  , nq_iv , vx_rv , vy_rv , x_rv , y_rv , z , w 
     .  ( oqp ) )

        wtbzm = 2d0 / nsp
C        if (preskp) then
C          call defdr(owtkp,nkp)
C          call dcopy(nkp,wtbzm*nsp,0,w(owtkp),1)
C        endif
      endif

      if (.not. ltet) call info2(20,1,1,' Integration by'//
     .' Methfessel-Paxton sampling,  order=%i  width=%d',norder,width)

C ... Open file containing bands and dos weights
      if (wtfn .eq. ' ') then
        wtfn = 'MOMS'
        if (tbdos) wtfn = 'BAND'
      endif
      call info0(20,1,0,
     .' ASADOS: reading weights from file '//wtfn//'%a')
      call info2(20,0,0,
     .'         expecting file to be resolved by '//
     .'%?#n==0#site##%-1j'//
     .'%?#(n==1|n==2)#l##%-1j'//
     .'%?#n==4#l and m##',ldos,0)
      nfilem = fopno(wtfn)

C ... Get dimensions (nsp,nspc,nband,nfstg); allocate memory for eband
      call iomomq(nfilem,0,nlf,nsp,nspc,ndum,nband,nfstg,j,ndum,
     .nchan,nchan,i,w,w,w,w,efermi,vmtz)
      call isanrg(nlf,nl,nl,' asados (warning)','file nl',.false.)
      nspd = nsp
C     Empirical TB stores stuff a little differently
      if (tbdos) then
        if (TBU) then
          nsp = 2
          nspc = 1
        else
          nband = nband*nsp
          nsp = 1
          nspc = 1
        endif
      endif
      nspx = nsp / nspc
      nbandx = nband*nspc
      allocate(eband_rv(nband*nsp*nkp))
Ctakao_ZeroClear_NotRequiered        if (nband*nsp*nkp<0) eband_rv(:)=0.0d0


C ... Determine number of DOS channels (nchan=0 -> total dos only)
      read (nfilem) nchan,ndum
      call info2(20,0,0,'         file has %i channel(s)',nchan,0)
      if (nchan .ne. 0) then
        if (mod(nfstg/10,10) .eq. 0)
     .  call rx(prgnam//'%a:  moments file missing DOS wts')
      else
        totdos = .true.
      endif
      if (nlstc .ne. 0) then
        nchds = nlstc*nl
      else
        nchds = nchan
      endif
C     No decomposition of dos; nl is no longer relevant
      if (totdos) then
        nchan = 1
        nchds = 1
        nl = nlf
      endif

C ... Read npts, emin, emax if not already available
      if (.not. bzmap) then
        emin = drange(1)
        emax = drange(2)
        if (npts .eq. 0 .or. .not. lrange) then
          if (npts .eq. 0) npts = 501
          i = awrite(' Enter npts (def=%i), emin and emax (def='//
     .    '%;3d,%;3d): ',outs,80,0,npts,emin,emax,0,0,0,0,0)
          call cwrite(outs,0,i-1,0)
          read (*,*) npts, emin, emax
        endif
      endif

C ... Allocate arrays for dos and moms if tetrahedron
      if (ltet) then
        allocate(doswt_rv(abs(-nchds*nbandx*nsp*(nkp+1))))
        if (-nchds*nbandx*nsp*(nkp+1)<0) doswt_rv(:)=0.0d0

        nfstg = 11
        if (nspc .eq. 2) nfstg = 21
C       Suppress projection of unity into parts for total dos
        if (totdos) then
          nfstg = 1
          call dvset ( doswt_rv , 1 , nchds * nbandx * nsp * ( nkp + 
     .    1 ) , 1d0 )

        endif
      else
        allocate(doswt_rv(nchds*nbandx*nspc))
Ctakao_ZeroClear_NotRequiered          if (nchds*nbandx*nspc<0) doswt_rv(:)=0.0d0

        nfstg = 1
      endif
      if (bzmap) then
        allocate(dos_rv(nqmx*nsp*nchds*nwmx*npln))
Ctakao_ZeroClear_NotRequiered          if (nqmx*nsp*nchds*nwmx*npln<0) dos_rv(:)=0.0d0

      else
        allocate(wk_rv(npts))
Ctakao_ZeroClear_NotRequiered          if (npts<0) wk_rv(:)=0.0d0

        allocate(dos_rv(npts*nsp*nchds))
Ctakao_ZeroClear_NotRequiered          if (npts*nsp*nchds<0) dos_rv(:)=0.0d0

      endif

C ... Make sure missing bands are high enough
      call dvset ( eband_rv , 1 , nband * nsp * nkp , 9d9 )

C ... Ensure nl,nsp,nspc,nkp,nband match; read in bands (& moms if tet.)
      if (nchds .eq. nchan) then
        call iomomq ( nfilem , 32 , nl , nsp , nspc , nkp , nband , nfstg 
     .  , j , nband , nchan , nchan , nevmx , eband_rv , doswt_rv 
     .  , doswt_rv , w , efermi , vmtz )

      else
        call iomomx ( nfilem , 32 , nl , nsp , nspc , nkp , nband , nfstg 
     .  , j , nband , nchan , nchan2 , nchds , nevmx , nlstc , lsite 
     .  , eband_rv , w , doswt_rv , w , efermi , vmtz )

      endif
C     call dfdump(w(oeband),nband*nsp*nkp,-66)
      if (j .ne. nkp) call rx(prgnam//': moments file missing qpts')
      if (efmax-.2d0 .lt. emax) call info2(10,1,0,' ASADOS'//
     .'%a (warning) efmax (=%;3d) is not >> emax (=%;3d)',efmax,emax)

C ... Printout switches
      call awrit4('%N ASADOS:  make dos for %i points from %i bands in'
     .//' window (%;3d,%;3d)',' ',80,stdo,npts,nevmx,emin,emax)
      outs = ' '
      call awrit6(' options '//
     .'%?#n# mode: '//modstr(icond)//'%a;##'//
     .'%?#(n==1|n==2)#%b,v=(%3;3d);#%j#'//
     .'%?#n# totdos;##'//
     .'%?#n# idos;##'//
     .'%?#n# tbdos;##'//
     .'%b %b',outs,80,0,
     .icond,icond,cvec(1),isw(totdos),lidos,tbdos)
      if (outs .ne. ' options') call awrit0('%a',outs,-80,-stdo)
      if (nlstc .ne. 0) call awrit2('%10fclasses: %n:1i',' ',
     .scrwid,stdo,nlstc,lsite)
      if (nblst .ne. 0) call awrit2('%10fbands: %n:1i',' ',
     .scrwid,stdo,nblst,iblst)

C ... Sanity checks
      if (nlstc .ne. 0 .and. totdos)
     .call rx(' class list is incompatible with totdos')
      do  i = 1, nblst
        if (iblst(i) .gt. nband .or. iblst(i) .le. 0)
     .  call fexit2(-1,111,' Exit -1  band %i '
     .  //'exceeds number of available bands (%i)',iblst(i),nband)
      enddo
C ... cvec(2) reverts to cvec(1) if it wasn't defined.
      if (ddot(3,cvec(4),1,cvec(4),1) .eq. 0)
     .call dcopy(3,cvec(1),1,cvec(4),1)
C ... Normalize cvec(1),cvec(2)
      if (icond .eq. 1 .or. icond .eq. 2) then
        if (ddot(3,cvec(1),1,cvec(1),1) .eq. 0) call
     .  rx(' direction vector required for this mode but not defined')
        call dscal(3,1/sqrt(ddot(3,cvec(1),1,cvec(1),1)),cvec(1),1)
        call dscal(3,1/sqrt(ddot(3,cvec(4),1,cvec(4),1)),cvec(4),1)
      endif

C --- Make BZ maps (integrated DOS vs k) ---
      if (bzmap) then
        call rxx(nchds.ne.nchan,'maps not set for subchannels')
        call nosbzm ( nfilem , nbandx , nsp , nspc , nchan , lmet , norder 
     .  , width , nkp , wtbzm , eband_rv , doswt_rv , npln , nwmx 
     .  , nqmx , nw_iv , ew_rv , nq_iv , dos_rv )

        call bzmio ( dos_rv , nchan * nsp , npln , nwmx , nqmx , nw_iv 
     .  , ew_rv , nq_iv , vx_rv , vy_rv , x_rv , y_rv 
     .  , z , nkp , nspd , lmet , norder , width , wtbzm , efermi , - 
     .  fopnn ( 'BZMP' ) )

C --- Ballistic conductivity ---
      elseif (lcond) then
        call rxx(nchds.ne.nchan,'cond not set for subchannels')
        if (.not. ltet) call
     .  rx('LMDOS: conductivity only with tetrahedron integration')
C   ... To save memory, allocate dos later
        if (allocated(wk_rv)) deallocate(wk_rv)
        if (allocated(dos_rv)) deallocate(dos_rv)

C   ... for debugging
C        call upack('bz nkp oqp',sbz,nkp,oqp,0,0,0)
C        iblst(1) = 0
C        nblst = 1
C        nevmx = 1
C        nband = 1
C        nchan = 1
C        nbandx = 1
C        nsp = 1
C        nspx = 1
C        call snot(w(oqp),nkp,w(oeband),plat)
C   ... Remake qp with reduced symmetry; assume no grp ops for now
        npgrp = 1

        slat%npgrp=npgrp

        call awrit2('          remake qp with %i symop%?#n>1#s##',
     .  ' ',80,stdo,npgrp,npgrp)
        call mkqp(sctrl,sbz,slat,ltet,F,0,2)

        nkp=sbz%nkp
        i_copy_size=size(sbz%nkabc)
        call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
        oqp=sbz%oqp
        owtkp=sbz%owtkp
        ogstar=sbz%ostar


        ntet=sbz%ntet
        oidtet=sbz%oidtet
        oipq=sbz%oipq

        allocate(dos_rv(npts*nsp*nchan))
Ctakao_ZeroClear_NotRequiered          if (npts*nsp*nchan<0) dos_rv(:)=0.0d0

        if (nblst .eq. 0) nblst = nevmx
        call contet ( icond , nbandx , nsp , nspx , nblst , nchan , nkxyz 
     .  ( 1 ) , nkxyz ( 2 ) , nkxyz ( 3 ) , ntet , w ( oidtet ) , w ( 
     .  oqp ) , w ( oipq ) , w ( ogstar ) , iblst , eband_rv , cvec 
     .  ( 1 ) , doswt_rv , npts , emin , emax , lidos , wk_rv , 
     .  dos_rv )

        del = 0d0
C --- make DOS or NOS ---
      else
        if (ltet) then
          call dostet ( nbandx , nsp , nspx , nevmx , nchds , nkxyz ( 1 
     .    ) , nkxyz ( 2 ) , nkxyz ( 3 ) , ntet , w ( oidtet ) , eband_rv 
     .    , doswt_rv , npts , emin , emax , lidos , wk_rv , dos_rv 
     .    )

          del = 0d0

C   ... dosspl reads bands and moms on the fly to save work space
        else
          call rxx(nchds.ne.nchan,'sampline not set for subchannels')
          call dosspl ( nfilem , nbandx , nsp , nspc , nchan , norder , 
     .    width , nkp , w ( owtkp ) , eband_rv , doswt_rv , npts 
     .    , emin , emax , lidos , wk_rv , dos_rv )

          del = norder+width
        endif
        if ( nspd .gt. nsp ) call dscal ( npts * nsp * nchan , .5d0 , 
     .  dos_rv , 1 )

      endif
      if (.not. bzmap) then
        i = nsp
        if (lncol .gt. 0) i = 1
        call dosio ( dos_rv , npts * ( 1 + nsp - i ) , nchds , npts 
     .  , nchds , emin , emax , i , efermi , del , 1 , - fopnn ( 'DOS' 
     .  ) )

      endif
      if (allocated(dos_rv)) deallocate(dos_rv)
      if (allocated(doswt_rv)) deallocate(doswt_rv)


C --- List channels in DOS file ---
      if (iprint() .ge. 30 .and. .not. totdos .and.
     ..not. cmdopt('--cls',5,0,strn)) then
        if (nsp .eq. 1) print 332
        if (nsp .eq. 2) print 332, 'spin-2'
  332   format(/' Channels in dos file generated by LMDOS:'
     .  /' site class label   spin-1':23x,a)
        if (nsite .ne. 0) then
          imax = nbas
        else
          imax = nclass
        endif

        il = 1
        ichds = 0
C     Loop over all sites (or classes)
        do  40  i = 1, imax
          if (lsite(1) .ne. -1) then
            if (i .ne. lsite(il)) goto 40
          endif

          if (nsite .ne. 0) then
            ib = i
            ic = ival(w(oipc),ib)
          else
            ic = i
            ib = iclbas(i,w(oipc))
          endif

C       Make a new nlo if channels are site-dependent
          if (lfp) then
C         Global upper bound to lmax
            nl = nll

            is=ssite(ib)%spec

            lmax = int(sspec(is)%lmxa)

            if ( cmdopt ( '--mull' , 6 , 0 , strn ) ) lmax = int(sspec(is)%lmxb)

            nll = min(lmax+1,nl)
            nlo = nll
            if (ldos .eq. 0) nlo = 1
            if (ldos .eq. 2) nlo = nll
            if (ldos .eq. 4) nlo = nll*nll
          endif

          do i_spacks=ib,ib
            call spacks_copy('u',ssite(i_spacks)%clabel,ib,ib,clabl,i_spacks)
          enddo

          outs = ' '
          write(outs,333) i, ic, clabl
  333     format(2i5,3x,a8)
          j = 0
          do  43  l = 0, nlo-1
            j = j+1
            iblst(j) = ichds+nsp*l+1
   43     continue
          call ilst2a(iblst,j,outs(21:))
C        do  44  l = 0, nlo-1
C   44   call awrit1('%a%i,',outs(21:),80,0,ichds+nsp*l+1)
C        call awrit0('%a%b ',outs(21:),80,0)
          if (nsp .eq. 2) then
            j = 0
            do  45 l = 0, nlo-1
              j = j+1
              iblst(j) = ichds+nsp*l+2
   45       continue
            call ilst2a(iblst,j,outs(50:))
C          do  46  l = 0, nlo-1
C   46     call awrit1('%a%i,',outs(50:),80,0,ichds+nsp*l+2)
C          call awrit0('%a%b ',outs(50:),80,0)
          endif
          call awrit0('%a',outs,80,-stdo)
          il = il+1

          ichds = ichds + nlo*nsp
   40   continue
      endif
Cki#error, have return with len(w_varlist)>0 at line 771
      if (allocated(eband_rv)) deallocate(eband_rv)
      if (allocated(ew_rv)) deallocate(ew_rv)
      if (allocated(z_rv)) deallocate(z_rv)
      if (allocated(y_rv)) deallocate(y_rv)
      if (allocated(x_rv)) deallocate(x_rv)
      if (allocated(vy_rv)) deallocate(vy_rv)
      if (allocated(vx_rv)) deallocate(vx_rv)
      if (allocated(nq_iv)) deallocate(nq_iv)
      if (allocated(nw_iv)) deallocate(nw_iv)
      return


  999 continue
      call rxs('asados: failed to parse dos options:',sopts)
      end subroutine asados

C      subroutine snot(qp,nkp,eband,plat)
C      implicit none
C      integer nkp
C      double precision qp(3,nkp),eband(nkp),plat(3,3)
C      integer iq
C      double precision q(3),qlat(3,3),xx
C
C      call dinv33(plat,1,qlat,xx)
C      do  iq = 1, nkp
C        call dcopy(3,qp(1,iq),1,q,1)
C        call shorbz(qp(1,iq),q,qlat,plat)
C        eband(iq) = q(1)**2+q(2)**2+q(3)**2
CC        print 333, iq, qp(1,iq),qp(2,iq),qp(3,iq), q, eband(iq)
CC  333   format(i5,6f12.6,2x,f12.6)
C      enddo
C      end


