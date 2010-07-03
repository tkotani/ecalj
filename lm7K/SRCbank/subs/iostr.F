      logical function iostr(lio,name,nl,nbas,nkap,kap2,itral,ckbas,
     .lmaxw,nitab,oalpha,oiax,onpr,os)
C- Real space structure constants file read or write
C ----------------------------------------------------------------------
Ci Inputs
Ci   lio   :specifies conditions for file i/o (condensed into digits)
Ci         :digit  specification
Ci         :1-10s   bits in 1-10s digit control program flow
Ci                  0 for file read
Ci                  1 for file write
Ci                  2 (read only) No header matching requirement;
Ci                    header parms returned into passed arguments
Ci                  4 read/write file header only and exit.
Ci                  8 Add this number to allocate arrays before reading
Ci                    or to release arrays after writing
Ci                 16 file is formatted (not implemented)
Ci         :100-1000 bits in these digits contain info about structure of s
Ci                  1 s is complex (otherwise s is real)
Ci                  2 s is stored in slices according to site.
Ci                    row ordering that of the cluster determined
Ci                    by iax table.
Ci                    Default storage (2's bit zero):
Ci                    s is stored by (nds,nds) blocks
Ci                  4 s is the structure matrix
Ci                    Default: s is the matrix of expansion coffs
Ci                  8 file contains energy derivative of s
Ci         :10000-100000 info about the functions s corresponds to
Ci                  0 2nd generation Hankel functions with Andersen's
Ci                    conventions
Ci                  1 NMTO Hankel functions with Tank and Andersen
Ci                    conventions: alpha=alpha(2,2,nl**2,nbas,nkap)
Ci                    and s=s(nl**2,nl**2,nsite,nkap)
Ci                  2 Conventional Hankel functions, Methfessel defs
Ci   name  :file name holding strux
Ci   nl    :1 + global-maximum-l for (nl**2,nl**2) format
Ci          total dimension of s for vector format
Ci   nbas  :number of atoms in the basis
Ci   nkap  :number of energies for which strux are calculated
Ci   kap2  :Hankel energies
Ci   itral :characterizes structure matrix transformations,
Ci          (Used by NMTO only; see mktra2.f)
Ci   ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Ci   lmaxw :maximum l for Watson-sphere
Cio Inputs/Outputs:
Cio  nitab :total number of inequivalent pairs in iax table (strscr.f)
Cio         Unused for vector format
Cio  oalpha:offset for tight-binding screening parameters
Cio  oiax  :offset for iax array
Cio  onpr  :offset for number of sites in each cluster
Cio  os    :offset for screened structure constant matrix
Cb Bugs
Cb   Only lio 10s digit only 0 for MPI
Cr Remarks
Cr   If file read, allocates memory for and reads alpha,iax,npr,s
Cr   Input parameters nl,nbas,bas,nkap,lmaxw,itral must match file
Cr   for iostr to complete file read.
Cu Updates
Cu   06 Aug 06 Redesigned to work with 2-kappa strux
Cu   8 Jun 04  (MPI) read from master node only and broadcast
Cu             (ccomp with MASTERIO.) (Warning: not checked.)
Cu   1 Aug 98  revised for 3rd generation lmto
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character name*(*)
      integer itral,lio,nitab,nl,nbas,nkap,lmaxw,oalpha,oiax,onpr,os
      double precision kap2(nkap),ckbas
C ... Local parameters
      character locnam*8
      logical iostr1,iostr2,iostr4,lerr
      integer nlr,nbr,ifi,ipr,nkr,nttab,lgunit,itralr,lmaxwr,niax,ncplx,
     .lio01,lio23,lio45,getdig,nds,lscoff,nmto,lior,nblk,nkaps,nkapn
      integer mpipid,procid,vsn
      double precision tiny,ckbasr
C     double precision sttime,entime,MPI_WTIME
      parameter (tiny=1d-8,niax=10,vsn=1)
C ... Heap
      integer w(1)
      common /w/ w

C --- Setup information ---
      call getpr(ipr)
      locnam = name
      lio01 = mod(lio/1,100)
      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      ncplx = 1+getdig(lio23,0,2)

      procid = mpipid(1)
C     Leading dimension of 'alpha' array depends on whether it
C     corresponds to 2nd gen alpha or tral matrix

C --- File write ---
      if (mod(lio,2) .ne. 0) then

C       Use nds=nl**2 for now, until nds is passed
        nds = nl*nl
        nmto = getdig(lio45,0,2)
        nkaps = nkap
        nkapn = 1
        if (nmto .eq. 1) then
          nkaps = 1
          nkapn = nkap
        endif
        lscoff = nds*nbas*nkaps**2*nkapn

C   ... Write header
        nttab = w(onpr+nbas)
        if (.not. iostr1(lio,name,itral,nttab,nitab,nl,nbas,
     .  nkap,lmaxw,ckbas,ifi)) goto 99
C   ... Write strux, block format
        if (getdig(lio23,1,2) .eq. 0) then
          nblk = nitab*nkaps**2*nkapn
          call rxx(getdig(lio23,0,2).eq.1,'iostr not ready for complex')
          if (.not. iostr2(-ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .    w(onpr),w(oiax),kap2,w(oalpha),w(os))) goto 99
C   ... Write strux, vector format
        else
          if (.not. iostr4(-ifi,nbas,nttab,ncplx,nl,nkap,
     .    w(onpr),w(oiax),kap2,w(os))) goto 99
        endif

        if (ipr .ge. 40) call awrit4(' IOSTR: wrote to file '
     .  //locnam//'%a %i sites'//
     .  '%?#n#, %i inequivalent sites#%j#'//
     .  '%?#n>1#, %-1j%i energies##',
     .  ' ',80,lgunit(1),nttab,nttab-nitab,nitab,nkap)
        call fclose(ifi)

C --- File read ---
      else
C   ... Open file, and read header
#if MASTERIO
        if (procid .eq. 0) then
#endif
          lior = 0
          lerr = .not. iostr1(lior,name,itralr,nttab,nitab,nlr,nbr,
     .    nkr,lmaxwr,ckbasr,ifi)
          if (lerr) goto 99
C       Copy file parameters to those in argument list
          if (getdig(lio01,1,2) .eq. 1) then
            nl = nlr
            nbas = nbr
            nkap = nkr
            itral = itralr
            ckbas = ckbasr
            lmaxw = lmaxwr
          endif

          call isanrg(nlr,nl,nl,'file: IOSTR:','nl',.true.)
          call isanrg(nbr,nbas,nbas,'file: IOSTR:','nbas',.true.)
          call isanrg(nkr,nkap,nkap,'file: IOSTR:','nkap',.true.)
          call fsanrg(ckbasr,ckbas,ckbas,tiny,'file: IOSTR:','cksum',
     .    .true.)
C   ... This is only a warning, not fatal error
          call isanrg(lmaxwr,lmaxw,lmaxw,'file: IOSTR:','lmaxw',.false.)

#if MASTERIO
        endif

        call mpibc1(lerr,1,0,.false.,'iostr','lerr')
        call mpibc1(nl,1,2,.false.,'iostr','nl')
        call mpibc1(nbas,1,2,.false.,'iostr','nbas')
        call mpibc1(nkap,1,2,.false.,'iostr','nkap')
        call mpibc1(itral,1,2,.false.,'iostr','itral')
        call mpibc1(lmaxw,1,2,.false.,'iostr','lmaxw')
        call mpibc1(nttab,1,2,.false.,'iostr','nttab')
        call mpibc1(nitab,1,2,.false.,'iostr','nitab')
#endif

C   ... If failed to read header, give up
        if (lerr) goto 99

C       Use nds=nl**2 for now, until nds is passed
        nds = nl*nl
        nmto = getdig(lio45,0,2)
        nkaps = nkap
        nkapn = 1
        if (nmto .eq. 1) then
          nkaps = 1
          nkapn = nkap
        endif

C   ... Allocate memory for arrays
        if (getdig(lio01,3,2) .eq. 1) then
          lscoff = nds*nbas*nkaps**2*nkapn
          call defdr(oalpha,lscoff)
          call defi(oiax,niax*nttab)
          call defi(onpr,nbas+1)
          call defdr(os,nds**2*nitab*nkaps**2*ncplx*nkapn)
        endif

C   ... If only header is to be read, exit
        if (getdig(lio01,2,2) .eq. 1) then
          iostr = .true.
Cki#error, have return with len(w_varlist)>0 at line 196
          return
        endif

C   ... Read strux, block format
        if (getdig(lio23,1,2) .eq. 0) then
          call rxx(getdig(lio23,0,2).eq.1,'iostr not ready for complex')
#if MASTERIO
          if (procid .eq. 0) then
#endif
            nblk = nitab*nkap**2
            if (nmto .eq. 1) nblk = nitab*nkap
            lerr = iostr2(ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .      w(onpr),w(oiax),kap2,w(oalpha),w(os))
#if MASTERIO
          endif
C         sttime = MPI_WTIME()
C         call info0(20,0,-1,' iostr: MPI broadcast strx data ...')
          call mpibc1(lerr,1,0,.false.,'iostr','lerr')
          call mpibc1(w(onpr),nbas+1,2,.false.,'iostr','ntab')
          call mpibc1(w(oiax),niax*nttab,2,.false.,'iostr','iax')
          call mpibc1(kap2,nkap,4,.false.,'iostr','kap2')
          call mpibc1(w(oalpha),lscoff,4,.false.,'iostr','alp')
          call mpibc1(w(os),nds**2*nitab*nkap,4,.false.,'iostr','s')
C         entime = MPI_WTIME()
C         call info2(20,0,0,'  took %;1d sec',(entime-sttime),0)
#endif
C   ... Read strux, vector format
        else
#if MASTERIO
          if (procid .eq. 0) then
#endif
            lerr = iostr4(ifi,nbas,nttab,ncplx,nl,nkap,
     .      w(onpr),w(oiax),kap2,w(os))
          endif
#if MASTERIO
          call mpibc1(lerr,1,0,.false.,'iostr','lerr')
          call mpibc1(w(onpr),nbas+1,2,.false.,'iostr','ntab')
          call mpibc1(w(oiax),niax*nttab,2,.false.,'iostr','iax')
          call mpibc1(kap2,nkap,4,.false.,'iostr','kap2')
          call mpibc1(w(os),nds**2*nitab*nkap*ncplx,4,.false.,'iostr',
     .    's')
        endif
#endif
        if (.not. lerr) goto 99
        if (ipr .ge. 100) call awrit3(' IOSTR: read %i sites'//
     .  '%?#n#, %i inequivalent sites##',' ',80,
     .  lgunit(1),nttab,nttab-nitab,nitab)
      endif
      iostr = .true.
Cki#error, have return with len(w_varlist)>0 at line 245
      return

C ... Error handling
   99 continue
      if ((lio/10) .eq. 0) call rxs('IOSTR: file mismatch, file ',name)
      if (ipr .ge. 10)
     .print *, 'IOSTR (warning): failed to read file ',name
      iostr = .false.

      end
      logical function iostr1(lio,name,itral,nttab,nitab,nl,nbas,
     .nkap,lmaxw,ckbas,ifi)
C- Low level structure constants file open and header I/O
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  lio   :1s digit 0 for read, otherwise write
Cio        :On file READ, returns with file contents of lio
Cio  name  :file name
Cio  itral :characterizes structure matrix transformations (mktral.f)
Cio  nttab :total number of pairs in neighbor and iax (pairs.f)
Cio  nitab :total number of inequivalent pairs in iax table (strscr.f)
Cio  nl    :global maximum l + 1
Cio  nbas  :number of atoms in the basis
Cio  nkap  :number of kinetic energies for which strux are calculated
Cio  lmaxw :maximum l for Watson-sphere
Cio  ckbas :checksum for basis vectors; make with cksum(bas,3*nbas)
Co Outputs
Co   ifi   :file logical unit number
Co   iostr1:T if read successful, false otherwise
Cr Remarks
Cr   Reads screened strux header data from file, versions 0 and 1.
Cr
Cr   Structure of file(including header) version 2
Cr     record    content
Cr       1         -99  vsn  lio
Cr       2         itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
Cr   Structure of file(including header) version 1
Cr     record    content
Cr       1         -99  vsn  ihold
Cr       2         itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
Cr       ... end of header
Cr       3         ntab,iax,kap2,X,s
Cr                 where X=alpha or tral, depending on ihold
Cr
Cr   Structure of file(including header) version 0
Cr     record    content
Cr       1         itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
Cr       ... end of header
Cr       2         ntab,iax,kap2,alpha,s
Cu Updates
Cu   8 Jun 04 Implement version 1.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character name*(*)
      integer lio,itral,nttab,nitab,nl,nbas,nkap,lmaxw
      double precision ckbas
C ... Local parameters
      integer ifi,fopnx
      integer i,vsn,rdvsn
      character strn*120
      parameter (vsn=2)

      ifi = fopnx(name,100,16+8+4+0,-1)
      rewind ifi

C ... File READ
      if (mod(lio,2) .eq. 0) then

C ... Check for version number.  If nonexistent, vsn=0
        read(ifi,end=99,err=99) i,rdvsn,lio
        if (i .eq. -99 .and. rdvsn .eq. 2) then
          read(ifi,end=99,err=99) itral,nttab,nitab,nl,nbas,nkap,lmaxw,
     .    ckbas
        elseif (i .eq. -99 .and. rdvsn .eq. 1) then
          read(ifi,end=99,err=99) itral,nttab,nitab,nl,nbas,nkap,lmaxw,
     .    ckbas
        else
          rdvsn = 0
          lio = 0
          backspace ifi
          read(ifi,end=99,err=99) itral,nttab,nitab,nl,nbas,nkap,lmaxw,
     .    ckbas
        endif

        strn = ' IOSTR : read file='//name//
     .  '%a lio=%i nkap=%i nl=%i nttab=%i'
        i = lio

C ... File WRITE
      else
C       This is version 0
C       write(ifi) itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
C       This is version 1
C        write(ifi) -99,vsn,ihold
C        write(ifi) itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas
C       This is version 2
        write(ifi) -99,vsn,lio-1
        write(ifi) itral,nttab,nitab,nl,nbas,nkap,lmaxw,ckbas

        strn = ' IOSTR : write file='//name//
     .  '%a lio=%i nkap=%i nl=%i nttab=%i'
        i = lio-1


      endif

      call info5(110,1,0,strn,i,nkap,nl,nttab,0)

      iostr1 = .true.
      return
   99 iostr1 = .false.
      end
      logical function iostr2(ifi,nbas,nttab,nblk,lscoff,nds,nkap,
     .ntab,iax,kap2,scoffs,s)
C- Low level structure constants I/O
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   nbas  :number of atoms in the basis
Ci   nttab :total number of pairs in neighbor and iax (pairs.f)
Ci   nblk  :number of blocks of s to write
Ci   lscoff:size of scoffs matrix
Ci   nds   :leading dimensions of s
Ci   nkap  :number of kinetic energies for which strux are calculated
Cio Inputs/Outputs
Cio  ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairs.f)
Cio  iax   :array of parameters containing info about each pair
Cio  kap2  :interstitial kinetic energies
Cio  scoffs :tight-binding screening parameters
Cio  s     :real-space structure constant matrix
C ----------------------------------------------------------------------
C     implicit none
      integer ifi,nbas,nttab,nblk,nds,nkap,niax,ntab(nbas+1)
      parameter (niax=10)
      integer iax(niax,nttab),lscoff
      double precision kap2(nkap),scoffs(lscoff),s(nds,nds,nblk)

      iostr2 = .true.
      if (ifi .gt. 0) then
        read(ifi,end=99,err=99) ntab,iax,kap2,scoffs,s
      else
        write(-ifi) ntab,iax,kap2,scoffs,s
      endif
      return
   99 continue
      iostr2 = .false.
      end
      logical function iostr4(ifi,nbas,nttab,ncplx,ndim,nkap,ntab,
     .iax,kap2,s)
C- Low level structure constants I/O, vector format
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file handle
Ci   nbas  :number of atoms in the basis
Ci   nttab :total number of pairs in neighbor and iax (pairs.f)
Ci   nl2   :leading dimensions of s
Ci   nkap  :number of kinetic energies for which strux are calculated
Cio Inputs/Outputs
Cio  ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairs.f)
Cio  iax   :array of parameters containing info about each pair
Cio  kap2  :interstitial kinetic energies
Cio  alpha :tight-binding screening parameters
Cio  s     :real-space structure constant matrix
C ----------------------------------------------------------------------
C     implicit none
      integer ifi,nbas,nttab,ncplx,ndim,nkap,niax,ntab(nbas+1)
      parameter (niax=10)
      integer iax(niax,nttab)
      double precision kap2(ncplx,nkap),s(ncplx,ndim,nkap)
      iostr4 = .true.
      if (ifi .gt. 0) then
        read(ifi,end=99,err=99) kap2,s
      else
        write(-ifi) kap2,s
      endif
      return
   99 continue
      iostr4 = .false.
      end

