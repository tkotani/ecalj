C#define MPI-SCHEDULED
      subroutine bndtb(sctrl,sarray,sbz,slat,sspec,nbas,nl,nspc,nsp,
     .  nsp1,lmx,idxdn,nclass,ips,ipc,nrclas,pv,force,nbmax,zval,mull,
     .  npln,nwmx,nqmx,nw,ew,nqpl,vxpl,vypl,xpl,ypl,zpl,wtbzm,nsite,iax,
     .  npr,xyzfrz,hrs,h0,dh,dhcf,vso,hso,srs,ds,dscf,pot0,iwk,eband,
     .  efermi,sumev,entrpy,f,thrpv,esite,rho,rholm,rhoc,rhon,zos,index,
     .  drhosl,ldim)
C- k-integration of tight-binding bands
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas, no. of atoms in the basis; nl, number of l's;
Ci   nspc (2 for coupled spins) nsp (2 for spin polarised TB)
Ci   nps1 is for dimensioning arrays (equivalent to nspx in LMASA)
Ci   lmx, max. l for classes;
Ci   if idxdn(l,i) > 1 then this L for ith class is excluded from basis;
Ci   nclass, no. of classes - atoms in the same class are
Ci   symmetry-related; the jth atom belongs to class ipc(j); 
Ci   nrclas(i), no. of atoms in the ith class;  force: TRUE for
Ci   forces; nbmax, max no. of bands;
Ci   zval, total number of valence electrons
Ci   mull, Mulliken DOS switch, see mkwttb; 
Ci   npln,nwmx,nqmx,nw,ew,nqpl,vxpl,vypl,xpl,ypl,zpl,wtbzm:
Ci   parameters for BZ maps (integrated DOS vs k), see getbzp and bzmio
Ci   nsite:   total number of neighbors in all clusters;
Ci   iax, neighbor lists; npr, see tbham;
Ci   hrs,h0,dh: real-space hamiltonian, and dh(x,y,z,r) its derivatives
Ci   vso:    table of spin-orbit parameters and hso is 
Ci   hso:    the spin-orbit hamiltonian
Ci   srs,ds: real-space overlap matrix, and ds(x,y,z,r) its derivatives
Ci   dhcf,dscf: the crystal field ham. and overlap derivatives
Ci   pot0: monopole potential at site R (eq. 7.81, Finnis)
Co Outputs
Co   iwk, number of orbitals for each atom; eband, energy bands; efermi;
Co   sumev, sum occupied e'vals; f(3,nbas), force on each atom from the
Co   entrpy, "entropy term" actually TS
Co   band structure; thrpv, 3PV from the band structure; esite, band
Co   energies for each atom; 
Co   rho, s,p, and d Mulliken charges;
Co   rholm, {lm} decomposed Mulliken charges;
Co   If STONER=T zos returns the local number of states and index points
Co   to the d-electron components.
Co   rho, rhoc and rhon are, respectively, local charges by l, 
Co   eigenvector products for TB-L and TB+U (see tbfrce and tbesel)
Co   drhosl: (c*_RL dS_RLR'L'/dR c_R'L' + cc) s-c TB+ovlp
Co          built in tbfrce
Co   ldim is passed back up for dimensioning drhosl in tbesel
Cf  Files:
Cf    BAND file has structure (compatability with ASA)
Cf    1.   nl 1 1 nkp ldim*nspc nfstg
Cf    For each qpt, the following records:
Cf    2.   nchan  nev (if nfstg nonzero)
Cf         eband      (if 1s   digit of nfstg nonzero)
Cf         doswt      (if 10s  digit of nfstg 2)
Cf    3.   efermi, 0d0
Cm MPI
Cm   Parallelise over the k-loop. Each process computes contributions
Cm   to efermi, sumev, f, thrpv, esite, rho. These need to be
Cm   reduced (merged) before returning. This has to be done using
Cm   additional buffers until MPI supports MPI_IN_PLACE. These buffers
Cm   are allocated locally. Thus FORTRAN90 is needed.
Cm   eband is assembled at different k by different processors. They
Cm   broadcast their bands to the other processors at the end of the
Cm   k loop.
Cm   It's best to give each process a contiguous set of k-points
Cm   to reduce the total number of broadcasts.
Cm
Cm   evm = T is set in line, which means the eigenvectors are kept in 
Cm   allocatable memory rather than written to disc if evd = T.
Cm   It is probably better to do this since the eigenvectors are
Cm   distributed across the processes so each process only needs
Cm   local memory for its own eigenvectors. On the other hand putting
Cm   them temporarily on disc as in the serial version may slow down
Cm   the application depending upon whether the processes have access
Cm   to local disc storage. evm may be set to F to achieve this if
Cm   required.  
Cu Updates
Cu   04 Jun 08 (ATP) Modificiations for gamma-only ; faster execution
Cu   10 Apr 02 Bug fix --- nfstg now set correctly
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   03 Feb 01 use ldos instead of lbzio
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
C#ifdefC MPI
C      include "mpif.h"
C      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*20 ext
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C      logical mlog
C#endif
      double precision sctrl(1),sarray(1),sbz(1),slat(1),sspec(1)
      integer nbas,nl,nsp,nspc,nsp1,nclass,nbmax,nx,mull,npln,nwmx,nqmx,
     .        nsite
      integer lmx(1),idxdn(0:nl-1,1),ips(1),ipc(1),nrclas(1),
     .  nw(npln),nqpl(2,npln),iax(1),npr(1),iwk(nbas),index(1)
      double precision zval,wtbzm,efermi,sumev,entrpy,thrpv
      double precision 
     .  ew(2,nwmx,npln),vxpl(3,npln),vypl(3,npln),xpl(2,npln),
     .  ypl(2,npln),zpl(npln),hrs(1),h0(1),dh(1),vso(1),hso(1),
     .  srs(1),ds(1),eband(nbmax,nsp1,1),f(3,nbas),esite(nbas,1),
     .  rho(nl,2,nbas),rholm(nl**2,2,nbas),rhoc(nl**2,nl**2,nbas),
     .  rhon(nl**2,nl**2,nbas,2),zos(1),dhcf(1),dscf(1),drhosl(1),
     .  pot0(nbas)
      logical pv,force,xyzfrz(3),evecs
C Heap:
      integer w(1)
      common /w/ w
C Local variables
C#ifdefC MPI
C      integer pid,i,jp,jkp,nkproc,ikproc,lgunit
C      integer,          dimension(:),       allocatable :: kpproc
C      double precision, dimension(:),       allocatable :: buffer
C      double precision, dimension(:,:,:,:), allocatable :: zll
C      double precision thrpv_buffer
C#endif
      integer ldos,lncol,lstnr(3),ltb,getef
      integer i,l,ikp,nfilet,nfileb,inxsh(0:3),ldim,idim,
     .  lihdim,ipr,iprx,nev,fopn,fopnT,iprint,icheck,ibas,indx,nchan,
     .  ic,ifi,nevhi,nevlo,nat0,nlm,mull0,mull1,fopno,i1mach,lm,ll,
     .  imfil,ndos0,ndos,nfstg,isp,ispx,ldimx,nbmaxx,ik0,
     .  nevmx0,j,iomoms
      integer lmet,nkabc(3),n1,n2,n3,nkp,ntet,oidtet,mpsord,norder,nevmx
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      integer nsgrp,oistab,oag,osymgr
      integer owtkb,odos,oidxsh,ozll,odoswt,oxm0,ozwk,owk,ofwk,osll,
     .  ozwk2,opdos,obc,obcos,opwk,obwk,ogd,ogd2,owkcf,oeband,oqp,owtkp,
     .  owklm,ozup,ozdn
C --- for reconstructing RS hamiltonian
      integer ohrs,oors,ontab,opos,oclabl,oalpha
      double precision dum,del,evlo,qp(3),dval,dosw(2),plat(3,3),alat,
     .  srnge,swidth,efmax,dosef(2),qval(2)
      logical metal,tetra,leig,doswt,fsym,evd,charge,rl,trh,lso,donly,
     .  lmol,lgamma,lgamsp,lov,bzmp,stoner,cryf,ocryf,rstrt,swnos,
     .  TBU,rn,bittst,evm,cmdopt
      character*80 outs
      character*45 multyp(0:3,0:1)
      data multyp
     . /' Partial DOS, resolved by class and l',
     .  ' Bond charge, resolved by class pairs',
     .  ' Bond charge, resolved by {class, l} pairs',
     .  ' Bond charge, resolved by {class, l, m} pairs',
     .  ' Partial DOS, resolved by atom, l, and m',
     .  ' Bond charge, resolved by atom pairs',
     .  ' Bond charge, resolved by {atom, l} pairs',
     .  ' Bond charge, resolved by {atom, l, m} pairs'/

C#ifdefC MPI
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      mlog = cmdopt('--mlog',6,0,outs)
C      evm = .true.
C#else
      evm = .false.
C#endif

C --- Setup ---
      call upack('array opos oclabl',sarray,opos,oclabl,0,0,0)
      call upack('ctrl ltb ldos lncol lstonr',sctrl,ltb,ldos,lncol,
     .  lstnr,0)
      call upack('bz nkabc nkp ntet oidtet lmet',sbz,nkabc,nkp,ntet,
     .  oidtet,lmet)
      call upack('bz ndos dosw oqp owtkp',sbz,ndos0,dosw,oqp,owtkp,0)
      call upack('bz n range w nevmx efmax',sbz,mpsord,srnge,
     .  swidth,nevmx,efmax)
      mpsord = isign(1,mpsord) * mod(iabs(mpsord),100)
      call upack('lat plat alat',slat,plat,alat,0,0,0)
      call upack('lat nsgrp oistab oag osymgr',slat,nsgrp,oistab,oag,
     .  osymgr,0)

      call rxx(bittst(ltb,2**12),'Old MRS theory no longer implemented')
      tetra  = ntet .gt. 0
      metal  = lmet .ne. 0
      stoner = lstnr(1) .ne. 0
C      doswt  = bittst(ldos,2) .or. bittst(ldos,4)
      doswt  = (mull .ge. 0)
      trh    = bittst(ltb,2**10)
      evd    = (bittst(ltb,2**8) .or. tetra) .and. metal
      charge = bittst(ltb,2**12) .or. bittst(ltb,2**11)
     .    .or. bittst(ltb,2**13) .or. bittst(ltb,2**15)
      rl     = bittst(ltb,2**15)
      lso    = bittst(lncol,4)
      donly  = bittst(ltb,2**6)
      cryf   = bittst(ltb,2)
      lov    = bittst(ltb,1)
      ocryf  = bittst(ltb,4)
      bzmp   = bittst(ldos,8)
C     T if just gamma point
      lgamma = bittst(ltb,2**17)
      lmol   = bittst(ltb,2**18)
      if (lmol) then
        lgamma = .true.
        evd = .false.
      endif

C spin polarised gamma only, or molecule:
      lgamsp = lgamma .and. nsp .eq. 2
      if (lgamsp) then
        metal = .true.
        evm = .true.
        evd = .false.
        call rxx(tetra,'BNDTB: for spin pol at gamma use TETRA=F')
      endif
      evecs  = (force .or. charge .or. trh .or. pv)
      TBU = bittst(ltb,2**13)
      rn = TBU
      if (rl .or. rn) trh = .true.
      if (stoner) doswt = .true.
      rstrt = cmdopt('-cont',5,0,outs)
      fsym = nsgrp .gt. 1
      if (lmol) fsym = .false.
      leig = force .or. doswt .or. charge .or. trh .or. pv
      swnos = ndos0 .lt. 0
      ndos = iabs(ndos0)
      nfstg = 1
      if (doswt) nfstg = 11
      call rxx(lso .and. nspc .ne. 2,'BNDTB: must set NSPIN=2 for SO=T')
C --- Set up permutation matrix for Bloch ---
      lihdim = nbas * nl**2
      call defi(oidxsh,lihdim)
      if (iprint() .ge. 30) print *
      call makidx(nl,1,1,nbas,0,sspec,ips,-1,w(oidxsh),inxsh)
      ldim = inxsh(0)
      idim = inxsh(1) - ldim
      call rxx(idim .ne. 0,'BNDTB: bad orbital switches')
      ldimx = ldim*nspc
      nbmaxx = nbmax*nspc
      nx = nbmax*nsp1

      if (rl .or. rn) then
        call dpzero(rhoc,nl**4*nbas)
        if (lov) then
          call dpzero(drhosl,3*ldim**2)
        endif
      endif
      if (rn) call dpzero(rhon,nl**4*nsp*nbas*2)

C#ifndef MPI
      if (lgamma) then
        if (lgamsp) then
C --- keep up and down spin eigenvectors in memory (cf evd for metals)
          call defdr(ozup,ldimx*ldimx)
          call defdr(ozdn,ldimx*ldimx)
        else
          call defdr(ozll,ldimx*ldimx)
        endif
      else
        call defdc(ozll,ldimx*ldimx)
      endif
C#endif
C --- Get number of orbitals on each atom ---
      if (evecs .and. .not. donly) then
        call iinit(iwk,nbas)
        icheck = 0
        do  20  ibas = 1, nbas
          do  10  l = 0, nl-1
            indx = idxdn(l,ipc(ibas))
            if (indx .gt. 1) goto 10
            icheck = icheck + 2*l + 1
            iwk(ibas) = iwk(ibas) + 2*l + 1
   10     continue
   20   continue
        call rxx(icheck .ne. ldim,'BNDTB: icheck .ne. ldim')
      endif

C --- Setup for Mulliken DOS ---
      if (stoner .and. mull .ne. 0)
     .  call rx('BNDTB: only MULL=0 allowed for STONER=T')
      if (doswt) then
        mull0 = mod(mull,10)
        mull1 = mod(mull,100) / 10
        imfil = mull / 100
        if (mull1 .gt. 1) imfil = 0
        call rxx(mull0 .ne. 0 .and. mull0 .ne. 1,
     .    'BNDTB: bad ones digit in mull switch')
        nevhi = 0
        nevlo = ldimx
        nat0 = nclass
        if (mull0 .eq. 1) nat0 = nbas

C --- Partial DOS ---
        if (mull1 .eq. 0) then
          nlm = nl
          if (mull0 .eq. 1) nlm = nl*nl
          if (imfil .ne. 0) goto 26
          nchan = 0
          do  25  ibas = 1, nat0
            do  25  lm = 1, nlm
              if (mull0 .eq. 0) then
                l = lm - 1
                ic = ibas
              else
                l = ll(lm)
                ic = ipc(ibas)
              endif
              if (stoner .and. l .eq. 2) then
                if (idxdn(l,ic) .gt. 1) call rx('STONER=T and IDXDN>1')
                index(ibas) = nchan
              endif
              if (idxdn(l,ic) .gt. 1) goto 25
              nchan = nchan + 1
   25     continue
          nchan = nchan*nspc
   26     nlm = nlm*nspc

C --- Bond charge by atom pairs ---
        elseif (mull1 .eq. 1) then
          nlm = 1
          if (imfil .ne. 0) goto 35
          nchan = nat0*(nat0+1)/2
          if (mull0 .eq. 0) then
            indx = 0
            do  30  ic = 1, nclass
              if (nrclas(ic) .gt. 1) indx = indx + 1
   30       continue
            nchan = nchan + indx
          endif
   35     continue

C --- Bond charge by orbital pairs ---
        elseif (mull1 .eq. 2) then
          nlm = nl*nspc
        elseif (mull1 .eq. 3) then
          nlm = nl*nl*nspc
        else
          call rx('BNDTB: bad tens digit in mull switch')
        endif

C --- Determine number of DOS channels from file MULL ---
        if (imfil .ne. 0 .or. mull1 .eq. 2 .or. mull1 .eq. 3) then
          ifi = fopno('MULL')
          rewind ifi
          nchan = 0
   40     read(ifi,*,err=50,end=50) indx
          nchan = nchan + 1
          goto 40
   50     continue
        endif

C --- Allocate DOS arrays ---
        opdos = 1
        opwk  = 1
        obc   = 1
        obcos = 1
        obwk  = 1
        if (mull1 .eq. 0) then
          call defdr(opdos,nlm*nat0)
          if (fsym .and. mull0 .eq. 1) call defdr(opwk,nlm*nat0)
        else
          call defdr(obc,nlm*nlm*nat0*nat0)
          if (fsym .and. mull0 .eq. 1 .and.
     .      (mull1 .eq. 1 .or. mull1 .eq. 2))
     .      call defdr(obwk,nlm*nlm*nat0*nat0)
          if (mull0 .eq. 0) call defdr(obcos,nlm*nlm*nat0)
        endif
        if (iprint() .ge. 20) 
     .    write(*,510) multyp(mull1,mull0),nchan
  510   format(/' BNDTB: Mulliken pop. analysis: ',a45
     .         /'        Number of Mulliken DOS channels = ',i4)

C --- Set up for DOS symmetrization ---
        ogd = 1
        if (fsym .and. (mod(mull,100) .eq. 1 .or. mull1 .eq. 3)) then
          if (nl .gt. 2 .and. mod(mull,100) .eq. 1) then
            call defdr(ogd, -5*5*nsgrp)
            call symtbd(nsgrp,w(osymgr),w(ogd))
          endif
          if (iprint() .ge. 10) then
            if (mull1 .eq. 3) then
              write(*,515)
            elseif (nl .gt. 3) then
              write(*,520)
            else
              write(*,530)
            endif
          endif
        endif
  515   format(/' ***WARNING*** DOS not symmetrized, use full BZ')
  520   format(/' ***WARNING*** DOS for L > 2 not symmetrized')
  530   format(/' DOS weights symmetrized up to L = 2')
      endif

C#ifdefC MPI
CC      if (evd .and. .not. evm .and. evecs) nfilet = fopn('TMP')
C      if (evd .and. .not. evm .and. evecs) nfilet = fopnT('TMP',-1,4,0)
C#else
C --- Open TMP and BAND files; look for fermi energy from last run ---
      nfileb = fopn('BAND')
      if (.not. lgamma .and. evecs .and. metal
     .  .and. .not. evd .and. .not. donly) then
        i = getef(nfileb,11,efermi)
        if (iprint() .gt. 10) write (*,540) efermi
  540   format(/' ***WARNING*** Fermi energy taken from BAND file : ',
     .    f10.6/' **** restart with EVDISC=T in control file **** ')
      elseif (donly) then
        i = getef(nfileb,10,efermi)
      endif
      if (rstrt) then
        call rxx(evecs .and. .not. donly,
     .    'BNDTB: continuation only allowed for e-vals and DOS weights')
        print *
        call iomomq(nfileb,22,nl,1,1,nkp,ldimx,1,ik0,nbmaxx,nchan,nchan,
     .    nev,eband,w,w,w,efermi,dum)
        if (ik0 .lt. 0) call rx('BNDTB: bad moments file')
        ik0 = ik0+1
      else
        rewind nfileb
        j = iomoms(-nfileb,nl,nsp1,nspc,nkp,ldimx,nfstg,1,0,0,0,0,0,0,0,
     .             0d0,0d0,0d0,0d0,0d0,0d0)
        ik0 = 1
      endif
C     if (evd .and. evecs .and. .not. donly) nfilet = fopn('TMP')
      if (evd .and. evecs .and. .not. donly) nfilet=fopnT('TMP',-1,4,0)
C#endif

C --- Determine verbosity for secmtb ---
      ipr=0
      if (iprint() .ge. 30) ipr=1
      if (iprint() .ge. 35) ipr=2
      if (iprint() .gt. 40) ipr=3

C --- Begin loop over spin and k-points ---
      nevmx0 = (int(zval) + 1)/2
      if (nevmx .eq. 0) then
        nevmx = nevmx0
        if (metal) nevmx = min(max(nevmx+nevmx/2,nl*nl),ldim)
        nevmx = nspc*min(nevmx,nbmax)
      endif
      nevmx0 = nspc*min(nevmx0,nbmax)
      if (iprint() .gt. 10 .and. nevmx .lt. nevmx0 .and. leig)
     .  call awrit2(' ***WARNING*** nevmx=%i < %i, '//
     .  'expect errors in integrals',' ',80,i1mach(2),nevmx,nevmx0)
      oxm0 = 1
      osll = 1
      ozwk2 = 1
      owkcf = 1
      thrpv = 0d0
      if (force) call dpzero(f,3*nbas)
      if (trh) call dpzero(esite,nbas*nsp1)
      if (trh .or. charge) call dpzero(rho,6*nbas)
      if (trh .or. charge) call dpzero(rholm,nl**2*2*nbas)
C#ifdefC MPI
C#ifdefC MPI-SCHEDULED
CC MPI: Distribute contiguous k points over processes
CC MPI: for nsp=2 process will do both spins
CC MPI: Find the number >= nkp divisible by numprocs
C      allocate (kpproc(0:numprocs), stat=ierr)
C      call dstrbp(nkp,numprocs,1,kpproc(0))
C      if (evecs) then
C        nkproc = 1
C        if (evd .and. evm) nkproc = kpproc(procid+1) - kpproc(procid)
C        allocate(zll(ldimx,ldimx,2,nkproc), stat=ierr)
C      endif
C      ikproc = 0
C      do  60  ikp = kpproc(procid), kpproc(procid+1)-1
C       if (ikp .eq. kpproc(procid)) then
C         if (mlog) then
C           call gettime(datim)
C           call awrit4(' MPI '//datim//' Process %i of %i on '
C     .       //shortname(procid)(1:namelen(procid))//
C     .       ' starting k-points %i to %i',' ',256,lgunit(3),
C     .       procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
C         endif
C       endif
C#elseC
C      if (evecs) then
C        nkproc = 1
C        if (evd .and. evm) then
CC How many eigenvectors will I do?
C          nkproc = 0
C          do ikp = 1+procid, nkp, numprocs
C            nkproc = nkproc + 1
C          enddo
C        endif
C        if (nkproc .gt. 0) allocate(zll(ldimx,ldimx,2,nkproc),stat=ierr)
C      endif
C      call dpzero(eband,nbmax*nkp)
C      ikproc = 0
C      do  60  ikp = 1+procid, nkp, numprocs
C#endifC
C       if (evd .and. evm) then
C         ikproc = ikproc + 1
C       else
C         ikproc = 1
C       endif
C#else
      ohrs = 1
      oors = 1
      if (cmdopt('--invbl',7,0,outs)) then
        call defdr(ohrs,-nl**4*nsite)
        if (lov) then
          call defdr(oors,-nl**4*nsite)
        endif
      endif
      do  60  ikp = ik0, nkp
C#endif
       if (lgamma .and. ikp .ne. 1) then
         if (iprint() .gt. 5) then
           print *,
     .       'BNDTB: for gamma only or mol set NKABC=1 and BZJOB=0'
         endif
         call fexit(0,119,'  ',0d0)
       endif
       call dpscop(w(oqp),qp,3,3*ikp-2,1,1d0)
       iprx = ipr
       if (ikp .ge. ik0+1) iprx = max(ipr-1,0)
       if (iprint() .ge. 30 .and. mod(ikp,100) .eq. 0) iprx=max(ipr,1)
       do  60  isp = 1, nsp
C --- Make bands and possibly eigenvectors ---
C#ifdefC MPI
C       if (TBU) call rx('BNDTB: secmtb call not yet set')
C       call secmtb(sctrl,plat,nbas,nl,nspc,nsp,ispx,lmx,ipc,w(oidxsh),
C     .   ldimx,nevmx,efmax,ikp,nkp,qp,nsite,iax,npr,hrs,vso,hso,srs,
C     .   pot0,rl,0,leig,nev,zll(1,1,1,ikproc),eband(1,isp,ikp),
C     .   w(ohrs),w(oors))
C#else
       if (nspc .eq. 2) then
C --- coupled spins:
         ispx = 1
       else
         ispx = isp
       endif
       if (lgamsp) then
         if (isp .eq. 1) then
           ozll = ozup
         else
           ozll = ozdn
         endif
       endif
       call secmtb(sctrl,plat,nbas,nl,nspc,nsp,ispx,lmx,ipc,w(oidxsh),
     .   ldimx,nevmx,efmax,ikp,nkp,qp,nsite,iax,npr,hrs,vso,hso,srs,
     .   pot0,rl,iprx,leig,nev,w(ozll),eband(1,isp,ikp),w(ohrs),w(oors))
     .   
C#endif
       if (cmdopt('--invbl',7,0,outs)) goto 60
C --- Make forces, etc. or write e-vec's to disc ---
       if (evecs .and. .not. donly) then
         if (evd) then
C#ifdefC MPI
C           if (.not. evm) call tiozll(lgamma,-nfilet,ldimx,nev,zll)
C#else
           call tiozll(lgamma,-nfilet,ldimx,nev,w(ozll))
C#endif
         elseif (lgamsp) then
C --- Do nothing: up and down eigenvectors stored in memory ---
         else
           if (lgamma) then
             call defdr(ozwk,ldimx*ldimx)
           else
             call defdc(ozwk,ldimx*ldimx)
           endif
           if (trh .or. charge) call defdr(oxm0,nbas*nspc*nl**2)
           if (lov) then
             if (lgamma) then
               call defdr(ozwk2,ldimx*ldimx)
             else
               call defdc(ozwk2,ldimx*ldimx)
             endif
           endif
           if (cryf .or. ocryf) call defdr(owkcf,nbas)
           call defdr(owtkb,-ldimx*nkp*nsp1)
C#ifdefC MPI
C           call tbfrce(0,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp,
C     .       lmx,ipc,nrclas,idxdn,w(oidxsh),iwk,ldimx,nev,zval,ikp,qp,
C     .       abs(dval(w(owtkp),ikp))/nsp1,
C     .       w(owtkb),mpsord,swidth,metal,efermi,nsite,iax,npr,xyzfrz,
C     .       eband(1,isp,ikp),zll,hrs,h0,dh,dhcf,vso,hso,srs,ds,dscf,
C     .       iprx,charge,rl,rn,trh,pv,force,lso,lov,cryf,ocryf,
C     .       w(oxm0),w(ozwk),w(ozwk2),w(owkcf),sumev,entrpy,f,thrpv,
C     .       esite,rho,rholm,rhoc,rhon,drhosl)
C#else
           call tbfrce(0,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp,
     .       lmx,ipc,nrclas,idxdn,w(oidxsh),iwk,ldimx,nev,zval,ikp,qp,
     .       abs(dval(w(owtkp),ikp))/nsp1,
     .       w(owtkb),mpsord,swidth,metal,efermi,nsite,iax,npr,xyzfrz,
     .       eband(1,isp,ikp),w(ozll),hrs,h0,dh,dhcf,vso,hso,srs,ds,
     .       dscf,iprx,charge,rl,rn,trh,pv,force,lso,lov,cryf,ocryf,
     .       w(oxm0),w(ozwk),w(ozwk2),w(owkcf),sumev,entrpy,f,thrpv,
     .       esite,rho,rholm,rhoc,rhon,drhosl)
C#endif
           call rlse(ozwk)
         endif
       endif

C --- Write Mulliken DOS weights to disc ---
C MPI: Can't do this in parallel 'cos you can't all access band file
C#ifndef MPI
       if (doswt) then
         call tcn('DOS weights')
         if (tetra) then
           nevhi = max0(nevhi,nev)
           nevlo = min0(nevlo,nev)
         endif
         if (iprint() .gt. 30 .and. ikp .eq. ik0 .and. isp .eq. 1) then
           print *
           print *,' BNDTB: Putting DOS weights onto disc ...'
         endif
         call defdr(odoswt,-nchan*nev)
         if (lov) then
           call defdc(osll,ldimx*ldimx)
           call defi(owk,nbas)
           call tbloch(.false.,qp,nl,nspc,1,1,nbas,plat,lmx,ipc,
     .       w(oidxsh),nsite,iax,npr,srs,vso,hso,.false.,ldimx,
     .       w(osll),w(owk))
         endif
         call mkwttb(nfileb,nbas,nl,nspc,nclass,ipc,nrclas,idxdn,
     .     w(oidxsh),nev,ldimx,mod(mull,100),imfil,nchan,nlm,nat0,ifi,
     .     nsgrp,w(oistab),w(osymgr),w(ogd),lov,w(ozll),w(osll),
     .     w(opdos),w(obc),w(obcos),w(opwk),w(obwk),eband(1,isp,ikp),
     .     w(odoswt))
         call rlse(odoswt)
         call tcx('DOS weights')
       else
         write (nfileb) 0, ldimx
         call dpdump(eband(1,isp,ikp),ldimx,-nfileb)
       endif
C#endif
   60 continue
      if (cmdopt('--invbl',7,0,outs)) then
        call awrit0('Real space hamiltonian:',' ',128,i1mach(2))
        call defdr(oalpha,-nbas*nl**2)
        call defi(ontab,nbas+1)
        call npr2tb(0,nbas,npr,w(ontab))
        if (iprint() .ge. 60) then
          call shostr(nl**2,nsite,nbas,plat,w(opos),0,w(oalpha),
     .                iax,w(ontab),w(ohrs),1,1,0d0,0,1d0)
          if (lov) then
            call awrit0('Real space overlap:',' ',128,i1mach(2))
            call shostr(nl**2,nsite,nbas,plat,w(opos),0,w(oalpha),
     .                 iax,w(ontab),w(oors),1,1,0d0,0,1d0)
          endif
        endif
        call wrirsh(.false.,.true.,lov,nl,nbas,plat,w(opos),alat,
     .              w(oclabl),ipc,nsite,iax,w(ohrs),w(oors),h0)
        call rlse(ohrs)
        call fexit(0,9,'Real space hamiltonian written to disc',0)
      endif
C#ifdefC MPI
CC MPI: Wait here till all the bands are made
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPI-SCHEDULED
CC MPI: Broadcast the bands at your k points
C      do  pid = 0, numprocs-1
C        ikp = kpproc(pid)
C        jkp = kpproc(pid+1) - ikp
C        if (procid .eq. pid) then
C          if (mlog) then
C            call gettime(datim)
C            call awrit5(' MPI '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' bcast k-pts %i to %i, %i evls',' ',256,lgunit(3),
C     .        procid,numprocs,kpproc(procid),kpproc(procid+1)-1,
C     .        jkp*nx)
C          endif
C        endif
C        call MPI_BCAST(eband(1,1,ikp),jkp*nx,MPI_DOUBLE_PRECISION,pid,
C     .                 MPI_COMM_WORLD,ierr)
C      enddo
C#elseC
CC MPI: Pool the bands by ALLREDUCE
C      allocate(buffer(nx*nkp), stat=ierr)
C      call MPI_ALLREDUCE(eband,buffer,nx*nkp,MPI_DOUBLE_PRECISION,
C     .                   MPI_SUM,MPI_COMM_WORLD,ierr)
C      call dcopy(nx*nkp,buffer,1,eband,1)
C      deallocate(buffer, stat=ierr)
C#endifC
C      if (procid .eq. master) then
C        do  ikp = 1, nkp
C        do  isp = 1, nsp
C          call awrit2(' MPI Bands at kpt %i of %i',' ',120,lgunit(1),
C     .                ikp, nkp)
C          write(lgunit(1),'(255(5x,9f8.4:/))')
C     .         (eband(i,isp,ikp), i=1,ldim)
C        enddo
C        enddo
C      endif
CC MPI: Merge f,esite,rho (if !evd)
C      if (evecs .and. .not. evd) then
C        if (force) then
C          allocate(buffer(3*nbas), stat=ierr)
C          call MPI_ALLREDUCE(f,buffer,3*nbas,MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(3*nbas,buffer,1,f,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (trh) then
C          allocate(buffer(nbas*nsp), stat=ierr)
C          call MPI_ALLREDUCE(esite,buffer,nbas*nsp,
C     .                       MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(nbas*nsp,buffer,1,esite,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (charge) then
C          allocate(buffer(nl*nsp*nbas), stat=ierr)
C          call MPI_ALLREDUCE(rho,buffer,nl*nsp*nbas,
C     .                       MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(nl*nsp*nbas,buffer,1,rho,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (rl) then
C          allocate(buffer(nl**4*nsp*nbas), stat=ierr)
C          call MPI_ALLREDUCE(rhol,buffer,nl**4*nsp*nbas,
C     .                       MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(nl**4*nsp*nbas,buffer,1,rhol,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (pv) then
C          call MPI_ALLREDUCE(thrpv,thrpv_buffer,1,MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          thrpv = thrpv_buffer
C        endif
C      endif
C#endif

C --- Check for partial bands in tetrahedron DOS ---
C#ifndef MPI
      if (doswt) then
        if (mull1 .eq. 0) then
          call rlse(opdos)
        else
          call rlse(obc)
        endif
        if (tetra) then
          call rxx(nevhi .lt. nevlo, 'BNDTB: programming error')
          if (nevhi .gt. nevlo .and. iprint() .ge. 10) then
            nevlo = nevlo+1
            evlo = eband(nevlo,1,1)
            do   ikp = 1, nkp
              do   isp = 1, nsp
                evlo = dmin1(evlo,eband(nevlo,isp,ikp))
              enddo
            enddo
            write(*,550) evlo
  550       format(/' BNDTB: Mulliken DOS only reliable below ',
     .        'E = ',f9.5)
          endif
        endif
      endif
C#endif

C --- BZ integration for fermi level, band sum and qp weights ---
C      if (.not. donly .and. (.not. lgamma .or. lgamsp)) then
      if (.not. donly) then
        call defdr(owtkb,ldim*nsp1*nkp)
        norder = mpsord
C --- enforce metal sampling even if filled bands encountered ---
C        norder = sign(1,mpsord)*(abs(mpsord) + 100)
C        call rxx(norder .lt. 0,
C     .    'F-D not implemented for spin pol. Restart with N>=0')
        call bzwts(nbmax,ldim,nsp1,nspc,n1,n2,n3,nkp,ntet,w(oidtet),
     .    zval,metal,tetra,norder,ndos,swidth,srnge,w(owtkp),eband,
     .    efmax,efermi,sumev,w(owtkb),dosef,qval,entrpy)
      endif
C#ifndef MPI
      write (nfileb) efermi, 0d0
      call fclose (nfileb)
      if (donly) goto 1020
C#endif

C --- Use e-vec's from disc to make forces, etc. ---
      if ((evd .and. evecs) .or. lgamsp) then
        if (.not. evm) rewind nfilet
C#ifdefC MPI
C        ikproc = 0
C#ifdefC MPI-SCHEDULED
C        do  130  ikp = kpproc(procid)+1, kpproc(procid+1)
C#elseC
C        do  130  ikp = 1+procid, nkp, numprocs
C#endifC
C          if (evd .and. evm) then
C            ikproc = ikproc + 1
C          else
C            ikproc = 1
C          endif
C#else
         do  ikp = 1, nkp
C#endif
           iprx = ipr
           if (ikp .ge. 2) iprx = max(ipr-1,0)
           call dpscop(w(oqp),qp,3,3*ikp-2,1,1d0)
           if (cmdopt('--LUMO',6,0,outs)) then
             call lumo(lmol,qval,ldim,nsp,nspc,nbmax,eband,w(owtkb))
           endif
           do  isp = 1, nsp
C#ifdefC MPI
C          if (.not. evm) call tiozll(lgamma,nfilet,ldimx,nev,zll)
C#else
           if (evd) then
             call tiozll(lgamma,nfilet,ldimx,nev,w(ozll))
           elseif (lgamsp) then
             if (isp .eq. 1) then
               ozll = ozup
             else
               ozll = ozdn
             endif
           else
             call rx('BNDTB error: what are you doing here?')
           endif
C#endif
           if (lgamma) then
             call defdr(ozwk,ldimx*ldimx)
           else
             call defdc(ozwk,ldimx*ldimx)
           endif 
           if (trh .or. charge) call defdr(oxm0,nbas*nspc*nl**2)
           if (lov) then
             if (lgamma) then
               call defdr(ozwk2,ldimx*ldimx)
             else
               call defdc(ozwk2,ldimx*ldimx)
             endif
           endif
           if (cryf .or. ocryf) call defdr(owkcf,nbas)
C#ifdefC MPI
C          call tbfrce(1,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp,lmx,
C     .      ipc,nrclas,idxdn,w(oidxsh),iwk,ldimx,nev,zval,ikp,qp,
C     .      abs(dval(w(owtkp),ikp))/nsp,
C     .      w(owtkb),mpsord,swidth,metal,efermi,nsite,iax,npr,xyzfrz,
C     .      eband(1,isp,ikp),zll(1,1,1,ikproc),hrs,h0,dh,dhcf,vso,hso,
C     .      srs,ds,dscf,iprx,charge,rl,rn,trh,pv,force,lso,lov,cryf,
C     .      ocryf,w(oxm0),w(ozwk),w(ozwk2),w(owkcf),sumev,entrpy,f,
C     .      thrpv,esite,rho,rholm,rhoc,rhon,drhosl)
C#else
           call tbfrce(1,lmol,lgamma,plat,nbas,nl,nspc,nsp,nsp1,isp,lmx,
     .       ipc,nrclas,idxdn,w(oidxsh),iwk,ldimx,nev,zval,ikp,qp,
     .       abs(dval(w(owtkp),ikp))/nsp1,
     .       w(owtkb),mpsord,swidth,metal,efermi,nsite,iax,npr,xyzfrz,
     .       eband(1,isp,ikp),w(ozll),hrs,h0,dh,dhcf,vso,hso,srs,ds,
     .       dscf,iprx,charge,rl,rn,trh,pv,force,lso,lov,cryf,ocryf,
     .       w(oxm0),w(ozwk),w(ozwk2),w(owkcf),sumev,entrpy,f,thrpv,
     .       esite,rho,rholm,rhoc,rhon,drhosl)
C#endif
           call rlse(ozwk)
         enddo
       enddo
C MPI: Merge f,esite,rho,rhol etc
        if (.not. evm) call fclose(nfilet)
      endif
C#ifdefC MPI
C      deallocate(kpproc, stat=ierr)
C      if (evecs) deallocate(zll, stat=ierr)
CC MPI: Merge f,esite,rho,rhol etc
C      if (evecs .and. evd) then
C        if (force) then
C          allocate(buffer(3*nbas), stat=ierr)
C          call MPI_ALLREDUCE(f,buffer,3*nbas,MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(3*nbas,buffer,1,f,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (trh) then
C          allocate(buffer(nbas*nsp), stat=ierr)
C          call MPI_ALLREDUCE(esite,buffer,nbas*nsp,
C     .                       MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(nbas*nsp,buffer,1,esite,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (charge) then
C          allocate(buffer(nl*nsp*nbas), stat=ierr)
C          call MPI_ALLREDUCE(rho,buffer,nl*nsp*nbas,
C     .                       MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(nl*nsp*nbas,buffer,1,rho,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (rl) then
C          allocate(buffer(nl**4*nsp*nbas), stat=ierr)
C          call MPI_ALLREDUCE(rhol,buffer,nl**4*nsp*nbas,
C     .                       MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          call dcopy(nl**4*nsp*nbas,buffer,1,rhol,1)
C          deallocate(buffer, stat=ierr)
C        endif
C        if (pv) then
C          call MPI_ALLREDUCE(thrpv,thrpv_buffer,1,MPI_DOUBLE_PRECISION,
C     .                       MPI_SUM,MPI_COMM_WORLD,ierr)
C          thrpv = thrpv_buffer
C        endif
C      endif
C#endif
      call rlse(oidxsh)

C --- Symmetrize forces ---
      if (fsym .and. force) then
        call defdr(ofwk,3*nbas)
        call symfor(nbas,1,w(osymgr),nsgrp,w(oistab),w(ofwk),f)
        call rlse(ofwk)
      endif

C --- Symmetrise rho, rholm and esite ---
      if (fsym .and. (charge .or. trh)) then
        call defdr(ofwk,nl*nbas*2)
        owklm = 1
        if (lov) call defdr(owklm,nl**2*nbas*2)
        owk = 1
        if (trh) call defdr(owk,nbas*2)
        call symre(lov,trh,nl,nsp1,nbas,nsgrp,w(oistab),w(ofwk),
     .             w(owklm),w(owk),rho,rholm,esite)
        call rlse(ofwk)
      endif

C --- Symmetrise rhoc, rhon ---
      if (fsym .and. (rl .or. rn)) then
        call defdr(ofwk,nl**4*nbas)
        call defdr(ogd2, -5*5*nsgrp)
        call symtbd(nsgrp,w(osymgr),w(ogd2))
        if (rn) then
          do  isp = 1, nsp
            call symrtb(nl,nsp,nbas,nsgrp,w(osymgr),w(ogd2),w(oistab),
     .                  w(ofwk),rhon(1,1,1,isp))
          enddo
        endif
        if (rl) then
          call symrtb(nl,1,nbas,nsgrp,w(osymgr),w(ogd2),w(oistab),
     .                w(ofwk),rhoc)
        endif
        call rlse(ofwk)
      endif

C#ifndef MPI
C --- Generate total DOS or BZ map on disk if requested ---
 1020 continue
      if (bittst(ldos,1)) then
        if (bzmp) then
C --- BZ map ---
          call defdr(odos,nqmx*nwmx*npln)
          if (iprint() .ge. 30) print 610
  610     format(/' ...Generating BZ map from bands...')
          call makbzm(nkp,ldimx,nbmaxx,nsp1,eband,metal,mpsord,swidth,
     .      npln,nwmx,nqmx,nw,ew,nqpl,wtbzm,w(odos))
          call bzmio(w(odos),1,npln,nwmx,nqmx,nw,ew,nqpl,vxpl,vypl,
     .      xpl,ypl,zpl,nkp,nsp1,metal,mpsord,swidth,wtbzm,efermi,
     .      -fopn('BZMP'))
        else
C --- Total DOS ---
          call tcn('total DOS')
          if (tetra) then
            call defdr(odos,(nsp1+1)*ndos)
            call bzints(n1,n2,n3,eband,dum,nkp,ldimx,nbmaxx,nsp1,
     .        dosw(1),dosw(2),w(odos),ndos,efermi,1,ntet,w(oidtet),
     .        dum,dum)
            if (swnos) then
              if (iprint() .ge. 30) print 620
  620         format(/' ... Generating integrated DOS ...')
            else
              call xxxdif(dosw(1),dosw(2),ndos,nsp1,0,w(odos))
              if (iprint() .ge. 30) print 630
  630         format(/' ... Generating total DOS ...')
            endif
            del = 0d0
          else
            call defdr(odos,nsp1*ndos)
            if (swnos) then
              call maknos(nkp,ldimx,nbmaxx,nsp1,w(owtkp),eband,mpsord,
     .          swidth,-srnge,dosw(1),dosw(2),ndos,w(odos))
              if (iprint() .ge. 30) print 620
            else
              call makdos(nkp,ldimx,nbmaxx,nsp1,w(owtkp),eband,mpsord,
     .          swidth,-srnge,dosw(1),dosw(2),ndos,w(odos))
              if (iprint() .ge. 30) print 630
            endif
            del = mpsord+swidth
          endif
          if (nspc .eq. 2) call dscal(ndos,.5d0,w(odos),1)
          call dosio(w(odos),ndos,nsp1,ndos,1,dosw(1),dosw(2),nsp1,
     .      efermi,del,1,-fopn('DOS'))
          call fclose(fopn('DOS'))
          call tcx('total DOS')
        endif
        call rlse(odos)
      endif

C --- Make local DOS for Stoner magnetism ---
      if (stoner .and. .not. donly) then
        call rxx (nl .lt. 3,'BNDTB: nl must be > 2 for Stoner')
        if (iprint() .ge. 30) call awrit3(
     .   ' BNDTB: Make number-of-states functions for Stoner model'//
     .   ', %i bins in (%d,%d)',' ',100,i1mach(2),ndos-1,dosw(1),
     .   dosw(2))
        nfileb = fopn('BAND')
        call defdr(owk,ndos)
        if (tetra) then
          call defdr(oeband,ldimx*nkp)
          call defdr(odoswt,-nchan*ldimx*nkp)
          nfstg = 11
          call iomomq(nfileb,22,nl,1,1,nkp,ldimx,nfstg,ik0,ldimx,nchan,
     .      nchan,nev,w(oeband),w,w(odoswt),w,dum,dum)
          call rxx(ik0 .ne. nkp,'BNDTB:  moments file missing qpts')
          call dostet(ldimx,1,1,nev,nchan,n1,n2,n3,ntet,w(oidtet),
     .      w(oeband),w(odoswt),ndos,dosw(1),dosw(2),
     .      .true.,w(owk),zos)
        else
          call defdr(oeband,ldimx)
          call defdr(odoswt,nchan*ldimx)
          nfstg = 0
          call iomomq(nfileb,12,nl,1,1,nkp,ldimx,nfstg,ik0,ldimx,nchan,
     .      nchan,nev,w,w,w,w,dum,dum)
          call rxx(ik0 .ne. nkp,'BNDTB:  moments file missing qpts')
          call dosspl(nfileb,ldimx,1,1,nchan,mpsord,swidth,nkp,w(owtkp),
     .      w(oeband),w(odoswt),ndos,dosw(1),dosw(2),.true.,
     .      w(owk),zos)
        endif
        call rlse(owk)
C --- scale zos to get per spin ---
        call dscal(ndos*nchan,0.5d0,zos,1)
        call fclose(nfileb)
      endif
C#endif

      end
      
