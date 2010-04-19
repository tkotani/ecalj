C#define BLAS3
      subroutine secmtb(sctrl,plat,nbas,nl,nsp,nspu,ispu,lmx,ipc,indxsh,
     .  ldim,nevmx,efmax,ikp,nkp,bk,nsite,iax,npr,hrs,vso,hso,srs,pot0,
     .  rl,ipr,leig,nev,z,eb,rhrs,rors)
C- Set up tight-binding Hamiltonian; diagonalize secular matrix
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat,nbas,nl,nsp,lmx,ipc,indxsh
Ci   nsp =  2 for coupled spins (empirical spin-orbit, lso=T)
Ci   nspu = 2 for TB+U, ispu is current spin
Ci   ldim: dimension of l - wave block of hamiltonian matrix
Ci   nevmx, max no. of e-vec's; efmax, highest eigenvalue to seek
Ci   ikp,nkp,bk: current k-point, total number of k-points, k-points
Ci   nsite, total number of neighbors in all clusters;
Ci   iax, neighbor lists; npr, see tbham;
Ci   hrs,srs: real-space hamiltonian and overlap matrices
Ci   vso,hso: table of spin-orbit parameters and the hamiltonian
Ci   ipr:  (verbosity) 0, print nothing; 1, print 1st 9 e'vals and
Ci         timings; 2, print all e'vals and timings
Ci   leig: true: calculate eigenvectors; false: do not
Ci   pot0: monopole potentials, from tbesel
Co Outputs:
Co   nev: number of eigenvectors found from diagno
Co   eigenvalues and eigenvectors are returned in eb, z (leig true)
Co   (with command line --invbl)
Co   rhrs, rsrs: H and O in real space formatted like strux.
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,nsp,nspu,ispu,ldim,nevmx,ikp,nkp,nsite,ipr,nev
      integer lmx(0:*),ipc(nbas),indxsh(nl**2*nbas),iax(1),npr(1)
      double precision sctrl(1),efmax,plat(3,3),bk(3),hrs(1),vso(1),
     .  hso(1),srs(1),z(ldim,ldim*2),eb(1),pot0(nbas),rhrs(1),rors(1),
     .  wtkp(nkp)
      logical leig,rl
C heap:
      integer w(1)
      common / w /   w
C local variables
      integer owk,odiawk,ohk,osk,bitand
      integer i,j,iprint,ii,lgunit,l2,nevec,ipr0,i1mach,linv,ltb,
     .  lncol,lov
      logical lgors,bittst,lx,lso,lgamma,invb,cmdopt
      character*80 outs

      invb = cmdopt('--invbl',7,0,outs)

      call tcn('secmtb')

C --- Set up dimensions and allocate arrays ---
      call upack('ctrl ltb lncol',sctrl,ltb,lncol,0,0,0)
      lov = bitand(ltb,1)
      lso = bittst(lncol,4)
      lgamma = bittst(ltb,2**17) .or. bittst(ltb,2**18)
      call rxx(invb.and.lgamma,'SECMTB: --invbl needs GAMMA=F')
      nevec = 0
      if (leig) then
        nevec = nevmx
        if (nevmx .eq. 0) nevec = ldim
      endif
      osk = 1
      l2 = ldim**2
      call defi(owk,nbas)
      if (lgamma) then
        call defdr(ohk,l2)
        if (lov .ne. 0) call defdr(osk,l2)
      else
        call defdc(ohk,l2)
        call dcopy(2*l2,0d0,0,w(ohk),1)
        if (lov .ne. 0) call defdc(osk,l2)
      endif

C --- Bloch-transform Hamiltonian and overlap ---
      call tcn('assemble H')
      call tbloch(lgamma,bk,nl,nsp,nspu,ispu,nbas,plat,lmx,ipc,indxsh,
     .  nsite,iax,npr,hrs,vso,hso,lso,ldim,w(ohk),w(owk))
      if (lov .ne. 0) then
C --- ispu = 1 here 'cos S is spin-independent
        call tbloch(lgamma,bk,nl,nsp,1,1,nbas,plat,lmx,ipc,indxsh,
     .    nsite,iax,npr,srs,vso,hso,.false.,ldim,w(osk),w(owk))
      endif
      call tcx('assemble H')

C --- add off-site self consistent shifts ---
      if (lov .ne. 0 .and. rl) then
        call addoff(nl,nbas,lmx,ipc,indxsh,ldim,pot0,w(osk),w(ohk),
     .              w(owk))
      endif

C --- Printout ---
      if (iprint() .gt. 100) then
C        call defi(oicl, nbas*nsp)
C        call icopy(nbas,ipc,1,w(oicl),1)
C        if (nsp .eq. 2) call icopy(nbas,ipc,1,w(oicl+nbas),1)
C        call prmsk2(lmx,lmx,w(oicl),2,ldim,ldim,w(ohk),'?','Ham',1,
C     .    ldim,indxsh)
C        if (lov .ne. 0) call prmsk2(lmx,lmx,w(oicl),2,ldim,ldim,w(osk),
C     .    '?','Overlap',1,ldim,indxsh)
        if (lgamma) then
          call yprm('ham',1,w(ohk),ldim*ldim,ldim,ldim,ldim)
          if (lov .ne. 0)
     .    call yprm('ovlp',1,w(osk),ldim*ldim,ldim,ldim,ldim)
        else
          call yprm('ham',12,w(ohk),ldim*ldim,ldim,ldim,ldim)
          if (lov .ne. 0)
     .    call yprm('ovlp',12,w(osk),ldim*ldim,ldim,ldim,ldim)
        endif
      endif

      if (invb) then
        call invbl(.false.,lov,nbas,nl,nsp,ldim,plat,nsite,iax,indxsh,
     .             ikp,nkp,bk,w(ohk),w(osk),rhrs,rors)
        goto 999
      endif

C --- Diagonalize Hamiltonian ---
C#ifdef BLAS3
      lx = .true.
C#elseC
C      lx = .false.
C#endif
      linv = 0
      if (nevec .gt. 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1
      if (linv .ne. 0) then
        call defdr(odiawk,ldim*11)
      else
C ...   Use 5*ldim for parallel implementations ...
        call defdr(odiawk,ldim*5)
      endif
      ipr0 = 0
      if (ipr .ge. 1 .and. ikp .eq. 1) ipr0 = ipr
      if (lgamma) then
        call dsev1 (ldim,w(ohk),w(osk),w(odiawk),ipr0,lx,lov.ne.0,linv,
     .    nevec,efmax,nev,z,eb)
      else
        call diagno(ldim,w(ohk),w(osk),w(odiawk),lx,lov,linv,nevec,
     .    efmax,nev,z,eb)
      endif
      if (iprint() .gt. 100) then
        if (lgamma) then
          call yprm('z',1,z,ldim*ldim,ldim,ldim,ldim)
        else
          call yprm('z',12,z,ldim*ldim,ldim,ldim,ldim)
        endif
      endif

C --- Printout ---
      if (ipr .ge. 1) then
        if (lgamma) then
          print *,' SECMTB: hamiltonian real '
        endif
        j = min(9*nsp,ldim)
        if (ipr .ge. 2) j = ldim
C#ifdefC LINUX_PGI
C        do  18  ii = 1, 1
C#else
        do  18  ii = 1, 2
C#endif
        write(lgunit(ii),'()')
        call awrit3(' SECMTB:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,bk)
   18   write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
        if (ipr .ge. 2) call awrit5
     .    (' nev, nevmx, nevec, ldim=  %i  %i  %i  %i  efmax= %1;5d',
     .    ' ',80,i1mach(2),nev,nevmx,nevec,ldim,efmax)
        call ftflsh(lgunit(1))
      endif
      call rlse(owk)

C --- exit ---
  999 continue 
      call tcx('secmtb')
      end
