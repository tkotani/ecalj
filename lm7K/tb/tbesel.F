C#define MPI-SCHEDULED
      subroutine tbesel(mixQ,ssite,ltb,nbas,nl,nsp,nclass,ipc,dclabl,
     .                  idxdn,z,symgr,nsgrp,istab,bas,awld,alat,vol,
     .                  dlat,nkd,glat,nkg,qpol,indxcg,jcg,cg,cy,gaunt,
     .                  ak,rho,mmom,drhosl,ldim,drhos,rhoc,rhon,qnu,
     .                  stni,idu,uh,jh,qmpol,vm,vu,vj,dh,pot0,ecorr,fwk,
     .                  f,fnou,fnom,ppdip,tpvq)
C- Make L-expanded electrostatic potential, multipole moments; force
C ----------------------------------------------------------------------
Ci Inputs: rhoc=c*_RL c_RL' + cc; 
Ci         qpol, 10 polarisation parameters;
Ci         rhon=c*_RL S c_RL' + cc (TB+U only);
Ci         drhosl=c*_RL dS_RLR'L'/dR c_R'L' + cc (ovlp only)
Ci         drhos: work array to make drhosl summed over all LL'
Ci         gaunt: Gaunt coefficients
Ci    stni: Stoner parameter (l=2)
Ci    idu: for TB+U, index to which channels get additional potential
Ci    uh, uj: Hubbard U's and J's for TB+U and spin pol TB-L
Ci            uh is the screened U for TB+U. Unscreened U is taken from
Ci            start parameters (qnu)
Ci    rho: s, p, d Mulliken charges, dimension nl,nsp,nbas (UL only)
Ci    qmpol: multipole moments on each site (from tbmpole)
Ci    mmom: magnetic moments on each site (from tbmpole)
Co Outputs:
Co    rho: s, p, d Mulliken charges, dimension nl,nsp,nbas
Co    vm: Madelung potential due to charge transfer
Co    vu: Hartree U part of the potential
Co    vj: Stoner J part of the potential
Co    dh: increment to on-site hamiltonian, \Delta H_{RLRL'},
Co    pot0: monopole potential at site R (eq. 7.81, Finnis)
Co        ---all these in Rydberg atomic units.
Co    f: force from e'static terms
Co    fnou: force from overlap dependence of monopoles (Hubbard part)
Co    fnom: force from overlap dependence of monopoles (Madelung part)
Co    ecorr: second order correction to total energy
Co           E_2 of Finnis' book, eq (7.75); or E^U in TB+U
Co    ppdip: accumulated dipole moment, RELAX=0 in ctrl is misused 
Co           to permit the dipole to be accumulated over a subset
Co           of the sites (see tbtote also)
Co    tpvq: contribution of multipoles to 3pV (pressure)
Cr Remarks
Cr    The multipoles Q_L and components of electrostatic potential V_L
Cr    are defined here in terms of the conventions of Stone:
Cr    "Theory of intermolecular forces." If the multipole moments of the
Cr    charge at sites i are Q_L then the components of the potential
Cr    at sites j are defined from 
Cr    V(r) = \Sum_L sqrt(4pi/(2l+1)) V_L r^l Y_L(r)
Cr         = \Sum_L sqrt(4pi/(2l+1))(2l+1)!! V_L J_L(r) where J_L are
Cr    Bessel functions as defined by Michael (see subroutine BESSL).
Cr    V_L = \Sum_L' sqrt((2l+1)*(2l'+1))/((2l+1)!!(2l'+1)!!) B_L'L Q_L'
Cr    where B_L'L are Michael's structure constants. Here the Q are
Cr    Stone's, not Jackson's multipole moments. They are related by
Cr    Q_L=sqrt(4pi/(2l+1)) Q_L(Jackson)
Cr    Jackson's are the ones used in the FP-LMTO and NFP programs.
Cr    That means that the conventional multipole moments (Jackson) are
Cr    these moments multiplied by sqrt{(2l+1)/4pi}. The V_L here are
Cr    defined such that the potential V(r) is expanded as
Cr    V(r)=\sum_L \sqrt(4pi/(2l+1)) V_L r^l Y_L(r), 
Cr    therefore they must be multiplied by \sqrt{4pi/(2l+1)} 
Cr    to get the V_L used in FP-LMTO, in which
Cr    V(r)=\sum_L V_L r^l Y_L(r).
Cr    The extra force from the electrostatics is \sum_L Q_L grad V_L
Cr    Note that in the Finnis 2nd order theory, multipoles are those
Cr    of the difference in charge with respect to free atoms. The old
Cr    MRS theory referred to total electronic and core charge is no
Cr    longer implemented. These theories are essentially equivalent
Cr    anyway.
Cr
Cr    If NOUAVG=F is set in CTRL, then U is averaged over all the
Cr    open channels on that site. Hence the Hubbard potential is
Cr    UdN where dN is the total charge transfer on that site. Otherwise
Cr    the Hubbard potential is V_l = \sum_l U_{ll'} dn_l' and dn_l is
Cr    the charge tansfer in the l-channel. This latter approach is
Cr    closer to the TB+U formula for E^U. In that case non diagonal
Cr    U_{ll'} must be defined. For now, we use U_{ll'}=min(U_l, U_l').
Cr    Other possibilities would be (U_l+U_l')/2 or \sqrt(U_l*U_l')
Cr    The Hubbard energy is then (1/2)\sum_{ll'}U_{ll'}dn_l dn_l'
Cr
Cr    In TB+U we determine on-site Coulomb integrals indirectly using
Cr    I=(U+2lJ)/(2l+1) in the d channels (l=2) and U is taken from
Cr    qnu and I is taken from the I= token in SPEC. In the s and p
Cr    channels we take U from qnu and set J=0. We then use F^0=U
Cr    and J=(5I-U)/4=(F^2+F^4)/14 if l=2. This furnishes us with the
Cr    Slater parameters F^0=U (all l); F^2=14J/2.6, F^4=1.6F^2 if l=2
Cr    F^2=F^4=0 if l<2.
Cu Upgrades
Cu    As of version 8.0, tbe uses (\delta q)*U rather than q*U as the
Cu    'Hubbard' potential. \delta q = q - q_0 and q_0 is computed in
Cu    each channel from qnu. The Hubbard U may now be the average 
Cu    over all channels. The older theory (as in the MRS paper) is
Cu    no longer implemented. 
Cu    TB+U retains all terms except the Hubbard potential which is
Cu    modified to be spin (s) and orbital dependent:
Cu    V^{s}_{LL'}=\sum_{L''L'''} V_{LL''L'''L'} dn^{-s}_{L''L'''} +
Cu                    (V_{LL''L'''L'}-V_{LL''L'L'''}) dn^{s}_{L''L'''}
Cu    n is the spin polarised density matrix (see tbfrce) and dn is
Cu    the difference given by subtracting the diagonal elements q_0,
Cu    the V are given in terms of the Slater radial integrals as 
Cu    V_{LL'L''L'''}=\sum_{k=0}^{2lmax} R^k(ll'l''l''') A^k(LL'L''L''')
Cu    (k an even number less that 2l) with
Cu    A^k(LL'L''L''')=4pi/2k+1 \sum_{p=-k}^{k} C_{LL'''K} C_{L'L''K}
Cu    where K={kp} as L={lm} and C are the Gaunt (CG) coefficients.
Cu    As a first approximation to the Slater integrals we will have
Cu    R^0(llll) = F^0 = U, in all channels
Cu    R^2(2222) = F^2, R^4(2222) = F^4 (ie l=2, d-channel)
Cu    all others zero.
Cu Updates
Cu   04 Jun 08 (ATP) Cleanup, multipole contribution to pressure
Cu   17 Mar 06 (ATP) Cleanup, strip out old MRS theory, start on 
Cu                   non orthogonal TB-L
Cu    6 Jun 05 (ATP) Added TB+U
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   03 Dec 01 (ATP) bug fix
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
C#ifdefC MPI
C      implicit none
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
C      integer lgunit
C#endif
      logical mixQ
      integer ltb,nbas,nl,nsp,nclass,nsgrp,istab(1),nkd,nkg,
     .        idxdn(0:nl-1,1),ldim
      integer ipc(1),indxcg(1),jcg(1),idu(4,nbas)
      double precision ecorr
      double precision ssite(1),bas(3,1),z(1),symgr(1),awld,vol,dlat(1),
     .  glat(1),qpol(10,nclass),rho(nl,2,nbas),
     .  rhoc(nl**2,nl**2,nbas),alat,mmom(nbas),
     .  qnu(3,0:nl-1,nsp,nclass),stni(nclass),uh(4,nbas),jh(4,nbas),
     .  dclabl(nclass),qmpol(9,nbas),cy(1),cg(1),gaunt(9,9,25),
     .  vm(25,nbas),vu(0:nl-1,nbas),vj(2,nbas),dh(nl**2,nl**2,nbas,nsp),
     .  pot0(nbas),rhon(nl**2,nl**2,nbas,2),tpvq,
     .  No(3,nbas,nsp),f(3,nbas),fnou(3,nbas),fnom(3,nbas),fwk(3,nbas),
     .  Ak(9,9,9,9,3),ppdip(3),drhosl(ldim,ldim,3),drhos(nbas,nbas,3)
C Local variables
C#ifdefC MPI
C      integer pid,nbproc,ibproc
C      integer, dimension(:), allocatable :: bproc
C      double precision, dimension(:), allocatable :: buffer
C#endif
      logical pv,lov,MOL,bittst,UL,Uav,TBU,diagn,point,field,cmdopt
      integer i,k,ic,jc,ib,jb,l,lp,ilm,ilmp,ilmpp,ip,i1mach,iprint,
     .        le,lmax,lmxf,lmxst,nlx,nlf,nlm,ilm1,ilm2,ilm3,ilm4,
     .        isp,l1,l2,l3,l4,ifrlx(3),opi,opj,nlmi,nlmj,
     .        ilmi,ilmj,it(3)
      integer ll,getnlm,parg
      double precision M,dq(0:2),v(25),dhd(25),qmp(9),
     .                 sumV,sumU,sumJ,dV1(2:4),dV2(5:9),tau(3),
     .                 hl(49),bl(49),strx(25,9),drstrx(25,9),U,J,F2,F4,
     .                 dQtot,drho1,drho2,drho3,q0,
     .                 VL(9,9,9,9),Rk(9,9,9,9,3),Efield(3)
      double precision pi,dsum,dsqrt,ddot,getavU,getavJ,d1mach,dmin1
C     double precision dasum,efg(5)
      character*8 clabl
      character*120 strn
C#ifdefC F90 | AUTO-ARRAY
C      double precision Ratm(nbas,3),A(nbas),B(nbas)
C#else
      integer nbmx
      parameter(nbmx=2000)
      double precision Ratm(nbmx,3),A(nbas),B(nbas)
      if (nbas .gt. nbmx) call rx(' TBESEL: nbas > nbmx')
C#endif

      call tcn('tbesel')
C#ifdefC MPI
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      mlog = cmdopt('--mlog',6,0,strn)
C#endif
      pi = 4d0*datan(1d0)
      call dpzero(No,3*nbas*nsp)
      call dpzero(ppdip,3)
      call dpzero(A,nbas)
      lov = bittst(ltb,1)
      call dpzero(Efield,3)
      if (iprint() .ge. 30) then
        print *
        print *, 'TBESEL (Stone''s definitions for Q and V)'
      endif
      field = .false.
C ... get external electric field from command line (could put in ctrl)
      if (cmdopt('-efield=',8,0,strn)) then
        ip = 8
        call skipbl(strn,len(strn),ip)
        i = parg(' ',4,strn,ip,len(strn),', ',2,3,it,efield)
        if (i .lt. 0) then
          call rxs2('TBESEL: failed to parse "',strn(1:ip+5),' ..."')
        else
          field = .true.
          if (iprint() .gt. 20) then
            call awrit1(' TBESEL external electric field E= %3:-2,4;4d',
     .        ' ',120,i1mach(2),efield)
          endif
        endif
      endif
C calculate pressure terms from e-statix
      pv = bittst(ltb,2**7)
C ... Molecule or cluster
      MOL =  bittst(ltb,2**18)
C ... 2nd order self consistent TB (TB-L)
      UL = bittst(ltb,2**15)
C ... use U averaged over each channel at a site
      Uav = (.not. bittst(ltb,2**14))
      call rxx(UL .and. .not. Uav,
     .  ' TBESEL: set NOUAVG=F in ctrl.'//
     .  ' For an orbital dependent potential, use TB+U')
C ... TB+U:
      TBU = bittst(ltb,2**13)
C ... ignore non diagonal density matrix element (mostly for debugging)
      diagn = cmdopt('--diagn',7,0,strn)
C ... keep point charges only ---
      if (bittst(ltb,2**9)) then
        point = .true.
      else
        point = cmdopt('--point',7,0,strn)
      endif
C --- sanity test ---
      if (TBU .and. UL) then
        call rx('TBESEL: cannot have UL and TB+U')
      endif

C --- loop to get dipole moments on selected atoms (misuse RELAX=) ---
      do  ib = 1, nbas
        ic = ipc(ib)
        call upack2('site relax',ssite,ib,ifrlx)
        if (ifrlx(1) .eq. 1) then
          do  ilm = 2, 4
            ppdip(ilm-1) = ppdip(ilm-1) + qmpol(ilm,ib)
          enddo
        endif
      enddo

C --- get components of the Madelung potential ---
      call dpzero(vm,25*nbas)
      le = 0
      lmax = 2
      lmxf = 4
      if (point) then
        if (iprint() .gt. 20) then
          call awrit0(' TBESEL using point charges only',' ',128,
     .                i1mach(2))
        endif
        lmax = 0
        lmxf = 1
      endif
      lmxst = lmax + lmxf
      nlx = (lmax + 1)**2
      nlf = (lmxf + 1)**2
      call tcn('Madelung potential')
C#ifdefC MPI
C#ifdefC MPI-SCHEDULED
CC MPI: Distribute atoms ib over processes
CC MPI: Find the number >= nbas divisible by numprocs
C      allocate (bproc(0:numprocs), stat=ierr)
C      call dstrbp(nbas,numprocs,1,bproc(0))
C      bproc(numprocs) = nbas+1
C      do  ib = bproc(procid), bproc(procid+1)-1
C        if (ib .eq. bproc(procid)) then
C          if (mlog) then
C            call gettime(datim)
C            call awrit4(' MPI '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' starting atoms %i to %i',' ',256,lgunit(3),
C     .        procid,numprocs,bproc(procid),bproc(procid+1)-1)
C          endif
C        endif
C#elseC
C      do  ib =  1+procid, nbas, numprocs
C#endifC
C#else
C --- loop to make Madelung potential ---
      tpvq = 0d0
      do  ib = 1, nbas
C#endif
        ic = ipc(ib)
        do  jb = 1, nbas
          jc = ipc(jb)
C .. make strx on the fly ..
C .. tau = R_j - R_i :
          call dmadd(bas(1,jb),3,1,1d0,bas(1,ib),3,1,-1d0,tau,3,1,3,1)
          if (.not. MOL) then
            call shortn(tau,tau,dlat,nkd)
            call rcnsl0(tau,awld,lmxst,alat,glat,nkg,dlat,nkd,vol,cy,hl)
            le = 0
          else
            call soldhj(tau,0d0,0,lmxf,hl,bl,cy)
            le = 1
          endif
          if (iprint() .gt. 60)
     .      call awrit2('tau=%3:1d, rstrx=%49:1d',' ',1028,i1mach(2),
     .                   tau,hl)
          call hstr(pv,le,strx,drstrx,nlf,nlx,nlf,hl,49,cg,indxcg,
     .              jcg,vol)
          call xxxind(strx)
          call xxxfac(strx)
          if (pv) then
            call xxxind(drstrx)
            call xxxfac(drstrx)
          endif
          do  ilmp = 1, 9
            do  ilm = 1, 25
              vm(ilm,ib) = vm(ilm,ib) +
     .                     2d0*strx(ilm,ilmp)*qmpol(ilmp,jb)
            enddo
            if (pv) then
              do  ilm = 1, 9
                tpvq = tpvq - qmpol(ilmp,jb) * qmpol(ilm,ib) *
     .                         2d0*drstrx(ilm,ilmp)
              enddo
            endif
          enddo
          if (lov) then
            A(ib) = A(ib) + 2d0*strx(1,1)*qmpol(1,jb)
          endif
        enddo
        if (lov) then
          if (.not. Uav) call rx('TBESEL: cannot yet have OVLP and U_l')
          B(ib) = getavU(nl,nsp,qnu,idxdn,ic) * qmpol(1,ib)
        endif
C .. add electric field
        if (field) then
          vm(1,ib) = vm(1,ib)
     .             + sqrt(2d0)*alat*ddot(3,efield,1,bas(1,ib),1)
          vm(2,ib) = vm(2,ib) + sqrt(2d0)*efield(1)
          vm(3,ib) = vm(3,ib) + sqrt(2d0)*efield(2)
          vm(4,ib) = vm(4,ib) + sqrt(2d0)*efield(3)
        endif
      enddo
      tpvq = tpvq / 2d0
      if (lov) then
C#ifdefC MPI
C      call rx('TBESEL: bug to fix, parallelise this loop. Fix A.')
C#endif
        if (.not. Uav) call rx('TBESEL: not ready for OVLP and U_l')
        do  ib = 1, nbas
          pot0(ib) = A(ib) + B(ib)
        enddo
      endif            
C#ifdefC MPI
CC MPI: Wait here till all the atoms are done
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPI-SCHEDULED
CC MPI: Broadcast the potential for your atoms
C      do  pid = 0, numprocs-1
C        ib = bproc(pid)
C        jb = bproc(pid+1) - ib
C        if (procid .eq. pid) then
C          if (mlog) then
C            call gettime(datim)
C            call awrit5(' tbesel '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' bcast atoms %i to %i, %i numbers',' ',256,lgunit(3),
C     .        procid,numprocs,bproc(procid),bproc(procid+1)-1,25*jb)
C          endif
C        endif
C        call MPI_BCAST(vm(1,ib),25*jb,MPI_DOUBLE_PRECISION,pid,
C     .                 MPI_COMM_WORLD,ierr)
C      enddo
C      deallocate (bproc, stat=ierr)
C#elseC
C      allocate(buffer(25*nbas), stat=ierr)
C      call MPI_ALLREDUCE(vm,buffer,25*nbas,MPI_DOUBLE_PRECISION,
C     .                   MPI_SUM,MPI_COMM_WORLD,ierr)
C      call dcopy(25*nbas,buffer,1,vm,1)
C      deallocate(buffer, stat=ierr)
C#endifC
C#endif
      call tcx('Madelung potential')

      nlm = nl**2
      call dcopy(nl**4*nbas*nsp,0d0,0,dh,1)
      call dcopy(3*nbas,0d0,0,Ratm,1)
C ... "U * dN^2" contribution to E_2 or E^U:
      sumU = 0d0
      sumJ = 0d0
C --- get components of the Hubbard potential and double counting ---
      if (TBU) then
        do  ib = 1, nbas
          ic = ipc(ib)
C --- assemble matrix of Slater integrals ---
          call dcopy(19683,0d0,0,Rk,1)
          do  ilm1 = 1, nlm
            do  ilm2 = 1, nlm
              do  ilm3 = 1, nlm
                do  ilm4 = 1, nlm
                  l1 = ll(ilm1)
                  l2 = ll(ilm2)
                  l3 = ll(ilm3)
                  l4 = ll(ilm4)
C --- first approximation: R^k(ll'l''l''')=R^k(llll)delta(ll'l''l''') 
                  if (l1.eq.l2 .and. l1.eq.l3 .and. l1.eq.l4) then
                    l = l1
                    U = qnu(3,l,1,ic)
                    Rk(ilm1,ilm2,ilm3,ilm4,1) = U
                    if (l .eq. 2) then
                      F2 = (14d0/(4d0*2.6d0))*(5d0*stni(ic) - U)
                      call rxx(F2 .lt. 0,
     .                         ' TBESEL: U, I parameter error. F2<0')
                      F4 = 1.6d0*F2
                      Rk(ilm1,ilm2,ilm3,ilm4,2) = F2
                      Rk(ilm1,ilm2,ilm3,ilm4,3) = F4
                      do  i = 1, 3
                        Ratm(ib,i) =  Rk(ilm1,ilm2,ilm3,ilm4,i)
                      enddo
                    else
                      Ratm(ib,l+1) = U
                    endif
                  endif
                enddo
              enddo
            enddo
          enddo
C --- assemble matrix of interaction parameters V = R^k * A^k ---
          if (iprint() .ge. 60) write (*,10)
          call dcopy(6561,0d0,0,VL,1)
          do  ilm1 = 1, nlm
            do  ilm2 = 1, nlm
              do  ilm3 = 1, nlm
                do  ilm4 = 1, nlm
                  if ((idxdn(ll(ilm1),ic) .eq. 1) .and.
     .                (idxdn(ll(ilm2),ic) .eq. 1) .and.
     .                (idxdn(ll(ilm3),ic) .eq. 1) .and.
     .                (idxdn(ll(ilm4),ic) .eq. 1)) then
                    do  i = 1, 3
                     VL(ilm1,ilm2,ilm3,ilm4) = VL(ilm1,ilm2,ilm3,ilm4) +
     .               Rk(ilm1,ilm2,ilm3,ilm4,i)*Ak(ilm1,ilm2,ilm3,ilm4,i)
                    enddo

C --- verbose output ---
                    if (iprint() .ge. 60 .and.
     .                Rk(ilm1,ilm2,ilm3,ilm4,1) .gt. d1mach(3))
     .                write(*,20)
     .                ilm1,ilm2,ilm3,ilm4,ll(ilm1),ll(ilm2),ll(ilm3),
     .                ll(ilm4),(Rk(ilm1,ilm2,ilm3,ilm4,k),k=1,3),
     .                (Ak(ilm1,ilm2,ilm3,ilm4,k),k=1,3),
     .                VL(ilm1,ilm2,ilm3,ilm4)
C ---------------------

                  endif
                enddo
              enddo
            enddo
          enddo
C --- put the on-site potential directly into dh ---
          if (iprint() .gt. 40) write (*,50)
          do  isp = 1, nsp
            do  ilm1 = 1, nlm
              do  ilm2 = 1, nlm
                do  ilm3 = 1, nlm
                  do  ilm4 = 1, nlm
                    if ((idxdn(ll(ilm1),ic) .eq. 1) .and.
     .                  (idxdn(ll(ilm2),ic) .eq. 1) .and.
     .                  (idxdn(ll(ilm3),ic) .eq. 1) .and.
     .                  (idxdn(ll(ilm4),ic) .eq. 1)) then
                      if (diagn) then
                        if (ilm1 .ne. ilm2 .or. ilm3 .ne. ilm4) goto 1
                      endif
                      U = VL(ilm1,ilm3,ilm4,ilm2)
                      J = VL(ilm1,ilm3,ilm2,ilm4)
                      drho1 = rhon(ilm3,ilm4,ib,2-isp/2)
                      drho2 = rhon(ilm3,ilm4,ib,isp)
                      if (ilm3 .eq. ilm4) then
                        l = ll(ilm3)
                        q0 = qnu(1,l,isp,ic) / (2*l + 1)
                        drho1 = drho1 - q0
                        drho2 = drho2 - q0
                      endif
                      drho3 = rhon(ilm1,ilm2,ib,isp)
                      if (ilm1 .eq. ilm2) then
                        l = ll(ilm1)
                        q0 = qnu(1,l,isp,ic) / (2*l + 1)
                        drho3 = drho3 - q0
                      endif
                      dh(ilm1,ilm2,ib,isp) = dh(ilm1,ilm2,ib,isp) +
     .                                       U*drho1 + (U-J)*drho2
                      sumU = sumU + (U*drho1 + (U-J)*drho2) * drho3

C --- verbose output ---
                      if ((iprint() .ge. 50 .and.
     .                (dabs(U)+dabs(J)) .gt. d1mach(3)) .or.
     .                (iprint() .gt. 40
     .                .and. ilm3 .eq. ilm4 .and. ilm1 .eq. ilm2))
     .                write (*,75) 
     .                ilm1,ilm2,ilm3,ilm4,ll(ilm1),ll(ilm2),ll(ilm3),
     .                ll(ilm4),isp,drho1,drho2,U,J,U-J
C ----------------------

    1                 continue 
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      else
C --- TB-L branch ---
        do  ib = 1, nbas
          ic = ipc(ib)
          if (Uav) then
            U = getavU(nl,nsp,qnu,idxdn,ic)
            dQtot = qmpol(1,ib)
            sumU = sumU +  U * dQtot*dQtot
          endif          
          do  l = 0, nl-1
            vu(l,ib) = 0d0
            if (Uav) then
C ... no real need for this:
C              if (idxdn(l,ic) .eq. 1) then
                vu(l,ib) = U * dQtot
C              endif
            else
              do  lp = 0, nl-1
                if (lp .eq. l) then
                  U = qnu(3,l,1,ic)
                else
                  U = dmin1(qnu(3,l,1,ic),qnu(3,lp,1,ic))
                endif
                vu(l,ib) = vu(l,ib) 
     .                      + U * (rho(lp+1,1,ib) - qnu(1,lp,1,ic))
                sumU = sumU + U * ((rho(l+1,1,ib) - qnu(1,l,1,ic)))
     .                          * ((rho(lp+1,1,ib) - qnu(1,lp,1,ic)))
              enddo
            endif
          enddo
C --- Spin polarised TB-L ---
          if (nsp .eq. 2) then
            J = getavJ(nl,jh(1,ib),idxdn,ic)
            vj(1,ib) = - 0.5d0 * J * (qmpol(1,ib) + mmom(ib))
            vj(2,ib) = - 0.5d0 * J * (qmpol(1,ib) - mmom(ib))
            sumJ = sumJ - J * (0.25d0 * (qmpol(1,ib) + mmom(ib))**2
     .                      +  0.25d0 * (qmpol(1,ib) - mmom(ib))**2)
          endif
        enddo
      endif

C --- increments to hamiltonian ---
      call tcn('delta H')
      do   ib = 1, nbas
        do  isp = 1, nsp
          if (iprint() .gt. 50 .and. isp .eq. 1) then
            if (UL) then
              write (*,150) 
            else
              write (*,175) 
            endif
          endif
          ic = ipc(ib)
          do  ilmp = 1, nl**2
            do  ilmpp = 1, nl**2
              if ((idxdn(ll(ilmp) ,ic) .eq. 1) .and.
     .            (idxdn(ll(ilmpp),ic) .eq. 1)) then
                if (.not. TBU) then
                  if (ilmp .eq. ilmpp) then
                    dh(ilmp,ilmpp,ib,isp) = dh(ilmp,ilmpp,ib,isp) +
     .                vu(ll(ilmp),ib) + vj(isp,ib)
                  endif
                endif
                do  ilm = 1, 25
                  call getM(ilm,ilmp,ilmpp,qpol(1,ic),M)
                  dh(ilmp,ilmpp,ib,isp) = dh(ilmp,ilmpp,ib,isp) +
     .                            vm(ilm,ib) * M * gaunt(ilmp,ilmpp,ilm)

C --- verbose output ---
                  if (iprint() .gt. 50 .and. M .ne. 0d0 .and. isp .eq. 1
     .                .and. gaunt(ilmp,ilmpp,ilm) .ne. 0d0) then 
                    if (ilmp .eq. ilmpp .and. UL) then
                      U = getavU(nl,nsp,qnu,idxdn,ic)
                      dQtot = dsum(nl,rho(1,1,ib),1)
     .                      - dsum(nl,qnu(1,0,1,ic),3)
                      write (*,250) 
     .                  ilmp,ilmpp,ilm,ll(ilmp),ll(ilmpp),ll(ilm),
     .                  M,gaunt(ilmp,ilmpp,ilm),vm(ilm,ib),
     .                  dQtot,U,vu(ll(ilmp),ib),dh(ilmp,ilmpp,ib,isp)
                    else
                      write (*,210) 
     .                  ilmp,ilmpp,ilm,ll(ilmp),ll(ilmpp),ll(ilm),
     .                  M,gaunt(ilmp,ilmpp,ilm),vm(ilm,ib),
     .                  dh(ilmp,ilmpp,ib,isp)
                    endif
                  endif
C ---------------------

                enddo
              endif
            enddo
          enddo
        enddo
      enddo
      call tcx('delta H')

C --- electrostatic energy and force ---
      if (bittst(ltb,2**4)) then
        call dpzero(f,3*nbas)
        if (lov) then
          call dpzero(fnou,3*nbas)
          call dpzero(fnom,3*nbas)
        endif
      endif
C ... dQ * V :
      sumV = 0d0
      do  ib = 1, nbas
        ic = ipc(ib)
        do  ilm = 1, 9
          sumV = sumV + qmpol(ilm,ib) * vm(ilm,ib)
        enddo
C ... Forces :
        if (bittst(ltb,2**4)) then
          do  i = 1, 3
            ilm = i+1
            f(i,ib) = f(i,ib) - vm(ilm,ib)*qmpol(1,ib)
C ... sum_m' Q_1m' * (d/dx_i)V_1m'
            do  ilmp = 2, 4
C ... make (d/dx_i)V_1m'
              call dcopy(3,0d0,0,dV1,1)                
              do  ilmpp = 5, 9
                dV1(ilmp) = dV1(ilmp) + 
     .               gaunt(ilm,ilmp,ilmpp)*vm(ilmpp,ib)
              enddo
              f(i,ib) = f(i,ib) -
     .                    dsqrt(20d0*pi)*qmpol(ilmp,ib)*dV1(ilmp)
            enddo
C ... sum_m' Q_2m' * (d/dx_i)V_2m'
            do  ilmp = 5, 9
C ... make (d/dx_i)V_2m'
              call dcopy(5,0d0,0,dV2,1)                
              do  ilmpp = 10, 16
                dV2(ilmp) = dV2(ilmp) + 
     .               gaunt(ilm,ilmp,ilmpp)*vm(ilmpp,ib)
              enddo
              f(i,ib) = f(i,ib) -
     .                    dsqrt(140d0*pi/3d0)*qmpol(ilmp,ib)*dV2(ilmp)
            enddo
          enddo
        endif
      enddo

      ecorr = 0.5d0*(sumV + sumU + sumJ)

C --- make terms in dervatives of the overlap for forces ---
      if (lov) then
        call dcopy(3*nbas**2,0d0,0,drhos,1)
        opi = 0
        do  ib = 1, nbas
          ic = ipc(ib)
          nlmi = getnlm(nl,idxdn,ic)
          opj = 0
          do  jb = 1, nbas
            jc = ipc(jb)
            nlmj = getnlm(nl,idxdn,jc)
            do  k = 1, 3
              drhos(ib,jb,k) = 0d0
              do  ilmi = 1, nlmi
                do  ilmj = 1, nlmj
                  drhos(ib,jb,k) = drhos(ib,jb,k)
     .              + drhosl(opi+ilmi,opj+ilmj,k)
                enddo
              enddo
            enddo
            opj = opj + nlmj
          enddo
          opi = opi + nlmi
        enddo
        call rxx(opi .ne. ldim,' bug in tbesel: opi .ne. ldim')
        call rxx(opj .ne. ldim,' bug in tbesel: opj .ne. ldim')
        if (iprint() .gt. 40) then
          print *,   ' TBESEL: drhos:'
          do  k = 1, 3
            if (k. eq. 1) print *, '         x :'
            if (k. eq. 2) print *, '         y :'
            if (k. eq. 3) print *, '         z :'
            do  ib = 1, nbas
              write(*,400) (drhos(ib,jb,k),jb=1,nbas)
            enddo
          enddo
        endif

C --- forces from overlap ---
        do  ib = 1, nbas
          do  jb = 1, nbas
            do  k = 1, 3
              fnou(k,ib) = fnou(k,ib) - (B(ib) + B(jb)) * drhos(ib,jb,k)
              fnom(k,ib) = fnom(k,ib) - (vm(1,ib) + vm(1,jb))
     .                                * drhos(ib,jb,k)
            enddo
          enddo
        enddo
      endif
      if (bittst(ltb,2**4) .and..not. MOL) then
        call symfor(nbas,1,symgr,nsgrp,istab,fwk,f)
        if (lov) then
          call symfor(nbas,1,symgr,nsgrp,istab,fwk,fnou)
          call symfor(nbas,1,symgr,nsgrp,istab,fwk,fnom)
        endif
      endif
      
      call tcx('tbesel')

C --- printout ---
      if (iprint() .lt. 30) return
      if (iprint() .eq. 30) goto 1000
      do ib = 1, nbas
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
        print *
        if (mixQ) then
          call awrit0('Atom '//clabl,' ',180,i1mach(2))
        else
          if (nsp .eq. 2) then
            if (nl .eq. 3) then
              call awrit8('Atom '//clabl//
     .        '%cN_s=%d+%d N_p=%d+%d N_d=%d+%d Q=%d mom=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(1,2,ib),rho(2,1,ib),
     .        rho(2,2,ib),rho(3,1,ib),rho(3,2,ib),
     .        dsum(nl,rho(1,1,ib),1)+dsum(nl,rho(1,2,ib),1),
     .        dsum(nl,rho(1,1,ib),1)-dsum(nl,rho(1,2,ib),1))
            elseif (nl .eq. 2) then
              call awrit6('Atom '//clabl//
     .        '%cN_s=%d+%d N_p=%d+%d Q=%d mom=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(1,2,ib),rho(2,1,ib),
     .        rho(2,2,ib),
     .        dsum(nl,rho(1,1,ib),1)+dsum(nl,rho(1,2,ib),1),
     .        dsum(nl,rho(1,1,ib),1)-dsum(nl,rho(1,2,ib),1))
            elseif (nl .eq. 1) then
              call awrit4('Atom '//clabl//
     .        '%cN_s=%d+%d Q=%d mom=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(1,2,ib),
     .        dsum(nl,rho(1,1,ib),1)+dsum(nl,rho(1,2,ib),1),
     .        dsum(nl,rho(1,1,ib),1)-dsum(nl,rho(1,2,ib),1))
            endif
          else
            if (nl .eq. 3) then
              call awrit4('Atom '//clabl//
     .        '%cN_s=%d N_p=%d N_d=%d Q=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(2,1,ib),
     .        rho(3,1,ib),dsum(nl,rho(1,1,ib),1))
            elseif (nl .eq. 2) then
              call awrit3('Atom '//clabl//
     .        '%cN_s=%d N_p=%d Q=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),rho(2,1,ib),
     .        dsum(nl,rho(1,1,ib),1))
            elseif (nl .eq. 1) then
              call awrit2('Atom '//clabl//
     .        '%cN_s=%d Q=%d',
     .        ' ',180,i1mach(2),rho(1,1,ib),
     .        dsum(nl,rho(1,1,ib),1))
            endif
          endif
        endif
        if (iprint() .gt. 31) then
          if (TBU) then
            print *, 'c*_RL S c_RL'':'
            print *, ' spin up:'
            do  i = 1, nl**2
              write (*,300) (rhon(i,k,ib,1),k=1,nl**2)
            enddo
            print *, ' spin down:'
            do  i = 1, nl**2
              write (*,300) (rhon(i,k,ib,2),k=1,nl**2)
            enddo
          endif
          print *, 'c*_RL c_RL'':'
          do  i = 1, nl**2
            write (*,300) (rhoc(i,k,ib),k=1,nl**2)
          enddo
        endif
        if (TBU) then
          if (nl .eq. 3) then
            call awrit5('        U=%d, F^2=%d, F^4=%d, J=%d, I=%d',' ',
     .                  180,i1mach(2),Ratm(ib,1),Ratm(ib,2),Ratm(ib,3),
     .                  (Ratm(ib,2)+Ratm(ib,3))/14d0,stni(ic))
          else
            call awrit2('        U_s=%d, U_p=%d, F^2=0',' ',180,
     .                  i1mach(2),Ratm(ib,1),Ratm(ib,2))
          endif
        endif
        call awrit7('        M=%d %d %d %d %d %d %d',' ',180,i1mach(2),
     .   sqrt(4*pi/3)*qpol(1,ic),sqrt(4*pi/5)*qpol(2,ic),
     .   sqrt(4*pi/5)*qpol(3,ic),sqrt(4*pi/3)*qpol(4,ic),
     .   sqrt(4*pi/5)*qpol(5,ic),sqrt(4*pi/7)*qpol(6,ic),
     .   sqrt(4*pi/9)*qpol(7,ic))
        call awrit3('        Q^e/e=%d, %3:1d, %5:1d',' ',180,
     .                i1mach(2),qmpol(1,ib),qmpol(2,ib),qmpol(5,ib))
        qmp(1) = qmpol(1,ib)/dsqrt(4d0*pi)
        call dpcopy(qmpol(2,ib),qmp(2),1,3,dsqrt(3d0/(4d0*pi)))
        call dpcopy(qmpol(5,ib),qmp(5),1,5,dsqrt(5d0/(4d0*pi)))
        call awrit3('        Qmpol=%d, %3:1d, %5:1d',' ',180,
     .                i1mach(2),qmp(1),qmp(2),qmp(5))
        call dcopy(25,vm(1,ib),1,v,1)
        call awrit3('        e dV=  %d, %3:1d, %5:1d',' ',180,
     .              i1mach(2),v,v(2),v(5))
        call awrit1('              %7:1d',' ',180,i1mach(2),v(10))
        call awrit1('              %9:1d',' ',180,i1mach(2),v(17))
        if (field) then
          call awrit1('               includes E=%3:-2,4;4d',
     .                ' ',120,i1mach(2),efield)
        endif
        if (.not. TBU) then
          U = getavU(nl,nsp,qnu,idxdn,ic)
          J = getavJ(nl,jh(1,ib),idxdn,ic)
          if (nl .eq. 3) then
            dq(0) = rho(1,1,ib) - qnu(1,0,1,ic)
            dq(1) = rho(2,1,ib) - qnu(1,1,1,ic)
            dq(2) = rho(3,1,ib) - qnu(1,2,1,ic)
            if (Uav) then
              call awrit3
     .        ('        Hubbard potential U=%d: '//
     .          '%1:1d, dq:%1:1d',' ',180,i1mach(2),
     .        U,vu(0,ib),qmpol(1,ib))
              if (nsp .eq. 2) then
                call awrit3
     .          ('        Stoner potentials J=%d: '//
     .            '%2:1d, moment: %d',' ',180,i1mach(2),
     .          J,vj(1,ib),mmom(ib))
              endif
            else
              call awrit6
     .        ('        Hubbard potentials U=%d %d %d: '//
     .          '%3:1d, dq:%3:1d (total=%d)',' ',180,i1mach(2),
     .        qnu(3,0,1,ic),qnu(3,1,1,ic),qnu(3,2,1,ic),
     .        vu(0,ib),dq,dsum(3,dq,1))
            endif
          endif
          if (nl .eq. 2) then
            dq(0) = rho(1,1,ib) - qnu(1,0,1,ic)
            dq(1) = rho(2,1,ib) - qnu(1,1,1,ic)
            if (Uav) then
              call awrit3
     .        ('        Hubbard potential U=%d: '//
     .          '%1:1d, dq:%1:1d ',' ',120,i1mach(2),U,
     .        vu(0,ib),qmpol(1,ib))
              if (nsp .eq. 2) then
                call awrit3
     .          ('        Stoner potentials J=%d: '//
     .            '%2:1d, moment: %d',' ',180,i1mach(2),
     .          J,vj(1,ib),mmom(ib))
              endif
            else
              call awrit5
     .        ('        Hubbard potentials U=%d %d: '//
     .          '%2:-2d, dq:%2:-2d (total=%d)',' ',180,i1mach(2),
     .        qnu(3,0,1,ic),qnu(3,1,1,ic),
     .        vu(0,ib),dq,dsum(2,dq,1))
            endif
          endif
          if (nl .eq. 1) then
            dq(0) = rho(1,1,ib) - qnu(1,0,1,ic)
            call awrit3
     .        ('        Hubbard potential U=%d: '//
     .        '%d, dq: %d ',' ',120,i1mach(2),U,vu(0,ib),dq)
          endif
          do isp = 1, nsp
            call dcopy(nl**2,dh(1,1,ib,isp),nl**2+1,dhd,1)
            if (isp .eq. 1) then
              if (nl .gt. 2) then
                call awrit3('        DH^e_LL=%1:-2d, %3:-2d, %5:-2d',
     .            ' ',180,i1mach(2),dhd,dhd(2),dhd(5))
              elseif (nl .gt. 1) then
                call awrit2('        DH^e_LL=%1:-2d, %3:-2d',' ',180,
     .            i1mach(2),dhd,dhd(2))
              endif
            else
              if (nl .gt. 2) then
                call awrit3('                %1:-2d, %3:-2d, %5:-2d',
     .            ' ',180,i1mach(2),dhd,dhd(2),dhd(5))
              elseif (nl .gt. 1) then
                call awrit2('                %1:-2d, %3:-2d',' ',180,
     .            i1mach(2),dhd,dhd(2))
              endif
            endif
          enddo
          if (lov .and. iprint() .ge. 40) then
            call awrit3('        Overlap terms: B=%d A=%d pot0=%d',' ',
     .                  120,i1mach(2),B(ib),A(ib),pot0(ib))
          endif
        endif
        if (iprint() .gt. 40 .or. (TBU .and. iprint() .gt. 31)) then
          print *, 'DH^e_LL'':'
          do  isp = 1, nsp
            if (nsp .eq. 2) then
              if (isp .eq. 1) print*, ' spin up:'
              if (isp .eq. 2) print*, ' spin down:'
            endif
            do  i = 1, nl**2
              write (*,300) (dh(i,k,ib,isp),k=1,nl**2)
            enddo
          enddo
        endif
C --- Electric field gradient ---
C        if (dasum(5,v(5),1) .gt. 10*dsqrt(d1mach(3))) then
C          efg(1) = v(5)*dsqrt(4d0*pi/5d0)
C          efg(2) = v(6)*dsqrt(4d0*pi/5d0)
C          efg(3) = v(9)*dsqrt(4d0*pi/5d0)
C          efg(4) = v(7)*dsqrt(4d0*pi/5d0)
C          efg(5) = v(8)*dsqrt(4d0*pi/5d0)
C          call efgrad(z(ic),efg)
C        endif
        if (bittst(ltb,2**4) .and. iprint() .ge. 20)
     .    call awrit1('        Force=%3:1d',' ',120,i1mach(2),f(1,ib))
      enddo
 1000 continue 
      print *
      if (iprint() .ge. 30) then
        if (TBU) then
          write(*,425) sumV/2d0,sumU/2d0
        else
          if (nsp .eq. 1) then
            write(*,450) sumV/2d0,sumU/2d0
          else
            write(*,475) sumV/2d0,sumU/2d0,sumJ/2d0
          endif
        endif
        write (*,500) ecorr
      endif

   10 format ('  L    L''   L''''  L'''''' l    l''   l''''  l''''''',
     .14x,' R^k',25x,' A^k',10x,'  V (accumulated)')
   20 format (8(i3,2x),7f10.6)
   50 format ('  L    L''   L''''  L'''''' l    l''   l''''  l'''''' s
     .      dn(-s)    dn(s)        U        J       U-J')
   75 format (9(i3,2x),5f10.6)
C 100 format ('  L''   L''''  L    l''   l''''  l      M         CG
C    .      rho_L''L''''')
  150 format ('  L''   L''''  L    l''   l''''  l      M         CG
     .   V_m        dQ      U     V_u       DH')
  175 format ('  L''   L''''  L    l''   l''''  l      M         CG
     .   V_m       DH')
C 120 format(1028g12.4)
C 200 format (6(i3,2x),2(2x,f6.2,2x),2f10.6)
  210 format (6(i3,2x),2(2x,f6.2,2x),1x,f10.6,25x,f10.6)
  250 format (6(i3,2x),2(2x,f6.2,2x),2f11.6,f4.1,2f10.6)
  300 format (5x,9f10.6)
  400 format (1028f10.6)
  425 format ('   (1/2) dQ dV             : ',f12.6/
     .        '   E_U                     : ',f12.6)
  450 format ('   (1/2) dQ dV             : ',f12.6/
     .        '   (1/2) U dN^2            : ',f12.6)
  475 format ('   (1/2) dQ dV             : ',f12.6/
     .        '   (1/2) U dN^2            : ',f12.6/
     .        '  -(1/2) J dN_u^2+dN_d^2   : ',f12.6)
  500 format ('   E_2                     : ',f12.6)
      end
      subroutine getM(ilm,ilmp,ilmpp,qpol,M)
C- 
C ----------------------------------------------------------------------
Ci Inputs:
Ci   
Co Outputs:
Co   
Cr Remarks
Cr  The tight binding parameters in qpol are as follows
Cr  qpol(1) = M(011) = M(101)
Cr  qpol(2) = M(112)
Cr  qpol(3) = M(022) = M(202)
Cr  qpol(4) = M(121) = M(211)
Cr  qpol(5) = M(222)
Cr  qpol(6) = M(123)
Cr  qpol(7) = M(224)
Cr  These are converted into Stone's definitions by multiplying the
Cr  values from ctrl by \sqrt{4\pi/(2l+1)}
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer ilm,ilmp,ilmpp
      double precision qpol(10),M
C Local Variables
      integer ll,l,lp,lpp
      double precision fourpi,sqr4pi,fac,dsqrt,datan
      fourpi = 16d0*datan(1d0)
      sqr4pi = dsqrt(fourpi)

      l = ll(ilm)
      lp = ll(ilmp)
      lpp = ll(ilmpp)
      fac = dsqrt(fourpi/(2*l + 1))
      M = 0d0
      if (lp .eq. 1 .and. lpp .eq. 0 .and. l .eq. 1) M = qpol(1)
      if (lp .eq. 0 .and. lpp .eq. 1 .and. l .eq. 1) M = qpol(1)
      if (lp .eq. 1 .and. lpp .eq. 1 .and. l .eq. 2) M = qpol(2)
      if (lp .eq. 2 .and. lpp .eq. 0 .and. l .eq. 2) M = qpol(3)
      if (lp .eq. 0 .and. lpp .eq. 2 .and. l .eq. 2) M = qpol(3)
      if (lp .eq. 2 .and. lpp .eq. 1 .and. l .eq. 1) M = qpol(4)
      if (lp .eq. 1 .and. lpp .eq. 2 .and. l .eq. 1) M = qpol(4)
      if (lp .eq. 2 .and. lpp .eq. 2 .and. l .eq. 2) M = qpol(5)
      if (lp .eq. 1 .and. lpp .eq. 2 .and. l .eq. 3) M = qpol(6)
      if (lp .eq. 2 .and. lpp .eq. 1 .and. l .eq. 3) M = qpol(6)
      if (lp .eq. 2 .and. lpp .eq. 2 .and. l .eq. 4) M = qpol(7)
      M = fac*M
      if (l .eq. 0) then
        if (lp .eq. lpp) then
          M = sqr4pi
        endif
      endif
      end

      double precision function getavU(nl,nsp,qnu,idxdn,ic)
C- Return average U over open channels
C     implicit none
C Passed parameters
      integer nl,nsp,ic,idxdn(nl,1)
      double precision qnu(3,nl,nsp,1)
C Local variables
      integer l,isp,nchan

      nchan = 0
      getavU = 0d0
      do  isp = 1, nsp
        do  l = 1, nl
          if (idxdn(l,ic) .eq. 1) then
            nchan = nchan + 1
            getavU = getavU + qnu(3,l,isp,ic)
          endif
        enddo
      enddo
      if (nchan. gt. 0) getavU = getavU / nchan
      end

      double precision function getavJ(nl,jh,idxdn,ic)
C- Return average J over open channels
C     implicit none
C Passed parameters
      integer nl,ic,idxdn(nl,1)
      double precision jh(4)
C Local variables
      integer l,nchan

      nchan = 0
      getavJ = 0d0
      do  l = 1, nl
        if (idxdn(l,ic) .eq. 1) then
          nchan = nchan + 1
          getavJ = getavJ + jh(l)
        endif
      enddo
      if (nchan. gt. 0) getavJ = getavJ / nchan
      end

      integer function getnlm(nl,idxdn,ic)
C- Return number of lm's at species ic
C     implicit none
      integer nl,ic,idxdn(0:nl-1,1),nlm,l
      nlm = 0
      do  l = 0, nl-1
       if (idxdn(l,ic) .eq. 1) nlm = nlm + (2*l + 1)
      enddo
      getnlm = nlm
      end

      subroutine xxxind(dm)
C- Shuffle 25X9 matrix from MSM lm indices to TBE lm indices
C ----------------------------------------------------------------------
Ci Inputs:
Ci   dm : 25X9 structure constant sub block
Co Outputs:
Co   dm reordered
Cr Remarks
Cr   Michael's structure constants are ordered according to the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     y     z     x    xy    yz 3z^2-r^2 zx   x^2-y^2
Cr   while the TBE programs use the scheme
Cr         1     2     3     4     5     6     7     8     9
Cr         1     x     y     z    xy    yz    zx  x^2-y^2  3z^2-r^2 
Cr   This widget rearranges the matrix made by Michael's HBLSTR into
Cr   the TBE order.
Cr   The l>2 ordering is unchanged
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      double precision dm(25,9)
C Local Variables
      double precision wk(25,9)
      integer i, j, ind(9)
      data ind /1,4,2,3,5,6,8,9,7/

      do  1  i = 1, 9
      do  1  j = 1, 9
        wk(i,j) = dm(ind(i),ind(j))
    1 continue 
      do  2  i = 10, 25
      do  2  j = 1, 9
        wk(i,j) = dm(i,ind(j))
    2 continue 
      call dcopy(25*9,wk,1,dm,1)
      end

      subroutine xxxfac(dm)
C- Include factor in 25X9 structure constant block
C ----------------------------------------------------------------------
Ci Inputs:
Ci   dm : 25X9 structure constant sub block
Co Outputs:
Co   dm * fac
Cr Remarks
Cr   To convert to Stone's convention the l_i,l_j element is multiplied
Cr   by sqrt((2l_j+1)*(2l_i+1))/((2l_j+1)!!(2l_i+1)!!)
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      double precision dm(25,9)
C Local Variables
      integer i,j,li,lj,df(0:4),ll
      double precision arg,fac
      external ll
C ... df(l) = (2l+1)!!
      data df /1,3,15,105,945/ 
      
      do  1  i = 1, 25
      do  1  j = 1, 9
        li = ll(i)
        lj = ll(j)
        arg = dble((2*lj+1))*dble((2*li+1))
        fac = dsqrt(arg)/(df(lj)*df(li))
        dm(i,j) = dm(i,j)*fac
    1 continue
      end
