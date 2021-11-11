
      subroutine chksgr(opt,ltrans,kcplx,plat,nsp,ndimh,ifis,hreal,
     .iprmb,nttabs,iaxs,hrs,ndhrs,tols,nqp,qplst,nbas,lrot,rotm,
     .delT)
      use m_lmfinit,only: nl,mxorb
      use m_hamindex,only : symops_af
C- Compare Bloch summed sigma against file sigma(k)
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 calculate sigma at list of qp in file,
Ci         :  compare against file sigma
Ci         :1 calculate sigma at list of qp, write to file sigma
Ci         :  NB: caller must write header to sigma file before
Ci         :  calling chksgr!
Ci         :2 like opt=0, but use qplst instead of qp
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ndimh :dimension of hamiltonian
Ci   ifis  :file logical unit for sigma
Ci   hreal :0 sigm(R) is complex
Ci         :1 sigm(R) is real
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nttabs:number of site pairs in sigm(R)
Ci   iaxs  :neighbor table for sigm(R)
Ci   hrs   :sigm(R)  <=================================== important input to be checked
Ci   ndhrs :leading dimensions of hrs
Ci   tols  :deviations in sigma > tols causes chksgr to abort
Ci         :tols=0 -> chksgr never aborts
Ci   qplst :(opt=0) not used
Ci         :(opt=1) list of qp at which to make sigma
Ci   lrot  :(opt=0) not used
Ci         :(opt=1) 0 -> do not rotate sigma matrix
Ci         :        1 -> rotate sigma matrix by rotm
Ci         :        2 -> additionally scale sigma matrix by phase
Ci         :             owing to change in shortening of basis
Ci         :             vectors
Ci   rotm  :(opt=0) not used
Ci         :(opt=1) list of qp at which to make sigma
Ci   delT  :(only used if lrot=2)
Ci         :change in lattice translation vectors on rotation
Ci         :to add phase shift to sigm.
Cio Inputs/Outputs
Cio  nqp   :Input, opt=0:  number of qp for which sigma was checked
Cio        :Output, opt=1: number of qp for which to make sigma
Cio  (opt=1) sigma written to logical unit ifis
Cl Local variables
Cl   sq    :(opt=0) file's contents of sigma, for comparison
Cl         :(opt=1) temporary array in rotating sigma
Cl   sq1   :bloch-summed sigma
Cr Remarks
Cr
Cu Updates
Cu   20 Dec 04 enable rotation of sigma before writing
Cu   11 Jul 03 enable write of sigma
Cu   09 May 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer opt,kcplx,hreal,ifis,ndimh,nsp,nttabs,ndhrs,nqp,lrot
      integer niax,n0,nkap0
      parameter (niax=10, n0=10, nkap0=3)
      integer iprmb(*),iaxs(niax,*),nbas
      double precision plat(3,3),tols,qplst(3,nqp),rotm(3,3),
     .delT(3,nbas)
      double complex hrs(ndhrs,ndhrs,nsp,nttabs)
      logical ltrans
C ... Local parameters
      logical lsplts
      integer i,j,isp,stdo,iq1,ipr,lgunit,isw,opt0,opt01,
     .nscnts
      integer ib,jb,is,norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
C     these are read by sigma header
      integer nk1,nk2,nk3,lshft(3),offh(nbas),ndim(nbas),offi,offj,
     .ndi,ndj
C     integer os
      double precision qp(3),errmxx,errmx,err,qpr(3),xx
      double complex hwk(ndhrs,ndhrs)
      complex(8),allocatable:: sq(:,:),sq1(:,:),uz(:,:)
      integer ::iwdummy
      integer ii,jj
      integer:: ndimhr !Read from sigm. This can be different from ndimh in the case of PMT, 
                       ! that is, ndimh< ndimhr=nlmto+napw !apr2012takao
      complex(8),allocatable:: sqr(:,:)
      integer:: lb1,lb2 !mar2014takao
      logical:: laf !jun2015takao

      call getpr(ipr)
      laf=allocated(symops_af) !jun2015takao
      stdo = lgunit(1)
      opt01 = mod(mod(opt,10),2)
      opt0  = mod(opt,10)
      nscnts = 0
      allocate(sq(ndimh,ndimh),sq1(ndimh,ndimh))
C --- Rotate hrs if switch set ---
      if (opt01 .ne. 0 .and. mod(lrot,2) .ne. 0) then
        allocate(uz(ndimh,ndimh))
        call rothrm(2,ndimh,iprmb,rotm,1,nbas,ndimh,uz,uz)
        do  ib = 1, nbas
          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,j,offl,ndim(ib))
          offh(ib) = j+1
        enddo
C       For each pair and spin, rotate hrs to U(ib) hrs U(jb)+
        do is = 1, nttabs
          ib = iaxs(1,is)
          jb = iaxs(2,is)
          offi = offh(ib)
          offj = offh(jb)
          ndi = ndim(ib)
          ndj = ndim(jb)
          do  isp = 1, nsp
C           U1 sigm(unrotated) U2+
            if (ltrans) then
C             ltrans: use (U1 h+ U2+)+ = (U1 (U2 h)+)+ = U2 h U1+
              call zgemm('N','N',ndj,ndi,ndj,(1D0,0D0),Uz(offj,offj),
     .        ndimh,hrs(1,1,isp,is),ndhrs,(0D0,0D0),hwk,ndhrs)
C             U sigm(unrotated) U+
              call zgemm('N','C',ndj,ndi,ndi,(1D0,0D0),hwk,ndhrs,
     .        Uz(offi,offi),ndimh,(0D0,0D0),hrs(1,1,isp,is),ndhrs)
            else
C             Make U1 h U2+
              call zgemm('N','N',ndi,ndj,ndi,(1D0,0D0),Uz(offi,offi),
     .        ndimh,hrs(1,1,isp,is),ndhrs,(0D0,0D0),hwk,ndhrs)
C             U sigm(unrotated) U+
              call zgemm('N','C',ndi,ndj,ndj,(1D0,0D0),hwk,ndhrs,
     .        Uz(offj,offj),ndimh,(0D0,0D0),hrs(1,1,isp,is),ndhrs)
            endif

          enddo
        enddo
      endif
C ... Read header data
      lsplts = .false.
      if (opt01 .eq. 0) then
        call iosigh(0,nscnts,i,ndimhr,nk1,nk2,nk3,nqp,
     .  lshft(1),lshft(2),lshft(3),ifis)
        if (nscnts .ne. 0 .and. nscnts .ne. 5) then
          call iosigh(2,0,i,ndimhr,nk1,nk2,nk3,nqp,
     .    lshft(1),lshft(2),lshft(3),ifis)
        endif
        if (i .eq. nsp) then
          lsplts = .false.
        elseif (i .gt. nsp) then
          call rx(
     .    'rdsigm: sigm file spin polarized but calculation is not')
        else
          lsplts = .true.
        endif
      endif
      do  isp = 1, nsp
        if(laf.and.isp==2) cycle !jun2015takao
        errmxx = 0
C     If sigma file not spin polarized, use sigma from spin 1
        if (isp .eq. 2 .and. lsplts) then
          call iosigh(1,nscnts,1,ndimhr,nk1,nk2,nk3,nqp,
     .    lshft(1),lshft(2),lshft(3),ifis)
        endif
C ... For each q, compare Bloch sum of hrs against file's s(q)
C     or write to file
        do  iq1 = 1, nqp
          if (opt01 .eq. 0) then
            read(ifis) qp
            if (opt0 .eq. 2) call dcopy(3,qplst(1,iq1),1,qp,1)
            if (nscnts .eq. 5) call dpdump(xx,1,ifis)
            allocate(sqr(ndimhr,ndimhr))
            if(ndimhr<ndimh) call rx('chksgr :ndimhr<ndimh')
            call dpdump(sqr,ndimhr**2*2,ifis)
            sq = sqr(1:ndimh,1:ndimh) !!! readin sq
            deallocate(sqr)
            call ztoyy(sq,ndimh,ndimh,ndimh,ndimh,kcplx,1)
          else
            call dcopy(3,qplst(1,iq1),1,qp,1)
          endif
C       Call bloch : 104010=perm orb, transpose, c*16
          i = 100000 + 4000 + 40*(1-hreal) + 10
          if (.not. ltrans) i = i - 4000
          call bloch ( i , qp , nl , plat , mxorb , iprmb , 1 , nttabs 
     .     , iaxs , hrs , ndhrs , isp , nsp , ndimh , ndimh , 0 , ndimh 
     .     , 0 , ndimh , 0 , sq1 , iwdummy , iwdummy )
          errmx = 0
C       Compare Bloch sum of hrs against file's s(q)
          if (opt01 .eq. 0) then
            ii = 1
            jj = 1
            do  j = 1, ndimh
              do  i = 1, ndimh
                err = abs(sq(i,j)-sq1(i,j))  !Comparison here.
                if (err .gt. errmx) then
                  ii = i
                  jj = j
                  if (err .gt. tols .and. tols .gt. 0) then
                    print 357, i,j, sngl(err),sq1(i,j),sq(i,j)
  357               format(/
     .              '   i   j',6x,'diff',14x,'bloch sum',16x,'file value'/
     .              2i4,f12.6,2x,2f12.6,2x,2f12.6)
                    call rx1('rdsigm: Bloch sum deviates more than'//
     .              ' allowed tolerance (tol=%g)',tols)
                  endif
                endif
                errmx = max(errmx,err)
              enddo
            enddo
            if (ipr .ge. 45) write(stdo,332) iq1,errmx,ii,jj
  332       format(' comparing s(iq) to interpolated sq for iq=',i4,
     .      '  errmx=',1pe8.1:' at i,j=',2i4)
            errmxx = max(errmxx,errmx)
          else
C ...  Here we had 'Rotation of sigma(k)" but no longer needed, because hrs(T) was rotated now.
C   ... Write to s(q)
C         Additional scaling by phase shifts
            if (lrot .ge. 2) then
              call dgemm('T','N',3,1,3,1d0,plat,3,qp,3,0d0,qpr,3)
              call rothph(02,qpr,delT,ndimh,iprmb,1,nbas,ndimh,sq1)
            endif
C         Rotate qp to new coordinate system
            if (mod(lrot,2) .ne. 0) then
C           In-line multiply avoids bug in DEC fort compiler
              do  i = 1, 3
                qpr(i) = rotm(i,1)*qp(1) +
     .          rotm(i,2)*qp(2) +
     .          rotm(i,3)*qp(3)
              enddo
              call dcopy(3,qpr,1,qp,1)
            endif
C         Write qp, sigm for this qp
            write(ifis) qp
            call dpdump(sq1,ndimh**2*2,-ifis)
C         call zprm('sigm',2,sq1,ndimh,ndimh,ndimh)
          endif
        enddo
        if (opt01 .eq. 0)
     .  call info5(20,0,0,' check FT s(R) against s(q) ... '//
     .  'maximum error = %;2g%?#n>0# < tol (%;2g)#%j#%?#n==2# spin 2##',
     .  errmxx,isw(tols.ne.0),tols,isp,0)
      enddo
      deallocate(sq,sq1)
      if (opt01 .ne. 0 .and. mod(lrot,2) .ne. 0) deallocate(uz)
      end subroutine chksgr
      
!!sssssssssssssssssssssssssssssssssssss
      subroutine rothrm(opt,ldima,iprmb,rotm,ib1,ib2,ldU,U,Uz)
C- Generate rotation matrix for hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :  0 use real harmonics; output U is real
Ci         :  1 use true spherical harmonics; output U is complex
Ci         :  2 use real harmonics, but return U as complex
Cixxx   ldmpa :offset to first orbital in downfolding block of h to rotate
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   rotm  :3x3 cartesian rotation matrix
Ci   ib1   :U  generated for (ib1..ib2) site subblock of h
Ci   ib2   : -//-
Ci   ldU  :leading dimension of U
Co Outputs
Co   U,Uz :Rotation matrix for (ib1..ib2) subblock of h
Co        :U is the real form (opt=0); Uz is the complex form.
Co        :Only of of (U,Uz) is returned.
Cl Local variables
Cl   rYL  :rotation matrix that rotates Y_lm
Cr Remarks
Cr   This routine generates rotation matrix U for the hamiltonian matrix
Cr   which consists of blocks which transform in the row and colum indices
Cr   as Y_lm.
Cr
Cr   The unrotated and rotated hamiltonians h and hr are related by :
Cr      hr = U h U+
Cr
Cr   See also rothph, below
Cb Bugs
Cb   This routine has only been checked for ib1=1 and ib2=nbas
Cu Updates
Cu   17 Dec 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer opt,ldU,iprmb(1),ldmpa,ldima,ib1,ib2
      double precision rotm(3,3),U(ldU,ldU)
      double complex Uz(ldU,ldU)
C ... Local parameters
      integer n0,nkap0,nlx,nl2,i,ibas,iorb,j,li,
     .nlm,nlmi,norb,offri,offi
      parameter (n0=10,nkap0=3,nlx=9,nl2=nlx*nlx)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      double precision xx,rYL(nl2,nl2)
      double complex rYLz(nl2,nl2)

      if (opt .ne. 0 .and. opt .ne. 2)
     .call rx('rothrm not ready for this opt')

C ... Rotation matrix for spherical harmonics
      call ylmrtg(nl2,rotm,rYL)
C     call prmx('rothrm: rYL',rYL,nl2,nl2,nl2)

      if (opt .eq. 1) then
        call ylmrtg(nl2,rotm,rYL)
        call s2sph(1,nlx,nlx,rYL,nl2,nl2,nl2,nl2,rYLz)
C       call yprm('rYL',3,rYL,nl2*nl2,nl2,nl2,nl2)
      endif

C ... Initialize U
      if (opt .eq. 0) then
        call dpzero(U,ldU*ldU)
      else
        call dpzero(Uz,ldU*ldU*2)
      endif

C --- For each (ibas,l) block, poke rotation matrix rotating Ylm(block)
      do  ibas = ib1, ib2

C   ... Row, column offsets to h
        call orbl(ibas,0,ldima,iprmb,norb,ltab,ktab,xx,offl,nlm)

C   --- For each (ibas,l) block, poke rYL(i,j) into U(i,j) ---
C   ... loop over pairs of orbital blocks
        do  iorb = 1, norb
          li = ltab(iorb)
          nlmi = 2*li + 1
          offi = offl(iorb) !- ldmpa
          offri = li**2
          if (li+1 .gt. nlx)
     .    call rxi('rothrm: increase nlx, need nlx=',(li+1)**2)

C     ... loop over pairs of m channels within this block
          if (opt .eq. 0) then
            do  j = 1, nlmi
              do  i = 1, nlmi
                U(i+offi,j+offi) = rYL(i+offri,j+offri)
              enddo
            enddo
          elseif (opt .eq. 1) then
            do  j = 1, nlmi
              do  i = 1, nlmi
                Uz(i+offi,j+offi) = rYLz(i+offri,j+offri)
              enddo
            enddo
          else
            do  j = 1, nlmi
              do  i = 1, nlmi
                Uz(i+offi,j+offi) = rYL(i+offri,j+offri)
              enddo
            enddo
          endif

        enddo

      enddo

C      if (opt .eq. 0) then
C        call prmx('U',U,ldU,ldU,ldU)
C      else
C        call zprm('Uz',2,Uz,ldU,ldU,ldU)
C      endif
      end

      subroutine rothph(opt,q,tbas,ldima,iprmb,ib1,ib2,ldH,H)
C- Extra scaling of hamiltonian by a phase
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0 scale on left by phase
Ci         :1 scale on right by conjugate of phase
Ci         :2 scale both left and right
Ci         :10s digit
Ci         :0 phase = exp(2*pi*tbas*q)
Ci         :1 phase = exp(-2*pi*tbas*q)
Ci   q     :k-point for current h
Ci         :Note: q and tbas must be in reciprocal units
Ci   tbas  :set of basis vector shifts; see Remarks
Ci         :Note: q and tbas must be in reciprocal units
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ib1   :H  generated for (ib1..ib2) site subblock of h
Ci   ib2   : -//-
Ci   ldh  :leading dimension of H
Co Outputs
Co   H    :H is scaled by a q-dependent phase; see Remarks
Cl Local variables
Cr Remarks
Cr   when the hamiltonian is constructed by Ewald sums, and a basis
Cr   vector is shortened by dpos, the hamiltonin is scaled by
Cr   a q-dependent phase.
Cr   hamiltonian rows belonging to that site are scaled by
Cr      phase  = exp(2*pi*i * q * dpos)
Cr   hamiltonian columns belonging to that site are scaled by
Cr      phase  = exp(-2*pi*i * q * dpos)
Cb Bugs
Cb   This routine has only been checked for ib1=1 and ib2=nbas
Cb   This routine assumes hamiltonian is contiguous block betw/ (ib1..ib2)
Cb   This routine uses off from orbl and assigns ldimp=0.
Cb     These assume only lower block is being rotated
Cu Updates
Cu   17 Dec 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer opt,ldH,iprmb(1),ldima,ib1,ib2
      double precision q(3),tbas(3,ib2)
      double complex H(ldH,ldH)
C ... Local parameters
      logical nphase
      integer n0,nkap0,i,ib,j,j1,j2,ndimi,norb,off,ipr,PRDBG
      double precision twopi,sp
      double complex phase
      parameter (n0=10,nkap0=3,PRDBG=80)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)

      if (q(1) .eq. 0 .and. q(2) .eq. 0 .and. q(3) .eq. 0) return

      twopi = 8d0*datan(1d0)
      nphase = .false.
      call getpr(ipr)
      if (mod(opt/10,10) .ne. 0) nphase = .true.

C ... Offsets to starting and ending regions of h(ib1..ib2)
      call orbl(ib1,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
      j1 = off+1
      call orbl(ib2,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)
      j2 = off + ndimi

C --- For each site ibas, construct poke rotation matrix rotating Ylm(ibas)
      do  ib = ib1, ib2

C   ... Row, column offsets to h
        call orbl(ib,0,ldima,iprmb,norb,ltab,ktab,off,offl,ndimi)

C   ... Phase
        sp = twopi*(tbas(1,ib)*q(1) + tbas(2,ib)*q(2) + tbas(3,ib)*q(3))
        if (nphase) sp = -sp
        phase = dcmplx(dcos(sp),dsin(sp))
        if (ipr .ge. PRDBG) then
          print 333, ib,sp*360/twopi,q
  333     format(' rothph: ib=',i4,'  phase=',f7.1,' deg.  q=',3f7.1)
        endif

C   --- For all orbitals belonging to site ib, h->phase*h and/or h->h*phase+
        if (mod(opt,10) .eq. 0 .or. mod(opt,10) .eq. 2) then
          do  j = j1, j2
            do  i = off+1, off+ndimi
              h(i,j) = phase*h(i,j)
            enddo
          enddo
        endif

        if (mod(opt,10) .eq. 1 .or. mod(opt,10) .eq. 2) then
          phase = dconjg(phase)
          do  i = off+1, off+ndimi
            do  j = j1, j2
              h(j,i) = h(j,i)*phase
            enddo
          enddo
        endif

      enddo

      end

