#define USEOPTIONAL
      module m_sigflg
      contains
      subroutine sigflg(dc,sio,strn,ledit,lwsig,lrot,lfbz,lnwmsh,lonesp,
     .qpoff,rots)

C- Reads flags for Sigma IO
C ----------------------------------------------------------------------
Ci Inputs
Ci   dc:   separator for switches
Ci  sio:   a string, which must be one of :
Ci         'rsig' --- parse parameters for reading sigma:
Ci                lfbz qpoff
Ci         'wsig' --- parse parameters for writing a new sigma:
Ci                ledit lwsig lrot lfbz lnwmsh lonesp qpoff rots
Ci   strn: string containing switches
Co Outputs
Co   ledit :(wsig only)
Co         :1  flag requesting sigma editor
Co   lwsig :(wsig only)
Co         :0 normal exit
Co         :1  transformation of sigm from orbital to LDA basis
Co         :-1 transformation of from LDA basis to orbital basis
Co         :2  Same as 1, but high,low parts of sigm replaced by diagonal
Co         :3  Write evals,evecs of LDA hamiltonian to file
Co         :4  Write evals,evecs of total hamiltonian to file
Co         : OLD; not used
Co         :4  Make QP hamiltonian from Z^LDA and U with E(qp)
Co         :5  Make orthogonalized QP hamiltonian from Z^LDA and U with E(qp)
Co   lrot  :(wsig only)
Co         :nonzero, flag to rotate sigma matrix
Co   lfbz  :(rdsig and wsig)
Co         :T, write sigma for qp in entire BZ
Co   lnwmsh:(wsig only)
Co         :T, write sigm for qp on new mesh
Co   lonesp:(wsig only)
Co         :nonzero, average up and down spins before writing
Co   qpoff :(rdsig and wsig)
Co         :k-point offset to be added to uniform mesh
Co   rots  :(wsig only)
Co         :string containing rotation matrix for sigma
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   21 Jun 08
Cu   15 Jan 07  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character dc*1,sio*(4)
#ifdef USEOPTIONAL
      character strn*(*)
      character,optional::rots*(*)
      integer,optional:: ledit,lwsig,lrot
      integer::lonesp
      logical lfbz
      logical,optional::lnwmsh
#else
      character strn*(*),rots*(*)
      integer ledit,lwsig,lrot,lonesp
      logical lfbz,lnwmsh
#endif
      double precision qpoff(3)
C ... Local parameters
      integer parg
      integer i,j,j1,j2,iv(10)
      logical lio

      logical:: isanrg, l_dummy_isanrg

      lio = sio .eq. 'wsig'
      if (.not. lio) then
        if (sio .ne. 'rsig') call rx('sigflg: bad argument '//sio)
      endif
      if (lio) then
        ledit = 0
        lwsig = 0
        lrot = 0
        lfbz = .false.
        lnwmsh = .false.
        lonesp = 0
        rots = ' '
        call dpzero(qpoff,3)
      else
        lfbz = .false.
        lonesp = 0
        call dpzero(qpoff,3)
      endif

      if (dc .ne. ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 0
   50   continue
        j2 = j2+1
        if (strn(j2:j2) .eq. dc) goto 50
        j1 = min(len(strn),j2)
        call nwordg(strn,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then
          elseif (strn(j1:j1+3) .eq. 'edit'.and. lio) then
            ledit = 1
          elseif (strn(j1:j1+4) .eq. 'onesp') then
            lonesp = 1
          elseif ((strn(j1:j1+3) .eq. 'rot=' .or.
     .    strn(j1:j1+3) .eq. 'ROT=') .and. lio) then
            rots = strn(j1+4:j2)
          elseif (strn(j1:j1+5) .eq. 'trans=' .and. lio) then
            j = 0
            if (parg('trans=',2,strn(j1:),j,len(strn(j1:)),
     .      dc//' ',1,1,i,lwsig) .lt. 0) goto 999
Ckino isanrg is logical function,             call isanrg(mod(lwsig,10),-1,5,'SIGFLG:','trans',.true.)
            l_dummy_isanrg=isanrg(mod(lwsig,10),-1,5,'SIGFLG:','trans',.true.)
            if (lwsig .ne. 0 .and. mod(lwsig,10) .eq. 0 .or.
     .      lwsig .ge. 40)
     .      call rxi('rdsigm: nonsensical trans',lwsig)
          elseif (strn(j1:j1+4) .eq. 'phase'.and. lio) then
            lrot = 2
          elseif (strn(j1:j1+2) .eq. 'fbz') then
            lfbz = .true.
          elseif (strn(j1:j1+4) .eq. 'newkp'.and. lio) then
            lnwmsh = .true.
          elseif (strn(j1:j1+5) .eq. 'shftq=') then
            j = 0
            if (parg('shftq=',4,strn(j1:),j,len(strn(j1:)),
     .      ', '//dc,2,3,iv,qpoff) .lt. 0) goto 999
          elseif (strn(j1:j1+4) .eq. 'shftq') then
            qpoff(1) = 4d0/4096
            qpoff(2) = 5d0/4096
            qpoff(3) = 7d0/4096
          else
            goto 59
          endif
          goto 50
   59     continue
  999     call rxs('rdsigm: failed to parse wsig option', strn)
        endif
      endif

      end subroutine sigflg
      end module m_sigflg


      subroutine iosigh(optio,mode,nsp,ndimh,nk1,nk2,nk3,nqp,
     .lshft1,lshft2,lshft3,ifi)
C- Read/write header information in sigma file
C ----------------------------------------------------------------------
Ci Inputs
Ci   optio :(read only) 1s digit
Ci         :0 read header data into
Ci            nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
Ci         :1 File data must passed matched value for
Ci            nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
Ci         :2 file mode must passed matched value
Ci         :These digits may be taken in any combination
Ci   ifi   :File logical unit: ifi>0 for read; ifi<0 for write
Cio Inputs/Outputs
Cio  mode  : switches with information about the contents of sigma
Cio        :1s digit:
Cio        :0 file contains sigma in the orbital basis
Cio        :1 file contains sigma in the lda basis, and LDA
Cio           eigenvalues
Cio        :2 file contains sigma in the lda basis, and LDA
Cio           eigenvalues
Cio        :3 file consists of LDA eigenvectors and eigenvalues
Cio        :4 file consists of LDA eigenvectors and eigenvalues
Cio        :5 file contains sigma in the orbital basis, and LDA
Cio           eigenvalues
Cio  nsp   :2 for spin-polarized case, otherwise 1
Cio  ndimh :dimension of hamiltonian
Cio  nk1..nk3 : number of divisions in 1st,2nd,3rd Qlat
Cio  nqp   :number of irr. k-points
Cio  lshft1:1 if offset BZ in Qlat(1), otherwise 0
Cio  lshft2:1 if offset BZ in Qlat(2), otherwise 0
Cio  lshft3:1 if offset BZ in Qlat(3), otherwise 0
Cl Local variables
Cr Remarks
Cr   Reads or writes header file for sigma, with sanity checks
Cr   sigm file contents: (1s digit mode = 0..3,5,6)
Cr    1.  header
Cr    2.  qp offset (version 2 and later)
Cr    ... 3 and 4 are repeated for nqp points (nqp is contained in header)
Cr    3.  qp
Cr    4.  sigm(qp) or ham(qp)
Cr   sigm file contents: (1s digit mode = 1 or 3)
Cr    1.  header
Cr    2.  qp offset (version 2 and later)
Cr    ... 3..5 are repeated for nqp points (nqp is contained in header)
Cr    3.  qp
Cr    4.  eval(LDA)
Cr    5.  sigm(LDA) (mode 1) or evec(LDA) (mode 3)
Cu Updates
Cu   27 Jun 08 Revised mode, and input checking.  Old optio=1 -> optio=3
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,ifi,optio,
     .mode
C ... Local parameters
      integer jfi,i1,i2,i3,i4,i5,i6,i7,i8,i9,ivsn,modef
      logical:: isanrg, l_dummy_isanrg
      if (ifi .lt. 0) then
        ivsn = 1
C       if (mode .eq. 0) ivsn = 0
        jfi = -ifi
        rewind jfi
C       Write 0 => file contains version number (versions 1 and later)
        if (ivsn .ne. 0) write(jfi) 0, ivsn, mode
        write(jfi) nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
      elseif (optio .eq. 0) then
        rewind ifi
C       Read version, checking whether header is pre- version 1
        read(ifi,err=99,end=99) i1,ivsn,mode
C       File header is pre-version 1; set version to zero
        if (i1 .ne. 0) then
          ivsn = 0
          mode = 0
          rewind ifi
        endif
        read(ifi,err=99,end=99)
     .  nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
      else
        rewind ifi
C       Read version, checking whether header is pre- version 1
        read(ifi,err=99,end=99) i1,ivsn,modef
C       File header is pre-version 1; set version, modef to zero
        if (i1 .ne. 0) then
          ivsn = 0
          modef = 0
          rewind ifi
        endif
C       Assign mode or check with file correspondence
        if (mod(optio/2,2) .ne. 0) then
Ckino isanrg is logical function,           call isanrg(modef,mode,mode,'iosigh:','file''s mode',.true.)
          l_dummy_isanrg=isanrg(modef,mode,mode,'iosigh:','file''s mode',.true.)
        else
          mode = modef
        endif
C       Assign header parameters or check with file correspondence
        if (mod(optio,2) .ne. 0) then
          read(ifi,err=99,end=99) i1,i2,i3,i4,i5,i6,i7,i8,i9
          l_dummy_isanrg=isanrg(i1,nsp,nsp,'iosigh:','file''s nsp',.true.)
          l_dummy_isanrg=isanrg(i2,ndimh,ndimh,'iosigh:','file''s ndimh',.true.)
Ckino isanrg is logical function,           call isanrg(i1,nsp,nsp,'iosigh:','file''s nsp',.true.)
Ckino isanrg is logical function,           call isanrg(i2,ndimh,ndimh,'iosigh:','file''s ndimh',.true.)
Ckino isanrg is logical function,           call isanrg(i3,nk1,nk1,'iosigh:','file''s nk1',.true.)
Ckino isanrg is logical function,           call isanrg(i4,nk2,nk2,'iosigh:','file''s nk2',.true.)
Ckino isanrg is logical function,           call isanrg(i5,nk3,nk3,'iosigh:','file''s nk3',.true.)
Ckino isanrg is logical function,           call isanrg(i6,nqp,nqp,'iosigh:','file''s nqp',.true.)
Ckino isanrg is logical function,           call isanrg(i7,lshft1,lshft1,'iosigh:','file''s lshft1',.true.)
Ckino isanrg is logical function,           call isanrg(i8,lshft2,lshft2,'iosigh:','file''s lshft2',.true.)
Ckino isanrg is logical function,           call isanrg(i9,lshft3,lshft3,'iosigh:','file''s lshft3',.true.)
          l_dummy_isanrg=isanrg(i3,nk1,nk1,'iosigh:','file''s nk1',.true.)
          l_dummy_isanrg=isanrg(i4,nk2,nk2,'iosigh:','file''s nk2',.true.)
          l_dummy_isanrg=isanrg(i5,nk3,nk3,'iosigh:','file''s nk3',.true.)
          l_dummy_isanrg=isanrg(i6,nqp,nqp,'iosigh:','file''s nqp',.true.)
          l_dummy_isanrg=isanrg(i7,lshft1,lshft1,'iosigh:','file''s lshft1',.true.)
          l_dummy_isanrg=isanrg(i8,lshft2,lshft2,'iosigh:','file''s lshft2',.true.)
          l_dummy_isanrg=isanrg(i9,lshft3,lshft3,'iosigh:','file''s lshft3',.true.)
        else
          read(ifi,err=99,end=99)
     .    nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3
        endif
      endif
      return
   99 continue
      call rx('iosigh: failed to read header from file')
      end subroutine iosigh


      subroutine chksgr(opt,ltrans,kcplx,plat,nsp,ndimh,ifis,hreal,
     .iprmb,nttabs,iaxs,hrs,ndhrs,tols,nqp,qplst,nbas,lrot,rotm,
     .delT)
      use m_globalvariables
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
Ci   hrs   :sigm(R)
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
      integer i,j,isp,nglob,stdo,iq1,nl,mxorb,ipr,lgunit,isw,opt0,opt01,
     .nscnts
      integer ib,jb,is,norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
C     these are read by sigma header
      integer nk1,nk2,nk3,lshft(3),offh(nbas),ndim(nbas),offi,offj,
     .ndi,ndj
C     integer os
      double precision qp(3),errmxx,errmx,err,qpr(3),xx
      double complex hwk(ndhrs,ndhrs)
      complex(8),allocatable:: sq(:,:),sq1(:,:),uz(:,:)
C ... Heap
      integer ::iwdummy

      integer ii,jj
      integer:: ndimhr !Read from sigm. This can be different from ndimh in the case of PMT, 
                       ! that is, ndimh< ndimhr=nlmto+napw !apr2012takao
      complex(8),allocatable:: sqr(:,:)

      call getpr(ipr)
      stdo = lgunit(1)
Changenglob      nl = nglob('nl')
      nl = globalvariables%nl
Changenglob      mxorb = nglob('mxorb')
      mxorb = globalvariables%mxorb
      opt01 = mod(mod(opt,10),2)
      opt0  = mod(opt,10)
      nscnts = 0

      allocate(sq(ndimh,ndimh),sq1(ndimh,ndimh))

C --- Rotate hrs if switch set ---
      if (opt01 .ne. 0 .and. mod(lrot,2) .ne. 0) then
        allocate(uz(ndimh,ndimh))
        call rothrm(2,0,ndimh,iprmb,rotm,1,nbas,ndimh,uz,uz)
C       call zprm('Uz',2,uz,ndimh,ndimh,ndimh)
C       Table of hamiltonian offsets and site dimensions
        do  ib = 1, nbas
          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,j,offl,ndim(ib))
          offh(ib) = j+1
        enddo

C       For each pair and spin, rotate hrs to U(ib) hrs U(jb)+
        do  is = 1, nttabs
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
C              call zprm('U S',2,hwk,ndhrs,ndj,ndi)
C             U sigm(unrotated) U+
              call zgemm('N','C',ndj,ndi,ndi,(1D0,0D0),hwk,ndhrs,
     .        Uz(offi,offi),ndimh,(0D0,0D0),hrs(1,1,isp,is),ndhrs)
C              call zprm('U S U+',2,hrs(1,1,isp,is),ndhrs,ndhrs,ndhrs)
            else
C             Make U1 h U2+
              call zgemm('N','N',ndi,ndj,ndi,(1D0,0D0),Uz(offi,offi),
     .        ndimh,hrs(1,1,isp,is),ndhrs,(0D0,0D0),hwk,ndhrs)
C             call zprm('U S',2,hwk,ndhrs,ndi,ndj)
C             U sigm(unrotated) U+
              call zgemm('N','C',ndi,ndj,ndj,(1D0,0D0),hwk,ndhrs,
     .        Uz(offj,offj),ndimh,(0D0,0D0),hrs(1,1,isp,is),ndhrs)
C             call zprm('U S U+',2,hrs(1,1,isp,is),ndhrs,ndhrs,ndhrs)
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
        errmxx = 0

C     If sigma file not spin polarized, use sigma from spin 1
        if (isp .eq. 2 .and. lsplts) then
          call iosigh(1,nscnts,1,ndimhr,nk1,nk2,nk3,nqp,
     .    lshft(1),lshft(2),lshft(3),ifis)
        endif

C ... For each q, compare Bloch sum of hrs against file's s(q)
C     or write to file
        do  iq1 = 1, nqp

C       print *, 'starting iq=',iq1

          if (opt01 .eq. 0) then
            read(ifis) qp
            if (opt0 .eq. 2) call dcopy(3,qplst(1,iq1),1,qp,1)
            if (nscnts .eq. 5) call dpdump(xx,1,ifis)
!! readin sq
c            call dpdump(sq,ndimh**2*2,ifis)
            allocate(sqr(ndimhr,ndimhr))
            if(ndimhr<ndimh) call rx('chksgr :ndimhr<ndimh')
            call dpdump(sqr,ndimhr**2*2,ifis)
            sq = sqr(1:ndimh,1:ndimh)
            deallocate(sqr)
!!
            call ztoyy(sq,ndimh,ndimh,ndimh,ndimh,kcplx,1)
C          print *, iq1
C          call zprm('sq',2,sq,ndimh,ndimh,ndimh)
          else
            call dcopy(3,qplst(1,iq1),1,qp,1)
          endif
C       Call bloch : 104010=perm orb, transpose, c*16
          i = 100000 + 4000 + 40*(1-hreal) + 10
          if (.not. ltrans) i = i - 4000
          call bloch ( i , qp , nl , plat , mxorb , iprmb , 1 , nttabs 
     .     , iaxs , hrs , ndhrs , isp , nsp , ndimh , ndimh , 0 , ndimh 
     .     , 0 , ndimh , 0 , sq1 , iwdummy , iwdummy )

C        print *, 'iq=',iq1,'qp=',sngl(qp)
C        call zprm('sigm',2,sq1,ndimh,ndimh,ndimh)
          errmx = 0

C       Compare Bloch sum of hrs against file's s(q)
          if (opt01 .eq. 0) then
            ii = 1
            jj = 1
            do  j = 1, ndimh
              do  i = 1, ndimh
                err = abs(sq(i,j)-sq1(i,j))
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

C   ... Write to s(q)
          else

C     ...  Rotation of sigma(k) (no longer needed, because hrs(T) was rotated)
C          if (mod(lrot,2) .ne. 0) then
C
CC           call zprm('unrotated S',2,sq1,ndimh,ndimh,ndimh)
C
CC           U sigm(unrotated)
C            call zgemm('N','N',ndimh,ndimh,ndimh,(1D0,0D0),
C     .        Uz,ndimh,sq1,ndimh,(0D0,0D0),sq,ndimh)
CC           call zprm('U S',2,sq,ndimh,ndimh,ndimh)
CC           U sigm(unrotated) U+
C            call zgemm('N','C',ndimh,ndimh,ndimh,(1D0,0D0),
C     .        sq,ndimh,Uz,ndimh,(0D0,0D0),sq1,ndimh)
C            call zprm('U S U+',2,sq1,ndimh,ndimh,ndimh)
C
C          endif

C         Additional scaling by phase shifts
            if (lrot .ge. 2) then
C           call zprm('unscaled S',2,sq1,ndimh,ndimh,ndimh)
              call dgemm('T','N',3,1,3,1d0,plat,3,qp,3,0d0,qpr,3)
C           call prmx('q in units of qlat',qpr,3,3,1)
              call rothph(02,qpr,delT,ndimh,iprmb,1,nbas,ndimh,sq1)
C           call zprm('scale S by phase',2,sq1,ndimh,ndimh,ndimh)
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
      if (opt01 .ne. 0 .and. mod(lrot,2) .ne. 0) then
        deallocate(uz)
      endif

      end subroutine chksgr


      subroutine siged(mode,nbas,nsp,ndhrs,plat,pos,ndimh,iprmb,hreal,
     .ntab,iax,hrs)
      use m_globalvariables
C- Sigma editor
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :compound of these digits:
Ci         :1s digit 1 : merge up and down spins, and exit
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ndhrs :dimensions hrs
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   ndimh :hamiltonian dimension
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   hreal :1 if hrs is real
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Cio Inputs/Outputs
Cio  hrs   :real-space hamiltonian, which may be modified on exit
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   01 May 06 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,nsp,ndhrs,ndimh,iprmb(1),hreal
      integer n0,nkap0,n0H,niax
      double precision plat(3,3),pos(3,nbas)
      parameter (n0=10,nkap0=3,n0H=5,niax=10)
      integer iax(niax,1),ntab(*)
      double precision hrs(2,ndhrs,ndhrs,nsp,1)
C ... Local parameters
      character outs*150, strn*80
      integer nw,ix(10),iwk(10),it1,it2,a2vec,i,j,stdo,nglob
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .blks(n0*nkap0),ntabh(n0*nkap0)
      integer it,ib,jb,iorb,i1,l1,ik1,norb,j1,j2,isp,isub(4)
      double precision xx(10),rr,drr2,vlat(3),rcut,dglob
      logical lin

C ... Default range: all pairs
      it1 = 1
      it2 = ntab(nbas+1)
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      rcut = 2
      if (hreal .ne. 0) call rx('siged not ready for hreal')

C ... Merge up, down spin potentials into average
C     print *, dsum(2*ndhrs*ndhrs*nsp*it2,hrs,1)
      if (mod(mode,10) .eq. 1) then
        do  i = it1, it2
C         print *, i
C         call zprm('h1,i',2,hrs(1,1,1,1,i),ndhrs,ndhrs,ndhrs)
C         call zprm('h2,i',2,hrs(1,1,1,2,i),ndhrs,ndhrs,ndhrs)
          call daxpy(2*ndhrs*ndhrs,1d0,hrs(1,1,1,2,i),1,
     .    hrs(1,1,1,1,i),1)
        enddo
        do  i = it1, it2
          call dcopy(2*ndhrs*ndhrs,hrs(1,1,1,1,i),1,hrs(1,1,1,i,1),1)
C         call zprm('hi',2,hrs(1,1,1,i,1),ndhrs,ndhrs,ndhrs)
        enddo
C       print *, dsum(2*ndhrs*ndhrs*it2,hrs,1)

        call dscal(2*ndhrs*ndhrs*it2,0.5d0,hrs,1)
C       print *, dsum(2*ndhrs*ndhrs*it2,hrs,1)
        nsp = 1
Changenglob        rr = dglob('nsp',1d0,1)
        globalvariables%nsp = 1d0; globalvariables%l_nsp = globalvariables%l_nsp +1; rr = 1d0
        return
      endif

C      print *, '!!'
C      call dvset(hrs,1,2*ndhrs*ndhrs*nsp*it2,-1d0)

   99 print 301
  301 format(//
     .' Welcome to the sigma file editor.')
  306 format(' Failed to parse string ',a,' ... try again.')
  100 continue
#if PRTNOCR
      print '('' Option (enter ? for menu) : '',$)'
#else
      print '('' Option (enter ? for menu) : '')'
#endif
      outs = ' '
      read(*,'(a150)') outs
      call locase(outs)
      if (.false.) then
      elseif (outs .eq. ' ') then
        print '('' Enter ''''q'''' to exit, ''''a'''' to abort'')'
        goto 100

C ... Specify pair
      elseif (outs(1:4) .eq. 'pair') then
        call words(outs,nw)
        if (nw .eq. 1) then
          it1 = 1
          it2 = ntab(nbas+1)
        else
          i = 5
          nw = a2vec(outs,len(outs),i,2,', ',2,-3,2,ix,iwk)
          if (nw .le. 0) then
            print 306, outs(4:)
            goto 100
          else if (nw .eq. 1) then
            iwk(2) = iwk(1)
          endif
          it1 = min(iwk(1),iwk(2))
          it2 = max(iwk(1),iwk(2))
        endif

        ib = iax(1,it1)
        jb = iax(2,it1)

        if (it1 .eq. it2) then
          print '('' edit pair'',i6,''   ib ='',i4,''   jb ='',i4)',
     .    it1,ib,jb

          ib = iax(1,it1)
          isub(1) = 1
          call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(2))
          jb = iax(2,it1)
          isub(3) = 1
          call orbl(jb,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(4))
        else
          print '('' edit range'',i5,''  to'',i5)', it1,it2
        endif
        goto 100

C ... Summarize hamiltonian for this pair
      elseif (outs(1:3) .eq. 'sum' .and. it1 .eq. it2) then

        ib = iax(1,it1)
  125   continue
        call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,xx)
        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntabh,blks)
        write(stdo,359) ib
  359   format(/
     .  ' Summary of hamiltonian subblocks for site',i4/
     .  '     kap  lstart  h(k)    h(rs)     size'/
     .  '                 start    range')

        do  iorb = 1, norb
          if (blks(iorb) .ne. 0) then
C         l1,ik1 = l and kaph indices, needed for sigma
            l1  = ltab(iorb)
            ik1 = ktab(iorb)
C         i1 = orbital index in iprmb order; ilm1 = augm. index in L order
            i1 = offl(iorb)
            outs = ' '
            call awrit2('%i,%i',outs,len(outs),0,1+i1-offl(1),
     .      i1-offl(1)+blks(iorb))
            print '(3i7,5x,a10,i3)', ik1, l1, 1+i1, outs(1:10), blks(iorb)
          endif
        enddo
        if (iax(2,it1) .ne. ib) then
          ib = iax(2,it1)
          goto 125
        endif
        print *, ' '

        goto 100

C ... Summarize pair table
      elseif (outs(1:3) .eq. 'sum' ) then

        write(stdo,356)
  356   format(/
     .  ' Summary of pair table for R.S. sigma'/
     .  ' site  start   end  size')
        do  ib = 1, nbas
          write(stdo,357) ib,ntab(ib)+1,ntab(ib+1),ntab(ib+1)-ntab(ib)
  357     format(i4,2i7,i5)
        enddo

        goto 100

C ... Get subblock of hrs
      elseif (outs(1:3) .eq. 'sub') then

        ib = iax(1,it1)
        isub(1) = 1
        call orbl(ib,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(2))
        jb = iax(2,it1)
        isub(3) = 1
        call orbl(jb,0,ndimh,iprmb,norb,ltab,ktab,xx,offl,isub(4))
        call words(outs,nw)
        if (nw .gt. 1) then
          call word(outs,2,j1,j2)
          i = j1-1
          nw = a2vec(outs,len(outs),i,2,', ',2,-3,4,ix,iwk)
          if (nw .ne. 4) then
            print 306, outs(j1:)
          else
            call icopy(4,iwk,1,isub,1)
          endif
          goto 100
        endif

C ... Mask subblock of hrs
      elseif (outs(1:4) .eq. 'mask' .and. it1 .eq. it2) then

        call words(outs,nw)
        if (nw .eq. 1) then
          print *, 'no subblocks specified ... nothing to mask'
          goto 100
        endif

        do  i = isub(1), isub(2)
          do  j = isub(3), isub(4)

C         See whether (i,j) is inside any mask
            lin = .false.
            do  i1 = 2, nw
              call word(outs,i1,j1,j2)
              l1 = j1-1
              j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,4,ix,iwk)
              if (j1 .ne. 4) then
                print 306, outs(i1:)
                goto 100
              endif
              lin = lin .or. (i .ge. iwk(1) .and. i .le. iwk(2) .and.
     .        j .ge. iwk(3) .and. j .le. iwk(4))
            enddo

C         If (i,j) is outside all masks, zero it out
            if (.not. lin) then
              hrs(1,i,j,1,it1) = 0
              hrs(2,i,j,1,it1) = 0
              hrs(1,i,j,nsp,it1) = 0
              hrs(2,i,j,nsp,it1) = 0
            endif


          enddo
        enddo

        goto 100

C ... Mask hrs not in range
      elseif (outs(1:4) .eq. 'mask') then

        call words(outs,nw)
        if (nw .eq. 1) then
          print *, 'no ranges specified ... nothing to mask'
          goto 100
        endif

        do  it = it1, it2

C         See whether (i,j) is inside any mask
          lin = .false.
          do  i1 = 2, nw
            call word(outs,i1,j1,j2)
            l1 = j1-1
            j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,2,ix,iwk)
            if (j1 .ne. 2) then
              print 306, outs(i1:)
              goto 100
            endif
            lin = lin .or. (it .ge. iwk(1) .and. it .le. iwk(2))
          enddo

C         If (it) is outside all masks, zero it out
          if (.not. lin) then
            call dpzero(hrs(1,1,1,1,it),2*ndhrs*ndhrs*nsp)
          endif

        enddo

        goto 100

C ... Show sigma for range of pairs
      elseif (outs(1:5) .eq. 'shows' .or. outs(1:6) .eq. 'showsd') then

        call words(outs,nw)

        do  i = 1, max(nw-1,1)
          if (nw .eq. 1) then
            call icopy(4,isub,1,iwk,1)
          else
            call word(outs,i+1,j1,j2)
            l1 = j1-1
            if (outs(1:6) .eq. 'showsd') then
              j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,2,ix,iwk)
              if (j1 .ne. 2) then
                print 306, outs(l1:)
                goto 100
              endif
            else
              j1 = a2vec(outs,len(outs),l1,2,', ',2,-3,4,ix,iwk)
              if (j1 .ne. 4) then
                print 306, outs(l1:)
                goto 100
              endif
            endif
          endif

          do  it = it1, it2
            ib = iax(1,it)
            jb = iax(2,it)
            rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .      iax(3,it),iax(4,it),iax(5,it),vlat))
            write(stdo,332)
            write(stdo,334) iax(1,it),iax(2,it),
     .      (vlat(j),j=1,3), rr, (iax(j,it), j=3,5), it

            do  isp = 1, nsp
              if (outs(1:6) .eq. 'showsd') then
                write(stdo,'(19f12.6)')
     .          (hrs(1,j1,j1,isp,it), j1=iwk(1),iwk(2))
                write(stdo,'(19f12.6)')
     .          (hrs(2,j1,j1,isp,it), j1=iwk(1),iwk(2))
              else
                call awrit5('rows (%i,%i) cols (%i,%i) spin %i',strn,80,0,
     .          iwk(1),iwk(2),iwk(3),iwk(4),isp)
                call ywrm(0,strn,3,stdo,'(19f12.6)',
     .          hrs(1,iwk(1),iwk(3),isp,it),1,ndhrs,
     .          iwk(2)-iwk(1)+1,iwk(4)-iwk(3)+1)
              endif
            enddo
          enddo
        enddo

        goto 100

      elseif (outs(1:5) .eq. 'showp' ) then

        i = 0
        do  it = it1, it2
          ib = iax(1,it)
          jb = iax(2,it)
          rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .    iax(3,it),iax(4,it),iax(5,it),vlat))
          if (ib .ne. i) then
            write(stdo,332)
  332       format(/'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,
     .      'd       -x-plat-   pair')
C            write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut
C  345       format(' pairc, ib=',i3,':',i4,' neighbors in range',f7.3)
          endif
          i = ib
          if (rr .lt. rcut) then
            write(stdo,334) iax(1,it),iax(2,it),
     .      (vlat(j),j=1,3), rr, (iax(j,it), j=3,5), it
  334       format(i4,i4,3f11.6,f9.4,3x,3i3,i6)
          endif
        enddo

        goto 100

C ... abort
      elseif (outs(1:1) .eq. 'a') then
        call rx0('aborting sigma editor ... no file written')

C ... quit
      elseif (outs(1:1) .eq. 'q') then
        goto 199

C ... help
      elseif (outs .eq. '?') then
        print 310
  310   format(/
     .  ' Editor works with a range of pairs,',
     .  ' or with a particular pair.'/
     .  ' In the latter case, you can edit a subblock',
     .  ' of h(rs) for that pair'//
     .  ' Select one of these options:'/
C    .  t4,'i',t10,'to import a portion of a sigma file'/
     .  t4,'pair',t10, '[#1 [#2]] specifies a pair or range of pairs',
     .  ' used by other options.'/
     .  ,t10, 'Use only #1 to select a specific pair.'/
     .  ,t10, 'Use no arguments to select all pairs.'/
     .  t4,'sum', t10, 'summarizes : pair table,',
     .  ' or sigma for given pair')
        print 311
  311   format(
     .  t4,'sub',  t10, '              work with entire block of hrs'/
     .  t4,'sub',  t10, 'r1,r2,c1,c2 : define subblock of hrs'/
     .  t10, 'other options use.'/
     .  t4,'mask', t10, 'range1 range2 ... zero hrs for pairs not in'/
     .  t10, 'one of supplied ranges'/
     .  t10, 'range* has structure #1,#2'/
     .  t10, 'Applies to range of pairs specified',
     .  'with "pair"'/
     .  t4,'mask', t10, 'subbl1 subbl2 ... zero hrs(current pair) for '/
     .  t10, 'any element not in one of supplied subblocks'/
     .  t10, 'subbl* has structure r1,r2,c1,c2'/
     .  t10, 'Applies to one pair only'/
     .  t4,'showp',t10, 'shows pairs within table'/
     .  t4,'shows' t10, '[r1,r2,c1,c2] shows subblock of sigma, for '/
     .  t10, 'list of pairs specified by ''pair'' command'/
     .  t4,'edit', t10, '(does nothing), or matrix element of sigma'/
     .  t4,'q',t10,'to quit the editor'/
     .  t4,'a',t10,'to abort (no sigma file is written)'/
     .  )
        goto 100

      else
        print '(1x,''unrecognized option: '',a10)', outs
        goto 100

      endif

C ... Cleanup and exit rst editor
  199 continue
      end subroutine siged




C      subroutine sigfn(io,lwsig,ifi)
CC- Opens file and returns file logical unit for sigm transformation modes
CC ----------------------------------------------------------------------
CCi Inputs
CCi   io    :0 => open input file ifi
CCi         :1 => open output file ifi
CCi   lwsig
CCo Outputs
CCo   ifi   :file logical unit
CCl Local variables
CCl         :
CCr Remarks
CCr
CCu Updates
CCu   20 Dec 07
CC ----------------------------------------------------------------------
C      implicit none
C      integer io,lwsig,ifi
C      integer fopna
C
C      if (lwsig .eq. 1 .or. lwsig .eq. 2) then
C        if (io .eq. 0) ifi = fopna('sigm',-1,4)
C        if (io .eq. 1) ifi = fopna('sigm2',-1,4)
C      elseif (lwsig .eq. -1) then
C        if (io .eq. 0) ifi = fopna('sigm2',-1,4)
C        if (io .eq. 1) ifi = fopna('sigm3',-1,4)
C      elseif (lwsig .eq. 3) then
C        if (io .eq. 0) ifi = fopna('sigm3',-1,4)
C        if (io .eq. 1) ifi = fopna('sigm2',-1,4)
C      else
C        ifi = fopna('sigm',-1,4)
C      endif
C
C      end




C ... titus
C      subroutine getoffH(offHo,offH,n0H,nkap0,nbas)
C      implicit none
C      integer n0H,nkap0,nbas
C      integer offHo(n0H,nkap0,nbas+1),offH(n0H,nkap0,nbas+1)
C      integer i,j,k
C      do 11 k = 1, (nbas+1)
C      do 11 j = 1, nkap0
C      do 11 i = 1, n0H
C
C        offH(i,j,k)= offHo(i,j,k)
C 11   continue
C
C      end
C ... titus

C      subroutine pvtrod(sfz,k1,k2,k3,ldim)
C      integer k1,k2,k3,ldim
C      double complex sfz(k1,k2,k3,ldim,ldim)
C      integer i1,i2
C
C      i1 = 38
C      i2 = 38
C      call zprm3('gf(i1,i2)',0,sfz(1,1,1,i1,i2),k1,k2,k3)
C      return
C
C      call zprm3('gf(2,1)',0,sfz(1,1,1,2,1),k1,k2,k3)
C      call zprm3('gf(1,2)',0,sfz(1,1,1,1,2),k1,k2,k3)
C      end


ctakaox
      subroutine divctakao(a,b,n)

      integer:: n,i,j
      complex(8):: a(n,n),b(n,n)
      print *,' sum check 1=',sum(abs(a))
      print *,' sum check 1=',sum(abs(b))
      do i=1,n
        do j=1,n
c        write(6,"(2i5,2x, 4d13.4,3x,2d13.5)")
c     &  i,j,atan(dimag(a(i,j)/dreal(a(i,j))))       if( abs(a(i,j))>1d-3 )
c     &   write(6,"(2i5,2x, 4d13.4,3x,2d13.5)")i,j,a(i,j)
c     &   ,atan(dimag(a(i,j)/dreal(a(i,j))))
c,b(i,j),
c     &  a(i,j)/( b(i,j)+b(i,j)*0.01/abs(b(i,j)) )
c        a(i,j) = a(i,j)/ ( b(i,j)+b(i,j)*0.01/abs(b(i,j)) )
c        a(i,j)=a(i,j)/b(i,j)
        enddo
      enddo
      print *,' sum check 2=',sum(abs(a))
      end subroutine divctakao

ckino:        subroutine rdsigm(mode,nbas,nsp,ndimh,slat,sham,sbz,ifis,rsrnge,
ckino:       .lwsig)
ckino:
ckino:        use m_struc_def  !Cgetarg
ckino:        use m_hft2rs
ckino:        use m_globalvariables
ckino:        use m_sigflg
ckino:  C- Read sigm(k) from file and generate sigm(R)
ckino:  C ----------------------------------------------------------------------
ckino:  Ci Inputs
ckino:  Ci   mode  :1s digit
ckino:  Ci         :0 do not read self-energy; just exit
ckino:  Ci         :1 read and FT sigma(k) to make sigma(T)
ckino:  Ci         :2 symmetrize sigma(T)
ckino:  Ci         :4 take the real part of sigma(T) only
ckino:  Ci         :10s digit (used here only for printout)
ckino:  Ci         :0 Simple bloch sum of sigma
ckino:  Ci         :1 approx high- and low-energy sigma with diagonal
ckino:  Ci         :2 perturbation approach, diagonalizing exactly only
ckino:  Ci         :  the first kappa block.
ckino:  Ci         :3 Linear interpolation of sigma from neighboring p
ckino:  Ci         :  sig(q) = sum_i wt_i A+_i sig(qpi) A_i
ckino:  Ci         :  where A_i = z^LDA_qi (z_q^LDA)^-1
ckino:  Ci         :100s digit specifies number of interpolation points
ckino:  Ci         :    for interpolation mode 3.
ckino:  Ci         :  0 -> use default (4 points)
ckino:  Ci         :1000s digit specifies ASA conventions:
ckino:  Ci         :      poke sigm(k)_RL,R'L' -> sigm(T)_RL,R'L'
ckino:  Ci         :      (10000s digit of hft2rs)
ckino:  Ci         :10000s digit
ckino:  Ci         :  1 specifies file sigm stored with no group operations
ckino:  Ci         :    sigm is stored at k-points in the full BZ.
ckino:  Ci         :    Equivalent to --rsig:fbz
ckino:  Ci         :  2 do not force R.S. neighbor table to conform to
ckino:  Ci         :    symmetry
ckino:  Ci         :  4 do not force file qp to match those generated by
ckino:  Ci              bzmesh.  Use qp from bzmesh.
ckino:  Ci   nbas  :size of basis
ckino:  Ci   nsp   :2 for spin-polarized case, otherwise 1
ckino:  Ci   ndimh :hamiltonian dimension
ckino:  Ci   slat  :struct for lattice information; see routine ulat
ckino:  Ci     Elts read: plat opos nsgrp oistab osymgr oag
ckino:  Ci     Stored:
ckino:  Ci     Passed to:
ckino:  Ci   sbz   :struct for the Brillouin Zone; see routine ubz
ckino:  Ci     Elts read: nkabc lshft nkp oqp
ckino:  Ci         : (used only when rdsigm writes sigma on a new mesh)
ckino:  Ci   ifis  :file logical unit for sigma
ckino:  Ci   rsrnge:maximum length for connecting vectors in sigma(R)
ckino:  Cio Inputs/Outputs:
ckino:  Cio  sham  :struct for parameters defining hamiltonian; see routine uham
ckino:  Cio    Elts read: ldham ooffH oindxo sigp rsstol
ckino:  Cio    Stored:    nqsig oqsig ndhrs ontabs oiaxs ohrs
ckino:  Cio Outputs:
ckino:  Co   lwsig :0  No special transformations of sigma
ckino:  Co         :1  Mode transforms sigma from orbital into LDA basis.
ckino:  Co         :   rdsigm reads sigma from file sigm in orbital basis;
ckino:  Co         :   Calling program generates and stores transformed sigma
ckino:  Co         :2  Similar to lwsig=1, except
ckino:  Co             low- and high- energy blocks replaced by diagonal parts
ckino:  Co         :-1 Mode transforms sigm from LDA to orbital basis
ckino:  Co         :   (Inverse operation of lwsig=1 or lwsig=2).
ckino:  Co         :   This mode requires both sigma and LDA eigenvectors go be
ckino:  Co         :   stored on disk in files 'sigm' and 'evec'
ckino:  Co         :   rdsigm generates and stores transformation in file 'sigm2'
ckino:  Co         :3  Returns lwsig=3 to flag calling program.  It should generate and
ckino:  Co         :   store LDA eigenvalues and eigenvectors. No sigm file is read.
ckino:  Co         :4  Returns lwsig=4 to flag calling program.  It should generate and
ckino:  Co         :   store eigenvalues and eigenvectors.
ckino:  Co         :5  Writes sigm as in lwsig=2 with low- and high- energy blocks
ckino:  Co         :   replaced by diagonal parts, but writes sigm in orbital bas.
ckino:  Cl Local variables
ckino:  Cl   sigp  :parameters for approximating self-energy sigma.  sigma
ckino:  Cl         :is approximated by its diagonal part sigii for energies
ckino:  Cl         :below a low-energy cutoff (specified nmin or emin) and
ckino:  Cl         :above a low-energy cutoff (specified nmax or emax).
ckino:  Cl         : arg 1: specifies how to set diagonal part sigii
ckino:  Cl         :        for states above the high-energy cutoff nmax or emax
ckino:  Cl         :        0 constrain sigii to be > asig+bsig*e
ckino:  Cl         :        1 constrain sigii to be = asig+bsig*e
ckino:  Cl         :        2 constrain sigii to be > asig and < bsig
ckino:  Cl         :        3 constraint same as case 1.
ckino:  Cl         :          arg1=3 differs in that the least-squares fit to
ckino:  Cl         :          sigii (for informational purposes only, to help
ckino:  Cl         :          estimate asig and bsig) is done for states between
ckino:  Cl         :          efit and nmax or emax
ckino:  Cl         : arg 2: nmin : sigma for states 1..nmin are approximated by sigii
ckino:  Cl         : arg 3: emin : (used only if nmin<0)
ckino:  Cl         :             : sigma for levels e<emin are approximated by sigii
ckino:  Cl         : arg 4: nmax : sigma for levels i>nmax are approximated by
ckino:  Cl         :               sigii AND constrained according to mode
ckino:  Cl         : arg 5: emax : (used only if nmax<=0)
ckino:  Cl         :             : sigma for levels e<emax are approximated by
ckino:  Cl         :               sigii AND constrained according to mode
ckino:  Cl         : arg 6: asig : constraint used to approximate
ckino:  Cl         :               sigii = asig + E * bsig  or
ckino:  Cl         :               asig < sigii < bsig
ckino:  Cl         : arg 7: bsig : constraint used to approximate
ckino:  Cl         :               sigii = asig + E * bsig  or
ckino:  Cl         :               asig < sigii < bsig
ckino:  Cl         : arg 8: efit : (mode 3) energy minimium
ckino:  Cl                         for fitting asig and bsig (not used here)
ckino:  Cl   Note  :sigp takes a similar, but slightly different meaning when
ckino:  Cl         :used in conjuction with sigm(LDA)->sigm(orbital), i.e. lwsig=-1
ckino:  Cl         :In this case, sigp applies to the new basis, whose
ckino:  Cl         :eigenvalues and eigenvectors are stored in 'evec.'
ckino:  Cl         :Arguments are used in rotevs, with the following relations:
ckino:  Cl         : arg1 (mode) has same meaning as before
ckino:  Cl         : arg2 (nmin) takes meaning of lcore.
ckino:  Cl         :       nmin>0 => sigm does not contain this block,  but the
ckino:  Cl         :                 new basis does.  Use in emin (arg3) in
ckino:  Cl         :                 place of sigm for this block .
ckino:  Cl         :       nmin<0 => new basis does not contain this block
ckino:  Cl         : arg3 (emin) If nmin>0, use emin for diag. sigma, this block
ckino:  Cl         : arg4, arg5 (nmax,emax) could play the role of lhigh in
ckino:  Cl         :       rotevs.  However, they are not used.  lhigh is
ckino:  Cl         :       determined from constraint ndims+lccore+lhigh = ndimz
ckino:  Cl
ckino:  Cl   hreal :0 allow sigm(R) to be complex
ckino:  Cl         :1 assume sigm(R) is real
ckino:  Cl  ltrans :specifies indexing of sigm(T)_RL,R'L'L
ckino:  Cl         : ltrans = F taken from sigm(k)_RL,R'L'
ckino:  Cl         : ltrans = T taken from sigm(k)_RL,R'L'
ckino:  Cl  kcplx  :0 sigm(k) has real, imaginary separated
ckino:  Cl         :1 sigm(k) is in complex*16 format:
ckino:  Cl         :2 sigm(k) has real, imaginary separated by columns
ckino:  Cl  lphase :T if phase convention phi = q * [(g R_j + a) - R_i]
ckino:  Cl         :  for rotations should be scaled by -1
ckino:  Cl  lssym  :10000s digit mode
ckino:  Cl  lfbzin :flags whether input self-energy file has suppressed symops
ckino:  Cl  nscnts :compound of switches containing file contents of sigma
ckino:  Cl  qoffi  :k-mesh offset for input sigma file
ckino:  Cl  qoffo  :k-mesh offset for generated sigma file
ckino:  Cr Remarks
ckino:  Cb Bugs
ckino:  Cb   Routine should be able to accomodate offset qp case (ifac<>0)
ckino:  Cb   Routine should be able to accomodate case when a different
ckino:  Cb   list of irreducible qp is used than the one generated by
ckino:  Cb   bzmesh.
ckino:  Cu Updates
ckino:  Cu   27 Jun 08 Redesigned transformation of sigma to new basis
ckino:  Cu   15 May 07 Parallelized symmetrizer
ckino:  Cu   24 Jan 07 New option --wsig:onesp
ckino:  Cu   20 Jan 07 Bug fix, --mixsig, spin polarized case
ckino:  Cu   16 Jan 07 New option lwsig=-1; updated iosigh
ckino:  Cu   24 Jul 06 MPI broadcast of sigma
ckino:  Cu   24 Feb 05 Rotation of sigma matrix now in hrs
ckino:  Cu             Switch lwsig to flag transformed sigm to be generated
ckino:  Cu   27 Jan 05 New 40000s digit mode (allow qp mismatch)
ckino:  Cu   20 Dec 04 New --wsig:rot and --wsig:phase switches
ckino:  Cu             New 10000s digit mode
ckino:  Cu   22 Nov 04 New --wsig:newkp switch
ckino:  Cu   25 Sep 04 rdsigm can read sigm(q) file for offset q-mesh
ckino:  Cu             rdsigm can write sigm(q) for mesh sbz->nkabc
ckino:  Cu    5 Apr 04 bug fix for --wsig switch.  New --mixsig switch
ckino:  Cu   15 Feb 04 Enable rdsigm to read a linear combination of
ckino:  Cu             data from file ifis and data from file 'sigm1'
ckino:  Cu    4 Jan 04 rdsigm checks bloch sum against ham->rsstol
ckino:  Cu   10 Oct 03 rdsigm can spin-split a non-spin-polarized sigma
ckino:  Cu   23 Sep 03 rdsigm can now read sigm fit to ASA hamiltonian
ckino:  Cu   11 Jul 03 enable write of sigma to full BZ
ckino:  Cu   24 May 03 Redesigned switches in accord with new interpolation
ckino:  Cu   09 May 03 Added symmetrization of sigma
ckino:  Cu   10 Jan 03 Some adaptations for Sergey's interpolation
ckino:  Cu   14 Aug 02 Added option to orthogonalize sigm.
ckino:  Cu   27 Jul 02 first created
ckino:  C ------------------------------------------------------------
ckino:        use m_psigd
ckino:  C     implicit none
ckino:  C ... Passed parameters
ckino:        integer mode,nbas,ifis,ndimh,lwsig
ckino:  C     integer osig(nbas),otau(nbas),oppi(nbas)
ckino:        type(s_lat)::slat
ckino:        type(s_ham)::sham
ckino:        type(s_bz)::sbz
ckino:
ckino:  C     double precision ssite(1),sspec(1)
ckino:        double precision rsrnge
ckino:  C ... Local parameters
ckino:        logical llshft(3),cmdopt,ltrans,lphase,lsplts,lnwmsh,
ckino:       .latvec,lfbzin,lfbzout
ckino:        integer parg,isw,lonesp
ckino:        character outs*80,out2*80,dc*1,rots*120
ckino:        integer i,j,ifis2,ifiz,isp,nsp,nglob,ival,ldham(16),hreal,lrsig,
ckino:       .moditp,nttabs,nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps,
ckino:       .ndhrs,j1,k1,k2,k3,iq1,nl,nspc,mxorb,nqsig,modsgp,nmin,nmax,
ckino:       .fopna,kcplx,lrot,iprint,stdo,lssym,ledit,niax,nscnts,
ckino:       .ndims,ndimz,iq,n123(4),lcore,lhigh
ckino:        parameter (niax=10)
ckino:        equivalence (nk1,nkxyz(1)),(nk2,nkxyz(2)),(nk3,nkxyz(3))
ckino:  ckino Dec.14.2011:                   integer,pointer :: iv_p_ontabs(:) =>NULL()
ckino:                   integer,allocatable :: iv_a_ontabs(:)
ckino:
ckino:  ckino Dec.14.2011:                  integer,pointer :: iv_p_oiaxs(:) =>NULL()
ckino:                  integer,allocatable :: iv_a_oiaxs(:)
ckino:
ckino:                 real(8),pointer :: rv_p_opos(:) =>NULL()
ckino:
ckino:                real(8),pointer :: rv_p_og(:) =>NULL()
ckino:
ckino:               real(8),pointer :: rv_p_oag(:) =>NULL()
ckino:
ckino:              real(8),pointer :: rv_p_ohrs(:) =>NULL()
ckino:
ckino:             integer,pointer :: iv_p_oistab(:) =>NULL()
ckino:
ckino:            integer,pointer :: iv_p_oiprmb(:) =>NULL()
ckino:
ckino:           real(8),pointer :: rv_p_oqsig(:) =>NULL()
ckino:
ckino:          real(8),pointer :: rv_p_oqp(:) =>NULL()
ckino:
ckino:         integer,pointer :: iv_p_ooffh(:) =>NULL()
ckino:
ckino:        integer ,allocatable :: gstar_iv(:)
ckino:        integer ,allocatable :: ipq_iv(:)
ckino:        real(8) ,allocatable :: wgt_rv(:)
ckino:        complex(8) ,allocatable :: sfz_zv(:)
ckino:        complex(8) ,allocatable :: wk_zv(:)
ckino:        complex(8) ,allocatable :: sigm_zv(:)
ckino:        complex(8) ,allocatable :: sigm2_zv(:)
ckino:        real(8) ,allocatable :: delt_rv(:)
ckino:        integer ,allocatable :: istb2_iv(:)
ckino:
ckino:        double precision tmp(3)
ckino:        double complex ztmp(3)
ckino:        real(8),allocatable:: evls(:),evlz(:),sigii(:)
ckino:        complex(8),allocatable:: siglda(:,:),z(:,:),sigo(:,:)
ckino:  C     For offset q mesh
ckino:        integer is(3),lshft(3),ifac(3),lqoffo
ckino:        double precision rb(3,3),qb(3,3),qoffi(3),qoffo(3)
ckino:  C MPI
ckino:        logical mlog
ckino:        integer procid,mpipid,master
ckino:  C ... for mixing sigma files: Use alf1*sigm + alf2*sigm1
ckino:        double precision alf(2)
ckino:  C     integer os
ckino:        double precision plat(3,3),qp(3),tolq,rsstol,sigp(10)
ckino:        double precision emin,emax,asig,bsig,efit,qoff(3),rotm(3,3),ddot
ckino:        equivalence (emin,sigp(3)),(emax,sigp(5))
ckino:        equivalence (asig,sigp(6)),(bsig,sigp(7)),(efit,sigp(8))
ckino:  C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
ckino:        double precision qk
ckino:        integer jj1,jj2,jj3,k
ckino:
ckino:        real(8):: rv_dummy(1)
ckino:        complex(8):: zv_dummy(1)
ckino:  C.....titus
ckino:  C      integer ib,ib1,ib2,nlmaa,lidim,ldh,osfz1,nsp1
ckino:  C      integer nkap0,n0H
ckino:  C      parameter  (nkap0=3,n0H=5)
ckino:  C      integer offH(n0H,nkap0,nbas)
ckino:  C.....titus
ckino:  C ... Heap
ckino:         integer ::iwdummy
ckino:
ckino:
ckino:        real(8):: eseavr
ckino:        integer:: debugmode,nhrss
ckino:        real(8),allocatable::hrss(:)
ckino:  ccccccccccccccccccccccccccccccc
ckino:  c      integer:: oss
ckino:  ccccccccccccccccccccccccccccccc
ckino:        logical:: isanrg, l_dummy_isanrg
ckino:
ckino:
ckino:
ckino:        qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
ckino:       .(jj2*ifac(2)-1)*qb(k,2) +
ckino:       .(jj3*ifac(3)-1)*qb(k,3)
ckino:
ckino:  C --- Setup and printout ---
ckino:        procid = mpipid(1)
ckino:        master = 0
ckino:        mlog = cmdopt('--mlog',6,0,outs)
ckino:
ckino:  C     tolq is allowed fuzz in qp
ckino:        tolq = 1d-6
ckino:  Changenglob      stdo = nglob('stdo')
ckino:        stdo = globalvariables%stdo
ckino:        ledit = 0
ckino:        lrsig = mod(mode,10)
ckino:        lwsig = 0
ckino:        moditp = mod(mode/10,10)
ckino:        alf(1) = 1
ckino:        alf(2) = 0
ckino:        nscnts = 0
ckino:        lqoffo = 0
ckino:        lfbzout = .false.
ckino:
ckino:  C     Switch flagging whether input file has symops or not
ckino:        hreal = 0
ckino:        if (lrsig .ge. 4) hreal = 1
ckino:        if (lrsig .eq. 0) return
ckino:
ckino:        ltrans = mod(mode/1000,10) .eq. 0
ckino:        kcplx = 1
ckino:        if (mod(mode/1000,10) .eq. 1) kcplx = 0
ckino:        lphase = mod(mode/1000,10) .eq. 0
ckino:        lssym = mod(mode/10000,10)
ckino:  Ckino isanrg is logical function,       call isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)
ckino:        l_dummy_isanrg=isanrg(lssym,0,7,'rdsigm:','10000s digit mode',.true.)
ckino:
ckino:  C ... Command-line switches for input sigma and defaults
ckino:  #ifdef USEOPTIONAL
ckino:         call sigflg ( dc=' ' , sio='rsig' ,strn= ' '
ckino:       . , lfbz=lfbzin  , lonesp=lonesp , qpoff=qoffi  )
ckino:  #else
ckino:         call sigflg ( ' ' , 'rsig' , ' ' , iwdummy , iwdummy , iwdummy
ckino:       . , lfbzin , iwdummy , lonesp , qoffi , iwdummy )
ckino:  #endif
ckino:
ckino:        if (cmdopt('--rsig',6,0,outs)) then
ckino:          out2 = outs(7:)
ckino:          dc = out2(1:1)
ckino:  #ifdef USEOPTIONAL
ckino:           call sigflg ( dc=dc , sio='rsig' ,strn= out2
ckino:       .   , lfbz=lfbzin  ,lonesp= lonesp , qpoff=qoffi  )
ckino:  #else
ckino:           call sigflg ( dc , 'rsig' , out2 , iwdummy , iwdummy , iwdummy
ckino:       .   , lfbzin , iwdummy , lonesp , qoffi , iwdummy )
ckino:  #endif
ckino:
ckino:        endif
ckino:  #if 1
ckino:        if ( lfbzin .or. mod(lssym,2).eq.1 ) then
ckino:          lfbzin = .true.
ckino:        else
ckino:          lfbzin = .false.
ckino:        endif
ckino:  #else
ckino:        lfbzin = lfbzin .or. mod(lssym,2)
ckino:  #endif
ckino:
ckino:        if (cmdopt('--mixsig=',9,0,outs)) then
ckino:          j = 0
ckino:          i = parg('--mixsig=',4,outs,j,len(outs),
ckino:       .  ', ',2,2,ifac,alf)
ckino:          if (i .ne. 1 .and. i .ne. 2) call rx
ckino:       .  ('rdsigm: failed to parse arguments to --mixsig=')
ckino:        endif
ckino:
ckino:        call info(10,1,0,' RDSIGM: read file sigm and create '//
ckino:       .'%?#n#REAL#COMPLEX# sigma(R) by FT ...',hreal,0)
ckino:
ckino:        i_copy_size=size(sham%ldham)
ckino:        call i8icopy(i_copy_size,sham%ldham,1,ldham,1)
ckino:         iv_p_ooffh => sham%iv_p_ooffh
ckino:
ckino:         iv_p_oiprmb => sham%iv_p_oindxo
ckino:
ckino:        i_copy_size=size(sham%sigp)
ckino:        call dcopy(i_copy_size,sham%sigp,1,sigp,1)
ckino:        rsstol=sham%rsstol
ckino:
ckino:        modsgp = nint(sigp(1))
ckino:        nmin   = nint(sigp(2))
ckino:        nmax   = nint(sigp(4))
ckino:        call info2(30,0,0,
ckino:       .'%9fSigm will be approximated by:  '//
ckino:       .'%?#(n==0)#Simple bloch sum##%-1j'//
ckino:       .'%?#(n==1)#diagonal Sigma for high and low states##%-1j'//
ckino:       .'%?#(n==2)#Perturbation for higher kappa blocks##%-1j'//
ckino:       .'%?#(n==3)#Interpolation from known points##%-1j'//
ckino:       .' ',moditp,0)
ckino:        if (moditp .eq. 1 .or. moditp .eq. 2) then
ckino:          call info5(30,0,0,
ckino:       .  '%9fApproximate sigma '//
ckino:       .  '%?#(n<0)#for energies E(lda)<%d; and %-2j#%-1j#'//
ckino:       .  '%?#(n>0)#for states up to %-1jn=%i; and ##%j'//
ckino:       .  '%?#(n<=0)#for energies E(lda)>%d%-2j#%-1j#'//
ckino:       .  '%?#(n>0)#for states above %-1jn=%i##%j',
ckino:       .  nmin,emin,nmax,emax,0)
ckino:          call info5(30,0,0,'%9fFor high states '//
ckino:       .  '%?#(n==0)#Sigii > %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
ckino:       .  '%?#(n==1)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
ckino:       .  '%?#(n==3)#Sigii = %,3;4d + %,3;4d * E(lda)%-2j##%-1j'//
ckino:       .  '%?#(n==2)#%,3;4d < Sigii < %,3;4d%-2j##%-1j'//
ckino:       .  ' ',modsgp,asig,bsig,0,0)
ckino:          if (nmax .gt. 0) call info5(30,0,0,
ckino:       .  '%9fLinear fit to diagonal sigma for E(lda)>%d'//
ckino:       .  '%?#(n==3)# and n<%i'//
ckino:       .  ' ',efit,modsgp,nmax,0,0)
ckino:          if (nmax .eq. 0) call info5(30,0,0,
ckino:       .  '%9fLinear fit to diagonal sigma for E(lda)>%d'//
ckino:       .  '%?#(n==3)# and E(lda)<%d'//
ckino:       .  ' ',efit,modsgp,emax,0,0)
ckino:  Ckino isanrg is logical function,         call isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
ckino:          l_dummy_isanrg=isanrg(modsgp,0,3,'rdsigm:','sig fit mode',.true.)
ckino:        endif
ckino:        if (lfbzin .or. ddot(3,qoffi,1,qoffi,1).ne.0) then
ckino:          call info5(30,0,0,
ckino:       .  '%?#(n==1)#%9fFile sigm in saved in FBZ: '//
ckino:       .  'symmetrization suppressed##'//
ckino:       .  '%?#(n==1)#%N%9fFile k-mesh is offset by:%3:2g#%j#'//
ckino:       .  ' ',isw(lfbzin),isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi,0,0)
ckino:        endif
ckino:
ckino:        if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
ckino:          call info5(30,0,0,
ckino:       .  '%9fUse for sigma: %d*(file sigm)%?#n# + '//
ckino:       .  '%d*(file sigm1)',alf(1),isw(alf(2).ne.0),alf(2),0,0)
ckino:        endif
ckino:
ckino:  C ... Setup for sigm transformation modes
ckino:        if (cmdopt('--wsig',6,0,outs) .or. cmdopt('-wsig',5,0,outs)) then
ckino:          if (procid .eq. master) then
ckino:            out2 = outs(7:)
ckino:            if (outs(1:5) .eq. '-wsig') out2 = outs(6:)
ckino:            dc = out2(1:1)
ckino:            call sigflg(dc,'wsig',out2,ledit,lwsig,lrot,lfbzout,lnwmsh,
ckino:       .    lonesp,qoffo,rots)
ckino:            lqoffo = isw(ddot(3,qoffo,1,qoffo,1) .ne. 0)
ckino:
ckino:  C       Setup, sigm transformation modes
ckino:            if (lwsig .ne. 0) then
ckino:
ckino:              call info5(10,0,0,
ckino:       .      '%9ftrans=%i:%19p'//
ckino:       .      '%?#(n==-1)#transform sigm, LDA -> orbital'//
ckino:       .      ' basis (file sigm3)##%-1j'//
ckino:       .      '%?#(n==1|n==2)#transform sigm, orbital ->'//
ckino:       .      ' LDA basis (file sigm2)##%-1j'//
ckino:       .      '%?#(n==5)#transform sigm, orbital ->'//
ckino:       .      ' orbital basis (file sigm2)##%-1j'//
ckino:       .      '%?#(n==2|n==5)#%N%19fhigh-energy block approximated by'//
ckino:       .      ' a diagonal matrix##%-1j'//
ckino:       .      '%?#(n==3)#make LDA eigenvalues, eigenvectors##%-1j'//
ckino:       .      '%j%?#(n==1)#%N%19fshift k-mesh by%3:2g#%j#'//
ckino:       .      '%?#(n==1)#%N%19fsigma in FBZ##'//
ckino:       .      ' ',lwsig,lwsig,lqoffo,qoffo,isw(lfbzout))
ckino:
ckino:  C     ... Transform sigm(lda)->sigm(orbital)
ckino:              if (lwsig .eq. -1) then
ckino:
ckino:  C         Read dimension of hamiltonian in evec file
ckino:  C         Check that evec file contains LDA eigenvectors
ckino:                call info0(30,0,0,'%9fchecking compatibility between '//
ckino:       .        'sigm2, evec files ...')
ckino:                ifiz = fopna('evec',-1,4)
ckino:                rewind ifiz
ckino:                call iosigh(2,3,i,ndimz,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
ckino:       .        lshft(3),ifiz)
ckino:  C         Verify that sigm2 file contains sigm(LDA); read dimension
ckino:                ifis = fopna('sigm2',-1,4)
ckino:                call iosigh(0,nscnts,i,ndims,nk1,nk2,nk3,nqp,lshft(1),
ckino:       .        lshft(2),lshft(3),ifis)
ckino:                if (nscnts .ne. 1 .and. nscnts .ne. 2)
ckino:       .        call rx1('File sigm2 (trans=%i) does not contain '//
ckino:       .        'sigm(LDA): create with trans=1 or 2 ',nscnts)
ckino:  C         Check that evec file matches (nsp,nqp)
ckino:                rewind ifiz
ckino:                if (lfbzin) nqp = nk1*nk2*nk3
ckino:                call iosigh(3,3,i,ndimz,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
ckino:       .        lshft(3),ifiz)
ckino:
ckino:                allocate(evls(ndims),evlz(ndimz),sigii(ndimz))
ckino:                allocate(siglda(ndims,ndims),z(ndimz,ndimz),sigo(ndimz,ndimz))
ckino:
ckino:  C         Write header for sigm3
ckino:                ifis2 = fopna('sigm3',-1,4)
ckino:                call iosigh(0,0,nsp,ndimz,nk1,nk2,nk3,nqp,
ckino:       .        lshft(1),lshft(2),lshft(3),-ifis2)
ckino:
ckino:  C         For each qp, rotate sigm(lda)->sigm(orb)
ckino:                do  iq = 1, nqp
ckino:  C         Check that qp match
ckino:                  read(ifis) qp, eseavr !add eseavr june2009
ckino:                  print *,' uuu1 eseavr=',iq,eseavr
ckino:                  read(ifiz) tmp
ckino:                  call daxpy(3,-1d0,qp,1,tmp,1)
ckino:                  if (ddot(3,tmp,1,tmp,1) .gt. 1d-9)
ckino:       .          call rx('rdsigm : qp mismatch, files sigm and evec')
ckino:                  call dpdump(evls,ndimh,ifis)
ckino:                  call dpdump(siglda,ndims**2*2,ifis)
ckino:                  call dpdump(evlz,ndimz,ifiz)
ckino:                  call dpdump(z,ndimz**2*2,ifiz)
ckino:                  call info5(30,0,0,' rdsigm:  kpt %i of %i, k=%3:2,5;5d',
ckino:       .          iq,nqp,qp,0,0)
ckino:                  write(stdo,542)
ckino:    542           format(8x,'Elda(evec)   Elda(sig)    diff       sig_ii',
ckino:       .          '      a+b*Elda')
ckino:                  do  i = 1, ndimz
ckino:                    write(stdo,543) i, evlz(i), evls(i),  evlz(i)-evls(i),
ckino:       .            dble(siglda(i,i)),asig+bsig*evlz(i)
ckino:    543             format(i5,5f12.6)
ckino:                  enddo
ckino:
ckino:  C          print *, '!!'; sigp(2) = 2
ckino:
ckino:  C         Determine lhigh = ndimz - ndims - lcore
ckino:  C         sigp(4) = nmax = ndimz - lhigh = ndims + lcore
ckino:                  sigp(4) = ndims + sigp(2)
ckino:
ckino:  C          print *, '!!'; ndimz = ndimz+5
ckino:  C          deallocate(evlz,sigii)
ckino:  C          allocate(evlz(ndimz),sigii(ndimz))
ckino:  C          evlz = evls(ndims)+1
ckino:  C          evlz(1:2) = -99
ckino:  C          evlz(1+2:ndims+2) = evls(1:ndims)
ckino:  C          evlz(ndimz) = 11
ckino:
ckino:                  call dpzero(sigii,ndimz)
ckino:                  call psigd(1,ndimz,ztmp,evlz,sigp,n123,sigii)
ckino:                  lcore = nint(sigp(2))
ckino:                  lhigh = ndimz - ndims - lcore
ckino:                   call rotevs ( 1 , ndimz , ndims , lcore , lhigh , iwdummy , siglda
ckino:       .           , sigii , z , sigo )
ckino:
ckino:                  ifis2 = fopna('sigm3',-1,4)
ckino:                  write(ifis2) qp
ckino:                  call dpdump(sigo,ndimh**2*2,-ifis2)
ckino:                enddo
ckino:
ckino:                deallocate(evls,evlz,sigii)
ckino:                deallocate(siglda,z,sigo)
ckino:
ckino:                call rx0('RDSIGM:  sigm(orbital basis) saved in file sigm3')
ckino:
ckino:              endif
ckino:
ckino:  C         No sigma for LDA eigenvectors; early exit
ckino:              if (lwsig .eq. 3 .or. lwsig .eq. 4) return
ckino:
ckino:            endif
ckino:          endif
ckino:        endif
ckino:
ckino:        call tcn('rdsigm')
ckino:
ckino:  C --- Read sigma(orbital basis) from file ---
ckino:
ckino:  C ... Read k-mesh parameters
ckino:        if (procid .eq. master) then
ckino:  C       Require file contains sigm(orbital basis)
ckino:          call iosigh(0,nscnts,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
ckino:       .  lshft(3),ifis)
ckino:          if (nscnts .ne. 0 .and. nscnts .ne. 5) then
ckino:            call iosigh(2,0,i,j,nk1,nk2,nk3,nqp,lshft(1),lshft(2),
ckino:       .    lshft(3),ifis)
ckino:          endif
ckino:        endif
ckino:
ckino:        call mpibc1(i,1,2,.false.,'rdsigm','nsp')
ckino:        call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
ckino:        call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
ckino:        call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
ckino:        if (i .eq. nsp) then
ckino:          lsplts = .false.
ckino:        elseif (i .gt. nsp) then
ckino:          call rx(
ckino:       .  'rdsigm: sigm file spin polarized but calculation is not')
ckino:        else
ckino:          lsplts = .true.
ckino:        endif
ckino:
ckino:  C     Check for consistency in remaining ndimh, nqp
ckino:  C     and set file sigma file pointer past header
ckino:        if (procid .eq. master) then
ckino:          rewind ifis
ckino:          if (lfbzin) nqp = nk1*nk2*nk3
ckino:          call iosigh(3,nscnts,i,ndimh,nk1,nk2,nk3,nqp,
ckino:       .  lshft(1),lshft(2),lshft(3),ifis)
ckino:        endif
ckino:        call mpibc1(nkxyz,3,2,mlog,'rdsigm','nkxyz')
ckino:        call mpibc1(nqp,1,2,mlog,'rdsigm','nqp')
ckino:        call mpibc1(lshft,3,2,mlog,'rdsigm','lshft')
ckino:        call info5(30,0,0,'%9fsigm file has %i irreducible QP: nk ='//
ckino:       .' ( %i %i %i )  shift=%3:1l',nqp,nk1,nk2,nk3,lshft)
ckino:  C     Sanity checks
ckino:        if (lsplts) call info0(30,0,0,
ckino:       .'%9f(warning) sigm file not spin pol .. splitting spins')
ckino:
ckino:        i_copy_size=size(slat%plat)
ckino:        call dcopy(i_copy_size,slat%plat,1,plat,1)
ckino:         rv_p_opos => slat%rv_p_opos
ckino:
ckino:
ckino:
ckino:        nsgrp=slat%nsgrp
ckino:         iv_p_oistab => slat%iv_p_oistab
ckino:
ckino:         rv_p_og => slat%rv_p_osymgr
ckino:
ckino:         rv_p_oag => slat%rv_p_oag
ckino:
ckino:
ckino:  C     Get leading dimension of hrs
ckino:         call offshp ( iv_p_oiprmb , 1 , nbas , 2 , 0 , ndimh , 0 , ndhrs
ckino:       . )
ckino:
ckino:  C     Setup for FFT
ckino:        call fftz30(nk1,nk2,nk3,k1,k2,k3)
ckino:        if (nk1.ne.k1 .or. nk2.ne.k2 .or. nk3.ne.k3) call rx(
ckino:       .'rdsigm: not ready for FFT w/ dimensions ne no div.')
ckino:  C ... Number of group operations for input file sigma
ckino:        if (mod(lssym,2) .eq. 0) nsgrps = nsgrp
ckino:        if (mod(lssym,2) .eq. 1) nsgrps = 1
ckino:        if (lfbzin) nsgrps = 1
ckino:  C ... Make is,ifac,qb,qlat,qoff
ckino:        do  8  i = 1, 3
ckino:          llshft(i) = lshft(i) .ne. 0
ckino:      8 continue
ckino:        call pshpr(0)
ckino:        call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
ckino:        call poppr
ckino:        do  i = 1, 3
ckino:          qoff(i) = qk(i,1,1,1) + qoffi(i)
ckino:        enddo
ckino:
ckino:  C ... Setup for r.s. hamiltonian: allocate and create iax,ntab; allocate hrs
ckino:        i = 1000 + 10*hreal + 0
ckino:  C      print *, '!!'; i = 1000*0 + 10*hreal + 0
ckino:        if (mod(lssym,4) .ge. 2) i = 10*hreal + 0
ckino:        if (ltrans) i = i + 10000
ckino:        if(debugmode()>0) print *,' go to hft2rs i=',i
ckino:            call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , 1 ,
ckino:       .    nsp , nbas , rv_p_og , rv_p_oag , nsgrp , rsrnge , iv_p_ooffh
ckino:  ckino Dec.14.2011:       .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , rv_p_opos , iv_p_ontabs ,
ckino:       .    , 1 , nbas , ndimh , ndimh , zv_dummy , plat , rv_p_opos , iv_a_ontabs ,
ckino:  ckino Dec.14.2011:       .    iv_p_oiaxs , ndhrs , rv_dummy )
ckino:       .    iv_a_oiaxs , ndhrs , rv_dummy )
ckino:
ckino:
ckino:
ckino:
ckino:  ckino Dec.14.2011:         nttabs = ival ( iv_p_ontabs , nbas + 1 )
ckino:         nttabs = ival ( iv_a_ontabs , nbas + 1 )
ckino:
ckino:  C     Allocate memory for hrs
ckino:        if (hreal .eq. 1) then
ckino:           allocate(rv_p_ohrs(ndhrs**2*nttabs*nsp))
ckino:           rv_p_ohrs(:)=0d0
ckino:        else
ckino:           allocate(rv_p_ohrs(2*ndhrs**2*nttabs*nsp))
ckino:           rv_p_ohrs(:)=0d0
ckino:        endif
ckino:
ckino:  c      stop 'xxxxxxxxxxxxxxxx takao test1 in rdsigm.F xxxxxxxxxx'
ckino:
ckino:
ckino:
ckino:  C --- Get list of irreducible k-points, and ipq and gstar arrays ---
ckino:        mxkp = nk1*nk2*nk3
ckino:         allocate(rv_p_oqsig(abs(3*mxkp)))
ckino:         if (3*mxkp<0) rv_p_oqsig(:)=0.0d0
ckino:
ckino:         allocate(rv_p_oqp(abs(3*mxkp)))
ckino:         if (3*mxkp<0) rv_p_oqp(:)=0.0d0
ckino:
ckino:        allocate(gstar_iv(abs(-mxkp-1)))
ckino:        if (-mxkp-1<0) gstar_iv(:)=0
ckino:
ckino:        gstar_iv = - 2
ckino:
ckino:        allocate(ipq_iv(mxkp))
ckino:        if (mxkp<0) ipq_iv(:)=0
ckino:
ckino:        allocate(wgt_rv(abs(-mxkp)))
ckino:        if (-mxkp<0) wgt_rv(:)=0.0d0
ckino:
ckino:        call info(20,1,0,
ckino:       .' q-points in full BZ where sigma calculable ...',0,0)
ckino:          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy
ckino:       .  , 0 , ipq_iv , rv_p_oqsig , wgt_rv , nqsig , mxkp , 0 , 0 )
ckino:
ckino:
ckino:
ckino:
ckino:        sham%nqsig=nqsig
ckino:         sham%rv_p_oqsig => rv_p_oqsig
ckino:
ckino:
ckino:        call dpzero ( wgt_rv , mxkp )
ckino:
ckino:        call info2(20,0,0,
ckino:       .' Irr. qp for which sigma is calculated ...',0,0)
ckino:          call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og
ckino:       .  , nsgrps , ipq_iv , rv_p_oqp , wgt_rv , nqps , mxkp , gstar_iv
ckino:       .  , 0 )
ckino:
ckino:
ckino:
ckino:  Ckino isanrg is logical function,       call isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
ckino:        l_dummy_isanrg=isanrg(nqp,nqps,nqps,'rdsigm','file sigm nqp',.true.)
ckino:  C     call prmx('q for which sigma is read',w(oqp),3,3,nqp)
ckino:        call info2(30,0,0,
ckino:       .'%?#(n==1)#%1fFile sigm has k-offset=%3:2g '//
ckino:       .'relative to the above#%j#'//
ckino:       .' ',isw(ddot(3,qoffi,1,qoffi,1).ne.0),qoffi)
ckino:
ckino:  C ... titus
ckino:  CC     List of sites
ckino:  C      call getoffH(w(ooffH),offH,n0H,nkap0,nbas)
ckino:  C      ib1 = 1
ckino:  C      ib2 = nbas
ckino:  C      lidim = offH(4,1,nbas+1)
ckino:  C
ckino:  C      print*,'lidim',lidim
ckino:  C
ckino:  C      do 100 ib = ib1, ib2
ckino:  C      rewind ifis
ckino:  C      call iosigh(0,nscnts,nsp1,ndimh,nk1,nk2,nk3,nqp,
ckino:  C     .  lshft(1),lshft(2),lshft(3),ifis)
ckino:  C      nlmaa = offH(4,1,ib+1) - offH(4,1,ib)
ckino:  C      print*,'nlmaa',nlmaa
ckino:  C      print*,'offH(4,1,ib)',offH(4,1,ib)
ckino:  C      print*,'offH(4,1,ib+1)',offH(4,1,ib+1)
ckino:  C ... titus
ckino:
ckino:  C --- Generate hrs = sigma(T) from file sigma(k) ---
ckino:  C ... Loop over spins, keeping sig in full BZ for only one spin
ckino:        if (procid .eq. master) then
ckino:          do  isp = 1, nsp
ckino:
ckino:  C     If sigma file not spin polarized, use sigma from spin 1
ckino:            if (isp .eq. 2 .and. lsplts) then
ckino:              call iosigh(3,0,1,ndimh,nk1,nk2,nk3,nqp,
ckino:       .      lshft(1),lshft(2),lshft(3),ifis)
ckino:            endif
ckino:
ckino:  C     If to mix in some other sigma file, set up ifis2
ckino:            if (alf(2) .ne. 0) then
ckino:              if (isp .eq. 1 .or. isp .eq. 2 .and. lsplts) then
ckino:                ifis2 = fopna('sigm1',-1,4)
ckino:                call iosigh(3,0,nsp,ndimh,nk1,nk2,nk3,nqp,
ckino:       .        lshft(1),lshft(2),lshft(3),ifis2)
ckino:              endif
ckino:            endif
ckino:
ckino:  C ... File sigma in irr BZ -> sigma in full BZ by rotations
ckino:            allocate(sfz_zv(k1*k2*k3*ndimh**2))
ckino:            if (k1*k2*k3*ndimh**2<0) sfz_zv(:)=0.0d0
ckino:
ckino:  C ... titus
ckino:  C      call defcc(osfz,k1*k2*k3*ndimh*nlmaa)
ckino:  C      call defcc(osfz1,k1*k2*k3*ndimh*nlmaa)
ckino:  C ... titus
ckino:            allocate(wk_zv(ndimh**2))
ckino:            if (ndimh**2<0) wk_zv(:)=0.0d0
ckino:
ckino:  Changenglob          nl = nglob('nl')
ckino:            nl = globalvariables%nl
ckino:            nspc = 1
ckino:            allocate(sigm_zv(ndimh**2))
ckino:            if (ndimh**2<0) sigm_zv(:)=0.0d0
ckino:
ckino:  ccccccccccccccccccccccccc
ckino:  ctakaox
ckino:  c      call defcc(oss,ndimh**2)
ckino:  ccccccccccccccccccccccccc
ckino:            do  iq1 = 1, nqps
ckino:              read(ifis) qp,  eseavr !june2009  add eseavr
ckino:              print *,' uuu2 iq1 eseavr=',iq1,eseavr
ckino:              call dpscop(qp,tmp,3,1,1,1d0)
ckino:               call dpsadd ( tmp , rv_p_oqp , 3 , 1 , 3 * iq1 - 2 , - 1d0 )
ckino:
ckino:              call dpsadd(tmp,qoffi,3,1,1,-1d0)
ckino:              if (.not. latvec(1,tolq,plat,tmp)) then
ckino:                 call dpscop ( rv_p_oqp , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 )
ckino:
ckino:                if (lssym .ge. 4) then
ckino:                  call info5(20,0,0,' rdsigm (warning) qp %i : '//
ckino:       .          'expected%3;8,4D read%3;8,4D',iq1,tmp,qp,0,0)
ckino:                  call dpscop(tmp,qp,3,1,1,1d0)
ckino:                else
ckino:                  print 456, iq1,tmp,qp
ckino:    456           format(' rdsigm: error on reading qp no',i4/
ckino:       .          ' Expected qp=',3f12.6/
ckino:       .          '     File qp=',3f12.6)
ckino:                  call rx(' incompatible q-mesh')
ckino:                endif
ckino:              endif
ckino:  C       File read sigma
ckino:              if (nscnts .eq. 5) call dpdump(tmp,1,ifis)
ckino:              call dpdump ( sigm_zv , ndimh * * 2 * 2 , ifis )
ckino:
ckino:  C        if (isp .eq. 2) then
ckino:  C          print *, iq1,isp
ckino:  C         call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
ckino:  C        endif
ckino:
ckino:  cccccccccccccccccccccccccccccccccccccc
ckino:  ctakaox !see also bndfp.F and hambls.F
ckino:  c        if(iq1==1) then
ckino:  c          print *,'reading ovlmat'
ckino:  c          ifio = fopna('ovlmat',-1,4)
ckino:  c        endif
ckino:  c        call dpdump(w(oss),ndimh**2*2,ifio)
ckino:  c        print *,'read ovlmat for iq1=',iq1
ckino:  c        call divctakao(w(osigm),w(oss),ndimh)
ckino:  ccccccccccccccccccccccccccccccccccccc
ckino:
ckino:              if (alf(1) .ne. 1d0) then
ckino:                call dscal ( ndimh * * 2 * 2 , alf ( 1 ) , sigm_zv , 1 )
ckino:
ckino:              endif
ckino:              if (alf(2) .ne. 0d0) then
ckino:                read(ifis2) qp
ckino:                call dpscop(qp,tmp,3,1,1,1d0)
ckino:                 call dpsadd ( tmp , rv_p_oqp , 3 , 1 , 3 * iq1 - 2 , - 1d0 )
ckino:
ckino:                call dpsadd(tmp,qoffi,3,1,1,-1d0)
ckino:                if (abs(tmp(1))+abs(tmp(2))+abs(tmp(3)) .gt. tolq) then
ckino:                   call dpscop ( rv_p_oqp , tmp , 3 , 3 * iq1 - 2 , 1 , 1d0 )
ckino:
ckino:                  if (lssym .ge. 4) then
ckino:                    call info5(20,0,0,' rdsigm (warning) qp %i : '//
ckino:       .            'expected%3;8,4D read%3;8,4D (file sigm1)',iq1,tmp,qp,
ckino:       .            0,0)
ckino:                  else
ckino:                    call rx(' file sigm1 has incompatible irr mesh')
ckino:                  endif
ckino:                endif
ckino:  C         call dpsadd(qp,w(oqp),3,1,3*iq1-2,1d0)
ckino:                allocate(sigm2_zv(ndimh**2))
ckino:                if (ndimh**2<0) sigm2_zv(:)=0.0d0
ckino:
ckino:                call dpdump ( sigm2_zv , ndimh * * 2 * 2 , ifis2 )
ckino:
ckino:                call daxpy ( ndimh * * 2 * 2 , alf ( 2 ) , sigm2_zv , 1 ,
ckino:       .        sigm_zv , 1 )
ckino:
ckino:                if (allocated(sigm2_zv)) deallocate(sigm2_zv)
ckino:
ckino:              endif
ckino:  C       call zprm('sigm',2,w(osigm),ndimh,ndimh,ndimh)
ckino:
ckino:              i = 100*kcplx + 00
ckino:              if (lphase) i = i+40
ckino:  C ... titus
ckino:  C        call hamfb2(nbas,nl,w(ooffH),w(oiprmb),ib,i,w(opos),iq1,nk1,
ckino:  C     .    nk2,nk3,k1,k2,k3,w(oipq),w(oistab),w(og),w(oag),
ckino:  C     .    w(ogstar),ifac,lidim,ndimh,nlmaa,ndimh,nlmaa,nspc,qb,
ckino:  C     .    w(osigm),w(owk),
ckino:  C     .    w(owk),w(osfz),w(osfz1))
ckino:  C ... titus
ckino:                    call hamfb3 ( nbas , nl , iv_p_ooffh , iv_p_oiprmb , i , rv_p_opos
ckino:       .            , iq1 , nk1 , nk2 , nk3 , k1 , k2 , k3 , ipq_iv , iv_p_oistab
ckino:       .            , rv_p_og , rv_p_oag , gstar_iv , ifac , ndimh , ndimh , ndimh
ckino:       .            , nspc , qb , sigm_zv , wk_zv , wk_zv , sfz_zv )
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:            enddo
ckino:
ckino:
ckino:  cccccccccccccccccccccccccccccccccccccccccccc
ckino:  ctakaox !see also bndfp.F and hambls.F
ckino:  c      print *,' end of iq1 loop'
ckino:  c      call fclose(ifio)
ckino:  c      call rlse(oss)
ckino:  c      print *,' end of iq1 loop xxxx'
ckino:  cccccccccccccccccccccccccccccccccccccccccccc
ckino:
ckino:
ckino:
ckino:  C      print *, 'one element of sigma(k)'
ckino:  C      call pvtrod(w(osfz),k1,k2,k3,ndimh)
ckino:
ckino:  C ... FT sfz and copy sfz(T) to hrs
ckino:            i = 10*hreal + 1
ckino:            if (ltrans) i = i + 10000
ckino:  C ... titus
ckino:  C      call hft2rs(i,nk1,nk2,nk3,k1,k2,k3,qoff,isp,nsp,nbas,w(og),w(oag)
ckino:  C     .  ,nsgrps,rsrnge,w(ooffH),ib,ib,nlmaa,ndimh,w(osfz),plat,
ckino:  C     .  w(opos),ontabs,
ckino:  C     .  oiaxs,ndhrs,w(ohrs))
ckino:  C ... titus
ckino:                 call hft2rs ( i , nk1 , nk2 , nk3 , k1 , k2 , k3 , qoff , isp
ckino:       .         , nsp , nbas , rv_p_og , rv_p_oag , nsgrps , rsrnge , iv_p_ooffh
ckino:  ckino Dec.14.2011:       .         , 1 , nbas , ndimh , ndimh , sfz_zv , plat , rv_p_opos , iv_p_ontabs
ckino:       .         , 1 , nbas , ndimh , ndimh , sfz_zv , plat , rv_p_opos , iv_a_ontabs
ckino:  ckino Dec.14.2011:       .         , iv_p_oiaxs , ndhrs , rv_p_ohrs )
ckino:       .         , iv_a_oiaxs , ndhrs , rv_p_ohrs )
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:  C      print *, 'print out sigma(T)'
ckino:  C      call pvtrof(w(ohrs),ndhrs,nsp,nk1,nk2,nk3,w(oiaxs),1,
ckino:  C     .  nttabs,w(osfz),ndimh)
ckino:  C      stop
ckino:
ckino:  C ... End loop over spins
ckino:            if (allocated(sigm_zv)) deallocate(sigm_zv)
ckino:            if (allocated(wk_zv)) deallocate(wk_zv)
ckino:            if (allocated(sfz_zv)) deallocate(sfz_zv)
ckino:
ckino:          enddo
ckino:  C ... titus
ckino:  C 100  continue
ckino:  C     exit of MPI master loop
ckino:        endif
ckino:        if (allocated(wgt_rv)) deallocate(wgt_rv)
ckino:        if (allocated(ipq_iv)) deallocate(ipq_iv)
ckino:        if (allocated(gstar_iv)) deallocate(gstar_iv)
ckino:
ckino:
ckino:
ckino:        sham%ndhrs=ndhrs
ckino:  ckino Dec.14.2011:         sham%iv_p_ontabs => iv_p_ontabs
ckino:  ckino Dec.14.2011:         sham%iv_a_ontabs => iv_a_ontabs
ckino:        call move_alloc(from= iv_a_ontabs ,to= sham%iv_a_ontabs )
ckino:
ckino:  ckino Dec.14.2011:         sham%iv_p_oiaxs => iv_p_oiaxs
ckino:  ckino Dec.14.2011:         sham%iv_a_oiaxs => iv_a_oiaxs
ckino:        call move_alloc(from= iv_a_oiaxs ,to= sham%iv_a_oiaxs )
ckino:
ckino:         sham%rv_p_ohrs => rv_p_ohrs
ckino:
ckino:
ckino:
ckino:  C ... Check that FT replicates sigma(k) at file qp to tolerance rsstol
ckino:  C     rsstol = 1d-6
ckino:  C     print *, '!!'
ckino:  C     if (alf(1) .ne. 1 .or. alf(2) .ne. 0 .or. .true.) then
ckino:        if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
ckino:          call info0(20,0,0,' skipping check of Bloch summed sigma'//
ckino:       .  ' (sigma scaled) ')
ckino:        else
ckino:          j1 = 0
ckino:          if (lssym .ge. 4) j1 = 2
ckino:          if (procid .eq. master) then
ckino:            call dpzero(tmp,3)
ckino:                 call chksgr ( j1 , ltrans , kcplx , plat , nsp , ndimh , ifis
ckino:  ckino Dec.14.2011:       .         , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs
ckino:       .         , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino:       .         , rsstol , i , rv_p_oqp , nbas , 0 , rotm , iwdummy )
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:          endif
ckino:        endif
ckino:
ckino:        if ( .not. ( cmdopt ( '--wsig' , 6 , 0 , outs ) .or.cmdopt (
ckino:       .'-wsig' , 5 , 0 , outs ) ) ) then
ckino:           if (associated(rv_p_oqp)) deallocate(rv_p_oqp)
ckino:
ckino:        endif
ckino:
ckino:
ckino:  C ... Broadcast R.S. sigma before symmetrization
ckino:        if (hreal .eq. 1) then
ckino:           call mpibc1 ( rv_p_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog
ckino:       .   , 'rdsigm' , 'hrs' )
ckino:
ckino:        else
ckino:           call mpibc1 ( rv_p_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4
ckino:       .   , mlog , 'rdsigm' , 'hrs' )
ckino:
ckino:        endif
ckino:
ckino:  C --- Symmetrize hrs ---
ckino:  C     Best to use nsgrp here instead of nsgrps.
ckino:  C     However, neighbor table was generated for nsgrps, and may be
ckino:  C     increase if nsgrp>nsgrps.  This branch will fail in that case.
ckino:        if (mod(lrsig,4) .ge. 2 .and. nsgrp .gt. 1) then
ckino:
ckino:  cC       Allocate memory for hrs
ckino:  c        if (hreal .eq. 1) then
ckino:  c          call defdr(ohrss,ndhrs**2*nttabs*nsp)
ckino:  cCchp1           call dcopy(ndhrs**2*nttabs*nsp*1,w(ohrs),1,w(ohrss),1)
ckino:  c           call dcopy ( ndhrs * * 2 * nttabs * nsp * 1 , rv_p_ohrs , 1 ,
ckino:  c     .     w ( ohrss ) , 1 )
ckino:  c        else
ckino:  c          call defdc(ohrss,ndhrs**2*nttabs*nsp)
ckino:  cCchp1           call dcopy(ndhrs**2*nttabs*nsp*2,w(ohrs),1,w(ohrss),1)
ckino:  c           call dcopy ( ndhrs * * 2 * nttabs * nsp * 2 , rv_p_ohrs , 1 ,
ckino:  c     .     w ( ohrss ) , 1 )
ckino:  c        endif
ckino:          nhrss= size(rv_p_ohrs)
ckino:          allocate(hrss(nhrss))
ckino:          hrss = rv_p_ohrs
ckino:
ckino:  Changenglob        mxorb = nglob('mxorb')
ckino:          mxorb = globalvariables%mxorb
ckino:  C       10s digit distributes ri-rj pairs to avg (ri-rj),(rj-ri) pairs
ckino:  C       Use 10*2 for hermitian matrices, 10*1 for symmetric ones
ckino:          i = 1-hreal + 10*2 + 100*(1-hreal) + 100000
ckino:          allocate(istb2_iv(nsgrp*nbas))
ckino:          if (nsgrp*nbas<0) istb2_iv(:)=0
ckino:
ckino:           call istbpm ( iv_p_oistab , nbas , nsgrp , istb2_iv )
ckino:
ckino:
ckino:  Changenglob        nl = nglob('nl')
ckino:          nl = globalvariables%nl
ckino:                call rsmsym ( i , plat , mxorb , iv_p_oiprmb , ndimh , nbas ,
ckino:  ckino Dec.14.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_p_ontabs , iv_p_oiaxs
ckino:  ckino Dec.14.2011:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_p_ontabs , iv_a_oiaxs
ckino:       .        rv_p_opos , nl , nsp , 1 , nttabs , iv_a_ontabs , iv_a_oiaxs
ckino:       .        , rv_p_og , istb2_iv , nsgrp , ndhrs , hrss , rv_p_ohrs )
ckino:
ckino:
ckino:
ckino:
ckino:          if (allocated(istb2_iv)) deallocate(istb2_iv)
ckino:  c        call rlse(ohrss)
ckino:          deallocate(hrss)
ckino:
ckino:
ckino:  C        Debugging
ckino:  C        ifis2 = fopna('out',-1,0)
ckino:  C        print *, 'hello',ndhrs,nttabs
ckino:  C        call ywrm(0,'sig',2,ifis2,'(9f20.10)',w(ohrs),1,ndhrs**2,
ckino:  C     .    ndhrs**2,nttabs*nsp)
ckino:  C        call rx0('done')
ckino:
ckino:  C   ... Check how well symmetrized FT replicates sigma(k) at file qp
ckino:          if (procid .eq. master) then
ckino:            rsstol = 0
ckino:            if (alf(1) .ne. 1 .or. alf(2) .ne. 0) then
ckino:              call info0(20,0,0,' skipping check of symmetrized sigma'//
ckino:       .      ' (sigma scaled) ')
ckino:            else
ckino:                  call chksgr ( 0 , ltrans , kcplx , plat , nsp , ndimh , ifis
ckino:  ckino Dec.14.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs
ckino:       .          , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino:       .          , rsstol , i , 0 , nbas , 0 , rotm , iwdummy )
ckino:
ckino:
ckino:
ckino:
ckino:            endif
ckino:          endif
ckino:
ckino:  C   ... Broadcast R.S. sigma after symmetrization
ckino:          if (hreal .eq. 1) then
ckino:             call mpibc1 ( rv_p_ohrs , ndhrs * * 2 * nttabs * nsp , 4 , mlog
ckino:       .     , 'rdsigm' , 'hrs' )
ckino:
ckino:          else
ckino:             call mpibc1 ( rv_p_ohrs , 2 * ndhrs * * 2 * nttabs * nsp , 4
ckino:       .     , mlog , 'rdsigm' , 'hrs' )
ckino:
ckino:          endif
ckino:
ckino:        endif
ckino:
ckino:  C ... lwsig=0 modes that internally transform sigma, save and exit
ckino:        if (lwsig .eq. 0) then
ckino:          if (cmdopt('--wsig',6,0,outs) .or. cmdopt('-wsig',5,0,outs)) then
ckino:            if (procid .eq. master) then
ckino:              ifis2 = fopna('sigm2',-1,4)
ckino:              rsstol = 0
ckino:              out2 = outs(7:)
ckino:              if (outs(1:5) .eq. '-wsig') out2 = outs(6:)
ckino:              dc = out2(1:1)
ckino:              call sigflg(dc,'wsig',out2,ledit,lwsig,lrot,lfbzout,lnwmsh,
ckino:       .      lonesp,qoffo,rots)
ckino:
ckino:              if (lonesp .eq. 1) then
ckino:                if (nsp .eq. 1) call rx(
ckino:       .        'onesp option nonsensical unless sigma file is spin pol')
ckino:                     call siged ( 1 , nbas , nsp , ndhrs , plat , rv_p_opos , ndimh
ckino:  ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_p_oiaxs , rv_p_ohrs
ckino:  ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino:       .             , iv_p_oiprmb , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino:       .             )
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:              endif
ckino:
ckino:              if (ledit .eq. 1) then
ckino:                     call siged ( 0 , nbas , nsp , ndhrs , plat , rv_p_opos , ndimh
ckino:  ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_p_oiaxs , rv_p_ohrs
ckino:  ckino Dec.14.2011:       .             , iv_p_oiprmb , hreal , iv_p_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino:       .             , iv_p_oiprmb , hreal , iv_a_ontabs , iv_a_oiaxs , rv_p_ohrs
ckino:       .             )
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:              endif
ckino:
ckino:  C       Handle case q-mesh for file output differs from file input
ckino:  C       ltmp = true in cases where file symops artificially reduced
ckino:  C        ltmp = mod(lssym,2) .eq. 1 .and. lwsig .ne. 3
ckino:  C        if (lnwmsh .or. ltmp) then
ckino:              if (lnwmsh) then
ckino:
ckino:                i_copy_size=size(sbz%nkabc)
ckino:                call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
ckino:                i_copy_size=size(sbz%lshft)
ckino:                call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)
ckino:                nqp=sbz%nkp
ckino:                 rv_p_oqp => sbz%rv_p_oqp
ckino:
ckino:
ckino:
ckino:                mxkp = nk1*nk2*nk3
ckino:                 allocate(rv_p_oqp(abs(3*mxkp)))
ckino:                 if (3*mxkp<0) rv_p_oqp(:)=0.0d0
ckino:
ckino:                allocate(gstar_iv(abs(-mxkp-1)))
ckino:                if (-mxkp-1<0) gstar_iv(:)=0
ckino:
ckino:                gstar_iv = - 2
ckino:
ckino:                allocate(ipq_iv(mxkp))
ckino:                if (mxkp<0) ipq_iv(:)=0
ckino:
ckino:                allocate(wgt_rv(abs(-mxkp)))
ckino:                if (-mxkp<0) wgt_rv(:)=0.0d0
ckino:
ckino:                do  i = 1, 3
ckino:                  llshft(i) = lshft(i) .ne. 0
ckino:                enddo
ckino:                call dpzero ( wgt_rv , mxkp )
ckino:
ckino:                if (lfbzout) then
ckino:                  call info0(20,1,0,
ckino:       .          ' Creating sigma for new mesh, full BZ (file sigm2) ...')
ckino:                    call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy
ckino:       .            , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 )
ckino:
ckino:
ckino:
ckino:                else
ckino:                  call info0(20,1,0,
ckino:       .          ' Creating sigma for new mesh (file sigm2) ...')
ckino:                    call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , rv_p_og
ckino:       .            , nsgrp , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , gstar_iv
ckino:       .            , 0 )
ckino:
ckino:
ckino:
ckino:
ckino:                endif
ckino:
ckino:              elseif (lfbzout) then
ckino:  cki              oqp = oqsig
ckino:                rv_p_oqp => rv_p_oqsig
ckino:                nqp = nqsig
ckino:                call info2(20,1,0,
ckino:       .        ' Creating sigma for full BZ, %i qp (file sigm2) ...',nqp,0)
ckino:              else
ckino:                call info2(20,0,0,
ckino:       .        ' Writing sigma for irr BZ, %i qp (file sigm2) ...',nqp,0)
ckino:              endif
ckino:
ckino:  C   ... Setup rotation matrix
ckino:              if (rots .ne. ' ') then
ckino:                call a2rotm(rots,.false.,0,rotm)
ckino:                if (iprint() .ge. 20) then
ckino:                  call info0(20,0,0,' Rotate sigma by rotation matrix:')
ckino:                  write (stdo,350) ((rotm(i,j),j=1,3),i=1,3)
ckino:    350           format(3f11.6)
ckino:                endif
ckino:                lrot = lrot + 1
ckino:              endif
ckino:
ckino:  C   ... Setup phase matrix
ckino:              if (lrot .ge. 2) then
ckino:                allocate(delt_rv(abs(-3*nbas)))
ckino:                if (-3*nbas<0) delt_rv(:)=0.0d0
ckino:
ckino:                call info0(0,0,0,' Phase shift of sigma.  '//
ckino:       .        'Translation vectors from shorps file:')
ckino:                call iopos ( .false. , - 1 , 'shorps' , nbas , delt_rv )
ckino:
ckino:              endif
ckino:  C       call prmx('q for which sigm is written',w(oqp),3,3,nqp)
ckino:
ckino:              call iosigh(0,nscnts,nsp,ndimh,nk1,nk2,nk3,nqp,
ckino:       .      lshft(1),lshft(2),lshft(3),-ifis2)
ckino:              if (lqoffo .ne. 0) then
ckino:  C         call prmx('qp',w(oqp),3,3,nqp)
ckino:                do  i = 1, nqp
ckino:                   call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 *
ckino:       .           i - 2 , 1 , 1d0 )
ckino:
ckino:                enddo
ckino:  C         call prmx('qp',w(oqp),3,3,nqp)
ckino:              endif
ckino:                  call chksgr ( 1 , ltrans , kcplx , plat , nsp , ndimh , ifis2
ckino:  ckino Dec.14.2011:       .          , hreal , iv_p_oiprmb , nttabs , iv_p_oiaxs , rv_p_ohrs , ndhrs
ckino:       .          , hreal , iv_p_oiprmb , nttabs , iv_a_oiaxs , rv_p_ohrs , ndhrs
ckino:       .          , rsstol , nqp , rv_p_oqp , nbas , lrot , rotm , delt_rv )
ckino:
ckino:
ckino:
ckino:
ckino:
ckino:              call fclose(ifis2)
ckino:              call rx0('done writing sigma, file sigm2')
ckino:  C      elseif (lwsig .eq. 1 .or. lwsig .eq. 2) then
ckino:  C        if (lfbzout) call info0(20,0,0,' ')
ckino:  C        if (lfbzout) goto 99
ckino:            endif
ckino:          endif
ckino:        endif
ckino:
ckino:  C     Broadcast entire sham structure (shouldn't be needed)
ckino:  C     call mpibc1(sham,i,4,mlog,'rdsigm','sham')
ckino:  C     Broadcast neighbor table (shouldn't be needed)
ckino:  C     call mpibc1(w(ontabs),nbas+1,2,mlog,'rdsigm','ntabs')
ckino:  C     call mpibc1(w(oiaxs),niax*nttabs,2,mlog,'rdsigm','iaxs')
ckino:
ckino:  C ... Transformations done by caller: write qp list to disk
ckino:     99 continue
ckino:        if (lwsig .ne. 0) then
ckino:
ckino:          call info0(30,0,0,' ')
ckino:
ckino:          if (lnwmsh) then
ckino:
ckino:            i_copy_size=size(sbz%nkabc)
ckino:            call i8icopy(i_copy_size,sbz%nkabc,1,nkxyz,1)
ckino:            i_copy_size=size(sbz%lshft)
ckino:            call i8icopy(i_copy_size,sbz%lshft,1,lshft,1)
ckino:            nqp=sbz%nkp
ckino:             rv_p_oqp => sbz%rv_p_oqp
ckino:
ckino:
ckino:            call info2(20,0,0,
ckino:       .    ' Transform sigma on new k-mesh:  '//
ckino:       .    'nk=%3:1i  shft=%3:1l',nkxyz,lshft)
ckino:          endif
ckino:
ckino:          if (lfbzout) then
ckino:            mxkp = nk1*nk2*nk3
ckino:            call info(20,0,0,' rdsigm: use full BZ (%i q-points)...',mxkp,0)
ckino:             allocate(rv_p_oqp(abs(3*mxkp)))
ckino:             if (3*mxkp<0) rv_p_oqp(:)=0.0d0
ckino:
ckino:            allocate(ipq_iv(mxkp))
ckino:            if (mxkp<0) ipq_iv(:)=0
ckino:
ckino:            allocate(wgt_rv(abs(-mxkp)))
ckino:            if (-mxkp<0) wgt_rv(:)=0.0d0
ckino:
ckino:
ckino:            i_copy_size=size(slat%plat)
ckino:            call dcopy(i_copy_size,slat%plat,1,plat,1)
ckino:
ckino:  C   ... Make is,ifac,qb,qlat
ckino:            do   i = 1, 3
ckino:              llshft(i) = lshft(i) .ne. 0
ckino:            enddo
ckino:            call pshpr(0)
ckino:            call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
ckino:            call poppr
ckino:              call bzmesh ( plat , qb , nk1 , nk2 , nk3 , llshft , iwdummy
ckino:       .      , 0 , ipq_iv , rv_p_oqp , wgt_rv , nqp , mxkp , 0 , 0 )
ckino:
ckino:
ckino:
ckino:  C       call pack2('bz nkp oqp',sbz,nqp,oqp)
ckino:            if (allocated(delt_rv)) deallocate(delt_rv)
ckino:            if (allocated(wgt_rv)) deallocate(wgt_rv)
ckino:
ckino:          endif
ckino:
ckino:  C ... Write points to disk, shifting by qoffo
ckino:          if (lqoffo .ne. 0) then
ckino:  C       call prmx('qp',w(oqp),3,3,nqp)
ckino:            do  i = 1, nqp
ckino:               call dmsadd ( rv_p_oqp , 1 , qoffo , 1 , 1 , 3 , 1 , 1 , 3 *
ckino:       .       i - 2 , 1 , 1d0 )
ckino:
ckino:            enddo
ckino:  C       call prmx('qp',w(oqp),3,3,nqp)
ckino:          endif
ckino:
ckino:          call info0(30,0,0,' RDSIGM: writing data to file QPTS ...')
ckino:          ifiz = fopna('qpts',-1,0)
ckino:            call getqp ( 1 , - ifiz , nqp , nkxyz , lshft , 0 , rv_p_oqp
ckino:       .    , iwdummy , iwdummy )
ckino:
ckino:
ckino:
ckino:          if (lfbzout) then
ckino:            if (allocated(ipq_iv)) deallocate(ipq_iv)
ckino:
ckino:          endif
ckino:
ckino:        endif
ckino:
ckino:        call tcx('rdsigm')
ckino:
ckino:        if (allocated(gstar_iv)) deallocate(gstar_iv)
ckino:
ckino:        end subroutine rdsigm
ckino:
