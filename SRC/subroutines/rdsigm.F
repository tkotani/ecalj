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

      subroutine iosigh(optio,mode,nsp,ndimh,nk1,nk2,nk3,nqp, lshft1,lshft2,lshft3,ifi)
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
      integer nsp,ndimh,nk1,nk2,nk3,nqp,lshft1,lshft2,lshft3,ifi,optio, mode
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
      use m_lmfinit,only: globalvariables
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

      integer:: lb1,lb2 !mar2014takao
      logical:: laf !jun2015takao

      call getpr(ipr)
      laf=allocated(symops_af) !jun2015takao
      stdo = lgunit(1)
      nl = globalvariables%nl
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

