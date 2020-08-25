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
      use m_globalvariables
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

