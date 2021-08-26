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
      character strn*(*)
      character,optional::rots*(*)
      integer,optional:: ledit,lwsig,lrot
      integer::lonesp
      logical lfbz
      logical,optional::lnwmsh
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
      implicit none
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



c$$$      subroutine divctakao(a,b,n)
c$$$      integer:: n,i,j
c$$$      complex(8):: a(n,n),b(n,n)
c$$$      print *,' sum check 1=',sum(abs(a))
c$$$      print *,' sum check 1=',sum(abs(b))
c$$$      do i=1,n
c$$$        do j=1,n
c$$$c        write(6,"(2i5,2x, 4d13.4,3x,2d13.5)")
c$$$c     &  i,j,atan(dimag(a(i,j)/dreal(a(i,j))))       if( abs(a(i,j))>1d-3 )
c$$$c     &   write(6,"(2i5,2x, 4d13.4,3x,2d13.5)")i,j,a(i,j)
c$$$c     &   ,atan(dimag(a(i,j)/dreal(a(i,j))))
c$$$c,b(i,j),
c$$$c     &  a(i,j)/( b(i,j)+b(i,j)*0.01/abs(b(i,j)) )
c$$$c        a(i,j) = a(i,j)/ ( b(i,j)+b(i,j)*0.01/abs(b(i,j)) )
c$$$c        a(i,j)=a(i,j)/b(i,j)
c$$$        enddo
c$$$      enddo
c$$$      print *,' sum check 2=',sum(abs(a))
c$$$      end subroutine divctakao
c$$$
