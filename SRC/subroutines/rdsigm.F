c$$$      module m_sigflg
c$$$      contains
c$$$      subroutine sigflg(dc,sio,strn,ledit,lwsig,lrot,lfbz,lnwmsh,lonesp,
c$$$     .qpoff,rots)
c$$$
c$$$C- Reads flags for Sigma IO
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   dc:   separator for switches
c$$$Ci  sio:   a string, which must be one of :
c$$$Ci         'rsig' --- parse parameters for reading sigma:
c$$$Ci                lfbz qpoff
c$$$Ci         'wsig' --- parse parameters for writing a new sigma:
c$$$Ci                ledit lwsig lrot lfbz lnwmsh lonesp qpoff rots
c$$$Ci   strn: string containing switches
c$$$Co Outputs
c$$$Co   ledit :(wsig only)
c$$$Co         :1  flag requesting sigma editor
c$$$Co   lwsig :(wsig only)
c$$$Co         :0 normal exit
c$$$Co         :1  transformation of sigm from orbital to LDA basis
c$$$Co         :-1 transformation of from LDA basis to orbital basis
c$$$Co         :2  Same as 1, but high,low parts of sigm replaced by diagonal
c$$$Co         :3  Write evals,evecs of LDA hamiltonian to file
c$$$Co         :4  Write evals,evecs of total hamiltonian to file
c$$$Co         : OLD; not used
c$$$Co         :4  Make QP hamiltonian from Z^LDA and U with E(qp)
c$$$Co         :5  Make orthogonalized QP hamiltonian from Z^LDA and U with E(qp)
c$$$Co   lrot  :(wsig only)
c$$$Co         :nonzero, flag to rotate sigma matrix
c$$$Co   lfbz  :(rdsig and wsig)
c$$$Co         :T, write sigma for qp in entire BZ
c$$$Co   lnwmsh:(wsig only)
c$$$Co         :T, write sigm for qp on new mesh
c$$$Co   lonesp:(wsig only)
c$$$Co         :nonzero, average up and down spins before writing
c$$$Co   qpoff :(rdsig and wsig)
c$$$Co         :k-point offset to be added to uniform mesh
c$$$Co   rots  :(wsig only)
c$$$Co         :string containing rotation matrix for sigma
c$$$Cl Local variables
c$$$Cl         :
c$$$Cr Remarks
c$$$Cr
c$$$Cu Updates
c$$$Cu   21 Jun 08
c$$$Cu   15 Jan 07  First created
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character dc*1,sio*(4)
c$$$      character strn*(*)
c$$$      character,optional::rots*(*)
c$$$      integer,optional:: ledit,lwsig,lrot
c$$$      integer::lonesp
c$$$      logical lfbz
c$$$      logical,optional::lnwmsh
c$$$      double precision qpoff(3)
c$$$C ... Local parameters
c$$$      integer parg
c$$$      integer i,j,j1,j2,iv(10)
c$$$      logical lio
c$$$
c$$$      logical:: isanrg, l_dummy_isanrg
c$$$
c$$$      lio = sio .eq. 'wsig'
c$$$      if (.not. lio) then
c$$$        if (sio .ne. 'rsig') call rx('sigflg: bad argument '//sio)
c$$$      endif
c$$$      if (lio) then
c$$$        ledit = 0
c$$$        lwsig = 0
c$$$        lrot = 0
c$$$        lfbz = .false.
c$$$        lnwmsh = .false.
c$$$        lonesp = 0
c$$$        rots = ' '
c$$$        call dpzero(qpoff,3)
c$$$      else
c$$$        lfbz = .false.
c$$$        lonesp = 0
c$$$        call dpzero(qpoff,3)
c$$$      endif
c$$$
c$$$      if (dc .ne. ' ') then
c$$$C   ... Return here to resume parsing for arguments
c$$$        j2 = 0
c$$$   50   continue
c$$$        j2 = j2+1
c$$$        if (strn(j2:j2) .eq. dc) goto 50
c$$$        j1 = min(len(strn),j2)
c$$$        call nwordg(strn,0,dc//' ',1,j1,j2)
c$$$        if (j2 .ge. j1) then
c$$$          if (.false.) then
c$$$          elseif (strn(j1:j1+3) .eq. 'edit'.and. lio) then
c$$$            ledit = 1
c$$$          elseif (strn(j1:j1+4) .eq. 'onesp') then
c$$$            lonesp = 1
c$$$          elseif ((strn(j1:j1+3) .eq. 'rot=' .or.
c$$$     .    strn(j1:j1+3) .eq. 'ROT=') .and. lio) then
c$$$            rots = strn(j1+4:j2)
c$$$          elseif (strn(j1:j1+5) .eq. 'trans=' .and. lio) then
c$$$            j = 0
c$$$            if (parg('trans=',2,strn(j1:),j,len(strn(j1:)),
c$$$     .      dc//' ',1,1,i,lwsig) .lt. 0) goto 999
c$$$Ckino isanrg is logical function,             call isanrg(mod(lwsig,10),-1,5,'SIGFLG:','trans',.true.)
c$$$            l_dummy_isanrg=isanrg(mod(lwsig,10),-1,5,'SIGFLG:','trans',.true.)
c$$$            if (lwsig .ne. 0 .and. mod(lwsig,10) .eq. 0 .or.
c$$$     .      lwsig .ge. 40)
c$$$     .      call rxi('rdsigm: nonsensical trans',lwsig)
c$$$          elseif (strn(j1:j1+4) .eq. 'phase'.and. lio) then
c$$$            lrot = 2
c$$$          elseif (strn(j1:j1+2) .eq. 'fbz') then
c$$$            lfbz = .true.
c$$$          elseif (strn(j1:j1+4) .eq. 'newkp'.and. lio) then
c$$$            lnwmsh = .true.
c$$$          elseif (strn(j1:j1+5) .eq. 'shftq=') then
c$$$            j = 0
c$$$            if (parg('shftq=',4,strn(j1:),j,len(strn(j1:)),
c$$$     .      ', '//dc,2,3,iv,qpoff) .lt. 0) goto 999
c$$$          elseif (strn(j1:j1+4) .eq. 'shftq') then
c$$$            qpoff(1) = 4d0/4096
c$$$            qpoff(2) = 5d0/4096
c$$$            qpoff(3) = 7d0/4096
c$$$          else
c$$$            goto 59
c$$$          endif
c$$$          goto 50
c$$$   59     continue
c$$$  999     call rxs('rdsigm: failed to parse wsig option', strn)
c$$$        endif
c$$$      endif
c$$$      end subroutine sigflg
c$$$      end module m_sigflg




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
