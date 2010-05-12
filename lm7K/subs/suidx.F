Cgetarg...info...           structure ['suidx', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['suidx', 'm_struc_def']
      subroutine suidx(nkaph,mode,nspec,sspec)

      use m_struc_def  !Cgetarg

C- Setup for idxdn
C ----------------------------------------------------------------------
Ci Inputs
Ci   nkaph :The maximum number of radial functions centered at
Ci         :particular R and l channel used in the lmto basis.
Ci         :This includes the number of envelope functions.
Ci         :and the number of local orbitals
Ci   mode  :1s digit
Ci         :  0 lmxb cutoff (e.g. ASA 2nd generation hamiltonian)
Ci         :    orbitals for lmxb<l<n0 are mapped to 'neglected'
Ci         :  1 lmxb and lmxa cutoffs
Ci         :    orbitals for lmxb<l<=lmxa are mapped into 'high'
Ci         :    orbitals for lmxa<l<n0 are mapped into 'neglected'
Ci         :  add 2 to replace any idxdn=0 with idxdn=1
Ci         :  NB: for floating orbitals, lmxa=-1.
Ci         :  Thus, use lmxa -> max(lmxa,lmxb)
Ci         :10s digit
Ci         :  1 lmf style hamiltonian:
Ci         :    orbitals with positive rsmh are defined as being
Ci         :    in the basis.
Ci   sspec :struct containing species-specific information
Cio Inputs/Outputs
Cio  sspec :struct for species-specific information; see routine uspec
Cio    Elts read: lmxb lmxa idxdn pz orbp
Cio    Stored:    idxdn pz
Cio        :Note: input pz serves both as a boundary condition for any
Cio        :local orbitals defined, and also the 10's digit flags how
Cio        :each local orbital is to be included in the hamiltonian.
Cio        :On output the 10s' digit is cleared.
Cio        :10s digit on input:
Cio        :0 value and slope of local orbital=0 at MT boundary
Cio        :1 a smooth Hankel tail is attached.  The orbital is included
Cio        :  in the basis
Cio        :2 a smooth Hankel tail is attached.  The orbital is not
Cio        :  included in the basis; instead the coupling between the
Cio        :  local orbital and remaining states in the system is
Cio        :  included perturbatively.
Cl Local variables
Cl    iloc :0 if basis has no local orbitals
Cl         :1 if basis has local orbitals
Cr Remarks
Cr  Each basis function is labelled by an l quantum number and a species
Cr  index (a function of this type is centered at every site corresponding
Cr  to the species index).  Also, there may be more than one kind of basis
Cr  function per site and l, but there can be at most nkaph of such kinds.
Cr  Thus the maximum possible number of function types associated with a
Cr  particular site is nkaph*lmxa.
Cr
Cr  Array idxdn(1..lmxb+1,1:nkaph) keeps track of how functions are used
Cr  in the construction of the basis.  For a particular 0<=l<=lmxb and
Cr  1<=ik<nkaph, idxdn takes one of the following:
Cr   value   Signifies
Cr     0     Role of this orbital has not yet been determined
Cr     1     Orbital is included as "active" orbitals, which means
Cr  	     they are included in the hamiltonian and diagonalized
Cr     2     Orbital is treated as an "intermediate" orbital, which
Cr  	     means it is downfolded and included in a perturbative way
Cr     3     Orbital is treated as an "high" orbital, which means
Cr  	     means it is included in augmentation for tails of other
Cr  	     orbitals, but is otherwise not part of the basis.
Cr     4     Orbital is neglected
Cr    10     Orbital is a local orbital whose value and slope are
Cr    	     constructed to be zero at the MT boundary.
Cr           It is included in the basis.
Cr    11     Orbital is a local orbital with a smooth Hankel tail
Cr           and it is included in the basis.
Cr    12     Orbital is a local orbital with a smooth Hankel tail
Cr  	     It is incorporated perturbatively to the hamiltonian
Cr  	     and is not assembled as part of the hamiltonian matrix
Cr
Cr   Routine makidx uses idxdn to construct the ordering  of orbitals
Cr   in the hamiltonian matrix.  See description of
Cr   offH and iprmb in makidx.
Cb Bugs
Cb   Perhaps iloc should be species-dependent.  Implies nkaph
Cb   channel is only local orbital for species which have them.
Cu Updates
Cu   24 Jun 04 Extended definition of local orbitals.  New arg list.
Cu   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
Cu   24 Aug 01 Extended to local orbitals.  Altered argument list.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nkaph,nspec
Cgetarg       double precision sspec(1)
      type(s_spec)::sspec(*)

C ... Local parameters
      logical lmf
      integer is,n0,nkap0,lmxb,lmxa,iloc
      parameter (n0=10,nkap0=3)
      integer idxdn(n0,nkap0),lp1,mode0,ik
      double precision orbp(n0,2,nkap0),pz(n0,2),dum

      mode0 = mod(mode,10)
      lmf = mod(mode/10,10) .eq. 1
      call isanrg(mode0,0,3,'suidx:','1s digit mode',.true.)
      call isanrg(mod(mode/10,10),0,1,'suidx:','10s digit mode',.true.)
C     Determine whether basis has any local orbitals
      call uspecb(0,0,sspec,1,nspec,dum,dum,dum,iloc)
      if (iloc/10 .ne. 0) then
        iloc = 1
      endif

      do  is = 1, nspec
Cgetarg...info...           integer :: lmxb , integer(8) :: sspec%lmxb 1
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idxdn(n0,nkap0) , integer(8) :: sspec%idxdn 30
Cgetarg...info...           real(8) :: pz(n0,2) , real(8) :: sspec%pz 20
Cgetarg       call upack('spec lmxb lmxa idxdn pz',sspec,is,lmxb,lmxa,idxdn,pz)

        lmxb=sspec(is)%lmxb
        lmxa=sspec(is)%lmxa
        i_copy_size=size(sspec(is)%idxdn)
        call i8icopy(i_copy_size,sspec(is)%idxdn,1,idxdn,1)
        i_copy_size=size(sspec(is)%pz)
        call dcopy(i_copy_size,sspec(is)%pz,1,pz,1)


        do  ik = 1, nkap0

C         Map all orbitals above max(lmxb,lmxa) into `neglected'
          call ivset(idxdn(1,ik),max(lmxb,lmxa)+2,n0,4)

C         mode0=0: Map any orbital above lmxb into `neglected'
          if (mod(mode0,2) .eq. 0) then
            call ivset(idxdn(1,ik),lmxb+2,n0,4)

C         Map any orbital between lmxb and lmxa into `high'
          else
            call ivset(idxdn(1,ik),lmxb+2,lmxa+1,3)
          endif

C         Mapping for local orbitals
          if (ik .eq. nkaph .and. iloc .eq. 1) then
            call ivset(idxdn(1,ik),1,n0,4)
            do  lp1  = 1, lmxb+1

              if (pz(lp1,1) .ne.  0) then
                call pz2idx(1,pz(lp1,1),idxdn(lp1,ik))
C              if (pz(lp1,1) .ne.  0) idxdn(lp1,ik) = 10
C              if (pz(lp1,1) .gt. 10) idxdn(lp1,ik) = 11
C              if (pz(lp1,1) .gt. 20) idxdn(lp1,ik) = 12
C              if (pz(lp1,1) .gt. 10) pz(lp1,1)=pz(lp1,1)-10
C              if (pz(lp1,1) .gt. 10) pz(lp1,1)=pz(lp1,1)-10
C              if (pz(lp1,1) .gt. 10)
C     .          call rx1('suidx: bad value for pz',pz(lp1,1)+20)
              endif
            enddo
          endif

C         Map all orbitals for ik>nkaph into `neglected'
          if (ik .gt. nkaph) call ivset(idxdn(1,ik),1,n0,4)

C     ... Automatic downfolding turned off ... set idxdn=0 to idxdn=1
          if (mode0 .ge. 2) then
            do  lp1  = 1, n0
              if (idxdn(lp1,ik) .eq. 0) idxdn(lp1,ik) = 1
            enddo
          endif

C     ... lmf-specific
          if (lmf .and. ik .le. nkaph-iloc) then
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: sspec%orbp 60
Cgetarg             call upack2('spec orbp',sspec,is,orbp)

            i_copy_size=size(sspec(is)%orbp)
            call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1)

            do  lp1 = 1, lmxb+1
              if (orbp(lp1,1,ik) .lt. 0) idxdn(lp1,ik) = 3
C             if (orbp(lp1,1,ik) .eq. 0) idxdn(lp1,ik) = max(idxdn(lp1,ik),2)
              if (orbp(lp1,1,ik) .eq. 0) idxdn(lp1,ik) = 3
            enddo
          endif

        enddo

Cgetarg...info...           integer :: idxdn(n0,nkap0) , integer(8) :: sspec%idxdn 30
Cgetarg         call pack2('spec idxdn',sspec,is,idxdn)

        i_copy_size=size(sspec(is)%idxdn)
        call ii8copy(i_copy_size,idxdn,1,sspec(is)%idxdn,1)

      enddo

Cgetarg       end
      end subroutine suidx


      subroutine pz2idx(mode,pz,idxdn)

C- Converts information about local orbital between pz and idxdn
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing; just return
Ci         :1 use 10's digit of pz to set idxdn
Ci         :2 use idxdn to set 10's digit of pz
Ci         :Add 10 clear 10's digit of pz before exit
Cio Inputs/Outputs
Cio   pz   :10s digit of pz contains the following information
Cio        :0 value and slope of local orbital=0 at MT boundary
Cio        :1 a smooth Hankel tail is attached.  The orbital is included
Cio        :  in the basis
Cio        :2 a smooth Hankel tail is attached.  The orbital is not
Cio        :  included in the basis; instead the coupling between the
Cio        :  local orbital and remaining states in the system is
Cio        :  included perturbatively.
Cio        :mode 2 or 12: pz input: sets idxdn from 10's digit of pz
Cio        :mode 1 pz output: 10's digit of pz from idxdn
Cio        :mode >=10 10's digit of pz cleared on exit
Cio   idxdn:contains the same kind information as the 10s digit of pz.
Cio        :<10 equivalent to pz=0
Cio        :10  equivalent to 10's digit pz=0
Cio        :11  equivalent to 10's digit pz=1
Cio        :12  equivalent to 10's digit pz=2
Cio        :mode 2 or 12: idxdn output: sets idxdn from 10's digit of pz
Cio        :mode 1 idxdn input: 10's digit of pz from idxdn
Cr Remarks
Cr   idxdn keeps information about how local orbitals are used because
Cr    it controls ordering of hamiltonian orbitals.  See makidx.
Cu Updates
Cu   25 Jun 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,idxdn
      double precision pz
C ... Local parameters

      call isanrg(mod(mode,10),0,2,'pz2idx','1s digit mode',.true.)
      call isanrg(mode/10,0,1,'pz2idx','10s digit mode',.true.)

      if (mod(mode,10) .eq. 1) then
        call fsanrg(pz,0d0,28d0,0d0,'pz2idx','pz',.true.)
        if (pz .ne.  0) idxdn = 10
        if (pz .ge. 10) idxdn = 11
        if (pz .ge. 20) idxdn = 12
      elseif (mod(mode,10) .eq. 2) then
        call isanrg(mod(idxdn,10),0,2,'pz2idx','1s digit idxdn',.true.)
        call isanrg(idxdn/10,0,1,'pz2idx','10s digit idxdn',.true.)
        if (idxdn .lt. 10) pz =  0
        if (idxdn .eq. 10) pz = mod(pz,10d0) +  0
        if (idxdn .eq. 11) pz = mod(pz,10d0) + 10
        if (idxdn .eq. 12) pz = mod(pz,10d0) + 20
      endif

      if (mode .ge. 10) pz = mod(pz,10d0)

Cgetarg       end
      end subroutine pz2idx


