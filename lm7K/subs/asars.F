      subroutine asars(mode,sctrl,ssite,sspec,slat,spot,sbz,
     .pnu,qnu,lbin,ifi)

      use m_struc_def  !Cgetarg

C- ASA file I/O to restart file
C ----------------------------------------------------------------------
Ci Inputs/Outputs
Ci   mode  :not used now
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: iors
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: pnu
Ci     Stored:    class pnu
Ci     Passed to: spackv asars1 iors
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: iors
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: iors
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: iors
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: iors
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
Cio Inputs/Outputs
Cio  pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Cio         pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Cio         class-based array; I/O as site-based
Cio  qnu   :energy-weighted moments of the sphere charges
Cio         class-based array; I/O as site-based
Cr Remarks
Cr   This is a bridging code to map class-resolved ASA parameters
Cr   into format suitable for iors.
Cl Local variables
Cl   lrs switches:  0=>use rst file data; 1=>ignore rst file data
Cl   lrs(1) site positions
Cl   lrs(2) starting fermi level
Cl   lrs(3) starting pnu's
Cu Updates
Cu   07 Jul 05 version -1.04 ... adapted to altered iors.
Cu   28 Feb 02 version -1.03 ... adapted to altered iors.
Cu   10 May 01 version -1.  Only writes pnu,qnu
C  ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lbin
      integer mode,ifi
      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_pot)::spot
      type(s_bz)::sbz

      double precision pnu(*),qnu(*)
C ... Local parameters
      character*32 jobid
      integer:: n0 , nbas , nl , nsp , nglob , i , iors , nit , nspec 
     ., nat
      integer ,allocatable :: ov_iv(:)
      integer ,allocatable :: ipc_iv(:)
      integer ,allocatable :: ipa_iv(:)
      real(8) ,allocatable :: pnu_rv(:)
      real(8) ,allocatable :: qnu_rv(:)

      parameter (n0=10)
C ... Heap
c      integer w(1)
c      common /w/ w

      jobid = 'ASA output'
      nbas = nglob('nbasp')
      nsp = nglob('nsp')
      nl = nglob('nl')
      nspec = nglob('nspec')
      nit = 1

C ... Get class list
      allocate(ipc_iv(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) ipc_iv(:)=0

      allocate(ipa_iv(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) ipa_iv(:)=0

      i_copy_size=1;
      do i_spackv=1,nbas
        call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %class 
     .  , i_copy_size , i_spackv + 1 - 1 , ipc_iv )

      enddo


C ... Memory for qnu, site resolved
      allocate(qnu_rv(abs(-3*n0*nsp*nbas)))
      if (-3*n0*nsp*nbas<0) qnu_rv(:)=0.0d0

      allocate(pnu_rv(abs(-n0*nsp*nbas)))
      if (-n0*nsp*nbas<0) pnu_rv(:)=0.0d0


C ... For now, just a placeholder.
      allocate(ov_iv(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) ov_iv(:)=0


C ... Poke class-based P,Q into site->pnu and oqnu
      if ( ifi .lt. 0 ) call asars1 ( 10 , ssite , nbas , nsp , nl 
     ., ov_iv , pnu_rv , qnu_rv , ipc_iv , ipa_iv , 
     .pnu , qnu )


C ... File I/O
C     for now
      nat = nbas
      i = iors ( 0 , sctrl , ssite , sspec , slat , spot , sbz , jobid 
     ., nbas , nat , nspec , w , qnu_rv , nit , lbin , ifi )


C ... Poke site struc and oqnu into P,Q
      if ( ifi .gt. 0 ) call asars1 ( 11 , ssite , nbas , nsp , nl 
     ., ov_iv , pnu_rv , qnu_rv , ipc_iv , ipa_iv , 
     .pnu , qnu )


      if (allocated(ov_iv)) deallocate(ov_iv)
      if (allocated(pnu_rv)) deallocate(pnu_rv)
      if (allocated(qnu_rv)) deallocate(qnu_rv)
      if (allocated(ipa_iv)) deallocate(ipa_iv)
      if (allocated(ipc_iv)) deallocate(ipc_iv)

      end subroutine asars

      subroutine asars1(mode,ssite,nbas,nsp,nl,ov,pnu,qnu,ipc,ipa,pl,ql)

      use m_struc_def  !Cgetarg

C- Kernel that copies ASA parameters into form readable by iors
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit :
Ci         : 0  Copies classed-based P,Q to site-based P,Q
Ci         : 1  Symmetrizes site-based P,Q and copy to class-based P,Q
Ci         :10s digit specifies how site P is packed/unpacked
Ci         : 1  Store or retrieve site P to/from site->pnu
Ci              Note: In storage mode, site->pnu is also copied to
Ci              array pnu as part of the symmetrization procedure.
Ci              Symmetrized pnu is repacked into site->pnu
Ci         : 2  Store or retrieve site P to/from array pnu
Ci         : 3  In class-to-site copy: combination of 1+2
Ci         :    In site-to-class copy, same as 2, but also
Ci         :    symmetrized pnu is packed into site->pnu.
Ci   mode  :0 copy ASA pl,ql to site-based pnu,qnu and also site->pnu
Ci         :1 symmetrize and reverse copy
Ci   ssite :struct containing site-specific information
Ci     Elts read: pnu (mode=1)
Ci     Stored:    pnu (mode=0)
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   ov    :not used now
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ipa   :integer work array (holds list of class members)
Cio Inputs/Outputs
Cio   pl   :sphere boundary conditions, by class
Cio        :Input if mode=0; output if mode=1
Cio   ql   :energy-weighted moments of the sphere charges, by class
Cio        :Input if mode=0; output if mode=1
Cio  pnu   :same as pl, but resolved by site (also different dim.)
Cio  qnu   :same as ql, but resolved by site (also different dim.)
Cio        :Output if mode=0; input if mode=1
Cb Bugs
Cb   ssite shouldn't be needed here ... should move copy to pnu out.
Cu Updates
Cu    1 Mar 02  Reworked mode options.
Cu   10 May 01  First attempt.  Writes P,Q only
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,mode,nl,nsp,ov(nbas),ipc(nbas),ipa(7),n0
      parameter (n0=10)
      real(8):: pl(nl,nsp,*) , ql(3,nl,nsp,*)
      type(s_site)::ssite(*)

      double precision pnu(n0,nsp,*),qnu(3,n0,nsp,*)
C ... Local parameters
      integer ib,ic,nrclas,mode0,mode1
C     character*8 clabl
      double precision xx,wk(3*n0*2),ploc(n0,2)
C     logical aiopot

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

      call isanrg(mode0,0,1,'asars1:','1s digit mode',.true.)
      call isanrg(mode1,1,3,'asars1:','10s digit mode',.true.)

C --- Create symmetrized pnu ---
      if (mode0 .eq. 1) then

C   ... Extract site->pnu into pnu array
        if (mode1 .eq. 1) then
          do  ib = 1, nbas

            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,ssite(ib)%pnu,1,ploc,1)

            call dcopy(n0*nsp,ploc,1,pnu(1,1,ib),1)
          enddo
        endif

C   ... Symmetrize pnu,qnu
        ic = 0
   10   continue
        ic = ic+1
        call psymr0(-2,-ic,nbas,ipc,xx,xx,ipa,nrclas)
        if (nrclas .gt. 0) then
          call psymq0(nrclas,nsp,ipa,wk,3*n0,qnu)
          call psymq0(nrclas,nsp,ipa,wk,n0,pnu)
          goto 10
        endif

C   ... Copy pnu to site->pnu
        if (mode1 .eq. 1 .or. mode1 .eq. 3) then
          call dpzero(ploc,n0*2)
          do  ib = 1, nbas
            call dcopy(n0*nsp,pnu(1,1,ib),1,ploc,1)

            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,ploc,1,ssite(ib)%pnu,1)

          enddo
        endif
      endif

C --- Poke P,Q (class) to/from P,Q (site) ---
      do  ib = 1, nbas
        ic = ipc(ib)

C   ... Pack pl into site->pnu and/or array pnu; and ql to qnu
        if (mode0 .eq. 0) then
          call dpzero(ploc,n0*2)
          call dcopy(nl,pl(1,1,ic),1,ploc,1)
          call dcopy(nl,pl(1,nsp,ic),1,ploc(1,nsp),1)
          if ( mod ( mode1 , 2 ) .eq. 1 ) then

            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,ploc,1,ssite(ib)%pnu,1)
          endif

          if (mode1 .ge. 2) then
            call dpzero(pnu(1,1,ib),n0*nsp)
            call dcopy(1*nl,pl(1,1,ic),  1,pnu(1,1,ib),  1)
            call dcopy(1*nl,pl(1,nsp,ic),1,pnu(1,nsp,ib),1)
          endif
          call dpzero(qnu(1,1,1,ib),3*n0*nsp)
          call dcopy(3*nl,ql(1,1,1,ic),  1,qnu(1,1,1,ib),  1)
          call dcopy(3*nl,ql(1,1,nsp,ic),1,qnu(1,1,nsp,ib),1)

C         call spacks(0,'site clabel',ssite,clabl,ib,ib)
C         ifi = fopna(clabl,30,0)
C         call defdr(ov(ib),   nr*nsp)
C         line = 'pot'
C         if (.not. aiopot(nr,nsp,a,rmx,-99d0,w(ov(ib)),ifi)) goto 999

C   ... Copy pnu to pl and qnu to ql
        else
C         call upack2('site pnu',ssite,ib,pnu(1,1,ib))
          call dcopy(nl,pnu(1,1,ib),1,pl(1,1,ic),1)
          call dcopy(nl,pnu(1,nsp,ib),1,pl(1,nsp,ic),1)
          call dcopy(3*nl,qnu(1,1,1,ib),  1,ql(1,1,1,ic),  1)
          call dcopy(3*nl,qnu(1,1,nsp,ib),1,ql(1,1,nsp,ic),1)
        endif

      enddo

      return

C ... Error exit
C  999 continue
C      call rx('iorsa, class'//clabl//': missing',line)


      end subroutine asars1


