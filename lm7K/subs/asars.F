Cdelw1 var ['oipc', '0', 'oipc', 'nbas', 'integer', 'iv_w_', 'asars']
Cdelw1 var ['oipa', '0', 'oipa', 'nbas', 'integer', 'iv_w_', 'asars']
Cdelw1 var ['oqnu', '0', 'oqnu', '-3*n0*nsp*nbas', 'real(8)', 'rv_w_', 'asars']
Cdelw1 var ['opnu', '0', 'opnu', '-n0*nsp*nbas', 'real(8)', 'rv_w_', 'asars']
Cdelw1 var ['oov', '0', 'oov', 'nbas', 'integer', 'iv_w_', 'asars']
Cgetarg...info...           structure ['asars', 'site', 'ssite']
Cgetarg...info...           structure ['asars1', 'site', 'ssite']
Cgetarg...info...           use_to_add ['asars', 'm_struc_def']
Cgetarg...info...           use_to_add ['asars1', 'm_struc_def']
      subroutine asars(mode,sctrl,ssite,sspec,slat,spot,sbz,
     .  pnu,qnu,lbin,ifi)
       
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
Cgetarg       double precision sctrl(1),ssite(1),sspec(1),slat(1),
Cgetarg      .  spot(1),sbz(1)
       type(s_ctrl)::sctrl
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_pot)::spot
       type(s_bz)::sbz

      double precision pnu(*),qnu(*)
C ... Local parameters
      character*32 jobid
Cdelw1       integer n0,nbas,nl,nsp,nglob,oov,i,iors,oipc,oipa,opnu,oqnu,nit,
Cdelw1      .  nspec,nat
       integer:: n0 , nbas , nl , nsp , nglob , i , iors , nit , nspec 
     . , nat 
       integer ,allocatable :: iv_w_oov(:)
       integer ,allocatable :: iv_w_oipc(:)
       integer ,allocatable :: iv_w_oipa(:)
       real(8) ,allocatable :: rv_w_opnu(:)
       real(8) ,allocatable :: rv_w_oqnu(:)

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
Cdelw1       call defi(oipc,nbas)
       allocate(iv_w_oipc(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_oipc(:)=0

Cdelw1       call defi(oipa,nbas)
       allocate(iv_w_oipa(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_oipa(:)=0

Cgetarg       call spackv(10,'site class',ssite,1,nbas,w(oipc))
       i_copy_size=1; 
       do i_spackv=1,nbas 
Cdelw1        call spackv_array_copy_i8_i('u',ssite(i_spackv)%class,i_copy_size,i_spackv+1-1,w(oipc))
        call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %class 
     .  , i_copy_size , i_spackv + 1 - 1 , iv_w_oipc ) 

       enddo


C ... Memory for qnu, site resolved
Cdelw1       call defrr(oqnu,-3*n0*nsp*nbas)
       allocate(rv_w_oqnu(abs(-3*n0*nsp*nbas)))
       if (-3*n0*nsp*nbas<0) rv_w_oqnu(:)=0.0d0

Cdelw1       call defrr(opnu,-n0*nsp*nbas)
       allocate(rv_w_opnu(abs(-n0*nsp*nbas)))
       if (-n0*nsp*nbas<0) rv_w_opnu(:)=0.0d0


C ... For now, just a placeholder.
Cdelw1       call defi(oov,nbas)
       allocate(iv_w_oov(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_oov(:)=0


C ... Poke class-based P,Q into site->pnu and oqnu
Cdelw1       if (ifi .lt. 0) call asars1(10,ssite,nbas,nsp,nl,w(oov),w(opnu),
Cdelw1      .  w(oqnu),w(oipc),w(oipa),pnu,qnu)
       if ( ifi .lt. 0 ) call asars1 ( 10 , ssite , nbas , nsp , nl 
     . , iv_w_oov , rv_w_opnu , rv_w_oqnu , iv_w_oipc , iv_w_oipa , 
     . pnu , qnu ) 


C ... File I/O
C     for now
      nat = nbas
Cdelw1       i = iors(0,sctrl,ssite,sspec,slat,spot,sbz,jobid,nbas,nat,nspec,w,
Cdelw1      .  w(oqnu),nit,lbin,ifi)
       i = iors ( 0 , sctrl , ssite , sspec , slat , spot , sbz , jobid 
     . , nbas , nat , nspec , w , rv_w_oqnu , nit , lbin , ifi ) 


C ... Poke site struc and oqnu into P,Q
Cdelw1       if (ifi .gt. 0) call asars1(11,ssite,nbas,nsp,nl,w(oov),w(opnu),
Cdelw1      .  w(oqnu),w(oipc),w(oipa),pnu,qnu)
       if ( ifi .gt. 0 ) call asars1 ( 11 , ssite , nbas , nsp , nl 
     . , iv_w_oov , rv_w_opnu , rv_w_oqnu , iv_w_oipc , iv_w_oipa , 
     . pnu , qnu ) 


Cdelw1 rlse name= oipc old_list= oipc oipa oqnu opnu oov 
Cdelw1 rlse name= oipc new_list= (None)
Cdelw1       call rlse(oipc)
       if (allocated(iv_w_oov)) deallocate(iv_w_oov)
       if (allocated(rv_w_opnu)) deallocate(rv_w_opnu)
       if (allocated(rv_w_oqnu)) deallocate(rv_w_oqnu)
       if (allocated(iv_w_oipa)) deallocate(iv_w_oipa)
       if (allocated(iv_w_oipc)) deallocate(iv_w_oipc)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cgetarg       double precision ssite(1),pl(nl,nsp,*),ql(3,nl,nsp,*)
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
Cgetarg...info...           real(8) :: ploc(n0,2) , real(8) :: ssite%pnu 20
Cgetarg             call upack2('site pnu',ssite,ib,ploc)
             
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
Cgetarg...info...           real(8) :: ploc(n0,2) , real(8) :: ssite%pnu 20
Cgetarg             call pack2('site pnu',ssite,ib,ploc)
             
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
Cgetarg...info...           real(8) :: ploc(n0,2) , real(8) :: ssite%pnu 20
Cgetarg           if (mod(mode1,2) .eq. 1) call pack2('site pnu',ssite,ib,ploc)
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

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

C ... Error exit
C  999 continue
C      call rx('iorsa, class'//clabl//': missing',line)


Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine asars1 


