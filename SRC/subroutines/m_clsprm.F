      module m_clsprm
      integer , allocatable ::  iv_a_oipqx(:) !for lmchk only
      real(8) , allocatable ::  rv_a_ormax(:)
      real(8) , allocatable ::  rv_a_obxc (:) !lmaux clsprm
      real(8) , allocatable ::  rv_a_ogrrme (:)
      real(8) , allocatable ::  rv_a_opmpol (:)
      real(8) , allocatable ::  rv_a_opnu (:)
      real(8) , allocatable ::  rv_a_opp (:)
      real(8) , allocatable ::  rv_a_opprel (:)
      real(8) , allocatable ::  rv_a_oqnu (:)
      real(8) , allocatable ::  rv_a_oqpp (:)
      real(8) , allocatable ::  rv_a_oqt (:)
      real(8) , allocatable ::  rv_a_orhrmx (:)
      real(8) , allocatable ::  rv_a_osop (:)
      real(8) , allocatable ::  rv_a_ovdif (:)
      real(8) , allocatable ::  rv_a_oves (:)
      real(8) , allocatable ::  rv_a_ovintr (:)
      real(8) , allocatable ::  rv_a_ovrmax (:)
      contains
      subroutine clsprm(mode,sspec)!,slat)!,sbz)!,sstr) !,sstrn)sctrl,sham,spot,,sarray
      use m_struc_def           !Cgetarg
      use m_supot,only: rv_a_omad 
      use m_mksym,only: iv_a_oics,ctrl_nclass
      use m_lmfinit,only: lso,  ctrl_nl,ctrl_nspin,ctrl_nl,ctrl_nspec,ctrl_nbas,
     &     ham_neula,ham_nbf,ham_nlibu
      use m_lattic,only: rv_a_opos
C     - Allocate memory for, and try to read, class-related data
C     ----------------------------------------------------------------------
C     i    mode :1s digit
C     i         : 1 allocates ASA arrays (see spot,sham,sarray)
C     i         : 2 allocates TBE arrays:  oifrlx
C     i   sctrl :struct for program flow parameters; see routine uctrl
C     i     Elts read: nbas nclass nl nspin nspec loptc nmap
C     i     Stored:    *
C     i     Passed to: lgors ioqpp aiocls clsprp
C     i   sarray:structure containing offsets to various arrays
C     i     Elts read: nclasp oics oipc oclabl opos
C     i     Stored:    ogroup ormax ohave oifrlx ozos
C     i     Passed to: aiocls clsprp
C     i   sham  :struct for parameters defining hamiltonian; see routine uham
C     i     Elts read: *
C     i     Stored:    osop oivso
C     i     Passed to: aiocls clsprp
C     i   spot  :struct for information about the potential; see routine upot
C     i     Elts read: ovrmax orhrmx opnu oqnu opp oves
C     i     Stored:    opnu oqnu oqc oqt oves ovdif ovrmax orhrmx opp orhos
C     i                osop ogrrme ovintr opmpol oqpp
C     i     Passed to: ioqpp aiocls clsprp
C     i   sspec :struct for species-specific information; see routine uspec
C     i     Elts read: rmt idmod z lmxa
C     i     Stored:    *
C     i     Passed to: aiocls clsprp
C     i   slat  :struct for lattice information; see routine ulat
C     i     Elts read: plat
C     i     Stored:    *
C     i     Passed to: aiocls clsprp
C     i   sbz   :struct for the Brillouin Zone; see routine ubz
C     i     Elts read: ndos
C     i     Stored:    *
C     i     Passed to: clsprp
C     i
C     i   sstr  :struct for parameters for screened strux; see routine ustr
C     i     Elts read: n map
C     i     Stored:    *
C     i     Passed to: *
C     i   sstrn :struct for global strings
C     i     Elts read: map
C     i     Stored:    *
C     i     Passed to: clsprp
C     u Updates
C     u   08 Nov 07 (J. Xu) LDA+U implementation; qpp is complex
C     u   30 Sep 04 Reads/writes relativistic ppar's
C     u   02 Apr 04 SO array redimensioned; make and pack pot->nrhos
C     u   18 Mar 03 Allocate space for relativistic potential parameters
C     u   07 Feb 03 SO array redimensioned
C     u   19 Feb 02 entry clsprp copies class data to doubly padded layers
C     u   28 Apr 98 allocate space for radial matrix elements of grad
C     ----------------------------------------------------------------------
      implicit none
      integer mode
      type(s_spec)::sspec(*)
c      character*(*) sstrn
c      type(s_ctrl)::sctrl
c      type(s_array)::sarray
c      type(s_ham)::sham
c      type(s_pot)::spot
c      type(s_lat)::slat
c      type(s_bz)::sbz
c      type(s_str)::sstr
C     Local variables
      logical lgors,cmdopt
      character*80 strn,clabl*8
      integer:: ogroup , oifrlx , oivso , oipc , opos , ozos , n0 
      parameter (n0=10)
      integer nbas,nclass,nspec,nl,nsp,nlspc,igets,nclasp,i1,i2,nmap,ic,
     .     is,ib,iclbsj,idmod(n0),iprint,lmx,nclspp,nqpp,lnsph,isw,neula,
     .     nbf,nrhos,nlibu
      double precision plat(9),z,xx,dval
      integer,allocatable:: ics(:)

      nbas=ctrl_nbas
      nclass=ctrl_nclass
      nl=ctrl_nl
      nsp=ctrl_nspin
      nspec=ctrl_nspec
      nlibu=ham_nlibu
      nclasp=ctrl_nclass !sarray%nclasp
      nclspp = max(2*nclasp-nclass,nspec)
      nlspc = nl*nsp*max(nclspp,nspec)

C     ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
         allocate(rv_a_opnu(nlspc))
         rv_a_opnu(:)=0.0d0
         allocate(rv_a_oqnu(3*nlspc))
         rv_a_oqnu(:)=0.0d0
         allocate(rv_a_opp(6*nlspc))
         rv_a_opp(:)=0.0d0
         allocate(rv_a_opprel(32*nl*nl*max(nclspp,nspec)))
         rv_a_opprel(:)=0.0d0
         allocate(rv_a_oqt(nclasp))
         rv_a_oqt(:)=0.0d0
         allocate(iv_a_oipqx(nclasp))
         iv_a_oipqx(:)=0
         allocate(rv_a_oves(nclspp))
         rv_a_oves(:)=0.0d0
         allocate(rv_a_ovdif(nclspp))
         rv_a_ovdif(:)=0.0d0
         allocate(rv_a_orhrmx(nclspp))
         rv_a_orhrmx(:)=0.0d0
         allocate(rv_a_ovrmax(2*nclspp))
         rv_a_ovrmax(:)=0.0d0
         allocate( ics(nclspp) )
         call icopy ( nclspp , iv_a_oics , 1 , ics , 1 ) 
         allocate(rv_a_ormax(nclspp))
         call dcopy ( nclspp , sspec ( ics ( 1:nclspp ) ) %rmt , 1 , rv_a_ormax, 1 ) 
C     Even in collinear case, bxc might be written to
         allocate(rv_a_obxc(abs(-3*nclasp)))
         rv_a_obxc(:)=0.0d0
         nrhos=0 
         if ( lso /=0 ) then 
            neula=ham_neula
            nbf=ham_nbf
            nrhos = max(nl,neula,nbf)
            if (allocated(rv_a_obxc)) deallocate(rv_a_obxc)
            allocate(rv_a_osop(nlspc*nsp*9))
            rv_a_osop(:)=0.0d0
            allocate(rv_a_obxc(3*nclasp))
            rv_a_obxc(:)=0.0d0
            nrhos=nrhos 
         endif
         if (allocated(rv_a_ogrrme)) deallocate(rv_a_ogrrme)
c         if ( iand(3,int(ctrl_lsx)) .ne.0 .or. iand(-1,int(ctrl_lscr)) .ne.0 ) then 
c            allocate(rv_a_ovintr(abs(-nclasp*(nl*nsp)**2)))
c            rv_a_ovintr(:)= 0d0
c         endif
         lnsph = isw ( .false.)
C     ... Allocate memory for, attempt to read qpp from file
         if (lnsph .ne. 0 .or. nlibu .ne. 0) then
            allocate(rv_a_opmpol(abs(-(2*nl-1)*nl**2*3*nsp*nclasp)))
            rv_a_opmpol(:)= 0d0
            i1 = nl**2
            nqpp = (i1*(i1+1))/2
            allocate(rv_a_oqpp(2*nqpp*4*nsp*nbas))
            call dvset ( rv_a_oqpp , 1 , 1 , - 1d0 ) 
            call ioqpp(.false.)!,spot)!,sctrl
         endif
      endif
      return

      entry clsprp(mode,sspec)!,slat)!,sbz),sctrl sham,,spot,sarray
c     .     sstr) !,sstrn)
C     - Copy class data to doubly padded case
      nclass=ctrl_nclass
      nl=ctrl_nl
      nsp=ctrl_nspin
      nspec=ctrl_nspec
      nclasp=ctrl_nclass !sarray%nclasp

      if (nclasp .eq. nclass) return
      nclspp = max(2*nclasp-nclass,nspec)
C     ... Allocate ASA arrays
      if (mod(mode,10) .eq. 1) then
         i1 = nclasp-nclass
         call dpscop ( rv_a_opnu , rv_a_opnu , nl * nsp * i1 , nl * nsp 
     .        * nclass + 1 , nl * nsp * nclasp + 1 , 1d0 ) 
         call dpscop ( rv_a_oqnu , rv_a_oqnu , 3 * nl * nsp * i1 , 3 * 
     .        nl * nsp * nclass + 1 , 3 * nl * nsp * nclasp + 1 , 1d0 ) 
         call dpscop ( rv_a_oves , rv_a_oves , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
         call dpscop ( rv_a_ovdif , rv_a_ovdif , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
         call dpscop ( rv_a_ovrmax , rv_a_ovrmax , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
         call dpscop ( rv_a_orhrmx , rv_a_orhrmx , i1 , nclass + 1 , nclasp + 1 , 1d0 ) 
      else
         call rxi('clsprp: bad mode',mode)
      endif
      end subroutine clsprm 
      end module m_clsprm
