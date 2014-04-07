CSFCPP#define F90 1

#define USEOPTIONAL 1

      module m_phmbls
      contains
      subroutine phmbls(mode,ndimh,eval,iprm,wk,h,zd,z,zhz)

C- Make (zd+)*h*z; order in ascending order of eigenvalues
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :consists of a sequence of bits that prescribed independent
Ci         :functions
Ci         :0  =mode0   do nothing; just return
Ci         :1  =mode1   make zhz = zd+ * h * z
Ci         :2  =mode2   poke diagonal of zhz into eval
Ci         :4  =mode4   permute eval in ascending eval order
Ci         :8  =mode8   permute z in ascending eval order
Ci         :16 =mode16  permute zhz in ascending eval order
Ci         :32 =mode32  Overwrite starting z with z^-1
Ci                      Thus, e.g. mode1 => (zd+)^-1 * h * (z)^-1
Ci         :64 =mode64  make zhz = zd+ * e * z
Ci                      If used in conjunction w/ mode=32, reconstructs
Ci                      hamiltonian from evals,evecs.
Ci   ndimh :dimension of zhz, zd, z
Ci   wk    :double complex work array of length ndimh**2
Ci   h     :hamiltonian
Ci         :Not used unless mode1 set
Cio Inputs/Outputs
Cio  iprm  :integer work of length ndimh.
Cio        :Space used if any permutations done (mode4,mode8,mode16)
Cio        :On output, holds permutations that order z and zhz in
Cio        :ascending eval order.
Cio  eval  :eigenvalues.
Cio        :If mode2, eval generated from diag (zd+ h z)
Cio        :Otherwise eval is input
Cio  zd    :(mode1,mode64 only) left evec to generate (zd+ h z).
Cio        :Usually zd and z point to the same address space;
Cio        :however, they may be different, e.g. in constructing
Cio        :off-diagonal spin blocks from spin diagonal z's
Cio  z     :approximate eigenvectors
Cio        :If mode16, z are permuted by ascending eval order
Cio  zhz   :zd+ * h * z, where z are evecs, h is hamiltonian
Cio        :If mode1, zhz is generated from z
Cio        :Otherwise zhz is input
Cio        :If mode8, output zhz is permuted by ascending eval order
Cio        :NB: zhz can occupy the same address space as h
Cio        :in which case h will be overwitten
Cr Remarks
Cr
Cu Updates
Cu   16 Jan 07 New mode 32
Cu   19 May 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#ifdef USEOPTIONAL
      integer mode,ndimh
      integer,optional::iprm(ndimh)
#else
      integer mode,ndimh,iprm(ndimh)
#endif
      double precision eval(ndimh)
#ifdef USEOPTIONAL
      complex(8),optional:: z(ndimh,ndimh),zd(ndimh,ndimh)
      complex(8)::zhz(ndimh,ndimh)
      complex(8),optional:: wk(ndimh,ndimh)
      complex(8),optional::h(ndimh,ndimh)
#else
      double complex z(ndimh,ndimh),zd(ndimh,ndimh),zhz(ndimh,ndimh)
      double complex wk(ndimh,ndimh),h(ndimh,ndimh)
#endif
C ... Local parameters
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      integer i,j,k
CSFCPP#if F90
      integer ipiv(ndimh)
CSFCPP#elif AUTO_ARRAY
CSFCPP      integer ipiv(ndimh)
CSFCPP#else
CSFCPP      integer ipiv
CSFCPP      call rx('phmbls: requires compiler with allocatable arrays')
CSFCPP#endif

      if (mode .eq. 0) return

C ... Invert z
      if (mod(mode/32,2) .ne. 0) then
        call zgetrf(ndimh,ndimh,z,ndimh,ipiv,j)
        call zgetri(ndimh,z,ndimh,ipiv,wk,ndimh**2,j)
C       call zprm('z^-1',2,z,ndimh,ndimh,ndimh)
      endif

C ... zhz <- zd+ * h * z
      if (mod(mode,2) .eq. 1) then
        call zgemm('C','N',ndimh,ndimh,ndimh,one,zd,ndimh,h,ndimh,
     .  zer,wk,ndimh)
        call zgemm('N','N',ndimh,ndimh,ndimh,one,wk,ndimh,z,ndimh,
     .  zer,zhz,ndimh)
C       call zprm('z+ h z',2,zhz,ndimh,ndimh,ndimh)
      endif

C ... Eigenvalue estimates
      if (mod(mode/2,2) .eq. 1) then
        do  i = 1, ndimh
          eval(i) = zhz(i,i)
        enddo
      endif

C ... iprm sorts eigenvalues in ascending order
      if (mod(mode/4,8) .ne. 0) then
        call dvheap(1,ndimh,eval,iprm,0d0,1)
      endif

C ... Permute evals
      if (mod(mode/4,8) .ne. 0) then
        call dvprm(1,ndimh,eval,wk,iprm,.true.)
      endif

C ... Permute eigenvectors
      if (mod(mode/8,2) .ne. 0) then
        call zcopy(ndimh**2,z,1,wk,1)
        do  i = 1, ndimh
          k = iprm(i)
          do  j = 1, ndimh
            z(j,i) = wk(j,k)
          enddo
        enddo
      endif

C ... Permute zhz
      if (mod(mode/16,2) .ne. 0) then
        call zcopy(ndimh**2,zhz,1,wk,1)
        do  i = 1, ndimh
          k = iprm(i)
          do  j = 1, ndimh
            zhz(j,i) = wk(iprm(j),k)
          enddo
        enddo
      endif

C ... zhz <- zd+ * e * z ... usu. used in conjunction w/ mode=32
      if (mod(mode/64,2) .ne. 0) then
        do  j = 1, ndimh
          do  i = 1, ndimh
            wk(i,j) = eval(i)*z(i,j)
          enddo
        enddo
C       call zprm('e * z',2,wk,ndimh,ndimh,ndimh)
        call zgemm('C','N',ndimh,ndimh,ndimh,one,zd,ndimh,wk,ndimh,
     .  zer,zhz,ndimh)
C       call zprm('z+ h z',2,zhz,ndimh,ndimh,ndimh)

      endif

      end subroutine phmbls
      end module m_phmbls


      subroutine hambls ( mode , nbas , ssite , sspec , slat , sham 
     .     , isp , ispc , q , k1 , k2 , k3 , qsig , nqsig , smpot , vconst 
     .     , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp , lso , alfa , ndimh 
     .     , napw , igvapw , h , s ,isize_h, ismidb  ,ndimsig)
      use m_struc_def           !Cgetarg
      use m_globalvariables
      use m_psigd
      use m_phmbls
!! NOTE mode is returned to be -1 when diagonalized within this routine.
!!
!! Make Hamiltonian in LDA case.
!! But in QSGW case lrsig/=0, this routine makes a interpolation 
!! (filter out noise) as for sigma-vxc, and do diagonalization of LDA+(sigma-vxc).
!! Confusing... Takao like to simplify this method (see newsigmasw() branch ---it is on the way half).
!!
!!
C     - Hamiltonian and overlap of LDA+(sigma-vxc) for one k-point.
C     ----------------------------------------------------------------------
C     i Inputs
cccc  mode is output. When sigma mode, mode=-1 is returned.
C     i   mode  :1s digit specifies LDA matrix elements
C     i         :  0 compute LDA hamiltonian and overlap
C     i         :    If the self-energy is not to be added, return with h,s
C     i         :  1 Assume h,s already computed.
C     i         :    Self-energy added to h, if 10s digit set
C     i         :10s digit specifies self-energy matrix elements
C     i         :  0 Self-energy is not added
C     i         :  1 Add Bloch sum of r.s. sigma; r.s. sigma is real
C     i         :  2 Add Bloch sum of r.s. sigma; r.s. sigma is complex
C     i         :100s digit specifies interpolation mode
C     i         :  0 Simple bloch sum of sigma
C     i         :  1 approx high- and low-energy sigma with diagonal
C     i         :    See Remarks
C     i         :  2 perturbation approach, diagonalizing exactly only
C     i         :    the first kappa block
C     i         :  3 Linear interpolation of sigma from neighboring points
C     i         :    sig(q) = sum_i wt_i A+_i sig(qpi) A_i
C     i         :    where A_i = z^LDA_qi (z_q^LDA)^-1
C     i         :1000s digit nonzero flags quantities other than h and s are
C     i         :      to be returned in arrays h and s
C     i         :    1 sigma not calculated.  Return with:
C     i         :      evecs of the LDA hamiltonian in array s;
C     i         :      evals of the LDA hamiltonian in array h, first column
C     i         :    2 return raw sigma, orbital basis, in s and
C     i         :      1st col of array h, LDA evals (real vector)
C     i         :    3 return in array s, sigma in LDA basis and
C     i         :      1st col of array h, LDA evals (real vector)
Ci         :    4 Same as 1000s digit mode 3, but
Ci         :      low- and high- energy blocks replaced by diagonal parts
Ci         :    5 sigm returned in orbital basis, but with
Ci         :      low- and high- energy blocks replaced by diagonal parts
Ci         :    6 Returns evecs of full hamiltonian in array s,
Ci         :      evals of the full hamiltonian in array h (NOT IMPLEMENTED)
Ci         :(Digit has a special meaning unique to interpolation mode 3:)
Ci         :    number of interpolation points
Ci         :    0 -> use default (4 points)
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: augmbl smhsbl hsibl
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: augmbl smhsbl hsibl
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc
Ci     Stored:    *
Ci     Passed to: augmbl smhsbl hsibl
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oindxo
Ci     Stored:    *
Ci     Passed to: *
Ci   isp   :spin index
Ci   ispc  :index where to store h,s. (for noncollinear mode)
Ci   q     :Bloch vector (k-point)
Ci   k1,k2,k3 dimensions of smpot
Ci   qsig  :q-points at which sigma can be computed
Ci   nqsig :number of qsig
Ci   smpot :smooth potential on uniform mesh (mkpot.f)
Ci   vconst:additional constant potential
Ci   lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   osig,otau,oppi  augmentation matrices
Ci   alfa  :add alfa * overlap to hamiltonian
Ci         :This is for stability in evals.  Preferably alfa=0
Ci   ndimh :dimension of hamiltonian and overlap h,s
Co Outputs
Co   mode  :if mode returns <0, hambls finds (approximate)
Co         :evals and evecs internally; see also description of h,s
Co   h     :Hamiltonian matrix
Co         :NB: if mode returned as -1, h has already been diagonalized
Co         :In that case, the diagonal of h holds eigenvalues and
Co         :s the eigenvectors.  The method for doing this, and the
Co         :motivation for it, are described in Remarks.
Co   s     :overlap matrix
Co         :NB: if mode is returned as -1, see description for h.
Co  ismidb :upper limit to middle block of sigma
Cl Local variables
Cl   lwsig :is 1000s digit of mode.  See Remarks for further description
Cl         :1 hambls does not generate hamiltonian, but returns
Cl         :  LDA eigenvectors in place of s
Cl         :  LDA eigenvalues in the first column of h
Cl         :  In this case, sigm is not generated or used
Cl         :2 hambls does not generate hamiltonian, but returns
Cl         :  sigm in orbital basis in place of h
Cl         :3 hambls does not generate hamiltonian, but returns
Cl         :  sigm rotated to lda basis in place of h
Cl         :4 Same as lwsig=3, but
Cl         :  low- and high- energy blocks replaced by diagonal parts
Cl         :5 sigm returned in orbital basis, but with
Cl         :  low- and high- energy blocks replaced by diagonal parts
Cl   mxitrp:maximum no qp to include in interpolation for self-energy
!!
!! takao think the following sigp mechanism is not so meaningful
!! I replaces it just by a number ESEAVER---
Cl   sigp  :parameters for approximating self-energy sigma.  sigma
Cl         :is approximated by its diagonal part sigii for energies
Cl         :below a low-energy cutoff (specified nmin or emin) and
Cl         :above a low-energy cutoff (specified nmax or emax).
Cl         : arg 1: mode : specifies how to set its diagonal part
Cl         :        for states above the high-energy cutoff nmax or emax
Cl         :        0 constrain sigii to be > asig+bsig*e
Cl         :        1 constrain sigii to be = asig+bsig*e
Cl         :        2 constrain sigii to be > asig and < bsig
Cl         :        3 constraint same as mode 1.
Cl         :          Mode 3 differs in that the least-squares fit to
Cl         :          sigii (for informational purposes only, to help
Cl         :          estimate asig and bsig) is done for states between
Cl         :          efit and nmax or emax
Cl         : arg 2: nmin : sigma for states 1..nmin are approximated
Cl         :               by sigii
Cl         : arg 3: emin : (used only if nmin<0)
Cl         :             : sigma for levels e<emin are approximated
Cl         :               by sigii
Cl         : arg 4: nmax : sigma for levels i>nmax are approximated by
Cl         :               sigii AND constrained according to mode
Cl         : arg 5: emax : (used only if nmax<=0)
Cl         :             : sigma for levels e<emax are approximated by
Cl         :               sigii AND constrained according to mode
Cl         : arg 6: asig : constraint used to approximate
Cl         :                 sigii = asig + E * bsig  or
Cl         :                 asig < sigii < bsig
Cl         : arg 7: bsig : constraint used to approximate
Cl         :                 sigii = asig + E * bsig  or
Cl         :                 asig < sigii < bsig
Cl         : arg 8: efit : fit sigii between efit and emax (mode 3)
!! Followings are still useful (but not so clean) method.
Cr Remarks
Cr   Adding the self-energy to the LDA hamiltonian is challenging,
Cr   because of problems with interpolation of the high-energy
Cr   parts of sigma to arbitrary k.
Cr
Cr   A standard approach has been developed, which approximates the
Cr   high-energy parts of sigma with a diagonal matrix.  This is what
Cr   is done when 100s digit of mode is set.  The procedure is:
Cr
Cr   A. generate evals and evecs of the LDA hamiltonian
Cr
Cr   B. Make a Bloch sum sigm(k) of sigm(RS)
Cr
Cr   C. Rotate (sigm-vxc)(k) to LDA basis:
Cr      sigm(Lda) = (Zlda)+ sigma-vxc (Zlda)
Cr
Cr   D. Overwrite sigm(k) with a diagonal array for high-lying
Cr      (and optionally low-lying) energies.  This step is the crucial
Cr      one that enables interpolation to proceed smoothly
Cr
Cr   E. Add Hlda = Elda to sigm generate total hamiltonian in LDA basis
Cr
Cr   F. Diagonalize Hlda+sigm, LDA basis
Cr
Cr   G. Rotate evecs from LDA basis to orbital basis
Cr
Cr   Other quantities may be sought, depending on lswig.
Cr   lwsig:
Cr     1  Return LDA evecs in s, evals in h (stop at step A above).
Cr        sigm is not used in this case.
Cr
Cr     2  Return with unmodified sigma (orbital basis) in s,
Cr        (stop at step B above) and LDA evals in h
Cr
Cr     3  Return with unmodified sigma (lda basis) in s
Cr        (stop at step C above) and LDA evals in h
Cr
Cr     4  Return with modified sigma (lda basis) in s
Cr        (stop at step D above) and LDA evals in h
Cu Updates
Cu   04 Jul 08 (T. Kotani) New PW addition to basis
Cu   27 Jun 08 Redesigned lwsig modes
Cu   26 May 07 Some preparation for rotation betw/ LDA, GW basis
Cu   16 Jan 07 Redefined functions of lwsig
Cu   20 Feb 05 Add switch (1000s digit mode) to return sigma in h
Cu    7 Feb 07 SO coupling made to with with self-energy
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   25 Jun 03 Update to incorporate sigp mode 3
Cu   24 May 03 New interpolation modes
Cu   10 Jan 03 First created
C ----------------------------------------------------------------------
      implicit none
      integer mode,nbas,isp,ispc,k1,k2,k3,lcplxp,ndimh,nqsig,lso,ismidb
      integer napw,igvapw(3,napw)
      type(s_rv1) :: sv_p_oppi(3,nbas)
      type(s_rv1) :: sv_p_otau(3,nbas)
      type(s_rv1) :: sv_p_osig(3,nbas)
      double precision qsig(3,nqsig),alfa
      real(8):: q(3) , vconst
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_ham)::sham
      type(s_lat)::slat
      double complex smpot(k1,k2,k3)
      integer:: isize_h(3)
      double complex h(ndimh,ndimh,isize_h(3)),s(ndimh,ndimh,isize_h(3))
      double precision qpi(3,4),tol,ddot,plat(3,3),qlat(3,3),norm(4,4),
     .     wt(4),swt,dq,emin,emax,fasig,fbsig,asig,bsig,xx,sigp(10)
      equivalence (emin,sigp(3)),(emax,sigp(5))
      equivalence (asig,sigp(6)),(bsig,sigp(7))
      integer i,j,lrsig,ipr,iprm(nqsig),ipiv(ndimh),dinv,nev,ndim1
      integer hreal,nl,nglob,mxorb,nttabs,ival,nsp,ndhrs,nsim,moditp,
     .     mxitrp,n123(4),off1,nspc,ndimhx,lwsig
      logical T,F
      parameter (T=.true.,F=.false.,tol=1d-7)
      complex(8),allocatable:: hqi(:,:,:),sigq(:,:),sqi(:,:,:)
      complex(8),allocatable:: zqi(:,:),z0(:,:)
      real(8),allocatable:: e(:),e0(:)
      integer,allocatable:: ipm(:)
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      double precision qloc(3,nqsig)
      double precision wk(ndimh*11),evl(ndimh)
      integer ::iwdummy
!     !takao
      complex(8),allocatable:: ss(:,:)
      integer:: sigswitch,debugmode
      logical:: pmtsig
      logical:: debug=.false.
      logical:: l_dummy_isanrg,isanrg,oncewrite
      integer:: ndimsig         !dimension of sigma
      complex(8),allocatable:: sigmat(:,:)
      complex(8),allocatable:: ovlmtoi(:,:),ovliovl(:,:)
      integer:: lb1,lb2,ilb
      call tcn('hambls')

!! === LDA hamiltonian ===
      if (mod(mode,10) .eq. 0) then
         j = 0
         if (lso .eq. 1) j = 10
         call hambl ( j , nbas , ssite , sspec , slat , sham , isp , q 
     .        , k1 , k2 , k3 , smpot , vconst , sv_p_osig , sv_p_otau , sv_p_oppi 
     .        , lcplxp , alfa , ndimh , napw , igvapw , h ( 1 , 1 , ispc ) 
     .        , s ( 1 , 1 , ispc ) , h ( 1 , 1 , 3 ) )
      endif
      lrsig = mod(mode/10,10)
      lwsig = mod(mode/1000,10)
      ismidb = 0
      l_dummy_isanrg=isanrg(lwsig,0,5,'hambls:','lwsig',T)
      if(debug) print *, ' mode napw lwsig',mode,napw,lwsig

c$$$c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccc             
c mar2014
c$$$      print *,'test ddddddddddddddddddddd'
c$$$      lb1=5                      
c$$$      lb2=9                      
c$$$      h(lb1:lb2,:,ispc) = 0d0    
c$$$      h(:,lb1:lb2,ispc) = 0d0    
c$$$      s(lb1:lb2,:,ispc) = 0d0    
c$$$      s(:,lb1:lb2,ispc) = 0d0    
c$$$      do ilb=lb1,lb2             
c$$$        h(ilb,ilb,ispc) = 99d0
c$$$        s(ilb,ilb,ispc) = 1d0    
c$$$      enddo                      
c$$$      lb1= 21                     
c$$$      lb2= 26
c$$$      h(lb1:lb2,:,ispc) = 0d0    
c$$$      h(:,lb1:lb2,ispc) = 0d0    
c$$$      s(lb1:lb2,:,ispc) = 0d0    
c$$$      s(:,lb1:lb2,ispc) = 0d0    
c$$$      do ilb=lb1,lb2             
c$$$        h(ilb,ilb,ispc) = 99d0
c$$$        s(ilb,ilb,ispc) = 1d0    
c$$$      enddo                      
c$$$      lb1= 26                     
c$$$      lb2= ndimh
c$$$      h(lb1:lb2,:,ispc) = 0d0    
c$$$      h(:,lb1:lb2,ispc) = 0d0    
c$$$      s(lb1:lb2,:,ispc) = 0d0    
c$$$      s(:,lb1:lb2,ispc) = 0d0    
c$$$      do ilb=lb1,lb2             
c$$$        h(ilb,ilb,ispc) = 99d0
c$$$        s(ilb,ilb,ispc) = 1d0    
c$$$      enddo                      
c$$$c mar2014

c$$$c----- Cu Fe 3d band extract test ---------------
c$$$      print *,'test ddddddddddddddddddddd'
c$$$      lb1=1                      
c$$$      lb2=4                      
c$$$      h(lb1:lb2,:,ispc) = 0d0    
c$$$      h(:,lb1:lb2,ispc) = 0d0    
c$$$      s(lb1:lb2,:,ispc) = 0d0    
c$$$      s(:,lb1:lb2,ispc) = 0d0    
c$$$      do ilb=lb1,lb2             
c$$$        h(ilb,ilb,ispc) = 99d0
c$$$        s(ilb,ilb,ispc) = 1d0    
c$$$      enddo                      
c$$$      lb1= 10                     
c$$$      lb2= 20
c$$$      h(lb1:lb2,:,ispc) = 0d0    
c$$$      h(:,lb1:lb2,ispc) = 0d0    
c$$$      s(lb1:lb2,:,ispc) = 0d0    
c$$$      s(:,lb1:lb2,ispc) = 0d0    
c$$$      do ilb=lb1,lb2             
c$$$        h(ilb,ilb,ispc) = 99d0
c$$$        s(ilb,ilb,ispc) = 1d0    
c$$$      enddo                      
c$$$      lb1= 26                     
c$$$      lb2= ndimh
c$$$      h(lb1:lb2,:,ispc) = 0d0    
c$$$      h(:,lb1:lb2,ispc) = 0d0    
c$$$      s(lb1:lb2,:,ispc) = 0d0    
c$$$      s(:,lb1:lb2,ispc) = 0d0    
c$$$      do ilb=lb1,lb2             
c$$$        h(ilb,ilb,ispc) = 99d0
c$$$        s(ilb,ilb,ispc) = 1d0    
c$$$      enddo                      
c$$$c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccc                

!! ... Store LDA eigenvectors, eigenvalues in s and h and exit
      if (lwsig .eq. 1) then
         allocate(zqi(ndimh,ndimh))
         allocate(z0(ndimh,ndimh))
         allocate(e(ndimh))
         call zhevx(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),1,T,ndimh,9d9,i,
     .        z0,F,e,ndimh,zqi)
         call dcopy(ndimh**2*2,zqi,1,s(1,1,ispc),1)
         call dcopy(ndimh,e,1,h,1)
         deallocate(z0,zqi,e)
         goto 999
      endif
      if (lrsig .eq. 0) goto 999
!! end of LDA part (lrsig=0)



!! === Start sigm=\Sigma-Vxc interpolation part ===
!! complecated. I like to skip to other way in future (see bndfp.F).
      nspc = globalvariables%nspc
      ndimhx = ndimh * nspc

C     ... Setup for Bloch transform
      if(debug) print *,' goto set up bloch transformation'
      ndhrs=sham%ndhrs
C     hreal=1 -> sig(T) is real
      hreal = 0
      if (lrsig .ge. 4) hreal = 1
      nl = globalvariables%nl
      mxorb = globalvariables%mxorb
      nttabs = ival ( sham%iv_a_ontabs , nbas + 1 )
      nsp  = globalvariables%nsp
      l_dummy_isanrg=isanrg(isp,1,nsp,'hambl:','isp',T)
      plat=slat%plat
      call dinv33(plat,1,qlat,swt)
      if(debug) print *,' hambles:eee 111222 sham%iv_a_ontabs=',sham%iv_a_ontabs
      sigp=sham%sigp
      if(debug) print *,' hambles: eeeee1111'



!! ===takao simplified case for the PMT method. ===
!! Not yet (pmtsig()=F). This does not work well.
      if (pmtsig()) then 
         print *,'hambles:  xxxxxxxxxx new pmtsig=',pmtsig()
C     In this branch, sigq is just a work array
         allocate(ss(ndimsig,ndimsig))
C     ... Bloch transform sigm(RS)-sigm(k).  This is step B in Remarks
C     Call to bloch : 104110=perm orb, transpose, no add, c*16
         i = 100000 + 4000 + 40*(1-hreal) + 000 + 10
c     kino Dec.9.2011:              call bloch ( i , q , nl , plat , mxorb , iv_p_oiprmb , 1 , nttabs
c     kino Dec.27.2011:              call bloch ( i , q , nl , plat , mxorb , sham%iv_p_oindxo , 1 , nttabs
         call bloch ( i , q , nl , plat , mxorb , sham%iv_a_oindxo , 1 , nttabs
c     kino Dec.9.2011:       .      , iv_p_oiaxs , rv_p_ohrs , ndhrs , isp , nsp , ndimh , ndimh
c     kino Dec.28.2011:       .      , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs , isp , nsp , ndimh , ndimh 
     .        , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs , isp , nsp , ndimsig , ndimsig 
     .        , 0 , ndimsig, 0 , ndimsig , 0 , ss ( 1 , 1 ) , iwdummy , iwdummy )
         allocate(ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
         ovlmtoi = s(1:ndimsig,1:ndimsig,ispc)

C cccccccccccccccccccccccccccccccccccccccccccc
C c mar2014 
C       lb1=137
C       lb2=168
C       ovlmtoi(lb1:lb2,:) = 0d0
C       ovlmtoi(:,lb1:lb2) = 0d0
C       do ilb=lb1,lb2
C          ovlmtoi(ilb,ilb)=1d0
C       enddo
C cccccccccccccccccccccccccccccccccccccccccccc

         call matcinv(ndimsig,ovlmtoi)
         ovliovl = matmul(ovlmtoi,s(1:ndimsig,1:ndimh,ispc))
         deallocate(ovlmtoi)
!! for sclaedsigma dec2012
         if(sham%scaledsigma/=1d0) then
            print *,' QSGW-LDA hybrid: ScaledSigma=alpha_Q=',sham%scaledsigma
            ss = sham%scaledsigma* ss
         endif   
         h(:,:,ispc) = h(:,:,ispc) + matmul(transpose(dconjg(ovliovl)), matmul(ss(:,:),ovliovl))
         deallocate(ss,ovliovl)
         call tcx('hambls')
         return
      endif


!! == Followings are original Mark's metod ==
      moditp = mod(mode/100,10)
      mxitrp = mod(mode/1000,10)
      if (mxitrp .eq. 0) mxitrp = 4
!! --- Interpolate sigma keeping higher energy LDA eigenvectors ---
      if(debugmode()>0) print *,' moditp mxitrp=',moditp,mxitrp
      if (moditp .eq. 1) then
C     In this branch, sigq is just a work array
         allocate(sigq(ndimh,ndimh))
         allocate(zqi(ndimh,ndimh))
         allocate(z0(ndimh,ndimh))
         allocate(e(ndimh))
         allocate(e0(ndimh))
!! === conversion matrix from sigma(1:ndimsig,1:ndimsig) to sigma(1:ndimh,1:ndimh). ===
         allocate(ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
         ovlmtoi = s(1:ndimsig,1:ndimsig,ispc)
         call matcinv(ndimsig,ovlmtoi)
         ovliovl = matmul(ovlmtoi,s(1:ndimsig,1:ndimh,ispc))
         deallocate(ovlmtoi)
!!  === zqi <- evecs, evals of LDA hamiltonian = Z.lda ===
C     This is step A in Remarks
*     call zprm('hamiltonian, LDA',2,h,ndimh,ndimh,ndimh)
*     call zprm('overlap',2,s,ndimh,ndimh,ndimh)
         call zhevx(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),1,T,ndimh,9d9,i,
     .        z0,F,e,ndimh,zqi)
C     ... Bloch transform from sham%rv_a_ohrs=sigm(RealSpace) to sigm(k).  This is step B in Remarks
C     Call to bloch : 104110=perm orb, transpose, no add, c*16
         i = 100000 + 4000 + 40*(1-hreal) + 000 + 10
c     call bloch ( i , q , nl , plat , mxorb , sham%iv_a_oindxo , 1 , nttabs
c     .      , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs , isp , nsp , ndimh , ndimh 
c     .      , 0 , ndimh , 0 , ndimh , 0 , s ( 1 , 1 , ispc ) , iwdummy , iwdummy )
c     print *,'goto bloch sigmat ndimsig=',ndimsig
         allocate(sigmat(ndimsig,ndimsig))
         call bloch ( i , q , nl , plat , mxorb , sham%iv_a_oindxo , 1 , nttabs
     .        , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs , isp , nsp , ndimsig , ndimsig
     .        , 0 , ndimsig , 0 , ndimsig , 0 , sigmat( 1 , 1  ) , iwdummy , iwdummy )
!! sigma (s) is moved back to s in ndimh.
         s(:,:,ispc)= matmul(transpose(dconjg(ovliovl)), matmul(sigmat(:,:),ovliovl))
         deallocate(sigmat,ovliovl)
!! for sclaedsigma dec2012takao
         if(sham%scaledsigma/=1d0) then
            if(oncewrite(9)) print *,' QSGW-LDA hybrid: ScaledSigma=alpha_Q=',sham%scaledsigma
            s(:,:,ispc) = sham%scaledsigma* s(:,:,ispc)
         endif   
C     LDA evals in h
         if (lwsig .ge. 2 .and. lwsig .le. 4) then
            call dcopy(ndimh,e,1,h,1)
         endif
C     Return with raw sigm (orbital basis) in s, LDA evals in h
         if (lwsig .eq. 2) goto 99
C     ... Rotate sigm to LDA basis; step C in Remarks.
C     s <- (Z.lda)+ sigma-vxc Z.lda = sigma-vxc in LDA basis
C     call zprm('sigma, orbital basis',2,s,ndimh,ndimh,ndimh)
#ifdef USEOPTIONAL
         call phmbls(mode=3,ndimh=ndimh,eval=e0,wk=z0,h=s(1,1,ispc),zd=zqi,z=zqi,zhz=s(1,1,ispc))
#else
         call phmbls(3,ndimh,e0,xx,z0,s(1,1,ispc),zqi,zqi,s(1,1,ispc))
#endif
C     Return with raw sigm (LDA basis) in s, LDA evals in h
         if (lwsig .eq. 3) goto 99


c     takao
c     C   ... Add to least-squares fit of sigma for high energies
c     call phmbl3(2,ndimh,0,e,e0,sigp,fasig,fbsig)
C     ... Overwrite z+ sigma z with approx for low,high-lying states

C     This is step D in Remarks
!! new psigd (defined in rdsigm.F)
         if(sigswitch() ==0) then
            call psigd(0,ndimh,s(1,1,ispc),e,sigp,n123,wk(1:ndimh),q=q,isp=isp,iout=0) 
!! previous version. kept for backword compatibility.
         else
            if(oncewrite(1)) print *,
     &           ' sigswitch()=1: a*e+b interpolation formula is not used.'//
     &           ' ---> Use asig=(average of ESEAVR), bsig=0'
            call psigd(0,ndimh,s(1,1,ispc),e,sigp,n123,wk, q=q,isp=isp,iout=1,iprx=1) 
!! this version read ESEAVR (constant. Not use a*e+b interpolation formula).
         endif
         ismidb = n123(3)
C     Return with modified sigm (LDA basis) in s, LDA evals in h
         if (lwsig .eq. 4) goto 99

C     ... Return with modified sigm (orb basis) in s, LDA evals in h
         if (lwsig .eq. 5) then
#ifdef USEOPTIONAL
            call phmbls(mode=33,ndimh=ndimh,eval=e0,wk=z0,h=s(1,1,ispc),zd=zqi,z=zqi,zhz=s(1,1,ispc))
#else
            call phmbls(33,ndimh,e0,xx,z0,s(1,1,ispc),zqi,zqi,s(1,1,ispc))
#endif
            call dcopy(ndimh,e,1,h,1)
            goto 99
         endif
C     ... s <- Z.lda+ (h(lda) + sigma) Z.lda = e(lda) + Z.lda+ sigma Z.lda
C     This is h+sigma(approx) in LDA basis, step E in Remarks
C     Initialize e0 with diagonal term of h+sigma_approx in LDA basis
C     For low and high-lying states, e0 is already the final answer
C     The middle eigenvalues e0 will be updated later.
         do  i = 1, ndimh
            s(i,i,ispc) = s(i,i,ispc) + e(i)
            e0(i)  = s(i,i,ispc)
         enddo
C     ... s <- evecs of h + sig in LDA basis, for states in middle block
C     This is step F described in Remarks
C     Middle block has range (off1+1:off1+ndim1,off1+1:off1+ndim1)
         off1  = n123(2)
         ndim1 = n123(3)-n123(2)
         call zhevx(ndim1,ndimh,s(1+off1,1+off1,ispc),xx,0,T,ndim1,9d9,i,
     .        sigq,F,e0(1+off1),ndimh,z0)
C     ... Rotate LDA evecs(orbital basis) to evecs of h+sig, orbital basis
C     This is step G described in Remarks
C     Given zqi = LDA evecs (orbital basis)
C     and   z0  = h+sig evecs, LDA basis
C     Then elda = zqi+  Hlda zqi  =>  e = z0+ zqi+ (Hlda+sig) zqi z0
C     Initialize with entire Z.lda; overwrite evecs off1+1..off1+ndim1
         call zcopy(ndimh**2,zqi,1,s(1,1,ispc),1)
         call zgemm('N','N',ndimh,ndim1,ndim1,one,zqi(1,off1+1),
     .        ndimh,z0,ndimh,zer,s(1,off1+1,ispc),ndimh)

c     takao jun2009
c     takao no re-ordering. So, some eigvenvalue in LDA(+constant) can be less than eigenvalues in qsgw block.
c     C   ... Reorder evecs,evals in ascending eval order
c     allocate(ipm(ndimh))
c     call phmbls(8+4,ndimh,e0,ipm,sigq,xx,s(1,1,ispc),s(1,1,ispc),xx)
c     deallocate(ipm)

C     ... Overwrite h with eigenvalues
         call dpzero(h(1,1,ispc),ndimh**2*2)
         do  i = 1, ndimh
            h(i,i,ispc) = e0(i)
         enddo

         if(debugmode()>0) then
            if(ispc==1) then
               do  i = 1, ndimh
                  write(6,"('ddd hhh:',2i4,12d13.5)") ispc,i,e0(i),sum(abs(s(:,i,1)))
               enddo
            elseif(ispc==2) then
               do  i = 1, ndimh
                  write(6,"('ddd hhh:',2i4,12d13.5)") ispc,i,e0(i),sum(abs(s(:,i,1))),sum(abs(s(:,i,2)))
               enddo
            endif  
         endif


!! === spin-orbit is added ===
!!   ... Add off-diagonal spin-orbit hamiltonian to h, and re-diagonalize
         if (lso .eq. 1 .and. ispc .eq. 2) then
c$$$  !! TEST of skip spin-off diagonal
c$$$  !! If you set this true, SO=1 becomes SO=2 because no LxSx+LySy is skiped.
c$$$  print *, '!! zero out SO ndimh ndimhx',ndimh,ndimhx
c$$$  h(1:ndimh,1:ndimh,3) = 0 
#ifdef USEOPTIONAL
            call phmbls(mode=1,ndimh=ndimh,eval=e0,wk=z0,
     .           h=h(1,1,3),zd=s(1,1,1),z=s(1,1,2),zhz=h(1,1,3))
#else
            call phmbls(1,ndimh,e0,xx,z0,h(1,1,3),s(1,1,1),s(1,1,2),h(1,1,3))
#endif
            deallocate(z0,zqi,e0)
            allocate(z0(ndimhx,ndimhx),zqi(ndimhx,ndimhx),e0(ndimhx))
!!    zqi <- evecs of H w/out SO in orbital basis, spin-coupled form 
!!    h <- full hamiltonian in zqi basis, spin-coupled form
            call sopert(1,ndimh,nspc,z0,s,zqi) !z0 work array.
            call sopert(0,ndimh,nspc,z0,h,h) !z0 work array
!!    Diagonalize h in the representation of zqi.
!!     h_ij= <zqi(:,i)|h|zqi(:,j)>, call zhevx with lov=0 means no overlap matrix.
            call zhevx(ndimhx,ndimhx,h,xx,0,T,ndimhx,9d9, 
     .           i,s,F,e0,ndimhx,z0) !z0 contains evec in the zqi basis.
C     Rotate evecs to orbital basis = z(so) z(orbital)
C     Given zqi = evecs w/out L.S (orbital basis)
C     and   z0  = h+sig evecs, LDA basis
C     Then e = z0+ zqi+ (H+SO) zqi z0
            call zgemm('N','N',ndimhx,ndimhx,ndimhx,one,zqi,
     .           ndimhx,z0,ndimhx,zer,s,ndimhx)
C     Replace h with evals on the diagonal
            call dpzero(h,ndimhx**2*2)
            call dcopy(ndimhx,e0,1,h,2*ndimhx+2)
*     call zprm('evals of h+sig+L.S.',2,h,ndimhx,ndimhx,ndimhx)
         endif
 99      mode = -1
         deallocate(sigq,zqi,z0,e,e0)
      else
         call rxi('unknown interpolation mode %i',moditp)
      endif
      call getpr(ipr)
  999 continue
      call tcx('hambls')
      return
      end subroutine hambls



      subroutine blsig(mode,nbas,sham,isp,nsp,nspc,plat,q,
     .lwtkb,zval,iq,wtkp,wtkb,ndimh,zq,sig)
      use m_struc_def  !Cgetarg
      use m_globalvariables
      use m_phmbls
C- sigma-vxc or sum_occ <sigma-vxc> for one k-point.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit specifies what to do:
Ci         :  0 Make sigma and store in sig
Ci         :  1 Make sum_states wtkb(iq) (zq+)*sig*zq
Ci         :10s digit specifies self-energy matrix elements
Ci         :  0 No self-energy to compute
Ci         :  1 r.s. sigma is real
Ci         :  2 r.s. sigma is complex
Ci   nbas  :size of basis
Ci     Elts read: ndhrs ontabs oiaxs ohrs oindxo
Ci     Stored:
Ci     Passed to:
Ci   isp   :generate self-energy for this spin
Ci         :Not used in the noncollinear case
Ci   plat  :primitive lattice vectors
Ci   q     :Bloch vector (k-point)
Ci   wtkp  :q-point weights from symmetry operations
Ci   lwtkb :(used for 1s digit mode=1)
Ci         :0 assign weight 1 for first zval states, 0 otherwise
Ci         :1 weights are required a priori, and are read from disk
Ci   wtkb  :integration weights, needed if lwtkb is 1
Ci   ndimh :dimension of hamiltonian and self-energy sig
Ci   zq    :evecs (used for 1s digit mode=1)
Cio Inputs/Outputs
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: ndhrs ontabs oiaxs ohrs oindxo eterms
Ci     Stored:    eterms  (1s digit mode=1)
Co Outputs
Co   sig   :(1s digit mode = 0) self-energy matrix
Co         :(1s digit mode = 1) <zq sig zq>
Co         :Note: in noncollinear case, passed sig MUST be
Co         :dimensioned sig(ndimh,nspc,ndimh,nspc) even though
Co         :dimension appears as sig(ndimh,ndimh) in this routine.
Cl Local variables
Cr Remarks
Cu Updates
Cu   11 Jan 05 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,iq,isp,ndimh,nsp,nspc,lwtkb
      real(8):: plat(3,3) , q(3) , wtkb(ndimh,nsp,1) , wtkp(iq) , zval
      type(s_ham)::sham

      double complex zq(ndimh,nspc,ndimh,nspc),sig(ndimh,ndimh,2)
C ... Local parameters
ckino Dec.9.2011:            integer,pointer :: iv_p_ontabs(:) =>NULL()

ckino Dec.9.2011:           integer,pointer :: iv_p_oiaxs(:) =>NULL()

ckino Dec.9.2011:          real(8),pointer :: rv_p_ohrs(:) =>NULL()

ckino Dec.9.2011:         integer,pointer :: iv_p_oiprmb(:) =>NULL()

      integer i,lrsig,ipr,jsp,ksp
      integer hreal,nl,nglob,mxorb,nttabs,ival,ndhrs,nevec,ndimhx
      logical T,F
      parameter (T=.true.,F=.false.)
      double precision eterms(20),rhosig,swgt
      equivalence (eterms(19),rhosig)
CSFCPP#if F90
      complex(8),allocatable:: wk(:,:,:,:)
CSFCPP#else
CSFCPP      double complex wk(ndimh,nspc,ndimh,nspc)
CSFCPP#endif
      double precision zsz(ndimh*nspc),ewgt(ndimh*nspc)

C ... Heap
      integer ::iwdummy

      logical:: l_dummy_isanrg,isanrg

      save swgt

      lrsig = mod(mode/10,10)
      if (lrsig .eq. 0) return
      call tcn('blsig')

C ... Setup for Bloch transform
      ndhrs=sham%ndhrs
C     hreal=1 -> sig(T) is real
      hreal = 0
      if (lrsig .ge. 4) hreal = 1
Changenglob      nl = nglob('nl')
      nl = globalvariables%nl
Changenglob      mxorb = nglob('mxorb')
      mxorb = globalvariables%mxorb
ckino Dec.9.2011:         nttabs = ival ( iv_p_ontabs , nbas + 1 )
ckino Dec.14.2011:         nttabs = ival ( sham%iv_p_ontabs , nbas + 1 )
      nttabs = ival ( sham%iv_a_ontabs , nbas + 1 )

Ckino isanrg is logical function,       call isanrg(isp,1,nsp,'blsig:','isp',T)
      l_dummy_isanrg=isanrg(isp,1,nsp,'blsig:','isp',T)
      ndimhx = ndimh*nspc

CSFCPP#if F90
      allocate(wk(ndimh,nspc,ndimh,nspc))
CSFCPP#endif

C --- Bloch sum sigma-vxc, orbital basis ---
C     Call to bloch : 104110=perm orb, transpose, no add, c*16
      i = 100000 + 4000 + 40*(1-hreal) + 000 + 10
      do  jsp = 1, nspc
        ksp = isp
        if (nspc .eq. 2) ksp = jsp
C     collinear case:     ksp = isp, jsp=1
C     noncollinear case:  ksp=1..2,  jsp=1..2
ckino Dec.9.2011:              call bloch ( i , q , nl , plat , mxorb , iv_p_oiprmb , 1 , nttabs
ckino Dec.27.2011:              call bloch ( i , q , nl , plat , mxorb , sham%iv_p_oindxo , 1 , nttabs
        call bloch ( i , q , nl , plat , mxorb , sham%iv_a_oindxo , 1 , nttabs
ckino Dec.9.2011:       .      , iv_p_oiaxs , rv_p_ohrs , ndhrs , ksp , nsp , ndimh , ndimh
ckino Dec.9.2011:       .      , sham%iv_p_oiaxs , rv_p_ohrs , ndhrs , ksp , nsp , ndimh , ndimh
ckino Dec.14.2011:       .      , sham%iv_p_oiaxs , sham%rv_p_ohrs , ndhrs , ksp , nsp , ndimh , ndimh
ckino Dec.28.2011:       .      , sham%iv_a_oiaxs , sham%rv_p_ohrs , ndhrs , ksp , nsp , ndimh , ndimh 
     .      , sham%iv_a_oiaxs , sham%rv_a_ohrs , ndhrs , ksp , nsp , ndimh , ndimh 
     .      , 0 , ndimh , 0 , ndimh , 0 , sig ( 1 , 1 , jsp ) , iwdummy , 
     .      iwdummy )

      enddo
C     Convert sig to full matrix with zero off-diagonal spin blocks
      if (nspc .eq. 2) then
        call sopert(1,ndimh,nspc,wk,sig,sig)
      endif

      call getpr(ipr)
      if (ipr .ge. 100 .and. mode .ge. 0) then
        call info(0,0,0,'sig for q=%3:1;6,6d',q,0)
        call zprm('sig',2,sig,ndimhx,ndimhx,ndimhx)
      endif

C ... ksp = current spin index (isp or 1, in noncollinear case)
      ksp = isp
      if (nspc .eq. 2) ksp = 1

C --- Compute double-counting contribution to total energy ---
      if (mod(mode,10) .eq. 1) then
C  ...  zsz <- diagonal part of <sigma-vxc>
#ifdef USEOPTIONAL
        call phmbls ( mode=3 , ndimh=ndimhx , eval=zsz  , wk=wk , h=sig ,zd= zq , z=zq 
     .   , zhz=sig )
#else
        call phmbls ( 3 , ndimhx , zsz , iwdummy , wk , sig , zq , zq 
     .   , sig )
#endif

C  ...  Make weights for d.c. sum
        call dpzero(ewgt,ndimhx)
        if (lwtkb .eq. 0) then
          nevec = zval/2d0 + 0.500001d0
          if (nspc .eq. 2) nevec = 2*nevec
          do  i = 1, nevec
            ewgt(i) = abs(wtkp(iq))/nsp
          enddo
        else
          do  i = ndimhx, 1, -1
            nevec = i
            if (abs(wtkb(i,ksp,iq)) .gt. 1d-8) goto 12
          enddo
   12     continue
          do  i = 1, nevec
            ewgt(i) = abs(wtkb(i,ksp,iq))
          enddo
        endif
C  ...  Contribution to d.c. sum for this qp

c        i_copy_size=size(sham%eterms)
c        call dcopy(i_copy_size,sham%eterms,1,eterms,1)
        eterms=sham%eterms

        if (iq .eq. 1 .and. ksp .eq. 1) then
          rhosig = 0
          swgt = 0
        endif
        do  i = 1, nevec
          swgt = swgt + ewgt(i)
          rhosig = rhosig + ewgt(i)*zsz(i)
        enddo

c        i_copy_size=size(sham%eterms)
c        call dcopy(i_copy_size,eterms,1,sham%eterms,1)
        sham%eterms=eterms
        

C  ...  Sanity check
C       print *, 'iq=',iq,'swgt=',swgt,'rhosig=',rhosig
      endif

CSFCPP#if F90
      deallocate(wk)
CSFCPP#endif

      call tcx('blsig')
      end subroutine blsig

      subroutine permhk(mode,nbas,ldh,offH,h,hp,lblock)
      use m_phmbls
C- Permutation of H by kappa blocks
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy permuted h into hp.  Both rows, cols permuted.
Ci         :1 copying in reverse permutation order hp to h.
Ci            mode=1 is inverse operation of mode 0.
Ci            Thus this sequence of calls:
Ci              call permhk(0,nbas,ldh,offH,h,hp,lblock)
Ci              h = 0
Ci              call permhk(1,nbas,ldh,offH,h,hp,lblock)
Ci            copies h to hp, permuting hp; then initializes h and
Ci            restores the original h from hp.
Ci         :2 copy permuted h into hp.  Only cols are permuted.
Ci         :3 copying in reverse permutation order hp to h.
Ci            mode=3 is inverse operation of mode 2.
Ci            Thus this sequence of calls:
Ci              call permhk(2,nbas,ldh,offH,h,hp,lblock)
Ci              h = 0
Ci              call permhk(3,nbas,ldh,offH,h,hp,lblock)
Ci            copies h to hp, permuting hp; then initializes h and
Ci            restores the original h from hp.
Ci         :4 copy permuted h into hp.  Only rows are permuted.
Ci         :5 copying in reverse permutation order hp to h.
Ci            mode=5 is inverse operation of mode 4
Ci            Thus this sequence of calls:
Ci              call permhk(4,nbas,ldh,offH,h,hp,lblock)
Ci              h = 0
Ci              call permhk(5,nbas,ldh,offH,h,hp,lblock)
Ci            copies h to hp, permuting hp; then initializes h and
Ci            restores the original h from hp.
Ci   nbas  :size of basis
Ci   ldh :dimension of h
Ci   offH  :table of offsets to unpermuted hamiltonian (makidx.f)
Ci         :offH(m,ik,ib) = offset to subblock (ik,ib), where
Ci         :m=1 for lower orbitals, m=2 for intermediate, etc.
Ci         :offH(m,ik+1,ib)-offH(m,ik,ib) = size of block
Cio Inputs/Outputs
Cio  h     :mode=0: (input) matrix to be permuted
Cio        :mode=2: (input) matrix to be permuted
Cio        :mode=1: (output) matrix to be restored from hp.
Cio        :        Both rows and columns are (reverse) permuted
Cio        :mode=3: (output) matrix to be restored from hp
Cio        :        Only columns are (reverse) permuted.
Cio  hp    :mode=0: (output) matrix copied from h where both
Cio                 rows and columns are permuted.
Cio        :mode=2: (output) matrix copied from h where only
Cio                 columns are permuted.
Cio        :mode=1: (input) matrix used to restore h
Cio        :mode=3: (input) matrix used to restore h
Co Outputs
Co   lblock:dimension of lower block
Cu Updates
Cu   17 May 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,ldh,lblock
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer offh(n0H,nkap0,nbas+1)
      double complex h(ldh,ldh),hp(ldh,ldh)
C ... Local parameters
      integer nblk,m,i,ib,ikap

      integer iprm(nbas,nkap0),offs(nkap0,nbas+1)
      double precision xx

      nblk = nkap0*nbas
      m = 1
      call icopy(nblk+1,offH(m,1,1),n0H,offs,1)
      i = 0
      lblock = 0
      do  ib = 1, nbas
        do  ikap = 1, nkap0
          i = i+1
          iprm(ib,ikap) = i
        enddo
        lblock = lblock + offs(2,ib) - offs(1,ib)
      enddo

      if (mode .eq. 0) then
        call zpmblk(nblk,iprm,offs,xx,ldh,h,8120,xx,ldh,ldh,hp,i)
      elseif (mode .eq. 1) then
        call zpmblk(nblk,iprm,xx,offs,ldh,hp,8202,xx,ldh,ldh,h,i)
      elseif (mode .eq. 2) then
        call zpmblk(nblk,iprm,offs,xx,ldh,h,120,xx,ldh,ldh,hp,i)
      elseif (mode .eq. 3) then
        call zpmblk(nblk,iprm,xx,offs,ldh,hp,202,xx,ldh,ldh,h,i)
      elseif (mode .eq. 4) then
        call zpmblk(nblk,iprm,offs,xx,ldh,h,4120,xx,ldh,ldh,hp,i)
      elseif (mode .eq. 5) then
        call zpmblk(nblk,iprm,xx,offs,ldh,hp,4202,xx,ldh,ldh,h,i)
      else
        call rxi('permhk: bad mode',mode)
      endif

      end subroutine permhk



C      subroutine psigd(mode,ndimh,sig,eval,sigp,n123,sigd)
CC- Approximate sigma for low,higher energies with diagonal part,
CC- and further add constraints for the higher energies
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :0 sig is a complex square matrix; poke sig(i,i)
CCi         :  sigd is not used
CCi         :1 sigd is a real diagonal matrix; poke sigd(i)
CCi         :  for those elements that satisfy constraints.
CCi         :  sig is not used
CCi   ndimh :dimension of hamiltonian, and number of energies
CCi   eval  :list of (LDA) eigenvalues
CCi   sigp  :parameters for approximating self-energy sigma.  sigma
CCi         :is approximated by its diagonal part sigii for energies
CCi         :below a low-energy cutoff (specified nmin or emin) and
CCi         :above a low-energy cutoff (specified nmax or emax).
CCi         : arg 1: specifies how to set diagonal part sigii
CCi         :        for states above the high-energy cutoff nmax or emax
CCi         :        0 constrain sigii to be > asig+bsig*e
CCi         :        1 constrain sigii to be = asig+bsig*e
CCi         :        2 constrain sigii to be > asig and < bsig
CCi         :        3 constraint same as case 1; for this routine, there
CCi         :          is no difference.  Elsewhere,
CCi         :          arg1=3 differs in that the least-squares fit to
CCi         :          sigii (for informational purposes only, to help
CCi         :          estimate asig and bsig) is done for states between
CCi         :          efit and nmax or emax
CCi         : arg 2: nmin : usage depends on mode above.
CCi         :               mode = 0: for states 1..nmin, off-diagonal
CCi         :               parts of sig(1:nmin,1:nmin) are zeroed out.
CCi         :               mode = 1: sigd(1..nmin) is filled with emin
CCi         : arg 3: emin : usage depends on mode above.
CCi         :               mode = 0: for states e_i<emin, off-diagonal
CCi         :               parts of sig(1:i,1:i) are zeroed out.
CCi         :               mode = 1: sigd(1..nmin) is filled with emin
CCi         : arg 4: nmax : sigma for levels i>nmax are approximated by
CCi         :               sigii AND constrained according to arg 1
CCi         : arg 5: emax : (used only if nmax<=0)
CCi         :             : sigma for levels e<emax are approximated by
CCi         :               sigii AND constrained according to arg 1
CCi         : arg 6: asig : constraint used to approximate
CCi         :               sigii = asig + E * bsig  or
CCi         :               asig < sigii < bsig
CCi         : arg 7: bsig : constraint used to approximate
CCi         :               sigii = asig + E * bsig  or
CCi         :               asig < sigii < bsig
CCi         : arg 8: efit : fit sigii between efit and emax
CCio Inputs/Outputs
CCio  sig   :sigma, in LDA representation
CCio        :On output:
CCio        : *high and low states are replaced by diagonal part
CCio        :  of sigma
CCio        : *diagonal part may be altered to satisfy constraints
CCo Outputs
CCo   n123  :blocks sigma into lower, middle, high parts
CCo         :n123(1) = 0
CCo         :n123(2) = index to highest orbital in 'low' block
CCo         :n123(3) = index to highest orbital in 'middle' block
CCo         :n123(4) = ndimh
CCl Local variables
CCl   llow  :T if this eigenvalue is below minimum cutoff
CCl   lhigh :T if this eigenvalue is below above max cutoff
CCr Remarks
CCr
CCu Updates
CCu   19 May 03 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,ndimh,n123(4)
C      double precision eval(ndimh),sigp(10),sigd(ndimh)
C      double complex sig(ndimh,ndimh)
CC ... Local parameters
C      logical llow,lhigh
C      integer i,stdo,lgunit,ipr,nmin,nmax,modsgp,PRTE
C      double precision emin,emax,asig,bsig,tol,siglin
C      double complex zer,sigii
C      parameter (zer=(0d0,0d0),tol=1d-7,PRTE=45)
C
C      call getpr(ipr)
C      stdo = lgunit(1)
C      call ivset(n123,1,3,0)
C      n123(4) = ndimh
C
C      modsgp = nint(sigp(1))
C      if (modsgp .eq. 3) modsgp = 1
C      nmin   = nint(sigp(2))
C      emin   = sigp(3)
C      nmax   = nint(sigp(4))
C      emax   = sigp(5)
C      asig   = sigp(6)
C      bsig   = sigp(7)
C
C      if (mode .eq. 0) then
C        call info5(PRTE,1,0,
C     .  ' hambls: approximate sigma'//
C     .  '%?#(n<0)# for energies E(lda)<%d; and%-2j#%-1j#'//
C     .  '%?#(n>0)# for states %-1jn=%i and below; and##%j'//
C     .  '%?#(n<=0)# for energies E(lda)>%d%-2j#%-1j#'//
C     .  '%?#(n>0)# for states above %-1jn=%i##%j'//
C     .  '%N  state    E(lda)%8fsig_ii%4fconstraint%6fuse',
C     .  nmin,emin,nmax,emax,0)
C      elseif (mode .eq. 1) then
C        call info5(PRTE,1,0,
C     .  ' hambls: new diagonal sigma for:'//
C     .  ' %?#(n>0)#%-1j %i DEEP states (E=%d)#%jno DEEP states# and'//
C     .  ' %?#(n>0)#%-1j %i#no# HIGH states'//
C     .    ' ',nmin,emin,ndimh-nmax,0,0)
C        if (nmin .le. 0 .and. nmax .ge. ndimh) return
C        call info0(PRTE,0,0,'  state    E(lda)%8fsig_ii')
C        stop 'for now'
C      else
C        call rxi('psigd: bad mode ',mode)
C      endif
C
C      do  i = 1, ndimh
C
CC       Require evals to be ordered
C        if (eval(i) .lt. eval(max(i-1,1))-tol)
C     .    call rxi('psigd: eval %i not ordered',i)
C
CC       Decide whether this eval is in low, middle, or high block
C        llow  = (nmin .lt. 0 .and. eval(i) .lt. emin) .or.
C     .          (nmin .ge. 0 .and. i .le. nmin)
C        lhigh = (nmax .le. 0 .and. eval(i) .gt. emax) .or.
C     .          (nmax .gt. 0 .and. i .gt. nmax)
C        if (mode .eq. 1) llow = i .le. nmin
C        if (llow) n123(2) = i
C        if (.not. lhigh) n123(3) = i
C
CC       Calculate new diagonal sigma that satisfies constraints
C        sigii = 0
C        if (mode .eq. 0) then
C          sigii = sig(i,i)
C        elseif (mode .eq. 1 .and. llow) then
C          sigii = emin
C        endif
C        siglin = asig + bsig*eval(i)
C        if (lhigh) then
C          if (modsgp .eq. 0) then
C            if (dble(sigii) .lt. siglin) sigii = siglin
C          elseif (modsgp .eq. 1) then
C            sigii = siglin
C          elseif (modsgp .eq. 2) then
C            if (dble(sigii) .lt. asig) sigii = asig
C            if (dble(sigii) .gt. bsig) sigii = bsig
C          else
Ckino isanrg is logical function, C            call isanrg(modsgp,0,2,'hambls:','sig fit mode',.true.)
C            l_dummy_isanrg=isanrg(modsgp,0,2,'hambls:','sig fit mode',.true.)
C          endif
C        endif
C
C
CC       Printout
C        if (mode .eq. 1) then
C        if (ipr .ge. PRTE .and. (llow.or.lhigh)) then
C          write(stdo,331) i,eval(i),dble(sigii)
C  331     format(i6,f12.6,2x,f12.6)
C          endif
C        else
C          if (ipr .ge. PRTE .and. lhigh .and. modsgp .eq. 2) then
C            write(stdo,332)
C     .        i,eval(i),dble(sig(i,i)),asig,bsig,dble(sigii)
C  332       format(i6,f12.6,2x,f12.6,f7.2,',',f5.2,f12.6)
C          elseif (ipr .ge. PRTE .and. lhigh) then
C            write(stdo,333) i,eval(i),dble(sig(i,i)),siglin,dble(sigii)
C  333       format(i6,f12.6,2x,2f12.6,f13.6)
C          elseif (ipr .ge. PRTE) then
C            write(stdo,334) i,eval(i),dble(sig(i,i)),dble(sigii)
C  334       format(i6,f12.6,2x,f12.6,12x,f13.6)
C          endif
C        endif
C
CC       Overwrite full sigma with diagonal matrix, or write to sigd
C        if (llow .or. lhigh) then
C          if (mode .eq. 0) then
C            call zscal(ndimh,zer,sig(1,i),1)
C            call zscal(ndimh,zer,sig(i,1),ndimh)
C            sig(i,i) = sigii
C          else
C            sigd(i) = sigii
C          endif
C        endif
C
C      enddo
C
C      end

      subroutine phmbl3(mode,ndimh,nsmidb,e,de,sigp,a,b)

C- Fit coffs a,b in linear fit a+b*e(lda) = de[e(lda)]
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 initialize normal matrix
Ci         :2 add points to normal matrix
Ci         :4 solve normal matrix for coefficients a and b
Ci         :8 print out nsmidb
Ci         :NB: these switches may be taken in combination
Ci   ndimh :number of eigenvalues
Ci   nsmidb:smallest value of nmax encountered
Ci         :(only used in printout, mode=8)
Ci   e     :list of (LDA) eigenvalues
Ci   de    :list of self-energies
Ci   sigp  :parameters for approximating self-energy sigma.  sigma
Ci         :is approximated by its diagonal part sigii for energies
Ci         :below a low-energy cutoff (specified nmin or emin) and
Ci         :above a low-energy cutoff (specified nmax or emax).
Ci         :NB: this routine uses parts of sigp only to determine
Ci         :    which states are to be included in fit
Ci         : arg 1: mode : specifies how to set its diagonal part
Ci         :        for states above the high-energy cutoff nmax or emax
Ci         :        0 constrain sigii to be > asig+bsig*e
Ci         :        1 constrain sigii to be = asig+bsig*e
Ci         :        2 constrain sigii to be > asig and < bsig
Ci         :        3 constraint same as mode 1.
Ci         :          Mode 3 differs in that the least-squares fit to
Ci         :          sigii (for informational purposes only, to help
Ci         :          estimate asig and bsig) is done for states between
Ci         :          efit and nmax or emax
Ci         : arg 2: nmin : sigma for states 1..nmin are approximated
Ci         :               by sigii
Ci         : arg 3: emin : (used only if nmin<0)
Ci         :             : sigma for levels e<emin are approximated
Ci         :               by sigii
Ci         : arg 4: nmax : sigma for levels i>nmax are approximated by
Ci         :               sigii AND constrained according to mode
Ci         : arg 5: emax : (used only if nmax<=0)
Ci         :             : sigma for levels e<emax are approximated by
Ci         :               sigii AND constrained according to mode
Ci         : arg 6: asig : constraint used to approximate
Ci         :               sigii = asig + E * bsig
Ci         : arg 7: bsig : constraint used to approximate
Ci         :               sigii = asig + E * bsig
Ci         : arg 8: efit : fit sigii between efit and emax
Ci         : arg 9: lwrite:write sigii to 'sigii' file
Co Outputs
Co   a     :(4s bit mode set) constant coefficient
Co   b     :(4s bit mode set) linear coefficient
Cr Remarks
Cr
Cu Updates
Cu    4 Nov 03 Write SE to sigii file also for states not in fit
Cu   25 Jun 03 Update to incorporate sigp mode 3
Cu   21 May 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ndimh,nsmidb
      double precision e(ndimh),de(ndimh),a,b,sigp(10)
C ... Local parameters
      logical lhigh,lfit
      integer i,nmin,nmax,modsgp,fopna,ifi,lwrite
      double precision norm(2,2),rhs(2),tnorm(2,2),emin,emax,elo,ehi,
     .efit
      save norm,rhs,elo,ehi
      data norm /4*0d0/,rhs/2*0d0/,elo/9d0/,ehi/-9d0/

      if (mod(mode,2) .eq. 1) then
        call dpzero(norm,4)
        call dpzero(rhs,2)
      endif
      if (mode .eq. 1) return

      modsgp = nint(sigp(1))
      nmin   = nint(sigp(2))
      emin   = sigp(3)
      nmax   = nint(sigp(4))
      emax   = sigp(5)
C     asig   = sigp(6)
C     bsig   = sigp(7)
      efit   = sigp(8)
      lwrite = nint(sigp(9))

      if (mod(mode/2,2) .eq. 1) then
        if (lwrite .eq. 1) then
          ifi = fopna('sigii',-1,0)
          write(ifi,'(''#  state    E(lda)      dE(sig)'')')
        endif
        do  i = 1, ndimh
          lhigh = (nmax .le. 0 .and. e(i) .gt. emax) .or.
     .    (nmax .gt. 0 .and. i .gt. nmax)
          if (modsgp .eq. 3) then
            lfit = e(i) .gt. efit .and. .not. lhigh
          else
            lfit = e(i) .gt. efit
          endif
          if (lfit) then
            norm(1,1) = norm(1,1) + 1
            norm(1,2) = norm(1,2) + e(i)
            norm(2,1) = norm(2,1) + e(i)
            norm(2,2) = norm(2,2) + e(i)**2
            rhs(1) = rhs(1) + de(i)
            rhs(2) = rhs(2) + e(i)*de(i)
            elo = min(elo,de(i))
            ehi = max(ehi,de(i))
            if (lwrite .eq. 1) then
              ifi = fopna('sigii',-1,0)
              write(ifi,'(i6,2f12.6)') i,e(i),de(i)
            endif
          else
            if (lwrite .eq. 1) then
              ifi = fopna('sigii',-1,0)
              write(ifi,'(''#'',i5,2f12.6)') i,e(i),de(i)
            endif
          endif
        enddo
      endif

      if (mod(mode/8,2) .eq. 1 .or. mod(mode/4,2) .eq. 1)
     .call info0(30,0,0,' ')


      if (mod(mode/8,2) .eq. 1) then
        call info2(30,0,0,' hambls: smallest nmax encountered for '//
     .  'sigm middle block = %i',nsmidb,0)
      endif

      if (mod(mode/4,2) .eq. 1) then
        if (norm(1,1) .ne. 0) then
          call dinv22(norm,tnorm)
          a = tnorm(1,1)*rhs(1) + tnorm(1,2)*rhs(2)
          b = tnorm(2,1)*rhs(1) + tnorm(2,2)*rhs(2)
          call info5(30,0,0,' hambls: sig(low,high) = %,4;4d,%,4;4d'//
     .    '  fit : %,4;4d + %,4;4d * E(lda) (%i points)',
     .    elo,ehi,a,b,nint(norm(1,1)))
          if (lwrite .eq. 1) then
            ifi = fopna('sigii',-1,0)
            write(ifi,345) modsgp,nmin,emin,nmax,emax,a,b,efit
  345       format('# mode nmin    emin   nmax    emax',8x,
     .      'asig      bsig',8x,'efit'/2i5,f10.4,i5,f10.4,2x,2f10.4,
     .      2x,f10.4)
          endif
        endif
      endif

      end subroutine phmbl3


