c$$$      module m_nbandgw
c$$$      real(8),allocatable,private:: qq(:,:,:)
c$$$      integer,allocatable,private:: ispp(:),intqxx(:,:)
c$$$      integer,private::nn
c$$$      contains
c$$$      subroutine getnbandgw(qin,isp,ginv,ntqxx) !june2009
c$$$c     obtain ntqxx for given q and isp
c$$$      implicit none
c$$$      real(8):: a,qin(3),ginv(*),qqq(3)
c$$$      integer:: isave=0,idummy,i,isp,ntqxx,ifi,ixx,iopen,iclose,indxq,iqindx,ispp,iii
c$$$      if(isave==0) then
c$$$         ifi= iopen ('NBANDGW',1,0,0)
c$$$         nn=0
c$$$         do
c$$$            read(ifi, *,end=1011) a
c$$$            nn=nn+1
c$$$         enddo
c$$$ 1011    continue
c$$$         rewind ifi
c$$$         allocate(qq(3,nn,2),intqxx(nn,2))
c$$$         do i=1,nn
c$$$            read(ifi, *,end=1011) qqq,ispp,idummy,iii
c$$$            qq(:,i,ispp)   = qqq
c$$$            intqxx(i,ispp) = iii
c$$$            write(6,*) 'xxx: ',i,qq(:,i,ispp),ispp,intqxx(i,ispp)
c$$$         enddo
c$$$         ixx= iclose('NBANDGW')
c$$$         isave=1
c$$$      endif
c$$$      print *,' qin indx=',qin
c$$$      ixx = iqindx(qin, ginv,qq(:,:,isp),nn)
c$$$      ntqxx = intqxx( ixx,isp)
c$$$      print *,' intqxx  ntqxx=',ntqxx,ixx
c$$$      end subroutine
c$$$      integer function indxq(qin,isp, q,ispp,nn)
c$$$      real(8):: qin(3),q(3,nn)
c$$$      integer:: nn,i,ispp(nn),isp
c$$$      do i=1,nn
c$$$         if(sum( abs(q(:,i)-qin(:)) ) <1d-8 .and. isp==ispp(i) ) then
c$$$            indxq=i
c$$$            return
c$$$         endif
c$$$      enddo
c$$$      stop 'indxq: sxcf_fac2.sc.F. can not find iddex for q'
c$$$      end function
c$$$      end module m_nbandgw

!! ----------------------------------------------------------------------------
c      subroutine sxcf_fal3_scz(wqfac,kount,qip,itq,ntq,ef,esmr,tiat,miat,
      subroutine sxcf_fal3_scz(kount,qip,itq,ntq,ef,esmr,tiat,miat,
c     i                  ifev, !removed at Jan2004
     i     nsp,isp,             !ifcphi, !12jan2004 ifrb,ifcb,ifrhb,ifchb,
c     i                  ifrcw,ifrcwi, 
     i     qbas,ginv,
     i     qibz,qbz,wk,nstbz,nstar,irkip,nrkip, !irkip,nrkip 25Sep2012
c     
     i     iclass,mdim,nlnmv,nlnmc,
     i     icore, ncore,imdim,
     i     ppb,                 !pdb,dpb,ddb,
     i     freqx,wx,dw,         ! !expa,ua,deltaw, !sf 23May02   no 'freq'
     i     ecore,
c     
     d     nlmto,nqibz,nqbz,nctot,
c     i          index_qbz, n_index_qbz, !12jan2004
     d     nl,nnc,nclass,natom,
     d     nlnmx,mdimx,nbloch,ngrp,niw,nq, !niwx removed
c     
c     &     nblochpmx ,ngpn,ngcni,ngpmx,ngcmx,geigB,ngvecpB,ngveccBr, !12jan2004
     &     nblochpmx ,ngpmx,ngcmx, !ngveccBr,
     &     wgt0,nq0i,q0i,symgg, alat, shtvg,nband, ifvcfpout,
c     &     shtw,  
     &     exchange,screen,cohtest, ifexsp,
c     tetra
c     tetra     &     wtet,wtetef,
c     tetra    &     ntqx,ibzx,tetra,
c     
c     i omega,iSigma_en, !iwini,iwend,  !sf 22May02
     i     iSigma_en,           !iwini,iwend,  !sf 22May02
     i     nbmx,ebmx,           !takao 18June2003 
     i     pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat 
     i     wklm,lxklm,          !mar2012 !wklm apr2012takao
     i     invg,il,in,im,nlnm,nn_,lx,nx_,nxx_,cgr,ppbrd, !sep 2012
     i    diagonly, eftrue,
     o     zsec,coh,nbandmx)
!! jan2013 eftrue is added.
!! ef=eftrue(true fermi energy) for valence exchange and correlation mode.
!! but ef is not the true fermi energy for core-exchange mode.
      use m_readqg
      use m_readeigen,only: readeval,readcphi,readgeig
      use keyvalue,only: getkeyvalue
c      use m_nbandgw,only: getnbandgw
      use m_mpi, only: mpi__rank 
      implicit none

!! == Calcualte full simga_ij(e_i)= <i|Re[S](e_i)|j> ==
!!  exchange=T : Calculate the exchange self-energy
!!            =F : Calculate correlated part of the self-energy
!! output
!!  zsec: S_ij= <i|Re[S](e_i)|j>
!!        Note that S_ij itself is not Hermite becasue it includes e_i.
!!        i and j are band indexes
!!
!! coh : this is dummy now
!! ---------------------------------------------------------------------
!---  from Jan2006!
!     "zsec from im-axis integral part"  had been symmetrized as
!     &        wtt*.5d0*(   sum(zwzi(:,itp,itpp))+ !S_{ij}(e_i)
!     &        dconjg( sum(zwzi(:,itpp,itp)) )   ) !S_{ji}^*(e_j)= S_{ij}(e_j)
!     However, I now do it just the 1st term.
!     &        wtt* sum(zwzi(:,itp,itpp))   !S_{ij}(e_i)
!     This is OK because the symmetrization is in hqpe.sc.F
!     Now zsec given in this routine is simply written as <i|Re[S](e_i)|j>.
!     ( In the version until Jan2006 (fpgw032f8), only the im-axis part was symmetrized.
!     But it was not necessary from the begining because it was done in hqpe.sc.F
!     
!     (Be careful as for the difference between
!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!     ---because e_i is included.
!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
!----------------------------------------------------------------------
C     r 2000 takao kotani. This sxcf is starting from sec.f F.Aryasetiawan.
C----------------------------------------------------------------------


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! original document for sce.f (correlation case) is kept. This is by F.Aryasetiawan.
!! Current GW code is very different now, but kept here as a reference.
c     92.02.24
c     93.10.18 from sec.f modified to take into account equivalent atoms
c     calculates the correlated part of the self-energy SE
c     SEc(q,t,t') = <psi(q,t) |SEc| psi(q,t'>
c     SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >

c     the zeroth order Green function
c     G(r,r';w)   = S[occ]   psi(kn,r) psi(kn,r')^* /(w-e(kn)-i*delta)
c     + S[unocc] psi(kn,r) psi(kn,r')^* /(w-e(kn)+i*delta)

c     the screened coulomb potential
c     Wc(r,r';w)  = W(r,r';w) - v(|r-r'|)
c     = < [r1,r2] v(|r-r1|) X(r1,r2;w) v(|r2-r'|) >
c     W(r,r';w)   = < [r''] ei(r,r'';w) v(|r''-r'| >
c     ei          = e^(-1), inverse dielectric matrix
c     = 1 + vX
c     e           = 1 - vX0 in RPA

c     expand Wc(r,r';w) in optimal product basis B
c     Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch]
c     B(k,i,r) Wc(k,w)(i,j) B(k,j,r')^*
c     Wc(k,w)(i,j) are  the matrix elements of Wc in B

c     SEc(q,t,t') = S[k=FBZ] S[n=occ]   S[i,j=1,nbloch]
c     <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c     (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)-i*delta)>
c     
c     + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c     <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c     (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)+i*delta)>

c     the analytic structure of GWc for w .le. ef
c     |
c     |   o = pole of G
c     ^   x = pole of Wc
c     |
c     |   ef-w
c     |----<-----
c     |          |
c     o  o  o  o  o |o  o  o   ^
c     x  x  x  x  x  x|          |
c     -----------------------------|---->------------------------------
c     |x  x  x  x  x  x  x  x
c     |              o  o  o  o  o
c     |       <----->
c     ^        gap in insulator
c     |
c     |

c     the analytic structure of GWc for w .gt. ef
c     |
c     |   o = pole of G
c     |   x = pole of Wc
c     |
c     gap in insulator      ^
c     <----->        |
c     o  o  o  o               |
c     x  x  x  x  x  x  x  x|
c     ------------------------>----|-----------------------------------
c     |           |x  x  x  x  x  x  x  x
c     ^   o  o  o  o  o  o  o
c     |           |
c     ------<----|
c     w-ef    |
c     ^
c     |

c     integration along the real axis from -inf to inf is equivalent to
c     the integration along the path shown

c------------------------------------------------------------
c     integration along the imaginary axis: wint (s. also wint.f) (takao ->wintz)
c     (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
c     the i*delta becomes irrelevant
c------------------------------------------------------------
c     
c     omit k and basis index for simplicity and denote e(q-k,n) = e
c     wint = (i/2pi) < [w'=-inf,inf] Wc(w')/(w+w'-e) >
c     
c     w' ==> iw', w' is now real
c     wint = - (1/pi) < [w'=0,inf] Wc(iw') (w-e)/{(w-e)^2 + w'^2} >
c     
c     transform: x = 1/(1+w')
c     this leads to a denser mesh in w' around 0 for equal mesh x
c     which is desirable since Wc and the lorentzian are peaked around w'=0
c     wint = - (1/pi) < [x=0,1] Wc(iw') (w-e)x^2/{(w-e)^2 + w'^2} >
c     
c     the integrand is peaked around w'=0 or x=1 when w=e
c     to handel the problem, add and substract the singular part as follows:
c     wint = - (1/pi) < [x=0,1] { Wc(iw') - Wc(0)exp(-a^2 w'^2) }
c     * (w-e)/{(w-e)^2 +w'^2}x^2 >
c     - (1/2) Wc(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c     
c     the second term of the integral can be done analytically, which
c     results in the last term
c     a is some constant
c     
c     when w = e, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w') and
c     the integral becomes -Wc(0)/2
c     this together with the contribution from the pole of G (s.u.)
c     gives the so called static screened exchange -Wc(0)

c--------------------------------------------
c     contribution from the poles of G: SEc(pole)
c--------------------------------------------
c     
c     for w .le. ef
c     SEc(pole) = - S[k=FBZ] S[n=occ] S[i,j=1,nbloch]
c     <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c     Wc(k,e(q-k,n)-w)(i,j) theta(e(q-k,n)-w)
c     
c     for w .gt. ef
c     SEc(pole) = + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c     <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c     Wc(k,w-e(q-k,n))(i,j) theta(w-e(q-k,n))
c     
c     theta(x)  = 1   if x > 0
c     = 1/2 if x = 0
c     = 0   if x < 0

c     FBZ = 1st BZ
c     NOTE: the routine only calculates the diagonal elements of the SE
c     i.e. SEc(q,t)

c     q       = q-vector in SEc(q,t)
c     itq     = states t at q
c     ntq     = no. states t
c     eq      = eigenvalues at q
c     ef      = fermi level in Rydberg
c     tr      = translational vectors in rot*R = R' + T
c     iatomp(R) = R'
c     ifrw,ifcw,ifrwi,ifcwi
c     = direct access unit files for Re and Im coulomb matrix
c     along real and imaginary axis
c     ifrb,ifcb,ifrhb,ifchb
c     = direct access unit files for Re and Im b,hb
c     qbas    = base reciprocal lattice vectors
c     ginv    = inverse of qbas s. indxrk.f
cxxxxxippb,ipdb,idpb,iddb = pointers to work array w for
c     ppb     = <phi(RLn) phi(RL'n') B(R,i)>
c     pdb     = <phi(RLn) phidot(RL'n') B(R,i)>
c     dpb     = <phidot(RLn) phi(RL'n') B(R,i)>
c     ddb     = <phidot(RLn) phidot(RL'n') B(R,i)>
c     freq    = frequencies along real axis
c     freqx   = gaussian frequencies x between (0,1)
c     freqw   = (1-freqx)/freqx
c     wx      = weights at gaussian points x between (0,1)
c     ua      = constant in exp(-ua^2 w'^2) s. wint.f
c     expa    = exp(-ua^2 w'^2) s. wint.f
c     dw      = frequency mesh along real axis
c     deltaw  = energy mesh in SEc(qt,w) ---Not used now
c     iclass  = given an atom, tells the class
c     wk      = weight for each k-point in the FBZ
c     indexk  = k-point index
c     qbz     = k-points in the 1st BZ
c     nstar   = no. stars for each k
c     irkip(k,R,nq) = gives index in the FBZ with k{IBZ, R=rotation
c     mdim    = dimension of B(R,i) for each atom R
c     work arrays:
c     rbq,cbq     = real and imaginary part of b(q)
c     rhbq,chbq   = real and imaginary part of hb(q)
c     rbkq,cbkq   = real and imaginary part of b(q-k)
c     rhbkq,chbkq = real and imaginary part of hb(q-k)
c     b is the eigenvector of the LMTO-Hamiltonian
c     ekq     = eigenvalues at q-k
c     rmel,cmel = real and imaginary part of
c     <psi(q,t') | psi(q-k,t) B(k,R,i)>
c     wr1 ... = work arrays
c     dimensions:
c     nqibz   = number of k-points in the irreducible BZ
c     n1,n2,n3= divisions along base reciprocal lattice vectors
c     natom   = number of atoms
c     nctot   = no. allowed core states
c     nbloch  = total number of Bloch basis functions
c     nlnmx   = maximum number of l,n,m
c     nlmto   = total number of LMTO basis functions
c     ngrp    = no. group elements (rotation matrices)
c     niw     = no. frequencies along the imaginary axis
c     nw      = no. frequencies along the real axis
c     
c     secq(t) = <psi(q,t) |SEc| psi(q,t)>
c----------------------------------------------------------------------
      integer :: ntq, natom,nqbz,nqibz,ngrp,nq,nw,niw,
     &     nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,
     &     irot,  iqisp,ikpisp,isp,nsp,  nlnmx,
c     &  iq, idxk,iwr1,iwr2,iwr3,iwr4,iwc1,iwc2,iwc3,iwc4
     &     ip, it,itp, ifcphi,  ! ifrb,ifcb,ifrhb,ifchb, 
     i     iiclass,mdim(*),
     i     ifrcw,ifrcwi,        !iindxk,
     &     ifvcfpout,ndummy1,ndummy2,kx,kr,ngc,ngb,nbloch,n1,n2,n3, 
     &     kp,nt0,nocc, nt0p,nt0m,irkp,i,j,nt0org,nmax,nt,
     &     nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,ixs,ixsmx,
     &     mdimx, 
     &     itq(ntq),iatomp(natom),miat(natom,ngrp),
     &     nstar(nqibz),irkip(nqibz,ngrp,nq),nrkip(nqibz,ngrp,nq),kount(nqibz,nq),  ntqxx
c     
      real(8) :: q(3),qbas(3*3),ginv(3*3),tiat(3,natom,ngrp),tr(3,natom),
     &     wk(nqbz),qibz(3,nqibz),qbz(3,nqbz),
     &     freqx(niw),wx(niw),expa_(niw),
     &     ekq(nband), ekc(nctot+nband),
     &     tpi,ef,esmr,efp,efm,wfac,wfacx,we,esmrx,ua_,
     &     dw,wtt,wexx,www , wfacx2,weavx2

      integer :: ngpmx, ngcmx, ngpn(nqbz), !ngcni(nqibz), 
     &     ngvecpB(3,ngpmx,nqbz), igc, !ngveccBr(3,ngcmx,nqibz),
     &     nadd(3)
      real(8) :: wgt0(nq0i,ngrp),qk(3),qfbz(3),
     &     qbasinv(3,3), det,qdiff(3),add(3),symgg(3,3,ngrp),symope(3,3),
     &     qxx(3),q0i(1:3,1:nq0i),shtvg(3,ngrp),shtv(3),alat,ecore(nctot),
     &     ppb(*),              !pdb(1),dpb(1),ddb(1), !*
     &     coh(ntq,nq)
      complex(8):: geigB(ngpmx,nband,nqbz),  alagr3z2,wintz,alagr3z
c     
      logical :: exchange,screen,cohtest
      real(8),allocatable:: freq_r(:),
     &     w1p(:,:,:),w2p(:,:,:) !,w3p(:,:)   
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
      complex(8),allocatable ::w3p(:,:,:)   
      logical :: debug, debugp,debug2=.false.,initp=.true.

      real(8),allocatable:: rmelt(:,:,:),cmelt(:,:,:),drealzzzmel(:,:,:),dimagzzzmel(:,:,:)
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzzmel(:,:,:),
     &     zw (:,:), zwz(:,:,:,:), zwz0(:,:,:),zwzi(:,:,:),zwm(:,:),zwix(:,:) ,
     &     zwzix(:,:,:) 
      complex(8) :: wintzav,wintzsg_npm

      integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen,npm
      integer(4),save::ifzwz=-999

      integer(4) :: ia,iSigma_en 
      real(8)    :: omg, omega(ntq),omega0    
      complex(8) :: zsec(ntq,ntq,nq)          
      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)

      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)
      integer(4) :: nt_max, igb1,igb2,iigb, itpp,ix0 
      complex(8),allocatable:: zmel_(:,:,:),zmel1(:),zmel3(:)
      complex(8), allocatable :: zw_(:,:),zzmel(:,:)
      complex(8), allocatable :: zwz2(:,:),zw2(:,:,:),zmel2(:,:)
      complex(8) ::  zz2, zwz3(3)
      complex(8), allocatable :: zmel1_(:,:,:),zz3(:), zw3(:,:,:),zw3x(:,:)
      complex(8), allocatable :: zwz4(:,:),zwz44(:,:)

      real(8)::weavx,wfaccut=1d-10

      complex(8)::zdotu
      complex(8), allocatable :: CC(:,:,:)
      complex(8) :: zdotc

      logical :: GaussSmear,gass !external
      real(8) :: ebmx(2),ddw
      integer(4):: nbmx(2),nbmxe,nstatetot

c     integer(4):: n_index_qbz
c     integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4):: bzcase,iqini,nstbz(nqbz),iqend
      real(8):: wgtq0p

      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      integer(4):: nrec
c--------------------------------------------------------------------
      integer(4)::verbose
      logical :: onlyimagaxis,ua_auto
      real(8):: freqw,ua2_(niw),ratio,freqw1
      real(8),allocatable:: uaa(:,:)

      integer(4) ::nbcut,nbcutc
C...  
      logical::smbasis
      integer(4):: nn,no,ifpomat,iclose,isx,iqx,ntp0
      complex(8),allocatable:: pomat(:,:), zwzs(:)
      real(8):: q_r(3)
      integer(4):: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
      complex(8):: pomatr(nnmx,nomx,nkpo),zwzs0
      real(8):: qrr(3,nkpo)

      logical:: timemix =.false., oncew
      real(8)   :: fffr(3)
      integer(4)::nwp,nw_i,ixsmin,iwp,iwm,isxmin,nwxi,nwx,iir

      integer(4)::  iix,ipp,  itini,itend
      logical :: test_symmetric_W,onceww

!! july2012
      integer:: lxklm,ivc,ifvcoud,idummy,iy,ngb0
      real(8):: wklm((lxklm+1)**2),pi,fpi,vc,qvv(3),aaaa
      complex(8)::zmelt1,zmelt0,zmelt2
      real(8)::voltot
      logical :: newaniso,testomitq0
      complex(8),allocatable:: ppovl(:,:),ppovlz(:,:),zcousq(:,:),zmeltt(:,:,:)
      real(8),allocatable::vcoud(:),vcousq(:)
      integer:: invg(ngrp),invrot,invr

      integer:: il(*),in(*),im(*),nlnm(*),nn_,lx(*),nx_(*),nxx_
      real(8)::  cgr(*),ppbrd(*)

      integer:: nprecx,mrecl,ifwd,nrot
      character(5):: charnum5

      real(8):: eq(nband),qip(3,*)

      logical ::eibz4sig
      integer::nbandmx(nq)      !sep2012 number of bands for calcualted self-energy.

      logical:: diagonly
c      real(8)::wqfac(nqbz)
      integer:: irot2
      real(8):: eftrue

      integer::ierr 

c      real(8):: emax_d_sigm
c      integer:: ntqxxd,if_ntqxxd,iopenxx
c      character(15):: fffn
c      integer:: inii,iendi
!$    integer:: omp_get_thread_num
c--------------------------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,'wklm=0 test xxxxxxxxxxxxxx'
c      wklm=0d0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!TIME0 
!TIME0 
c      call getkeyvalue("GWinput","emax_d_sigm", emax_d_sigm, default=1d10)
c      print *,'emax_d_sigm=',emax_d_sigm

      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
      if(verbose()>=90) debug=.true.
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      nbcutc=nctot+nbcut
      if(debug) then
         print *, ' sxcf_fal3_scz:'
         if(GaussSmear()) then
            print *,' new version: gaussian smearing!!!'
            print *,' wintzsg wfacx weavx are related parts...'
            gass=GaussSmear()
            write(6,*)' sxcf_fal3_scz: new version: gaussian smearing=',gass
         endif
      endif
      tpi         = 8d0*datan(1d0)
      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
      call dinv33(qbas,0,qbasinv,det)

      if(debug) print *, ' sxcf: 1'
      allocate(expikt(natom))
cccccccccccccccccccccc
c      print *, 'sssssssss',sizeof(expikt)
ccccccccccccccccccccccc

!TIME1 "before 1001"

!!== ip loop to spedify external q ==
      do 1001 ip = 1,nq         
         if(sum(irkip(:,:,ip))==0) cycle ! next ip
         write (*,*) ip,'  out of ',nq,'  k-points(extrnal q) ' ! call cputid  (0)
c     call readeval(q(1,ip),isp,eq)
!TIME0 
         q(1:3)= qip(1:3,ip)
         call readeval(q,isp,eq)
         do i  = 1,ntq
            omega(i) = eq(itq(i)) !sf 21May02 + 2d0*(dble(iw)-shtw)*deltaw
         enddo

c     if(nctot==0) then !this is to override -Mbounds check
c     deallocate(ecore)
c     allocate(ecore(1)) !this is dummry
c     endif
c     if(debug) write(6,*) "goto sxcf_fal3_scz"
c     print *, ' sxcf_fal3_scz:'
c     iq          = idxk (q,qbz,nqbz) ! index for q

c     direct access read e,b and hb at q
c     iqisp =isp + nsp*(iq-1)
c     call rwdd1   (ifcphi, iqisp, 2*nlmtobnd, cphikq)
         call  readcphi(q, nlmto,isp, quu, cphikq)
         if(debug) print *, ' sxcf: 2'

c     b and hb of states {t} at q
         do      it = 1,ntq
            itp        = itq(it)
            cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
         enddo
         if(debug) print *, ' sxcf: 3'

         if(bzcase()==1) then
            if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
            if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
         endif

c-----
         if(exchange.and.(.not.newaniso())) then
            rewind  ifvcfpout
            read(ifvcfpout) ndummy1, ndummy2
         endif

c     --- gauss_img : interpolation gaussion for W(i \omega).
         call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
         if(ua_<=0d0) then
            ua_auto =.true.
            write(6,"(' ua_auto=T')")
         else
            ua_auto =.false.
            do ix = 1,niw
               freqw     = (1d0 - freqx(ix))/ freqx(ix)
               expa_(ix) = exp(-(ua_*freqw)**2)
            enddo
         endif
c     print *,' ua_auto ua_=',ua_auto,ua_
c     if(newaniso()) ifvcoud = iopen('Vcoud',0,0,0)

         if(debug) then
            do ix = 1,niw
               write(6,"(' expa=',10d13.6)") freqx(ix),expa_(ix)
            enddo
c     print *, 'test1 bbb'
         endif

c===============================
c     loop over irreducible k-points
c===============================
c     if(newaniso()) then
c     continue
c     else
c     if(bzcase()==1) then
c     c         print *, 'test1 bbb2'
c     c           kx = 1  ! qibz(:,1)=0 contribution for kcount
c     c          print *, 'test1 bbb3'
c     c          if(irk(kx,irot)/=0) kount(kx,ip)= kount(kx,ip) + 1
c     endif
c     endif

c     --- main loop start
         iqini=2
         if(bzcase()==2) iqini=1
         iqend= nqibz+nq0i
c     print *, 'test1 bbb4'
         if(newaniso()) then    !takao2012apr
            iqini=1
            iqend=nqibz         !no sum for offset-Gamma points.
c     if(exchange) iqend=nqibz
         endif

         if(.not.exchange) then
            ifwd = iopen('WV.d',1,-1,0)
            read (ifwd,*) nprecx,mrecl
            ifwd = iclose('WV.d')
         endif 


!TIME1 "before 1100"

!! == Big loop for kx ==
         do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
            if(sum(irkip(kx,:,ip))==0) cycle ! next kx
!TIME0 
            print *, ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
            if( kx <= nqibz ) then
               qibz_k= qibz(:,kx)
            else
               qibz_k= 0d0
            endif
            if(verbose()>=40)  print *, ' sxcf_fal2_scz: loop 1100 kx=',kx
            if(timemix) call timeshow("11111 k-cycle")

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c test
c            if(kx/=iqini) then
c             cycle !Gamma only test.
c             print *,'Gamma onlytest skip kx=',kx
c            endif
c
c            if(abs(qibz_k(3))>1d-3) then
c               cycle
c               print *,' z-axis only skip kx=',kx,qibz_k
c            endif
c
c
c            if(abs(qibz_k(1))+abs(qibz_k(2))<1d-3) then
c               print *,' z-axis skip kx=',kx,qibz_k
c               cycle
c            endif
c
c            if(abs(qibz_k(1))+abs(qibz_k(2))><1d-3) then
c               print *,' z-axis skip kx=',kx,qibz_k
c               cycle
c            endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
c     ngc = ngcni(k)  ! k-points in IBZ
            ngb = nbloch + ngc
            if(debug) print *, ' sxcf: ngb=',ngb,nbloch

!! ===Readin diagonalized Coulomb interaction===
!! note sep102012takao
!!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
!!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
!!  Vould file is written in hvccfp0.m.F.
!! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
!! These can be also separeted into WVR.ID and WVI.ID files.
            if(newaniso()) then
               if(kx<=nqibz) qxx=qibz_k
               if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
               ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
               do
                  read(ifvcoud) ngb0
                  read(ifvcoud) qvv
c     print *,'readin qvv=',qvv
                  if(allocated(vcoud)) deallocate(vcoud)
                  allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
                  read(ifvcoud) vcoud


c$$$cccccccccccccccccccccccccccccccccccccccccccccc
c$$$                  if(sum(abs(qxx(1:2)))<1d-4) then
c$$$                  do irot2 = 1,ngrp
c$$$                     kr = irkip(kx,irot2,ip) ! index for rotated kr in the FBZ
c$$$                     if(kr==0) cycle ! next irot
c$$$                     vcoud(1) =vcoud(1)*wqfac(kr)
c$$$                     exit
c$$$                  enddo
c$$$                  print *,'wwwwwwwwwww why here ?wwwwwwwwww'
c$$$                  endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccc



                  read(ifvcoud) zcousq
                  if(sum(abs(qvv-qxx))<1d-6) goto 1133
               enddo
               if(sum(abs(qvv-qxx))>1d-6) then
                  print *,'qvv =',qvv
                  print *,'qxx=',qxx,kx
                  call rx( 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent')
               endif
 1133          continue
               if( ngb0/=ngb ) then !sanity check
                  print *,' qxx ngb0 ngb=',qxx,ngb0,ngb
                  call rx( 'hsfp0.m.f:ngb0/=ngb')
               endif
c$$$  if(sum(abs(qibz_k))<1d-6) then
c$$$  idummy  = iclose('Vcoud') !close and open again. This is because first data in Voud is for q=0
c$$$  ifvcoud = iopen('Vcoud',0,0,0)
c$$$  endif
c$$$  read(ifvcoud) ngb0
c$$$  if( ngb0/=ngb ) then
c$$$  print *,' qibz_k=',qibz_k,ngb0,ngb
c$$$  print *,' qibz_k=',qibz_k
c$$$  stop 'hsfp0.m.f:ngb0/=ngb'
c$$$  endif
c$$$  read(ifvcoud) qvv
c$$$  if(sum(abs(qvv-qibz_k))>1d-6) then
c$$$  print *,'qvv =',qvv
c$$$  print *,'qibz_k=',qibz_k,kx
c$$$  stop 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent'
c$$$  endif
c$$$  if(allocated(zcousq)) deallocate(zcousq,vcousq,vcoud)
c$$$  allocate( zcousq(ngb0,ngb0),vcousq(ngb0),vcoud(ngb0) )
c$$$  read(ifvcoud) vcoud
c$$$  read(ifvcoud) zcousq
c$$$  vcousq=sqrt(vcoud) !

!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!! ppovlz= O Zcousq
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
               if(allocated(ppovlz)) deallocate(ppovlz)
               allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
c               print *,'sssssssssssssss',sizeof(ppovlz)/1000000.
c               print *,'sssssssssssssss',sizeof(ppovlz(1,1))
c               print *,'sssssssssssssss',sizeof(ppovlz(1,1)),sizeof(reshape(ppovlz,(/1/)))
               call readppovl0(qibz_k,ngc,ppovl)
c nnnnnnnnnnn comment out mulipled by ppovl. 
               ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
               ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:)) !may2013takao not ppovl inlcuded.
c               ppovlz(1:nbloch+ngc,:) = zcousq(1:nbloch+ngc,:) !may2013takao not ppovl inlcuded.
               deallocate(zcousq,ppovl)
            endif

!! read Coulomb interaction for exchange for newaniso=F
            if(exchange) then
c     if(testomitq0()) then
c     print *,'xxxxxxxxxxx takao test omit q0i aaaaa'
c     if(allocated(vcoul)) deallocate(vcoul)
c     endif
               if(.not.newaniso()) then 
                  read(ifvcfpout) nn !oct2005
                  allocate(vcoul(nn,nn))
                  read(ifvcfpout) vcoul(1:nn,1:nn)
               endif  
            endif

!! === open WVR,WVI ===
            if(.not.exchange) then
               ifrcw  = iopen('WVR.'//charnum5(kx),0,-1,mrecl)
               ifrcwi = iopen('WVI.'//charnum5(kx),0,-1,mrecl)
            endif

            nrot=0
            do irot = 1,ngrp
               if( kx <= nqibz) then
                  kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
                  if(kr==0) cycle ! next irot
                  qbz_kr= qbz (:,kr) 
               else
                  kr=-99999     !for sanity check
                  qbz_kr= 0d0
                  if( wgt0(kx-nqibz,irot)==0d0 ) cycle ! next irot
               endif
               nrot=nrot+1
            enddo  
!TIME1 "before 1000"
            
!! === loop 1000 over rotations irot ===
            do 1000 irot = 1,ngrp
!TIME0
               if( kx <= nqibz) then
                  kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
                  if(kr==0) then
                     cycle      ! next irot
                  end if
                  qbz_kr= qbz (:,kr) 
               else
                  kr=-99999     !for sanity check
                  qbz_kr= 0d0
                  if( wgt0(kx-nqibz,irot)==0d0 ) then
                     cycle      ! next irot
                  end if
               endif
c     write(*,"('ip,kx irot=',3i5, 'out of',2i4)") ip,kx,irot, iqend,ngrp
c               write(*,"('ip,kx irot=',3i5, ' out of(ipend,nrot)',3i4)") ip,kx,irot, iqend,nrot,ngrp

c$$$  if( kx <= nqibz ) then
c$$$  kr = irk(kx,irot) ! index for rotated k in the FBZ
c$$$  c          qibz_k= qibz(:,kx)
c$$$  if(kr/=0) qbz_kr= qbz (:,kr) !feb2006
c$$$  else
c$$$  kr=-99999 !for sanity check
c$$$  qibz_k= 0d0
c$$$  qbz_kr= 0d0
c$$$  endif

c     if( sum(abs( irk(:,irot) )) ==0 .and.
c     &           sum(abs( wgt0(:,irot))) == 0d0 ) then
c     cycle
c     endif
c     call cputid (0)
c     if (irot == 1 .or. irot == ngrp) then
c     call cputid(0); write(*,*)' ppba '
c     endif

!!    rotate atomic positions invrot*R = R' + T
               invr  = invrot (irot,invg,ngrp)
               tr    = tiat(:,:,invr)
               iatomp= miat(:,invr)
               symope= symgg(:,:,irot)
               shtv  = matmul(symope,shtvg(:,invr))
!TIME1 "before ppbafp_v2"
!TIME0 
c     -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
               call ppbafp_v2 (irot,ngrp,isp,nsp,
     i              il,in,im,nlnm, !w(i_mnl),
     d              nl,nn_,nclass,nlnmx,
     i              mdimx,lx,nx_,nxx_, !Bloch wave    
     i              cgr, nl-1,  !rotated CG
     i              ppbrd,      !radial integrals
     o              ppb)

!! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
               qk =  q - qbz_kr !qbz(:,kr)
!! Read eigenvalue and Cphi at qk
               call readeval(qk, isp, ekq) !jan2004
               call readcphi(qk, nlmto,isp, quu, cphikq)

!! no. occupied (core+valence) and unoccupied states at q-rk
               ekc(nctot+1:nctot+nband) = ekq (1:nband)
               nt0 = nocc (ekc,ef,.true.,nstatetot)
               ddw= .5d0
               if(GaussSmear()) ddw= 10d0
               efp= ef+ddw*esmr
               efm= ef-ddw*esmr
               nt0p = nocc (ekc,efp,.true.,nstatetot)
               nt0m = nocc (ekc,efm,.true.,nstatetot)
c     tetra
c     if(tetra) then
c     irkp = ibzx(kp)  ! IBZ corresponding to qk
c     do i=1,nband
c     if( wtetef(i,irkp) == 0d0 ) exit
c     enddo
c     if(sum(abs(wtetef(i:nband,irkp)))/=0d0)
c     &   stop ' sxcf: wtetef sum err1'
c     nt0org = nocc (ekc,ef,.true.,nstatetot)
c     nt0 = i-1 + nctot
c     endif

!! nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!! nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W
               if(exchange) then
                  nbmax = nt0p-nctot
               else
                  nbmax = nband
                  nbmxe = nocc (ekc,ebmx(2),.true.,nstatetot)-nctot
                  nbmax  = min(nband,nbmx(2),nbmxe)
                  if(initp) then
                     print *,' nbmax=',nbmax
                     initp=.false.
                  endif
               endif
               if(debug) print *, ' sxcf: 1100c'
!! ntqxx is number of bands for <i|sigma|j>.
               ntqxx = nocc (omega-eftrue,ebmx(1),.true.,ntq) 
                             !bug -ef is added jan2013
                             !previous version do not give wrong results, but inefficient.
               ntqxx = min(ntqxx, nbmx(1))
               if(ntqxx<nband) then
                  do i=ntqxx,1,-1 !redudce ntqxx when band tops are degenerated. !sep2012
                     if(omega(i+1)-omega(i)<1d-2) then
                        ntqxx=i-1
                     else
                        exit
                     endif
                  enddo
               endif  
               nbandmx(ip)=ntqxx !number of bands to be calculated Sep2012.
c$$$ccccccccccccccccccccccccccccccccccccccccccccc
c$$$c sep2013t
c$$$               ntqxxd = nocc (omega-eftrue,emax_d_sigm,.true.,ntq) 
c$$$                             !bug -ef is added jan2013
c$$$                             !previous version do not give wrong results, but inefficient.
c$$$               ntqxxd = min(ntqxxd, nbmx(1))
c$$$               if(ntqxxd<nband) then
c$$$                  do i=ntqxxd,1,-1 !redudce ntqxx when band tops are degenerated. !sep2012
c$$$                     if(omega(i+1)-omega(i)<1d-2) then
c$$$                        ntqxxd=i-1
c$$$                     else
c$$$                        exit
c$$$                     endif
c$$$                  enddo
c$$$               endif  
c$$$               if(isp==1) fffn='NTQXXD.up.'//charnum5(ip)
c$$$               if(isp==2) fffn='NTQXXD.dn.'//charnum5(ip)
c$$$               if_ntqxxd=iopen(fffn,1,-1,0)
c$$$               write(if_ntqxxd,"(3f13.6,i5)") q(1:3),ntqxxd
c$$$               if_ntqxxd=iclose(fffn)
c$$$               write(6,*)' ntqxxd=',ip,isp,q(1:3),ntqxxd
c$$$ccccccccccccccccccccccccccccccccccccccccccccc
               if(debug) print *,' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
               nstate = nctot + nbmax ! = nstate for the case of correlation
!TIME1 "before expikt"
!TIME0

!!====================================================
!!     matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
!!     including the phase factor exp(ik.T)
!!     B(rot*k,r) = B(k,invrot*r)
c====================================================

!! phase factors expikt(ia) is for exp(ik.T(R))
               do ia = 1,natom
                  expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
               end do
               if(debug) print *, ' sxcf: tr=',tr
               if(debug) print *, ' sxcf: goto psicb2'
!! matrix elements
!!   core
               nt   = nctot + nbmax ! = nstate for the case of correlation
               allocate( zzzmel(nbloch,nt,ntqxx))
               call psicb_v2  (icore,ncore,ntqxx,iclass,
     i              dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i              cphiq, 
     i              ppb,  
     i              nlnmv,nlnmc,mdim,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,nt,ntqxx,natom,nclass,
     d              nl,nnc,     
     o              zzzmel)
               if(debug)  print *, ' sxcf_fal2sc: goto psi2bc1'
!!   valence
               call psi2b_v2  (nbmax, ntqxx,iclass,
     i              dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i              cphikq,     !occ    q-rk nband
     i              cphiq,      !unocc  q    ntq
     i              ppb,  
     i              nlnmv,nlnmc,mdim,nctot,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx, nband, nt,ntqxx,
     d              natom,nclass,
     o              zzzmel)
               if(verbose()>50) call timeshow("4 after psi2bc1")
               if(debug2) then
                  write(6,"('sum of zmel abszmel=',4d23.16)") sum(zzzmel),sum(abs(zzzmel) )
               end if
!TIME1 "bfore psi2b_v2"
!TIME0
!! -- IPW part.
               if(debug) print *, ' sxcf_fal1: goto drvmelp2 xxx111'
               allocate(drealzzzmel(nbloch,nt,ntqxx),dimagzzzmel(nbloch,nt,ntqxx))
               drealzzzmel=dreal(zzzmel)
               dimagzzzmel=dimag(zzzmel)
               deallocate(zzzmel)
               allocate( rmelt(ngb, nctot+nbmax, ntqxx), ! nstate= nctot+nband
     &              cmelt(ngb, nctot+nbmax, ntqxx))
               call drvmelp2( q,             ntqxx, ! q in FBZ
     i              q-qbz_kr,  nbmax, ! q-rk
     i              qibz_k,     ! k in IBZ for mixed product basis. rk = symope(qibz_k) 
     i              isp,ginv,
     i              ngc,ngcmx, ngpmx,nband,itq, 
     i              symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
c     i       dreal(zzzmel), dimag(zzzmel), nbloch, nt,nctot,
     i              drealzzzmel, dimagzzzmel, nbloch, nt,nctot,
     o              rmelt,cmelt)
               if(debug) print *, ' sxcf_fal1: end of drvmelp2'
               deallocate(drealzzzmel,dimagzzzmel)
               if(verbose()>50) call timeshow("5 after drvmelp")
               if(nbcut/=0.and.(.not.exchange)) then
                  do it= nctot+1,nctot+min(nbcut,nbmax)
                     rmelt(:, it,:) =0d0
                     cmelt(:, it,:) =0d0
                  enddo
               endif
!TIME1 "after drvmelp2"

!! NOTE:=========================================
!! zmelt = rmelt(igb(qbz_kr), iocc(q), iunocc(q-qbz_kr)) + i* cmelt
!! iunocc: band index at target  q.
!! iocc:   band index at intermediate vector qk = q - qbz_kr 
!! igb: index of mixed product basis       at qbz_kr (or written as rk)
!!   igb=1,ngb 
!!   ngb=nbloch+ngc  ngb: # of mixed product basis
!!                   nbloch: # of product basis (within MTs)
!!                   ngc: # of IPW for the Screened Coulomb interaction.
!!                   igc is for given 
!! See readgeig in drvmelp2.
!! ==================================================

!TIME0
!! smbasis 
               if(smbasis()) then !
                  ntp0= ntqxx
                  nn= nnr(kx)
                  no= nor(kx)
                  allocate( pomat(nn,no) )
                  pomat= pomatr(1:nn,1:no,kx)
                  if( sum(abs(qibz_k-qrr(:,kx)))>1d-10 .and.kx <= nqibz ) then
                     call rx( 'qibz/= qrr')
                  endif
                  if(no /= ngb.and.kx <= nqibz) then
!!     A bit sloppy check only for kx<nqibz because qibze is not supplied...
                     write(6,"(' q  ngb  ',3d13.5,3i5)")  qibz_k,ngb
                     write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
                     call rx( 'x0kf_v2h: POmat err no/=ngb')
                  endif
ccccccccccccccccccccccccccc
c     write(6,*)
c     write(6,*)'kx=',kx
c     write(6,"('ngb nn=',3i4)") ngb,nn,no
c     write(6,"(' q  ngb  ',3d13.5,3i5)")  qibz_k,ngb
c     write(6,"(' q_r  nn no',3d13.5,3i5)") q_r,nn,no
c     if(ngb/=nn) stop ' sxcf_fal2 ngb/=nn'
c     goto 1113
ccccccccccccccccccccccccccc
                  if(timemix) call timeshow("xxx2222 k-cycle")
                  ngb = nn      ! Renew ngb !!!
                  allocate ( zmel  (nn, nctot+nbmax, ntp0) )
                  call matm( pomat, dcmplx(rmelt,cmelt), zmel, 
     &                 nn, no, (nctot+nbmax)*ntp0 )
                  deallocate(rmelt, cmelt)
                  allocate( rmelt(ngb, nctot+nbmax, ntp0), !ngb is reduced.
     &                 cmelt(ngb, nctot+nbmax, ntp0) )
                  rmelt = dreal(zmel)
                  cmelt = dimag(zmel)
                  deallocate(zmel,pomat)
c     1113 continue
               else
                  nn=ngb
                  no=ngb
               endif
               if( oncew() ) then
                  write(6,"('ngb nn no=',3i6)") ngb,nn,no
               endif
               if(timemix) call timeshow("22222 k-cycle")
!! === End of zmelt ; we now have matrix element zmelt= rmelt + img* cmelt === 
               if(allocated(zzzmel))deallocate(zzzmel) !rmel,cmel)
               if(debug) print *, ' sxcf: goto wtt'
               if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)

cccccccccccccccccccccccccccccccccccccccccccccccccc
c     A matrix element check in Si. See the corresponding section of x0kf.f
c     write(6,"(3f13.5)") q         !ntq
c     write(6,"(3f13.5)") qbz(:,kr) !product basis
c     print *,' sxcf: sumcheck rmelt cmelt='
c     &   ,sum(rmelt(1:ngb,4,1))
c     &   ,sum(cmelt(1:ngb,4,1))
c     &   ,sum(abs(rmelt(1:ngb,4,1)))
c     &   ,sum(abs(cmelt(1:ngb,4,1)))
c     stop ' a matrix ele check test end'
ccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccc
c     testmelt0
c     takao
c     if(nctot==0) then
c     write(6,*)"testmode: clear rmelt cmelt(:,nctot+1:nctot+5,:) "
c     rmelt(:,nctot+1:nctot+5,:)=0d0
c     cmelt(:,nctot+1:nctot+5,:)=0d0
c     endif
ccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccc
c     do  it  = 1,nstate
c     itp=1
c     it =1
c     do  itp = 1,ntq
c     print *,' ngb nbloch ngc =',ngb,nbloch,ngcni(k)
c     do  ibl  = 1,ngb
c     write(6,"(' rmelt cmelt =',3i4,4d14.6)")it,itp,ibl
c     &  ,rmelt(ibl,it,itp) ,cmelt(ibl,it,itp)
c     &  ,rmel (ibl,it,itp) ,cmel (ibl,it,itp)
c     enddo
cccccccccccccccccccccccccccccccccccccccccccccc
c
c     if(kx<= nqibz) then
c     wtx = 1d0
c     wtt = wk(kr)
c     else
c     wtx = wgt0(kx-nqibz,irot)
c     wtt = wk(1)*wgt0(kx-nqibz,irot)
c     endif
               if(bzcase()==2)then
                  if(kx<=nqibz) then
                     wtt = wk(kr)
                     if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
                  elseif(kx>nqibz) then !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
                     wtt= wgt0(kx-nqibz,irot)
                  endif
               else
                  if(kx<= nqibz) then
c     wtx = 1d0
                     wtt = wk(kr)
                  else
c     wtx = wgt0(kx-nqibz,irot)
                     wtt = wk(1)*wgt0(kx-nqibz,irot)
                     if(abs(wk(1)-1d0/dble(nqbz))>1d-10) call rx( 'sxcf:wk(1) inconsistent')
                  endif
               endif
!!
               if(eibz4sig()) then
                  wtt=wtt*nrkip(kx,irot,ip) !25sep2012
               endif   

!!--------------------------------------------------------
!! --- exchange section ---
!!--------------------------------------------------------
!TIME1 "start exchange section"
!TIME0
               if(exchange) then
c need to check following comments ----     
c     S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c     v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
c     
c     > z1p(j,n,t) = S[i=1,nbloch] <psi(q,t) | psi(q-rk,n) B(rk,i)> v(k)(i,j)
c     
c     write(6,*)' vcoulsum=',sum(vcoul)
c     if(debug) write(6,*)'  sumz=',dcmplx(rmelt,cmelt),sum(vcoul)
c     
c     --- screened exchange case
c     if(screen) then
c     allocate( zw (nblochpmx,nblochpmx))
c     ix = 1
c     ! write(*,*)(kx-2)*(nw_w+1)+ix
c     read(ifrcw,rec=((kx-2)*nw+ix)) zw  ! Readin W(0) - v         !sf 22May02
c     !nw is number of frequency points in general mesh: freq_r(nw), freq_r(1)=0
c     vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
c     deallocate(zw)
c     endif
                  allocate( zmel (ngb, nctot+nbmax, ntqxx),w3p(nctot+nbmax,ntqxx,ntqxx))
                  zmel  = dcmplx (rmelt,cmelt)
cccccccccccccccccccccccccccccc
                  if(debug) then
                     do it=1,nctot+nbmax
                        write(6,"('wwwwwsc ',i5,2f10.4)") it,sum(abs(zmel(:,it,1)))
                     enddo
                     write(6,*)'sumcheck ppovlz=',sum(abs(ppovlz(:,:)))
                     write(6,*)'sumcheck vcoud=',sum(abs(vcoud(:)))
                  endif
ccccccccccccccccccccccccccccc
                  if(.not.newaniso() ) then
                     allocate(vcoult(1:ngb,1:ngb),z1p(ngb,nctot+nbmax,ntqxx))
                     vcoult= transpose(vcoul)
                     call matm( vcoult, zmel, z1p, ngb,ngb,(nctot+nbmax)*ntqxx )
                     deallocate(vcoult) !, zmel)
                  endif  
!TIME1 "end of matm vcoult"
!TIME0

                  if(newaniso()) then
c                    write(*,'(a,5I10)')'kino: ntqxx,nctot+nbmax,ngb=',ntqxx,nctot+nbmax,ngb
                     allocate(zmeltt(nctot+nbmax,ntqxx,ngb))
!$OMP parallel  private(vc) 
!$OMP do 
                     do itp= 1,ntqxx
                        do it = 1,nctot+nbmax
                           do ivc=1,ngb
                              zmeltt(it,itp,ivc) =  sum( zmel(:,it,itp)* ppovlz(:,ivc) )
                           enddo
                        enddo
                     enddo
!$OMP end do 
!$OMP do 
                     do 992 itpp= 1,ntqxx
                        do 993 itp = 1,ntqxx
                           if(diagonly.and.(itpp/=itp)) cycle
!! sep2013t a test:c           if(itpp>ntqxxd .and.itp/=itpp) cycle
                           do 994 it  = 1,nctot+nbmax
                              w3p(it,itp,itpp) = 0d0
                              do ivc=1,ngb
                                 if(ivc==1.and.kx==iqini) then
                                    vc= wklm(1)* fpi*sqrt(fpi) /wk(kx)
c     print *,'wklm(1) vc=',wklm(1),vc
                                 else
                                    vc= vcoud(ivc)
                                 endif
c     zmelt1 =  sum( zmel(:,it,itp)  *ppovlz(:,ivc) )
c     zmelt2 =  sum( zmel(:,it,itpp) *ppovlz(:,ivc) )
                                 w3p(it,itp,itpp) = w3p(it,itp,itpp)
     &                                + vc * zmeltt(it,itp,ivc)*dconjg(zmeltt(it,itpp,ivc))
                              enddo
 994                       continue
 993                    continue
 992                 continue
!$OMP end do 
!$OMP end parallel 
                     deallocate(zmeltt)
                  else
!$OMP parallel do 
                     do itpp= 1,ntqxx
                        do itp = 1,ntqxx
                           if(diagonly.and.(itpp/=itp)) cycle
c sep2013t a test:c  if(itpp>ntqxxd .and.itp/=itpp) cycle
                           do it  = 1,nctot+nbmax
                              w3p(it,itp,itpp) =dcmplx(
     &                             sum ( dreal(z1p(:,it,itp))*rmelt(:,it,itpp)
     &                             +   dimag(z1p(:,it,itp))*cmelt(:,it,itpp) ) ,
     &                             sum ( dimag(z1p(:,it,itp))*rmelt(:,it,itpp) 
     &                             -   dreal(z1p(:,it,itp))*cmelt(:,it,itpp) ) )
                           enddo
                        enddo
                     enddo
!$OMP end parallel do 
                     deallocate(z1p)
                  endif  
                  deallocate(zmel)
!!-- Write the Spectrum function for exchange May. 2001
                  if(ifexsp/=0) then
                     do it  = 1, nctot+nbmax
                        do itp = 1,ntqxx
                           write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")
     &                          ip,itp,it, qbz_kr, ekc(it), -wtt*dreal(w3p(it,itp,itp))
                        enddo
                     enddo
                  endif
!TIME1 "end of write ifsexsp"
!TIME0 
                  if(debug) then
                     do  it  = 1,nctot+nbmax
                        do  itp = 1,ntqxx
                           write(6,"(' w3p =',2i4,2d14.6)") it,itp,w3p(it,itp,itp)
                        enddo
                     enddo
                  endif
!! --- Correct weigts wfac for valence by esmr
                  do it = nctot+1, nctot+nbmax
c     gaussian
                     wfac = wfacx(-1d99, ef, ekc(it), esmr)
                     do itp=1,ntqxx
                        w3p(it,itp,1:ntqxx) = wfac * w3p(it,itp,1:ntqxx)
                     enddo
                  enddo

                  do itpp=1,ntqxx
                     do itp = 1,ntqxx !S[j=1,nbloch]  z1p(j,t,n) <B(rk,j) psi(q-rk,n) |psi(q,t')>
                        if(diagonly.and.(itpp/=itp)) cycle
c sep2013t a test:c  if(itpp>ntqxxd .and.itp/=itpp) cycle
                        zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
     &                       - wtt * sum( w3p(:,itp,itpp) )
                     enddo
                  enddo
                  deallocate( w3p,rmelt,cmelt)
                  if(.not.newaniso()) deallocate(vcoul)
!TIME1 "enddo of zsec wtt sum"
                  cycle         ! next irot

               endif            !exchange
!! End of exchange section --------------
               if(timemix) call timeshow("33333 k-cycle")
!TIME1 "end of exchange section"
!TIME0

!!----------------------------------------------------------
!!---  correlation section ---------------------------------
!!----------------------------------------------------------
!! reading general energy mesh from file 'freq_r'
               open(UNIT=3111,file='freq_r')
               read(3111,*)nwp,nw_i !nwp=nw+1
               nw=nwp-1
               allocate(freq_r(nw_i:nw)) !freq_r(0)=0d0 from feb2006 (it WAS freq_r(1)=0d0)
               do iw= nw_i,nw
                  read(3111,*) freq_r(iw)
               enddo
               close(3111)
!TIME1 "end reading freq_r 3111"
!TIME0 
!!  Check freq_r for nw_i case.
!!  NOTE: We assume freq_r(i) == -freq_r(-i) in this code. feb2006
               if(nw_i/=0) then
                  if(nw/= -nw_i)        call rx( "sxcf_fal3_scz: nw/=-nw_i")
                  if(freq_r(0)/=0d0)    call rx( "sxcf_fal3_scz: freq_r(0)/=0")
                  if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0) 
     &                 call rx( "sxcf_fal3_scz: freq_r /= -freq_r")
               endif
!!  NOTE: this program assumes freq_r(iw)=freq_r(-iw). feb2006
!!        So freq_r(iw <0) is a little redundant because of the above relation.
!TIME1 "end check freq_r"

!!  The matrix elements zmel.
cKINO              write(*,'(a,9i7)')'kino: ngb,nstate,ntqxx,size(ppovlz)=',
cKINO&             ngb,nstate,ntqxx,size(ppovlz,dim=1),size(ppovlz,dim=2)
!TIME0 
               allocate( zmel (ngb, nstate, ntqxx) )
               zmel = dcmplx (rmelt,-cmelt)
               if(newaniso()) then
!$OMP parallel 
!$OMP do 
                  do itp=1,ntqxx
                     do it=1,nstate
                        zmel(:,it,itp) =  matmul(zmel(:,it,itp),dconjg(ppovlz(:,:)))
                     enddo
                  enddo
!$OMP end do nowait 
!$OMP end parallel
               endif  
               deallocate(rmelt,cmelt)
!TIME1 "matmul zmel ppovlz"
!TIME0
               if(debug) print *,' end of zmel'
!!================================================================
!! need to check the following notes.
!!     The correlated part of the self-energy:
!!     S[n=all] S[i,j=1,nbloch]
!!     <psi(q,t) |psi(q-rk,n) B(rk,i)>
!!     < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
!!     <B(rk,j) psi(q-rk,n) |psi(q,t)>
!!     e = e(q-rk,n), w' is real, Wc = W-v
!!================================================================

!!====================================================================
!!     contribution to SEc(qt,w) from integration along the imaginary axis
!!====================================================================

!!------------------------------------------------
!!     loop over w' = (1-x)/x, frequencies in Wc(k,w')
!!     {x} are gaussian-integration points between (0,1)
!!------------------------------------------------
               allocate( zw (nblochpmx,nblochpmx), zwz0(nstate,ntqxx,ntqxx) )
c$$$  cccccccccccccccccccccccccccccccccccccccccccc
c$$$  c write test Apr2006
c$$$  c---------------------------------------
c$$$  write(6,"(' qibz_k  ngb  ',3d13.5,3i5)")  qibz_k,ngb
c$$$  allocate( zwm (nblochpmx,nblochpmx) )
c$$$  do ix= 0,10
c$$$  iix =  ix + 1 - nw_i
c$$$  nrec=(kx-2)*(nw-nw_i+1) + iix
c$$$  read(ifrcw,rec=nrec) zw
c$$$  iix = -ix + 1 - nw_i
c$$$  nrec=(kx-2)*(nw-nw_i+1) + iix
c$$$  read(ifrcw,rec=nrec) zwm
c$$$  do ipp=1,5
c$$$  write(6,"('vvvvvvvvvvvvvvvvvvvvv ipp ix=',2i3)") ipp,ix
c$$$  do i = 200,220
c$$$  if(ipp==1) j = i
c$$$  if(ipp>2) j = ipp+i
c$$$  write(6,"('zw= ', i2,2x, 2i5,2d13.5)") ix,i,j,zw (i,j)
c$$$  write(6,"('zwm=', i2,2x, 2i5,2d13.5)") ix,i,j,zwm(i,j)
c$$$  write(6,*)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  stop 'test end------'
c$$$  ccccccccccccccccccccccccccccccccccccccccccc 
               allocate(zmel1(ngb),zmel1_(ntqxx,ngb,nstate))
               do itpp= 1,ntqxx
                  do it  = 1,nstate
                     zmel1_(itpp,1:ngb,it) = zmel(1:ngb,it,itpp)
                  enddo
               enddo
               if(timemix) call timeshow("bbb33333 k-cycle")
               ix = 1 - nw_i    !at omega=0
               nrec=ix
c     nrec=(kx-iqini)*(nw-nw_i+1) + ix !2 -->iqini newaniso
c     if(bzcase()==2) nrec= (kx-1)*(nw-nw_i+1) + ix
c     if(debug) print *,' kx, nw,nw_i,ix,nrec=',kx,nw,nw_i,ix,nrec
!TIME1 "before read frcw zw"
!TIME0 
               read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v

c$$$cccccccccccccccccccccccccccccccccccccccccccccc
c$$$                  if(sum(abs(qxx(1:2)))<1d-4) then
c$$$                  do irot2 = 1,ngrp
c$$$                     kr = irkip(kx,irot2,ip) ! index for rotated kr in the FBZ
c$$$                     if(kr==0) cycle ! next irot
c$$$                     zw(1,1) = zw(1,1)*wqfac(kr)
c$$$                     exit
c$$$                  enddo
c$$$                  print *,'wwwwwwwwwww why here ?wwwwwwwwww11111111'
c$$$                  endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccc

!TIME1 "read ifrwc zw before matzwz2"
!TIME0
c     zw=(zw+dconjg(transpose(zw)))/2d0

c$$$ sep2013test
c$$$               call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxxd, nstate,ngb,   
c$$$     o              zwz0(1:nstate,1:ntqxxd,1:ntqxxd))       ! zwz0 = zmel*(W(0)-v)*zmel
c$$$!! mode 1 of matzwz2 is only for diagonal part                
c$$$               call matzwz2(2,zw(1:ngb,1:ngb),zmel(1:ngb,1:nstate,ntqxxd+1:ntqxx),ntqxx-ntqxxd, nstate,ngb,
c$$$     o            zwz0(1:nstate,ntqxxd+1:ntqxx,ntqxxd+1:ntqxx)) ! zwz = zmel*(W(0)-v)*zmel
               call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
     o              zwz0)       ! zwz0 = zmel*(W(0)-v)*zmel
!TIME1 "after matzwz2"
!TIME0
c     call matzwzs(zw(1:ngb,1:ngb), zmel, ntq, nstate,ngb,
c     o    zwz0) ! zwz0 = zmel*(W(0)-v)*zmel
c     print *,'sum check2---',sum(abs(zwz0))


ccccccccccccccccccccccccccccccccccccccccccc
!!     COH term test ----- The sum of the all states for zwz00 gives the delta function.
cccccccccccccccccccccccccccccccfaleev 22May02 we do not use cohtest
c     if(cohtest) then
c     do itp = 1,ntq
c     coh(itp,ip)  = coh(itp,ip)
c     &       + .5d0*wtt*sum(dreal(zwz0(1:nstate,itp,itp)))
c     enddo
c     deallocate(zw,zwz0,zmel)
c     cycle
c     endif
ccccccccccccccccccccccccccccccccccccccccccc     
               nx  = niw
               if(niw <1) call rx( " sxcf:niw <1")
               npm = 1
               if(nw_i/=0) npm = 2 ! TimeReversal off case.
               allocate(zwz(niw,nstate,ntqxx,ntqxx), zwzi(nstate,ntqxx,ntqxx)) !sf 22may02
               zwzi=0d0
c     if(screen) allocate(zwz00(nstate,ntq,ntq)) !sf 22may02
               if(verbose()>50) call timeshow("6 before matzwz in ix cycle")

C     ******************************* BEGIN 1 ******************************
C     if(.false.) then !This is wrong ---do loop is not paralell
C     allocate(zz3(ntq))
C     do ix = 1,nx            ! imaginary frequency w'-loop
C     read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis
C     c         read(ifrcwi,rec=((kx-2)*niw+ix)) zw4  ! Readin W-v on imag axis
C     c             do i=1,ngb; do j=1,ngb
C     c             zw(j,i)=zw4(j,i) !zw4 is complex(4) to reduce size WVR and WVI files
C     c             enddo;      enddo
C     c         zw(1:ngb,1:ngb)=zw4(1:ngb,1:ngb) !zw4 is complex(4) to reduce size WVR and WVI files
C     do itp=1,ntq
C     do it =1,nstate
C     zmel1(:)=dconjg(zmel(:,it,itp))
C     zz3 = 0d0
C     do igb2=2,ngb !sum over igb1=1,igb2-1
C     ccccccccccccccccccc
C     c               zz2=sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) ) +
C     c    &              .5d0*zmel1(igb2)*zw(igb2,igb2) !second time consuming part
C     cccccccccccccccc
C     zz2= zdotu(igb2-1,zmel1(1),1,zw(1,igb2),1) +
C     &              .5d0*zmel1(igb2)*zw(igb2,igb2) !second time consuming part
C     cccccccccccccccccc
C     do itpp = 1,ntq
C     if(iSigma_en==5.and.itpp/=itp) cycle
C     zz3(itpp)=zz3(itpp)+zz2*zmel1_(itpp,igb2,it)
C     enddo ! itpp
C     enddo !igb2
C     do itpp = 1,ntq
C     if(iSigma_en==5.and.itpp/=itp) cycle
C     zwz(ix,it,itp,itpp)=zz3(itpp)+
C     &	       .5d0*dconjg(zmel(1,it,itp))*zw(1,1)*zmel1_(itpp,1,it)
C     enddo !itpp
C     enddo !it
C     enddo !itp
C     c          endif
C     do itpp=1,ntq
C     do itp =1,ntq
C     if(iSigma_en==5.and.itpp/=itp) cycle
C     do it=1,nstate
C     zwz(ix,it,itp,itpp)=zwz(ix,it,itp,itpp)+
C     &                          dconjg(zwz(ix,it,itpp,itp)) !!This is wrong ---do loop is not paralell
C     ccccccccc here we effectively summed  over the rest: igb1=igb2,ngb
C     ccccccccc we used W(i*omg)_{IJ}=conjg[W(i*omg)_{JI}]
C     enddo !it
C     enddo !itp
C     enddo !itpp
C     enddo !ix
C     deallocate(zz3)
C     else
C     c---instead
C     allocate(CC(ngb,nstate,ntq) )
C     do ix = 1,nx     ! imaginary frequency w'-loop
C     read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis
C     call matm(zw(1:ngb,1:ngb),zmel,cc, ngb, ngb, nstate*ntq)
C     do itp  = 1,ntq
C     do itpp = 1,ntq
C     do  it  = 1,nstate
C     zwz(ix,it,itp,itpp)
C     &     = zdotc(ngb,zmel(1,it,itp),1,CC(1,it,itpp),1)
C     enddo
C     enddo
C     enddo
C     print *,'sum check1---',sum(abs(zwz(ix,:,:,:)))
C     enddo
C     deallocate(CC)
C     endif
C     ********************** END 1 ***********************************************

               if(timemix) call timeshow("ccc33333 k-cycle")
!TIME1 "before ccc333 k-cycle"
ccccccccccccccccccccccccvariant ccccccccccccccccccccccccccccccccc
c     zwz=0d0
c!KINO               write(*,'(a,7i7)')'kino: nx,ngb,nstate,ntqxx=',nx,ngb,nstate,ntqxx
!TIME0 
              allocate( zwix(nblochpmx,nblochpmx),stat=ierr )
              if (ierr.ne.0) then 
                 write(6,'(i3,a,i3,a,i5)')
     &           mpi__rank,'failed to allocate zwix',ierr,'size=',nblochpmx**2; call flush(6) 
                 call mpi_abort()
                 call rx( '')
              endif
!$OMP parallel  private(zwzix) 
               write(6,*)mpi__rank,'zwzix allocate'; call flush(6)
               allocate(zwzix(1:nstate,1:ntqxx,1:ntqxx),stat=ierr)
               if (ierr.ne.0) then 
                  write(6,'(i5,a,i5,a,i7)') mpi__rank,
     &            'failed to allocate zwzix',ierr,'size=',nstate*ntqxx**2
!$                 write(6,'(i5)') 'omp_get_thread=',omp_get_thread_num()
                  call flush(6)
                  call mpi_abort()
                  call rx( '')
               endif
!$OMP do               
               do ix=1,nx
                  nrec=ix
                  read(ifrcwi,rec=nrec) zwix(:,:)  ! direct access read Wc(0) = W(0) - v
c     read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! direct access read Wc(0) = W(0) - v
c$$$cccccccccccccccccccccccccccccccccccccccccccccc
c$$$                  if(sum(abs(qxx(1:2)))<1d-4) then
c$$$                  do irot2 = 1,ngrp
c$$$                     kr = irkip(kx,irot2,ip) ! index for rotated kr in the FBZ
c$$$                     if(kr==0) cycle ! next irot
c$$$                     zw(1,1) = zw(1,1)*wqfac(kr)
c$$$                     exit
c$$$                  enddo
c$$$                  print *,'wwwwwwwwwww why here ?wwwwwwwwww22222'
c$$$                  endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccc
c!KINO                   write(*,'(i3,a,20i7)') mpi__rank,'before matzwz2',
c!KINO     &       ngb,ngb,size(zmel,dim=1),size(zmel,dim=2),size(zmel,dim=3),nstate,ntqxx,ntqxx
c                  call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
c$$$ccccccccccccccccccccccccccccccccccccccccccc
c$$$c!! sep2013t
c$$$                  zwzix=0d0
c$$$                  call matzwz2(2, zwix(1:ngb,1:ngb,ix), zmel, ntqxxd, nstate,ngb,   
c$$$     o                 zwzix(1:nstate,1:ntqxxd,1:ntqxxd)) ! zwz = zmel*(W(0)-v)*zmel
c$$$                  if(ntqxx>ntqxxd) then
c$$$!! mode 1 of matzwz2 is only for diagonal part                
c$$$                  call matzwz2(1,zwix(1:ngb,1:ngb,ix),zmel(1:ngb,1:nstate,ntqxxd+1:ntqxx),ntqxx-ntqxxd, nstate,ngb,
c$$$     o            zwzix(1:nstate,ntqxxd+1:ntqxx,ntqxxd+1:ntqxx)) ! zwz = zmel*(W(0)-v)*zmel
c$$$                  endif
c$$$ccccccccccccccccccccccccccccccccccccccccccc
                  call matzwz2(2, zwix(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
c     o                  zwzix(1:nstate,1:ntqxx,1:ntqxx)) ! zwz = zmel*(W(0)-v)*zmel
     o                 zwzix) ! zwz = zmel*(W(0)-v)*zmel
c!KINO                   write(6,*)mpi__rank,'matzwz2 called',ix;call flush(6)
                  zwz(ix,1:nstate,1:ntqxx,1:ntqxx)= zwzix(1:nstate,1:ntqxx,1:ntqxx) 
c     call matzwzs(zw(1:ngb,1:ngb), zmel, ntq, nstate,ngb,
c     o    zwz(ix,1:nstate,1:ntq,1:ntq)) ! zwz = zmel*(W(0)-v)*zmel
ccc   print *,'sum check2---',sum(abs(zwz(ix,:,:,:)))
               enddo ! ix
!$OMP end  do 
               deallocate(zwzix) 
!$OMP end parallel 
               deallocate(zwix) 
!TIME1 "end do xxx33333 k-cycle"
               if(timemix) call timeshow("xxx33333 k-cycle")
!TIME0
cccccccccccccccccccccccccccccccccccccccccccccc
c     do  it  = 1,nstate
c     itp=1
c     ix =1
c     do  itp = 1,ntq
c     do  ix = 1,nx
c     write(6,"(' zwz on img =',3i4,2d14.6)")ix,it,itp,zwz(ix,it,itp)
c     enddo
c     enddo
c     enddo
c     stop " *** TEST END *** "
cccccccccccccccccccccccccccccccccccccccccccccc

c--------------------------------------------------------------
c     S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c     Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
c     needed to take care of the singularity in the w' integration
c     when w-e(q-rk,n) is small
c--------------------------------------------------------------
c     if(screen) then
c     zwz00 = zwz0
c     zwz0  = 0d0
c     do ix = 1,nx
c     zwz(ix,:,:,:)=zwz(ix,:,:,:) - zwz00 !sf 22May02
c     enddo
c     endif

!!-------------------------
!!     loop over w in SEc(qt,w)
!! -------------------------
c     iSigma_en==0 SE_nn'(ef)+imag integral:delta_nn'(SE_nn(e_n)-SE_nn(ef))
c     iSigma_en==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
c     iSigma_en==2 SE_nn'((e_n+e_n')/2)
c     iSigma_en==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
c     iSigma_en==4 SE_nn'(ef)
c     iSigma_en==5 delta_nn' SE_nn(e_n)
c     output file in hsfp0  should contain hermitean part of SE
c     ( hermitean of SE_nn'(e_n) means SE_n'n(e_n')^* )
c
c     we use that zwz(itp,itpp)=dconjg( zwz(itpp,itp) )
               if(ua_auto) then
                  allocate(uaa(nstate,ntqxx))
                  do itp = 1,ntqxx
                     do  it = 1,nstate
                        ratio = abs(zwz(niw,it,itp,itp)/zwz0(it,itp,itp))
                        call gen_uaa(ratio,freqx(niw),  uaa(it,itp))
                        if(verbose()>45) then
                           write(6,"(' it itp uaa=',2i4,12f8.4)")it,itp,uaa(it,itp)
                        elseif(verbose()>40.and.mod(it,10)==1.and.mod(itp,10)==1) then
                           write(6,"(' it itp uaa=', 2i4,12f8.4)")it,itp,uaa(it,itp)
                        endif
                     enddo
                  enddo
               endif
c----------------------------------------------------------
               if(iSigma_en==5) goto 2005
c----------------------------------------------------------
!TIME1 "before 1385"
!TIME0

!!  Frequency integration along the imaginary axis, s. wint.f
!!     for each e(q-rk,n) and w in SEc(qt,w)
ccc!$OMP parallel do private(omega0,we,esmrx,ratio,freqw1,ua2_,ua_)
ccc!$OMP& private(zwzs,zwzs0)
               allocate(zwzs(npm*nx))
               do 1385   itpp= 1,ntqxx
c$$$c sep2013test cccccccccccccccccccccc
c$$$                  if(itpp>ntqxxd) then
c$$$                     inii=itpp; iendi=itpp
c$$$                  else
c$$$                     inii=1; iendi=ntqxxd
c$$$                  endif   
c$$$                  do 1386    itp = inii,iendi
                  do 1386    itp = 1,ntqxx
                     if(diagonly.and.(itpp/=itp)) cycle

                     do 1387  it = 1,nstate
                        if (iSigma_en==0 .or. iSigma_en==1.or.iSigma_en==4) then
                           omega0=ef
                        elseif (iSigma_en==2) then
                           omega0=.5d0*( omega(itp)+omega(itpp) )
                        elseif (iSigma_en==3) then
                           omega0=omega(itp)
                        else
                           call rx( "sxcf: iSigma_en /=0,1,2, or 3")
                        endif
                        we =.5d0*( omega0 -ekc(it))
                        esmrx = esmr
                        if(it <= nctot) esmrx = 0d0

                        if(verbose()>50) then
                           do  ix = 1,niw
                              ratio  = abs(zwz(ix,it,itp,itpp)/zwz0(it,itp,itpp))
                              freqw1 = (1d0 - freqx(ix))/ freqx(ix)
                              ua2_(ix) = ratio !sqrt(- 1d0/freqw1*log(ratio))
c     write(6,"(' ix freqw=',i4,f10.4)") ix,freqw1
                           enddo
                           write(6,"(' sxcf_fal2sc:',3i4,' abs=',d12.4,
     &                          ' ratio= ',12f8.2)") it,itp,itpp,
     &                          abs(zwz0(it,itp,itpp)),ua2_(1:niw)
                        endif
c     if(ua_auto) then
c     ratio = .5d0 *( abs(zwz(niw,it,itp,itp  )/zwz0(it,itp,itp  ))
c     &                     +abs(zwz(niw,it,itpp,itpp)/zwz0(it,itpp,itpp)) )
c     call gen_ua(ratio,niw,freqx
c     o      ,expa_,ua_)
c     if(verbose()>45) then
c     write(6,"(' it itp itpp ua_=',3i4,12f8.4)")it,itp,itpp,ua_
c     elseif(verbose()>40.and.mod(it,10)==1.and.itp==itpp.and.mod(itp,10)==1) then
c     write(6,"(' it itp itpp ua_=',3i4,12f8.4)")it,itp,itpp,ua_
c     endif
c     endif
                        if(ua_auto) then
                           ua_ = 0.5d0*(uaa(it,itp)+uaa(it,itpp))
                           call gen_expa(niw,freqx,ua_,  expa_)
                        endif

                        do ix=1,nx
                           zwzs(ix   ) = (          zwz(ix,it,itp,itpp) ! w(iw) + w(-iw) symmetric part
     &                          + dconjg( zwz(ix,it,itpp,itp) ) )/2d0 ! w(iw) + w(-iw) symmetric part
                           if(npm==2) then
                              zwzs(ix+nx) = (          zwz(ix,it,itp,itpp) 
     &                             - dconjg( zwz(ix,it,itpp,itp) ) )/2d0/img
                           endif
                        enddo

                        if(GaussSmear()) then
                           zwzs0= (zwz0(it,itp,itpp) +dconjg(zwz0(it,itpp,itp)))/2d0 !22May2006 but not necessary.
                           zwzi(it,itp,itpp) =
c     Gaussian smearing
c     &      wintzsg_npm (npm, zwz(1,it,itp,itpp),zwz0(it,itp,itpp)
     &                          wintzsg_npm (npm, zwzs, zwzs0 !zwz0(it,itp,itpp)
     &                          ,freqx,wx,ua_,expa_,we,nx,esmrx)
                        else
                           if(npm==2) 
     &                          call rx( ' ###Not impliment wintzav for npm=2. Use Gausssmear.')
                           zwzi(it,itp,itpp) =  
     &                          wintzav ( zwz(1,it,itp,itpp),zwz0(it,itp,itpp)
     &                          ,freqx,wx,ua_,expa_,we,nx, esmrx)
                        endif
 1387                continue   !it
 1386             continue      !itp
 1385          continue         !itpp
               deallocate(zwzs)
!TIME1 "end 1385 loop"
!TIME0
c     sum over both occupied and unoccupied states and multiply by weight
               do     itpp= 1,ntqxx
c$$$c sep2013test cccccccccccccccccccccc
c$$$                  if(itpp>ntqxxd) then
c$$$                     inii=itpp; iendi=itpp
c$$$                  else
c$$$                     inii=1; iendi=ntqxxd
c$$$                  endif   
c$$$                  do itp = inii,iendi
                  do     itp = 1,ntqxx
                     if(diagonly.and.(itpp/=itp)) cycle

                     if ((iSigma_en==1 .or. iSigma_en==0) .and. itp==itpp) then
                     else
                        zsec(itp,itpp,ip)  = zsec(itp,itpp,ip) + !sf 22may02
!     new from Jan2006! I think this should be OK.  ----------------------------
!     The output of sxcf_fal2 is  <i|Re[S](e_i)|j> ------------
!     Im-axis integral gives Hermitian part of S.
!     (Be careful as for the difference between
!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!     ---because e_i is included.
!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
     &                       wtt* sum(zwzi(:,itp,itpp)) !S_{ij}(e_i)
!     old befor Jan2006
!     &        wtt*.5d0*(   sum(zwzi(:,itp,itpp))+ !S_{ij}(e_i)
!     &        dconjg( sum(zwzi(:,itpp,itp)) )   ) !S_{ji}^*(e_j)= S_{ij}(e_j)
!-----------------------------------------------------------------------------

                     endif      !this zsec contribution  is hermitean automatically
c     we use that zwz(itp,itpp)=dconjg( zwz(itpp,itp) )
!     This contribution zwzi to zsec is hermitian.
!     So this zsec is
                  enddo         !itp
               enddo            !itpp
               if(timemix) call timeshow("yyy33333 k-cycle")
!TIME1 "end yyy333 k-cycle"
!TIME0
cccccccccccccccccccccccccccccccccccccc
               if(iSigma_en==4) then
                  deallocate(zwz,zwz0,zwzi)
                  deallocate(zmel,zmel1,zmel1_,zw,freq_r) !zw4,
                  cycle ! next irot
               endif
cccccccccccccccccccccccccccccccccccccc
 2005          continue
               if (.not.allocated(zwzs)) allocate(zwzs(npm*nx))
               if (iSigma_en==1 .or. iSigma_en==0 .or.iSigma_en==5) then !..................
                  do     itp = 1,ntqxx
                     do      it = 1,nstate

                        we =.5d0*( omega(itp) -ekc(it))
                        esmrx = esmr
                        if(it <= nctot) esmrx = 0d0

                        if(verbose()>50) then
                           do  ix = 1,niw
                              ratio  = abs(zwz(ix,it,itp,itp)/zwz0(it,itp,itp))
                              freqw1 = (1d0 - freqx(ix))/ freqx(ix)
                              ua2_(ix) = ratio !sqrt(- 1d0/freqw1*log(ratio))
c     write(6,"(' ix freqw=',i4,f10.4)") ix,freqw1
                           enddo
                           write(6,"(' sxcf_fal2sc:',3i4,' abs=',d12.4,
     &                          ' ratio= ',12f8.2)") it,itp,itp,
     &                          abs(zwz0(it,itp,itp)),ua2_(1:niw)
                        endif

                        if(ua_auto) then
                           ua_= uaa(it,itp)
                           call gen_expa(niw,freqx, ua_,  expa_)
                        endif

                        do ix=1,nx
                           zwzs(ix   ) = dreal( zwz(ix,it,itp,itp)) ! w(iw) + w(-iw) symmetric part
                           if(npm==2) then
                              zwzs(ix+nx) = dimag( zwz(ix,it,itp,itp)) ! w(iw) - w(-iw)
                           endif
                        enddo

                        if(GaussSmear()) then
                           zwzi(it,itp,itp) =
c     Gaussian smearing
     &                          wintzsg_npm (npm, zwzs, zwz0(it,itp,itp)
     &                          ,freqx,wx,ua_,expa_,we,nx,esmrx)
                        else
                           if(npm==2) 
     &                          call rx( ' ###Not impliment wintzav for npm=2. Use Gausssmear.')
                           zwzi(it,itp,itp) =  
     &                          wintzav ( zwz(1,it,itp,itp),zwz0(it,itp,itp),
     &                          freqx,wx,ua_,expa_,we,nx, esmrx)
                        endif
c     .    wintz (zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx)
                     enddo      !it
                  enddo         !itp

c     sum over both occupied and unoccupied states and multiply by weight

                  do     itp = 1,ntqxx
                     zsec(itp,itp,ip)  = zsec(itp,itp,ip) + !sf 22may02
     &                    wtt*sum(zwzi(:,itp,itp)) !this  contribution to zsec is not hermitean
c     &     + dconjg( sum(zwzi(:,itp,itp)) )   )*.5d0          !sf 22May02
                  enddo         !itp
               endif            ! iSigma_en==1 .or. iSigma_en==0!...................
c     end of SEc w-loop
               if(ua_auto) deallocate(uaa)
               if(debug2) then
                  print *,' ntqxx nstate sum(zwzi)=',ntqxx,nstate,sum(zwzi)
                  print *,' ntq   nstate sum(zwz )=',ntq,  nstate,sum(zwz)
                  do itp = 1,ntq
                     write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip) !sf 22May02
                  enddo
                  call rx( 'xxxxxxxxxxxxx test end ecal xxxxxxxxxxxxxx')
               endif
               deallocate(zwz,zwz0,zwzi)
               if (iSigma_en==0) then
                  deallocate(zmel,zmel1,zmel1_,zw,freq_r,zwzs) !zw4,
                  cycle ! next irot        !no pole contribution for SE(e_f)  sf 23may02
               endif
               deallocate(zwzs)
!TIME1 "end of deallocate zmel zmel1"
!TIME0
ccccccccccccccccccccccccccccccccccccc
c     if(onlyimagaxis()) then
c     print *,' sxcf_fal2: onlyimagaxis cycle'
c     deallocate(freq_r,zmel,zw,zmel1)
c     goto 1100
c     endif
ccccccccccccccccccccccccccccccccccccc

c==============================================
c     contribution to SEc(qt,w) from the poles of G
c==============================================
               if(debug) print *,' go to poles'
c---------------------------------------
c     maximum ixs reqired.
c---------------------------------------
               ixsmx =0
               ixsmin=0
               do 3001 itp = 1,ntqxx
                  omg  = omega(itp) !sf 23 May02
                  if (omg < ef) then
                     itini= 1
                     itend= nt0p
                  else
                     itini= nt0m+1
                     itend= nstate
                  endif
                  do 3011 it=itini,itend
                     esmrx = esmr
                     if(it<=nctot) esmrx = 0d0
                     wfac = wfacx2(omg,ef, ekc(it),esmrx)
                     if(GaussSmear()) then
                        if(wfac<wfaccut) cycle ! next it
                        we = .5d0*(weavx2(omg,ef,ekc(it),esmr)-omg)
                     else
                        if(wfac==0d0) cycle ! next it
                        if(omg>=ef) we = max( .5d0*(omg-ekc(it)), 0d0) ! positive
                        if(omg< ef) we = min( .5d0*(omg-ekc(it)), 0d0) ! negative
                     endif
                     if(it<=nctot) then
                        if(wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
                     endif
                     do iwp = 1,nw
                        ixs=iwp
                        if(freq_r(iwp)>abs(we)) exit
                     enddo
c     This change is because G(omega-omg') W(omg') !may2006
c     if(ixs>ixsmx  .and. omg<=ef ) ixsmx  = ixs
c     if(ixs>ixsmin .and. omg> ef ) ixsmin = ixs
                     if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
                     if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
                     wexx  = we
                     if(ixs+1 > nw) then
                        write (*,*) ' nw_i ixsmin',nw_i, ixsmin
                        write (*,*) ' wexx, dw ',wexx,dw
                        write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
                        call rx( ' sxcf 222: |w-e| out of range')
                     endif
 3011             continue
 3001          continue         !end of SEc w and qt -loop
!TIME1 "end of 3001"
!TIME0
               if(nw_i==0) then
                  nwxi = 0
                  nwx  = max(ixsmx+1,ixsmin+1)
               else
                  nwxi = -ixsmin-1
                  nwx  =  ixsmx+1
               endif
               if (nwx > nw   ) then
                  call rx( ' sxcf_fal3_sc nwx check : |w-e| > max(w)')
               endif
               if (nwxi < nw_i) then
                  call rx( ' sxcf_fal3_sc nwxi check: |w-e| > max(w)')
               endif
               if(debug) print *,'nw, nwx=',nw,nwx
               if(verbose()>50)call timeshow("10before alagr3z iw,itp,it ")

C...  Find nt_max ------------------------------------
               nt_max=nt0p      !initial nt_max
               do 4001 itp = 1,ntqxx
                  omg     = omega(itp)
                  if (omg > ef) then
                     do  it = nt0m+1,nstate ! nt0m corresponds to efm
                        wfac = wfacx2 (ef,omg, ekc(it),esmr)
                        if( (GaussSmear().and.wfac>wfaccut)
     &                       .or.(.not.GaussSmear().and.wfac/=0d0)) then
                           if (it > nt_max) nt_max=it ! nt_max is  unocc. state
                        endif   ! that ekc(it>nt_max)-omega > 0
                     enddo      ! so it > nt_max does not contribute to omega pole integral
                  endif
 4001          continue         !end of  w and qt -loop

!TIME1 "end of 4001"
!TIME0
****************************iSigma_en =1,3  begin ************************
               if(timemix) call timeshow("444444 k-cycle")
               if (iSigma_en == 1 .or. iSigma_en == 3.or.iSigma_en==5) then
               else
                  call rx( 'sxcf_fal2_sc: iSigma_en /= 1 3 5')
               endif
               allocate( zw3(ngb,ngb,nwxi:nwx))
               do      ix = nwxi,nwx
                  nrec= ix-nw_i+1
                  read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v
c     read(ifrcw,rec=((kx-2)*nw+ix)) zw
c     nrec=(kx-iqini)*(nw+1-nw_i) + ix-nw_i+1
c     if(bzcase()==2) nrec=(kx-1)*(nw+1-nw_i) + ix-nw_i+1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c!!!!!!!!!!!!!! for W0test !!!!!!!!!!!!!!!!
c                  nrec=1 - nw_i !at omega=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cccccccccccccccccccccccccccccccccccccccccccccc
c$$$                  if(sum(abs(qxx(1:2)))<1d-4) then
c$$$                  do irot2 = 1,ngrp
c$$$                     kr = irkip(kx,irot2,ip) ! index for rotated kr in the FBZ
c$$$                     if(kr==0) cycle ! next irot
c$$$                     zw(1,1) = zw(1,1)*wqfac(kr)
c$$$                     exit
c$$$                  enddo
c$$$                  print *,'wwwwwwwwwww why here ?wwwwwwwwww333333333333'
c$$$                  endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccc
                  do igb2=1,ngb
                     do igb1=1,ngb
                        zw3(igb1,igb2,ix)=.5d0*(zw(igb1,igb2)+dconjg(zw(igb2,igb1)) )
                     enddo      !making potential hermitean, because
                  enddo         !we need only hermitean part of Sigma_nn'
               enddo
               deallocate(zw)   !,zw4
               if(timemix) call timeshow("55555 k-cycle")
               if(test_symmetric_W().and.npm==2) then
                  if(onceww(4)) print *,' test_symmetric_W()=',test_symmetric_W(),nwxi,nwx
                  allocate(zw3x(ngb,ngb))
                  do ix= 1,min(abs(nwxi),nwx)
                     zw3x = 0.5d0* (zw3(:,:,ix) + zw3(:,:,-ix))
                     zw3(:,:, ix)=zw3x
                     zw3(:,:,-ix)=zw3x
                  enddo
                  deallocate(zw3x)
               endif
!TIME1 "before 2001"
!KINO               write(*,'(a,5i6)')'kino: ntqxx,itini,itend,ngb=',ntqxx,itini,itend,ngb
!TIME0 
!!     May2006 simplified version ===================================================
               allocate(zwz44(3,ntqxx),zwz4(ntqxx,3))
               do 2001 itp = 1,ntqxx ! loop over states (q-k,n)
                  omg = omega(itp)
                  if (omg >= ef) then
                     itini= nt0m+1
                     itend= nt_max
                     iii=  1
                  else
                     itini= 1
                     itend= nt0p
                     iii= -1
                  endif

                  do 2011 it = itini,itend ! nt0p corresponds to efp
                     esmrx = esmr
                     if(it<=nctot) esmrx = 0d0
                     wfac = wfacx2(omg,ef, ekc(it),esmrx)
                     if(GaussSmear()) then
                        if(wfac<wfaccut) cycle ! next it
c     we = .5d0* abs( weavx2(omg,ef, ekc(it),esmr)- omg )
                        we = .5d0* abs( omg-weavx2(omg,ef, ekc(it),esmr) )
                        if(it<=nctot) then !faleev
                           if(wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
                        endif
                     else
                        if( wfac==0d0) cycle ! next it
                        if( omg >= ef) we = 0.5d0* abs( max(omg-ekc(it), 0d0) )
                        if( omg <  ef) we = 0.5d0* abs( min(omg-ekc(it), 0d0) )
                        if( it<=nctot) then !faleev
                           if(wfac/=0) call rx( "sxcf:  it<=nctot.and.wfac/=0")
                        endif
                     endif
                     if(debug)write(6,"( ' xxx1',10d13.6)") omg,ef, ekc(it),wfac
                     wfac= iii* wfac*wtt
                                !ixs  = idint (we/dw) + 1

                     do iwp = 1,nw !sf 23May02
                        ixs=iwp
                        if(freq_r(iwp)>we) exit
                     enddo      !sf 23May02
                     if(nw_i==0) then
                        if(ixs+1>nwx) call rx( ' sxcf: ixs+1>nwx xxx2')
                     else
c     write(6,*) "pppp ixs nwxi=",ixs,nwxi,freq_r(ixs-1),we,freq_r(ixs)
c     if(ixs+1>abs(nwxi)) stop ' sxcf: ixs+1>nwx yyy2'
                        if(omg >=ef .and. ixs+1> nwx ) then
                           write(6,*)'ixs+1 nwx=',ixs+1,nwx
                           call rx( ' sxcf: ixs+1>nwx yyy2a')
                        endif
                        if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
                           write(6,*)'ixs+1 nwxi=',ixs+1,nwxi
                           call rx( ' sxcf: ixs-1<nwi yyy2b')
                        endif
                     endif
                     iir = 1
                     if(omg < ef .and. nw_i/=0) iir = -1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write(6,"('qqqq iir ixs we freq=',2i4,d13.6,3d13.6)") iir,ixs,we,freq_r(ixs-1:ixs+1)
c     if(we<freq_r(ixs-1).or. freq_r(ixs)<we) stop 'qqqq xxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                     zmel1(:)=dconjg(zmel(:,it,itp))
                     if (iSigma_en == 1.or.iSigma_en==5) then
                        zwz3=(0d0,0d0)
!$OMP parallel do private(ix,zz2)
                        do ix0=1,3
                           ix=ixs+ix0-2
                           do igb2=2,ngb !**** most time consuming part for iSigma_en=1 ******
                              zz2=sum(zmel1(1:igb2-1)*zw3(1:igb2-1,igb2,iir*ix)  ) +
     &                             .5d0* zmel1(igb2)*zw3(igb2,igb2,iir*ix)
                              zwz3(ix0)=zwz3(ix0)+zz2*zmel(igb2,it,itp)
                           enddo !igb2
                           zwz3(ix0)=2d0*dreal(zwz3(ix0))+
     &                          zmel1(1)*zw3(1,1, iir*ix)*zmel(1,it,itp)
                        enddo   !ix
!$OMP end parallel do 
                        if(npm==1) then
                           zsec(itp,itp,ip) = zsec(itp,itp,ip) 
     .                          + wfac*alagr3z2(we,freq_r(ixs-1),zwz3,itp,itp)
                        else
                           zsec(itp,itp,ip) = zsec(itp,itp,ip) 
     .                          + wfac*alagr3z(we,freq_r(ixs-1),zwz3)
                        endif
c     this contribution to zsec_nn is real (hermitean)
                     elseif(iSigma_en == 3) then
c                        allocate(zwz4(ntqxx,3)) !itpp
                        zwz4=(0d0,0d0)
CYY!$OMP parallel do private(ix,zz2)
                        do ix0=1,3
                           ix=ixs+ix0-2
                           do igb2=1,ngb !**** most time consuming part for iSigma_en=3 ******
                             zz2=sum(zmel1(1:ngb)*zw3(1:ngb,igb2, iir*ix)  )
cc sep2013t cccccccccccccccccccccccccccccccccccccccc
c                              do itpp=1,ntqxx
c                                 zwz4(ix0,itpp)=zwz4(ix0,itpp)+ zz2*zmel1_(itpp,igb2,it)
c                              enddo !itpp
ccccccccccccc
c$$$                             zwz4(1:ntqxxd,ix0)=zwz4(1:ntqxxd,ix0)+ zz2*zmel1_(1:ntqxxd,igb2,it)
c$$$c                            call zaxpy(ntqxxd,zz2,zmel1_(1,igb2,it),1,zwz4(1,ix0),1)
c$$$                             if(itp>ntqxxd) zwz4(itp,ix0)=zwz4(itp,ix0)+ zz2*zmel1_(itp,igb2,it)
cccccccccccccccccccccccccccccccccccccccccccccccccc
c                            zwz4(1:ntqxx,ix0)=zwz4(1:ntqxx,ix0)+ zz2*zmel1_(1:ntqxx,igb2,it)
                             call zaxpy(ntqxx,zz2,zmel1_(1,igb2,it),1,zwz4(1,ix0),1)
                           enddo !igb2
                        enddo   !ix
c                        allocate(zwz44(3,ntqxx))
                        zwz44 = transpose(zwz4)
c                        deallocate(zwz4)
c$$$c sep2013t
c$$$                        if(itp>ntqxxd) then
c$$$                           inii=itp; iendi=itp
c$$$                        else
c$$$                           inii=1; iendi=ntqxxd
c$$$                        endif
CYY!$OMP end parallel do 
c$$$                        do 2021 itpp=inii,iendi
                        do 2021 itpp=1,ntqxx
                           if(diagonly.and.(itpp/=itp)) cycle
                           if(npm==1) then
                              zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
     .                             + wfac*alagr3z2(we,freq_r(ixs-1),zwz44(1,itpp),itp,itpp)
                           else
                              zsec(itp,itpp,ip) = zsec(itp,itpp,ip) 
     .                             + wfac*alagr3z(we,freq_r(ixs-1),zwz44(1,itpp))
                           endif
 2021                   continue !itpp
c                        deallocate(zwz44)
                     endif      ! inner iSigma_en=1 or 3
c     this contribution to zsec_nn' is not hermitean because W(e_n)
c     and must be made hermitean when zsec will be written on disc
 2011             continue
 2001          continue         !itp
!TIME1 "end 2001 loop"
!TIME0 
               deallocate(zw3, zmel, zwz44,zwz4)
               if(verbose()>50)call timeshow("11after alagr3z iw,itp,it cycles")
               if(debug) print *,' end of do 2001 or 2002 '
 1050          continue
               if(debug) then
                  do itp = 1,ntq
                     write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip)
                  enddo
               endif
               deallocate(zmel1,zmel1_)
               deallocate(freq_r)
 1000       continue ! end do irot

c            print *, 'end kx=',kx,' out of iqend=',iqend,' exchange=',exchange
            if(newaniso()) ifvcoud =iclose('Vcoud.'//charnum5(kx))
            if(.not.exchange) then
               ifrcw  = iclose('WVR.'//charnum5(kx))
               ifrcwi = iclose('WVI.'//charnum5(kx))
            endif   
!TIME1 "end of 1100 loop"
 1100    continue               ! end of kx-loop
         if(newaniso()) ifvcoud =iclose('Vcoud')
         if(irot==1) write(6,"('  sum(abs(zsec))=',d23.15)") sum(abs(zsec))
ccccccccccccccccccccccccccccc
c     print *," end kount ip irot=", irot, ip, sum(kount(:,ip))
c     print *," sumkount 1=",irot,ip,sum(kount(:,:))
cccccccccccccccccccccccccccc
         if (allocated(vcoul))deallocate(vcoul)
 1001 continue ! end do ip
      if (allocated(expikt))deallocate(expikt)
!TIME1 "end of sxcf_fal2_scz"
      end subroutine sxcf_fal3_scz




