      subroutine sxcf_fal3_scz(kount,qip,itq,ntq,ef,esmr,tiat,miat,
     i     nsp,isp,             
     i     qbas,ginv,
     i     qibz,qbz,wk,nstbz,nstar,irkip,nrkip, 
     i     iclass,mdim,nlnmv,nlnmc,
     i     icore, ncore,imdim,
     i     ppb,                 
     i     freqx,wx,dw,         
     i     ecore,
     d     nlmto,nqibz,nqbz,nctot,
     d     nl,nnc,nclass,natom,
     d     nlnmx,mdimx,nbloch,ngrp,niw,nq, 
     &     nblochpmx ,ngpmx,ngcmx, 
     &     wgt0,nq0i,q0i,symgg, alat, shtvg,nband, ifvcfpout,
     &     exchange,screen,cohtest, ifexsp,
     i     iSigma_en,           
     i     nbmx,ebmx,           
     i     pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, 
     i     wklm,lxklm,          
     i     invg,il,in,im,nlnm,nn_,lx,nx_,nxx_,cgr,ppbrd, 
     i    diagonly, eftrue,
     o     zsec,coh,nbandmx)
      use m_readqg
      use m_readeigen,only: readeval,readcphi,readgeig
      use keyvalue,only: getkeyvalue
      use m_mpi, only: mpi__rank 
      implicit none
!> \brief
!! Calcualte full simga_ij(e_i)= <i|Re[Sigma](e_i)|j> 
!! ---------------------
!! \param exchange 
!!   - T : Calculate the exchange self-energy
!!   - F : Calculate correlated part of the self-energy
!! \param zsec
!!   - S_ij= <i|Re[S](e_i)|j>
!!   - Note that S_ij itself is not Hermite becasue it includes e_i.
!!     i and j are band indexes
!! \param coh dummy 
!! \param screen dummy 
!!
!!
!! \remark
!! \verbatim
!! Jan2013: eftrue is added.
!!   ef=eftrue(true fermi energy) for valence exchange and correlation mode.
!!   but ef is not the true fermi energy for core-exchange mode.
!!
!! Jan2006
!!     "zsec from im-axis integral part"  had been symmetrized as
!!     &        wtt*.5d0*(   sum(zwzi(:,itp,itpp))+ !S_{ij}(e_i)
!!     &        dconjg( sum(zwzi(:,itpp,itp)) )   ) !S_{ji}^*(e_j)= S_{ij}(e_j)
!!     However, I now do it just the 1st term.
!!     &        wtt* sum(zwzi(:,itp,itpp))   !S_{ij}(e_i)
!!     This is OK because the symmetrization is in hqpe.sc.F
!!     Now zsec given in this routine is simply written as <i|Re[S](e_i)|j>.
!!     ( In the version until Jan2006 (fpgw032f8), only the im-axis part was symmetrized.
!!     But it was not necessary from the begining because it was done in hqpe.sc.F
!!     
!!     (Be careful as for the difference between
!!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!!     ---because e_i is included.
!!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
!! \endverbatim

c original document for sce.f (correlation case) is kept. This is by F.Aryasetiawan.
c Current GW code is very away from his original, but kept here as a reference.
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
      integer:: dummy4doxygen

! input variables
      logical, intent(in) :: exchange,screen,cohtest,diagonly
      integer, intent(in) :: ntq,natom,nqbz,nqibz,ngrp,nq,niw
      integer, intent(in) :: nband,nlmto,nq0i,nctot,isp,nsp,nlnmx,mdim(*)
      integer, intent(in) :: ifvcfpout,nbloch,nclass,nl,nnc, nblochpmx
      integer, intent(in) :: mdimx,itq(ntq),miat(natom,ngrp),nstar(nqibz)
      integer, intent(in) :: irkip(nqibz,ngrp,nq),nrkip(nqibz,ngrp,nq)
      integer, intent(in) :: kount(nqibz,nq),ngpmx,ngcmx,ifexsp,iSigma_en
      integer, intent(in) :: nbmx(2),nlnmv(*),nlnmc(*),iclass(*),icore(*)
      integer, intent(in) :: ncore(*),imdim(*),nstbz(nqbz),nnmx,nomx,nkpo
      integer, intent(in) :: nnr(nkpo),nor(nkpo),lxklm,invg(ngrp)
      integer, intent(in) :: il(*),in(*),im(*),nlnm(*),nn_,lx(*),nx_(*),nxx_
      real(8), intent(in) :: wgt0(nq0i,ngrp),symgg(3,3,ngrp)
      real(8), intent(in) :: q0i(1:3,1:nq0i),shtvg(3,ngrp),alat,ecore(nctot)
      real(8), intent(in) :: coh(ntq,nq),qbas(3*3),ginv(3*3)
      real(8), intent(in) :: tiat(3,natom,ngrp),wk(nqbz),qibz(3,nqibz)
      real(8), intent(in) :: qbz(3,nqbz),freqx(niw),wx(niw),ef,esmr,dw
      real(8), intent(in) :: ebmx(2),qrr(3,nkpo),wklm((lxklm+1)**2)
      real(8), intent(in) :: cgr(*),ppbrd(*),qip(3,*),eftrue
      complex(8), intent(in) :: pomatr(nnmx,nomx,nkpo)
      real(8), intent(in) :: ppb(*)

! output variables
      integer, intent(out) ::nbandmx(nq)
      complex(8), intent(out) :: zsec(ntq,ntq,nq)          

! local variables
      integer :: ifrcw,ifrcwi
      logical :: initp=.true.
c$$$      logical :: ua_auto !fixed to be .false.

      integer :: ip, it, itp, i, ix, kx, ngc, ngb, irot, kr
      integer :: iatomp(natom), nt0p, nt0m, nbmax, ntqxx, nstate
      integer :: nt,nw,nx,ixsmx,ixs,iw,ivc,ifvcoud,ngb0
      integer :: nprecx,mrecl,ifwd,nrot,nwp,nw_i,ierr 
      integer :: nbcut,nstatetot,iqini,iqend
      integer :: invr,nbmxe,ia,nn,ntp0,no,itpp,nrec,npm,ixsmin,itini,itend
      integer :: iwp,nwxi,nwx,iir,nt_max, igb1,igb2,ix0,iii

      real(8) :: tpi, ekc(nctot+nband), qbasinv(3,3), det, q(3),ua_
      real(8) :: expa_(niw), qxx(3), tr(3,natom), symope(3,3),shtv(3)
      real(8) :: qk(3), ekq(nband),efp,efm,wtt,wfac,we,esmrx
      real(8) :: qvv(3),pi,fpi,eq(nband),omega(ntq),quu(3),freqw,ratio
      real(8) :: qibz_k(3),qbz_kr(3),ddw,vc,omega0,omg

      complex(8) :: cphiq(nlmto,nband), cphikq(nlmto,nband)
      complex(8) :: zwzs0,zz2,zwz3(3)


! local arrays
      real(8),allocatable :: freq_r(:),rmelt(:,:,:),cmelt(:,:,:),vcoud(:)
      real(8),allocatable :: drealzzzmel(:,:,:), dimagzzzmel(:,:,:),uaa(:,:)
      complex(8),allocatable :: vcoul(:,:),w3p(:,:,:)   
      complex(8),allocatable :: zmel(:,:,:),zzzmel(:,:,:),zw (:,:)
      complex(8),allocatable :: zwz(:,:,:,:), zwz0(:,:,:),zwzi(:,:,:)
      complex(8),allocatable :: zwix(:,:),zwzix(:,:,:),expikt(:),zmel1(:)
      complex(8), allocatable :: zmel1_(:,:,:), zw3(:,:,:),zw3x(:,:)
      complex(8), allocatable :: zwz4(:,:),zwz44(:,:),pomat(:,:), zwzs(:)
      complex(8),allocatable :: ppovl(:,:),ppovlz(:,:),zcousq(:,:)
      complex(8),allocatable :: zmeltt(:,:,:),z1r(:,:),z2r(:,:),w3pi(:,:)

      real(8), parameter :: wfaccut=1d-10
      complex(8), parameter :: img=(0d0,1d0)

! external function
      logical :: smbasis
c      logical :: test_symmetric_W
c      logical :: GaussSmear !fixed to be T
c      logical :: newaniso !fixed to be T
      character(5) :: charnum5
c      integer(4) :: bzcase !fixed to be 1
      integer(4) :: iopen,iclose
      integer :: invrot
      complex(8) :: wintzav, wintzsg_npm
      integer :: nocc
      real(8) :: wfacx
      real(8) :: wfacx2
      real(8) :: weavx2
      complex(8) :: alagr3z
      complex(8) :: alagr3z2

      integer:: ndummy1,ndummy2,nlmtobnd,nt0
      real(8):: wexx
c      complex(8),allocatable :: z1p(:,:,:),vcoult(:,:)
      logical :: debug, debugp,debug2=.false.
      logical :: gass !external
c      real(8):: wgtq0p
      integer(4)::verbose
      real(8):: ua2_(niw),freqw1
      integer(4) ::nbcutc
      real(8):: q_r(3)
      logical:: timemix =.false., oncew
      logical :: onceww
      logical ::eibz4sig

!TIME0 
!TIME0 
!TIME0 
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
      if(verbose()>=90) debug=.true.
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      nbcutc=nctot+nbcut
      if(debug) then
         write(6,*) ' sxcf_fal3_scz: Gausssmear=T'
c         if(GaussSmear()) then
c            write(6,*)' new version: gaussian smearing!!!'
c            write(6,*)' wintzsg wfacx weavx are related parts...'
c            gass=GaussSmear()
c            write(6,*)' sxcf_fal3_scz: new version: gaussian smearing' ,gass
c         endif
      endif
      tpi         = 8d0*datan(1d0)
      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
      call dinv33(qbas,0,qbasinv,det)
      if(debug) write(6,*) ' sxcf: 1'
      allocate(expikt(natom))

!!== ip loop to spedify external q ==
!TIME1 "before 1001"
      do 1001 ip = 1,nq         
         if(sum(irkip(:,:,ip))==0) cycle ! next ip
         write (6,*) ip,'  out of ',nq,'  k-points(extrnal q) '
!TIME0 
         q(1:3)= qip(1:3,ip)
         call readeval(q,isp,eq)
         do i  = 1,ntq
            omega(i) = eq(itq(i)) 
         enddo

c     if(nctot==0) then !this is to override -Mbounds check
c     deallocate(ecore)
c     allocate(ecore(1)) !this is dummry
c     endif
c     if(debug) write(6,*) "goto sxcf_fal3_scz"
c     write(6,*) ' sxcf_fal3_scz:'
c     iq          = idxk (q,qbz,nqbz) ! index for q

c     direct access read e,b and hb at q
c     iqisp =isp + nsp*(iq-1)
c     call rwdd1   (ifcphi, iqisp, 2*nlmtobnd, cphikq)
         call  readcphi(q, nlmto,isp, quu, cphikq)
         if(debug) write(6,*) ' sxcf: 2'

!!  b and hb of states {t} at q
         do      it = 1,ntq
            itp        = itq(it)
            cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
         enddo

!! we only consider bzcase()==1
         if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
         if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')

!! gauss_img : interpolation gaussion for W(i \omega).
         call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
c         if(ua_<=0d0) then
c            ua_auto =.true.
c            write(6,"(' ua_auto=T')")
c         else
c            ua_auto =.false.
         do ix = 1,niw
            freqw     = (1d0 - freqx(ix))/ freqx(ix)
            expa_(ix) = exp(-(ua_*freqw)**2)
         enddo
c         endif

!!  Loop over irreducible k-points 
         iqini=1
         iqend=nqibz            !no sum for offset-Gamma points.
         if(.not.exchange) then
            ifwd = iopen('WV.d',1,-1,0)
            read (ifwd,*) nprecx,mrecl
            ifwd = iclose('WV.d')
         endif 
!TIME1 "before 1100"
         do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
            if(sum(irkip(kx,:,ip))==0) cycle ! next kx
!TIME0 
            write(6,*) ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
            if( kx <= nqibz ) then
               qibz_k= qibz(:,kx)
            else
               qibz_k= 0d0
            endif
            if(verbose()>=40)  write(6,*) ' sxcf_fal2_scz: loop 1100 kx=',kx
            if(timemix) call timeshow("11111 k-cycle")
            call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
            ngb = nbloch + ngc
            if(debug) write(6,*) ' sxcf: ngb=',ngb,nbloch

!! ===Readin diagonalized Coulomb interaction===
!! note sep102012takao
!!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
!!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
!!  Vould file is written in hvccfp0.m.F.
!! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
!! These can be also separeted into WVR.ID and WVI.ID files.
c$$$            if(newaniso()) then
            if(kx<=nqibz) qxx=qibz_k
            if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
            ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
            do
               read(ifvcoud) ngb0
               read(ifvcoud) qvv
               if(allocated(vcoud)) deallocate(vcoud)
               allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
               read(ifvcoud) vcoud
               read(ifvcoud) zcousq
               if(sum(abs(qvv-qxx))<1d-6) goto 1133
            enddo
            if(sum(abs(qvv-qxx))>1d-6) then
               write(6,*)'qvv =',qvv
               write(6,*)'qxx=',qxx,kx
               call rx( 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent')
            endif
 1133       continue
            if( ngb0/=ngb ) then !sanity check
               write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
               call rx( 'hsfp0.m.f:ngb0/=ngb')
            endif
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!! ppovlz= O Zcousq
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
            if(allocated(ppovlz)) deallocate(ppovlz)
            allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
            call readppovl0(qibz_k,ngc,ppovl)
            ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
            ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
            deallocate(zcousq,ppovl)
c$$$            endif


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
c               if(GaussSmear()) ddw= 10d0
               ddw= 10d0
               efp= ef+ddw*esmr
               efm= ef-ddw*esmr
               nt0p = nocc (ekc,efp,.true.,nstatetot)
               nt0m = nocc (ekc,efm,.true.,nstatetot)

!! nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!! nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W
               if(exchange) then
                  nbmax = nt0p-nctot
               else
                  nbmax = nband
                  nbmxe = nocc (ekc,ebmx(2),.true.,nstatetot)-nctot
                  nbmax  = min(nband,nbmx(2),nbmxe)
                  if(initp) then
                     write(6,*)' nbmax=',nbmax
                     initp=.false.
                  endif
               endif
               if(debug) write(6,*) ' sxcf: 1100c'
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
               if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
               nstate = nctot + nbmax ! = nstate for the case of correlation
!TIME1 "before expikt"
!TIME0

!! ====================================================
!!     matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
!!     including the phase factor exp(ik.T)
!!     B(rot*k,r) = B(k,invrot*r)
!! ====================================================
!! phase factors expikt(ia) is for exp(ik.T(R))
               do ia = 1,natom
                  expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) 
               end do
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
               if(debug)  write(6,*) ' sxcf_fal2sc: goto psi2bc1'
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
               if(debug) write(6,*) ' sxcf_fal1: goto drvmelp2 xxx111'
               allocate(drealzzzmel(nbloch,nt,ntqxx),dimagzzzmel(nbloch,nt,ntqxx))
               drealzzzmel=dreal(zzzmel)
               dimagzzzmel=dimag(zzzmel)
               deallocate(zzzmel)
               allocate( rmelt(ngb, nctot+nbmax, ntqxx), ! nstate= nctot+nband
     &              cmelt(ngb, nctot+nbmax, ntqxx))
               call drvmelp2( q,             ntqxx, ! q in FBZ
     i              q-qbz_kr,  nbmax, ! q-rk
     i              qibz_k,  ! k in IBZ for mixed product basis. rk = symope(qibz_k) 
     i              isp,ginv,
     i              ngc,ngcmx, ngpmx,nband,itq, 
     i              symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i              drealzzzmel, dimagzzzmel, nbloch, nt,nctot,
     o              rmelt,cmelt)
               if(debug) write(6,*) ' sxcf_fal1: end of drvmelp2'
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
               if(debug) write(6,*) ' sxcf: goto wtt'
               if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)")sum(rmelt),sum(cmelt)

c$$$               if(bzcase()==2)then
c$$$                  if(kx<=nqibz) then
c$$$                     wtt = wk(kr)
c$$$                     if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
c$$$                  elseif(kx>nqibz) then !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
c$$$                     wtt= wgt0(kx-nqibz,irot)
c$$$                  endif
c$$$               else
                  if(kx<= nqibz) then
c     wtx = 1d0
                     wtt = wk(kr)
                  else
c     wtx = wgt0(kx-nqibz,irot)
                     wtt = wk(1)*wgt0(kx-nqibz,irot)
                     if(abs(wk(1)-1d0/dble(nqbz))>1d-10) call rx( 'sxcf:wk(1) inconsistent')
                  endif
c$$$               endif
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
c                  if(.not.newaniso() ) then
c                     allocate(vcoult(1:ngb,1:ngb),z1p(ngb,nctot+nbmax,ntqxx))
c                     vcoult= transpose(vcoul)
c                     call matm( vcoult, zmel, z1p, ngb,ngb,(nctot+nbmax)*ntqxx )
c                     deallocate(vcoult) !, zmel)
c                  endif  
!TIME1 "end of matm vcoult"
!TIME0

c$$$                  if(newaniso()) then
c                    write(*,'(a,5I10)')'kino: ntqxx,nctot+nbmax,ngb=',ntqxx,nctot+nbmax,ngb
                     allocate(zmeltt(nctot+nbmax,ntqxx,ngb))
c#ifdef USE_GEMM_FOR_SUM
                     if(verbose()>39) then
                        write(*,*)'info: USE GEMM FOR SUM (zmeltt=zmel*ppovlz) in sxcf_fal2.sc.F'
                        write(*,*)'zgemmsize',ntqxx*(nctot+nbmax),ngb,ngb
                        write(*,*)'size ,zmel',size(zmel,dim=1),size(zmel,dim=2),size(zmel,dim=3)
                        write(*,*)'size ,ppovlz',size(ppovlz,dim=1),size(ppovlz,dim=2)
                        write(*,*)'size ,zmeltt',size(zmeltt,dim=1),size(zmeltt,dim=2),size(zmeltt,dim=3)
                     endif
		     call flush(6)	
                     call zgemm('T','N',ntqxx*(nctot+nbmax),ngb,ngb,(1d0,0d0),
     .                 zmel,ngb,ppovlz,ngb,(0d0,0d0),zmeltt,ntqxx*(nctot+nbmax) )
                     if(verbose()>39) write(*,*)'phase1'
                     vc = vcoud(1) ! save vcoud(1)
                     if (kx == iqini) vcoud(1) = wklm(1)* fpi*sqrt(fpi) /wk(kx)
                     allocate(z1r(ntqxx,ngb),z2r(ntqxx,ngb),w3pi(ntqxx,ntqxx))
                     do  it = 1, nctot+nbmax
                        do  ivc = 1, ngb
                           do  itp = 1, ntqxx
                              z1r(itp,ivc) = zmeltt(it,itp,ivc) * vcoud(ivc)
                              z2r(itp,ivc) = zmeltt(it,itp,ivc)
                           enddo ! ivc
                        enddo   ! it
                        call zgemm('N','C',ntqxx,ntqxx,ngb,(1d0,0d0),z1r,ntqxx,
     .                    z2r,ntqxx,(0d0,0d0),w3pi,ntqxx)
                        if(verbose()>39) write(*,*)'phase2'
C              call zprm('w3pi',w3p,ntqxx,ntqxx,ntqxx)
C             Faster, but harder to parallelize
!                        call zqsmpy(11,'N','C',ntqxx,ngb,z1r,ntqxx,z2r,ntqxx,
!     .                       (0d0,0d0),w3pi,ntqxx)
C             call zprm('w3pi',w3p,ntqxx,ntqxx,ntqxx)
                        do  itp = 1, ntqxx
                           do itpp = 1, ntqxx
                              w3p(it,itp,itpp) = w3pi(itp,itpp)
                           enddo
                        enddo
                        if(verbose()>39)write(*,*)'phase3'
                     enddo
                     vcoud(1) = vc
                     deallocate(z1r,z2r,w3pi)
                     if(verbose()>=30) call cputid2(' complete w3p',0)

c$$$#else
c$$$!kino 2014-08-13  !$OMP parallel  private(vc) 
c$$$!kino 2014-08-13  !$OMP do 
c$$$                     do itp= 1,ntqxx
c$$$                        do it = 1,nctot+nbmax
c$$$                           do ivc=1,ngb
c$$$                              zmeltt(it,itp,ivc) =  sum( zmel(:,it,itp)* ppovlz(:,ivc) )
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$!kino 2014-08-13  !$OMP end do 
c$$$!kino 2014-08-13  !$OMP do 
c$$$                     do 992 itpp= 1,ntqxx
c$$$                        do 993 itp = 1,ntqxx
c$$$                           if(diagonly.and.(itpp/=itp)) cycle
c$$$!! sep2013t a test:c           if(itpp>ntqxxd .and.itp/=itpp) cycle
c$$$                           do 994 it  = 1,nctot+nbmax
c$$$                              w3p(it,itp,itpp) = 0d0
c$$$                              do ivc=1,ngb
c$$$                                 if(ivc==1.and.kx==iqini) then
c$$$                                    vc= wklm(1)* fpi*sqrt(fpi) /wk(kx)
c$$$c     write(6,*)'wklm(1) vc=',wklm(1),vc
c$$$                                 else
c$$$                                    vc= vcoud(ivc)
c$$$                                 endif
c$$$c     zmelt1 =  sum( zmel(:,it,itp)  *ppovlz(:,ivc) )
c$$$c     zmelt2 =  sum( zmel(:,it,itpp) *ppovlz(:,ivc) )
c$$$                                 w3p(it,itp,itpp) = w3p(it,itp,itpp)
c$$$     &                                + vc * zmeltt(it,itp,ivc)*dconjg(zmeltt(it,itpp,ivc))
c$$$                              enddo
c$$$ 994                       continue
c$$$ 993                    continue
c$$$ 992                 continue
c$$$!kino 2014-08-13  !$OMP end do 
c$$$!kino 2014-08-13  !$OMP end parallel 
c$$$#endif
!KINO                     write(*,*)'kino: w3p checksum=',sum(w3p)
                     deallocate(zmeltt)
c$$$                  else
c$$$!kino 2014-08-13  !$OMP parallel do 
c$$$                     do itpp= 1,ntqxx
c$$$                        do itp = 1,ntqxx
c$$$                           if(diagonly.and.(itpp/=itp)) cycle
c$$$c sep2013t a test:c  if(itpp>ntqxxd .and.itp/=itpp) cycle
c$$$                           do it  = 1,nctot+nbmax
c$$$                              w3p(it,itp,itpp) =dcmplx(
c$$$     &                          sum ( dreal(z1p(:,it,itp))*rmelt(:,it,itpp)
c$$$     &                          +   dimag(z1p(:,it,itp))*cmelt(:,it,itpp) ) ,
c$$$     &                          sum ( dimag(z1p(:,it,itp))*rmelt(:,it,itpp) 
c$$$     &                          -   dreal(z1p(:,it,itp))*cmelt(:,it,itpp) ) )
c$$$                           enddo
c$$$                        enddo
c$$$                     enddo
c$$$!kino 2014-08-13  !$OMP end parallel do 
c$$$                     deallocate(z1p)
c$$$                  endif  
                  deallocate(zmel)
c$$$!!-- Write the Spectrum function for exchange May. 2001
c$$$                  if(ifexsp/=0) then
c$$$                     do it  = 1, nctot+nbmax
c$$$                        do itp = 1,ntqxx
c$$$                           write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")
c$$$     &                          ip,itp,it, qbz_kr, ekc(it), -wtt*dreal(w3p(it,itp,itp))
c$$$                        enddo
c$$$                     enddo
c$$$                  endif
c$$$!TIME1 "end of write ifsexsp"
!TIME1 "end of w3p newaniso"
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
     &                    - wtt * sum( w3p(:,itp,itpp) )
                     enddo
                  enddo
                  deallocate( w3p,rmelt,cmelt)
c$$$                  if(.not.newaniso()) deallocate(vcoul)
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
!TIME0 
               allocate( zmel (ngb, nstate, ntqxx) )
               zmel = dcmplx (rmelt,-cmelt)
c$$$               if(newaniso()) then
!kino 2014-08-13  !$OMP parallel 
!kino 2014-08-13  !$OMP do 
c#ifdef USE_GEMM_FOR_SUM
                  if(verbose()>39)write(*,*)'info: USE GEMM FOR SUM (zmel=zmel*ppovlz) in sxcf_fal2.sc.F'
                  allocate(zmeltt(ngb, nstate, ntqxx))
                  zmeltt= dcmplx (rmelt,-cmelt)
                  call zgemm('C','N',ngb,nstate*ntqxx,ngb,(1d0,0d0),
     .                 ppovlz,ngb,zmeltt,ngb,(0d0,0d0),zmel,ngb)
                  deallocate(zmeltt)
c$$$#else
c$$$              do itp=1,ntqxx
c$$$              do it=1,nstate
c$$$                zmel(:,it,itp) =  matmul(zmel(:,it,itp),dconjg(ppovlz(:,:)))
c$$$              enddo
c$$$              enddo
c$$$#endif
!kino 2014-08-13  !$OMP end do nowait 
!kino 2014-08-13  !$OMP end parallel
c$$$               endif  
               deallocate(rmelt,cmelt)
!TIME1 "matmul zmel ppovlz"
!TIME0
               if(debug) write(6,*)' end of zmel'
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
               allocate(zmel1(ngb),zmel1_(ntqxx,ngb,nstate))
               do itpp= 1,ntqxx
                  do it  = 1,nstate
                     zmel1_(itpp,1:ngb,it) = zmel(1:ngb,it,itpp)
                  enddo
               enddo
               if(timemix) call timeshow("bbb33333 k-cycle")
               ix = 1 - nw_i    !at omega=0
               nrec=ix
c     if(bzcase()==2) nrec= (kx-1)*(nw-nw_i+1) + ix
!TIME1 "before read frcw zw"
!TIME0 
               read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v

!TIME1 "read ifrwc zw before matzwz2"
!TIME0
c     zw=(zw+dconjg(transpose(zw)))/2d0

               call matzwz2(2, zw(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
     o              zwz0)       ! zwz0 = zmel*(W(0)-v)*zmel
!TIME1 "after matzwz2"
!TIME0

ccccccccccccccccccccccccccccccccccccccccccc
!!     COH term test ----- The sum of the all states for zwz00 gives the delta function.
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
               npm = 1             ! npm=1 Timeveversal on
               if(nw_i/=0) npm = 2 ! npm=2 TimeReversal off 
               allocate(zwz(niw,nstate,ntqxx,ntqxx), zwzi(nstate,ntqxx,ntqxx)) 
               zwzi=0d0
               if(verbose()>50) call timeshow("6 before matzwz in ix cycle")

!TIME1 "before ccc333 k-cycle"
!TIME0 
              allocate( zwix(nblochpmx,nblochpmx),stat=ierr )
              if (ierr.ne.0) then 
                 write(6,'(i3,a,i3,a,i5)')
     &           mpi__rank,'failed to allocate zwix',ierr,'size=',nblochpmx**2; call flush(6) 
                 call rx( 'failed to allocate zwix')
              endif
!kino 2014-08-13  !$OMP parallel  private(zwzix) 
               write(6,*)mpi__rank,'zwzix allocate'; call flush(6)
               allocate(zwzix(1:nstate,1:ntqxx,1:ntqxx),stat=ierr)
               if (ierr.ne.0) then 
                  write(6,'(i5,a,i5,a,i7)') mpi__rank,
     &            'failed to allocate zwzix',ierr,'size=',nstate*ntqxx**2
                  call flush(6)
!kino 2014-09-18                  call mpi_abort()
                  call rx( 'failed to allocate zwzix')
               endif
!kino 2014-08-13  !$OMP do               
               do ix=1,nx
                  nrec=ix
                  read(ifrcwi,rec=nrec) zwix(:,:)  !direct access read Wc(0)=W(0)-v
                  call matzwz2(2, zwix(1:ngb,1:ngb), zmel, ntqxx, nstate,ngb,   
     o                 zwzix) ! zwz = zmel*(W(0)-v)*zmel
                  zwz(ix,1:nstate,1:ntqxx,1:ntqxx)= zwzix(1:nstate,1:ntqxx,1:ntqxx) 
               enddo ! ix
!kino 2014-08-13  !$OMP end  do 
               deallocate(zwzix) 
!kino 2014-08-13  !$OMP end parallel 
               deallocate(zwix) 
!TIME1 "end do xxx33333 k-cycle"
!TIME0

c--------------------------------------------------------------
c     S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c     Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
c     needed to take care of the singularity in the w' integration
c     when w-e(q-rk,n) is small
c--------------------------------------------------------------

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

c$$$               if(ua_auto) then
c$$$                  allocate(uaa(nstate,ntqxx))
c$$$                  do itp = 1,ntqxx
c$$$                     do  it = 1,nstate
c$$$                        ratio = abs(zwz(niw,it,itp,itp)/zwz0(it,itp,itp))
c$$$                        call gen_uaa(ratio,freqx(niw),  uaa(it,itp))
c$$$                        if(verbose()>45) then
c$$$                           write(6,"(' it itp uaa=',2i4,12f8.4)")it,itp,uaa(it,itp)
c$$$                        elseif(verbose()>40.and.mod(it,10)==1.and.mod(itp,10)==1) then
c$$$                           write(6,"(' it itp uaa=', 2i4,12f8.4)")it,itp,uaa(it,itp)
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$               endif
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
                           enddo
                           write(6,"(' sxcf_fal2sc:',3i4,' abs=',d12.4,
     &                          ' ratio= ',12f8.2)") it,itp,itpp,
     &                          abs(zwz0(it,itp,itpp)),ua2_(1:niw)
                        endif

c$$$                        if(ua_auto) then
c$$$                           ua_ = 0.5d0*(uaa(it,itp)+uaa(it,itpp))
c$$$                           call gen_expa(niw,freqx,ua_,  expa_)
c$$$                        endif

                        do ix=1,nx
                           zwzs(ix) = ( zwz(ix,it,itp,itpp)       ! w(iw) + w(-iw) symmetric part
     &                      + dconjg( zwz(ix,it,itpp,itp) ) )/2d0 ! w(iw) + w(-iw) symmetric part
                           if(npm==2) then
                              zwzs(ix+nx) = (          zwz(ix,it,itp,itpp) 
     &                             - dconjg( zwz(ix,it,itpp,itp) ) )/2d0/img
                           endif
                        enddo

c$$$                        if(GaussSmear()) then
                           zwzs0= (zwz0(it,itp,itpp) +dconjg(zwz0(it,itpp,itp)))/2d0 !22May2006 but not necessary.
                           zwzi(it,itp,itpp) =
c     Gaussian smearing
c     &      wintzsg_npm (npm, zwz(1,it,itp,itpp),zwz0(it,itp,itpp)
     &                          wintzsg_npm (npm, zwzs, zwzs0 !zwz0(it,itp,itpp)
     &                          ,freqx,wx,ua_,expa_,we,nx,esmrx)
c$$$                        else
c$$$                           if(npm==2) 
c$$$     &                          call rx( ' ###Not impliment wintzav for npm=2. Use Gausssmear.')
c$$$                           zwzi(it,itp,itpp) =  
c$$$     &                          wintzav ( zwz(1,it,itp,itpp),zwz0(it,itp,itpp)
c$$$     &                          ,freqx,wx,ua_,expa_,we,nx, esmrx)
c$$$                        endif
 1387                continue   !it
 1386             continue      !itp
 1385          continue         !itpp
               deallocate(zwzs)
!TIME1 "end 1385 loop"
!TIME0

!! sum over both occupied and unoccupied states and multiply by weight
               do     itpp= 1,ntqxx
                  do     itp = 1,ntqxx
                     if(diagonly.and.(itpp/=itp)) cycle

                     if ((iSigma_en==1 .or. iSigma_en==0) .and. itp==itpp) then
                     else
                        zsec(itp,itpp,ip)  = zsec(itp,itpp,ip) + 
     &                       wtt* sum(zwzi(:,itp,itpp)) !S_{ij}(e_i)
!     new from Jan2006! I think this should be OK.  ----------------------------
!     The output of sxcf_fal2 is  <i|Re[S](e_i)|j> ------------
!     Im-axis integral gives Hermitian part of S.
!     (Be careful as for the difference between
!     <i|Re[S](e_i)|j> and transpose(dconjg(<i|Re[S](e_i)|j>)).
!     ---because e_i is included.
!     The symmetrization (hermitian) procedure is inlucded in hqpe.sc.F
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

               if(iSigma_en==4) then
                  deallocate(zwz,zwz0,zwzi)
                  deallocate(zmel,zmel1,zmel1_,zw,freq_r) !zw4,
                  cycle ! next irot
               endif


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
                           enddo
                           write(6,"(' sxcf_fal2sc:',3i4,' abs=',d12.4,
     &                          ' ratio= ',12f8.2)") it,itp,itp,
     &                          abs(zwz0(it,itp,itp)),ua2_(1:niw)
                        endif

c$$$                        if(ua_auto) then
c$$$                           ua_= uaa(it,itp)
c$$$                           call gen_expa(niw,freqx, ua_,  expa_)
c$$$                        endif

                        do ix=1,nx
                           zwzs(ix   ) = dreal( zwz(ix,it,itp,itp)) ! w(iw) + w(-iw) symmetric part
                           if(npm==2) then
                              zwzs(ix+nx) = dimag( zwz(ix,it,itp,itp)) ! w(iw) - w(-iw)
                           endif
                        enddo

c$$$                        if(GaussSmear()) then
!! Gaussian smearing
                           zwzi(it,itp,itp) =
     &                          wintzsg_npm (npm, zwzs, zwz0(it,itp,itp)
     &                          ,freqx,wx,ua_,expa_,we,nx,esmrx)
c$$$                        else
c$$$                           if(npm==2) 
c$$$     &                          call rx( ' ###Not impliment wintzav for npm=2. Use Gausssmear.')
c$$$                           zwzi(it,itp,itp) =  
c$$$     &                          wintzav ( zwz(1,it,itp,itp),zwz0(it,itp,itp),
c$$$     &                          freqx,wx,ua_,expa_,we,nx, esmrx)
c$$$                        endif
c     .    wintz (zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx)
                     enddo      !it
                  enddo         !itp

!! sum over both occupied and unoccupied states and multiply by weight
                  do     itp = 1,ntqxx
                     zsec(itp,itp,ip)  = zsec(itp,itp,ip) + 
     &                    wtt*sum(zwzi(:,itp,itp)) !this  contribution to zsec is not hermitean
                  enddo         !itp
               endif            ! iSigma_en==1 .or. iSigma_en==0!...................
c     end of SEc w-loop

c$$$               if(ua_auto) deallocate(uaa)
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
c     write(6,*)' sxcf_fal2: onlyimagaxis cycle'
c     deallocate(freq_r,zmel,zw,zmel1)
c     goto 1100
c     endif
ccccccccccccccccccccccccccccccccccccc

c==============================================
c     contribution to SEc(qt,w) from the poles of G
c==============================================
               if(debug) write(6,*)' go to poles'
c---------------------------------------
c     maximum ixs reqired.
c---------------------------------------
               ixsmx =0
               ixsmin=0
               do 3001 itp = 1,ntqxx
                  omg  = omega(itp) 
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
c$$$                     if(GaussSmear()) then
                        if(wfac<wfaccut) cycle ! next it
                        we = .5d0*(weavx2(omg,ef,ekc(it),esmr)-omg)
c$$$                     else
c$$$                        if(wfac==0d0) cycle ! next it
c$$$                        if(omg>=ef) we = max( .5d0*(omg-ekc(it)), 0d0) ! positive
c$$$                        if(omg< ef) we = min( .5d0*(omg-ekc(it)), 0d0) ! negative
c$$$                     endif
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
               if(debug) write(6,*)'nw, nwx=',nw,nwx
               if(verbose()>50)call timeshow("10before alagr3z iw,itp,it ")

C...  Find nt_max ------------------------------------
               nt_max=nt0p      !initial nt_max
               do 4001 itp = 1,ntqxx
                  omg     = omega(itp)
                  if (omg > ef) then
                     do  it = nt0m+1,nstate ! nt0m corresponds to efm
                        wfac = wfacx2 (ef,omg, ekc(it),esmr)
c                        if( (GaussSmear().and.wfac>wfaccut)
c     &                       .or.(.not.GaussSmear().and.wfac/=0d0)) then
                        if(wfac>wfaccut) then
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
                  do igb2=1,ngb
                     do igb1=1,ngb
                        zw3(igb1,igb2,ix)=.5d0*(zw(igb1,igb2)+dconjg(zw(igb2,igb1)) )
                     enddo      !making potential hermitean, because
                  enddo         !we need only hermitean part of Sigma_nn'
               enddo
               deallocate(zw)   !,zw4
               if(timemix) call timeshow("55555 k-cycle")
c$$$               if(test_symmetric_W().and.npm==2) then
c$$$                  if(onceww(4)) write(6,*)' test_symmetric_W()=',test_symmetric_W(),nwxi,nwx
c$$$                  allocate(zw3x(ngb,ngb))
c$$$                  do ix= 1,min(abs(nwxi),nwx)
c$$$                     zw3x = 0.5d0* (zw3(:,:,ix) + zw3(:,:,-ix))
c$$$                     zw3(:,:, ix)=zw3x
c$$$                     zw3(:,:,-ix)=zw3x
c$$$                  enddo
c$$$                  deallocate(zw3x)
c$$$               endif

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
c$$$                     if(GaussSmear()) then
                        if(wfac<wfaccut) cycle ! next it
c     we = .5d0* abs( weavx2(omg,ef, ekc(it),esmr)- omg )
                        we = .5d0* abs( omg-weavx2(omg,ef, ekc(it),esmr) )
                        if(it<=nctot) then !faleev
                           if(wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
                        endif
c$$$                     else
c$$$                        if( wfac==0d0) cycle ! next it
c$$$                        if( omg >= ef) we = 0.5d0* abs( max(omg-ekc(it), 0d0) )
c$$$                        if( omg <  ef) we = 0.5d0* abs( min(omg-ekc(it), 0d0) )
c$$$                        if( it<=nctot) then !faleev
c$$$                           if(wfac/=0) call rx( "sxcf:  it<=nctot.and.wfac/=0")
c$$$                        endif
c$$$                     endif
                     if(debug)write(6,"( ' xxx1',10d13.6)") omg,ef, ekc(it),wfac
                     wfac= iii* wfac*wtt
                     do iwp = 1,nw 
                        ixs=iwp
                        if(freq_r(iwp)>we) exit
                     enddo      
                     if(nw_i==0) then
                        if(ixs+1>nwx) call rx( ' sxcf: ixs+1>nwx xxx2')
                     else
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
                     zmel1(:)=dconjg(zmel(:,it,itp))
                     if (iSigma_en == 1.or.iSigma_en==5) then
                        zwz3=(0d0,0d0)
!kino 2014-08-13  !$OMP parallel do private(ix,zz2)
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
!kino 2014-08-13  !$OMP end parallel do 
                        if(npm==1) then
                           zsec(itp,itp,ip) = zsec(itp,itp,ip) 
     .                          + wfac*alagr3z2(we,freq_r(ixs-1),zwz3,itp,itp)
                        else
                           zsec(itp,itp,ip) = zsec(itp,itp,ip) 
     .                          + wfac*alagr3z(we,freq_r(ixs-1),zwz3)
                        endif
!!  this contribution to zsec_nn is real (hermitean)
                     elseif(iSigma_en == 3) then
                        zwz4=(0d0,0d0)
!$OMP parallel private(ix,zz2)
                        do ix0=1,3
                           ix=ixs+ix0-2
!$OMP do reduction(+:zwz4)
                           do igb2=1,ngb !**** most time consuming part for iSigma_en=3 ******
                             zz2=sum(zmel1(1:ngb)*zw3(1:ngb,igb2, iir*ix)  )
                             call zaxpy(ntqxx,zz2,zmel1_(1,igb2,it),1,zwz4(1,ix0),1)
                           enddo !igb2
                        enddo   !ix
!$OMP end parallel
                        zwz44 = transpose(zwz4)
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
                     endif      ! inner iSigma_en=1 or 3
!!     this contribution to zsec_nn' is not hermitean because W(e_n)
!!     and must be made hermitean when zsec will be written on disc
 2011             continue
 2001          continue         !itp
!TIME1 "end 2001 loop"
!TIME0 
               deallocate(zw3, zmel, zwz44,zwz4)
               if(verbose()>50)call timeshow("11after alagr3z iw,itp,it cycles")
               if(debug) write(6,*)' end of do 2001 or 2002 '
 1050          continue
               if(debug) then
                  do itp = 1,ntq
                     write(6,'(" zsec=",i3,2d15.7)') itp,zsec(itp,itp,ip)
                  enddo
               endif
               deallocate(zmel1,zmel1_)
               deallocate(freq_r)
 1000       continue ! end do irot
c            if(newaniso()) ifvcoud =iclose('Vcoud.'//charnum5(kx))
            ifvcoud =iclose('Vcoud.'//charnum5(kx))
            if(.not.exchange) then
               ifrcw  = iclose('WVR.'//charnum5(kx))
               ifrcwi = iclose('WVI.'//charnum5(kx))
            endif   
!TIME1 "end of 1100 loop"
 1100    continue               ! end of kx-loop
c        if(newaniso()) ifvcoud =iclose('Vcoud')
         ifvcoud =iclose('Vcoud')
         if(irot==1) write(6,"('  sum(abs(zsec))=',d23.15)") sum(abs(zsec))
         if (allocated(vcoul))deallocate(vcoul)
 1001 continue ! end do ip
      if (allocated(expikt))deallocate(expikt)
!TIME1 "end of sxcf_fal2_scz"
      end subroutine sxcf_fal3_scz




