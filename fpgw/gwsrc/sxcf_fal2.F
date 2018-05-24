      subroutine sxcf_fal3z(kount,ixc,deltaw,shtw,qip,itq, ntq,ef,ef2,esmr,esmr2,
     i nsp,isp,                 !tiat,miat, 
     i qbas,ginv,
     i qibz,qbz,wk,nstbz,wik,  
     i nstar,irkip,  
     i freq_r,freqx,wx, 
     i dwdummy,ecore, 
     d nlmto,nqibz,nqbz,nctot,
c     d nl,nnc,nclass,natom,mdimx,
     d nbloch,ngrp, nw_i,nw ,niw,niwx,nq, !nlnmx,
     & nblochpmx ,ngpmx,ngcmx, 
     & wgt0,nq0i,q0i,symgg,alat, nband, ifvcfpout, !shtvg,
     & exchange,tote,screen,cohtest, ifexsp,
     i iwini,iwend,
     i nbmx,ebmx, 
     i wklm,lxklm, 
     i dwplot,
     o zsec,coh,exx)
      use m_readqg
      use m_readeigen,only: readeval
      use m_keyvalue,only: getkeyvalue
      use m_zmel,only: get_zmelt,
     o ppovlz, zmel,zmeltt 
      implicit none
!! TimeReversal off. when nw_i is not zero.
!! Calcualte diagonal part only version of simga_ii(e_i)= <i|Re[S](e)|i> 
!! Similar with sxcf_fal2.sc.F
Co zsec: S_ij= <i|Re[S](e)|i> where e=e_i and e_i \pm deltaw
Co
Cr  exchange=T : Calculate the exchange self-energy
Cr          =F : Calculate correlated part of the self-energy
Cr
Cr
Cr---- 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F
Cr 2000 takao kotani. This sxcf is starting from sec.f F.Aryasetiawan.
C---------------------------------------------------------------


c---- original document for sce.f (correlation case) by ferdi.Aryasetiawan.
c 92.02.24
c 93.10.18 from sec.f modified to take into account equivalent atoms
c calculates the correlated part of the self-energy SE
c SEc(q,t,t') = <psi(q,t) |SEc| psi(q,t'>
c SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >

c the zeroth order Green function
c G(r,r';w)   = S[occ]   psi(kn,r) psi(kn,r')^* /(w-e(kn)-i*delta)
c             + S[unocc] psi(kn,r) psi(kn,r')^* /(w-e(kn)+i*delta)

c the screened coulomb potential
c Wc(r,r';w)  = W(r,r';w) - v(|r-r'|)
c             = < [r1,r2] v(|r-r1|) X(r1,r2;w) v(|r2-r'|) >
c W(r,r';w)   = < [r''] ei(r,r'';w) v(|r''-r'| >
c ei          = e^(-1), inverse dielectric matrix
c             = 1 + vX
c e           = 1 - vX0 in RPA

c expand Wc(r,r';w) in optimal product basis B
c Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch]
c               B(k,i,r) Wc(k,w)(i,j) B(k,j,r')^*
c Wc(k,w)(i,j) are  the matrix elements of Wc in B

c SEc(q,t,t') = S[k=FBZ] S[n=occ]   S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)-i*delta)>
c
c             + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)+i*delta)>

c the analytic structure of GWc for w .le. ef
c                               |
c                               |   o = pole of G
c                               ^   x = pole of Wc
c                               |
c                               |   ef-w
c                               |----<-----
c                               |          |
c                 o  o  o  o  o |o  o  o   ^
c               x  x  x  x  x  x|          |
c  -----------------------------|---->------------------------------
c                               |x  x  x  x  x  x  x  x
c                               |              o  o  o  o  o
c                               |       <----->
c                               ^        gap in insulator
c                               |
c                               |

c the analytic structure of GWc for w .gt. ef
c                               |
c                               |   o = pole of G
c                               |   x = pole of Wc
c                               |
c         gap in insulator      ^
c                <----->        |
c      o  o  o  o               |
c         x  x  x  x  x  x  x  x|
c  ------------------------>----|-----------------------------------
c                   |           |x  x  x  x  x  x  x  x
c                   ^   o  o  o  o  o  o  o
c                   |           |
c                    ------<----|
c                       w-ef    |
c                               ^
c                               |
c integration along the real axis from -inf to inf is equivalent to
c the integration along the path shown
c------------------------------------------------------------
c integration along the imaginary axis: wint (s. also wint.f) (takao ->wintz)
c   (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
c the i*delta becomes irrelevant
c------------------------------------------------------------
c
c omit k and basis index for simplicity and denote e(q-k,n) = e
c wint = (i/2pi) < [w'=-inf,inf] Wc(w')/(w+w'-e) >
c
c w' ==> iw', w' is now real
c wint = - (1/pi) < [w'=0,inf] Wc(iw') (w-e)/{(w-e)^2 + w'^2} >
c
c transform: x = 1/(1+w')
c this leads to a denser mesh in w' around 0 for equal mesh x
c which is desirable since Wc and the lorentzian are peaked around w'=0
c wint = - (1/pi) < [x=0,1] Wc(iw') (w-e)x^2/{(w-e)^2 + w'^2} >
c
c the integrand is peaked around w'=0 or x=1 when w=e
c to handel the problem, add and substract the singular part as follows:
c wint = - (1/pi) < [x=0,1] { Wc(iw') - Wc(0)exp(-a^2 w'^2) }
c                          * (w-e)/{(w-e)^2 +w'^2}x^2 >
c        - (1/2) Wc(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c
c the second term of the integral can be done analytically, which
c results in the last term
c a is some constant
c
c when w = e, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w') and
c the integral becomes -Wc(0)/2
c this together with the contribution from the pole of G (s.u.)
c gives the so called static screened exchange -Wc(0)

c--------------------------------------------
c contribution from the poles of G: SEc(pole)
c--------------------------------------------
c
c for w .le. ef
c SEc(pole) = - S[k=FBZ] S[n=occ] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c             Wc(k,e(q-k,n)-w)(i,j) theta(e(q-k,n)-w)
c
c for w .gt. ef
c SEc(pole) = + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
c             Wc(k,w-e(q-k,n))(i,j) theta(w-e(q-k,n))
c
c theta(x)  = 1   if x > 0
c           = 1/2 if x = 0
c           = 0   if x < 0

c FBZ = 1st BZ
c NOTE: the routine only calculates the diagonal elements of the SE
c       i.e. SEc(q,t)

c q       = q-vector in SEc(q,t)
c itq     = states t at q
c ntq     = no. states t
c eq      = eigenvalues at q
c ef      = fermi level in Rydberg
c tr      = translational vectors in rot*R = R' + T
c iatomp(R) = R'
c ifrw,ifcw,ifrwi,ifcwi
c   = direct access unit files for Re and Im coulomb matrix
c     along real and imaginary axis
c ifrb,ifcb,ifrhb,ifchb
c         = direct access unit files for Re and Im b,hb
c qbas    = base reciprocal lattice vectors
c ginv    = inverse of qbas s. indxrk.f
cxxxxx ippb,ipdb,idpb,iddb = pointers to work array w for
c  ppb     = <phi(RLn) phi(RL'n') B(R,i)>
c  pdb     = <phi(RLn) phidot(RL'n') B(R,i)>
c  dpb     = <phidot(RLn) phi(RL'n') B(R,i)>
c  ddb     = <phidot(RLn) phidot(RL'n') B(R,i)>
c freq    = frequencies along real axis
c freqx   = gaussian frequencies x between (0,1)
c freqw   = (1-freqx)/freqx
c wx      = weights at gaussian points x between (0,1)
c ua      = constant in exp(-ua^2 w'^2) s. wint.f
c expa    = exp(-ua^2 w'^2) s. wint.f
c dw      = frequency mesh along real axis
c deltaw  = energy mesh in SEc(qt,w) ---Not used now
c iclass  = given an atom, tells the class
c wk      = weight for each k-point in the FBZ
c indexk  = k-point index
c qbz     = k-points in the 1st BZ
c nstar   = no. stars for each k
c irk(k,R,nq) = gives index in the FBZ with k{IBZ, R=rotation
c mdim    = dimension of B(R,i) for each atom R
c work arrays:
c rbq,cbq     = real and imaginary part of b(q)
c rhbq,chbq   = real and imaginary part of hb(q)
c rbkq,cbkq   = real and imaginary part of b(q-k)
c rhbkq,chbkq = real and imaginary part of hb(q-k)
c   b is the eigenvector of the LMTO-Hamiltonian
c ekq     = eigenvalues at q-k
c rmel,cmel = real and imaginary part of
c             <psi(q,t') | psi(q-k,t) B(k,R,i)>
c wr1 ... = work arrays
c dimensions:
c nqibz   = number of k-points in the irreducible BZ
c n1,n2,n3= divisions along base reciprocal lattice vectors
c natom   = number of atoms
c nctot   = no. allowed core states
c nbloch  = total number of Bloch basis functions
c nlnmx   = maximum number of l,n,m
c nlmto   = total number of LMTO basis functions
c ngrp    = no. group elements (rotation matrices)
c niw     = no. frequencies along the imaginary axis
c nw      = no. frequencies along the real axis
c niwx    = max(niw,nw)
c
c secq(t) = <psi(q,t) |SEc| psi(q,t)>
c----------------------------------------------------------------------
      intent(in)::
     i kount,ixc,deltaw,shtw,qip,itq, ntq,ef,ef2,esmr,esmr2,
     i nsp,isp,                 !tiat,miat, 
     i qbas,ginv,
     i qibz,qbz,wk,nstbz,wik,  
     i nstar,irkip,  
c     i iclass,mdim,nlnmv,nlnmc,
c     i icore,ncore,imdim,
c     i ppb, 
     i freq_r,freqx,wx, 
     i dwdummy,ecore, 
     d nlmto,nqibz,nqbz,nctot,
c     d nl,nnc,nclass,natom,mdimx,
     d nbloch,ngrp, nw_i,nw ,niw,niwx,nq, !nlnmx,
     & nblochpmx ,ngpmx,ngcmx, 
     & wgt0,nq0i,q0i,symgg,alat, nband, ifvcfpout, !shtvg,
     & exchange,tote,screen,cohtest, ifexsp,
     i iwini,iwend,
     i nbmx,ebmx, 
     i wklm,lxklm 
c     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, 
c     i invg,!il,in,im,nn_, lx,nx_,nxx_,dwplot !ppbrd, !cgr,,nlnm

      integer :: ntq, nqbz,nqibz,ngrp,nq,nw,niw, !natom,
     & nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,nstatex,
     & irot,  iqisp,ikpisp,isp,nsp, !nlnmx, !iq, idxk,
c     &  iwr1,iwr2,iwr3,iwr4,iwc1,iwc2,iwc3,iwc4
     & ip, it,itp,              !ifcphi,    ! ifrb,ifcb,ifrhb,ifchb, 
c     i iiclass,                 !mdim(*),
     i ifrcw,ifrcwi,            !iindxk,
     & ifvcfpout,ndummy1,ndummy2,kx,kr,ngc,ngb,nbloch, !n1,n2,n3, k,
     & kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0,
     & nbmax,nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, !nclass,nl,nnc, 
     &  nwx,niwx,
     & itq(ntq), !,iatomp(natom),  !,miat(natom,ngrp),
     & nstar(nqibz),irkip(nqibz,ngrp,nq),kount(nqibz,nq)
c
      real(8) :: q(3),qbas(3*3),ginv(3*3), !tr(3,natom), !,tiat(3,natom,ngrp)
     & wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),
     & freqx(niw),wx(niw),      !expa(niw),
     & eq(nband,nq), 
     & ekq(nband), ekc(nctot+nband),
     & tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx, !ua,
     & dwdummy,wtt,wexx,www,exx,exxq ,wfacx2,weavx2,wex
c      complex(8) :: zsec(-1:1,ntq,nq)
c      real(8)    ::  shtw
c                       ! This shft is  to avoid some artificial resonance effects.
c                       ! shtw can be zero for esmr/=0 given by takao.
c
      integer(4):: ngpmx, ngcmx, !ngcni(nqibz), !ngpn(nqbz),
     & igc,                     !ngvecpB(3,ngpmx,nqbz),ngveccBr(3,ngcmx,nqibz),
     & nadd(3)
      real(8) :: wgt0(nq0i,ngrp),qk(3), !qfbz(3),
     & qdiff(3),add(3),symgg(3,3,ngrp),symope(3,3), !qbasinv(3,3), det,
     & qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot), !shtvg(3,ngrp),
c     &   ppb(1), !pdb(1),dpb(1),ddb(1), !*
     & coh(ntq,nq)              !, pos(3,natom)
      complex(8)::   alagr3zz,wintz !geigB  (ngpmx,nband,nqbz),

c
c      real(8),allocatable:: !rmel(:,:,:),cmel(:,:,:),
c     &                      rmelt(:,:,:),cmelt(:,:,:)
      complex(8),allocatable :: zz(:),zzmel(:,:,:),
     & zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:)
c for exchange --------------------
      logical :: exchange,screen,cohtest,tote
      real(8),allocatable:: 
     & w1p(:,:,:),w2p(:,:,:),w3p(:,:)
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
      integer:: invrot,invr
c      integer:: invg(ngrp),il(*),in(*),im(*),nn_,lx(*),nx_(*),nxx_ !nlnm(*),
c      real(8)::  cgr(*),ppbrd(*)

c- debugwrite ---------------------
      logical :: debug=.false. ,onceww

cccccccccccccc
c tetra
c      integer(4) :: ntqx
c      integer(4) :: ibzx(nqbz)
c      real(8)    :: wtet  (nband,nqibz,1:ntqx), wtetef(nband,nqibz)
c                ! where the last index is 3*itq+iw-1,itq=1,ntq,iw=-1,1
c     logical    :: tetraex
cccccczzcccccc

      complex(8) :: wintzav,wintzsg_npm,wintzsg

      integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen
      integer(4),save::ifzwz=-999

      integer(4) :: iwini, iwend, ia
      real(8)    :: esec, omega(ntq, iwini:iwend)
      complex(8) :: zsec(iwini:iwend,ntq,nq)
c      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)
ctakao
c      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)

      integer(4) :: nt_max, igb1,igb2,iigb,  nw_i !nw_i is at feb2006 TimeReversal off case
      complex(8),allocatable:: zmel3(:) !zmel1(:),
      complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)
      complex(8), allocatable :: zwz2(:,:),zw2(:,:),zmel2(:,:) !0 variant
      complex(8) ::  zz2 ,zwz3(3) ,zwz3x
      real(8) :: dd,omg_c,dw2,omg
      real(8) :: freq_r(nw_i:nw)
      complex(8), allocatable :: zw3(:,:,:)


      real(8)::weavx,wfaccut=1d-10,qqqq

      logical :: GaussSmear=.true.,gass
      real(8) :: ebmx,ddw
      integer(4):: nbmx,nbmxe,nstatetot

c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

c      integer(4)::icore(*),ncore(*),imdim(*) !,iclass(*),nlnmv(*),nlnmc(*),

      integer(4)::verbose,nstbz(nqbz),bzcase=1,iqini,iqend
      real(8):: wgtq0p

      integer(4):: nrec,kxx
      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      logical ::  onlyimagaxis 

      logical ::zwz3mode


      real(8):: ua_,expa_(niw),ua2,freqw,freqw1,ratio,ua2_(niw)
c$$$      logical :: ua_auto !fixed to be .false.
      integer(4):: icc=0
      real(8),allocatable:: uaa(:,:)

c      logical ::testimx=.false.
ccccc zvz test cccccccccccccccccccccccccc
      integer(4):: ngbx
c      complex(8):: vcoul(ngbx,ngbx)
      complex(8),allocatable:: vzz(:,:,:),aaa(:), zwzs(:)
      complex(8):: zvz,zvz1
      integer(4):: ib1,ib2,ifix
ccccccccccccccccccccccccccccccccccc
      logical ::iww2=.true., oncew


C...
c      logical::smbasis
      integer(4):: iclose,isx,iqx !nn,no,ifpomat,
c      complex(8),allocatable:: pomat(:,:)
      real(8):: q_r(3)
c      integer(4):: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
c      complex(8):: pomatr(nnmx,nomx,nkpo)
c      real(8):: qrr(3,nkpo)

      real(8):: elxx,ehxx,ekxx,efxx
      integer(4):: ixsmin,iwm,iir,nwxi, itini,itend, npm
      real(8)   :: fffr(3),ppp
      complex(8):: zwzz(3)

      real(8),allocatable:: ebb(:)
      integer(4):: ii,iq
      logical ::evaltest        !, imgonly

      integer:: lxklm,ivc,ifvcoud,idummy,iy,ngb0
      real(8):: wklm((lxklm+1)**2),pi,fpi,vc,qvv(3),aaaa
      complex(8)::zmelt1,zmelt0
      real(8)::voltot
c      logical :: newaniso !fixed to be T

      complex(8),allocatable:: ppovl(:,:),zcousq(:,:) !,ppovlz(:,:)
      real(8),allocatable::vcoud(:),vcousq(:)
      integer:: mrecl,nprecx,ifwd
      character(5):: charnum5

      integer:: ixc
      real(8):: qip(3,*),deltaw,shtw,eqx(nband),dwplot,tolq=1d-8
      complex(8),allocatable:: zmelt(:,:)
      integer:: ntqxx,nrot
c--------------------------------------------------------------------
      write(6,*)'sxcf_fal3z'
c      timemix=.false.
      pi  = 4d0*datan(1d0)
      fpi = 4d0*pi
      debug=.false.
      if(verbose()>=90) debug=.true.
!!
      if(.not.exchange) then
        ifwd = iopen('WV.d',1,-1,0)
        read (ifwd,*) nprecx,mrecl
        ifwd = iclose('WV.d')
c$$$!! --- gauss_img : interpolation gaussion for W(i \omega).
c$$$      call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
c$$$      if(ua_<=0d0) then
c$$$        ua_auto =.true.
c$$$        write(6,"(' ua_auto=T')")
c$$$      else
c$$$        ua_auto =.false.
c$$$        do ix = 1,niw
c$$$          freqw     = (1d0 - freqx(ix))/ freqx(ix)
c$$$          expa_(ix) = exp(-(ua_*freqw)**2)
c$$$        enddo
c$$$      endif
        call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
        do ix = 1,niw           !! Energy mesh; along im axis.
          freqw     = (1d0 - freqx(ix))/ freqx(ix)
          expa_(ix) = exp(-(ua_*freqw)**2)
        enddo
        npm = 1                 ! npm=1    Timeveversal case
        if(nw_i/=0) npm = 2     ! npm=2 No TimeReversal case. Need negative energy part of W(omega)
      endif 

      tpi         = 8d0*datan(1.d0)
      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
c      call dinv33(qbas,0,qbasinv,det)
c      allocate(expikt(natom))


!! == ip loop to spedify external q ==
      do 1001 ip = 1,nq   
        if(sum(irkip(:,:,ip))==0) cycle
        q = qip(1:3,ip) 
        write (*,*) ip,'  out of ',nq,'  k-points ' ! call cputid  (0)
        if(ixc==2) then
          call readeval(q,isp,eqx)
          do iw = iwini,iwend
            do i  = 1,ntq
              omega(i,iw) = eqx(itq(i)) + 2d0*(dble(iw)-shtw)*deltaw
            enddo
          enddo
        endif
!!
        if(ixc==4) then
c        dwplot=0.01
        do iw = iwini,iwend
          omega(1:ntq,iw) =  dwplot* iw + ef
        enddo
        endif

        call readeval(q, isp, eq(1,ip))
!! we only consider bzcase()==1
        if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
        if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
        If (tote) exxq = 0.d0

!! == Big loop for kx ==
!! kx is for irreducible k points, kr=irk(kx,irot) runs all k points in the full BZ.
        iqini=1
        iqend=nqibz             !no sum for offset-Gamma points.
        do 1100 kx = iqini,iqend
          if(sum(irkip(kx,:,ip))==0) cycle
          write(6,*) ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
c          if( kx <= nqibz ) then
            qibz_k= qibz(:,kx)
c          else
c            qibz_k= 0d0
c         endif
          if(verbose()>=40)  write(6,*) ' sxcf_fal3z: loop 1100 kx=',kx
          call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
          ngb = nbloch + ngc    !oct2005
          if(debug) write(6,*) ' sxcf: ngb=',ngb,nbloch

!! ===Readin diagonalized Coulomb interaction===
!!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
!!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
!!  Vould file is written in hvccfp0.m.F.
!! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
!! These can be also separeted into WVR.ID and WVI.ID files.
!! NOTE: vcoud and zcousq are in module m_zmelt.
c          if(kx<=nqibz) qxx=qibz_k
c          if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
          qxx=qibz_k
          ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
          do
            read(ifvcoud) ngb0
            read(ifvcoud) qvv
            if(allocated(vcoud)) deallocate(vcoud)
            if(allocated(zcousq)) deallocate(zcousq)
            allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
            read(ifvcoud) vcoud
            read(ifvcoud) zcousq
            if(sum(abs(qvv-qxx))<tolq) goto 1133
          enddo
          if(sum(abs(qvv-qxx))>tolq) then
            write(6,*)'qvv =',qvv
            write(6,*)'qxx=',qxx,kx
            call rx( 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent')
          endif
 1133     continue
          if( ngb0/=ngb ) then  !sanity check
            write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
            call rx( 'hsfp0.m.f:ngb0/=ngb')
          endif
!! used in get_zmel
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
!! ppovlz= O Zcousq
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
          if(allocated(ppovl)) deallocate(ppovl,ppovlz)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
          call readppovl0(qibz_k,ngc,ppovl)
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
          ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
          deallocate(zcousq)
!! === open WVR,WVI for correlation mode ===
          if(.not.exchange) then
            ifrcw  = iopen('WVR.'//charnum5(kx),0,-1,mrecl)
            ifrcwi = iopen('WVI.'//charnum5(kx),0,-1,mrecl)
          endif
          nrot=0
          do irot = 1,ngrp
c            if( kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
              if(kr==0) cycle   ! next irot
              qbz_kr= qbz (:,kr) 
c            else
c              kr=-99999         !for sanity check
c              qbz_kr= 0d0
c              if( wgt0(kx-nqibz,irot)==0d0 ) cycle ! next irot
c            endif
            nrot=nrot+1
          enddo  

!! == loop over rotations ==
!! We may extend 
          do 1000 irot = 1,ngrp
c            if( kx <= nqibz) then
              kr = irkip(kx,irot,ip) ! index for rotated k in the FBZ
              if(kr==0) cycle
              qbz_kr= qbz (:,kr) 
c            else
c              kr=-99999         !for sanity check
c              qbz_kr= 0d0
c              if( wgt0(kx-nqibz,irot)==0d0 ) cycle
c            endif
            write(*,"('ip,kx irot=',3i5, ' out of',2i4)") ip,kx,irot, iqend,ngrp

c qk = q - rk, rk is inside 1st BZ, not restricted to the irreducible BZ
            qk =  q - qbz_kr    ! qbz(:,kr)
            call readeval(qk, isp, ekq)
            ekc(nctot+1:nctot+nband) = ekq (1:nband)
            nt0 = nocc (ekc,ef,.true.,nstatetot)
            ddw= .5d0
c        if(GaussSmear) ddw= 10d0
            ddw=10d0
            efp= ef+ddw*esmr
            efm= ef-ddw*esmr
            nt0p = nocc (ekc,efp,.true.,nstatetot)
            nt0m = nocc (ekc,efm,.true.,nstatetot)
!! nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!! nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W
            if(exchange) then
              nbmax = nt0p-nctot
              if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
            else
              nbmax = nband
              nbmxe = nocc (ekc,ebmx,.true.,nstatetot)-nctot
              nbmax  = min(nband,nbmx,nbmxe)
              if(onceww(3)) write(6,*)' nbmax=',nbmax
            endif
            nstate = nctot + nbmax ! = nstate for the case of correlation

!! all are identical.
            ntp0 = ntq
            ntqxx= ntp0

!! Get matrix element zmelt= rmelt + img*cmelt, defined in m_zmel.F---
            if(debug) write(6,*)'zzBBB ppovlz =',sum(abs(ppovlz(:,:))),kx,irot
            if(allocated(zmel)) deallocate(zmel)
            if(allocated(zmeltt)) deallocate(zmeltt)
! this return zmeltt (for exchange), or zmel (for correlation)
            call get_zmelt(exchange,q,kx,qibz_k,irot,qbz_kr,kr,isp,
     &       ngc,ngb,nbmax,ntqxx,nctot,ncc=0)
            if(kx<= nqibz) then
              wtt = wk(kr)      !         wtx = 1d0
            else
              wtt = wk(1)*wgt0(kx-nqibz,irot) !       wtx = wgt0(kx-nqibz,irot)
              if(abs(wk(1)-1d0/dble(nqbz))>1d-10)call rx( 'sxcf:wk(1)inconsistent')
            endif
            if(debug) write(6,*) 'ssssssss',size(zmel),ntqxx*nstate*ngb
            if(debug) write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
            if(debug) write(6,*)' 000 sumzmel=',ngb, nstate, ntp0,sum(abs(real(zmel))),sum(abs(imag(zmel)))

!!--------------------------------------------------------
!! === exchange section ===
!!--------------------------------------------------------
c
c S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c                        v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
c
c> z1p(j,t,t') = S[i=1,nbloch] <psi(q,t') | psi(q-rk,t) B(rk,i)> v(k)(i,j)
!! NOTE: zmel(igb, nctot+nbmax, ntp0) ---> <phi phi  |igb>

c --- screened exchange case ----
c          if(screen) then
c            ix = 1
c            nrec=(kx-iqini)*nw+ix
c            if(bzcase()==2) nrec= (kx-1)*nw+ix
c            read(ifrcw,rec=nrec) zw  ! Readin W(0) - v
c            vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
c          endif

c          allocate( zmel(ngb, nctot+nbmax, ntp0), w3p( nctot+nbmax,ntp0))
c          zmel  = dcmplx (rmelt,cmelt)
            if(exchange) then
              allocate( w3p( nctot+nbmax,ntp0))
              do 992 itp = 1,ntp0
                do 993 it  = 1,nctot+nbmax
                  w3p(it,itp) = 0d0
                  do 994 ivc=1,ngb
                    if(ivc==1.and.kx==1) then
                      vc= wklm(1)* fpi*sqrt(fpi) /wk(kx)
c                    write(6,*)'wklm(1) vc=',wklm(1),vc
                    else
                      vc= vcoud(ivc)
                    endif
                    w3p(it,itp) = w3p(it,itp)+ vc * abs(zmeltt(it,itp,ivc))**2
 994              continue
 993            continue
 992          continue
              if(debug) then
                do  it  = 1,nctot+nbmax
                  do  itp = 1,ntp0
                    write(6,"(' w3p =',2i4,2d14.6)") it,itp,w3p(it,itp)
                  enddo
                enddo
              endif

!! Write the Spectrum function for exchange May. 2001. 
!!!!!! Probably, Need to fix this....
              if(ifexsp/=0) then
                do it  = 1, nctot+nbmax
                  do itp = 1,ntp0
                    write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")
     &               ip,itp,it, qbz_kr, ekc(it), -wtt*w3p(it,itp)
                  enddo
                enddo
              endif

!! --- Correct weigts wfac for valence by esmr
              do it = nctot+1, nctot+nbmax
                wfac = wfacx(-1d99, ef, ekc(it), esmr) !gaussian
                w3p(it,1:ntp0) = wfac * w3p(it,1:ntp0)
              enddo

              if (.not.tote) then !total energy mode tote
                do itp = 1,ntp0 !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
                  zsec(iwini,itp,ip) = zsec(iwini,itp,ip) 
     &             - wtt * sum( w3p(:,itp) )
                enddo
              else
                do itp = 1,ntp0
                  wfac = wfacx(-1d99, ef2, eq(itq(itp),ip), esmr2) !gaussian
                  w3p(1:nctot+nbmax,itp) = wfac * w3p(1:nctot+nbmax,itp)
                  exxq = exxq - wtt * sum( w3p(:,itp) )
                enddo
              endif
              deallocate( w3p)  !,rmelt,cmelt)
              cycle
            endif
c-- End of exchange section --------------



c--------------------------------------------------------------------------
c--- correlation section --------------------------------------------------
c--------------------------------------------------------------------------
c$$$c--- The matrix elements zmel.
c$$$c        allocate( zmel (ngb, nstate, ntp0) )
c$$$c        zmel = dcmplx (rmelt,-cmelt)
c$$$c        if(newaniso) then
c$$$c#ifdef USE_GEMM_FOR_SUM
c$$$          if(verbose()>39)write(*,*)'info: USE GEMM FOR SUM (zmel=zmel*ppovlz), in sxcf_fal2.F'
c$$$          allocate( zmelt (ngb, nstate) )
c$$$          do itp=1,ntp0
c$$$          zmelt = dcmplx(rmelt(:,:,itp),-cmelt(:,:,itp))
c$$$          call zgemm('C','N',ngb,nstate,ngb,(1d0,0d0),
c$$$     .      ppovlz,ngb,zmelt,ngb,(0d0,0d0),zmel(1,1,itp),ngb)
c$$$          enddo
c$$$          deallocate(zmelt)
c$$$#else
c$$$          do itp=1,ntp0
c$$$            do it=1,nstate
c$$$              zmel(:,it,itp) =  matmul(zmel(:,it,itp),dconjg(ppovlz(:,:)))
c$$$            enddo
c$$$          enddo
c$$$#endif
c$$$c        endif
c        deallocate(rmelt,cmelt)
c        if(debug) write(6,*)' end of zmel'

c================================================================
c The correlated part of the self-energy:
c S[n=all] S[i,j=1,nbloch]
c <psi(q,t) |psi(q-rk,n) B(rk,i)>
c  < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
c                                <B(rk,j) psi(q-rk,n) |psi(q,t)>
c e = e(q-rk,n), w' is real, Wc = W-v
c================================================================
            allocate( zw (nblochpmx,nblochpmx) )
c====================================================================
c contribution to SEc(qt,w) from integration along the imaginary axis
c====================================================================
c------------------------------------------------
c loop over w' = (1-x)/x, frequencies in Wc(k,w')
c {x} are gaussian points between (0,1)
c------------------------------------------------
            allocate( zwz0(nstate,ntp0) )
            ix = 1  - nw_i      !at omega=0
c        nrec=(kx-iqini)*(nw-nw_i+1) +ix ! 2---> iqini
c        if(bzcase()==2) nrec= (kx-1)*(nw-nw_i+1) +ix
            nrec=ix 
            if(debug) write(6,*)' wvr nrec kx nw nw_i ix=',nrec,kx,nw,nw_i,ix
            read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v
            zwz0=0d0
!! this loop looks complicated but just in order to get zwz0=zmel*zwz0*zmel
!! Is this really efficient???
CCC!$OMP parallel do private(itp,it,igb2,zz2)
            do itp=1,ntp0
              do it=1,nstate
                do igb2=2,ngb
                  zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )
                  zwz0(it,itp) = zwz0(it,itp)+zz2*zmel(igb2,it,itp)*2d0+
     &             dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp)
                enddo           !igb2
                zwz0(it,itp) = zwz0(it,itp)+
     &           dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
              enddo             !it
            enddo               !itp
            zwz0 = dreal(zwz0)
c COH term test ----- The sum of the all states for zwz00 gives the delta function.
            if(cohtest) then
              do itp = 1,ntq
                coh(itp,ip)  = coh(itp,ip) 
     &           + .5d0*wtt*sum(dreal(zwz0(1:nstate,itp)))
              enddo
              deallocate(zw,zwz0,zmel)
              cycle
            endif
c
            nx  = niw
            if(niw <1) call rx( " sxcf:niw <1")
            if(allocated(zwz)) deallocate(zwz)
            if(allocated(zwzi)) deallocate(zwzi)
            allocate( zwz(niw*npm, nstate,ntp0),  zwzi(nstate,ntp0) )
            if(screen) allocate( zwz00(nstate,ntp0) )
            if(verbose()>50) write(*,'("6 before matzwz in ix cycle ",$)')
            if(verbose()>50) call cputid(0)

            zwz=0d0       
            do ix = 1,nx        !*npm            ! imaginary frequency w'-loop
              nrec= ix
              if(debug) write(6,*)' wvi nrec=',nrec
              read(ifrcwi,rec=nrec) zw ! Readin W-v on imag axis
              if(npm==1) then   !then zwz is real so, we can use mode c2.
                do itp= 1,ntp0
                  do it = 1,nstate
                    ppp=0d0
                    do igb2 = 2,ngb
                      zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )
! only take real part
                      ppp = ppp + dreal(zz2*zmel(igb2,it,itp)) * 2d0
     &                 + dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp)
                    enddo       !igb2
                    zwz(ix,it,itp) = ppp +
     &               dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
                  enddo         !it
                enddo           !itp
              else              !we need to use mode2 because zwz is not real now.
                call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb, 
     o           zwz(ix,1:nstate,1:ntp0))
              endif
              if(debug) write(6,*)' sumzw=',sum(abs(zw))
            enddo               !ix
            if(verbose()>50) write(*,'("xxx:6.1 before matzwz in ix cycle ",$)')
            if(verbose()>50) call cputid(0)
            if(debug) write(6,*)' sumzmel=',ngb, nstate, ntp0,sum(abs(real(zmel))),sum(abs(imag(zmel)))

c--------------------------------------------------------------
c S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
c needed to take care of the singularity in the w' integration
c when w-e(q-rk,n) is small
c--------------------------------------------------------------
            if(screen) then
              zwz00 = zwz0
              zwz0  = 0d0
              do ix = 1,nx
                zwz(ix,:,:)=zwz(ix,:,:) - zwz00
              enddo
            endif

c------------------------------------------------
c loop over w in SEc(qt,w)
c------------------------------------------------
c$$$        if(ua_auto) then
c$$$          allocate(uaa(nstate,ntq))
c$$$          do itp = 1,ntq
c$$$            do  it = 1,nstate
c$$$              ratio = abs(zwz(niw,it,itp)/zwz0(it,itp))
c$$$              call gen_uaa(ratio,freqx(niw),  uaa(it,itp))
c$$$              if(verbose()>45) then
c$$$                write(6,"(' it itp uaa=',2i4,12f8.4)") it,itp,uaa(it,itp)
c$$$              elseif(verbose()>40.and.mod(it,10)==1.and.mod(itp,10)==1) then
c$$$                write(6,"(' it itp uaa=', 2i4,12f8.4)") it,itp,uaa(it,itp)
c$$$              endif
c$$$            enddo
c$$$          enddo
c$$$        endif
            allocate(zwzs(npm*nx))
            do iw = iwini,iwend
c frequency integration along the imaginary axis, s. wint.f
c for each e(q-rk,n) and w in SEc(qt,w)
              do 1385  itp = 1,ntq
                do 1387 it = 1,nstate
                  we =.5d0*( omega(itp,iw) -ekc(it)) != .5d0*( eq(itq(itp),ip)+2d0*(dble(iw)-shtw)*deltaw-ekc(it))
                  if(verbose()>50) then
                    do  ix = 1,niw
                      ratio  = abs(zwz(ix,it,itp)/zwz0(it,itp)) 
                      freqw1 = (1d0 - freqx(ix))/ freqx(ix)
                      ua2_(ix) = sqrt(- 1d0/freqw1*log(ratio))
                    enddo
                    write(6,"(' sxcf_fal2: ua=sqrt(1/w1*log(v0/v1))=',12f8.4)") ua2_(1:niw)
                  endif
c          if(ua_auto) then
c            call gen_ua(abs(zwz(niw,it,itp)/zwz0(it,itp)), niw,freqx, expa_,ua_)
c            if(iw==ini) then
c            if(verbose()>45) then
c              write(6,"(' it itp ua_=',2i4,12f8.4)")it,itp,ua_
c            elseif(verbose()>40.and.mod(it,20)==1.and.mod(itp,20)==1) then
c              write(6,"(' it itp ua_=',3i4,12f8.4)")it,itp,ua_
c            elseif(irot==1.and.mod(it,10)==1.and.itp==it) then
c              write(6,"(' it itp ua_=',3i4,12f8.4)")it,itp,ua_
c            endif
c            endif
c          endif
c$$$              if(ua_auto) then
c$$$                ua_ = .5d0*uaa(it,itp)
c$$$                call gen_expa(niw,freqx,ua_,  expa_)
c$$$              endif
                  esmrx = esmr
                  if(it <= nctot) esmrx = 0d0
                  do ix=1,nx
                    zwzs(ix   ) = dreal( zwz(ix,it,itp)) ! w(iw) + w(-iw) symmetric part
                    if(npm==2) then
                      zwzs(ix+nx) = dimag( zwz(ix,it,itp)) ! w(iw) - w(-iw)
                    endif
                  enddo
c                  if(GaussSmear) then
                    zwzi(it,itp) = 
     &               wintzsg_npm(npm, zwzs,zwz0(it,itp),freqx,wx,ua_,expa_,we,nx, esmrx)
c                  else
c                    if(npm==2) 
c     &               call rx( ' ###Not impliment wintzav for npm=2. Use Gausssmear.')
c                    zwzi(it,itp) = 
c     &               wintzav( zwzs,zwz0(it,itp),freqx,wx,ua_,expa_,we,nx, esmrx)
c                  endif
c    .    wintz (zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx)
ccccccccccccccccccccccccccccccc
c          if(verbose()>45) then
c          if(it==50.and.itp==1) then
c          write(6,"(' it itp abs(zwzi)=',2i4,12d13.5)")it,itp,abs( zwzi(it,itp))
c          icc=icc+1
c          if(icc==10) stop 'test end'
c          endif
c          endif
cccccccccccccccccccccccccccccc
 1387           continue
 1385         continue
c sum over both occupied and unoccupied states and multiply by weight
              do  itp = 1,ntq
                zsec(iw,itp,ip)  = zsec(iw,itp,ip) + wtt*sum(zwzi(:,itp))
              enddo
c end of SEc w-loop
            enddo
            deallocate(zwzs)
            if(debug) then
              write(6,*)' ntq nstate sum(zwzi)=',ntq,nstate,sum(zwzi)
              write(6,*)' ntq nstate sum(zwz )=',ntq,nstate,sum(zwz)
              do itp = 1,ntq
                write(6,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
              enddo
            endif
            deallocate(zwz,zwz0,zwzi)

c===============================================================================
c contribution to SEc(qt,w) from the poles of G
c===============================================================================
!    We assume freq_r(i) == -freq_r(-i) in this code. feb2006
c---------------------------------------
c maximum ixs finder
c---------------------------------------
c      write(6,*)' ekc at nt0p nt0m+1=', ekc(nt0p),ekc(nt0m+1)
c      write(6,*)'  nt0p nt0m+1=', nt0p, nt0m+1
            ixsmx =0
            ixsmin=0
            do 3001 iw  = iwini,iwend
              do 3002 itp = 1,ntq
                omg = omega(itp,iw)
                if (omg < ef) then
                  itini= 1
                  itend= nt0p
                else
                  itini= nt0m+1
                  itend= nstate
                endif
                do 3011 it= itini,itend
                  esmrx = esmr
                  if(it<=nctot) esmrx = 0d0
                  wfac = wfacx2(omg,ef, ekc(it),esmrx)
                  if(GaussSmear) then
                    if(wfac<wfaccut) cycle
                    we = .5d0*(omg-weavx2(omg,ef,ekc(it),esmr))
                  else
                    if(wfac==0d0) cycle
                    if(omg>=ef) we = max( .5d0*(omg-ekc(it)), 0d0) ! positive
                    if(omg< ef) we = min( .5d0*(omg-ekc(it)), 0d0) ! negative
                  endif
                  do iwp  = 1,nw ! may2006
                    ixs = iwp   ! ixs = iwp= iw+1
c                write (*,*) 'xxx freq we=',freq_r(iwp),abs(we)
                    if(freq_r(iwp) > abs(we)) exit
                  enddo
c This change is because G(omega-omg') W(omg') !may2006
c             if(ixs>ixsmx  .and. omg<=ef ) ixsmx  = ixs
c             if(ixs>ixsmin .and. omg> ef ) ixsmin = ixs
                  if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
                  if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
                  wexx  = we
                  if(ixs+1 > nw) then
                    write (*,*) ' nw_i ixsmin',nw_i, ixsmin
c                    write (*,*) ' wexx, dw ',wexx,dw
                    write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
Cstop2rx 2013.08.09 kino                stop ' sxcf 222: |w-e| out of range'
                    call rx( ' sxcf 222: |w-e| out of range')
                  endif
 3011           continue
 3002         continue          !end of SEc w and qt -loop
 3001       continue            !end of SEc w and qt -loop
            if(nw_i==0) then
              nwxi = 0
              nwx  = max(ixsmx+1,ixsmin+1)
            else
              nwxi = -ixsmin-1
              nwx  =  ixsmx+1
            endif
            if (nwx > nw   ) then
              call rx( ' sxcf nwx check : |w-e| > max(w)')
            endif
            if (nwxi < nw_i) then
              call rx( ' sxcf nwxi check: |w-e| > max(w)')
            endif
            if(debug) write(6,*)' nwxi nwx nw=',nwxi,nwx,nw

C... Find nt_max ------------------------------------
            nt_max=nt0p         !initial nt_max
            do 4001 iw  = iwini,iwend
              do 4002 itp = 1,ntq
                omg     = omega(itp,iw)
                if (omg > ef) then
                  do  it = nt0m+1,nstate ! nt0m corresponds to efm
                    wfac = wfacx2 (ef,omg, ekc(it),esmr)
                    if( (GaussSmear.and.wfac>wfaccut)
     &               .or.(.not.GaussSmear.and.wfac/=0d0)) then
                      if (it > nt_max) nt_max=it ! nt_max is  unocc. state
                    endif       ! that ekc(it>nt_max)-omega > 0
                  enddo
                endif
 4002         continue
 4001       continue

C... Set zw3 or zwz -----------------------------------
            zwz3mode=.true.
            if(iwend-iwini>2) then
              zwz3mode=.false.
            endif
            if(zwz3mode) then
              allocate( zw3(ngb,ngb,nwxi:nwx))
              do ix = nwxi,nwx  ! real frequency w'-loop
                nrec=ix-nw_i+1
                if(debug) write(6,*)' wvr3 nrec=',nrec,nblochpmx,kx,ix,nw
                read(ifrcw,rec=nrec) zw
                zw3(1:ngb,1:ngb,ix) = zw(1:ngb,1:ngb)
                if(evaltest()) then
                  write(6,"('iii --- EigenValues for zw --------')")
                  allocate(ebb(ngb))
                  call diagcvh2((zw(1:ngb,1:ngb)-transpose(dconjg(zw(1:ngb,1:ngb))))/2d0/img, 
     &             ngb, ebb)
                  do ii=1,ngb
                    if(abs(ebb(ii))>1d-8.and.ebb(ii)>0) then
                      write(6,"('iii1xxx:  iw ii eb=',2i4,d13.5)") ix,ii,ebb(ii)
                    else
                      write(6,"('iii1:  iw ii eb=',2i4,d13.5)") ix,ii,ebb(ii)
                    endif
                  enddo
                  deallocate(ebb)
                endif
              enddo
              deallocate(zw)
            else
              nstatex= max(ntp0,nt_max)
              if(allocated(zwz)) deallocate(zwz)
              allocate( zwz(nwxi:nwx,1:nstatex,ntp0) )
              do      ix = nwxi,nwx
                nrec= ix-nw_i+1
                read(ifrcw,rec=nrec) zw ! Readin (W-v)(k,w')(i,j) at k and w' on imag axis
c zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
                call matzwz(zw(1:ngb,1:ngb), zmel(1:ngb,1:nstatex,1:ntp0), ntp0,nstatex,ngb,   
     o           zwz(ix,1:nstatex,1:ntp0))
! zmel (ngb, nstate, ntp0)
              enddo
              deallocate(zmel)
              deallocate(zw)
            endif
c---------------------------------------------
            if(screen) then
              if(zwz3mode) call rx( ' this mode is not implimented')
              do ix = nw_i,nwx
                zwz(ix,:,:)=zwz(ix,:,:) - zwz00
              enddo
              deallocate(zwz00)
            endif

c-------------------------------
c loop over w and t in SEc(qt,w)
c-------------------------------
            if(debug) write(6,*)' sss ngb, nstate, ntp0=',ngb,nstate,ntp0
            if(debug) write(6,*)' sss zmel=',sum(abs(zmel(:,:,:)))

            if(verbose()>50) write(*,'("10 wfacx  iw,itp,it cycles ",$)')
            if(verbose()>50) call cputid(0)
            do 2001 iw  = iwini,iwend
              do 2002 itp = 1,ntq
                if(debug) write(6,*)'2011 0 zmel=',sum(abs(zmel(:,:,:)))
                omg = omega(itp,iw)
                if (omg >= ef) then
                  itini= nt0m+1
                  itend= nt_max
                  iii=  1
                else
                  itini= 1
                  itend= nt0p
                  iii= -1
                endif

                do 2011 it= itini,itend
                  if(debug) write(6,*)'2011 1 loop--- it=',iw,itp,it,sum(abs(zmel(:,:,:)))
                  esmrx = esmr
                  if(it<=nctot) esmrx = 0d0
                  wfac = wfacx2(omg,ef, ekc(it),esmrx)
                  if(GaussSmear) then
                    if(wfac<wfaccut) cycle
                    we = .5d0*abs(omg-weavx2(omg,ef, ekc(it),esmr))
                  else
                    if(wfac==0d0) cycle
                    if(omg>=ef) we = 0.5d0* abs(max(omg-ekc(it), 0d0)) ! positive
                    if(omg< ef) we = 0.5d0* abs(min(omg-ekc(it), 0d0)) ! negative
                  endif

                  wfac= iii* wfac*wtt
c three-point interpolation for Wc(we)
                  do iwp = 1,nw
                    ixs=iwp
                    if(freq_r(iwp)>we) exit
                  enddo
                  if(nw_i==0) then
                    if(ixs+1>nwx) then
                      write(6,*)' ixs,nwx, we =',ixs,nwx,we
                      call rx( ' sxcf: ixs+1>nwx xxx2')
                    endif   
                  else          !   write(6,*)" ixs nwxi=",ixs,nwxi,freq_r(ixs-1),we,freq_r(ixs)
                    if(omg >=ef .and. ixs+1> nwx ) then
                      write(6,*)'ixs+1 nwx=',ixs+1,nwx
                      call rx( ' sxcf: ixs+1>nwx yyy2a')
                    endif
                    if(omg < ef .and. abs(ixs+1)> abs(nwxi) ) then
                      write(6,*)'ixs+1 nwxi=',ixs+1,nwxi
                      call rx( ' sxcf: ixs-1<nwi yyy2b')
                    endif
                  endif

                  iir=1
                  if(omg < ef .and. nw_i/=0) iir = -1 !May2006 because of \int d omega' G(omega-omega') W(omega')
                  if(zwz3mode) then
                    zwz3=(0d0,0d0)
                    if(debug) write(6,"('wwwwwww ixs=',10i4)"),ixs,igb2,it,itp
                    if(debug) write(6,*)'2011 www zmel aaa=',sum(abs(zmel(:,:,:)))
                    do ix = ixs, ixs+2
                      do igb2=1,ngb
                        zz2 = sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,iir*(ix-1)) )
                        zwz3(ix-ixs+1) = zwz3(ix-ixs+1)+zz2 *zmel(igb2,it,itp)
                      enddo     !igb2
                    enddo       !ix
                    if(debug) write(6,"('w xxxxxxxxxxxxx ixs loopend=',i4)"),ixs
                    if(debug) write(6,*)zwz3(1:3) !,freq_r(ixs-1),zwz3(1:3)
                    if(debug) write(6,*)'we frez zwz3=', we,ixs,freq_r(ixs-1:ixs+1)
                    if(debug) write(6,*)'2011 bbb www zmel=',sum(abs(zmel(:,:,:)))

                    zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     &               + wfac *alagr3zz(we,freq_r(ixs-1),zwz3) !faleev

                    if(debug) write(6,*)'2011 ccc www zmel=',sum(abs(zmel(:,:,:)))
                    if(debug) write(6,"('wwwwwww eo zsecsum')")
                  else
                    zwzz(1:3) = zwz(iir*(ixs-1):iir*(ixs+1):iir, it,itp)
                    zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     &               + wfac*alagr3zz(we,freq_r(ixs-1),zwzz)
                  endif
 2011           continue
 2002         continue          !end of SEc w and qt -loop
 2001       continue            !end of SEc w and qt -loop
            if(debug) write(6,*)' end of do 2001'
            if(verbose()>50) then
              write(*,'("11 after alagr3zz iw,itp,it cycles ",$)')
              call cputid(0)
            endif
            if(debug) then
              do itp = 1,ntq
                write(6,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
              enddo
            endif
            if(zwz3mode) then
              deallocate(zmel,zw3)
            else
              deallocate(zwz)
            endif
 1000     continue
c      if(newaniso) ifvcoud =iclose('Vcoud.'//charnum5(kx))
          ifvcoud =iclose('Vcoud.'//charnum5(kx))
          if(.not.exchange) then
            ifrcw  = iclose('WVR.'//charnum5(kx))
            ifrcwi = iclose('WVI.'//charnum5(kx))
          endif
 1100   continue                ! end of k-loop
        if (tote) then
          exx = exx + wik(ip) * exxq * 0.25d0
        endif
        if (allocated(zz)) deallocate(zz)
        if (allocated(zmel)) deallocate(zmel)
        if (allocated(zzmel))deallocate(zzmel)
        if (allocated(zw)) deallocate(zw)
        if (allocated(zwz)) deallocate(zwz)
        if (allocated(zwz0)) deallocate(zwz0)
        if (allocated(zwzi)) deallocate(zwzi)
        if (allocated(zwz00)) deallocate(zwz00)
        if (allocated(w1p)) deallocate(w1p)
        if (allocated(w2p)) deallocate(w2p)
        if (allocated(w3p)) deallocate(w3p)
        if (allocated(z1p)) deallocate(w1p)
        if (allocated(vcoul)) deallocate(vcoul)
        if (allocated(vcoult)) deallocate(vcoul)
c      if (allocated(zmel1)) deallocate(zmel1)
        if (allocated(zmel3)) deallocate(zmel3)
        if (allocated(zw_)) deallocate(zw_)
        if (allocated(zwz2)) deallocate(zwz2)
c      if (allocated(zw2)) deallocate(zw2)
        if (allocated(zmel2)) deallocate(zmel2)
        if (allocated(zw3)) deallocate(zw3)
        if (allocated(uaa)) deallocate(uaa)
 1001 continue
c      if (allocated(expikt)) deallocate(expikt)
      end
