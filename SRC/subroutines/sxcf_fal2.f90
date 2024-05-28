subroutine sxcf_fal3z(&
     kount,ixc,deltaw,shtw,qip,itq, ntq,ef,ef2,esmr,esmr2,&
     nsp,isp,  qbas,ginv,qibz,qbz,wk,nstbz,wik,  &
     nstar,irkip,  freq_r,freqx,wx, &
     dwdummy,ecore, nlmto,nqibz,nqbz,nctot,&
     nbloch,ngrp, nw_i,nw ,niw,niwx,nq, &
     nblochpmx ,ngpmx,ngcmx, &
     wgt0,nq0i,q0i,symgg,alat, nband, ifvcfpout, &
     exchange,tote,screen,cohtest, ifexsp,&
     iwini,iwend, nbmx,ebmx, wklm,lxklm, dwplot,&
     zsec,coh,exx)
  use m_readqg,only:readqg0
  use m_readeigen,only: readeval
  use m_keyvalue,only: getkeyvalue
  use m_zmel,only: get_zmel_init,setppovlz, zmel
  use m_readVcoud,only:   Readvcoud, vcoud,vcousq,zcousq,ngb,ngc
  use m_wfac,only:wfacx2,weavx2
  implicit none
  intent(in)::&
       kount,ixc,deltaw,shtw,qip,itq, ntq,ef,ef2,esmr,esmr2, nsp,isp,  &
       qbas,ginv,   qibz,qbz,wk,nstbz,wik,  nstar,irkip,  freq_r,freqx,wx, &
       dwdummy,ecore, nlmto,nqibz,nqbz,nctot, nbloch,ngrp, nw_i,nw ,niw,niwx,nq, &
       nblochpmx ,ngpmx,ngcmx, wgt0,nq0i,q0i,symgg,alat, nband, ifvcfpout,&
       exchange,tote,screen,cohtest, ifexsp, iwini,iwend,&
       nbmx,ebmx, wklm,lxklm 
  !! TimeReversal off. when nw_i is not zero.
  !! Calcualte diagonal part only version of simga_ii(e_i)= <i|Re[S](e)|i> 
  !! Similar with sxcf_fal2.sc.F
  !o zsec: S_ij= <i|Re[S](e)|i> where e=e_i and e_i \pm deltaw
  !o
  !r  exchange=T : Calculate the exchange self-energy
  !r          =F : Calculate correlated part of the self-energy
  !r
  !r
  !r---- 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F
  !r 2000 takao kotani. This sxcf is starting from sec.f F.Aryasetiawan.
  !---------------------------------------------------------------


  !---- original document for sce.f (correlation case) by ferdi.Aryasetiawan.
  ! 92.02.24
  ! 93.10.18 from sec.f modified to take into account equivalent atoms
  ! calculates the correlated part of the self-energy SE
  ! SEc(q,t,t') = <psi(q,t) |SEc| psi(q,t'>
  ! SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >

  ! the zeroth order Green function
  ! G(r,r';w)   = S[occ]   psi(kn,r) psi(kn,r')^* /(w-e(kn)-i*delta)
  !             + S[unocc] psi(kn,r) psi(kn,r')^* /(w-e(kn)+i*delta)

  ! the screened coulomb potential
  ! Wc(r,r';w)  = W(r,r';w) - v(|r-r'|)
  !             = < [r1,r2] v(|r-r1|) X(r1,r2;w) v(|r2-r'|) >
  ! W(r,r';w)   = < [r''] ei(r,r'';w) v(|r''-r'| >
  ! ei          = e^(-1), inverse dielectric matrix
  !             = 1 + vX
  ! e           = 1 - vX0 in RPA

  ! expand Wc(r,r';w) in optimal product basis B
  ! Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch]
  !               B(k,i,r) Wc(k,w)(i,j) B(k,j,r')^*
  ! Wc(k,w)(i,j) are  the matrix elements of Wc in B

  ! SEc(q,t,t') = S[k=FBZ] S[n=occ]   S[i,j=1,nbloch]
  !        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
  !        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)-i*delta)>
  !
  !             + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
  !        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
  !        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)+i*delta)>

  ! the analytic structure of GWc for w .le. ef
  !                               |
  !                               |   o = pole of G
  !                               ^   x = pole of Wc
  !                               |
  !                               |   ef-w
  !                               |----<-----
  !                               |          |
  !                 o  o  o  o  o |o  o  o   ^
  !               x  x  x  x  x  x|          |
  !  -----------------------------|---->------------------------------
  !                               |x  x  x  x  x  x  x  x
  !                               |              o  o  o  o  o
  !                               |       <----->
  !                               ^        gap in insulator
  !                               |
  !                               |

  ! the analytic structure of GWc for w .gt. ef
  !                               |
  !                               |   o = pole of G
  !                               |   x = pole of Wc
  !                               |
  !         gap in insulator      ^
  !                <----->        |
  !      o  o  o  o               |
  !         x  x  x  x  x  x  x  x|
  !  ------------------------>----|-----------------------------------
  !                   |           |x  x  x  x  x  x  x  x
  !                   ^   o  o  o  o  o  o  o
  !                   |           |
  !                    ------<----|
  !                       w-ef    |
  !                               ^
  !                               |
  ! integration along the real axis from -inf to inf is equivalent to
  ! the integration along the path shown
  !------------------------------------------------------------
  ! integration along the imaginary axis: wint (s. also wint.f) (takao ->wintz)
  !   (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
  ! the i*delta becomes irrelevant
  !------------------------------------------------------------
  !
  ! omit k and basis index for simplicity and denote e(q-k,n) = e
  ! wint = (i/2pi) < [w'=-inf,inf] Wc(w')/(w+w'-e) >
  !
  ! w' ==> iw', w' is now real
  ! wint = - (1/pi) < [w'=0,inf] Wc(iw') (w-e)/{(w-e)^2 + w'^2} >
  !
  ! transform: x = 1/(1+w')
  ! this leads to a denser mesh in w' around 0 for equal mesh x
  ! which is desirable since Wc and the lorentzian are peaked around w'=0
  ! wint = - (1/pi) < [x=0,1] Wc(iw') (w-e)x^2/{(w-e)^2 + w'^2} >
  !
  ! the integrand is peaked around w'=0 or x=1 when w=e
  ! to handel the problem, add and substract the singular part as follows:
  ! wint = - (1/pi) < [x=0,1] { Wc(iw') - Wc(0)exp(-a^2 w'^2) }
  !                          * (w-e)/{(w-e)^2 +w'^2}x^2 >
  !        - (1/2) Wc(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
  !
  ! the second term of the integral can be done analytically, which
  ! results in the last term
  ! a is some constant
  !
  ! when w = e, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w') and
  ! the integral becomes -Wc(0)/2
  ! this together with the contribution from the pole of G (s.u.)
  ! gives the so called static screened exchange -Wc(0)

  !--------------------------------------------
  ! contribution from the poles of G: SEc(pole)
  !--------------------------------------------
  !
  ! for w .le. ef
  ! SEc(pole) = - S[k=FBZ] S[n=occ] S[i,j=1,nbloch]
  !        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
  !             Wc(k,e(q-k,n)-w)(i,j) theta(e(q-k,n)-w)
  !
  ! for w .gt. ef
  ! SEc(pole) = + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
  !        <psi(q,t) |psi(q-k,n) B(k,i)> <B(k,j) psi(q-k,n) |psi(q,t')>
  !             Wc(k,w-e(q-k,n))(i,j) theta(w-e(q-k,n))
  !
  ! theta(x)  = 1   if x > 0
  !           = 1/2 if x = 0
  !           = 0   if x < 0

  ! FBZ = 1st BZ
  ! NOTE: the routine only calculates the diagonal elements of the SE
  !       i.e. SEc(q,t)

  ! q       = q-vector in SEc(q,t)
  ! itq     = states t at q
  ! ntq     = no. states t
  ! eq      = eigenvalues at q
  ! ef      = fermi level in Rydberg
  ! tr      = translational vectors in rot*R = R' + T
  ! iatomp(R) = R'
  ! ifrw,ifcw,ifrwi,ifcwi
  !   = direct access unit files for Re and Im coulomb matrix
  !     along real and imaginary axis
  ! ifrb,ifcb,ifrhb,ifchb
  !         = direct access unit files for Re and Im b,hb
  ! qbas    = base reciprocal lattice vectors
  ! ginv    = inverse of qbas s. indxrk.f
  !xxxxx ippb,ipdb,idpb,iddb = pointers to work array w for
  !  ppb     = <phi(RLn) phi(RL'n') B(R,i)>
  !  pdb     = <phi(RLn) phidot(RL'n') B(R,i)>
  !  dpb     = <phidot(RLn) phi(RL'n') B(R,i)>
  !  ddb     = <phidot(RLn) phidot(RL'n') B(R,i)>
  ! freq    = frequencies along real axis
  ! freqx   = gaussian frequencies x between (0,1)
  ! freqw   = (1-freqx)/freqx
  ! wx      = weights at gaussian points x between (0,1)
  ! ua      = constant in exp(-ua^2 w'^2) s. wint.f
  ! expa    = exp(-ua^2 w'^2) s. wint.f
  ! dw      = frequency mesh along real axis
  ! deltaw  = energy mesh in SEc(qt,w) ---Not used now
  ! iclass  = given an atom, tells the class
  ! wk      = weight for each k-point in the FBZ
  ! indexk  = k-point index
  ! qbz     = k-points in the 1st BZ
  ! nstar   = no. stars for each k
  ! irk(k,R,nq) = gives index in the FBZ with k{IBZ, R=rotation
  ! mdim    = dimension of B(R,i) for each atom R
  ! work arrays:
  ! rbq,cbq     = real and imaginary part of b(q)
  ! rhbq,chbq   = real and imaginary part of hb(q)
  ! rbkq,cbkq   = real and imaginary part of b(q-k)
  ! rhbkq,chbkq = real and imaginary part of hb(q-k)
  !   b is the eigenvector of the LMTO-Hamiltonian
  ! ekq     = eigenvalues at q-k
  ! rmel,cmel = real and imaginary part of
  !             <psi(q,t') | psi(q-k,t) B(k,R,i)>
  ! wr1 ... = work arrays
  ! dimensions:
  ! nqibz   = number of k-points in the irreducible BZ
  ! n1,n2,n3= divisions along base reciprocal lattice vectors
  ! natom   = number of atoms
  ! nctot   = no. allowed core states
  ! nbloch  = total number of Bloch basis functions
  ! nlnmx   = maximum number of l,n,m
  ! nlmto   = total number of LMTO basis functions
  ! ngrp    = no. group elements (rotation matrices)
  ! niw     = no. frequencies along the imaginary axis
  ! nw      = no. frequencies along the real axis
  ! niwx    = max(niw,nw)
  !
  ! secq(t) = <psi(q,t) |SEc| psi(q,t)>
  !----------------------------------------------------------------------
  integer :: ntq, nqbz,nqibz,ngrp,nq,nw,niw, nband,  &
       nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,&
       nstatex, irot,  iqisp,ikpisp,isp,nsp, ip, it,itp,&
       ifrcw,ifrcwi,          ifvcfpout,ndummy1,&
       ndummy2,kx,kr,nbloch, kp,nt0,nocc, nt0p,nt0m,irkp,i,&
       nt0org,nmax,nt,ntp0, nbmax,nblochpmx,ix,nx,iw,&
       iwp,ixs,ixsmx, nwx,niwx, itq(ntq),  &
       nstar(nqibz),irkip(nqibz,ngrp,nq),kount(nqibz,nq)
  real(8) :: q(3),qbas(3*3),ginv(3*3), wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),&
       &freqx(niw),wx(niw),     eq(nband,nq), ekq(nband), ekc(nctot+nband), &
       &tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx, dwdummy,&
       &wtt,wexx,www,exx,exxq,wex
  integer:: ngpmx, ngcmx, igc,                     nadd(3)
  real(8) :: wgt0(nq0i,ngrp),qk(3), qdiff(3),add(3),symgg(3,3,ngrp),symope(3,3),&
       & qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot), coh(ntq,nq)
  complex(8)::   alagr3zz,wintz
  complex(8),allocatable :: zz(:),zzmel(:,:,:),&
       & zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:)
  ! for exchange --------------------
  logical :: exchange,screen,cohtest,tote
  real(8),allocatable:: w1p(:,:,:),w2p(:,:,:),w3p(:,:)
  complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
  integer:: invrot,invr
  logical :: debug=.false. ,onceww
  complex(8) :: wintzsg_npm
  integer :: ibl,iii,ivsumxxx,ifexsp 
  integer,save::ifzwz=-999
  integer :: iwini, iwend, ia
  real(8)    :: esec, omega(ntq, iwini:iwend)
  complex(8) :: zsec(iwini:iwend,ntq,nq)
  complex(8):: img=(0d0,1d0)
  integer :: nt_max, igb1,igb2,iigb,  nw_i !nw_i is at feb2006 TimeReversal off case
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
  integer:: nbmx,nbmxe,nstatetot
  integer::verbose,nstbz(nqbz),bzcase=1,iqini,iqend
  real(8):: wgtq0p
  integer:: nrec,kxx
  real(8)::quu(3),qibz_k(3),qbz_kr(3)
  logical ::  onlyimagaxis 
  logical ::zwz3mode
  real(8):: ua_,expa_(niw),ua2,freqw,freqw1,ratio,ua2_(niw)
  integer:: icc=0
  real(8),allocatable:: uaa(:,:)
  !cccc zvz test cccccccccccccccccccccccccc
  integer:: ngbx
  !      complex(8):: vcoul(ngbx,ngbx)
  complex(8),allocatable:: vzz(:,:,:),aaa(:), zwzs(:)
  complex(8):: zvz,zvz1
  integer:: ib1,ib2,ifix
  !cccccccccccccccccccccccccccccccccc
  logical ::iww2=.true., oncew
  !...
  !      logical::smbasis
  integer:: iclose,isx,iqx !nn,no,ifpomat,
  !      complex(8),allocatable:: pomat(:,:)
  real(8):: q_r(3)
  !      integer:: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
  !      complex(8):: pomatr(nnmx,nomx,nkpo)
  !      real(8):: qrr(3,nkpo)

  real(8):: elxx,ehxx,ekxx,efxx
  integer:: ixsmin,iwm,iir,nwxi, itini,itend, npm
  real(8)   :: fffr(3),ppp
  complex(8):: zwzz(3)

  real(8),allocatable:: ebb(:)
  integer:: ii,iq
  logical ::evaltest        !, imgonly

  integer:: lxklm,ivc,ifvcoud,idummy,iy,ngb0
  real(8):: wklm((lxklm+1)**2),pi,fpi,vc,qvv(3),aaaa
  complex(8)::zmelt1,zmelt0
  real(8)::voltot
  !      logical :: newaniso !fixed to be T

  !      complex(8),allocatable:: ppovl(:,:)!,zcousq(:,:) !,ppovlz(:,:)
  !      real(8),allocatable::vcoud(:)!,vcousq(:)
  integer:: mrecl,nprecx,ifwd
  character(10):: i2char

  integer:: ixc
  real(8):: qip(3,*),deltaw,shtw,eqx(nband),dwplot,tolq=1d-8
  complex(8),allocatable:: zmelt(:,:)
  integer:: ntqxx,nrot
  !--------------------------------------------------------------------
  write(6,*)'sxcf_fal3z'
  !      timemix=.false.
  pi  = 4d0*datan(1d0)
  fpi = 4d0*pi
  debug=.false.
  if(verbose()>=90) debug=.true.
  !!
  if(.not.exchange) then
     open(newunit=ifwd,file='WV.d')
     read (ifwd,*) nprecx,mrecl
     close(ifwd)
     !$$$!! --- gauss_img : interpolation gaussion for W(i \omega).
     !$$$      call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
     !$$$      if(ua_<=0d0) then
     !$$$        ua_auto =.true.
     !$$$        write(6,"(' ua_auto=T')")
     !$$$      else
     !$$$        ua_auto =.false.
     !$$$        do ix = 1,niw
     !$$$          freqw     = (1d0 - freqx(ix))/ freqx(ix)
     !$$$          expa_(ix) = exp(-(ua_*freqw)**2)
     !$$$        enddo
     !$$$      endif
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
  !      call dinv33(qbas,0,qbasinv,det)
  !      allocate(expikt(natom))


  !! == ip loop to spedify external q ==
  do 1001 ip = 1,nq   
     if(sum(irkip(:,:,ip))==0) cycle
     q = qip(1:3,ip) 
     write (*,*) ip,'  out of ',nq,'  k-points ' ! call cputid  (0)
     if(ixc==2) then
        eqx= readeval(q,isp)
        do iw = iwini,iwend
           do i  = 1,ntq
              omega(i,iw) = eqx(itq(i)) + 2d0*(dble(iw)-shtw)*deltaw
           enddo
        enddo
     endif
     !!
     if(ixc==4) then
        !        dwplot=0.01
        do iw = iwini,iwend
           omega(1:ntq,iw) =  dwplot* iw + ef
        enddo
     endif

     eq(:,ip)= readeval(q, isp) 
     !! we only consider bzcase()==1
     if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
     if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
     If (tote) exxq = 0.d0

     !! == Big loop for kx ==
     !! kx is for irreducible k points, kr=irk(kx,irot) runs all k points in the full BZ.
     iqini=1
     iqend=nqibz             !no sum for offset-Gamma points.
     kxloop: do 1100 kx = iqini,iqend
        if(sum(irkip(kx,:,ip))==0) cycle
        write(6,*) ' ### do 1100 start kx=',kx,' from ',iqini,' through', iqend
        !          if( kx <= nqibz ) then
        qibz_k= qibz(:,kx)
        !          else
        !            qibz_k= 0d0
        !         endif
        !          if(verbose()>=40)  write(6,*) ' sxcf_fal3z: loop 1100 kx=',kx
        !          call readqg0('QGcou',qibz_k,  quu,ngc)
        !          ngb = nbloch + ngc    !oct2005
        !          write(6,*) ' xxxxxxxxxxxx sxcf: ngb=',ngb,nbloch

        !! ===Readin diagonalized Coulomb interaction===
        !!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
        !!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
        !!  Vould file is written in hvccfp0.m.F.
        !! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
        !! These can be also separeted into WVR.ID and WVI.ID files.
        !! NOTE: vcoud and zcousq are in module m_zmelt.
        !          if(kx<=nqibz) qxx=qibz_k
        !          if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
        !          qxx=qibz_k
        !          write(6,*) ' xxxxxxxxxxxx2222222222 sxcf: ngb=',ngb,nbloch
        !          ifvcoud = iopen('Vcoud.'//i2char(kx),0,0,0)
        !          do
        !            read(ifvcoud) ngb0
        !            read(ifvcoud) qvv
        !            if(allocated(vcoud)) deallocate(vcoud)
        !            if(allocated(zcousq)) deallocate(zcousq)
        !            allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
        !            read(ifvcoud) vcoud
        !            read(ifvcoud) zcousq
        !            if(sum(abs(qvv-qxx))<tolq) goto 1133
        !          enddo
        !          if(sum(abs(qvv-qxx))>tolq) then
        !            write(6,*)'qvv =',qvv
        !            write(6,*)'qxx=',qxx,kx
        !            call rx( 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent')
        !          endif
        ! 1133     continue
        !          if( ngb0/=ngb ) then  !sanity check
        !            write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
        !            call rx( 'hsfp0.m.f:ngb0/=ngb')
        !          endif
        !! used in get_zmel
        !! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
        !! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi phi|I> O^-1_IJ
        !! ppovlz= O Zcousq
        !! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
        !$$$          if(allocated(ppovl)) deallocate(ppovl,ppovlz)
        !$$$          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
        !$$$          call readppovl0(qibz_k,ngc,ppovl)
        !$$$          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
        !$$$          ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))

        call Readvcoud(qibz_k,kx,NoVcou=.false.) ! Readin vcousq,zcousq !for the Coulomb matrix
        call Setppovlz(qibz_k,matz=.true.,npr=ngb) !add npr at 2024-5-23 obata
        if(debug) write(6,*) ' sxcf_fal2: ngb ngc nbloch=',ngb,ngc,nbloch

        !! === open WVR,WVI for correlation mode ===
        if(.not.exchange) then
           open(newunit=ifrcw, file='WVR.'//i2char(kx),form='unformatted',access='direct',recl=mrecl)
           open(newunit=ifrcwi,file='WVI.'//i2char(kx),form='unformatted',access='direct',recl=mrecl)
        endif
        nrot=0
        do irot = 1,ngrp
           !            if( kx <= nqibz) then
           kr = irkip(kx,irot,ip) ! index for rotated kr in the FBZ
           if(kr==0) cycle   ! next irot
           qbz_kr= qbz (:,kr) 
           !            else
           !              kr=-99999         !for sanity check
           !              qbz_kr= 0d0
           !              if( wgt0(kx-nqibz,irot)==0d0 ) cycle ! next irot
           !            endif
           nrot=nrot+1
        enddo

        !! == loop over rotations ==
        !! We may extend 
        do 1000 irot = 1,ngrp
           !            if( kx <= nqibz) then
           kr = irkip(kx,irot,ip) ! index for rotated k in the FBZ
           if(kr==0) cycle
           qbz_kr= qbz (:,kr) 
           !            else
           !              kr=-99999         !for sanity check
           !              qbz_kr= 0d0
           !              if( wgt0(kx-nqibz,irot)==0d0 ) cycle
           !            endif
           write(*,"('ip,kx irot=',3i5, ' out of',2i4)") ip,kx,irot, iqend,ngrp

           ! qk = q - rk, rk is inside 1st BZ, not restricted to the irreducible BZ
           qk =  q - qbz_kr    ! qbz(:,kr)
           ekq = readeval(qk, isp)
           ekc(nctot+1:nctot+nband) = ekq (1:nband)
           nt0 = count(ekc<ef) !nocc(ekc,ef,nstatetot)!
           ddw= .5d0
           !        if(GaussSmear) ddw= 10d0
           ddw=10d0
           efp= ef+ddw*esmr
           efm= ef-ddw*esmr
           nt0p = count(ekc<efp) !nocc(ekc,efp,nstatetot)!
           nt0m = count(ekc<efm) !nocc(ekc,efm,nstatetot)!
           !! nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
           !! nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W
           if(exchange) then
              nbmax = nt0p-nctot
              if(debug) write(6,*)' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
           else
              nbmax = nband
              nbmxe = count(ekc<ebmx)-nctot !nocc(ekc,ebmx,nstatetot)-nctot!
              nbmax  = min(nband,nbmx,nbmxe)
              if(onceww(3)) write(6,*)' nbmax=',nbmax
           endif
           nstate = nctot + nbmax ! = nstate for the case of correlation
           !! all are identical.
           ntp0 = ntq
           ntqxx= ntp0
           !! Get matrix element zmelt= rmelt + img*cmelt, defined in m_zmel.F---
           call Get_zmel_init(q,qibz_k,irot,qbz_kr, 1,nbmax+nctot,isp, 1,ntqxx,isp, nctot,ncc=0,iprx=.false.,zmelconjg=.false.)
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
           if(exchange) then
              allocate( w3p( nctot+nbmax,ntp0))
              do 992 itp = 1,ntp0
                 do 993 it  = 1,nctot+nbmax
                    w3p(it,itp) = 0d0
                    do 994 ivc=1,ngb
                       if(ivc==1.and.kx==1) then
                          vc= wklm(1)* fpi*sqrt(fpi) /wk(kx)
                       else
                          vc= vcoud(ivc)
                       endif
                       w3p(it,itp) = w3p(it,itp)+ vc * abs(zmel(ivc,it,itp))**2
994                 enddo!continue
993              enddo!continue
992           enddo!continue
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
                       write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")&
                            ip,itp,it, qbz_kr, ekc(it), -wtt*w3p(it,itp)
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
                    zsec(iwini,itp,ip) = zsec(iwini,itp,ip) &
                         - wtt * sum( w3p(:,itp) )
                 enddo
              else
                 do itp = 1,ntp0
                    wfac = wfacx(-1d99, ef2, eq(itq(itp),ip), esmr2) !gaussian
                    w3p(1:nctot+nbmax,itp) = wfac * w3p(1:nctot+nbmax,itp)
                    exxq = exxq - wtt * sum( w3p(:,itp) )
                 enddo
              endif
              if(debug) then
                 do itp = 1,ntp0
                    write(6,'(" exchange zsec=",i3,6d15.7)') itp,zsec(iwini,itp,ip)
                 enddo
              endif

              deallocate( w3p)  !,rmelt,cmelt)
              cycle
           endif
           !-- End of exchange section --------------



           !--------------------------------------------------------------------------
           !--- correlation section --------------------------------------------------
           !--------------------------------------------------------------------------
           !$$$c--- The matrix elements zmel.
           !$$$c        allocate( zmel (ngb, nstate, ntp0) )
           !$$$c        zmel = dcmplx (rmelt,-cmelt)
           !$$$c        if(newaniso) then
           !$$$c#ifdef USE_GEMM_FOR_SUM
           !$$$          if(verbose()>39)write(*,*)'info: USE GEMM FOR SUM (zmel=zmel*ppovlz), in sxcf_fal2.F'
           !$$$          allocate( zmelt (ngb, nstate) )
           !$$$          do itp=1,ntp0
           !$$$          zmelt = dcmplx(rmelt(:,:,itp),-cmelt(:,:,itp))
           !$$$          call zgemm('C','N',ngb,nstate,ngb,(1d0,0d0),
           !$$$     .      ppovlz,ngb,zmelt,ngb,(0d0,0d0),zmel(1,1,itp),ngb)
           !$$$          enddo
           !$$$          deallocate(zmelt)
           !$$$#else
           !$$$          do itp=1,ntp0
           !$$$            do it=1,nstate
           !$$$              zmel(:,it,itp) =  matmul(zmel(:,it,itp),dconjg(ppovlz(:,:)))
           !$$$            enddo
           !$$$          enddo
           !$$$#endif
           !$$$c        endif
           !        deallocate(rmelt,cmelt)
           !        if(debug) write(6,*)' end of zmel'

           !================================================================
           ! The correlated part of the self-energy:
           ! S[n=all] S[i,j=1,nbloch]
           ! <psi(q,t) |psi(q-rk,n) B(rk,i)>
           !  < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
           !                                <B(rk,j) psi(q-rk,n) |psi(q,t)>
           ! e = e(q-rk,n), w' is real, Wc = W-v
           !================================================================
           allocate( zw (nblochpmx,nblochpmx) )
           !====================================================================
           ! contribution to SEc(qt,w) from integration along the imaginary axis
           !====================================================================
           !------------------------------------------------
           ! loop over w' = (1-x)/x, frequencies in Wc(k,w')
           ! {x} are gaussian points between (0,1)
           !------------------------------------------------
           allocate( zwz0(nstate,ntp0) )
           ix = 1  - nw_i      !at omega=0
           !        nrec=(kx-iqini)*(nw-nw_i+1) +ix ! 2---> iqini
           !        if(bzcase()==2) nrec= (kx-1)*(nw-nw_i+1) +ix
           nrec=ix 
           if(debug) write(6,*)' wvr nrec kx nw nw_i ix=',nrec,kx,nw,nw_i,ix
           read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v
           zwz0=0d0
           !! this loop looks complicated but just in order to get zwz0=zmel*zwz0*zmel 
           !! Is this really efficient???
           do itp=1,ntp0
              do it=1,nstate
                 do igb2=2,ngb
                    zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )
                    zwz0(it,itp) = zwz0(it,itp)+zz2*zmel(igb2,it,itp)*2d0+&
                         &             dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp)
                 enddo           !igb2
                 zwz0(it,itp) = zwz0(it,itp)+     dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
              enddo             !it
           enddo               !itp
           zwz0 = dreal(zwz0)
           ! COH term test ----- The sum of the all states for zwz00 gives the delta function.
           if(cohtest) then
              do itp = 1,ntq
                 coh(itp,ip)  = coh(itp,ip)  + .5d0*wtt*sum(dreal(zwz0(1:nstate,itp)))
              enddo
              deallocate(zw,zwz0)
              cycle
           endif
           nx  = niw
           if(niw <1) call rx( " sxcf:niw <1")
           if(allocated(zwz)) deallocate(zwz)
           if(allocated(zwzi)) deallocate(zwzi)
           allocate( zwz(niw*npm, nstate,ntp0),  zwzi(nstate,ntp0) )
           if(screen) allocate( zwz00(nstate,ntp0) )
           if(verbose()>50) write(*,'("6 before matzwz in ix cycle ",$)')
           if(verbose()>50) call cputid(0)
           zwz=0d0       
           do ix = 1,nx          ! imaginary frequency w'-loop
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
                          ppp = ppp + dreal(zz2*zmel(igb2,it,itp)) * 2d0&
                               &                 + dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp)
                       enddo       !igb2
                       zwz(ix,it,itp) = ppp +      dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
                    enddo         !it
                 enddo           !itp
              else              !we need to use mode2 because zwz is not real now.
                 call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb, &
                      zwz(ix,1:nstate,1:ntp0))
              endif
              if(debug) write(6,*)' sumzw=',sum(abs(zw))
           enddo               !ix
           if(verbose()>50) write(*,'("xxx:6.1 before matzwz in ix cycle ",$)')
           if(verbose()>50) call cputid(0)
           if(debug) write(6,*)' sumzmel=',ngb, nstate, ntp0,sum(abs(real(zmel))),sum(abs(imag(zmel)))
           !--------------------------------------------------------------
           ! S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
           !                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
           ! needed to take care of the singularity in the w' integration
           ! when w-e(q-rk,n) is small
           !--------------------------------------------------------------
           if(screen) then
              zwz00 = zwz0
              zwz0  = 0d0
              do ix = 1,nx
                 zwz(ix,:,:)=zwz(ix,:,:) - zwz00
              enddo
           endif
           !------------------------------------------------
           ! loop over w in SEc(qt,w)
           !------------------------------------------------
           !$$$        if(ua_auto) then
           !$$$          allocate(uaa(nstate,ntq))
           !$$$          do itp = 1,ntq
           !$$$            do  it = 1,nstate
           !$$$              ratio = abs(zwz(niw,it,itp)/zwz0(it,itp))
           !$$$              call gen_uaa(ratio,freqx(niw),  uaa(it,itp))
           !$$$              if(verbose()>45) then
           !$$$                write(6,"(' it itp uaa=',2i4,12f8.4)") it,itp,uaa(it,itp)
           !$$$              elseif(verbose()>40.and.mod(it,10)==1.and.mod(itp,10)==1) then
           !$$$                write(6,"(' it itp uaa=', 2i4,12f8.4)") it,itp,uaa(it,itp)
           !$$$              endif
           !$$$            enddo
           !$$$          enddo
           !$$$        endif
           allocate(zwzs(npm*nx))
           do iw = iwini,iwend
              ! frequency integration along the imaginary axis, s. wint.f
              ! for each e(q-rk,n) and w in SEc(qt,w)
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
                    !          if(ua_auto) then
                    !            call gen_ua(abs(zwz(niw,it,itp)/zwz0(it,itp)), niw,freqx, expa_,ua_)
                    !            if(iw==ini) then
                    !            if(verbose()>45) then
                    !              write(6,"(' it itp ua_=',2i4,12f8.4)")it,itp,ua_
                    !            elseif(verbose()>40.and.mod(it,20)==1.and.mod(itp,20)==1) then
                    !              write(6,"(' it itp ua_=',3i4,12f8.4)")it,itp,ua_
                    !            elseif(irot==1.and.mod(it,10)==1.and.itp==it) then
                    !              write(6,"(' it itp ua_=',3i4,12f8.4)")it,itp,ua_
                    !            endif
                    !            endif
                    !          endif
                    !$$$              if(ua_auto) then
                    !$$$                ua_ = .5d0*uaa(it,itp)
                    !$$$                call gen_expa(niw,freqx,ua_,  expa_)
                    !$$$              endif
                    esmrx = esmr
                    if(it <= nctot) esmrx = 0d0
                    do ix=1,nx
                       zwzs(ix   ) = dreal( zwz(ix,it,itp)) ! w(iw) + w(-iw) symmetric part
                       if(npm==2) then
                          zwzs(ix+nx) = dimag( zwz(ix,it,itp)) ! w(iw) - w(-iw)
                       endif
                    enddo
                    !                  if(GaussSmear) then
                    zwzi(it,itp) =& 
                         &               wintzsg_npm(npm, zwzs,zwz0(it,itp),freqx,wx,ua_,expa_,we,nx, esmrx)
                    !                  else
                    !                    if(npm==2) 
                    !     &               call rx( ' ###Not impliment wintzav for npm=2. Use Gausssmear.')
                    !                    zwzi(it,itp) = 
                    !     &               wintzav( zwzs,zwz0(it,itp),freqx,wx,ua_,expa_,we,nx, esmrx)
                    !                  endif
                    !    .    wintz (zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx)
                    !cccccccccccccccccccccccccccccc
                    !          if(verbose()>45) then
                    !          if(it==50.and.itp==1) then
                    !          write(6,"(' it itp abs(zwzi)=',2i4,12d13.5)")it,itp,abs( zwzi(it,itp))
                    !          icc=icc+1
                    !          if(icc==10) stop 'test end'
                    !          endif
                    !          endif
                    !ccccccccccccccccccccccccccccc
1387             enddo!continue
1385          enddo!continue
              ! sum over both occupied and unoccupied states and multiply by weight
              do  itp = 1,ntq
                 zsec(iw,itp,ip)  = zsec(iw,itp,ip) + wtt*sum(zwzi(:,itp))
              enddo
              ! end of SEc w-loop
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
           !===============================================================================
           ! contribution to SEc(qt,w) from the poles of G
           !===============================================================================
           !    We assume freq_r(i) == -freq_r(-i) in this code. feb2006
           !---------------------------------------
           ! maximum ixs finder
           !---------------------------------------
           !      write(6,*)' ekc at nt0p nt0m+1=', ekc(nt0p),ekc(nt0m+1)      write(6,*)'  nt0p nt0m+1=', nt0p, nt0m+1
           ixsmx =0
           ixsmin=0
           iwloop3: do 3001 iw  = iwini,iwend
             itploop3: do 3002 itp = 1,ntq
               omg = omega(itp,iw)
               if (omg < ef) then
                 itini= 1
                 itend= nt0p
               else
                 itini= nt0m+1
                 itend= nstate
               endif
               itloop3: do 3011 it= itini,itend
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
                   !                write (*,*) 'xxx freq we=',freq_r(iwp),abs(we)
                   if(freq_r(iwp) > abs(we)) exit
                 enddo
                 ! This change is because G(omega-omg') W(omg') !may2006
                 !             if(ixs>ixsmx  .and. omg<=ef ) ixsmx  = ixs
                 !             if(ixs>ixsmin .and. omg> ef ) ixsmin = ixs
                 if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
                 if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
                 wexx  = we
                 if(ixs+1 > nw) then
                   write (*,*) ' nw_i ixsmin',nw_i, ixsmin
                   !                    write (*,*) ' wexx, dw ',wexx,dw
                   write (*,*) ' omg ekc(it) ef ', omg,ekc(it),ef
                   !stop2rx 2013.08.09 kino                stop ' sxcf 222: |w-e| out of range'
                   call rx( ' sxcf 222: |w-e| out of range')
                 endif
3011           enddo itloop3
3002         enddo itploop3          !end of SEc w and qt -loop
3001       enddo iwloop3           !end of SEc w and qt -loop
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
           !... Find nt_max ------------------------------------
           nt_max=nt0p         !initial nt_max
           do 4001 iw  = iwini,iwend
              do 4002 itp = 1,ntq
                 omg     = omega(itp,iw)
                 if (omg > ef) then
                    do  it = nt0m+1,nstate ! nt0m corresponds to efm
                       wfac = wfacx2 (ef,omg, ekc(it),esmr)
                       if( (GaussSmear.and.wfac>wfaccut) &
                            &               .or.(.not.GaussSmear.and.wfac/=0d0)) then
                          if (it > nt_max) nt_max=it ! nt_max is  unocc. state
                       endif       ! that ekc(it>nt_max)-omega > 0
                    enddo
                 endif
4002          enddo!continue
4001       enddo!continue
           !... Set zw3 or zwz -----------------------------------
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
                    call diagcvh2((zw(1:ngb,1:ngb)-transpose(dconjg(zw(1:ngb,1:ngb))))/2d0/img,& 
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
                 ! zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
                 call matzwz(zw(1:ngb,1:ngb), zmel(1:ngb,1:nstatex,1:ntp0), ntp0,nstatex,ngb,   &
                      zwz(ix,1:nstatex,1:ntp0))
                 ! zmel (ngb, nstate, ntp0)
              enddo
              !              deallocate(zmel)
              !call Deallocate_zmel()
              deallocate(zw)
           endif
           !---------------------------------------------
           if(screen) then
              if(zwz3mode) call rx( ' this mode is not implimented')
              do ix = nw_i,nwx
                 zwz(ix,:,:)=zwz(ix,:,:) - zwz00
              enddo
              deallocate(zwz00)
           endif
           !-------------------------------
           ! loop over w and t in SEc(qt,w)
           !-------------------------------
           if(debug) write(6,*)' sss ngb, nstate, ntp0=',ngb,nstate,ntp0
           if(debug) write(6,*)' sss zmel=',sum(abs(zmel(:,:,:)))
           if(verbose()>50) write(*,'("10 wfacx  iw,itp,it cycles ",$)')
           if(verbose()>50) call cputid(0)
           iwloop: do 2001 iw  = iwini,iwend
             itploop: do 2002 itp = 1,ntq
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
               itloop: do 2011 it= itini,itend
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
                 ! three-point interpolation for Wc(we)
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
                   if(debug) write(6,"('wwwwwww ixs=',10i4)")ixs,igb2,it,itp
                   if(debug) write(6,*)'2011 www zmel aaa=',sum(abs(zmel(:,:,:)))
                   do ix = ixs, ixs+2
                     do igb2=1,ngb
                       zz2 = sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,iir*(ix-1)) )
                       zwz3(ix-ixs+1) = zwz3(ix-ixs+1)+zz2 *zmel(igb2,it,itp)
                     enddo     !igb2
                   enddo       !ix
                   if(debug) write(6,"('w xxxxxxxxxxxxx ixs loopend=',i4)")ixs
                   if(debug) write(6,*)zwz3(1:3) !,freq_r(ixs-1),zwz3(1:3)
                   if(debug) write(6,*)'we frez zwz3=', we,ixs,freq_r(ixs-1:ixs+1)
                   if(debug) write(6,*)'2011 bbb www zmel=',sum(abs(zmel(:,:,:)))

                   zsec(iw,itp,ip) = zsec(iw,itp,ip) &
                        &               + wfac *alagr3zz(we,freq_r(ixs-1),zwz3) !faleev

                   if(debug) write(6,*)'2011 ccc www zmel=',sum(abs(zmel(:,:,:)))
                   if(debug) write(6,"('wwwwwww eo zsecsum')")
                 else
                   zwzz(1:3) = zwz(iir*(ixs-1):iir*(ixs+1):iir, it,itp)
                   zsec(iw,itp,ip) = zsec(iw,itp,ip) &
                        &               + wfac*alagr3zz(we,freq_r(ixs-1),zwzz)
                 endif
2011           enddo itloop
2002         enddo itploop      !end of SEc w and qt -loop
2001       enddo iwloop         !end of SEc w and qt -loop
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
              deallocate(zw3)
           else
              deallocate(zwz)
           endif
1000    enddo
        if(.not.exchange) then
           close(ifrcw) 
           close(ifrcwi)
        endif
1100 enddo kxloop
     if (tote) exx = exx + wik(ip) * exxq * 0.25d0
     if (allocated(zz)) deallocate(zz)
     if (allocated(zzmel)) deallocate(zzmel)
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
     if (allocated(zmel3)) deallocate(zmel3)
     if (allocated(zw_)) deallocate(zw_)
     if (allocated(zwz2)) deallocate(zwz2)
     if (allocated(zmel2)) deallocate(zmel2)
     if (allocated(zw3)) deallocate(zw3)
     if (allocated(uaa)) deallocate(uaa)
1001 enddo !continue
end subroutine sxcf_fal3z
