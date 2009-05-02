      subroutine sxcf_v2z(kount,ip,q,irot,itq, ntq,ef,ef2,esmr,esmr2,tr,
     i          iatomp,nsp,isp, ! ,ifcphi jan2004   ifrb,ifcb,ifrhb,ifchb,
     i                  ifrcw,ifrcwi, 
     i                  qbas,ginv,
     i         qibz,qbz,wk,nstbz,wik,nstar,irk,  !kount,,iindxk
c
c     i                  iiclass,mdim,inlnmv,inlnmc,
c     i                  iicore,incore,iimdim,
     i                  iclass,mdim,nlnmv,nlnmc,
     i                  icore,ncore,imdim,
     i                  ppb, !pdb,dpb,ddb,
     i                  freq,freqx,wx,expa,ua,dw, !deltaw,
     i                  ecore,
c
     d            nlmto,nqibz,nqbz,nctot,
c     i            index_qbz,n_index_qbz, !jan2004
     d                  nl,nnc,nclass,natom,
     d                  nlnmx,mdimx,nbloch,ngrp,nw,niw,niwx,nq,
c
c     &     nblochpmx ,ngpn,ngcni,ngpmx,ngcmx,geigB,ngvecpB,ngveccBr,
     &     nblochpmx ,ngpmx,ngcmx,!ngveccBr,!jan2004
     &     wgt0,nq0i,q0i,symope,alat, shtv,nband, ifvcfpout,
c     &     shtw,  
     &     exchange,tote,screen,cohtest, ifexsp,
ctetra
     &     wtetef,
     &     ibzx,tetraex, 
     &     ixc,
c
     i omega,iwini,iwend,
     i     nbmx,ebmx, !takao 18June2003 
     o zsec,coh,exx)
c 2002 May Takashi Miyake, Exx for total energy calc.
c 2002 Apr takao. Multiple number of argumentation waves in each MT.
c 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F
c 2000 takao kotani.
c This code is originally from sec.f by ferdi.Aryasetiawan.
c
c  exchange=T : Calculate the exchange self-energy
c  exchange=F : Calculate correlated part of the self-energy
c

c---- correlation case documents.  Originally by Ferdi and takao edited.
c Read GW manual.
c
c FBZ = 1st BZ
c NOTE: the routine only calculates the diagonal elements of the SE
c       i.e. SEc(q,t)
c
c calculates the correlated part of the self-energy SE
c SEc(q,t,t') = <psi(q,t) |SEc| psi(q,t'>
c SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
c
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
!
c expand Wc(r,r';w) in the Mixed  basis M
c Wc(r,r';w)  = S[k=FBZ] S[i,j=1,nbloch+ngc]
c               M(k,i,r) Wc(k,w)(i,j) M(k,j,r')^*
c Wc(k,w)(i,j) are  the matrix elements of Wc in M
!
c SEc(q,t,t') = S[k=FBZ] S[n=occ]   S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) M(k,i)> <M(k,j) psi(q-k,n) |psi(q,t')>
c        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)-i*delta)>
c
c             + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch]
c        <psi(q,t) |psi(q-k,n) M(k,i)> <M(k,j) psi(q-k,n) |psi(q,t')>
c        (i/2pi) <[w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n)+i*delta)>


c=========================================================================
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
c
c=========================================================================
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
c See the GW manual, which contains clear figure.
c
c integration along the real axis from -inf to inf is equivalent to
c the integration along the path shown
!
!
c------------------------------------------------------------
c integration along the imaginary axis: wint (s. also wint.f) (takao ->wintzav)
c   (i/2pi) < [w'=-inf,inf] Wc(k,w')(i,j)/(w'+w-e(q-k,n) >
c the i*delta becomes irrelevant
c------------------------------------------------------------
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
c
Cr!!!  takao : I modified this method so that all the pole has some width esmr (subroutine wintzav).
!
c--------------------------------------------
c contribution from the poles of G: SEc(pole)
c--------------------------------------------
c for w .le. ef
c SEc(pole) = - S[k=FBZ] S[n=occ] S[i,j=1,nbloch+ngc]
c        <psi(q,t) |psi(q-k,n) M(k,i)> <M(k,j) psi(q-k,n) |psi(q,t')>
c             Wc(k,e(q-k,n)-w)(i,j) theta(e(q-k,n)-w)
c
c for w .gt. ef
c SEc(pole) = + S[k=FBZ] S[n=unocc] S[i,j=1,nbloch+ngc]
c        <psi(q,t) |psi(q-k,n) M(k,i)> <M(k,j) psi(q-k,n) |psi(q,t')>
c             Wc(k,w-e(q-k,n))(i,j) theta(w-e(q-k,n))
c
c theta(x)  = 1   if x > 0
c           = 1/2 if x = 0
c           = 0   if x < 0
c takao: As I said, each pole is now assumed to have a width esmr.
!
!
c------------- kotani rewrited  at Apr 2002  ------------------------------
c q       = q-vector in SEc(q,t)
c itq     = states t at q
c ntq     = no. states t
c eq      = eigenvalues at q
c ef      = fermi level in Rydberg
c tr      = translational vectors in rot*R = R' + T
c iatomp(R) = R'
c ifvcfpout = file handle of Coulomb matrix
c
c qbas    = base reciprocal lattice vectors
c ginv    = inverse of qbas s. indxrk.f
c
c  ppb     = <phi(RLn) phi(RL'n') B(R,i)>
c
c freq    = frequencies along real axis
c freqx   = gaussian frequencies x between (0,1)
c freqw   = (1-freqx)/freqx
c wx      = weights at gaussian points x between (0,1)
c ua      = constant in exp(-ua^2 w'^2)  wint.f or wintzav
c expa    = exp(-ua^2 w'^2) s. wint.f
c dw      = frequency mesh along real axis
c iclass  = given an atom, tells the class
c wk      = weight for each k-point in the FBZ
c indexk  = k-point index
c qbz     = k-points in the 1st BZ
c nstar   = no. stars for each k
c irk(k,R) = gives index in the FBZ with k{IBZ, R=rotation
c mdim     = dimension of B(R,i) for each atom R
c work arrays:
c cphiq       = real and imaginary part of b(q)
c cphikq      = real and imaginary part of b(q-k)
c           b is the coefficients of eigenfunctions for argumentation wave in each MT.
c ekq     = eigenvalues at q-k
c zzmel  =    <psi(q,t') | psi(q-k,t) B(k,R,i)>
c zmel   =    <psi(q,t') | psi(q-k,t) M(k,R,i)>
c dimensions:
c nqibz   = number of k-points in the irreducible BZ
c natom   = number of atoms
c nctot   = no. allowed core states
c nbloch  = total number of Product basis functions
c nlnmx   = maximum number of l,n,m
c nlmto   = the total number of the argumentation wave.
c ngrp    = no. group elements (rotation matrices)
c niw     = no. frequencies along the imaginary axis
c nw      = no. frequencies along the real axis
c
c zsec(iw, t, iq) = <psi(q,t) |SEc (iw)| psi(q,t)>
c----------------------------------------------------------------------
      use m_readqg
      use m_readeigen
      implicit none
      integer(4) :: ntq, natom,nqbz,nqibz,ngrp,nq,nw,niw,
     &  nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,
     &   irot,        iqisp,ikpisp,isp,nsp,  nlnmx, !iq, idxk,
c     &  iwr1,iwr2,iwr3,iwr4,iwc1,iwc2,iwc3,iwc4,
     &  ip, it,itp,   !ifcphi,  ! ifrb,ifcb,ifrhb,ifchb, 
     i  iiclass,mdim(*),
     i  ifrcw,ifrcwi, !iindxk,
     &  ifvcfpout,ndummy1,ndummy2,kx,kr,ngc,ngb,nbloch, !n1,n2,n3, k,
     &  kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0,
     &  nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,ixs,ixsmx,
     &  mdimx, nwx,niwx,
     &  itq(ntq),iatomp(natom),      
     &  nstar(nqibz),irk(nqibz,ngrp),kount(nqibz,nq),
     &  ntq0,ntq0p,ntq0m
c
      real(8) :: q(3),qbas(3*3),ginv(3*3),tr(3,natom),
     &  wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),
     &  freq(nw),freqx(niw),wx(niw),expa(niw),
     &  eq(nband,nq), ekq(nband), ekc(nctot+nband),
c     &  rbq (nlmto,ntq),rbkq (nlmto,nband),
c     &  cbq (nlmto,ntq),cbkq (nlmto,nband),
c     &  rhbq(nlmto,ntq),rhbkq(nlmto,nband),
c     &  chbq(nlmto,ntq),chbkq(nlmto,nband),
     &  tpi,ef,ef2,esmr,esmr2,efp,efm,wfac,wfacx,we,esmrx,ua,
     &  dw,wtt,wexx,www,exx,exxq,efptmp,efmtmp
c      complex(8) :: zsec(-1:1,ntq,nq)
c      real(8)    ::  shtw
c                       ! This shft by Ferdi was to avoid some artificial resonance effects.
c                       ! shtw can be zero for esmr/=0 given by takao.
c
      integer(4):: ngpmx, ngcmx, ngpn(nqbz), !ngcni(nqibz), 
     &  ngvecpB(3,ngpmx,nqbz), igc, !ngveccBr(3,ngcmx,nqibz),
     &  nadd(3)
      real(8) :: wgt0(nq0i,ngrp),qk(3),qfbz(3),
     &  qbasinv(3,3), det,qdiff(3),add(3),symope(3,3),
     &  qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot),
     &   ppb(1),   !,pdb(1),dpb(1),ddb(1), !
     &  coh(ntq,nq)
      complex(8)::   alagr3z,wintz !geigB  (ngpmx,nband,nqbz),
c
      real(8),allocatable::  !rmel(:,:,:),cmel(:,:,:),
     &                 rmelt(:,:,:),cmelt(:,:,:)
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:),
     &     zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:)

c for exchange --------------------
      logical :: exchange,screen,cohtest,tote
      real(8),allocatable:: 
     &   w1p(:,:,:),w2p(:,:,:),w3p(:,:)
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)

c- debug write ---------------------
      logical :: debug=.false.

cccccccccccccc
c tetra
c      integer(4) :: ntqx
c      real(8)    :: wtet  (nband,nqibz,1:ntqx)
c                ! where the last index is 3*itq+iw-1,itq=1,ntq,iw=-1,1
      integer(4) :: ibzx(nqbz)
      real(8)    :: wtetef(nband,nqibz)
      logical    :: tetraex
cccccczzcccccc

      complex(8) :: wintzav,wintzsg

c worksize in megabytes (1 word = 4 bytes)
c     parameter (mbytes=60)
c     parameter (iwksize=mbytes * 1000 *1000/4)
c      integer w
c      common /w/ w(1)

      integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen
      integer(4),save::ifzwz=-999

      integer(4) :: iwini, iwend, ia
      real(8)    :: esec, omega(ntq, iwini:iwend)
      complex(8) :: zsec(iwini:iwend,ntq,nq)
      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)

      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)

      integer(4)::ixc

      real(8)::weavx,wfaccut=1d-10

      logical :: GaussSmear
      real(8) :: ebmx,ddw
      integer(4):: nbmx,nbmxe,nstatetot


c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4)::verbose,nstbz(nqbz),bzcase,iqini,iqend
      real(8):: wgtq0p

      integer(4):: iqindx
      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      integer(4):: nrec
      logical :: test_omitq0p
c--------------------------------------------------------------------
ccccc zvz test cccccccccccccccccccccccccc
      integer(4):: ngbx
c      complex(8):: vcoul(ngbx,ngbx)
      complex(8),allocatable:: vzz(:,:,:)
      complex(8):: zvz,zvz1
      integer(4):: ib1,ib2,ifix
ccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      kount=1
c      return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write(6,*)'iwini iwend=',iwini,iwend

c      inquire(file='DEBUG', exist=debug)
      debug=.false.
      if(verbose()>=100) debug=.true.

c
      if(debug) then
        write(6,"(' sxcf_v2z: q  =',3f8.4)")q
        write(6,*)' sxcf_v2z: new version: gaussian smearing=',GaussSmear()
        print *,' wintzsg wfacx weavx are related parts...'
      endif
      tpi = 8d0*datan(1.d0)
c      iq  = idxk (q,qbz,nqbz) ! index for q
c      write(6,"(' iq q  =',i4,3f8.4)")iq,q

cccc      iq          = idxk (q,qbze,nqbze) ! index for q
      ekc(1:nctot)= ecore(1:nctot)  ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
      call dinv33(qbas,0,qbasinv,det)

c work arrays for psi2br.f
      if(debug) print *, ' sxcf: 1'
      allocate(expikt(natom))

c direct access read e,b and hb at q
c      call rwdd1   (ifev, iq, nband,eq(1,ip))
      call readeval(q,isp, eq(1,ip))
c      iqisp =isp + nsp*(iq-1)
c      call rwdd1   (ifcphi, iqisp, 2*nlmtobnd, cphikq)
      call  readcphi(q, nlmto,isp, quu, cphikq)
      if(debug) print *, ' sxcf: 2'

c b and hb of states {t} at q
      do      it = 1,ntq
        itp        = itq(it)
        cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
      enddo
      if(debug) print *, ' sxcf: 3'

      if(bzcase()==1) then
        if(abs(sum(qibz(:,1)**2))/=0d0) stop ' sxcf assumes 1st qibz/=0 '
        if(abs(sum( qbz(:,1)**2))/=0d0) stop ' sxcf assumes 1st qbz /=0 '
      endif

c-----
      if(exchange) then
        rewind  ifvcfpout
        read(ifvcfpout) ndummy1, ndummy2
      endif

c===============================
c loop over irreducible k-points
c===============================
ccccccccccccccccccccccccccccccccc
c      iii = ivsumxxx(irk,nqibz*ngrp)
c      print *,' sxcf:sum non-zero irk=',iii
c      stop "sss"
c
ccccccccccccccccccccccccccccccccc
      if(bzcase()==1) then
        kx = 1  ! qibz(:,1)=0 contribution for kcount
        if(irk(kx,irot)/=0) kount(kx,ip)= kount(kx,ip) + 1
      endif

      if (tote) exxq = 0.d0

c === main loop start ====================================================
      iqini=2
      if(bzcase()==2) iqini=1
      iqend=nqibz+nq0i
cccccccccccccccccccccccccccccc
      if(test_omitq0p()) then
        iqend=nqibz
        print *,'iqend=',iqend
      endif
cccccccccccccccccccccccccccccc
      do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
        if(debug)  print *, ' sxcf: goto loop kx=',kx

        if( kx <= nqibz ) then
c          k  = kx
          kr = irk(kx,irot) ! index for rotated k in the FBZ
          qibz_k= qibz(:,kx)
          qbz_kr= qbz (:,kr)
        else
c          k = 1  ! corresponds to q=0
c          kr= 1  ! corresponds to q=0
c          k = iqindx((/0d0,0d0,0d0/), ginv, qibz,nqibz)
c          kr = iqindx((/0d0,0d0,0d0/), ginv, qbz,  nqbz)
          kr=-99999 !for sanity check
          qibz_k= 0d0
          qbz_kr= 0d0
        endif

        call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
c        ngc = ngcni(k)  ! k-points in IBZ
        ngb = nbloch + ngc

c --- exchange
        if(exchange) then
          allocate(vcoul(ngb,ngb))
          read(ifvcfpout) vcoul(1:ngb,1:ngb)
        endif

c- weight check for cycle or not.
        if( kx <= nqibz ) then
          if (kr == 0)    then
            if(exchange) deallocate(vcoul)
            cycle
          endif
          kount(kx,ip)= kount(kx,ip) + 1     ! count the no. times k
          ! appears in the 1st BZ
cccccccccccccccccccccccccccccccccccccccccccccccccc
c         print *,' irot,ip, k, kount in  =',irot, ip, k, kount(k,ip)
c          deallocate(vcoul)
c          cycle
c          print *,' kount out =',kount(k,ip)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
          if (kount(kx,ip) > nstar(kx)) stop 'sexc: too many stars'
        else
          if( wgt0(kx-nqibz,irot) == 0d0 ) then
            if(exchange) deallocate(vcoul)
            cycle
          endif
        endif

c qk = q - rk, rk is inside 1st BZ, not restricted to
c the irreducible BZ
c      call vminv   (q,qbz(1,kr),3,qk)
        qk =  q - qbz_kr !qbz(:,kr)

c bring q-rk into the 1st BZ ==> qfbz labelled by kp
c        call fbz     (qk,n1,n2,n3,qbas,ginv,w(iindxk),
c     o               qfbz,kp)
c         call fbz2(qk,  ginv,index_qbz,n_index_qbz, qbz,nqbz, !qbz=rk
c     o                qfbz,kp)


c        if(debug)  print *,'xxx2',kp,qk,n1,n2,n3,qbas,ginv
c read eigenvalues,b and hb at q-rk
c        call rwdd1   (ifev, kp,nband,ekq)
        call readeval (qk,isp, ekq)
c        ikpisp=isp + nsp*(kp-1)
c        call rwdd1  (ifcphi, ikpisp,2*nlmtobnd, cphikq)
        call  readcphi(qk, nlmto,isp, quu, cphikq)


c no. occupied (core+valence) and unoccupied states at q-rk
c      call dcopy   (nband,ekq,1,ekc(nctot+1),1)
        ekc(nctot+1:nctot+nband) = ekq (1:nband)
        nt0 = nocc (ekc,ef,.true.,nstatetot)
c
        ddw= .5d0
        if(GaussSmear()) ddw= 10d0
        efp= ef+ddw*esmr
        efm= ef-ddw*esmr
        nt0p = nocc (ekc,efp,.true.,nstatetot)
        nt0m = nocc (ekc,efm,.true.,nstatetot)

ctetra
c        if(tetraex) then
ccccccccccccccccccccccccccccccccccccc
c           stop 'tetraex is not implimented yet ...'
ccccccccccccccccccccccccccccccccccccc
c         irkp = ibzx(kp) kp???  ! IBZ corresponding to qk
c         do i=1,nband
c          if( wtetef(i,irkp) == 0d0 ) exit
c         enddo
c         if(sum(abs(wtetef(i:nband,irkp)))/=0d0)
c     &     stop ' sxcf: wtetef sum err1'
c         nt0p = i-1 + nctot
c        endif

        if(exchange) then
          nbmax = nt0p-nctot
          if(debug) print *,' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
        else
c-------------
          nbmax = nband
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          nbmxe = nocc (ekc,ebmx,.true.,nstatetot)-nctot
c          nbmax  = min(nband,nbmx,nbmxe)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          print *,' nbmax=',nbmax
c-------------
        endif
        nstate = nctot + nbmax ! = nstate for the case of correlation

c====================================================
c matrix elements <psi(q,t') | psi(q-rk,t) B(rk,R,i)>
c including the phase factor exp(ik.T)
c B(rot*k,r) = B(k,invrot*r)
c====================================================

c> phase factors exp(ik.T(R))
c        call cv      (tpi,qibz(1,k),3,qfbz)
c        call expkt   (qfbz,tr,natom,
c     o              w(icoskt),w(isinkt))
c equivalent with above phase
        do ia = 1,natom
          expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
        end do
c
        if(debug) print *, ' sxcf: tr=',tr
        if(debug) print *, ' sxcf: goto psicb2'
c> matrix elements
c>> core
        nt   = nctot + nbmax ! = nstate for the case of correlation
        ntp0 = ntq
        allocate( zzmel(nbloch,nt,ntq)) !     &      rmel    cmel(nbloch,nt,ntq))
        call psicb_v2  (icore,ncore,ntp0,iclass,
     i              dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i              cphiq, !rbq,cbq,rhbq,chbq,
     i              ppb,   ! pdb,dpb,ddb,
     i              nlnmv,nlnmc,mdim,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,nt,ntp0,natom,nclass,
     d              nl,nnc,                ! nt,ntp0 denotes dimenstion of rmel
     o              zzmel) !rmel,cmel)
        if(debug)  print *, ' sxcf: goto psi2bc1'
c>> valence
        call psi2b_v2  (nbmax, ntp0,iclass,
     i    dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i              cphikq, !rbkq,cbkq,rhbkq,chbkq, !occ    q-rk nband
     i              cphiq,  !rbq,cbq,rhbq,chbq,     !unocc  q    ntp0
     i              ppb,  !pdb,dpb,ddb,
     i              nlnmv,nlnmc,mdim,nctot,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0,
     d              natom,nclass,
     o             zzmel) ! rmel,cmel)

ccccccccccccccccccccccccccccccccccccccccccc
c test rmel cmel=0
c      rmel=0d0
c      cmel=0d0
ccccccccccccccccccccccccccccccccccccccccccc

        if(debug) write(6,"('sum of zmel abszmel=',4d23.16)")
     &                 sum(zzmel),sum(abs(zzmel) )

cccccccccccccccccccccccccccccccccccccccccccccc
c      do  it  = 1,nstate
c      it=32
c      itp=1
c     do  itp = 1,ntp0
c      do  ibl  = 1,nbloch
c        write(6,"(' rmel cmel  =',3i4,2d14.6)")it,itp,ibl
c     & ,rmel(ibl,it,itp),cmel(ibl,it,itp)
c      enddo
c      enddo
c     enddo
c      stop " *** TEST END *** "
cccccccccccccccccccccccccccccccccccccccccccccc

c-plane
        allocate( rmelt(ngb, nctot+nbmax, ntp0), ! nstate= nctot+nband
     &            cmelt(ngb, nctot+nbmax, ntp0))
        if(debug) print *, ' sxcf_v2z: goto drvmelp2'
c        call drvmelp( iq,   ntp0, ! q in FBZ
c     i                kp, nbmax,  ! q-rk
c     i                k,          ! k in IBZ for e-product basis
c     i       ngc,ngcmx,ngpn, ngpmx,ngvecpB,ngveccBr,geigB,nband,itq,
c     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
c     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
c     o       rmelt,cmelt)
        call drvmelp2( q,   ntp0, ! q in FBZ
     i                q - qbz_kr, nbmax,  ! q-rk
     i                qibz_k,           ! k in IBZ for e-product basis
     i           isp,ginv,
c     i       ngc,ngcmx,ngpmx,ngveccBr,nband,itq, 
     i       ngc,ngcmx,ngpmx,nband,itq, 
     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
     o       rmelt,cmelt)

ccccccccccccccccccccccccccccccccccccccccccc
c A matrix element check in Si. See the corresponding section of x0kf.f
c      write(6,"(3f13.5)") q         !ntp0
c      write(6,"(3f13.5)") qbz(:,kr) !product basis
c      print *,' sxcf: sumcheck rmelt cmelt='
c     &   ,sum(rmelt(1:ngb,4,1))
c     &   ,sum(cmelt(1:ngb,4,1))
c     &   ,sum(abs(rmelt(1:ngb,4,1)))
c     &   ,sum(abs(cmelt(1:ngb,4,1)))
c      stop ' a matrix ele check test end'
ccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
c testmelt0
c      rmelt(1: nbloch, :,:)=0d0
c      cmelt(1: nbloch, :,:)=0d0
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccc
c      do  it  = 1,nstate
c      itp=1
c      it =1
c     do  itp = 1,ntp0
c      print *,' ngb nbloch ngc =',ngb,nbloch,ngcni(k)
c      do  ibl  = 1,ngb
c        write(6,"(' rmelt cmelt =',3i4,4d14.6)")it,itp,ibl
c     &  ,rmelt(ibl,it,itp) ,cmelt(ibl,it,itp)
c     &  ,rmel (ibl,it,itp) ,cmel (ibl,it,itp)
c      enddo
cccccccccccccccccccccccccccccccccccccccccc

        deallocate(zzmel) !rmel,cmel)
c
        if(debug) print *, ' sxcf: goto wtt'
        if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum( rmelt),sum(cmelt )

c        if(kx<= nqibz) then
c         wtx = 1d0
c         wtt = wk(kr)
c        else
c         wtx = wgt0(kx-nqibz,irot)
c         wtt = wk(1)*wgt0(kx-nqibz,irot)
c        endif
        if(bzcase()==2)then
          if(kx<=nqibz) then
            wtt = wk(kr)
            if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
          elseif(kx>nqibz) then
c            wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
            wtt= wgt0(kx-nqibz,irot)
          endif
        else
          if(kx<= nqibz) then !          wtx = 1d0
            wtt = wk(kr)
          else
            wtt = wk(1)*wgt0(kx-nqibz,irot)  !       wtx = wgt0(kx-nqibz,irot)
            if(abs(wk(1)-1d0/dble(nqbz))>1d-10) stop 'sxcf:wk(1) inconsistent'
          endif
        endif

        if(debug) then
          write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
        endif

c--------------------------------------------------------
c --- exchange section ---
c--------------------------------------------------------
c
c S[i,j=1,nbloch] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c                        v(k)(i,j) <B(rk,j) psi(q-rk,n) |psi(q,t')>
c
c> z1p(j,t,t') = S[i=1,nbloch] <psi(q,t') | psi(q-rk,t) B(rk,i)> v(k)(i,j)
c
c
c      write(6,*)' vcoulsum=',sum(vcoul)
c      if(debug) write(6,*)'  sumz=',dcmplx(rmelt,cmelt),sum(vcoul)
c
        if(exchange) then
c --- screened exchange case
          if(screen) then
            allocate( zw (nblochpmx,nblochpmx))
            ix = 1
            nrec=(kx-2)*nw+ix
            if(bzcase()==2) nrec= (kx-1)*nw+ix
            read(ifrcw,rec=nrec) zw  ! Readin W(0) - v
            vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
            deallocate(zw)
          endif
c
          allocate( zmel (ngb, nctot+nbmax, ntp0),vcoult(1:ngb,1:ngb))
          zmel  = dcmplx (rmelt,cmelt)
          vcoult= transpose(vcoul)
c
          allocate( z1p(ngb,nctot+nbmax,ntp0), w3p( nctot+nbmax,ntp0))
          call matm( vcoult, zmel,
     &           z1p, ngb,ngb,(nctot+nbmax)*ntp0 )
          deallocate(vcoult, zmel)

c         print *,' nctot nbmax=',nctot,nbmax

          do itp = 1,ntp0
            do it  = 1,nctot+nbmax
              w3p(it,itp) = sum ( dreal(z1p(:,it,itp))*rmelt(:,it,itp)
     .                   +      dimag(z1p(:,it,itp))*cmelt(:,it,itp) )
            enddo
          enddo

c Write the Spectrum function for exchange May. 2001
          if(ifexsp/=0) then
            do it  = 1, nctot+nbmax
              do itp = 1,ntp0
                write(ifexsp,"(3i4, 3f12.4, ' ',d23.15,'  ',d23.15)")
     &      ip,itp,it, qbz_kr, ekc(it), -wtt*w3p(it,itp)
              enddo
            enddo
          endif

c DEBUG write
          if(debug) then
            do  it  = 1,nctot+nbmax
              do  itp = 1,ntp0
                write(6,"(' w3p =',2i4,2d14.6)") it,itp,w3p(it,itp)
              enddo
            enddo
          endif

c--- tetraex
c          if(tetraex) then
cccccccccccccccccccccccccccccccccccccc
cc           stop 'tetraex is not implimented yet ...'
cccccccccccccccccccccccccccccccccccccc
c           if (.not.tote) then
c            do itp = 1,ntp0  !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
c            zsec(-1,itp,ip)  = zsec(-1,itp,ip)
c     &      - wtt*sum( w3p(1:nctot,itp))
c     &      - wtx*sum(wtetef(1:nbmax,irkp)*w3p(nctot+1:nctot+nbmax,itp))
ccccc electron number check
cccc     &        - wtt*sum(wtetef(1:nbmax,irkp))/rydberg()
ccccc
c            enddo
c           else
c            do itp = 1,ntp0
c              wfac = wtetef(itq(itp),ip)
c              if(ixc==5) then
c                if(nctot/=0) stop ' sxcf_v2z: nctot/=0'
c                exxq = exxq  -wfac*wtx*
c     &          sum( wtetef(1:nbmax,irkp)*w3p(nctot+1:nctot+nbmax,itp))
c              elseif(ixc==6) then
c                exxq = exxq - wfac*wtt*sum( w3p(1:nctot,itp))
c              endif
c            enddo
c           endif
c           deallocate( z1p,w3p,rmelt,cmelt,vcoul)
c           cycle
c          endif

c--- Correct weigts wfac for valence by esmr
          do it = nctot+1, nctot+nbmax
cgaussian
            wfac = wfacx(-1d99, ef, ekc(it), esmr)
            if(debug) print *, ' it wfac=',it, wfac,ef, ekc(it), esmr
            w3p(it,1:ntp0) = wfac * w3p(it,1:ntp0)
          enddo

          if (.not.tote) then
            do itp = 1,ntp0  !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
              zsec(iwini,itp,ip) = zsec(iwini,itp,ip) 
     &        - wtt * sum( w3p(:,itp) )
            enddo
          else
            do itp = 1,ntp0
cgaussian
              wfac = wfacx(-1d99, ef2, eq(itq(itp),ip), esmr2) !june-2002 takao
              w3p(1:nctot+nbmax,itp) = wfac * w3p(1:nctot+nbmax,itp)
              exxq = exxq - wtt * sum( w3p(:,itp) )
            enddo
          endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctestzzz
c          do itp = 1,ntp0  !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
c            if(kx ==2)
c     &      zsec(iwini,itp,ip) = zsec(iwini,itp,ip)
c     &      - wtt * sum( w3p(1:1,itp) )
c          enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccc
c      if(kx>nqibz) then
c         do ix= 1,ngb
c           write(666,*)ix,rmelt(ix, 1,1)
c         enddo
c         do itp = 1,ntp0  !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
c           zsec(iwini,itp,ip) = zsec(iwini,itp,ip)
c     &     - wtt *sum(abs(dcmplx(rmelt(:,1,itp),cmelt(:,1,itp))))
c         enddo
c      endif
cccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccc
c      do itp=1,ntp0
c      write(6,"(3i3,3d13.6)")
c     &    kx,itp,ip,zsec(iwini,itp,ip),sum( w3p(:,itp))
c      enddo
cccccccccccccccccccccccccccccccc

          deallocate( z1p,w3p,rmelt,cmelt,vcoul)
          cycle
        endif
c-- End of exchange section --------------




c--------------------------------------------------------------------------
c--- correlation section --------------------------------------------------
c--------------------------------------------------------------------------

c--- The matrix elements zmel.
        allocate( zmel (ngb, nstate, ntp0) )
        zmel = dcmplx (rmelt,-cmelt)
        deallocate(rmelt,cmelt)


        if(debug) print *,' end of zmel'
c================================================================
c The correlated part of the self-energy:
c S[n=all] S[i,j=1,nbloch]
c <psi(q,t) |psi(q-rk,n) B(rk,i)>
c  < [w'=0,inf] (1/pi) (w-e)/{(w-e)^2 + w'^2} Wc(k,iw')(i,j) >
c                                <B(rk,j) psi(q-rk,n) |psi(q,t)>
c e = e(q-rk,n), w' is real, Wc = W-v
c================================================================

c====================================================================
c contribution to SEc(qt,w) from integration along the imaginary axis
c====================================================================

c------------------------------------------------
c loop over w' = (1-x)/x, frequencies in Wc(k,w')
c {x} are gaussian points between (0,1)
c------------------------------------------------
        allocate( zw (nblochpmx,nblochpmx), zwz0(nstate,ntp0) )
        ix = 1
        nrec=(kx-2)*nw+ix
        if(bzcase()==2) nrec= (kx-1)*nw+ix
        read(ifrcw,rec=nrec) zw  ! direct access read Wc(0) = W(0) - v
        call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,   
     o    zwz0) ! zwz0 = zmel*(W(0)-v)*zmel

c COH term test ----- The sum of the all states for zwz00 gives the delta function.
        if(cohtest) then
          do itp = 1,ntq
            coh(itp,ip)  = coh(itp,ip) 
     &       + .5d0*wtt*sum(dreal(zwz0(1:nstate,itp)))
          enddo
          deallocate(zw,zwz0,zmel)
          cycle
        endif
c
        nx  = niw
        if(niw <1) stop " sxcf:niw <1"
        allocate( zwz(niw,nstate,ntp0),  zwzi(nstate,ntp0) )
        if(screen) allocate(zwz00(nstate,ntp0))
c
        do ix = 1,nx     ! imaginary frequency w'-loop
          nrec=(kx-2)*niw+ix
          if(bzcase()==2) nrec= (kx-1)*niw+ix
          read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis

c zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
c                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
c        do itp = 1,ntp0
c        do  it = 1,nstate
c          zwz(ix,it,itp) = sum(
c     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
c        enddo
c        enddo
          call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb, 
     o      zwz(ix,1:nstate,1:ntp0))

ccccccccccccccccccccccccccccccccccccccc
c zwz is the diagonal term and real.
c      write(6,*) ' img sumzwz=',sum(abs(imag(zwz(ix,1:nstate,1:ntp0))))
ccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccc
c      write(6,*) ' sumzwz=',sum(zmel)
c     & ,sum(zw(1:ngb,1:ngb)),sum(zwz(ix,1:nstate,1:ntp0))
ccccccccccccccccccccccccccccccccc
        enddo



cccccccccccccccccccccccccccccccccccccccccccccc
c      do  it  = 1,nstate
c      itp=1
c      ix =1
c     do  itp = 1,ntp0
c      do  ix = 1,nx
c        write(6,"(' zwz on img =',3i4,2d14.6)")ix,it,itp,zwz(ix,it,itp)
c      enddo
c      enddo
c     enddo
c      stop " *** TEST END *** "
cccccccccccccccccccccccccccccccccccccccccccccc




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

ccccccccccccccccccccccccccc
c       if(ifzwz==-999) ifzwz = iopen('zwz.hsfp0',1,-1,0)
c        do itp = 1,ntq
c        do it  = 1,nstate
c           do iw=1,niw
c             write(ifzwz,"(4i5,2d24.16)") iw,it,itp,kx,zwz(iw,it,itp)
c           enddo
c           iw = 0
c           write  (ifzwz,"(4i5,2d24.16)") iw,it,itp,kx,zwz0(it,itp)
c           write(ifzwz,*)
c        enddo
c        enddo
ccccccccccccccccccccccccccc


c-------------------------
c loop over w in SEc(qt,w)
c-------------------------
        do      iw = iwini,iwend
c frequency integration along the imaginary axis, s. wint.f
c for each e(q-rk,n) and w in SEc(qt,w)
          do     itp = 1,ntq
            do      it = 1,nstate
              we =.5d0*( omega(itp,iw) -ekc(it)) != .5d0*( eq(itq(itp),ip)+2d0*(dble(iw)-shtw)*deltaw-ekc(it))
              esmrx = esmr
              if(it <= nctot) esmrx = 0d0
              if(GaussSmear()) then
                zwzi(it,itp) =
cGaussian smearing
     &      wintzsg ( zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx,
     &              esmrx)
              else
                zwzi(it,itp) =  
     &      wintzav ( zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx,
     &              esmrx)
              endif
cccccccccccccccccccccccccc
cc     .    wintz (zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua,expa,we,nx)
ccccccccccccccccccccccccccc
            enddo
          enddo

c sum over both occupied and unoccupied states and multiply by weight
          do     itp = 1,ntq
            zsec(iw,itp,ip)  = zsec(iw,itp,ip) + wtt*sum(zwzi(:,itp))
          enddo

c end of SEc w-loop
        enddo
        if(debug) then
          print *,' ntq nstate sum(zwzi)=',ntq,nstate,sum(zwzi)
          print *,' ntq nstate sum(zwz )=',ntq,nstate,sum(zwz)
          do itp = 1,ntq
            write(6,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
          enddo
        endif
        deallocate(zwz,zwz0,zwzi)

c==============================================
c contribution to SEc(qt,w) from the poles of G
c==============================================

ccccccccccccccccccccccccccccccccccccc
c      cycle
ccccccccccccccccccccccccccccccccccccc


        if(debug) print *,' go to poles'


c---------------------------------------
c maximum ixs reqired.
c---------------------------------------
c      print *,' ekc at nt0p nt0m+1=', ekc(nt0p),ekc(nt0m+1)
c      print *,'  nt0p nt0m+1=', nt0p, nt0m+1
        ixsmx = 1
        do 3001 iw  = iwini,iwend
        do 3001 itp = 1,ntq
          esec  = omega(itp, iw)  ! eq(itq(itp),ip) + 2.d0*(dble(iw)-shtw)*deltaw
          if (esec <= ef) then
            do  it = 1,nt0p            ! nt0p corresponds to efp
c             we   = abs( 0.5d0*(ekc(it) - esec) )
              esmrx = esmr
              if(it<=nctot) esmrx = 0d0
cgaussian
              wfac = wfacx(esec,ef, ekc(it),esmrx)
              if(GaussSmear()) then
                if(wfac<wfaccut) cycle
                we   = abs( .5d0*(weavx(esec,ef,ekc(it),esmr)-esec))
              else
                if(wfac==0d0) cycle
                we   = abs( 0.5d0*(ekc(it) - esec) )
              endif

              if(it<=nctot) then
                if(wfac/=0) stop " it<=nctot.and.wfac/=0"
              endif

              if(wfac<wfaccut) cycle
              ixs  = idint (we/dw) + 1
              if(ixs>ixsmx) then
                wexx  = we
                ixsmx = ixs
c------
                if (ixsmx+2 > nw) then
                  write (*,*) ' nw',nw
                  write (*,*) ' dw ',wexx,dw
                  write (*,*) ' esec ekc(it) ef ', esec,ekc(it),ef,wfac
                  stop ' sxcf 1112: |w-e| > max(w)'
                endif
c------
              endif
            enddo
          else
            do  it = nt0m+1,nstate    ! nt0m corresponds to efm
c             we   = abs( 0.5d0*(esec - ekc(it)) )
cgaussian
              wfac = wfacx(ef,esec, ekc(it),esmr)
              if(GaussSmear()) then
                if(wfac<wfaccut) cycle
                we = abs( .5d0*(esec - weavx(ef,esec, ekc(it),esmr)))
              else
                if(wfac==0d0) cycle
                we = abs( 0.5d0*(esec - ekc(it)) )
              endif

              ixs  = idint (we/dw) + 1
              if(ixs>ixsmx) then
                wexx  = we
                ixsmx = ixs
c------
                if (ixsmx+2 > nw) then
                  write (*,*) ' ef esec ekc esmr',ef,esec,ekc(it),esmr
                  write (*,*) ' nw',nw
                  write (*,*) ' dw ',wexx,dw
                  write (*,*) ' esec ekc(it) ef ', esec,ekc(it),ef
                  stop ' sxcf 222: |w-e| > max(w)'
                endif
c------
              endif
            enddo
          endif
 3001   continue  !end of SEc w and qt -loop
        nwx = ixsmx+2
        if (nwx > nw) then
          stop ' sxcf xxx: |w-e| > max(w)'
        endif
        if(debug) print *,' nwx=',nwx


c--------------------------------------------
c loop over real frequencies of Wc(k,w')(i,j)
c--------------------------------------------
        allocate( zwz(nw,nstate,ntp0) )
        do      ix = 1,nwx                   ! imaginary frequency w'-loop
          nrec=(kx-2)*nw+ix
          if(bzcase()==2) nrec= (kx-1)*nw+ix
          read(ifrcw,rec=nrec) zw  ! Readin (W-v)(k,w')(i,j) at k and w' on imag axis
c zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
          call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,   
     o    zwz(ix,1:nstate,1:ntp0))
        enddo
        deallocate(zmel)

cccccccccccccccccccccccccccccccc
c      do  it  = 1,nstate
c     do  itp = 1,ntp0
c      do  ix  = 1,nwx
c        write(6,"(' zwz on real =',3i4,2d14.6)")ix,it,itp,zwz(ix,it,itp)
c      enddo
c     enddo
c     enddo
cccccccccccccccccccccccccccccccc

c
        if(screen) then
          do ix = 1,nwx
            zwz(ix,:,:)=zwz(ix,:,:) - zwz00
          enddo
          deallocate(zwz00)
        endif


c-------------------------------
c loop over w and t in SEc(qt,w)
c-------------------------------
        do 2001 iw  = iwini,iwend
        do 2001 itp = 1,ntq
c...........................................
c loop over occupied states (q-k,n); w.le.ef
c...........................................
          esec     = omega(itp,iw)              !eq(itq(itp),ip) + 2.d0*(dble(iw)-shtw)*deltaw
          if (esec <= ef) then
            do  it = 1,nt0p            ! nt0p corresponds to efp
! energy window is   esec<  e  <ef
! each ekc has the width, \pm 0.5*esmr.
c           we   = abs( 0.5d0*(ekc(it) - esec) )
ccccccccccccccccccccccccccccccccccc
ctest
c           we=0d0
ccccccccccccccccccccccccccccccccccc
              esmrx = esmr
              if(it<=nctot) esmrx = 0d0
cgaussian
              wfac = wfacx(esec,ef, ekc(it),esmrx)
              if(GaussSmear()) then
                if(wfac<wfaccut) cycle
                we = abs( .5d0*(weavx(esec,ef, ekc(it),esmr)- esec))
                if(it<=nctot) then                                     !faleev
                  if(wfac>wfaccut) stop "sxcf: it<=nctot.and.wfac/=0"
                endif
              else
                if(wfac==0d0) cycle
                we   = abs( 0.5d0*(ekc(it) - esec) )
                if(it<=nctot) then                                     !faleev
                  if(wfac/=0) stop "sxcf:  it<=nctot.and.wfac/=0"
                endif
              endif

              if(debug)write(6,"( ' xxx1',10d13.6)") esec,ef, ekc(it),wfac
c three-point interpolation for Wc(we)
              ixs  = idint (we/dw) + 1
              if(ixs+2>nwx) then
                print *,' ixs+2=',ixs,we,nwx,dw
                stop ' sxcf: ixs+2>nw xxx1'
              endif
              zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     .      - wfac*wtt*alagr3z(we,freq(ixs),zwz(ixs,it,itp))
            enddo
          else
c.............................................
c loop over unoccupied states (q-k,n); w.gt.ef
c.............................................
            do  it = nt0m+1,nstate    ! nt0m corresponds to efm
c           we   = abs( 0.5d0*(esec - ekc(it)) )
cccccccccccccccccccccccccccccccccc
ctest
c           we=0d0
cccccccccccccccccccccccccccccccccc
c          wfac = wfacx(ef,esec, ekc(it),esmr)
c          if(wfac<wfaccut) cycle
              wfac = wfacx(ef,esec, ekc(it),esmr)
              if(GaussSmear()) then
                if(wfac<wfaccut) cycle
                we   = abs( .5d0*(esec-weavx(ef,esec, ekc(it),esmr)))
              else
                if(wfac==0d0) cycle
                we = abs( .5d0*(esec- ekc(it)) )
              endif

              if(debug)write(6,"(' xxx2',10d13.6)") esec,ef, ekc(it),we,wfac
c three-point interpolation for Wc(we)
              ixs  = idint (we/dw) + 1
              if(ixs+2>nwx) then
                print *,'debug dat=',nwx,ixs
                stop ' sxcf: ixs+2>nwx xxxx2'
              endif
              zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     .      + wfac*wtt*alagr3z(we,freq(ixs),zwz(ixs,it,itp))
            enddo
          endif
 2001   continue  !end of SEc w and qt -loop
        if(debug) print *,' end of do 2001'


        if(debug) then
          do itp = 1,ntq
            write(6,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
          enddo
        endif
        deallocate(zwz,zw)
 1100 continue  ! end of k-loop

c 0.25 = 0.5 * 0.5
c 0.5: normalization factor of wik
c 0.5: EEX = (Sum SEx) * 0.5
      if (tote) then
c         if (tetraex) then
ccccccccccccccccccccccccccccccccccccccc
ccc           stop 'tetraex is not implimented yet ...'
cccccccccccccccccccccccccccccccccccccccc
c            exx = exx + nstar(ip) * exxq * 0.5d0
c         else
        exx = exx + wik(ip) * exxq * 0.25d0
c         endif
      endif

c      call rlse (iwr1)
ccccccccccccccccccccccccccccc
c      print *," end kount ip irot=", irot, ip, sum(kount(:,ip))
c       print *," sumkount 1=",irot,ip,sum(kount(:,:))
cccccccccccccccccccccccccccc
      return
      end


c -------------------------------------------------------------------
c --------------------------------------------------------------------
      complex(8) function wintzav(v,v0,x,wt,a,expa,we, 
     i                                nx,esmr)
C- Smeared-pole version of wintz. takao developed from wintz by ferdi.
C Assume thaat each eigenvale, and w-e as a result, has the width of esmr( Ry).
c
      implicit none
      integer(4) :: nx,i
      real(8) :: x(nx),wt(nx),expa(nx),a,we,we2,weh,wel,weh2,wel2,cons,
     & omd,omd2,rup,rdn,sss
      complex(8) ::v(nx),v0,sum,   wintzav1  !, wintzav2,sum2
      complex(8) :: wwz,wintz
      real(8)    :: esmr,wintsf,errsum=0d0
      integer(4) :: ie,nav = 2000
      real(8)    :: pi=3.1415926535897932d0, rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
     &
c---------------------------------------------------
c      print *,' wintzav:'
      if(esmr==0d0) then
        wintzav=wintz(v,v0,x,wt,a,expa,we,nx)
        return
      endif
c-----------------------------------------------------
c- check part by direct integration by wintz.
c- this procedure is equivalent to whole of this routine. But time-consuming.
c      wwz =0d0
c      do ie = 1, nav
c        wex = we + 0.5d0*esmr*(-0.5d0 + dble(ie-1)/(nav-1))
c        www=1d0
c        if(ie==1.or.ie==nav) www=0.5d0
c        wwz = wwz+ www*wintz(v,v0,x,wt,a,expa,wex,nx)
cxx       write(1112,"(4d24.16)") wex,wintz(v,v0,x,wt,a,expa,wex,nx)
c      enddo
c      wintzav1=wwz/(nav-1)
c      wintzav=wintzav1
c      return

c frequency integral
      we2  = we**2
      weh  = we + 0.25d0*esmr !because we is in a.u and esmr is in Ry.
      wel  = we - 0.25d0*esmr
      weh2 = weh**2
      wel2 = wel**2
      sum  = 0d0
c
c      sum2 = 0d0
c      wintzav2 = 0d0

ccase1 --------------------
      if ( max(dabs(weh),dabs(wel)) < rmax/a) then
ccase2 --------------------
ccase2 if (.false.) then
        do  i = 1,nx
          omd  = 1d0/x(i) - 1d0
          omd2 = omd**2
c   x = 1/(1+w')
c   We smear each pole with a pole distribution between wel ans weh,
c   So  it means \int_wel^weh dwe [ wt(i)/(x(i)**2) * we/(we2 + omd**2)]
          cons  = 0.5d0*log( (weh2+omd2)/(wel2+ omd2)) /(x(i)**2)
          sum   = sum + wt(i)*(v(i) - v0*expa(i))*cons
ccccccccccccccccccccccccccccccccccccccccccccccccc
c        sum2  = sum2 + wt(i)*v(i) *cons
c      write(6,"(d15.6,2x,2d15.6,2x,2d15.6)")
c     & x(i), v(i)*cons,v0*expa(i)*cons
c      write(6,"(d15.6,2x,2d16.6,2x,2d16.6)")
c     & x(i), v(i),v0*expa(i)
ccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo

        if(abs(weh)>abs(wel)) then
          sss = 1d0
          rup = abs(weh)
          rdn = abs(wel) !;    iix=1
        else
          sss = -1d0
          rup = abs(wel)
          rdn = abs(weh) !;   iix=11
        endif
        wintzav  = - sum/pi - 0.5d0*v0* sss* wintsf(a*rdn,a*rup)/a
c       wintzav2 = - sum2/pi
ccccccccccccccccccccccccccccccccccccccccccccccccc
c      xxx=0d0
c      write(6,"(d15.6,2x,2d16.6,2x,2d16.6)")
c     & xxx, v0 ,v0
c      write(6,"(' - sum/pi =',2d13.6)") -sum/pi
c      write(6,"(' wintzav  =',2d13.6)")   wintzav
c      write(6,"(' - sum2/pi =',2d13.6)") -sum2/pi
c      stop ' ==== test end of wintz ========'
ccccccccccccccccccccccccccccccccccccccccccccccccc
      else
c simple integration scheme.
        do  i = 1,nx
          omd  = 1d0/x(i) - 1d0
          omd2 = omd**2
          cons =  0.5d0*log( (weh2+omd2)/(wel2+ omd2)) /(x(i)**2)
          sum  = sum + wt(i)*v(i)*cons
        enddo
        wintzav = -sum/pi  !;   iix =2
      endif
      wintzav = wintzav/(0.5d0*esmr) !because esmr is in Ry.

cccccccccccccccccccccccccccc
      write(1116, "(' wel weh wintzav =',2f8.3,4f14.6)" )wel,weh,wintzav
ccccccccccccccccccccccccccc

c      wintzav2 = wintzav2/(0.5d0*esmr)
c      write(1116, "(' wintzav diff =',4d13.6,f10.6)" ) wintzav,wintzav2
c     &   ,abs(wintzav2-wintzav)/abs(wintzav)
cccccccc
c      stop ' ==== test end of wintz ========'
c      if(wintzav/=0d0 ) then
c        err = abs(wintzav-wintzav1)
c        if(err>1d-6) then
c          write(3302,"(3d24.16, i3,' err= ',2d13.5)")
c     &    we, wintzav, iix, err
c        endif
c     endif
cccccccccccccccccccccccccccc
      end

c----------------------------------------------------
      real(8) function wintsf(rdn,rup)
c integral \int^rup_rdn exp(x**2)*derfc(x)
c up to x <rmax<15 which is limited by the range of derfc(x).
      implicit real*8 (a-h,o-z)
      integer(4),parameter:: nd=2001
      real(8) :: polinta, rmax=14.9d0
      real(8),save :: eeint(nd),xval(nd)
      logical :: init=.true.
c--------------
      if(rup>rmax) stop 'wintsf: too large rup'
      if(rup<rdn ) stop 'wintsf: rup<rdn'
      if(init) then
        eeint(1) = 0d0
        xval(1)  = 0d0
        do ix= 2,nd
          x0= rmax*dble(ix-2)/(nd-1)
          x = rmax*dble(ix-1)/(nd-1)
          xm= .5d0*(x + x0)
          xval(ix) = x
          sum  = 1/6d0 * exp(x0**2)*derfc(x0)
     &         + 4/6d0 * exp(xm**2)*derfc(xm)
     &         + 1/6d0 * exp(x **2)*derfc(x )
          eeint(ix) = eeint(ix-1) + rmax/(nd-1) * sum
ccccccccccccccccccccccc
c        write(1113,"(4d24.16)") x,derfc(x),exp(x**2)*derfc(x),eeint(ix)
ccccccccccccccccccccccc
        enddo
        init=.false.
      endif
      wintsf = polinta(rup,xval,eeint,nd) - polinta(rdn,xval,eeint,nd)
      end


c--------------------------------------------------------------------
      complex(8) function wintz(v,v0,x,wt,a,expa,we,
     d                                nx)
c takao complex version of wint by ferdi
c 92.02.24
c wintz = < [w'=-inf,inf] (i/2pi) v(w')/(w+w'-e) >

c the integration is done along the imaginary axis
c w' ==> iw', w' is now real
c wintz = - < [w'=0,inf] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2} >

c transform: x = 1/(1+w')
c wintz = - < [x=0,1] v(iw') (1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >

c the integrand is peak around w'=0 or x=1 when we=w-e=0
c to handel the problem, add and substract the singular part as follows:
c wintz = - < [x=0,1] {v(iw') - v(0)exp(-a^2 w'^2)}
c                   *(1/pi) (w-e)/{(w-e)^2 + w'^2}x^2 >
c        - (1/2) v(0) sgn(w-e) exp(a^2 (w-e)^2) erfc(a|w-e|)
c the second term of the integral can be done analytically, which
c results in the last term

c when we=w-e ==> 0, (1/pi) (w-e)/{(w-e)^2 + w'^2} ==> delta(w')
c the integral becomes -v(0)/2

c v       = v(iw')
c v0      = v(0)
c x       = s.o.
c wt      = weights for integration
c a       = a constant determining the range of exp(-a^2 w'^2)
c expa(x) = exp(-a^2 w'^2)
c we      = w - e
c nx      = no. x points

      implicit real*8 (a-h,o-z)
      dimension x(nx),wt(nx),expa(nx)
      complex(8) ::v(nx),v0,sum
      real(8) :: rmax=2d0 !rmax =2 is by ferdi. Is it optimum? See wintz
      data pi/3.1415926535897932d0/,
cmod 93.06.23
c    d     tol/1.d-6/
     d     tol/1.d-8/

c if w = e the integral = -v(0)/2
cmod 93.06.23
c     if (dabs(we) .lt. tol) then
c     wintz       = -0.5d0*v0
c     return
c     endif
      if (dabs(we) .lt. tol) then
        print *, ' we=',we
        stop 'wintz: |w-e| < tol'
      endif
c frequency integral
      we2        = we*we
      sum        = 0.d0

      if (dabs(we) < rmax/a) then
        do       i = 1,nx
          onemx      = 1.d0 - x(i)
          cons       = we/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + cons*(v(i) - v0*expa(i))*wt(i)
        enddo

        wintz      = -sum/pi 
     .           - 0.5d0*v0*dsign(1.d0,we)*dexp(we2*a*a)
     .                                    *derfc(a*dabs(we))
      else
        do       i = 1,nx
          onemx      = 1.d0 - x(i)
          cons       = we/(we2*x(i)*x(i) + onemx*onemx)
          sum        = sum + cons*v(i)*wt(i)
        enddo

        wintz      = -sum/pi
      endif

      return
      end

c--------------------------------------------------------------------
      complex(8) function alagr3z (x,xi,fi)
      implicit real*8 (a-h,o-z)
c evne function version of alagr3z
c return the interpolated value on x for fi(xi).
c We assume even function.
c Imag part is corrected to be >0
      dimension  xi(3)
      complex(8) :: fi(3)
      real(8) :: amatinv(3,3),amat(3,3),ratio

c--- Assure that dimag(fi) is positive.
      if(x<0d0) stop ' alagr3z: x<0d0'
      do i=1,3
        if(dimag(fi(i))>1d-8) then
          ratio= abs(dimag(fi(i))/dreal(fi(i)))
          if( ratio >1d-6) then !jan2004. I added a bit of alloance
            print *,' alagr3z: error! dimag(fi(i)) is positive',i,xi,fi
            stop    ' alagr3z: error! dimag(fi(i)) is positive'
          endif
        endif
      enddo
c
      amat(1:3,1) = 1d0
      amat(1:3,2) = xi(1:3)**2
      amat(1:3,3) = xi(1:3)**4
      call dinv33(amat,0,amatinv,detxx)
      alagr3z = dcmplx (
     &          sum (matmul(amatinv,dreal(fi)) * (/1d0,x**2,x**4/) ),
     &          sum (matmul(amatinv,dimag(fi)) * (/1d0,x**2,x**4/) ) )

c--- imag part of alagr3z is positive because interpolation could give
c    negative values even if all fi(i) <0.
c
      if(dimag(alagr3z)>0d0 ) then
        alagr3z = dcmplx( dreal(alagr3z),0d0)
      endif
      end

c
c -----------------------------------

      complex(8) function alagr3z_old (x,xi,fi)
c complex version of fi
c 92.03.02
c three-point interpolation with arbitrary mesh
c given a function fi at three points xi, the routine interpolates
c the function at x
c f(x) = [ { (x-x2)(x-x3) } / { (x1-x2)(x1-x3) } ] f1
c      + [ { (x-x1)(x-x3) } / { (x2-x1)(x2-x3) } ] f2
c      + [ { (x-x1)(x-x2) } / { (x3-x1)(x3-x2) } ] f3
c
c x  = the point at which the function is to be interpolated
c xi(3) = points where the function is given
c fi(3) = the function at xi
      implicit real*8 (a-h,o-z)
      dimension xi(3)
      complex(8) :: fi(3)
      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x23        = xi(2)-xi(3)
      alagr3z_old = xx2*xx3*fi(1)/(x12*x13)
     .           - xx1*xx3*fi(2)/(x12*x23)
     .           + xx1*xx2*fi(3)/(x13*x23)
      return
      end

c-------------------------------------
      subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntp0) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
      do itp = 1,ntp0
        do  it = 1,nstate
          zwz(it,itp) = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it,itp))
        enddo
      enddo
      deallocate(CC)
      end

c-------------------------------------
c   This does not work.
c      subroutine matzwz(zw,zmel, ntp0,nstate,ngb, zwz)
c      implicit none
c      integer(4) :: nstate,ntp0,itp,it,ngb
c      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
c     & ,wx(ngb)
c --- This routine is really time consuming! ---
c        do itp = 1,ntp0
c        do  it = 1,nstate
c        call mml(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp), wx(1:ngb))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp))*wx )
c           wx = matmul(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp)) * wx )
c        enddo
c        enddo
c      end

c-------------------------------------
      subroutine matzwz_notusedx(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),zwz(nstate,ntp0)
c     & ,wx(ngb)
c --- This routine is really time consuming! ---
      do itp = 1,ntp0
        do  it = 1,nstate
          zwz(it,itp) = sum( 
     &   dconjg(zmel(:,it,itp))*matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)))
c           wx = matmul(zw(1:ngb,1:ngb),zmel(1:ngb,it,itp))
c         zwz(it,itp) = sum( dconjg(zmel(:,it,itp)) * wx )
        enddo
      enddo
      end

