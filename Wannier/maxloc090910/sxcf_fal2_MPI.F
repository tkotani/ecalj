      subroutine sxcf_fal3z_RSMPI(kount,ip,q,irot,itq, ntq,
     i          itini,itfin,ef,ef2,esmr,esmr2,tr,
     i          iatomp,nsp,isp, !ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     i                  ifrcw,ifrcwi, 
     i                  qbas,ginv,
     i         qibz,qbz,wk,nstbz,wik,nstar,irk,  !koun,,iindxk
c
     i                  iclass,mdim,nlnmv,nlnmc,
     i                  icore,ncore,imdim,
     i                  ppb, !pdb,dpb,ddb,
     i                  freq_r,freqx,wx, !expa,ua, !deltaw,freq
     i                  dw,ecore, 
     d                nlmto,nqibz,nqbz,nctot,
c     i          index_qbz, n_index_qbz,  !jan2004
     d                  nl,nnc,nclass,natom,
     d                  nlnmx,mdimx,nbloch,ngrp, nw_i,nw ,niw,niwx,nq, 
                             ! nw_i,nw feb2006
                             ! before  feb2006, nw has the meaning of nw+1
! See freq_r
c     &     nblochpmx ,ngpn,ngcni,ngpmx,ngcmx,geigB,ngvecpB,ngveccBr,
     &     nblochpmx ,ngpmx,ngcmx, !ngveccBr,!Jan2004
     &     wgt0,nq0i,q0i,symope,alat, shtv,nband, ifvcfpout,
c     &     shtw,  
     &     exchange,tote,screen,cohtest, ifexsp, 
ctetra
ctetra     &     wtet,wtetef,
ctetra    &     ntqx,ibzx,tetraex, 
c
     i omega,iwini,iwend,
     i     nbmx,ebmx, !takao 18June2003 
     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat 
c     o zsec,coh,exx)
ctm 070621
     i  loffd,
ctm 081021
     i  lwang,iko_ix,iko_fx,nwf,umn,
     o  zsec,zsec2,coh,exx) ! acuumulation variable
      use m_readqg
      use m_readeigen
      use keyvalue
c RS: modules for MPI 
      use rsmpi
      use rsmpi_qkgroup
      use rsmpi_rotkindex
C- Calcualte diagonal part simga_ii(e_i)= <i|Re[S](e)|i> -------
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
c irk(k,R) = gives index in the FBZ with k{IBZ, R=rotation
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
c      implicit real*8(a-h,o-z)
      implicit none
      integer(4) :: ntq, natom,nqbz,nqibz,ngrp,nq,nw,niw,
     &  nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,nstatex,
     &   irot,  iqisp,ikpisp,isp,nsp,  nlnmx, !iq, idxk,
c     &  iwr1,iwr2,iwr3,iwr4,iwc1,iwc2,iwc3,iwc4
     &  ip, it,itp, !ifcphi,    ! ifrb,ifcb,ifrhb,ifchb, 
     i  iiclass,mdim(*),
     i  ifrcw,ifrcwi, !iindxk,
     &  ifvcfpout,ndummy1,ndummy2,kx,kr,ngc,ngb,nbloch, !n1,n2,n3, k,
     &  kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0,
     &  nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, !iwp sf 21May02
     &  mdimx, nwx,niwx,
     &  itq(ntq),iatomp(natom),      
     &  nstar(nqibz),irk(nqibz,ngrp),kount(nqibz,nq)
     &  ,itini,itfin
c
      real(8) :: q(3),qbas(3*3),ginv(3*3),tr(3,natom),
     &  wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),
     &  freqx(niw),wx(niw),!expa(niw),
     &  eq(nband,nq), 
     &   ekq(nband), ekc(nctot+nband),
     &  tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx, !ua,
     &  dw,wtt,wexx,www,exx,exxq
c      complex(8) :: zsec(-1:1,ntq,nq)
c      real(8)    ::  shtw  
c                       ! This shft is  to avoid some artificial resonance effects.
c                       ! shtw can be zero for esmr/=0 given by takao.                 
c
      integer(4):: ngpmx, ngcmx,  !ngcni(nqibz), !ngpn(nqbz),
     &   igc, !ngvecpB(3,ngpmx,nqbz),ngveccBr(3,ngcmx,nqibz),
     &  nadd(3)
      real(8) :: wgt0(nq0i,ngrp),qk(3), !qfbz(3),
     &  qbasinv(3,3), det,qdiff(3),add(3),symope(3,3),
     &  qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot),
     &   ppb(1), !pdb(1),dpb(1),ddb(1), !*
     &  coh(ntq,nq)
      complex(8)::   alagr3z,alagr3z2,wintz !geigB  (ngpmx,nband,nqbz),

c
      real(8),allocatable:: !rmel(:,:,:),cmel(:,:,:),
     &                      rmelt(:,:,:),cmelt(:,:,:) 
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:),
     &     zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:)
c for exchange --------------------
      logical :: exchange,screen,cohtest,tote
      real(8),allocatable:: 
     &   w1p(:,:,:),w2p(:,:,:),w3p(:,:)
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
    
c- debug write ---------------------
      logical :: debug=.true.

cccccccccccccc
c tetra
c      integer(4) :: ntqx
c      integer(4) :: ibzx(nqbz)
c      real(8)    :: wtet  (nband,nqibz,1:ntqx), wtetef(nband,nqibz) 
c                ! where the last index is 3*itq+iw-1,itq=1,ntq,iw=-1,1
c     logical    :: tetraex
cccccczzcccccc

      complex(8) :: wintzav,wintzsg

      integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen
      integer(4),save::ifzwz=-999

      integer(4) :: iwini, iwend, ia
      real(8)    :: esec, omega(ntq, iwini:iwend)
      complex(8) :: zsec(iwini:iwend,ntq,nq)
      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)
ctakao
      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)

      integer(4) :: nt_max, igb1,igb2,iigb,  nw_i  !nw_i is at feb2006
      complex(8),allocatable:: zmel1(:),zmel3(:)
      complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)  
      complex(8), allocatable :: zwz2(:,:),zw2(:,:,:),zmel2(:,:) !0 variant
      complex(8) ::  zz2 ,zwz3(3) 
      real(8) :: dd,omg_c,dw2
      real(8) :: freq_r(nw_i:nw)
      complex(8), allocatable :: zw3(:,:,:)
        

      real(8)::weavx,wfaccut=1d-10

      logical :: GaussSmear,gass
      real(8) :: ebmx,ddw
      integer(4):: nbmx,nbmxe,nstatetot

c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4)::verbose,nstbz(nqbz),bzcase,iqini,iqend
      real(8):: wgtq0p

      integer(4):: iqindx,nrec,kxx
      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      logical :: onlyQ0P, onlyimagaxis ,noq0p !,test_omitq0p,

      logical ::zwz3mode


      real(8):: ua_,expa_(niw),ua2,freqw,freqw1,ratio,ua2_(niw)
      logical :: ua_auto
      integer(4):: icc=0
      real(8),allocatable:: uaa(:,:)

c      logical ::testimx=.false.
ccccc zvz test cccccccccccccccccccccccccc
      integer(4):: ngbx
c      complex(8):: vcoul(ngbx,ngbx)
      complex(8),allocatable:: vzz(:,:,:),aaa(:)
      complex(8):: zvz,zvz1
      integer(4):: ib1,ib2,ifix
ccccccccccccccccccccccccccccccccccc
      integer(4) ::nbcut,nbcutc
      logical ::iww2=.true., oncew


C... 
      logical::smbasis
      integer(4):: nn,no,ifpomat,iclose,isx,iqx
      complex(8),allocatable:: pomat(:,:)
      real(8):: q_r(3)
      integer(4):: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo) 
      complex(8):: pomatr(nnmx,nomx,nkpo)
      real(8):: qrr(3,nkpo)

      real(8):: elxx,ehxx,ekxx,efxx
      integer(4):: ixsmin,iwm,iir,nwxi
      real(8)   :: fffr(3)
      complex(8):: zwzz(3)
ctm, 070621 for off-diagonal elements
      logical:: loffd
      complex(8) :: zsec2(ntq,ntq,nq,iwini:iwend),ctmp
      integer(4) :: itp2
      complex(8),allocatable:: w4p(:,:,:)
     &       , zzwz(:,:,:,:), zzwz0(:,:,:),zzwzi(:,:,:),zzwz00(:,:,:)
ctm, 081021, for Wannier gauge
      logical::  lwang
      integer(4):: iko_ix,iko_fx,nwf
      complex(8):: umn(iko_ix:iko_fx,nwf)
c--------------------------------------------------------------------
c RS: MPI parameters
      integer :: kx_local

c--------------------------------------------------------------------
      debug=.false.
c      if(verbose()>=90) debug=.true.
      if(verbose()>=90) then
        call RSMPI_Stop("sxcf_fal2_RSMPI: debug-mode is disabled")
      endif

c oct2005
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      nbcutc=nctot+nbcut

c      if(debug) then
c RS: ifile_rsmpi is defined in gwsrc/RSMPI_mod.F
      gass=GaussSmear() 
      write(ifile_rsmpi,"(' sxcf_fal3z_RSMPI: q  =',3f8.4)")q
      write(ifile_rsmpi,*)
     &     ' sxcf_fal3z_RSMPI: new version: gaussian smearing=',gass
      write(ifile_rsmpi,*) ' wintzsg wfacx weavx are related parts...'
c      endif

      tpi         = 8d0*datan(1.d0)   
c      iq         = idxk (q,qbz,nqbz) ! index for q
c      write(6,"(' iq q  =',i4,3f8.4)")iq,q
cccc      iq          = idxk (q,qbze,nqbze) ! index for q
      ekc(1:nctot)= ecore(1:nctot)  ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
      call dinv33(qbas,0,qbasinv,det) 

c work arrays for psi2br.f
c      if(debug) print *, ' sxcf: 1'
      allocate(expikt(natom))

c read eval and eigen
      call readeval (q, isp, eq(1,ip))
      if (lwang) then
        call  readcphi_mlw(q, nlmto,isp, quu, cphikq(:,1:nwf))
      else ! lwang
      call  readcphi(q, nlmto,isp, quu, cphikq)
      endif ! lwang
c      if(debug) print *, ' sxcf: 2'

c eigenfunction for itp at q
      do      it = 1,ntq
        itp        = itq(it)
        cphiq(1:nlmto,it) = cphikq(1:nlmto,itp)
      enddo
c     if(debug) print *, ' sxcf: 3'
       
      if(bzcase()==1) then
        if(abs(sum(qibz(:,1)**2))/=0d0) then
          call RSMPI_Stop(' sxcf assumes 1st qibz/=0 ')
        endif
        if(abs(sum( qbz(:,1)**2))/=0d0) then
          call RSMPI_Stop(' sxcf assumes 1st qbz /=0 ')
        endif
      endif

c-----
c RS: this if- block is no longer necessary because 
c            subroutine "set_vcoul_rsmpi" is called later
c      if(exchange) then
c        rewind  ifvcfpout
c        read(ifvcfpout) ndummy1, ndummy2
c      endif


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
        kx = 1                  ! qibz(:,1)=0 contribution for kcount
        if(irk(kx,irot)/=0) kount(kx,ip)= kount(kx,ip) + 1 
      endif
      
      if (tote) exxq = 0.d0

c --- gauss_img : interpolation gaussion for W(i \omega).
      call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)
      if(ua_<=0d0) then
        ua_auto =.true.
        write(ifile_rsmpi,"(' ua_auto=T')")
      else
        ua_auto =.false.
        do ix = 1,niw
        freqw     = (1d0 - freqx(ix))/ freqx(ix)
        expa_(ix) = exp(-(ua_*freqw)**2)
        enddo
      endif
      
c --- main loop start
      iqini=2
      if(bzcase()==2) iqini=1
      iqend=nqibz+nq0i

c RS:
      write(ifile_rsmpi,*) "RS: irot = ", irot
      write(ifile_rsmpi,*) "RS: main loop start"

cccccccccccccccccccccccccccccc
      call getkeyvalue("GWinput","TestOnlyQ0P",onlyq0p,default=.false.)
      if (onlyq0p) then
        call RSMPI_Stop("TestOnlyQ0P=T is disabled!")
      endif
cccccccccccccccc
      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      if (noq0p) then
        call RSMPI_Stop("TestNoQ0P=T is disabled!")
      endif

c      if(noq0p) iqend=nqibz
      if(noq0p) call RSMPI_Stop("noq0p==T")
cccccccccccccccc
c      iqend=nqibz
c      if(test_omitq0p()) then
c        iqend=nqibz
c        print *,'iqend=',iqend
c      endif  
cccccccccccccccccccccccccccc

c      if(smbasis()) ifpomat = iopen('POmat',0,-1,0) !oct2005

c      do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
c RS: local index
      do 1100 kx_local = 1,nk_local_rotk(irot)
c RS: true index
        kx = ik_index_rotk(irot,kx_local)

        write(ifile_rsmpi,*) ' sxcf_RSMPI: goto loop kx=',kx

        
       ! write(*,'("1  begin k-cycle",$)') 
       !  call cputid(0)
c          write(*,*)'kx, ip, irot=',kx, ip,irot

        if( kx <= nqibz ) then
c          k  = kx
          kr = irk(kx,irot)     ! index for rotated k in the FBZ
          qibz_k= qibz(:,kx)
          if(kr/=0) qbz_kr= qbz (:,kr) !feb2006 
        else
c          k = 1  ! corresponds to q=0
c          kr= 1  ! corresponds to q=0
c          k = iqindx((/0d0,0d0,0d0/), ginv, qibz,nqibz)
c          kr= iqindx((/0d0,0d0,0d0/), ginv, qbz,  nqbz)
          kr=-99999 !for sanity check
          qibz_k= 0d0
          qbz_kr= 0d0
        endif
c        ngc = ngcni(k)  ! k-points in IBZ
c        print *, ' k ngc=',k,ngc
c        ngb = nbloch + ngcni(k)

        call readqg0('QGcou',qibz_k,ginv,  quu,ngc)
c        ngc = ngcni(k)  ! k-points in IBZ
        ngb = nbloch + ngc !oct2005

c --- exchange
        if(exchange) then
c RS: set pointer to the right place
          call set_vcoul_rsmpi(ifvcfpout,kx-iqini+1)

          read(ifvcfpout) nn    !oct2005
          allocate(vcoul(nn,nn))
          read(ifvcfpout) vcoul(1:nn,1:nn)
        endif



c- weight check for cycle or not.
        if( kx <= nqibz ) then
          if (kr == 0)    then
            if(exchange) deallocate(vcoul)

c RS: this cannot happen..
            call RSMPI_Stop("in sxcf_RSPMI,cycle(kx<=nqibz)! this cannot happen!")
            cycle 
          endif
          kount(kx,ip)= kount(kx,ip) + 1 ! count the no. times k
                                             ! appears in the 1st BZ
cccccccccccccccccccccccccccccccccccccccccccccccccc
c         print *,' irot,ip, k, kount in  =',irot, ip, k, kount(k,ip)
c          deallocate(vcoul)
c          cycle
c          print *,' kount out =',kount(k,ip)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
          if (kount(kx,ip) > nstar(kx)) then
            call RSMPI_Stop('sexc: too many stars')
          endif
        else
          if( wgt0(kx-nqibz,irot) == 0d0 ) then 
            if(exchange) deallocate(vcoul)
c RS: this cannot happen..
            call RSMPI_Stop("in sxcf_RSPMI,cycle(kx>nqibz)! this cannot happen!")
            cycle
          endif
        endif

c---test
        if(OnlyQ0P.and.kx<=nqibz) then
          call RSMPI_Stop("this cannot happen")
        endif
    
c qk = q - rk, rk is inside 1st BZ, not restricted to 
c the irreducible BZ
c      call vminv   (q,qbz(1,kr),3,qk)
        qk =  q - qbz_kr ! qbz(:,kr)

c bring q-rk into the 1st BZ ==> qfbz labelled by kp
c         call fbz2(qk, ginv,index_qbz,n_index_qbz, qbz,nqbz, !qbz=rk
c     o                qfbz,kp)

c        if(debug)  print *,'xxx2',kp,qk,n1,n2,n3,qbas,ginv
c read eigenvalues,b and hb at q-rk
c        call rwdd1   (ifev, kp,nband,ekq)
        call readeval(qk, isp, ekq)

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
c         irkp = ibzx(kp) kp??? ! IBZ corresponding to qk
c         do i=1,nband
c          if( wtetef(i,irkp) == 0d0 ) exit
c         enddo
c        if(sum(abs(wtetef(i:nband,irkp)))/=0d0) 
c     &   stop ' sxcf: wtetef sum err1'
c         nt0org = nocc (ekc,ef,.true.,nstatetot)
c         nt0 = i-1 + nctot
c        endif

        if(exchange) then 
          nbmax = nt0p-nctot
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          nbmax = 5
c          write (6,*)' test xxxx=',nbmax
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          write(ifile_rsmpi,*) ' sxcf: nbmax nctot nt0p =',nbmax,nctot,nt0p
        else
c-------------
          nbmax = nband
          nbmxe = nocc (ekc,ebmx,.true.,nstatetot)-nctot
          nbmax  = min(nband,nbmx,nbmxe)
c          write(ifile_rsmpi,*) ' nbmax=',nbmax
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
        nt   = nctot + nbmax    ! = nstate for the case of correlation
        ntp0 = ntq
        allocate( zzmel(nbloch,nt,ntq)) !     &      rmel    cmel(nbloch,nt,ntq))
        call psicb_v2 (icore,ncore,ntp0,iclass,
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
        call psi2b_v2 (nbmax, ntp0, iclass,
     i    dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i              cphikq, !rbkq,cbkq,rhbkq,chbkq, !occ    q-rk nband
     i              cphiq,  !rbq,cbq,rhbq,chbq,     !unocc  q    ntp0
     i              ppb,  !pdb,dpb,ddb,
     i              nlnmv,nlnmc,mdim,nctot,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0,
     d              natom,nclass,
     o              zzmel) ! rmel,cmel)

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
        if(debug) print *, ' sxcf_fal2: goto drvmelp'
c        call drvmelp( iq,   ntp0, ! q in FBZ
c     i                kp, nbmax,  ! q-rk
c     i                k,          ! k in IBZ for e-product basis



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do iq=1,nqibz
c        write(6,"('--goto drvmelp2 qibz ngc=',3f9.4,'  ',i5)")
c     &       qibz(1:3,iq),ngcni(iq)
c      enddo  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if (lwang) then
        call drvmelp2_wan(iko_ix,iko_fx,nwf,umn,
     i                 q,   ntp0,              ! q in FBZ
     i                 q-qbz_kr, nbmax,  ! q-rk
     i                 qibz_k,           ! k in IBZ for e-product basis
     i                 isp,ginv,   
     i       ngc,ngcmx,ngpmx,nband,itq, 
     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
     o       rmelt,cmelt) 
        else ! lwang
        call drvmelp2( q,   ntp0, ! q in FBZ
     i                 q-qbz_kr, nbmax,  ! q-rk
     i                 qibz_k,           ! k in IBZ for e-product basis
     i                 isp,ginv,   
     i       ngc,ngcmx,ngpmx,nband,itq, 
     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
     o       rmelt,cmelt) 
        endif ! lwang

c zmel(I,ij) = dcmplx(rmelt, cmelt) = < psi_i psi^*_j | M_I>
c oct2005
        if(nbcut/=0.and.(.not.exchange)) then  
          do it= nctot+1,nctot+min(nbcut,nbmax)
            rmelt(:, it,:) =0d0
            cmelt(:, it,:) =0d0
          enddo  
        endif    

c oct2005 This block is on the fly reading of Pomat from disk.
c      if(smbasis()) then !
cC- smoothed mixed basis !oct2005
cC This replace original zmelt with new zmelt based on smoothed mixed basis.
c       do
c         read(ifpomat) q_r,nn,no,iqx !readin reduction matrix pomat
c         allocate( pomat(nn,no) )
c         read(ifpomat) pomat
c         if( sum(abs(qibz_k-q_r))<1d-10 .and.kx <= nqibz ) then
c           write(6,*) 'ok find the section for give qibz_k'
c           exit
c         elseif (iqx == kx ) then 
c! It may be better to check agreement with q_r itself
c! But now I did not generated qibze in hsfp0.
c! iqx ordering is by hsfp0
c           exit 
c         endif
c         deallocate(pomat)
c       enddo
c     endif
      
! need pomatr, nnr, nor 
!oct2005
        if(smbasis()) then      !
c RS: smbasis=T is disabled..
          call RSMPI_Stop("sxcf(),smbasis==T")
c 1113 continue
        else
          nn=ngb
          no=ngb
        endif

        if( oncew() ) then
          write(ifile_rsmpi,"('ngb nn=',3i4)") ngb,nn,no
        endif

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
c      allocate( zmel (ngb, nctot+nbmax, ntp0), aaa(ngb))
c      zmel = dcmplx (rmelt,-cmelt)
c      ifix=9090
c      open (ifix,file='zvz_test')
c      print *,' ======= goto zvz test  ========='
c      do ib1=1,nctot+nbmax
c      do ib2=1,ntp0
c! zmelt*vcoul
c        zvz1= sum( dconjg(zmel(1:nbloch,ib1,ib2))
c     &          * matmul(vcoul(1:nbloch,1:nbloch),zmel(1:nbloch,ib1,ib2)) )
c        aaa= matmul(vcoul,zmel(:,ib1,ib2) )
c        zvz= sum( dconjg(zmel(:,ib1,ib2))*aaa(:) )
c        write(ifix,"(2i4,  2d13.5,'  ',2d13.5,'  ',2d13.5)") ib1,ib2,zvz ,zvz1, zvz-zvz1
c      enddo
c      enddo
c      close(ifix)
c      stop '--- zvz test end------------------'
ccccccccccccccccccccccccccccccccccccccccccc



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
cccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccc
        deallocate(zzmel)       !rmel,cmel)
c
        write(ifile_rsmpi,*) ' sxcf: goto wtt'
        write(ifile_rsmpi,"('sum of rmelt cmelt=',4d23.16)") sum( rmelt),sum(cmelt)

        if(bzcase()==2)then 
          if(kx<=nqibz) then
            wtt = wk(kr)
            if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
          elseif(kx>nqibz) then !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
            wtt= wgt0(kx-nqibz,irot)
          endif   
        else
          if(kx<= nqibz) then 
            wtt = wk(kr)        !         wtx = 1d0
          else  
            wtt = wk(1)*wgt0(kx-nqibz,irot) !       wtx = wgt0(kx-nqibz,irot)
            if(abs(wk(1)-1d0/dble(nqbz))>1d-10) then
              call RSMPI_Stop('sxcf:wk(1) inconsistent')
            endif
          endif
        endif  

c      if(debug) then
        write(ifile_rsmpi,"(' kx wtt=',i4,f12.8)") kx,wtt
c      endif

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
ccc --- screened exchange case
c          if(screen) then
c            ix = 1
c            nrec=(kx-2)*nw+ix
c            if(bzcase()==2) nrec= (kx-1)*nw+ix
c            read(ifrcw,rec=nrec) zw  ! Readin W(0) - v 
c            vcoul = vcoul + zw(1:ngb,1:ngb) !c  screen test
c          endif
c
          allocate( zmel (ngb, nctot+nbmax, ntp0),vcoult(1:ngb,1:ngb))
          zmel  = dcmplx (rmelt,cmelt)
          vcoult= transpose(vcoul)
c     
          allocate( z1p(ngb,nctot+nbmax,ntp0), w3p( nctot+nbmax,ntp0))
ctm, 070621
          if (loffd) allocate( w4p(nctot+nbmax,ntp0,ntp0) )
          call matm( vcoult, zmel,
     &         z1p, ngb,ngb,(nctot+nbmax)*ntp0 )
          deallocate(vcoult, zmel)

c         print *,' nctot nbmax=',nctot,nbmax
          do itp = 1,ntp0
            do it  = 1,nctot+nbmax
              w3p(it,itp) = sum ( dreal(z1p(:,it,itp))*rmelt(:,it,itp)
     .             +      dimag(z1p(:,it,itp))*cmelt(:,it,itp) )
            enddo
          enddo
ctm, 070621
          if (loffd) then
            allocate( zmel (ngb, nctot+nbmax, ntp0) )
            zmel = dcmplx (rmelt,-cmelt)
            do itp = 1,ntp0
            do itp2= 1,ntp0
            do it  = 1,nctot+nbmax
              w4p(it,itp,itp2) = sum( z1p(:,it,itp)*zmel(:,it,itp2) )
            enddo
            enddo
            enddo
            deallocate(zmel)
          endif

cccccccccccccccccccccccccccccccccccc
c test
c          print *,'nctot=',nctot
c          do it=1,nctot+nbmax
c            if(ekc(it)<-6d0) w3p(it,:)=0d0
c          enddo
c          do itp=1,ntp0
c             write (6,*)' ic eq-', it,eq(itq(itp),ip) 
c             if(eq(itq(itp),ip)<-6d0)  w3p(:,itp)=0d0
c          enddo
ccccccccccccccccccccccccccccccccccccc


c DEBUG write
          if(debug) then
            do  it  = 1,nctot+nbmax
              do  itp = 1,ntp0
                write(ifile_rsmpi,"(' w3p =',2i4,2d14.6)") it,itp,w3p(it,itp)
              enddo
            enddo
          endif

ctetra
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

c Write the Spectrum function for exchange May. 2001
c RS: this is disabled..
          if(ifexsp/=0) then
            call RSMPI_Stop("sxcf, ifexsp!=0")
          endif

c--- Correct weigts wfac for valence by esmr
          do it = nctot+1, nctot+nbmax
cgaussian
            wfac = wfacx(-1d99, ef, ekc(it), esmr)
c          if(debug) print *, ' it wfac=',it, wfac,ef, ekc(it), esmr
            write(ifile_rsmpi,*) ' it wfac=',it, wfac,ef, ekc(it), esmr
            w3p(it,1:ntp0) = wfac * w3p(it,1:ntp0)
ctm, 070621
            if (loffd) w4p(it,:,:) = wfac * w4p(it,:,:)
          enddo

          if (.not.tote) then
            do itp = 1,ntp0     !S[j=1,nbloch]  z1p(j,t,t') <B(rk,j) psi(q-rk,n) |psi(q,t')>
              zsec(iwini,itp,ip) = zsec(iwini,itp,ip) 
     &           - wtt * sum( w3p(:,itp) )
            enddo
ctm, 070621
            if (loffd) then
             do itp = 1,ntp0
             do itp2= 1,ntp0
              zsec2(itp,itp2,ip,iwini) = zsec2(itp,itp2,ip,iwini)
     &        - wtt * sum( w4p(:,itp,itp2) )
             enddo
             enddo
             deallocate(w4p)
            endif   
          else  
            do itp = 1,ntp0
cgaussian
              wfac = wfacx(-1d99, ef2, eq(itq(itp),ip), esmr2) !june-2002 takao
              w3p(1:nctot+nbmax,itp) = wfac * w3p(1:nctot+nbmax,itp)
              exxq = exxq - wtt * sum( w3p(:,itp) )
            enddo
          endif
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
cccccccccccccccccccccccc
c      do itp=1,ntp0
c        write(6,"(2i3,3d13.6)") itp,ip,zsec(iwini,itp,ip),sum( w3p(:,itp) )
c      enddo
cccccccccccccccccccccccc
          write(ifile_rsmpi,*) "goto dealloc z1p.."
          deallocate( z1p,w3p,rmelt,cmelt,vcoul)
c        if(debug) write(6,*) "end of ex section"
          write(ifile_rsmpi,*) "end of ex section"
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

c      if(debug) print *,' end of zmel'
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
ctm, 070621
        if (loffd) allocate(zzwz0(nstate,ntp0,ntp0))
        ix = 1  - nw_i

cccccccccccccccccccccccccccccccccc
c        kxx=34
c        ix=22
c        nrec= (kxx-1)*nw+ix
c        print *,' xxx wvr3 nrec=',nrec,nblochpmx,kxx,ix,nw
c        read(ifrcw,rec=nrec) zw 
c        stop 'test end xxxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccc 


ccccccccccccccccccccccccccccccccccccccccccc faleev cccccccccccccccccccccc        
        nrec=(kx-2)*(nw-nw_i+1) +ix
        if(bzcase()==2) nrec= (kx-1)*(nw-nw_i+1) +ix
        write(ifile_rsmpi,*) ' wvr nrec=',nrec
        read(ifrcw,rec=nrec) zw ! direct access read Wc(0) = W(0) - v  
c1
c1        call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,   
c1     o    zwz0) ! zwz0 = zmel*(W(0)-v)*zmel
c2
        zwz0=0d0                
c        allocate(zmel1(ngb))
!$OMP parallel do private(itp,it,igb2,zz2)
        do itp=1,ntp0
          do it=1,nstate 
c          zmel1(:)=dconjg(zmel(:,it,itp))
            do igb2=2,ngb
c            zz2 = sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) )
              zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )

              zwz0(it,itp) = zwz0(it,itp)+zz2*zmel(igb2,it,itp)*2+
     &           dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp) 
            enddo               !igb2
            zwz0(it,itp) = zwz0(it,itp)+
     &         dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
          enddo                 !it
        enddo                   !itp
!$OMP end parallel do
        zwz0 = dreal(zwz0)
c2end

          !deallocate(zmel1)
ctm, 070621
        if (loffd) then
        zzwz0=0d0                
!$OMP parallel do private(itp,it,igb2,zz2)
        do itp=1,ntp0
        do itp2=1,ntp0
        do it=1,nstate 
c          zmel1(:)=dconjg(zmel(:,it,itp))
          do igb2=2,ngb
c            zz2 = sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) )
            zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )

            zzwz0(it,itp,itp2) = zzwz0(it,itp,itp2)
     &       +zz2*zmel(igb2,it,itp2)*2+
     &       dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp2) 
          enddo !igb2
          zzwz0(it,itp,itp2) = zzwz0(it,itp,itp2)+
     &    dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp2)
        enddo !it
        enddo !itp2
        enddo !itp
!$OMP end parallel do
        endif ! loffd
cccccccccccccccccccccccc old variant ccccccccccccccccccccccccccccccccc
c        read(ifrcw,rec=((kx-2)*nw+ix)) zw  ! direct access read Wc(0) = W(0) - v  
c        call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb,   
c     o    zwz0) ! zwz0 = zmel*(W(0)-v)*zmel
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c COH term test ----- The sum of the all states for zwz00 gives the delta function.
c RS: disabled..
        if(cohtest) then
          call RSMPI_Stop("sxcf, cohtet=T")
        endif
c
        nx  = niw
        if(niw <1) call RSMPI_Stop(" sxcf:niw <1")
        allocate( zwz(niw,nstate,ntp0),  zwzi(nstate,ntp0) )
ctm, 070621
        if (loffd)
     &  allocate( zzwz(niw,nstate,ntp0,ntp0),  zzwzi(nstate,ntp0,ntp0) )
        if(screen) allocate(zwz00(nstate,ntp0))
ctm, 070621
        if (screen.and.loffd) allocate(zzwz00(nstate,ntp0,ntp0))
c   
       
        if(verbose()>50) write(ifile_rsmpi,'("6 before matzwz in ix cycle ",$)') 
        if(verbose()>50) call cputid(ifile_rsmpi)

**************************** BEGIN 1 ************************************************
        zwz=0d0                 !...........................faleev, temp
ctm, 070621
        if (loffd) zzwz = 0d0
       ! allocate(zmel1(ngb))
        do ix = 1,nx            ! imaginary frequency w'-loop
          nrec=(kx-2)*niw+ix
          if(bzcase()==2) nrec= (kx-1)*niw+ix
          write(ifile_rsmpi,*) ' wvi nrec=',nrec
          read(ifrcwi,rec=nrec) zw ! Readin W-v on imag axis 
c1          call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb, 
c1     o      zwz(ix,1:nstate,1:ntp0))
c2
!$OMP parallel do private(itp,it,igb2,zz2)
          do itp= 1,ntp0
            do it = 1,nstate 
c             zmel1(:)=dconjg(zmel(:,it,itp))
              do igb2 = 2,ngb
c               zz2 = sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) )
                zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )

                zwz(ix,it,itp)=zwz(ix,it,itp)+zz2*zmel(igb2,it,itp)*2+
     &             dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp) 
              enddo             !igb2
              zwz(ix,it,itp) = zwz(ix,it,itp)+
     &           dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp)
c2end          
            enddo               !it
          enddo                 !itp
!$OMP end parallel do
ctm, 070621
          if (loffd) then
!$OMP parallel do private(itp,it,igb2,zz2)
          do itp= 1,ntp0
          do itp2= 1,ntp0
          do it = 1,nstate 
c             zmel1(:)=dconjg(zmel(:,it,itp))
             do igb2 = 2,ngb
c               zz2 = sum( zmel1(1:igb2-1)*zw(1:igb2-1,igb2) )
               zz2 = sum( dconjg(zmel(1:igb2-1,it,itp))*zw(1:igb2-1,igb2) )

               zzwz(ix,it,itp,itp2)=zzwz(ix,it,itp,itp2)
     &        +zz2*zmel(igb2,it,itp2)*2+
     &         dconjg(zmel(igb2,it,itp))*zw(igb2,igb2)*zmel(igb2,it,itp2) 
             enddo !igb2
             zzwz(ix,it,itp,itp2) = zzwz(ix,it,itp,itp2)+
     &       dconjg(zmel(1,it,itp))*zw(1,1)*zmel(1,it,itp2)
c2end          
           enddo !it
           enddo!itp2
           enddo!itp
!$OMP end parallel do
         endif ! loffd
!
        enddo                   !ix  
        if(verbose()>50) write(ifile_rsmpi,'("6.1 before matzwz in ix cycle ",$)') 
        if(verbose()>50) call cputid(ifile_rsmpi)
         !deallocate(zmel1)
******************************* END 1 ***********************************************
        zwz=dreal(zwz)



c********************************* BEGIN 3 ***************************************
c       if(.false.) then
c        do ix = 1,nx     ! imaginary frequency w'-loop         !********* faleev
c         nrec=(kx-2)*niw+ix
c         if(bzcase()==2) nrec= (kx-1)*niw+ix
c         if(debug) print *,' wvi2 nrec=',nrec
c         read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis 
c
cc zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> 
cc                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
cc        do itp = 1,ntp0
cc        do  it = 1,nstate
cc          zwz(ix,it,itp) = sum( 
cc     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
cc        enddo
cc        enddo
c          call matzwz( zw(1:ngb,1:ngb), zmel, ntp0,nstate,ngb, !********* faleev
c     o      zwz(ix,1:nstate,1:ntp0))
c        enddo !ix   
c        endif                        !********* faleev
cccccccccccccccccccccccccccccccccccccccc
cc zwz is the diagonal term and real.
cc      write(6,*) ' img sumzwz=',sum(abs(imag(zwz(ix,1:nstate,1:ntp0))))
cccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccc
cc      write(6,*) ' sumzwz=',sum(zmel)
cc     & ,sum(zw(1:ngb,1:ngb)),sum(zwz(ix,1:nstate,1:ntp0))
cccccccccccccccccccccccccccccccccc
c************************************ END 3 *************************************
c        
c        if(verbose()>50) write(*,'("7 after matzwz in ix cycle ",$)') 
c        if(verbose()>50) call cputid(0)

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
ctm, 070621
          if (loffd) then
           zzwz00 = zzwz0
           zzwz0  = 0d0      
           do ix = 1,nx
           zzwz(ix,:,:,:)=zzwz(ix,:,:,:) - zzwz00
           enddo 
          endif
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
ccccccccccccccccccccccccccccccccccccc
ctest imx
c        testimx=.true.
c        if(testimx) zwzi=0d0
ccccccccccccccccccccccccccccccccccccc


        if(ua_auto) then 
          allocate(uaa(nstate,ntq))
          do itp = 1,ntq
            do  it = 1,nstate
              ratio = abs(zwz(niw,it,itp)/zwz0(it,itp))
              call gen_uaa(ratio,freqx(niw),  uaa(it,itp))
              if(verbose()>45) then
                write(ifile_rsmpi,"(' it itp uaa=',2i4,12f8.4)") it,itp,uaa(it,itp)
              elseif(verbose()>40.and.mod(it,10)==1.and.mod(itp,10)==1) then
                write(ifile_rsmpi,"(' it itp uaa=', 2i4,12f8.4)") it,itp,uaa(it,itp)
              endif
            enddo
          enddo
        endif

        do      iw = iwini,iwend
c frequency integration along the imaginary axis, s. wint.f
c for each e(q-rk,n) and w in SEc(qt,w)
c        do     itp = 1,ntq
        do     itp = itini,itfin
           if (itini .gt. itfin) cycle
            do      it = 1,nstate
              we =.5d0*( omega(itp,iw) -ekc(it)) != .5d0*( eq(itq(itp),ip)+2d0*(dble(iw)-shtw)*deltaw-ekc(it))
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ctest imx
c         if(testimx) then 
c           if(ef-omega(itp,iw)>0) then
c             if(  +0.001<-we .and. ekc(it)< ef) cycle
c           endif
c           if(ef-omega(itp,iw)<=0) then
c             if(  -0.001>-we .and. ekc(it)> ef) cycle
c           endif
c         endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccc
              if(verbose()>50) then
                do  ix = 1,niw
                  ratio  = abs(zwz(ix,it,itp)/zwz0(it,itp))
                  freqw1 = (1d0 - freqx(ix))/ freqx(ix)
                  ua2_(ix) = sqrt(- 1d0/freqw1*log(ratio))
c            write(6,"(' ix freqw1=',i4,f10.4)") ix,freqw1
                enddo
                write(ifile_rsmpi,"(' sxcf_fal2:
     & ua=sqrt(1/w1*log(v0/v1))= ',12f8.4)")ua2_(1:niw)
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
              if(ua_auto) then
                ua_ = .5d0*uaa(it,itp)
                call gen_expa(niw,freqx,ua_,  expa_)
              endif

c----------------------
              esmrx = esmr
              if(it <= nctot) esmrx = 0d0
              if(GaussSmear()) then
                zwzi(it,itp) = 
c Gaussian smearing
     &               wintzsg ( zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua_,expa_,we,nx, esmrx)
              else
                zwzi(it,itp) = 
     &              wintzav ( zwz(1,it,itp),zwz0(it,itp),freqx,wx,ua_,expa_,we,nx, esmrx)
              endif 
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
            enddo
          enddo

c sum over both occupied and unoccupied states and multiply by weight
c         do     itp = 1,ntq 
        do     itp = itini,itfin
           if (itini .gt. itfin) cycle
             zsec(iw,itp,ip)  = zsec(iw,itp,ip) + wtt*sum(zwzi(:,itp))  
        enddo

c end of SEc w-loop
        enddo

ctm, 070621
      if (loffd) then
       do iw = iwini,iwend
c        do     itp = 1,ntq
        do     itp = itini,itfin
           if (itini .gt. itfin) cycle
        do     itp2 = itini,itfin
        do      it = 1,nstate
          we =.5d0*( omega(itp,iw) -ekc(it)) != .5d0*( eq(itq(itp),ip)+2d0*(dble(iw)-shtw)*deltaw-ekc(it))
!          we =.5d0*( ef - 1.0d-4 -ekc(it)) ! omega=Ef for off-diagonal SE(w)
ctm, 070702
c          if (itp.eq.itp2) we =.5d0*( omega(itp,0) -ekc(it))
          esmrx = esmr
          if(it <= nctot) esmrx = 0d0
          if(GaussSmear()) then
            zzwzi(it,itp,itp2) = 
c Gaussian smearing
     &      wintzsg ( zzwz(1,it,itp,itp2),zzwz0(it,itp,itp2),freqx,wx,
     &                ua_,expa_,we,nx, esmrx)
          else
            zzwzi(it,itp,itp2) = 
     &      wintzav ( zzwz(1,it,itp,itp2),zzwz0(it,itp,itp2),freqx,wx,
     &      ua_,expa_,we,nx, esmrx)
          endif 
        enddo
        enddo
        enddo

c sum over both occupied and unoccupied states and multiply by weight
c        do     itp = 1,ntq 
        do     itp = itini,itfin
           if (itini .gt. itfin) cycle
        do     itp2 = itini,itfin
          zsec2(itp,itp2,ip,iw)  = zsec2(itp,itp2,ip,iw)
     &                    + wtt*sum(zzwzi(:,itp,itp2))  
        enddo
        enddo
       enddo ! iw 
      endif ! loffd

c       if(debug) then
        write(ifile_rsmpi,*)' ntq nstate sum(zwzi)=',ntq,nstate,sum(zwzi)
        write(ifile_rsmpi,*) ' ntq nstate sum(zwz )=',ntq,nstate,sum(zwz)
        do itp = 1,ntq
          write(ifile_rsmpi,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
        enddo
c     endif
        deallocate(zwz,zwz0,zwzi)
ctm, 070621
        if (loffd) deallocate(zzwz,zzwz0,zzwzi)

        if(ua_auto) deallocate(uaa)
ccccccccccccccccccccccccccccccccccccc
c        if(onlyimagaxis()) then
c          print *,' sxcf_fal2: onlyimagaxis cycle'
c          deallocate(freq_r,zmel,zw,zmel1)
c         goto 1100
c        endif
c        if(testimx) then
c          print *,' sxcf_fal2: testimx cycle'
c          deallocate(freq_r,zmel,zw,zmel1)
c         goto 1100
c        endif
ccccccccccccccccccccccccccccccccccccc
c        do iw=nw_i,nw
c        write(6,*) iw,freq_r(iw)
c        enddo
c        stop 'cccccccccccccccccccccccccccccccccccccccccccc'




c===============================================================================
c contribution to SEc(qt,w) from the poles of G
c===============================================================================
        if(debug) print *,' go to poles'

C... Check freq_r for nw_i case. 
!    We assume freq_r(i) == -freq_r(-i) in this code. feb2006
        if(nw_i/=0) then
          if(nw/= -nw_i) call RSMPI_Stop("sxcf_fal2z: nw/=-nw_i")
          if(freq_r(0)/=0d0) call RSMPI_Stop("sxcf_fal2z: freq_r(0)/=0")
          if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0) 
     &         call RSMPI_Stop("sxcf_fal2z: freq_r /= -freq_r")
        endif

c---------------------------------------
c maximum ixs reqired.
c---------------------------------------
c      print *,' ekc at nt0p nt0m+1=', ekc(nt0p),ekc(nt0m+1)
c      print *,'  nt0p nt0m+1=', nt0p, nt0m+1
        ixsmx =0
        ixsmin=0
        do 3001 iw  = iwini,iwend
c        do 3001 itp = 1,ntq
        do 3001 itp = itini,itfin
            if (itini .gt. itfin) cycle
            esec  = omega(itp, iw) ! eq(itq(itp),ip) + 2.d0*(dble(iw)-shtw)*deltaw
            if (esec <= ef) then
              do  it = 1,nt0p   ! nt0p corresponds to efp 
c     we   = abs( 0.5d0*(ekc(it) - esec) )
                esmrx = esmr
                if(it<=nctot) esmrx = 0d0
cgaussian
                wfac = wfacx(esec,ef, ekc(it),esmrx)
                if(GaussSmear()) then
                  if(wfac<wfaccut) cycle
                  we = ( .5d0*(weavx(esec,ef,ekc(it),esmr)-esec))
                                ! we should be positive W(we)
                else
                  if(wfac==0d0) cycle
                  we = max( .5d0*(ekc(it) - esec), 0d0) !feb2006 ---enforcee is positive
!               we = ( .5d0*(ekc(it) - esec) )       !this should give little effects 
                endif

c             if(it<=nctot) then
c               if(wfac>wfaccut) then
c                 print *,' debug info aaa1:',wfac,wfaccut,it,nctot
c                 stop "sxcf: it<=nctot.and.wfac/=0"
c               endif  
c             endif

c RS:
c                do iwp  = 1,nw+1 ! feb2006
                do iwp  = 1,nw ! feb2006
                  ixs=iwp       ! ixs = iwp= iw+1
                  if(freq_r(iwp) > we) exit
                enddo           !sf 23May02

                if(ixs>ixsmx) then
                  wexx  = we
                  ixsmx = ixs
c------
                  if (ixsmx+1 > nw) then 
                    write (ifile_rsmpi,*) ' nw',nw
                    write (ifile_rsmpi,*) ' wexx, dw ',wexx,dw
                    write (ifile_rsmpi,*) ' esec ekc(it) ef ', esec,ekc(it),ef
                    call RSMPI_Stop(' sxcf 1112: |w-e| > max(w)')
                  endif
c------
                endif  
              enddo
            else                ! esec > ef
              do  it = nt0m+1,nstate ! nt0m corresponds to efm 
c             we   = abs( 0.5d0*(esec - ekc(it)) ) 
cgaussian 
                wfac = wfacx(ef,esec, ekc(it),esmr)
                if(GaussSmear()) then
                  if(wfac<wfaccut) cycle
                  we = .5d0*(weavx(ef,esec, ekc(it),esmr)-esec)
               !this we is negtive definite. 
                else            !abs is not necessary feb2006
                  if(wfac==0d0) cycle 
                  we = min( .5d0*(ekc(it)-esec), 0d0 ) 
               ! this we is negative definite
                endif

c                do iwm = 1,nw+1 !feb2006
                do iwm = 1,nw !feb2006
                  ixs=iwm       !ixs=iwp
c              write(6,"(' iwp freq_r we=',i5,f13.5,f13.5)") iwp,freq_r(iwp),we
                  if( abs(freq_r(iwm))>abs(we) ) exit
                enddo           !we assume freq_r(iw) = -freq_r(-iw) feb2006
                if( ixs>ixsmin) then
                  wexx  = we
                  ixsmin = abs(ixs)
                  if (ixsmin+1 > nw) then 
                    write (ifile_rsmpi,*) ' nw_i ixsmin',nw_i, ixsmin
                    write (ifile_rsmpi,*) ' wexx, dw ',wexx,dw
                    write (ifile_rsmpi,*) ' esec ekc(it) ef ', esec,ekc(it),ef
                    call RSMPI_Stop(' sxcf 222: |w-e| < min(w)')
                  endif
                endif
              enddo     
            endif
 3001     continue              !end of SEc w and qt -loop 
          if(nw_i==0) then
            nwxi = 0
            nwx  = max(ixsmx+1,ixsmin+1)
          else  
            nwxi = -ixsmin-1
            nwx  =  ixsmx+1
          endif
          if (nwx > nw   ) then 
            call RSMPI_Stop(' sxcf nwx check : |w-e| > max(w)')
          endif
          if (nwxi < nw_i) then 
            call RSMPI_Stop(' sxcf nwxi check: |w-e| > max(w)')
          endif
c      if(debug) print *,'nw, nwx=',nw,nwx
          write(ifile_rsmpi,*) 'nw, nwx=',nw,nwx


C... Find nt_max ------------------------------------
          nt_max=nt0p           !initial nt_max
          do 4001 iw  = iwini,iwend
c          do 4001 itp = 1,ntq
          do 4001 itp = itini,itfin
           if (itini .gt. itfin) cycle
c...........................................
c loop over occupied states (q-k,n); w.le.ef
c...........................................
              esec = omega(itp,iw) !eq(itq(itp),ip) + 2.d0*(dble(iw)-shtw)*deltaw
              if (esec > ef) then
c.............................................
c loop over unoccupied states (q-k,n); w.gt.ef
c.............................................
                do  it = nt0m+1,nstate ! nt0m corresponds to efm 
cgaussian
c           we   = abs( 0.5d0*(esec - ekc(it)) )
                  wfac = wfacx (ef,esec, ekc(it),esmr)
                  if( (GaussSmear().and.wfac>wfaccut)
     &                 .or.(.not.GaussSmear().and.wfac/=0d0)) then           
                    if (it > nt_max) nt_max=it ! nt_max is  unocc. state 
                  endif         ! that ekc(it>nt_max)-omega > 0
                enddo
              endif
 4001       continue            !end of  w and qt -loop 
c        if(nt0p>nt_max) then
c           print *, ' nt0p, nt_max=',nt0p,nt_max
c        endif
   
C... Set zw3 or zwz -----------------------------------
          zwz3mode=.true.
          if(iwend-iwini>2) then
            zwz3mode=.false.
          endif

          if(zwz3mode) then
            allocate( zw3(ngb,ngb,nwxi:nwx))
            do ix = nwxi,nwx    ! real frequency w'-loop
          ! write(*,*)'r_om', ix,nwx,nw_w
              nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
              if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
              write(ifile_rsmpi,*) ' wvr3 nrec=',nrec,nblochpmx,kx,ix,nw
              read(ifrcw,rec=nrec) zw 
              zw3(1:ngb,1:ngb,ix) = zw(1:ngb,1:ngb)  
            enddo
            deallocate(zw)
          else
            nstatex= max(ntp0,nt_max)
            allocate( zwz(nwxi:nwx,1:nstatex,ntp0) )
ctm, 070621
            if (loffd) allocate( zzwz(nwxi:nwx,1:nstatex,ntp0,ntp0) )
            do      ix = nwxi,nwx        
              nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
              if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
              read(ifrcw,rec=nrec) zw ! Readin (W-v)(k,w')(i,j) at k and w' on imag axis 
c zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)> 
              call matzwz(zw(1:ngb,1:ngb), zmel(1:ngb,1:nstatex,1:ntp0), ntp0,nstatex,ngb,   
     o             zwz(ix,1:nstatex,1:ntp0))
! zmel (ngb, nstate, ntp0)
ctm, 070621
              if (loffd)
     &        call matzwz2(2,zw(1:ngb,1:ngb), zmel(1:ngb,1:nstatex,1:ntp0), 
     &        ntp0,nstatex,ngb,   
     o        zzwz(ix,1:nstatex,1:ntp0,1:ntp0))
             enddo
            deallocate(zmel)
            deallocate(zw)
          endif  

c---------------------------------------------
          if(screen) then
            if(zwz3mode) call RSMPI_Stop(' this mode is not implimented')
            do ix = nw_i,nwx
              zwz(ix,:,:)=zwz(ix,:,:) - zwz00      
            enddo
            deallocate(zwz00)
ctm, 070621
            if (loffd) then
             do ix = nw_i,nwx
               zzwz(ix,:,:,:)=zzwz(ix,:,:,:) - zzwz00      
             enddo
             deallocate(zzwz00)
            endif   
          endif

c--------------------------------------------
c loop over real frequencies of Wc(k,w')(i,j)
c--------------------------------------------

c-------------------------------
c loop over w and t in SEc(qt,w)
c-------------------------------
          if(verbose()>50) write(ifile_rsmpi,'("10 wfacx
     &  iw,itp,it cycles ",$)') 
          if(verbose()>50) call cputid(ifile_rsmpi)

          do 2001 iw  = iwini,iwend
c          do 2001 itp = 1,ntq
          do 2001 itp = itini,itfin
             if (itini .gt. itfin) cycle
c...........................................
c loop over occupied states (q-k,n); w.le.ef
c...........................................
              esec     = omega(itp,iw) !eq(itq(itp),ip) + 2.d0*(dble(iw)-shtw)*deltaw
              if (esec <= ef) then
                do  it = 1,nt0p ! nt0p corresponds to efp 
! energy window is   esec<  e  <ef
! each ekc has the width, \pm 0.5*esmr.
c           we   = abs( 0.5d0*(ekc(it) - esec) )
                  esmrx = esmr
                  if(it<=nctot) esmrx = 0d0
cgaussian
                  wfac = wfacx(esec,ef, ekc(it),esmrx)
                  if(GaussSmear()) then
                    if(wfac<wfaccut) cycle
                    we = ( .5d0*(weavx(esec,ef, ekc(it),esmr)- esec))
            !this should be positive
c            if(it<=nctot) then                                     !faleev
c               if(wfac>wfaccut) stop "sxcf: it<=nctot.and.wfac/=0"
c            endif
                  else
                    if(wfac==0d0) cycle
                    we   = max( 0.5d0*(ekc(it) - esec), 0d0 )
            ! this should be positive
c            if(it<=nctot) then                                     !faleev
c               if(wfac/=0) stop "sxcf:  it<=nctot.and.wfac/=0"
c            endif
                  endif
ccccccccccccccccccccccccccccccccccccc
ctest imx
c           if(testimx.and. we>0.001) cycle
ccccccccccccccccccccccccccccccccccccc
c                  write(ifile_rsmpi,"( ' xxx1',10d13.6)") esec,ef, ekc(it),wfac
c three-point interpolation for Wc(we)
                  do iwp = 1,nw 
                    ixs=iwp
                    if(freq_r(iwp)>we) exit
                  enddo    
                  if(ixs+1>nwx) then
                    write(ifile_rsmpi,*)' ixs+1=',ixs,we,nwx,dw
                    call RSMPI_Stop(' sxcf: ixs+1>nw xxx1')
                  endif
cccccccccccccccccccccccccc
c           write(6,*)' bbb22222222111 ixs=',ixs,freq_r(ixs),we
cccccccccccccccccccccccccc
                  if(zwz3mode) then
c            zmel1(:)=dconjg(zmel(:,it,itp)) 
                    zwz3=(0d0,0d0)
!$OMP parallel do private( ix,igb2,zz2)
                    do ix = ixs, ixs+2
                      do igb2=1,ngb
c              zz2 = sum(zmel1(1:ngb)*zw3(1:ngb,igb2,ix)  ) 
                        zz2 = sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,ix-1) )
                        zwz3(ix-ixs+1) = zwz3(ix-ixs+1)+zz2*zmel(igb2,it,itp)
                      enddo     !igb2
                    enddo       !ix  
!$OMP end parallel do
                    zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     .                   - wfac*wtt*alagr3z(we,freq_r(ixs-1),zwz3) !faleev
                  else
                    zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     .                   - wfac*wtt*alagr3z(we,freq_r(ixs-1),zwz(ixs-1,it,itp))
                  endif
c------------
                enddo
              else
c.............................................
c loop over unoccupied states (q-k,n); w.gt.ef
c.............................................
                do it = nt0m+1,nt_max ! nt0m corresponds to efm         !faleev,nstates->nt_max  
c           we   = abs( 0.5d0*(esec - ekc(it)) )
cgaussian
                  wfac = wfacx(ef,esec, ekc(it),esmr)
                  if(GaussSmear()) then
                    if(wfac<wfaccut) cycle
                    we  = abs( .5d0*(weavx(ef,esec, ekc(it),esmr)-esec) )
                                !    this should be negative   
                  else
                    if(wfac==0d0) cycle
                    we  = abs( .5d0*min(ekc(it)-esec,0d0) )
                  endif         !this should be negative  

ccccccccccccccccccccccccccccccccccccc
ctest imx
c           if(testimx.and. we > 0.001) cycle
ccccccccccccccccccccccccccccccccccccc

c                  write(ifile_rsmpi,"(' xxx2',10d13.6)") esec,ef, ekc(it),we,wfac
c three-point interpolation for Wc(we)
c          if (we .le. dw/2d0+dw**2/omg_c/4d0) then
c           ixs=0
c          else 
c           ixs=idint( omg_c/dw*( 
c     &           dsqrt(1d0+2d0/omg_c*(we-dw**2/8d0/omg_c))-1d0 )+.5d0 )
c          endif  
                  do iwp = 1,nw 
                    ixs=iwp
                    if(freq_r(iwp)>we) exit
                  enddo                
                  if(nw_i==0) then
                    if(ixs+1>nwx) call RSMPI_Stop(' sxcf: ixs+1>nwx xxx2')
                  else
                    write(ifile_rsmpi,*)" ixs nwxi=",ixs,nwxi,freq_r(ixs-1),we,freq_r(ixs)
                    if(ixs+1>abs(nwxi)) call RSMPI_Stop(' sxcf: ixs+1>nwx yyy2')
                  endif

         ! dd=0d0 !(we-freq(ixs))/dw
cccccccccccccccccccccccccc
c           write(6,*)'bbb 2222222233 ixs=',ixs,freq_r(ixs),we
cccccccccccccccccccccccccc


c--------------------------------------------------------
c if nw_i==-nw, use W(-we) ---> iir=-1    :Feb2006
                  iir=1
                  if(nw_i/=0) iir=-1
                  fffr = freq_r( iir*(ixs-1):iir*(ixs+1):iir )
ccccccccccccccccccccccccccccc
c          write(6,*) 'fffr=',fffr
ccccccccccccccccccccccccccccc
                  if(zwz3mode) then
c            zmel1(:)=dconjg(zmel(:,it,itp))
                    zwz3=(0d0,0d0)
!$OMP parallel do private(ix,igb2,zz2)
                    do ix=ixs, ixs+2
                      do igb2=1,ngb
c              zz2=sum(zmel1(1:ngb)*zw3(1:ngb,igb2,ix)  ) 
                        zz2=sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,iir*(ix-1))  )
                        zwz3(ix-ixs+1)=zwz3(ix-ixs+1)+zz2*zmel(igb2,it,itp)
                      enddo     !igb2
                    enddo       !ix  
!$OMP end parallel do 
                    zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     .                   + wfac*wtt*alagr3z(we,fffr,zwz3) !faleev
                  else
                    zwzz(1:3) = zwz(iir*(ixs-1):iir*(ixs+1):iir, it,itp)
                    zsec(iw,itp,ip) = zsec(iw,itp,ip) 
     &                 + wfac*wtt*alagr3z(we,fffr,zwzz)
                  endif
c--------------
                enddo
              endif
 2001     continue            !end of SEc w and qt -loop 
          write(ifile_rsmpi,*) ' end of do 2001'
!
ctm, 070621
      if (loffd) then
        do 2501 iw = iwini,iwend  
        do 2501 itp = itini,itfin
           if (itini .gt. itfin) cycle
        do 2501 itp2 = itini,itfin
c...........................................
c loop over occupied states (q-k,n); w.le.ef
c...........................................
         esec     = omega(itp,iw)              !eq(itq(itp),ip) + 2.d0*(dble(iw)-shtw)*deltaw
!         esec = ef-1.d-4
ctm, 070702
c         if (itp.eq.itp2) esec     = omega(itp,0)
         if (esec <= ef) then
          do  it = 1,nt0p            ! nt0p corresponds to efp 
! energy window is   esec<  e  <ef
! each ekc has the width, \pm 0.5*esmr.
c           we   = abs( 0.5d0*(ekc(it) - esec) )
           esmrx = esmr
           if(it<=nctot) esmrx = 0d0
cgaussian
           wfac = wfacx(esec,ef, ekc(it),esmrx)
           if(GaussSmear()) then
            if(wfac<wfaccut) cycle
            we = ( .5d0*(weavx(esec,ef, ekc(it),esmr)- esec))
            !this should be positive
c            if(it<=nctot) then                                     !faleev
c               if(wfac>wfaccut) stop "sxcf: it<=nctot.and.wfac/=0"
c            endif
           else
            if(wfac==0d0) cycle
            we   = max( 0.5d0*(ekc(it) - esec), 0d0 )
            ! this should be positive
c            if(it<=nctot) then                                     !faleev
c               if(wfac/=0) stop "sxcf:  it<=nctot.and.wfac/=0"
c            endif
           endif
ccccccccccccccccccccccccccccccccccccc
ctest imx
c           if(testimx.and. we>0.001) cycle
ccccccccccccccccccccccccccccccccccccc
           if(debug.and.Is_IO_Root_RSMPI())
     &      write(6,"( ' xxx1',10d13.6)") esec,ef, ekc(it),wfac
c three-point interpolation for Wc(we)
           do iwp = 1,nw 
            ixs=iwp
            if(freq_r(iwp)>we) exit
           enddo    
           if(ixs+1>nwx) then
            if (Is_IO_Root_RSMPI())print *,' ixs+1=',ixs,we,nwx,dw
            call RSMPI_Stop(' sxcf: ixs+1>nw xxx1')
           endif
cccccccccccccccccccccccccc
c           write(6,*)' bbb22222222111 ixs=',ixs,freq_r(ixs),we
cccccccccccccccccccccccccc
           if(zwz3mode) then
c            zmel1(:)=dconjg(zmel(:,it,itp)) 
            zwz3=(0d0,0d0)
!$OMP parallel do private( ix,igb2,zz2)
            do ix = ixs, ixs+2
            do igb2=1,ngb
c              zz2 = sum(zmel1(1:ngb)*zw3(1:ngb,igb2,ix)  ) 
              zz2 = sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,ix-1) )
              zwz3(ix-ixs+1) = zwz3(ix-ixs+1)+zz2*zmel(igb2,it,itp2)
            enddo !igb2
            enddo !ix  
ctm, 070811
            if (wfac*wtt.gt.1d-4) then
!$OMP end parallel do
            zsec2(itp,itp2,ip,iw) = zsec2(itp,itp2,ip,iw) 
     .      - wfac*wtt*alagr3z2(we,freq_r(ixs-1),zwz3,itp,itp2) !faleev
            endif
           else
            zsec2(itp,itp2,ip,iw) = zsec2(itp,itp2,ip,iw) 
     .      - wfac*wtt*alagr3z2(we,freq_r(ixs-1),zzwz(ixs-1,it,itp,itp2),itp,itp2)
           endif
c------------
          enddo
         else
c.............................................
c loop over unoccupied states (q-k,n); w.gt.ef
c.............................................
         do it = nt0m+1,nt_max  ! nt0m corresponds to efm         !faleev,nstates->nt_max  
c           we   = abs( 0.5d0*(esec - ekc(it)) )
cgaussian
           wfac = wfacx(ef,esec, ekc(it),esmr)
           if(GaussSmear()) then
             if(wfac<wfaccut) cycle
             we  = abs( .5d0*(weavx(ef,esec, ekc(it),esmr)-esec) )
                   !    this should be negative   
           else
             if(wfac==0d0) cycle
             we  = abs( .5d0*min(ekc(it)-esec,0d0) )
           endif   !this should be negative  

ccccccccccccccccccccccccccccccccccccc
ctest imx
c           if(testimx.and. we > 0.001) cycle
ccccccccccccccccccccccccccccccccccccc

c          write(6,"(' xxx2',10d13.6)") esec,ef, ekc(it),we,wfac
c three-point interpolation for Wc(we)
c          if (we .le. dw/2d0+dw**2/omg_c/4d0) then
c           ixs=0
c          else 
c           ixs=idint( omg_c/dw*( 
c     &           dsqrt(1d0+2d0/omg_c*(we-dw**2/8d0/omg_c))-1d0 )+.5d0 )
c          endif  
          do iwp = 1,nw 
             ixs=iwp
             if(freq_r(iwp)>we) exit
          enddo                
          if(nw_i==0) then
             if(ixs+1>nwx) call RSMPI_Stop(' sxcf: ixs+1>nwx xxx2')
          else
             if (Is_IO_Root_RSMPI())write(6,*)
     &        " ixs nwxi=",ixs,nwxi,freq_r(ixs-1),we,freq_r(ixs)
             if(ixs+1>abs(nwxi))call RSMPI_Stop(' sxcf: ixs+1>nwx yyy2')
          endif

         ! dd=0d0 !(we-freq(ixs))/dw
cccccccccccccccccccccccccc
c           write(6,*)'bbb 2222222233 ixs=',ixs,freq_r(ixs),we
cccccccccccccccccccccccccc


c--------------------------------------------------------
c if nw_i==-nw, use W(-we) ---> iir=-1    :Feb2006
          iir=1
          if(nw_i/=0) iir=-1
          fffr = freq_r( iir*(ixs-1):iir*(ixs+1):iir )
ccccccccccccccccccccccccccccc
c          write(6,*) 'fffr=',fffr
ccccccccccccccccccccccccccccc
          if(zwz3mode) then
c            zmel1(:)=dconjg(zmel(:,it,itp))
            zwz3=(0d0,0d0)
!$OMP parallel do private(ix,igb2,zz2)
            do ix=ixs, ixs+2
            do igb2=1,ngb
c              zz2=sum(zmel1(1:ngb)*zw3(1:ngb,igb2,ix)  ) 
              zz2=sum(dconjg(zmel(1:ngb,it,itp))*zw3(1:ngb,igb2,iir*(ix-1))  )
              zwz3(ix-ixs+1)=zwz3(ix-ixs+1)+zz2*zmel(igb2,it,itp2)
            enddo !igb2
            enddo !ix  
!$OMP end parallel do 
            zsec2(itp,itp2,ip,iw) = zsec2(itp,itp2,ip,iw) 
     .      + wfac*wtt*alagr3z2(we,fffr,zwz3,itp,itp2) !faleev
          else
            zwzz(1:3) = zzwz(iir*(ixs-1):iir*(ixs+1):iir, it,itp,itp2)
            zsec2(itp,itp2,ip,iw) = zsec2(itp,itp2,ip,iw) 
     &      + wfac*wtt*alagr3z2(we,fffr,zwzz,itp,itp2)
          endif
c--------------
         enddo
        endif
 2501   continue  !end of SEc w and qt -loop 
        if(debug.and.Is_IO_Root_RSMPI()) print *,' end of do 2501'
      endif !loffd
!
          write(ifile_rsmpi,'("11 after alagr3z iw,itp,it cycles ",$)') 
          call cputid(ifile_rsmpi)

          if(debug) then
            do itp = 1,ntq
              write(ifile_rsmpi,'(" zsec=",i3,6d15.7)') itp,zsec(iwini:iwini+2,itp,ip)
            enddo
          endif
          if (allocated(zmel1))  deallocate(zmel1) !,zwz)
          if(zwz3mode) then
            deallocate(zmel,zw3)
          else  
            deallocate(zwz)
            if (loffd) deallocate(zzwz)
         endif
 1100 continue                ! end of k-loop
c     if(smbasis()) isx = iclose("POmat")
      write(ifile_rsmpi,"(' sum(abs(zsec))=',d23.15)") sum(abs(zsec))

c      call rlse (iwr1)
ccccccccccccccccccccccccccccc
c      print *," end kount ip irot=", irot, ip, sum(kount(:,ip)) 
c       print *," sumkount 1=",irot,ip,sum(kount(:,:)) 
cccccccccccccccccccccccccccc
      if (tote) then 
        call RSMPI_Stop("sxcf_RSMPI,tote==T")
c         if (tetraex) then
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
      if (allocated(rmelt)) deallocate(rmelt)
      if (allocated(cmelt)) deallocate(cmelt)
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
      if (allocated(expikt)) deallocate(expikt)
      if (allocated(zmel1)) deallocate(zmel1)
      if (allocated(zmel3)) deallocate(zmel3)
      if (allocated(zw_)) deallocate(zw_)
      if (allocated(zwz2)) deallocate(zwz2)
      if (allocated(zw2)) deallocate(zw2)
      if (allocated(zmel2)) deallocate(zmel2) 
      if (allocated(zw3)) deallocate(zw3)
      if (allocated(uaa)) deallocate(uaa)
      end

c -------------------------------------------------------------------
      subroutine drvmelp2_wan(iko_ix,iko_fx,nwf,umn,
     i        q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx,nband,itp0, 
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       rmelt,cmelt) 
      use m_readqg
      use m_readeigen
c ! q1=q    ---> iq 1:ntp0 q-point
c ! q2=q-rk ---> kp 1:nt0  occupied                        
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c  zmelpl
c from drvmelp2, 081021 by TM
c-----------------------------------------------------
      implicit none
      real(8):: q(3),q_rk(3),qik(3),ginv(3,3)
      integer(4):: ngp1, ngp2, ngpmx,nqbz,nqibz,  !iq,kp,k,
     &      ngcmx ,nctot,nband,
     &      ntp0,nt0,nbloch,noccx,  itx, ngc,nnum,inum,ig1,ig2,igc,
     &             ngvecpB1(3,ngpmx),  
     &             ngvecpB2(3,ngpmx),  
     &             ngveccBr(3,ngcmx), itp0(ntp0),
     &     nadd(3),isp  !,ngpn(nqbz)
      complex(8),allocatable::  zmelpl(:,:,:),geigq(:,:)
      real(8):: qbas(3,3),shtv(3),qdiff(3),add(3)
     &         ,qibz(3,nqibz),qbz(3,nqbz),qbasinv(3,3),symope(3,3)
     &         ,rmel(nbloch,noccx,ntp0) 
     &         ,cmel(nbloch,noccx,ntp0) 
     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     & ,pi=3.1415926535897932D0 
!      complex(8) :: geigB  (ngpmx,nband,nqbz)
      complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
ctm, 081021
      integer(4)::iko_ix,iko_fx,nwf
      complex(8):: umn(iko_ix:iko_fx,nwf)
c
c      print *,' drvmelp',ngpn(iq),ngpn(kp),ngc
      if(verbose()>50) print *,' drvmelp2:'
c      ngp1 = ngpn(iq)  ! q        in FBZ 
c      ngp2 = ngpn(kp)  ! q-rk     in FBZ
      call readqg('QGpsi', q,   ginv, qt,   ngp1, ngvecpB1)
      call readqg('QGpsi', q_rk,ginv, q_rkt,ngp2, ngvecpB2)
ctm, 081021
!      call readgeig(q,    ngpmx, isp, qu1, geig1)
      call readgeig_mlw(q,    ngpmx, isp, qu1, geig1(:,1:nwf))
      call readgeig(q_rk, ngpmx, isp, qu2, geig2)
      if(sum(abs(qt-qu1))>1d-10) stop 'drvmelp2;qu1/=qu1x'
      if(sum(abs(q_rkt-qu2))>1d-10) stop 'drvmelp2;qu2/=qu2x'
      if(verbose()>=100) print *,' end of read geig '

c     ngc  = ngcni(k)  ! k-points in IBZ
c      allocate( geigq(ngp1,ntp0) )
c      do itx = 1,ntp0
c        geigq(1:ngp1,itx) = geigB(1:ngp1,itp0(itx), iq)
c      enddo

c     print *,' drvmelp x1',symope,qibz(1:3,k),qbz(:,iq),qbz(:,kp)
c      qdiff = matmul(symope,qibz(1:3,k)) - qbz(:,iq) + qbz(:,kp) 
      qdiff = matmul(symope,qik)  - qt + q_rkt 
                           ! rk    -q  +(q-rk) is not zero.
      add  = matmul(qbasinv,qdiff)
      nadd = idint( add + dsign(.5d0,add))
      if(sum(abs(add-nadd))>1d-7) then
        print *,' drvmelp qdiff',qdiff
        write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)') 
     &  add,nadd,abs(add-nadd)
        stop "drvmlp: abs(add-nadd))>1d-10"
      endif 

      allocate( zmelpl(ngc,nt0,ntp0) )

cccccccccccccccccccccccccccccccccccccccccccccccc
c      if(.false.) then !old version
c      call melpln( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           , ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           , ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
c     i       shtv, matmul(symope,qibz(1:3,k)), qbas, 
c     o       zmelpl)
c      else             !new version from Sep-10 2000
c       call melpln2( ngp1, ngvecpB  (:,1:ngp1, iq)
c     &           ,  ngp2, ngvecpB  (:,1:ngp2, kp)
c     &           ,  ngc,  ngveccBr (:,1:ngc , k ), nadd,
c     &       geigq, ntp0,                   ! q1=q    ---> iq ngp1 1:ntp0 q-point
c     &       geigB(1:ngp2,1:nt0, kp), nt0,  ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
c     i       shtv, matmul(symope,qibz(1:3,k)),qibz(1:3,k), symope, qbas, 
c     o       zmelpl)
c       endif
       call melpln2( ngp1, ngvecpB1
     &           ,   ngp2, ngvecpB2 
c     &           ,   ngc,  ngveccBr, nadd,
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,  ! q1=q    ---> iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,1:nt0), nt0,     ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
     i       shtv, matmul(symope,qik),qik, symope, qbas, 
     o       zmelpl)
c      deallocate(geig1,geig2)
      rmelt=0d0
      cmelt=0d0
      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = 
     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)  
      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =  
     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
c extensiton to nbloch+ngc   
      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )  
      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )  
      deallocate(zmelpl) 
      end



c------------------------------------------------------
c --------------------------------------------------------------------
