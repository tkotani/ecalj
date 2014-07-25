      subroutine wmatq(kount,irot,ef,ef2,esmr,esmr2,tr,
ctm, 070501
c     i          iatomp,rsite,nsp,isp, !ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     i          iatomp,
c     i          rws,irws,nrws,
     i          rws1,rws2,nrws1,nrws2,nrws,
     i          nsp,isp, !ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     i                  ifrcw,ifrcwi, 
     i                  qbas,ginv,
     i         qibz,qbz,wk,nstbz,wik,nstar,irk,  !koun,,iindxk
c
     i                  iclass,mdim,nlnmv,nlnmc,
     i                  icore,ncore,imdim,
     i                  ppb, !pdb,dpb,ddb,
     i                  freq_r,freqx,wx,expa,ua,dw, !deltaw,freq
     i                  ecore, 
c
     d                  nlmto,nqibz,nqbz,nctot,
c     i          index_qbz, n_index_qbz,  !jan2004
     d                  nl,nnc,nclass,natom,
     d                  nlnmx,mdimx,nbloch,ngrp,nw_i,nw,nrw,niw,niwx,nq,
c
c     &     nblochpmx ,ngpn,ngcni,ngpmx,ngcmx,geigB,ngvecpB,ngveccBr,
     &     nblochpmx ,ngpmx,ngcmx, !ngveccBr,!Jan2004
     &     wgt0,wqt,nq0i,q0i,symope,alat, shtv,nband, ifvcfpout,
c     &     shtw,  
     &     exchange, !tote,screen,cohtest, ifexsp, 
ctetra
ctetra     &     wtet,wtetef,
ctetra    &     ntqx,ibzx,tetraex, 
c
c     i omega,iwini,iwend,
     i     nbmx,ebmx, !takao 18June2003 
     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat 
     i     nwf,
     o          rw_w,cw_w,rw_iw,cw_iw)
c
c 2006 May Takashi Miyake, updated for new fpgw
c 2004 Sep Takashi Miyake, off-site W
c 2004 Jul Takashi Miyake,
c 2004 Apr Takashi Miyake, from sxcf_fal2.f
c
c 2001 Sep. esec=omega(itp,iw). Genral iw mode for exchange =F
c
c 2000 takao kotani. This sxcf is modified from sec.f F.Aryasetiawan.
c
c  exchange=T : Calculate the exchange self-energy
c  exchange=F : Calculate correlated part of the self-energy
c
c---- correlation case documents by ferdi.Aryasetiawan.  -----------------
c 92.02.24
c 93.10.18 from sec.f modified to take into account equivalent atoms
c
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
c
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
c----------------------------------------------------------------------
c      implicit real*8(a-h,o-z)
      use m_readqg
      use m_readeigen
      use keyvalue
      implicit none
      integer(4) :: ntq, natom,nqbz,nqibz,ngrp,nq,nw_i,nw,niw,
     &  nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,nstatex,
     &   irot,  iqisp,ikpisp,isp,nsp,  nlnmx, iq, ! idxk,
c     &  iwr1,iwr2,iwr3,iwr4,iwc1,iwc2,iwc3,iwc4
     &  ip, it,itp, it2, itp2, !ifcphi,    ! ifrb,ifcb,ifrhb,ifchb, 
     i  iiclass,mdim(*),
     i  ifrcw,ifrcwi, !iindxk,
     &  ifvcfpout,ndummy1,ndummy2,kx,kr,kr2,kr3,ngc,ngb,nbloch, !n1,n2,n3, k,
     &  kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0,
     &  nbmax,nclass,nl,nnc, nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, !iwp sf 21May02
     &  mdimx, nwx,niwx,
     &  iatomp(natom),      
     &  nstar(nqibz),irk(nqibz,ngrp),kount(nqibz),nwf !,iclose
c
      real(8) :: q(3),qbas(3*3),ginv(3*3),tr(3,natom),
     &  wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz),
     &  freqx(niw),wx(niw),expa(niw),
     &  eq(nband), 
     &   ekq(nband), ekc(nctot+nband),
     &  tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx,ua,
     &  dw,wtt,wexx,www,exx,exxq,weight
c      complex(8) :: zsec(-1:1,ntq,nq)
c      real(8)    ::  shtw  
c                       ! This shft is  to avoid some artificial resonance effects.
c                       ! shtw can be zero for esmr/=0 given by takao.                 
c
      integer(4):: ngpmx, ngcmx,  !ngcni(nqibz), !ngpn(nqbz),
     &   igc, !ngvecpB(3,ngpmx,nqbz),ngveccBr(3,ngcmx,nqibz),
     &  nadd(3)
      real(8) :: wgt0(nq0i,ngrp),wqt(nq0i),qk(3), !qfbz(3),
     &  qbasinv(3,3), det,qdiff(3),add(3),symope(3,3),
     &  qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot),
     &   ppb(1) !pdb(1),dpb(1),ddb(1)
      complex(8)::   alagr3z,wintz !geigB  (ngpmx,nband,nqbz),

c
      real(8),allocatable:: !rmel(:,:,:),cmel(:,:,:),
     &                      rmelt(:,:,:),cmelt(:,:,:), 
     &                      rmelt2(:,:,:),cmelt2(:,:,:), 
c     &                      rmelt3(:,:,:),cmelt3(:,:,:) 
     &                      rmelt3(:,:,:,:),cmelt3(:,:,:,:) 
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:),
     &     zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:),
     &     zmelt(:,:,:),zmelc(:,:,:,:)
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
c      integer(4) :: ibzx(nqbz)
c      real(8)    :: wtet  (nband,nqibz,1:ntqx), wtetef(nband,nqibz) 
c                ! where the last index is 3*itq+iw-1,itq=1,ntq,iw=-1,1
c     logical    :: tetraex
cccccczzcccccc

      complex(8) :: wintzav,wintzsg

      integer(4) :: ibl,iii,ivsumxxx,ifexsp ,iopen
      integer(4),save::ifzwz=-999

      integer(4) :: iwini, iwend, ia
c      real(8)    :: esec, omega(ntq, iwini:iwend)
      real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw),
     &           cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw),
     &           rw_iw(nwf,nwf,nwf,nwf,nrws,niw),
     &           cw_iw(nwf,nwf,nwf,nwf,nrws,niw)
      complex(8),allocatable:: expikt(:)
      complex(8):: img=(0d0,1d0)
ctakao
      complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)
     .           , cphiqtmp(nlmto,nband)

cccccccccccccccccccccccccccccccccccccccccccccccc faleev 2002
      integer(4) :: nt_max, igb1,igb2,iigb, nw_w
      complex(8),allocatable:: zmel1(:)
      complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)  
      complex(8), allocatable :: zwz2(:,:),zw2(:,:,:,:) !0 variant
      complex(8) ::  zz2 ,zwz3(3) 
      real(8) :: dd,omg_c,dw2
      real(8) :: freq_r(nw_i:nw)
      complex(8), allocatable :: zw3(:,:,:)
        

      real(8)::weavx,wfaccut=1d-10

      logical :: GaussSmear
      real(8) :: ebmx,ddw
      integer(4):: nbmx,nbmxe,nstatetot

c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4)::verbose,nstbz(nqbz),bzcase,iqini,iqend
      real(8):: wgtq0p

      integer(4):: iqindx,nrec,kxx
      real(8)::quu(3),qibz_k(3),qbz_kr(3)
      logical :: onlyQ0P, onlyimagaxis ,noq0p !,noq0p,test_omitq0p,

      logical ::zwz3mode
c      logical ::testimx=.false.

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

ctm
      integer(4) :: nqbz2,nwf2,iko_ix,iko_fx,iqtmp,ifmlw,nko,iqk
     .             ,ifi,in1,in2,imp,ilp,ii,jj,nrws,nrws1,nrws2
     .             ,ir1,ir2,ir3,ir,nrw
      real(8) :: norm2,qtmp(3),rws1(3,nrws1),rws2(3,nrws2),tmp
      complex(8) :: ztmp,expiqR1(nrws1),expiqR2
      complex(8),allocatable :: cnk(:,:,:),zmel2(:,:,:),zmel3(:,:,:)
      integer(4) :: itq(nwf)
      complex(8) :: weightc(nrws1)

c--------------------------------------------------------------------
      debug=.false.
      if(verbose()>=90) debug=.true.

c oct2005
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      nbcutc=nctot+nbcut

      tpi         = 8d0*datan(1.d0)   
c      iq         = idxk (q,qbz,nqbz) ! index for q
c      write(6,"(' iq q  =',i4,3f8.4)")iq,q
cccc      iq          = idxk (q,qbze,nqbze) ! index for q
      ekc(1:nctot)= ecore(1:nctot)  ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
      call dinv33(qbas,0,qbasinv,det) 

c work arrays for psi2br.f
      if(debug) print *, ' sxcf: 1'
      allocate(expikt(natom))

      if(bzcase()==1) then
      if(abs(sum(qibz(:,1)**2))/=0d0) stop ' sxcf assumes 1st qibz/=0 '
      if(abs(sum( qbz(:,1)**2))/=0d0) stop ' sxcf assumes 1st qbz /=0 '
      endif

      do it = 1,nwf
         itq(it) = it
      enddo

ctm debug
c        write(*,*)'isp,',isp,rw_w(1,1,1,1,1,0),rw_w(2,2,2,2,1,0)

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
      if(irk(kx,irot)/=0) kount(kx)= kount(kx) + 1 
c      kount(kx)= kount(kx) + 1 
      endif
      
c --- main loop start
      iqini=2
      if(bzcase()==2) iqini=1
      iqend=nqibz+nq0i

cccccccccccccccccccccccccccccc
      call getkeyvalue("GWinput","TestOnlyQ0P",onlyq0p,default=.false.)
      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      if (.not.noq0p)
     &  call getkeyvalue("GWinput","NoQ0P",noq0p,default=.false.)
      if (noq0p)write(*,*)'noq0p mode'
      if(noq0p) iqend=nqibz
c      iqend=nqibz
c      if(test_omitq0p()) then
c        iqend=nqibz
c        print *,'iqend=',iqend
c      endif  
cccccccccccccccccccccccccccc
      do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
c debug:
c      do 1100 kx = iqini,iqini !kx=1 corresponds to q=0 is omitted.
c        write (6,"(i3,'  out of ',i3,$)") kx,iqend
        if(debug)  print *, ' sxcf: goto loop kx=',kx
        
       ! write(*,'("1  begin k-cycle",$)') 
c         call cputid(0)
c          write(*,*)'kx, ip, irot=',kx, ip,irot

        if( kx <= nqibz ) then
c          k  = kx
          kr = irk(kx,irot) ! index for rotated k in the FBZ
          qibz_k= qibz(:,kx)
c          qbz_kr= qbz (:,kr)
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
        ngb = nbloch + ngc

c --- exchange
        if(exchange) then
c          allocate(vcoul(ngb,ngb))
c          read(ifvcfpout) vcoul(1:ngb,1:ngb)
          read(ifvcfpout) nn !oct2005
          allocate(vcoul(nn,nn))
          read(ifvcfpout) vcoul(1:nn,1:nn)
        endif

c- weight check for cycle or not.
        if( kx <= nqibz ) then
          if (kr == 0)    then
c            stop 'wmat: kr=0'
            if(exchange) deallocate(vcoul)
            cycle 
          endif
          kount(kx)= kount(kx) + 1     ! count the no. times k
                                             ! appears in the 1st BZ
cccccccccccccccccccccccccccccccccccccccccccccccccc
c         print *,' irot,ip, k, kount in  =',irot, ip, k, kount(k,ip)
c          deallocate(vcoul)
c          cycle
c          print *,' kount out =',kount(k)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c          if (kount(kx) > nstar(kx)) stop 'sexc: too many stars'
          if (kount(kx) > nstar(kx)) stop 'wmat: kount > nkstar'
c          if (kount(kx) > 1) stop 'wmat: kount > 1'
        else
          if( wgt0(kx-nqibz,irot) == 0d0 ) then 
            if(exchange) deallocate(vcoul)
            cycle
          endif
        endif

c---test
      if(OnlyQ0P.and.kx<=nqibz) then
        if(exchange) deallocate(vcoul)
        cycle
      endif
    
c phase factor for off-site W
      do ir1=1,nrws1
         expiqR1(ir1) = exp(-img*tpi* sum(qbz_kr(:)*rws1(:,ir1)))
      enddo
c      tmp =  sum(qbz_kr*rsite)
c      write(*,*)'tmp',tmp
c      write(*,*)'qbz_kr',qbz_kr
c      write(*,*)'rsite',rsite

c loop over FBZ
      allocate( rmelt3(ngb,nwf,nwf,nrws2),cmelt3(ngb,nwf,nwf,nrws2))
      rmelt3 = 0d0
      cmelt3 = 0d0

      do iq = 1,nqbz
ctm debug      do iq = 1,0

      q(:) = qbz(:,iq)

c read eval and eigen
      call readeval (qbz(:,iq), isp, eq)
      call readcphi (qbz(:,iq), nlmto,isp, quu, cphiq)
ctm debug
c      write(7400,*)iq,nlmto,nwf
c      write(7400,*)cphiq
      if(debug) print *, ' sxcf: 2'

c qk = q - rk, rk is inside 1st BZ, not restricted to 
c the irreducible BZ
c      call vminv   (q,qbz(1,kr),3,qk)
        qk =  q - qbz_kr ! qbz(:,kr)
c        iqk = iqindx(qk, ginv,qbz,nqbz)

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
        nbmax = nwf
        nt   = nctot + nbmax ! = nstate for the case of correlation
        ntp0 = nwf
        allocate( zzmel(nbloch,nt,ntp0))
        zzmel = 0d0
c>> valence        
        call psi2b_v2 (nbmax, ntp0, iclass,
     i    dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i              cphikq, !rbkq,cbkq,rhbkq,chbkq, !  q-rk nko
     i              cphiq,  !rbq,cbq,rhbq,chbq,     !  q    nko
     i              ppb,  !pdb,dpb,ddb,
     i              nlnmv,nlnmc,mdim,nctot,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0,
     d              natom,nclass,
     o              zzmel) ! rmel,cmel)

        if(debug) write(6,"('sum of zmel abszmel=',4d23.16)")
     &                 sum(zzmel),sum(abs(zzmel) )

c-plane
c        allocate( rmelt(ngb, nctot+nbmax, ntp0), ! nstate= nctot+nband
c     &            cmelt(ngb, nctot+nbmax, ntp0))
        allocate( zmelt(ngb, nctot+nbmax, ntp0) )
        if(debug) print *, ' sxcf_fal2: goto drvmelp'
c        call drvmelp( iq,   ntp0, ! q in FBZ
c     i                kp, nbmax,  ! q-rk
c     i                k,          ! k in IBZ for e-product basis


c        call drvmelp2( q,   ntp0, ! q in FBZ
c     i                 q-qbz_kr, nbmax,  ! q-rk
c     i                 qibz_k,           ! k in IBZ for e-product basis
c     i                 isp,ginv,   
c     i       ngc,ngcmx,ngpmx,nband,itq, 
c     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
c     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
c     o       rmelt,cmelt) 

        call drvmelp3( q,   ntp0, ! q in FBZ
     i                 q-qbz_kr, nbmax,  ! q-rk
     i                 qibz_k,           ! k in IBZ for e-product basis
     i                 isp,ginv,   
     i       ngc,ngcmx,ngpmx,nband,itq, 
     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
     o       zmelt) 

        deallocate(zzmel) !rmel,cmel)
c
        if(debug) print *, ' sxcf: goto wtt'
        if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum( rmelt),sum(cmelt)

        do ir2=1,nrws2
         expiqR2 = exp( img*tpi* sum(q(:)*rws2(:,ir2)))
c rmelt3,cmelt3
         rmelt3(:,:,:,ir2) = rmelt3(:,:,:,ir2)
     &                     + wk(iq) * dreal(zmelt(:,:,:)*expiqR2)
         cmelt3(:,:,:,ir2) = cmelt3(:,:,:,ir2)
     &                     + wk(iq) * dimag(zmelt(:,:,:)*expiqR2)
        enddo  ! ir2
        deallocate(zmelt)

c end of loop over FBZ
      end do

ctm debug
c      do jj = 1,nwf
c        norm2 = 0d0
c        do ii = 1,ngb
c           norm2 = norm2 + (rmelt3(ii,jj,jj)**2+cmelt3(ii,jj,jj)**2)
c        enddo
c        write(*,*)'norm check',jj,norm2
c      enddo

        if(bzcase()==2)then 
          if(kx<=nqibz) then
            wtt = wk(kr)
            if(nstbz(kr)/=0) wtt = wk(kr)*(1d0-wgtq0p()/nstbz(kr))
ctm nqibz            wtt = wik(kx)*0.5d0
ctm nqibz            if(nstbz(kr)/=0) wtt = wtt*(1d0-wgtq0p()/nstbz(kr))
          elseif(kx>nqibz) then  !   wtx= wgt0(kx-nqibz,irot)/dble(nqbz)
            wtt= wgt0(kx-nqibz,irot)
ctm nqibz            wtt= wqt(kx-nqibz)*wgtq0p()/dble(nqbz)
          endif   
        else
          if(kx<= nqibz) then 
            wtt = wk(kr)                   !         wtx = 1d0
ctm nqibz            wtt = wik(kx)*0.5d0             !
          else  
            wtt = wk(1)*wgt0(kx-nqibz,irot) ! wtx = wgt0(kx-nqibz,irot)
ctm nqibz            wtt = wk(1)*wqt(kx-nqibz)  !
            if(abs(wk(1)-1d0/dble(nqbz))>1d-10) stop 'sxcf:wk(1) inconsistent'
          endif
        endif  

        weight = wtt
ctm        weight = weight * wik(kx) * 0.5d0 / wk(kr)

        if(debug) then
           write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
        endif

ctm, phase factor arising from off-site W
c        weight = weight
c        weightc = weight
        do ir1=1,nrws1
           weightc(ir1) = weight*expiqR1(ir1)
        enddo

c--------------------------------------------------------
c --- bare Coulomb section ---
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

         if (debug) write(*,*) 'bare coulomb section begins'
c
         allocate(z1p(ngb,nwf,nwf))
         allocate(zmel1(ngb))
         allocate(zmel(ngb, nwf, nwf),vcoult(1:ngb,1:ngb))
         do ir2=1,nrws2
c> zmel(j,t,t') = <psi(q,t') | psi(q-rk,t) B(rk,j)>
          zmel  = dcmplx (rmelt3(:,:,:,ir2),cmelt3(:,:,:,ir2))
          vcoult= transpose(vcoul)
c
c          deallocate(vcoul)

c          allocate( z1p(ngb,nwf,nwf))
          call matm( vcoult, zmel,
     &           z1p, ngb,ngb,nwf*nwf )
c          deallocate(vcoult, zmel)
c          allocate(zmel1(ngb))

c         print *,' nctot nbmax=',nctot,nbmax


          do ir3=1,nrws2
c
          do itp2 = 1,nwf
          do it2  = 1,nwf
          do it   = 1,nwf
          do itp  = 1,nwf

            zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),-cmelt3(:,it,itp,ir3))
            ztmp = sum ( z1p(:,it2,itp2)*zmel1 )
            do ir1=1,nrws1
            ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
            rw_w(itp2,it2,it,itp,ir,0) = rw_w(itp2,it2,it,itp,ir,0)
     &                                 + real(ztmp*weightc(ir1))
            cw_w(itp2,it2,it,itp,ir,0) = cw_w(itp2,it2,it,itp,ir,0)
     &                                 + imag(ztmp*weightc(ir1))
            enddo ! ir1
          enddo
          enddo
          enddo
          enddo
c
          enddo ! ir3 
         enddo ! ir2
         deallocate(vcoul)
         deallocate(vcoult, zmel)
         deallocate(z1p,rmelt3,cmelt3,zmel1)

         if (debug) write(*,*) 'bare coulomb section finished'

c-- End of bare-Coulomb section --------------

      else
c--------------------------------------------------------------------------
c--- screening effect section----------------------------------------------
c--------------------------------------------------------------------------

c--------------------------------------------------------------
c S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> 
c                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n') |psi(q,t')>
c--------------------------------------------------------------


c--- The matrix elements zmel.
        allocate( zmelc (ngb, nwf, nwf,nrws2),
     c            zw (nblochpmx,nblochpmx),
     c            zw2(nwf,nwf,nwf,nwf) )
c> zmelc(j,t,t') = < B(rk,j) psi(q-rk,t) | psi(q,t') >
        zmelc = dcmplx (rmelt3,-cmelt3)
        deallocate(rmelt3,cmelt3)
        if(debug) print *,' end of zmel'


c====================================================================
c Wc(qt,w) along the imaginary axis
c====================================================================

c------------------------------------------------
c loop over w' = (1-x)/x, frequencies in Wc(k,w')
c {x} are gaussian points between (0,1)
c------------------------------------------------
c
        nx  = niw
c        if(niw <1) stop " sxcf:niw <1"
c
        do ix = 1,nx     ! imaginary frequency w'-loop
          nrec=(kx-2)*niw+ix
          if(bzcase()==2) nrec= (kx-1)*niw+ix
          read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis 
c          read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis 

c zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> 
c                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
c        do itp = 1,ntp0
c        do  it = 1,nstate
c          zwz(ix,it,itp) = sum( 
c     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
c        enddo
c        enddo

          do ir3=1,nrws2
          do ir2=1,nrws2
          call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), 
     i                                   zmelc(:,:,:,ir3), 
     i      nwf,nwf,ngb, 
     o      zw2)

          do ir1=1,nrws1
            ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
          rw_iw(:,:,:,:,ir,ix)  
     &  = rw_iw(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
          cw_iw(:,:,:,:,ir,ix)  
     &  = cw_iw(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
          enddo ! ir1
          enddo ! ir2
          enddo ! ir3

        enddo


c====================================================================
c Wc(qt,w) along the real axis
c====================================================================

        if(debug) print *,' go to poles'

c        do      ix = 0,nw                    ! real frequency w'-loop
        do      ix = 0,nrw                    ! real frequency w'-loop
          nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
          if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
          read(ifrcw,rec=nrec) zw 

c zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)> 
          do ir3=1,nrws2
          do ir2=1,nrws2
          call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), 
     i                                   zmelc(:,:,:,ir3), 
     i      nwf,nwf,ngb, 
     o      zw2)

          do ir1=1,nrws1
            ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
          rw_w(:,:,:,:,ir,ix)  
     &  = rw_w(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
          cw_w(:,:,:,:,ir,ix)  
     &  = cw_w(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
          enddo ! ir1
          enddo ! ir2
          enddo ! ir3

c debug: >
c          write(*,*)'nw,nrw=',nw,nrw
c          write(*,*)'ix=',ix
c          write(*,*)'nrws=',nrws
c          do i=1,5
c          write(*,*)i,rw_w(i,i,i,i,1,ix)
c          enddo
c          stop 'debug: 1'
c debug: <             

        enddo

        deallocate(zmelc,zw,zw2)

        if(debug) print *,' end of screen-if'

c end of if (exchange)
      endif


 1100 continue  ! end of k-loop


      return
      end

c -------------------------------------------------------------------
      subroutine matzwz2(zw,zmel, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),
     c              zwz(ntp0,nstate,nstate,ntp0)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntp0) )
      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
      do itp2 = 1,ntp0
      do  it2 = 1,nstate
      do  it  = 1,nstate
      do itp  = 1,ntp0
        zwz(itp,it,it2,itp2) 
     & = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
      enddo
      enddo
      enddo
      enddo
      deallocate(CC)
      end
c -------------------------------------------------------------------
      subroutine matzwz3(zw,zmel1,zmel2, ntp0,nstate,ngb, zwz)
      implicit none
      integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
      complex(8) :: zw(ngb,ngb),zmel1(ngb,nstate,ntp0),
     c                          zmel2(ngb,nstate,ntp0),
     c              zwz(ntp0,nstate,nstate,ntp0)
      complex(8), allocatable :: CC(:,:,:)
      allocate(CC(ngb,nstate,ntp0) )
      call matm(zw,zmel2,cc, ngb, ngb, nstate*ntp0)
      do itp2 = 1,ntp0
      do  it2 = 1,nstate
      do  it  = 1,nstate
      do itp  = 1,ntp0
        zwz(itp,it,it2,itp2) 
     & = sum( dconjg(zmel1(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
      enddo
      enddo
      enddo
      enddo
      deallocate(CC)
      end
c -------------------------------------------------------------------
      subroutine drvmelp3( q, ntp0, q_rk,nt0, qik, isp,ginv,
     i        ngc,ngcmx,ngpmx,nband,itp0, 
     i        symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       rmel, cmel, nbloch,noccx,nctot,
     o       zmelt) 
      use m_readqg
      use m_readeigen
c ! q1=q    ---> iq 1:ntp0 q-point
c ! q2=q-rk ---> kp 1:nt0  occupied                        
c- Contributions to matrix elements <psi |psi Plane> from interstitial plane wave.
c  zmelpl
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
c     &         ,rmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
c     &         ,cmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
     & ,pi=3.1415926535897932D0 
!      complex(8) :: geigB  (ngpmx,nband,nqbz)
      complex(8):: geig1(ngpmx,nband),geig2(ngpmx,nband)
     &            ,zmelt(1:nbloch+ngc,1:nctot+nt0,1:ntp0)
      real(8):: q_rkt(3),qt(3),qu1(3),qu2(3)
      integer(4)::verbose
c
c      print *,' drvmelp',ngpn(iq),ngpn(kp),ngc
      if(verbose()>50) print *,' drvmelp2:'
c      ngp1 = ngpn(iq)  ! q        in FBZ 
c      ngp2 = ngpn(kp)  ! q-rk     in FBZ
      call readqg('QGpsi', q,   ginv, qt,   ngp1, ngvecpB1)
      call readqg('QGpsi', q_rk,ginv, q_rkt,ngp2, ngvecpB2)
      call readgeig(q,    ngpmx, isp, qu1, geig1)
      call readgeig(q_rk, ngpmx, isp, qu2, geig2)
c debug
c      geig1=0d0
c      geig2=0d0

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
      if(sum(abs(add-nadd))>1d-7) stop "drvmlp: abs(add-nadd))>1d-10"
c     print *,' drvmelp qdiff',qdiff
c      write(6,'(" add nadd diff =",3f12.5,3i12,d13.6)') 
c     &  add,nadd,abs(add-nadd)

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
      if (nt0 .ne. ntp0) stop 'drvmelp3: nt0 .ne. ntp0'
       call melpln2( ngp1, ngvecpB1
     &           ,   ngp2, ngvecpB2 
c     &           ,   ngc,  ngveccBr, nadd,
     &           ,   ngc,  nadd,
     &       geig1(1:ngp1,itp0(1:ntp0)), ntp0,  ! q1=q    ---> iq ngp1 1:ntp0 q-point
     &       geig2(1:ngp2,itp0(1:nt0)), nt0,     ! q2=q-rk ---> kp ngp2 1:nt0  occupied                        
     i       shtv, matmul(symope,qik),qik, symope, qbas, 
     o       zmelpl)
c      deallocate(geig1,geig2)
c      rmelt=0d0
c      cmelt=0d0
      zmelt=0d0
c      rmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = 
c     &rmel (1:nbloch, 1:nctot+nt0, 1:ntp0)  
c      cmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) =  
c     &cmel (1:nbloch, 1:nctot+nt0, 1:ntp0)
       zmelt(1:nbloch, 1:nctot+nt0, 1:ntp0) = dcmplx
     &(rmel (1:nbloch, 1:nctot+nt0, 1:ntp0),  
     & cmel (1:nbloch, 1:nctot+nt0, 1:ntp0))
c extensiton to nbloch+ngc   
c      rmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c     & = dreal( zmelpl(1:ngc,         1:nt0,        1:ntp0) )  
c      cmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
c     & = dimag( zmelpl(1:ngc,         1:nt0,        1:ntp0) )  
       zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0)
     &=zmelpl(1:ngc,               1:nt0,            1:ntp0) 
      deallocate(zmelpl) 
      end



c------------------------------------------------------
