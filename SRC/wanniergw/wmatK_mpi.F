      subroutine wmatqk_mpi(kount,irot,ef,ef2,esmr,esmr2,tr,
!takao to fit to new hvccfp0.
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
c     i     nbmx,ebmx, !takao 18June2003 
     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat 
     i     nwf,
     o          rw_w,cw_w,rw_iw,cw_iw)
      use m_zmel,only: drvmelp3
      
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
      use m_readqg,only: readqg0
      use m_readeigen,only:readcphiw
      use m_keyvalue,only: getkeyvalue
      use m_read_bzdata,only: wklm

c RS: modules for MPI    
      use rsmpi
      use rsmpi_qkgroup
      use rsmpi_rotkindex
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
c     &   ekq(nband), ekc(nctot+nband),
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
     &  qbasinv(3,3),qdiff(3),add(3),symope(3,3),
     &  qxx(3),q0i(1:3,1:nq0i),shtv(3),alat,ecore(nctot),
     &   ppb(1) !pdb(1),dpb(1),ddb(1)
      real(8),allocatable:: !rmel(:,:,:),cmel(:,:,:),
     &                      rmelt(:,:,:),cmelt(:,:,:), 
     &                      rmelt2(:,:,:),cmelt2(:,:,:), 
c     &                      rmelt3(:,:,:),cmelt3(:,:,:) 
     &                      rmelt3(:,:,:,:),cmelt3(:,:,:,:) 
      complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:),
     &     zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:),
     &     zmelt(:,:,:),zmelc(:,:,:,:),zmelcc(:,:,:,:)
c for exchange --------------------
      logical :: exchange,screen,cohtest,tote
      real(8),allocatable:: 
     &   w1p(:,:,:),w2p(:,:,:)
      complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
   
c- debug write ---------------------
      logical :: debug=.false.

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
      real(8) :: ddw !ebmx,
      integer(4):: nbmxe,nstatetot !nbmx,

c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

      integer(4)::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

      integer(4)::verbose,nstbz(nqbz),iqini,iqend !bzcase,
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
      complex(8) :: weightc(nrws1),zmeltt1

      complex(8),allocatable:: ppovl(:,:),ppovlz(:,:),zcousq(:,:),zmeltt(:,:,:)
      real(8),allocatable::vcoud(:),vcousq(:)
      real(8),parameter:: pi=4d0*atan(1d0),fpi=4d0*4d0*atan(1d0)

      integer:: ifvcoud,ivc,ngb0
      real(8)::qvv(3),vc
      logical :: newansisoW
      character(5):: charnum5

c      real(8),allocatable:: wklm(:),dmlx(:,:),epinvq0i(:,:)
c      integer:: ifidmlx,lxklm,nq0ix
      complex(8)::w3p

      integer:: mrecl,nprecx,nwordr,il

c--------------------------------------------------------------------
c RS: MPI parameters
      integer :: kx_local
c--------------------------------------------------------------------
      debug=.false.
      if(verbose()>=90) debug=.true.

!     ! === readin Vcoud and EPSwklm for newansisoW()=T ===
c$$$      if(newansisoW()) then
c$$$         ifidmlx = iopen('EPSwklm',0,0,0)
c$$$         read(ifidmlx) nq0ix,lxklm
c$$$         if(nq0i/=nq0ix) then
c$$$            write(6,*)'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
c$$$            call rx( 'nq0i from EPSwklm /= nq0i')
c$$$         endif
c$$$         allocate( dmlx(nq0i,9))
c$$$         allocate( epinvq0i(nq0i,nq0i) )
c$$$         allocate( wklm((lxklm+1)**2))
c$$$         read(ifidmlx) dmlx, epinvq0i
c$$$         read(ifidmlx) wklm
c$$$c         do il=1,(lxklm+1)**2
c$$$c         write(6,"('EPSwklm=', i5,3f12.5)") il,wklm(il)
c$$$c         enddo
c$$$         ifidmlx = iclose('EPSwklm')
c$$$      endif


c oct2005
      call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
      nbcutc=nctot+nbcut

      tpi         = 8d0*datan(1.d0)   
c      iq         = idxk (q,qbz,nqbz) ! index for q
c      write(6,"(' iq q  =',i4,3f8.4)")iq,q
cccc      iq          = idxk (q,qbze,nqbze) ! index for q
c      if(nctot/=0) ekc(1:nctot)= ecore(1:nctot) ! core
      nlmtobnd    = nlmto*nband
      nstatetot      = nctot + nband
      call minv33(qbas,qbasinv) 

c work arrays for psi2br.f
      if(debug) write(6,*) ' sxcf: 1'
      allocate(expikt(natom))

c      if(bzcase()==1) then
      if(abs(sum(qibz(:,1)**2))/=0d0) call RSMPI_Stop( ' sxcf assumes 1st qibz/=0 ')
      if(abs(sum( qbz(:,1)**2))/=0d0) call RSMPI_Stop( ' sxcf assumes 1st qbz /=0 ')
c      endif

      do it = 1,nwf
         itq(it) = it
      enddo

ctm debug
c        write(*,*)'isp,',isp,rw_w(1,1,1,1,1,0),rw_w(2,2,2,2,1,0)

c-----
      if(exchange.and.(.not.newansisoW())) then
         rewind  ifvcfpout
         read(ifvcfpout) ndummy1, ndummy2
      endif

c===============================
c loop over irreducible k-points
c===============================
ccccccccccccccccccccccccccccccccc
c      iii = ivsumxxx(irk,nqibz*ngrp)
c      write(6,*)' sxcf:sum non-zero irk=',iii
c      stop "sss"
c     
ccccccccccccccccccccccccccccccccc 

c      if(bzcase()==1) then
      kx = 1  ! qibz(:,1)=0 contribution for kcount
      if(irk(kx,irot)/=0) kount(kx)= kount(kx) + 1 
c      kount(kx)= kount(kx) + 1 
c      endif
      
c --- main loop start
      iqini=2
c      if(bzcase()==2) iqini=1
      iqend=nqibz+nq0i

      if(newansisoW()) then       !takao2012apr
         iqini=1
         iqend=nqibz            !no sum for offset-Gamma points.
c     if(exchange) iqend=nqibz
      endif
ccccccccccccccccccccccccccc
c      iqini=2
ccccccccccccccccccc

c RS: ifile_rsmpi is defined in gwsrc/RSMPI_mod.F
      write(ifile_rsmpi,*) "RS: irot = ", irot
      write(ifile_rsmpi,*) "RS: main loop start"

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
c        write(6,*)'iqend=',iqend
c      endif  
cccccccccccccccccccccccccccc
      do 1100 kx_local = 1,nk_local_rotk(irot)
        kx = ik_index_rotk(irot,kx_local)
c      do 1100 kx = iqini,iqend !kx=1 corresponds to q=0 is omitted.
c debug:
c      do 1100 kx = iqini,iqini !kx=1 corresponds to q=0 is omitted.
c        write (6,"(i3,'  out of ',i3,$)") kx,iqend
        write(ifile_rsmpi,*) ' wmat_MPI: goto loop kx=',kx
        if(debug)  write(6,*) ' sxcf: goto loop kx=',kx
        
!        write(*,'("1  begin k-cycle",$)') 
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
c        write(6,*) ' k ngc=',k,ngc
c        ngb = nbloch + ngcni(k)

        call readqg0('QGcou',qibz_k,  quu,ngc)
c        ngc = ngcni(k)  ! k-points in IBZ
        ngb = nbloch + ngc

!! ===Readin diagonalized Coulomb interaction===
!! note sep102012takao
!!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
!!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
!!  Vould file is written in hvccfp0.m.F.
!! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
!! These can be also separeted into WVR.ID and WVI.ID files.
        if(newansisoW()) then
           qxx=qibz_k
c           if(kx<=nqibz) qxx=qibz_k
c           if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
           ifvcoud = iopen('Vcoud.'//charnum5(kx),0,0,0)
           do
              read(ifvcoud) ngb0
              read(ifvcoud) qvv
c              write(6,"('readin qvv ngb0=',3f9.4,i5)")qvv,ngb0
c              write(6,"('readin qxx ngb0=',3f9.4,i5)")qxx
              if(allocated(vcoud)) deallocate(vcoud)
              allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
              read(ifvcoud) vcoud
c$$$  cccccccccccccccccccccccccccccccccccccccccccccc
c$$$  if(sum(abs(qxx(1:2)))<1d-4) then
c$$$  do irot2 = 1,ngrp
c$$$  kr = irkip(kx,irot2,ip) ! index for rotated kr in the FBZ
c$$$  if(kr==0) cycle ! next irot
c$$$  vcoud(1) =vcoud(1)*wqfac(kr)
c$$$  exit
c$$$  enddo
c$$$  write(6,*)'wwwwwwwwwww why here ?wwwwwwwwww'
c$$$  endif
c$$$  cccccccccccccccccccccccccccccccccccccccccccccc
              read(ifvcoud) zcousq
              if(sum(abs(qvv-qxx))<1d-6) goto 1133
           enddo
           if(sum(abs(qvv-qxx))>1d-6) then
              write(6,*)'qvv =',qvv
              write(6,*)'qxx=',qxx,kx
              call rx( 'wmatK: qvv/=qibz(:,kx) hvcc is not consistent')
           endif
 1133      continue
           ifvcoud = iclose('Vcoud.'//charnum5(kx))
           if( ngb0/=ngb ) then !sanity check
              write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
              call rx( 'hsfp0.m.f:ngb0/=ngb')
           endif
c$$$  if(sum(abs(qibz_k))<1d-6) then
c$$$  idummy  = iclose('Vcoud') !close and open again. This is because first data in Voud is for q=0
c$$$  ifvcoud = iopen('Vcoud',0,0,0)
c$$$  endif
c$$$  read(ifvcoud) ngb0
c$$$  if( ngb0/=ngb ) then
c$$$  write(6,*)' qibz_k=',qibz_k,ngb0,ngb
c$$$  write(6,*)' qibz_k=',qibz_k
c$$$  stop 'hsfp0.m.f:ngb0/=ngb'
c$$$  endif
c$$$  read(ifvcoud) qvv
c$$$  if(sum(abs(qvv-qibz_k))>1d-6) then
c$$$  write(6,*)'qvv =',qvv
c$$$  write(6,*)'qibz_k=',qibz_k,kx
c$$$  stop 'sxcf_fal2: qvv/=qibz(:,kx) hvcc is not consistent'
c$$$  endif
c$$$  if(allocated(zcousq)) deallocate(zcousq,vcousq,vcoud)
c$$$  allocate( zcousq(ngb0,ngb0),vcousq(ngb0),vcoud(ngb0) )
c$$$  read(ifvcoud) vcoud
c$$$  read(ifvcoud) zcousq
c$$$  vcousq=sqrt(vcoud) !

!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi_q,itp |phi_q-rk,it B_rk,I> O^-1_IJ
!! ppovlz= O Zcousq
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
           if(allocated(ppovlz)) deallocate(ppovlz)
           allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
           call readppovl0(qibz_k,ngc,ppovl)
           ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
           ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
           deallocate(zcousq,ppovl)
        endif

        if(.not.newansisoW()) then 
           if(exchange) then
c RS: set pointer to the right place
          call set_vcoul_rsmpi(ifvcfpout,kx-iqini+1)
c     allocate(vcoul(ngb,ngb))
c     read(ifvcfpout) vcoul(1:ngb,1:ngb)
              read(ifvcfpout) nn !oct2005
              allocate(vcoul(nn,nn))
              read(ifvcfpout) vcoul(1:nn,1:nn)
           endif
!     ! weight check for cycle or not.
           if(kx <= nqibz ) then
              if (kr == 0)    then
c     stop 'wmat: kr=0'
                 if(exchange) deallocate(vcoul)
                 cycle 
              endif
              kount(kx)= kount(kx) + 1 ! count the no. times k
! appears in the 1st BZ
cccccccccccccccccccccccccccccccccccccccccccccccccc
c     write(6,*)' irot,ip, k, kount in  =',irot, ip, k, kount(k,ip)
c     deallocate(vcoul)
c     cycle
c     write(6,*)' kount out =',kount(k)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     if (kount(kx) > nstar(kx)) stop 'sexc: too many stars'
              if (kount(kx) > nstar(kx)) stop 'wmat: kount > nkstar'
c     if (kount(kx) > 1) stop 'wmat: kount > 1'
           else
              if( wgt0(kx-nqibz,irot) == 0d0 ) then 
                 if(exchange) deallocate(vcoul)
                 cycle
              endif
           endif
        else
           if (kr == 0) cycle 
        endif
c---test
      if(OnlyQ0P.and.kx<=nqibz) then
        if(exchange) deallocate(vcoul)
        cycle
      endif
!! phase factor for off-site W
      do ir1=1,nrws1
         expiqR1(ir1) = exp(-img*tpi* sum(qbz_kr(:)*rws1(:,ir1)))
      enddo

      
!! ===================================================================
      allocate( rmelt3(ngb,nwf,nwf,nrws2),cmelt3(ngb,nwf,nwf,nrws2))
      rmelt3 = 0d0
      cmelt3 = 0d0
!! loop over FBZ
      do iq = 1,nqbz 
        q(:) = qbz(:,iq)
        call readcphiW (qbz(:,iq), nlmto,isp, quu, cphiq)
        qk =  q - qbz_kr          ! qbz(:,kr)
        call  readcphiW(qk, nlmto,isp, quu, cphikq)
        do ia = 1,natom
         expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
        end do
        if(debug) write(6,*) ' sxcf: tr=',tr
        if(debug) write(6,*) ' sxcf: goto psicb2'
        nbmax = nwf
        nt   = nctot + nbmax      ! = nstate for the case of correlation
        ntp0 = nwf
        allocate( zzmel(nbloch,nt,ntp0)) ! rk,ibloch  q-rk,it  q,itp
        zzmel = 0d0
        call psi2b_v2 (nbmax, ntp0, iclass,
     i     dreal(expikt(1:natom)),dimag(expikt(1:natom)),
     i     cphikq,              !rbkq,cbkq,rhbkq,chbkq, !  q-rk nko
     i     cphiq,               !rbq,cbq,rhbq,chbq,     !  q    nko
     i     ppb,                 !pdb,dpb,ddb,
     i     nlnmv,nlnmc,mdim,nctot,
     i     imdim,iatomp,
     d     mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0,
     d     natom,nclass,
     o     zzmel)               ! rmel,cmel)
        if(debug) write(6,"('sum of zmel abszmel=',4d23.16)")
     &                 sum(zzmel),sum(abs(zzmel) )
        allocate( zmelt(ngb, nctot+nbmax, ntp0) ) ! rk,ibloch  q-rk,it  q,itp
        if(debug) write(6,*) ' sxcf_fal2: goto drvmelp xxxxx1',ngb,nctot,nbmax,ntp0
        call drvmelp3( q,   ntp0, ! q in FBZ  q,itp
     i                 q-qbz_kr, nbmax,  ! q-rk,it
     i                 qibz_k,           ! k in IBZ for e-product basis
     i                 isp,ginv,   
     i       ngc,ngcmx,ngpmx,nband,itq, 
     i       symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz,
     i       dreal(zzmel), dimag(zzmel), nbloch, nt,nctot,
     o       zmelt) 
        deallocate(zzmel) !rmel,cmel)
        if(debug) write(6,*) ' sxcf: goto wtt'
        if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum(zmelt)
        do ir2=1,nrws2
         expiqR2 = exp( img*tpi* sum(q(:)*rws2(:,ir2)))
         rmelt3(:,:,:,ir2) = rmelt3(:,:,:,ir2) + wk(iq) * dreal(zmelt(:,:,:)*expiqR2)
         cmelt3(:,:,:,ir2) = cmelt3(:,:,:,ir2) + wk(iq) * dimag(zmelt(:,:,:)*expiqR2)
        enddo  ! ir2
        deallocate(zmelt)
      enddo
!! ===================================================================


      
      if(kx<= nqibz) then 
         wtt = wk(kr)           !         wtx = 1d0
      else  
         wtt = wk(1)*wgt0(kx-nqibz,irot) ! wtx = wgt0(kx-nqibz,irot)
         if(abs(wk(1)-1d0/dble(nqbz))>1d-10) stop 'sxcf:wk(1) inconsistent'
      endif
      weight = wtt
      if(debug) then
         write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
      endif
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
         allocate(zmel1(ngb))
         allocate(zmel(ngb, nwf, nwf))
         if(.not.newansisoW()) allocate(vcoult(1:ngb,1:ngb),z1p(ngb,nwf,nwf))
         do ir2=1,nrws2
!!  (rmelt3,cmelt3)  !rk,ibloch  q-rk,it  q,itp
          zmel  = dcmplx (rmelt3(:,:,:,ir2),cmelt3(:,:,:,ir2)) !<psi_itp|psi_it B>
          if(.not.newansisoW() ) then
            vcoult= transpose(vcoul)
            call matm( vcoult, zmel, z1p,  ngb,ngb,nwf*nwf ) !z1p= <psi_itp|psi_it B> * voul
c            deallocate(vcoult)
            do ir3=1,nrws2
            do itp2 = 1,nwf
            do it2  = 1,nwf
            do it   = 1,nwf
            do itp  = 1,nwf
            zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),-cmelt3(:,it,itp,ir3)) ! <B psi_it|psi_itp>
            ztmp = sum ( z1p(:,it2,itp2)*zmel1 ) ! <psi_itp2|psi_it2 B>*vcoul*<B psi_it|psi_itp>
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
            enddo ! ir3 
          else
! based on E_I basis. See Christoph's paper
            allocate(zmeltt(nwf,nwf,ngb))
            do itp= 1,nwf
            do it = 1,nwf
            do ivc=1,ngb
            zmeltt(it,itp,ivc)=sum(zmel(:,it,itp)*ppovlz(:,ivc)) ! <psi_itp|psi_it E_I> (I=ivc) 
            enddo
            enddo
            enddo
            do ir3=1,nrws2
            do itp2 = 1,nwf
            do it2  = 1,nwf
            do it   = 1,nwf
            do itp  = 1,nwf
              zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),cmelt3(:,it,itp,ir3)) ! <psi_itp|psi_it B_I>
              w3p=0d0
              do ivc=1,ngb
                zmeltt1 =  sum( zmel1(:)*ppovlz(:,ivc) ) !<psi_itp|psi_it E_I>
                if(ivc==1.and.kx==iqini) then
                   vc= wklm(1)* fpi*sqrt(fpi) /wk(kx) !kx right?
                else
                   vc= vcoud(ivc)
                endif
                w3p=w3p + zmeltt(it2,itp2,ivc)  *vc* dconjg(zmeltt1) 
                       ! <psi_itp2|psi_it2 E_I> *vc* <E_I psi_it|psi_itp>
              enddo
              ztmp= w3p
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
            enddo ! ir3 
            deallocate(zmeltt)
          endif
         enddo ! ir2
         if(allocated(vcoul)) deallocate(vcoul,vcoult)
         if(allocated(z1p))   deallocate(z1p)
         if(allocated(rmelt3)) deallocate(rmelt3,cmelt3,zmel1, zmel)
         if (debug) write(*,*) 'bare coulomb section finished'
!! --- End of bare-Coulomb section --------------
      else
c--------------------------------------------------------------------------
c--- screening effect section----------------------------------------------
c--------------------------------------------------------------------------
c S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> 
c                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n') |psi(q,t')>
c t-->itp   n -->it  
c t'-->itp2 n'-->it2
c--------------------------------------------------------------
!!--- The matrix elements zmelc.
!! zmelc  = < E(rk,j) psi(q-rk,it) | psi(q,itp) >
!! E basis is ghe Christoph's basis diagonalize the Coulomb inteaction
        allocate( zmelc (ngb, nwf, nwf,nrws2),zmelcc (ngb, nwf, nwf,nrws2),
     c            zw (nblochpmx,nblochpmx),
     c            zw2(nwf,nwf,nwf,nwf) )
        zmelcc = dcmplx (rmelt3,-cmelt3) !zmelcc = < B(rk,j) psi(q-rk,it) | psi(q,itp) >
        do it =1,nwf
        do itp=1,nwf
        do ir2=1,nrws2
           zmelc(:,it,itp,ir2) =  matmul(zmelcc(:,it,itp,ir2),dconjg(ppovlz(:,:)))
        enddo
        enddo
        enddo
        deallocate(rmelt3,cmelt3,zmelcc)
        if(debug) write(6,*)' end of zmel'
c====================================================================
c Wc(qt,w) along the imaginary axis
c====================================================================
c------------------------------------------------
c loop over w' = (1-x)/x, frequencies in Wc(k,w')
c {x} are gaussian points between (0,1)
c------------------------------------------------
        nx  = niw
        nprecx=8
        mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
        ifrcwi = iopen('WVI.'//charnum5(kx),0,-1,mrecl)
        do ix = 1,nx     ! imaginary frequency w'-loop
ccccccccccccccccccccccccccccccccc
c          nrec=(kx-2)*niw+ix
c          if(bzcase()==2) nrec= (kx-1)*niw+ix
          nrec=ix
          read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis 
c          read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis 
ccccccccccccccccccccccccccccccccc

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
          call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), zmelc(:,:,:,ir3), 
     i      nwf,nwf,ngb, 
     o      zw2)
          do ir1=1,nrws1
            ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
          rw_iw(:,:,:,:,ir,ix) = rw_iw(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
          cw_iw(:,:,:,:,ir,ix) = cw_iw(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
          enddo ! ir1
          enddo ! ir2
          enddo ! ir3
        enddo
        ifrcwi = iclose('WVI.'//charnum5(kx))

c====================================================================
c Wc(qt,w) along the real axis
c====================================================================
        if(debug) write(6,*)' go to poles'
        nprecx=8
        mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
        ifrcw = iopen('WVR.'//charnum5(kx),0,-1,mrecl)
        do      ix = 0,nrw                    ! real frequency w'-loop
ccccccccccccccccccccccccccc
c          nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
c          if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
          nrec=ix-nw_i+1
          read(ifrcw,rec=nrec) zw 
cccccccccccccccccccccccccc
c zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)> 
          do ir3=1,nrws2
          do ir2=1,nrws2
          call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), zmelc(:,:,:,ir3), 
     i      nwf,nwf,ngb, 
     o      zw2)
          do ir1=1,nrws1
            ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
          rw_w(:,:,:,:,ir,ix)  = rw_w(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
          cw_w(:,:,:,:,ir,ix)  = cw_w(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
          enddo ! ir1
          enddo ! ir2
          enddo ! ir3
        enddo
        ifrcw = iclose('WVR.'//charnum5(kx))
        deallocate(zmelc,zw,zw2)
        if(debug) write(6,*)' end of screen-if'
c end of if (exchange)
      endif
 1100 continue  ! end of k-loop
      return
c     end subroutine wmatqk_mpi

!!-------------------------------------------------------------------------------      
      contains
     
!! psi2b_v2 and psicb_v2 are older versions of psi2b_v3 and psicb_v3 in m_zmel
      subroutine psi2b_v2(nt0,ntp0,iclass,coskt,sinkt,
     i                   cphik, !rbk,cbk,rhbk,chbk,
     i                   cphikq,!rbkq,cbkq,rhbkq,chbkq,
     i                   ppb, !pdb,dpb,ddb,
     i                   nlnmv,nlnmc,mdim,nctot,imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,noccxv,nt,ntp,
     d                   natom,nclass,
     o zpsi2b)
c originaly 92.03.17 by Ferdi.
c takao modified at Apr 2002
c calculates <psi(k',t') | psi(k,t) B(R,i)>
c for all R
c psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
c B(R,i)   = Bloch orthonormal product basis for atom R
c psi(k,t) is stored after nctot
c
c nt0        = no. t
c ntp0       = no. t'
c coskt,sinkt= exp(ik.T)
c cphik b(k)
c cphikq b(k')
c
c ppb        = <phi(RLn) phi(RL'n') B(R,i)>
c
c ddb        = <phidot(RLn) phidot(RL'n') B(R,i)>, s. ppbl.f
c nlnmv      = number of l,n,m for valence
c nlnmc      = number of n,l,m for core states
c mdim       = number of optimal product basis functions for each class
c nctot      = total no. allowed core states
c nbloch     = total no. optimal product basis
c nlnmx      = maximum number of l,n,m
c noccxv     = max. no. occupied valence states
c nt         = maximum number of occupied states
c ntp        = maximum number of unoccupied states
c
c zpsi2b     =  the matrix elements
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
      integer::nt0,ntp0,natom,nlmto,nclass,nlnmc,mdim,nctot
      complex(8):: cphik(nlmto,noccxv),cphikq(nlmto,ntp0)
     &            ,zpsi2b(nbloch,nt,ntp),phase
      dimension 
     i          ppb(nlnmx,nlnmx,mdimx,nclass),
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i          coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)

      integer(4),allocatable::iasx(:)
      integer(4) :: nzppb1,nzppb2
      complex(8),allocatable :: zz(:,:), zppb(:,:)
      complex(8) :: alpha,beta

c zppb is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
      nzppb1=max(nt0,nlnmx)
      nzppb2=max(ntp0,nlnmx)
      allocate( zz(nlnmx,ntp) )
      allocate( zppb(nzppb1,nzppb2) )
      beta=0d0  ; alpha=1d0

c check dimensions
      if (ntp0 .gt. ntp) call rx( 'psi2bc: ntp exceeded')
      if (mdimx .ne. maxval(mdim)) call rx( 'psi2bc: wrong mdimx')
      if (nctot+nt0 .gt. nt) call rx( 'psi2bc: nt exceeded')
      if (nt0 .gt. noccxv) call rx( 'psi2bc: noccxv exceeded')
      if ( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v2: wrong nbloch')
      allocate(iasx(natom))
      ias = 1
      do ia = 1,natom
        iasx(ia) = ias
        ias = ias + nlnmv(iclass(ia))
      enddo
      if(ias-1/=nlmto) call rx( ' psi2b_v2:sum(nlnmv)/= nlmto')
c loop over atoms
      do  ia = 1,natom
        ic   = iclass(ia)
        nc   = nlnmc(ic)
        nv   = nlnmv(ic)
        nc1  = nc + 1
        if (nc+ nlnmv(ic) .gt. nlnmx) call rx( 'psi2b_v2: nlnmx exceeded')
        phase= dcmplx(coskt(ia),sinkt(ia))
        ias  = iasx(ia)
        iap  = iatomp(ia)
        icp  = iclass(iap)
        do   i = 1,mdim(icp) ! loop over optimal product basis
c---------------------------------------------------
cc sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)>
cc for a given i, for all L'n' and t'
cc bkq is complex but < > is real
c1      do     itp = 1,ntp0
c1      do      jp = 1,nlnmv(ic)
c1      zz (jp,itp)=dconjg(
c1     &    sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,nc+jp,i,icp)) )
c1      end do
c1      end do
c
c2      zz(1:nv,1:ntp0) =dconjg(
c2     & matmul(  transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp))
c2     &         ,cphikq(ias:ias+nv-1,1:ntp0)) )
c
c3        call dgemm('T','N',nv,ntp0,nv,
c3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
c3     &          dreal(cphikq(ias:ias+nv-1,1:ntp0)), nv,
c3     &          0d0,
c3     &   rr, nlnmx )
c3        call dgemm('T','N',nv,ntp0,nv,
c3     &   1d0, ppb(nc1:nc+nv,nc1:nc+nv,i,icp),     nv,
c3     &          dimag(cphikq(ias:ias+nv-1,1:ntp0)), nv,
c3     &          0d0,
c3     &   cc, nlnmx )
          zppb(1:nv,1:nv) = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp)
          call zgemm('T','N',nv,ntp0,nv,
     &   alpha, zppb,nzppb1, cphikq(ias,1), nlmto,  beta,
     &   zz,  nlnmx )
          do itp = 1,ntp0
            do jp = 1,nv
              zz(jp,itp)= dconjg(zz(jp,itp) )
            enddo
          enddo
c----------------------------------------------------
c <psi(k+q,t') | psi(k,t) B(i)>
c1      do      it = 1,nt0
c1      do     itp = 1,ntp0
c1       zpsi2b(ib,nctot+it,itp)=
c1     &   phase * sum( zz(1:nv,itp)*cphik(ias:ias+nv-1,it) )
c1c end of t'(unoccupied)-loop
c1      end do
c1c end of t(occupied)-loop
c1      end do
c3        call zgemm('T','N',nt0,ntp0,nv,
c3     &   phase, cphik(ias:ias+nv-1,1:nt0),  nv,
c3     &          dcmplx(rr(1:nv,1:ntp0),-cc(1:nv,1:ntp0)),  nv,
c3     &          0d0,
c3    &   zpsi2b(imdim(iap)-1+i,nctot+1:nctot+nt0,1:ntp0), nt0)
          call zgemm('T','N', nt0,ntp0,nv,
     &  phase, cphik(ias,1),nlmto, zz,nlnmx, beta,
     &  zppb, nzppb1 )
          ib = imdim(iap)-1+i
          zpsi2b(ib,nctot+1:nctot+nt0,1:ntp0)=zppb(1:nt0,1:ntp0)
c------------------------------------------------------
        end do !end of optimal product basis-loop
      end do !end of atom-loop
c      deallocate(rr,cc,iasx)
      deallocate(zz,zppb,iasx)
      end

c------------------------------------------------------------------------------------
      subroutine psicb_v2 (icore,ncore,ntp0,iclass,coskt,sinkt,
     i                   cphikq, !rbkq,cbkq,rhbkq,chbkq,
     i                   ppb,!pdb,dpb,ddb
     i                   nlnmv,nlnmc,mdim,
     i                   imdim,iatomp,
     d                   mdimx,nlmto,nbloch,nlnmx,nt,ntp,natom,nclass,
     d                   nl,nnc,
     o zpsi2b)! rpsi2b,cpsi2b)
c written by Ferdi  92.03.17
c takao modified at Apr 2002
c
c calculates <psi(k+q,t') | core(k,t) B(R,i)>
c for all R
c psi(k,t) = S[RLn] b(RLn,k,t)*X(RLn,k)
c          = S[RLn] b(RLn,k,t)*Phi(RLn,k) + hb(RLn,k,t)*Phidot(RLn,k)
c core(k,t)= core states
c B(R,i)   = Bloch orthonormal product basis for atom R
c
c <psi(k+q,t') | core(k,t) B(R,i)>
c = S[RLn]  b(RLn,k+q,t')^* <Phi(RLn)    |core(k,t) B(R,i)>
c + S[RLn] hb(RLn,k+q,t')^* <Phidot(RLn) |core(k,t) B(R,i)>

c ntp0       = no. unoccupied states
c coskt,sinkt= exp(ik.T)
c cphikq  = real and imaginary part of b(k+q).
c            coefficients of eigenfunctions for argumentationwaves in each MT
c
c icore      = index for core states
c ncore      = no. core states in each class
c ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>
c
c nlnmv      = number of l,n,m for valence
c nlnmc      = number of l,n,m for core states
c mdim       = number of optimal product basis functions for each class
c nbloch     = total no. optimal product basis
c nlnmx      = maximum number of l,n,m
c nt         = maximum number of occupied states
c ntp        = maximum number of unoccupied states
c
c zpsi2b     =  the matrix elements
c
      implicit real*8(a-h,o-z)
      implicit integer (i-n)
      complex(8):: cphikq(nlmto,ntp0),zpsi2b(nbloch,nt,ntp),phase
      dimension
c                rbkq(nlmto,ntp0),cbkq(nlmto,ntp0),
c     i          rhbkq(nlmto,ntp0),chbkq(nlmto,ntp0),
     i          icore(nl*nl*nnc,nclass),ncore(nclass),
     i          ppb(nlnmx,nlnmx,mdimx,nclass),
c     i          pdb(nlnmx,nlnmx,mdimx,nclass),
c     i          dpb(nlnmx,nlnmx,mdimx,nclass),
c     i          ddb(nlnmx,nlnmx,mdimx,nclass),
     i          nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),
     i          coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)
c      dimension rpsi2b(nbloch,nt,ntp),
c     o          cpsi2b(nbloch,nt,ntp)
c initialise matrix elements
c      call dinit   (rpsi2b,nbloch*nt*ntp)
c      call dinit   (cpsi2b,nbloch*nt*ntp)

      zpsi2b = 0d0
c loop over atoms
      ib         = 0
      ias        = 1
      ics        = 0
      do      ia = 1,natom
        ic         = iclass(ia)
        nc         = nlnmc(ic)
        nv         = nlnmv(ic)
        nc1        = nc + 1
        phase  =  dcmplx(coskt(ia),sinkt(ia))
c loop over optimal product basis
        iap        = iatomp(ia)
        icp        = iclass(iap)
        ib         = imdim(iap)-1
        do       i = 1,mdim(icp)
          ib         = ib + 1

c S[Ln] bkq(Ln,t')^(*) * <Phi(Ln) core(L'n') B(i)> etc.
c for a given i, for all L'n' and t'
c bkq is complex but < > is real
          do     itp = 1,ntp0
            do      it = 1,ncore(ic)
              icr        = icore(it,ic)

cccccccccccccc
c      write(6,*),it,ic,icore(it,ic)
cccccccccccccc
c
c      rs1        = vdv ( rbkq(ias,itp),ppb(nc1,icr,i,icp),nv)
c     .           + vdv (rhbkq(ias,itp),dpb(nc1,icr,i,icp),nv)
c      cs1        = vdv ( cbkq(ias,itp),ppb(nc1,icr,i,icp),nv)
c     .           + vdv (chbkq(ias,itp),dpb(nc1,icr,i,icp),nv)
c
c      rpsi2b(ib,ics+it,itp) = rs1*coskt(ia) + cs1*sinkt(ia)
c      cpsi2b(ib,ics+it,itp) = rs1*sinkt(ia) - cs1*coskt(ia)
cccccccccccccccccccccccccccccccc
c      if(abs(rs1)>1d8.or.abs(cs1)>1d8) then
c        print *,'  psicb2*:'
c        print *, nc1,icr,i,icp, ppb(nc1,icr,i,icp)
c        print *,    dpb(nc1,icr,i,icp)
c       stop
c      endif
cccccccccccccccccccccccccccccccc

              zpsi2b(ib,ics+it,itp) = phase* 
     & dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))

c end of t'(unoccupied)-loop
            end do
c end of t(occupied)-loop
          end do

c end of optimal product basis-loop
        end do

c end of atom-loop
        ias        = ias + nlnmv(ic)
        ics        = ics + ncore(ic)
      end do

      return
      end


c$$$c -------------------------------------------------------------------
c$$$      subroutine matzwz2(zw,zmel, ntp0,nstate,ngb, zwz)
c$$$      implicit none
c$$$      integer(4) :: nstate,ntp0,itp,it,itp2,it2,ngb
c$$$      complex(8) :: zw(ngb,ngb),zmel(ngb,nstate,ntp0),
c$$$     c              zwz(ntp0,nstate,nstate,ntp0)
c$$$      complex(8), allocatable :: CC(:,:,:)
c$$$      allocate(CC(ngb,nstate,ntp0) )
c$$$      call matm(zw,zmel,cc, ngb, ngb, nstate*ntp0)
c$$$      do itp2 = 1,ntp0
c$$$      do  it2 = 1,nstate
c$$$      do  it  = 1,nstate
c$$$      do itp  = 1,ntp0
c$$$        zwz(itp,it,it2,itp2) 
c$$$     & = sum( dconjg(zmel(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      enddo
c$$$      deallocate(CC)
c$$$      end
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

      end subroutine wmatqk_mpi
      
      logical function newansisoW()
      newansisoW=.true.
      end



