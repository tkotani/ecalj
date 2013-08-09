      subroutine x0kf_v2z (
     i              wwk,n1b,n2b,nbnbx, nbnb, tetra,  ! for tetra
c     i                  q,ifev1,ifev2,ifcphi,nsp,isp,    !ifrb,ifcb,ifrhb,ifchb,
     i                  q,isp1,isp2,nsp, !isp,   ! ifcphi,ifrb,ifcb,ifrhb,ifchb,
     i                  qbas,ginv,ef,rk,wk,      !iindxk,
c     i                index_qbz,n_index_qbz,
cc     i                mdim,inlnm,inlnmv,inlnmc,iiclass,
     i                  mdim, nlnm, nlnmv, nlnmc, iclass, ! not the w(pointer)
     i                  ppb, !pdb,dpb,ddb,
     i                  freqr,freqi,delta,realomega,imagomega,
cc     i                  iicore,incore,
     i                  icore,ncore, ! not the w(pointer)
     &                  ecore,
c
     d                  nlmto,nqbz,noccx,noccxv,
     d                  nl,nclass,natom,nnc,
     d                  nlnmx,mdimx,nbloch,nwt,niwt,
c
c     i      iq,  ngpn,ngc, ngpmx_in,ngcmx, !geigB,
c     i       ngvecpB,ngveccB,
     i      iq, ngc, ngpmx_in,ngcmx,
c     i      ngveccB,
     i      nqbze, nband,nqibz,
c
     o      zxq,zxqi,
     i      nolfco,zzr,
     o      x0mean, !)
     i   vcoul,ngbx)
      use m_readqg, only : readqg
      use m_readeigen, only : readeval,readcphi,readgeig
      use keyvalue
c takao kotani Aug 2002   This originated from Ferdi's x0k.
c
Ci realomega: =T -> calculate x0 along real axis
Ci imagomega: =T -> calculate x0 along imag axis
Ci
c zadd and zaddr dominates the cpu time.
c
c calculates x0(i,j)(q,w) = <M(i)|x0(q,w)|M(j)> for a given spin
c direct access read
c
c x0(i,j)(q,w) = S[k=FBZ] S[t=occ] S[t'=unocc]
c  <M(q,i) psi(k,t) |psi(k+q,t')> <psi(k+q,t')| psi(k,t) M(q,j)>
c  { 1/[w-e(k+q,t')+e(k,t)+i*delta] - 1/[w+e(k+q,t')-e(k,t)-i*delta] }
c  ; w is real or imaginary
c
c zzmel =  <psi(k+q,t') | psi(k,t) B(R,i)>
c zmel  =  <psi(k+q,t') | psi(k,t) M(R,i)>
c zxq zxqi=  zeroth order response function along
c            the positive real or imaginary axis
c
c q       = q-vector in x(q,iw)
c ifchi   = direct access unit file for cphi, the coefficient of eigenfunction for argumentation wave.
c qbas    = base reciprocal lattice vectors
c ginv    = inverse of qbas s. indxrk.f
c ef      = fermi level in Rydberg
c
c ppb     = <phi(RLn) phi(RL'n') B(R,i)>
c
c iclass  = given an atom, tells the class
c iindxk  = index for k-points in the FBZ
c rk      = k-points in the 1st BZ
c wk      = weight for each k-point in the 1st BZ
c freq    = frequency points along positive imaginary axis
c delta   = broadening for delta function
c
c mdim    = dimension of B(R,i) for each atom R
c nlnmx   = maximum number of l,n,m
c nlmto   = total number of LMTO basis functions
c nqbz    = number of k-points in the 1st BZ
c n1,n2,n3= divisions along base reciprocal lattice vectors
c natom   = number of atoms
c noccx   = maximum number of occupied states
c noccxv  = maximum number of occupied valence states
c nbloch  = total number of Bloch basis functions
c niw     = number of frequency points along
c           the positive imaginary axis
c
c cphik cphikq:  b(k) and b(k+q)
c  : coefficients of eigenfunctions for argumentation waves in each MT
c
c ek      = eigenvalues at k-points in the 1st BZ
c ekq     = eigenvalues at k+q, k in the 1st BZ
c
c
c
      implicit real*8(a-h,o-z)
      dimension q(3),qbas(3,3),ginv(3,3),rk(3,nqbz),wk(nqbz),
     i          freqr(nwt),freqi(niwt),
c     w          rmel(nbloch,noccx,nband),
c     w          cmel(nbloch,noccx,nband),
c     w          rbk (nlmto,noccxv),rbkq (nlmto,nband),
c     w          cbk (nlmto,noccxv),cbkq (nlmto,nband),
c     w          rhbk(nlmto,noccxv),rhbkq(nlmto,nband),
c     w          chbk(nlmto,noccxv),chbkq(nlmto,nband),
     w          ek(nband),ekq(nband),ekc(noccx), ecore(noccx - noccxv)
c automatic
c     w          w1(nbloch+ngcmx,nbloch+ngcmx)
c     &         ,w2(nbloch+ngcmx,nbloch+ngcmx)
c
c      real(8):: rxq (nblochpmx,nblochpmx, nwt),
c     o          cxq (nblochpmx,nblochpmx, nwt),
c     o          rxqi(nblochpmx,nblochpmx,niwt),
c     o          cxqi(nblochpmx,nblochpmx,niwt)
      complex(8):: zxq (nbloch+ngc,nbloch+ngc, nwt),
     o             zxqi(nbloch+ngc,nbloch+ngc ,niwt)

c      dimension qk(3),qbkp(3)
      complex(8) :: imag=(0d0,1d0),trc

      complex(8),allocatable:: cphik(:,:),cphikq(:,:)

c      integer w
c      common /w/ w(1)

ctakao -----------------------------------------------------------
      integer(4):: ngpmx_in, ngcmx, nqbze, nband,
     &             ngc,nadd(3), !ngpn(nqbze),ngvecpB(3,ngpmx_in,nqbze),
     &             ngvecpB1(3,ngpmx_in),ngvecpB2(3,ngpmx_in) !, ngveccB(3,ngcmx),igc
c      complex(8):: geigB  (ngpmx_in,nband,nqbze)
      complex(8):: geig1(ngpmx_in,nband),geig2(ngpmx_in,nband)
c      complex(8),allocatable:: zmelpl(:,:,:)
c      real(8),allocatable::  rmelt(:,:,:),cmelt(:,:,:)
      complex(8),allocatable :: zmelt(:,:,:)
c
      real(8) :: qbasinv(3,3), det,qdiff(3),add(3),symope(3,3),shtv(3)
c     &  symope(3,3)=(/1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/),
c     &  shtv(3)=(/0d0,0d0,0d0/)
      data symope /1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0/
      data shtv   /0d0,0d0,0d0/


      real(8) :: ppb(*) !,pdb(*),dpb(*),ddb(*)
c      real(8)::     ppb(nlnmx,nlnmx,mdimx,nclass),

      real(8) :: coskt(natom), sinkt(natom),qx(3)
      integer(4) :: mdim(natom)
      complex(8),allocatable :: ttx(:,:),z1p(:,:)
      complex(8) :: zviw
      logical :: realomega, imagomega

c tetra ---
      logical :: tetra
      complex(8) ::  wwk(nbnbx,nqbz,niwt+nwt)
      integer(4) ::  nbnb(nqbz), 
     &  n1b(nbnbx,nqbz), n2b(nbnbx,nqbz)

#ifdef USE_X0KBLAS
      complex(8),allocatable,save:: wwkr(:,:),wwki(:,:), z1pp(:,:,:)
      complex(8):: zalpha,zbeta
#endif

      complex(8),allocatable:: zzmel(:,:,:)
      integer(4)::imdim(natom),iatomp(natom)

c for ixc==2 or ixc==22 mode
c      integer(4):: ixc
      logical :: nolfco
      complex(8):: x0mean(nwt),zq0,zzr(nbloch+ngc)
      real(8)   :: zq0zq0


c      integer(4):: n_index_qbz
c      integer(4):: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

c-----
      integer(4):: nlnm(*),nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*)

      integer(4):: isp1,isp2,verbose
      real(8):: qu1(3),qu2(3),qu1x(3),qu2x(3),quu(3)


      integer(4),allocatable:: m_in(:),n_in(:),l_in(:),ibas_in(:)


      real(8):: eband(nband)
      logical :: iww1=.true.,iww2=.true.

      real(8):: eee(2),ecut,ecuts,wcut
      logical ::ipr=.false.

ccccc zvz test cccccccccccccccccccccccccc
      integer(4):: zvztest
      integer(4):: ngbx
      complex(8):: vcoul(ngbx,ngbx)
      complex(8),allocatable:: vzz(:,:,:),ppovl(:,:),ppovlinv(:,:),aaa(:)
      complex(8):: zvz,zvz1
      integer(4):: ib1,ib2,ifix,iqx
      integer(4):: ippovlinv=501,ngcx
      logical :: ippovlinv_test=.true.
      if(verbose()>35) ipr=.true.
      if(zvztest()==1) then
Cstop2rx 2013.08.09 kino        if(ngbx /= ngc+nbloch) stop 'x0kf_v2z:ngbx /= ngc+nbloch'
        if(ngbx /= ngc+nbloch) call rx( 'x0kf_v2z:ngbx /= ngc+nbloch')
        if(ippovlinv_test) then
          open (ippovlinv,form='unformatted',file='PPOVLinv')
          do iqx=1,nqibz
            read(ippovlinv) qx,ngcx
            allocate( ppovl(1:ngcx,1:ngcx),ppovlinv(1:ngcx,1:ngcx))
            print *,' ngc ngcx=',ngc,ngcx
            read(ippovlinv) ppovl(1:ngcx,1:ngcx)
            read(ippovlinv) ppovlinv(1:ngcx,1:ngcx)
            if(sum (abs(q-qx))<1d-10) then
              print *, ' iqx=',iqx
              print *, qx
Cstop2rx 2013.08.09 kino              if(ngc/=ngcx) stop 'x0kf: ngx/=ngcx'
              if(ngc/=ngcx) call rx( 'x0kf: ngx/=ngcx')
              exit
            endif
            deallocate(ppovlinv,ppovl)
          enddo
        endif
      endif
ccccccccccccccccccccccccccccccccccccccccc

      allocate(cphik (nlmto,nband),cphikq(nlmto,nband))

c ------------------------------------------------
c      print *,' x0kf_v2z ---------'
      if(nolfco) then !ixc==2.or.ixc==22) then
        if(.not.tetra)
Cstop2rx 2013.08.09 kino     &    stop'  x0kf_v2z: tetra=t should be for this mode'
     &    call rx( '  x0kf_v2z: tetra=t should be for this mode')
        if(niwt/=0) 
Cstop2rx 2013.08.09 kino     &    stop ' x0kf_v2z: niwt==0 should be for this mode'
     &    call rx( ' x0kf_v2z: niwt==0 should be for this mode')
        if(.not.realomega) 
Cstop2rx 2013.08.09 kino     &    stop ' x0kf_v2z: realomega=t should be for this mode'
     &    call rx( ' x0kf_v2z: realomega=t should be for this mode')
      endif

      imdim(1) = 1
      do iatom = 1,natom
        iatomp(iatom) = iatom
        if(iatom<natom) imdim(iatom+1)=imdim(iatom)+mdim(iatom)
      enddo

c dimensions and constants
      nlmtoband  = nlmto*nband
      nctot      = noccx - noccxv
      call dinv33(qbas,0,qbasinv,det)
      coskt=1d0
      sinkt=0d0

ccccccccccccccccccccccccc
c      print *,' wksum=',sum(wk(:))
ccccccccccccccccccccccccc

C ---- overlaid ntp0 ---------- kq 21Jun2003 takao
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      call getnemx(nbmx,ebmx,7,.true.)

      call getkeyvalue("GWinput","ecut_pol",eee,2, default=(/1d10,1d10/) )
c        call getkeyvalue("GWinput","ecut_high",eee,2, default=(/1d10,1d10/) )
      ecut = eee(1)
      ecuts= eee(2)
      print *,' ecut ecuts=',ecut,ecuts

c loop over k-points
      do  1000  k = 1,nqbz
        if(nbnb(k)==0) cycle !May 2002.
c qk = k + q, k is inside 1st BZ, not restricted to

c      if(iq <=nqibz) then
cc the irreducible BZ
c        qk  = rk(:,k) + q
cc bring qk into the 1st BZ ==> qbz labelled by kp
cc      print *, qk,n1,n2,n3,qbas,ginv
ccc        call fbz   (qk,n1,n2,n3,qbas,ginv,w(iindxk),
ccc     o              qbkp,kp)
c         call fbz2(qk,  ginv,index_qbz,n_index_qbz,rk,nqbz, !qbz=rk
c     o              qbkp,kp)
ccc
c      else
c        kp   = nqbz *(iq - nqibz) + k
c        qbkp = rk(:,k) + q
c      endif

        if(k<=5) then
          write(6,"(' x0kf: k rk=',i4,3f10.4,$)")k, rk(:,k)
          call cputid(0)
        elseif(mod(k,10)==1.or.k==nqbz) then
          write(6,"(' x0kf: k rk=',i4,3f10.4)")k, rk(:,k)
        endif

c direct access read eigenvalues at k and k+q
c      call rwdd1   (ifev1,k, nband, ek)
c      call rwdd1   (ifev2,kp,nband,ekq)
        call readeval( rk(:,k)    ,isp1,  ek)
        call readeval( rk(:,k) + q,isp2, ekq)

        if(ipr) print *, '111111111111111111'
c direct access read b,hb at k and k+q
c      kisp= isp + nsp*(k-1)
c      call rwdd1   (ifcphi, kisp, 2*nlmto*noccxv, cphik)
        call readcphi(rk(:,k),   nlmto,isp1, quu, cphik)
        call readcphi(q+rk(:,k), nlmto,isp2, quu, cphikq)
        if(ipr) print *, '222222222222222'

c no. occupied (valence) and unoccupied states
c      nt0        = nocc (ek,ef,.true.,nlmto)
c      ntp0       = nocc (ekq,ef,.false.,nlmto)
        nt0        = nocc (ek, ef,.true. ,nband) !occ   (OUT) q2
        ntp0       = nocc (ekq,ef,.false.,nband) !unocc (IN)  q1
        ekc(1:nctot)           = ecore(1:nctot)
        ekc(nctot+1:nctot+nt0) = ek   (1:nt0)

c> starting index for unoccupied state
c      itps       = nlmto - ntp0 + 1
        itps       = nband - ntp0 + 1


c--- tetra ------ override nt0 itps ntp0 ---------------
        if(tetra) then                                  ! Tetrahedron case
          n1bmax  = 0; n2bmin = nband + nctot +1000000  ! for security
c        do iw   = 1, niwt+nwt
          do ibib = 1, nbnb(k)
            if( n1b(ibib,k)<=nband) n1bmax = max(n1b(ibib,k),n1bmax)
            n2bmin = min(n2b(ibib,k),n2bmin)
          enddo
c        enddo
          nt0  = n1bmax           ! n1bmax = the num of maxmum occ
          itps = n2bmin           ! n2bmin = the num of min unocc
          ntp0 = nband - (n2bmin -1)
c        if(mod(k,10)==1) print *,' k n1bmax n2bmin =', k, n1bmax,n2bmin
!       print *,' nt0 itps...=',nt0,ntp0,itps
          call readeval(q+rk(:,k),isp2,eband)
          n2bmax = nocc (eband,ebmx,.true.,nband)
          ntp0   = n2bmax - (itps -1)
          if(ntp0 <=0 .or.ipr) write(6,"('  n1bmax n2bmin n2bmax=',3i10)")n1bmax,n2bmin,n2bmax
          if(ntp0 <=0) cycle
        endif

        if(ipr) print *, ' 3333333333333333333333'

c calculate the matrix elements <psi(k+q,t') | psi(k,t) B(R,i)>

c> core
c      write (*,*) 'matrix elements'
c     call cputid(0)
        allocate( zzmel(nbloch,noccx,nband) )
c---

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(zvztest()==2.and.k==2) then
c        call readcphi(rk(:,k), nlmto,isp, quu, cphikq)
          cphik= cphikq
          zzmel=0d0
          itps=1
          nt0 =nband
          ntp0=nband
          noccxv = nband
          noccx  = nband +nctot
          deallocate(zzmel)
          allocate( zzmel(nbloch,noccx,nband) )
          print *,' zvztest=2 nctot=',nctot
        endif

cccccccccccccccc
        if(zvztest()/=0) then
          if(k==1) allocate(m_in(nbloch),n_in(nbloch),l_in(nbloch),ibas_in(nbloch))
          ifmn=4301
          open(ifmn,file='@MNLA_CPHI')
          read(ifmn,*)
          ix =1
          do
            read(ifmn,*,err=4302,end=4302) m_in(ix),n_in(ix),l_in(ix),ibas_in(ix)
            write(6,"(5i5)") m_in(ix),n_in(ix),l_in(ix),ibas_in(ix),ix
            ix=ix+1
          enddo
 4302     continue
          close(ifmn)
Cstop2rx 2013.08.09 kino          if(ix-1/=nlmto) stop 'readin @MNLA_CPHI but nlmto /= readin num'
          if(ix-1/=nlmto) call rx( 'readin @MNLA_CPHI but nlmto /= readin num')
          do ix=1,nlmto
            if(l_in(ix)>1) cphik (ix,1:nband)=0d0
            if(l_in(ix)>1) cphikq(ix,1:nband)=0d0
          enddo
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! The difference between psicb_x0 and psicb is imdim,iatomp for rotation.
c      call psicb_x0k_v2  (w(iicore),w(incore),ntp0,w(iiclass),
c     i              coskt,sinkt,
c     i              cphikq(1,itps), !rbkq(1,itps),cbkq(1,itps),rhbkq(1,itps),chbkq(1,itps),
c     i              ppb,            !,pdb,dpb,ddb,
c     i              w(inlnmv),w(inlnmc),mdim,
c     d              mdimx,nlmto,nbloch,nlnmx,noccx,ntp0,natom,nclass,
c     d              nl,nnc,
c     o              zzmel)
c> valence
c      call psi2b_x0k_v2  (nt0,ntp0,w(iiclass),coskt,sinkt,
c     i              cphik,         !rbk,cbk,rhbk,chbk,
c     i              cphikq(1,itps),!rbkq(1,itps),cbkq(1,itps),
c     i                             ! rhbkq(1,itps),chbkq(1,itps),
c     i              ppb,           !,pdb,dpb,ddb,
c     i              w(inlnmv),w(inlnmc),mdim,nctot,
cc     w              w(iwr1),w(iwr2),w(iwr3),w(iwr4),
cc     w              w(iwc1),w(iwc2),w(iwc3),w(iwc4),
c     d              mdimx,nlmto,nbloch,nlnmx,noccxv,noccx,ntp0,
c     d              natom,nclass,
c     o              zzmel)

        call psicb_v2  (icore,ncore,ntp0,iclass,
     i              coskt,sinkt,
     i              cphikq(1,itps), !rbkq(1,itps),cbkq(1,itps),rhbkq(1,itps),chbkq(1,itps),
     i              ppb,            !,pdb,dpb,ddb,
     i              nlnmv, nlnmc, mdim, 
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,noccx,ntp0,natom,nclass,
     d              nl,nnc,
     o              zzmel)
c> valence
        call psi2b_v2  (nt0,ntp0,iclass,coskt,sinkt,
     i              cphik,         ! rbk,cbk,rhbk,chbk,
     i              cphikq(1,itps),! rbkq(1,itps),cbkq(1,itps),rhbkq(1,itps),chbkq(1,itps),
     i              ppb,           !,pdb,dpb,ddb,
     i              nlnmv,nlnmc,mdim,nctot,
     i              imdim,iatomp,
     d              mdimx,nlmto,nbloch,nlnmx,noccxv,noccx,ntp0,
     d              natom,nclass,
     o              zzmel)
        if(ipr) print *, '4444444444444444444'


cccccccccccccccccccccccccccccccccccccccc
        if(zvztest()==2) then
c      allocate( vzz(ngb, nctot+nt0,ntp0) ,aaa(ngc))
c      call matm(vcoul,zmelt,  vzz,ngb,ngb, (nctot+nt0)*ntp0)
          ifix=9090
          open (ifix,  file='zvz_test_Bphiphi')
c      write(ifix,
c     & "(' ib1 ib2     <B psi_ib1 psi_ib2>')")
          write(6,"(' q+rk q+rk=',3f8.5,2x,3f8.5)") q+rk(:,k),q+rk(:,k)
c      write(6,"('noccx nctot nt0=',3i5)") noccx,nctot,nt0
c        if(noccx< nctot+nt0) then
c         stop 'x0kf:zvztest2:noccx < nctot+nt0'
c        endif
ccc     zzmel(1:nbloch, 1:nctot+nt0, 1:ntp0)
          do ib1=1,noccxv
            ib2=ib1
            zvz = zzmel(1,ib1,ib2)*zzmel(1,1,1)
            write(ifix,  "(2i4,  2f13.8,'  ',2f13.8,'  ',2f13.8)") ib1,ib2,zvz
          enddo
          close(ifix)
c      deallocate(vzz,aaa)
Cstop2rx 2013.08.09 kino          if( k==2) stop '--- zvz test mode=2 end------------------' !Takeo k==2 case
          if( k==2) call rx( '--- zvz test mode=2 end------------------')
        endif


cccccccccccccccccccccccccccccccccccccccc
c      q =(/0.05d0,0d0,0d0/)
c      rk(:,k)=(/0.125d0,0d0,-0.125d0/)
cccccccccccccccccccccccccccccccccccccc
        if(ipr) print *, 'goto readqg'
c plane ---------------------------------------------------------------
        call readqg('QGpsi',q+rk(:,k),ginv, qu1x,ngp1, ngvecpB1)
        call readqg('QGpsi',  rk(:,k),ginv, qu2x,ngp2, ngvecpB2)
        if(ipr) print *, 'end readqg'

c      ngp1 = ngpn(kp)  ! q+k   ntp0 in FBZ
c      ngp2 = ngpn(k)   ! k     np0  in FBZ
!     ngc              ! q          in IBZ

        ngb  = nbloch + ngc
        allocate( zmelt(ngb,nctot+nt0,ntp0))

#ifndef USE_X0KBLAS
        allocate(z1p(ngb,ngb) )
#endif

c     z1p =1d99
c zmelpl(ngc,nt0,ntp0),
c     &          rmelt(ngb,1:nctot+nt0,1:ntp0),
c     &          cmelt(ngb,1:nctot+nt0,1:ntp0) )
        if(ipr) print *, '55555555555555'

        call readgeig(q+rk(:,k), ngpmx, isp2, qu1, geig1)
cccccccccccccc
c      write(6,"('uuu sumtest g1a=',6f8.4,2i4,2d13.5)")
c     &  q+rk(:,k), qu1,itps,itps+ntp0-1, sum(geig1(1:ngp1, itps:itps+ntp0-1))
cccccccccccc
        call readgeig(  rk(:,k), ngpmx, isp1, qu2, geig2)
        if(sum(abs(qu1-qu1x))>1d-10) then
          print *,' qu1=',qu1
          print *,' qu1x=',qu1x
Cstop2rx 2013.08.09 kino          stop 'x0kf_v2z;qu1/=qu1x'
          call rx( 'x0kf_v2z;qu1/=qu1x')
        endif
        if(sum(abs(qu2-qu2x))>1d-10) then
          print *,' qu2=',qu2
          print *,' qu2x=',qu2x
Cstop2rx 2013.08.09 kino          stop 'x0kf_v2z;qu2/=qu2x'
          call rx( 'x0kf_v2z;qu2/=qu2x')
        endif
cccccccccccccccccccccccc
c         print *,' qu1=',qu1
c         print *,' qu1x=',qu1x
c         stop 'test end -------------------'
cccccccccccccccccccccccc
        qdiff= q+qu2 - qu1
c      qdiff = q   - qbkp(:) + rk(:,k)
c            ! q   - (q+k)   + k  is not zero.
c            ! qc  -  q1     + q2
        add   = matmul(qbasinv, qdiff)
        nadd  = idint( add + dsign(.5d0,add))
Cstop2rx 2013.08.09 kino        if(sum(abs(add-nadd))>1d-10) stop "sexc: abs(add-nadd))>1d-10"
        if(sum(abs(add-nadd))>1d-10) call rx( "sexc: abs(add-nadd))>1d-10")
c     print *,' qdif=',qdiff,qbkp(:),rk(:,k)

        if(ipr) print *, '6666666666666'

        zmelt = 0d0
c      call melpln2( ngp1, ngvecpB(:,1:ngp1, kp)  ! q1=q+k  ; kp ngp1 1:ntp0 q-point
c     &           ,  ngp2, ngvecpB(:,1:ngp2, k)    ! q2=k    ; k  ngp2 1:nt0  occupied
c     &           ,  ngc,  ngveccB,  nadd,
c     &       geigB(1:ngp1, itps:nband,kp), ntp0, ! q1=q+k  ; kp
c     &       geigB(1:ngp2, 1:nt0,    k ) ,  nt0, ! q2=k    ; k
        call melpln2( ngp1,  ngvecpB1  ! q1=q+k  ; kp ngp1 1:ntp0 q-point
     &           ,  ngp2,  ngvecpB2  ! q2=k    ; k  ngp2 1:nt0  occupied
c     &           ,  ngc,   ngveccB,  nadd,
     &           ,  ngc,   nadd,
     &       geig1(1:ngp1, itps:nband), ntp0, ! q1=q+k  ; kp
     &       geig2(1:ngp2, 1:nt0     ),  nt0, ! q2=k    ; k
     i       shtv,  q, q,  symope,qbas,
c extensiton to nbloch+ngc       
     o       zmelt (nbloch+1:nbloch+ngc, nctot+1:nctot+nt0,1:ntp0))


        if(ipr) print *, '77777777777777'

!                                              k            q+k
cccccccccccccccccccccccccccccccccccc
c      zmelt=0d0
cccccccccccccccccccccccccccccccccccc
        zmelt        (1:nbloch, 1:nctot+nt0, 1:ntp0) = 
     &        zzmel(1:nbloch, 1:nctot+nt0, 1:ntp0)
c     &dcmplx(  rmel(1:nbloch, 1:nctot+nt0, 1:ntp0),
c     &         cmel(1:nbloch, 1:nctot+nt0, 1:ntp0) )
!                               k            q+k
!      deallocate(zmelpl)
        deallocate(zzmel)

cccccccccccccccccccccccccccccccc
c      print *,' sumchk nnn zmelt=',nctot+nt0, ntp0,sum(zmelt(1:ngb,1:nctot+nt0, 1:ntp0))
ccccccccccccccccccccccccccccccccccccc


c matrix element test case
        if(zvztest()==1) then
          allocate( vzz(ngb, nctot+nt0,ntp0) ,aaa(ngc))
          call matm(vcoul,zmelt,  vzz,ngb,ngb, (nctot+nt0)*ntp0)
          ifix=9090
          open (ifix,  file='zvz_test')
          write(ifix,  
     & "(' ib1 ib2     Total  =   MT  + (MT -Total)  :Coulomb mat')")
          open (ifix+1,file='zvz_test2')
          write(ifix+1,
     & "(' ib1 ib2     Total  =   MT  + Interstitial :overlap mat')")
          print *,' ======= goto zvz test  ========= rk '
          write(6,"(' rk q+rk=',3f8.5,2x,3f8.5)") rk(:,k),q+rk(:,k)
          write(ifix+1,"(' rk q+rk=',3f8.5,2x,3f8.5)") rk(:,k),q+rk(:,k)
          do ib1=1,nctot+nt0
            do ib2=1,ntp0
! zmelt*vcoul
              zvz1= sum(zmelt(1:nbloch,ib1,ib2)*matmul(vcoul(1:nbloch,1:nbloch),dconjg(zmelt(1:nbloch,ib1,ib2))) )
              zvz= sum((zmelt(:,ib1,ib2))*matmul(vcoul,dconjg(zmelt(:,ib1,ib2))) )
              write(ifix,  "(2i4,  2f13.8,'  ',2f13.8,'  ',2f13.8)") ib1,ib2+itps-1,zvz ,zvz1, zvz-zvz1

              zvz1= sum(zmelt(1:nbloch,ib1,ib2)*dconjg(zmelt(1:nbloch,ib1,ib2)))
              aaa = matmul(ppovl, dconjg(zmelt(nbloch+1:ngb,ib1,ib2)) )
              zvz = sum(zmelt(nbloch+1:ngb,ib1,ib2)*aaa(1:ngc))
              write(ifix+1,"(2i4,  2f13.8,'  ',2f13.8,'  ',2f13.8)") ib1,ib2+itps-1, zvz+zvz1, zvz1, zvz
            enddo
          enddo
          close(ifix)
          close(ifix+1)
          deallocate(vzz,aaa)
Cstop2rx 2013.08.09 kino          if(k==2) stop '--- zvz test mode=1 end------------------' !Takeo this is k==2 case
          if(k==2) call rx( '--- zvz test mode=1 end------------------')
        endif
c-end of matrix elements --------------------------------------------------


ccccccccccccccccccccccccccccccccccccccccccc
c A matrix element check in Si. See the corresponding section of sxcf.f
c      write(6,"(3f13.5)") qbkp(:) !ntp0
c      write(6,"(3f13.5)") q       !product basis
c      print *,' x0kf:sumcheck rmelt cmelt='
c     &   ,sum(zmelt(1:ngb,4,1))
c     &   ,sum(abs(dreal(zmelt(1:ngb,4,1))))
c     &   ,sum(abs(dimag(zmelt(1:ngb,4,1))))
c      stop 'test end'
c      stop ' a matrix ele check test end'
ccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccc
!test to masure CPU time !
! This test shows the lower part from here dominates CPU time.
c      deallocate(zmelt,z1p)
c      cycle
cccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' nt0 ntp0 itps nctot =', nt0,ntp0,itps,nctot
c      do  itp = 1,nband
c        print *, ' ekq x=',ekq(itp)
c      enddo
c      do  itp = 1,ntp0
c        print *, ' ekq=',ekq(itps+itp-1)
c      enddo
c      do  it = 1,nt0+nctot
c        print *, ' ekc=',ekc(it)
c     enddo
c     write (*,*) 'sum over states'
c     call cputid(0)
cccccccccccccccccccccccccccccccccccccccccccccccc

        if(tetra) then
          if(nolfco) then !ixc==2.or.ixc==22) then ! for <exp(iq r)|x0(q,\omega)|exp(iq r)>
            do ibib = 1, nbnb(k) !---  ibib loop
c            if(ipr) print *,'tetra nolfco: xxx1 ibib =',ibib
c--- nbcut
              if( n2b(ibib,k) > nbmx) then
                if(iww1) then
                  print *,' nband_chi0 nbmx=',nbmx
                  iww1=.false.
                endif
                cycle
              endif
              if( n1b(ibib,k) <= nbcut.and. n2b(ibib,k)>nbcut2) then
                if(iww2) then
                  write(6,"(' nbcut nbcut2 n2b n1b=',4i6)") nbcut,nbcut2,n2b(ibib,k),n1b(ibib,k)
                  iww2=.false.
                endif
                cycle
              endif
c--- ecut Jan2005
              wcut = 1d0
              if(ecut<1d9) then
c            wcut = wcutef(ekq(n2b(ibib,k))-ef, ecut,ecuts)
                wcut = wcutef(ekq(n2b(ibib,k))- ek(n1b(ibib,k)), ecut,ecuts)
                if(wcut <1d-6) cycle
              endif
              if(n1b(ibib,k) <= nband) then
                it = nctot + n1b(ibib,k) !valence
              else
                it = n1b(ibib,k) - nband !core
              endif
              itp = n2b(ibib,k) - itps + 1
              if(itp> n2bmax-itps+1) cycle !takao July2005
ccccccccccccccccccccc
c            if(ipr) then
c               print *
c               print *,' ttt tetra nolfco: xxx2 ibib =',ibib
c               print *,' ttt sumchk zzr=',sum(zzr(1:ngb)),it,itp
c               print *, 'ttt ngb,nctot+nt0,ntp0=',ngb,nctot+nt0,ntp0
c               print *,' ttt sumchk zmelt=',sum(zmelt(1:ngb,it,itp))
c            endif
ccccccccccccccccccccc
              zq0 = sum(zmelt(1:ngb,it,itp)*zzr(1:ngb))
              zq0zq0 = dconjg(zq0)*zq0
c$$$cccccccccccccccccccccccccccc
c$$$            if(k/=10) cycle
c$$$            write(6,"('uuu',14i5)")k,ibib,n1b(ibib,k),n2b(ibib,k),it,itp
c$$$            write(6,"('uuu zq0zq0 sumzzr=',2d24.16)")zq0zq0,sum(abs(zzr(1:ngb)))
c$$$            write(6,"('uuu sumzmelt=',2d24.16)") sum(zmelt(1:ngb,it,itp))
c$$$            write(6,"('uuu sumzmelt1=',2d24.16)") sum(zmelt(1:nbloch,it,itp))
c$$$            write(6,"('uuu sumzmelt2=',2d24.16)") sum(zmelt(nbloch+1:ngb,it,itp))
c$$$           write(6,"('uuu ng chk=',18i5)")
c$$$     &       ngp1,  sum(abs(ngvecpB1(:,1:ngp1))), ngp2,sum(abs(ngvecpB2(:,1:ngp2)))
c$$$     &       ,  ngc,   nadd,  ntp0, nt0
c$$$           write(6,"('uuu q chk=',3f9.4,9f9.4)") q(:), qbas(:,:)
c$$$           write(6,"('uuu itps ntp0=',18i5)") itps,itps+ntp0-1,nband
c$$$           write(6,"('uuu geig=',d23.16)")
c$$$     &        sum(abs(geig1(1:ngp1, itps:itps+ntp0-1))),
c$$$     &        sum(abs(geig2(1:ngp2, 1:nt0)))
c$$$cccccccccccccccccccccccccccc
              do iw  = 1, nwt  !----- iw loop
                x0mean(iw) = x0mean(iw) + wwk(ibib,k,iw)*zq0zq0
              enddo ! ----- iw   loop
            enddo ! ----- ibib loop
            deallocate(zmelt)
            cycle
          endif

c  z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
c  zxq(iw,ibg1,igb2) = sum_ibib wwk(iw,ibib)* z1p(ibib, igb1,igb2)
c
c
#ifdef USE_X0KBLAS
          allocate( wwkr( nbnb(k), nwt) , wwki( nbnb(k), niwt) )
#ifdef X0KBLAS_DIV
c       ngbblk = 2
c Mbytes_X0KBLAS_DIV is the available memory for allocating ngbblock.
          ngbdiv = 16d0* dble(nbnb(k))*dble(ngb)**2
     &         /(1024d0**2*Mbytes_X0KBLAS_DIV)+1
          ngbblk = min( ngb/ngbdiv+1, ngb )
          allocate( z1pp( nbnb(k), ngb,ngbblk) )
          if(k==1) write(*,'(a,I7,a,I7,2i7,3d14.6)')
     &     'x0kf: tetrahedron+blas3 mode, ngbblk/ngb ngbdiv=',
     &     ngbblk,'/', ngb, ngbdiv,nbnb(k)
#else
          allocate( z1pp( nbnb(k), ngb,ngb) )
#endif

          zalpha=1.0d0
          zbeta = 1.0d0
          do iw=1,nwt
            wwkr(:,iw) = wwk(:,k,iw+niwt)
          enddo
          do iw=1,niwt
            wwki(:,iw) = dreal(wwk(:,k,iw))
          enddo

#ifdef X0KBLAS_DIV
          do igb2=1,ngb, ngbblk
            igb2end = min(ngb, igb2+ngbblk-1)
            do ibib = 1,nbnb(k)
c--- nbcut
              if( n2b(ibib,k) > nbmx) then
                if(iww1) then
                  print *,' nband_chi0 nbmx=',nbmx
                  iww1=.false.
                endif
                cycle
              endif
              if( n1b(ibib,k) <= nbcut.and. n2b(ibib,k)>nbcut2) then
                if(iww2) then
                  write(6,"(' nbcut nbcut2 n2b n1b=',4i6)") nbcut,nbcut2,n2b(ibib,k),n1b(ibib,k)
                  iww2=.false.
                endif
                cycle
              endif
c--- ecut Jan2005
              wcut = 1d0
              if(ecut<1d9) then
c            wcut = wcutef(ekq(n2b(ibib,k))-ef, ecut,ecuts)
                wcut = wcutef(ekq(n2b(ibib,k))- ek(n1b(ibib,k)), ecut,ecuts)
                if(wcut <1d-6) cycle
              endif
c---
              if(n1b(ibib,k) <= nband) then
                it = nctot + n1b(ibib,k) !valence
              else
                it = n1b(ibib,k) - nband !core
              endif
              itp  = n2b(ibib,k) - itps + 1
              if(itp> n2bmax-itps+1) cycle !takao Nov2004


              do igb2start = igb2, igb2end
                do igb1=1,ngb
                  z1pp(ibib,igb1,igb2start-igb2+1) = 
     &        dconjg(zmelt(igb1,it,itp)) * zmelt(igb2start,it,itp)*wcut !Jan2005
                enddo
              enddo
            enddo
            if ( imagomega )
     &    call zgemm( 't','n', ngb*(igb2end-igb2+1), niwt, nbnb(k) ,
     &     zalpha, z1pp, nbnb(k), wwki, nbnb(k),
     &     zbeta,  zxqi(1,igb2,1), ngb**2 )
            if (realomega)
     &    call zgemm( 't','n', ngb*(igb2end-igb2+1), nwt, nbnb(k) ,
     &     zalpha, z1pp, nbnb(k), wwkr, nbnb(k),
     &     zbeta,  zxq(1,igb2,1), ngb**2 )
          enddo

#else
          do ibib=1,nbnb(k)
c--- nbcut Oct2004
            if( n2b(ibib,k) > nbmx) then
              if(iww1) then
                print *,' nband_chi0 nbmx=',nbmx
                iww1=.false.
              endif
              cycle
            endif
            if( n1b(ibib,k) <= nbcut.and. n2b(ibib,k)>nbcut2) then
              if(iww2) then
                write(6,"(' nbcut nbcut2 n2b n1b=',4i6)") nbcut,nbcut2,n2b(ibib,k),n1b(ibib,k)
                iww2=.false.
              endif
              cycle
            endif
c--- ecut Jan2005
            wcut = 1d0
            if(ecut<1d9) then
c            wcut = wcutef(ekq(n2b(ibib,k))-ef, ecut,ecuts)
              wcut = wcutef(ekq(n2b(ibib,k))- ekq(n1b(ibib,k)), ecut,ecuts)
              if(wcut <1d-6) cycle
            endif
c
            if(n1b(ibib,k) <= nband) then
              it = nctot + n1b(ibib,k) !valence
            else
              it = n1b(ibib,k) - nband !core
            endif
            itp  = n2b(ibib,k) - itps + 1

            do igb2=1,ngb
              do igb1=1,ngb
                z1pp(ibib,igb1,igb2) = dconjg(zmelt(igb1,it,itp)) 
     &            * zmelt(igb2,it,itp)*wcut !jan2005
              enddo
            enddo

          enddo

          if ( imagomega ) 
     &   call zgemm( 't','n',   ngb**2  , niwt, nbnb(k) ,
     &    zalpha, z1pp, nbnb(k), wwki, nbnb(k),  
     &    zbeta,  zxqi, ngb**2 )

          if (realomega) 
     &   call zgemm( 't','n',    ngb**2  , nwt, nbnb(k) ,
     &    zalpha, z1pp, nbnb(k), wwkr, nbnb(k),  
     &    zbeta,  zxq, ngb**2 )
#endif

          deallocate( z1pp )
          deallocate( wwkr, wwki )


#else
          do 30 ibib = 1, nbnb(k) !---  ibib loop
c--- nbcut Oct2004
            if( n2b(ibib,k) > nbmx) then
              if(iww1) then
                print *,' nband_chi0 nbmx=',nbmx
                iww1=.false.
              endif
              cycle
            endif
            if( n1b(ibib,k) <= nbcut.and. n2b(ibib,k)>nbcut2) then
              if(iww2) then
                write(6,"(' nbcut nbcut2 n2b n1b=',4i6)") nbcut,nbcut2,n2b(ibib,k),n1b(ibib,k)
                iww2=.false.
              endif
              cycle
            endif
c--- ecut Jan2005
            wcut = 1d0
            if(ecut<1d9) then
c            wcut = wcutef(ekq(n2b(ibib,k))-ef, ecut,ecuts)
              wcut = wcutef(ekq(n2b(ibib,k))- ekq(n1b(ibib,k)), ecut,ecuts)
              if(wcut <1d-6) cycle
            endif

          do 30 iw   = 1, niwt + nwt  !----- iw loop
            rviw  = dreal(wwk(ibib,k,iw))
            zviw  = wwk(ibib,k,iw)

            if(n1b(ibib,k) <= nband) then
              it = nctot + n1b(ibib,k) !valence
            else
              it = n1b(ibib,k) - nband !core
            endif
            itp  = n2b(ibib,k) - itps + 1
c ---
            if(iw==1) then !Assume all it itp are the same for all iw. ! See chknbnb in hx0fp0
              do ic = 1,ngb
                z1p(1:ngb,ic) =  
     &      zmelt(ic,it,itp)*dconjg(zmelt(1:ngb,it,itp))*wcut !jan2005
              end do
            endif
c
!!! this part dominates the cpu time !!!!!!!!!!!!!!!!!!!!!
            if(iw <= niwt.and.imagomega) then  !imag axis
              ! zxqi(1:ngb,1:ngb,iw) = zxqi(1:ngb,1:ngb,iw)+ rviw*z1p(1:ngb,1:ngb)
              call zaddr(zxqi(1,1,iw),rviw,z1p,ngb**2)
            endif
            if(iw > niwt.and.realomega) then  !real axis
              ! zxq(1:ngb,1:ngb,iw-niwt) = zxq(1:ngb,1:ngb,iw-niwt) + zviw*z1p(1:ngb,1:ngb)
              call zadd(zxq(1,1,iw-niwt),zviw,z1p,ngb**2)
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   30     continue ! ----- iw ibib loop
#endif
        else
! occ and unocc loop
          do 10 it  = 1,nt0+nctot   ! loop over occupied states
          do 10 itp = 1,ntp0        ! loop over unoccupied states
            de  = .5d0*( ekq(itps+itp-1) - ekc(it) ) ! energy in Hartree

cccccccccccccccccccccccccccccccccccccccccccccccccc
c test
c      de=0.5d0*0.5d0
cccccccccccccccccccccccccccccccccccccccccccccccccc

c frequency independent part
c> loop over basis
            do ic = 1,ngb
              z1p(1:ngb,ic) =  
     &    wk(k)* zmelt(ic,it,itp)*dconjg(zmelt(1:ngb,it,itp))
            end do
c
            if (realomega) then  !real axis
              do iw = 1,nwt
                call vw    (freqr(iw),de,delta,
     o              rviw,cviw)
                zviw = dcmplx(rviw,cviw)
!!! this part dominates the cpu time !!!!!!!!!!!!!!!!!!!!!
c1          zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)
c1     &        + zviw*z1p(1:ngb,1:ngb)
c3          zxq(:,:,iw) = zxq(:,:,iw) + zviw*z1p
                call zadd(zxq(1,1,iw),zviw,z1p,ngb**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              enddo
            endif
            if (imagomega) then  ! img axis
              do  iw = 1,niwt
                rviw     = -2d0*pilrntz(freqi(iw),de)
!!! this part dominates the cpu time !!!!!!!!!!!!!!!!!!!!!
c1          zxqi(1:ngb,1:ngb,iw) = zxqi(1:ngb,1:ngb,iw)
c1     &       + rviw*z1p(1:ngb,1:ngb)
c2          call daxpy(2*ngb**2,rviw,z1p,1,zxqi(1,1,iw),1)
c3          zxqi(:,:,iw) = zxqi(:,:,iw) + rviw*z1p
                call zaddr(zxqi(1,1,iw),rviw,z1p,ngb**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              enddo
            endif
   10     continue !!end of unoccupied-loop! end of occupied-loo
        endif


c-----------------------------------------
        if(.false.) then
          print *,' ---- trace chck ---- '
c      ippovl=221
c      open(ippovl,file= "ppovl",form='unformatted')
c     allocate(ttx(ngc,ngc))
c      do
c        read(ippovl) qx,ngcx
c       allocate(ppovl(ngcx,ngcx))
c        read(ippovl) ppovl
c        print *, ' qx=',qx
c        print *, ' q =',q
c       if (sum(abs(q-qx))<1d-10) exit
c       deallocate(ppovl)
c     enddo
c      close(ippovl)
c      print *,' ---- readed ppovl ---- '
c      if(ngcx/=ngc) then
c       print *, ngc ngcx
c        stop 'ngcx/=ngc'
c      endif
c       nn1= nbloch+1
c       nn2= ngb
          do iw=1,nwt
            trc = 0d0
            do ir=1,nbloch
              trc = trc + zxq(ir,ir,iw)
            enddo
c       ttx = matmul( zxq(nn1:nn2,nn1:nn2,iw), ppovl)
c        do ir=1,ngc
c        trc = trc + ttx(ir,ir)
c        enddo
            write(6,"(' iw  trc=',i4,2d14.6)")iw,trc
          enddo
c
          do iw=1,niwt
            trc = 0d0
            do ir=1,nbloch
              trc = trc + zxqi(ir,ir,iw)
            enddo
c       ttx = matmul(
c     &   zxqi(nn1:nn2,nn1:nn2,iw), ppovl)
c        do ir=1,ngc
c        trc = trc + ttx(ir,ir)
c        enddo
            write(6,"(' iwi trc=',i4,2d14.6)")iw,trc
          enddo
c     deallocate(ttx,ppovl)
c      stop ' ---- trace chck ---- end test '
        endif
c-----------------------------------------
        deallocate(zmelt)
#ifndef USE_X0KBLAS
        deallocate(z1p)
#endif
 1000 continue
c      call rlse (iwr1)
      end

c --------------------------------------------
      subroutine zadd(a,fac,b,ldim)
      complex(8) :: a(ldim),fac,b(ldim)
      do l=1,ldim
        a(l) = a(l) + fac*b(l)
      enddo
      end
      subroutine zaddr(a,fac,b,ldim)
      complex(8) :: a(ldim),b(ldim)
      real(8) :: fac
      do l=1,ldim
        a(l) = a(l) + fac*b(l)
      enddo
      end
c--------------------------------------------------------------------
ctaken from x0.f
      subroutine vw (w,de,delta,
     o               rviw,cviw)

c calculates rviw + i*cviw
c = 1/(w-de+ i*delta)  -  1/(w+de - i*delta)

      implicit double precision (a-h,o-z)

      wde        = w - de
      delta2     = delta*delta
      denom      = 1.d0/( wde*wde + delta2 )
      rviw       = denom*wde
      cviw       =-denom*delta

      wde        = w + de
      denom      = 1.d0/( wde*wde + delta2 )
      rviw       = rviw - denom*wde
      cviw       = cviw - denom*delta

      return
      end

      real(8) function wcutef(e,ecut,ecuts)
      real(8):: e,ecut,ecuts
c      wcutef = 1d0/( exp((e-ecut)/ecuts)+ 1d0)
      wcutef = exp( -(e/ecut)**2 ) ! ecuts is not used in this case
      end
