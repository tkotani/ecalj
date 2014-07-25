      program hnocc_mlw
c-------------------------------------------------------------
c Calculates the Fermi level and # occ. el. in maxloc states
c
c Aug 19, 2004 Takashi Miyake
c------------------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
      implicit none
c------------------------------------
      real(8)    :: esmr2,shtw
c      integer(4) :: mxclass,ngnmax,mbytes,mwords,iwksize,
c     &   natom,nclass,ipos,igrp,
c     &   iqibz,
c     &   iqbz,
c     &   iinvg,
c     o   nspin,nl,nn,nnv,nnc,
c     o   inindx,inindxv,inindxc,iiclass,
c     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o   iz, 
c     o   iil,iin,iim,iilnm,i_mnl,
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,
c     o   incwf,iecore,ikonf,iicore,incore,nctot,
c     o   imagw,niw,nw,ifreq,  
      integer(4)::
     &   ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,ifmloc,
     &   nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &   mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2)
c     &    ifev(2),
     &    ifwmat(2) !,ifcphi
     &   ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &   ifsecomg(2),ifexx,ndble=8
c
c      real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     & qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     & zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
c
c     &   lowesteval !defined in readeigen
ccc   qbas(3,3),ginv(3,3) plat(3,3),
      logical lqall,laf
c      character*120 symgrp

c class parameters
c      parameter (mxclass=100)
c      character*6 clabl(mxclass)
c symmetry group
c      parameter (ngnmax=10)
c      real(8) :: gen(9,ngnmax)
c
      integer(4),allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c takao
      integer(4),allocatable :: ngvecpB(:,:,:),!ngveccB(:,:,:),
     &  ngvecp(:,:), ngvecc(:,:),iqib(:), !,ngpn(:)ngcni(:)
     &  kount(:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:), !symope(:,:,:)=symgg, ! qibz(:,:),
     &  ppb(:), eq(:), !,pdb(:),dpb(:),ddb(:)
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
     &         , rw_w(:,:,:,:,:),cw_w(:,:,:,:,:),
     &           rw_iw(:,:,:,:,:),cw_iw(:,:,:,:,:)
      complex(8),allocatable:: geigB(:,:,:,:)
c 
      logical :: screen, exchange, cohtest, legas, tote 
      real(8) ::  rydberg,hartree 
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c space group infermation 
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c tetra 
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &      ,iene(:,:,:),ibzx(:) ! ,nstar(:)
c      real(8) :: qbasmc(3,3)

c worksize in megabytes (1 word = 4 bytes)
c      integer(4) :: mbytes,mwords,iwksize
c      parameter (mbytes=60)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)

      integer(4) ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum
c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
      integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     & ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
      character(3)  :: charnum3,sss
      character(12) :: filenameex
      logical :: exspwrite=.false.
      character*8 xt

      integer(4) :: iwini,iwend
c      real(8),allocatable:: omega(:)
      real(8) ::  omegamax,dwplot,omegamaxin
c      logical :: sergeys

      integer(4)::nqbze,ini,nq0it,idummy
      real(8),allocatable:: qbze(:,:)

      real(8)   :: ebmx
      integer(4):: nbmx
      
      real(8):: volwgt

      integer(4)::nwin, incwfin
      real(8)::efin,ddw
      integer(4),allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:),expa(:)

      logical:: GaussSmear !readgwinput, 
      integer(4)::ret
      character*(150):: ddd


      integer(4):: bzcase,  ngpn1,mrecg,verbose,ngcn1,nwxx
      real(8)   :: wgtq0p,quu(3)

      real(8),allocatable:: freq_r(:)

ctm
      real(8)   :: noccmlw1,noccmlw2
      integer(4):: nwf,nqbze2
c---------------------------------------
      hartree=2d0*rydberg()

      iii=verbose()
      print *,' verbose=',iii

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
ccccc    ! &   ,n_index_qbz  
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:)          !,index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      print *,' nqbz  =',nqbz
c      print *,  qbz 
      print *,' nqibz ngrp=',nqibz,ngrp
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

c set up work array
c      call wkinit (iwksize)
      call pshprt(60)


C--- readin GWIN and LMTO, then allocate and set datas.
      nwin =-999    !not readin NW file
      efin =-999d0  !not readin EFERMI
      incwfin= -1  !use 7th colmn for core at the end section of GWIN
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c---  These are allocated and setted.
c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw, !not readin nw
c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr   
c      character(120):: symgrp
c      character(6),allocatable :: clabl(:)
c      integer(4),allocatable:: iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)   
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c      real(8), allocatable:: 
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:),  symgg(:,:,:) ! symgg=w(igrp),freq(:)
c-----------------------------------------------------------------------

c--- Get maximums takao 18June03
      call getnemx(nbmx,ebmx,8,.true.) !8+1 th line of GWIN0

c-------------------------------------------------------------------
c      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      print *,' hsfp0: end of genallcf2'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

      shtw = 0d0
      if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)

      legas = .false.


c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen

c
      call getkeyvalue("GWinput","<MLWF>",unit=ifmloc,status=ret) 
      read(ifmloc,*)nwf
      close(ifmloc)

c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
! When esmr is negative, esmr is geven automatically by efsimplef.
        call efsimplef2a_mlw(nspin,wibz,qibz,ginv, 
     i        nband,nwf,nqibz
     i       ,konf,z,nl,natom,iclass,nclass
     i       ,valn, legas, esmr,! valn is input for legas=T, output otherwise. 
     i        qbz,nqbz !index_qbz, n_index_qbz,
     o       ,efnew,noccmlw1,noccmlw2)
c
c        print *,' end of efsimple'
c        ef = efnew
c- check total ele number -------
c        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
        ntot  = nocctotg2(nspin, efnew,esmr, qbz,wbz, nband,nqbz) !wbz
c        print *,' ef    =',ef
        print *,' ef    =',efnew
        print *,' esmr  =',esmr
        print *,' valn  =',valn
        print *,' ntot  =',ntot

c      ifcphi  = iopen('CPHI',0,0,mrecb)

      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen


      call wefn(ef,efnew,esmr,valn,noccmlw1,noccmlw2)

c------------
c close files
c------------
      isx = iclose ('hbe.d')

c-----------------------------------------------------------------------
      call cputid(0)
      stop ' OK! hnocc_mlw'
      end

c-----------------------------------------------------------------------
      subroutine efsimplef2a_mlw (ispin,wibz,qibz,ginv,
     i             nband,nwf,nqibz,
     i             konfig,z,nl,natom,iclass,nclass, 
     i             valn, legas, esmr, !input only for empty case
     i             qbz,nqbz, !index_qbz, n_index_qbz,
     o             ef,noccmlw1,noccmlw2)
      use m_readeigen
C- Calculate efermi for discrete sum. (not for tetrahedron method)
Cr user readeval (readeigen.f) to get eigenvalues.
Cr You need to call init_reaeigen before you call this.
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level 
c noccmlw = # occ. maxloc Wannier orbitals
c nband   = no. states
c nqbz    = no. k-points
c valn    = number of valence electron.
c
c ------------------- 
c      e(iband) < efm : occupation is one 
c efm< e(iband) < efp : occupation is wfacef. 
c efp< e(iband)       : occupation is zero 
c
      implicit none
      integer(4):: is,iq,ispin,nqibz,ik,isig,kpx,nband,ifev(2)
      integer(4):: ieaord(nband*nqibz*ispin),mbytes,mwords,iwksize,
     &             iqibz  !,iindxk 
      real(8)   :: ekt(nband, nqibz,ispin), ektx(nband*nqibz*ispin)
      real(8)   :: wgt(nband, nqibz,ispin), wgtx(nband*nqibz*ispin)
      real(8)   :: qbzx(3), qx(3),qbas(3,3),ginv(3,3),wwg
c
      integer(4):: nclass,natom,nl,ncore,l,ia,ic !,indxk(*)
      real(8)   :: wibz(nqibz),valn,ef, z(nclass),qibz(3,nqibz)
      integer(4):: iclass(natom),konfig(0:nl-1,nclass),ierr
c
      integer(4) :: nbnqnsp,ix,ikx,ikini,nne
      real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini
     &            ,eee2,wwg2 ,enumef
c      real(8) :: efp,efm,wwgo,wfacef
      logical :: legas,autoew,GaussSmear !is external

      integer(4):: nqbz
c      integer(4):: n_index_qbz,index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
      real(8)   :: qbz(3,nqbz)

ctm
      complex(8),allocatable:: evec(:,:,:,:),dnk(:,:,:,:)
      real(8),allocatable:: eval(:,:,:)
      real(8)   :: noccmlw1,noccmlw2
     &            ,enumef_gauss_mlw1,enumef_mlw1
     &            ,enumef_gauss_mlw2,enumef_mlw2
      integer(4):: nwf,ifenk,ifcnk,ib,
     &             nqbz2,nwf2,iko_ix,iko_fx,
     &             nqbze,nqbze2
C--------------------------------------------------------------------
      autoew =.false.
      if(GaussSmear()) then
      print *,' efsimplef2(gaussian mode):start'
      else
      print *,' efsimplef2(rectangular mode):start'
      endif
      if(esmr<=0d0) autoew=.true.
c total valence charge
      if(legas) then
        print *,' efsimplef2: legas=T use given valn = ',valn
      else 
        valn    = 0d0
        do ia   = 1,natom
          ic    = iclass(ia)
          valn  = valn + z(ic)
          print *,' ia z(ic)=',ia, z(ic)
        do    l = 0,nl-1
          print *,' l (konfig(l,ic)-l-1) 2*(2l+1)=',l,(konfig(l,ic)-l-1),( 2*l +1)*2
          valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
        end do
        end do
        print *,' valn=',valn
      endif
c
      do is = 1,ispin
      do iq = 1,nqibz
        call readeval(qibz(:,iq),is, ekt(1,iq,is) )
c        call catch1 (qibz,iq,3,nqibz,qx)
c         call fbz2(qx,  ginv,index_qbz,n_index_qbz,qbz,nqbz,
c     o              qbzx,kpx)
c        call rwdd1( ifev(is),   !read eigenvalues
c     d              kpx,nband,
c     o              ekt(1,iq,is) )
      enddo
      enddo
c
      if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
         print *, 'sum (wibz)=', sum(wibz(1:nqibz))
         stop 'efsimplef2: wibzsumerr'
      endif
      do is = 1,ispin
      do iq = 1,nqibz
        wgt(1:nband,iq,is) = wibz(iq)
        if(ispin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
      enddo
      enddo


c ekt and wgt 
      call dcopy ( nband*nqibz*ispin, ekt,1, ektx ,1)
      call dcopy ( nband*nqibz*ispin, wgt,1, wgtx ,1)
      call sortea( ektx,ieaord,nband*nqibz*ispin,isig)
c
c      print *,nband,nqibz,ispin,nband*nqibz*ispin
c      do ik  = 1, nband*nqibz*ispin
c        write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wgtx(ieaord(ik))
c      enddo
c
      open(8301,file = "DOSACC.lda")
      open(8302,file = "DOSACC2.lda")
      wwg = 0d0
      wwg2= 0d0
      eee2= -1d99
      ikini= 1
      ierr= 1
      nne = nband*nqibz*ispin
      do ik = 1, nne
c---
        if(eee2 +1d-4< ektx(ieaord(ik)).or.ik==nne ) then    
! degeneracy check 
          if(ik/=1) write(8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
          wwg2 = wgtx(ieaord(ik))
          eee2 = ektx(ieaord(ik))
          ikini =ik
        else
          wwg2= wwg2 + wgtx(ieaord(ik))
        endif  
c--- 
        wwg = wwg + wgtx(ieaord(ik))
        if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg
c
        write(8301,"(2i6,3d23.15)")
     &  ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))
c
        if( wwg>valn-1d-8.and.ierr==1 ) then
          write(6,*)
          efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
          if(autoew) then
             if(ik<3) stop ' efsimplef2: ik<3'
             esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
          endif
          ierr=0
        endif
        if( wwg > valn+1d0) ikx=ik
      enddo
c
      if(ierr==1) stop ' efsimplef2: ierr=1 given nval is too large'
c
      nbnqnsp = nband*nqibz*ispin
ccc gaussian
      if(GaussSmear()) then
      valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,efini,esmr,nbnqnsp)
      else
      valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,efini,esmr,nbnqnsp)
      endif
      print *, 'valx at efini=',efini,valx
      if(abs(valx-valn)<1d-8) then
        ef=efini
        goto 8891
      endif

      ew1= ektx(ieaord(1))-0.01d0
      ew2= ektx(ieaord(ikx))
      nbnqnsp =nband*nqibz*ispin
      do ix = 1,100
      ein = 0.5d0*(ew1+ew2)
      if(GaussSmear()) then
        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp))
     &          ,ektx(ieaord(1:nbnqnsp)),ein,esmr,nbnqnsp)  
      else
        valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))
     &               ,ein,esmr,nbnqnsp)  
      endif
        if(valx>valn) ew2=ein
        if(valx<valn) ew1=ein
        if(abs(ew1-ew2)<1d-15) exit
c        write(*,*)'valx at',ix,valx,ein
      enddo    
      ef = 0.5d0*(ew1+ew2)
c
 8891 write(8301,*) " ef=",ef
      close(8301)
      write(8302,*) " ef=",ef
      close(8302)

c debug:
c        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp))
c     &          ,ektx(ieaord(1:nbnqnsp)),ef,esmr,nbnqnsp)  
c      write(*,*)'valx, conv',valx,ef

ctm, for maxloc
c read data for maxloc states
      ifcnk = 8311
      ifenk = 8312
      do is = 1,ispin
         if (is.eq.1) then
            open(ifcnk,file='MLWU',form='unformatted')
            open(ifenk,file='MLWEU',form='unformatted')
         else
            open(ifcnk,file='MLWD',form='unformatted')
            open(ifenk,file='MLWED',form='unformatted')
         endif

         read(ifcnk)nqbz2,nwf2,iko_ix,iko_fx
         if (nqbz .ne. nqbz2) stop 'efsimple2a_mlw: nqbz error'
         if (nwf .ne. nwf2) stop 'efsimple2a_mlw: nwf error'

         if (is.eq.1) 
     &   allocate(dnk(iko_ix:iko_fx,nwf,nqibz,ispin))

         call r_dnk(is,ifcnk,iko_ix,iko_fx,
     i              ginv,qbz,qibz,
     d              nwf,nqbz,nqibz,ispin,
     o              dnk)

         read(ifcnk)nqbze,nwf2
         if (nwf2 .ne. nwf) stop 'efsimple2a_mlw: nwf error 2'
         read(ifenk)nqbze2,nwf2
         if (nqbze2 .ne. nqbze) stop 'efsimple2a_mlw: nqbze error'
         if (nwf2 .ne. nwf) stop 'efsimple2a_mlw: nwf error 3'
         if (is.eq.1) 
     &   allocate(evec(iko_ix:iko_fx,nwf,nqibz,ispin),
     &            eval(nwf,nqibz,ispin))

         call r_evalevec(is,ifenk,ifcnk,
     i                 ginv,qbz,qibz,
     d                 nwf,nqbz,nqibz,ispin,
     o                 eval,evec)

         close(ifcnk)
         close(ifenk)
      enddo

c debug:
c      do ik=1,nqibz
c      do ib=4,nwf+3
c         write(500,*)ib,ik,ekt(ib,ik,1)
c      enddo   
c      enddo   
c      do ik=1,nqibz
c      do ib=1,nwf
c         write(600,*)ib,ik,eval(ib,ik,1)
c      enddo   
c      enddo   

c noccmlw1: projected occupation number
      if(GaussSmear()) then
        noccmlw1= enumef_gauss_mlw1(dnk,iko_ix,iko_fx
     &          ,wgt,ekt,ef,esmr
     &          ,nband,nwf,nqibz,ispin)  
      else
        noccmlw1= enumef_mlw1(dnk,iko_ix,iko_fx
     &          ,wgt,ekt,ef,esmr
     &          ,nband,nwf,nqibz,ispin)  
      endif



c noccmlw2: occupation number in maxloc bands
      if(GaussSmear()) then
        noccmlw2= enumef_gauss_mlw2(eval,iko_ix,iko_fx
     &          ,wgt,ekt,ef,esmr
     &          ,nband,nwf,nqibz,ispin)  
      else
        noccmlw2= enumef_mlw2(eval,iko_ix,iko_fx
     &          ,wgt,ekt,ef,esmr
     &          ,nband,nwf,nqibz,ispin)  
      endif


      deallocate(dnk,eval,evec)


c      print *,' esmr        =',esmr
c      print *,' determined ef =',ef
c----------------------------------
c      wwg = 0d0
c      do ik  = 1, nband*nqibz*ispin
c       wwgo = wwg 
c       wwg  = wwg + wgtx(ieaord(ik))
c       if( abs(wwg-valn)<1d-6) then 
c         ef = 0.5d0*( ektx(ieaord(ik))+ektx(ieaord(ik+1)) )
c         efp = ef + 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
c         efm = ef - 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efm is just below the fermi
c        elseif(wwg>valn) then
c         ef      = ektx(ieaord(ik))
c         wfacef  = (valn-wwgo)/wgtx(ieaord(ik))
c         efp = ef + 0.5d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
c         efm = ef - 0.5d0*(ektx(ieaord(ik))-ektx(ieaord(ik-1)))  !efm is just below the fermi
c          print *,' determined ef    =',ef
c          print *,'            efp   =',efp
c          print *,'            efm   =',efm
c          print *,'           wfacef =',wfacef
c          return
c        endif 
c     enddo
      if(GaussSmear()) then
      print *,' efsimplef2(gauss):end'
      else
      print *,' efsimplef2:end'
      endif
      end
c----------------------------------------------------------------------------
      subroutine r_dnk(is,ifcnk,iko_ix,iko_fx,
     i                 ginv,qbz,qibz,
     d                 nwf,nqbz,nqibz,ispin,
     o                 dnk)

      implicit real*8(a-h,o-z)
 
      parameter (eps = 1d-8)
      complex(8), allocatable :: dnk_in(:,:,:)
      complex(8) :: dnk(iko_ix:iko_fx,nwf,nqibz,ispin)
      real(8) :: ginv(3,3)
      real(8) :: qbz(3,nqbz),qibz(3,nqibz),q(3)

      allocate(dnk_in(iko_ix:iko_fx,nwf,nqbz))

c read cnk
      do iqbz = 1,nqbz
         read(ifcnk)iqbz2,q(1:3)
         if (iqbz2 .ne. iqbz) stop 'r_dnk: iqbz error'
         read(ifcnk)dnk_in(iko_ix:iko_fx,1:nwf,iqbz)
      enddo

c
      do iqibz = 1,nqibz
         iqbz = iqindx(qibz(:,iqibz), ginv,qbz,nqbz)
         dnk(   iko_ix:iko_fx,1:nwf,iqibz,is)
     & = dnk_in(iko_ix:iko_fx,1:nwf,iqbz)
      enddo

      deallocate(dnk_in)

      end
c----------------------------------------------------------------------------
      subroutine r_evalevec(is,ifenk,ifcnk,
     i              ginv,qbz,qibz,
     d              nwf,nqbz,nqibz,ispin,
     o              eval,evec)

      implicit real*8(a-h,o-z)
 
      parameter (eps = 1d-8)
      complex(8), allocatable :: evec_in(:,:,:)
      complex(8) :: evec(nwf,nwf,nqibz,ispin)
      real(8), allocatable :: eval_in(:,:)
      real(8) :: ginv(3,3),eval(nwf,nqibz,ispin)
      real(8) :: qbz(3,nqbz),qibz(3,nqibz),q(3)

      allocate(evec_in(nwf,nwf,nqbz),eval_in(nwf,nqbz))

      do iqbz = 1,nqbz
         read(ifcnk)iqbz2,q(1:3)
         if (iqbz2 .ne. iqbz) stop 'r_evalevec: iqbz error'
         read(ifcnk)evec_in(1:nwf,1:nwf,iqbz)

         read(ifenk)iqbz2,q(1:3)
         if (iqbz2 .ne. iqbz) stop 'r_evalevec: iqbz error 2'
         read(ifenk)eval_in(1:nwf,iqbz)
      enddo

      do iqibz = 1,nqibz
         iqbz = iqindx(qibz(:,iqibz), ginv,qbz,nqbz)
         evec(:,:,iqibz,is) = evec_in(:,:,iqbz)
         eval(:,iqibz,is) = eval_in(:,iqbz)
      enddo

      deallocate(eval_in,evec_in)

      end
c----------------------------------------------------------------------------
      real(8) function enumef_gauss_mlw1(cnk,iko_ix,iko_fx
     &          ,wgt,ekt,ein,esmr
     &          ,nband,nwf,nqibz,ispin)
      implicit none
      integer(4):: nbnqnsp,ik,iko_ix,iko_fx,nband,nwf,nqibz,ispin,
     &           ibqs,ib,iq,is,iwf
      real(8) :: ekt(nband,nqibz,ispin),wgt(nband,nqibz,ispin),wwg,
     & derfcx,ein,esmr,weight
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqibz,ispin)


      wwg = 0d0
      do is = 1, ispin
      do iq = 1, nqibz
      do ib = iko_ix,iko_fx
         weight = 0d0
         do iwf = 1,nwf
            weight = weight + abs(cnk(ib,iwf,iq,is))**2
         enddo
         wwg= wwg + wgt(ib,iq,is)
     &       *0.5d0* derfcx( -(ein-ekt(ib,iq,is))/sqrt(2d0)/esmr )
     &       *weight
      enddo
      enddo
      enddo
      enumef_gauss_mlw1 = wwg
c     print *,' ein enumef=', ein, enumef

      end
c---------------------------------------------------------------
      real(8) function enumef_mlw1(cnk,iko_ix,iko_fx
     &               ,wgt,ekt,ein,esmr
     &               ,nband,nwf,nqibz,ispin)
      implicit real*8(a-h,o-z)
      integer(4):: iko_ix,iko_fx,nband,nwf,nqibz,ispin,
     &           ib,iq,is,iwf
      real(8) :: ekt(nband,nqibz,ispin),wgt(nband,nqibz,ispin)
      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqibz,ispin)


c     print *, esmr
      wwg = 0d0
      do is = 1, ispin
      do iq = 1, nqibz
      do ib = iko_ix,iko_fx
         weight = 0d0
         do iwf = 1,nwf
            weight = weight + abs(cnk(ib,iwf,iq,is))**2
         enddo

        if    (  ekt(ib,iq,is) + 0.5d0*esmr < ein ) then
          wwg  = wwg + wgt(ib,iq,is)*weight
        elseif(  ekt(ib,iq,is) - 0.5d0*esmr < ein ) then
          wwg  = wwg + wgt(ib,iq,is)
     &        *(ein- (ekt(ib,iq,is)-0.5d0*esmr))/esmr *weight
        endif
      enddo
      enddo
      enddo
      enumef_mlw1 = wwg
c     print *,' ein enumef=', ein, enumef

      end
c---------------------------------------------------------------
      real(8) function enumef_gauss_mlw2(enk,iko_ix,iko_fx
     &          ,wgt,ekt,ein,esmr
     &          ,nband,nwf,nqibz,ispin)
      implicit none
      integer(4):: nbnqnsp,ik,iko_ix,iko_fx,nband,nwf,nqibz,ispin,
     &           ibqs,ib,iq,is,iwf
      real(8) :: ekt(nband,nqibz,ispin),wgt(nband,nqibz,ispin),wwg,
     &           derfcx,ein,esmr,weight
     &          ,enk(nwf,nqibz,ispin)


      wwg = 0d0
      do is = 1, ispin
      do iq = 1, nqibz
      do ib = 1, nwf
         wwg= wwg + wgt(ib,iq,is)
     &       *0.5d0* derfcx( -(ein-enk(ib,iq,is))/sqrt(2d0)/esmr )
      enddo
      enddo
      enddo
c debug:
c      wwg = 0d0
c      do is = 1, ispin
c      do iq = 1, nqibz
c      do ib = 1, nband
c         wwg= wwg + wgt(ib,iq,is)
c     &       *0.5d0* derfcx( -(ein-ekt(ib,iq,is))/sqrt(2d0)/esmr )
c      enddo
c      enddo
c      enddo
      enumef_gauss_mlw2 = wwg
c     print *,' ein enumef=', ein, enumef

c debug:
c      do ik=1,nqibz
c      do ib=4,nwf+3
c         write(500,*)ib,ik,ekt(ib,ik,1)
c      enddo   
c      enddo   
c      do ik=1,nqibz
c      do ib=1,nwf
c         write(600,*)ib,ik,enk(ib,ik,1)
c      enddo   
c      enddo   

      end
c---------------------------------------------------------------
      real(8) function enumef_mlw2(enk,iko_ix,iko_fx
     &               ,wgt,ekt,ein,esmr
     &               ,nband,nwf,nqibz,ispin)
      implicit real*8(a-h,o-z)
      integer(4):: iko_ix,iko_fx,nband,nwf,nqibz,ispin,
     &           ib,iq,is,iwf
      real(8) :: ekt(nband,nqibz,ispin),wgt(nband,nqibz,ispin)
     &          ,enk(nwf,nqibz,ispin)


c     print *, esmr
      wwg = 0d0
      do is = 1, ispin
      do iq = 1, nqibz
      do ib = 1, nwf
        if    (  enk(ib,iq,is) + 0.5d0*esmr < ein ) then
          wwg  = wwg + wgt(ib,iq,is)
        elseif(  enk(ib,iq,is) - 0.5d0*esmr < ein ) then
          wwg  = wwg + wgt(ib,iq,is)*(ein- (enk(ib,iq,is)-0.5d0*esmr))/esmr
        endif
      enddo
      enddo
      enddo
      enumef_mlw2 = wwg
c     print *,' ein enumef=', ein, enumef

      end
c---------------------------------------------------------------
      subroutine wefn(ef,efnew,esmr,valn,noccmlw1,noccmlw2)
      implicit real*8(a-h,o-z)

      real(8) :: noccmlw1,noccmlw2

      open(8321,file='NOCC_MLW')

      write(8321,*)'*** efin'
      write(8321,*)ef
      write(8321,*)'*** ef, esmr'
      write(8321,*)efnew,esmr
      write(8321,*)'*** nocc(all)'
      write(8321,*)valn
      write(8321,*)'*** nocc(mlw)'
      write(8321,*)noccmlw1,noccmlw2

      close(8321)

      end
c---------------------------------------------------------------
