      program heftet
C - Calculates the Fermi energy by tetrahedron method.
c This routine is from hsfp0
c So it contains many unused variables.
c-------------------------------------------------------------
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use m_readeigen,only: init_readeigen,readeval
      implicit none
      real(8) :: shtw
      integer(4) :: mxclass,ngnmax,
     &   iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &   nprecx,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim,
     &   ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2),ifev(2),ifsec(2)
c     &   ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &    nq2,ntq2,
     &   ndble=8,ifev(2)
c
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &   qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     &   zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,efxx2,
c
     &   lowesteb

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
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:,:),
     &  ngvecp(:,:), ngvecc(:,:),ngpn(:),ngcni(:),iqib(:),
     &  kount(:,:), nx(:,:),nblocha(:),lx(:),ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:,:),
     &  ppb(:),pdb(:),dpb(:),ddb(:), eq(:,:),
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
c
      logical :: exchange, legas
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c space group infermation
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:),eband2(:,:,:), ene(:)  !ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &      ,iene(:,:,:),ibzx(:)
c      real(8) :: qbasmc(3,3)

c worksize in megabytes (1 word = 4 bytes)
c      parameter (mbytes=60)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)


      integer(4) ::ib,igp,iii,ivsumxxx,isx
c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,qq(3),
     &              elda, eqp01, eqp02
      integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     & ,ifexspxx , itx, iqx
      character(12) :: filenameex
      logical :: exspwrite=.false.

      integer(4) :: nptdos=101 !,iflegas

      logical :: metal      ,qbzreg
      integer(4) ::ldim,nbmax,nsp,nkp,
     &    itt,ifief ,imode,iftote,it
      real(8) :: zval,elo,ehi,e1,e2,efermi,dum,dosef,efxx,rydberg,dum111(1,1,1)
     &    ,tol=1d-10

      real(8)::volwgt
      real(8),allocatable:: bot(:),top(:),dos(:)
      integer(4):: nwin,incwfin
      real(8)::efin ,ddq(3),bandgap

      integer(4):: bzcase=1

c---------------------------------------
      write(6,*) 'mode=(1-4)?'
      read(5,*) imode
      call headver('heftet',imode)
      write(6,*) '--- heftet: calculation mode =',imode
Cstop2rx 2013.08.09 kino      if(imode<1.or.imode>4) stop 'mode out of range(1-4)'
      if(imode<1.or.imode>4) call rx( 'mode out of range(1-4)')
c---------------------------------------

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz
c       real(8):: qbas(3,3),ginv(3,3),dq_bzcase2
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      write(6,*)' nqbz qbz =',nqbz
c      write(6,*)  qbz
      write(6,*)' nqibz ngrp=',nqibz,ngrp2
c      write(6,*)' irk=',irk
      write(6,*)' #### idtetf: ####'
c      write(6,*) idtetf


cccc set up work array
cccc      call wkinit (iwksize)

c-------------------------------
c generate all quantities needed
c-------------------------------
      incwfin= -1  !use 7th colmn for core at the end section of GWIN
      nwin   = -999    !Not read NW
c      ifphiv(1)=-999 ! not readin ifphi
      efin =-999d0     ! not readin ef
C--- readin GWIN and LMTO, then allocate and set datas.
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
Cstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
c------------------------------
c      integer(4)::
c     o  nclass,natom,nspin,nl,nn,nnv,nnc,
c     o  ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o  nctot,
c     o  niw,nw
c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
c      character(120):: symgrp
c      character(6),pointer :: clabl(:)
c      integer(4),pointer::
c     &   iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:)
c     &  ,ncwf(:,:,:) ,
c     o    invg(:),
c     o    il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c      real(8), pointer::
c     o  plat(:,:),pos(:,:),z(:), symgg(:,:,:), !w(igrp)
c     o  ecore(:,:), freq(:)
c--------------------------------

c      allocate(ecore(nctot,nspin)) !core energies
c      do is = 1,nspin
c        if (nctot > 0) call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c      enddo

c-------------------------------------------------------------------
c      if (nclass > mxclass) stop ' heftet:increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
Cstop2rx 2013.08.09 kino      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      if (nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
      write(6,*)' hsfp0: end of genallcf2'
c
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat
c
      shtw = 0d0
      if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

c--- ef is taken as rs for the empty-sphere test case of legas=T case -------------
      legas = .false.
c$$$      INQUIRE (FILE = 'LEGAS', EXIST = legas)
c$$$      if(legas) then !!! test for electron gas case.
c$$$        write(6,*)' find LEGAS. legas =',legas
c$$$        iflegas = 2101
c$$$        open (iflegas,file='LEGAS')
c$$$        read(iflegas,*)rs
c$$$        close(iflegas)
c$$$        alpha = (9*pi/4d0)**(1d0/3d0)
c$$$        qfermi = alpha/rs
c$$$        efx  = qfermi**2
c$$$        valn = efx**1.5d0*voltot/3d0/pi**2
c$$$        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
c$$$        write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
c$$$        write (6,*)' egas  Exact Fermi energy    Ef  =', efx
c$$$      endif

c--------------------
      write(6, *) ' --- computational conditions --- '
      write(6,'("    deltaw  =",f13.6)') deltaw
      write(6,'("    esmr    =",f13.6)') esmr
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
c      write(6,'("    niw nw dw   =",2i3,f13.6)') niw,nw,dw

c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband
Cstop2rx 2013.08.09 kino      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'
      if (nprecb == 4) call rx( 'hsfp0: b,hb in single precision')
      isx = iclose('hbe.d')

C --- Locate band limits ---
c$$$      if(legas) then
c$$$        write(6,*)' efsimplef2: legas=T use given valn = ',valn
c$$$      else
        call getvaln(
     i       konf,z,nl,natom,iclass,nclass
     o      ,valn)
c$$$      endif
cccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccc
c      valn=29d0
c      write(6,*)' valn=',valn
cccccccccccccccccccccccccccccccccccccccccc

      if(imode==1) then
        call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
        ! EValue is readin in init_readeigen
        nsp   = nspin !26jan2006
        allocate( eband2(nband,nsp,nqbz) )
        ddq = 0d0
c        if(bzcase()==2) ddq= dq_bzcase2
        if(.not.qbzreg()) ddq= ddq- dq_bzcase2
        write(6,*)' nsp nqbz =', nsp,nqbz
c        write(6,*)' dq_bzcase2 =',dq_bzcase2 ,nsp,nqbz
        do is  = 1,nsp !Readin eband
          ix = 0
          do iqi = 1,nqbz
            ix = ix+1
c          write(6,*)'is ix=',is,ix
c          write(6,"(' --- qbz=',2i5,3f9.4)") iqi,is!,qbz(:,iqi)+ddq
            call readeval(qbz(:,iqi)+ddq, is, eband2(:,is,ix))
            write(6,"(' qbz band=',3f9.4,' ',10f9.4)") qbz(:,iqi)+ddq,eband2(1:5,is,ix)
          enddo
        enddo
        ldim  = nband
        nbmax = nband
        zval  = valn
        nkp   = nqibz
        allocate( bot(nband*nsp), top(nband*nsp))
        call efrang3(nsp,nqbz,nband,nband,zval,e1,e2,eband2,
     &               bot,top,bot,top,elo,ehi,bandgap)
        deallocate(bot,top)
c
        ifief=201
        open(ifief,file='EFERMI')
c --- ef smering method with esmr
c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
! When esmr is negative, esmr is geven automatically by efsimplef.
c$$$        call efsimplef2a(nspin,wibz,qibz,ginv,
c$$$     i        nband,nqibz
c$$$     i       ,konf,z,nl,natom,iclass,nclass
c$$$     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
c$$$     i        qbz,nqbz !index_qbz, n_index_qbz,
c$$$     o       ,efnew)
c$$$        ef = efnew
c$$$c ... check total ele number -------
c$$$        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !require init_readeigen
c$$$        write(6,*)' ef by smearing   =',ef
c$$$        write(6,*)' esmr             =',esmr
c$$$        write(6,*)' valn             =',valn
c$$$        write(6,*)' ntot             =',ntot

c ... read eband
        if(qbzreg() .and. bzcase==1 ) then
          allocate( eband(nband,nspin,nqibz))
          do is  = 1,nspin !Readin eband
            do iqi = 1,nqibz
              call readeval(qibz(:,iqi),is,  eband(:,is,iqi))
              write(6,"(' -------- qibz=',2i5,3f9.4)")iqi,is,qibz(:,iqi)
c           write(6,*)eband(:,is,iqi)
            enddo
          enddo
        endif
      elseif(imode==2.or.imode==3.or.imode==4) then
        ifief=201
        if(imode==2) open(ifief,file='EFERMI.check')
        if(imode==3) open(ifief,file='EFERMI.QP')
        if(imode==4) open(ifief,file='EFERMI.QPz=1')
        iftote=202
        open(iftote,file='TOTE.UP')
        read(iftote,*) nq,ntq, efxx
Cstop2rx 2013.08.09 kino        if(nq/= nqibz) stop ' heftet: nq TOTE/= nqibz'
        if(nq/= nqibz) call rx( ' heftet: nq TOTE/= nqibz')
        nband = ntq
        allocate( eband(nband,nspin,nqibz), qz(3,nqibz) )
        call dcopy (3*nqibz,qibz,1,qz,1)
        do  is  = 1,nspin !Readin eband
          if(is==2) then
            open(iftote,file='TOTE.DN')
            read(iftote,*) nq2,ntq2, efxx2
Cstop2rx 2013.08.09 kino            if(nq2/=nq)   stop 'nq dirrerent TOTE.UP TOTE.DN'
            if(nq2/=nq)   call rx( 'nq dirrerent TOTE.UP TOTE.DN')
Cstop2rx 2013.08.09 kino            if(ntq2/=ntq) stop 'ntq dirrerent TOTE.UP TOTE.DN'
            if(ntq2/=ntq) call rx( 'ntq dirrerent TOTE.UP TOTE.DN')
Cstop2rx 2013.08.09 kino            if(efxx/=efxx)stop 'efxx dirrerent TOTE.UP TOTE.DN'
            if(efxx/=efxx)call rx( 'efxx dirrerent TOTE.UP TOTE.DN')
          endif
          do  iqi = 1,nqibz
            do   it = 1,nband
              read(iftote,*) qq(1:3), itx, iqx, elda, eqp01, eqp02, zfac
Cstop2rx 2013.08.09 kino              if(it /= itx ) stop ' heftet: it /=itx '
              if(it /= itx ) call rx( ' heftet: it /=itx ')
              if(sum(abs(qq(1:3)-qz(1:3,iqi)))>1d-6 ) then
                write(6,*) 'iqi it=',iqi,it
                write(6,*) qq(1:3)
                write(6,*) qz(1:3,iqi)
Cstop2rx 2013.08.09 kino                stop ' heftet: q in TOTE /=qz (iqi) '
                call rx( ' heftet: q in TOTE /=qz (iqi) ')
              endif
              if(imode==2) eband(it,is,iqi) = elda /rydberg() + efxx
              if(imode==3) eband(it,is,iqi) = eqp01/rydberg() + efxx
              if(imode==4) eband(it,is,iqi) = eqp02/rydberg() + efxx
            enddo
          enddo
          close(iftote)
        enddo

caaaa jan2006 bug fix: this bug was between 032f5 to 032f8
        ldim  = nband
        nbmax = nband
        zval  = valn
        nsp   = nspin
        nkp =nqibz
caaaa
        allocate(bot(ldim*nsp),top(ldim*nsp))
        call efrang3(nsp,nkp,nbmax,ldim,zval,e1,e2,eband,
     .            bot,top,bot,top,elo,ehi,bandgap)
        deallocate(bot,top)
      endif
ccccccccccccccccccccccccccccccccccc
c        do  is  = 1,nspin !Readin eband
c        do  iqi = 1,nqibz
c        do   it = 1,nband
c          write(6,'("is q ib eband =",3i5,d13.6)')is,iqi,it,eband(it,is,iqi)
c        enddo
c        enddo
c        enddo
c        stop 'test end xxxxxxxxxx'
cccccccccccccccccccccccccccccccccccc
      write(6,"(' end of efrang2: e1 e2 nteti=',2d15.7,i6)") e1,e2,nteti
      write(6,"(' elo ehi = ',2d15.7,i6)") elo,ehi
      write(6,*)'e1=',e1
      write(6,*)'e2=',e2
cccccccccccccccccccccccccccccccccccccccccccc
      if(e1 /= e2) then
        if(.not. ( qbzreg().and.bzcase==1 ) ) then
          write(6,*)"Fermi energy finder for off-regular mesh"//
     &          "(Chi_RegQbz off)or bzcase==2 is not implimented yet..."
Cstop2rx 2013.08.09 kino          stop 
          call rx( ''//
     &     "Ef finder for off-regular mesh. NotYet implimented")
        endif
      endif
cccccccccccccccccccccccccccccccccccccccccccc

      if(e1 == e2) then ! --- Find band gap  e1=e2 is at the middle of gap---
        efermi = e1
        write(6,*)'heftet: only filled or empty bands encountered.'
        write(ifief,"(2d23.15,a)") efermi,bandgap,
     &  ' ! This efermi bandgap are obtained by heftet: '
        write(ifief,"(a)")'heftet: only filled or empty bands encountered.'
      else ! --- Find the Fermi energy by the Tetrahedron method ---
        write(6,*)' goto bzints2 for dos'
        allocate(dos(nptdos))
        volwgt = (3.d0 - nsp) / ntetf ! ntetf was =6*n1*n2*n3
        write (6,"(1x,'heftet : ---- Tetrahedron Integration ---- ')")
        do itt = 1, 10
          call bzints2x(volwgt,eband,dum111,nkp,ldim,nbmax,
     .              nsp,e1,e2,dos,nptdos,efermi,1,nteti,idteti)
          call fermi(zval,dos,nptdos,e1,e2,efermi,e1,e2,dosef)
c          write(*,"(i3,7x,5(d13.6,1x))")itt, efermi, e1, e2, e2-e1, dosef
          if(e2 - e1 < tol) goto 4
        enddo
        write(*,*) 'heftet:bug in fermi level finder or tol too small'
Cstop2rx 2013.08.09 kino        stop       'heftet:bug in fermi level finder or tol too small'
        call rx( 'heftet:bug in fermi level finder or tol too small')
c        call
c     &  fexit(-1,'heftet:bug in fermi level finder or tol too small')
    4   continue
        deallocate(dos)
        bandgap=0d0
        write(ifief,"(2d23.15,a)") efermi,bandgap," ! This efermi bandgap are obtained by heftet"
      endif
      close(201)

c- eband check
      if(imode/=1) then
        do  is  = 1,nspin
          do  iqi = 1,nqibz
            if( eband(nband,is,iqi) <efermi) then
              write(6,*)' heftet: WARNING! eband(maxband) is less than efermi:
     &    you might have to enlarge nband to get efermi by tetrahedron'
              goto 666
            endif
          enddo
        enddo
        write(6,*) ' check OK! eband(nband,:,:) is greater than efermi.'
 666    continue
      endif
c
      if(imode==1) then
        write(6,"(' Tet EFERMI gap = ',2f24.15)") efermi,bandgap
Cstop2rx 2013.08.09 kino        stop ' OK! heftet mode=1 EFERMI generated '
        call rx0( ' OK! heftet mode=1 EFERMI generated ')
      elseif(imode==2) then
        write(6,"(' Tet EFERMI.check gap= ',2f24.15)") efermi,bandgap
Cstop2rx 2013.08.09 kino        stop ' OK! heftet mode=2 EFERMI.check generated '
        call rx0( ' OK! heftet mode=2 EFERMI.check generated ')
      elseif(imode==3) then
        write(6,"(' Tet EFERMI.QP gap   = ',2f24.15)") efermi,bandgap
Cstop2rx 2013.08.09 kino        stop ' OK! heftet mode=3 EFERMI.QP generated '
        call rx0( ' OK! heftet mode=3 EFERMI.QP generated ')
      elseif(imode==4) then
        write(6,"(' Tet EFERMI.QPz=1 gap= ',2f24.15)") efermi,bandgap
Cstop2rx 2013.08.09 kino        stop ' OK! heftet mode=4 EFERMI.QPz=1 generated '
        call rx0( ' OK! heftet mode=4 EFERMI.QPz=1 generated ')
      endif
      end
