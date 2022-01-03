      program hsfp0_sc
!> Calculates the self-energy \Sigma in GW approximation,  checked 2020jul
!!  including Off-diagonal components.
!!  (hsfp0.F is for diagonal part only).
!! ----------------------------------------
!!    SEx(q,itp,itpp) = <psi(q,itp) |SEx| psi(q,itpp)>
!!    SEc(q,itp,itpp) = <psi(q,itp) |SEc| psi(q,itpp)>
!!    Here SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
!!
!! ----------------------------------------
!! See papers;
!! [1]T. Kotani and M. van Schilfgaarde, Quasiparticle self-consistent GW method: 
!!     A basis for the independent-particle approximation, Phys. Rev. B, vol. 76, no. 16, 
!!     p. 165106[24pages], Oct. 2007.
!! [2]T. Kotani, Quasiparticle Self-Consistent GW Method Based on the Augmented Plane-Wave 
!!     and Muffin-Tin Orbital Method, J. Phys. Soc. Jpn., vol. 83, no. 9, p. 094711 [11 Pages], Sep. 2014.
!!
!! EIBZ symmetrization;
!! See [3] C. Friedrich, S. Bl?gel, and A. Schindlmayr, 
!!   Efficient implementation of the GW approximation within the all-electron FLAPW method, 
!!   Physical Review B, vol. 81, no. 12, Mar. 2010.
!!
!! Usage: This routine is called from a script for QSGW, ecalj/fpgw/exec/gwsc, 
!! which calls as "echo 2|../exec/hsfp0_sc >lsc" when mode=2.
!! In the gwsc script, we call hsfp0_sc three times with mode=1,2,3.    
!!   mode= 1: exchange    mode SEx, the exchange part of the self-energy
!!   mode= 2: correlation mode SEc, the correlated part of the self-energy
!!   mode= 3: core exchange mode SEXcore
!! (unused now) mode= 4: plot spectrum function ---See manual ---> this is performed by echo 4|hsfp0 
!!     
!! iSigMode parameter, which determines approximation for self-energy, is given by GWinput file as iSigMode.
!!     iSigMode==0 SE_nn'(ef)+image integr:delta_nn'(SE_nn(e_n)-SE_nn(ef))
!!     iSigMode==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
!!       xxx not support this mode now ... iSigMode==2 SE_nn'((e_n+e_n')/2)
!! --> iSigMode==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2  <--- this is used as default.
!!     iSigMode==5 delta_nn' SE_nn(e_n)
!!     Output file contain hermitean part of SE for energies to be real
!!    (for example, hermitean conjunction of SE_nn'(e_n) means SE_n'n(e_n')^* )
!!     
!!     History: We learned so much from GW-LMTO-ASA codeds developed by F.Aryasetiawan.
!!      
!! Module prograing style
!!   For writing codes, we recommend our rule of module programing. Minimum are
!!   1)Classification of variables are in modules
!!      We define variables in modules. For example, we use natom (number of atoms in the cell)
!!      as well as pos(3,natom) in the module m_genallcf_v3.
!!      Anywhere in the prograom, we have to read natom and pos contained in m_genallcf_v3.
!!      Not make copies of them. Not declear pos(3,natom) except very simple subroutines.
!!   2) Module variables are 'protected,public' or 'private'. Never use unprotected public variables.
!!      With 'private' as defaults, it is better to show subrouitnes as public. (default as private).
!!   3) subroutines which are not in the modules should be very simple mathematical ones.
!!   4) The intent should be written like this
!!       subroutine foobar (nctot,ncc,nt0,ntp0,  iclass, phase, 
!!       implicit none
!!      i                   icore,ncore,nl,nnc,  
!!      o                   zpsi2b)
!!       intent(in)::       nctot,ncc,nt0,ntp0,  iclass, phase, 
!!      i                   icore,ncore,nl,nnc  
!!       intent(out)::      zpsi2b
!!           ...
!!!! ----------------------------------------
      use m_readefermi,only: READEFERMI,ef_read=>ef
      use m_readqg,only: READQG0,READNGMX2, ngpmx,ngcmx
      use m_readeigen,only: INIT_READEIGEN,INIT_READEIGEN2,READEVAL,LOWESTEVAL
      use m_READ_BZDATA,only: READ_BZDATA,
     &     nqbz,nqibz,nqbzw,nteti,ntetf
     &     ,n1,n2,n3,qbas,ginv,qbz,wbz,qibz,wibz,qbzw,idtetf,ib1bz,idteti
     &     ,nstar,irk,nstbz,ngrp2=>ngrp 
      use m_genallcf_v3,only: GENALLCF_V3,Setesmr,
     &     nclass,natom,nspin,nl,nn, ngrp,
     &     nlmto,nlnmx, nctot,niw, !nw_input=>nw,
     &     alat, deltaw,symgrp,clabl,iclass,esmr_in=>esmr, !diw,dw, delta,
     &     invg, il,in,im,nlnm, 
     &     plat, pos,z,ecore, symgg, konf,nlnx, iantiferro
c      use m_keyvalue,only: GETKEYVALUE
!! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
      use m_rdpp,only: Rdpp,    !"call rdpp" generate following data.
     &     nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
!     ! Generate matrix element for "call get_zmelt".
      use m_zmel,only:          ! folloiwng data set are stored in this module in the main routin, 
     &      Mptauof_zmel ! and used when call get_zmelt, get_zmelt2. Ppbafp_v2_zmel,
c     &                          !ppbir,miat,tiat,shtvg,
      use m_itq,only:
     &     Setitq_hsfp0sc,!Setnband,nband
     &     itq,nbandmx,
     &     ntq
!     ! antiferro condition. only laf is used, after 'call anfcond()'
      use m_anf,only: Anfcond,
     &     laf
!     ! subroutine only
      use m_sxcfsc,only: Sxcf_scz
!     ! MPI
      use m_mpi,only: 
     &     MPI__Initialize,MPI__real8send,MPI__real8recv,MPI__send_iv,MPI__recv_iv,MPI__sxcf_rankdivider,
     &     MPI__Finalize,MPI__root,MPI__Broadcast,MPI__rank,MPI__size,MPI__allreducesum,
     &     MPI__consoleout,
     &     MPI__barrier, MPI__reduceSum
      use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      use m_readfreq_r,only: Readfreq_r, nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i,nw, freq_r
      use m_readq0p,only: ReadQ0P, nq0i, wqt,q0i
      use m_selectqp,only:Getqpoint,q,nq
      use m_eibzhs,only: Seteibzhs,  eibzsym,irkip_all,nrkip_all,tiii
      use m_readgwinput,only: ReadGwinputKeys,SetIsigMode, ebmx_sig,nbmx_sig,ISigMode
      use m_readepswklm,only: Readepswklm!,wklm,lxklm
      implicit none
!! -----------------------------------------------------------------------------------
!!     real(8),parameter :: ua  = 1d0 ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
c------------------------------------
!!!   test switches to calculate the self-energy based on an another separation of \Sigma.
!!!     \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
!!!     I found COH term has inevitably poor accuracy.
c
!\Sigma_{img axis} + \Sigma_{pole} for mode 2
c     & cohtest= .false.         ! \Sigma_{coh}. mode swich is not required.
c     &  , tetra  = .false.  ! test switch for tetrahedron method test.
c     ! tetra=T is only effective for exchange=T case.
c     ! Tetrahedron mehod for correlation is a bit
! difficult and I gave up for a while.
! If you want to calculate with tetra=T for exchange, you
! have to uncomment tetra related part in
! sxcf.f, and a part calling sxcf in this routine. Note wtet wtetef!
! They sometimes cause array destruction if you run tetra=T without comment them.
c
      integer::
     & ixc, noccxv, maxocc, noccx, iaf, ip, is, nspinmx, i, ifoutsex, ix, ifoutsec,
     & ifsec(2), ifxc(2),ifsex(2), ifsex2(2),ifsec2(2),  ifsecomg(2),
     & nq0ix, incwfin, ngpn1,ngcn1,
     &  iqxend,iqxini, timevalues(8) , ret,dest ,irank
      real(8) :: voltot,valn,efnew, rydberg,hartree,qreal(3), nocctotg2,tripl !rs,alpha,ntot,
      real(8) ::  wgtq0p,quu(3), eftrue,esmref,esmr,ef !,ecorem
      character(128) :: ixcc
      logical :: legas, exonly, iprintx, eibz4sig, selectqp=.false.,diagonly=.false.
      logical :: exchange, hermitianW, screen = .false. 
      integer,allocatable:: irkip(:,:,:,:), nrkip(:,:,:,:)
      real(8),allocatable:: vxcfp(:,:,:), eqt(:), eq(:), eqx(:,:,:),eqx0(:,:,:)
      complex(8),allocatable,target:: zsecall(:,:,:,:)
      complex(8),pointer::zsec(:,:,:)
c---------------------------------------
c      if(cohtest) then          !currently not used (may need fixing if necessary)
c        screen = .true.
c        ixc = 2                 ! ; nz=0
c        if671=ifile_hanlde()
c        open(if671,file='COH')
      call MPI__Initialize() 
      call date_and_time(values=timevalues)
      write(6,"('mpirank=',i5,' YYYY.MM.DD.HH.MM.msec=',9i4)")mpi__rank,timevalues(1:3),timevalues(5:8)
      hartree=2d0*rydberg()
      hermitianW=.true.
      if(MPI__root) then
        write(6,*) ' --- Choose modes below ------------'
        write(6,*) '  Sx(1) Sc(2) ScoreX(3) '
        write(6,*) '  [option --- (+ QPNT.{number} ?)] '
        write(6,*) '  Add 1000, eg, 1001 is diagonal only mode for one-shot Z=1'
        write(6,*) ' --- Put number above ! ------------'
        read(5,*) ixc 
        write(6,*) ' ixc=', ixc !computational mode index
      endif
      call MPI__Broadcast(ixc)
      if(ixc>1000) then         !selected QP
        ixc=mod(ixc,1000)
        selectqp  =.true.
        diagonly  =.true.
        hermitianW=.false.
        write(6,*) "--- Diagonal-only mode. jobsw=5; see description at the top of sxcf_fal2.sc.F."
        write(6,*) "--- This is the same as one-shot calculaiton with iSigMode5 in GWinput."
      endif
      write(ixcc,"('.mode=',i4.4)")ixc
      call MPI__consoleout('hsfp0_sc'//trim(ixcc)) !Open console output stdout.irank.hsfp0_sc.mode 
      call pshprt(60)
!! === Readin GWIN and LMTO, then allocate and set datas. ===
!! See use m_genallcf_v3 at the top of this routine
      if(ixc==3) then; incwfin= -2 !core exchange mode
      else           ; incwfin= -1 !use 7th colmn for core at the end section of GWIN
      endif
      call GENALLCF_V3(incwfin)    ! readin basic data
      call READ_BZDATA() ! Readin BZ data. See gwsrc/rwbzdata.f ===
!!CAUTION!!  WE currently ASSUME iclass(iatom)= iatom (because of historical reason)
      if(nclass /= natom ) call rx('hsfp0_sc: nclass /= natom ')
      if(ngrp/= ngrp2)     call rx('ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      write(6,"(' nqbz nqibz ngrp=',3i12)") nqbz,nqibz,ngrp
      esmr  = esmr_in !read from GWinput
      esmref= esmr
      call ReadGwinputKeys()
      if(diagonly) call SetIsigMode(5)
      call pshprt(30)
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3))) ! primitive cell volume
c$$$!! ef is taken as rs for the empty-sphere test case of legas=T case 
c$$$!! HOMOGENIOUS GAS code. Usually not used. Need fixing if necessary.
c$$$!! Keep this just as a memo.
c$$$      legas = .false.
c$$$      if(.false.) then
c$$$        INQUIRE (FILE = 'LEGAS', EXIST = legas)
c$$$        if(legas) then          !!! test for electron gas case.
c$$$          write(6,*)' find LEGAS. legas =',legas
c$$$          iflegas = 2101
c$$$          open (iflegas,file='LEGAS')
c$$$          read(iflegas,*)rs
c$$$          close(iflegas)
c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
c$$$          qfermi = alpha/rs
c$$$          efx  = qfermi**2
c$$$          valn = efx**1.5d0*voltot/3d0/pi**2
c$$$          write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
c$$$          write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
c$$$          write (6,*)' egas  Exact Fermi energy    Ef  =', efx
c$$$          if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
c$$$        endif
c$$$      endif
!!
      if(ixc==1) then
        exchange = .true.
        write(6,*) ' --- Exchange mode --- '
      elseif(ixc==2) then
        exchange=.false.
        write(6,*) ' --- Correlation mode --- '
      elseif(ixc==3) then
         esmr=0d0
         exchange = .true.
         write(6,*) ' --- CORE Exchange mode --- '
      else
        call rx(' hsfp0_sc: Need input (std input) 1(Sx) 2(Sc) or 3(ScoreX)!')
      endif
      call setesmr(esmr_in=esmr) !set esmr back in genalloc_v3
      call Readhbe()        ! Read dimensions in m_readhbe
      call INIT_READEIGEN() ! initialization for readeigen readchpi readgeig.
      call INIT_READEIGEN2()! initialize m_readeigen
      call Mptauof_zmel(symgg,ngrp) ! Put space-group transformation information to m_zmel
      call Readngmx2() !return ngpmx and ngcmx in m_readqg
      write(6,*)' max number of G for QGpsi and QGcou: ngcmx ngpmx=',ngcmx,ngpmx
      call pshprt(60)
      if(.not.exchange.or.(exchange.and.screen)) then !screen means screened exchange case
         call readfreq_r() !Readin WV.d and freq_r 
      endif
c      call READEFERMI() !! efermi by tetrahedron. 
c      ef = ef_read
!! Determine Fermi energy ef for given valn (legas=T), or Get charge given by z and konf.==
      legas=.false. ! if legas=T, homogenius electron gas test case.
      call efsimplef2ax(legas, esmref, 
     o     valn, ef)
      eftrue = ef
      if(ixc==3) then
         ef = LOWESTEVAL() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
      endif
!!--------------      
      write(6,'(" --- computational conditions --- ")')
      write(6,'("    deltaw  =",f13.6)') deltaw
      write(6,'("    esmr    =",f13.6)') esmr
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
      write(6,*)' ef    =',ef
      write(6,*)' esmr  =',esmr
      write(6,*)' valn  =',valn
      
!! Core-exchange case. We set Ef just below the valence eigenvalue (to pick up only cores).==
! exonly removed
      if(ixc==3) then
c        ef = LOWESTEVAL() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
c        if(wex==0d0) then
c          exonly=.false.
c        else !we need to check wheter exorem trick is really useful or not. jul2020
c           exonly=.true.
c           ecorem = 1000d0 !write(6,*)' exonly=T ecore shift: ecore---> ecore-ecorem'
c        endif
        write(6,"(a)")' CoreEx mode: We change ef as ef=LOWESTEVAL-1d-3, slightly below the bottom of valence.'
        write(6,"(a,f13.5,i5,i5)")' CoreEx mode: ef nspin nctot=',ef,nspin,nctot
        do ix=1,nctot
           write(6,"(i4,x,d13.5,x,d13.5)") ix,(ecore(ix,is),is=1,nspin)
        enddo
      endif
      call Getqpoint(selectqp,nqibz,qibz) !Get q(3,nq) from <QPNT> or qibz
      call anfcond()
      nspinmx = nspin
      if(laf) nspinmx=1         !!! Antiferro case. Only calculate up spin
!! itq.
!!  We calculate <ki|\sigma|kj> for i \in itq (and j \in itq).
!!  During iteration, we use NTQXX file, to keep itq set.
!!  If we already have NTQXX (in your pre gwsc calculion), the NTQXX is read.
      if( mpi__root .and. mpi__rank/=0) call rx('mpi__root .and. mpi__rank/=0')
      do irank = 0,mpi__size-1  ! irank=0 may write NTQXX if it exists. irank>0 is reading mode.
        if(mpi__rank==irank) call Setitq_hsfp0sc(qibz,nqibz,nq,nspin,nbmx_sig,ebmx_sig,eftrue,nspinmx) !read NTQXX
        call MPI__barrier() !mpi barrier is only for root ?
      enddo
!! Read eigenvalues given by lmf for writing out files.
      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
      do is = 1,nspin
        do ip = 1,nq
          eqt= READEVAL(q(1,ip),is)
          eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
          eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- eftrue)
        enddo
      enddo
      deallocate(eqt)
!!-----------------------------------------------------------------------
      call Hswriteinit()        !internal subroutine.  Write initial part of output files
!!      
      if(ixc==3.and.nctot==0) then !make dummy SEXcore
         allocate(zsecall(ntq,ntq,nq,1)) !, coh(ntq,nq) ) kount(nqibz,nq),
         zsecall = 0d0
         do is=1,nspinmx
            zsec => zsecall(:,:,:,1)
            call HsWriteResult() !internal subroutine.
         enddo
         goto 9999
      endif
      
!! Offset Gamma point Q0P (slightly shifted from Gamma)
      call Readq0p()
      call Readepswklm()
!! == irkip:  paralellization is contrlloed by irkip ==
!! We have to distribute non-zero irkip to all ranks (irkip is dependent on rank).
!! When irkip(nqibz,ngrp,nq,nspinmx)/=0, we expect grain-size 
!! on each job of (iqibz,igrp,iq,isp) is almost the same.
!! Our pupose is to calculate the self-energy zsec(itp,itpp,iq).
      call Seteibzhs(nspinmx,nq,q,iprintx=MPI__root)
      allocate(irkip(nspinmx,nqibz,ngrp,nq)) ! local
      allocate(nrkip(nspinmx,nqibz,ngrp,nq)) ! nrkip is weight correspoinding to irkip for a node.
      call MPI__sxcf_rankdivider(irkip,irkip_all,nspinmx,nqibz,ngrp,nq) ! MIZUHO-IR
      nrkip = nrkip_all         ! we don't need to change this for MPI case. 
! ----> We need to distribute non-zero irkip to ranks.
!! Main loop to calculate the the self-energy SEx(ip) or SEc(ip)
      if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
      if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
      allocate(zsecall(ntq,ntq,nq,nspinmx)) !, coh(ntq,nq) ) kount(nqibz,nq),
      zsecall = 0d0
      
      do 2000 is = 1,nspinmx
        zsec => zsecall(:,:,:,is)
!! == ip loop to spedify external q ==
        call sxcf_scz (qip= q, ef= ef,esmr= esmr, isp= is,   
     i   irkip= irkip(is,:,:,:),nrkip= nrkip(is,:,:,:),nq= nq,
     i   exchange= exchange, 
     i   hermitianW= hermitianW,
     i   jobsw= iSigMode, nbandmx= nbandmx(1:nq,is), !nbandmx is input mar2015
     o   zsec= zsec) 
c$$$!!  electron gas bare exchange (exact)
c$$$        if (legas.and.exchange) then
c$$$          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
c$$$          pi         = 4.d0*datan(1.d0)
c$$$          tpia       = 2.d0*pi/alat
c$$$          qfermi= dsqrt(efz)
c$$$          alpha = (9*pi/4d0)**(1d0/3d0)
c$$$          write (6,*)' --- exact electron gas bare exchange --- '
c$$$          write (6,*)' density parameter rs= ', alpha/qfermi
c$$$          write (6,*)' kf= ',qfermi
c$$$          do      ip = 1,nq
c$$$            qreal =  tpia*q(1:3,ip)
c$$$            qm    = dsqrt ( sum(qreal**2) )
c$$$            xsex  = hartree * egex (qm,efz)
c$$$            write (6,*)
c$$$            write (6,"(' True qm-ef Sx=',2f14.6,' q/qf=',f14.6)")
c$$$     &       rydberg()*(qm**2-efz), xsex, qm/qfermi
c$$$            write (6,"(' Num  qm-ef Sx=',2f14.6)") 
c$$$     &       eqx(1,ip,is),        hartree*dreal(zsec(1,1,ip)) !sf 21May02
c$$$            write (6,"(' === diff     =',2f14.6)") 
c$$$     &       rydberg()*(qm**2-efz)-eqx(1,ip,is)
c$$$     &       , xsex - hartree*dreal(zsec(1,1,ip)) !sf 21May02
c$$$            write (661,"(' qm True qm-ef Sx=',3f14.6)") 
c$$$     &       qm,rydberg()*(qm**2-efz), xsex
c$$$            write (662,"(' qm Num  qm-ef Sx=',3f14.6)") 
c$$$     &       qm,eqx(1,ip,is),     hartree*dreal(zsec(1,1,ip)) !sf 21May02
c$$$ccc   write (ifsex(is),6600) qreal(1),qreal(2),qreal(3),xsex
c$$$ccc   write (6,6600) qreal(1),qreal(2),qreal(3),xsex
c$$$ccc   6600   format (' qreal =',3f8.4,'   SEx(q) =',d13.5)
c$$$            write (663,"(2f14.6)") qm/qfermi, qfermi
c$$$          end do
c$$$        endif
        if(eibz4sig())call zsecsym(zsec,ntq,nq,nband,nbandmx,nspinmx, nspin,eibzsym,ngrp,tiii,q,is)
 2000 continue   !end of spin-loop
!!      
      call MPI__reduceSum(root=0, data=zsecall, sizex=ntq*ntq*nq*nspinmx )
      if(MPI__root) then
        do is=1,nspinmx
        zsec => zsecall(:,:,:,is)
        call HsWriteResult() !internal subroutine.
        enddo
      endif
!!      
 9999 continue
      call cputid(0)
      if(ixc==1 ) call rx0( ' OK! hsfp0_sc: Exchange mode')
      if(ixc==2 ) call rx0( ' OK! hsfp0_sc: Correlation mode')
      if(ixc==3 ) call rx0( ' OK! hsfp0_sc: Core-exchange mode')
      stop 
!-----------------------------------------------------------------------------
      contains
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      subroutine Hswriteinit() !only write out files
      implicit none
      write (6,*)' ***'
      call READQG0('QGpsi',qibz(1:3,1), quu,ngpn1)
      call READQG0('QGcou',qibz(1:3,1), quu,ngcn1)
      write (6,6700) nspin,nq,ntq
 6700 format (1x,3i4,'  nspin  nq  ntq')
      write (6,6501) is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,ef,esmr
 6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4
     & ,'  nqbz =',i6,'  nqibz =',i6,'   ef=', f10.4,' Rydberg'
     & ,/,d23.16,' <= deltaw(Hartree)'
     & ,/,d23.16,' <= alat'
     & ,/,d23.16,' <= ef '
     &     ,/,d23.16,' <= esmr')
!! Print LDA exchange-correlation files XCU and XCD
      if(ixc==1) then
        allocate(  vxcfp(ntq,nq,nspin) )
        call rsexx(nspin,itq,q,ntq,nq, ginv, vxcfp) !add ginv july2011
        if(MPI__root) then
           do is = 1,nspinmx
            if(is==1) open(newunit=ifxc(1),file='XCU')!//xt(nz))
            if(is==2) open(newunit=ifxc(2),file='XCD')!//xt(nz))
            write (ifxc(is),*) '==================================='
            write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
            write (ifxc(is),*) '==================================='
            call winfo(ifxc(is),nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
            write (ifxc(is),*)' ***'
            write (ifxc(is),"(a)") ' jband   iq ispin                  qvec eigen-Ef (in eV)     LDA XC (in eV)'
            ifoutsex = ifxc(is)
            write(6,*)
            do ip = 1,nq
            do i  = 1,ntq
               write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), vxcfp(i,ip,is)
               if(eqx(i,ip,is) <1d20.and.vxcfp(i,ip,is)/=0d0) then
!              !takao june2009. See lmf2gw (evl_d=1d20; in Ry.. but eqx is in eV. no problem for inequality).
                  write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f10.4,'  Sxc(LDA)=',f10.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), vxcfp(i,ip,is)
               endif
            enddo
            enddo
            close(ifxc(is))
          enddo                 !     end of spin-loop
        endif                   !MPI__root
        deallocate(vxcfp)
      endif
      end subroutine
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      subroutine HsWriteResult() !only write out files
      implicit none
      integer:: iii,ixx,iyy
      if(exchange) then
         if(is==1.and.ixc==1) then
            open(newunit=ifsex(1),file='SEXU') !//xt(nz))
            open(newunit=ifsex2(1),file='SEX2U',form='unformatted') !out SEX_nn'
         elseif(is==2.and.ixc==1) then
            open(newunit=ifsex(2),file='SEXD') !//xt(nz))
            open(newunit=ifsex2(2),file='SEX2D',form='unformatted') !out SEX_nn'
         elseif(is==1.and.ixc==3) then
            open(newunit=ifsex(1),file='SEXcoreU') !//xt(nz))
            open(newunit=ifsex2(1),file='SEXcore2U',form='unformatted') !out SEXcore_nn'
cccccccccccccccccccccccccc            
c            open(newunit=iii,file='SEXxxx') !out SEX_nn'
cccccccccccccccccccccccccc            
         elseif(is == 2.and.ixc==3) then
            open(newunit=ifsex(2),file='SEXcoreD') !//xt(nz))
            open(newunit=ifsex2(2),file='SEXcore2D',form='unformatted') !out SEXcore_nn'
         endif
         write(ifsex(is),*) '======================================='
         write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
         write(ifsex(is),*) '======================================='
         call winfo(ifsex(is),nspin,nq,ntq,is,nbloch,ngpn1,
     &        ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
         write (ifsex(is),*)' *** '
         write (ifsex(is),"(a)") ' jband   iq ispin                  '//
     &        '           qvec            eigen-Ef (in eV)           exchange (in eV)'
         write(ifsex2(is)) nspin, nq, ntq,nqbz,nqibz, n1,n2,n3
         ifoutsex=ifsex(is)
         write(6,*)
         do ip = 1,nq
            do i  = 1,ntq
               write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &              itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &              hartree*dreal(zsec(i,i,ip)) !sf 21May02
               if( eqx(i,ip,is)<1d20.and.abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                  write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,' eig=',f10.4,'  Sx=',f10.4)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 hartree*dreal(zsec(i,i,ip)) !sf 21May02
               endif
            enddo
            write(ifsex2(is)) is, q(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
ccccccccccccccccccccccccccccc
c            if(is==1) then
c            do ixx=1,ntq
c            do iyy=1,ntq
c               write(iii,"('zsec:',i2, 3f6.3,x,2i4,2x,9(2f8.3,x))") is, q(1:3,ip),ixx,iyy,zsec(ixx,iyy,ip)
c            enddo
c            enddo
c            endif
ccccccccccccccccccccccccccccc            
         enddo
         close(ifsex(is))
         close(ifsex2(is))
      elseif(ixc==2) then
         if(is == 1.and.ixc==2) then
            open(newunit=ifsec(1),file='SECU') !//xt(nz))
            open(newunit=ifsec2(1),file='SEC2U',form='unformatted') !out SEC_nn'
         elseif(is == 2.and.ixc==2) then
            open(newunit=ifsec(2),file='SECD') !//xt(nz))
            open(newunit=ifsec2(2),file='SEC2D',form='unformatted') !out SEC_nn'
         endif
         write(ifsec2(is)) nspin, nq, ntq ,nqbz,nqibz  ,n1,n2,n3
         write(ifsec(is),*) '=========================================='
         write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)") is
         write(ifsec(is),*) '=========================================='
         call winfo(ifsec(is),nspin,nq,ntq,is,nbloch,ngpn1,
     &        ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
         write (ifsec(is),*)' *** '
         write (ifsec(is),"(a)") ' jband   iq ispin                  '//
     &        '           qvec            eigen-Ef (in eV)           '//
     &        'Re(Sc) 3-points (in eV)                        '//
     &        '           In(Sc) 3-points (in eV)                Zfactor(=1)'
         ifoutsec = ifsec(is)
         do ip = 1,nq
            do i  = 1,ntq
               if( eqx(i,ip,is)<1d20.and.abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                  write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f8.4,'  Re(Sc) =',f8.4,'  Img(Sc) =',f8.4 )") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zsec(i,i,ip)), 
     &                 hartree*dimag(zsec(i,i,ip))  
               endif
               write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,d24.16, 3x,d24.16)")
     &              itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &              hartree*dreal(zsec(i,i,ip)), 
     &              hartree*dimag(zsec(i,i,ip))  
            end do
            write(ifsec2(is)) is, q(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
         end do
         close(ifsec(is))
         close(ifsec2(is))
      endif                     !ixc
      end subroutine
      end program hsfp0_sc
      
