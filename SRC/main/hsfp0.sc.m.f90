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
  use m_hamindex,only:   Readhamindex, symgg=>symops
  use m_readeigen,only: INIT_READEIGEN,INIT_READEIGEN2,READEVAL,LOWESTEVAL
  use m_READ_BZDATA,only: READ_BZDATA, &
       nqbz,nqibz,nqbzw,nteti,ntetf &
       ,n1,n2,n3,ginv,qbz,wbz,qibz,wibz,qbzw,idtetf,ib1bz,idteti &
       ,nstar,irk,nstbz, nq0i, wqt=>wt, q0i
  use m_genallcf_v3,only: GENALLCF_V3,Setesmr, nclass,natom,nspin,nl,nn, &
       nlmto,nlnmx, nctot,niw, alat, deltaw,clabl,iclass,esmr_in=>esmr, il,in,im,nlnm, &
       plat, pos,z,ecore, konf,nlnx
  !! Base data to generate matrix elements zmel*. Used in "call get_zmelt".
  use m_rdpp,only: Rdpp, nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
  !     ! Generate matrix element for "call get_zmelt".
  use m_zmel,only: Mptauof_zmel
  use m_itq,only: setitq_hsfp0sc,itq,nbandmx, ntq
  use m_anf,only: Anfcond, laf
  use m_sxcf_count,only: sxcf_scz_count
  use m_sxcf_main,only: sxcf_scz_main,zsecall
  use m_mpi,only: &
       MPI__Initialize,MPI__real8send,MPI__real8recv,MPI__send_iv,MPI__recv_iv, &
       MPI__Finalize,MPI__root,MPI__Broadcast,MPI__rank,MPI__size,MPI__allreducesum, &
       MPI__consoleout, &
       MPI__barrier, MPI__reduceSum
  use m_readhbe,only: Readhbe, nband !nprecb,mrecb,mrece,nlmtot,nqbzt,,mrecg
  use m_readfreq_r,only: Readfreq_r, nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i,nw, freq_r
!  use m_selectqp,only:Getqpoint,qvec,nq
!  use m_eibzhs,only: Seteibzhs, eibzsym,tiii
  use m_readgwinput,only: ReadGwinputKeys,SetIsigMode, ebmx_sig,nbmx_sig,ISigMode
  use m_hamindex,only:ngrp
  use m_lgunit,only:m_lgunit_init,stdo
  implicit none
  ! real(8),parameter :: ua  = 1d0 ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
!!!   test switches to calculate the self-energy based on an another separation of \Sigma.
!!!     \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
!!!     I found COH term has inevitably poor accuracy.

  ! Sigma_{img axis} + \Sigma_{pole} for mode 2
  !     & cohtest= .false.         ! \Sigma_{coh}. mode swich is not required.
  !     &  , tetra  = .false.  ! test switch for tetrahedron method test.
  !     ! tetra=T is only effective for exchange=T case.
  !     ! Tetrahedron mehod for correlation is a bit
  ! difficult and I gave up for a while.
  ! If you want to calculate with tetra=T for exchange, you
  ! have to uncomment tetra related part in
  ! sxcf.f, and a part calling sxcf in this routine. Note wtet wtetef!
  ! They sometimes cause array destruction if you run tetra=T without comment them.

  integer:: &
       ixc,  maxocc, iaf, ip, is, nspinmx, i, ifoutsex, ix, ifoutsec, &
       ifsec(2), ifxc(2),ifsex(2), ifsex2(2),ifsec2(2),  ifsecomg(2), &
       nq0ix, incwfin, ngpn1,ngcn1, &
       iqxend,iqxini, timevalues(8) , ret,dest ,irank,isp,nq
  real(8) :: voltot,valn,efnew, rydberg,hartree,qreal(3), tripl !rs,alpha,ntot,
  real(8) ::  wgtq0p,quu(3), eftrue,esmref,esmr,ef !,ecorem
  character(128) :: ixcc
  logical :: legas, exonly, iprintx,diagonly=.false. !, selectqp=.false.
  logical :: exchange, hermitianW=.true.!, screen = .false.
  integer,allocatable:: irkip(:,:,:,:), nrkip(:,:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), eqt(:), eq(:), eqx(:,:,:),eqx0(:,:,:),qvec(:,:)
  complex(8),pointer::zsec(:,:,:)
  call MPI__Initialize()
  call M_lgunit_init()
  call date_and_time(values=timevalues)
  write(6,"('mpirank=',i5,' YYYY.MM.DD.HH.MM.msec=',9i4)")mpi__rank,timevalues(1:3),timevalues(5:8)
  if(MPI__root) then
     write(6,*) ' --- Choose modes below ------------'
     write(6,*) '  Sx(1) Sc(2) ScoreX(3) '
     write(6,*) '  [option --- (+ QPNT.{number} ?)] '
     write(6,*) ' --- Put number above ! ------------'
     read(5,*) ixc
     write(6,*) ' ixc=', ixc !computational mode index
  endif
  call MPI__Broadcast(ixc)
  write(ixcc,"('.mode=',i4.4)")ixc
  call MPI__consoleout('hsfp0_sc'//trim(ixcc)) !Open console output stdout.irank.hsfp0_sc.mode
  call pshpr(60)
  if(ixc==3) then; incwfin= -2 !core exchange mode
  else           ; incwfin= -1 !use 7th colmn for core at the end section of GWIN
  endif
  call GENALLCF_V3(incwfin)    ! readin basic data
  call READ_BZDATA() ! Readin BZ data. See gwsrc/rwbzdata.f ===
  if(nclass /= natom ) call rx('hsfp0_sc: sanitiy check nclass /= natom ')! CAUTION. ASSUME iclass(iatom)= iatom (because of historical reason)
  write(6,"(' nqbz nqibz ngrp=',3i12)") nqbz,nqibz,ngrp
  esmr  = esmr_in !read from GWinput
  esmref= esmr
  call ReadGwinputKeys()
  call pshpr(30)
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3))) ! primitive cell volume
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
  call Readhamindex()
  call INIT_READEIGEN() ! initialization for readeigen readchpi readgeig.
  call INIT_READEIGEN2()! initialize m_readeigen
  call Mptauof_zmel(symgg,ngrp) ! Put space-group transformation information to m_zmel
  call Readngmx2() !return ngpmx and ngcmx in m_readqg
  write(stdo,"(*(g0))")' max number of G for QGpsi and QGcou: ngcmx ngpmx=',ngcmx,ngpmx
  call pshpr(60)
  if(.NOT.exchange) call readfreq_r()  !Readin WV.d and freq_r
  legas=.false. ! if legas=T, homogenius electron gas test case.
  call efsimplef2ax(legas,esmref, valn,ef)!Get num of val electron valn and Fermi energy ef. legas=T give ef for given valn.
  eftrue = ef
  if(ixc==3)ef = LOWESTEVAL() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
  write(6,'(" --- computational conditions --- ")')
  write(6,'("    deltaw  =",f13.6)') deltaw
  write(6,'("    esmr    =",f13.6)') esmr
  write(6,'("    alat voltot =",2f13.6)') alat, voltot
  write(6,*)' ef    =',ef
  write(6,*)' esmr  =',esmr
  write(6,*)' valn  =',valn
  nspinmx = nspin
  call anfcond()
  if(laf) nspinmx=1         !!! Antiferro case. Only calculate up spin
  !!  We calculate <ki|\sigma|kj> for i \in itq (and j \in itq).
  !!  During iteration, we use NTQXX file, to keep itq set.
  !!  If we already have NTQXX (in your pre gwsc calculion), the NTQXX is read.
  if( mpi__root .AND. mpi__rank/=0) call rx('mpi__root .AND. mpi__rank/=0')
  do irank = 0,mpi__size-1  ! irank=0 may write NTQXX if it exists. irank>0 is reading mode.
     if(mpi__rank==irank) call Setitq_hsfp0sc(qibz,nqibz,nspin,nbmx_sig,ebmx_sig,eftrue,nspinmx) !read NTQXX
  enddo
  call MPI__barrier() !mpi barrier is only for root ?
  SchedulingSelfEnergyCalculation: Block
    if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
    if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
    if ( .NOT. iSigMode==3) call rx('sxcf_scz: only for jobsw=3') !nbandmx is input mar2015 jobsw= iSigMode=3 only
    call sxcf_scz_count(ef,esmr,exchange,nbandmx,ixc,nspinmx) !initializatition quick
  endblock SchedulingSelfEnergyCalculation
  
  if(ixc==3) then ! Core-exchange mode. We set Ef just below the valence eigenvalue (to pick up only cores)
     write(6,"(a)")'CoreEx mode: We change ef as ef=LOWESTEVAL-1d-3, slightly below the bottom of valence.'
     write(6,"(a,f13.5,i5,i5)")' CoreEx mode: ef nspin nctot=',ef,nspin,nctot
     do ix=1,nctot
        write(6,"(i4,x,d13.5,x,d13.5)") ix,(ecore(ix,is),is=1,nspin)
     enddo
  endif
  !! Read eigenvalues given by lmf for writing out files.
  allocate(eqx(ntq,nqibz,nspin),eqt(nband),qvec(3,nqibz)) !,eqx0(ntq,nq,nspin)
  qvec=qibz
  nq  =nqibz
  do is = 1,nspin
     do ip = 1,nqibz
        eqt= READEVAL(qibz(1,ip),is)
        eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- eftrue)
     enddo
  enddo
  deallocate(eqt)
  hartree=2d0*rydberg()
  call Hswriteinit()        !internal subroutine.  Write initial part of output files
  !! Currently, irkip, involved in eibz (extended BZ) mode is removed.
  !! == irkip:  parallelization is controled by irkip ==
  !! We have to distribute non-zero irkip to all ranks (irkip is dependent on rank).
  !! When irkip(nqibz,ngrp,nq,nspinmx)/=0, we expect grain-size
  !! on each job of (iqibz,igrp,iq,isp) is almost the same.
  !! Our pupose is to calculate the self-energy zsec(itp,itpp,iq).
  ! Remove eibzmode symmetrizer 2023Jan22
  !call Seteibzhs(nspinmx,nq,qvec,iprintx=MPI__root)
  
  Main4SelfEnergy: Block !time-consuming part Need highly paralellized
    call sxcf_scz_main(qibz,ef,esmr,nqibz,exchange,nbandmx,ixc,nspinmx) !main part of job
  EndBlock Main4SelfEnergy
! Remove eibzmode symmetrizer 2023Jan22 (extended irreducibel BZ mode)
!  SymmetrizeZsec :Block
!    logical:: eibz4sig
!    if(eibz4sig())then
!       do isp=1,nspinmx
!          call zsecsym(zsecall(:,:,:,isp),ntq,nq,nband,nbandmx,nspinmx,nspin,eibzsym,ngrp,tiii,qvec,isp)
!       enddo
!    endif
!  EndBlock SymmetrizeZsec
  call MPI__reduceSum(root=0, data=zsecall, sizex=ntq*ntq*nqibz*nspinmx )
  if(MPI__root) then
     do is=1,nspinmx
        zsec => zsecall(:,:,:,is)
        call HsWriteResult() !internal subroutine.
     enddo
  endif
  call cputid(0)
  if(ixc==1 ) call rx0( ' OK! hsfp0_sc: Exchange mode')
  if(ixc==2 ) call rx0( ' OK! hsfp0_sc: Correlation mode')
  if(ixc==3 ) call rx0( ' OK! hsfp0_sc: Core-exchange mode')
  stop
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
contains 
  subroutine Hswriteinit() !contained in hsfp0_sc. Only write out files, no side effect
    implicit none
    write (6,*)' ***'
    call READQG0('QGpsi',qibz(1:3,1), quu,ngpn1)
    call READQG0('QGcou',qibz(1:3,1), quu,ngcn1)
    write (6,6700) nspin,nq,ntq
6700 format (1x,3i4,'  nspin  nq  ntq')
    write (6,6501) is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,ef,esmr
6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4 &
         ,'  nqbz =',i6,'  nqibz =',i6,'   ef=', f10.4,' Rydberg' &
         ,/,d23.16,' <= deltaw(Hartree)',/,d23.16,' <= alat' &
         ,/,d23.16,' <= ef ',/,d23.16,' <= esmr')
    !! Print LDA exchange-correlation files XCU and XCD
    if(ixc==1) then
       allocate(  vxcfp(ntq,nq,nspin) )
       call rsexx(nspin,itq,qvec,ntq,nq, ginv, vxcfp) !add ginv july2011
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
                   write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                        itq(i),ip,is, qvec(1:3,ip), eqx(i,ip,is), vxcfp(i,ip,is)
                   if(eqx(i,ip,is) <1d20 .AND. vxcfp(i,ip,is)/=0d0) then
                      !    !takao june2009. See lmf2gw (evl_d=1d20; in Ry.. but eqx is in eV. no problem for inequality).
                      write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f10.4,'  Sxc(LDA)=',f10.4)") &
                           itq(i),ip,is, qvec(1:3,ip), eqx(i,ip,is), vxcfp(i,ip,is)
                   endif
                enddo
             enddo
             close(ifxc(is))
          enddo                 !     end of spin-loop
       endif                   !MPI__root
       deallocate(vxcfp)
    endif
  end subroutine Hswriteinit
  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine HsWriteResult() !contained in hsfp0_sc. Only write out files, no side effect
    implicit none
    integer:: iii,ixx,iyy
    if(exchange) then
       if(is==1 .AND. ixc==1) then
          open(newunit=ifsex(1),file='SEXU') !//xt(nz))
          open(newunit=ifsex2(1),file='SEX2U',form='unformatted') !out SEX_nn'
       elseif(is==2 .AND. ixc==1) then
          open(newunit=ifsex(2),file='SEXD') !//xt(nz))
          open(newunit=ifsex2(2),file='SEX2D',form='unformatted') !out SEX_nn'
       elseif(is==1 .AND. ixc==3) then
          open(newunit=ifsex(1),file='SEXcoreU') !//xt(nz))
          open(newunit=ifsex2(1),file='SEXcore2U',form='unformatted') !out SEXcore_nn'
       elseif(is == 2 .AND. ixc==3) then
          open(newunit=ifsex(2),file='SEXcoreD') !//xt(nz))
          open(newunit=ifsex2(2),file='SEXcore2D',form='unformatted') !out SEXcore_nn'
       endif
       write(ifsex(is),*) '======================================='
       write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
       write(ifsex(is),*) '======================================='
       call winfo(ifsex(is),nspin,nq,ntq,is,nbloch,ngpn1, &
            ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
       write (ifsex(is),*)' *** '
       write (ifsex(is),"(a)") ' jband   iq ispin                  '// &
            '           qvec            eigen-Ef (in eV)           exchange (in eV)'
       write(ifsex2(is)) nspin, nq, ntq,nqbz,nqibz, n1,n2,n3
       ifoutsex=ifsex(is)
       write(6,*)
       do ip = 1,nq
          do i  = 1,ntq
             write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                  itq(i),ip,is, qvec(1:3,ip), eqx(i,ip,is), &
                  hartree*dreal(zsec(i,i,ip)) !sf 21May02
             if( eqx(i,ip,is)<1d20 .AND. abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,' eig=',f10.4,'  Sx=',f10.4)") &
                     itq(i),ip,is, qvec(1:3,ip), eqx(i,ip,is), &
                     hartree*dreal(zsec(i,i,ip)) !sf 21May02
             endif
          enddo
          write(ifsex2(is)) is, qvec(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
       enddo
       close(ifsex(is))
       close(ifsex2(is))
    elseif(ixc==2) then
       if(is == 1 .AND. ixc==2) then
          open(newunit=ifsec(1),file='SECU') !//xt(nz))
          open(newunit=ifsec2(1),file='SEC2U',form='unformatted') !out SEC_nn'
       elseif(is == 2 .AND. ixc==2) then
          open(newunit=ifsec(2),file='SECD') !//xt(nz))
          open(newunit=ifsec2(2),file='SEC2D',form='unformatted') !out SEC_nn'
       endif
       write(ifsec2(is)) nspin, nq, ntq ,nqbz,nqibz  ,n1,n2,n3
       write(ifsec(is),*) '=========================================='
       write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)") is
       write(ifsec(is),*) '=========================================='
       call winfo(ifsec(is),nspin,nq,ntq,is,nbloch,ngpn1, &
            ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
       write (ifsec(is),*)' *** '
       write (ifsec(is),"(a)") ' jband   iq ispin                  '// &
            '           qvec            eigen-Ef (in eV)           '// &
            'Re(Sc) 3-points (in eV)                        '// &
            '           In(Sc) 3-points (in eV)                Zfactor(=1)'
       ifoutsec = ifsec(is)
       do ip = 1,nq
          do i  = 1,ntq
             if( eqx(i,ip,is)<1d20 .AND. abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f8.4,'  Re(Sc) =',f8.4,'  Img(Sc) =',f8.4 )") &
                     itq(i),ip,is, qvec(1:3,ip), eqx(i,ip,is), &
                     hartree*dreal(zsec(i,i,ip)), &
                     hartree*dimag(zsec(i,i,ip))
             endif
             write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,d24.16, 3x,d24.16)") &
                  itq(i),ip,is, qvec(1:3,ip), eqx(i,ip,is), &
                  hartree*dreal(zsec(i,i,ip)), &
                  hartree*dimag(zsec(i,i,ip))
          end do
          write(ifsec2(is)) is, qvec(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
       end do
       close(ifsec(is))
       close(ifsec2(is))
    endif                     !ixc
  end subroutine HsWriteResult
end program hsfp0_sc

subroutine rsexx (nspin, itq, q, ntq,nq,ginv, vxco)
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  dimension vxco(ntq,nq,nspin),q(3,nq),itq(ntq) !itq is not dependent on q, right?
  real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
  logical ::nocore,lfind
  real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3)
  integer:: ikpx=999999
  write(6,*)' OPEN VXCFP '
  open(newunit=ifvxcfp,file='VXCFP',form='unformatted')
  read(ifvxcfp) ldim,nqbz
  write(6,*)' rsexx ldim,nqbz',ldim,nqbz
  allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
  do ikp = 1,nqbz
     read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
     write(6,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
  enddo
  close(ifvxcfp)
  do iq=1,nq
     do ikp=1,nqbz
        lfind=.false.
        if(sum( (qqq(1:3,ikp)-q(1:3,iq))**2) <tolq) then
           lfind=.true.
        else
           call rangedq( matmul(ginv,q(1:3,iq)-qqq(:,ikp)), qx)
           if(sum(abs(qx))< tolq) lfind= .TRUE. 
        endif
        if(lfind) then
           ikpx=ikp
           goto 100
        endif
     enddo
     call rx( ' rsexx: not find ikp')
100  continue
     vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(itq(1:ntq),ikpx,1:nspin)
  enddo
end subroutine rsexx
