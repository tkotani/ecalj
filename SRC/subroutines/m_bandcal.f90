module m_bandcal
  use m_struc_def,only: s_rv1
  use m_suham,  only: ndhamx=>ham_ndhamx,nspx=>ham_nspx
  use m_qplist, only: nkp
  use m_lmfinit,only: nsp,nlibu,lmaxu,nbas,nl!,nlmto
  use m_mkqp,only: ntet=> bz_ntet
  ! ccccccccccc
  use m_mkqp,only: bz_nabc
  ! ccccccccccc
  use m_qplist,only: qplist
  use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq
  use m_lmfinit,only: lrsig=>ham_lsig, lso,ham_scaledsigma,lekkl, &
       lmet=>bz_lmet,stdo,nbas,epsovl=>ham_oveps,nspc,plbnd,lfrce=>ctrl_lfrce, &
       pwmode=>ham_pwmode,pwemax,stdl,iprmb
  use m_MPItk,only: mlog, master_mpi, procid,strprocid, numprocs=>nsize, mlog_MPIiq
  use m_subzi, only: nevmx,lswtk,rv_a_owtkb
  use m_supot, only: k1,k2,k3
  use m_mkpot,only: m_Mkpot_init,m_Mkpot_deallocate, osmpot,vconst, &
       sv_p_osig, sv_p_otau, sv_p_oppi,ohsozz,ohsopm
  use m_rdsigm2,only: senex,sene,getsenex,dsene,ndimsig
  use m_procar,only: m_procar_init,m_procar_closeprocar
  use m_clsmode,only: m_clsmode_set1
  use m_addrbl,only: Addrbl,swtk,Swtkzero
  !! outputs ---------------------------
  integer,allocatable,protected::     ndimhx_(:),nev_(:),nevls(:,:)
  real(8),allocatable,protected::     evlall(:,:,:),frcband(:,:), orbtm_rv(:,:,:)
  complex(8),allocatable,protected::  smrho_out(:),dmatu(:,:,:,:)
  type(s_rv1),allocatable,protected:: sv_p_oeqkkl(:,:), sv_p_oqkkl(:,:)
  !! ------------------------------------------------
  logical,private:: debug,sigmamode,call_m_bandcal_2nd,procaron,writeham
  logical,private:: dmatuinit=.true.
  real(8),private::  sumqv(3,2),sumev(3,2)

contains
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  ! Set up Hamiltonian, diagonalization
  subroutine m_bandcal_init(iqini,iqend,ispini,ispend,lrout,ef0,ifih,lwtkb)
    implicit none
    complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),hammhso(:,:,:) !Hamiltonian,Overlapmatrix
    integer:: iq,nmx,ispinit,isp,jsp,nev,ifih,lwtkb,iqini,iqend,lrout,ifig,ispini,ispend,ispendx
    real(8):: qp(3),ef0,def=0d0,xv(3),q(3)
    real(8),allocatable:: evl(:,:) !eigenvalue
    complex(8),allocatable :: evec(:,:) !eigenvector
    integer::  iprint,i,ibas,iwsene
    logical:: ltet,cmdopt0,dmatuinit=.true.,wsene
    character(3):: charnum3
    call tcn('m_bandcal_init')
    if(master_mpi) write(stdo,*)' m_bandcal_init: start'
    sigmamode = mod(lrsig,10)/=0
    writeham= cmdopt0('--writeham')
    PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
    debug    = cmdopt0('--debugbndfp')
    !! pdos mode = PROCARon .and. fullmesh (--mkprocar and --fullmesh. See job_pdos).
    ltet = ntet>0
    if(plbnd==0 .AND. lso/=0 .AND. lwtkb==0) call rx('metal weights required to get orb.moment')
    if(lso/=0) allocate(orbtm_rv(nl,nsp,nbas)) !for spin-orbit coupling
    if(lso/=0) orbtm_rv=0d0
    allocate( evlall(ndhamx,nspx,nkp))
    if(lfrce>0) then
       allocate( frcband(3,1:nbas)) !force for band
       frcband  = 0d0
    endif
    allocate( ndimhx_(nkp),nev_(nkp),nevls(nkp,nspx))
    ndimhx_=0
    nev_   =0
    nevls  =0

    if(nlibu>0 .AND. dmatuinit) then
       allocate( dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
       dmatuinit=.false.
    endif

    if(lrout/=0) then
       allocate( sv_p_oeqkkl(3,nbas), sv_p_oqkkl(3,nbas))
       call dfqkkl( sv_p_oqkkl ) !zero clear
       if(lekkl==1) call dfqkkl( sv_p_oeqkkl )!zero clear
       allocate( smrho_out(k1*k2*k3*nsp) )
       smrho_out = 0d0
    endif

    call_m_bandcal_2nd =.false.
    if(plbnd==0) call_m_bandcal_2nd= (lmet>=0 .AND. lrout>0 )
    if(call_m_bandcal_2nd) open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
    allocate( evl(ndhamx,nspx))

    !! These are accumulation varivables
    if(nlibu>0)  dmatu=0d0    !density matrix for U initialization
    sumev = 0d0
    sumqv = 0d0
    if (lswtk==1)  call swtkzero()
    do 2010 iq = iqini, iqend !This is a big iq loop
       qp = qplist(:,iq)
       if(debug) print *,' do 2010 iq procid=',iq,procid,iq,iqini,iqend
       if(iq==iqini) call mlog_MPIiq(iq,iqini,iqend)
       call m_Igv2x_setiq(iq)    ! Get napw and so on for given qp
       if(lso==1) then
          allocate(hamm(ndimh,ndimh,4),ovlm(ndimh,ndimh,4)) !spin offdiagonal included
       else
          allocate(hamm(ndimh,ndimh,1),ovlm(ndimh,ndimh,1)) !only for one spin
       endif
       ispendx = nsp
       if(lso==1) ispendx=1
       do 2005 isp = 1,ispendx
          ! ccccccccccccc
          if(iq==iqini .AND. ispini==2 .AND. isp==1) cycle
          if(iq==iqend .AND. ispend==1 .AND. isp==2) cycle
          ! ccccccccccccc
          jsp = isp
          if(lso==1) jsp = 1
          ! Hambl calls augmbl.
          ! See Appendix C in http://dx.doi.org/10.7566/JPSJ.84.034702
          !! finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
          !! Note that F2=Hankel head at a site + Hankel tail contributions from the other site.

          !! == Set up Hamiltonian by hambl. ==============
          !!    Hamiltonian: hamm(1:ndimh,1:ndimh,3) means off-diagonal section when SO=1.
          !!    Overlap matrix: ovlm
          !!    senex:  Sigma-Vxc
          !! ==========================================
          !! Generate sene(Sigma-Vxc) for given sfz.
          !! Determine interpolated self-energy sene at qp from sfz.
          !! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only)
          !!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
          !! Main input  => ham_iv_a_oiaxs,ham_rv_a_ohrs
          !     ! Main output sene. See m_seneinterp

          !! See Eq.(36) and appendix in http://dx.doi.org/10.7566/JPSJ.84.034702
          !! Hamm and ovlm are made from smooth part and augmentation part.

          ! SOC Hamiltonian hammhso is calculated.
          if(lso/=0 .AND. ( .NOT. allocated(hammhso))) then
             allocate(hammhso(ndimh,ndimh,3))
             call aughsoc(qp, ohsozz,ohsopm,ndimh, hammhso)!aug2021 Hso
          endif
          !!
          wsene = cmdopt0('--writesene')
          if(wsene) then
             open(newunit=iwsene,file='sene.isp:'//charnum3(isp)//'_iq:'//charnum3(iq) &
                  ,form='unformatted')
             if(iq==1 .AND. isp==1) write(iwsene) nsp,ndimsig,bz_nabc,nkp,0,0,0
          endif

          hamm=0d0
          ovlm=0d0
          if(lso==1) then !L.S case nspc=2
             call hambl(1,qp,osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi, hamm(1,1,1),ovlm(1,1,1))
             call hambl(2,qp,osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi, hamm(1,1,2),ovlm(1,1,2))
             hamm(:,:,1:2)= hamm(:,:,1:2)+hammhso(:,:,1:2) !spin-diag SOC elements (1,1), (2,2) added
             hamm(:,:, 3) = hammhso(:,:,3) !spin-offdiagonal SOC elements (1,2) added
             if(sigmamode) then
                call getsenex(qp, 1, ndimh,ovlm(1,1, 1))
                hamm(:,:,1) = hamm(:,:,1) + ham_scaledsigma*senex !senex_up= Vxc(QSGW)-Vxc(LDA)
                if(wsene) write(iwsene) qp,1
                if(wsene) write(iwsene) sene
                call dsene()
                call getsenex(qp, 2, ndimh,ovlm(1,1, 2))
                hamm(:,:,2) = hamm(:,:,2) + ham_scaledsigma*senex !senex_dn= Vxc(QSGW)-Vxc(LDA)
                if(wsene) write(iwsene) qp,2
                if(wsene) write(iwsene) sene
                call dsene()
             endif
             call sopert2 ( hamm , hamm ) !L.S case. re-ordered to be 2x2 spin matrix.
             call sopert2 ( ovlm , ovlm )
          else ! lso=0 (No SO) or lso=2(Lz.Sz)  Spin Diagonal case.spin diagonal)
             call hambl(isp,qp,osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,hamm(1,1,1),ovlm(1,1,1))
             if(lso==2) hamm(:,:, 1) = hamm(:,:, 1) + hammhso(:,:,isp)
             if(sigmamode) then !!Add  Vxc(QSGW)-Vxc
                call getsenex(qp,isp,ndimh,ovlm(1,1,1))
                hamm(:,:, 1) = hamm(:,:, 1) + ham_scaledsigma*senex !senex= Vxc(QSGW)-Vxc(LDA)
                if(wsene) write(iwsene) qp,isp
                if(wsene) write(iwsene) sene
                call dsene()
             endif
          endif
          if(wsene) close(iwsene)
          nmx=min(nevmx,ndimhx)!nmx: max number of eigenfunctions we will obtain. Smaller is faster.

          if(iprint()>=30) write(stdo,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4, &
               " ndimh = nmto+napw = ",3i5,f13.5)') iq,nkp,qp,ndimh,ndimh-napw,napw
          if(writeham) then
             !              write(stdo,"(a,3f9.5)") "Hamiltonian: Writing hamm and ovlm for qp= ",qp
             write(ifih) qp,ndimhx,lso,epsovl,jsp
             if(lso==1) then  !L.S case ndimhx=ndimh*nspc nspc=2
                write(ifih) ovlm !Note sopert2. When you read, use ovlm(1:ndimhx, 1:ndimhx)
                write(ifih) hamm
             else             !spin diagonal case nspc=1 ndimhx=ndimh
                write(ifih) ovlm ! When you read, use ovlm(1:ndimhx, 1:ndimhx)
                write(ifih) hamm
             endif
          endif

          allocate(evec(ndimhx,nmx))
          !! == Diagonalize Hamiltonian ==
          !! ndimhx: dimension of Hamitonian
          !! hamm:Hamiltonian, ovlm: overlap matrix
          !! evec:eigenfunciton. evl: eigenvalue.
          !!
          !! nmx: input, number of requested eigenvalues(functions).
          !!      If nmx=0, no eigenfunctions but all eigenvalues. <== WARNNNNNNNNNN!
          !! nev: out number of obtained eigenfvalues(funcitons)
          call zhev_tk4(ndimhx, hamm, ovlm, nmx, nev, evl(1, jsp ), evec, epsovl)
          if(writeham .AND. master_mpi) call prtev(evec, ndimhx , evl(1, jsp ) , nmx , nev )
          if(call_m_bandcal_2nd) then
             write(ifig) nev,nmx,ndimhx
             write(ifig) evl(1:nev,jsp)
             write(ifig) evec(1:ndimhx,1:nmx)
          endif

          !! nev: number of eigenvalues
          evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
          nevls(iq,jsp) = nev  !nov2014 isp and jsp is confusing...
          evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
          if(master_mpi .AND. epsovl>=1.000001d-14 .AND. plbnd/=0) then
             write(stdo,"(' : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
          endif

          if(plbnd==0 .AND. lwtkb/=-1) then !lwtkb=-1,0,1
             if(nlibu>0 .AND. nev>0 .AND. lwtkb==0) call rx('metal weights required for LDA+U calculation')
             if(lso/=0)  call mkorbm(jsp, nev, iq,qp, evec,  orbtm_rv) !Orbital magnetic moment
             if(nlibu>0 .AND. nev>0) call mkdmtu(jsp, iq, qp, nev, evec,  dmatu) !density matrix
             if(cmdopt0('--cls'))  call m_clsmode_set1(nmx,jsp,iq,qp,nev,evec) !all inputs
          endif

          if(lrout/=0 .AND. lwtkb>=0) then
             !               call readindensitymodesetting() !dummy
             ! ccumulate output density and sampling DOS.
             call addrbl (jsp, qp &
                  , iq , lfrce,  osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi &
                  , evec,evl,nev, smrho_out, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frcband)
          endif
          if(PROCARon) call m_procar_init(iq,isp,ef0,evl,ndimh,jsp,qp,nev,evec,ndimhx,nmx)
          if(allocated(evec)) deallocate(evec)
          continue !== end loop over isp (main loop in parallel mode)==
2005   enddo
       if(allocated(hammhso)) deallocate(hammhso)
       if(allocated(hamm)) deallocate(hamm,ovlm)
       ndimhx_(iq)= ndimhx
       nev_(iq)   = nev
       continue                  !end of iq loop
2010 enddo
    !! ... Average forces so net force on system is zero (APW case)
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       do i = 1, 3
          xv(i)=sum(frcband(i,1:nbas))/nbas
       enddo
       do ibas= 1, nbas
          frcband(:,ibas) = frcband(:,ibas) - xv(:)
       enddo
    endif
    if(PROCARon) call m_procar_closeprocar()
    if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid
    if(call_m_bandcal_2nd) close(ifig)
    deallocate(evl)
    call tcx('m_bandcal_init')
  end subroutine m_bandcal_init
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine m_bandcal_2nd(iqini,iqend,ispini,ispend,lrout)!,ef0) !,emin,emax)!,ndos)
    ! ccumulation by addrbl
    implicit none
    integer:: iq,nmx,ispinit,isp,jsp,nev,iqini,iqend,lrout,ifig,i,ibas,ispini,ispend,ispendx
    real(8):: qp(3),ef0,def=0d0,xv(3)
    real(8),allocatable:: evl(:,:)
    complex(8),allocatable :: evec(:,:)
    logical:: cmdopt0
    call tcn('m_bandcal_2nd')
    if(master_mpi) write(stdo,*)' mmmmm m_bandcal_2nd'
    call dfqkkl( sv_p_oqkkl ) !zero clear
    if(lekkl==1) call dfqkkl( sv_p_oeqkkl ) !zero clear
    if (lfrce>0)  frcband  = 0d0
    if (lswtk==1) call swtkzero()
    if(lso/=0) orbtm_rv=0d0
    if(allocated(smrho_out)) deallocate(smrho_out)
    allocate( smrho_out(k1*k2*k3*nsp) )
    smrho_out = 0d0
    sumev = 0d0
    sumqv = 0d0
    allocate(evl(ndhamx,nspx))
    open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
    do 12010 iq = iqini, iqend !This is a big iq loop
       qp = qplist(:,iq)
       call m_Igv2x_setiq(iq) !qp)   ! Get napw and so on for given qp
       ispendx = nsp
       if(lso==1) ispendx=1
       do 12005 isp = 1,ispendx
          ! ccccccccccccc
          if(iq==iqini .AND. ispini==2 .AND. isp==1) cycle
          if(iq==iqend .AND. ispend==1 .AND. isp==2) cycle
          ! ccccccccccccc
          jsp = isp
          if(lso==1) jsp = 1
          read(ifig) nev,nmx  !ndimhx <---supplied by m_Igv2x_set
          if (allocated(evec)) deallocate(evec)
          allocate(evec(ndimhx,nmx))
          read(ifig) evl(1:nev,jsp)
          read(ifig) evec(1:ndimhx,1:nmx)
          evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
          if( lso/=0)            call mkorbm(jsp, nev, iq,qp, evec,  orbtm_rv)
          if( nlibu>0 .AND. nev>0) call mkdmtu(jsp, iq,qp, nev, evec,  dmatu)
          if( cmdopt0('--cls'))  call m_clsmode_set1(nmx,jsp,iq,qp,nev,evec) !all inputs
          call addrbl (jsp, qp &
               , iq , lfrce,  osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi &
               , evec,evl,nev, smrho_out, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frcband)
          if(allocated(evec)) deallocate(evec)
12005  enddo
12010 enddo
    !! ... Average forces so net force on system is zero (APW case)
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       do i = 1, 3
          xv(i)=sum(frcband(i,1:nbas))/nbas
       enddo
       do  ibas= 1, nbas
          frcband(:,ibas) = frcband(:,ibas) - xv(:)
       enddo
    endif
    close(ifig)
    deallocate(evl)
    call tcx('m_bandcal_2nd')
  end subroutine m_bandcal_2nd
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine m_bandcal_clean()
    if (allocated(orbtm_rv)) deallocate(orbtm_rv)
    !      if (allocated(dos_rv))   deallocate(dos_rv)
    if (allocated(smrho_out)) deallocate(smrho_out)
    if (allocated(frcband))  deallocate(frcband)
    if (allocated(ndimhx_))  deallocate(ndimhx_,nev_,nevls)
    if(allocated(sv_p_oqkkl)) deallocate( sv_p_oqkkl)
    if(allocated(sv_p_oeqkkl))deallocate( sv_p_oeqkkl)
    deallocate(evlall)
  end subroutine m_bandcal_clean
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  !     !  Allreduce density-related quantities
  subroutine m_bandcal_allreduce(lwtkb)
    integer:: nnn,ib,i,lwtkb
    if(debug) print *,'goto m_bandcal_allreduce'
    !      if (lrout .ne. 0) then
    nnn=size(sumqv);    call mpibc2_real(sumqv,nnn,'bndfp_sumqv')
    nnn=size(sumev);    call mpibc2_real(sumev,nnn,'bndfp_sumev')
    nnn=size(smrho_out); call mpibc2_complex(smrho_out,nnn,'bndfp_smrho')
    if (lswtk==1) then
       nnn=size(swtk);  call mpibc2_complex(swtk,'bndfp_swtk')
    endif
    do  ib = 1, nbas
       do  i = 1, 3
          if(allocated(sv_p_oqkkl(i,ib)%v)) then
             nnn = size(sv_p_oqkkl(i,ib)%v)
             if(nnn>0) call mpibc2_real(sv_p_oqkkl(i,ib)%v,nnn,'bndfp_qkkl')
          endif
          if(lekkl==1 .AND. allocated(sv_p_oeqkkl(i,ib)%v)) then
             nnn = size(sv_p_oeqkkl(i,ib)%v)
             if(nnn>0) call mpibc2_real(sv_p_oeqkkl(i,ib)%v,nnn,'bndfp_eqkkl')
          endif
       enddo
    enddo
    !         if( ndos>0 ) nnn=size(dos_rv)
    !         if( ndos>0 ) call mpibc2_real(dos_rv,nnn,'bndfp_dos')
    if(lfrce/=0) nnn=size(frcband)
    if(lfrce/=0) call mpibc2_real(frcband,nnn,'bndfp_frcband')
    if(nlibu>0)  nnn=size(dmatu)
    if(nlibu>0)  call mpibc2_complex(dmatu,nnn,'bndfp_dmatu')
    if(lso/=0 .AND. lwtkb/=-1) nnn=size(orbtm_rv)
    if(lso/=0 .AND. lwtkb/=-1) call mpibc2_real(orbtm_rv,nnn,'bndfp_orbtm')
    !      endif
  end subroutine m_bandcal_allreduce

  !!------------------------------
  !$$$      subroutine m_bandcal_dosw(lwtkb,lrout,  dosw,evtop,ecbot) !goto99)
  !$$$      ! lwtkb,lrout,def are used only for
  !$$$      use m_lmfinit,only:ctrl_zbak,bz_w
  !$$$      intent(in)   ::           lwtkb,lrout
  !$$$      intent(inout)::                             dosw !input is just initial guess
  !$$$      intent(out)  ::                                  evtop,ecbot
  !$$$c      intent(out)::                                  goto99
  !$$$!!  ===   Repeat loop for printout and goto 99 ===
  !$$$!!   jsp=isp in the collinear case; jsp=1 in the noncollinear
  !$$$!!     Thus jsp should be used in place of isp
  !$$$!     !     isp serves as a flag for the noncollinear case
  !$$$!     ! block10
  !$$$      logical:: ltet !goto99,
  !$$$      real(8):: ef00,ef0,dosw(2),qp(3),qbg,evl(ndhamx,nspx),dum,ebot,ecbot,evtop
  !$$$      integer:: iq,ipr,lwtkb,isp,jsp,nev_iq,i,lrout,iprint
  !$$$      character(10):: i2char
  !$$$c      goto99=.false.
  !$$$      ltet = ntet>0
  !$$$      qbg = ctrl_zbak(1) !homogenious background charge
  !$$$      ipr=iprint()
  !$$$      evtop = -99999
  !$$$      ecbot = -evtop
  !$$$      ebot = 1000d0
  !$$$      do  iq = 1, nkp
  !$$$         qp=qplist(:,iq)
  !$$$         do isp = 1, nspx
  !$$$            jsp = isp
  !$$$            evl(1:ndhamx,jsp) = evlall(1:ndhamx,jsp,iq)
  !$$$            nev_iq    = nev_(iq)
  !$$$            if(plbnd==0.and.ipr>=10 .and. iq==1) write (stdl,"('fp evl',8f8.4)") (evl(i,jsp),i=1,nev_iq)
  !$$$            ebot = dmin1(ebot,evl(1,jsp))
  !$$$            i = max(1,nint(qval-qbg)/(3-nspc))
  !$$$            evtop = max(evtop,evl(i,jsp))
  !$$$            ecbot = min(ecbot,evl(i+1,jsp))
  !$$$            if (lmet==0 .and. iq==1 .and. jsp==1) ef0 = evtop
  !$$$            if(debug) print *,'eeeee44444444444 plbnd=',plbnd
  !$$$            if(plbnd==0.and.lwtkb/=-1) then
  !$$$                  if (iq .eq. 1 .and. jsp .eq. nsp ) then !
  !$$$                     ef00 = ef0
  !$$$                     call fixef0(qval-qbg,jsp,1,nev_iq,ndhamx,evl,dosw,ef0)
  !$$$!!        :on output, dosw is revised if ebot<dosw(1) or dosw(2)<ef0
  !$$$                     if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
  !$$$     .                    lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
  !$$$                          if (master_mpi) write(stdo,"(a)")
  !$$$     .                       ' ... Fermi level reset in second spin channel ... restart band pass'
  !$$$                          call rx0('tk think ecalj is going though not maintained branch).')
  !$$$                         !goto99=.true.
  !$$$                        exit    !this was the case of make co test at ecalj/TestInstall/ (did I remove this?)
  !$$$                     endif
  !$$$                  endif
  !$$$            endif
  !$$$!!     Check for cases when nevmx is too small : i=2 => fatal error
  !$$$            if(plbnd==0.and.lwtkb/=-1) then
  !$$$                  i = 0
  !$$$                  if (nevmx.ge.0 .and. lmet .ne. 0) then
  !$$$                     dum = evl(max(nev_iq,1),jsp)
  !$$$                     if (.not. ltet .and. ef0+5*bz_w .gt. dum) i=2
  !$$$c                     if (lmet.eq.4 .and. ef0+def+5*bz_w .gt.dum)i=2
  !$$$                  endif
  !$$$                  if(i .eq. 2) then
  !$$$                     write(stdo,"(a,f13.5,a, f13.5)")
  !$$$     &                 'evl(nev='//trim(i2char(nev_iq))//')=',
  !$$$     &                 evl(max(nev_iq,1),jsp),' but ef0=',ef0
  !$$$                     call rx('bndfp:... restart with larger nevmx: bndfp')
  !$$$                  endif
  !$$$
  !$$$            endif
  !$$$         enddo                  ! end second loop over isp
  !$$$      enddo                     !end second loop over iq
  !$$$      end subroutine

  subroutine m_bandcal_symsmrho()
    call tcn('m_bandcal_symsmrho')
    call symsmrho(smrho_out)
    call tcx('m_bandcal_symsmrho')
  end subroutine m_bandcal_symsmrho

end module m_bandcal