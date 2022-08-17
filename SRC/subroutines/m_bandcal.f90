module m_bandcal !band structure calculation
  use m_struc_def,only: s_rv1
  use m_suham,  only: ndhamx=>ham_ndhamx,nspx=>ham_nspx
  use m_qplist, only: nkp
  use m_lmfinit,only: nsp,nlibu,lmaxu,nbas,nl
  use m_mkqp,only: ntet=> bz_ntet, bz_nabc
  use m_qplist,only: qplist
  use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq
  use m_lmfinit,only: lrsig=>ham_lsig, lso,ham_scaledsigma,lekkl, &
       lmet=>bz_lmet,nbas,epsovl=>ham_oveps,nspc,plbnd,lfrce=>ctrl_lfrce, &
       pwmode=>ham_pwmode,pwemax,stdl
  use m_MPItk,only: mlog, master_mpi, procid,strprocid, numprocs=>nsize, mlog_MPIiq
  use m_subzi, only: nevmx,lswtk,rv_a_owtkb
  use m_supot, only: k1,k2,k3
  use m_mkpot,only: m_Mkpot_init,m_Mkpot_deallocate, osmpot,vconst, osig, otau, oppi,ohsozz,ohsopm
  use m_rdsigm2,only: senex,sene,getsenex,dsene,ndimsig
  use m_procar,only: m_procar_init,m_procar_closeprocar
  use m_clsmode,only: m_clsmode_set1
  use m_addrbl,only: Addrbl,swtk,Swtkzero
  use m_lgunit,only:stdo
  use m_augmbl,only: aughsoc
  use m_makusq,only: makusq
  use m_ftox
  !! outputs ---------------------------
  public m_bandcal_init,m_bandcal_2nd,m_bandcal_clean,m_bandcal_allreduce,m_bandcal_symsmrho
  integer,allocatable,protected,public::     ndimhx_(:),nev_(:),nevls(:,:)
  real(8),allocatable,protected,public::     evlall(:,:,:),frcband(:,:), orbtm_rv(:,:,:)
  complex(8),allocatable,protected,public::  smrho_out(:),dmatu(:,:,:,:)
  type(s_rv1),allocatable,protected,public:: sv_p_oeqkkl(:,:), sv_p_oqkkl(:,:)
  !! ------------------------------------------------
  logical,private:: debug,sigmamode,call_m_bandcal_2nd,procaron,writeham
  logical,private:: dmatuinit=.true.
  real(8),private::  sumqv(3,2),sumev(3,2)
  private
contains
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
    if(nlibu>0)  dmatu=0d0    !density matrix for U initialization
    sumev = 0d0
    sumqv = 0d0
    if (lswtk==1)  call swtkzero()
    bandcalculation_q: do 2010 iq = iqini, iqend 
       qp = qplist(:,iq)
       !write(stdo,ftox)'m_bandcal_init: procid iq=',procid,iq,ftof(qp)
       if(iq==iqini) call mlog_MPIiq(iq,iqini,iqend)
       call m_Igv2x_setiq(iq)    ! Get napw and so on for given qp
       if(lso==1) then
          allocate(hamm(ndimh,ndimh,4),ovlm(ndimh,ndimh,4)) !spin offdiagonal included
       else
          allocate(hamm(ndimh,ndimh,1),ovlm(ndimh,ndimh,1)) !only for one spin
       endif
       ispendx = nsp
       if(lso==1) ispendx=1
       bandcalculation_spin: do 2005 isp = 1,ispendx
          if(iq==iqini .AND. ispini==2 .AND. isp==1) cycle
          if(iq==iqend .AND. ispend==1 .AND. isp==2) cycle
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

          Setup_hamiltonian_and_diagonalize : block
            integer:: iprint
            character:: charnum3
            if(lso/=0 .AND. ( .NOT. allocated(hammhso))) then
               allocate(hammhso(ndimh,ndimh,3))
               call aughsoc(qp, ohsozz,ohsopm,ndimh, hammhso)! SOC part of Hamiltonian hammhso is calculated.
            endif
            wsene = cmdopt0('--writesene')
            if(wsene) then
               open(newunit=iwsene,file='sene.isp:'//charnum3(isp)//'_iq:'//charnum3(iq) &
                    ,form='unformatted')
               if(iq==1 .AND. isp==1) write(iwsene) nsp,ndimsig,bz_nabc,nkp,0,0,0
            endif
            hamm=0d0
            ovlm=0d0
            if(lso==1) then !L.S case nspc=2
               call hambl(1,qp,osmpot,vconst,osig,otau,oppi, hamm(1,1,1),ovlm(1,1,1))
               call hambl(2,qp,osmpot,vconst,osig,otau,oppi, hamm(1,1,2),ovlm(1,1,2))
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
               call hambl(isp,qp,osmpot,vconst,osig,otau,oppi,hamm(1,1,1),ovlm(1,1,1))
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
            !! nmx: input, number of requested eigenvalues(functions).
            !!      If nmx=0, no eigenfunctions but all eigenvalues. <== WARNNNNNNNNNN!
            !! nev: out number of obtained eigenfvalues(funcitons)
            diagonalize_hamilatonian: block 
              call zhev_tk4(ndimhx, hamm, ovlm, nmx, nev, evl(1, jsp ), evec, epsovl)
            endblock diagonalize_hamilatonian
          endblock Setup_hamiltonian_and_diagonalize
          if(writeham.AND.master_mpi) call prtev(evec, ndimhx , evl(1, jsp ) , nmx , nev )
          if(call_m_bandcal_2nd) then
             write(ifig) nev,nmx,ndimhx !nev: number of eigenvalues
             write(ifig) evl(1:nev,jsp)
             write(ifig) evec(1:ndimhx,1:nmx)
          endif
          evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
          nevls(iq,jsp) = nev   !nov2014 isp and jsp is confusing...
          evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
          if(master_mpi .AND. epsovl>=1.000001d-14 .AND. plbnd/=0) then
             write(stdo,"(' : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
          endif
          if(plbnd==0 .AND. lwtkb/=-1) then !lwtkb=-1,0,1
             if(nlibu>0.AND.nev>0.AND.lwtkb==0)call rx('metal weights required for LDA+U calculation')
             if(lso/=0)  call mkorbm(jsp, nev, iq,qp, evec,  orbtm_rv) !Orbital magnetic moment
             if(nlibu>0 .AND. nev>0) call mkdmtu(jsp, iq, qp, nev, evec,  dmatu) !density matrix
             if(cmdopt0('--cls'))  call m_clsmode_set1(nmx,jsp,iq,qp,nev,evec) !all inputs
          endif
          if(lrout/=0 .AND. lwtkb>=0) then ! accumulate output density and sampling DOS.
             call addrbl (jsp, qp &
                  , iq , lfrce,  osmpot,vconst,osig,otau,oppi &
                  , evec,evl,nev, smrho_out, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frcband)
          endif
          if(PROCARon) call m_procar_init(iq,isp,ef0,evl,ndimh,jsp,qp,nev,evec,ndimhx,nmx)
          if(allocated(evec)) deallocate(evec)
          continue !== end loop over isp (main loop in parallel mode)==
2005   enddo bandcalculation_spin
       if(allocated(hammhso)) deallocate(hammhso)
       if(allocated(hamm)) deallocate(hamm,ovlm)
       ndimhx_(iq)= ndimhx     !Ham dimenstion
       nev_(iq)   = nev        !calculated number of bands  
2010 enddo bandcalculation_q
    !! ... Average forces so net force on system is zero (APW case)
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       xv(:)=[(sum(frcband(i,1:nbas))/nbas,i=1,3)]
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
  subroutine m_bandcal_2nd(iqini,iqend,ispini,ispend,lrout)! accumule evec things by addrbl
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
    iqloop: do 12010 iq = iqini, iqend !This is a big iq loop
       qp = qplist(:,iq)
       call m_Igv2x_setiq(iq) !qp)   ! Get napw and so on for given qp
       ispendx = nsp
       if(lso==1) ispendx=1
       isploop: do 12005 isp = 1,ispendx
          if(iq==iqini .AND. ispini==2 .AND. isp==1) cycle
          if(iq==iqend .AND. ispend==1 .AND. isp==2) cycle
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
               , iq , lfrce,  osmpot,vconst,osig,otau,oppi &
               , evec,evl,nev, smrho_out, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frcband)
          if(allocated(evec)) deallocate(evec)
12005  enddo isploop
12010 enddo iqloop
    !! ... Average forces so net force on system is zero (APW case)
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       xv(:)=[(sum(frcband(i,1:nbas))/nbas,i=1,3)]
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
    call mpibc2_real(sumqv,size(sumqv),'bndfp_sumqv')
    call mpibc2_real(sumev,size(sumev),'bndfp_sumev')
    call mpibc2_complex(smrho_out,size(smrho_out),'bndfp_smrho')
    if(lswtk==1) call mpibc2_complex(swtk,size(swtk),'bndfp_swtk')
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
    if(lfrce/=0) nnn=size(frcband)
    if(lfrce/=0) call mpibc2_real(frcband,nnn,'bndfp_frcband')
    if(nlibu>0)  nnn=size(dmatu)
    if(nlibu>0)  call mpibc2_complex(dmatu,nnn,'bndfp_dmatu')
    if(lso/=0 .AND. lwtkb/=-1) call mpibc2_real(orbtm_rv,size(orbtm_rv),'bndfp_orbtm')
  end subroutine m_bandcal_allreduce

  subroutine m_bandcal_symsmrho()
    call tcn('m_bandcal_symsmrho')
    call symsmrho(smrho_out)
    call tcx('m_bandcal_symsmrho')
  end subroutine m_bandcal_symsmrho

  subroutine mkorbm(isp,nev,iq,qp,evec, orbtm)
    use m_lmfinit,only: ispec,sspec=>v_sspec,nbas,nlmax,nsp,nspc,nl,n0,nppn,nab
    use m_igv2x,only: napw,ndimh,ndimhx,igvapw=>igv2x
    use m_mkpot,only: ppnl=>ppnl_rv
    use m_subzi, only: wtkb=>rv_a_owtkb
    use m_qplist,only: nkp
    use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx
    !- Decomposition of norm from projection of w.f into MT sphere
    ! ----------------------------------------------------------------------
    !i   isp   :current spin channel (1 or 2)
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   nlmax :leading dimension of aus
    !i   ndham :dimensions aus,wtkp
    !i   nev   :number of eigenvectors to accumulate orbital moment
    !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
    !i   iq    :current k-point
    !i   nbas  :size of basis
    !i   ppnl  :NMTO potential parameters; see eg potpus.f
    !i   aus   :values of (phi,phidot) MT sphere boundary; see makusq
    !i   nl    :(global maximum l) + 1
    !i   nkp   :number of irreducible k-points (bzmesh.f)
    !o Outputs
    !o   orbtm :orbital moments accumulated for this qp
    !l Local variables
    !l   ispc  :the current spin index in the coupled spins case.
    !l         :Some quantities have no separate address space for each
    !l         :spin in the indepedent-spins case (evec,evl,ewgt) but do
    !l         :in the coupled-spins case.  A separate loop ispc=1..nspc
    !l         :must be added for the latter case
    !l         :ispc is the appropriate index for objects which distinguish
    !l         :spins in the spin-coupled case only
    !l   isp   :isp  is the appropriate index for objects which distinguish
    !l         :spins in the spin-uncoupled case only
    !l   ksp   :the current spin index in both independent and coupled
    !l         :spins cases.
    !l         :ksp is appropriate spin index for quantities that have
    !l         :separate address space for each spin in every case
    !l         :(potential- and density-like objects).
    !u Updates
    !u   25 Apr 05 (A. Chantis) extended to local orbitals
    !u   24 Dec 04 Extended to spin-coupled case
    !u   30 Aug 04 (A. Chantis) first written, adapted from mkpdos
    ! ----------------------------------------------------------------------
    implicit none
    integer :: isp,nev,iq
    integer :: lmxa,lmxax,lmdim,ichan,ib,is,igetss,iv,ilm,l,m,nlma, ll,lc,em,ispc,ksp
    real(8):: qp(3)
    real(8):: suml(11),s11,s22,s12,s33,s31,s32,s13,s23, suma,rmt,sab(nab,n0,2),orbtm(nl,nsp,*) 
    complex(8):: au,as,az,iot=(0d0,1d0),evec(ndimh,nsp,nev)
    complex(8),allocatable ::aus(:,:,:,:,:)
    allocate(aus(nlmax,ndham*nspc,3,nsp,nbas))
    aus=0d0
    call makusq(0 , nbas,[-999], nev, isp,1,qp,evec, aus )
    lmxax = ll(nlmax)
    iot = dcmplx(0d0,1d0)
    ichan = 0
    ibloop: do  ib = 1, nbas
       is = ispec(ib)
       lmxa=sspec(is)%lmxa
       rmt= sspec(is)%rmt
       lmxa = min(lmxa,lmxax)
       if (lmxa == -1) cycle 
       nlma = (lmxa+1)**2
       lmdim = nlma
       call phvsfp(nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)
       !       In noncollinear case, isp=1 always => need internal ispc=1..2
       !       ksp is the current spin index in both cases:
       !       ksp = isp  in the collinear case
       !           = ispc in the noncollinear case
       !       ispc=1 for independent spins, and spin index when nspc=2
       ispcloop: do  ispc = 1, nspc
          ksp = max(ispc,isp)
          ivloop: do  iv = 1, nev
             suml=0d0
             suma=0d0
             ilm = 0
             !  ....  Rotate from real to spherical harmonics (order assumed: m,...,-m).
             !        |Psi>_l = \Sum_{m}(A_l,m * u_l + B_l,m * s_l)*R_l,m --->
             !        |Psi>_l = \Sum_{m}(C_l,m * u_l + D_l,m * s_l)*Y_l,m
             !        R_l,m and Y_l,m are the real and spherical harmonics respectively.

             !              | (-1)^m/sqrt(2)*A_l,-m + i*(-1)^m/sqrt(2)*A_l,m , m>0
             !        C_l,m=|  A_l,m                                         , m=0
             !              |  1/sqrt(2)*A_l,-m -  i*1/sqrt(2)*A_l,m         , m<0

             !       Same relationships are valid between D and B.
             lloop: do  l = 0, lmxa
                lc = (l+1)**2 - l
                mloop: do  m = -l, l
                   em = abs(m)
                   ilm = ilm+1
                   !    ...    m,...,-m order
                   !            if (m .lt. 0) then
                   !            au = 1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) -
                   !     .           iot/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
                   !            as = 1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) -
                   !     .           iot/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
                   !            az = 1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) -
                   !     .           iot/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
                   !            else if (m .gt. 0) then
                   !            au = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) +
                   !     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
                   !            as = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) +
                   !     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
                   !            az = (-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) +
                   !     .           iot*(-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
                   !            else
                   !    ...   -m,...,m order
                   if (m < 0) then
                      au = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,1,ksp,ib) + &
                           1d0/dsqrt(2d0)*aus(lc+em,iv,1,ksp,ib)
                      as = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,2,ksp,ib) + &
                           1d0/dsqrt(2d0)*aus(lc+em,iv,2,ksp,ib)
                      az = iot*1d0/dsqrt(2d0)*aus(lc-em,iv,3,ksp,ib) + &
                           1d0/dsqrt(2d0)*aus(lc+em,iv,3,ksp,ib)
                   else if (m > 0) then
                      au = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,1,ksp,ib) + &
                           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,1,ksp,ib)
                      as = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,2,ksp,ib) + &
                           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,2,ksp,ib)
                      az = -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,3,ksp,ib) + &
                           (-1)**m/dsqrt(2d0)*aus(lc+m,iv,3,ksp,ib)
                   else
                      au = aus(ilm,iv,1,ksp,ib)
                      as = aus(ilm,iv,2,ksp,ib)
                      az = aus(ilm,iv,3,ksp,ib)
                   endif
                   !           If (au,as) are coefficients to (u,s), use this
                   s11 = dconjg(au)*au*sab(1,l+1,ksp)
                   s12 = 2*dconjg(au)*as*sab(2,l+1,ksp)
                   s22 = dconjg(as)*as*sab(4,l+1,ksp)
                   s33 = dconjg(az)*az*sab(5,l+1,ksp)
                   s31 = dconjg(az)*au*sab(6,l+1,ksp)
                   s32 = dconjg(az)*as*sab(7,l+1,ksp)
                   s13 = dconjg(au)*az*sab(6,l+1,ksp)
                   s23 = dconjg(as)*az*sab(7,l+1,ksp)
                   orbtm(l+1,ksp,ib)=orbtm(l+1,ksp,ib)+m*(s11+s12+s22+ &
                        s33+s32+s23+s31+s13)*wtkb(iv,isp,iq) !corrected. it should be wtkb
                enddo mloop
             enddo lloop
          enddo ivloop
       enddo ispcloop    !      print*, isp, ib, 'ORB.MOMNT=',orbtm(isp,ib)
    enddo ibloop
    deallocate(aus)
  end subroutine mkorbm
end module m_bandcal
