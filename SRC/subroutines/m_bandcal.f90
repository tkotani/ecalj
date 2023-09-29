!>band structure calculation
module m_bandcal 
  use m_struc_def,only: s_rv1,s_rv5
  use m_suham,  only: ndhamx=>ham_ndhamx,nspx=>ham_nspx
  use m_qplist, only: nkp
  use m_mkqp,only: ntet=> bz_ntet, bz_nabc
  use m_qplist,only: qplist,niqisp,iqproc,isproc
  use m_igv2x,only: m_igv2x_setiq, napw,ndimh,ndimhx,igv2x
  use m_lmfinit,only: lrsig=>ham_lsig, lso,ham_scaledsigma,lmet=>bz_lmet,nbas,epsovl=>ham_oveps,nspc,plbnd,lfrce,&
       pwmode=>ham_pwmode,pwemax,stdl,nsp,nlibu,lmaxu,lmxax
  use m_MPItk,only: mlog, master_mpi, procid,strprocid, numprocs=>nsize, mlog_mpiiq
  use m_subzi, only: nevmx,rv_a_owtkb
  use m_supot, only: n1,n2,n3
  use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate, osmpot,vconst, osig, otau, oppi,ohsozz,ohsopm
  use m_rdsigm2,only: senex,sene,getsenex,dsene,ndimsig
  use m_procar,only: m_procar_init,m_procar_closeprocar
  use m_clsmode,only: m_clsmode_set1
  use m_addrbl,only: addrbl!,swtk,Swtkzero
  use m_lgunit,only:stdo
  use m_augmbl,only: aughsoc
  use m_makusq,only: makusq
  use m_zhev,only: zhev_tk4
  use m_ftox
  use m_hambl,only: hambl
  ! outputs ---------------------------
  public m_bandcal_init, m_bandcal_2nd, m_bandcal_clean, m_bandcal_allreduce, m_bandcal_symsmrho
  integer,allocatable,protected,public::     ndimhx_(:,:),nevls(:,:) 
  real(8),allocatable,protected,public::     frcband(:,:), orbtm_rv(:,:,:),evlall(:,:,:)
  complex(8),allocatable,protected,public::  smrho_out(:),dmatu(:,:,:,:)
  type(s_rv5),allocatable,protected,public:: oeqkkl(:,:), oqkkl(:,:)
  !------------------------------------------------
  logical,private:: debug,sigmamode,call_m_bandcal_2nd,procaron,writeham
  logical,private:: dmatuinit=.true.
  real(8),private:: sumqv(3,2),sumev(3,2)
  integer,allocatable,private::neviqis(:),ndimhxiqis(:)
  real(8),allocatable,private::evliqis(:,:)
  complex(8),allocatable,private:: eveciqis(:,:,:)
  private
contains
  subroutine m_bandcal_init(lrout,ef0,vmag0,ifih) ! Set up Hamiltonian, diagonalization
    implicit none
    intent(in)::            lrout,ef0,vmag0,ifih
    complex(8),allocatable:: hamm(:,:,:,:),ovlm(:,:,:,:),hammhso(:,:,:) !Hamiltonian,Overlapmatrix
    integer:: iq,nmx,ispinit,isp,nev,ifih,lwtkb,lrout,ifig,i,ibas,iwsene,idat,ikp
    real(8):: qp(3),ef0,def=0d0,xv(3),q(3),vmag0
    real(8),allocatable    :: evl(:,:)  !eigenvalue (nband,nspin)
    complex(8),allocatable :: evec(:,:) !eigenvector( :,nband)
    logical:: ltet,cmdopt0,dmatuinit=.true.,wsene
    character(3):: charnum3
!    real(8)::  evlall(ndhamx,nspx,nkp)
    call tcn('m_bandcal_init')
    if(master_mpi) write(stdo,ftox)' m_bandcal_init: start'
    sigmamode = mod(lrsig,10)/=0
    writeham = cmdopt0('--writeham')
    PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
    debug    = cmdopt0('--debugbndfp')
    ltet = ntet>0
!    nspx=nsp/nspc
    if(plbnd==0 .AND. lso/=0 .AND. lmet==0 ) call rx('metal weights required to get orb.moment')
    if(lso/=0) allocate(orbtm_rv(lmxax+1,nsp,nbas),source=0d0) !for spin-orbit coupling
    if(lfrce>0) allocate( frcband(3,1:nbas),source=0d0) !force for band
    allocate( ndimhx_(nkp,nspx),nevls(nkp,nspx),source=0) 
    allocate( evlall(ndhamx,nspx,nkp),source=0d0)
    if(nlibu>0 .AND. dmatuinit) then
       allocate( dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
       dmatuinit=.false.
    endif
    if(nlibu>0)  dmatu=0d0    !density matrix for U initialization
    if(lrout/=0) then
       allocate( oeqkkl(3,nbas), oqkkl(3,nbas))
       call dfqkkl( oqkkl  )!allocate and zero clear
       call dfqkkl( oeqkkl )!
       allocate( smrho_out(n1*n2*n3*nsp),source=(0d0,0d0) )
    endif
    call_m_bandcal_2nd =.false.
    if(plbnd==0) call_m_bandcal_2nd= (lmet>=0 .AND. lrout>0 )
!    if(call_m_bandcal_2nd) open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
    if(call_m_bandcal_2nd) then
       if(allocated(neviqis))deallocate(neviqis,ndimhxiqis,evliqis,eveciqis)
       allocate(neviqis(niqisp),ndimhxiqis(niqisp),evliqis(nevmx,niqisp),eveciqis(ndhamx,nevmx,niqisp))
    endif
    allocate( evl(ndhamx,nspx))
    sumev = 0d0
    sumqv = 0d0
!    if (lswtk==1)  call swtkzero() !write(stdo,*)'iiiiqqq procid iqini iqend=',procid,iqini,iqend
    bandcalculation_q: do 2010 idat=1,niqisp
       iq = iqproc(idat)
       qp = qplist(:,iq) !write(stdo,ftox)'m_bandcal_init: procid iq=',procid,iq,ftof(qp)
       isp= isproc(idat) !NOTE: isp=1:nspx=nsp/nspc
       if(cmdopt0('--afsym').and.isp==2) cycle
       call m_Igv2x_setiq(iq) ! Get napw,ndimh,ndimhx, and igv2x
       allocate(hamm(ndimh,nspc,ndimh,nspc),ovlm(ndimh,nspc,ndimh,nspc)) !spin offdiagonal included
       ! hambl calls augmbl. See Appendix C in http://dx.doi.org/10.7566/JPSJ.84.034702
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
         integer:: iprint,ispc
         character:: charnum3
         if(lso/=0 .AND. ( .NOT. allocated(hammhso))) then
            allocate(hammhso(ndimh,ndimh,3))
            call aughsoc(qp, ohsozz,ohsopm,ndimh, hammhso)! SOC part of Hamiltonian hammhso is calculated.
         endif
         wsene = cmdopt0('--writesene')
         if(wsene) then
            open(newunit=iwsene,file='sene.isp:'//charnum3(isp)//'_iq:'//charnum3(iq),form='unformatted')
            if(iq==1 .AND. isp==1) write(iwsene) nsp,ndimsig,bz_nabc,nkp,0,0,0
         endif
         ovlm=0d0
         if(lso==1) then !L.S case nspc=2
            do ispc=1,nspc
               call hambl(ispc,qp,osmpot,vconst,osig,otau,oppi, hamm(:,ispc,:,ispc),ovlm(:,ispc,:,ispc))
               hamm(:,ispc,:,ispc)= hamm(:,ispc,:,ispc) + hammhso(:,:,ispc) !spin-diag SOC elements (1,1), (2,2) added
            enddo
            hamm(:,1,:,2)= hammhso(:,:,3)                    !spin-offdiagonal SOC elements (1,2) added
            hamm(:,2,:,1)= transpose(dconjg(hammhso(:,:,3)))
            if(sigmamode) then
               do ispc=1,nspc
                  call getsenex(qp, 1, ndimh, ovlm(:,ispc,:,ispc))
                  hamm(:,ispc,:,ispc) = hamm(:,ispc,:,ispc) + ham_scaledsigma*senex !senex_up= Vxc(QSGW)-Vxc(LDA)
                  if(wsene) write(iwsene) qp,ispc
                  if(wsene) write(iwsene) sene
                  call dsene()
               enddo
            endif
         else ! lso=0 (No SO) or lso=2(Lz.Sz)  Spin Diagonal case.spin diagonal)
            call hambl(isp,qp,osmpot,vconst,osig,otau,oppi,hamm(:,1,:,1), ovlm(:,1,:,1))
            if(lso==2) hamm(:,1,:,1) = hamm(:,1,:,1) + hammhso(:,:,isp)
            if(sigmamode) then !!Add  Vxc(QSGW)-Vxc
               call getsenex(qp,isp,ndimh,ovlm(:,1,:,1))
               hamm(:,1,:, 1) = hamm(:,1,:,1) + ham_scaledsigma*senex !senex= Vxc(QSGW)-Vxc(LDA)
               if(wsene) write(iwsene) qp,isp
               if(wsene) write(iwsene) sene
               call dsene()
            endif
         endif
         if(wsene) close(iwsene)
         nmx=min(nevmx,ndimhx)!nmx:maximum number of eigenfunctions we will obtain. Smaller is faster.
         if(iprint()>=30) write(stdo,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4, &
              " ndimh = nmto+napw = ",3i5,f13.5)') iq,nkp,qp,ndimh,ndimh-napw,napw
         if(writeham) then
            write(ifih) qp,ndimhx,lso,epsovl,isp ! ndimhx=ndimh*nspc 
            write(ifih) ovlm ! When you read, use ovlm(1:ndimhx, 1:ndimhx)
            write(ifih) hamm
         endif
         allocate(evec(ndimhx,nmx))
         !== Diagonalize Hamiltonian ==
         ! ndimhx: dimension of Hamitonian
         ! hamm:Hamiltonian, ovlm: overlap matrix
         ! evec:eigenfunciton. evl: eigenvalue.
         ! nmx: input, number of requested eigenvalues(functions).
         !      If nmx=0, no eigenfunctions but all eigenvalues. <== WARNNNNNNNNNN!
         ! nev: out number of obtained eigenfvalues(funcitons)
         diagonalize_hamilatonian: block 
           call zhev_tk4(ndimhx, hamm, ovlm, nmx, nev, evl(1, isp ), evec, epsovl)
         endblock diagonalize_hamilatonian
       endblock Setup_hamiltonian_and_diagonalize
       if(writeham.AND.master_mpi) write(stdo,"(9f8.4)") (evl(i,isp), i=1,nev)
       if(call_m_bandcal_2nd) then
          neviqis(idat)   =nev
          ndimhxiqis(idat)=ndimhx
          evliqis(1:nev,idat)=evl(1:nev,isp)
          eveciqis(1:ndimhx,1:nev,idat)=evec(1:ndimhx,1:nev)
          !write(ifig) nev,ndimhx !nev: number of eigenvalues
          !write(ifig) evl(1:nev,isp)
          !write(ifig) evec(1:ndimhx,1:nev)
       endif
       evl(nev+1:ndhamx,isp)=1d99  !padding. flag to skip these data
       nevls(iq,isp)  = nev        !nov2014 isp and isp is confusing...
       ndimhx_(iq,isp)= ndimhx     !Hamiltonian dimension
       evlall(1:ndhamx,isp,iq) = evl(1:ndhamx,isp)
       if(cmdopt0('--afsym')) then
          evlall(1:ndhamx,2,iq) = evl(1:ndhamx,1)
          nevls(iq,2)  = nev        
          ndimhx_(iq,2)= ndimhx     !Hamiltonian dimension
       endif   
       if(master_mpi.AND.epsovl>=1d-14.AND.plbnd/=0) write(stdo,&
            "(' : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
       if(PROCARon) call m_procar_init(iq,isp,ef0,vmag0,evl,ndimh,qp,nev,evec,ndimhx)
       if(allocated(evec)) deallocate(evec)
       if(allocated(hammhso)) deallocate(hammhso)
       if(allocated(hamm)) deallocate(hamm,ovlm)
2010 enddo bandcalculation_q
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       xv(:)=[(sum(frcband(i,1:nbas))/nbas,i=1,3)]
       do ibas= 1, nbas
          frcband(:,ibas) = frcband(:,ibas) - xv(:) ! Average forces so net force on system is zero (APW case)
       enddo
    endif
    if(PROCARon) call m_procar_closeprocar()
    if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid !if(call_m_bandcal_2nd) close(ifig)
    deallocate(evl)
    call tcx('m_bandcal_init')
  end subroutine m_bandcal_init
  subroutine m_bandcal_2nd(lrout)! accumule evec things by addrbl
    implicit none
    integer:: iq,ispinit,isp,nev,lrout,ifig,i,ibas,idat
    real(8):: qp(3),def=0d0,xv(3)
    real(8),allocatable:: evl(:,:)
    complex(8),allocatable :: evec(:,:)!,evecbackup(:,:)
    logical:: cmdopt0
    call tcn('m_bandcal_2nd')
    if(master_mpi) write(stdo,ftox)
    if(master_mpi) write(stdo,ftox)' m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi'
    call dfqkkl( oqkkl ) !zero clear
    call dfqkkl( oeqkkl ) !zero clear if(lekkl==1) 
    if (lfrce>0)  frcband  = 0d0
    if(lso/=0) orbtm_rv=0d0
    if(allocated(smrho_out)) deallocate(smrho_out)
    allocate( smrho_out(n1*n2*n3*nsp) )
    smrho_out = 0d0
    sumev = 0d0
    sumqv = 0d0
    allocate(evl(ndhamx,nspx))
!    open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
    iqloop: do 12010 idat=1,niqisp !iq = iqini, iqend !This is a big iq loop
       iq = iqproc(idat)
       qp = qplist(:,iq)  !write(stdo,ftox)'m_bandcal_init: procid iq=',procid,iq,ftof(qp)
       isp= isproc(idat)
       if(cmdopt0('--afsym').and.isp==2) cycle
       call m_Igv2x_setiq(iq) ! Get napw and so on for given qp
       !read(ifig) nev,ndimhx  !ndimhx <---supplied by m_Igv2x_set
       !allocate(evec(ndimhx,nev))
       !read(ifig) evl(1:nev,isp)
       !read(ifig) evec(1:ndimhx,1:nev)
       !evl(nev+1:ndhamx,isp)=1d99 !to skip these data

       nev   = neviqis(idat)
       ndimhx= ndimhxiqis(idat)
       allocate(evec(ndimhx,nev))
       evl(1:nev,isp)=evliqis(1:nev,idat)
       evl(nev+1:ndhamx,isp)=1d99 !padding 
       evec(1:ndimhx,1:nev)=eveciqis(1:ndimhx,1:nev,idat)
       if(lso/=0)              call mkorbm(isp, nev, iq,qp, evec,  orbtm_rv)
       if(nlibu>0 .AND. nev>0) call mkdmtu(isp, iq,qp, nev, evec,  dmatu)
       if(cmdopt0('--cls'))    call m_clsmode_set1(nev,isp,iq,qp,nev,evec) !all inputs
       call addrbl(isp,qp,iq, osmpot,vconst,osig,otau,oppi,evec,evl,nev, smrho_out, sumqv, sumev, oqkkl,oeqkkl, frcband)

       afsymifffffffffffffffff:  if(cmdopt0('--afsym')) then
          if(idat==1) write(stdo,ftox)'--afsym:'
          afsymblock: block !isp2 is given by isp=1
            use m_rotwave,only:  rotevec
            use m_mksym,only: symops,ngrp,ngrpAF,ag
            use m_qplist,only: qplist,nkp
            use m_lattic,only: plat=>lat_plat
            use m_ftox
            logical:: cmdopt0
            integer:: igrp,isp2,ikp,iev,ndeltaG(3),ikpx
            real(8):: qtarget(3),platt(3,3),diffq(3),tol=1d-4,qpr(3)
            complex(8):: evecrot(ndimhx,nev)
            platt=transpose(plat)
            evl(1:nev,2)=evl(1:nev,1)
            do igrp = ngrp + 1, ngrp+ngrpAF !AF symmetry
               do ikp=1,nkp ! write(stdo,ftox)'ssssym',iq,igrp,ikp,ftof(matmul(symops(:,:,igrp),qplist(:,ikp)),3)
                  diffq = matmul(platt, (qplist(:,ikp)-matmul(symops(:,:,igrp),qp)) )
                  if(sum(abs(diffq-nint(diffq)))<tol) goto 1018
               enddo
            enddo
            DebugWrite: block
              write(stdo,ftox)'error afsymmode'
              write(stdo,ftox)' igrp ikp',igrp,ikp,'qp=',ftof(qp,3),'qplist=',ftof(qplist(:,ikp),3) !,'deltaG=',ndeltaG
              write(stdo,ftox)'qp=',ftof(qp,3)
              do ikpx=1,nkp
                 write(stdo,ftox)'qplist=',ikpx,ftof(qplist(:,ikpx),3)
              enddo
              call rx('can not find qtarget by afsym')
            endblock DebugWrite
1018        continue!write(stdo,ftox)'ikp qp=',ikp,ftof(qp,3),'is mapped to',ftof(matmul(symops(:,:,igrp),qp),3),' by symops igp=',igrp
            qpr = qplist(:,ikp)
            isp2 = 2
            call rotevec(igrp,qp, qpr,ndimhx,napw,nev,evec(:,1:nev), evecrot(:,1:nev))! evec at qp is roteted to be evecrot at qpr by symops(:,:,igrp)
            if( lso/=0)              call mkorbm(isp2, nev, ikp,qpr, evecrot,  orbtm_rv)
            if( nlibu>0 .AND. nev>0) call mkdmtu(isp2,      ikp,qpr, nev, evecrot,  dmatu)
            if( cmdopt0('--cls'))    call m_clsmode_set1(nev,isp2,ikp,qpr,nev,evecrot) 
            call addrbl(isp2,qpr,ikp, osmpot,vconst,osig,otau,oppi,evecrot,evl,nev, smrho_out, sumqv, sumev, oqkkl,oeqkkl, frcband)
          endblock afsymblock
       endif afsymifffffffffffffffff
       deallocate(evec)
12010 enddo iqloop
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       xv(:)=[(sum(frcband(i,1:nbas))/nbas,i=1,3)]
       do  ibas= 1, nbas
          frcband(:,ibas) = frcband(:,ibas) - xv(:) ! Average forces so net force on system is zero (APW case)
       enddo
    endif
!    close(ifig)
    deallocate(evl)
    call tcx('m_bandcal_2nd')
  end subroutine m_bandcal_2nd
  subroutine m_bandcal_clean() !cleaning allocation
    if (allocated(orbtm_rv)) deallocate(orbtm_rv)
    if (allocated(smrho_out)) deallocate(smrho_out)
    if (allocated(frcband))  deallocate(frcband)
    if (allocated(ndimhx_))  deallocate(ndimhx_,nevls,evlall)
    if(allocated(oqkkl)) deallocate( oqkkl)
    if(allocated(oeqkkl))deallocate( oeqkkl)
  end subroutine m_bandcal_clean
  subroutine m_bandcal_allreduce()!lwtkb) !  Allreduce density-related quantities
    integer:: nnn,ib,i
    if(debug) print *,'goto m_bandcal_allreduce'
    call mpibc2_real(sumqv,size(sumqv),'bndfp_sumqv')
    call mpibc2_real(sumev,size(sumev),'bndfp_sumev')
    call mpibc2_complex(smrho_out,size(smrho_out),'bndfp_smrho')
!    if(lswtk==1) call mpibc2_complex(swtk,size(swtk),'bndfp_swtk')
    do  ib = 1, nbas
       do  i = 1, 3
          if(allocated(oqkkl(i,ib)%v)) then
             nnn = size(oqkkl(i,ib)%v)
             if(nnn>0) call mpibc2_real(oqkkl(i,ib)%v,nnn,'bndfp_qkkl')
          endif
          if(allocated(oeqkkl(i,ib)%v)) then !lekkl==1 
             nnn = size(oeqkkl(i,ib)%v)
             if(nnn>0) call mpibc2_real(oeqkkl(i,ib)%v,nnn,'bndfp_eqkkl')
          endif
       enddo
    enddo
    if(lfrce/=0) nnn=size(frcband)
    if(lfrce/=0) call mpibc2_real(frcband,nnn,'bndfp_frcband')
    if(nlibu>0)  nnn=size(dmatu)
    if(nlibu>0)  call mpibc2_complex(dmatu,nnn,'bndfp_dmatu')
    if(lso/=0) call mpibc2_real(orbtm_rv,size(orbtm_rv),'bndfp_orbtm')
  end subroutine m_bandcal_allreduce
  subroutine m_bandcal_symsmrho()
    call tcn('m_bandcal_symsmrho')
    call symsmrho(smrho_out)
    call tcx('m_bandcal_symsmrho')
  end subroutine m_bandcal_symsmrho
  subroutine mkorbm(isp,nev,iq,qp,evec, orbtm) !decomposed orbital moments within MT
    use m_ll,only:ll
    use m_lmfinit,only: ispec,sspec=>v_sspec,nbas,nlmax,nsp,nspc,n0,nppn,lmxax
    use m_igv2x,only: napw,ndimh,ndimhx,igvapw=>igv2x
    use m_mkpot,only: sab_rv
    use m_subzi, only: wtkb=>rv_a_owtkb
    use m_qplist,only: nkp
    use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx
    !i   isp   :current spin channel (1 or 2)
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   nlmax :leading dimension of aus
    !i   ndham :dimensions aus,wtkp
    !i   nev   :number of eigenvectors to accumulate orbital moment
    !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
    !i   iq    :current k-point
    !i   nbas  :size of basis
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
    integer :: lmxa,lmdim,ichan,ib,is,igetss,iv,ilm,l,m,nlma, lc,em,ispc,ksp
    real(8):: qp(3),diff
    real(8):: suml(11),s11,s22,s12,s33,s31,s32,s13,s23, suma,rmt,orbtm(lmxax+1,nsp,*) 
    complex(8):: au,as,az,iot=(0d0,1d0),evec(ndimh,nsp,nev),auasaz(3)
    complex(8),allocatable ::aus(:,:,:,:,:)
!    real(8):: sab(nab,n0,2)
    allocate(aus(nlmax,ndham*nspc,3,nsp,nbas))
    call makusq(nbas,[-999], nev, isp,1,qp,evec, aus )
!    lmxax = ll(nlmax)
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
                   if (m < 0) then
                      auasaz=  iot*1d0/dsqrt(2d0)    *aus(lc-em,iv,:,ksp,ib) + 1d0/dsqrt(2d0)   *aus(lc+em,iv,:,ksp,ib)
                   elseif (m > 0) then
                      auasaz= -iot*(-1)**m/dsqrt(2d0)*aus(lc-m,iv,:,ksp,ib)  +(-1)**m/dsqrt(2d0)*aus(lc+m,iv,:,ksp,ib)
                   else
                      auasaz= aus(ilm,iv,:,ksp,ib)
                   endif ! (au as az) are for (u,s,gz) functions where gz=gz'=0 at MT
                   orbtm(l+1,ksp,ib)= orbtm(l+1,ksp,ib) +m*wtkb(iv,isp,iq)*sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,ksp,ib),auasaz))
                enddo mloop
             enddo lloop
          enddo ivloop ! print*, l, ksp,ib,'ORB.MOMNT=',(orbtm(l+1,ksp,ib),l=0,lmxa)
       enddo ispcloop
    enddo ibloop
    deallocate(aus)
  end subroutine mkorbm
  subroutine mkdmtu(isp,iq,qp,nev,evec,dmatu) !Get density matrix dmatu for LDA+U (phi-projected density matrix)
    use m_lmfinit,only: ispec,sspec=>v_sspec,nbas,nlmax,nsp,nspc,n0,nppn,nlibu,lmaxu,nlibu,lldau,idu
    use m_mkpot,only: phzdphz
    use m_subzi, only: wtkb=>rv_a_owtkb
    use m_igv2x,only: ndimh
    use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx
    use m_makusq,only: makusq
    use m_locpot,only: rotp
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   sspec : lmxa idu
    !i   wtkb  :eigenvalue weights for BZ integration of occupied states
    !i   isp   :current spin channel (1 or 2)
    !i   iq    :qp index, used only to address element in wtkb
    !i         :NB: aus is stored only for current qp
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   ndham :dimensions wtkb,aus
    !i   nlmax :1st dimension of aus (maximum nlma over all sites)
    !i   nbas  :size of basis
    !i   nev   :actual number of eigenvectors generated
    !i   phzdphz  :phz dphz
    !i   aus   :coefficients to phi and phidot made previously by makusqldau
    !i  lldau  :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat beginning at dmats(*,lldau(ib))
    !o Outputs
    !o   dmatu :density matrix for specified LDA+U channels
    !b Bugs
    !b   Never checked for noncollinear case
    !u Updates
    !u   09 Nov 05 Convert dmat to complex form
    !u   28 Jun 05 bug fix for nspc=2
    !u   09 Jun 05 (MvS) extended to local orbitals
    !u   30 Apr 05 (WRL) first created
    ! ----------------------------------------------------------------------
    implicit none
    integer :: isp,iq,nev
    double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    double complex add,au,as,az,ap1,ap2
    double precision :: dlphi,rmt,dlphip,phi,phip,dphi,dphip,r(2,2),det,phz,dphz
    integer :: lmxa,ilm1,ilm2,l,iv,m1,m2,ib,is,igetss,iblu,ispc, ksp
    complex(8),allocatable ::aus(:,:,:,:,:)
    double complex evec(ndimh,nsp,nev)
    real(8)::qp(3)
    complex(8):: auas(2)
    allocate(aus(nlmax,ndham*nspc,3,nsp,nbas))! aus(2*nlmax*ndhamx*3*nsp*nbas))
    aus=0d0
    call makusq(nbas,[0] , nev,  isp, 1, qp, evec, aus )
    iblu = 0
    do  ib = 1, nbas
       if (lldau(ib) == 0) goto 10
       is = ispec(ib)
       lmxa=sspec(is)%lmxa
       rmt= sspec(is)%rmt
       do  l = 0, min(lmxa,3)
          if (idu(l+1,is) /= 0) then
             iblu = iblu+1
             !           In noncollinear case, isp=1 always => need internal ispc=1..2
             !           ksp is the current spin index in both cases:
             !           ksp = isp  in the collinear case
             !               = ispc in the noncollinear case
             !           ispc=1 for independent spins, and spin index when nspc=2
             do  ispc = 1, nspc
                ksp = max(ispc,isp)
                phz    = phzdphz(1,l+1,ksp,ib)
                dphz   = phzdphz(2,l+1,ksp,ib)
                ilm1 = l*l
                do  m1 = -l, l
                   ilm1 = ilm1+1
                   ilm2 = l*l
                   do  m2 = -l, l
                      ilm2 = ilm2+1
                      add = (0d0,0d0)
                      !  Since (au,as,az) are coefficients to (u,s,gz), (gz is local orbital with val=slo=0 at MT)
                      !  Local orbital contribution adds to u,s
                      !  deltau = -phi(rmax) * az   deltas = -dphi(rmax) * az
                      do  iv = 1, nev
                         az = aus(ilm1,iv,3,ksp,ib)
                         au = aus(ilm1,iv,1,ksp,ib) - phz*az
                         as = aus(ilm1,iv,2,ksp,ib) - dphz*az !u,s components
                         auas= matmul([au,as],rotp(l,ksp,:,:,ib)) ! rotation (u,s) to (phi,phidot) comp.
                         ap1 = auas(1) !au*r(1,1) + as*r(2,1) !projection to phi components.
                         az = aus(ilm2,iv,3,ksp,ib)
                         au = aus(ilm2,iv,1,ksp,ib) - phz*az
                         as = aus(ilm2,iv,2,ksp,ib) - dphz*az
                         auas= matmul([au,as],rotp(l,ksp,:,:,ib))
                         ap2 = auas(1) !au*r(1,1) + as*r(2,1)
                         add = add + ap1*dconjg(ap2)*wtkb(iv,isp,iq)
                      enddo
                      dmatu(m1,m2,ksp,iblu) = dmatu(m1,m2,ksp,iblu) + add !dmatu is for phi-projected density matrix
                   enddo
                enddo
             enddo
          endif
       enddo
10     continue
    enddo
    deallocate(aus)
  end subroutine mkdmtu
end module m_bandcal
subroutine dfqkkl( oqkkl ) !Allocates arrays to accumulate output site density
  use m_lmfinit,only: nsp,nbas,ispec,sspec=>v_sspec,nkaphh
  use m_struc_def,only:s_rv5   !o oqkkl : memory is allocated for qkkl
  implicit none
  type(s_rv5) :: oqkkl(3,nbas)
  integer :: ib,is,kmax,lmxa,lmxh,nlma,nlmh ,nkaph
  do  ib = 1, nbas
     is = ispec(ib) 
     nkaph=nkaphh(is)
     lmxa=sspec(is)%lmxa
     if (lmxa == -1) cycle
     nlma = (lmxa+1)**2
     nlmh = (sspec(is)%lmxb+1)**2
     kmax =  sspec(is)%kmxt
     if(allocated(oqkkl(1,ib)%v)) deallocate(oqkkl(1,ib)%v,oqkkl(2,ib)%v,oqkkl(3,ib)%v)
     allocate(oqkkl(1,ib)%v(0:kmax,0:kmax, nlma,nlma ,nsp), source=0d0)! Pkl*Pkl
     allocate(oqkkl(2,ib)%v(nkaph, 0:kmax, nlmh,nlma ,nsp), source=0d0)! Pkl*Hsm
     allocate(oqkkl(3,ib)%v(nkaph,  nkaph, nlmh,nlmh ,nsp), source=0d0)! Hsm*Hsm
  enddo
end subroutine dfqkkl
