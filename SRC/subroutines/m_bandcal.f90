!>band structure calculation !How to learng this? Instead of reading all sources, understand I/O. See bndfp.f90 and 'use m_bandcal'.
module m_bandcal 
  use m_ftox
  use m_lgunit,only:stdo,stdl
  use m_lmfinit,only: lmxa_i=>lmxa,rmt_i=>rmt,afsym,nspx
  use m_struc_def,only: s_rv1,s_rv2,s_rv5
  use m_qplist, only: nkp
  use m_mkqp,only: ntet=> bz_ntet, bz_nabc
  use m_qplist,only: qplist,niqisp,iqproc,isproc
  use m_igv2x,only: m_igv2x_setiq, napw,ndimh,ndimhx,igv2x,nbandmx
  use m_lmfinit,only: lrsig=>ham_lsig, lso,ham_scaledsigma,lmet=>bz_lmet,nbas,epsovl=>ham_oveps,nspc,plbnd,lfrce
  use m_lmfinit,only: pwmode=>ham_pwmode,pwemax,nsp,nlibu,lmaxu,lmxax
  use m_MPItk,only: master_mpi, procid,strprocid, numprocs=>nsize, comm
  use m_subzi, only: nevmx
  use m_supot, only: n1,n2,n3
  use m_rdsigm2,only: senex,sene,getsenex,dsene,ndimsig
  use m_procar,only: m_procar_add,m_procar_closeprocar
  use m_clsmode,only: m_clsmode_set1
  use m_addrbl,only: addrbl
  use m_augmbl,only: aughsoc
  use m_makusq,only: makusq
  use m_zhev,only: zhev_tk4, zhev_gpu_cleanup
  use m_hambl,only : hambl
  use m_mkpot,only : m_mkpot_init,  osmpot,vconst                !main inputs for potential
  use m_locpot,only:                                   osig,otau,oppi,ohsozz,ohsopm !main inputs
  use m_lmfinit,only: ispec,nkaphh,kmxt_i=>kmxt,lmxb_i=>lmxb
  use m_lmfinit,only: nlmax,nspc,n0,lldau,idu
  use m_struc_def,only:s_rv5   !o oqkkl : memory is allocated for qkkl
  use m_mpiio, only: writem_d, openm, closem, openedm, writem_struct, record_item, record_item_from
  ! outputs ---------------------------
  public m_bandcal_init, m_bandcal_2nd, m_bandcal_clean, m_bandcal_allreduce, m_bandcal_symsmrho
  public :: m_bandcal_gather_evlall, m_bandcal_gather_spinweightall
  integer,allocatable,protected,public::     ndimhx_(:,:),nevls(:,:) 
  real(8),allocatable,protected,public::     frcband(:,:), orbtm_rv(:,:,:),evlall(:,:,:), spinweightall(:,:,:) !all data is accumulated in mpi master 
  complex(8),allocatable,protected,public::  smrho_out(:,:,:,:),dmatu(:,:,:,:)
  type(s_rv5),allocatable,protected,public:: oeqkkl(:,:), oqkkl(:,:)
  type(s_rv1), allocatable, protected, public :: t_evl(:,:)
  type(s_rv2), allocatable, protected, public :: t_spinweight(:)
  !------------------------------------------------
  logical,private:: debug,sigmamode,call_m_bandcal_2nd,procaron,dmatuinit=.true.!,writeham
  real(8),private:: sumqv(3,2),sumev(3,2)
  integer,allocatable,private::neviqis(:),ndimhxiqis(:)
  complex(8),allocatable,private:: eveciqis(:,:,:)
  private
contains
  subroutine m_bandcal_init(lrout,ef0,vmag,writeham) ! Set up Hamiltonian, diagonalization
    use m_gpu, only: use_gpu, ngpu_ranks
    implicit none
    intent(in)::            lrout,ef0,vmag,writeham
    complex(8),allocatable:: hamm(:,:,:,:),ovlm(:,:,:,:),hammhso(:,:,:),ovlms(:,:,:,:) !Hamiltonian,Overlapmatrix
    integer:: iq,nmx,ispinit,isp,nev,ifih,lwtkb,lrout,ifig,i,ibas,iwsene,idat,ikp,istat
    real(8):: qp(3),ef0,def=0d0,xv(3),q(3),vmag
    real(8),allocatable    :: evl(:,:), spinweight(:,:)
    complex(8),allocatable :: evec(:,:) !eigenvector( :,nband)
    logical:: ltet,cmdopt0,dmatuinit=.true.,wsene,magexist,writeham
    character(3):: charnum3
#ifdef __GPU
    ! Batch storage for GPU diagonalization
    complex(8),allocatable :: hamm_batch(:,:,:), ovlm_batch(:,:,:), ovlm_save(:,:,:,:,:)
    integer,allocatable :: ndimhx_batch(:), isp_batch(:), iq_batch(:), nmx_batch(:)
#endif
    call tcn('m_bandcal_init')
    if(master_mpi) write(stdo,ftox)'m_bandcal_init: start'
    sigmamode = mod(lrsig,10)/=0
    ! writeham = cmdopt0('--writeham')
    PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
    debug    = cmdopt0('--debugbndfp')
    ltet = ntet>0 !   nspx=nsp/nspc !nspc=1 only for so=1 
    if(plbnd==0 .AND. lso/=0 .AND. lmet==0 ) call rx('metal weights required to get orb.moment')
    if(lso/=0) allocate(orbtm_rv(lmxax+1,nsp,nbas),source=0d0) !for spin-orbit coupling
    if(lfrce>0) allocate( frcband(3,1:nbas),source=0d0) !force for band
    if(master_mpi) write(stdo,"('MagField added to Hailtonian -vmag/2 for isp=1, +vmag/2 for isp=2: vmag(Ry)=',d13.6)") vmag
    magexist= abs(vmag)>1d-6
    allocate( ndimhx_(nkp,nspx),nevls(nkp,nspx),source=0) 
    allocate( t_evl(nspx,nkp))
    if(lso==1) allocate(t_spinweight(nkp))
    if(master_mpi) allocate( evlall(nbandmx,nspx,nkp),source=0d0)
    if(master_mpi .and. lso==1) allocate( spinweightall(nbandmx,nsp,nkp),source=0d0) !nsp=2 for lso==1
    if(nlibu>0 .AND. dmatuinit) then
       allocate( dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
       dmatuinit=.false.
    endif
    if(nlibu>0)  dmatu=0d0    !density matrix for U initialization
    if(lrout/=0) then
       allocate( oeqkkl(3,nbas), oqkkl(3,nbas)) !pointer arrays. energy and charge for each ibas. 3 is for three radial channel
       call dfqkkl( oqkkl  )! allocate and zero clear
       call dfqkkl( oeqkkl )!
       allocate( smrho_out(n1,n2,n3,nsp),source=(0d0,0d0) )
    endif
    call_m_bandcal_2nd =.false.
    if(plbnd==0) call_m_bandcal_2nd= (lmet>=0 .AND. lrout>0 )
    if(call_m_bandcal_2nd) then 
       if(allocated(neviqis))deallocate(neviqis,ndimhxiqis,eveciqis)
       allocate(neviqis(niqisp),ndimhxiqis(niqisp),eveciqis(nbandmx,nevmx,niqisp)) 
    endif
    allocate( evl(nbandmx,nspx), spinweight(nbandmx,nsp))
    sumev = 0d0
    sumqv = 0d0
    if(writeham) then
      PrepWriteHamiltonianPMT:block
      integer :: ifihh_info, mrech
        mrech = 8*3+4+16*nbandmx*nbandmx*2
        if(master_mpi) then
          open(newunit=ifihh_info, file='__HamiltonianPMT.info', form='unformatted')
          write(ifihh_info) nbandmx, mrech
          close(ifihh_info)
        endif
        istat = openm(newunit=ifih,file='__HamiltonianPMT',recl=mrech, comm=comm)
        write(stdo,ftox) 'xxxx',nbandmx, mrech, ifih
      endblock PrepWriteHamiltonianPMT
    endif
#ifdef __GPU
    if(use_gpu) then
    gpumem_before: block
      use cudafor
      integer(8) :: free_mem, total_mem
      integer :: ierr_mem
      ierr_mem = cudaMemGetInfo(free_mem, total_mem)
      write(6,'(a,2f10.1,a)') ' GPU mem before k-loop: free/total(MB)=', &
           free_mem/1d6, total_mem/1d6, ' MB'
    endblock gpumem_before
    endif
#endif
    bandcalculation_q: do 2010 idat=1,niqisp
       iq = iqproc(idat)
       qp = qplist(:,iq) !write(stdo,ftox)'m_bandcal_init: procid iq=',procid,iq,ftof(qp)
       isp= isproc(idat) !NOTE: isp=1:nspx=nsp/nspc
       if((.not.writeham).and.(afsym.and.isp==2)) cycle 
       call m_Igv2x_setiq(iq) ! NOTE: we get napw,ndimh,ndimhx, and igv2x here.
       allocate(hamm(ndimh,nspc,ndimh,nspc),ovlm(ndimh,nspc,ndimh,nspc)) !Spin-offdiagonal block included since nspc=2 for lso=1.
       !!  hambl calls augmbl. See Appendix C in http://dx.doi.org/10.7566/JPSJ.84.034702
       !! We finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
       !! Note that F2=Hankel head at a site + Hankel tail contributions from the other site.
       Setup_hamiltonian_and_diagonalize : block !write(6,*) 'goto Setup_hamiltonian_and_diagonalize : block lso=',lso
         integer:: iprint,ispc
         character:: charnum3
         !! == Set up Hamiltonian by hambl. ==============
         !!    Hamiltonian: hamm(1:ndimh,1:ndimh,3) means off-diagonal section when SO=1.
         !!    Overlap matrix: ovlm
         !!    senex:  Sigma-Vxc
         !! ==========================================
         !! Generate senex=(Sigma-Vxc) for given sfz.
         !! Determine interpolated self-energy senex at qp from sfz.
         !! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only)
         !!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
         !! Main input  => ham_iv_a_oiaxs,ham_rv_a_ohrs
         !     ! Main output sene. See m_seneinterp

         !! See Eq.(36) and appendix in http://dx.doi.org/10.7566/JPSJ.84.034702
         !! Hamm and ovlm are made from smooth part and augmentation part.
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
            do ispc=1,2  ! nspc==2
               call hambl(ispc,qp,osmpot,vconst,osig,otau,oppi, hamm(:,ispc,:,ispc),ovlm(:,ispc,:,ispc))
               hamm(:,ispc,:,ispc)= hamm(:,ispc,:,ispc) + hammhso(:,:,ispc) !spin-diag SOC elements (1,1), (2,2) added
            enddo
            if (cmdopt0('--testso')) then !this is for AHC test
               hamm(:,1,:,2) = 0d0
               hamm(:,2,:,1) = 0d0
            else
              hamm(:,1,:,2)= hammhso(:,:,3)                    !spin-offdiagonal SOC elements (1,2) added
              hamm(:,2,:,1)= transpose(dconjg(hammhso(:,:,3)))
            endif   
            if(sigmamode) then
               do ispc=1,nspc
                  call getsenex(qp, ispc, ndimh, ovlm(:,ispc,:,ispc)) !bugfix at 2024-4-24 obata: ispc was 1 when 2023-9-20
                  hamm(:,ispc,:,ispc) = hamm(:,ispc,:,ispc) + ham_scaledsigma*senex !sene= Vxc(QSGW)-Vxc(LDA)
                  if(wsene) write(iwsene) qp,ispc
                  if(wsene) write(iwsene) sene
                  call dsene()
               enddo
            endif
            allocate(ovlms,source=ovlm)
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
         if(iprint()>=30) write(stdo,'(" bndfp: kpt ",i5," of ",i7, " k=",3f8.4, &
              " ndimh = nmto+napw = ",3i5,f13.5)') iq,nkp,qp,ndimh,ndimh-napw,napw
         if(writeham) then
           WriteHamiltonianPMT: block
              type(record_item), allocatable :: items(:)
              integer :: iqqisp
              complex(8) :: ovlm_(nbandmx, nbandmx), hamm_(nbandmx, nbandmx)
              ovlm_(1:ndimhx,1:ndimhx) = reshape(ovlm, shape=[ndimhx,ndimhx])
              hamm_(1:ndimhx,1:ndimhx) = reshape(hamm, shape=[ndimhx,ndimhx])
              iqqisp= isp + nspx*(iq-1)
              items = [record_item_from(qp), record_item_from(ndimhx), record_item_from(ovlm_), record_item_from(hamm_)]
              istat = writem_struct(ifih, rec=iqqisp, items=items)
            ! write(ifih) qp,ndimhx,lso,epsovl,isp ! ndimhx=ndimh*nspc 
            ! write(ifih) ovlm ! When you read, use ovlm(1:ndimhx, 1:ndimhx)
            ! write(ifih) hamm
          endblock WriteHamiltonianPMT
         endif
         nmx=min(nevmx,ndimhx)! nmx:maximum number of eigenfunctions we will obtain. Smaller is faster.
         if(magexist) then
            if(nspc==2) then
               hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
               hamm(:,2,:,2)= hamm(:,2,:,2) + vmag/2d0*ovlm(:,2,:,2)
            else
               if(isp==1) hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
               if(isp==2) hamm(:,1,:,1)= hamm(:,1,:,1) + vmag/2d0*ovlm(:,1,:,1)
            endif
         endif
#ifdef __GPU
         ! --- GPU batched diag: store H, S for batch diagonalization after k-loop ---
         StoreForBatch: block
           complex(8), allocatable :: hamm_flat(:,:), ovlm_flat(:,:)
           if(.not.allocated(hamm_batch)) then
             allocate(hamm_batch(nbandmx,nbandmx,niqisp), ovlm_batch(nbandmx,nbandmx,niqisp))
             allocate(ndimhx_batch(niqisp), isp_batch(niqisp), iq_batch(niqisp), nmx_batch(niqisp))
             if(lso==1) allocate(ovlm_save(ndimh,nspc,ndimh,nspc,niqisp))
           endif
           allocate(hamm_flat(ndimhx,ndimhx), ovlm_flat(ndimhx,ndimhx))
           hamm_flat = reshape(hamm, shape=[ndimhx,ndimhx])
           ovlm_flat = reshape(ovlm, shape=[ndimhx,ndimhx])
           hamm_batch(1:ndimhx,1:ndimhx,idat) = hamm_flat
           ovlm_batch(1:ndimhx,1:ndimhx,idat) = ovlm_flat
           deallocate(hamm_flat, ovlm_flat)
           ndimhx_batch(idat) = ndimhx
           isp_batch(idat) = isp
           iq_batch(idat) = iq
           nmx_batch(idat) = nmx
           if(lso==1) ovlm_save(:,:,:,:,idat) = ovlm
         endblock StoreForBatch
       endblock Setup_hamiltonian_and_diagonalize
#else
         ! --- CPU path: diagonalize immediately ---
         allocate(evec(ndimhx,nmx))
         Diagonalize_hamilatonian: block
           call zhev_tk4(ndimhx, hamm, ovlm, nmx, nev, evl(1, isp ), evec, epsovl)
         endblock Diagonalize_hamilatonian
       endblock Setup_hamiltonian_and_diagonalize
       if(writeham.AND.master_mpi) write(stdo,"(9f8.4)") (evl(i,isp), i=1,nev)
       if(call_m_bandcal_2nd) then
          neviqis(idat)   =nev
          ndimhxiqis(idat)=ndimhx
          if(nmx/=0) eveciqis(1:ndimhx,1:nev,idat)=evec(1:ndimhx,1:nev)
       endif
       evl(nev+1:nbandmx,isp)=1d99
       nevls(iq,isp)  = nev
       ndimhx_(iq,isp)= ndimhx
       GetSpinWeightSOC1: if(lso==1.and.nmx/=0) then
          associate(nd=>ndimh)
            spinweight(:,:) = 0d0
            spinweight(1:nev,1)= [(sum(dconjg(evec(1:nd,i))*matmul(ovlms(:,1,:,1),evec(1:nd,i))),i=1,nev)]
            spinweight(1:nev,2)= [(sum(dconjg(evec(nd+1:nd+nd,i))*matmul(ovlms(:,2,:,2),evec(nd+1:nd+nd,i))),i=1,nev)]
            NormalizationcheckFORspinweightSOC: if(any([(abs(sum(spinweight(i,:))-1d0)>1d-6,i=1,nev-10)])) then
               do i=1,nev
                  write(stdo,ftox)'spinweightsoc=',i,ftof(spinweight(i,1:2)),sum(spinweight(i,1:2))
               enddo
               call rx('m_bandcal: error of SpinWeight sumcheck')
            endif NormalizationcheckFORspinweightSOC
          endassociate
       endif GetSpinWeightSOC1
       if(allocated(t_evl(isp,iq)%v)) deallocate(t_evl(isp,iq)%v)
       allocate(t_evl(isp,iq)%v(nbandmx), source = evl(:,isp))
       if(lso==1) allocate(t_spinweight(iq)%v(nbandmx,nsp), source = spinweight)
       if(afsym) then
          if(allocated(t_evl(2,iq)%v)) deallocate(t_evl(2,iq)%v)
          allocate(t_evl(2,iq)%v(nbandmx), source = evl(:,1))
          nevls(iq,2)  = nev
          ndimhx_(iq,2)= ndimhx
       endif
       if(master_mpi.AND.epsovl>=1d-14.AND.plbnd/=0) write(stdo,&
            "(' : ndimhx=',i5,' --> nev=',i5,' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
       if(PROCARon) call m_procar_add(iq,isp,ef0,evl,qp,nev,evec,ndimhx)
       if(allocated(evec)) deallocate(evec)
#endif
       if(allocated(hammhso)) deallocate(hammhso)
       if(allocated(hamm)) deallocate(hamm,ovlm)
       if(allocated(ovlms)) deallocate(ovlms)
2010 enddo bandcalculation_q
#ifdef __GPU
    ! === GPU batched diagonalization: all k-points on GPU ===
    if(niqisp > 0 .and. use_gpu) then
    BatchDiag: block
      use cusolverdn
      use cublas_v2
      use m_zhev_gpu_handles
      use cudafor
      complex(8), device, allocatable :: omat_d(:,:), h_d(:,:), zz_d(:,:), hhm_d(:,:), hh_d(:,:), z_d(:,:)
      real(8), device, allocatable :: eo_d(:), e_d(:)
      complex(8), device, allocatable :: work_d(:)
      integer, device, allocatable :: devinfo
      complex(8), allocatable :: zz_h(:,:), evec_tmp(:,:)
      real(8) :: eo(nbandmx)
      integer :: jdat, nd, nmx_j, nev_j, iq_j, isp_j, ix_j, ni_j, nm_j, lwork_j, m_out, istat_g
      ! Create handles once for all k-points
      if(.not. zhev_gpu_handles_init) then
        istat_g = cusolverDnCreate(zhev_cusolver_handle)
        istat_g = cublasCreate(zhev_cublas_handle)
        zhev_gpu_handles_init = .true.
      endif
      allocate(devinfo)
      do jdat = 1, niqisp
        nd = ndimhx_batch(jdat)
        nmx_j = nmx_batch(jdat)
        isp_j = isp_batch(jdat)
        iq_j = iq_batch(jdat)
        ! Step 1: Diag overlap on GPU
        allocate(omat_d(nd,nd), eo_d(nd))
        omat_d = ovlm_batch(1:nd,1:nd,jdat)
        istat_g = cusolverDnZheevdx_bufferSize(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
             CUSOLVER_EIG_RANGE_ALL, CUBLAS_FILL_MODE_UPPER, nd, omat_d, nd, &
             0d0, 0d0, 1, nd, m_out, eo_d, lwork_j)
        allocate(work_d(lwork_j))
        istat_g = cusolverDnZheevdx(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
             CUSOLVER_EIG_RANGE_ALL, CUBLAS_FILL_MODE_UPPER, nd, omat_d, nd, &
             0d0, 0d0, 1, nd, m_out, eo_d, work_d, lwork_j, devinfo)
        deallocate(work_d)
        eo(1:nd) = eo_d(1:nd)
        ! Step 2: Build projection zz (overlap reduction by epsovl)
        ni_j = 1
        do ix_j = 1, nd
          if(eo(ix_j) > epsovl) then; ni_j = ix_j; exit; endif
        enddo
        nm_j = nd - ni_j + 1
        allocate(zz_h(nd, nm_j))
        zz_h = omat_d(:, ni_j:nd)  ! copy eigenvectors to host
        do ix_j = ni_j, nd
          zz_h(:, ix_j-ni_j+1) = zz_h(:, ix_j-ni_j+1) / sqrt(eo(ix_j))
        enddo
        allocate(zz_d(nd, nm_j))
        zz_d = zz_h
        deallocate(zz_h, omat_d, eo_d)
        ! Step 3: Project H: hh = zz^H * H * zz
        allocate(h_d(nd,nd), hhm_d(nm_j,nd), hh_d(nm_j,nm_j))
        h_d = hamm_batch(1:nd,1:nd,jdat)
        istat_g = cublasZgemm_v2(zhev_cublas_handle, CUBLAS_OP_C, CUBLAS_OP_N, nm_j, nd, nd, &
             (1d0,0d0), zz_d, nd, h_d, nd, (0d0,0d0), hhm_d, nm_j)
        istat_g = cublasZgemm_v2(zhev_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, nm_j, nm_j, nd, &
             (1d0,0d0), hhm_d, nm_j, zz_d, nd, (0d0,0d0), hh_d, nm_j)
        deallocate(hhm_d, h_d)
        ! Step 4: Diag reduced H
        nev_j = min(nmx_j, nm_j)
        allocate(e_d(nm_j))
        istat_g = cusolverDnZheevdx_bufferSize(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
             CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, nm_j, hh_d, nm_j, &
             0d0, 0d0, 1, nev_j, m_out, e_d, lwork_j)
        allocate(work_d(lwork_j))
        istat_g = cusolverDnZheevdx(zhev_cusolver_handle, CUSOLVER_EIG_MODE_VECTOR, &
             CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, nm_j, hh_d, nm_j, &
             0d0, 0d0, 1, nev_j, m_out, e_d, work_d, lwork_j, devinfo)
        deallocate(work_d)
        evl(1:nev_j, isp_j) = e_d(1:nev_j)
        nev_j = m_out
        deallocate(e_d)
        ! Step 5: Back-transform z = zz * hh_d(:,1:nev)
        allocate(z_d(nd, nev_j))
        istat_g = cublasZgemm_v2(zhev_cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, nd, nev_j, nm_j, &
             (1d0,0d0), zz_d, nd, hh_d, nm_j, (0d0,0d0), z_d, nd)
        allocate(evec_tmp(nd, nev_j))
        evec_tmp = z_d
        deallocate(zz_d, hh_d, z_d)
        nev = nev_j
        ! Store results
        if(writeham.AND.master_mpi) write(stdo,"(9f8.4)") (evl(i,isp_j), i=1,nev)
        if(call_m_bandcal_2nd) then
          neviqis(jdat) = nev
          ndimhxiqis(jdat) = nd
          if(nmx_j/=0) eveciqis(1:nd,1:nev,jdat) = evec_tmp(1:nd,1:nev)
        endif
        evl(nev+1:nbandmx,isp_j) = 1d99
        nevls(iq_j,isp_j) = nev
        ndimhx_(iq_j,isp_j) = nd
        if(lso==1 .and. nmx_j/=0) then
          associate(ndd => nd/nspc)
            spinweight(:,:) = 0d0
            spinweight(1:nev,1) = [(sum(dconjg(evec_tmp(1:ndd,i))*matmul(ovlm_save(:,1,:,1,jdat),evec_tmp(1:ndd,i))),i=1,nev)]
            spinweight(1:nev,2) = [(sum(dconjg(evec_tmp(ndd+1:ndd+ndd,i))*matmul(ovlm_save(:,2,:,2,jdat),evec_tmp(ndd+1:ndd+ndd,i))),i=1,nev)]
          endassociate
        endif
        if(allocated(t_evl(isp_j,iq_j)%v)) deallocate(t_evl(isp_j,iq_j)%v)
        allocate(t_evl(isp_j,iq_j)%v(nbandmx), source = evl(:,isp_j))
        if(lso==1) allocate(t_spinweight(iq_j)%v(nbandmx,nsp), source = spinweight)
        if(afsym) then
          if(allocated(t_evl(2,iq_j)%v)) deallocate(t_evl(2,iq_j)%v)
          allocate(t_evl(2,iq_j)%v(nbandmx), source = evl(:,isp_j))
          nevls(iq_j,2) = nev
          ndimhx_(iq_j,2) = nd
        endif
        if(PROCARon) then
          allocate(evec(nd,nev))
          evec(1:nd,1:nev) = evec_tmp(1:nd,1:nev)
          call m_procar_add(iq_j,isp_j,ef0,evl,qplist(:,iq_j),nev,evec,nd)
          deallocate(evec)
        endif
        deallocate(evec_tmp)
      enddo
      deallocate(devinfo)
      deallocate(hamm_batch, ovlm_batch, ndimhx_batch, isp_batch, iq_batch, nmx_batch)
      if(allocated(ovlm_save)) deallocate(ovlm_save)
    endblock BatchDiag
    endif ! niqisp > 0 .and. use_gpu
    if(use_gpu) then
    gpumem_after: block
      use cudafor
      integer(8) :: free_mem, total_mem
      integer :: ierr_mem
      ierr_mem = cudaMemGetInfo(free_mem, total_mem)
      write(6,'(a,2f10.1,a)') ' GPU mem after k-loop: free/total(MB)=', &
           free_mem/1d6, total_mem/1d6, ' MB'
    endblock gpumem_after
    ! Release GPU handles to free memory before bandcal_2nd
    call zhev_gpu_cleanup()
    gpumem_freed: block
      use cudafor
      integer(8) :: free_mem, total_mem
      integer :: ierr_mem
      ierr_mem = cudaMemGetInfo(free_mem, total_mem)
      write(6,'(a,2f10.1,a)') ' GPU mem after cleanup: free/total(MB)=', &
           free_mem/1d6, total_mem/1d6, ' MB'
    endblock gpumem_freed
    endif ! use_gpu
    ! === Redistribute evec from GPU ranks to all ranks for bandcal_2nd ===
    if(ngpu_ranks > 0 .and. ngpu_ranks < numprocs .and. call_m_bandcal_2nd) then
    RedistributeEvec: block
      use mpi
      use m_dstrbp, only: dstrbp
      use m_gpu, only: ngpu_ranks
      use m_qplist, only: m_qplist_redistribute_all, kpproc, owner
      integer :: ndata_all, ierr_r, nspxx_r
      integer, allocatable :: kpproc_gpu(:), kpproc_new(:)
      integer :: niqisp_new, jdat, iq_r, isp_r, itag_r
      integer :: nev_r, nd_r, src_r, dest_r
      integer :: mpi_status(mpi_status_size)
      complex(8), allocatable :: eveciqis_old(:,:,:)
      integer, allocatable :: neviqis_old(:), ndimhxiqis_old(:)
      nspxx_r = nspx
      if((.not.cmdopt0('--jobgw')).and.(.not.cmdopt0('--writeham')).and.afsym) nspxx_r=1
      ndata_all = nkp * nspxx_r
      ! Broadcast nevls and ndimhx_ from GPU ranks to all ranks (needed for redistribution)
      call mpi_bcast(nevls, size(nevls), mpi_integer, 0, comm, ierr_r)
      call mpi_bcast(ndimhx_, size(ndimhx_), mpi_integer, 0, comm, ierr_r)
      ! Save GPU distribution
      allocate(kpproc_gpu(0:numprocs))
      kpproc_gpu = kpproc
      ! Save old eveciqis
      if(niqisp > 0 .and. allocated(eveciqis)) then
        allocate(eveciqis_old, source=eveciqis)
        allocate(neviqis_old, source=neviqis)
        allocate(ndimhxiqis_old, source=ndimhxiqis)
      endif
      ! Switch to standard all-rank distribution
      call m_qplist_redistribute_all()
      ! Compute new distribution for looking up new owners
      allocate(kpproc_new(0:numprocs))
      kpproc_new = kpproc
      niqisp_new = niqisp
      ! Reallocate eveciqis for new distribution
      if(allocated(neviqis)) deallocate(neviqis, ndimhxiqis, eveciqis)
      allocate(neviqis(max(1,niqisp_new)), ndimhxiqis(max(1,niqisp_new)))
      allocate(eveciqis(nbandmx, nevmx, max(1,niqisp_new)))
      ! Transfer: all ranks loop over all k-points; sender sends, receiver recvs
      do itag_r = 1, ndata_all
        iq_r = (itag_r - 1) / nspxx_r + 1
        isp_r = mod(itag_r - 1, nspxx_r) + 1
        nev_r = nevls(iq_r, isp_r)
        nd_r = ndimhx_(iq_r, isp_r)
        ! Old owner (GPU rank)
        src_r = -1
        do jdat = 0, numprocs - 1
          if(itag_r >= kpproc_gpu(jdat) .and. itag_r < kpproc_gpu(jdat+1)) then
            src_r = jdat; exit
          endif
        enddo
        ! New owner (all-rank distribution)
        dest_r = owner(isp_r, iq_r)  ! already updated by redistribute_all
        if(src_r < 0 .or. dest_r < 0) cycle
        if(src_r == dest_r .and. procid == src_r) then
          block
            integer :: old_idx, new_idx
            old_idx = itag_r - kpproc_gpu(procid) + 1
            new_idx = itag_r - kpproc_new(procid) + 1
            neviqis(new_idx) = nev_r
            ndimhxiqis(new_idx) = nd_r
            eveciqis(1:nd_r, 1:nev_r, new_idx) = eveciqis_old(1:nd_r, 1:nev_r, old_idx)
          endblock
        else
          if(procid == src_r) then
            block
              integer :: old_idx
              complex(8), allocatable :: sendbuf(:,:)
              old_idx = itag_r - kpproc_gpu(procid) + 1
              allocate(sendbuf(nd_r, nev_r))
              sendbuf = eveciqis_old(1:nd_r, 1:nev_r, old_idx)
              call mpi_send(sendbuf, nd_r*nev_r, mpi_double_complex, dest_r, itag_r, comm, ierr_r)
              deallocate(sendbuf)
            endblock
          endif
          if(procid == dest_r) then
            block
              integer :: new_idx
              complex(8), allocatable :: recvbuf(:,:)
              new_idx = itag_r - kpproc_new(procid) + 1
              neviqis(new_idx) = nev_r
              ndimhxiqis(new_idx) = nd_r
              allocate(recvbuf(nd_r, nev_r))
              call mpi_recv(recvbuf, nd_r*nev_r, mpi_double_complex, src_r, itag_r, comm, mpi_status, ierr_r)
              eveciqis(1:nd_r, 1:nev_r, new_idx) = recvbuf
              deallocate(recvbuf)
            endblock
          endif
        endif
      enddo
      if(allocated(eveciqis_old)) deallocate(eveciqis_old, neviqis_old, ndimhxiqis_old)
      deallocate(kpproc_gpu, kpproc_new)
    endblock RedistributeEvec
    endif
#endif
    if(writeham) istat = closem(ifih)
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       xv(:)=[(sum(frcband(i,1:nbas))/nbas,i=1,3)]
       forall(ibas= 1:nbas) frcband(:,ibas) = frcband(:,ibas) - xv(:) ! Average forces so net force on system is zero (APW case)
    endif
    if(PROCARon) call m_procar_closeprocar()
    if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid !if(call_m_bandcal_2nd) close(ifig)
    deallocate(evl)
    if(allocated(spinweight))deallocate(spinweight)
    call tcx('m_bandcal_init')
  end subroutine m_bandcal_init
  subroutine m_bandcal_2nd()! accumulate eval,evec-related quantities by addrbl
    implicit none
    integer:: iq,ispinit,isp,nev,ifig,i,ibas,idat
    real(8):: qp(3),def=0d0,xv(3)
    real(8),allocatable:: evl(:,:)
    complex(8),allocatable :: evec(:,:)!,evecbackup(:,:)
    logical:: cmdopt0
    call tcn('m_bandcal_2nd')
    if(master_mpi) write(stdo,ftox)'m_bandcal_2nd: to fill eigenfunctions**2 up to Efermi'
    call dfqkkl( oqkkl ) !zero clear
    call dfqkkl( oeqkkl ) !zero clear if(lekkl==1) 
    if (lfrce>0)  frcband  = 0d0
    if(lso/=0) orbtm_rv=0d0
    if(allocated(smrho_out)) deallocate(smrho_out)
    allocate( smrho_out(n1,n2,n3,nsp) )
    smrho_out = 0d0
    sumev = 0d0
    sumqv = 0d0
    allocate(evl(nbandmx,nspx))
    iqloop: do 12010 idat=1,niqisp !iq = iqini, iqend !This is a big iq loop
       iq = iqproc(idat)
       qp = qplist(:,iq)  !write(stdo,ftox)'m_bandcal_init: procid iq=',procid,iq,ftof(qp)
       isp= isproc(idat)
       if(afsym.and.isp==2) cycle !cmdopt0('--afsym').and.isp==2) cycle
       call m_Igv2x_setiq(iq) ! Get napw and so on for given qp
       nev   = neviqis(idat)
       !write(6,*)'nnnnnnnnnn iq nev=',iq,nev
       ndimhx= ndimhxiqis(idat)
       allocate(evec(ndimhx,nev))
       evl(1:nev,isp)= t_evl(isp,iq)%v(1:nev)
       evl(nev+1:nbandmx,isp)=1d99 !padding 
       evec(1:ndimhx,1:nev)=eveciqis(1:ndimhx,1:nev,idat)
       if(lso/=0)              call mkorbm(isp, nev, iq,qp, evec,  orbtm_rv)
       if(nlibu>0 .AND. nev>0) call mkdmtu(isp, iq,qp, nev, evec,  dmatu)
       if(cmdopt0('--cls'))    call m_clsmode_set1(nev,isp,iq,qp,nev,evec) !all inputs
       call addrbl(isp,qp,iq, osmpot,vconst,osig,otau,oppi,evec,evl,nev, smrho_out, sumqv, sumev, oqkkl,oeqkkl, frcband)
       afsymGETevecFROMisponeANDaccumulate:  if(afsym) then !cmdopt0('--afsym')) then
          if(idat==1.and.master_mpi) write(stdo,ftox)'m_bandcal: afsymblock'
          afsymblock: block !isp2 is given by isp=1
            use m_rotwave,only:  rotevec
            use m_mksym,only: symops,ngrp,ngrpAF,ag
            use m_qplist,only: qplist,nkp
            use m_lattic,only: plat=>lat_plat
            use m_ftox
            use m_subzi, only: m_subzi_copy_wtkb
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
              write(stdo,ftox)'igrp ikp',igrp,ikp,'qp=',ftof(qp,3),'qplist=',ftof(qplist(:,ikp),3) !,'deltaG=',ndeltaG
              write(stdo,ftox)'qp=',ftof(qp,3)
              do ikpx=1,nkp; write(stdo,ftox)'qplist=',ikpx,ftof(qplist(:,ikpx),3); enddo
              call rx('sygmrpAF afsym mode: can not find qtarget by afsym')
            endblock DebugWrite
1018        continue!write(stdo,ftox)'ikp qp=',ikp,ftof(qp,3),'is mapped to',ftof(matmul(symops(:,:,igrp),qp),3),' by symops igp=',igrp
            qpr = qplist(:,ikp)
            isp2 = 2
            call m_Igv2x_setiq(ikp) ! Get napw and so on for given qp !Bug fix at 2023-10-29. This set igv.
            call rotevec(igrp,qp, qpr,ndimhx,napw,nev,evec(:,1:nev), evecrot(:,1:nev))! evec at qp is roteted to be evecrot at qpr by symops(:,:,igrp)
            call m_subzi_copy_wtkb(isp, iq, isp2, ikp) !copy t_wtkb(isp,iq)%v to t_wtkb(isp2,ikp)%v
            if( lso/=0)              call mkorbm(isp2, nev, ikp,qpr, evecrot,  orbtm_rv)
            if( nlibu>0 .AND. nev>0) call mkdmtu(isp2,      ikp,qpr, nev, evecrot,  dmatu)
            if( cmdopt0('--cls'))    call m_clsmode_set1(nev,isp2,ikp,qpr,nev,evecrot) 
            call addrbl(isp2,qpr,ikp, osmpot,vconst,osig,otau,oppi,evecrot,evl,nev, smrho_out, sumqv, sumev, oqkkl,oeqkkl, frcband)
          endblock afsymblock
       endif afsymGETevecFROMisponeANDaccumulate
       deallocate(evec)
12010 enddo iqloop
    if (pwemax>0 .AND. mod(pwmode,10)>0 .AND. lfrce/=0) then
       xv(:)=[(sum(frcband(i,1:nbas))/nbas,i=1,3)]
       do  ibas= 1, nbas
          frcband(:,ibas) = frcband(:,ibas) - xv(:) ! Average forces so net force on system is zero (APW case)
       enddo
    endif
    deallocate(evl)
    call tcx('m_bandcal_2nd')
  end subroutine m_bandcal_2nd
  subroutine m_bandcal_gather_evlall()
    use m_qplist, only: owner
    use m_MPItk, only: master, master_mpi, procid, comm
    use mpi, only: mpi_double_precision, mpi_status_size
    implicit none
    integer:: status(MPI_Status_size), iq, isp, itag, ierr
    do iq = 1, nkp
      do isp = 1, nspx
        itag = (iq-1)*nspx + isp
        if(master_mpi) then
          if(owner(isp,iq) == master) evlall(:,isp,iq) = t_evl(isp,iq)%v
          if(owner(isp,iq) /= master) call mpi_recv(evlall(:,isp,iq), nbandmx, mpi_double_precision, &
                                                    owner(isp,iq), itag, comm, status, ierr)
        else
          if(owner(isp,iq) == procid) call mpi_send(t_evl(isp,iq)%v, nbandmx, mpi_double_precision, &
                                                    master, itag, comm, ierr)
        endif
      enddo
    enddo
  end subroutine
  subroutine m_bandcal_gather_spinweightall()
    use m_qplist, only: owner
    use m_MPItk, only: master, master_mpi, procid, comm
    use mpi, only: mpi_double_precision, mpi_status_size
    implicit none
    integer:: status(MPI_Status_size), iq, isp, itag, ierr
    do iq = 1, nkp
      do isp = 1, nspx !nspx should be 1
        itag = (iq-1)*nspx + isp
        if(master_mpi) then
          if(owner(isp,iq) == master) spinweightall(:,:,iq) = t_spinweight(iq)%v
          if(owner(isp,iq) /= master) call mpi_recv(spinweightall(:,:,iq), nbandmx*2, mpi_double_precision, &
                                                    owner(isp,iq), itag, comm, status, ierr)
        else
          if(owner(isp,iq) == procid) call mpi_send(t_spinweight(iq)%v, nbandmx*2, mpi_double_precision, &
                                                    master, itag, comm, ierr)
        endif
      enddo
    enddo
  end subroutine
  subroutine m_bandcal_clean() !cleaning allocation
    if(allocated(orbtm_rv)) deallocate(orbtm_rv)
    if(allocated(smrho_out)) deallocate(smrho_out)
    if(allocated(frcband))  deallocate(frcband)
    if(allocated(ndimhx_))  deallocate(ndimhx_,nevls)
    if(allocated(evlall)) deallocate(evlall)
    if(allocated(spinweightall)) deallocate(spinweightall)
    if(allocated(oqkkl)) deallocate( oqkkl)
    if(allocated(oeqkkl))deallocate( oeqkkl)
    if(allocated(t_evl))  deallocate(t_evl)
    if(allocated(t_spinweight)) deallocate(t_spinweight)
  end subroutine m_bandcal_clean
  subroutine m_bandcal_allreduce()!  Allreduce density-related quantities
    integer:: nnn,ib,i
    if(debug) print *,'goto m_bandcal_allreduce'
    call mpibc2_real(sumqv,size(sumqv),'bndfp_sumqv')
    call mpibc2_real(sumev,size(sumev),'bndfp_sumev')
    call mpibc2_complex(smrho_out,size(smrho_out),'bndfp_smrho')
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
    use m_igv2x,only: napw,ndimh,ndimhx,igvapw=>igv2x
    use m_locpot,only: sab_rv=>sab
    use m_subzi, only: t_wtkb
    use m_qplist,only: nkp
    !i   isp   :current spin channel (1 or 2)
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 for so=1 (SOC), 1 otherwise.
    !i   nlmax :leading dimension of aus
    !i   nev   :number of eigenvectors to accumulate orbital moment
    !i   iq    :current k-point
    !i   aus   :values of (phi,phidot,pz) MT sphere boundary; see makusq
    !i   nkp   :number of irreducible k-points
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
    ! ----------------------------------------------------------------------
    implicit none
    integer :: isp,nev,iq,ispx
    integer :: lmxa,lmdim,ichan,ib,is,igetss,iv,ilm,l,m,nlma, lc,em,ispc,ksp
    real(8):: qp(3),diff
    real(8):: suml(11),s11,s22,s12,s33,s31,s32,s13,s23, suma,rmt,orbtm(lmxax+1,nsp,*) 
    complex(8):: au,as,az,iot=(0d0,1d0),evec(ndimh,nsp,nev),auasaz(3)
    complex(8),allocatable ::aus(:,:,:,:,:)
    allocate(aus(nlmax,nbandmx,3,nsp,nbas))
    call makusq(nbas,[-999], nev, isp,1,qp,evec, aus )
    ichan = 0
    ibloop: do  ib = 1, nbas
       is = ispec(ib)
       lmxa=lmxa_i(is)
       rmt= rmt_i(is)
       lmxa = min(lmxa,lmxax)
       if (lmxa == -1) cycle 
       nlma = (lmxa+1)**2
       lmdim = nlma
       !nspc=2 if two spins are coupled(lso=1), nspc=1 otherwize 
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
                   ispx=merge(1,isp,lso==1)
                   ! orbtm(l+1,ksp,ib)= orbtm(l+1,ksp,ib) +m*wtkb(iv,ispx,iq)&
                   orbtm(l+1,ksp,ib)= orbtm(l+1,ksp,ib) +m*t_wtkb(ispx,iq)%v(iv)&
                        *sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,ksp,ib),auasaz))
                enddo mloop
             enddo lloop
          enddo ivloop ! print*, l, ksp,ib,'ORB.MOMNT=',(orbtm(l+1,ksp,ib),l=0,lmxa)
       enddo ispcloop
    enddo ibloop
    deallocate(aus)
  end subroutine mkorbm
  subroutine mkdmtu(isp,iq,qp,nev,evec,dmatu) !Get density matrix dmatu for LDA+U (phi-projected density matrix)
    use m_locpot,only: phzdphz
    use m_subzi, only: t_wtkb
    use m_igv2x,only: ndimh
    use m_makusq,only: makusq
    use m_locpot,only: rotp
    !i   wtkb  :eigenvalue weights for BZ integration of occupied states
    !i   isp   :current spin channel (1 or 2)
    !i   iq    :qp index, used only to address element in wtkb
    !i         :NB: aus is stored only for current qp
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   nlmax :1st dimension of aus (maximum nlma over all sites)
    !i   nbas  :size of basis
    !i   nev   :actual number of eigenvectors generated
    !i   phzdphz  :phz dphz
    !i   aus   :coefficients to phi and phidot made previously by makusqldau
    !i  lldau  :lldau(ib)=0 => no U on this site otherwise
    !i         :U on site ib with dmat beginning at dmats(*,lldau(ib))
    !o Outputs
    !o   dmatu :density matrix for specified LDA+U channels
    implicit none
    integer :: isp,iq,nev
    double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
    double complex add,au,as,az,ap1,ap2
    double precision :: dlphi,rmt,dlphip,phi,phip,dphi,dphip,r(2,2),det,phz,dphz
    integer :: lmxa,ilm1,ilm2,l,iv,m1,m2,ib,is,igetss,iblu,ispc, ksp
    complex(8) ::aus(nlmax,nbandmx,3,nsp,nbas), evec(ndimh,nsp,nev)
    real(8)::qp(3)
    complex(8):: auas(2)
    call makusq(nbas,[0] , nev,  isp, 1, qp, evec, aus )
    iblu = 0
    do  ib = 1, nbas
       if(lldau(ib) == 0) cycle
       is  = ispec(ib)
       lmxa=lmxa_i(is) 
       rmt = rmt_i(is)
       do l = 0, min(lmxa,3)
          if (idu(l+1,is) ==0) cycle
          iblu = iblu+1
          do  ispc = 1, nspc 
             ksp = max(ispc,isp) !! ksp is the current spin index in both cases:  ksp = isp  in the collinear case, = ispc in the noncollinear case
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
                      ! add = add + ap1*dconjg(ap2)*wtkb(iv,isp,iq)
                      add = add + ap1*dconjg(ap2)*t_wtkb(isp,iq)%v(iv)
                   enddo
                   dmatu(m1,m2,ksp,iblu) = dmatu(m1,m2,ksp,iblu) + add !dmatu is for phi-projected density matrix
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine mkdmtu
  subroutine dfqkkl( oqkkl ) !Allocates arrays to accumulate output site density
    implicit none
    type(s_rv5) :: oqkkl(3,nbas)
    integer :: ib,is,kmax,lmxa,lmxh,nlma,nlmh ,nkaph
    do  ib = 1, nbas
      is = ispec(ib) 
      nkaph=nkaphh(is)
      lmxa=lmxa_i(is)
      if (lmxa == -1) cycle
      nlma = (lmxa+1)**2
      nlmh = (lmxb_i(is)+1)**2
      kmax =  kmxt_i(is)
      if(allocated(oqkkl(1,ib)%v)) deallocate(oqkkl(1,ib)%v,oqkkl(2,ib)%v,oqkkl(3,ib)%v)
      allocate(oqkkl(1,ib)%v(0:kmax,0:kmax, nlma,nlma ,nsp), source=0d0)! Pkl*Pkl
      allocate(oqkkl(2,ib)%v(nkaph, 0:kmax, nlmh,nlma ,nsp), source=0d0)! Pkl*Hsm
      allocate(oqkkl(3,ib)%v(nkaph,  nkaph, nlmh,nlmh ,nsp), source=0d0)! Hsm*Hsm
    enddo
  end subroutine dfqkkl
end module m_bandcal
