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
#ifdef __GPU
    use m_gpu, only: use_gpu
#endif
    implicit none
    intent(in)::            lrout,ef0,vmag,writeham
    complex(8),allocatable:: hamm(:,:,:,:),ovlm(:,:,:,:),hammhso(:,:,:),ovlms(:,:,:,:) !Hamiltonian,Overlapmatrix
    integer:: iq,nmx,ispinit,isp,nev,ifih,lwtkb,lrout,ifig,i,ibas,iwsene,idat,ikp,istat
    real(8):: qp(3),ef0,def=0d0,xv(3),q(3),vmag
    real(8),allocatable    :: evl(:,:), spinweight(:,:)
    complex(8),allocatable :: evec(:,:) !eigenvector( :,nband)
    logical:: ltet,cmdopt0,dmatuinit=.true.,wsene,magexist,writeham
    character(3):: charnum3
    ! Batch storage for hambl results (all ranks) and GPU diagonalization
    complex(8),allocatable :: hamm_batch(:,:,:), ovlm_batch(:,:,:)
    integer,allocatable :: ndimhx_batch(:), isp_batch(:), iq_batch(:), nmx_batch(:)
#ifdef __GPU
    complex(8),allocatable :: ovlm_save(:,:,:,:,:)
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
         ! --- Store H, S for batch diagonalization after k-loop (all ranks) ---
         StoreForBatch: block
           complex(8), allocatable :: hamm_flat(:,:), ovlm_flat(:,:)
           if(.not.allocated(hamm_batch)) then
             allocate(hamm_batch(nbandmx,nbandmx,niqisp), ovlm_batch(nbandmx,nbandmx,niqisp))
             allocate(ndimhx_batch(niqisp), isp_batch(niqisp), iq_batch(niqisp), nmx_batch(niqisp))
#ifdef __GPU
             if(lso==1) allocate(ovlm_save(ndimh,nspc,ndimh,nspc,niqisp))
#endif
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
#ifdef __GPU
           if(lso==1) ovlm_save(:,:,:,:,idat) = ovlm
#endif
         endblock StoreForBatch
       endblock Setup_hamiltonian_and_diagonalize
       if(allocated(hammhso)) deallocate(hammhso)
       if(allocated(hamm)) deallocate(hamm,ovlm)
       if(allocated(ovlms)) deallocate(ovlms)
2010 enddo bandcalculation_q
#ifdef __GPU
    ! === Phase 2: Gather hamm/ovlm from CPU ranks to GPU ranks ===
    GatherToGPU: block
      use m_gpu, only: use_gpu, ngpu_ranks
      use mpi, only: MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_STATUS_SIZE
      integer :: src, dst_gpu, nk_recv, nk_total, jdat_ext, ierr_g, niqisp_local
      integer :: status_g(MPI_STATUS_SIZE)
      integer :: nk_all(0:numprocs-1)
      complex(8), allocatable :: hamm_all(:,:,:), ovlm_all(:,:,:)
      integer, allocatable :: ndimhx_all(:), isp_all(:), iq_all(:), nmx_all(:)
      ! Exchange k-point counts
      call mpi_allgather(niqisp, 1, MPI_INTEGER, nk_all, 1, MPI_INTEGER, comm, ierr_g)
      niqisp_local = niqisp
      if(use_gpu) then
        ! GPU rank: collect from assigned CPU ranks (round-robin: src%ngpu_ranks == procid)
        nk_total = 0
        do src = 0, numprocs-1
          if(mod(src, ngpu_ranks) == procid) nk_total = nk_total + nk_all(src)
        enddo
        allocate(hamm_all(nbandmx,nbandmx,nk_total), ovlm_all(nbandmx,nbandmx,nk_total))
        allocate(ndimhx_all(nk_total), isp_all(nk_total), iq_all(nk_total), nmx_all(nk_total))
        ! Copy my own data first
        jdat_ext = 0
        do jdat_ext = 1, niqisp_local
          hamm_all(:,:,jdat_ext) = hamm_batch(:,:,jdat_ext)
          ovlm_all(:,:,jdat_ext) = ovlm_batch(:,:,jdat_ext)
          ndimhx_all(jdat_ext) = ndimhx_batch(jdat_ext)
          isp_all(jdat_ext) = isp_batch(jdat_ext)
          iq_all(jdat_ext) = iq_batch(jdat_ext)
          nmx_all(jdat_ext) = nmx_batch(jdat_ext)
        enddo
        jdat_ext = niqisp_local
        ! Receive from CPU ranks in ascending order
        do src = 0, numprocs-1
          if(src == procid) cycle
          if(mod(src, ngpu_ranks) /= procid) cycle
          nk_recv = nk_all(src)
          if(nk_recv == 0) cycle
          call mpi_recv(hamm_all(1,1,jdat_ext+1), nbandmx*nbandmx*nk_recv, MPI_DOUBLE_COMPLEX, src, 100, comm, status_g, ierr_g)
          call mpi_recv(ovlm_all(1,1,jdat_ext+1), nbandmx*nbandmx*nk_recv, MPI_DOUBLE_COMPLEX, src, 101, comm, status_g, ierr_g)
          call mpi_recv(ndimhx_all(jdat_ext+1), nk_recv, MPI_INTEGER, src, 102, comm, status_g, ierr_g)
          call mpi_recv(isp_all(jdat_ext+1), nk_recv, MPI_INTEGER, src, 103, comm, status_g, ierr_g)
          call mpi_recv(iq_all(jdat_ext+1), nk_recv, MPI_INTEGER, src, 104, comm, status_g, ierr_g)
          call mpi_recv(nmx_all(jdat_ext+1), nk_recv, MPI_INTEGER, src, 105, comm, status_g, ierr_g)
          jdat_ext = jdat_ext + nk_recv
        enddo
        ! Replace batch with gathered data
        deallocate(hamm_batch, ovlm_batch, ndimhx_batch, isp_batch, iq_batch, nmx_batch)
        call move_alloc(hamm_all, hamm_batch)
        call move_alloc(ovlm_all, ovlm_batch)
        call move_alloc(ndimhx_all, ndimhx_batch)
        call move_alloc(isp_all, isp_batch)
        call move_alloc(iq_all, iq_batch)
        call move_alloc(nmx_all, nmx_batch)
      else
        ! CPU rank: send data to assigned GPU rank
        dst_gpu = mod(procid, ngpu_ranks)
        if(niqisp > 0) then
          call mpi_send(hamm_batch, nbandmx*nbandmx*niqisp, MPI_DOUBLE_COMPLEX, dst_gpu, 100, comm, ierr_g)
          call mpi_send(ovlm_batch, nbandmx*nbandmx*niqisp, MPI_DOUBLE_COMPLEX, dst_gpu, 101, comm, ierr_g)
          call mpi_send(ndimhx_batch, niqisp, MPI_INTEGER, dst_gpu, 102, comm, ierr_g)
          call mpi_send(isp_batch, niqisp, MPI_INTEGER, dst_gpu, 103, comm, ierr_g)
          call mpi_send(iq_batch, niqisp, MPI_INTEGER, dst_gpu, 104, comm, ierr_g)
          call mpi_send(nmx_batch, niqisp, MPI_INTEGER, dst_gpu, 105, comm, ierr_g)
        endif
        ! CPU ranks keep hamm_batch/ovlm_batch for CPU fallback diag in ScatterResults
      endif
    endblock GatherToGPU
    ! === Phase 3: Pad matrices to uniform size + GPU batched diagonalization ===
    if(use_gpu) then
    PadAndBatchDiag: block
      use cusolverdn
      use m_blas, only: zmm => zmm_d, m_op_C
      use cudafor
      ! GPU device arrays
      complex(8), device, allocatable :: omat_d(:,:), zz_d(:,:), h_d(:,:), hhm_d(:,:), hh_d(:,:)
      complex(8), device, allocatable :: evec_d(:,:), z_d(:,:), work_d(:)
      real(8), device, allocatable :: eo_d(:)
      integer, device, allocatable :: devinfo
      ! Host arrays
      complex(8), allocatable :: zz_h(:,:), evec_tmp(:,:)
      real(8) :: eo(nbandmx)
      integer :: jdat, nd, nd_max, nmx_j, nev_j, iq_j, isp_j, ix_j, ni_j, nm_j
      integer :: istat_g, ndiag_total, jp, lwork_j, m_out
      type(cusolverDnHandle) :: cs_h
      logical :: skip_ovldiag
      ndiag_total = size(iq_batch)
      skip_ovldiag = (epsovl < 1d-14)
      ! === Pad all matrices to nd_max ===
      nd_max = maxval(ndimhx_batch(1:ndiag_total))
      do jdat = 1, ndiag_total
        nd = ndimhx_batch(jdat)
        if(nd < nd_max) then
          hamm_batch(1:nd, nd+1:nd_max, jdat) = (0d0, 0d0)
          hamm_batch(nd+1:nd_max, 1:nd, jdat) = (0d0, 0d0)
          hamm_batch(nd+1:nd_max, nd+1:nd_max, jdat) = (0d0, 0d0)
          do jp = nd+1, nd_max; hamm_batch(jp, jp, jdat) = (9999d0, 0d0); enddo
          ovlm_batch(1:nd, nd+1:nd_max, jdat) = (0d0, 0d0)
          ovlm_batch(nd+1:nd_max, 1:nd, jdat) = (0d0, 0d0)
          ovlm_batch(nd+1:nd_max, nd+1:nd_max, jdat) = (0d0, 0d0)
          do jp = nd+1, nd_max; ovlm_batch(jp, jp, jdat) = (1d0, 0d0); enddo
        endif
      enddo
      if(master_mpi) write(stdo,'(a,i5,a,i5,a,l1)') &
           ' BatchDiag(Ozaki+cuSOLVER): ndiag=',ndiag_total,' nd_max=',nd_max,' skip_ovldiag=',skip_ovldiag
      ! Create cuSOLVER handle (for eigensolve only, GEMM uses Ozaki)
      istat_g = cusolverDnCreate(cs_h)
      allocate(devinfo)
      if(cmdopt0('--chefsi') .or. cmdopt0('--tridiag') .or. cmdopt0('--qdwh')) then
        ! === Alternative eigensolver path ===
        alt_diag: block
          use mpi, only: MPI_WTIME
          complex(8), device, allocatable :: evecs_dd(:,:)
          integer :: jd, nd_j, nmx_jj, isp_jj, iq_jj, info_cf
          real(8) :: t0cf
          t0cf = MPI_WTIME()
          do jd = 1, ndiag_total
            nd_j = ndimhx_batch(jd); nmx_jj = nmx_batch(jd)
            isp_jj = isp_batch(jd); iq_jj = iq_batch(jd)
            allocate(evecs_dd(nd_j, nmx_jj))
            if(cmdopt0('--qdwh')) then
              use_qdwh: block
                use m_qdwh_eigensolver, only: qdwh_batched_zhegv
                ! QDWH: process single k-point (batch=1 for now, full batch later)
                call qdwh_batched_zhegv(nd_j, nmx_jj, 1, hamm_batch(1,1,jd), ovlm_batch(1,1,jd), &
                     nbandmx, evl(1,isp_jj), evecs_dd, info_cf)
              endblock use_qdwh
            elseif(cmdopt0('--tridiag')) then
              use_tridiag: block
                use m_ozaki_tridiag, only: ozaki_tridiag_zhegv
                call ozaki_tridiag_zhegv(nd_j, nmx_jj, hamm_batch(1,1,jd), nbandmx, &
                     ovlm_batch(1,1,jd), evl(1,isp_jj), evecs_dd, epsovl, info_cf)
              endblock use_tridiag
            else
              use_chefsi: block
                use m_chefsi, only: chefsi_zhegv
                call chefsi_zhegv(nd_j, nmx_jj, hamm_batch(1,1,jd), nbandmx, ovlm_batch(1,1,jd), &
                     evl(1,isp_jj), evecs_dd, 15, 20, 1d-6, info_cf)
              endblock use_chefsi
            endif
            nev = nmx_jj
            if(call_m_bandcal_2nd .and. jd <= niqisp) then
              neviqis(jd) = nev; ndimhxiqis(jd) = nd_j
              block; complex(8),allocatable::etmp(:,:); allocate(etmp(nd_j,nev))
              etmp = evecs_dd; eveciqis(1:nd_j,1:nev,jd) = etmp; deallocate(etmp); endblock
            endif
            evl(nev+1:nbandmx,isp_jj) = 1d99
            nevls(iq_jj,isp_jj) = nev; ndimhx_(iq_jj,isp_jj) = nd_j
            if(allocated(t_evl(isp_jj,iq_jj)%v)) deallocate(t_evl(isp_jj,iq_jj)%v)
            allocate(t_evl(isp_jj,iq_jj)%v(nbandmx), source = evl(:,isp_jj))
            if(afsym) then
              if(allocated(t_evl(2,iq_jj)%v)) deallocate(t_evl(2,iq_jj)%v)
              allocate(t_evl(2,iq_jj)%v(nbandmx), source = evl(:,isp_jj))
              nevls(iq_jj,2) = nev; ndimhx_(iq_jj,2) = nd_j
            endif
            deallocate(evecs_dd)
          enddo
          if(master_mpi) write(stdo,'(a,f8.3,a)') ' Alt eigensolver total: ',MPI_WTIME()-t0cf,' s'
        endblock alt_diag
        deallocate(hamm_batch, ovlm_batch, ndimhx_batch, isp_batch, iq_batch, nmx_batch)
        if(allocated(ovlm_save)) deallocate(ovlm_save)
      elseif(cmdopt0('--batch')) then
      ! === Batched cuSOLVER+Ozaki with CUDA streams (epsovl=0, padded) ===
      ! === Batched cuSOLVER+Ozaki: fully pipelined CUDA streams ===
      ! === Batched cuSOLVER+Ozaki: all-GPU pipelined ===
      batched_diag: block
        use mpi, only: MPI_WTIME
        integer, parameter :: NSTR = 4
        integer(kind=8) :: strms(NSTR)
        type(cusolverDnHandle) :: cs_str(NSTR)
        complex(8), device, allocatable :: omat_all(:,:,:), Linv_all(:,:,:)
        real(8), device, allocatable :: eo_all(:,:)
        complex(8), device, allocatable :: S_d(:,:), work_dpf(:)
        complex(8), allocatable :: Lh(:,:)
        real(8), allocatable :: eo_h(:,:)
        real(8) :: t0b, t1b, t2b, t3b
        integer :: is, jd2, nd2, lw_pf, lw_ev, info_l2, jj2
        t0b = MPI_WTIME()
        nd2 = nd_max
        do is = 1, NSTR
          istat_g = cudaStreamCreate(strms(is))
          istat_g = cusolverDnCreate(cs_str(is))
          istat_g = cusolverDnSetStream(cs_str(is), strms(is))
        enddo
        allocate(omat_all(nd2,nd2,ndiag_total), Linv_all(nd2,nd2,ndiag_total))
        allocate(eo_all(nd2,ndiag_total))
        ! === Phase A1: GPU zpotrf on streams (pipelined) ===
        allocate(S_d(nd2,nd2))  ! temp for bufferSize
        block; real(8),device,allocatable::ed(:); allocate(ed(1))
        istat_g = cusolverDnZpotrf_bufferSize(cs_str(1), CUBLAS_FILL_MODE_LOWER, nd2, S_d, nd2, lw_pf)
        deallocate(ed); endblock
        deallocate(S_d)
        allocate(work_dpf(lw_pf * NSTR))
        do jd2 = 1, ndiag_total
          is = mod(jd2-1, NSTR) + 1
          allocate(S_d(nd2,nd2))
          S_d = ovlm_batch(1:nd2, 1:nd2, jd2)
          omat_all(:,:,jd2) = S_d  ! save S for zpotrf result
          deallocate(S_d)
          istat_g = cusolverDnZpotrf(cs_str(is), CUBLAS_FILL_MODE_LOWER, nd2, &
               omat_all(1,1,jd2), nd2, work_dpf((is-1)*lw_pf+1), lw_pf, devinfo)
        enddo
        do is = 1, NSTR; istat_g = cudaStreamSynchronize(strms(is)); enddo
        deallocate(work_dpf)
        t1b = MPI_WTIME()
        ! === Phase A2a: L^{-1} via GPU cuBLAS ztrsm on streams ===
        ! Solve L * X = I → X = L^{-1} for each k-point
        block
          use cublas_v2
          type(cublasHandle) :: cb_str(NSTR)
          complex(8), device, allocatable :: L_dk(:,:), I_dk(:,:)
          complex(8), allocatable :: Ih(:,:)
          integer :: is2
          ! Create per-stream cuBLAS handles
          do is2 = 1, NSTR
            istat_g = cublasCreate(cb_str(is2))
            istat_g = cublasSetStream(cb_str(is2), strms(is2))
          enddo
          ! Prepare identity on host (reuse)
          allocate(Ih(nd2,nd2)); Ih = (0d0,0d0)
          do jj2 = 1, nd2; Ih(jj2,jj2) = (1d0,0d0); enddo
          ! Submit all ztrsm on streams (pipelined, no sync between k-points)
          do jd2 = 1, ndiag_total
            is2 = mod(jd2-1, NSTR) + 1
            allocate(L_dk(nd2,nd2), I_dk(nd2,nd2))
            L_dk = omat_all(:,:,jd2)  ! L from zpotrf (lower triangle)
            I_dk = Ih                  ! identity
            istat_g = cublasZtrsm_v2(cb_str(is2), CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, &
                 CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, nd2, nd2, (1d0,0d0), L_dk, nd2, I_dk, nd2)
            Linv_all(:,:,jd2) = I_dk  ! L^{-1} (will sync via device-to-device copy)
            deallocate(L_dk, I_dk)
          enddo
          ! Sync all streams
          do is2 = 1, NSTR; istat_g = cudaStreamSynchronize(strms(is2)); enddo
          do is2 = 1, NSTR; istat_g = cublasDestroy(cb_str(is2)); enddo
          deallocate(Ih)
        endblock
        ! === Phase A2b: All H' projection on GPU (pipelined on streams) ===
        do jd2 = 1, ndiag_total
          is = mod(jd2-1, NSTR) + 1
          block
            complex(8), device, allocatable :: Li_d(:,:), H_d2(:,:), T_d2(:,:), Hp_d2(:,:)
            allocate(Li_d(nd2,nd2), H_d2(nd2,nd2), T_d2(nd2,nd2), Hp_d2(nd2,nd2))
            Li_d = Linv_all(:,:,jd2)
            H_d2 = hamm_batch(1:nd2, 1:nd2, jd2)
            istat_g = zmm(Li_d, H_d2, T_d2, m=nd2, n=nd2, k=nd2)
            istat_g = zmm(T_d2, Li_d, Hp_d2, m=nd2, n=nd2, k=nd2, opB=m_op_C)
            omat_all(:,:,jd2) = Hp_d2
            deallocate(Li_d, H_d2, T_d2, Hp_d2)
          endblock
        enddo
        t2b = MPI_WTIME()
        ! === Phase B: ALL Zheevd on streams (NO per-k sync) ===
        block
          complex(8), device, allocatable :: dummy_d(:,:)
          real(8), device, allocatable :: dummy_e(:)
          allocate(dummy_d(nd2,nd2), dummy_e(nd2))
          istat_g = cusolverDnZheevd_bufferSize(cs_str(1), CUSOLVER_EIG_MODE_VECTOR, &
               CUBLAS_FILL_MODE_UPPER, nd2, dummy_d, nd2, dummy_e, lw_ev)
          deallocate(dummy_d, dummy_e)
        endblock
        allocate(work_dpf(lw_ev * NSTR))
        do jd2 = 1, ndiag_total
          is = mod(jd2-1, NSTR) + 1
          istat_g = cusolverDnZheevd(cs_str(is), CUSOLVER_EIG_MODE_VECTOR, &
               CUBLAS_FILL_MODE_UPPER, nd2, omat_all(1,1,jd2), nd2, eo_all(1,jd2), &
               work_dpf((is-1)*lw_ev+1), lw_ev, devinfo)
        enddo
        do is = 1, NSTR; istat_g = cudaStreamSynchronize(strms(is)); enddo
        deallocate(work_dpf)
        t3b = MPI_WTIME()
        ! === Phase D: D2H eigenvalues + back-transform + store ===
        allocate(eo_h(nd2, ndiag_total)); eo_h = eo_all; deallocate(eo_all)
        do jd2 = 1, ndiag_total
          nd = ndimhx_batch(jd2); nev = nmx_batch(jd2)
          evl(1:nev, isp_batch(jd2)) = eo_h(1:nev, jd2)
          evl(nev+1:nbandmx, isp_batch(jd2)) = 1d99
          block
            complex(8), device, allocatable :: Y_d2(:,:), evec_d2(:,:), Li_d(:,:)
            allocate(Li_d(nd2,nd2), Y_d2(nd2,nev), evec_d2(nd2,nev))
            Li_d = Linv_all(:,:,jd2)
            Y_d2 = omat_all(1:nd2, 1:nev, jd2)
            istat_g = zmm(Li_d, Y_d2, evec_d2, m=nd2, n=nev, k=nd2, opA=m_op_C)
            if(call_m_bandcal_2nd .and. jd2 <= niqisp) then
              neviqis(jd2) = nev; ndimhxiqis(jd2) = nd
              block; complex(8),allocatable::et(:,:); allocate(et(nd,nev))
              et = evec_d2(1:nd,1:nev); eveciqis(1:nd,1:nev,jd2) = et; deallocate(et); endblock
            endif
            deallocate(Li_d, Y_d2, evec_d2)
          endblock
          nevls(iq_batch(jd2),isp_batch(jd2)) = nev
          ndimhx_(iq_batch(jd2),isp_batch(jd2)) = nd
          if(allocated(t_evl(isp_batch(jd2),iq_batch(jd2))%v)) &
               deallocate(t_evl(isp_batch(jd2),iq_batch(jd2))%v)
          allocate(t_evl(isp_batch(jd2),iq_batch(jd2))%v(nbandmx), source=evl(:,isp_batch(jd2)))
          if(afsym) then
            if(allocated(t_evl(2,iq_batch(jd2))%v)) deallocate(t_evl(2,iq_batch(jd2))%v)
            allocate(t_evl(2,iq_batch(jd2))%v(nbandmx), source=evl(:,isp_batch(jd2)))
            nevls(iq_batch(jd2),2) = nev; ndimhx_(iq_batch(jd2),2) = nd
          endif
        enddo
        deallocate(eo_h, omat_all, Linv_all)
        do is = 1, NSTR
          istat_g = cusolverDnDestroy(cs_str(is))
          istat_g = cudaStreamDestroy(strms(is))
        enddo
        if(master_mpi) write(stdo,'(a,f8.3,a,3(a,f6.3))') ' Batched diag: ',MPI_WTIME()-t0b,' s', &
             ' zpotrf=',t1b-t0b,' proj+inv=',t2b-t1b,' zheevd=',t3b-t2b
        deallocate(hamm_batch, ovlm_batch, ndimhx_batch, isp_batch, iq_batch, nmx_batch)
        if(allocated(ovlm_save)) deallocate(ovlm_save)
      endblock batched_diag
      else
      ! === cuSOLVER path (original sequential) ===
      timing: block
      use mpi, only: MPI_WTIME
      real(8) :: t_ovl, t_zz, t_gemm, t_hdiag, t_bt, t_store, t0
      logical :: use_chefsi_dummy
      t_ovl=0; t_zz=0; t_gemm=0; t_hdiag=0; t_bt=0; t_store=0
      do jdat = 1, ndiag_total
        nd = ndimhx_batch(jdat)
        nmx_j = nmx_batch(jdat)
        isp_j = isp_batch(jdat)
        iq_j = iq_batch(jdat)
        if(skip_ovldiag) then
        ! === epsovl=0: Cholesky path (no eigensolve for overlap) ===
        ! Step 1: Cholesky S = L*L^H on CPU, then L^{-1}
        t0 = MPI_WTIME()
        nm_j = nd_max
        cholesky_path: block
          complex(8), allocatable :: linv_h(:,:)
          complex(8), device, allocatable :: linv_d(:,:)
          integer :: info_l, ii, jj
          allocate(linv_h(nd_max, nd_max))
          linv_h = ovlm_batch(1:nd_max, 1:nd_max, jdat)
          call zpotrf('L', nd_max, linv_h, nd_max, info_l)  ! S = L*L^H
          call ztrtri('L', 'N', nd_max, linv_h, nd_max, info_l)  ! L^{-1}
          ! Zero upper triangle (ztrtri only modifies lower)
          do jj = 2, nd_max; linv_h(1:jj-1, jj) = (0d0, 0d0); enddo
          t_ovl = t_ovl + MPI_WTIME() - t0
          ! Step 3: H' = L^{-1} * H * L^{-H} on GPU with Ozaki GEMM
          t0 = MPI_WTIME()
          allocate(linv_d(nd_max,nd_max), h_d(nd_max,nd_max), hhm_d(nd_max,nd_max), hh_d(nd_max,nd_max))
          linv_d = linv_h
          h_d = hamm_batch(1:nd_max, 1:nd_max, jdat)
          istat_g = zmm(linv_d, h_d, hhm_d, m=nd_max, n=nd_max, k=nd_max)  ! L^{-1} * H
          istat_g = zmm(hhm_d, linv_d, hh_d, m=nd_max, n=nd_max, k=nd_max, opB=m_op_C)  ! (L^{-1}*H) * L^{-H}
          deallocate(hhm_d, h_d)
          ! Store linv_d for back-transform (Step 5 uses L^{-H})
          allocate(zz_d(nd_max, nd_max))
          zz_d = linv_d  ! zz_d = L^{-1}, back-transform uses L^{-H} via opA=C
          deallocate(linv_d, linv_h)
        endblock cholesky_path
        allocate(zz_h(1,1))  ! dummy for dealloc later
        t_gemm = t_gemm + MPI_WTIME() - t0
        t_zz = 0
        else
        ! === epsovl>0: Full eigensolve path ===
        ! Step 1: Overlap eigensolve on GPU
        t0 = MPI_WTIME()
        allocate(omat_d(nd_max,nd_max), eo_d(nd_max))
        omat_d = ovlm_batch(1:nd_max,1:nd_max,jdat)
        istat_g = cusolverDnZheevd_bufferSize(cs_h, CUSOLVER_EIG_MODE_VECTOR, &
             CUBLAS_FILL_MODE_UPPER, nd_max, omat_d, nd_max, eo_d, lwork_j)
        allocate(work_d(lwork_j))
        istat_g = cusolverDnZheevd(cs_h, CUSOLVER_EIG_MODE_VECTOR, &
             CUBLAS_FILL_MODE_UPPER, nd_max, omat_d, nd_max, eo_d, work_d, lwork_j, devinfo)
        deallocate(work_d)
        eo(1:nd_max) = eo_d
        deallocate(eo_d)
        t_ovl = t_ovl + MPI_WTIME() - t0
        ! Step 2: Build zz = S^{-1/2} on host
        t0 = MPI_WTIME()
        ni_j = 1
        do ix_j = 1, nd_max
          if(eo(ix_j) > epsovl) then; ni_j = ix_j; exit; endif
        enddo
        nm_j = nd_max - ni_j + 1
        allocate(zz_h(nd_max, nm_j))
        zz_h = omat_d(:, ni_j:nd_max)
        deallocate(omat_d)
        do ix_j = ni_j, nd_max
          zz_h(:, ix_j-ni_j+1) = zz_h(:, ix_j-ni_j+1) / sqrt(max(eo(ix_j), 1d-30))
        enddo
        t_zz = t_zz + MPI_WTIME() - t0
        ! Step 3: Project H on GPU with Ozaki GEMM
        t0 = MPI_WTIME()
        allocate(zz_d(nd_max, nm_j), h_d(nd_max, nd_max), hhm_d(nm_j, nd_max), hh_d(nm_j, nm_j))
        zz_d = zz_h
        h_d = hamm_batch(1:nd_max, 1:nd_max, jdat)
        istat_g = zmm(zz_d, h_d, hhm_d, m=nm_j, n=nd_max, k=nd_max, opA=m_op_C)
        istat_g = zmm(hhm_d, zz_d, hh_d, m=nm_j, n=nm_j, k=nd_max)
        deallocate(hhm_d, h_d)
        t_gemm = t_gemm + MPI_WTIME() - t0
        endif
        ! === Step 4: H eigensolve on GPU (cuSOLVER Zheevd) ===
        t0 = MPI_WTIME()
        allocate(eo_d(nm_j))
        istat_g = cusolverDnZheevd_bufferSize(cs_h, CUSOLVER_EIG_MODE_VECTOR, &
             CUBLAS_FILL_MODE_UPPER, nm_j, hh_d, nm_j, eo_d, lwork_j)
        allocate(work_d(lwork_j))
        istat_g = cusolverDnZheevd(cs_h, CUSOLVER_EIG_MODE_VECTOR, &
             CUBLAS_FILL_MODE_UPPER, nm_j, hh_d, nm_j, eo_d, work_d, lwork_j, devinfo)
        deallocate(work_d)
        ! Take lowest nev_j eigenvalues
        nev_j = min(nmx_j, nm_j)
        evl(1:nev_j, isp_j) = eo_d(1:nev_j)  ! D2H: eigenvalues (lowest nev_j)
        deallocate(eo_d)
        t_hdiag = t_hdiag + MPI_WTIME() - t0
        ! hh_d now contains all eigenvectors; take first nev_j columns
        ! === Step 5: Back-transform on GPU with Ozaki GEMM ===
        t0 = MPI_WTIME()
        allocate(evec_d(nm_j, nev_j), z_d(nd_max, nev_j))
        evec_d = hh_d(1:nm_j, 1:nev_j)  ! first nev_j eigenvectors
        if(skip_ovldiag) then
          ! Cholesky path: x = L^{-H} * y = (L^{-1})^H * y; zz_d = L^{-1}
          istat_g = zmm(zz_d, evec_d, z_d, m=nd_max, n=nev_j, k=nm_j, opA=m_op_C)
        else
          ! Eigensolve path: x = zz * y
          istat_g = zmm(zz_d, evec_d, z_d, m=nd_max, n=nev_j, k=nm_j)
        endif
        deallocate(evec_d, hh_d, zz_d)
        allocate(evec_tmp(nd_max, nev_j))
        evec_tmp = z_d  ! D2H
        deallocate(z_d)
        t_bt = t_bt + MPI_WTIME() - t0
        nev = nev_j
        ! Store results (use original nd)
        t0 = MPI_WTIME()
        if(call_m_bandcal_2nd .and. jdat <= niqisp) then
          neviqis(jdat) = nev
          ndimhxiqis(jdat) = nd
          if(nmx_j/=0) eveciqis(1:nd,1:nev,jdat) = evec_tmp(1:nd,1:nev)
        endif
        evl(nev+1:nbandmx,isp_j) = 1d99
        nevls(iq_j,isp_j) = nev
        ndimhx_(iq_j,isp_j) = nd
        if(allocated(t_evl(isp_j,iq_j)%v)) deallocate(t_evl(isp_j,iq_j)%v)
        allocate(t_evl(isp_j,iq_j)%v(nbandmx), source = evl(:,isp_j))
        if(afsym) then
          if(allocated(t_evl(2,iq_j)%v)) deallocate(t_evl(2,iq_j)%v)
          allocate(t_evl(2,iq_j)%v(nbandmx), source = evl(:,isp_j))
          nevls(iq_j,2) = nev
          ndimhx_(iq_j,2) = nd
        endif
        t_store = t_store + MPI_WTIME() - t0
        deallocate(evec_tmp, zz_h)
      enddo
      if(master_mpi) write(stdo,'(a,6(a,f8.3))') ' BatchDiag timing(s):', &
           ' ovl=',t_ovl,' zz=',t_zz,' gemm=',t_gemm,' Hdiag=',t_hdiag,' bt=',t_bt,' store=',t_store
      endblock timing
      deallocate(devinfo)
      istat_g = cusolverDnDestroy(cs_h)
      deallocate(hamm_batch, ovlm_batch, ndimhx_batch, isp_batch, iq_batch, nmx_batch)
      if(allocated(ovlm_save)) deallocate(ovlm_save)
      endif ! --chefsi / cuSOLVER
    endblock PadAndBatchDiag
    endif ! use_gpu
    ! === Phase 4: Scatter results back to CPU ranks ===
    ScatterResults: block
      use m_gpu, only: use_gpu, ngpu_ranks
      use mpi, only: MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE, MPI_MAX, MPI_SUM, MPI_IN_PLACE
      integer :: ierr_s
      ! Combine nevls and ndimhx_ from all ranks
      nevls_reduce: block
        integer, allocatable :: nevls_recv(:,:), ndimhx_recv(:,:)
        allocate(nevls_recv(size(nevls,1),size(nevls,2)), ndimhx_recv(size(ndimhx_,1),size(ndimhx_,2)))
        call mpi_allreduce(nevls, nevls_recv, size(nevls), MPI_INTEGER, MPI_MAX, comm, ierr_s)
        nevls = nevls_recv
        call mpi_allreduce(ndimhx_, ndimhx_recv, size(ndimhx_), MPI_INTEGER, MPI_MAX, comm, ierr_s)
        ndimhx_ = ndimhx_recv
        deallocate(nevls_recv, ndimhx_recv)
      endblock nevls_reduce
      ! Broadcast t_evl: pack into contiguous array, allreduce, unpack
      evl_bcast: block
        use m_qplist, only: nkp
        use m_gpu, only: use_gpu
        use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
        real(8), allocatable :: evl_all(:,:,:)
        integer :: jq, jsp
        allocate(evl_all(nbandmx, nspx, nkp), source=0d0)
        if(use_gpu) then
          do jq = 1, nkp
            do jsp = 1, nspx
              if(allocated(t_evl(jsp,jq)%v)) evl_all(:,jsp,jq) = t_evl(jsp,jq)%v
            enddo
          enddo
        endif
        evl_sum: block
          real(8), allocatable :: evl_recv(:,:,:)
          allocate(evl_recv(nbandmx, nspx, nkp))
          call mpi_allreduce(evl_all, evl_recv, size(evl_all), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr_s)
          evl_all = evl_recv
          deallocate(evl_recv)
        endblock evl_sum
        do jq = 1, nkp
          do jsp = 1, nspx
            if(.not.allocated(t_evl(jsp,jq)%v)) allocate(t_evl(jsp,jq)%v(nbandmx))
            t_evl(jsp,jq)%v = evl_all(:,jsp,jq)
          enddo
        enddo
        deallocate(evl_all)
      endblock evl_bcast
      ! Scatter eveciqis: GPU ranks send diag results back to CPU ranks
      if(call_m_bandcal_2nd) then
      scatter_evec: block
        use mpi, only: MPI_STATUS_SIZE
        integer :: src_gpu, dst, nk_send, jdat_off, ierr_e, jd
        integer :: status_e(MPI_STATUS_SIZE)
        integer :: nk_all(0:numprocs-1)
        integer, allocatable :: nev_send(:), ndimhx_send(:)
        complex(8), allocatable :: evec_send(:,:,:)
        call mpi_allgather(niqisp, 1, MPI_INTEGER, nk_all, 1, MPI_INTEGER, comm, ierr_e)
        if(use_gpu) then
          ! GPU rank: send eveciqis for remote k-points back to CPU ranks
          ! BatchDiag stored: jdat 1..niqisp = local, jdat niqisp+1..ndiag_total = remote
          ! Remote k-points need eveciqis sent back. We stored them in a temp buffer during BatchDiag.
          ! But current BatchDiag doesn't store eveciqis for remote k-points (jdat > niqisp).
          ! We need to re-run diag or store them. For simplicity, send neviqis/ndimhxiqis/eveciqis.
          ! Actually, BatchDiag skipped storing for jdat > niqisp. Fix: store in temp arrays.
          ! For now, send nevls/ndimhx_ (already in allreduce) and let CPU ranks re-diag for eveciqis.
          continue  ! GPU rank done - results already in nevls/ndimhx_/t_evl via allreduce
        else
          ! CPU rank: need eveciqis for bandcal_2nd. Diag locally with CPU LAPACK.
          cpu_diag: block
            complex(8), allocatable :: hamm_loc(:,:), ovlm_loc(:,:), evec_loc(:,:)
            integer :: nd_loc, nmx_loc, nev_loc, iq_loc, isp_loc
            do jd = 1, niqisp
              iq_loc = iqproc(jd)
              isp_loc = isproc(jd)
              nd_loc = ndimhx_batch(jd)
              nmx_loc = nmx_batch(jd)
              allocate(hamm_loc(nd_loc,nd_loc), ovlm_loc(nd_loc,nd_loc), evec_loc(nd_loc,nmx_loc))
              hamm_loc = hamm_batch(1:nd_loc,1:nd_loc,jd)
              ovlm_loc = ovlm_batch(1:nd_loc,1:nd_loc,jd)
              call zhev_tk4(nd_loc, hamm_loc, ovlm_loc, nmx_loc, nev_loc, evl(1,isp_loc), evec_loc, epsovl)
              neviqis(jd) = nev_loc
              ndimhxiqis(jd) = nd_loc
              nevls(iq_loc,isp_loc) = nev_loc
              ndimhx_(iq_loc,isp_loc) = nd_loc
              if(nmx_loc/=0) eveciqis(1:nd_loc,1:nev_loc,jd) = evec_loc(1:nd_loc,1:nev_loc)
              deallocate(hamm_loc, ovlm_loc, evec_loc)
            enddo
          endblock cpu_diag
          if(allocated(hamm_batch)) deallocate(hamm_batch, ovlm_batch, ndimhx_batch, isp_batch, iq_batch, nmx_batch)
        endif
      endblock scatter_evec
      endif
    endblock ScatterResults
    if(use_gpu) then
    gpumem_after: block
      use cudafor
      integer(8) :: free_mem, total_mem
      integer :: ierr_mem
      ierr_mem = cudaMemGetInfo(free_mem, total_mem)
      write(6,'(a,2f10.1,a)') ' GPU mem after k-loop: free/total(MB)=', &
           free_mem/1d6, total_mem/1d6, ' MB'
    endblock gpumem_after
    call zhev_gpu_cleanup()
    gpumem_freed: block
      use cudafor
      integer(8) :: free_mem, total_mem
      integer :: ierr_mem
      ierr_mem = cudaMemGetInfo(free_mem, total_mem)
      write(6,'(a,2f10.1,a)') ' GPU mem after cleanup: free/total(MB)=', &
           free_mem/1d6, total_mem/1d6, ' MB'
    endblock gpumem_freed
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
#ifdef __GPU
  subroutine blocked_ztrtri_ozaki(n, L, nb)
    !! Blocked lower-triangular inversion. Off-diagonal updates via Ozaki GEMM.
    !! LAPACK ztrtri algorithm: ztrmm→Ozaki GEMM, ztrsm→small_inv+Ozaki GEMM.
    use m_blas, only: zmm_oz => zmm_d, m_op_C_oz => m_op_C
    integer, intent(in) :: n, nb
    complex(8), intent(inout) :: L(n, n)
    complex(8), device, allocatable :: Lblk_d(:,:), B_d(:,:), T_d(:,:), Dinv_d(:,:)
    complex(8), allocatable :: Dinv_h(:,:)
    integer :: j, jb, nn, istat_oz, info_oz, k2
    do j = 1, n, nb
      jb = min(nb, n - j + 1)
      nn = n - j - jb + 1
      if(nn > 0) then
        ! Step 1: T = -L(j+jb:n, j+jb:n) * L(j+jb:n, j:j+jb-1)  [Ozaki GEMM]
        ! L(j+jb:n, j+jb:n) is already inverted from previous iterations
        allocate(Lblk_d(nn,nn), B_d(nn,jb), T_d(nn,jb))
        Lblk_d = L(j+jb:n, j+jb:n)
        B_d = L(j+jb:n, j:j+jb-1)
        istat_oz = zmm_oz(Lblk_d, B_d, T_d, m=nn, n=jb, k=nn)
        ! Negate
        block; complex(8),allocatable::th(:,:); allocate(th(nn,jb)); th=T_d; th=-th; T_d=th; deallocate(th); endblock
        ! Step 2: Invert diagonal block jb×jb on CPU
        allocate(Dinv_h(jb,jb))
        Dinv_h = L(j:j+jb-1, j:j+jb-1)
        call ztrtri('L', 'N', jb, Dinv_h, jb, info_oz)
        L(j:j+jb-1, j:j+jb-1) = Dinv_h
        ! Step 3: T = T * Dinv  [Ozaki GEMM]
        allocate(Dinv_d(jb,jb))
        Dinv_d = Dinv_h
        istat_oz = zmm_oz(T_d, Dinv_d, B_d, m=nn, n=jb, k=jb)
        L(j+jb:n, j:j+jb-1) = B_d
        deallocate(Lblk_d, B_d, T_d, Dinv_d, Dinv_h)
      else
        ! Last block: invert on CPU
        call ztrtri('L', 'N', jb, L(j,j), n, info_oz)
      endif
    enddo
  end subroutine
#endif
end module m_bandcal
