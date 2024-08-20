!> Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz(inverse of overlap matrix) !  "call get_zmel_init" return zmel 
!  All dependencies (use foobar below ) are inputs (must be protected).
module m_zmel
  use m_genallcf_v3,only: nclass,natom,nspin,nl,nn,nnv,nnc, nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw
  use m_genallcf_v3,only: alat,delta,deltaw,esmr,iclass,nlnmv,nlnmc,icore,ncore,plat,pos,z,ecore,mnl=>nlnm,nl,nn,nlnmx,il,in,im
  use m_hamindex,only: ngrp, symgg=>symops,invg=>invgx
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
  use m_read_bzdata,only: nqbz,nqibz,  qlat,ginv,qbz,qibz,wbz, done_read_bzdata
  use m_readhbe,only: nband
  use m_readQG,only: ngpmx,ngcmx,Readqg
  use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
  use m_ftox
  use m_lgunit,only:stdo
  use m_MPItk,only:master_mpi
  use m_mem,only: memused,writemem
  public:: get_zmel_init_gemm, Mptauof_zmel,Setppovlz,Setppovlz_chipm ! Call mptauof_zmel and setppovlz in advance to get_zmel_init
  complex(8),allocatable,protected,public :: zmel(:,:,:) ! OUTPUT: zmel(nbb,nmtot, nqtot) ,nbb:mixproductbasis, nmtot:middlestate, nqtot:endstate
  complex(8),allocatable,protected,public :: ppovlz(:,:)
  real(8),allocatable,protected,public :: tiat(:,:,:),shtvg(:,:), ppbir(:,:,:,:,:,:)
  integer,protected,public:: nbb 
  integer,allocatable,protected,public :: miat(:,:)
  private
  real(8),parameter:: kk=1000
contains
  subroutine setppovlz(q,matz,npr) ! Set ppovlz for given q
    intent(in)::       q,matz,npr
    logical:: matz
    ! 2024-5-24; add npr
    ! Set ppovlz(ngb,npr), where
    !    ngb: the size of MPB
    !   npr: the degree of freedom to calculate polarization funciton.
    !          npr=1 for eps mode nolfc., npr=ngb for no nolfc (epsPP0 mode)
    !
    !    ppolvz(igb,ivcou)= (1    0 ) \times  zcousq(igb, ivcou)
    !                       (0 ppovl)
    !    If matz=F, no multiplication by ivcou.  Thus we have ppolz(igb,igb)
    !     
    real(8) :: q(3)
    complex(8),allocatable :: ppovl_(:,:),ppovl(:,:)!,ppovlzinv(:,:) !    logical:: eibz4x0
    integer:: i,ngc_r,ippovl0,npr
    real(8):: qx(3),tolq=1d-8
    if(allocated(ppovlz)) deallocate(ppovlz)
    if(allocated(ppovl))  deallocate(ppovl)
    allocate( ppovl(ngc,ngc),ppovlz(ngb,npr))
    open(newunit=ippovl0,file='PPOVL0',form='unformatted') !inefficient search for PPOVLO for given q
    do 
       read(ippovl0) qx,ngc_r
       if(sum(abs(qx-q))<tolq) then
          if(ngc_r/=ngc) call rx( 'readin ppovl: ngc_r/=ngc')
          read(ippovl0) ppovl
          exit
       endif
    enddo
    close(ippovl0)
    ppovlz(1:nbloch,1:npr) = zcousq(1:nbloch,1:npr)
    ppovlz(nbloch+1:nbloch+ngc,1:npr)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,1:npr))
    deallocate(ppovl)
    nbb=npr    ! ngb obatabugfix 2025-5-23. We had set nbb=ngb every time. Thus we had memory(and computational) loss for nolfc case.
  end subroutine setppovlz
  subroutine setppovlz_chipm(zzr,nmbas1) !Set ppovlz for chipm case
    intent(in)::             zzr,nmbas1
    integer::nmbas1
    complex(8):: zzr(ngb,nmbas1)
    if(allocated(ppovlz)) deallocate(ppovlz)
    allocate(ppovlz(ngb,nmbas1))
    ppovlz= zzr
    nbb   = nmbas1
  end subroutine setppovlz_chipm
  subroutine mptauof_zmel(symops,ng)! Set miat,tiat,invgx,shtvg, and then get ppbir
    use m_mksym_util,only:mptauof
    use m_hamindex0,only: Readhamindex0,iclasst
    intent(in)::          symops,ng
    integer:: ng
    real(8):: symops(9,ng)
    integer,allocatable ::  invgx(:)
    call readhamindex0()
    allocate(invgx(ng),miat(natom,ng),tiat(3,natom,ng),shtvg(3,ng))
    call mptauof(symops,ng,plat,natom,pos,iclasst,miat,tiat,invgx,shtvg ) !Set miat,tiat,shtvg for space group.
    deallocate(invgx)
    call rdpp(ng,symops)  !return ppbrd:radial integrals and cgr:rotated cg coeffecients. 
    ppbafp_v2_zmel: block 
      integer :: is,irot,lmxax, ic, i,lb,nb,mb,lmb,i1,ibas,i2, np,lp,mp,lmp,n,l,m,lm
      lmxax=nl-1
      allocate(ppbir(nlnmx,nlnmx,mdimx,nclass,ng,nspin)) ! ppbir is rotated <Phi(SLn,r) Phi(SL'n',r) B(S,i,rot^{-1}(r))> by rotated cg coefficients cgr
      do irot = 1,ng
         do is = 1,nspin 
            do concurrent (ic=1:nclass)
               ibas = ic
               i = 0 !i = product basis index.
               do lb  = 0, lx (ibas)
                  do nb  = 1, nx (lb,ibas)
                     do mb  = -lb, lb
                        i    = i+1           !The number of product basis is  =(i at the end of loop).
                        lmb  = lb*lb + lb + mb + 1
                        do concurrent (i2 = 1:mnl(ic),i1 = 1:mnl(ic)) !phi1 phi2 index
                           np= in(i2,ic);  lp= il(i2,ic);  mp= im(i2,ic);  lmp=lp*lp + lp + mp + 1
                           n = in(i1,ic);  l = il(i1,ic);  m = im(i1,ic);  lm = l*l + l + m + 1
                           ppbir(i1,i2,i,ic,irot,is)=cgr(lm,lmp,lmb,irot)*ppbrd(l,n,lp,np,lb,nb,is+nspin*(ic-1))
                           !note ppbir is not necesssairy if i1,i2 are different spins
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
    endblock ppbafp_v2_zmel
  end subroutine mptauof_zmel

  subroutine get_zmel_init_gemm(q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc,  & ! get_zmel_init for blas/cuBLAS by M. Obata 2024-05-18
       iprx,zmelconjg,comm,maxmem)
    use m_readeigen,only: readcphif 
    use m_readeigen,only: readgeigf
    use m_itq,only: itq, ntq
    use m_kind, only: kp => kindgw
    use m_blas, only: m_op_c, m_op_n, m_op_t, int_split
    use m_blas, only: gemm => zmm, gemm_batch => zmm_batch
#ifdef __GPU
    use openacc, only: acc_is_present
#endif
    implicit none
    include "mpif.h"
    intent(in)::           q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg
    integer, optional, intent(in) :: comm
    real(8), parameter::tolq=1d-8
    complex(8), parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
    integer:: ns1, ns2, nqmax, irot, ispq, ispm, nqini, nctot, ncc, ncnv, ncorec, nccc, mdim, it, ia
    integer:: ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
    integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime
    real(8):: quu(3),q(3), kvec(3),rkvec(3),qkt(3),qt(3), qdiff(3)
    real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass) ! ppb= <Phi(SLn,r-R)_q,isp1 |Phi(SL'n',r-R)_qk,isp2 B_k(S,i,rot^{-1}(r-R))>
    logical:: iprx,zmelconjg,debug,cmdopt0
    integer,allocatable:: ngveccR(:,:),igcgp2i_(:,:)
    complex(8)::cphiq(nlmto,nband), cphim(nlmto,nband)
    complex(8),allocatable:: geigq(:,:),dgeigqk(:,:)
    integer:: invr,nt0,ntp0,nmtot,nqtot
    integer:: iasx(natom),icsx(natom),iatomp(natom),imdim(natom)
    real(8)::tr(3,natom),qk(3),symope(3,3),shtv(3)
    integer :: ierr, nqini_rank, nqmax_rank, ntp0_rank ,nm1,nm2,nm1c,nm2c,nm1v,nm2v,nm1cc,nm2cc
    integer :: mpi_rank, mpi_size, ini_index, end_index, num_index, mpi_info, irank
    real(8),optional::maxmem
    character(8),external:: charext
    complex(kind=kp), parameter :: CONE = (1_kp, 0_kp), CZERO = (0_kp, 0_kp)
    complex(8),allocatable:: zmelp0(:,:,:),ggitp_(:,:) 
    complex(8),allocatable:: ggitp(:,:), gggmat(:,:), ggitp_work(:,:)
    complex(8), allocatable :: zmel_buf(:,:,:)
    complex(8), allocatable :: zmelt_d(:,:,:), zmelt(:,:,:)
    complex(8), allocatable :: ppbvphiq_d(:,:,:), cphim_d(:,:), cphiq_d(:,:), ppbc_d(:,:,:), ppbv_d(:,:,:)
    debug=cmdopt0('--debugzmel')
#ifdef __GPU
    attributes(device) :: ppb
#endif
    if(allocated(zmel)) then
#ifdef __GPU
       if(acc_is_present(zmel, size(zmel))) then
         !$acc exit data delete(zmel)
       endif
#endif
      deallocate(zmel)
    endif
    nm1=ns1
    nm2=ns2
    nmtot = nm2-nm1+1 
    SetRangeOfCoreAndValence :if(nm1>nctot) then ! Core nm1c:nm2c, Valence nm1v:nm2v
       nm1c=0
       nm2c=-1
       nm1v=nm1
       nm2v=nm2
    elseif(nm2>nctot) then
       nm1c=nm1
       nm2c=nctot
       nm1v=nctot+1
       nm2v=nm2
    else !skip valence
       nm1c=nm1
       nm2c=nm2
       nm1v=0
       nm2v=-1
    endif SetRangeOfCoreAndValence
    !! For given range of band index nm1;nm2, We now get nm1v:nmv2 and nm1c:nm2c, which are range of valence and core index.
    !! Note that band index is nctot+nvalence order.
    ntp0   = nqmax-nqini+1
    nqtot  = ncc + ntp0    ! = phi_end
    ini_index = 1
    end_index = ntp0
    nqini_rank = nqini
    nqmax_rank = nqmax
    SetByCPU :block
      qk =  q - rkvec ! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      associate(cphitemp=> readcphif(q,ispq))
        cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq)) 
      endassociate
      cphim = readcphif(qk, ispm) 
      symope= symgg(:,:,irot)
      allocate(geigq(ngpmx,nband),dgeigqk(ngpmx,nband))
      if(ngc/=0) then
        call readqg('QGpsi',q,     qt, ngp1, ngvecpB1) !q is mapped to qt in BZ
        call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2)
        qdiff = matmul(symope,kvec)-qt+qkt ! rkvec + qkt - qt is not zero. <M(rkvec) Phi(qk) |Phi(q)>
        nadd = nint(matmul(transpose(plat), qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
        block
          use m_read_ppovl,only: getppx2, ngvecc,ngcread
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
          call getppx2(qlat,kvec) ! read and allocate ppovlinv
          if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')
          allocate(ngveccR(1:3,1:ngc))
          call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
        endblock
        geigq   = readgeigf(q, ispq) !read IPW part at q   !G1 for ngp1
        dgeigqk = readgeigf(qk,ispm) !read IPW part at qk  !G2 for ngp2
        dgeigqk = dconjg(dgeigqk)
      endif
      !$acc data copyin(ppbir(1:nlnmx,1:nlnmx,1:mdimx,1:nclass,irot,ispq))
      !$acc kernels
      ppb = ppbir(:,:,:,:,irot,ispq)           !MPB has no spin dependence
      !$acc end kernels
      !$acc end data 
      invr  = invg(irot)       !invrot (irot,invg,ngrp) ! Rotate atomic positions invrot*R = R' + T
      tr    = tiat(:,:,invr)
      iatomp= miat(:,invr)
      shtv  = matmul(symope,shtvg(:,invr))
      imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
      iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
      icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
      if(present(comm)) then
        call mpi_comm_rank(comm, mpi_rank, mpi_info)
        call mpi_comm_size(comm, mpi_size, mpi_info)
        call int_split(ntp0, mpi_size, mpi_rank, ini_index, end_index, num_index)
        ntp0 = num_index
        nqini_rank = nqini + ini_index - 1
        nqmax_rank = nqini + end_index - 1
      endif
    end block SetByCPU
    if(debug) write(stdo,ftox)'zmel_init gpu',nbloch,ngc,nm1,nm2,nqtot
    ZmelBlock:block
#ifdef __GPU
      attributes(device) :: zmelt, zmelt_d
#endif
      call writemem('mmmmm_zmel000: zmelsize='//ftof((nbloch+ngc)*(nm2-nm1+1)*nqtot*16/kk**3)//' GB')
      allocate(zmelt(1:nbloch+ngc,nm1:nm2,1:nqtot))
!$acc kernels
      zmelt(1:nbloch+ngc,nm1:nm2,1:nqtot) = czero
!$acc end kernels
      ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
        complex(8):: phasea(natom) 
#ifdef __GPU
        attributes(device) :: ppbvphiq_d, cphim_d, cphiq_d, ppbc_d, ppbv_d
#endif
        phasea = [(exp(-img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)]
        iatomloop: do ia = 1, natom
          ic    = iclass(ia)
          nc    = nlnmc(ic) !nlnmc      = number of l,n,m for core states
          nv    = nlnmv(ic) !nlnmv      = number of l,n,m for valence
          nc1   = nc + 1
          ncnv  = nc+nv  !ncore + nvalence
          iap   = iatomp(ia)  
          icp   = iclass(iap)
          ias   = iasx(ia)  !start of nlnmv for ia
          iae   = ias+nv-1  !end  
          ics   = icsx(ia)   
          ims   = imdim(iap) !start of PB for ia
          ime   = ims-1+nblocha(icp)
          mdim  = nblocha(icp)
          ncorec= merge(ncore(ic),0,nctot>0)
          if(ncc/=0) call rx('ncc/=0 get_zmel_init: Not implemented nccc=0 yet')
!          nccc  = merge(ncore(ic),0,ncc>0) !ncc is expected to be zero
          if (nm2v>=nm1v) allocate(cphim_d(nv,nm1v:nm2v), source= cphim(ias:iae,nm1v:nm2v)) !copy from CPU to GPU
          allocate(                cphiq_d(nv,ntp0),      source= cphiq(ias:iae,nqini_rank:nqmax_rank))
          ValenceValence: if (nm2v>=nm1v) then 
            allocate(ppbvphiq_d(nv,ntp0,mdim))
            allocate(   zmelt_d(nm1v:nm2v,ntp0,mdim))
            allocate(    ppbv_d(nv,nv,mdim))
            !$acc kernels
            ppbv_d(1:nv,1:nv,1:mdim) = dcmplx(ppb(nc1:ncnv,nc1:ncnv,1:mdim,icp)) 
            !$acc end kernels
            ierr=gemm_batch(ppbv_d, cphiq_d,ppbvphiq_d,M=nv,N=ntp0,K=nv,         NBATCH=mdim,                  opA=m_op_T,&
                 sameB=.true.)
            ierr=gemm_batch(cphim_d,ppbvphiq_d,zmelt_d,M=nm2v-nm1v+1,N=ntp0,K=nv,NBATCH=mdim,alpha=phasea(ia), opA=m_op_C,&
                 sameA=.true.)
            !$acc kernels
            do i = 1, mdim
              zmelt(i-1+ims,nm1v:nm2v,ncc+1:ncc+ntp0) = zmelt_d(nm1v:nm2v,1:ntp0,i)
            enddo
            !$acc end kernels
            deallocate(zmelt_d, ppbvphiq_d, ppbv_d)
          endif ValenceValence
          nm1cc = max(nm1c,ics+1)        !core index range between [nm1c,nm2c]. This corresponds to  nm1cc:nm2cc for atom ia.
          nm2cc=  min(nm2c,ics+ncorec)  !       write(6,*)'ia nm1cc nm2cc=',ia, nm1cc,nm2cc,ntp0,ncc,mdim
          CoreValence: if(nm2cc>=nm1cc) then 
            allocate(zmelt_d(mdim,ntp0,nm1cc:nm2cc), ppbc_d(mdim,nv,nm1cc:nm2cc))
            !$acc kernels
            do i = 1, mdim
              ppbc_d(i,1:nv,nm1cc:nm2cc) = dcmplx(ppb(nc1:ncnv,nm1cc-ics:nm2cc-ics,i,icp))
            enddo
            !$acc end kernels
            ierr = gemm_batch(ppbc_d,cphiq_d,zmelt_d, M=mdim,N=ntp0,K=nv,NBATCH=nm2cc-nm1cc+1, alpha=phasea(ia), sameB = .true.)
            !$acc kernels
            do it = nm1cc,nm2cc
               zmelt(ims:ime,it,ncc+1:ncc+ntp0) = zmelt_d(1:mdim,1:ntp0,it)
            enddo
            !$acc end kernels
            deallocate(zmelt_d) ! write(6,*)'xxxxxxx ia nm1cc nm2ccxxx =',ia, nm1cc,nm2cc,sum(abs(zmelt(ims:ime,nm1cc:nm2cc,ncc+1:ncc+ntp0)))
          endif CoreValence
          ! !(MO) valence-core part NOT TESTED
          if(ncc>0) call rx('m_zmel_init: not yet ncc>0')
         !  ValenceCore: if(nt0 > 0 .and. nccc > 0) then
         !    allocate(zmelt_d(mdim,nt0,nccc))
         !    if (.not.allocated(cphim_d)) allocate(cphim_d(nv,nt0), source = cphim(ias:iae,nmini:nmmax)) ! this may be unnecessary
         !    if (allocated(ppbc_d)) deallocate(ppbc_d)
         !    allocate(ppbc_d(mdim,nv,nccc))
         !    !$acc kernels
         !    do i = 1, mdim
         !      ppbc_d(i,1:nv,1:nccc) = dcmplx(ppb(nc1:ncnv,1:nccc,i,icp))
         !    enddo
         !    !$acc end kernels
         !    ierr = gemm_batch(ppbc_d, cphim_d, zmelt_d, mdim, nt0, nv, nccc, alpha = phasea(ia), sameB = .true.)
         !    !$acc kernels
         !    do itp = 1, nccc
         !      zmelt(ims:ime,nctot+1:nctot+nt0,ics+itp) = zmelt_d(1:mdim,1:nt0,itp)
         !    enddo
         !    !$acc end kernels
         !    deallocate(zmelt_d)
         !  endif ValenceCore
          if (allocated(cphim_d)) deallocate(cphim_d)
          if (allocated(ppbc_d)) deallocate(ppbc_d)
          deallocate(cphiq_d)
        enddo iatomloop
      endblock ZmelWithinMT
      call writemem('mmmmm_zmel111 ngc= '//trim(charext(ngc))//' nm1v nm2v= '//trim(charext(nm1v))//' '//trim(charext(nm2v)))
      flush(stdo)
      ZmelIPWif: if(ngc/=0 .and. nm1v<=nm2v) then
        ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
          use m_read_ppovl,only:igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv,nnxi,nnxe,nnyi,nnye,nnzi,nnze,nggg
          integer:: igcgp2,nn(3), iggg, igp1, itp, igc, igp2 !, igcgp2i_(ngc,ngp2)
          complex(8):: phase(ngc)!zmelp0(ngc,nm1v:nm2v,ntp0)
          allocate( ggitp(ngcgp,ntp0), gggmat(ngcgp,ngp1), ggitp_work(ngc, ngp2),igcgp2i_(ngc,ngp2))
          if(debug) call writemem('mmmmm_zmel111aaa')
          allocate(zmelp0(ngc,nm1v:nm2v,ntp0))
          if(debug) call writemem('mmmmm_zmel111bbb')
#ifdef __GPU
          attributes(device) :: zmelp0, ggitp, gggmat, ggitp_work, igcgp2i_, nn
#endif
          phase(:)=[(exp( -img*tpi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]  !prepared by CPU
          !$acc data copyin(dgeigqk, geigq, phase, ngvecpB1, ngvecpB2, ngveccR, nadd, ggg(1:nggg), nvgcgp2(1:3,1:ngcgp), &
          !$acc             igggi(nxi:nxe,nyi:nye,nzi:nze), igcgp2i(nnxi:nnxe,nnyi:nnye,nnzi:nnze), ppovlinv(1:ngc,1:ngc))
          !$acc kernels loop independent collapse(2)
          if(debug) call writemem('mmmmm_zmel111ccc')
          do igcgp2 = 1, ngcgp
            do igp1 =1, ngp1
              nn(1:3) = ngvecpB1(1:3,igp1) - nvgcgp2(1:3,igcgp2) - nadd(1:3)
              if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
              iggg = igggi(nn(1),nn(2),nn(3))
              if(iggg>=0) gggmat(igcgp2,igp1)=ggg(iggg)
            enddo
          enddo
          !$acc end kernels
          if(debug) call writemem('mmmmm_zmel111ddd')
          !$acc kernels loop independent collapse(2)
          do igc = 1, ngc
           do igp2 = 1, ngp2
              nn(1:3) = ngveccR(1:3,igc) + ngvecpB2(1:3,igp2)
              igcgp2i_(igc,igp2) = igcgp2i(nn(1),nn(2),nn(3))
            enddo
          enddo
          !$acc end kernels
          ierr = gemm(gggmat, geigq(:,itq(nqini_rank:nqmax_rank)), ggitp, ngcgp, ntp0, ngp1, ldB = ngpmx)
          !      2024-8-20 bugfix for sxcf when itq in not contiguous. geigq(1,itq(nqini_rank)) --> geigq(:,itq(nqini_rank:nqmax_rank))
          deallocate(gggmat)
          if(debug) call writemem('mmmmm_zmel111eee')
          do itp = 1, ntp0
            !$acc kernels loop independent collapse(2)
            do igp2 = 1, ngp2
              do igc = 1, ngc
                ggitp_work(igc, igp2) = ggitp(igcgp2i_(igc,igp2), itp)
              enddo
            enddo
            !$acc end kernels
            ierr = gemm(ggitp_work, dgeigqk(1,nm1v), zmelp0(1,nm1v,itp), ngc, nm2v-nm1v+1, ngp2, ldB = ngpmx)
          enddo
          deallocate(ggitp_work,ggitp,igcgp2i_)
          if(debug) call writemem('mmmmm_zmel111fff')
          !$acc kernels
          do igc = 1, ngc
            zmelp0(igc,nm1v:nm2v,1:ntp0) = phase(igc)*zmelp0(igc,nm1v:nm2v,1:ntp0)
          enddo
          !$acc end kernels
          allocate(zmelt_d(ngc,nm1v:nm2v,ntp0))
          if(debug) call writemem('mmmmm_zmel111hhh')
          ierr = gemm(ppovlinv, zmelp0, zmelt_d, ngc, ntp0*nmtot, ngc) 
          !$acc kernels
          zmelt(nbloch+1:nbloch+ngc,nm1v:nm2v,ncc+1:ncc+ntp0) = zmelt_d(1:ngc,nm1v:nm2v,1:ntp0)
          !$acc end kernels
          !$acc end data
          if(debug) call writemem('mmmmm_zmel333')
          deallocate(zmelp0)
          deallocate(zmelt_d)
        endblock ZmelIPW
      endif ZmelIPWif
      if(debug) call writemem('mmmmm_zmel endof ZmelIPWif')
      allocate(zmel(nbb,ns1:ns2,nqtot))
      if(debug) call writemem('mmmmm_zmel deallocate zmel')
      !$acc enter data create(zmel)
      !$acc host_data use_device(zmel)
      !$acc data copyin(ppovlz(1:ngb,1:nbb))
      ierr = gemm(ppovlz, zmelt, zmel(1,ns1,1), nbb, nmtot*ncc, ngb, opA = m_op_C)
      ierr = gemm(ppovlz, zmelt, zmel(1,ns1,ncc+ini_index), nbb, nmtot*ntp0, ngb, opA = m_op_C)
      !$acc end data
      !$acc end host_data
      if (present(comm)) then
        block
          integer, allocatable :: data_disp(:), data_size(:)
          integer :: ini, num, end
          allocate(zmel_buf, mold = zmel)
          allocate(data_size(0:mpi_size-1), data_disp(0:mpi_size-1))
          do irank = 0, mpi_size - 1
            call int_split(nqmax-nqini+1, mpi_size, irank, ini, end, num)
            data_size(irank) = nmtot*nbb*num
            data_disp(irank) = nmtot*nbb*(ini-1)
          enddo
          !GPU aware MPI style does not work ...why??
          !!$acc data create(zmel_buf) present(zmel)
          !!$acc host_data use_device(zmel, zmel_buf)
          !call mpi_allgatherv(zmel(1,ns1,ncc+ini_index), data_size(mpi_rank), mpi_complex16, zmel_buf, data_size, data_disp, &
          !          &  mpi_complex16, comm, mpi_info)
          !!$acc kernels
          !zmel = zmel_buf
          !!$acc end kernels
          !!$acc end host_data
          !!$acc end data
          ! data copy from GPU to CPU for MPI routine. this would be slow
          !$acc update host(zmel)
          call mpi_allgatherv(zmel(1,ns1,ncc+ini_index), data_size(mpi_rank), mpi_complex16, zmel_buf, data_size, data_disp, &
                    &  mpi_complex16, comm, mpi_info)
          zmel = zmel_buf
          !$acc update device(zmel)
          if(debug) call writemem('mmmmm_zmel after mpi=allgatherv')
          deallocate(zmel_buf, data_size, data_disp)
        end block
      endif
      if(zmelconjg) then
        !$acc kernels
        zmel = dconjg(zmel)
        !$acc end kernels
      endif
      if(present(maxmem)) maxmem=memused() ! MaxUsed memory in GB 
      deallocate(zmelt)
    endblock ZmelBlock
  end subroutine get_zmel_init_gemm
end module m_zmel


!This is original of get_zmel_init. Removed at 2024-8-20. 
!   subroutine get_zmel_init(q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg,maxmem) !maxmem is optional Maxused memeory in GB
!     use m_readeigen,only: readcphif 
!     use m_readeigen,only: readgeigf
!     use m_itq,only: itq,ntq
!     intent(in)::           q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg
!     !note: ncc can be =nctot or 0 
!     ! nm1:nm2 is the batch range of middle states (count including core index (i=1,nctot). Thus kth valence is at nctot+k.
!     ! kvec is in the IBZ, rk = Rot_irot(kvec)
!     !! \parameter all inputs
!     !! \parameter matrix <MPB psi|psi>
!     !! igb: index of mixed product basis       at rkvec (or written as rk)
!     !!   igb=1,ngb
!     !!   ngb=nbloch+ngc  ngb: # of mixed product basis
!     !!                   nbloch: # of product basis (within MTs)
!     !!                   ngc: # of IPW for the Screened Coulomb interaction.
!     !!                   igc is for given
!     ! zmelt(itp|it,ib)= ZO^-1 <MPB(rkvec,ngb) phim(q-rkvec, nm1:nm2, ispm)|phiq(q,ncc+nqmax,ispq)> , or dconjg( ) when zmelconjg=T
!     ! zmel = transpose(ppovlz,zmelt(:,it,itp))
!     !
!     ! matmul(symgg(:,:,irot),kvec)-rkvec can have difference of reciprocal vectors.
!     !
!     real(8),parameter::tolq=1d-8
!     complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
!     integer:: isp,nm1,nm2,nqmax,irot,ispq,ispm,nqini, nctot,ncc,ncnv,ncorec,nccc,mdim
!     integer:: it,ia,kx,verbose,nstate,ispqk,nm1c,nm2c,nm1v,nm2v
!     integer:: ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
!     integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime
!     real(8):: quu(3),q(3), kvec(3),rkvec(3),qkt(3),qt(3), qdiff(3)
!     real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass) ! ppb= <Phi(SLn,r-R)_q,isp1 |Phi(SL'n',r-R)_qk,isp2 B_k(S,i,rot^{-1}(r-R))>
!     logical:: iprx,debug=.false.,cmdopt0
!     logical:: zmelconjg
!     integer,allocatable:: ngveccR(:,:)
!     complex(8)::cphiq(nlmto,nband), cphim(nlmto,nband)
!     complex(8):: geigq(ngpmx,nband),dgeigqk(ngpmx,nband)
!     integer:: invr,nt0,ntp0,nmtot,nqtot
!     integer:: iasx(natom),icsx(natom),iatomp(natom),imdim(natom)
!     real(8)::tr(3,natom)
!     real(8)::qk(3),symope(3,3),shtv(3)
!     real(8),optional::maxmem
!     integer::ns1,ns2
!     ! nblocha     = number of optimal product basis functions for each class
!     ! nlnmx     = maximum number of l,n,m
!     ! nctot      = total no. allowed core states
!     ! nbloch     = total no. product basis within MT
!     !           if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
!     !           if(sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
!     !           if(sum(nlnmv(iclass(1:natom)))/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')
!     !           if(sum(ncore(iclass(1:natom)))/= nctot) call rx( "psicb_v3:sum(ncore) wrong")
!     !           if(ncc/=0 .AND. ncc/=nctot) call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
!     !! zmelp(igc(qi),it(qk),itp(q)) = <igc it(for q2+G2) |itp(for q1+G1)>
!     !! NOTE: shtv = g(delta_{g^-1})
!     !!-- zmelp(igc,it,itp) = <igc it(for G2)|itp(for G1)|> matrix element.
!     !!   zmelp0(igc,it,itp)= <Gc' G2|G1> geig(G1,itp) geig^*(G2,it)
!     !!   zmelp(igc,it,itp) = <Gc|Gc'>^-1 zmelp0(Gc',it,itp) 
!     !!   (<Gc|Gc'>^-1 = ppovlinv)
!     !! New ggg matrix <Gc |G1 G2> is introduced.
!     !!
!     !!    <Gc G2|G1> is equivalent to <-Gc+G1-G2>; described by ggg
!     !! Readin input
!     !!    ggg(1:nggg) = <Gc+G2-G1>
!     !!    nvggg(3,1:nggg)   for Gc+G2-G1
!     !!    nvgcgp2(3,ngcgp2) for Gc+G2
!     !!    ppovlinv(ngc,ngc) <Gc|Gc> matrix
!     !
!     ! ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shtv)*exp(-i(q-Gadd)*shtv)
!     ! NOTE: nvgcgp2(:,igcgp2) means symope(Gc)+ G2
!     !========================================================
!     ! NOTE \bfr'= g (\bfr) +\delta_g. Then mapping is ROT[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
!     ! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <Gc'G2|G1> geigq(G1,itp) geigqk*(G2,it) = <Gc' itp(G2)|it(G1)>
!     ! core + valence ordered index: nm1,nm2
!     ! nm1 :starting index of middle state  (nctot+nvalence order)
!     ! nm2 :end      index of middle state  (nctot+nvalence order)
!     nm1=ns1
!     nm2=ns2
!     debug=cmdopt0('--debugzmel')
!     if(allocated(zmel)) deallocate(zmel)
!     nt0  = nm2-nm1+1
!     ntp0 = nqmax-nqini+1
!     nmtot  = nt0 !nctot + nt0     ! = phi_middle nmtot=nm2-nm1+1
!     nqtot  = ncc + ntp0    ! = phi_end
!     if(nm1>nctot) then !skip core 
!        nm1c=0
!        nm2c=-1
!        nm1v=nm1
!        nm2v=nm2
!     elseif(nm2>nctot) then
!        nm1c=nm1
!        nm2c=nctot
!        nm1v=nctot+1
!        nm2v=nm2
!     else !skip valence
!        nm1c=nm1
!        nm2c=nm2
!        nm1v=0
!        nm2v=-1
!     endif   
!     if(debug) write(stdo,ftox)'mmmmmmmmmmmmmm ncc nqtot',ncc, nqtot,'ntp0 nt0 nbb=',ntp0,nt0,nbb
!     if(nmtot<=0.or.nqtot<=0) return
!     qk =  q - rkvec ! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
!     associate(cphitemp=> readcphif(q,ispq))    
!       cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq)) 
!     endassociate
!     cphim = readcphif(qk, ispm) 
!     symope= symgg(:,:,irot)
!     if(ngc/=0) then
!       call readqg('QGpsi',q,     qt, ngp1, ngvecpB1) !q is mapped to qt in BZ
!       call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2)
!       qdiff = matmul(symope,kvec)-qt+qkt ! rkvec + qkt - qt is not zero. <M(rkvec) Phi(qk) |Phi(q)>
!       nadd = nint(matmul(transpose(plat), qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
!       block
!         use m_read_ppovl,only: getppx2, ngvecc,ngcread
!         integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
!         call getppx2(qlat,kvec) ! read and allocate ppovlinv
!         if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')
!         allocate(ngveccR(1:3,1:ngc))
!         call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
!       endblock
!       geigq   = readgeigf(q, ispq) !read IPW part at q   !G1 for ngp1
!       dgeigqk = readgeigf(qk,ispm) !read IPW part at qk  !G2 for ngp2
!       dgeigqk = dconjg(dgeigqk)
!    endif

!     ppb = ppbir(:,:,:,:,irot,ispq)           !MPB has no spin dependence
!     invr  = invg(irot)       !invrot (irot,invg,ngrp) ! Rotate atomic positions invrot*R = R' + T
!     tr    = tiat(:,:,invr)
!     iatomp= miat(:,invr)
!     shtv  = matmul(symope,shtvg(:,invr))
!     imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
!     iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)] !offset counter
!     icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
!     if(debug) write(stdo,ftox)'qqqqqqqq org',nbloch,ngc,nm1,nm2,nqtot
!     ZmelBlock:block
!       complex(8):: zmelt(1:nbloch+ngc,nm1:nm2,nqtot)
!       logical ::oncewrite
! !      real(8)::memused
!       character(8):: charext
!       zmelt=0d0
!       ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
!         complex(8):: phasea(natom) 
!         phasea = [(exp(-img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)] 
!         iatomloop: do concurrent(ia = 1:natom)
!           ic    = iclass(ia)
!           nc    = nlnmc(ic) !nlnmc      = number of l,n,m for core states
!           nv    = nlnmv(ic) !nlnmv      = number of l,n,m for valence
!           nc1   = nc + 1
!           ncnv  = nc+nv  !ncore + nvalence
!           iap   = iatomp(ia)  
!           icp   = iclass(iap)
!           ias   = iasx(ia)  !start of nlnmv for ia
!           iae   = ias+nv-1  !end  
!           ics   = icsx(ia)   
!           ims   = imdim(iap) !start of PB for ia
!           ime   = ims-1+nblocha(icp)
!           mdim  = nblocha(icp)
!           ncorec= merge(ncore(ic),0,nctot>0)
!           nccc  = merge(ncore(ic),0,ncc>0)
!           associate(      &
!                dcphiqk => dconjg(cphim(ias:iae,nm1v:nm2v)),&
!                tdcphiqk=> transpose(dconjg(cphim(ias:iae,nm1v:nm2v))),&
!                cphiq   => cphiq(ias:iae,nqini:nqmax),&
!                ppbc    => ppb(nc1:ncnv,icore(1:ncorec,ic),:,icp)) !all readin but we need only bands nm1: or nqini:
!             ValenceValence:if(nm2v>=nm1v) then           
!                do concurrent(i=1:mdim)    !valence-valence  !=== this may be time-consuming block ==================
!                   zmelt(i-1+ims, nm1v:nm2v,ncc+1:ncc+ntp0)=phasea(ia) &
!                        *matmul(tdcphiqk,matmul(transpose(ppb(nc1:ncnv,nc1:ncnv,i,icp)),cphiq))
!                enddo
!             endif ValenceValence
!             CoreValence: if(nm2c>=nm1c) then
!                do concurrent(it=1:ncorec) !core-valence
!                   if(nm1c<=ics+it.and.ics+it<=nm2c) then
!                      zmelt(ims:ime, ics+it, ncc+1:ncc+ntp0)=phasea(ia)* matmul(transpose(ppbc(:,it,1:mdim)),cphiq(:,1:ntp0))
!                   endif
!                enddo
!                do concurrent(itp=1:nccc)  !valence-core
!                   zmelt(ims:ime, nm1v:nm2v, ics+itp)=phasea(ia)*matmul(transpose(ppbc(:,itp,1:mdim)),dcphiqk(:,nm1v:nm2v)) !1:nt0))
!                enddo                              !^^^^^^^^^phasea(ia) right? 2024-1-7 or phasea(ia) or dcongj(phasea(ia))?
!             endif CoreValence
!           endassociate
!         enddo iatomloop
!       endblock ZmelWithinMT
!       if(debug) call writemem('mmmmm_zmel111 ngc nm1v nm2v '//charext(ngc)//' '//charext(nm1v)//' '//charext(nm2v))
!       if(debug) flush(stdo)
!       if(ngc/=0.and.nm1v <= nm2v)then
!         ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
!           use m_read_ppovl,only: igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv
!           integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2),itqx,inq
!           complex(8):: phase(ngc),ggitp(ngcgp,ntp0),gggp(ngpmx)!NOTE: ngcgp is the index for (QpGcou+QpGphi). See sugw.f90
!           complex(8),allocatable:: zmelp0(:,:,:),ggitp_(:,:) !,gggmat(:,:)
!           do concurrent(igcgp2=1:ngcgp)
!              gggp=0d0
!              do igp1=1,ngp1 !G synthesized 
!                 nn = ngvecpB1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd.  Note that -Gadd= -rk + qt -qkt
!                 if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
!                 iggg = igggi(nn(1),nn(2),nn(3)) !inversion table
!                 if(iggg>=0) gggp(igp1)=ggg(iggg)
!              enddo
!              ggitp(igcgp2,:)= [(sum(gggp(1:ngp1)*geigq(1:ngp1,itq(inq))),inq=nqini,nqmax)]
!           enddo
!           do concurrent (igc=1:ngc,igp2=1:ngp2) !igc for B !G synthesized 
!             nn = ngveccR(:,igc) + ngvecpB2(:,igp2)
!             igcgp2i_(igc,igp2)=igcgp2i(nn(1),nn(2),nn(3))
!           enddo
!           if(debug) call writemem('mmmmm_zmel222 Mems '//&
!                ftof(4*size(igcgp2i_)/kk**3+4*size(ngvecpB2)/kk**3+4*size(ngvecpB1)/kk**3,3)//' '//ftof(16*size(ggitp)/kk**3,3))
!           if(debug) flush(stdo)
!           phase(:)=[(exp( -img*tpi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]
!           !               geigq   => readgeigf(q, ispq),&         !read IPW part at q   !G1 for ngp1
!           !               dgeigqk => dconjg(readgeigf(qk,ispqk))) !read IPW part at qk  !G2 for ngp2
!           allocate(zmelp0(ngc,nm1v:nm2v,ntp0),ggitp_(ngc,ngp2))
!           do concurrent (itp= 1:ntp0) !=== this may be time-consuming block (or maynot)==================
!              forall(igc=1:ngc,igp2=1:ngp2) ggitp_(igc,igp2)= ggitp(igcgp2i_(igc,igp2),itp)
!              zmelp0(:,:,itp)= matmul(ggitp_,dgeigqk(1:ngp2,nm1:nm2))
!           enddo
!           deallocate(ggitp_)
!           forall(igc=1:ngc) zmelp0(igc,:,:)=phase(igc)*zmelp0(igc,:,:) 
!           call matm(ppovlinv,zmelp0,zmelt(nbloch+1:nbloch+ngc,nm1:nm2,ncc+1:ncc+ntp0),ngc,ngc,ntp0*nt0)
!           if(debug) call writemem('mmmmm_zmel333')
!           deallocate(zmelp0)
!         endblock ZmelIPW
!       endif
!       allocate(zmel(nbb,nm1:nm2, nqtot))
!       call matm(dconjg(transpose(ppovlz)), zmelt, zmel,nbb,ngb,nmtot*nqtot) !MultiplePPOVLZ
!       if(present(maxmem)) maxmem=memused() ! MaxUsed memory in GB 
!       if(zmelconjg) zmel=dconjg(zmel)
!     endblock ZmelBlock
!   end subroutine get_zmel_init

