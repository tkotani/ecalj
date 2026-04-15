!>  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
module m_uumat
  public uumatrix
  private
contains
subroutine uumatrix()
  ! ixc=2: <u(k) | u(k+b)>
  ! ixc=3: <u(k) | u(k+q0)>
  ! ixc=4: <u(k) | u(k+q)>
  ! Takashi Miyake, Mar 2008, parallelized.  originally written by Takao Kotani, April, 2004
  use m_readqg,only: readngmx, readqg0, readqg
  use m_hamindex,only: Readhamindex, ngrp, symops
  use m_readeigen,only: init_readeigen, init_readeigen2, readeval
  use m_readeigen, only: readcphif_mpi, readgeigf_mpi
  use m_read_bzdata,only: read_bzdata, nqbz, nqibz, qbas=>qlat, qibz, qbz, nq0i_=> nq0i, nq0i=>nq0ix, q0i
  use m_genallcf_v3,only: genallcf_v3, natom, nspin, nl, nn, plat, pos, alat, nindx, ndima, nqbzt, nband, nspc, nspx, lmxa
  use m_keyvalue,only: getkeyvalue
  use m_pwmat,only: mkppovl2
  use m_ll,only: ll
  use m_mpi,only: mpi__broadcast, mpi__root, mpi__size, mpi__rank, comm, mpi__allreducesum, mpi__reducesum
  use m_lgunit,only: m_lgunit_init, stdo
  use m_setqibz_lmfham,only: set_qibz, irotg
  use m_mlo_ham, only: mlo_read_hma_rs => read_ham_rs, mlo_nwf => ndimMTO
  use m_mlo_scrw, only: mlo_nnwf_init => nnwf_init, mlo_nnwf => nnwf
  use m_mlo_wfs, only: cmlo_init, get_geig_cmlo, get_cphi_cmlo
  use m_mpiio, only: openm, writem, closem
  use m_lmfinit,only: m_lmfinit_init
  use m_lattic,only: m_lattic_init
  use m_mksym,only: m_mksym_init
  use m_mpitk, only: m_mpitk_init
  use,intrinsic :: ieee_arithmetic
  use m_ftox
  implicit none
  integer:: i,ix,ngrpx ,is, nxx ,ibas ,ibas1, ngpmx, ifphi, nbas, nradmx, ncoremx, &
            nrx, ic, icx, isp, l, n, irad, ifoc, ldim2, ixx, ngp1, ngp2, &
            ixc, nbbloop, ifuu(2), ifbb, nbb, iko_ixs(2), iko_fxs(2), &
            iqbz, ibb, itmp, iqb, ibb2, iqtmp, ibbtmp, ndg1(3),ndg2(3), &
            j1, j2, j1max, j1min, ispin ,l1, l2, lm1, lm2, ibas2, lm3, ir, ia1, ia2, m2, l3, m1, lxx, &
            ico, lx, ierr, n1, n2, ii, timevalues(8), ib,  ie, ioc, iog, ispc
  integer, allocatable :: ngvecpB(:,:,:),ngveccB(:,:), ngvecpf1(:,:), ngvecpf2(:,:),nx(:,:),nblocha(:),ifppb(:)
  integer, allocatable :: ncindx(:,:), lcindx(:,:), nrad(:), nindx_r(:,:), lindx_r(:,:), nc_max(:,:), &
                          m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:), ikbidx(:,:), ncore(:)
  real(8), allocatable :: ppbrd (:,:,:,:,:,:,:), cg(:,:,:), symope(:,:), phij(:), psij(:), rprodx(:,:), rphiphi(:)
  real(8), allocatable :: phitoto(:,:,:,:,:), aa(:), rr(:,:), phitotr(:,:,:,:,:), bb(:), zz(:), rmax(:), cy(:), yl(:)
  real(8), allocatable :: bbv(:,:), eval1(:), eval2(:)
  real(8) :: q1(3), q2(3), dq(3), absqg2, absdq, r2s, absqg, ylk
  real(8) :: ef, q(3), q1x(3),q2x(3)
  real(8),parameter ::  pi = 4d0*atan(1d0),fpi = 4d0*pi
  complex(8),parameter :: img=(0d0,1d0)
  complex(8),allocatable :: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:), uum(:,:,:), ppovl(:,:), ppj(:,:,:,:)
  complex(8) :: phaseatom
  logical :: cmdopt2, cmdopt0
  logical :: use_bbvec_file, use_mlo, spinflip
  complex(8), allocatable :: uumq(:,:,:,:)
  character(8) :: head(2:3,2)
  character(4) charnum4
  character*7:: charnum7
  character(20):: outs=''
  procedure(readgeigf_mpi), pointer :: get_geig => readgeigf_mpi
  procedure(readcphif_mpi), pointer :: get_cphi => readcphif_mpi
  call M_lgunit_init()
  call m_MPItk_init(comm)
  !for rotMTO
  call m_lmfinit_init('uumat',comm)! Read ctrlp into module m_lmfinit.
  call m_lattic_init()       ! lattice setup (for ewald sum)
  call m_mksym_init()  !symmetry go into m_lattic and m_mksym

  call date_and_time(values=timevalues)
  write(stdo,"('mpirank=',i5,' YYYY.MM.DD.HH.MM.msec=',9i4)")mpi__rank,timevalues(1:3),timevalues(5:8)
  if(mpi__root) then
    if(cmdopt2('--job=',outs)) then
      read(outs,*) ixc
    else
      write(stdo,*) ' --- Choose modes below -------------------'
      write(stdo,*) '  (2) (q,q+b), (3) (q,q+q0), (4) sum_k (k,k+q) with spinflip'
      write(stdo,*) ' --- Put number above ! ------------'
      read(5,*) ixc
      write(stdo,*) ' ixc=', ixc !computational mode index
    endif
  endif
  call MPI__Broadcast(ixc)
  if(.not.(ixc == 2.or. ixc==3 .or. ixc==4))call rx('main_huumat_MPI: ixc error')
  use_bbvec_file = (ixc /= 4)
  use_mlo         = (ixc == 4)
  spinflip      = (ixc == 4)
  call read_BZDATA()
  if (mpi__root) write(stdo,*)' ======== nqbz nqibz ngrp, nq0i, nq0i_=',nqbz,nqibz,ngrp, nq0i, nq0i_
  call genallcf_v3(incwfx=0) !readin condition. use ForX0 for core in GWIN !  call Readhbe()    !Read dimensions of h,hb
  call getsrdpp2(natom,nl,nxx)    ! --- read by rdpp ; Radial integrals ppbrd and plane wave part
  call readngmx('QGpsi',ngpmx)

  ReadPHIVC: block
    open(newunit=ifphi,file='__PHIVC',form='unformatted')     ! PHIV+PHIC augmentation wave and core
    read(ifphi) nbas, nradmx, ncoremx,nrx
    if(ixc == 4 .and. nspin /=2) call rx('ixc == 4 works only nspin=2')
    if(ixc == 4 .and. mpi__root) write(stdo,ftox) '!!WARNING ixc ==4 requires --phispinsym'
    if(nqbz  /= nqbzt) call rx( ' hx0fp0: nqbz /=nqbzt  in hbe.d')
    if(nbas  /= natom) call rx(' nbas(PHIVC) /= natom ')
    allocate(  ncindx(ncoremx,nbas), lcindx(ncoremx,nbas), &
      nrad(nbas), nindx_r(1:nradmx,1:nbas), lindx_r(1:nradmx,1:nbas), &
      aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) , &
      phitoto(nrx,0:nl-1,nn,nbas,nspin), &
      phitotr(nrx,0:nl-1,nn,nbas,nspin), &
      nc_max(0:nl-1,nbas),ncore(nbas),rmax(nbas) )
    read(ifphi) nrad(1:nbas)
    read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
    nc_max=0
    do ibas=1,nbas
      ic = ibas
      read(ifphi) ncore(ic), ncoremx                            !core
      read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
      read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
      if(ic/=icx) call rx(' h_uu: ic/=icx')
      read(ifphi) rr(1:nrofi(ic),ic)
      rmax(ic) = rr(nrofi(ic),ic)
      do isp = 1, nspin
        if (mpi__root)  write(stdo,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
        do ico = 1, ncore(ic) !core
          l =  lcindx(ico,ic)
          n =  ncindx(ico,ic)
          read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
          phitotr(1:nrofi(ic),l,n, ic,isp)=  phitoto(1:nrofi(ic),l,n, ic,isp) ! core raw= core orthgonal
          if(n>nc_max(l,ic)) nc_max(l,ic)=n
        enddo
        do irad = 1, nrad(ic)   !valence
          l = lindx_r (irad,ic)
          n = nindx_r (irad,ic) + nc_max(l,ic)
          read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
          read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
        enddo
      enddo
    enddo
    close(ifphi)
  endblock ReadPHIVC

  ngrpx=1
  allocate( cg(nl**2,nl**2,(2*nl-1)**2),source=0d0)
  allocate( symope(3,3),source=reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],shape=[3,3])) !symops identicay matrix with ng=1
  call rotcg(nl-1,symope,ngrpx,cg) !CG coefficient
  call Readhamindex()
  call init_readeigen()   !Initialization for readeigen
  call init_readeigen2()
  if(use_mlo) then
    call mlo_read_hma_rs() !set nwf
    call mlo_nnwf_init(nnwf_size_reduction=.true.)
    get_geig => get_geig_cmlo
    get_cphi => get_cphi_cmlo
    call cmlo_init()
  endif
  call readngmx('QGpsi',ngpmx) !max number of the set q+G
  if(use_mlo) then
    allocate(geig1(ngpmx*nspc,mlo_nwf),geig2(ngpmx*nspc,mlo_nwf))
    allocate(cphi1(ndima*nspc,mlo_nwf),cphi2(ndima*nspc,mlo_nwf))
  else
    allocate(geig1 (ngpmx*nspc,nband),geig2(ngpmx*nspc,nband),eval1(nband),eval2(nband))
    allocate(cphi1 (ndima*nspc,nband),cphi2(ndima*nspc,nband) )
  endif

  open(newunit=ifoc,file='@MNLA_CPHI')
  ldim2 = ndima
  read(ifoc,*)
  allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
  do ix =1,ndima
    read(ifoc,*) m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx !m,m,l,ibas index 
    if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
  enddo
  close(ifoc)
  if(mpi__root) then
    if(ixc==4) then
    write(stdo,*) ' Used k number in Q0P =', nq0i_
    write(stdo,"(i3,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i_)
    else
    write(stdo,*) ' Used k number in Q0P =', nq0i
    write(stdo,"(i3,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
    endif
  endif

  if(use_bbvec_file) then
    Readbbvec: block
      integer :: nspin2, nqbz2, itmp2
      open(newunit=ifbb,file='BBVEC')
      read(ifbb,*)
      read(ifbb,*)nbb, nqbz2
      if (nqbz /= nqbz2) call rx('readbb: nqbz is wrong!')
      allocate(bbv(3,nbb),ikbidx(nbb,nqbz))    !call readbb(ifbb,nqbz,nspin,nbb, bbv, ikbidx, iko_ixs,iko_fxs,noxs)
      do i = 1,nbb
        read(ifbb,*) bbv(1:3,i)
      enddo
      do iqbz = 1,nqbz
        read(ifbb,*) !itmp,u(1:3)
        do ib = 1,nbb
          read(ifbb,*)itmp,itmp2,ikbidx(ib,iqbz) !,u(1:3)
        enddo
      enddo
      read(ifbb,*)
      read(ifbb,*)nspin2
      if(nspx /= nspin2) call rx('nspin is wrong!')
      do is = 1,nspx
        read(ifbb,*)iko_ixs(is),iko_fxs(is)
      enddo
      close(ifbb)
    endblock Readbbvec
  else !bbv and nbb are not used
    iko_ixs(1:nspx) = 1
    iko_fxs(1:nspx) = mlo_nwf
  endif

  if(ixc ==2 .or. ixc == 3) then
    head(2,1:2)=['UUU.','UUD.']
    head(3,1:2)=['UUq0U.','UUq0D.']
    if(mpi__root) then
      do isp=1,nspx
        if(cmdopt0('--ahc')) then
          open(newunit=ifuu(isp),file=trim(head(ixc,isp))//charnum7(0),form='unformatted')
        else
          open(newunit=ifuu(isp),file=trim(head(ixc,isp))//charnum4(0),form='unformatted')
        endif
        if(ixc==2)then
          write(ifuu(isp))'nqbz,nbb,iko_ixs(isp),iko_fxs(isp)',isp
          write(ifuu(isp))nqbz,nbb,iko_ixs(isp),iko_fxs(isp)
        elseif(ixc==3) then
          write(ifuu(isp))'nqbz,nq0i,iko_ixs(isp),iko_fxs(isp)',isp
          write(ifuu(isp))nqbz,nq0i,iko_ixs(isp),iko_fxs(isp)
        endif
        close(ifuu(isp))
      enddo
    endif
  endif
  ! --- Set q1(j1range) q2(j2range); Note that the true q when we generate eigenfunctions are q1x and q2x.
  ! q1-q1x should be a G vector.  So you may need to take into account the phase shift to <u|u> vectors.
  j1min = minval(iko_ixs(1:nspx)) !starting band index
  j1max = maxval(iko_fxs(1:nspx))
  allocate( uum(j1min:j1max, j1min:j1max,nspx) ) ! uumatrix allocated
  if(cmdopt0('--qibzonly')) call set_qibz(plat,qbz,nqbz,symops,ngrp) !If only at qibz, we need to set irotg
  if (ixc == 2) nbbloop = nbb
  if (ixc == 3) nbbloop = nq0i
  if (ixc == 4) nbbloop = nq0i_ ! same with nq0i

  lxx=2*(nl-1)
  allocate(ppj(ndima,ndima,nspin,nbbloop), source = (0d0,0d0))
  allocate(ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas), rprodx(nrx,0:lxx), phij(0:lxx), psij(0:lxx), rphiphi(nrx))
  allocate(cy((lxx+1)**2),yl((lxx+1)**2))

  ibbloop0: do ibb = 1,nbbloop
    if(ixc == 2) dq=-bbv(:,ibb)  !q1(:) = qbz(:,iqbz)        !q2(:) = qbz(:,iqbz) + bbv(:,ibb)
    if(ixc == 3) dq=-q0i(:,ibb) !q1(:) = qbz(:,iqbz)         !q2(:) = qbz(:,iqbz) + q0i(:,ibb)
    if(ixc == 4) dq=-q0i(:,ibb) !q1(:) = qbz(:,iqbz)         !q2(:) = qbz(:,iqbz) + qbz(:,ibb)
    if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
    if(cmdopt0('--q2q1test')) dq=1d-10
    absdq = sqrt(sum(dq**2))
    absqg2 = (2*pi/alat)**2 *sum(dq**2)
    absqg =sqrt(absqg2)
    call sylmnc(cy,lxx)
    call sylm(dq/absdq,yl,lxx,r2s) !spherical factor Y(dq)

    ppbrd=0d0
    ibasloop0: do ibas = 1,nbas ! radial integral  ppbrd = <phi phi j_l>
      ic = ibas
      do ir =1,nrofi(ic)
        call bessl(absqg2*rr(ir,ibas)**2,lxx,phij,psij) ! phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
        do lx = 0, lxx
          rprodx(ir,lx) = merge(0d0,rr(ir,ibas)* phij(lx)* (absqg*rr(ir,ibas))**lx,rr(ir,ibas)==0d0)
          !             = r \times j_l(|dq|r)  !bessel function for the expansion of exp(i(q1-q2) r)
        enddo
      enddo
      ispinloop00: do isp = 1,nspin
        do  l1 = 0, lmxa(ic)
          do  l2 = 0, lmxa(ic)
            do  n1 = 1, nindx(l1+1,ic)
              do  n2 = 1, nindx(l2+1,ic)
                rphiphi(1)       = 0d0
                rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)*phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:nrofi(ic),ic) ! phi = u = r \phi
                do lx = 0, 2*(nl-1)
                  if(lx <abs(l1-l2) .OR. l1+l2<lx) cycle
                  call gintxx( rprodx(1,lx), rphiphi,aa(ic),bb(ic),nrofi(ic), ppbrd(l1, n1,l2, n2, lx, isp,ibas) )
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo ispinloop00
   enddo ibasloop0
 
    ! Calcuate <u{q1x j1} | u_{q2x j2}> = < psi^*{q1x j1} exp(i(q1x-q2x)r) psi_{q2x j2} >
    ! Note that exp(i(q1x-q2x)r) is expanded in the spherical bessel function within MT.
    ! MT part ldim2=ndima; n_indx(1;ldim2):n(phi=1 phidot=2 localorbital=3); l_indx(1:ldim2):l index ; ibas_indx(1:ldim2):ibas index.
    ispinloop02: do ispin=1,nspin
      ia1loop: do 10201 ia1 = 1,ndima
        ibas1= ibas_indx(ia1)
        ia2loop: do 10101 ia2 = 1,ndima
          ibas2= ibas_indx(ia2)
          if(ibas2/=ibas1) cycle
          l1=l_indx(ia1); m1=m_indx(ia1); n1=n_indx(ia1)+ nc_max(l1,ibas1); lm1= l1**2+l1+1+ m1
          l2=l_indx(ia2); m2=m_indx(ia2); n2=n_indx(ia2)+ nc_max(l2,ibas2); lm2= l2**2+l2+1+ m2
          phaseatom = exp( img* 2d0*pi*sum(dq*pos(:,ibas1)) )
          do lm3= (l1-l2)**2+1, (l1+l2+1)**2 ! l3 takes |l1-l2|,...l1+l2
            l3 = ll(lm3)
            ylk= cy(lm3)*yl(lm3)
            ppj(ia1,ia2,ispin,ibb) = ppj(ia1,ia2,ispin,ibb)&
                 + ppbrd(l1,n1,l2,n2,l3,ispin,ibas1) *cg(lm1,lm2, lm3) * fpi* img**l3* phaseatom * ylk
            ! cg(lm1,lm2,lm3)= \int Y_lm3(\hat(r)) Y_lm2(\hat(r)) Y_lm1(\hat(r)) \frac{d \Omega}{4\pi}.
            ! This is based on inverse expansion. See the book of angular momentum book of Rose.Eq.3.8.
         enddo
10101    enddo ia2loop
10201  enddo ia1loop
    enddo ispinloop02
  enddo ibbloop0
  deallocate(ppbrd, rprodx, phij, psij, rphiphi, cy, yl)

  if(allocated(uumq)) deallocate(uumq)
  if(ixc == 4) allocate(uumq(mlo_nwf,mlo_nwf,nq0i_,nspx), source = (0d0,0d0))

  iqbz4uum: do 1070 iqbz = 1,nqbz  !qibzonly need to be improved to balance load in ranks.
    if(mod(iqbz-1,mpi__size)/=mpi__rank) cycle !MPI
    if (cmdopt0('--qibzonly')) then
      if(irotg(iqbz)/=1)  cycle !only irreducible q point
    endif
    if(ixc == 2 .or. ixc == 3)  then
    do isp=1,nspx
      if(cmdopt0('--ahc')) then
        open(newunit=ifuu(isp),file=trim(head(ixc,isp))//charnum7(iqbz),form='unformatted')
      else
        open(newunit=ifuu(isp),file=trim(head(ixc,isp))//charnum4(iqbz),form='unformatted')
      endif
    enddo
    endif
    ibbloop: do 1080 ibb = 1,nbbloop
      if(ixc == 2) then
        iqb = ikbidx(ibb,iqbz)
        q1(:) = qbz(:,iqbz)
        if(cmdopt0('--q2q1test')) then
          q2(:) = qbz(:,iqbz)
        else
          q2(:) = qbz(:,iqbz) + bbv(:,ibb)
        endif
        if (iqb < iqbz) then
          iqtmp = iqb
          do ibb2 = 1,nbb
            itmp = ikbidx(ibb2,iqtmp)
            if (itmp == iqbz) then
              ibbtmp = ibb2
              goto 1200
            endif
          enddo
          call rx('huumat: (iq,ib) error')
1200      continue
          do ispin = 1,nspx
            write(ifuu(ispin))-20
            write(ifuu(ispin))iqbz,ibb,iqtmp,ibbtmp
          enddo
          cycle !we may quit ibbloop here
        endif
      elseif (ixc == 3 .or. ixc == 4) then
        q1(:) = qbz(:,iqbz)
        q2(:) = qbz(:,iqbz) + q0i(:,ibb)
      endif
      call readqg0('QGpsi',q1,  q1x, ngp1) ! write(stdo,"('uuuiq q1 q1x=',3f9.4,3x,3f9.4,i5)") q1,q1x,ngp1
      call readqg0('QGpsi',q2,  q2x, ngp2) ! write(stdo,"('uuuiq q2 q2x=',3f9.4,3x,3f9.4,i5)") q2,q2x,ngp2
      allocate( ngvecpf1(3,ngp1), ngvecpf2(3,ngp2), ppovl(ngp1,ngp2) )
      call readqg('QGpsi',q1, q1x, ngp1, ngvecpf1)
      call readqg('QGpsi',q2, q2x, ngp2, ngvecpf2)
      ndg1= nint(matmul((q1x-q1),plat(:,:)))
      ndg2= nint(matmul((q2x-q2),plat(:,:)))
      do i = 1,3
        ngvecpf1(i,1:ngp1) = ngvecpf1(i,1:ngp1) + ndg1(i)
        ngvecpf2(i,1:ngp2) = ngvecpf2(i,1:ngp2) + ndg2(i)
      enddo
      call mkppovl2(alat,plat,qbas, ngp1,ngvecpf1, ngp2,ngvecpf2, nbas,rmax,pos, ppovl) !--- ppovl= <P_{q1+G1}|P_{q2+G2}>
      ispinloop2: do 1050 ispin=1,nspx !note that nspx=nsp/nspc where nspc=2 for lso=1 (nspc=1 for lso=0,2)
        ii = iko_ixs(ispin)
        ie = iko_fxs(ispin)
        cphi1 = get_cphi(q1,ispin) ! MT part of eigenfunctions
        cphi2 = get_cphi(q2, merge(3-ispin, ispin, spinflip))
        geig1 = get_geig(q1,ispin) ! IPW part of eigenfunctions
        geig2 = get_geig(q2, merge(3-ispin, ispin, spinflip))
        eval1 = readeval(q1,ispin) !eigenvalue at q1
        eval2 = readeval(q2,ispin)
        uum(:,:,ispin) = 0d0
        do ispc=1,nspc ! For lso=0 or 2,ispin=1,nsp. For lso=1, ispin=1 ispc=1,2 nspc=2 
          ioc=(ispc-1)*ndima
          iog=(ispc-1)*ngpmx
          uum(ii:ie,ii:ie,ispin) =uum(ii:ie,ii:ie,ispin) &
               + matmul(transpose(dconjg(cphi1(ioc+1:ioc+ndima,ii:ie))),&
               matmul(ppj(1:ndima,1:ndima,ispc,ibb),cphi2(ioc+1:ioc+ndima,ii:ie))) &! MT part
               + matmul(dconjg(transpose(geig1(iog+1:iog+ngp1, ii:ie))), matmul(ppovl,geig2(iog+1:iog+ngp2,ii:ie))) !IPW part
        enddo
        if(ixc/=4) write(ifuu(ispin)) -10 !dummy
        if(ixc==2) write(ifuu(ispin)) iqbz,ibb,ikbidx(ibb,iqbz)
        if(ixc==3) write(ifuu(ispin)) iqbz,ibb
        if(ixc==2 .or. ixc == 3) write(ifuu(ispin)) ((uum(j1,j2,ispin),j1=ii,ie),j2=ii,ie)
        if(ixc==4) uumq(:,:,ibb,ispin) = uumq(:,:,ibb,ispin) + uum(:,:,ispin)
        if(ixc==4) cycle
        checkwirte: block
          do j1=ii,ie
            do j2=ii,ie !; do j2=j2min,j2max !checkwrite  !if(j1==j2)
              if(eval1(j1)>1d10.or.eval2(j2)>1d10) cycle ! see sugw.f90. padding by huge number
!            write(stdo,ftox)'uumatrix: iq isp j1 j2 q1 q2 <uu>= ',iqbz,ispin,j1,j2,ftof(q1,4),ftof(q1-q2,4),ftof(uum(j1,j2),4),ftof(abs(uum(j1,j2)))
            ! write(stdo,ftox)'uumatrix: iq isp j1 j2 q1 q2 <uu>= ',iqbz,ispin,j1,j2,ftof(q1,4),ftof(q1-q2,4),ftof(uum(j1,j2,ispin),4),'abs',ftof(abs(uum(j1,j2,ispin)))
            write(stdo,ftox)'uumatrix: iq isp j1 j2 q1 q2 <uu>= ',iqbz,ispin,j1,j2,ftof(q1,4),ftof(q1-q2,4),uum(j1,j2,ispin),'abs',abs(uum(j1,j2,ispin))
          enddo ; enddo
        endblock checkwirte
1050  enddo ispinloop2
      deallocate(ngvecpf1, ngvecpf2, ppovl)
      ! write(stdo,*) !'============ result --- diagonal --- ==============',nspx,j1min,j1max,j2min,j2max
1080 enddo ibbloop
    if(ixc/=4) then
    close(ifuu(1))
    if(nspin==2) close(ifuu(2))
    endif
1070 enddo iqbz4uum

  if(ixc==4) then
    call mpi__allreducesum(uumq, size(uumq))
    uumq(:,:,:,:)  = uumq(:,:,:,:)/dble(nqbz)
    if(mpi__root) then
      block
        character(64) :: datfile
        integer :: ifile(2), iq, ifile_handle, recl
        recl = 16*mlo_nwf**2
        do iq=1, nq0i_
          write(stdo,ftox) 'q, diag sum uumq(updw,dwup)/nwf:', q0i(:,iq), &
          (sum([(uumq(i,i,iq,isp),i=1,mlo_nwf)])/dble(mlo_nwf),isp=1,nspin)
        enddo
        ! Info file (sequential): header + q-point list
        open(newunit=ifile_handle, file='__MLOFormFactorQ.info', form='unformatted', status='replace')
        write(ifile_handle) mlo_nwf, nqbz, nspin, nq0i_
        write(ifile_handle) q0i(:,1:nq0i_)
        close(ifile_handle)
        ! Data files (direct-access): one record per q-point
        do isp=1, nspin
          if(isp==1 .and. (.not.spinflip)) datfile = '__MLOFormFactorQ.UP'
          if(isp==2 .and. (.not.spinflip)) datfile = '__MLOFormFactorQ.DN'
          if(isp==1 .and. spinflip) datfile = '__MLOFormFactorQ.UPDN'
          if(isp==2 .and. spinflip) datfile = '__MLOFormFactorQ.DNUP'
          open(newunit=ifile_handle, file=trim(datfile), form='unformatted', access='direct', recl=recl, status='replace')
          do iq=1, nq0i_
            write(ifile_handle, rec=iq) uumq(:,:,iq,isp)
          enddo
          close(ifile_handle)
        enddo
        write(stdo,ftox) 'uumq written to __MLOFormFactor.info / __MLOFormFactorQ.UP / .DN'
      endblock
    endif
  endif
  deallocate(uum)
  if (mpi__root) write(stdo,*)'====== end ========================================'
  call mpi_finalize(ierr)
end subroutine uumatrix

end module m_uumat
