!>  Calculate <u|u> matrix . u_kj(r) is the perodic part of eigencuntion.
subroutine h_uumatrix()
  ! ixc=2: <u(k) | u(k+b)>
  ! ixc=3: <u(k) | u(k+q0)>
  ! Takashi Miyake, Mar 2008, parallelized.  originally written by Takao Kotani, April, 2004
  use m_readqg,only: readngmx,ngcmx,readqg0,readqg
  use m_hamindex,only:   Readhamindex,ngrp,symops
  use m_readeigen,only:init_readeigen,init_readeigen2,readcphif,readgeigf,readeval
  use m_read_bzdata,only: read_bzdata, nqbz,nqibz,nqbzw,nteti,ntetf,qbas=>qlat, ginv, &
    dq_,wbz,qibz,wibz,qbzw, qbz, idtetf,ib1bz,idteti, nstar,irk,nstbz,  nq0i=>nq0ix,q0i
  use m_genallcf_v3,only: genallcf_v3, ncore2=>ncore,nrxx=>nrx, natom,nclass,nspin,nl,nn,nnv,nnc, &
       nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,plat,pos,alat,nindx,& !nlmto,
       nprecb,mrecb,mrece,ndima,nqbzt,nband,mrecg,nspc
  use m_keyvalue,only: getkeyvalue
  use m_pwmat,only: mkppovl2
  use m_ll,only: ll
!  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_mpi,only: mpi__broadcast,mpi__root, mpi__size,mpi__rank,mpi__initialize
  use m_lgunit,only: m_lgunit_init,stdo
  use m_setqibz_lmfham,only: set_qibz,irotg
  use m_ftox
  implicit none
  integer:: nw_input, i,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx &
    ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp, &
    ifepscond,nxx,ifvxcpout,ifgb0vec,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx &
    ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas &
    ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0 &
    ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode,ngpmx
  integer:: nwin,incwfin,verbose,ifphi,nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc, ldim2,ixx,ngp1,ngp2,nq0it
  integer:: iqindx,nqbandx,nqband,j1min_c(2),j1max_c(2),nbmin,nbmax, nmin,nmax,iq2,ntmp,if99,ifile_handle
  integer:: ixc,idummy,idummy2,i1,i2,i3,nbbloop, ifq0p,ifuu(2), ifbb,nbb,iko_ixs(2),iko_fxs(2), &
    iqibz,iqbz,ibb,itmp,itmp2,nqibz2,nqbz2,iqb,ibb2,iqtmp,ibbtmp,ndg(3),ndg1(3),ndg2(3), &
    nb1d,iq0i,j1,j2,j1max,j2max,j1min,j2min,ispin ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx &
    ,ico,lxd,lx, ierr,iclose,input3(3),n1,n2,ig, nproc1,nproc2,nq_proc,ii,jj,kk,iftmp,if101,&
    timevalues(8),ib,nspin2,ie,ioc,iog,ispc
  integer,allocatable :: ngvecpB(:,:,:),ngveccB(:,:), ngvecpf1(:,:), ngvecpf2(:,:),nx(:,:),nblocha(:),ifppb(:)
  integer,allocatable:: ncindx(:,:), lcindx(:,:), nrad(:), nindx_r(:,:), lindx_r(:,:), nc_max(:,:), &
    m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
  integer,allocatable:: ikidx(:),ikbidx(:,:), ibidx(:,:),ibidxs(:,:),ibidx0(:,:,:), ij1idx(:),ij2idx(:)
  integer,allocatable:: ncore(:),iq_proc(:)
  real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:), phij(:),psij(:),rprodx(:,:),rphiphi(:)
  real(8):: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg, ylk
  real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini, dwry,dwh,omg_c,omg2,xxx
  real(8):: ef,q(3),  qgbin(3),qx(3), qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
  real(8):: uunorm,dqx(3),dqx0(3),dq0(3),dg(3),dqmin(3),adq0, q0wf(3),wgt, emin,emax,rydberg
  real(8),parameter::  pi = 4d0*atan(1d0),fpi =    4d0*pi
  real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:) ,phitotr(:,:,:,:,:), bb(:),zz(:),rmax(:),cy(:),yl(:)
  real(8),allocatable :: bbv(:,:), qbandx(:,:),qband(:,:),eband(:)
  complex(8),parameter:: img=(0d0,1d0)
  complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:) ,uum(:,:,:), ppovl(:,:),ppj(:,:,:,:)
  complex(8):: phaseatom
  logical:: qbzreg, lbnds,cmdopt2,cmdopt0
  character(8) :: xt,head(2:3,2)
  character(4) charnum4
  character(20):: outs=''
  call MPI__Initialize()
  call M_lgunit_init()
  call date_and_time(values=timevalues)
  write(stdo,"('mpirank=',i5,' YYYY.MM.DD.HH.MM.msec=',9i4)")mpi__rank,timevalues(1:3),timevalues(5:8)
  if(mpi__root) then
    if(cmdopt2('--job=',outs)) then
      read(outs,*) ixc
    else
      write(stdo,*) ' --- Choose modes below -------------------'
      write(stdo,*) '  (2) (q,q+b), (3) (q,q+q0)'
      write(stdo,*) ' --- Put number above ! ------------'
      read(5,*) ixc
      write(stdo,*) ' ixc=', ixc !computational mode index
    endif
  endif
  call MPI__Broadcast(ixc)
  if(.not.(ixc == 2.or. ixc==3))call rx('main_huumat_MPI: ixc error')
  call read_BZDATA()
  if (mpi__root) write(stdo,*)' ======== nqbz nqibz ngrp=',nqbz,nqibz,ngrp
  call genallcf_v3(incwfx=0) !readin condition. use ForX0 for core in GWIN !  call Readhbe()    !Read dimensions of h,hb
  call getsrdpp2(nclass,nl,nxx)    ! --- read by rdpp ; Radial integrals ppbrd and plane wave part
  call readngmx('QGpsi',ngpmx)
  open(newunit=ifphi,file='PHIVC',form='unformatted')     ! PHIV+PHIC augmentation wave and core
  read(ifphi) nbas, nradmx, ncoremx,nrx
  if(nclass/= natom) call rx(' nclass /= natom ') !WE ASSUME iclass(iatom)= iatom
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
  ngrpx=1
  allocate( cg(nl**2,nl**2,(2*nl-1)**2),source=0d0)
  allocate( symope(3,3),source=reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],shape=[3,3])) !symops identicay matrix with ng=1
  call rotcg(nl-1,symope,ngrpx,cg) !CG coefficient
  call Readhamindex()
  call init_readeigen()   !Initialization for readeigen
  call init_readeigen2()
  call readngmx('QGpsi',ngpmx) !max number of the set q+G
  allocate(geig1 (ngpmx*nspc,nband),geig2(ngpmx*nspc,nband))
  allocate(cphi1 (ndima*nspc,nband),cphi2(ndima*nspc,nband) )
  open(newunit=ifoc,file='@MNLA_CPHI')
  ldim2 = ndima
  read(ifoc,*)
  allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
  do ix =1,ldim2
    read(ifoc,*) m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx !m,m,l,ibas index 
    if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
  enddo
  close(ifoc)
  if(mpi__root) then
    write(stdo,*) ' Used k number in Q0P =', nq0i
    write(stdo,"(i3,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
  endif
  readbbvec: block
    open(newunit=ifbb,file='BBVEC')
    read(ifbb,*)
    read(ifbb,*)nbb,nqbz2
    if (nqbz /= nqbz2) call rx('readbb: nqbz is wrong!')
    allocate(bbv(3,nbb),ikbidx(nbb,nqbz))    !call readbb(ifbb,nqbz,nspin,nbb, bbv, ikbidx, iko_ixs,iko_fxs,noxs)
    do i = 1,nbb
      read(ifbb,*) bbv(1:3,i)
    enddo
    do iq = 1,nqbz
      read(ifbb,*) !itmp,u(1:3)
      do ib = 1,nbb
        read(ifbb,*)itmp,itmp2,ikbidx(ib,iq) !,u(1:3)
      enddo
    enddo
    read(ifbb,*)
    read(ifbb,*)nspin2
    if(nspin /= nspin2) call rx('nspin is wrong!')
    do is = 1,nspin
      read(ifbb,*)iko_ixs(is),iko_fxs(is)
    enddo
    close(ifbb)
  end block readbbvec
  head(2,1:2)=['UUU.','UUD.']
  head(3,1:2)=['UUq0U.','UUq0D.']
  if(mpi__root) then
    do isp=1,nspin
      open(newunit=ifuu(isp),file=trim(head(ixc,isp))//charnum4(0),form='unformatted')
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
  ! --- Set q1(j1range) q2(j2range); Note that the true q when we generate eigenfunctions are q1x and q2x.
  ! q1-q1x should be a G vector.  So you may need to take into account the phase shift to <u|u> vectors.
  j1min = minval(iko_ixs(1:nspin)) !starting band index
  j1max = maxval(iko_fxs(1:nspin))
  j2min = j1min
  j2max = j1max
  allocate( uum(j1min:j1max, j2min:j2max,nspin) ) ! uumatrix allocated
  if(cmdopt0('--qibzonly')) call set_qibz(plat,qbz,nqbz,symops,ngrp) !If only at qibz, we need to set irotg
  if (ixc == 2) nbbloop = nbb
  if (ixc == 3) nbbloop = nq0i
  lxx=2*(nl-1)
  allocate(ppj(ndima,ndima,nspin,nbbloop),source=(0d0,0d0)) ! ppj: ovalap matrix within MT
  allocate(ppbrd(0:nl-1,nn,0:nl-1,nn,0:2*(nl-1),nspin,nbas), rprodx(nrx,0:lxx), phij(0:lxx),psij(0:lxx),rphiphi(nrx))
  allocate(cy((lxx+1)**2),yl((lxx+1)**2))
  ibbloop0: do ibb = 1,nbbloop
    if(ixc == 2) dq=-bbv(:,ibb)  !q1(:) = qbz(:,iqbz)        !q2(:) = qbz(:,iqbz) + bbv(:,ibb)
    if(ixc == 3) dq=-q0i(:,ibb) !q1(:) = qbz(:,iqbz)         !q2(:) = qbz(:,iqbz) + q0i(:,ibb)
    if(sum(abs(dq))<1d-8) dq=(/1d-10,0d0,0d0/)
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
        do  l1 = 0, nl-1
          do  l2 = 0, nl-1
            do  n1 = 1, nindx(l1+1,ic)
              do  n2 = 1, nindx(l2+1,ic)
                rphiphi(1)       = 0d0
                rphiphi(2:nrofi(ic)) = phitoto(2:nrofi(ic),l1,n1,ic,isp)*phitoto(2:nrofi(ic),l2,n2,ic,isp)/rr(2:,ic) ! phi = u = r \phi
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
      ii = iko_ixs(ispin)
      ie = iko_fxs(ispin)
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
10201  enddo ia1loop !              !if(cmdopt0('--fpmt')) then; geig1 = readgeigf0(q1,ispin); geig2 = readgeigf0(q2,ispin);else
    enddo ispinloop02
  enddo ibbloop0
  deallocate(ppbrd, rprodx, phij, psij, rphiphi, cy, yl)
  iqbz4uum: do 1070 iqbz = 1,nqbz  !qibzonly need to be improved to balance load in ranks.
    if(mod(iqbz-1,mpi__size)/=mpi__rank) cycle !MPI
    if (cmdopt0('--qibzonly').and. irotg(iqbz)/=1)  cycle !only irreducible q point
    write(stdo,ftox)'iq =',iqbz, 'out of',nqbz
    do isp=1,nspin
      open(newunit=ifuu(isp),file=trim(head(ixc,isp))//charnum4(iqbz),form='unformatted')
    enddo
    if (ixc == 2) nbbloop = nbb
    if (ixc == 3) nbbloop = nq0i
    ibbloop: do 1080 ibb = 1,nbbloop
      write(stdo,ftox)'  ibbloop iq =',ibb,iq
      if(ixc == 2) then
        iqb = ikbidx(ibb,iqbz)
        q1(:) = qbz(:,iqbz)
        q2(:) = qbz(:,iqbz) + bbv(:,ibb)
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
          do ispin = 1,nspin
            write(ifuu(ispin))-20
            write(ifuu(ispin))iqbz,ibb,iqtmp,ibbtmp
          enddo
          cycle
        endif
      elseif (ixc == 3) then
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
      ispinloop2: do 1050 ispin=1,nspin
        ii = iko_ixs(ispin)
        ie = iko_fxs(ispin)
        cphi1 = readcphif(q1,ispin) ! MT part of eigenfunctions 
        cphi2 = readcphif(q2,ispin) 
        geig1 = readgeigf(q1,ispin) ! IPW part of eigenfunctions
        geig2 = readgeigf(q2,ispin) 
        !Since 1d20 padding in sugw.f90, uum(i,j) can be huge number for unuvailabe eigenfunctions.
        !uum(ii:ie,ii:ie,ispin) = matmul(transpose(dconjg(cphi1(:,ii:ie))),     matmul(ppj(:,:,ispin,ibb),cphi2(:,ii:ie))) &
        !     +                   matmul(transpose(dconjg(geig1(1:ngp1,ii:ie))),matmul(ppovl,geig2(1:ngp2,ii:ie)))
        uum=0d0
        do ispc=1,nspc ! For lso=0 or 2,ispin=1,nsp. For lso=1, ispin=1 ispc=1,2 nspc=2 
           ioc=(ispc-1)*ndima
           iog=(ispc-1)*ngpmx
           uum(ii:ie,ii:ie,ispin) = uum(ii:ie,ii:ie,ispin) &
             +matmul(transpose(dconjg(cphi1(ioc+1:ioc+ndima,ii:ie))),matmul(ppj(:,:,ispin,ibb),cphi2(ioc+1:ioc+ndima,ii:ie))) &
             +matmul(transpose(dconjg(geig1(iog+1:iog+ngp1, ii:ie))),matmul(ppovl,geig2(iog+1:iog+ngp2,ii:ie)))
        enddo   
        write(ifuu(ispin)) -10 !dummy
        if(ixc==2) write(ifuu(ispin)) iqbz,ibb,ikbidx(ibb,iqbz)
        if(ixc==3) write(ifuu(ispin)) iqbz,ibb
        write(ifuu(ispin)) ((uum(j1,j2,ispin),j1=ii,ie),j2=ii,ie)
        checkwirte: block
          do j1=ii,ie; j2=j1 !; do j2=j2min,j2max !checkwrite  !if(j1==j2)
            write(stdo,"('uuuiq isp=',i5,i2,' j1j2=',2i2,' q1 q2-q1=',3f8.4,x,3f8.4,' <u|u>=',2f9.4,x,f9.3)") &
                 iqbz,ispin,j1,j2,q1,q1-q2,uum(j1,j2,ispin),abs(uum(j1,j2,ispin))
          enddo !; enddo
        endblock checkwirte
1050  enddo ispinloop2
      write(stdo,*)' ============ result --- diagonal --- ==============',nspin,j1min,j1max,j2min,j2max
      deallocate(ngvecpf1, ngvecpf2, ppovl)
      close(ifuu(ispin))
1080 enddo ibbloop
    close(ifuu(1))
    if(nspin==2) close(ifuu(2))
1070 enddo iqbz4uum
  deallocate(uum)
  if (mpi__root) write(stdo,*) ' ====== end ========================================'
  call mpi_finalize(ierr)
end subroutine h_uumatrix
