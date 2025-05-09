module m_rdpp !Read PPBRDV2_*, radial integerals <p|p b> and rotated cg coefficients cgr.
  !note nbloch is the total number of ProductBasis (within MTs).
  use m_genallcf_v3,only: nl,nn,natom,nspin 
  use m_readqg,only: ngcmx
  use m_kind, only: kp => kindrcxq
  use m_mpi,only:ipr
  public:: Rdpp
  integer,protected,public:: mdimx,nbloch,nxx,nblochpmx, mrecl
  integer,allocatable,protected,public:: nblocha(:) ,lx(:), nx(:,:)
  real(8),allocatable,protected,public:: ppbrd (:,:,:,:,:,:,:), cgr(:,:,:,:)
  logical,protected,public:: done_rdpp=.false.
  integer, protected, public :: nprecx = kp
contains
  subroutine Rdpp( ngrp, symope) 
    implicit none
    integer :: is,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1,nwordr
    integer:: ngrp,ngpmx,nqbz,nqibz, nband, n1,n2,n3,iq0, ifppb(natom)
    real(8) ::  symope(3,3,ngrp),  pi
    if(ipr) write(6,*)" rdpp: natom=",natom
    if(done_rdpp) call rx('rdpp is already called')
    allocate( nblocha(natom) ,lx(natom), nx(0:2*(nl-1),natom))
    do ic = 1,natom
       open(newunit=ifppb(ic),file='__PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10)),form='unformatted')
       read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
    enddo
    nxx = maxval( nx )
    allocate( ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*natom), cgr(nl**2,nl**2,(2*nl-1)**2,ngrp) ) 
    if(ipr) write(6,*)' ppbrd size',nl,nn,nxx,natom,nspin
    do ic = 1,natom
       do isp= 1,nspin
          nxic = maxval( nx(0:2*(nl-1),ic) )
          read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1)) !  Radial integrals ppbrd
       enddo
       close(ifppb(ic))
    enddo
    ! Belows overide the values given by genallc.
    mdimx  = maxval(nblocha) 
    nbloch = sum(nblocha)
    nblochpmx = nbloch + ngcmx ! Maximum of MPB = PBpart +  IPWpartforMPB     !! ---------- WV.d
    mrecl  = nprecx*2*nblochpmx*nblochpmx !/nwordr()!record size
    if(ipr) write(6,*)' rdpp mdimx=',mdimx
    cgr=1d99
    call rotcg(nl-1,symope,ngrp,cgr) ! --- rotated CG setup
    done_rdpp=.true.
!    write(6,*)' rdpp:end '
  end subroutine rdpp
end module m_rdpp
subroutine rdpp_v3(nxx, nl,ngrp, nn, natom, nspin,symope, &
  nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr)
  implicit none
  integer :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp, &
       natom,nspin,nn,nblochpmx,nbloch,mdimx, &
       n1,n2,n3,iq0, &
       nblocha(natom) ,lx(natom),ifppb(natom)
  real(8)    ::  symope(3,3,ngrp) !, pi !qbas(3,3)
  integer :: is,iqi,iq,ic,isp,ip1,ip2,ioff,nxic, &
       ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
  integer:: nx(0:2*(nl-1),natom)
  real(8):: ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*natom), &
       cgr(nl**2,nl**2,(2*nl-1)**2,ngrp)
  write(6,*)" rdpp_v3: "
  !!  Radial integrals ppbrd
  do ic = 1,natom
     open(newunit=ifppb(ic),file='__PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10)),form='unformatted')
     read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
  enddo
  write(6,*)' ppbrd size',nl,nn,nxx,natom,nspin
  do ic = 1,natom
     do isp= 1,nspin
        nxic = maxval( nx(0:2*(nl-1),ic) )
        read(ifppb(ic)) ppbrd(:,:,:,:,:,1:nxic, isp+nspin*(ic-1))
     enddo
     close(ifppb(ic))
  enddo
  mdimx  = maxval(nblocha)
  nbloch = sum(nblocha)
  cgr=1d99
  call rotcg(nl-1,symope,ngrp,cgr)
  write(6,*)' rdpp_v3:end '
end subroutine rdpp_v3
subroutine Getsrdpp2(natom,nl,nxx)
  integer,intent(in):: natom,nl
  integer :: nx(0:2*(nl-1),natom),nxx,ifppb,ic,lxx,nblocha
  do ic = 1,natom
     open(newunit=ifppb, file='__PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10)), action='read',form='unformatted')
     read(ifppb) nblocha,lxx, nx(0:2*(nl-1),ic)
     close(ifppb)
  enddo
  nxx   = maxval( nx )
end subroutine Getsrdpp2
