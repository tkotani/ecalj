module m_rdpp !Read PPBRDV2_*, radial integerals <p|p b> and rotated cg coefficients cgr.
  use m_genallcf_v3,only: nl,nn,nclass,nspin 
  use m_readqg,only: ngcmx
  public:: Rdpp
  integer,protected,public:: mdimx,nbloch,nxx,nblochpmx, nprecx, mrecl
  integer,allocatable,protected,public:: nblocha(:) ,lx(:), nx(:,:)
  real(8),allocatable,protected,public:: ppbrd (:,:,:,:,:,:,:), cgr(:,:,:,:)
  logical,protected,public:: done_rdpp=.false.
  private
  integer:: ndble=8
contains
  subroutine Rdpp( ngrp, symope) 
    implicit none
    integer :: is,iqi,iq,ic,isp,ip1,ip2,ioff,nxic,ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1,nwordr
    integer:: ngrp,ngpmx,nqbz,nqibz, nband, n1,n2,n3,iq0, ifppb(nclass)
    real(8) ::  symope(3,3,ngrp),  pi
    character(11) :: filename(nclass)
    write(6,*)" rdpp: nclass=",nclass
    if(done_rdpp) call rx('rdpp is already called')
    allocate( nblocha(nclass) ,lx(nclass), nx(0:2*(nl-1),nclass))
    do ic = 1,nclass
       filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
       open(newunit=ifppb(ic),file=trim(filename(ic)),form='unformatted')
       read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
    enddo
    nxx = maxval( nx )
    allocate( ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass), cgr(nl**2,nl**2,(2*nl-1)**2,ngrp) ) 
    write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
    do ic = 1,nclass
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
    nprecx = ndble          ! We use double precision arrays only.
    mrecl  = nprecx*2*nblochpmx*nblochpmx !/nwordr()!record size
    write(6,*)' rdpp mdimx=',mdimx
    cgr=1d99
    call rotcg(nl-1,symope,ngrp,cgr) ! --- rotated CG setup
    done_rdpp=.true.
    write(6,*)' rdpp:end '
  end subroutine rdpp
end module m_rdpp

subroutine rdpp_v3(nxx, nl,ngrp, nn, nclass, nspin,symope, &
  nblocha, lx, nx,  ppbrd , mdimx,nbloch, cgr)
  implicit none
  integer(4) :: ngpmx,ngcmx,nxx,  nqbz,nqibz, nband,nl,ngrp, &
       nclass,nspin,nn,nblochpmx,nbloch,mdimx, &
       n1,n2,n3,iq0, &
       nblocha(nclass) ,lx(nclass),ifppb(nclass)
  real(8)    ::  symope(3,3,ngrp) !, pi !qbas(3,3)
  integer(4) :: is,iqi,iq,ic,isp,ip1,ip2,ioff,nxic, &
       ifplane ,ngpmx_dum, ngcmx_dum,iqbzx,idxk,ngp,ngc,ig1
  character(11) :: filename(nclass)
  integer:: nx(0:2*(nl-1),nclass)
  real(8):: ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass), &
       cgr(nl**2,nl**2,(2*nl-1)**2,ngrp)
  write(6,*)" rdpp_v3: "
  !!  Radial integrals ppbrd
  do ic = 1,nclass
     filename(ic)='PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
     open(newunit=ifppb(ic),file=trim(filename(ic)),form='unformatted')
     read(ifppb(ic)) nblocha(ic),lx(ic),nx(0:2*(nl-1),ic)
  enddo
  write(6,*)' ppbrd size',nl,nn,nxx,nclass,nspin
  do ic = 1,nclass
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
subroutine Getsrdpp2(nclass,nl,nxx)
  integer(4),intent(in):: nclass,nl
  integer(4) :: nx(0:2*(nl-1),nclass),nxx,ifppb,ic,lxx,nblocha
  character(20) :: filename
  do ic = 1,nclass
     filename = 'PPBRD_V2_'//char( 48+ic/10 )//char( 48+mod(ic,10))
     open(newunit=ifppb, file=trim(filename), action='read',form='unformatted')
     read(ifppb) nblocha,lxx, nx(0:2*(nl-1),ic)
     close(ifppb)
  enddo
  nxx   = maxval( nx )
end subroutine Getsrdpp2
