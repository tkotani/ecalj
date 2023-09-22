!> Get nlat npair for FFT of H(R+T i, R'j)
module m_gennlat  
  use m_mkqp,only: bz_nabc
  use m_lattic,only: qlat=>lat_qlat,pos=>rv_a_opos,plat=>lat_plat
  use m_lmfinit,only: nbas,  alat=>lat_alat
  integer,protected:: npairmx
  integer,allocatable,protected:: npair(:,:),nlat(:,:,:,:),nqwgt(:,:,:)
  integer,allocatable,protected:: nlatS(:,:,:,:,:),nlatE(:,:,:,:,:)
contains
  subroutine m_gennlat_init(nkk)
    use m_shortn3,only:gennlat
    implicit none
    integer:: iq,ik1,ik2,ik3,nkk(3)
    integer:: nkk1,nkk2,nkk3
    logical:: ok
    call tcn('m_gennlat_init')
    nkk1=nkk(1)
    nkk2=nkk(2)
    nkk3=nkk(3)
    allocate(npair(nbas,nbas))
    npairmx=nkk1*nkk2*nkk3*2  !initial size of npairmx
    do ! to get reasonable npairmx satisfying npairmx >npair
       npairmx = npairmx + (nkk1*nkk2*nkk3+1)*.5
       allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
       allocate(nlatS(0:nkk1-1,0:nkk2-1,0:nkk3-1,nbas,nbas))
       allocate(nlatE(0:nkk1-1,0:nkk2-1,0:nkk3-1,nbas,nbas))
       call gennlat(pos,nbas,plat,nkk1,nkk2,nkk3,npairmx,npair,nlat,nqwgt,nlatS,nlatE,ok)
       if(ok) exit
       deallocate( nlat, nqwgt, nlatS,nlatE )
    enddo
    !      print *,'m_gennlat_init: We got nlat and nqwgt'
    call tcx('m_gennlat_init')
  end subroutine m_gennlat_init
end module m_gennlat
!> this is essentially a copy of m_gennlat but for sigm. nlat for FFT of H(R+T i, R'j) for sigm
module m_gennlat_sig 
  use m_lattic,only: qlat=>lat_qlat,pos=>rv_a_opos,plat=>lat_plat
  use m_lmfinit,only: nbas, alat=>lat_alat
  integer,protected:: npairmx
  integer,allocatable,protected:: npair(:,:),nlat(:,:,:,:),nqwgt(:,:,:)
  integer,allocatable,protected:: nlatS(:,:,:,:,:),nlatE(:,:,:,:,:)
contains
  subroutine m_gennlat_init_sig(nkk)
    use m_shortn3,only:gennlat
    implicit none
    integer:: iq,ik1,ik2,ik3,nkk(3),nkk1,nkk2,nkk3
    logical:: ok
    call tcn('m_gennlat_sig')
    nkk1=nkk(1)
    nkk2=nkk(2)
    nkk3=nkk(3)
    allocate(npair(nbas,nbas))
    npairmx=nkk1*nkk2*nkk3*2  !initial size of npairmx
    do ! to get reasonable npairmx satisfying npairmx >npair
       npairmx = npairmx + (nkk1*nkk2*nkk3+1)*.5
       allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
       allocate(nlatS(0:nkk1-1,0:nkk2-1,0:nkk3-1,nbas,nbas))
       allocate(nlatE(0:nkk1-1,0:nkk2-1,0:nkk3-1,nbas,nbas))
       call gennlat(pos,nbas,plat,nkk1,nkk2,nkk3,npairmx,npair,nlat,nqwgt,nlatS,nlatE,ok)
       if(ok) exit
       deallocate( nlat, nqwgt, nlatS,nlatE )
    enddo
    call tcx('m_gennlat_sig')
  end subroutine m_gennlat_init_sig
end module m_gennlat_sig
