module m_gennlat  ! Get nlat npair for FFT of H(R+T i, R'j)
  use m_mkqp,only: bz_nabc
  use m_lattic,only: qlat=>lat_qlat,pos=>rv_a_opos,plat=>lat_plat
  use m_lmfinit,only: nbas,  alat=>lat_alat
  integer,protected:: npairmx
  integer,allocatable,protected:: npair(:,:),nlat(:,:,:,:),nqwgt(:,:,:)
  integer,allocatable,protected:: nlatS(:,:,:,:,:),nlatE(:,:,:,:,:)
contains
  subroutine m_gennlat_init(nkk)
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
! ssssssssssssssssssssssssssssssssss
module m_gennlat_sig !this is essentially a copy of m_gennlat
  ! nlat for FFT of H(R+T i, R'j) sigm
  use m_lattic,only: qlat=>lat_qlat,pos=>rv_a_opos,plat=>lat_plat
  use m_lmfinit,only: nbas, alat=>lat_alat
  integer,protected:: npairmx
  integer,allocatable,protected:: npair(:,:),nlat(:,:,:,:),nqwgt(:,:,:)
  integer,allocatable,protected:: nlatS(:,:,:,:,:),nlatE(:,:,:,:,:)
contains
  subroutine m_gennlat_init_sig(nkk)
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
! ssssssssssssssssssssssssssssssssss
subroutine gennlat(pos,nbas,plat,nk1,nk2,nk3,npairmx,npair,nlat,nqwgt,nlatS,nlatE,ok)
  use m_shortn3,only: shortn3_initialize,shortn3,nout,nlatout
  implicit none
  intent(in)::       pos,nbas,plat,nk1,nk2,nk3,npairmx
  intent(out)::                                        ok!,npair,  t,nqwgt
  integer:: nbas,nk1,nk2,nk3,npairmx ,npair(nbas,nbas)
  integer:: ib1,ib2,nnn,ik,nmax(3),ix,i,j,ik1,ik2,ik3,ni
  integer:: nd(3),nlat(3,npairmx,nbas,nbas)
  integer:: debug=0,nqwgt(npairmx,nbas,nbas)
  real(8):: pos(3,nbas),plat(3,3), pi,qlat(3,3),q(3),pin(3),rmax2 !,nrmax(3)
  real(8):: eps=1d-8
  real(8):: posp(3),rxlat(3,3),rxprod(3,3),dummy,rrr,rmax
  logical :: ok
  real(8):: rlatp(3,3),xmx2(3)
  integer:: nlatS(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas)
  integer:: nlatE(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas)
  call tcn('gennlat')
  !      print *,'gennlat:'
  ok=.true.
  pi=4d0*datan(1d0)
  call dinv33(plat,1,qlat,dummy)
  ! periodicity for R mesh for nk1 nk2 nk3 division.
  rxlat(:,1)=plat(:,1)*nk1
  rxlat(:,2)=plat(:,2)*nk2
  rxlat(:,3)=plat(:,3)*nk3
  call shortn3_initialize(rxlat)
  !      do i=1,3
  !      do j=1,3
  !        rxprod(i,j) = sum(rxlat(:,i)*rxlat(:,j))
  !      enddo
  !      enddo
  !      print *, 'goto ellips',rxprod
  !      call ellipsoidxmax(rxprod,xmx2)
  !      print *, 'end of ellips xmx2=',xmx2
  ! et Maxmum x_i**2 for ellipsoid for 1d0 = sum_{i,j} x_i rxprod(i,j) x_j
  ! Maximum value for x_i for ellipsoid
  ! Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.

  !      allocate(nlatS(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas))
  !      allocate(nlatE(0:nk1-1,0:nk2-1,0:nk3-1,nbas,nbas))
  npair=0
  do ib1=1,nbas
     do ib2=1,nbas
        if(debug>0) print *
        if(debug>0) print *,' ---- ib1 ib2=',ib1,ib2
        posp = matmul( pos(:,ib1)-pos(:,ib2), qlat(:,:)) !posp is on plat-coodinate
        nnn=0
        do ik1=0,nk1-1
           do ik2=0,nk2-1
              do ik3=0,nk3-1
                 ! rint *,' ik1,ik2,ik3=',ik1,ik2,ik3
                 nd=[ik1,ik2,ik3]
                 if( ik1 >nk1/2) nd= nd-[nk1, 0,  0]
                 if( ik2 >nk2/2) nd= nd-[0, nk2,  0]
                 if( ik3 >nk3/2) nd= nd-[0, 0,  nk3]
                 ! above procedures are not necessary, but these may accelarate shortn3 a little.
                 pin(1)= (posp(1)+nd(1))/nk1
                 pin(2)= (posp(2)+nd(2))/nk2
                 pin(3)= (posp(3)+nd(3))/nk3
                 call shortn3(pin) ! return nout,nlatout
                 if(nnn+nout>npairmx) then
                    ok=.false.
                    return
                 endif
                 nlatS(ik1,ik2,ik3,ib1,ib2)= nnn+1
                 nlatE(ik1,ik2,ik3,ib1,ib2)= nnn+nout
                 do ix=1,nout
                    nlat(:,nnn+ix,ib1,ib2)= &
                         nd + [nk1*nlatout(1,ix), nk2*nlatout(2,ix), nk3*nlatout(3,ix)]
                    nqwgt(nnn+ix,ib1,ib2)= nout
                 enddo
                 do ix=1,nout
                    nqwgt(nnn+ix,ib1,ib2)= nout
                 enddo
                 nnn = nnn+nout
              enddo
           enddo
        enddo
        npair(ib1,ib2)= nnn
        !$$$            write(6,"(a,2i4,2x,i6)") ' ib1 ib2 npair=',ib1,ib2,npair(ib1,ib2)
        !$$$            do ni = 1,npair(ib1,ib2)
        !$$$              posp =  pos(:,ib1)-pos(:,ib2) + matmul(plat,nlat(:,ni,ib1,ib2))
        !$$$              rrr = sqrt(sum(posp**2))
        !$$$              write(6,"(i6,3x,3i3,f8.3,i5)") ni,nlat(1:3,ni,ib1,ib2),rrr,nqwgt(ni,ib1,ib2)
        !$$$            enddo
     enddo
  enddo
  do ib1=1,nbas
     do ib2=1,nbas
        if(abs(sum(1d0/nqwgt(1:npair(ib1,ib2),ib1,ib2))-nk1*nk2*nk3)>1d-8) call rx('bug:nqwgt sum is not unity')
     enddo
  enddo
  call tcx('gennlat')
end subroutine gennlat


