!>shortest p (real space) vector
module m_shortn3_plat 
! usage 2022-jul  
!  pp=matmul(transpose(qlat),p)
!  call shortn3_plat(pp)
!  p1 = matmul(plat,pp+nlatout(:,1)) !this is one of shortest vector
  implicit none
  public shortn3_plat,qlatx,platx
  integer,private,parameter:: noutmx=48
  integer,public:: nout,nlatout(3,noutmx)
  real(8):: platx(3,3),qlatx(3,3)
  
  real(8),private:: rlatp(3,3),xmx2(3)
  logical,private:: init=.true.
contains
  subroutine shortn3_plat(pin)!return shortest plat as p = qin + matmul(qlat(:,:),nlatout(:,i)) for i=1,nout
    use m_shortvec,only: shortvec, shortvecinitialize
    use m_lattic,only:   lattic_init, lat_plat,lat_qlat,lattic_init !lmf mode
    use m_hamindex,only: plat,qlat,readhamindex_init !gw mode
    real(8):: pin(3)
    if(init) then
       if(lattic_init) then !lmf case
          platx=lat_plat
          qlatx=lat_qlat
       elseif(readhamindex_init) then !gw case
          platx=plat
          qlatx=qlat
       else
          call rx('m_shortn3_plat: can not read plat')
       endif
       call shortvecinitialize(platx,rlatp,xmx2)
       init=.false.
    endif   
    call shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
  end subroutine shortn3_plat
end module m_shortn3_plat

!>shortest q 
module m_shortn3_qlat 
! usage 2022-jul  
!  pp=matmul(transpose(qlat),p)
!  call shortn3_plat(pp)
!  p1 = matmul(plat,pp+nlatout(:,1)) !this is one of shortest vector
  implicit none
  public shortn3_qlat
  integer,private,parameter:: noutmx=48
  integer,public:: nout,nlatout(3,noutmx)
  real(8):: platx(3,3),qlatx(3,3)
  
  real(8),private:: rlatp(3,3),xmx2(3)
  logical,private:: init=.true.
contains
  subroutine shortn3_qlat(qin) !return shortest qlat as q = qin + matmul(qlat(:,:),nlatout(:,i)) for i=1,nout
    use m_shortvec,only: shortvec, shortvecinitialize
    use m_lattic,only:   lattic_init, lat_plat,lat_qlat,lattic_init !lmf mode
    use m_hamindex,only: plat,qlat,readhamindex_init !gw mode
    real(8):: qin(3)
    if(init) then
       init=.false.
       if(lattic_init) then !lmf case
          platx=lat_plat
          qlatx=lat_qlat
       elseif(readhamindex_init) then !gw case
          platx=plat
          qlatx=qlat
       else
          call rx('m_shortn3_plat: can not read plat')
       endif
       call shortvecinitialize(qlatx, rlatp,xmx2)
    endif   
    call shortvec(qin,rlatp,xmx2,noutmx,  nout,nlatout)
  end subroutine shortn3_qlat
end module m_shortn3_qlat

