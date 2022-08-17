subroutine mkppovl2(alat,plat,qlat,ng1,ngvec1,ng2,ngvec2,nbas,rmax,bas, ppovl)
  !- < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
  implicit none
  integer::  nbas, ng1,ng2,nx(3),ig1,ig2, ngvec1(3,ng1),ngvec2(3,ng2), &
       n1x,n2x,n3x,n1m,n2m,n3m
  real(8) :: tripl,rmax(nbas),pi = 4d0*datan(1d0)
  real(8) :: plat(3,3),qlat(3,3), &
       alat,tpibaqlat(3,3),voltot, bas(3,nbas)
  complex(8) :: ppovl(ng1,ng2)
  complex(8),allocatable :: ppox(:,:,:)
  logical:: debug=.false.
  voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  tpibaqlat =  2*pi/alat *qlat
  ! < G1 | G2 >
  n1x = maxval( ngvec2(1,:)) - minval( ngvec1(1,:))
  n1m = minval( ngvec2(1,:)) - maxval( ngvec1(1,:))
  n2x = maxval( ngvec2(2,:)) - minval( ngvec1(2,:))
  n2m = minval( ngvec2(2,:)) - maxval( ngvec1(2,:))
  n3x = maxval( ngvec2(3,:)) - minval( ngvec1(3,:))
  n3m = minval( ngvec2(3,:)) - maxval( ngvec1(3,:))
  if(debug) print *,' mkppovl2: 1 ',n1x,n1m,n2x,n2m,n3x,n3m
  allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
  ppox = 1d99
  do ig1  = 1, ng1
     do ig2  = 1, ng2
        nx(1:3) = ngvec2(1:3,ig2) - ngvec1(1:3,ig1) ! G2-G1
        if( ppox(nx(1),nx(2),nx(3)) == 1d99 ) then
           call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat, &
                nx(1:3),ppox( nx(1),nx(2),nx(3)))
        endif
     enddo
  enddo
  if(debug) print *,' mkppovl2: 2 ',n1x,n1m,n2x,n2m,n3x,n3m
  do ig1 = 1,ng1
     do ig2 = 1,ng2
        nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
        ppovl(ig1,ig2) = ppox( nx(1),nx(2),nx(3) )
     enddo
  enddo
  deallocate(ppox)
end subroutine mkppovl2

