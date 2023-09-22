subroutine mkppmt(alat,plat,qlat, q, ng1,ngvec1, rmax,  nbas,  bas, lmxa, lmxax,  ppmt)
  !- ppmt value and slope at MT bounadaries for q+G plane waves.
  !o  ppmt(1,lm,nbas): value for (lm,ibas)
  !o  ppmt(2,lm,nbas):
  use m_ll,only: ll
  implicit none
  integer(4) ::  nbas,  ng1,lmxax,ig1,ibas,la,lma,lmh, ngvec1(3,ng1),lmxa(nbas)
  real(8) :: absqg, rmax(nbas),pi4,r2s, q(3),plat(3,3),qlat(3,3),qgg(3),pi,alat,tpiba, bas(3,nbas) ,facl
  complex(8) :: ppmt(2,(lmxax+1)**2,nbas,ng1)
  integer(4) :: verbose,it,ip
  real(8),allocatable:: cy(:),yl(:),ylr(:)
  real(8) :: ak1(200),aj1(200), dk1(200),dj1(200),tpi,absqg2,qqq(3)
  complex(8) :: img =(0d0,1d0),fac,phase
  real(8):: theta,phi,rdir(3)
  complex(8)::valx1,valx2
  logical ::debug=.false.
  if(debug) allocate(ylr((lmxax+1)**2) )
  if(verbose()>50) print *,' mkppmt:'
  pi  = 4d0*datan(1d0)
  tpi = 2d0*pi
  pi4 = 4d0*pi
  tpiba = tpi/alat
  allocate(cy((lmxax+1)**2),yl((lmxax+1)**2) )
  call sylmnc(cy,lmxax)
  ppmt = 0d0
  do ig1 =1,ng1
     qgg(1:3)=  tpiba*( q(1:3)+ matmul(qlat, ngvec1(1:3,ig1)) )
     !   exp(i qgg *\bfr)

     absqg2  = sum(qgg(1:3)**2)
     absqg = sqrt(absqg2)
     if(absqg==0d0) then
        qqq = (/0d0,0d0,1d0/)
     else
        qqq = qgg/absqg
     endif
     call sylm(qqq,yl,lmxax,r2s) !spherical factor Y( q+G )

     do ibas = 1,nbas
        phase = exp( img*sum( qgg*bas(:,ibas) )*alat  )
        call radkj(absqg2, rmax(ibas),lmxa(ibas),ak1,aj1,dk1,dj1,0)
        do lma = 1,(lmxa(ibas)+1)**2
           la = ll(lma)
           fac = pi4* img**la * cy(lma)*yl(lma) * facl(absqg,la)
           ppmt(1,lma,ibas,ig1) = fac*phase* aj1(la+1)
           ppmt(2,lma,ibas,ig1) = fac*phase* dj1(la+1)
        enddo
     enddo
     if(debug) then!--- debug check ! val at mt in two kinds of calculations.
        do it= 0,11
           do ip= 0,13
              do ibas = 1,nbas
                 phase = exp( img*sum( qgg*bas(:,ibas) )*alat  )
                 theta = pi/2d0* it/19d0
                 phi   = 2*pi  * ip/17d0
                 rdir  = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)
                 call sylm(rdir,ylr,lmxax,r2s) !spherical factor Y( q+G )
                 rdir  = rdir*rmax(ibas)
                 valx1 = exp(img* sum( qgg*( bas(:,ibas)*alat+rdir) ) )
                 valx2 = 0d0
                 do lma = 1,(lmxa(ibas)+1)**2
                    valx2= valx2 + ppmt(1,lma,ibas,ig1)*cy(lma) *ylr(lma)
                 enddo
                 if(abs(valx2 - valx1)/abs(valx1)>0) & 
                      write(6,"(' ig itip ibas valx1 err=',4i3,2d11.3,3x,d11.2)") &
                      ig1,it,ip,ibas, valx1, abs(valx2 - valx1)/abs(valx1)
              enddo
           enddo
        enddo
        write(6,*)
     endif
  enddo
  if(debug) call rx( 'test end ----------')
  deallocate(yl,cy)
end subroutine mkppmt
real(8) function facl(a,l)
  integer :: l
  real(8) a
  if(l==0) then
     facl=1d0
  else
     facl=a**l
  endif
END function facl
subroutine mkppovl2test(alat,plat,qlat,ng1,ngvec1,ng2,ngvec2,nbas,rmax,bas, ppovl) ! < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
  implicit none
  integer ::  nbas, ng1,ng2,nx(3),ig1,ig2,ibas, ngvec1(3,ng1),ngvec2(3,ng2), n1x,n2x,n3x,n1m,n2m,n3m
  real(8) :: absqg1,absqg2,tripl,rmax(nbas),pi
  real(8) :: plat(3,3),qlat(3,3), alat,tpiba,tpibaqlat(3,3),voltot, bas(3,nbas)
  complex(8) :: img =(0d0,1d0)
  complex(8) :: ppovl(ng1,ng2)
  complex(8),allocatable :: ppox(:,:,:)
  print *,' mkppovl2test:'
  ppovl = 0d0
  do ig1 = 1,ng1
     do ig2 = 1,ng2
        nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
        if(sum(nx**2)==0)  ppovl(ig1,ig2) = 1d0
     enddo
  enddo
end subroutine mkppovl2test
