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
subroutine mkppovl(alat,plat,qlat,q,ng1,ngvec1,ng2,ngvec2,rmax,nbas,bas,lmxa,lmxax, ppovl)!<P^q1_G1 | P^q2_G2 > matrix where P^q1_G1 denotes IPW, zero within sphere.
  implicit none
  integer ::  nbas,  ng1,ng2,lmxax,ig1,ig2,ibas,la,lma,lmh,ngvec1(3,ng1),ngvec2(3,ng2),lmxa(nbas),ll
  real(8) :: absqg1,absqg2,tripl,rmax(nbas),pi4,r2s
  real(8) :: q(3),plat(3,3),qlat(3,3),qgg1(3,ng1),qgg2(3,ng2), facl &
       ,pi,alat,tpiba,cost,voltot,plegn,qg1(3),qg2(3),qqq(3), &
       bas(3,nbas), fkk(0:lmxax),fkj(0:lmxax),fjk(0:lmxax),fjj(0:lmxax) &
       ,absqg1x(ng1),absqg2x(ng2)
  real(8),allocatable::cy(:),yl(:),yl1(:),yl2(:)
  complex(8) :: img =(0d0,1d0),phase
  complex(8) :: ppovl(ng1,ng2)
  integer:: verbose
  if(verbose()>50) print *,' mkppovl:'
  pi = 4d0*datan(1d0)
  pi4=16d0*datan(1d0)
  tpiba=2*pi/alat
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  allocate(cy((lmxax+1)**2),yl((lmxax+1)**2) )
  call sylmnc(cy,lmxax)
  ! <P^q_G1 | P^q_G2 > matrix
  do ig1 = 1,ng1
     qgg1(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvec1(1:3,ig1)))
     absqg1x(ig1)  = sqrt(sum(qgg1(1:3,ig1)**2))
  enddo
  do ig2 = 1,ng2
     qgg2(1:3,ig2) = tpiba * (q(1:3)+ matmul(qlat, ngvec2(1:3,ig2)))
     absqg2x(ig2)  = sqrt(sum(qgg2(1:3,ig2)**2))
  enddo
  ppovl = 0d0
  do ig1 =1,ng1
     qg1(1:3) = qgg1(1:3,ig1)
     absqg1   = absqg1x(ig1)
     do ig2 =1,ng2
        qg2(1:3) = qgg2(1:3,ig2)
        absqg2   = absqg2x(ig2)
        if(sum(abs(ngvec1(:,ig1)-ngvec2(:,ig2)))==0) ppovl(ig1,ig2)= voltot
        if( absqg1*absqg2 == 0d0) then
           cost = 1d0
        else
           cost = sum(qg1*qg2)/absqg1/absqg2
        endif
        do ibas = 1,nbas
           call wronkj( absqg1**2, absqg2**2, rmax(ibas), lmxa(ibas),fkk,fkj,fjk,fjj)
           do la  = 0,lmxa(ibas)
              ppovl(ig1,ig2) = ppovl(ig1,ig2) - exp( img* sum((qg2-qg1)*bas(1:3,ibas))*alat )   & 
                   * pi4 *(2*la+1d0) * plegn(la,cost) * (-fjj(la))  *facl(absqg1*absqg2,la)
           enddo
        enddo
     enddo
  enddo
  if (allocated(cy)) deallocate(cy)
  if (allocated(yl)) deallocate(yl)
end subroutine mkppovl
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
  use m_lldata,only: ll
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
