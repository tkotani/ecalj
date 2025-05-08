subroutine GramSchmidt2(nspc,n,nv1,nv2,nv2mx, omat1,omat2, zmel1,zmel2)!Modified GramSchmidt. MTpart+IPW. Originally for AHC branch.
  implicit none
  intent(in)::          nspc,n,nv1,nv2,nv2mx, omat1,omat2
  intent(inout)::                                          zmel1,zmel2
  integer:: it,itt,n,nv1,nv2,nv2mx,nspc
  complex(8):: dnorm2,ov(n),                  vec1(nv1*nspc),   vec2(nv2mx*nspc)
  complex(8):: omat1(nv1,nv1),omat2(nv2,nv2), zmel1(nv1*nspc,n),zmel2(nv2mx*nspc,n)
  do it = 1,n
     if(abs(zmel1(1,it))>1d6) exit !To skip padding parts.
     vec1= zmel1(:,it)
     vec2= zmel2(:,it)
     do itt = 1,it-1
        ov(itt) = sum( dconjg(zmel1(1:nv1,itt))*matmul(omat1,vec1(1:nv1))) \
        +         sum( dconjg(zmel2(1:nv2,itt))*matmul(omat2,vec2(1:nv2)))
        if(nspc==2) ov(itt) = ov(itt) \
        +         sum( dconjg(zmel1(nv1+1:2*nv1,itt))*matmul(omat1,vec1(nv1+1:2*nv1))) \
        +         sum( dconjg(zmel2(nv2mx+1:nv2mx+nv2,itt))*matmul(omat2,vec2(nv2mx+1:nv2mx+nv2)))
        vec1= vec1-zmel1(:,itt)*ov(itt) !Modified GramSchmidt. 'modified is not so effective, probably because almost orthogonal already'.
        vec2= vec2-zmel2(:,itt)*ov(itt) !
     enddo
!     vec1 = vec1 - matmul(zmel1(:,1:it-1),ov(1:it-1)) !Original GramSchmidt
!     vec2 = vec2 - matmul(zmel2(:,1:it-1),ov(1:it-1))
     dnorm2 = sum( dconjg(vec1(1:nv1))*matmul(omat1,vec1(1:nv1)))   +sum( dconjg(vec2(1:nv2))*matmul(omat2,vec2(1:nv2)))
     if(nspc==2) dnorm2=dnorm2 \
     + sum( dconjg(vec1(nv1+1:2*nv1))*matmul(omat1,vec1(nv1+1:2*nv1))) +sum( dconjg(vec2(nv2mx+1:nv2mx+nv2))*matmul(omat2,vec2(nv2mx+1:nv2mx+nv2)))
     zmel1(:,it) = vec1/dnorm2**.5
     zmel2(:,it) = vec2/dnorm2**.5
  enddo
end subroutine GramSchmidt2

! subroutine GramSchmidt2(nspc,n,nv1,nv2,nv2mx, omat1,omat2, zmel1,zmel2)!MTpart+IPW part. Originally in main_huumat in AHC branch
!   implicit none
!   intent(in)::          nspc,n,nv1,nv2,nv2mx, omat1,omat2
!   intent(inout)::                                          zmel1,zmel2
!   integer:: it,itt,n,nv1,nv2,nv2mx,nspc
!   complex(8):: dnorm2,ov(n),                  vec1(nv1*nspc),   vec2(nv2mx*nspc)
!   complex(8):: omat1(nv1,nv1),omat2(nv2,nv2), zmel1(nv1*nspc,n),zmel2(nv2mx*nspc,n)
!   do it = 1,n
!      if(abs(zmel1(1,it))>1d6) exit
!      vec1= zmel1(:,it)
!      vec2= zmel2(:,it)
!      do itt = 1,it-1
!         ov(itt) = sum( dconjg(zmel1(1:nv1,itt))*matmul(omat1,vec1(1:nv1))) \
!         +         sum( dconjg(zmel2(1:nv2,itt))*matmul(omat2,vec2(1:nv2)))
!         if(nspc==2) ov(itt) = ov(itt) \
!         +         sum( dconjg(zmel1(nv1+1:2*nv1,itt))*matmul(omat1,vec1(nv1+1:2*nv1))) \
!         +         sum( dconjg(zmel2(nv2mx+1:nv2mx+nv2,itt))*matmul(omat2,vec2(nv2mx+1:nv2mx+nv2)))
!      enddo
!      vec1 = vec1 - matmul(zmel1(:,1:it-1),ov(1:it-1))
!      vec2 = vec2 - matmul(zmel2(:,1:it-1),ov(1:it-1))
!      dnorm2 = sum( dconjg(vec1(1:nv1))*matmul(omat1,vec1(1:nv1)))   +sum( dconjg(vec2(1:nv2))*matmul(omat2,vec2(1:nv2)))
!      if(nspc==2) dnorm2=dnorm2 \
!      + sum( dconjg(vec1(nv1+1:2*nv1))*matmul(omat1,vec1(nv1+1:2*nv1))) +sum( dconjg(vec2(nv2mx+1:nv2mx+nv2))*matmul(omat2,vec2(nv2mx+1:nv2mx+nv2)))
!      zmel1(:,it) = vec1/dnorm2**.5
!      zmel2(:,it) = vec2/dnorm2**.5
!   enddo
! end subroutine GramSchmidt2
