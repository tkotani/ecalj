module m_GramSchmidt
  use m_blas,only: zmm => zmm_h, m_op_C, zmv => zmv_h, zvv => zvv_h
  public GramSchmidt2,CB_GramSchmidt
contains
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
        ov(itt) = sum( dconjg(zmel1(1:nv1,itt))*matmul(omat1,vec1(1:nv1))) &
             +         sum( dconjg(zmel2(1:nv2,itt))*matmul(omat2,vec2(1:nv2)))
        if(nspc==2) ov(itt) = ov(itt) &
             +         sum( dconjg(zmel1(nv1+1:2*nv1,itt))*matmul(omat1,vec1(nv1+1:2*nv1))) &
             +         sum( dconjg(zmel2(nv2mx+1:nv2mx+nv2,itt))*matmul(omat2,vec2(nv2mx+1:nv2mx+nv2)))
        vec1= vec1-zmel1(:,itt)*ov(itt) !Modified GramSchmidt. 'modified is not so effective, probably because almost orthogonal already'.
        vec2= vec2-zmel2(:,itt)*ov(itt) !
      enddo
      !     vec1 = vec1 - matmul(zmel1(:,1:it-1),ov(1:it-1)) !Original GramSchmidt
      !     vec2 = vec2 - matmul(zmel2(:,1:it-1),ov(1:it-1))
      dnorm2 = sum( dconjg(vec1(1:nv1))*matmul(omat1,vec1(1:nv1)))   +sum( dconjg(vec2(1:nv2))*matmul(omat2,vec2(1:nv2)))
      if(nspc==2) dnorm2=dnorm2 &
           + sum( dconjg(vec1(nv1+1:2*nv1))*matmul(omat1,vec1(nv1+1:2*nv1))) +sum( dconjg(vec2(nv2mx+1:nv2mx+nv2))*matmul(omat2,vec2(nv2mx+1:nv2mx+nv2)))
      zmel1(:,it) = vec1/dnorm2**.5
      zmel2(:,it) = vec2/dnorm2**.5
    enddo
  end subroutine GramSchmidt2

  subroutine CB_GramSchmidt(nspc,n,nv1,nv2,nv2mx, omat1,omat2, zmel1,zmel2)!Original GramSchmidt NOT Modifed version. Obata 2025-6-22
    ! Column-wise Blocking Classical Gram-Schmidt (CBCGS)
    ! https://www.ccs.tsukuba.ac.jp/wp-content/uploads/sites/14/2016/12/comp_takahashi.pdf
    use m_lgunit,only:stdo
    use m_ftox
    implicit none
    intent(in)::            nspc,n,nv1,nv2,nv2mx, omat1,omat2
    intent(inout)::                                            zmel1,zmel2
    integer:: it,n,nv1,nv2,nv2mx,nspc
    complex(8):: dnorm2, vec1(nv1*nspc), vec2(nv2mx*nspc)
    complex(8):: omat1(nv1,nv1), omat2(nv2,nv2), zmel1(nv1*nspc,n),zmel2(nv2mx*nspc,n)
    integer :: istat, cblock, jt, h, nrows, mcb = 256 !blocking size
    logical :: debug = .false.
    complex(8), allocatable :: w(:), oz1(:,:), oz2(:,:), wmat(:,:), q1(:,:), q2(:,:)
    if(debug) write(stdo,ftox) 'CB_GramSchmidt: Original GramSchmidt, NOT Modified version, with Column-wise Blocking'
    allocate(q1, source = zmel1)
    allocate(q2, source = zmel2)

    cblock_loop: do jt = 1, n, mcb
      h = min(mcb, n-jt+1)
      allocate(w(h))
      if(debug) write(06,*) 'CB_GramSchmidt: = ', n, mcb, jt, h
      do it = jt, jt + h - 1
        if(abs(zmel1(1,it))>1d6) exit cblock_loop !To skip padding parts.
        if(it /= jt) then
          vec1(:) = (0d0,0d0)
          vec2(:) = (0d0,0d0)
          istat = zmv(omat1, zmel1(:,it), vec1, m=nv1, n=nv1)
          istat = zmv(omat2, zmel2(:,it), vec2, m=nv2, n=nv2)
          if(nspc==2) then
            istat = zmv(omat1, zmel1(nv1+1  ,it), vec1(nv1+1),   m=nv1, n=nv1)
            istat = zmv(omat2, zmel2(nv2mx+1,it), vec2(nv2mx+1), m=nv2, n=nv2)
          endif
          istat = zmv(q1(1,jt), vec1, w, m=nv1*nspc  , n=it-jt, opA=m_op_C)
          istat = zmv(q2(1,jt), vec2, w, m=nv2mx*nspc, n=it-jt, opA=m_op_C, beta=(1d0,0d0))
          istat = zmv(q1(1,jt), w, vec1, m=nv1*nspc  , n=it-jt)
          istat = zmv(q2(1,jt), w, vec2, m=nv2mx*nspc, n=it-jt)
          q1(:,it) = q1(:,it) - vec1(:)
          q2(:,it) = q2(:,it) - vec2(:)
        endif
        dnorm2 = get_vxv(q1(1,it), omat1, nv1) + get_vxv(q2(1,it), omat2, nv2)
        if(nspc==2) dnorm2 = dnorm2 + get_vxv(q1(nv1+1,it), omat1, nv1) + get_vxv(q2(nv2mx+1,it), omat2, nv2)
        q1(:,it) = q1(:,it)/dnorm2**.5
        q2(:,it) = q2(:,it)/dnorm2**.5
        if(it == n) exit cblock_loop
      enddo
      deallocate(w)
      nrows = n - (jt + h) + 1
      allocate(oz1(nv1*nspc  ,nrows), source=(0d0,0d0))
      allocate(oz2(nv2mx*nspc,nrows), source=(0d0,0d0))
      allocate(wmat(h,nrows))

      istat = zmm(omat1, zmel1(1,jt+h), oz1, m=nv1, n=nrows, k=nv1, ldb=nv1*nspc  , ldc=nv1*nspc  )
      istat = zmm(omat2, zmel2(1,jt+h), oz2, m=nv2, n=nrows, k=nv2, ldb=nv2mx*nspc, ldc=nv2mx*nspc)
      if(nspc==2) then
        istat = zmm(omat1, zmel1(nv1+1,  jt+h), oz1(nv1+1,  1), m=nv1, n=nrows, k=nv1, ldb=nv1*nspc  , ldc=nv1*nspc  )
        istat = zmm(omat2, zmel2(nv2mx+1,jt+h), oz2(nv2mx+1,1), m=nv2, n=nrows, k=nv2, ldb=nv2mx*nspc, ldc=nv2mx*nspc)
      endif
      istat = zmm(q1(1,jt), oz1, wmat, m=h, n=nrows, k=nv1*nspc  , opA=m_op_C)
      istat = zmm(q2(1,jt), oz2, wmat, m=h, n=nrows, k=nv2mx*nspc, opA=m_op_C, beta=(1d0,0d0))
      istat = zmm(q1(1,jt), wmat, oz1, m=nv1*nspc  , n=nrows, k=h)
      istat = zmm(q2(1,jt), wmat, oz2, m=nv2mx*nspc, n=nrows, k=h)
      q1(:,jt+h:n) = q1(:,jt+h:n) - oz1(:,1:nrows)
      q2(:,jt+h:n) = q2(:,jt+h:n) - oz2(:,1:nrows)
      deallocate(oz1, oz2, wmat)
    enddo cblock_loop
    zmel1(:,:) = q1(:,:)
    zmel2(:,:) = q2(:,:)
    deallocate(q1, q2)
  contains
  end subroutine CB_GramSchmidt
  complex(8) function get_vxv(vvec, xmat, nsize) result(vxv)
    integer, intent(in) :: nsize
    complex(8), intent(in) :: vvec(nsize), xmat(nsize,nsize)
    complex(8) :: xv(nsize)
    integer:: istat
    istat = zmv(xmat, vvec, xv, nsize, nsize) 
    istat = zvv(vvec, xv, nsize, vxv) 
  end function get_vxv
end module m_gramschmidt
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
