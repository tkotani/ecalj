subroutine dspbmm(mode,i1,i2,j1,j2,nc,a,lda,ija,offs,x,ldx,b,ldb)
  !- Sparse-block-matrix dense-matrix multiply
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i          0  multiply contents into b
  !i          1  add product  a x to existing contents of b
  !i          2  add product -a x to existing contents of b
  !i   i1,i2 :calculate a(i,j)*x(j) for row subblocks i = i1..i2
  !i   j1,j2 :calculate a(i,j)*x(j) for col subblocks j = j1..j2
  !i   nc    :number of columns in x and the result matrix b
  !i   a     :sparse matrix, stored in block form by rows.
  !i          a consists of a vector of matrix subblocks:
  !i          a(*,*,i) = matrix subblock i
  !i   lda   :leading dimension of a
  !i   ija   :column index packed array pointer data to array a
  !i         ija(1,*) follows essentially the same conventions
  !i         as for scalar packed arrays (see da2spr)
  !i         except that indices now refer to matrix subblocks.
  !i         ija(1,1)= n+2, where n = max(number of rows, number of cols)
  !i         ija(1,i), i = 1,..n+1 = points to first entry in a for row i
  !i         ija(1,i), i = n+2... column index element a(i).  Thus
  !o                   for row i, k ranges from ija(i) <= k < ija(i+1) and
  !o                   sum_j a_ij x_j -> sum_k a_(ija(2,k)) x_(ija(1,k))
  !i         ija(2,*)  pointers to the matrix subblocks blocks in a:
  !i         ija(2,i), i=1..n  pointers to blocks on the diagonal of a
  !i         ija(2,i), i=n+2.. pointers to elements of a, grouped by rows
  !i   offs  :offsets to first entries in matrix subblocks
  !i          offs(i,i=1..n) offset to first row in x and b for subblock i
  !i          Thus the dimension of row i = offs(i+1) - offs(i)
  !i          If a consists of scalar subblocks, offs(i) = i-1.
  !i   x     :dense matrix, and second operand
  !i   ldx   :leading dimension of x
  !o Outputs
  !o   b     :result matrix
  !o   ldb   :leading dimension of b
  !r Remarks
  !r   This routine multiplies a sparse matrix whose elements
  !r   are matrix subblocks, by a dense matrix.
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: i1,i2,j1,j2,nc,lda,ldb,ldx,ija(2,*),offs(i2),mode
  double precision :: a(lda,lda,1),x(ldx,nc),b(ldb,nc)
  ! ... Local parameters
  integer :: ir,k,pa,ofxb,ofbb,nra,nca,ic,i,j,n
  double precision :: alp
  !     integer ofx,ofb
  logical:: l_dummy_isanrg,isanrg

  ! --- Setup ---
  ! ino isanrg is logical function,       call isanrg(mode,0,2,'dspbmm:','mode',.true.)
  l_dummy_isanrg=isanrg(mode,0,2,'dspbmm:','mode',.true.)
  alp = 1
  if (mode == 2) alp = -1
  ! --- Initialize contents of b ---
  if (mode == 0) then
     do    j = 1, nc
        do   i = offs(i1)+1, offs(i2+1)
           b(i,j) = 0
        enddo
     enddo
  endif

  ! --- For each row ir, multiply a(ir,ic) x(ic) ---
  do  10  ir = i1, i2
     !       offset to b for this subblock
     !       ofbb = offs(ir) - ofb
     ofbb = offs(ir)
     !       row dimension of this subblock
     nra = offs(ir+1) - offs(ir)
     !       pointer to diagonal subblock in a
     pa  = ija(2,ir)

     !   ... b(ir) += a(ir)*x(ir).  Skip if missing diagonal element
     if (pa /= 0 .AND. ir >= j1 .AND. ir <= j2) then
        !         offset to x for this subblock
        !         ofxb = offs(ir) - ofx
        ofxb = offs(ir)
        call dgemm('N','N',nra,nc,nra,alp,a(1,1,pa),lda, &
             x(1+ofxb,1),ldx,1d0,b(ofbb+1,1),ldb)
     endif

     !  ...  b(ir) = b(ir) + a(ir,ija(k))*x(ija(k)) for nonzero blocks in a
     do  11  k = ija(1,ir), ija(1,ir+1)-1
        !         column index to a and row index to x
        ic  = ija(1,k)
        if (ic < j1 .OR. ic > j2) goto 11
        !         offset to row x for this subblock
        !         ofxb = offs(ic) - ofx
        ofxb = offs(ic)
        !         col dimension of subblock a and row dimension of x
        nca = offs(ic+1) - offs(ic)
        !         pointer to subblock in a
        pa = ija(2,k)
        !         b(ir) += a(ir,ija(k))*x(ija(k))
        !          do j = 1, nc
        !              do i = 1, nra
        !            do kk = 1, nca
        !                b(i+ofbb,j) = b(i+ofbb,j) + a(i,kk,pa) * x(kk+ofxb,j)
        !              end do
        !            end do
        !          end do
        call dgemm('N','N',nra,nc,nca,alp,a(1,1,pa),lda, &
             x(1+ofxb,1),ldx,1d0,b(1+ofbb,1),ldb)
11   enddo
10 enddo
end subroutine dspbmm
subroutine dmspbm(mode,i1,i2,j1,j2,nr,a,lda,ija,offs,x,ldx,b,ldb)
  !- Dense-matrix sparse-block-matrix multiply
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i          0  multiply contents into b
  !i          1  add product  x a to existing contents of b
  !i          2  add product -x a to existing contents of b
  !i   i1,i2 :calculate x(i)*a(i,j) for row subblocks i = i1..i2
  !i   j1,j2 :calculate x(i)*a(i,j) for col subblocks j = j1..j2
  !i   nr    :number of rows in x and the result matrix b
  !i   a     :sparse matrix, stored in block form by rows.
  !i          a consists of a vector of matrix subblocks:
  !i          a(*,*,i) = matrix subblock i
  !i   lda   :leading dimension of a
  !i   ija   :column index packed array pointer data to array a
  !i         ija(1,*) follows essentially the same conventions
  !i         as for scalar packed arrays (see da2spr)
  !i         except that indices now refer to matrix subblocks.
  !i         ija(1,1)= n+2, where n = max(number of rows, number of cols)
  !i         ija(1,i), i = 1,..n+1 = points to first entry in a for row i
  !i         ija(1,i), i = n+2... column index element a(i).  Thus
  !o                   for row i, k ranges from ija(i) <= k < ija(i+1) and
  !o                   sum_j a_ij x_j -> sum_k a_(ija(2,k)) x_(ija(1,k))
  !i         ija(2,*)  pointers to the matrix subblocks blocks in a:
  !i         ija(2,i), i=1..n  pointers to blocks on the diagonal of a
  !i         ija(2,i), i=n+2.. pointers to elements of a, grouped by rows
  !i   offs  :offsets to first entries in matrix subblocks
  !i          offs(i,i=1..n) offset to first row in x and b for subblock i
  !i          Thus the dimension of row i = offs(i+1) - offs(i)
  !i          If a consists of scalar subblocks, offs(i) = i-1.
  !i   x     :dense matrix, and first operand
  !i   ldx   :leading dimension of x
  !o Outputs
  !o   b     :result matrix
  !o   ldb   :leading dimension of b
  !r Remarks
  !r   This routine multiplies x a, with a=sparse matrix whose elements
  !r   are matrix subblocks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: i1,i2,j1,j2,nr,lda,ldb,ldx,ija(2,*),offs(i2),mode
  double precision :: a(lda,lda,1),x(ldx,nr),b(ldb,nr)
  ! ... Local parameters
  integer :: ir,k,pa,ofxb,ofbb,nra,nca,ic,i,j,k1,k2
  double precision :: alp
  !     integer ofx,ofb
  logical:: isanrg, l_dummy_isanrg

  ! --- Setup ---
  ! ino isanrg is logical function,       call isanrg(mode,0,2,'dmspbm:','mode',.true.)
  l_dummy_isanrg=isanrg(mode,0,2,'dmspbm:','mode',.true.)
  alp = 1
  if (mode == 2) alp = -1
  !     offsets shifting origin of x and b.
  !     ofx = 0
  !     ofb = 0


  ! ... Initialize contents of b
  if (mode == 0) then
     do    j = offs(j1)+1, offs(j2+1)
        do   i = 1, nr
           b(i,j) = 0
        enddo
     enddo

  endif

  ! ... Product of diagonal elements x(ir)*a(ir). Skip if missing.
  k1 = max(i1,j1)
  k2 = min(i2,j2)
  do  8  ir = k1, k2
     !       Pointer to diagonal subblock in a
     pa  = ija(2,ir)
     !       offset to x and b for this subblock
     !       ofxb = offs(ir) - ofx
     !       ofbb = offs(ir) - ofb
     ofxb = offs(ir)
     ofbb = offs(ir)
     !       row and column dimension of this subblock
     nra = offs(ir+1) - offs(ir)
     !   ... b(ir) <- x(ir)*a(ir).
     if (pa /= 0) then
        call dgemm('N','N',nr,nra,nra,alp,x(1,ofxb+1),ldx,a(1,1,pa), &
             lda,1d0,b(1,ofbb+1),ldb)
     endif
8 enddo

  ! ... b(ija(k)) += x(ir) * a(ir,ija(k)) for nonzero blocks ija
  do  10  ir = i1, i2
     !       Pointer to diagonal subblock in a
     pa  = ija(2,ir)
     !       offset to x for this subblock
     !       ofxb = offs(ir) - ofx
     ofxb = offs(ir)
     !       col dimension of x and row dimension of a in this subblock
     nra = offs(ir+1) - offs(ir)
     do  11  k = ija(1,ir), ija(1,ir+1)-1
        !         column index to a and row index to x
        ic  = ija(1,k)
        if (ic < j1 .OR. ic > j2) goto 11
        !         offset to row x for this subblock
        !         ofbb = offs(ic) - ofb
        ofbb = offs(ic)
        !         row dimension of subblocks a and col dimension of b
        nca = offs(ic+1) - offs(ic)
        !         pointer to subblock in a
        pa = ija(2,k)
        !         b(ija(k)) += x(ir) * a(ir,ija(k))
        call dgemm('N','N',nr,nca,nra,alp,x(1,1+ofxb),ldx, &
             a(1,1,pa),lda,1d0,b(1,1+ofbb),ldb)
11   enddo
10 enddo
end subroutine dmspbm

! #if TEST
! !     Test dspbmm
! subroutine fmain
!   implicit none
!   integer :: ldap,na,nmax,mda
!   parameter(ldap=3,na=5,mda=15,nmax=na*na+1)
!   double precision :: ap(ldap,ldap,na),b(mda,2),bb(mda,2),x(mda,2)
!   double precision :: a(mda,mda),xt(2,mda),bt(2,mda),bbt(2,mda)
!   integer :: offs(na+1),ija(2,nmax)
!   integer :: ipiax(na,na),i,j,k,n,ii,ip,jj,m,ni,nj,offi,offj
!   integer :: i1,i2,j1,j2

!   data ap/-3d0,1d0,2d0, &
!        1d0,-4d0,-7d0, &
!        -1d0,2d0,9d0, &
!        1d0,-4d0,-6d0, &
!        -3d0,8d0,1d0, &
!        5d0,-2d0,-5d0, &
!        -1d0,9d0,8d0, &
!        4d0,-2d0,-5d0, &
!        -1d0,3d0,2d0, &
!        3d0,-2d0,-1d0, &
!        -6d0,7d0,6d0, &
!        2d0,-3d0,-4d0, &
!        -1d0,9d0,8d0, &
!        2d0,-5d0,-2d0, &
!        -2d0,3d0,1d0/

!   data x/-3d0,2d0, &
!        -7d0,0d0, &
!        -1d0,9d0, &
!        -6d0,1d0, &
!        -3d0,1d0, &
!        -5d0,5d0, &
!        -9d0,8d0, &
!        -5d0,2d0, &
!        2d0,3d0, &
!        3d0,-2d0, &
!        -6d0,7d0, &
!        -3d0,2d0, &
!        -0d0,9d0, &
!        -5d0,0d0, &
!        0d0,3d0/

!   data ipiax / 1,0,3,0,0, &
!        5,0,2,4,0, &
!        4,5,3,2,1, &
!        0,0,0,0,4, &
!        1,0,2,0,3/

!   data offs / 0,1,2,5,8,9/

!   ! ... Unpack a
!   call dpzero(a,mda*mda)
!   n = na
!   m = na
!   do  10  i = 1, n
!      do  12  j = 1, m
!         ip = ipiax(i,j)
!         offi = offs(i)
!         offj = offs(j)
!         ni   = offs(i+1) - offi
!         nj   = offs(j+1) - offj
!         if (ip == 0 .AND. i /= j) goto 12
!         print 332, i,j,ip,ni,nj
! 332     format(' block',2i2,' ip =',i2,' dimensions',2i2)
!         if (ip == 0) goto 12
!         do  20  ii = 1, ni
!            do  20  jj = 1, nj
!               a(offi+ii,offj+jj) = ap(ii,jj,ip)
!               !          print *, ii,jj,offi+ii,offj+jj,ip,ap(ii,jj,ip)
! 20         enddo
! 12      enddo
! 10   enddo
!      ni = offs(n+1)
!      nj = offs(m+1) &
!           all yprm('a',1,a,mda*mda,mda,ni,nj) &
!           all yprm('x',1,x,mda*mda,mda,ni,2)

!      ! ... Assemble ija
!      ija(1,1) = n+2
!      k = ija(1,1)-1
!      do  30  i = 1, n
!         do  32  j = 1, m
!            ip = ipiax(i,j)
!            if (ip /= 0 .AND. i /= j) then
!               k = k+1
!               ija(2,k) = ip
!               ija(1,k) = j
!            endif
! 32      enddo
!         ija(2,i) = ipiax(i,i)
!         ija(1,i+1) = k+1
! 30   enddo

!      call dvset(b,1,2*mda,-99d0)

!      print 331, offs
! 331  format(' offsets =',6i3)
!      i1 = 2
!      i2 = n
!      j1 = 1 + 1
!      j2 = m-1

!      ! ... Test a x
!      call dspbmm(0,i1,i2,j1,j2,2,ap,ldap,ija,offs,x,mda,b,mda)
!      !     call yprm('b',1,b,mda*2,mda,ni,2)

!      print 345, ni,nj, nj,2,i1,i2,j1,j2
! 345  format(/' Multiply a x, dimensioned (',i2,',',i2,') and (', &
!           i2,',',i2,').'/' Multiply subblocks',4i2)

!      call dgemm('N','N',offs(i2+1)-offs(i1),2,offs(j2+1)-offs(j1),1d0, &
!           a(offs(i1)+1,offs(j1)+1),mda,x(offs(j1)+1,1),mda,0d0, &
!           bb(1+offs(i1),1),mda)
!      write(*,'(t2,a,t9,a,t22,a,t42,a)') &
!           'row','reference','dspbmm result','diff'
!      do  40  i = offs(i1)+1,offs(i2+1)
!         write(*,333) i,bb(i,1),bb(i,2),b(i,1),b(i,2), &
!              b(i,1)-bb(i,1),b(i,2)-bb(i,2)
! 333     format(i3,t5,2f7.1,t20,2f7.1,t35,2f7.1)
! 40   enddo

!      ! ... Test xT a
!      call dvset(xt,1,2*mda,-99d0)
!      call dvset(bt,1,2*mda,-99d0)
!      do  100  i = 1, mda
!         xt(1,i) = x(i,1)
!         xt(2,i) = x(i,2)
! 100  enddo
!      !      call yprm('a',1,a,mda*mda,mda,ni,nj)
!      !      call yprm('xt',1,xt,2*mda,2,2,ni)

!      call dmspbm(0,i1,i2,j1,j2,2,ap,ldap,ija,offs,xt,2,bt,2)
!      !     call yprm('bt',1,bt,2*mda,2,2,ni)

!      print 355, 2,ni, ni,nj, i1,i2,j1,j2
! 355  format(/' Multiply x a, dimensioned (',i2,',',i2,') and (', &
!           i2,',',i2,').'/' Multiply subblocks',4i2)

!      call dgemm('N','N',2,offs(j2+1)-offs(j1),offs(i2+1)-offs(i1),1d0, &
!           xt(1,offs(i1)+1),2,a(offs(i1)+1,offs(j1)+1),mda,0d0, &
!           bbt(1,1+offs(j1)),2)

!      write(*,'(t2,a,t9,a,t22,a,t42,a)') &
!           'row','reference','dspbmm result','diff'
!      do  50  i = offs(j1)+1, offs(j2+1)
!         ii = i
!         write(*,333) ii,bbt(1,i),bbt(2,i),bt(1,i),bt(2,i), &
!              bt(1,i)-bbt(1,i),bt(2,i)-bbt(2,i)
! 50   enddo

!    end subroutine fmain
! #endif

