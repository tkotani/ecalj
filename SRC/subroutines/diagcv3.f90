!     For diagonalizing Non-Hermite matrix (2017/09/22, Okumura)
!     For Hermite matrix, use 'diagcv2.F'
subroutine diagcvuh3(uh,n,ev)
  ! uh:Non-Hermite matrix, n:dimenstion, ev:eigenvalue(COMPLEX)
  implicit none
  integer(4):: nmx,nev,i,n
  complex(8):: uh(n,n),hhx(n,n),oo(n,n),zz(n,n)
  complex(8):: ev(n) !complex eigenvalue
  hhx=uh
  nmx=n
  oo = 0d0
  do i=1,n
     oo(i,i) = 1d0
  enddo
  call diagcvz(oo,hhx,zz,n, ev,nmx,1d99,nev)
end subroutine diagcvuh3
! c
subroutine diagcvz(s,uh,t,n,evl,nmx,emx,nev)
  !  diagonalizes and returns the nev lowest eigenstates
  !  eigenvecs are returned for i.lt.nmx. emx is dummy
  implicit none
  integer :: n,nmx,nev
  double precision :: emx,abstol,dlamch
  complex(8)::evl(n)
  complex*16 s(n,n),uh(n,n),t(n,n)
  complex(8),allocatable:: WORK(:)
  real(8),allocatable:: RWORK(:)
  !      complex*16:: vl,vr       !bugfix (whis was not decleared) jan2013
  complex*16:: vl(1),vr(1)
  integer:: LWORK,info

  abstol= 2d0*DLAMCH('S')

  !---find optimum work size
  allocate( WORK(1))
  lwork = -1
  call ZGEEV( 'N', 'N', n, uh, n, evl, vl, 1,  vr, 1, &
       WORK, LWORK, RWORK, INFO)
  LWORK = WORK(1) + 1
  deallocate(WORK)
  !      LWORK = max(1,2*n-1)    ! This caused a  problem---why?
  !      LWORK = max(1,2*n-1)+1  ! This caused no problem---why?

  allocate( WORK(LWORK),RWORK(7*n))
!!! not return eigenvector (vl and vr)
  call ZGEEV( 'N', 'N', n, uh, n, evl, vl, 1,  vr, 1, &
       WORK, LWORK, RWORK, INFO)
  if(INFO/=0) then
     print *, 'See http://www.netlib.org/cgi-bin/netlibget.pl', &
          '/lapack/complex16/zgeev.f'
     print *, 'ZGEEV error info=',info
  endif
  deallocate(WORK,RWORK)
  return
end subroutine diagcvz

! ccc
!     For Hermite matrix, use 'diagcv2.F'
subroutine diagcvuh3_vec(uh,n,ev,vr)
  ! uh:Non-Hermite matrix, n:dimenstion, ev:eigenvalue(COMPLEX), vr:right-eigenvecter
  implicit none
  integer(4):: nmx,nev,i,n
  complex(8):: uh(n,n),hhx(n,n),oo(n,n),zz(n,n)
  complex(8):: ev(n) !complex eigenvalue
  complex(8):: vr(n,n)
  hhx=uh
  nmx=n
  oo = 0d0
  do i=1,n
     oo(i,i) = 1d0
  enddo
  call diagcvz_r(oo,hhx,zz,n, ev,nmx,1d99,nev,vr)
end subroutine diagcvuh3_vec
! cc
subroutine diagcvz_r(s,uh,t,n,evl,nmx,emx,nev,vr)
  !  diagonalizes and returns the nev lowest eigenstates
  !  eigenvecs are returned for i.lt.nmx. emx is dummy
  implicit none
  integer :: n,nmx,nev
  double precision :: emx,abstol,dlamch
  complex(8)::evl(n)
  complex*16 s(n,n),uh(n,n),t(n,n)
  complex(8),allocatable:: WORK(:)
  real(8),allocatable:: RWORK(:)
  !      complex*16:: vl,vr       !bugfix (whis was not decleared) jan2013
  complex*16:: vl(1) !vr(1)
  complex(8):: vr(n,n)
  integer:: LWORK,info

  abstol= 2d0*DLAMCH('S')

  !---find optimum work size
  allocate( WORK(1))
  lwork = -1
  call ZGEEV( 'N', 'V', n, uh, n, evl, vl, 1,  vr, n, &
       WORK, LWORK, RWORK, INFO)
  LWORK = WORK(1) + 1
  deallocate(WORK)
  !      LWORK = max(1,2*n-1)    ! This caused a  problem---why?
  !      LWORK = max(1,2*n-1)+1  ! This caused no problem---why?

  allocate( WORK(LWORK),RWORK(7*n))
!!! return eigenvector (vl and vr)
  call ZGEEV( 'N', 'V', n, uh, n, evl, vl, 1,  vr, n, &
       WORK, LWORK, RWORK, INFO)
  if(INFO/=0) then
     print *, 'See http://www.netlib.org/cgi-bin/netlibget.pl', &
          '/lapack/complex16/zgeev.f'
     print *, 'ZGEEV error info=',info
  endif
  deallocate(WORK,RWORK)
  return
end subroutine diagcvz_r

!-----------------------------------------------------------------------------
! lapack driver routine
! http://www.netlib.org/lapack/complex16/zgeev.f

SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
     WORK, LWORK, RWORK, INFO )

  !  -- LAPACK driver routine (version 3.1) --
  !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  !     November 2006

  !     .. Scalar Arguments ..
  CHARACTER          JOBVL, JOBVR
  INTEGER ::            INFO, LDA, LDVL, LDVR, LWORK, N
  !     ..
  !     .. Array Arguments ..
  DOUBLE PRECISION ::   RWORK( * )
  COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
       W( * ), WORK( * )
  !     ..

  !  Purpose
  !  =======

  !  ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
  !  eigenvalues and, optionally, the left and/or right eigenvectors.

  !  The right eigenvector v(j) of A satisfies
  !                   A * v(j) = lambda(j) * v(j)
  !  where lambda(j) is its eigenvalue.
  !  The left eigenvector u(j) of A satisfies
  !                u(j)**H * A = lambda(j) * u(j)**H
  !  where u(j)**H denotes the conjugate transpose of u(j).

  !  The computed eigenvectors are normalized to have Euclidean norm
  !  equal to 1 and largest component real.

  !  Arguments
  !  =========

  !  JOBVL   (input) CHARACTER*1
  !          = 'N': left eigenvectors of A are not computed;
  !          = 'V': left eigenvectors of are computed.

  !  JOBVR   (input) CHARACTER*1
  !          = 'N': right eigenvectors of A are not computed;
  !          = 'V': right eigenvectors of A are computed.

  !  N       (input) INTEGER
  !          The order of the matrix A. N >= 0.

  !  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
  !          On entry, the N-by-N matrix A.
  !          On exit, A has been overwritten.

  !  LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= max(1,N).

  !  W       (output) COMPLEX*16 array, dimension (N)
  !          W contains the computed eigenvalues.

  !  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
  !          If JOBVL = 'V', the left eigenvectors u(j) are stored one
  !          after another in the columns of VL, in the same order
  !          as their eigenvalues.
  !          If JOBVL = 'N', VL is not referenced.
  !          u(j) = VL(:,j), the j-th column of VL.

  !  LDVL    (input) INTEGER
  !          The leading dimension of the array VL.  LDVL >= 1; if
  !          JOBVL = 'V', LDVL >= N.

  !  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
  !          If JOBVR = 'V', the right eigenvectors v(j) are stored one
  !          after another in the columns of VR, in the same order
  !          as their eigenvalues.
  !          If JOBVR = 'N', VR is not referenced.
  !          v(j) = VR(:,j), the j-th column of VR.

  !  LDVR    (input) INTEGER
  !          The leading dimension of the array VR.  LDVR >= 1; if
  !          JOBVR = 'V', LDVR >= N.

  !  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
  !          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.

  !  LWORK   (input) INTEGER
  !          The dimension of the array WORK.  LWORK >= max(1,2*N).
  !          For good performance, LWORK must generally be larger.

  !          If LWORK = -1, then a workspace query is assumed; the routine
  !          only calculates the optimal size of the WORK array, returns
  !          this value as the first entry of the WORK array, and no error
  !          message related to LWORK is issued by XERBLA.

  !  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)

  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value.
  !          > 0:  if INFO = i, the QR algorithm failed to compute all the
  !                eigenvalues, and no eigenvectors have been computed;
  !                elements and i+1:N of W contain eigenvalues which have
  !                converged.

  !  =====================================================================

  !     .. Parameters ..
  DOUBLE PRECISION ::   ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
  !     ..
  !     .. Local Scalars ..
  LOGICAL ::            LQUERY, SCALEA, WANTVL, WANTVR
  CHARACTER          SIDE
  INTEGER ::            HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU, &
       IWRK, K, MAXWRK, MINWRK, NOUT
  DOUBLE PRECISION ::   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
  COMPLEX*16         TMP
  !     ..
  !     .. Local Arrays ..
  LOGICAL ::            SELECT( 1 )
  DOUBLE PRECISION ::   DUM( 1 )
  !     ..
  !     .. External Subroutines ..
  EXTERNAL           DLABAD, XERBLA, ZDSCAL, ZGEBAK, ZGEBAL, ZGEHRD, &
       ZHSEQR, ZLACPY, ZLASCL, ZSCAL, ZTREVC, ZUNGHR
  !     ..
  !     .. External Functions ..
  LOGICAL ::            LSAME
  INTEGER ::            IDAMAX, ILAENV
  DOUBLE PRECISION ::   DLAMCH, DZNRM2, ZLANGE
  EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DZNRM2, ZLANGE
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          DBLE, DCMPLX, DCONJG, DIMAG, MAX, SQRT
  !     ..
  !     .. Executable Statements ..

  !     Test the input arguments


  INFO = 0
  LQUERY = ( LWORK.EQ.-1 )
  WANTVL = LSAME( JOBVL, 'V' )
  WANTVR = LSAME( JOBVR, 'V' )
  IF( ( .NOT. WANTVL ) .AND. ( .NOT. LSAME( JOBVL, 'N' ) ) ) THEN
     INFO = -1
  ELSE IF( ( .NOT. WANTVR ) .AND. ( .NOT. LSAME( JOBVR, 'N' ) ) ) THEN
     INFO = -2
  ELSE IF( N < 0 ) THEN
     INFO = -3
  ELSE IF( LDA < MAX( 1, N ) ) THEN
     INFO = -5
  ELSE IF( LDVL < 1 .OR. ( WANTVL .AND. LDVL < N ) ) THEN
     INFO = -8
  ELSE IF( LDVR < 1 .OR. ( WANTVR .AND. LDVR < N ) ) THEN
     INFO = -10
  END IF

  !     Compute workspace
  !      (Note: Comments in the code beginning "Workspace:" describe the
  !       minimal amount of workspace needed at that point in the code,
  !       as well as the preferred amount for good performance.
  !       CWorkspace refers to complex workspace, and RWorkspace to real
  !       workspace. NB refers to the optimal block size for the
  !       immediately following subroutine, as returned by ILAENV.
  !       HSWORK refers to the workspace preferred by ZHSEQR, as
  !       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
  !       the worst case.)

  IF( INFO == 0 ) THEN
     IF( N == 0 ) THEN
        MINWRK = 1
        MAXWRK = 1
     ELSE
        MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
        MINWRK = 2*N
        IF( WANTVL ) THEN
           MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', &
                ' ', N, 1, N, -1 ) )
           CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VL, LDVL, &
                WORK, -1, INFO )
        ELSE IF( WANTVR ) THEN
           MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', &
                ' ', N, 1, N, -1 ) )
           CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VR, LDVR, &
                WORK, -1, INFO )
        ELSE
           CALL ZHSEQR( 'E', 'N', N, 1, N, A, LDA, W, VR, LDVR, &
                WORK, -1, INFO )
        END IF
        HSWORK = WORK( 1 )
        MAXWRK = MAX( MAXWRK, HSWORK, MINWRK )
     END IF
     WORK( 1 ) = MAXWRK

     IF( LWORK < MINWRK .AND. .NOT. LQUERY ) THEN
        INFO = -12
     END IF
  END IF

  IF( INFO /= 0 ) THEN
     CALL XERBLA( 'ZGEEV ', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF

  !     Quick return if possible

  IF( N == 0 ) &
       RETURN

  !     Get machine constants

  EPS = DLAMCH( 'P' )
  SMLNUM = DLAMCH( 'S' )
  BIGNUM = ONE / SMLNUM
  CALL DLABAD( SMLNUM, BIGNUM )
  SMLNUM = SQRT( SMLNUM ) / EPS
  BIGNUM = ONE / SMLNUM

  !     Scale A if max element outside range [SMLNUM,BIGNUM]

  ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
  SCALEA = .FALSE.
  IF( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN
     SCALEA = .TRUE.
     CSCALE = SMLNUM
  ELSE IF( ANRM > BIGNUM ) THEN
     SCALEA = .TRUE.
     CSCALE = BIGNUM
  END IF
  IF( SCALEA ) &
       CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )

  !     Balance the matrix
  !     (CWorkspace: none)
  !     (RWorkspace: need N)

  IBAL = 1
  CALL ZGEBAL( 'B', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )

  !     Reduce to upper Hessenberg form
  !     (CWorkspace: need 2*N, prefer N+N*NB)
  !     (RWorkspace: none)

  ITAU = 1
  IWRK = ITAU + N
  CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
       LWORK-IWRK+1, IERR )

  IF( WANTVL ) THEN

     !        Want left eigenvectors
     !        Copy Householder vectors to VL

     SIDE = 'L'
     CALL ZLACPY( 'L', N, N, A, LDA, VL, LDVL )

     !        Generate unitary matrix in VL
     !        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
     !        (RWorkspace: none)

     CALL ZUNGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
          LWORK-IWRK+1, IERR )

     !        Perform QR iteration, accumulating Schur vectors in VL
     !        (CWorkspace: need 1, prefer HSWORK (see comments) )
     !        (RWorkspace: none)

     IWRK = ITAU
     CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, &
          WORK( IWRK ), LWORK-IWRK+1, INFO )

     IF( WANTVR ) THEN

        !           Want left and right eigenvectors
        !           Copy Schur vectors to VR

        SIDE = 'B'
        CALL ZLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
     END IF

  ELSE IF( WANTVR ) THEN

     !        Want right eigenvectors
     !        Copy Householder vectors to VR

     SIDE = 'R'
     CALL ZLACPY( 'L', N, N, A, LDA, VR, LDVR )

     !        Generate unitary matrix in VR
     !        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
     !        (RWorkspace: none)

     CALL ZUNGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
          LWORK-IWRK+1, IERR )

     !        Perform QR iteration, accumulating Schur vectors in VR
     !        (CWorkspace: need 1, prefer HSWORK (see comments) )
     !        (RWorkspace: none)

     IWRK = ITAU
     CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, &
          WORK( IWRK ), LWORK-IWRK+1, INFO )

  ELSE

     !        Compute eigenvalues only
     !        (CWorkspace: need 1, prefer HSWORK (see comments) )
     !        (RWorkspace: none)

     IWRK = ITAU
     CALL ZHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, &
          WORK( IWRK ), LWORK-IWRK+1, INFO )
  END IF

  !     If INFO > 0 from ZHSEQR, then quit

  IF( INFO > 0 ) &
       GO TO 50

  IF( WANTVL .OR. WANTVR ) THEN

     !        Compute left and/or right eigenvectors
     !        (CWorkspace: need 2*N)
     !        (RWorkspace: need 2*N)

     IRWORK = IBAL + N
     CALL ZTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
          N, NOUT, WORK( IWRK ), RWORK( IRWORK ), IERR )
  END IF

  IF( WANTVL ) THEN

     !        Undo balancing of left eigenvectors
     !        (CWorkspace: none)
     !        (RWorkspace: need N)

     CALL ZGEBAK( 'B', 'L', N, ILO, IHI, RWORK( IBAL ), N, VL, LDVL, &
          IERR )

     !        Normalize left eigenvectors and make largest component real

     DO 20 I = 1, N
        SCL = ONE / DZNRM2( N, VL( 1, I ), 1 )
        CALL ZDSCAL( N, SCL, VL( 1, I ), 1 )
        DO 10 K = 1, N
           RWORK( IRWORK+K-1 ) = DBLE( VL( K, I ) )**2 + &
                DIMAG( VL( K, I ) )**2
10      enddo
        K = IDAMAX( N, RWORK( IRWORK ), 1 )
        TMP = DCONJG( VL( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
        CALL ZSCAL( N, TMP, VL( 1, I ), 1 )
        VL( K, I ) = DCMPLX( DBLE( VL( K, I ) ), ZERO )
20   enddo
  END IF

  IF( WANTVR ) THEN

     !        Undo balancing of right eigenvectors
     !        (CWorkspace: none)
     !        (RWorkspace: need N)

     CALL ZGEBAK( 'B', 'R', N, ILO, IHI, RWORK( IBAL ), N, VR, LDVR, &
          IERR )

     !        Normalize right eigenvectors and make largest component real

     DO 40 I = 1, N
        SCL = ONE / DZNRM2( N, VR( 1, I ), 1 )
        CALL ZDSCAL( N, SCL, VR( 1, I ), 1 )
        DO 30 K = 1, N
           RWORK( IRWORK+K-1 ) = DBLE( VR( K, I ) )**2 + &
                DIMAG( VR( K, I ) )**2
30      enddo
        K = IDAMAX( N, RWORK( IRWORK ), 1 )
        TMP = DCONJG( VR( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
        CALL ZSCAL( N, TMP, VR( 1, I ), 1 )
        VR( K, I ) = DCMPLX( DBLE( VR( K, I ) ), ZERO )
40   enddo
  END IF

  !     Undo scaling if necessary

50 CONTINUE
  IF( SCALEA ) THEN
     CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), &
          MAX( N-INFO, 1 ), IERR )
     IF( INFO > 0 ) THEN
        CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR )
     END IF
  END IF

  WORK( 1 ) = MAXWRK
  RETURN

  !     End of ZGEEV

END SUBROUTINE ZGEEV