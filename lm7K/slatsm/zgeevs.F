c takao modified this based on delw_2.2.
      subroutine zgeevs( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     .LDVL, VR, LDVR, ABNRM, INFO )
C- Eigenvalues, eigenvectors of a complex nonsymmetric matrix A
C  This is a driver for LAPACK zgeevx, with the following arguments
C  allocated internally: ILO,IHI,SCALE,RCONDE,RCONDV,RWORK,WORK.
C  The remaining arguments are as in zgeevx, which see for description.
C
Cu Updates
Cu   04 Mar 02 Zero out work to avoid bug in DEC cxml library
C     implicit none
*     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            INFO, LDA, LDVL, LDVR, N
      DOUBLE PRECISION   ABNRM
*     ..
*     .. Array Arguments ..
C      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
C     $                   W( * )
      DOUBLE PRECISION    A(2, LDA, * ), VL(2, LDVL, * ),
     .VR(2, LDVR, * ), W(2, * )

      INTEGER LWORK,ILO,IHI
      DOUBLE PRECISION SCALE(N),RCONDE(N),RCONDV(N),RWORK(2*N)
      REAL(8),ALLOCATABLE:: WORK(:,:)

      call tcn('zgeevs')
      LWORK = 2*N*N+2*N
      allocate(WORK(2,LWORK))
c      call dpzero(work,2*lwork)
      work=0d0
      CALL ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     .LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE,
     .RCONDV, WORK, LWORK, RWORK, INFO )
      deallocate(WORK)
      call tcx('zgeevs')
      END

