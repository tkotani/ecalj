#define F90 1
      subroutine zgeevs( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     .                   LDVL, VR, LDVR, ABNRM, INFO )
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
     .                    VR(2, LDVR, * ), W(2, * )


      INTEGER LWORK,ILO,IHI
#if AUTO_ARRAY | F90
      DOUBLE PRECISION SCALE(N),RCONDE(N),RCONDV(N),RWORK(2*N)
#if F90
      REAL(8),ALLOCATABLE:: WORK(:,:)
#elif AUTO_ARRAY
      DOUBLE PRECISION WORK(2,2*N*N+2*N)
#endif

      call tcn('zgeevs')
      LWORK = 2*N*N+2*N
#if F90
      allocate(WORK(2,LWORK))
#endif
      call dpzero(work,2*lwork)
      CALL ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     .             LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE,
     .             RCONDV, WORK, LWORK, RWORK, INFO )
#if F90
      deallocate(WORK)
#endif

      call tcx('zgeevs')
#else
      integer oscale,orcnde,orcndv,orwork,owork
C ... Heap
      integer ww(1)
      common /w/ ww

      call tcn('zgeevs')
      call defrr(oscale,n)
      call defrr(orcnde,n)
      call defrr(orcndv,n)
      call defrr(orwork,2*n)
      lwork = 2*n*n+2*n
C     Initialize to zero to avoid bug in dec cxml library
      call defcc(owork,-lwork)

      CALL ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
     .  LDVL, VR, LDVR, ILO, IHI, ww(oSCALE), ABNRM, ww(oRCNDE),
     .  ww(oRCNDV), ww(oWORK), LWORK, ww(oRWORK), INFO )

      call rlse(oscale)
      call tcx('zgeevs')
#endif

      END

