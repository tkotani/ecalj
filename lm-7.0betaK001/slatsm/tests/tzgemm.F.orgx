      subroutine fmain
      implicit none
      INTEGER   lda, ldb, ldc, ldw, n, m, k, l, i, ic, count
      double precision mflops
      PARAMETER (lda = 901, ldb = 901, ldc = 901, ldw=901)
C     PARAMETER (lda = 401, ldb = 401, ldc = 401, ldw=401)
C      PARAMETER (n = 600, m = 600, l = 600)
*
      double complex sa( lda, lda ), sb( ldb, ldb ), sc( ldc, ldc )
      double complex pa( lda, lda ), pb( ldb, ldb ), pc( ldc, ldc ),
     .                 pw( ldw, ldw)
      double precision s_time, p_time, temp, second, s_mflops, p_mflops,
     .  q_time, q_mflops
      double precision cpusec
      common /static/ sa, sb, sc, pa, pb, pc

C --- Setup ---
      n = 300
      m = 300
      l = 300
      print *, 'n,m,l?'
      read(*,*) n,m,l
      mflops = 2.0d-6 * n * m * l * 4
      call finits(2,0,0,i)
      call nada
      WRITE (6,100) n, m, n, l, l, m

C --- In-line matrix multiplication ---
      CALL init ( sc,sa,sb,ldc,lda,ldb,n,m,l )
C     call ywrm(0,'a',1,6,'(4f16.10)',sa,1,lda,lda,n)
C     call ywrm(0,'b',1,6,'(4f16.10)',sb,1,ldb,l,m)
      temp = cpusec()
      count = 0
      s_time = 0
    5 CALL scalar_matrix_multiply( sc, sa, sb, ldc, lda, ldb, n, m, l )
      count = count+1
      if (cpusec() - temp .lt. 1d0 .and. count .lt. 10) goto 5
      temp = cpusec()
      print *, 'execute multiply',count,' times'
      do  6  ic = 1, count
    6 CALL scalar_matrix_multiply( sc, sa, sb, ldc, lda, ldb, n, m, l )
      s_time = (cpusec() - temp)/count
      if (s_time .eq. 0d0) s_time = 1
      s_mflops = mflops / s_time
*     call ywrm(0,'c',1,6,'(4f16.10)',sc,1,ldc,n,m)

C --- zgemm matrix multiplication ---
      CALL init ( pc,pa,pb,ldc,lda,ldb,n,m,l )
      temp = cpusec()
      do  7  ic = 1, count
      call zgemm('N','N',n,m,l,(1d0,0d0),pa,lda,pb,ldb,(0d0,0d0),pc,ldc)
c      CALL fast_matrix_multiply  ( pc, pa, pb, ldc, lda, ldb, n, m, l)
    7 continue
c      call mmul(n, m, l, pa, lda, pb, ldb, pc, ldc)
      p_time = (cpusec() - temp)/count
      if (p_time .eq. 0d0) p_time = 1
      p_mflops = mflops / p_time
      CALL diff('compare in-line to zgemm', sc, pc, ldc, n, m )

C --- qmpy matrix multiplication ---
C      CALL init ( pc,pa,pb,ldc,lda,ldb,n,m,l )
C      temp = cpusec()
C      do  9  ic = 1, count
C    9 call qmpy('N','N',n,m,l,pw,ldw,pa,lda,pb,ldb,pc,ldc)
C      q_time = (cpusec() - temp)/count
C      CALL diff('compare in-line to qmpy ', sc, pc, ldc, n, m )
C      if (q_time .eq. 0d0) q_time = 1
C      q_mflops = mflops / q_time

      print 110, s_time, s_mflops, p_time, p_mflops

  100 format( 1X,'Matrix Multiply: ',I3,'x',I3,
     .        ' <- ',I3,'x',I3,' x ',I3,'x',I3)
  110 format(/1X,'Serial time: ',F7.2,'  Serial MFlops: ',F7.2,
     .       /1X,' zgemm Time: ',F7.2,'   zgemm MFlops: ',F7.2)
C     .       /1X,'  qmpy time: ',F7.2,'    qmpy MFlops: ',F7.2)

      end

      SUBROUTINE init ( c, a, b, ldc, lda, ldb, n, m, l )
C- Initialize arrays
      integer lda, ldb, ldc, n, m, l
      double complex a( lda, * ), b( ldb, * ), c( ldc, * )
*
      do i = 1, n
         do j = 1, m
            c(i,j) = dble( i+j ) * 0d0
            end do
         end do
*
      do i = 1, n
         do j = 1, l
            a(i,j) = dble( i-j )
            end do
         end do
*
      do i = 1, l
         do j = 1, m
            b(i,j) = dble( i )/j
            end do
         end do
      end

      SUBROUTINE diff (strn, sc, pc, ldc, n, m )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, i, j
      double complex sc( ldc, * )
      double complex pc( ldc, * )
      DO i = 1, n
         DO j = 1, m
            IF ( ABS( sc(i,j) - pc(i,j) ) .GT. 1.0d-6 ) THEN
               WRITE (6,100) strn
               RETURN
               END IF
            END DO
        END DO
      WRITE (6,110) strn
      RETURN
  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results')
  110 FORMAT(1X,a,'... arrays have the same results')
      END
      SUBROUTINE fast_matrix_multiply(c, a, b, ldc, lda, ldb, n, m, l)
C- Matrix Multiplication with zgemm
      integer lda, ldb, ldc, n, m, l
      double complex a (lda,*), b(ldb,*), c(ldc,*)
      call zgemm ('N','N',n,m,l,(1d0,0d0),a(1,1),lda,b(1,1),ldb,
     .  (0d0,0d0),c(1,1),ldc)
      end

      SUBROUTINE scalar_matrix_multiply(c, a, b, ldc, lda, ldb, n, m, l)
C- In-line matrix multipliation
      integer lda, ldb, ldc, n, m, l
      double complex a (lda,*), b(ldb,*), c(ldc,*)
      do j = 1, m
        do i = 1, n
          c(i,j) = 0
        end do
      end do
      do j = 1, m
        do k = 1, l
          do i = 1, n
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do
      end
C      subroutine qmpy(transa,transb,m,n,k,w,ldw,a,lda,b,ldb,c,ldc)
CC- Quick matrix multiply
CC  See Numerical Recipes, 2.11.
CC  In practice, at best 10% faster than dgemm for SGI, HP, AIX
CC  The same result on 4x4 matrices a,b
CC  mc -f4f16.10 a -split a 1,3,5 1,3,5 b -split b 1,3,5 1,3,5 -x\
CC  a12 a22 -- b21 b22 -+ -x -a q7  a11 a22 -+ b11 b22 -+ -x -a q1\
CC  a21 a22 -+ b11 -x -a q2 a11 b12 b22 -- -x  -a q3\
CC  a22 b21 b11 -- -x -a q4  a11 a12 -+ b22 -x -a q5\
CC  a21 a11 -- b11 b12 -+ -x -a q6  q1 q4 -+ q5 -- q7 -+ -a c11\
CC  q2 q4 -+ -a c21 q3 q5 -+ -a c12 q1 q3 -+ q2 -- q6 -+ -a c22\
CC  c11 c12 -ccat c21 c22 -ccat -rcat
C      implicit none
C      character*1 transa, transb
C      integer m,n,k,lda,ldb,ldc,ldw
C      double complex a(lda,1),b(ldb,1),c(ldc,1),w(ldw,1)
C      integer mmin,npass,m1,m2,i,j
C      parameter(mmin=4)
CC ... debugging
C      character*10 fmt
C      integer ip
C      data ip /0/
C      data fmt /'(4f16.10)'/
C      save ip
C
C      if (m .ne. n .or. m .ne. k)  then
C        if (ip .eq. 0) print *, 'qmpy m must be n,k for now'
C        ip = 1
C        return
C      endif
C
C      if ((transa .ne. 'n' .and. transa .ne. 'N') .or.
C     .    (transb .ne. 'n' .and. transb .ne. 'N'))
C     .  stop 'qmpy: not implemented trans'
C
C      npass = dlog(dble(m)/mmin)/dlog(2d0) + 1d-9
C      if (npass .eq. 0) then
C        call dgemm(transa,transb,m,n,k,1d0,a,lda,b,ldb,0d0,c,ldc)
C        return
C      endif
C
C      m1 = m/2
C      m2 = 2*m1
C
CC --- Q6 in c22 ---
C      do  10  j = 1, m1
C      do  10  i = 1, m1
C      w(i,j)    = a(m1+i,j) - a(i,j)
C   10 w(i+m1,j) = b(i,j) + b(i,j+m1)
C*     call ywrm(0,'a21-a11',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b11+b12',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
C     .  c(1+m1,1+m1),ldc)
C*     call ywrm(0,'Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
C
CC --- Q7 in c11 ---
C      do  20  j = 1, m1
C      do  20  i = 1, m1
C      w(i,j)    = a(i,m1+j) - a(m1+i,m1+j)
C   20 w(i+m1,j) = b(m1+i,j) + b(m1+i,j+m1)
C*     call ywrm(0,'a12-a22',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b21+b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
C     .  c,ldc)
C*     call ywrm(0,'Q7',1,6,fmt,c,1,ldc,m1,m1)
C
CC --- Q1 in c12 ---
C      do  30  j = 1, m1
C      do  30  i = 1, m1
C      w(i,j)    = a(i,j) + a(i+m1,j+m1)
C   30 w(i+m1,j) = b(i,j) + b(i+m1,j+m1)
C*     call ywrm(0,'a11+a22',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b11+b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,w(1+m1,1),ldw,0d0,
C     .  c(1,1+m1),ldc)
C*     call ywrm(0,'Q1',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
C
CC ... c11 += Q1 and c22 += Q1
C      do  40  j = 1, m1
C      do  40  i = 1, m1
C      c(i,j)       = c(i,j)       + c(i,j+m1)
C   40 c(i+m1,j+m1) = c(i+m1,j+m1) + c(i,j+m1)
C*     call ywrm(0,'Q1+Q7',1,6,fmt,c(1,1),1,ldc,m1,m1)
C*     call ywrm(0,'Q1+Q6',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
C
CC --- Q2 in c21, Q4 in c12 ---
C      do  50  j = 1, m1
C      do  50  i = 1, m1
C      w(i,j)    = a(i+m1,j) + a(i+m1,j+m1)
C   50 w(i+m1,j) = b(i+m1,j) - b(i,j)
C*     call ywrm(0,'a21+a22',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b21-b11',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,b,ldb,0d0,
C     .  c(1+m1,1),ldc)
C      call dgemm(transa,transb,m1,m1,m1,1d0,a(1+m1,1+m1),lda,
C     .  w(1+m1,1),ldw,0d0,c(1,1+m1),ldc)
C*     call ywrm(0,'Q2',1,6,fmt,c(1+m1,1),1,ldc,m1,m1)
C*     call ywrm(0,'Q4',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
C
CC ... c11 += Q4, c22 -= Q2, c21 = Q2+Q4 (watch optimizer here)
C      do  60  j = 1, m1
C      do  60  i = 1, m1
C      c(i,j)       = c(i,j)       + c(i,j+m1)
C      c(i+m1,j+m1) = c(i+m1,j+m1) - c(i+m1,j)
C      c(i+m1,j) = c(i+m1,j) + c(i,j+m1)
C   60 continue
C*     call ywrm(0,'Q1+Q7+Q4',1,6,fmt,c(1,1),1,ldc,m1,m1)
C*     call ywrm(0,'Q1+Q6-Q2',1,6,fmt,c(1+m1,1+m1),1,ldc,m1,m1)
C*     call ywrm(0,'c21=Q2+Q4',1,6,fmt,c(1+m1,1),1,ldc,m1,m1)
C
CC --- Q5 in c12, Q3 in w ---
C      do  70  j = 1, m1
C      do  70  i = 1, m1
C      w(i,j)    = a(i,j) + a(i,j+m1)
C   70 w(i+m1,j) = b(i,j+m1) - b(i+m1,j+m1)
C*     call ywrm(0,'a11+a12',1,6,fmt,w,1,ldw,m1,m1)
C*     call ywrm(0,'b12-b22',1,6,fmt,w(1+m1,1),1,ldw,m1,m1)
C      call dgemm(transa,transb,m1,m1,m1,1d0,w,ldw,b(1+m1,1+m1),ldb,0d0,
C     .  c(1,1+m1),ldc)
C      call dgemm(transa,transb,m1,m1,m1,1d0,a,lda,
C     .  w(1+m1,1),ldw,0d0,w,ldw)
C*     call ywrm(0,'Q5',1,6,fmt,c(1,1+m1),1,ldc,m1,m1)
C*     call ywrm(0,'Q3',1,6,fmt,w,1,ldw,m1,m1)
C
CC ... c11 -= Q5, c22 += Q3, c12 = Q3+Q5
C      do  80  j = 1, m1
C      do  80  i = 1, m1
C      c(i,j)       = c(i,j)       - c(i,j+m1)
C      c(i+m1,j+m1) = c(i+m1,j+m1) + w(i,j)
C   80 c(i,j+m1) = c(i,j+m1) + w(i,j)
C
C*     call ywrm(0,'c',1,6,fmt,c,1,ldc,m2,m2)
C      end
