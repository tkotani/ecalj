      SUBROUTINE GINTxx(G1,G2,A,B,NR,SUM)
C- Integrate product of two wave functions, Simpson rule
C ----------------------------------------------------------------
c takao \sum_i g(i)*g(i) drdi
      implicit none
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer:: ir,nr
      real(8):: G1(NR),G2(NR),a,b,sum,dedi,ea2,ea4,rmpb,drdi
      EA2 = DEXP(A + A)
      EA4 = EA2*EA2
      SUM = 0D0
      DRDI = (A*B)*DEXP(A)
      DO  10  IR = 2, NR-1, 2
        SUM = SUM + G1(IR)*G2(IR)*DRDI
        DRDI = DRDI*EA2
  10  CONTINUE
      SUM = SUM + SUM
      DRDI = (A*B)*EA2
      DO  11  IR = 3, NR-2, 2
        SUM = SUM + G1(IR)*G2(IR)*DRDI
        DRDI = DRDI*EA2
  11  CONTINUE
      RMPB = B*DEXP(A*(NR-1))
      SUM = (2*SUM+ G1(1)*G2(1)*(A*B) + G1(NR)*G2(NR)*(A*RMPB))/3d0
      END

      function gintegral(G1,G2,A,B,NR) result(sum)
C- Integrate product of two wave functions, Simpson rule
C ----------------------------------------------------------------
c takao \sum_i g(i)*g(i) drdi
      implicit none
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer:: ir,nr
      real(8):: G1(NR),G2(NR),a,b,sum,dedi,ea2,ea4,rmpb,drdi
      EA2 = DEXP(A + A)
      EA4 = EA2*EA2
      SUM = 0D0
      DRDI = (A*B)*DEXP(A)
      DO  10  IR = 2, NR-1, 2
        SUM = SUM + G1(IR)*G2(IR)*DRDI
        DRDI = DRDI*EA2
  10  CONTINUE
      SUM = SUM + SUM
      DRDI = (A*B)*EA2
      DO  11  IR = 3, NR-2, 2
        SUM = SUM + G1(IR)*G2(IR)*DRDI
        DRDI = DRDI*EA2
  11  CONTINUE
      RMPB = B*DEXP(A*(NR-1))
      SUM = (2*SUM+ G1(1)*G2(1)*(A*B) + G1(NR)*G2(NR)*(A*RMPB))/3d0
      END
      
