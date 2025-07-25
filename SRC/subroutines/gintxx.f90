SUBROUTINE gintxx(G1,G2,A,B,NR,SUM) !- Integrate product of two wave functions, Simpson rule
  ! takao \sum_i g(i)*g(i) drdi
  implicit none
  integer:: ir,nr
  real(8):: G1(NR),G2(NR),a,b,sum,dedi,ea2,ea4,rmpb,drdi
  EA2 = DEXP(A + A)
  EA4 = EA2*EA2
  SUM = 0D0
  DRDI = (A*B)*DEXP(A)
  DO  10  IR = 2, NR-1, 2
     SUM = SUM + G1(IR)*G2(IR)*DRDI
     DRDI = DRDI*EA2
10 enddo
  SUM = SUM + SUM
  DRDI = (A*B)*EA2
  DO  11  IR = 3, NR-2, 2
     SUM = SUM + G1(IR)*G2(IR)*DRDI
     DRDI = DRDI*EA2
11 enddo
  RMPB = B*DEXP(A*(NR-1))
  SUM = (2*SUM+ G1(1)*G2(1)*(A*B) + G1(NR)*G2(NR)*(A*RMPB))/3d0
end SUBROUTINE gintxx
function gintegral(G1,G2,A,B,NR) result(sum) !- Integrate product of two wave functions, Simpson rule
  ! takao \sum_i g(i)*g(i) drdi
  implicit none
  integer:: ir,nr
  real(8):: G1(NR),G2(NR),a,b,sum,dedi,ea2,ea4,rmpb,drdi
  EA2 = DEXP(A + A)
  EA4 = EA2*EA2
  SUM = 0D0
  DRDI = (A*B)*DEXP(A)
  DO  10  IR = 2, NR-1, 2
     SUM = SUM + G1(IR)*G2(IR)*DRDI
     DRDI = DRDI*EA2
10 enddo
  SUM = SUM + SUM
  DRDI = (A*B)*EA2
  DO  11  IR = 3, NR-2, 2
     SUM = SUM + G1(IR)*G2(IR)*DRDI
     DRDI = DRDI*EA2
11 enddo
  RMPB = B*DEXP(A*(NR-1))
  SUM = (2*SUM+ G1(1)*G2(1)*(A*B) + G1(NR)*G2(NR)*(A*RMPB))/3d0
end function gintegral

module m_gint
  public:: gintpp
contains
  pure function gintpp(G1,G2,A,B,NR) result(sum) !- Integrate product of two wave functions, Simpson rule
    ! takao \sum_i g(i)*g(i) drdi
    implicit none
    intent(in):: g1,g2,a,b,nr
    integer:: ir,nr
    real(8):: G1(NR),G2(NR),a,b,sum,dedi,ea2,ea4,rmpb,drdi
    EA2 = DEXP(A + A)
    EA4 = EA2*EA2
    SUM = 0D0
    DRDI = (A*B)*DEXP(A)
    DO  10  IR = 2, NR-1, 2
      SUM = SUM + G1(IR)*G2(IR)*DRDI
      DRDI = DRDI*EA2
10  enddo
    SUM = SUM + SUM
    DRDI = (A*B)*EA2
    DO  11  IR = 3, NR-2, 2
      SUM = SUM + G1(IR)*G2(IR)*DRDI
      DRDI = DRDI*EA2
11  enddo
    RMPB = B*DEXP(A*(NR-1))
    SUM = (2*SUM+ G1(1)*G2(1)*(A*B) + G1(NR)*G2(NR)*(A*RMPB))/3d0
  end function gintpp
endmodule m_gint
