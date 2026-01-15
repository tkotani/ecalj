module m_bz_integ
  implicit none
  private
  public :: ibz_integ
contains
  pure real(8) function tetrahedron_volume(p1, p2, p3, p4) result(volume)
    real(8), intent(in) :: p1(3), p2(3), p3(3), p4(3)
    real(8) :: v1(3), v2(3), v3(3)
    v1 = p1 - p4
    v2 = p2 - p4
    v3 = p3 - p4
    volume = v1(1)*(v2(2)*v3(3) - v2(3)*v3(2)) &
           + v1(2)*(v2(3)*v3(1) - v2(1)*v3(3)) &
           + v1(3)*(v2(1)*v3(2) - v2(2)*v3(1))
    volume = abs(volume)/6.0d0
  endfunction tetrahedron_volume
  pure complex(8) function ibz_integ(num_tetra, tetra_nodes, tetra_weight, qibz, nqibz, func_values) result(integral_value)
    integer, intent(in) :: num_tetra, nqibz
    integer, intent(in) :: tetra_nodes(4, num_tetra), tetra_weight(num_tetra)
    real(8), intent(in) :: qibz(3,nqibz)
    complex(8), intent(in) :: func_values(nqibz)
    integer :: i, n1, n2, n3, n4
    complex(8) :: avg_func_val
    real(8) :: p1(3), p2(3), p3(3), p4(3), tet_vol
    integral_value = (0.0d0, 0d0)
    do i = 1, num_tetra
      n1 = tetra_nodes(1,i)
      n2 = tetra_nodes(2,i)
      n3 = tetra_nodes(3,i)
      n4 = tetra_nodes(4,i)
      p1 = qibz(:,n1)
      p2 = qibz(:,n2)
      p3 = qibz(:,n3)
      p4 = qibz(:,n4)
      tet_vol = tetrahedron_volume(p1, p2, p3, p4)
      avg_func_val = (func_values(n1) + func_values(n2) + &
                      func_values(n3) + func_values(n4))/4.0d0
      integral_value = integral_value + avg_func_val*tet_vol*tetra_weight(i)
    enddo
  endfunction
endmodule m_bz_integ
