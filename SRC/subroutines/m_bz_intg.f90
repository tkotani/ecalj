module m_bz_integ
  implicit none
  private
  public :: ibz_integ
  interface ibz_integ
    module procedure ibz_zinteg, ibz_dinteg
  end interface
contains
  pure complex(8) function ibz_zinteg(tetra_vol, num_tetra, tetra_nodes, tetra_weight, qibz, nqibz, func_values) result(integral_value)
    integer, intent(in) :: num_tetra, nqibz
    integer, intent(in) :: tetra_nodes(4, num_tetra), tetra_weight(num_tetra)
    real(8), intent(in) :: qibz(3,nqibz), tetra_vol
    complex(8), intent(in) :: func_values(nqibz)
    integer :: i, n1, n2, n3, n4
    complex(8) :: avg_func_val
    integral_value = (0.0d0, 0d0)
    do i = 1, num_tetra
      n1 = tetra_nodes(1,i)
      n2 = tetra_nodes(2,i)
      n3 = tetra_nodes(3,i)
      n4 = tetra_nodes(4,i)
      avg_func_val = (func_values(n1) + func_values(n2) + &
                      func_values(n3) + func_values(n4))/4.0d0
      integral_value = integral_value + avg_func_val*tetra_vol*tetra_weight(i)
    enddo
  endfunction
  pure real(8) function ibz_dinteg(tetra_vol, num_tetra, tetra_nodes, tetra_weight, qibz, nqibz, func_values) result(integral_value)
    integer, intent(in) :: num_tetra, nqibz
    integer, intent(in) :: tetra_nodes(4, num_tetra), tetra_weight(num_tetra)
    real(8), intent(in) :: qibz(3,nqibz), tetra_vol
    real(8), intent(in) :: func_values(nqibz)
    integer :: i, n1, n2, n3, n4
    real(8) :: avg_func_val
    integral_value = 0.0d0
    do i = 1, num_tetra
      n1 = tetra_nodes(1,i)
      n2 = tetra_nodes(2,i)
      n3 = tetra_nodes(3,i)
      n4 = tetra_nodes(4,i)
      avg_func_val = (func_values(n1) + func_values(n2) + &
                      func_values(n3) + func_values(n4))/4.0d0
      integral_value = integral_value + avg_func_val*tetra_vol*tetra_weight(i)
    enddo
  endfunction
endmodule m_bz_integ
