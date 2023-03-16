module m_ll
  contains
  pure integer function ll(ilm) ! Returns l, given lm index
    intent(in):: ilm
    integer :: m,l,ilm
    integer,parameter:: lmaxx=16
    integer,parameter :: lla(1:(lmaxx+1)**2) = [((l,m=1,2*l+1),l=0,16)]
    ll = lla(ilm)
  end function ll
endmodule
