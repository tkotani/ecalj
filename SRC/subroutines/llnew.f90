module m_ll
  contains
  pure integer function ll(ilm) ! Returns l, given lm index
    intent(in)::           ilm
    integer ::  ilm,m,l
    integer,parameter:: lmaxx=24, lla(1:(lmaxx+1)**2) = [((l,m=1,2*l+1),l=0,lmaxx)]
    ll = lla(ilm)
  end function ll
endmodule
