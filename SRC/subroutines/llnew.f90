module m_ll
  contains
  pure integer function ll(ilm) ! Returns l, given lm index
    intent(in)::           ilm
    integer ::  ilm,m,l,ioff
    integer,parameter:: lmaxx=24
    integer,parameter:: lla(1:(lmaxx+1)**2) = [((l,m=1,2*l+1),l=0,lmaxx)] !nvfortran24.1 did not work
    !integer,save:: lla(1:(lmaxx+1)**2)
    !logical,save:: init=.true.
    !if(init) then
    !   do l=0,lmaxx
    !      ioff=l**2
    !      lla(ioff+1:ioff+2*l+1)=l
    !   enddo
    !   init=.false.
    !endif
    ll = lla(ilm)
  end function ll
endmodule
