subroutine slinz(volwgt,ec,emin,emax,dosi,nr)! Adds to number-of-states for one tetrahedron
  !i Inputs
  !i   volwgt:weight on tetrahedron
  !i   ec    :energies at corners of tetrahedron
  !i   emin, emax, energy window
  !i   nr    :number of bins + 1
  !o Outputs
  !o   dosi(k), integrated density in kth bin from tethdn.
  !o   ec:   sorted on output
  implicit none
  integer :: nr, i,i01,i02,i03,i04,i1,i2,i3,i4,j,i2p,i3p,i4p
  real(8) :: c0,c1,c2,c3,cc,de,e,e1,e2,e3,e4,x,mxmin,xx,cc4,ec(4),dosi(nr),volwgt,emin,emax
!  parameter (fuzz=1d-8)
!  mxmin(xx) = min(max(xx,-2*de),nr*de)/de+1.9999999d0 !     Guard against overflow, underflow. Make -1 < i0[1..4] < nr+1
  mxmin(xx) = xx/de+1.9999999d0 
  ! --- Sort the ec ---
  do   i = 1, 3
     do   j = 1, 4-i
        if (ec(j) > ec(j+1)) then
           e = ec(j)
           ec(j) = ec(j+1)
           ec(j+1) = e
        endif
     enddo
  enddo
  e1 = ec(1) !lowest
  e2 = ec(2)
  e3 = ec(3)
  e4 = ec(4) !highest
  de = (emax-emin)/(nr-1)
  i01 = mxmin(e1-emin)
  i02 = mxmin(e2-emin)
  i03 = mxmin(e3-emin)
  i04 = mxmin(e4-emin)
  i1 = max0(i01,1)
  i2 = min0(i02-1,nr)
  i2p = max0(i02,1)
  i3 = min0(i03-1,nr)
  i3p = max0(i03,1)
  i4 = min0(i04-1,nr)
  i4p = max0(i04,1)
  if (i1 <= i2) then
     cc = volwgt/((e2-e1)*(e3-e1)*(e4-e1))
     dosi(i1:i2) = dosi(i1:i2) + [(cc*(emin-e1+(i-1)*de)**3,i=i1,i2)]
  endif
  if (i2p <= i3) then !     print *,'eeeeee=',e1,e2,e3,e4
     c3 = volwgt*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
     c2 = volwgt*3d0/((e3-e1)*(e4-e1))
     c1 = c2*(e2-e1)
     c0 = c1*(e2-e1)/3d0
     do i = i2p, i3
        x = emin - e2 + (i-1)*de
        dosi(i) = dosi(i) + c0 + x*(c1 + x*(c2 + x*c3))
     enddo
  endif
  if (i3p <= i4) then
     cc4 = volwgt/((e3-e4)*(e2-e4)*(e1-e4))
     do i = i3p, i4
        x = emin - e4 + (i-1)*de
        dosi(i) = dosi(i) + volwgt - cc4*x**3
     enddo
  endif
  dosi(i4p:nr) = dosi(i4p:nr) + volwgt
end subroutine slinz
