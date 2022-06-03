subroutine fswgts(volwgt,e,ef,etop,w)
  !- Makes weights for integration up to Ef for one tetrahedron.
  ! ----------------------------------------------------------------
  !i Inputs
  !i
  !o Outputs
  !o   w
  !r Remarks
  !r   w(i,1): normal weights.  w(i,2): Bloechl-correction.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  double precision :: e(4),ef,volwgt,etop,w(4,2)
  ! Local parameters
  integer :: i,i00,j,m,n,isort(4)
  double precision :: ex(4),w1(4),w2(4),ec(4), &
       a14,a21,a24,a31,a32,a34,a41,a42,df=0d0,e1,e2,e3,e4,v1,v2,v3,vw4,xxx
  vw4 = volwgt/4d0
  w=0d0
  if(ef >= etop) then
     w(:,1)=vw4
     return
  endif
  ! --- Sort energies into ec ---
  ex = e
  do  3  i = 1, 4
     i00 = 1
     do  4  j = 2, 4
        if (ex(j) < ex(i00)) i00 = j
4    enddo
     ec(i) = ex(i00)
     isort(i) = i00
     ex(i00) = etop + 1d0
3 enddo
  e1 = ec(1)
  e2 = ec(2)
  e3 = ec(3)
  e4 = ec(4)
  ! --- Case Ef between e2,e3 ---
  if (e2 < ef .AND. ef <= e3) then
     a31 = (ef-e1)/(e3-e1)
     a41 = (ef-e1)/(e4-e1)
     a32 = (ef-e2)/(e3-e2)
     a42 = (ef-e2)/(e4-e2)
     v1 = a31*a41
     v2 = a31*a42*(1-a41)
     v3 = a42*a32*(1-a31)
     w1(1) = (v1*(3-a31-a41) + v2*(2-a31-a41) + v3*(1-a31))*vw4
     w1(2) = (v1 + v2*(2-a42) + v3*(3-a32-a42))*vw4
     w1(3) = (v1*a31 + v2*a31 + v3*(a31+a32))*vw4
     w1(4) = (v1*a41 + v2*(a41+a42) + v3*a42)*vw4
     df = ((e1+e2-e3-e4)*a32*a42 + 2*ef - e1 - e2)/((e3-e1)*(e4-e1))
     df = 3*volwgt*df
     ! --- Case Ef between e1,e2 ---
  else if (e1 < ef .AND. ef <= e2) then
     a21 = (ef-e1)/(e2-e1)
     a31 = (ef-e1)/(e3-e1)
     a41 = (ef-e1)/(e4-e1)
     xxx = a21*a31*a41*vw4
     w1(1) = xxx*(4-a21-a31-a41)
     w1(2) = xxx*a21
     w1(3) = xxx*a31
     w1(4) = xxx*a41
     df = 3*volwgt*a31*a41/(e2-e1)
     ! --- Case Ef between e3,e4 ---
  else if (e3 < ef .AND. ef <= e4) then
     a14 = (ef-e4)/(e1-e4)
     a24 = (ef-e4)/(e2-e4)
     a34 = (ef-e4)/(e3-e4)
     xxx = a14*a24*a34*vw4
     w1(1) = vw4 - xxx*a14
     w1(2) = vw4 - xxx*a24
     w1(3) = vw4 - xxx*a34
     w1(4) = vw4 - xxx*(4-a14-a24-a34)
     df = -3*volwgt*a14*a24/(e3-e4)
  endif
  ! --- Bloechl correction ---
  w2=0d0
  do   m = 1, 4
     do  n = 1, 4
        w2(m) = w2(m) + (ec(n)-ec(m))*df/40
     enddo
  enddo
  ! ----------------------------------------------
  do  35  i = 1, 4
     j = isort(i)
     w(j,1) = w1(i)
     w(j,2) = w2(i)
35 enddo
end subroutine fswgts

