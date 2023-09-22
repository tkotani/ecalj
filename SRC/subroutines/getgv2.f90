subroutine getgv2(alat,plat,qlat,q, QpGcut,job, ng,ngvec) ! Set up a list of recip vectors within cutoff |q+G| < QpGcut a.u.
  ! job==1 -> return ng (number of G ) only
  ! job==2 -> return ng and ngvec
  ! True G is given by  G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
  implicit none
  integer :: ng
  real(8):: s_lat(1),q(3),plat(3,3),qlat(3,3),qpg(3),enor(3),pi,alat,tpiba,QpGmax,QpGmax2,QpGcut,Qenor
  integer :: n1max,n1min,n2max,n2min,n3max,n3min, i1,i2,i3,ig,job,imx,ngvec(3,ng)
  pi=4d0*datan(1d0)
  tpiba=2*pi/alat
  QpGmax   = QpGcut/tpiba  ! QpGcut in a.u.= tpiba*Qcut
  QpGmax2  = QpGmax**2
  call eprod(qlat(1:3,2),qlat(1:3,3),enor);
  Qenor = sum(qlat(1:3,1)*enor); n1max = QpGmax/abs(Qenor) -sum(q*enor)/Qenor +1;n1min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1
  call eprod(qlat(1:3,3),qlat(1:3,1),enor)
  Qenor = sum(qlat(1:3,2)*enor); n2max = QpGmax/abs(Qenor) -sum(q*enor)/Qenor +1;n2min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1
  call eprod(qlat(1:3,1),qlat(1:3,2),enor)
  Qenor = sum(qlat(1:3,3)*enor); n3max = QpGmax/abs(Qenor) -sum(q*enor)/Qenor +1;n3min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1
  ig=0
  imx=-9999
  GETngvecWithinTheLimit: do i1 = n1min, n1max
     do i2 = n2min, n2max
        do i3 = n3min, n3max
           qpg(1:3)= q(1:3) + matmul(qlat,[i1,i2,i3]) 
           if( sum(qpg(1:3)**2) < QpGmax2) then
              ig = ig+1
              if(job==2) ngvec(1:3,ig) = [i1,i2,i3]
              if(job==1) imx=max(imx,abs(i1),abs(i2),abs(i3))
           endif
        enddo
     enddo
  enddo GETngvecWithinTheLimit
  ng = ig
  if(job==1) ngvec(1,1)=imx 
end subroutine getgv2
