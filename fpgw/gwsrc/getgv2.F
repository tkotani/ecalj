      subroutine getgv2(alat,plat,qlat,q, QpGcut,job,
     o                 ng, ngvec)
!! == Set up a list of recip vectors within cutoff |Q+G| < QpGcut a.u. ==
!! job==1 -> return ng (number of G ) and imx(as ngvec(1,1));mar2012takao add imx.
!! job==2 -> return ng and ngvec
!! True G is given as
!!    G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
!! NOTE: we need some geometorial consideration for this routine.
!!   Consiser ellipsoid. Takao need to give more detailed explanation...
!! -------------------------------------------------------------
      implicit none
      integer(4) :: ng
      real(8):: s_lat(1),q(3),plat(3,3),qlat(3,3),qpg(3),enor(3)
     & ,pi,alat,tpiba,QpGmax,QpGmax2,QpGcut,Qenor
      integer(4) ::  ngvec(3,ng),
     &    n1max,n1min,n2max,n2min,n3max,n3min,
     &    i1,i2,i3,ig,job,imx
      pi=4d0*datan(1d0)
c      call u_lat_vecs(s_lat,alat,plat,qlat)
      tpiba=2*pi/alat
      QpGmax   = QpGcut/tpiba  ! QpGcut in a.u.= tpiba*Qcut
      QpGmax2  = QpGmax**2
c search range of n1 n2 n3
      call eprod(qlat(1:3,2),qlat(1:3,3),enor)
      Qenor = sum(qlat(1:3,1)*enor)
      n1max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
      n1min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1

      call eprod(qlat(1:3,3),qlat(1:3,1),enor)
      Qenor = sum(qlat(1:3,2)*enor)
      n2max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
      n2min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1

      call eprod(qlat(1:3,1),qlat(1:3,2),enor)
      Qenor = sum(qlat(1:3,3)*enor)
      n3max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
      n3min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1

c get ngvec within the limit.
      ig=0
      imx=-9999
      do i1 = n1min, n1max
        do i2 = n2min, n2max
          do i3 = n3min, n3max
ccccccccc TEST
c      do i1 = n1min-3, n1max+3
c      do i2 = n2min-3, n2max+3
c      do i3 = n3min-3, n3max+3
            qpg(1:3)= q(1:3) + 
     &    qlat(1:3,1)*i1 +qlat(1:3,2)*i2 +qlat(1:3,3)*i3
            if( sum(qpg(1:3)**2) < QpGmax2) then
              ig = ig+1
              if(job==2) ngvec(1:3,ig) = (/i1,i2,i3/)
              if(job==1) imx=max(imx,abs(i1),abs(i2),abs(i3))
c          if(job==2) write(1116,'(f8.4,3i3)')
c     &      tpiba*sqrt(sum(qpg(1:3)**2)),ngvec(1:3,ig) ! check write
            endif
          enddo
        enddo
      enddo
      ng = ig
      if(job==1) ngvec(1,1)=imx !mar2012takao
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine eprod(a,b, c)
c- c gives normalized normal vector for a and b.
      real(8) :: a(3),b(3),c(3),cnorm !anorm,bnorm,
      c(1)= a(2)*b(3)-a(3)*b(2)
      c(2)= a(3)*b(1)-a(1)*b(3)
      c(3)= a(1)*b(2)-a(2)*b(1)
c      anorm = sqrt(sum(a(1:3)**2))
c      bnorm = sqrt(sum(b(1:3)**2))
c      c = c/anorm/bnorm
      cnorm = sqrt(sum(c(1:3)**2))
      c = c/cnorm
      end

