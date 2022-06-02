! Taken from Ferdi's GW  -------------------------------------------------------------------
subroutine genqbz (qbas,n1,n2,n3, &
     qbz,wbz, nstbz, nadd,half)
  ! generates the k-points in the 1BZ
  ! the 1BZ is a parallepiped formed by G1,G2,G3 (qbas(3,3))
  ! this is divided into microcells defined by G1/n1,G2/n2,G3/n3
  ! the k-points may be thought of as being centred at each microcell
  ! the sampling weight for each k-point is the same (1/n1*n2*n3)

  ! qbas = base reciprocal vectors G1,G2,G3
  ! n1,n2,n3 = divisions along G1,G2,G3

  ! qbz  = k-points in the 1BZ
  ! wbz  = sampling weight for qbz

  implicit real*8 (a-h,o-z)
  real(8):: qbas(3,3)
  real(8):: qbz(3,*),wbz(*) !wbz(n1*n2*n3)
  real(8):: qmic(3,3),w1(3),w2(3),w3(3),half(3)
  integer:: nstbz(*),nadd,n1,n2,n3,i1,i2,i3,kount,nnn
  !! icase=1
  !! vectors forming microcells
  qmic(:,1)= qbas(:,1)/dble(n1)
  qmic(:,2)= qbas(:,2)/dble(n2)
  qmic(:,3)= qbas(:,3)/dble(n3)
  nnn=(n1+nadd)*(n2+nadd)*(n3+nadd)
  nstbz(1:nnn)=0
  !! sampling weight
  weight     = 1.d0/dble(n1*n2*n3)
  kount      = 0
  do      i1 = 1,n1+nadd
     do      i2 = 1,n2+nadd
        do      i3 = 1,n3+nadd
           kount      = kount + 1
           qbz(:,kount) = (i1-1+half(1))*qmic(:,1) + (i2-1+half(2))*qmic(:,2) + (i3-1+half(3))*qmic(:,3)
           wbz(kount) = weight
           !        if(icase==2.and.(i1==1.or.i1==n1).and.(i2==1.or.i2==n2).and.(i3==1.or.i3==n3)) then
           !          nstbz(kount) = 2*2*2
           !        endif
           !        write(6,"(' iq qbz =',i5,3f9.4)")kount,qbz(:,kount)
        enddo
     enddo
  enddo
  if(kount /= (n1+nadd)*(n2+nadd)*(n3+nadd))call rx( 'genqbz: wrong no. k-points')
end subroutine genqbz
!--------------------
!      subroutine cv (c,v,n,
!     o w )
! forms w(i) = c * v(i)

!      implicit real*8(a-h,o-z)
!      dimension v(n)
!      dimension w(n)

!      do       i = 1,n
!      w(i)       = c*v(i)
!      end do

!      return
!      end
!----------------

