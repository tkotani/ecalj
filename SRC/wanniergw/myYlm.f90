
! Calculate rotation matrix D_{MK}^{J}(\alpha,\beta,\gamma)
! alpha, beta, and gamma are Euler angles

! R(a,b,g) YLM = \sum_{M'} D_{M' M}^{J}(a,b,g) YLM'

subroutine calc_D(J,alpha,beta,gamma,D)
  implicit none
  integer,intent(in) :: J
  double precision,intent(in) :: alpha,beta,gamma
  double complex :: D(2*J+1,2*J+1)
  ! local
  double precision :: DJbeta(2*J+1,2*J+1)
  integer :: i1,i2,M,K
  double complex :: I
  parameter (I=dcmplx(0.0d0,1.0d0))
  call calc_DJbeta(J,beta,DJbeta)
  do i2=1,2*J+1
     K=i2-J-1
     do i1=1,2*J+1
        M=i1-J-1
        D(i1,i2)=exp(-I*M*alpha)*DJbeta(i1,i2) &
             *exp(-I*K*gamma)
        !          write(*,"(a,i3,3f5.2,2i3,2f10.5)")
        !     &         'J,alpha,beta,gamma,M,K,D_{MK}=',
        !     &         J,alpha,beta,gamma,M,K,D(i1,i2)
     enddo
  enddo

end subroutine calc_D

! Calculate rotation matrix D_{MK}^{J}(\beta)

subroutine calc_DJbeta(J,beta,DJbeta)
  implicit none
  integer,intent(in) :: J
  double precision,intent(in) :: beta
  double precision :: DJbeta(2*J+1,2*J+1)

  double precision :: zeta
  integer :: i1,i2,M,K,l
  double precision :: tmp1,tmp2,tmp3
  ! function
  integer :: myfact,mycomb
  double precision :: mydif
  !c
  zeta=sin(0.5d0*beta)**2
  do i2=1,2*J+1
     K=i2-J-1
     do i1=1,2*J+1
        M=i1-J-1
        tmp1=dble(myfact(J+M))
        tmp2=dble(myfact(J-M)*myfact(J+K)*myfact(J-K))* &
             zeta**(M-K)*(1.0d0-zeta)**(M+K)

        tmp3=0.0d0
        do l=0,J+K
           tmp3=tmp3+(-1.0d0)**(J+K-l)*mycomb(J+K,l)* &
                mydif(zeta,2*J-l,J-M)
        enddo
        !          DJbeta(i1,i2)=(-1.0d0)**(M-K)*sqrt(tmp1/tmp2)*tmp3
        DJbeta(i1,i2)=(-1.0d0)**(M-K)*sqrt(tmp1/(tmp2+1.0d-20))*tmp3
     enddo
  enddo
end subroutine calc_DJbeta
! cccccccccccccccccccc
! Calculate (complex) Ylm(theta,phi)
! theta=beta,phi=alpha
subroutine calc_Ylm(L,beta,alpha,Y,Yreal)
  implicit none
  integer,intent(in) :: L
  double precision ,intent(in) :: beta,alpha
  double complex :: Y(2*L+1)
  double precision :: Yreal(2*L+1)
  ! local
  double precision :: pi
  double complex :: D(2*L+1,2*L+1)
  integer :: i1,m
  pi=4.0d0*atan(1.0d0)
  call calc_D(L,alpha,beta,0.0d0,D)
  do i1=1,2*L+1
     Y(i1)=sqrt(dble(2*L+1)/(4.0d0*pi))*dconjg(D(i1,L+1))
  enddo
  do m=-L,L
     if (m < 0) then
        Yreal(m+L+1)=(-1.0d0)**(-m)*sqrt(2.0d0)*dimag(Y(-m+L+1))
     elseif (m == 0) then
        Yreal(m+L+1)=dble(Y(m+L+1))
     else
        Yreal(m+L+1)=(-1.0d0)**(m)*sqrt(2.0d0)*dble(Y(m+L+1))
     endif
  enddo
end subroutine calc_Ylm

! cccccccccccccccccccc
subroutine calc_B(L,B,Binv)
  implicit none
  integer,intent(in) :: L
  double complex :: B(2*L+1,2*L+1),Binv(2*L+1,2*L+1)
  ! local
  integer :: i1,i2,m1,m2,i3
  double precision :: sq2
  double complex :: I
  parameter (I=dcmplx(0.0d0,1.0d0))

  double complex :: BBtmp

  sq2=sqrt(2.0d0)
  B(1:2*L+1,1:2*L+1)=0.0d0
  Binv(1:2*L+1,1:2*L+1)=0.0d0

  do i1=1,2*L+1
     m1=i1-L-1
     if (m1 < 0) then
        B(i1,i1)=-1.0d0/(sq2*I)
        B(-m1+L+1,i1)=(-1.0d0)**(-m1)/(sq2*I)

        Binv(i1,i1)=-I/(sq2)
        Binv(-m1+L+1,i1)=1.0d0/(sq2)

     elseif (m1 == 0) then
        B(i1,i1)=1.0d0
        Binv(i1,i1)=1.0d0
     else
        B(i1,i1)=(-1.0d0)**(m1)/sq2
        B(-m1+L+1,i1)=1.0d0/sq2

        Binv(i1,i1)=(-1.0d0)**(m1)/sq2
        Binv(-m1+L+1,i1)=(-1.0d0)**(m1)/sq2*I
     endif
  enddo
  do i1=1,2*L+1
     do i2=1,2*L+1
        BBtmp=0.0d0
        do i3=1,2*L+1
           BBtmp=BBtmp+B(i3,i1)*Binv(i2,i3)
        enddo
        write(*,"(a,2i3,2f10.5)") 'i1,i2,B*Binv=',i1,i2,BBtmp
     enddo
  enddo
end subroutine calc_B
! cccccccccccccccccccc
