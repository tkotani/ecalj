subroutine Veecomp(Vee,l,U_H,J_H)
  !- Calculate Vee from U, J, Slater integrals, and Clebsch-Gordon coefficients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   l     :l for which U_H and J_H are defined
  !i   U_H   :Screened direct coulomb integral (Hubbard U)
  !i   J_H   :Exchange integral (Hubbard J)
  !o Outputs
  !o   Vee   :Vee(m,m2,m3,m4) = <m,m2|Vee|m1,m3>
  !l Local variables
  !r Remarks
  !u Updates
  !u   02 Jun 05 Lambrecht bug fixes
  !u   27 Apr 05 Larson first created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: l
  double precision :: U_H,J_H
  double precision :: Vee(-3:3,-3:3,-3:3,-3:3)
  ! ... Local parameters
  integer :: m,m1,m2,m3,k,jj,sg
  double precision :: F(0:2*3),r42,r62,t1,t2,t3,ak
  !     Slater integrals from U_H, J_H (this will be updated later)

  call dpzero(vee,7**4)
  call dpzero(F,7)

  F(0) = U_H
  !     d-states -> r42 = F(4)/F(2) = 0.625
  if (l == 2) then
     r42 = 0.625d0
     F(2) = 14d0*J_H/(1d0+r42)
     F(4) = F(2)*r42
  endif
  !     f-states
  if (l == 3) then
     r42 = 451d0/675d0
     r62 = 1001d0/2025d0
     F(2) = 3d0*J_H/(2d0/15d0+1d0/11d0*r42+50d0/429d0*r62)
     F(4) = F(2)*r42
     F(6) = F(2)*r62
  endif
  !     Calculate <m,m2|Vee|m1,m3>
  do  m = -l, l
     do  m2 = -l, l
        do  m1 = -l, l
           do  m3 = -l, l
              do  k = 0, 2*l,2
                 !                     print *,'k=',k,'F=',F(k)
                 ak = 0d0
                 call t3j_all(l,k,l,0,0,0,t1)
                 !                     print *,'l=',l,'k=',k,'F=',F(k),'t1=',t1
                 do  jj = -k, k
                    call t3j_all(l,k,l,-m,jj,m1,t2)
                    call t3j_all(l,k,l,-m2,-jj,m3,t3)
                    sg = (-1)**(m+jj+m2)
                    ak = ak + sg*t2*t3
                 enddo
                 ak = (2*l+1)**2*t1**2*ak
                 !                     print '(5i3,3f8.2)',
                 !     .               m,m2,m1,m3,k,Vee(m,m2,m1,m3),ak,F(k)
                 Vee(m,m2,m1,m3) = Vee(m,m2,m1,m3) + ak*F(k)
                 !                     print '(5i3,3f9.4)',
                 !     .               m,m2,m1,m3,k,Vee(m,m2,m1,m3),ak,F(k)
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine Veecomp

subroutine t3j_all(j1,j2,j3,m1,m2,m3,t3j)
  ! ----------------------------------------------------------------
  !- Calculate 3-j symbols.
  !     implicit none
  integer :: j1,j2,j3,m1,m2,m3
  integer :: J,L1,L2,K0,K1,K2,KMAX,KMIN,trr,k
  double precision :: sgJ,sgk,C,ab,t3j
  double precision :: a1,a2,a3,a4,a5,a6,a7,a8,a9,b
  double precision :: a1a,a2a,a3a,a4a,a5a,a6a,a7a,a8a,a9a


  t3j = 0d0
  if ((m1+m2+m3) /= 0) goto 142
  call tri_rule(j1,j2,j3,trr)
  if (trr < 0) goto 142


  call fctor1( j1 + j2 - j3,a1)
  call fctor1( j1 - j2 + j3,a2)
  call fctor1(-j1 + j2 + j3,a3)
  call fctor1( j1 + m1,a4)
  call fctor1( j1 - m1,a5)
  call fctor1( j2 + m2,a6)
  call fctor1( j2 - m2,a7)
  call fctor1( j3 + m3,a8)
  call fctor1( j3 - m3,a9)
  call fctor1( j1 + j2 + j3 + 1,b)
  ab = b
  a1a = a1
  a2a = a2
  a3a = a3
  a4a = a4
  a5a = a5
  a6a = a6
  a7a = a7
  a8a = a8
  a9a = a9
  J   = j1 - j2 - m3
  sgJ = (-1)**(J)
  C   = sgJ*DSQRT(a1a*a2a*a3a*a4a*a5a*a6a*a7a*a8a*a9a/ab)

  K0   = j1 + j2 - j3
  K1   = j1 - m1
  K2   = j2 + m2
  L1   = j2 - j3 - m1
  L2   = j1 - j3 + m2
  KMAX = min0(K1,K2,K0)
  KMIN = max0(0,L1,L2)

  t3j = 0d0
  do  k = KMIN, KMAX
     call fctor1(k,a1)
     call fctor1(j1 +j2 - j3 - k,a2)
     call fctor1(j1 -m1 - k,a3)
     call fctor1(j2 +m2 - k,a4)
     call fctor1(j3 -j2 + m1 + k,a5)
     call fctor1(j3 -j1 - m2 + k,a6)
     sgk = (-1d0)**k
     a1a = a1
     a2a = a2
     a3a = a3
     a4a = a4
     a5a = a5
     a6a = a6
     t3j = t3j + sgk/(a1a*a2a*a3a*a4a*a5a*a6a)
  enddo

  t3j = C*t3j

142 continue
end subroutine t3j_all
subroutine fctor1(nfac,fctor2)
  ! ----------------------------------------------------------------
  !     implicit none
  double precision :: f,f2(0:100),fctor2
  integer :: nfac,i
  f2(0) = 1d0
  f = 0d0
  do  i = 1,nfac
     f = f+1d0
     f2(i) = f2(i-1)*f
  enddo
  fctor2 = f2(nfac)
end subroutine fctor1

subroutine tri_rule(j1,j2,j3,tri_ru)
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: j1,j2,j3,tri_ru
  !     check the triangular rule
  tri_ru = 1
  if ((j1+j2-j3) < 0)tri_ru = -1
  if ((j1-j2+j3) < 0)tri_ru = -1
  if ((-j1+j2+j3) < 0)tri_ru = -1
end subroutine tri_rule

