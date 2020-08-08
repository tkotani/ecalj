      subroutine Veecomp(Vee,l,U_H,J_H)
C- Calculate Vee from U, J, Slater integrals, and Clebsch-Gordon coefficients.
C ----------------------------------------------------------------------
Ci Inputs
Ci   l     :l for which U_H and J_H are defined
Ci   U_H   :Screened direct coulomb integral (Hubbard U)
Ci   J_H   :Exchange integral (Hubbard J)
Co Outputs
Co   Vee   :Vee(m,m2,m3,m4) = <m,m2|Vee|m1,m3>
Cl Local variables
Cr Remarks
Cu Updates
Cu   02 Jun 05 Lambrecht bug fixes
Cu   27 Apr 05 Larson first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer l
      double precision U_H,J_H
      double precision Vee(-3:3,-3:3,-3:3,-3:3)
C ... Local parameters
      integer m,m1,m2,m3,k,jj,sg
      double precision F(0:2*3),r42,r62,t1,t2,t3,ak
C     Slater integrals from U_H, J_H (this will be updated later)

      call dpzero(vee,7**4)
      call dpzero(F,7)

      F(0) = U_H
C     d-states -> r42 = F(4)/F(2) = 0.625
      if (l .eq. 2) then
        r42 = 0.625d0
        F(2) = 14d0*J_H/(1d0+r42)
        F(4) = F(2)*r42
      endif
C     f-states
      if (l .eq. 3) then
        r42 = 451d0/675d0
        r62 = 1001d0/2025d0
        F(2) = 3d0*J_H/(2d0/15d0+1d0/11d0*r42+50d0/429d0*r62)
        F(4) = F(2)*r42
        F(6) = F(2)*r62
      endif
C     Calculate <m,m2|Vee|m1,m3>
      do  m = -l, l
        do  m2 = -l, l
          do  m1 = -l, l
            do  m3 = -l, l
              do  k = 0, 2*l,2
c                     print *,'k=',k,'F=',F(k)
                ak = 0d0
                call t3j_all(l,k,l,0,0,0,t1)
c                     print *,'l=',l,'k=',k,'F=',F(k),'t1=',t1
                do  jj = -k, k
                  call t3j_all(l,k,l,-m,jj,m1,t2)
                  call t3j_all(l,k,l,-m2,-jj,m3,t3)
                  sg = (-1)**(m+jj+m2)
                  ak = ak + sg*t2*t3
                enddo
                ak = (2*l+1)**2*t1**2*ak
c                     print '(5i3,3f8.2)',
c     .               m,m2,m1,m3,k,Vee(m,m2,m1,m3),ak,F(k)
                Vee(m,m2,m1,m3) = Vee(m,m2,m1,m3) + ak*F(k)
c                     print '(5i3,3f9.4)',
c     .               m,m2,m1,m3,k,Vee(m,m2,m1,m3),ak,F(k)
              enddo
            enddo
          enddo
        enddo
      enddo
      end

      subroutine t3j_all(j1,j2,j3,m1,m2,m3,t3j)
C ----------------------------------------------------------------
C- Calculate 3-j symbols.
C     implicit none
      integer j1,j2,j3,m1,m2,m3
      integer J,L1,L2,K0,K1,K2,KMAX,KMIN,trr,k
      double precision sgJ,sgk,C,ab,t3j
      double precision a1,a2,a3,a4,a5,a6,a7,a8,a9,b
      double precision a1a,a2a,a3a,a4a,a5a,a6a,a7a,a8a,a9a
C
C
      t3j = 0d0
      if ((m1+m2+m3) .ne. 0) goto 142
      call tri_rule(j1,j2,j3,trr)
      if (trr .lt. 0) goto 142
C
C
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

 142  continue
      end
      subroutine fctor1(nfac,fctor2)
C ----------------------------------------------------------------
C     implicit none
      double precision f,f2(0:100),fctor2
      integer nfac,i
      f2(0) = 1d0
      f = 0d0
      do  i = 1,nfac
        f = f+1d0
        f2(i) = f2(i-1)*f
      enddo
      fctor2 = f2(nfac)
      end

      subroutine tri_rule(j1,j2,j3,tri_ru)
C ----------------------------------------------------------------
C     implicit none
      integer j1,j2,j3,tri_ru
C     check the triangular rule
      tri_ru = 1
      if ((j1+j2-j3).lt.0)tri_ru = -1
      if ((j1-j2+j3).lt.0)tri_ru = -1
      if ((-j1+j2+j3).lt.0)tri_ru = -1
      end

